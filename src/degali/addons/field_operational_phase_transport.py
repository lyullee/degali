"""Fail-safe operational aggregation for phase-routed pool-vapour transport.

This layer is intentionally after the physical phase-to-field handoff.  It
does not use one nominal pool history to waive source, weather, surface,
stability, detector or numerical uncertainty: every physical corner and every
declared sensor calibration corner receives its own refinement/decision.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
import math
from typing import Literal

from .field_contracts import BoundedValue, FieldApplicability, SensorModel
from .field_decision import (
    FieldOperationalScreeningDecision,
    _aggregate_gate_codes,
    _validate_aggregate_operational_decision,
    _validate_case_operational_decision,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from .field_phase_routing import FieldPhaseRoutingConfig, FieldPhaseRoutingUncertainty
from .field_phase_transport import (
    FieldPhaseRoutingTransportEnvelope,
    FieldPhaseRoutingTransportEnvelopeCase,
    field_phase_routing_transport_envelope_report,
    run_field_phase_routing_transport_envelope,
)
from .field_pool_launch import PoolVapourLaunchBoundary
from .field_workflow import (
    FieldSensorDeployment,
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    run_field_semi_fv_refinement_study,
    run_field_semi_fv_screening,
)


FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA = (
    "degali.field-operational-phase-routing-transport-envelope.v1"
)


def _validate_selection(
    selection: object,
    *,
    name: str,
    allow_strings: bool,
) -> None:
    if not isinstance(selection, tuple):
        raise TypeError(f"{name} must be a tuple of key/value pairs")
    seen: set[str] = set()
    for item in selection:
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError(f"{name} entries must be 2-item tuples")
        key, value = item
        if not isinstance(key, str) or not key.strip():
            raise ValueError(f"{name} keys must be non-empty strings")
        if key in seen:
            raise ValueError(f"{name} keys must be unique")
        seen.add(key)
        if isinstance(value, bool):
            raise TypeError(f"{name} values cannot be boolean")
        if isinstance(value, (int, float)):
            if not math.isfinite(float(value)):
                raise ValueError(f"{name} numeric values must be finite")
        elif not allow_strings or not isinstance(value, str) or not value.strip():
            raise TypeError(
                f"{name} values must be finite numbers"
                + (" or non-empty strings" if allow_strings else "")
            )


def _declared_sensors(request: FieldSemiFVRequest) -> tuple[tuple[str, SensorModel], ...]:
    primary = () if request.scenario.sensor is None else (("field_sensor", request.scenario.sensor),)
    return primary + tuple((item.label, item.sensor) for item in request.sensor_deployments)


def _sensor_selections(request: FieldSemiFVRequest) -> tuple[tuple[tuple[str, float], ...], ...]:
    choices = {
        f"{label}.{field}": (
            (value.nominal,) if value.is_exact else value.corners()
        )
        for label, sensor in _declared_sensors(request)
        for field, value in sensor.uncertainty_fields().items()
    }
    names = tuple(choices)
    return tuple(
        tuple((name, float(value)) for name, value in zip(names, combination))
        for combination in product(*(choices[name] for name in names))
    )


def _sensor_at_selection(
    sensor: SensorModel,
    label: str,
    values: dict[str, float],
) -> SensorModel:
    fields = sensor.uncertainty_fields()

    def exact(name: str) -> BoundedValue:
        original = fields[name]
        return BoundedValue(
            values[f"{label}.{name}"], unit=original.unit,
            source=f"{original.source}; phase-transport sensor corner",
        )

    return replace(
        sensor,
        response_time_s=exact("sensor_response_time_s"),
        gain=exact("sensor_gain"),
        bias_mole_fraction=exact("sensor_bias_mole_fraction"),
    )


def _calibrated_phase_request(
    physical_request: FieldSemiFVRequest,
    original_request: FieldSemiFVRequest,
    selection: tuple[tuple[str, float], ...],
) -> FieldSemiFVRequest:
    """Apply only detector calibration values to an otherwise fixed phase corner."""
    values = dict(selection)
    original_sensors = dict(_declared_sensors(original_request))
    primary = physical_request.scenario.sensor
    if primary is not None:
        primary = _sensor_at_selection(original_sensors["field_sensor"], "field_sensor", values)
    deployments = tuple(
        FieldSensorDeployment(
            deployment.label,
            _sensor_at_selection(original_sensors[deployment.label], deployment.label, values),
        )
        for deployment in physical_request.sensor_deployments
    )
    return replace(
        physical_request,
        scenario=replace(physical_request.scenario, sensor=primary),
        sensor_deployments=deployments,
        stability_uncertainty_resolved=(physical_request.stability_alternatives is not None),
        phase_routing_uncertainty_resolved=any(
            source.source_kind == "pool_vapour"
            for source in physical_request.distributed_vapour_sources
        ),
    )


def _with_phase_warnings(
    screening: FieldSemiFVScreeningResult,
    physical_case: FieldPhaseRoutingTransportEnvelopeCase,
) -> FieldSemiFVScreeningResult:
    """Keep the phase/pool scope in the decision-facing local result."""
    applicability = screening.applicability
    warnings = tuple(dict.fromkeys(applicability.warnings + physical_case.warnings))
    if applicability.status == "blocked":
        merged = FieldApplicability(
            "blocked", applicability.reasons, warnings, False,
        )
    else:
        merged = FieldApplicability(
            "conditional", (), warnings, applicability.uncertainty_complete,
        )
    return replace(screening, applicability=merged)


def _withheld_decision(reason: str) -> FieldOperationalScreeningDecision:
    reason_lower = reason.lower()
    detail_code = (
        "phase_routing_incomplete"
        if "phase-routing ledger" in reason_lower
        else "pool_transport_withheld"
        if "pool-vapour field transport withheld" in reason_lower
        else None
    )
    gate_codes = ("transport_incomplete",) if detail_code is None else (
        "transport_incomplete", detail_code,
    )
    return FieldOperationalScreeningDecision(
        "withheld", False, False, (reason,), (),
        refinement_required=True, refinement_available=False,
        uncertainty_resolution_required=True, uncertainty_resolved=False,
        gate_codes=gate_codes,
    )


@dataclass(frozen=True)
class FieldOperationalPhaseRoutingTransportEnvelopeCase:
    """One fully specified phase/pool and detector-calibration disposition."""

    phase_selection: tuple[tuple[str, float | str], ...]
    sensor_selection: tuple[tuple[str, float], ...]
    screening: FieldSemiFVScreeningResult | None
    refinement: FieldSemiFVRefinementResult | None
    decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        _validate_selection(
            self.phase_selection, name="phase_selection", allow_strings=True,
        )
        _validate_selection(
            self.sensor_selection, name="sensor_selection", allow_strings=False,
        )
        if self.screening is not None and not isinstance(
            self.screening, FieldSemiFVScreeningResult
        ):
            raise TypeError(
                "phase operational screening must be a FieldSemiFVScreeningResult or None"
            )
        if self.refinement is not None and not isinstance(
            self.refinement, FieldSemiFVRefinementResult
        ):
            raise TypeError(
                "phase operational refinement must be a FieldSemiFVRefinementResult or None"
            )
        if self.screening is None and self.refinement is not None:
            raise ValueError("phase operational refinement requires a screening result")
        if not isinstance(self.decision, FieldOperationalScreeningDecision):
            raise TypeError(
                "phase operational decision must be a FieldOperationalScreeningDecision"
            )
        if self.screening is None and self.decision.screening_allowed:
            raise ValueError(
                "phase operational decision cannot allow a missing screening result"
            )
        if self.refinement is None and self.decision.screening_allowed:
            raise ValueError(
                "phase operational decision cannot allow a missing refinement result"
            )
        _validate_case_operational_decision(
            self.decision, self.screening, self.refinement,
            context="phase operational case",
        )

    @property
    def label(self) -> str:
        return f"phase={dict(self.phase_selection)!r}, sensor={dict(self.sensor_selection)!r}"


@dataclass(frozen=True)
class FieldOperationalPhaseRoutingTransportEnvelope:
    """Complete phase/pool/sensor envelope with a fail-safe aggregate decision."""

    physical_envelope: FieldPhaseRoutingTransportEnvelope
    sensor_selections: tuple[tuple[tuple[str, float], ...], ...]
    cases: tuple[FieldOperationalPhaseRoutingTransportEnvelopeCase, ...]
    operational_decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        if not isinstance(self.physical_envelope, FieldPhaseRoutingTransportEnvelope):
            raise TypeError("physical_envelope must be a FieldPhaseRoutingTransportEnvelope")
        if not isinstance(self.sensor_selections, tuple):
            raise TypeError("phase operational sensor_selections must be a tuple")
        if not self.sensor_selections:
            raise ValueError("phase operational envelope requires at least one sensor selection")
        if not isinstance(self.cases, tuple):
            raise TypeError("phase operational cases must be a tuple")
        if not self.cases:
            raise ValueError("phase operational envelope requires at least one case")
        for selection in self.sensor_selections:
            _validate_selection(
                selection, name="sensor_selections", allow_strings=False,
            )
        if len(set(self.sensor_selections)) != len(self.sensor_selections):
            raise ValueError("phase operational sensor selections must be unique")
        if any(
            not isinstance(item, FieldOperationalPhaseRoutingTransportEnvelopeCase)
            for item in self.cases
        ):
            raise TypeError(
                "phase operational envelope cases must contain only "
                "FieldOperationalPhaseRoutingTransportEnvelopeCase values"
            )
        if not isinstance(self.operational_decision, FieldOperationalScreeningDecision):
            raise TypeError(
                "phase operational operational_decision must be a "
                "FieldOperationalScreeningDecision"
            )
        expected = {
            (physical_case.selection, sensor_selection)
            for physical_case in self.physical_envelope.cases
            for sensor_selection in self.sensor_selections
        }
        actual = {(case.phase_selection, case.sensor_selection) for case in self.cases}
        if actual != expected:
            raise ValueError("operational phase envelope must retain every phase and sensor corner")
        _validate_aggregate_operational_decision(
            self.operational_decision, self.cases,
            context="operational phase envelope",
        )


def _aggregate(
    cases: tuple[FieldOperationalPhaseRoutingTransportEnvelopeCase, ...],
    *,
    allow_conditional: bool,
) -> FieldOperationalScreeningDecision:
    withheld = tuple(case for case in cases if not case.decision.screening_allowed)
    refinement_available = all(case.decision.refinement_available for case in cases)
    uncertainty_resolved = all(case.decision.uncertainty_resolved for case in cases)
    if withheld:
        reasons = []
        actions = []
        for case in withheld:
            reasons.extend(f"corner {case.label}: {value}" for value in case.decision.reasons)
            actions.extend(f"corner {case.label}: {value}" for value in case.decision.required_actions)
        return FieldOperationalScreeningDecision(
            "withheld", False, False, tuple(dict.fromkeys(reasons)),
            tuple(dict.fromkeys(actions)), True, refinement_available,
            True, uncertainty_resolved,
            _aggregate_gate_codes(tuple(case.decision for case in cases)),
        )
    conditional = any(case.decision.status == "conditional_allowed" for case in cases)
    if conditional and not allow_conditional:
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            ("one or more phase-routing transport corners remain conditional without accountable review opt-in",),
            ("set allow_conditional=True only after reviewing every phase/pool/sensor corner",),
            True, refinement_available, True, uncertainty_resolved,
            _aggregate_gate_codes(
                tuple(case.decision for case in cases),
                "conditional_review_required",
            ),
        )
    reasons = []
    actions = []
    for case in cases:
        reasons.extend(f"corner {case.label}: {value}" for value in case.decision.reasons)
        actions.extend(f"corner {case.label}: {value}" for value in case.decision.required_actions)
    status: Literal["screening_allowed", "conditional_allowed"] = (
        "conditional_allowed" if conditional else "screening_allowed"
    )
    return FieldOperationalScreeningDecision(
        status, True, False, tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(actions)),
        True, refinement_available, True, uncertainty_resolved,
        _aggregate_gate_codes(tuple(case.decision for case in cases)),
    )


def run_field_operational_phase_routing_transport_envelope(
    request: FieldSemiFVRequest,
    phase_config: FieldPhaseRoutingConfig,
    pool_launch: PoolVapourLaunchBoundary,
    *,
    pool_vertical_sigma_m: float,
    pool_vertical_sigma_uncertainty: BoundedValue | None = None,
    phase_uncertainty: FieldPhaseRoutingUncertainty | None = None,
    max_cases: int = 128,
    table_nodes: int = 161,
    include_refinement: bool = True,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    max_cell_steps: int = 20_000_000,
    allow_conditional: bool = False,
) -> FieldOperationalPhaseRoutingTransportEnvelope:
    """Refine and decide every phase/pool physical and sensor corner."""
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(include_refinement, bool) or not isinstance(allow_conditional, bool):
        raise TypeError("include_refinement and allow_conditional must be boolean")
    if phase_uncertainty is not None and not isinstance(
        phase_uncertainty, FieldPhaseRoutingUncertainty
    ):
        raise TypeError("phase_uncertainty must be FieldPhaseRoutingUncertainty or None")
    sensor_selections = _sensor_selections(request)
    if len(sensor_selections) < 1:
        raise ValueError("phase operational envelope has no sensor selection")
    physical_envelope = run_field_phase_routing_transport_envelope(
        request, phase_config, pool_launch, pool_vertical_sigma_m=pool_vertical_sigma_m,
        pool_vertical_sigma_uncertainty=pool_vertical_sigma_uncertainty,
        phase_uncertainty=phase_uncertainty,
        max_cases=max_cases, table_nodes=table_nodes,
    )
    total_cases = len(physical_envelope.cases) * len(sensor_selections)
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} operational phase-routing transport corners exceed max_cases={max_cases}"
        )
    cases = []
    for physical_case in physical_envelope.cases:
        for sensor_selection in sensor_selections:
            if physical_case.screening is None:
                reason = "; ".join(physical_case.warnings) or "phase-to-field transport is unavailable"
                cases.append(FieldOperationalPhaseRoutingTransportEnvelopeCase(
                    physical_case.selection, sensor_selection, None, None,
                    _withheld_decision(reason),
                ))
                continue
            calibrated_request = _calibrated_phase_request(
                physical_case.screening.request, request, sensor_selection,
            )
            screening = _with_phase_warnings(
                run_field_semi_fv_screening(calibrated_request), physical_case,
            )
            refinement = (
                run_field_semi_fv_refinement_study(
                    calibrated_request,
                    refinement_factors=refinement_factors,
                    relative_tolerance=relative_tolerance,
                    max_cell_steps=max_cell_steps,
                )
                if include_refinement and screening.completed else None
            )
            decision = evaluate_field_operational_screening(
                screening if refinement is None else refinement,
                require_refinement=True,
                require_in_plane_sensor=True,
                require_resolved_uncertainty=True,
                allow_conditional=allow_conditional,
            )
            cases.append(FieldOperationalPhaseRoutingTransportEnvelopeCase(
                physical_case.selection, sensor_selection, screening, refinement, decision,
            ))
    complete_cases = tuple(cases)
    return FieldOperationalPhaseRoutingTransportEnvelope(
        physical_envelope, sensor_selections, complete_cases,
        _aggregate(complete_cases, allow_conditional=allow_conditional),
    )


def field_operational_phase_routing_transport_envelope_report(
    result: FieldOperationalPhaseRoutingTransportEnvelope,
) -> dict[str, object]:
    """Serialise all phase/pool/sensor decisions and the aggregate holdback."""
    if not isinstance(result, FieldOperationalPhaseRoutingTransportEnvelope):
        raise TypeError("result must be a FieldOperationalPhaseRoutingTransportEnvelope")
    from .field_report import field_refinement_report, field_screening_report
    from .field_operational_envelope import deterministic_sensor_envelope_record

    return {
        "schema": FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA,
        "physical_phase_transport_envelope": field_phase_routing_transport_envelope_report(
            result.physical_envelope
        ),
        "operational_screening": field_operational_screening_decision_record(
            result.operational_decision
        ),
        "deterministic_sensor_envelope": deterministic_sensor_envelope_record(
            result.cases, request=result.physical_envelope.request,
        ),
        "cases": [
            {
                "phase_selection": dict(case.phase_selection),
                "sensor_selection": dict(case.sensor_selection),
                "field_result": (
                    None if case.screening is None else (
                        field_screening_report(case.screening)
                        if case.refinement is None else field_refinement_report(case.refinement)
                    )
                ),
                "operational_screening": field_operational_screening_decision_record(
                    case.decision
                ),
            }
            for case in result.cases
        ],
    }


__all__ = [
    "FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA",
    "FieldOperationalPhaseRoutingTransportEnvelopeCase",
    "FieldOperationalPhaseRoutingTransportEnvelope",
    "run_field_operational_phase_routing_transport_envelope",
    "field_operational_phase_routing_transport_envelope_report",
]
