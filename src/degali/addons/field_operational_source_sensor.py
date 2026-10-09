"""Joint operational envelope for field inputs and supplemental sensors.

The physical source/weather/geometry corners are solved once each with
supplemental detector calibration fixed at nominal values.  Every detector
calibration corner is then applied to the unchanged true receptor traces.
This closes the Cartesian uncertainty path without treating detector
calibration as transport physics or waiving any physical uncertainty gate.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
import math
from typing import Literal

from .field_contracts import BoundedValue
from .field_decision import (
    FieldOperationalScreeningDecision,
    _aggregate_gate_codes,
    _validate_aggregate_operational_decision,
    _validate_case_operational_decision,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from .field_observation import apply_sensor_model
from .field_report import (
    field_refinement_report,
    field_screening_report,
    field_semi_fv_envelope_report,
)
from .field_source_io import FieldAtmosphericSourceSchedule
from .field_workflow import (
    FieldSensorArrayUncertaintyCase,
    FieldSensorDeployment,
    FieldSensorDeploymentResult,
    FieldSemiFVEnvelope,
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    run_field_semi_fv_envelope,
    run_field_semi_fv_refinement_study,
)


FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA = (
    "degali.field-operational-source-sensor-envelope.v1"
)


def _validate_selection(
    selection: object,
    *,
    name: str,
    numeric_only: bool,
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
        elif numeric_only or not isinstance(value, str) or not value.strip():
            raise TypeError(
                f"{name} values must be finite numbers"
                + (" or non-empty strings" if not numeric_only else "")
            )


def _sensor_at_corner(
    sensor,
    label: str,
    values: dict[str, float],
):
    fields = sensor.uncertainty_fields()

    def exact(name: str) -> BoundedValue:
        original = fields[name]
        return BoundedValue(
            values[f"{label}.{name}"],
            unit=original.unit,
            source=f"{original.source}; joint source-sensor corner",
        )

    return replace(
        sensor,
        response_time_s=exact("sensor_response_time_s"),
        gain=exact("sensor_gain"),
        bias_mole_fraction=exact("sensor_bias_mole_fraction"),
    )


def _nominal_sensor(sensor):
    fields = sensor.uncertainty_fields()

    def exact(name: str) -> BoundedValue:
        original = fields[name]
        return BoundedValue(
            original.nominal,
            unit=original.unit,
            source=f"{original.source}; joint source-sensor nominal transport",
        )

    return replace(
        sensor,
        response_time_s=exact("sensor_response_time_s"),
        gain=exact("sensor_gain"),
        bias_mole_fraction=exact("sensor_bias_mole_fraction"),
    )


def _sensor_selections(
    request: FieldSemiFVRequest,
    *,
    max_cases: int,
) -> tuple[tuple[tuple[str, float], ...], ...]:
    fields = {
        f"{deployment.label}.{name}": value
        for deployment in request.sensor_deployments
        for name, value in deployment.sensor.uncertainty_fields().items()
    }
    if not fields:
        return ((),)
    choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in fields.items()
    }
    count = math.prod(len(values) for values in choices.values())
    if count * 1 > max_cases:
        raise ValueError(
            f"{count} supplemental sensor calibration corners exceed max_cases={max_cases}"
        )
    names = tuple(choices)
    return tuple(
        tuple((name, float(value)) for name, value in zip(names, combination))
        for combination in product(*(choices[name] for name in names))
    )


def _screening_at_sensor_corner(
    base: FieldSemiFVScreeningResult,
    original_request: FieldSemiFVRequest,
    sensor_case: FieldSensorArrayUncertaintyCase,
) -> FieldSemiFVScreeningResult:
    values = dict(sensor_case.values)
    exact_deployments = tuple(
        FieldSensorDeployment(
            deployment.label,
            _sensor_at_corner(deployment.sensor, deployment.label, values),
        )
        for deployment in original_request.sensor_deployments
    )
    exact_request = replace(base.request, sensor_deployments=exact_deployments)
    traces_by_label = (
        {} if base.transport is None
        else {trace.receptor.label: trace for trace in base.transport.receptor_traces}
    )
    existing = {item.label: item for item in base.sensor_results}
    sensor_results: list[FieldSensorDeploymentResult] = []
    primary_trace = base.sensor_trace
    if base.request.scenario.sensor is not None:
        primary = existing.get("field_sensor")
        if primary is None:
            sensor_results.append(FieldSensorDeploymentResult(
                "field_sensor", base.request.scenario.sensor, None,
                ("primary field sensor result was absent from the physical corner",),
            ))
            primary_trace = None
        else:
            sensor_results.append(primary)
            primary_trace = primary.trace
    for deployment in exact_deployments:
        true_trace = traces_by_label.get(deployment.label)
        warnings = existing.get(deployment.label)
        if true_trace is None:
            sensor_results.append(FieldSensorDeploymentResult(
                deployment.label, deployment.sensor, None,
                (() if warnings is None else warnings.warnings)
                or (f"sensor {deployment.label!r} receptor trace was unavailable",),
            ))
            continue
        trace = apply_sensor_model(
            true_trace.time_s,
            true_trace.concentration_kg_m3,
            deployment.sensor,
            ambient_air_density_kg_m3=base.request.ambient_air_density_kg_m3,
        )
        sensor_results.append(FieldSensorDeploymentResult(
            deployment.label, deployment.sensor, trace,
            () if warnings is None else warnings.warnings,
        ))
    return replace(
        base,
        request=exact_request,
        sensor_trace=primary_trace,
        sensor_results=tuple(sensor_results),
    )


@dataclass(frozen=True)
class FieldOperationalSourceSensorEnvelopeCase:
    """One physical field corner crossed with one sensor-calibration corner."""

    physical_selection: tuple[tuple[str, float | str], ...]
    sensor_case: FieldSensorArrayUncertaintyCase
    screening: FieldSemiFVScreeningResult
    refinement: FieldSemiFVRefinementResult | None
    decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        _validate_selection(
            self.physical_selection,
            name="joint source-sensor physical_selection",
            numeric_only=False,
        )
        if not isinstance(self.sensor_case, FieldSensorArrayUncertaintyCase):
            raise TypeError(
                "joint source-sensor sensor_case must be a "
                "FieldSensorArrayUncertaintyCase"
            )
        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError("joint source-sensor screening must be a FieldSemiFVScreeningResult")
        if self.refinement is not None and not isinstance(
            self.refinement, FieldSemiFVRefinementResult
        ):
            raise TypeError(
                "joint source-sensor refinement must be a "
                "FieldSemiFVRefinementResult or None"
            )
        if not isinstance(self.decision, FieldOperationalScreeningDecision):
            raise TypeError(
                "joint source-sensor decision must be a "
                "FieldOperationalScreeningDecision"
            )
        if self.refinement is None and self.decision.screening_allowed:
            raise ValueError(
                "joint source-sensor decision cannot allow a missing refinement result"
            )
        _validate_case_operational_decision(
            self.decision, self.screening, self.refinement,
            context="joint source-sensor case",
        )

    @property
    def sensor_selection(self) -> tuple[tuple[str, float], ...]:
        return self.sensor_case.values


@dataclass(frozen=True)
class FieldOperationalSourceSensorEnvelope:
    """Complete physical-plus-sensor Cartesian envelope."""

    physical_envelope: FieldSemiFVEnvelope
    sensor_selections: tuple[tuple[tuple[str, float], ...], ...]
    cases: tuple[FieldOperationalSourceSensorEnvelopeCase, ...]
    operational_decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        if not isinstance(self.physical_envelope, FieldSemiFVEnvelope):
            raise TypeError("physical_envelope must be a FieldSemiFVEnvelope")
        if not isinstance(self.sensor_selections, tuple) or not self.sensor_selections:
            raise ValueError("joint source-sensor envelope needs sensor selections")
        if not isinstance(self.cases, tuple):
            raise TypeError("joint source-sensor cases must be a tuple")
        for selection in self.sensor_selections:
            _validate_selection(
                selection,
                name="joint source-sensor sensor_selection",
                numeric_only=True,
            )
        if len(set(self.sensor_selections)) != len(self.sensor_selections):
            raise ValueError("joint source-sensor sensor selections must be unique")
        if any(
            not isinstance(item, FieldOperationalSourceSensorEnvelopeCase)
            for item in self.cases
        ):
            raise TypeError(
                "joint source-sensor cases must contain only "
                "FieldOperationalSourceSensorEnvelopeCase values"
            )
        expected = {
            (physical.values, sensor_selection)
            for physical in self.physical_envelope.cases
            for sensor_selection in self.sensor_selections
        }
        actual = {
            (case.physical_selection, case.sensor_selection)
            for case in self.cases
        }
        if actual != expected:
            raise ValueError(
                "joint source-sensor envelope must retain every physical and sensor corner"
            )
        if not isinstance(self.operational_decision, FieldOperationalScreeningDecision):
            raise TypeError(
                "joint source-sensor operational_decision must be a "
                "FieldOperationalScreeningDecision"
            )
        _validate_aggregate_operational_decision(
            self.operational_decision, self.cases,
            context="joint source-sensor envelope",
        )


def _aggregate(
    cases: tuple[FieldOperationalSourceSensorEnvelopeCase, ...],
    *,
    allow_conditional: bool,
) -> FieldOperationalScreeningDecision:
    decisions = tuple(case.decision for case in cases)
    if not decisions:
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            ("joint source-sensor envelope produced no executable corners",),
            ("resolve the physical field envelope before operational screening",),
            refinement_required=True, refinement_available=False,
            uncertainty_resolution_required=True, uncertainty_resolved=False,
            gate_codes=("corner_withheld",),
        )
    refinement_available = all(item.refinement_available for item in decisions)
    uncertainty_resolved = all(item.uncertainty_resolved for item in decisions)
    withheld = tuple(case for case in cases if not case.decision.screening_allowed)
    if withheld:
        reasons = []
        actions = []
        for case in withheld:
            label = f"physical={dict(case.physical_selection)!r}, sensor={dict(case.sensor_selection)!r}"
            reasons.extend(f"corner {label}: {reason}" for reason in case.decision.reasons)
            actions.extend(f"corner {label}: {action}" for action in case.decision.required_actions)
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(actions)),
            refinement_required=True, refinement_available=refinement_available,
            uncertainty_resolution_required=True,
            uncertainty_resolved=uncertainty_resolved,
            gate_codes=_aggregate_gate_codes(decisions),
        )
    conditional = any(item.status == "conditional_allowed" for item in decisions)
    if conditional and not allow_conditional:
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            ("one or more joint source-sensor corners remain conditional without accountable review opt-in",),
            ("set allow_conditional=True only after reviewing every joint corner's scope limits",),
            refinement_required=True, refinement_available=refinement_available,
            uncertainty_resolution_required=True,
            uncertainty_resolved=uncertainty_resolved,
            gate_codes=_aggregate_gate_codes(
                decisions, "conditional_review_required",
            ),
        )
    reasons = []
    actions = []
    for case in cases:
        label = f"physical={dict(case.physical_selection)!r}, sensor={dict(case.sensor_selection)!r}"
        reasons.extend(f"corner {label}: {reason}" for reason in case.decision.reasons)
        actions.extend(f"corner {label}: {action}" for action in case.decision.required_actions)
    status: Literal["screening_allowed", "conditional_allowed"] = (
        "conditional_allowed" if conditional else "screening_allowed"
    )
    return FieldOperationalScreeningDecision(
        status, True, False,
        tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(actions)),
        refinement_required=True, refinement_available=refinement_available,
        uncertainty_resolution_required=True,
        uncertainty_resolved=uncertainty_resolved,
        gate_codes=_aggregate_gate_codes(decisions),
    )


def run_field_operational_source_sensor_envelope(
    request: FieldSemiFVRequest,
    *,
    atmospheric_source_schedule: FieldAtmosphericSourceSchedule | None = None,
    max_cases: int = 64,
    table_nodes: int = 161,
    include_refinement: bool = True,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    max_cell_steps: int = 20_000_000,
    allow_conditional: bool = False,
) -> FieldOperationalSourceSensorEnvelope:
    """Run every physical corner crossed with supplemental sensor corners."""
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(include_refinement, bool) or not isinstance(allow_conditional, bool):
        raise TypeError("include_refinement and allow_conditional must be boolean")
    if request.scenario.sensor is None and not request.sensor_deployments:
        raise ValueError(
            "joint source-sensor envelope requires at least one declared sensor"
        )
    sensor_selections = _sensor_selections(request, max_cases=max_cases)
    nominal_deployments = tuple(
        FieldSensorDeployment(deployment.label, _nominal_sensor(deployment.sensor))
        for deployment in request.sensor_deployments
    )
    physical_request = replace(request, sensor_deployments=nominal_deployments)
    physical_envelope = run_field_semi_fv_envelope(
        physical_request,
        max_cases=max_cases,
        table_nodes=table_nodes,
        atmospheric_source_schedule=atmospheric_source_schedule,
    )
    total = len(physical_envelope.cases) * len(sensor_selections)
    if total > max_cases:
        raise ValueError(
            f"{total} joint source-sensor uncertainty corners exceed max_cases={max_cases}"
        )
    cases = []
    for physical_case in physical_envelope.cases:
        resolved_request = replace(
            physical_case.result.request,
            stability_uncertainty_resolved=(
                physical_case.result.request.stability_alternatives is not None
            ),
        )
        base_screening = replace(physical_case.result, request=resolved_request)
        refinement = (
            run_field_semi_fv_refinement_study(
                resolved_request,
                refinement_factors=refinement_factors,
                relative_tolerance=relative_tolerance,
                max_cell_steps=max_cell_steps,
            )
            if include_refinement else None
        )
        for sensor_selection in sensor_selections:
            values = dict(sensor_selection)
            calibration_results = []
            traces_by_label = (
                {} if base_screening.transport is None else {
                    trace.receptor.label: trace
                    for trace in base_screening.transport.receptor_traces
                }
            )
            existing = {item.label: item for item in base_screening.sensor_results}
            for deployment in request.sensor_deployments:
                sensor = _sensor_at_corner(deployment.sensor, deployment.label, values)
                true_trace = traces_by_label.get(deployment.label)
                prior = existing.get(deployment.label)
                if true_trace is None:
                    calibration_results.append(FieldSensorDeploymentResult(
                        deployment.label, sensor, None,
                        (() if prior is None else prior.warnings)
                        or (f"sensor {deployment.label!r} receptor trace was unavailable",),
                    ))
                else:
                    calibration_results.append(FieldSensorDeploymentResult(
                        deployment.label,
                        sensor,
                        apply_sensor_model(
                            true_trace.time_s,
                            true_trace.concentration_kg_m3,
                            sensor,
                            ambient_air_density_kg_m3=base_screening.request.ambient_air_density_kg_m3,
                        ),
                        () if prior is None else prior.warnings,
                    ))
            primary = existing.get("field_sensor")
            if primary is not None:
                calibration_results.insert(0, primary)
            sensor_case = FieldSensorArrayUncertaintyCase(
                sensor_selection, tuple(calibration_results),
            )
            screening = _screening_at_sensor_corner(
                base_screening, request, sensor_case,
            )
            corner_refinement = (
                None if refinement is None else replace(refinement, screening=screening)
            )
            decision = evaluate_field_operational_screening(
                screening if corner_refinement is None else corner_refinement,
                require_refinement=True,
                require_in_plane_sensor=True,
                require_resolved_uncertainty=True,
                allow_conditional=allow_conditional,
            )
            cases.append(FieldOperationalSourceSensorEnvelopeCase(
                physical_selection=physical_case.values,
                sensor_case=sensor_case,
                screening=screening,
                refinement=corner_refinement,
                decision=decision,
            ))
    complete_cases = tuple(cases)
    return FieldOperationalSourceSensorEnvelope(
        physical_envelope=physical_envelope,
        sensor_selections=sensor_selections,
        cases=complete_cases,
        operational_decision=_aggregate(
            complete_cases, allow_conditional=allow_conditional,
        ),
    )


def field_operational_source_sensor_envelope_report(
    result: FieldOperationalSourceSensorEnvelope,
) -> dict[str, object]:
    """Serialise the physical envelope and every joint decision corner."""
    if not isinstance(result, FieldOperationalSourceSensorEnvelope):
        raise TypeError("result must be a FieldOperationalSourceSensorEnvelope")
    from .field_operational_envelope import deterministic_sensor_envelope_record

    return {
        "schema": FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA,
        "physical_field_envelope": field_semi_fv_envelope_report(
            result.physical_envelope
        ),
        "sensor_selections": [dict(selection) for selection in result.sensor_selections],
        "operational_screening": field_operational_screening_decision_record(
            result.operational_decision
        ),
        "deterministic_sensor_envelope": deterministic_sensor_envelope_record(
            result.cases, request=result.physical_envelope.request,
        ),
        "cases": [
            {
                "physical_selection": dict(case.physical_selection),
                "sensor_selection": dict(case.sensor_selection),
                "field_result": (
                    field_screening_report(case.screening)
                    if case.refinement is None
                    else field_refinement_report(case.refinement)
                ),
                "operational_screening": field_operational_screening_decision_record(
                    case.decision
                ),
            }
            for case in result.cases
        ],
    }


__all__ = [
    "FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA",
    "FieldOperationalSourceSensorEnvelopeCase",
    "FieldOperationalSourceSensorEnvelope",
    "run_field_operational_source_sensor_envelope",
    "field_operational_source_sensor_envelope_report",
]
