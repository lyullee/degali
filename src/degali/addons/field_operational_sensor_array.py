"""Operational aggregation for detector-calibration uncertainty corners.

The scalar field is solved once because response time, gain and bias are
observation operators.  Each calibration corner still receives its own
decision-facing screening record, while unresolved source, weather, geometry
or numerical uncertainty remains withheld by the ordinary operational gate.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal

from .field_decision import (
    FieldOperationalScreeningDecision,
    _aggregate_gate_codes,
    _validate_aggregate_operational_decision,
    _validate_case_operational_decision,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from .field_report import (
    field_refinement_report,
    field_screening_report,
    field_sensor_array_uncertainty_envelope_report,
)
from .field_workflow import (
    FieldSensorArrayUncertaintyCase,
    FieldSensorArrayUncertaintyEnvelope,
    FieldSensorDeployment,
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    run_field_sensor_array_uncertainty_envelope,
    run_field_semi_fv_refinement_study,
)


FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA = (
    "degali.field-operational-sensor-array-envelope.v1"
)


@dataclass(frozen=True)
class FieldOperationalSensorArrayEnvelopeCase:
    """One detector-calibration corner and its operational evidence."""

    calibration_case: FieldSensorArrayUncertaintyCase
    screening: FieldSemiFVScreeningResult
    refinement: FieldSemiFVRefinementResult | None
    decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        if not isinstance(self.calibration_case, FieldSensorArrayUncertaintyCase):
            raise TypeError(
                "operational sensor-array calibration_case must be a "
                "FieldSensorArrayUncertaintyCase"
            )
        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError(
                "operational sensor-array screening must be a "
                "FieldSemiFVScreeningResult"
            )
        if self.refinement is not None and not isinstance(
            self.refinement, FieldSemiFVRefinementResult
        ):
            raise TypeError(
                "operational sensor-array refinement must be a "
                "FieldSemiFVRefinementResult or None"
            )
        if not isinstance(self.decision, FieldOperationalScreeningDecision):
            raise TypeError(
                "operational sensor-array decision must be a "
                "FieldOperationalScreeningDecision"
            )
        if self.refinement is None and self.decision.screening_allowed:
            raise ValueError(
                "operational sensor-array decision cannot allow a missing refinement result"
            )
        _validate_case_operational_decision(
            self.decision, self.screening, self.refinement,
            context="operational sensor-array case",
        )

    @property
    def selection(self) -> tuple[tuple[str, float], ...]:
        return self.calibration_case.values


@dataclass(frozen=True)
class FieldOperationalSensorArrayEnvelope:
    """Complete detector-calibration envelope with one aggregate decision."""

    sensor_envelope: FieldSensorArrayUncertaintyEnvelope
    cases: tuple[FieldOperationalSensorArrayEnvelopeCase, ...]
    operational_decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        if not isinstance(self.sensor_envelope, FieldSensorArrayUncertaintyEnvelope):
            raise TypeError(
                "sensor_envelope must be a FieldSensorArrayUncertaintyEnvelope"
            )
        if not isinstance(self.cases, tuple):
            raise TypeError("operational sensor-array cases must be a tuple")
        if any(
            not isinstance(item, FieldOperationalSensorArrayEnvelopeCase)
            for item in self.cases
        ):
            raise TypeError(
                "operational sensor-array cases must contain only "
                "FieldOperationalSensorArrayEnvelopeCase values"
            )
        if len(self.cases) != len(self.sensor_envelope.cases):
            raise ValueError(
                "operational sensor-array envelope must retain every calibration corner"
            )
        expected = {tuple(case.values) for case in self.sensor_envelope.cases}
        actual = {tuple(case.selection) for case in self.cases}
        if actual != expected:
            raise ValueError(
                "operational sensor-array selections do not cover every calibration corner"
            )
        if not isinstance(self.operational_decision, FieldOperationalScreeningDecision):
            raise TypeError(
                "operational sensor-array operational_decision must be a "
                "FieldOperationalScreeningDecision"
            )
        if not self.cases:
            if self.operational_decision.status != "withheld":
                raise ValueError(
                    "empty operational sensor-array envelope must be withheld"
                )
            if self.operational_decision.refinement_available or self.operational_decision.uncertainty_resolved:
                raise ValueError(
                    "empty operational sensor-array envelope cannot claim resolved evidence"
                )
        _validate_aggregate_operational_decision(
            self.operational_decision, self.cases,
            context="operational sensor-array envelope",
        )


def _screening_at_calibration_corner(
    request: FieldSemiFVRequest,
    base: FieldSemiFVScreeningResult,
    calibration_case: FieldSensorArrayUncertaintyCase,
) -> FieldSemiFVScreeningResult:
    """Attach one exact detector calibration to the shared transport audit."""
    by_label = {item.label: item for item in calibration_case.sensor_results}
    primary = request.scenario.sensor
    if primary is not None:
        primary_result = by_label.get("field_sensor")
        if primary_result is None:
            raise ValueError(
                "sensor-array calibration corner is missing the primary field_sensor result"
            )
        primary = primary_result.sensor
    deployments = tuple(
        FieldSensorDeployment(
            deployment.label,
            by_label[deployment.label].sensor,
        )
        for deployment in request.sensor_deployments
    )
    exact_request = replace(
        request,
        scenario=replace(request.scenario, sensor=primary),
        sensor_deployments=deployments,
    )
    primary_trace = None
    if request.scenario.sensor is not None:
        primary_trace = by_label["field_sensor"].trace
    return replace(
        base,
        request=exact_request,
        sensor_trace=primary_trace,
        sensor_results=calibration_case.sensor_results,
    )


def _aggregate(
    cases: tuple[FieldOperationalSensorArrayEnvelopeCase, ...],
    *,
    allow_conditional: bool,
) -> FieldOperationalScreeningDecision:
    decisions = tuple(case.decision for case in cases)
    if not decisions:
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            ("sensor calibration uncertainty envelope produced no executable corners",),
            ("resolve the base field screening failure before operational screening",),
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
            label = repr(dict(case.selection))
            reasons.extend(
                f"corner {label}: {reason}" for reason in case.decision.reasons
            )
            actions.extend(
                f"corner {label}: {action}"
                for action in case.decision.required_actions
            )
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
            ("one or more sensor calibration corners remain conditional without accountable review opt-in",),
            ("set allow_conditional=True only after reviewing every calibration corner's scope limits",),
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
        label = repr(dict(case.selection))
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


def run_field_operational_sensor_array_envelope(
    request: FieldSemiFVRequest,
    *,
    max_cases: int = 128,
    include_refinement: bool = True,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    max_cell_steps: int = 20_000_000,
    allow_conditional: bool = False,
) -> FieldOperationalSensorArrayEnvelope:
    """Refine and decide every declared detector-calibration corner.

    Calibration corners reuse the exact true scalar transport field. This
    helper resolves only detector calibration uncertainty; source, weather,
    obstacle, stability and transport uncertainty remain explicit gates.
    """
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(include_refinement, bool) or not isinstance(allow_conditional, bool):
        raise TypeError("include_refinement and allow_conditional must be boolean")
    sensor_envelope = run_field_sensor_array_uncertainty_envelope(
        request, max_cases=max_cases,
    )
    if not sensor_envelope.cases:
        return FieldOperationalSensorArrayEnvelope(
            sensor_envelope=sensor_envelope,
            cases=(),
            operational_decision=_aggregate((), allow_conditional=allow_conditional),
        )
    base_refinement = (
        run_field_semi_fv_refinement_study(
            request,
            refinement_factors=refinement_factors,
            relative_tolerance=relative_tolerance,
            max_cell_steps=max_cell_steps,
        )
        if include_refinement else None
    )
    cases = []
    for calibration_case in sensor_envelope.cases:
        screening = _screening_at_calibration_corner(
            request, sensor_envelope.screening, calibration_case,
        )
        refinement = (
            None if base_refinement is None
            else replace(base_refinement, screening=screening)
        )
        decision = evaluate_field_operational_screening(
            screening if refinement is None else refinement,
            require_refinement=True,
            require_in_plane_sensor=True,
            require_resolved_uncertainty=True,
            allow_conditional=allow_conditional,
        )
        cases.append(FieldOperationalSensorArrayEnvelopeCase(
            calibration_case=calibration_case,
            screening=screening,
            refinement=refinement,
            decision=decision,
        ))
    complete_cases = tuple(cases)
    return FieldOperationalSensorArrayEnvelope(
        sensor_envelope=sensor_envelope,
        cases=complete_cases,
        operational_decision=_aggregate(
            complete_cases, allow_conditional=allow_conditional,
        ),
    )


def field_operational_sensor_array_envelope_report(
    result: FieldOperationalSensorArrayEnvelope,
) -> dict[str, object]:
    """Serialise calibration evidence and every operational corner decision."""
    if not isinstance(result, FieldOperationalSensorArrayEnvelope):
        raise TypeError("result must be a FieldOperationalSensorArrayEnvelope")
    from .field_operational_envelope import deterministic_sensor_envelope_record

    return {
        "schema": FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA,
        "sensor_calibration_envelope": field_sensor_array_uncertainty_envelope_report(
            result.sensor_envelope
        ),
        "operational_screening": field_operational_screening_decision_record(
            result.operational_decision
        ),
        "deterministic_sensor_envelope": deterministic_sensor_envelope_record(
            result.cases, request=result.sensor_envelope.screening.request,
        ),
        "cases": [
            {
                "selection": dict(case.selection),
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
    "FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA",
    "FieldOperationalSensorArrayEnvelopeCase",
    "FieldOperationalSensorArrayEnvelope",
    "run_field_operational_sensor_array_envelope",
    "field_operational_sensor_array_envelope_report",
]
