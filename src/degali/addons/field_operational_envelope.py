"""Operational disposition across a complete deterministic field envelope.

A nominal field result with bounded inputs cannot stand in for its uncertainty
corners.  This module runs every declared direct-flash source/weather/primary-
sensor corner, optionally refines every one, and aggregates the existing
fail-safe decisions without labelling the deterministic sensitivity set a
probability interval.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Literal

from .field_decision import (
    FieldOperationalScreeningDecision,
    _aggregate_gate_codes,
    _validate_aggregate_operational_decision,
    _validate_case_operational_decision,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from .field_workflow import (
    FieldSemiFVEnvelope,
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    run_field_semi_fv_envelope,
    run_field_semi_fv_refinement_study,
)
from .field_source_io import FieldAtmosphericSourceSchedule


FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA = (
    "degali.field-operational-uncertainty-envelope.v1"
)

_OPERATIONAL_SENSOR_METRICS = (
    "peak_true_mole_fraction",
    "peak_indicated_mole_fraction",
    "final_true_mole_fraction",
    "final_indicated_mole_fraction",
    "time_average_true_mole_fraction",
    "time_average_indicated_mole_fraction",
)


def _trace_metric(trace: object, metric: str) -> float:
    """Evaluate one finite summary operator without resampling a trace."""
    if metric.startswith("peak_"):
        signal_name = metric.removeprefix("peak_")
    elif metric.startswith("final_"):
        signal_name = metric.removeprefix("final_")
    elif metric.startswith("time_average_"):
        signal_name = metric.removeprefix("time_average_")
    else:  # pragma: no cover - the constant tuple owns the vocabulary
        raise ValueError(f"unsupported operational sensor metric {metric!r}")
    values = tuple(float(value) for value in getattr(trace, signal_name))
    time = tuple(float(value) for value in getattr(trace, "time_s"))
    if not values or len(values) != len(time):
        raise ValueError("operational sensor trace arrays must be non-empty and aligned")
    if any(not math.isfinite(value) for value in values + time):
        raise ValueError("operational sensor trace arrays must be finite")
    if metric.startswith("peak_"):
        result = max(values)
    elif metric.startswith("final_"):
        result = values[-1]
    else:
        duration = time[-1] - time[0]
        if duration <= 0.0:
            result = values[-1]
        else:
            integral = math.fsum(
                0.5 * (left + right) * (right_time - left_time)
                for left, right, left_time, right_time in zip(
                    values[:-1], values[1:], time[:-1], time[1:]
                )
            )
            result = integral / duration
    if not math.isfinite(result):
        raise ValueError("operational sensor metric must be finite")
    return float(result)


def deterministic_sensor_envelope_record(
    cases: tuple[object, ...],
    *,
    request: FieldSemiFVRequest | None = None,
) -> dict[str, object]:
    """Summarise deterministic sensor extrema while preserving withheld cases.

    The returned bounds are extrema over completed deterministic corners, not a
    confidence interval.  A missing or withheld trace is counted explicitly;
    its value is never replaced with zero or an interpolated estimate.
    """
    if not isinstance(cases, tuple) or not cases:
        raise ValueError("sensor envelope requires at least one operational case")
    if request is not None and not isinstance(request, FieldSemiFVRequest):
        raise TypeError("sensor envelope request must be a FieldSemiFVRequest or None")
    expected_case_count = len(cases)
    buckets: dict[str, dict[str, object]] = {}

    def bucket_for(label: str, position: tuple[float, float, float]) -> dict[str, object]:
        if label not in buckets:
            buckets[label] = {
                "label": label,
                "position_m": [float(value) for value in position],
                "available_case_count": 0,
                "withheld_case_count": 0,
                "withheld_reasons": [],
                "values": {metric: [] for metric in _OPERATIONAL_SENSOR_METRICS},
            }
        return buckets[label]

    if request is not None:
        declared = {}
        if request.scenario.sensor is not None:
            declared["field_sensor"] = request.scenario.sensor
        declared.update({item.label: item.sensor for item in request.sensor_deployments})
        for label, sensor in declared.items():
            bucket_for(label, tuple(float(value) for value in sensor.position_m))

    for case in cases:
        screening = getattr(case, "screening", None)
        refinement = getattr(case, "refinement", None)
        if refinement is not None:
            screening = refinement.screening
        if screening is None:
            continue
        declared: dict[str, object] = {}
        if screening.request.scenario.sensor is not None:
            declared["field_sensor"] = screening.request.scenario.sensor
        declared.update({item.label: item.sensor for item in screening.request.sensor_deployments})
        results = {item.label: item for item in screening.sensor_results}
        legacy_primary_trace = (
            screening.sensor_trace
            if screening.request.scenario.sensor is not None
            else None
        )
        for label, sensor in declared.items():
            bucket = bucket_for(label, tuple(float(value) for value in sensor.position_m))
            item = results.get(label)
            trace = (
                legacy_primary_trace
                if item is None and label == "field_sensor"
                else None if item is None else item.trace
            )
            if trace is None:
                bucket["withheld_case_count"] = int(bucket["withheld_case_count"]) + 1
                reasons = bucket["withheld_reasons"]
                assert isinstance(reasons, list)
                reasons.extend(
                    list(item.warnings) if item is not None and item.warnings
                    else [f"sensor {label!r} has no completed trace in this corner"]
                )
                continue
            bucket["available_case_count"] = int(bucket["available_case_count"]) + 1
            values = bucket["values"]
            assert isinstance(values, dict)
            for metric in _OPERATIONAL_SENSOR_METRICS:
                values[metric].append(_trace_metric(trace, metric))

    sensor_records = []
    for label, bucket in buckets.items():
        values = bucket["values"]
        assert isinstance(values, dict)
        bounds = {
            metric: {
                "minimum": min(metric_values),
                "maximum": max(metric_values),
            }
            for metric, metric_values in values.items()
            if metric_values
        }
        available = int(bucket["available_case_count"])
        recorded_withheld = int(bucket["withheld_case_count"])
        withheld = expected_case_count - available
        if withheld > recorded_withheld:
            reasons_list = bucket["withheld_reasons"]
            assert isinstance(reasons_list, list)
            reasons_list.append(
                f"sensor {label!r} had no screening trace in "
                f"{withheld - recorded_withheld} corner(s)"
            )
        reasons = tuple(dict.fromkeys(str(reason) for reason in bucket["withheld_reasons"]))
        sensor_records.append({
            "label": label,
            "position_m": bucket["position_m"],
            "available_case_count": available,
            "withheld_case_count": withheld,
            "complete": available == expected_case_count and withheld == 0,
            "withheld_reasons": list(reasons),
            "bounds": bounds,
        })
    complete = bool(sensor_records) and all(item["complete"] for item in sensor_records)
    return {
        "status": "complete" if complete else "withheld",
        "case_count": expected_case_count,
        "sensor_count": len(sensor_records),
        "interpretation": (
            "deterministic extrema over completed uncertainty corners; "
            "not a probability interval or confidence interval"
        ),
        "sensors": sensor_records,
    }


def _validate_selection(selection: object, *, name: str) -> None:
    """Validate a JSON/report-safe deterministic corner selection."""
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
        elif not isinstance(value, str) or not value.strip():
            raise TypeError(
                f"{name} values must be finite numbers or non-empty strings"
            )


@dataclass(frozen=True)
class FieldOperationalUncertaintyEnvelopeCase:
    """One fully specified input corner and its refinement/decision evidence."""

    selection: tuple[tuple[str, float | str], ...]
    screening: FieldSemiFVScreeningResult
    refinement: FieldSemiFVRefinementResult | None
    decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        _validate_selection(self.selection, name="uncertainty-envelope selection")
        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError("uncertainty-envelope screening must be FieldSemiFVScreeningResult")
        if self.refinement is not None and not isinstance(
            self.refinement, FieldSemiFVRefinementResult
        ):
            raise TypeError("uncertainty-envelope refinement must be FieldSemiFVRefinementResult or None")
        if not isinstance(self.decision, FieldOperationalScreeningDecision):
            raise TypeError("uncertainty-envelope decision must be FieldOperationalScreeningDecision")
        if self.refinement is None and self.decision.screening_allowed:
            raise ValueError(
                "uncertainty-envelope decision cannot allow a missing refinement result"
            )
        _validate_case_operational_decision(
            self.decision, self.screening, self.refinement,
            context="uncertainty-envelope case",
        )


@dataclass(frozen=True)
class FieldOperationalUncertaintyEnvelope:
    """A complete deterministic field envelope with one aggregate decision."""

    envelope: FieldSemiFVEnvelope
    cases: tuple[FieldOperationalUncertaintyEnvelopeCase, ...]
    operational_decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        if not isinstance(self.envelope, FieldSemiFVEnvelope):
            raise TypeError("envelope must be FieldSemiFVEnvelope")
        if not isinstance(self.cases, tuple):
            raise TypeError("operational uncertainty envelope cases must be a tuple")
        if not self.cases:
            raise ValueError("operational uncertainty envelope needs at least one case")
        if len(self.cases) != len(self.envelope.cases):
            raise ValueError("operational uncertainty envelope must retain every field corner")
        if any(
            not isinstance(item, FieldOperationalUncertaintyEnvelopeCase)
            for item in self.cases
        ):
            raise TypeError(
                "operational uncertainty envelope cases must contain only "
                "FieldOperationalUncertaintyEnvelopeCase values"
            )
        if not isinstance(self.operational_decision, FieldOperationalScreeningDecision):
            raise TypeError("operational_decision must be FieldOperationalScreeningDecision")
        expected = {tuple(sorted(case.values)) for case in self.envelope.cases}
        actual = {tuple(sorted(case.selection)) for case in self.cases}
        if actual != expected:
            raise ValueError("operational uncertainty envelope selections do not cover every declared corner")
        _validate_aggregate_operational_decision(
            self.operational_decision, self.cases,
            context="operational uncertainty envelope",
        )


def _aggregate_decision(
    cases: tuple[FieldOperationalUncertaintyEnvelopeCase, ...],
    *,
    allow_conditional: bool,
) -> FieldOperationalScreeningDecision:
    """Allow screening only if every deterministic corner passed its own gate."""
    decisions = tuple(case.decision for case in cases)
    withheld = tuple(case for case in cases if not case.decision.screening_allowed)
    refinement_available = all(item.refinement_available for item in decisions)
    if withheld:
        reasons = []
        actions = []
        for case in withheld:
            label = repr(dict(case.selection))
            reasons.extend(f"corner {label}: {reason}" for reason in case.decision.reasons)
            actions.extend(f"corner {label}: {action}" for action in case.decision.required_actions)
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(actions)),
            refinement_required=True, refinement_available=refinement_available,
            uncertainty_resolution_required=True,
            uncertainty_resolved=all(
                item.decision.uncertainty_resolved for item in cases
            ),
            gate_codes=_aggregate_gate_codes(decisions),
        )
    conditional = any(item.status == "conditional_allowed" for item in decisions)
    if conditional and not allow_conditional:
        # This should normally have been withheld by each individual decision,
        # but retaining the check keeps the aggregate gate self-contained.
        return FieldOperationalScreeningDecision(
            "withheld", False, False,
            ("one or more field uncertainty corners remain conditional without accountable review opt-in",),
            ("set allow_conditional=True only after reviewing every corner's recorded scope limits",),
            refinement_required=True, refinement_available=refinement_available,
            uncertainty_resolution_required=True, uncertainty_resolved=True,
            gate_codes=_aggregate_gate_codes(
                decisions, "conditional_review_required",
            ),
        )
    status: Literal["screening_allowed", "conditional_allowed"] = (
        "conditional_allowed" if conditional else "screening_allowed"
    )
    reasons = []
    actions = []
    for case in cases:
        label = repr(dict(case.selection))
        reasons.extend(f"corner {label}: {reason}" for reason in case.decision.reasons)
        actions.extend(f"corner {label}: {action}" for action in case.decision.required_actions)
    return FieldOperationalScreeningDecision(
        status, True, False,
        tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(actions)),
        refinement_required=True, refinement_available=refinement_available,
        uncertainty_resolution_required=True, uncertainty_resolved=True,
        gate_codes=_aggregate_gate_codes(decisions),
    )


def run_field_operational_uncertainty_envelope(
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
) -> FieldOperationalUncertaintyEnvelope:
    """Run and aggregate every declared field-input uncertainty corner.

    The underlying field envelope refuses direct-vapour schedules and
    supplemental detector uncertainty because those require their dedicated
    time-aligned or sensor-array envelopes.  An already-atmospheric source
    schedule may be supplied through ``atmospheric_source_schedule``; its
    nominal/lower/upper deterministic corners are then refined and decided
    exactly like the weather and surface corners.  When refinement is
    disabled, the diagnostic corner reports are still produced but the
    aggregate decision remains withheld by the existing refinement requirement.
    """
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be FieldSemiFVRequest")
    if not isinstance(include_refinement, bool) or not isinstance(allow_conditional, bool):
        raise TypeError("include_refinement and allow_conditional must be boolean")
    envelope = run_field_semi_fv_envelope(
        request,
        max_cases=max_cases,
        table_nodes=table_nodes,
        atmospheric_source_schedule=atmospheric_source_schedule,
    )
    cases = []
    for field_case in envelope.cases:
        resolved_request = replace(
            field_case.result.request,
            stability_uncertainty_resolved=(
                field_case.result.request.stability_alternatives is not None
            ),
        )
        screening = replace(field_case.result, request=resolved_request)
        refinement = (
            run_field_semi_fv_refinement_study(
                resolved_request,
                refinement_factors=refinement_factors,
                relative_tolerance=relative_tolerance,
                max_cell_steps=max_cell_steps,
            )
            if include_refinement else None
        )
        decision = evaluate_field_operational_screening(
            screening if refinement is None else refinement,
            require_refinement=True,
            require_in_plane_sensor=True,
            require_resolved_uncertainty=True,
            allow_conditional=allow_conditional,
        )
        cases.append(FieldOperationalUncertaintyEnvelopeCase(
            selection=field_case.values,
            screening=screening,
            refinement=refinement,
            decision=decision,
        ))
    complete_cases = tuple(cases)
    return FieldOperationalUncertaintyEnvelope(
        envelope=envelope,
        cases=complete_cases,
        operational_decision=_aggregate_decision(
            complete_cases, allow_conditional=allow_conditional,
        ),
    )


def field_operational_uncertainty_envelope_report(
    result: FieldOperationalUncertaintyEnvelope,
) -> dict[str, object]:
    """Return the aggregate decision with every corner's exact evidence."""
    if not isinstance(result, FieldOperationalUncertaintyEnvelope):
        raise TypeError("result must be FieldOperationalUncertaintyEnvelope")
    from .field_report import (
        field_refinement_report,
        field_semi_fv_envelope_report,
        field_screening_report,
    )

    return {
        "schema": FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
        "field_uncertainty_envelope": field_semi_fv_envelope_report(result.envelope),
        "operational_screening": field_operational_screening_decision_record(
            result.operational_decision
        ),
        "deterministic_sensor_envelope": deterministic_sensor_envelope_record(
            result.cases, request=result.envelope.request,
        ),
        "cases": [
            {
                "selection": dict(case.selection),
                "field_result": (
                    field_screening_report(case.screening)
                    if case.refinement is None else field_refinement_report(case.refinement)
                ),
                "operational_screening": field_operational_screening_decision_record(
                    case.decision
                ),
            }
            for case in result.cases
        ],
    }


__all__ = [
    "FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA",
    "FieldOperationalUncertaintyEnvelopeCase", "FieldOperationalUncertaintyEnvelope",
    "run_field_operational_uncertainty_envelope",
    "field_operational_uncertainty_envelope_report",
    "deterministic_sensor_envelope_record",
]
