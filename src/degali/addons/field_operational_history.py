"""Operational aggregate for the joint measured-history field envelope.

This is the source-history counterpart to the nominal field uncertainty
envelope. It preserves the existing global lower/upper historian selection,
weather and detector corner construction, then requires numerical refinement
and the ordinary fail-safe decision for every resulting case.
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
from .field_history import (
    FieldJointMeasuredHistoryScreeningEnvelope,
    MeasuredHistoryQualityCriteria,
    MeasuredReleaseHistory,
    PressureDrivenMeasuredHistory,
    run_field_joint_pressure_driven_history_envelope,
    run_field_joint_measured_history_envelope,
)
from .field_workflow import (
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    run_field_semi_fv_refinement_study,
)


FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA = (
    "degali.field-operational-measured-history-envelope.v1"
)


def _validate_source_selection(selection: object) -> None:
    if not isinstance(selection, tuple):
        raise TypeError("joint-history source_selection must be a tuple of key/value pairs")
    seen: set[str] = set()
    for item in selection:
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError("joint-history source_selection entries must be 2-item tuples")
        key, value = item
        if not isinstance(key, str) or not key.strip():
            raise ValueError("joint-history source_selection keys must be non-empty strings")
        if key in seen:
            raise ValueError("joint-history source_selection keys must be unique")
        seen.add(key)
        if not isinstance(value, str) or not value.strip():
            raise TypeError(
                "joint-history source_selection values must be non-empty strings"
            )


def _validate_field_selection(selection: object) -> None:
    if not isinstance(selection, tuple):
        raise TypeError("joint-history field_selection must be a tuple of key/value pairs")
    seen: set[str] = set()
    for item in selection:
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError("joint-history field_selection entries must be 2-item tuples")
        key, value = item
        if not isinstance(key, str) or not key.strip():
            raise ValueError("joint-history field_selection keys must be non-empty strings")
        if key in seen:
            raise ValueError("joint-history field_selection keys must be unique")
        seen.add(key)
        if isinstance(value, bool):
            raise TypeError("joint-history field_selection values cannot be boolean")
        if isinstance(value, (int, float)):
            if not math.isfinite(float(value)):
                raise ValueError("joint-history field_selection numeric values must be finite")
        elif not isinstance(value, str) or not value.strip():
            raise TypeError(
                "joint-history field_selection values must be finite numbers or non-empty strings"
            )


@dataclass(frozen=True)
class FieldOperationalMeasuredHistoryEnvelopeCase:
    """One source-history/weather/detector corner and its decision evidence."""

    source_selection: tuple[tuple[str, str], ...]
    field_selection: tuple[tuple[str, float | str], ...]
    screening: FieldSemiFVScreeningResult
    refinement: FieldSemiFVRefinementResult | None
    decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        _validate_source_selection(self.source_selection)
        _validate_field_selection(self.field_selection)
        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError("joint-history operational case requires FieldSemiFVScreeningResult")
        if self.refinement is not None and not isinstance(
            self.refinement, FieldSemiFVRefinementResult
        ):
            raise TypeError("joint-history refinement must be FieldSemiFVRefinementResult or None")
        if not isinstance(self.decision, FieldOperationalScreeningDecision):
            raise TypeError("joint-history decision must be FieldOperationalScreeningDecision")
        if self.refinement is None and self.decision.screening_allowed:
            raise ValueError(
                "joint-history operational decision cannot allow a missing refinement result"
            )
        _validate_case_operational_decision(
            self.decision, self.screening, self.refinement,
            context="joint-history operational case",
        )

    @property
    def label(self) -> str:
        return f"source={dict(self.source_selection)!r}, field={dict(self.field_selection)!r}"


@dataclass(frozen=True)
class FieldOperationalMeasuredHistoryEnvelope:
    """Complete joint source-history envelope and aggregate screening gate."""

    envelope: FieldJointMeasuredHistoryScreeningEnvelope
    cases: tuple[FieldOperationalMeasuredHistoryEnvelopeCase, ...]
    operational_decision: FieldOperationalScreeningDecision

    def __post_init__(self) -> None:
        if not isinstance(self.envelope, FieldJointMeasuredHistoryScreeningEnvelope):
            raise TypeError("envelope must be FieldJointMeasuredHistoryScreeningEnvelope")
        if not isinstance(self.cases, tuple):
            raise TypeError("operational measured-history envelope cases must be a tuple")
        if len(self.cases) != len(self.envelope.cases) or not self.cases:
            raise ValueError("operational measured-history envelope must retain every joint corner")
        if any(
            not isinstance(item, FieldOperationalMeasuredHistoryEnvelopeCase)
            for item in self.cases
        ):
            raise TypeError(
                "operational measured-history envelope cases must contain only "
                "FieldOperationalMeasuredHistoryEnvelopeCase values"
            )
        if not isinstance(self.operational_decision, FieldOperationalScreeningDecision):
            raise TypeError("operational_decision must be FieldOperationalScreeningDecision")
        expected = {
            (case.source_selection, case.field_values)
            for case in self.envelope.cases
        }
        actual = {
            (case.source_selection, case.field_selection)
            for case in self.cases
        }
        if actual != expected:
            raise ValueError("operational measured-history envelope selections do not cover every joint corner")
        _validate_aggregate_operational_decision(
            self.operational_decision, self.cases,
            context="operational measured-history envelope",
        )


def _aggregate(
    cases: tuple[FieldOperationalMeasuredHistoryEnvelopeCase, ...],
    *,
    allow_conditional: bool,
) -> FieldOperationalScreeningDecision:
    """Allow only when every history/field corner passed its own gate."""
    decisions = tuple(case.decision for case in cases)
    refinement_available = all(item.refinement_available for item in decisions)
    uncertainty_resolved = all(item.uncertainty_resolved for item in decisions)
    withheld = tuple(case for case in cases if not case.decision.screening_allowed)
    if withheld:
        reasons = []
        actions = []
        for case in withheld:
            reasons.extend(f"corner {case.label}: {item}" for item in case.decision.reasons)
            actions.extend(f"corner {case.label}: {item}" for item in case.decision.required_actions)
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
            ("one or more measured-history field corners remain conditional without accountable review opt-in",),
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
        reasons.extend(f"corner {case.label}: {item}" for item in case.decision.reasons)
        actions.extend(f"corner {case.label}: {item}" for item in case.decision.required_actions)
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


def _operationalize_joint_history_envelope(
    envelope: FieldJointMeasuredHistoryScreeningEnvelope,
    *,
    history_provenance: tuple[tuple[str, str], ...],
    include_refinement: bool,
    refinement_factors: tuple[int, ...],
    relative_tolerance: float,
    max_cell_steps: int,
    allow_conditional: bool,
) -> FieldOperationalMeasuredHistoryEnvelope:
    """Apply the common refinement and fail-safe decision to a joint envelope."""
    if not isinstance(envelope, FieldJointMeasuredHistoryScreeningEnvelope):
        raise TypeError("envelope must be FieldJointMeasuredHistoryScreeningEnvelope")
    provenance = tuple(history_provenance)
    if any(
        not isinstance(key, str) or not key.strip()
        or not isinstance(value, str) or not value.strip()
        for key, value in provenance
    ):
        raise ValueError("history_provenance must contain non-empty string key/value pairs")
    if len({key for key, _value in provenance}) != len(provenance):
        raise ValueError("history_provenance keys must be unique")
    # Keep the nested screening envelope's serialized requests consistent with
    # the refined reports.  This is especially important for pressure-driven
    # imported histories: the schedule carries derivation metadata, while the
    # execution must bind every nested report to the fingerprinted CSV record.
    envelope = replace(
        envelope,
        cases=tuple(
            replace(
                joint_case,
                result=replace(
                    joint_case.result,
                    request=replace(
                        joint_case.result.request,
                        measured_history_provenance=provenance,
                        measured_history_source_uncertainty_resolved=True,
                    ),
                ),
            )
            for joint_case in envelope.cases
        ),
    )
    cases = []
    for joint_case in envelope.cases:
        resolved_request = replace(
            joint_case.result.request,
            measured_history_provenance=provenance,
            measured_history_source_uncertainty_resolved=True,
            stability_uncertainty_resolved=(
                joint_case.result.request.stability_alternatives is not None
            ),
        )
        screening = replace(joint_case.result, request=resolved_request)
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
        cases.append(FieldOperationalMeasuredHistoryEnvelopeCase(
            source_selection=joint_case.source_selection,
            field_selection=joint_case.field_values,
            screening=screening,
            refinement=refinement,
            decision=decision,
        ))
    complete_cases = tuple(cases)
    return FieldOperationalMeasuredHistoryEnvelope(
        envelope=envelope,
        cases=complete_cases,
        operational_decision=_aggregate(complete_cases, allow_conditional=allow_conditional),
    )


def run_field_operational_joint_measured_history_envelope(
    request: FieldSemiFVRequest,
    history: MeasuredReleaseHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria,
    history_provenance: tuple[tuple[str, str], ...] = (),
    max_cases: int = 64,
    include_refinement: bool = True,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    max_cell_steps: int = 20_000_000,
    allow_conditional: bool = False,
) -> FieldOperationalMeasuredHistoryEnvelope:
    """Refine and decide every approved joint measured-history corner.

    ``request`` must be the pre-schedule field request. The helper itself
    constructs the time-aligned direct-vapour schedules, so a caller cannot
    pretend that one nominal historian selection has resolved the whole source
    uncertainty. Optional provenance is copied into every refined field report
    after the complete joint envelope has been constructed.
    """
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be FieldSemiFVRequest")
    if not isinstance(history, MeasuredReleaseHistory):
        raise TypeError("history must be MeasuredReleaseHistory")
    if not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria")
    if request.direct_vapour_schedule is not None:
        raise ValueError("joint measured-history operational envelope requires a pre-schedule request")
    if not isinstance(include_refinement, bool) or not isinstance(allow_conditional, bool):
        raise TypeError("include_refinement and allow_conditional must be boolean")
    envelope = run_field_joint_measured_history_envelope(
        request, history, quality_criteria=quality_criteria, max_cases=max_cases,
    )
    return _operationalize_joint_history_envelope(
        envelope,
        history_provenance=history_provenance,
        include_refinement=include_refinement,
        refinement_factors=refinement_factors,
        relative_tolerance=relative_tolerance,
        max_cell_steps=max_cell_steps,
        allow_conditional=allow_conditional,
    )


def run_field_operational_joint_pressure_driven_history_envelope(
    request: FieldSemiFVRequest,
    history: PressureDrivenMeasuredHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria,
    history_provenance: tuple[tuple[str, str], ...] = (),
    max_cases: int = 64,
    include_refinement: bool = True,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    max_cell_steps: int = 20_000_000,
    allow_conditional: bool = False,
) -> FieldOperationalMeasuredHistoryEnvelope:
    """Refine and gate every pressure-driven source/field corner."""
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be FieldSemiFVRequest")
    if not isinstance(history, PressureDrivenMeasuredHistory):
        raise TypeError("history must be PressureDrivenMeasuredHistory")
    if not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria")
    if request.direct_vapour_schedule is not None:
        raise ValueError(
            "joint pressure-driven-history operational envelope requires a pre-schedule request"
        )
    if not isinstance(include_refinement, bool) or not isinstance(allow_conditional, bool):
        raise TypeError("include_refinement and allow_conditional must be boolean")
    envelope = run_field_joint_pressure_driven_history_envelope(
        request, history, quality_criteria=quality_criteria, max_cases=max_cases,
    )
    return _operationalize_joint_history_envelope(
        envelope,
        history_provenance=history_provenance,
        include_refinement=include_refinement,
        refinement_factors=refinement_factors,
        relative_tolerance=relative_tolerance,
        max_cell_steps=max_cell_steps,
        allow_conditional=allow_conditional,
    )


def field_operational_joint_measured_history_envelope_report(
    result: FieldOperationalMeasuredHistoryEnvelope,
) -> dict[str, object]:
    """Serialise complete historian provenance, corners and the aggregate gate."""
    if not isinstance(result, FieldOperationalMeasuredHistoryEnvelope):
        raise TypeError("result must be FieldOperationalMeasuredHistoryEnvelope")
    from .field_report import (
        field_joint_measured_history_envelope_report,
        field_refinement_report,
        field_screening_report,
    )
    from .field_operational_envelope import deterministic_sensor_envelope_record

    return {
        "schema": FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA,
        "joint_measured_history_envelope": field_joint_measured_history_envelope_report(
            result.envelope
        ),
        "operational_screening": field_operational_screening_decision_record(
            result.operational_decision
        ),
        "deterministic_sensor_envelope": deterministic_sensor_envelope_record(
            result.cases, request=result.envelope.request,
        ),
        "cases": [
            {
                "source_selection": dict(case.source_selection),
                "field_selection": dict(case.field_selection),
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
    "FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA",
    "FieldOperationalMeasuredHistoryEnvelopeCase",
    "FieldOperationalMeasuredHistoryEnvelope",
    "run_field_operational_joint_measured_history_envelope",
    "run_field_operational_joint_pressure_driven_history_envelope",
    "field_operational_joint_measured_history_envelope_report",
]
