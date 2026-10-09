"""Fail-safe operational disposition for reduced-order field screenings.

The semi-FV solver's applicability record describes a calculation.  This
module adds the separate question of whether that calculation may be consumed
as an operational *screening* input.  It never grants a design-basis or
approval decision.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Literal

from .field_workflow import FieldSemiFVRefinementResult, FieldSemiFVScreeningResult


FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA = "degali.field-operational-screening-decision.v1"
FIELD_OPERATIONAL_GATE_CODES = frozenset({
    "physical_applicability_blocked",
    "physical_applicability_conditional",
    "transport_incomplete",
    "refinement_missing",
    "refinement_failed",
    "sensor_missing",
    "sensor_trace_missing",
    "obstacle_unrepresented",
    "uncertainty_unresolved",
    "conditional_review_required",
    "corner_withheld",
    "phase_routing_incomplete",
    "pool_transport_withheld",
    "transport_mass_residual_exceeded",
    "source_ledger_mismatch",
    "source_schedule_mass_mismatch",
    "obstacle_transport_conditional",
})


@dataclass(frozen=True)
class FieldConditionalReviewAuthorization:
    """Auditable authorization for conditional screening use only."""

    review_id: str
    reviewer_id: str
    reviewer_role: str
    reviewed_at_utc: str
    evidence_id: str
    allowed_scope: str = "conditional_screening_only"

    def __post_init__(self) -> None:
        declared = {
            "review_id": self.review_id,
            "reviewer_id": self.reviewer_id,
            "reviewer_role": self.reviewer_role,
            "evidence_id": self.evidence_id,
        }
        for name, value in declared.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip().lower() == "unspecified"
            ):
                raise ValueError(f"conditional review {name} must be explicitly declared")
        if self.allowed_scope != "conditional_screening_only":
            raise ValueError(
                "conditional review allowed_scope must be "
                "'conditional_screening_only'"
            )
        if (
            not isinstance(self.reviewed_at_utc, str)
            or "T" not in self.reviewed_at_utc
            or not self.reviewed_at_utc.endswith("Z")
        ):
            raise ValueError(
                "conditional review reviewed_at_utc must be an ISO-8601 UTC "
                "timestamp ending in Z"
            )
        try:
            timestamp = datetime.fromisoformat(self.reviewed_at_utc[:-1] + "+00:00")
        except ValueError as error:
            raise ValueError(
                "conditional review reviewed_at_utc must be a valid ISO-8601 UTC timestamp"
            ) from error
        if timestamp.tzinfo is None or timestamp.utcoffset() != timezone.utc.utcoffset(timestamp):
            raise ValueError("conditional review reviewed_at_utc must use UTC")

    def as_record(self) -> dict[str, str]:
        return {
            "review_id": self.review_id,
            "reviewer_id": self.reviewer_id,
            "reviewer_role": self.reviewer_role,
            "reviewed_at_utc": self.reviewed_at_utc,
            "evidence_id": self.evidence_id,
            "allowed_scope": self.allowed_scope,
        }


@dataclass(frozen=True)
class FieldOperationalScreeningDecision:
    """Auditable, deliberately limited disposition of one field result."""

    status: Literal["screening_allowed", "conditional_allowed", "withheld"]
    screening_allowed: bool
    design_basis_allowed: bool
    reasons: tuple[str, ...]
    required_actions: tuple[str, ...]
    refinement_required: bool
    refinement_available: bool
    uncertainty_resolution_required: bool = True
    uncertainty_resolved: bool = False
    gate_codes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.status not in {"screening_allowed", "conditional_allowed", "withheld"}:
            raise ValueError("field operational screening status is unsupported")
        if not all(isinstance(value, bool) for value in (
            self.screening_allowed,
            self.design_basis_allowed,
            self.refinement_required,
            self.refinement_available,
            self.uncertainty_resolution_required,
            self.uncertainty_resolved,
        )):
            raise TypeError("field operational screening decision flags must be boolean")
        if not isinstance(self.reasons, tuple):
            raise TypeError("field operational screening reasons must be a tuple")
        if not isinstance(self.required_actions, tuple):
            raise TypeError("field operational screening actions must be a tuple")
        if not isinstance(self.gate_codes, tuple):
            raise TypeError("field operational screening gate_codes must be a tuple")
        if any(code not in FIELD_OPERATIONAL_GATE_CODES for code in self.gate_codes):
            raise ValueError(
                "field operational screening gate_codes contain an unsupported code"
            )
        if len(set(self.gate_codes)) != len(self.gate_codes):
            raise ValueError("field operational screening gate_codes must be unique")
        if any(not isinstance(value, str) or not value.strip() for value in self.reasons):
            raise ValueError("field operational screening reasons must be non-empty strings")
        if any(not isinstance(value, str) or not value.strip() for value in self.required_actions):
            raise ValueError("field operational screening actions must be non-empty strings")
        if self.design_basis_allowed:
            raise ValueError("reduced-order field screening cannot grant a design-basis decision")
        if self.status == "withheld" and self.screening_allowed:
            raise ValueError("withheld field screening cannot be allowed")
        if self.status == "withheld" and not self.reasons:
            raise ValueError("withheld field screening requires at least one reason")
        if self.status != "withheld" and not self.screening_allowed:
            raise ValueError("an allowed status requires screening_allowed=True")
        if self.screening_allowed and self.refinement_required and not self.refinement_available:
            raise ValueError(
                "screening_allowed cannot be true without required refinement evidence"
            )
        if self.screening_allowed and self.uncertainty_resolution_required and not self.uncertainty_resolved:
            raise ValueError(
                "screening_allowed cannot be true while required uncertainty remains unresolved"
            )
        if self.status == "conditional_allowed" and not self.reasons:
            raise ValueError("conditional_allowed screening requires a recorded reason")
        if self.status == "conditional_allowed" and not self.required_actions:
            raise ValueError(
                "conditional_allowed screening requires a recorded review action"
            )

    @property
    def approval_allowed(self) -> bool:
        """Always false: this model is not an approval or permitting authority."""
        return False

    def require_screening(self) -> None:
        """Raise instead of allowing a caller to consume a withheld result."""
        if not self.screening_allowed:
            detail = "; ".join(self.reasons + self.required_actions)
            raise ValueError("field operational screening is withheld: " + detail)


def _aggregate_gate_codes(
    decisions: tuple[FieldOperationalScreeningDecision, ...],
    *extra: str,
) -> tuple[str, ...]:
    """Preserve corner gate diagnostics on a Cartesian aggregate decision.

    Aggregate envelopes must not reduce a machine-readable corner failure to a
    status string.  Keep the first-seen order so JSON records remain stable,
    and add a neutral fallback for legacy/manual corner decisions that predate
    the typed gate-code field.
    """
    if not isinstance(decisions, tuple):
        raise TypeError("aggregate decision inputs must be a tuple")
    if any(not isinstance(item, FieldOperationalScreeningDecision) for item in decisions):
        raise TypeError("aggregate decision inputs must be FieldOperationalScreeningDecision values")
    codes: list[str] = []
    for decision in decisions:
        codes.extend(decision.gate_codes)
    if any(not decision.screening_allowed for decision in decisions) and not codes:
        codes.append("corner_withheld")
    codes.extend(extra)
    if any(code not in FIELD_OPERATIONAL_GATE_CODES for code in codes):
        raise ValueError("aggregate decision gate_codes contain an unsupported code")
    return tuple(dict.fromkeys(codes))


def _validate_aggregate_operational_decision(
    decision: FieldOperationalScreeningDecision,
    cases: tuple[object, ...],
    *,
    context: str,
) -> None:
    """Check that a Cartesian-envelope decision matches its corner decisions.

    Aggregate constructors are public dataclasses, so checking only the type of
    their ``operational_decision`` is insufficient: a caller could otherwise
    pair a withheld corner set with a fabricated allowed aggregate (or reverse
    the status) before report serialization.
    """
    if not cases:
        return
    decisions = tuple(getattr(case, "decision", None) for case in cases)
    if any(not isinstance(item, FieldOperationalScreeningDecision) for item in decisions):
        raise TypeError(f"{context} cases must expose FieldOperationalScreeningDecision values")
    if any(not item.refinement_required for item in decisions):
        raise ValueError(f"{context} cases must require numerical refinement")
    if any(not item.uncertainty_resolution_required for item in decisions):
        raise ValueError(f"{context} cases must require uncertainty resolution")
    expected_refinement_available = all(item.refinement_available for item in decisions)
    expected_uncertainty_resolved = all(item.uncertainty_resolved for item in decisions)
    if decision.refinement_required is not True:
        raise ValueError(f"{context} aggregate must require refinement")
    if decision.uncertainty_resolution_required is not True:
        raise ValueError(f"{context} aggregate must require uncertainty resolution")
    if decision.refinement_available != expected_refinement_available:
        raise ValueError(f"{context} aggregate refinement availability does not match its corners")
    if decision.uncertainty_resolved != expected_uncertainty_resolved:
        raise ValueError(f"{context} aggregate uncertainty resolution does not match its corners")
    all_allowed = all(item.screening_allowed for item in decisions)
    has_conditional = any(item.status == "conditional_allowed" for item in decisions)
    if not all_allowed:
        expected_statuses = {"withheld"}
    elif has_conditional:
        # Without review opt-in the aggregate is deliberately withheld; with
        # opt-in it may be conditional_allowed.
        expected_statuses = {"withheld", "conditional_allowed"}
    else:
        expected_statuses = {"screening_allowed"}
    if decision.status not in expected_statuses:
        raise ValueError(
            f"{context} aggregate status does not match its corner decisions"
        )


def _validate_case_operational_decision(
    decision: FieldOperationalScreeningDecision,
    screening: object | None,
    refinement: object | None,
    *,
    context: str,
) -> None:
    """Prevent an allowed corner decision from hiding an incomplete result."""
    if not decision.screening_allowed:
        return
    if screening is None or not bool(getattr(screening, "completed", False)):
        raise ValueError(f"{context} allowed decision requires a completed screening result")
    if refinement is None or not bool(getattr(refinement, "completed", False)):
        raise ValueError(f"{context} allowed decision requires a completed refinement result")
    evidence = refinement if refinement is not None else screening
    applicability = getattr(evidence, "applicability", None)
    applicability_status = getattr(applicability, "status", None)
    if applicability_status not in {"accepted", "conditional", "blocked"}:
        raise TypeError(
            f"{context} evidence must expose a valid FieldApplicability status"
        )
    if decision.status == "screening_allowed" and applicability_status != "accepted":
        raise ValueError(
            f"{context} screening_allowed decision requires accepted physical applicability"
        )
    if decision.status == "conditional_allowed" and applicability_status != "conditional":
        raise ValueError(
            f"{context} conditional_allowed decision requires conditional physical applicability"
        )
    if applicability_status == "blocked":
        raise ValueError(
            f"{context} blocked physical applicability cannot carry an allowed decision"
        )


def _declared_sensor_labels(result: FieldSemiFVScreeningResult) -> tuple[str, ...]:
    request = result.request
    primary = () if request.scenario.sensor is None else ("field_sensor",)
    return primary + tuple(item.label for item in request.sensor_deployments)


def _unrepresented_obstacle_labels(
    result: FieldSemiFVScreeningResult,
) -> tuple[str, ...]:
    """Return declared cuboids absent from the local wind-plane solve.

    A cuboid that does not cross the selected plane, or lies wholly upwind,
    may still generate a three-dimensional wake.  The local x-z operator has
    no basis for treating that as harmless.  This is distinct from a cuboid
    intersecting the plane, which is represented as a deliberately conditional
    full-plane obstruction by the conservative semi-FV solver.
    """
    request = result.request
    declared = ((request.obstacle,) if request.obstacle is not None else request.obstacles)
    projections = result.obstacle_projections
    if not projections and result.obstacle_projection is not None:
        projections = (result.obstacle_projection,)
    if len(projections) != len(declared):
        # A projection/result cardinality mismatch is itself unsupported: do
        # not let ``zip`` silently drop a declared obstacle from the gate.
        return tuple(cuboid.label for cuboid in declared)
    return tuple(
        cuboid.label
        for cuboid, projection in zip(declared, projections)
        if projection.obstacle is None
    )


def _unpropagated_field_uncertainty(
    result: FieldSemiFVScreeningResult,
) -> tuple[str, ...]:
    """Return bounded inputs that affect this nominal field result.

    Surface heat-transfer values are intentionally absent: the direct-vapour
    scalar path records but does not use them while post-flash liquid remains
    unrouted. Ambient temperature, pressure and air-density bounds are
    included because they affect the flash and sensor conversion. A declared
    direct-vapour schedule replaces the static source thermodynamic/rate
    boundary, so those upstream bounds are omitted; release geometry and
    timing remain in this check because the local source location, direction
    contract, and finite-duration window are still part of the field case.
    Dedicated measured-history, pool and droplet envelopes retain their own
    source-boundary gates. Generic distributed vapour source geometry and
    vertical-width bounds, plus declared global obstacle geometry bounds, are
    included explicitly as field uncertainty. Evidence-backed lower/upper rate
    histories on distributed sources are also unresolved until their common-
    clock schedule corners have been executed.
    """
    request = result.request
    fields = dict(request.scenario.weather.uncertainty_fields())
    fields.update(request.ambient_uncertainty_fields())
    schedule_uncertainty_labels: list[str] = []
    source_fields = request.scenario.source.uncertainty_fields()
    if request.direct_vapour_schedule is None:
        fields = source_fields | fields
    else:
        # A post-flash schedule supplies the atmospheric source history, but it
        # does not relocate the release or make an unresolved source timing or
        # direction bound disappear.  Keep those fields in the operational
        # gate while deliberately omitting upstream thermodynamic/rate bounds
        # that the schedule replaces.
        fields.update({
            name: value
            for name, value in source_fields.items()
            if name == "duration_s"
            or name.startswith("source_location_")
            or name.startswith("source_direction_")
        })
    for source in request.distributed_vapour_sources:
        fields.update(source.uncertainty_fields())
        if source.has_schedule_uncertainty:
            schedule_uncertainty_labels.append(
                f"distributed source {source.label!r} atmospheric rate schedule bounds"
            )
    for obstacle in request.obstacle_geometry_uncertainty:
        fields.update(obstacle.uncertainty_fields())
    sensor = request.scenario.sensor
    if sensor is not None:
        fields.update({f"field_sensor.{name}": value for name, value in sensor.uncertainty_fields().items()})
    for deployment in request.sensor_deployments:
        fields.update({
            f"{deployment.label}.{name}": value
            for name, value in deployment.sensor.uncertainty_fields().items()
        })
    unresolved = [name for name, value in fields.items() if not value.is_exact]
    unresolved.extend(schedule_uncertainty_labels)
    if (
        request.measured_history_quality is not None
        and not request.measured_history_source_uncertainty_resolved
    ):
        unresolved.append("measured-history source P/T/flow/phase bounds")
    if (
        request.stability_alternatives is not None
        and not request.stability_uncertainty_resolved
    ):
        unresolved.append("declared atmospheric stability alternatives")
    if (
        any(source.source_kind == "pool_vapour" for source in request.distributed_vapour_sources)
        and not request.phase_routing_uncertainty_resolved
    ):
        unresolved.append("phase-routed source/weather/surface pool-vapour bounds")
    return tuple(unresolved)


def evaluate_field_operational_screening(
    result: FieldSemiFVScreeningResult | FieldSemiFVRefinementResult,
    *,
    require_refinement: bool = True,
    require_in_plane_sensor: bool = True,
    require_resolved_uncertainty: bool = True,
    allow_conditional: bool = False,
) -> FieldOperationalScreeningDecision:
    """Decide whether a field result may be used for screening only.

    A conditional physics result needs an explicit caller opt-in. A complete
    numerical refinement study is required by default. Any blocked transport,
    missing local sensor trace, or missing refinement evidence withholds the
    output. The function preserves the model's distinction between screening,
    design basis and approval.
    """
    if not isinstance(result, (FieldSemiFVScreeningResult, FieldSemiFVRefinementResult)):
        raise TypeError("result must be a FieldSemiFVScreeningResult or FieldSemiFVRefinementResult")
    if not all(isinstance(value, bool) for value in (
        require_refinement, require_in_plane_sensor, require_resolved_uncertainty,
    )):
        raise TypeError(
            "require_refinement, require_in_plane_sensor and "
            "require_resolved_uncertainty must be boolean"
        )
    if not isinstance(allow_conditional, bool):
        raise TypeError("allow_conditional must be boolean")

    refinement = result if isinstance(result, FieldSemiFVRefinementResult) else None
    screening = result.screening if refinement is not None else result
    applicability = result.applicability if refinement is not None else screening.applicability
    reasons: list[str] = list(applicability.reasons)
    actions: list[str] = []
    gate_codes: list[str] = []
    if applicability.status == "blocked":
        gate_codes.append("physical_applicability_blocked")
    elif applicability.status == "conditional":
        gate_codes.append("physical_applicability_conditional")
    diagnostic_reasons = applicability.reasons + applicability.warnings
    if any("mass-conservation residual" in item for item in diagnostic_reasons):
        gate_codes.append("transport_mass_residual_exceeded")
    if any("source-mass ledger residual" in item for item in diagnostic_reasons):
        gate_codes.append("source_ledger_mismatch")
    if any("source-mass schedule residual" in item for item in diagnostic_reasons):
        gate_codes.append("source_schedule_mass_mismatch")
    if any("obstacle" in item.lower() for item in diagnostic_reasons):
        gate_codes.append("obstacle_transport_conditional")
    unpropagated_uncertainty = _unpropagated_field_uncertainty(screening)
    uncertainty_resolved = not unpropagated_uncertainty
    unrepresented_obstacles = _unrepresented_obstacle_labels(screening)
    declared_sensor_labels = _declared_sensor_labels(screening)
    sensor_results_by_label = {item.label: item for item in screening.sensor_results}
    missing_sensor_traces = tuple(
        label for label in declared_sensor_labels
        if label not in sensor_results_by_label
        or sensor_results_by_label[label].trace is None
    )

    if not screening.completed:
        reasons.append("field transport did not complete")
        gate_codes.append("transport_incomplete")
    if require_refinement:
        if refinement is None:
            reasons.append("numerical refinement evidence was not supplied")
            actions.append("run run_field_semi_fv_refinement_study before operational screening")
            gate_codes.append("refinement_missing")
        elif not refinement.completed:
            reasons.append("numerical refinement did not complete or did not pass its declared gate")
            actions.append("resolve the refinement applicability reasons and rerun the declared grid/time study")
            gate_codes.append("refinement_failed")
    if require_in_plane_sensor:
        if not declared_sensor_labels:
            reasons.append("no declared in-plane sensor/receptor supports an operational concentration screen")
            actions.append("declare at least one in-plane sensor with position and response calibration")
            gate_codes.append("sensor_missing")
        else:
            if missing_sensor_traces:
                reasons.append(
                    "declared sensor has no local wind-plane trace: " + ", ".join(missing_sensor_traces)
                )
                actions.append("move the sensor into the supported wind plane or use a validated three-dimensional model")
                gate_codes.append("sensor_trace_missing")
    if unrepresented_obstacles:
        reasons.append(
            "declared obstacle is not represented in the local wind-plane transport: "
            + ", ".join(unrepresented_obstacles)
        )
        actions.append(
            "use a validated three-dimensional obstacle/wake model or demonstrate that each omitted obstacle is outside the decision domain"
        )
        gate_codes.append("obstacle_unrepresented")
    if require_resolved_uncertainty and unpropagated_uncertainty:
        reasons.append(
            "declared field uncertainty was not propagated through this nominal result: "
            + ", ".join(unpropagated_uncertainty)
        )
        actions.append(
            "run the matching deterministic source/weather/surface/stability/sensor uncertainty envelope "
            "before operational screening; direct scalar results are not probability intervals"
        )
        gate_codes.append("uncertainty_unresolved")

    reasons.extend(applicability.warnings)
    reasons = list(dict.fromkeys(reasons))
    actions = list(dict.fromkeys(actions))
    gate_codes = list(dict.fromkeys(gate_codes))
    refinement_available = refinement is not None and refinement.completed
    if reasons and applicability.status == "blocked":
        return FieldOperationalScreeningDecision(
            "withheld", False, False, tuple(reasons), tuple(actions),
            require_refinement, refinement_available,
            require_resolved_uncertainty, uncertainty_resolved,
            tuple(gate_codes),
        )
    # Requirements added above are non-negotiable, independently of whether
    # the physical applicability record is conditional or accepted.
    nonnegotiable = (
        not screening.completed
        or (require_refinement and not refinement_available)
        or (require_in_plane_sensor and bool(missing_sensor_traces))
        or (require_in_plane_sensor and not declared_sensor_labels)
        or bool(unrepresented_obstacles)
        or (require_resolved_uncertainty and not uncertainty_resolved)
    )
    if nonnegotiable:
        return FieldOperationalScreeningDecision(
            "withheld", False, False, tuple(reasons), tuple(actions),
            require_refinement, refinement_available,
            require_resolved_uncertainty, uncertainty_resolved,
            tuple(gate_codes),
        )
    if applicability.status == "conditional":
        gate_codes.append("conditional_review_required")
        if allow_conditional:
            actions.append("review every conditional warning before acting on the screening result")
            return FieldOperationalScreeningDecision(
                "conditional_allowed", True, False, tuple(reasons), tuple(dict.fromkeys(actions)),
                require_refinement, refinement_available,
                require_resolved_uncertainty, uncertainty_resolved,
                tuple(dict.fromkeys(gate_codes)),
            )
        actions.append("set allow_conditional=True only after an accountable reviewer accepts the recorded scope limits")
        return FieldOperationalScreeningDecision(
            "withheld", False, False, tuple(reasons), tuple(dict.fromkeys(actions)),
            require_refinement, refinement_available,
            require_resolved_uncertainty, uncertainty_resolved,
            tuple(dict.fromkeys(gate_codes)),
        )
    return FieldOperationalScreeningDecision(
        "screening_allowed", True, False, tuple(reasons), tuple(actions),
        require_refinement, refinement_available,
        require_resolved_uncertainty, uncertainty_resolved,
        tuple(gate_codes),
    )


def field_operational_screening_decision_record(
    decision: FieldOperationalScreeningDecision,
) -> dict[str, object]:
    """Return a JSON-native decision record suitable for a batch manifest."""
    if not isinstance(decision, FieldOperationalScreeningDecision):
        raise TypeError("decision must be a FieldOperationalScreeningDecision")
    return {
        "schema": FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA,
        "status": decision.status,
        "screening_allowed": decision.screening_allowed,
        "design_basis_allowed": decision.design_basis_allowed,
        "approval_allowed": decision.approval_allowed,
        "reasons": list(decision.reasons),
        "required_actions": list(decision.required_actions),
        "refinement_required": decision.refinement_required,
        "refinement_available": decision.refinement_available,
        "uncertainty_resolution_required": decision.uncertainty_resolution_required,
        "uncertainty_resolved": decision.uncertainty_resolved,
        "gate_codes": list(decision.gate_codes),
    }


__all__ = [
    "FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA", "FIELD_OPERATIONAL_GATE_CODES",
    "FieldConditionalReviewAuthorization",
    "FieldOperationalScreeningDecision",
    "evaluate_field_operational_screening", "field_operational_screening_decision_record",
]
