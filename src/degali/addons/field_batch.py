"""Fail-safe export of auditable field-screening case batches."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from collections import Counter
from typing import Mapping

from .field_decision import (
    FieldConditionalReviewAuthorization,
    FieldOperationalScreeningDecision,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from .field_report import (
    FIELD_SCREENING_REPORT_SCHEMA,
    field_refinement_report,
    field_screening_report,
)
from .field_operational_envelope import (
    FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
    FieldOperationalUncertaintyEnvelope,
    field_operational_uncertainty_envelope_report,
    run_field_operational_uncertainty_envelope,
)
from .field_source_io import FieldAtmosphericSourceSchedule
from .field_json import strict_json_loads
from .field_workflow import (
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    run_field_semi_fv_refinement_study,
    run_field_semi_fv_screening,
)
FIELD_BATCH_MANIFEST_SCHEMA = "degali.field-batch.v3"
FIELD_BATCH_SUMMARY_SCHEMA = "degali.field-batch-summary.v1"
_SAFE_CASE_LABEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")
_FIELD_BATCH_REPORT_SCHEMAS = frozenset({
    FIELD_SCREENING_REPORT_SCHEMA,
    FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
})


def _write_json_exclusive(path: Path, payload: str, *, artifact: str) -> None:
    """Create one batch JSON artifact without replacing an existing path."""
    try:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing field batch {artifact}: {path}"
        ) from error


def _read_batch_json(path: Path, *, artifact: str) -> dict[str, object]:
    """Read one emitted JSON artifact before accepting its typed wrapper."""
    try:
        value = strict_json_loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"batch {artifact} must be valid UTF-8 JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"batch {artifact} must contain a JSON object")
    return value


def _sha256_file(path: Path) -> str:
    """Return the content fingerprint used by batch artifact provenance."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha256_payload(payload: str) -> str:
    """Return the UTF-8 fingerprint of a payload before it is written."""
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def field_batch_decision_summary(
    cases: tuple["FieldBatchCase", ...] | list["FieldBatchCase"],
) -> dict[str, object]:
    """Return a fail-safe aggregate of the already evaluated case decisions.

    This summary is intentionally derived from the typed operational decisions
    rather than from physical applicability labels.  A conditional physical
    case that was not authorized therefore remains counted as withheld, and a
    batch cannot look fully usable merely because every numerical report was
    written successfully.
    """
    if not isinstance(cases, (tuple, list)) or not cases:
        raise ValueError("field batch decision summary requires at least one case")
    if any(not isinstance(case, FieldBatchCase) for case in cases):
        raise TypeError("field batch decision summary cases must be FieldBatchCase values")
    status_counts = {
        "screening_allowed": 0,
        "conditional_allowed": 0,
        "withheld": 0,
    }
    for case in cases:
        status_counts[case.operational_decision.status] += 1
    screening_allowed = sum(
        case.operational_decision.screening_allowed for case in cases
    )
    gate_code_counts = dict(sorted(
        Counter(
            code
            for case in cases
            for code in case.operational_decision.gate_codes
        ).items()
    ))
    return {
        "schema": FIELD_BATCH_SUMMARY_SCHEMA,
        "case_count": len(cases),
        "decision_status_counts": status_counts,
        "screening_allowed_case_count": screening_allowed,
        "conditional_allowed_case_count": status_counts["conditional_allowed"],
        "withheld_case_count": status_counts["withheld"],
        "all_screening_allowed": screening_allowed == len(cases),
        "gate_code_counts": gate_code_counts,
    }


@dataclass(frozen=True)
class FieldBatchCase:
    """One written case report and its explicit applicability status."""

    label: str
    report_path: Path
    applicability_status: str
    refinement_requested: bool
    refinement_completed: bool
    operational_decision: FieldOperationalScreeningDecision
    conditional_review: FieldConditionalReviewAuthorization | None = None
    conditional_review_applied: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or _SAFE_CASE_LABEL.fullmatch(self.label) is None:
            raise ValueError("batch case label must match the safe report-name pattern")
        if not isinstance(self.report_path, Path):
            raise TypeError("batch case report_path must be a Path")
        if self.applicability_status not in {"accepted", "conditional", "blocked"}:
            raise ValueError("batch case applicability_status must be accepted, conditional or blocked")
        for name, value in {
            "refinement_requested": self.refinement_requested,
            "refinement_completed": self.refinement_completed,
            "conditional_review_applied": self.conditional_review_applied,
        }.items():
            if not isinstance(value, bool):
                raise TypeError(f"batch case {name} must be boolean")
        if not isinstance(self.operational_decision, FieldOperationalScreeningDecision):
            raise TypeError("batch case operational_decision must be a FieldOperationalScreeningDecision")
        if self.conditional_review is not None and not isinstance(
            self.conditional_review, FieldConditionalReviewAuthorization
        ):
            raise TypeError(
                "batch case conditional_review must be a "
                "FieldConditionalReviewAuthorization or None"
            )
        if self.conditional_review_applied and self.conditional_review is None:
            raise ValueError("batch case cannot apply a missing conditional review")
        if self.operational_decision.status == "conditional_allowed" and (
            self.conditional_review is None or not self.conditional_review_applied
        ):
            raise ValueError(
                "conditional_allowed batch case requires an applied conditional review"
            )
        decision_status = self.operational_decision.status
        if decision_status == "screening_allowed" and self.applicability_status != "accepted":
            raise ValueError(
                "screening_allowed batch case requires accepted physical applicability"
            )
        if decision_status == "conditional_allowed" and self.applicability_status != "conditional":
            raise ValueError(
                "conditional_allowed batch case requires conditional physical applicability"
            )
        if self.applicability_status == "blocked" and self.operational_decision.screening_allowed:
            raise ValueError(
                "blocked physical applicability cannot carry an allowed batch decision"
            )
        if self.refinement_completed and not self.refinement_requested:
            raise ValueError(
                "batch case cannot report completed refinement when refinement was not requested"
            )
        if self.refinement_completed != self.operational_decision.refinement_available:
            raise ValueError(
                "batch case refinement_completed must match operational decision availability"
            )


@dataclass(frozen=True)
class FieldBatchExport:
    """Locations and statuses emitted by one non-overwriting batch export."""

    output_directory: Path
    manifest_path: Path
    cases: tuple[FieldBatchCase, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.output_directory, Path) or not isinstance(self.manifest_path, Path):
            raise TypeError("batch export paths must be Path values")
        if not isinstance(self.cases, tuple) or not self.cases:
            raise ValueError("batch export requires at least one case")
        if any(not isinstance(case, FieldBatchCase) for case in self.cases):
            raise TypeError("batch export cases must contain only FieldBatchCase values")
        if len({case.label for case in self.cases}) != len(self.cases):
            raise ValueError("batch export case labels must be unique")
        if len({case.report_path for case in self.cases}) != len(self.cases):
            raise ValueError("batch export report paths must be unique")
        if not self.output_directory.is_dir():
            raise ValueError("batch export output_directory must exist and be a directory")
        if self.manifest_path.name != "manifest.json":
            raise ValueError("batch export manifest_path must be manifest.json")
        if self.manifest_path.resolve().parent != self.output_directory.resolve():
            raise ValueError("batch export manifest_path must be inside output_directory")
        if not self.manifest_path.is_file():
            raise ValueError("batch export manifest_path must point to an existing file")
        manifest = _read_batch_json(self.manifest_path, artifact="manifest")
        if manifest.get("schema") != FIELD_BATCH_MANIFEST_SCHEMA:
            raise ValueError("batch manifest schema does not match the typed export")
        raw_cases = manifest.get("cases")
        if not isinstance(raw_cases, list) or not raw_cases:
            raise ValueError("batch manifest cases must be a non-empty array")
        manifest_labels = []
        manifest_by_label: dict[str, dict[str, object]] = {}
        for raw_case in raw_cases:
            if not isinstance(raw_case, dict):
                raise ValueError("batch manifest cases must contain JSON objects")
            label = raw_case.get("label")
            if not isinstance(label, str) or _SAFE_CASE_LABEL.fullmatch(label) is None:
                raise ValueError("batch manifest case labels must match the safe report-name pattern")
            if label in manifest_by_label:
                raise ValueError("batch manifest case labels must be unique")
            manifest_labels.append(label)
            manifest_by_label[label] = raw_case
        expected_labels = [case.label for case in self.cases]
        if manifest_labels != expected_labels:
            raise ValueError("batch manifest cases do not match the typed export cases")
        summary = manifest.get("summary")
        if summary != field_batch_decision_summary(self.cases):
            raise ValueError("batch manifest summary does not match the typed export decisions")
        for case in self.cases:
            if case.report_path.name != f"{case.label}.json":
                raise ValueError(
                    "batch export report path must match its case label"
                )
            if case.report_path.resolve().parent != self.output_directory.resolve():
                raise ValueError(
                    "batch export report paths must be inside output_directory"
                )
            if not case.report_path.is_file():
                raise ValueError(
                    "batch export report paths must point to existing files"
                )
            if case.report_path.resolve() == self.manifest_path.resolve():
                raise ValueError("batch export report path cannot be manifest.json")
            manifest_case = manifest_by_label[case.label]
            if manifest_case.get("report_file") != case.report_path.name:
                raise ValueError(
                    f"batch manifest report_file does not match case {case.label!r}"
                )
            if manifest_case.get("report_sha256") != _sha256_file(case.report_path):
                raise ValueError(
                    f"batch manifest report_sha256 does not match case {case.label!r}"
                )
            if manifest_case.get("applicability_status") != case.applicability_status:
                raise ValueError(
                    f"batch manifest applicability does not match case {case.label!r}"
                )
            expected_review = (
                None if case.conditional_review is None else {
                    **case.conditional_review.as_record(),
                    "applied_to_this_execution": case.conditional_review_applied,
                }
            )
            if manifest_case.get("conditional_review") != expected_review:
                raise ValueError(
                    f"batch manifest conditional review does not match case {case.label!r}"
                )
            for name, expected in (
                ("refinement_requested", case.refinement_requested),
                ("refinement_completed", case.refinement_completed),
            ):
                if manifest_case.get(name) is not expected:
                    raise ValueError(
                        f"batch manifest {name} does not match case {case.label!r}"
                    )
            decision = manifest_case.get("operational_decision")
            expected_decision = field_operational_screening_decision_record(
                case.operational_decision
            )
            if decision != expected_decision:
                raise ValueError(
                    f"batch manifest operational decision does not match case {case.label!r}"
                )
            report = _read_batch_json(case.report_path, artifact=f"report {case.label!r}")
            report_schema = report.get("schema")
            if report_schema not in _FIELD_BATCH_REPORT_SCHEMAS:
                raise ValueError(
                    f"batch report {case.label!r} has an unsupported schema"
                )
            execution_kind = manifest_case.get("execution_kind")
            expected_execution_kind = (
                "operational_uncertainty_envelope"
                if report_schema == FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA
                else "single_field_screen"
            )
            if execution_kind != expected_execution_kind:
                raise ValueError(
                    f"batch report/execution kind does not match case {case.label!r}"
                )


def export_field_screening_batch(
    cases: Mapping[str, FieldSemiFVRequest],
    output_directory: str | Path,
    *,
    include_refinement: bool = True,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    max_cell_steps: int = 20_000_000,
    allow_conditional_operational_screening: bool = False,
    conditional_review_authorizations: Mapping[
        str, FieldConditionalReviewAuthorization
    ] | None = None,
    require_refinement_for_operational_screening: bool = True,
    require_in_plane_sensor_for_operational_screening: bool = True,
    require_resolved_uncertainty_for_operational_screening: bool = True,
    atmospheric_source_schedules: Mapping[str, FieldAtmosphericSourceSchedule] | None = None,
    max_uncertainty_cases: int = 64,
    table_nodes: int = 161,
) -> FieldBatchExport:
    """Run named cases and emit one fail-safe JSON record per case.

    Every requested case is computed before any file is created. This prevents
    a model exception from being confused with a complete-looking partial
    export. The target directory may be new or empty, but is never overwritten.
    With the default ``include_refinement=True``, each case report also carries
    its explicit grid/time convergence result or an explicit blocked reason.
    The manifest separately records a fail-safe operational-screening
    decision. Skipping refinement does not silently waive the default
    refinement requirement for that decision. Conditional opt-in requires a
    separate review authorization for every named case; each record and
    whether it was applied are retained in manifest v3.
    """
    if not cases:
        raise ValueError("at least one named field case is required")
    for name, value in {
        "include_refinement": include_refinement,
        "allow_conditional_operational_screening": allow_conditional_operational_screening,
        "require_refinement_for_operational_screening": require_refinement_for_operational_screening,
        "require_in_plane_sensor_for_operational_screening": require_in_plane_sensor_for_operational_screening,
        "require_resolved_uncertainty_for_operational_screening": require_resolved_uncertainty_for_operational_screening,
    }.items():
        if not isinstance(value, bool):
            raise TypeError(f"{name} must be boolean")
    schedules = dict(atmospheric_source_schedules or {})
    if any(not isinstance(label, str) for label in schedules):
        raise TypeError("atmospheric source schedule case labels must be strings")
    unknown_schedules = sorted(set(schedules) - set(cases))
    if unknown_schedules:
        raise ValueError(
            "atmospheric source schedule supplied for unknown batch case: "
            + ", ".join(unknown_schedules)
        )
    if not isinstance(max_uncertainty_cases, int) or isinstance(max_uncertainty_cases, bool) or max_uncertainty_cases < 1:
        raise ValueError("max_uncertainty_cases must be a positive integer")
    if not isinstance(table_nodes, int) or isinstance(table_nodes, bool) or table_nodes < 2:
        raise ValueError("table_nodes must be an integer of at least two")
    if schedules and not (
        require_refinement_for_operational_screening
        and require_in_plane_sensor_for_operational_screening
        and require_resolved_uncertainty_for_operational_screening
    ):
        raise ValueError(
            "atmospheric source schedule batch cases require refinement, in-plane "
            "sensor, and resolved-uncertainty operational gates"
        )
    for label, schedule in schedules.items():
        if not isinstance(schedule, FieldAtmosphericSourceSchedule):
            raise TypeError(
                f"atmospheric source schedule for batch case {label!r} must be "
                "FieldAtmosphericSourceSchedule"
            )
    reviews = dict(conditional_review_authorizations or {})
    if any(not isinstance(label, str) for label in reviews):
        raise TypeError("conditional review case labels must be strings")
    unknown_reviews = sorted(set(reviews) - set(cases))
    if unknown_reviews:
        raise ValueError(
            "conditional review supplied for unknown batch case: "
            + ", ".join(unknown_reviews)
        )
    for label, review in reviews.items():
        if not isinstance(review, FieldConditionalReviewAuthorization):
            raise TypeError(
                f"conditional review for batch case {label!r} must be "
                "FieldConditionalReviewAuthorization"
            )
    if allow_conditional_operational_screening:
        missing_reviews = sorted(set(cases) - set(reviews))
        if missing_reviews:
            raise ValueError(
                "conditional operational screening requires a review authorization "
                "for every batch case; missing: " + ", ".join(missing_reviews)
            )
    target = Path(output_directory)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty field batch directory: {target}")
    prepared: list[tuple[
        str,
        FieldSemiFVScreeningResult | FieldSemiFVRefinementResult | FieldOperationalUncertaintyEnvelope,
        FieldOperationalScreeningDecision,
        FieldConditionalReviewAuthorization | None,
        bool,
    ]] = []
    for label, request in cases.items():
        if not isinstance(label, str) or not _SAFE_CASE_LABEL.fullmatch(label):
            raise ValueError("case labels must be safe file names of 1 to 128 letters, digits, . _ or -")
        if not isinstance(request, FieldSemiFVRequest):
            raise TypeError("every batch case must be a FieldSemiFVRequest")
        source_schedule = schedules.get(label)
        if source_schedule is not None:
            result = run_field_operational_uncertainty_envelope(
                request,
                atmospheric_source_schedule=source_schedule,
                max_cases=max_uncertainty_cases,
                table_nodes=table_nodes,
                include_refinement=include_refinement,
                refinement_factors=refinement_factors,
                relative_tolerance=relative_tolerance,
                max_cell_steps=max_cell_steps,
                allow_conditional=allow_conditional_operational_screening,
            )
            decision = result.operational_decision
            refinement_completed = all(
                case.refinement is not None and case.refinement.completed
                for case in result.cases
            )
        else:
            if include_refinement:
                result = run_field_semi_fv_refinement_study(
                    request,
                    refinement_factors=refinement_factors,
                    relative_tolerance=relative_tolerance,
                    max_cell_steps=max_cell_steps,
                )
            else:
                result = run_field_semi_fv_screening(request)
            decision = evaluate_field_operational_screening(
                result,
                require_refinement=require_refinement_for_operational_screening,
                require_in_plane_sensor=require_in_plane_sensor_for_operational_screening,
                require_resolved_uncertainty=require_resolved_uncertainty_for_operational_screening,
                allow_conditional=allow_conditional_operational_screening,
            )
            refinement_completed = isinstance(result, FieldSemiFVRefinementResult) and result.completed
        prepared.append((label, result, decision, reviews.get(label), refinement_completed))

    # Render and validate every report payload before creating the output
    # directory.  A late JSON/finite-value error must not leave an earlier
    # case report that looks like a complete batch.
    rendered_reports = []
    for label, result, decision, conditional_review, refinement_completed in prepared:
        if isinstance(result, FieldOperationalUncertaintyEnvelope):
            report_record = field_operational_uncertainty_envelope_report(result)
            applicability_status = (
                "blocked" if any(
                    (
                        case.refinement.applicability
                        if case.refinement is not None else case.screening.applicability
                    ).status == "blocked"
                    for case in result.cases
                ) else "conditional" if any(
                    (
                        case.refinement.applicability
                        if case.refinement is not None else case.screening.applicability
                    ).status == "conditional"
                    for case in result.cases
                ) else "accepted"
            )
        elif isinstance(result, FieldSemiFVRefinementResult):
            report_record = field_refinement_report(result)
            applicability_status = result.applicability.status
        else:
            report_record = field_screening_report(result)
            applicability_status = result.applicability.status
        report_payload = json.dumps(
            report_record, indent=2, ensure_ascii=False, allow_nan=False,
        ) + "\n"
        rendered_reports.append((
            label, result, decision, conditional_review, refinement_completed,
            applicability_status, report_payload,
        ))

    target.mkdir(parents=True, exist_ok=True)
    exported = []
    execution_kind_by_label = {
        label: (
            "operational_uncertainty_envelope"
            if isinstance(result, FieldOperationalUncertaintyEnvelope)
            else "single_field_screen"
        )
        for label, result, _decision, _review, _completed, _status, _payload in rendered_reports
    }
    report_payload_by_label = {}
    for (
        label, result, decision, conditional_review, refinement_completed,
        applicability_status, report_payload,
    ) in rendered_reports:
        report_path = target / f"{label}.json"
        _write_json_exclusive(report_path, report_payload, artifact="case report")
        report_payload_by_label[label] = report_payload
        exported.append(FieldBatchCase(
            label=label,
            report_path=report_path,
            applicability_status=applicability_status,
            refinement_requested=include_refinement,
            refinement_completed=refinement_completed,
            operational_decision=decision,
            conditional_review=conditional_review,
            conditional_review_applied=(
                allow_conditional_operational_screening
                and conditional_review is not None
            ),
        ))
    manifest_path = target / "manifest.json"
    manifest_payload = json.dumps({
            "schema": FIELD_BATCH_MANIFEST_SCHEMA,
            "include_refinement": include_refinement,
            "refinement_factors": list(refinement_factors) if include_refinement else None,
            "relative_tolerance": relative_tolerance if include_refinement else None,
            "operational_screening_policy": {
                "allow_conditional": allow_conditional_operational_screening,
                "require_refinement": require_refinement_for_operational_screening,
                "require_in_plane_sensor": require_in_plane_sensor_for_operational_screening,
                "require_resolved_uncertainty": require_resolved_uncertainty_for_operational_screening,
                "conditional_review_required": allow_conditional_operational_screening,
            },
            "summary": field_batch_decision_summary(exported),
            "cases": [
                {
                    "label": item.label,
                    "report_file": item.report_path.name,
                    "report_sha256": _sha256_payload(report_payload_by_label[item.label]),
                    "applicability_status": item.applicability_status,
                    "refinement_requested": item.refinement_requested,
                    "refinement_completed": item.refinement_completed,
                    "execution_kind": execution_kind_by_label[item.label],
                    "conditional_review": (
                        None if item.conditional_review is None else {
                            **item.conditional_review.as_record(),
                            "applied_to_this_execution": item.conditional_review_applied,
                        }
                    ),
                    "operational_decision": field_operational_screening_decision_record(
                        item.operational_decision
                    ),
                }
                for item in exported
            ],
        }, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    _write_json_exclusive(manifest_path, manifest_payload, artifact="manifest")
    return FieldBatchExport(target, manifest_path, tuple(exported))


__all__ = [
    "FIELD_BATCH_MANIFEST_SCHEMA", "FIELD_BATCH_SUMMARY_SCHEMA",
    "FieldBatchCase", "FieldBatchExport",
    "field_batch_decision_summary",
    "export_field_screening_batch",
]
