"""Post-execution integrity checks for field artifacts.

The field CLI writes execution JSON with an input path and digest, but a saved
execution record is not evidence that its inputs remain unchanged. This module
reopens supported execution records, verifies the recorded input bytes, re-runs
the strict source/model/validation boundary, and (where deterministic) compares
the stored report with a fresh report. It never promotes a result or changes an
operational decision.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .field_comparison import (
    compare_field_model_sensor_sets,
    field_model_comparison_report,
)
from .field_comparison_io import read_field_model_comparison_case_json
from .field_json import strict_json_loads
from .field_source_io import (
    field_source_schedule_record,
    read_field_source_schedule_json,
)
from .field_case_io import read_field_screening_case_json
from .field_historian_io import ImportedPressureDrivenMeasuredHistory
from .field_batch_io import read_field_batch_input_json
from .field_decision import (
    FIELD_OPERATIONAL_GATE_CODES,
    FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from .field_operational_envelope import (
    FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
    field_operational_uncertainty_envelope_report,
    run_field_operational_uncertainty_envelope,
)
from .field_operational_history import (
    FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA,
    field_operational_joint_measured_history_envelope_report,
    run_field_operational_joint_measured_history_envelope,
    run_field_operational_joint_pressure_driven_history_envelope,
)
from .field_operational_phase_transport import (
    FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA,
    field_operational_phase_routing_transport_envelope_report,
    run_field_operational_phase_routing_transport_envelope,
)
from .field_operational_sensor_array import (
    FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA,
    field_operational_sensor_array_envelope_report,
    run_field_operational_sensor_array_envelope,
)
from .field_operational_source_sensor import (
    FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA,
    field_operational_source_sensor_envelope_report,
    run_field_operational_source_sensor_envelope,
)
from .field_report import (
    FIELD_DROPLET_HANDOFF_REFINEMENT_REPORT_SCHEMA,
    FIELD_JOINT_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA,
    FIELD_SCREENING_REPORT_SCHEMA,
    FIELD_SENSOR_ARRAY_UNCERTAINTY_ENVELOPE_REPORT_SCHEMA,
    FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA,
    FIELD_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA,
    FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA,
    field_refinement_report,
    field_screening_report,
)
from .field_validation import (
    field_validation_score_record,
    score_field_model_against_validation,
)
from .field_validation_io import read_field_validation_case_json
from .field_workflow import (
    run_field_semi_fv_refinement_study,
    run_field_semi_fv_screening,
)


FIELD_EXECUTION_VERIFICATION_SCHEMA = "degali.field-execution-verification.v1"

_SUPPORTED_EXECUTION_SCHEMAS = frozenset({
    "degali.field-atmospheric-schedule-execution.v1",
    "degali.field-model-comparison-execution.v1",
    "degali.field-validation-execution.v1",
    "degali.field-screening-execution.v1",
    "degali.field-batch-execution.v1",
})

_SUPPORTED_BATCH_REPORT_SCHEMAS = frozenset({
    FIELD_SCREENING_REPORT_SCHEMA,
    FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
})

_SUPPORTED_FIELD_REPORT_SCHEMAS = frozenset({
    FIELD_SCREENING_REPORT_SCHEMA,
    FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA,
    FIELD_DROPLET_HANDOFF_REFINEMENT_REPORT_SCHEMA,
    FIELD_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA,
    FIELD_JOINT_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA,
    FIELD_SENSOR_ARRAY_UNCERTAINTY_ENVELOPE_REPORT_SCHEMA,
    FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA,
    FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
    FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA,
    FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA,
    FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA,
    FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA,
})

# Keep the post-execution gate aligned with the field workflow's numerical
# applicability gate.  The verifier checks only the serialized contract for
# reports whose full inputs/options are not available (for example, legacy or
# option-incomplete envelope records); deterministic phase/history/sensor
# options are replayed when present.
_FIELD_MASS_RESIDUAL_REL_TOL = 1.0e-8
_FIELD_MASS_RESIDUAL_ABS_TOL_KG = 1.0e-12


def _mapping(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _source_schedule_record_with_defaults(
    value: object,
    name: str,
) -> dict[str, Any]:
    """Normalize legacy schedule records before provenance comparison.

    ``piecewise_constant`` was the only schedule operator in the v1 artifact
    schema.  A missing operator therefore has an unambiguous meaning and is
    upgraded in-memory; an explicitly supplied operator is never rewritten.
    """
    record = dict(_mapping(value, name))
    schedule_value = record.get("schedule")
    if isinstance(schedule_value, Mapping) and "rate_operator" not in schedule_value:
        schedule = dict(schedule_value)
        schedule["rate_operator"] = "piecewise_constant"
        record["schedule"] = schedule
    return record


def _validate_batch_summary_against_cases(
    summary: Mapping[str, Any],
    manifest_cases: list[Mapping[str, Any]],
) -> None:
    """Recompute the fail-safe batch summary from manifest decisions."""
    status_counts = {
        "screening_allowed": 0,
        "conditional_allowed": 0,
        "withheld": 0,
    }
    gate_code_counts: dict[str, int] = {}
    saw_gate_codes = False
    screening_allowed = 0
    for index, item in enumerate(manifest_cases):
        decision = _mapping(item.get("operational_decision"), f"manifest case {index}.operational_decision")
        status = decision.get("status")
        allowed = decision.get("screening_allowed")
        if status not in status_counts or not isinstance(allowed, bool):
            raise ValueError(f"manifest case {index} has an invalid operational decision")
        status_counts[status] += 1
        screening_allowed += int(allowed)
        codes = decision.get("gate_codes")
        if codes is None:
            continue
        saw_gate_codes = True
        if (
            not isinstance(codes, list)
            or any(
                not isinstance(code, str)
                or code not in FIELD_OPERATIONAL_GATE_CODES
                for code in codes
            )
            or len(set(codes)) != len(codes)
        ):
            raise ValueError(f"manifest case {index} has invalid operational gate_codes")
        for code in codes:
            gate_code_counts[code] = gate_code_counts.get(code, 0) + 1
    expected = {
        "schema": "degali.field-batch-summary.v1",
        "case_count": len(manifest_cases),
        "decision_status_counts": status_counts,
        "screening_allowed_case_count": screening_allowed,
        "conditional_allowed_case_count": status_counts["conditional_allowed"],
        "withheld_case_count": status_counts["withheld"],
        "all_screening_allowed": screening_allowed == len(manifest_cases),
    }
    if saw_gate_codes or "gate_code_counts" in summary:
        expected["gate_code_counts"] = dict(sorted(gate_code_counts.items()))
    if dict(summary) != expected:
        raise ValueError("batch manifest summary does not match its case decisions")


def _validate_operational_decision_record(
    decision: Mapping[str, Any], *, context: str,
) -> None:
    """Validate the fail-safe invariants of a serialized operational decision."""
    required = {
        "status", "screening_allowed", "reasons", "required_actions",
        "refinement_required", "refinement_available",
        "uncertainty_resolution_required", "uncertainty_resolved",
    }
    missing = sorted(required.difference(decision))
    if missing:
        raise ValueError(f"{context} is missing {', '.join(missing)}")
    status = decision["status"]
    if status not in {"screening_allowed", "conditional_allowed", "withheld"}:
        raise ValueError(f"{context}.status is invalid")
    for key in (
        "screening_allowed", "refinement_required", "refinement_available",
        "uncertainty_resolution_required", "uncertainty_resolved",
    ):
        if not isinstance(decision[key], bool):
            raise ValueError(f"{context}.{key} must be boolean")
    if "design_basis_allowed" in decision and decision["design_basis_allowed"] is not False:
        raise ValueError(f"{context}.design_basis_allowed must be false")
    if "approval_allowed" in decision and decision["approval_allowed"] is not False:
        raise ValueError(f"{context}.approval_allowed must be false")
    reasons = decision["reasons"]
    actions = decision["required_actions"]
    if not isinstance(reasons, list) or any(
        not isinstance(value, str) or not value.strip() for value in reasons
    ):
        raise ValueError(f"{context}.reasons must be an array of non-empty strings")
    if not isinstance(actions, list) or any(
        not isinstance(value, str) or not value.strip() for value in actions
    ):
        raise ValueError(f"{context}.required_actions must be an array of non-empty strings")
    if status == "withheld":
        if decision["screening_allowed"] or not reasons:
            raise ValueError(f"{context} withheld decision must be disallowed and explained")
    elif not decision["screening_allowed"]:
        raise ValueError(f"{context} allowed status must set screening_allowed=true")
    if status == "conditional_allowed" and (not reasons or not actions):
        raise ValueError(f"{context} conditional decision requires reasons and actions")
    gate_codes = decision.get("gate_codes")
    if gate_codes is not None:
        if (
            not isinstance(gate_codes, list)
            or any(
                not isinstance(code, str)
                or code not in FIELD_OPERATIONAL_GATE_CODES
                for code in gate_codes
            )
            or len(set(gate_codes)) != len(gate_codes)
        ):
            raise ValueError(f"{context}.gate_codes is invalid")


def _validate_report_decision_applicability(
    report: Mapping[str, Any],
    decision: Mapping[str, Any],
    *,
    context: str,
) -> None:
    """Keep serialized physical applicability and operational status aligned.

    Integrity-only verification cannot rely on a fresh typed result.  A saved
    report must therefore not be able to claim an allowed operational decision
    after its physical applicability was edited to ``blocked`` (or to a
    different conditional/accepted status).  Envelope reports carry one
    decision per case, so recurse only through explicit ``field_result`` /
    ``operational_screening`` pairs and never compare an aggregate decision to
    a nested corner's applicability.
    """
    schema = report.get("schema")
    if schema == FIELD_SCREENING_REPORT_SCHEMA:
        raw_applicability = report.get("refinement_applicability")
        if raw_applicability is None:
            raw_applicability = report.get("applicability")
        applicability = _mapping(raw_applicability, f"{context}.applicability")
        status = applicability.get("status")
        if status not in {"accepted", "conditional", "blocked"}:
            raise ValueError(f"{context}.applicability.status is invalid")
        decision_status = decision.get("status")
        allowed = decision.get("screening_allowed")
        if decision_status == "screening_allowed" and status != "accepted":
            raise ValueError(
                f"{context} screening_allowed decision requires accepted physical applicability"
            )
        if decision_status == "conditional_allowed" and status != "conditional":
            raise ValueError(
                f"{context} conditional_allowed decision requires conditional physical applicability"
            )
        if allowed is True and status == "blocked":
            raise ValueError(
                f"{context} blocked physical applicability cannot carry an allowed decision"
            )
        return
    cases = report.get("cases")
    if not isinstance(cases, list):
        return
    for index, raw_case in enumerate(cases):
        if not isinstance(raw_case, Mapping):
            continue
        field_result = raw_case.get("field_result")
        operational = raw_case.get("operational_screening")
        if isinstance(field_result, Mapping) and isinstance(operational, Mapping):
            _validate_operational_decision_record(
                operational,
                context=f"{context}.cases[{index}].operational_screening",
            )
            _validate_report_decision_applicability(
                field_result,
                operational,
                context=f"{context}.cases[{index}].field_result",
            )


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _input_file(
    payload: Mapping[str, Any], *, artifact: Path, allow_extra: bool = False,
) -> Path:
    data = _mapping(payload.get("input"), "execution input")
    required = {"path", "sha256"}
    if not required.issubset(data):
        raise ValueError("execution input must contain path and sha256")
    if not allow_extra and set(data) != required:
        raise ValueError("execution input must contain exactly path and sha256")
    raw_path = data["path"]
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise ValueError("execution input.path must be a non-empty string")
    raw_digest = data["sha256"]
    if (
        not isinstance(raw_digest, str)
        or len(raw_digest) != 64
        or any(character not in "0123456789abcdef" for character in raw_digest)
    ):
        raise ValueError("execution input.sha256 must be a lowercase SHA-256 digest")
    path = Path(raw_path)
    if not path.is_absolute():
        path = artifact.parent / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"execution input file is missing: {path}")
    observed = _sha256_file(path)
    if observed != raw_digest:
        raise ValueError(f"execution input changed since execution: {path}")
    return path


def _resolve_digest_file(
    value: object, *, name: str, base_directory: Path,
) -> Path:
    data = _mapping(value, name)
    if set(data) != {"path", "sha256"}:
        raise ValueError(f"{name} must contain exactly path and sha256")
    raw_path = data["path"]
    raw_digest = data["sha256"]
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise ValueError(f"{name}.path must be a non-empty string")
    if (
        not isinstance(raw_digest, str)
        or len(raw_digest) != 64
        or any(character not in "0123456789abcdef" for character in raw_digest)
    ):
        raise ValueError(f"{name}.sha256 must be a lowercase SHA-256 digest")
    path = Path(raw_path)
    if not path.is_absolute():
        path = base_directory / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{name}.path is missing: {path}")
    if _sha256_file(path) != raw_digest:
        raise ValueError(f"{name}.path changed since execution: {path}")
    return path


def _verify_batch_artifacts(wrapper: Mapping[str, Any], *, artifact: Path) -> None:
    manifest_path = _resolve_digest_file(
        {
            "path": wrapper.get("manifest_path"),
            "sha256": wrapper.get("manifest_sha256"),
        },
        name="manifest",
        base_directory=artifact.parent,
    )
    manifest = _mapping(
        strict_json_loads(manifest_path.read_text(encoding="utf-8-sig")),
        "batch manifest",
    )
    if manifest.get("schema") != "degali.field-batch.v3":
        raise ValueError("batch manifest schema is not degali.field-batch.v3")
    manifest_cases = manifest.get("cases")
    execution_cases = wrapper.get("cases")
    if (
        not isinstance(manifest_cases, list)
        or not isinstance(execution_cases, list)
        or not manifest_cases
    ):
        raise ValueError("batch execution and manifest cases must be arrays")
    manifest_summary = _mapping(manifest.get("summary"), "batch manifest summary")
    execution_summary = _mapping(wrapper.get("summary"), "batch execution summary")
    if dict(execution_summary) != dict(manifest_summary):
        raise ValueError("batch execution summary does not match the manifest summary")
    if any(not isinstance(item, Mapping) for item in manifest_cases):
        raise ValueError("batch manifest cases must be objects")
    _validate_batch_summary_against_cases(manifest_summary, manifest_cases)
    manifest_by_label = {
        item.get("label"): item for item in manifest_cases
        if isinstance(item, Mapping)
    }
    if len(manifest_by_label) != len(manifest_cases) or any(
        not isinstance(label, str) or not label for label in manifest_by_label
    ):
        raise ValueError("batch manifest case labels must be unique non-empty strings")
    execution_labels = [
        item.get("label") for item in execution_cases
        if isinstance(item, Mapping)
    ]
    if (
        len(execution_labels) != len(execution_cases)
        or any(not isinstance(label, str) or not label for label in execution_labels)
        or len(set(execution_labels)) != len(execution_labels)
    ):
        raise ValueError("batch execution case labels must be unique strings")
    for item in execution_cases:
        entry = _mapping(item, "batch execution case")
        label = entry.get("label")
        if not isinstance(label, str) or label not in manifest_by_label:
            raise ValueError("batch execution case is missing from the manifest")
        manifest_entry = _mapping(manifest_by_label[label], f"manifest case {label!r}")
        manifest_decision = _mapping(
            manifest_entry.get("operational_decision"),
            f"manifest case {label!r}.operational_decision",
        )
        _validate_operational_decision_record(
            manifest_decision,
            context=f"manifest case {label!r}.operational_decision",
        )
        expected_status = manifest_decision.get("status")
        expected_allowed = manifest_decision.get("screening_allowed")
        if not isinstance(expected_status, str) or not expected_status:
            raise ValueError(f"manifest case {label!r} has an invalid operational status")
        if not isinstance(expected_allowed, bool):
            raise ValueError(f"manifest case {label!r} has an invalid screening_allowed flag")
        if entry.get("operational_decision") != expected_status:
            raise ValueError(f"batch operational status mismatch for case {label!r}")
        if entry.get("screening_allowed") is not expected_allowed:
            raise ValueError(f"batch screening_allowed mismatch for case {label!r}")
        expected_gate_codes = manifest_decision.get("gate_codes")
        actual_gate_codes = entry.get("gate_codes")
        if expected_gate_codes is not None:
            if (
                not isinstance(expected_gate_codes, list)
                or any(
                    not isinstance(code, str)
                    or code not in FIELD_OPERATIONAL_GATE_CODES
                    for code in expected_gate_codes
                )
                or len(set(expected_gate_codes)) != len(expected_gate_codes)
            ):
                raise ValueError(
                    f"manifest case {label!r} has invalid operational gate_codes"
                )
            if actual_gate_codes != expected_gate_codes:
                raise ValueError(f"batch gate_codes mismatch for case {label!r}")
        elif actual_gate_codes is not None:
            raise ValueError(
                f"batch execution case {label!r} has gate_codes absent from the manifest"
            )
        report_file = entry.get("report_file")
        report_sha = entry.get("report_sha256")
        if not isinstance(report_file, str) or not report_file.strip():
            raise ValueError(f"batch report file is invalid for case {label!r}")
        if (
            not isinstance(report_sha, str)
            or len(report_sha) != 64
            or any(character not in "0123456789abcdef" for character in report_sha)
        ):
            raise ValueError(f"batch report digest is invalid for case {label!r}")
        if report_file != manifest_entry.get("report_file"):
            raise ValueError(f"batch report file mismatch for case {label!r}")
        if report_sha != manifest_entry.get("report_sha256"):
            raise ValueError(f"batch report digest mismatch for case {label!r}")
        report_path = (manifest_path.parent / str(report_file)).resolve()
        try:
            report_path.relative_to(manifest_path.parent.resolve())
        except ValueError as error:
            raise ValueError(f"batch report escapes the manifest directory: {report_path}") from error
        if not report_path.is_file():
            raise FileNotFoundError(f"batch report is missing: {report_path}")
        if _sha256_file(report_path) != report_sha:
            raise ValueError(f"batch report changed since execution: {report_path}")
        try:
            report = strict_json_loads(report_path.read_text(encoding="utf-8-sig"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError(f"batch report is not valid UTF-8 JSON: {report_path}") from error
        report_data = _mapping(report, f"batch report {label!r}")
        report_schema = report_data.get("schema")
        if report_schema not in _SUPPORTED_BATCH_REPORT_SCHEMAS:
            allowed = ", ".join(sorted(_SUPPORTED_BATCH_REPORT_SCHEMAS))
            raise ValueError(
                f"batch report {label!r} has unsupported schema {report_schema!r}; "
                f"expected one of {allowed}"
            )
        execution_kind = manifest_entry.get("execution_kind")
        expected_kind = (
            "operational_uncertainty_envelope"
            if report_schema == FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA
            else "single_field_screen"
        )
        if execution_kind != expected_kind:
            raise ValueError(f"batch report/execution kind does not match case {label!r}")
        _validate_report_transport_semantics(
            report_data,
            context=f"batch report {label!r}",
            decision=manifest_decision,
        )
        _validate_report_decision_applicability(
            report_data,
            manifest_decision,
            context=f"batch report {label!r}",
        )
    if set(manifest_by_label) != {
        item.get("label") for item in execution_cases if isinstance(item, Mapping)
    }:
        raise ValueError("batch execution case set does not match the manifest")


def _require_report(payload: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = payload.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"execution artifact is missing object {key!r}")
    return value


def _require_typed_report(
    payload: Mapping[str, Any], key: str, allowed_schemas: frozenset[str],
) -> Mapping[str, Any]:
    """Require a report object whose discriminator is in the supported set."""
    report = _require_report(payload, key)
    schema = report.get("schema")
    if schema not in allowed_schemas:
        allowed = ", ".join(sorted(allowed_schemas))
        raise ValueError(
            f"execution report {key!r} has unsupported schema {schema!r}; "
            f"expected one of {allowed}"
        )
    return report


def _report_number(value: object, name: str) -> float:
    """Read one finite JSON number used by a serialized transport diagnostic."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return value


def _validate_sensor_report_semantics(
    report: Mapping[str, Any], *, context: str,
) -> None:
    """Validate the compact sensor extrema emitted by a field report."""
    if report.get("schema") == FIELD_SCREENING_REPORT_SCHEMA:
        for key in ("sensor_result", "sensor_results"):
            if key not in report:
                raise ValueError(f"{context} is missing {key}")
    sensor = report.get("sensor_result")
    if sensor is not None:
        sensor_data = _mapping(sensor, f"{context}.sensor_result")
        required = {
            "time_start_s", "time_end_s", "maximum_true_mole_fraction",
            "maximum_indicated_mole_fraction", "final_indicated_mole_fraction",
            "ambient_air_density_kg_m3",
        }
        missing = sorted(required.difference(sensor_data))
        if missing:
            raise ValueError(
                f"{context}.sensor_result is missing {', '.join(missing)}"
            )
        start = _report_number(sensor_data["time_start_s"], f"{context}.sensor_result.time_start_s")
        end = _report_number(sensor_data["time_end_s"], f"{context}.sensor_result.time_end_s")
        if end < start:
            raise ValueError(f"{context}.sensor_result time bounds are reversed")
        true_peak = _report_number(
            sensor_data["maximum_true_mole_fraction"],
            f"{context}.sensor_result.maximum_true_mole_fraction",
        )
        if not 0.0 <= true_peak <= 1.0:
            raise ValueError(f"{context}.sensor_result true mole fraction is outside [0, 1]")
        _report_number(
            sensor_data["maximum_indicated_mole_fraction"],
            f"{context}.sensor_result.maximum_indicated_mole_fraction",
        )
        _report_number(
            sensor_data["final_indicated_mole_fraction"],
            f"{context}.sensor_result.final_indicated_mole_fraction",
        )
        density = _report_number(
            sensor_data["ambient_air_density_kg_m3"],
            f"{context}.sensor_result.ambient_air_density_kg_m3",
        )
        if density <= 0.0:
            raise ValueError(f"{context}.sensor_result ambient density must be positive")

    sensor_results = report.get("sensor_results")
    if sensor_results is None:
        return
    if not isinstance(sensor_results, list):
        raise ValueError(f"{context}.sensor_results must be an array")
    labels: set[str] = set()
    for index, item in enumerate(sensor_results):
        item_context = f"{context}.sensor_results[{index}]"
        data = _mapping(item, item_context)
        label = data.get("label")
        if not isinstance(label, str) or not label.strip() or label in labels:
            raise ValueError(f"{context}.sensor_results labels must be unique non-empty strings")
        labels.add(label)
        position = data.get("position_m")
        if (
            not isinstance(position, list)
            or len(position) != 3
            or any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
                for value in position
            )
        ):
            raise ValueError(f"{item_context}.position_m must contain three finite numbers")
        withheld = data.get("withheld")
        if not isinstance(withheld, bool):
            raise ValueError(f"{item_context}.withheld must be boolean")
        warnings = data.get("warnings")
        if not isinstance(warnings, list) or any(
            not isinstance(value, str) or not value.strip() for value in warnings
        ):
            raise ValueError(f"{item_context}.warnings must be an array of non-empty strings")
        result = data.get("result")
        if withheld:
            if result is not None or not warnings:
                raise ValueError(f"{item_context} withheld sensor must retain a warning and no result")
            continue
        result_data = _mapping(result, f"{item_context}.result")
        result_required = {
            "time_start_s", "time_end_s", "maximum_true_mole_fraction",
            "maximum_indicated_mole_fraction", "final_indicated_mole_fraction",
        }
        missing = sorted(result_required.difference(result_data))
        if missing:
            raise ValueError(f"{item_context}.result is missing {', '.join(missing)}")
        start = _report_number(result_data["time_start_s"], f"{item_context}.result.time_start_s")
        end = _report_number(result_data["time_end_s"], f"{item_context}.result.time_end_s")
        if end < start:
            raise ValueError(f"{item_context}.result time bounds are reversed")
        true_peak = _report_number(
            result_data["maximum_true_mole_fraction"],
            f"{item_context}.result.maximum_true_mole_fraction",
        )
        if not 0.0 <= true_peak <= 1.0:
            raise ValueError(f"{item_context}.result true mole fraction is outside [0, 1]")
        _report_number(
            result_data["maximum_indicated_mole_fraction"],
            f"{item_context}.result.maximum_indicated_mole_fraction",
        )
        _report_number(
            result_data["final_indicated_mole_fraction"],
            f"{item_context}.result.final_indicated_mole_fraction",
        )


def _validate_transport_source_provenance(
    report: Mapping[str, Any],
    diagnostics: Mapping[str, Any],
    *,
    context: str,
    tolerance: float,
) -> None:
    """Tie serialized source ledgers back to the declared transport inputs.

    A report can be verified without replaying the full solver when an input
    option is intentionally omitted from the execution record.  The transport
    input still contains immutable schedule provenance, however, so the
    source-level injected ledger must agree with every declared schedule mass.
    This catches a report that keeps a small global residual while changing a
    single source mass or source label.
    """
    transport_input = report.get("transport_input")
    if transport_input is None:
        return
    data = _mapping(transport_input, f"{context}.transport_input")
    current_schema = report.get("schema") == FIELD_SCREENING_REPORT_SCHEMA
    if current_schema and "declared_direct_vapour_schedule" not in data:
        raise ValueError(
            f"{context}.transport_input is missing declared_direct_vapour_schedule"
        )
    if current_schema and "distributed_vapour_sources" not in data:
        raise ValueError(
            f"{context}.transport_input is missing distributed_vapour_sources"
        )
    declared = data.get("declared_direct_vapour_schedule")
    known_masses: dict[str, float] = {}
    if declared is not None:
        direct = _mapping(
            declared,
            f"{context}.transport_input.declared_direct_vapour_schedule",
        )
        known_masses["primary"] = _report_number(
            direct.get("released_mass_kg"),
            f"{context}.transport_input.declared_direct_vapour_schedule.released_mass_kg",
        )
        if known_masses["primary"] < 0.0:
            raise ValueError(
                f"{context}.transport_input direct-vapour released mass must be non-negative"
            )

    distributed = data.get("distributed_vapour_sources")
    if distributed is None:
        return
    if not isinstance(distributed, list):
        raise ValueError(f"{context}.transport_input.distributed_vapour_sources must be an array")
    for index, item in enumerate(distributed):
        item_context = f"{context}.transport_input.distributed_vapour_sources[{index}]"
        source = _mapping(item, item_context)
        label = source.get("label")
        if not isinstance(label, str) or not label.strip():
            raise ValueError(f"{item_context}.label must be a non-empty string")
        ledger_label = f"distributed:{label}"
        if ledger_label in known_masses:
            raise ValueError(f"{context} transport source labels are not unique")
        released_mass = _report_number(
            source.get("released_mass_kg"),
            f"{item_context}.released_mass_kg",
        )
        if released_mass < 0.0:
            raise ValueError(f"{item_context}.released_mass_kg must be non-negative")
        known_masses[ledger_label] = released_mass

    injected = diagnostics.get("source_mass_injected_kg")
    if injected is None or not known_masses:
        return
    if not isinstance(injected, list):
        raise ValueError(f"{context}.transport_result.diagnostics.source_mass_injected_kg must be an array")
    observed: dict[str, float] = {}
    for index, item in enumerate(injected):
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError(
                f"{context}.transport_result.diagnostics.source_mass_injected_kg[{index}] must be [label, value]"
            )
        label = item[0]
        if not isinstance(label, str) or not label.strip() or label in observed:
            raise ValueError(
                f"{context}.transport_result.diagnostics.source_mass_injected_kg labels must be unique non-empty strings"
            )
        observed[label] = _report_number(
            item[1],
            f"{context}.transport_result.diagnostics.source_mass_injected_kg[{index}][1]",
        )
        if observed[label] < 0.0:
            raise ValueError(
                f"{context}.transport_result.diagnostics.source_mass_injected_kg values must be non-negative"
            )
    expected_labels = {"primary", *known_masses}
    if set(observed) != expected_labels:
        raise ValueError(
            f"{context} transport source labels do not match declared transport inputs"
        )
    for label, expected in known_masses.items():
        if abs(observed[label] - expected) > tolerance:
            raise ValueError(
                f"{context} injected mass for {label!r} does not match its declared schedule mass"
            )


def _validate_transport_report_semantics(
    report: Mapping[str, Any], *, context: str,
    decision: Mapping[str, Any] | None = None,
) -> None:
    """Fail closed on serialized semi-FV conservation diagnostics.

    Complex field envelopes are intentionally not numerically replayed by the
    execution verifier.  Their nested field reports must nevertheless retain
    the same semantic relationship enforced during execution: residuals are
    finite, the per-source schedule residual maximum is honest, and any
    residual beyond the numerical gate is accompanied by a blocked
    applicability reason.  This prevents an integrity-only verification from
    turning a tampered transport diagnostic into an operationally usable
    artifact.
    """
    if (
        report.get("schema") == FIELD_SCREENING_REPORT_SCHEMA
        and "transport_result" not in report
    ):
        raise ValueError(f"{context} is missing transport_result")
    transport = report.get("transport_result")
    if transport is None:
        return
    transport_data = _mapping(transport, f"{context}.transport_result")
    maximum_concentration = _report_number(
        transport_data.get("maximum_concentration_kg_m3"),
        f"{context}.transport_result.maximum_concentration_kg_m3",
    )
    if maximum_concentration < 0.0:
        raise ValueError(f"{context}.transport_result maximum concentration must be non-negative")
    receptor_count = transport_data.get("receptor_count")
    if (
        isinstance(receptor_count, bool)
        or not isinstance(receptor_count, int)
        or receptor_count < 0
    ):
        raise ValueError(f"{context}.transport_result.receptor_count must be a non-negative integer")
    _validate_sensor_report_semantics(report, context=context)
    diagnostics = _mapping(
        transport_data.get("diagnostics"),
        f"{context}.transport_result.diagnostics",
    )

    base_gate_keys = {
        "mass_injected_kg",
        "mass_domain_kg",
        "mass_outflow_kg",
        "maximum_mass_residual_kg",
        "source_mass_ledger_residual_kg",
    }
    if report.get("schema") == FIELD_SCREENING_REPORT_SCHEMA:
        missing = sorted(base_gate_keys.difference(diagnostics))
        if missing:
            raise ValueError(
                f"{context}.transport_result.diagnostics is missing {', '.join(missing)}"
            )

    gate_keys = {
        "maximum_mass_residual_kg",
        "final_mass_residual_kg",
        "source_mass_ledger_residual_kg",
        "source_mass_schedule_residual_kg",
        "maximum_source_mass_schedule_residual_kg",
    }
    if not gate_keys.intersection(diagnostics):
        # Reports written before the typed conservation fields were added are
        # still structurally supported; there is no semantic gate to check.
        return

    mass_keys = ("mass_injected_kg", "mass_domain_kg", "mass_outflow_kg")
    masses = {
        key: _report_number(diagnostics.get(key), f"{context}.transport_result.diagnostics.{key}")
        for key in mass_keys
    }
    negative = sorted(key for key, value in masses.items() if value < 0.0)
    if negative:
        raise ValueError(
            f"{context} transport mass totals must be non-negative: "
            + ", ".join(negative)
        )
    scale = max(
        abs(masses["mass_injected_kg"]),
        abs(masses["mass_domain_kg"]) + abs(masses["mass_outflow_kg"]),
        _FIELD_MASS_RESIDUAL_ABS_TOL_KG,
    )
    tolerance = max(
        _FIELD_MASS_RESIDUAL_ABS_TOL_KG,
        _FIELD_MASS_RESIDUAL_REL_TOL * scale,
    )
    final_inventory_residual = (
        masses["mass_injected_kg"]
        - masses["mass_domain_kg"]
        - masses["mass_outflow_kg"]
    )
    if abs(final_inventory_residual) > tolerance:
        raise ValueError(
            f"{context} transport mass totals do not close: "
            f"final inventory residual {final_inventory_residual:.6g} kg exceeds "
            f"the numerical gate {tolerance:.6g} kg"
        )
    reported_final_residual: float | None = None
    if "final_mass_residual_kg" in diagnostics:
        reported_final_residual = _report_number(
            diagnostics["final_mass_residual_kg"],
            f"{context}.transport_result.diagnostics.final_mass_residual_kg",
        )
        if abs(reported_final_residual - final_inventory_residual) > tolerance:
            raise ValueError(
                f"{context} final mass residual does not match mass totals"
            )
    _validate_transport_source_provenance(
        report,
        diagnostics,
        context=context,
        tolerance=tolerance,
    )

    applicability = _mapping(report.get("applicability"), f"{context}.applicability")
    status = applicability.get("status")
    reasons = applicability.get("reasons")
    if not isinstance(status, str) or status not in {"accepted", "conditional", "blocked"}:
        raise ValueError(f"{context}.applicability.status is invalid")
    if not isinstance(reasons, list) or any(not isinstance(reason, str) for reason in reasons):
        raise ValueError(f"{context}.applicability.reasons must be a string array")
    required_gate_codes: set[str] = set()

    def require_blocked(
        residual: float,
        reason_fragment: str,
        label: str,
        gate_code: str,
    ) -> None:
        if residual <= tolerance:
            return
        if status != "blocked" or not any(
            reason_fragment in reason.lower() for reason in reasons
        ):
            raise ValueError(
                f"{context} {label} exceeds the numerical gate without a matching "
                "blocked applicability reason"
            )
        required_gate_codes.add(gate_code)

    def validate_decision_gate_codes() -> None:
        if decision is None or not required_gate_codes:
            return
        decision_status = decision.get("status")
        if decision_status != "withheld" or decision.get("screening_allowed") is not False:
            raise ValueError(
                f"{context} transport gate failure is not reflected in the operational decision"
            )
        decision_codes = decision.get("gate_codes")
        if (
            not isinstance(decision_codes, list)
            or any(not isinstance(code, str) for code in decision_codes)
            or not required_gate_codes.issubset(decision_codes)
        ):
            raise ValueError(
                f"{context} operational decision is missing transport gate_codes"
            )

    if "maximum_mass_residual_kg" in diagnostics:
        residual = _report_number(
            diagnostics["maximum_mass_residual_kg"],
            f"{context}.transport_result.diagnostics.maximum_mass_residual_kg",
        )
        if residual < 0.0:
            raise ValueError(f"{context} maximum mass residual must be non-negative")
        require_blocked(
            residual,
            "mass-conservation residual",
            "maximum mass residual",
            "transport_mass_residual_exceeded",
        )
    if reported_final_residual is not None:
        require_blocked(
            abs(reported_final_residual),
            "final mass-conservation residual",
            "final mass residual",
            "transport_mass_residual_exceeded",
        )

    injected_labels: set[str] | None = None
    if "source_mass_ledger_residual_kg" in diagnostics:
        ledger_residual = _report_number(
            diagnostics["source_mass_ledger_residual_kg"],
            f"{context}.transport_result.diagnostics.source_mass_ledger_residual_kg",
        )
        residual = abs(ledger_residual)
        injected = diagnostics.get("source_mass_injected_kg")
        injected_total: float | None = None
        if injected is not None:
            if not isinstance(injected, list):
                raise ValueError(
                    f"{context}.source_mass_injected_kg must be an array"
                )
            injected_labels = set()
            injected_values: list[float] = []
            for index, item in enumerate(injected):
                if not isinstance(item, list) or len(item) != 2:
                    raise ValueError(
                        f"{context}.source_mass_injected_kg[{index}] must be [label, value]"
                    )
                label = item[0]
                if not isinstance(label, str) or not label.strip() or label in injected_labels:
                    raise ValueError(
                        f"{context}.source_mass_injected_kg labels must be unique non-empty strings"
                    )
                injected_labels.add(label)
                value = _report_number(
                    item[1],
                    f"{context}.source_mass_injected_kg[{index}][1]",
                )
                if value < 0.0:
                    raise ValueError(
                        f"{context}.source_mass_injected_kg values must be non-negative"
                    )
                injected_values.append(value)
            injected_total = math.fsum(injected_values)
            if abs(masses["mass_injected_kg"] - injected_total - ledger_residual) > tolerance:
                raise ValueError(
                    f"{context} source mass ledger residual does not match injected source values"
                )
        require_blocked(
            residual,
            "source-mass ledger residual",
            "source-mass ledger residual",
            "source_ledger_mismatch",
        )

    schedule_residuals = diagnostics.get("source_mass_schedule_residual_kg")
    schedule_max = diagnostics.get("maximum_source_mass_schedule_residual_kg")
    if (schedule_residuals is None) != (schedule_max is None):
        raise ValueError(
            f"{context} source schedule residual map and maximum must be present together"
        )
    if schedule_residuals is None:
        validate_decision_gate_codes()
        return
    if not isinstance(schedule_residuals, list):
        raise ValueError(f"{context} source schedule residuals must be an array")
    if injected_labels is not None:
        schedule_labels = {
            item[0] for item in schedule_residuals
            if isinstance(item, list) and len(item) == 2 and isinstance(item[0], str)
        }
        if schedule_labels != injected_labels:
            raise ValueError(
                f"{context} source schedule and injected-mass labels do not match"
            )
    if injected_labels is None:
        raise ValueError(
            f"{context} source schedule residuals require source_mass_injected_kg"
        )
    parsed_residuals: list[tuple[str, float]] = []
    labels: set[str] = set()
    for index, item in enumerate(schedule_residuals):
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError(f"{context} source schedule residual {index} must be [label, value]")
        label = item[0]
        if not isinstance(label, str) or not label.strip() or label in labels:
            raise ValueError(f"{context} source schedule residual labels must be unique strings")
        labels.add(label)
        parsed_residuals.append((label, _report_number(
            item[1], f"{context}.source_mass_schedule_residual_kg[{index}][1]",
        )))
    maximum = _report_number(
        schedule_max,
        f"{context}.transport_result.diagnostics.maximum_source_mass_schedule_residual_kg",
    )
    if maximum < 0.0:
        raise ValueError(f"{context} maximum source schedule residual must be non-negative")
    observed_maximum = max((abs(value) for _label, value in parsed_residuals), default=0.0)
    comparison_tolerance = max(1.0e-12, 1.0e-12 * max(1.0, observed_maximum, maximum))
    if abs(maximum - observed_maximum) > comparison_tolerance:
        raise ValueError(
            f"{context} maximum source schedule residual does not match its per-source values"
        )
    require_blocked(
        maximum,
        "source-mass schedule residual",
        "maximum source schedule residual",
        "source_schedule_mass_mismatch",
    )

    validate_decision_gate_codes()


def _validate_report_transport_semantics(
    value: object, *, context: str,
    decision: Mapping[str, Any] | None = None,
) -> None:
    """Validate every nested ``field_screening`` transport report."""
    if isinstance(value, Mapping):
        if value.get("schema") == FIELD_SCREENING_REPORT_SCHEMA or "transport_result" in value:
            _validate_transport_report_semantics(
                value,
                context=context,
                decision=decision,
            )
        for key, child in value.items():
            _validate_report_transport_semantics(
                child,
                context=f"{context}.{key}",
                decision=decision,
            )
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _validate_report_transport_semantics(
                child,
                context=f"{context}[{index}]",
                decision=decision,
            )


def _validate_pressure_driven_source_report(
    value: Mapping[str, Any], *, scenario: object, context: str,
    allow_corner: bool = False,
) -> None:
    """Bind a typed pressure-derived source to an integrity-only report.

    History/phase/sensor envelope reports cannot always be replayed from the
    saved execution options. A direct pressure-derived source is still
    replayable at the strict case boundary, so its typed ambient boundary,
    uncertainty label, and derivation metadata must not be editable in the
    stored report without detection.
    """
    source = getattr(scenario, "source", None)
    boundary = getattr(source, "pressure_driven_mass_flow", None)
    if boundary is None:
        return
    report_scenario = _mapping(value.get("scenario"), f"{context}.scenario")
    report_source = _mapping(
        report_scenario.get("source"), f"{context}.scenario.source"
    )
    expected_source = source.as_dict()
    expected_boundary = expected_source.get("pressure_driven_mass_flow")
    actual_boundary = _mapping(
        report_source.get("pressure_driven_mass_flow"),
        f"{context}.scenario.source.pressure_driven_mass_flow",
    )
    if not allow_corner:
        if actual_boundary != expected_boundary:
            raise ValueError(
                f"{context} pressure-driven source boundary does not match strict input"
            )
    else:
        expected_boundary_data = _mapping(
            expected_boundary, f"{context}.expected pressure-driven boundary"
        )
        if (
            actual_boundary.get("source_id") != expected_boundary_data.get("source_id")
            or actual_boundary.get("allow_supercritical_gas")
            != expected_boundary_data.get("allow_supercritical_gas")
        ):
            raise ValueError(
                f"{context} pressure-driven source identity does not match strict input"
            )
        expected_ambient = _mapping(
            expected_boundary_data.get("ambient_pressure_pa"),
            f"{context}.expected pressure-driven ambient pressure",
        )
        actual_ambient = _mapping(
            actual_boundary.get("ambient_pressure_pa"),
            f"{context}.scenario.source.pressure_driven_mass_flow.ambient_pressure_pa",
        )
        if (
            actual_ambient.get("unit") != expected_ambient.get("unit")
            or actual_ambient.get("source") != expected_ambient.get("source")
        ):
            raise ValueError(
                f"{context} pressure-driven ambient provenance does not match strict input"
            )
        expected_lower = _report_number(
            expected_ambient.get("lower"),
            f"{context}.expected pressure-driven ambient lower",
        )
        expected_upper = _report_number(
            expected_ambient.get("upper"),
            f"{context}.expected pressure-driven ambient upper",
        )
        actual_lower = _report_number(
            actual_ambient.get("lower"),
            f"{context}.pressure-driven ambient lower",
        )
        actual_nominal = _report_number(
            actual_ambient.get("nominal"),
            f"{context}.pressure-driven ambient nominal",
        )
        actual_upper = _report_number(
            actual_ambient.get("upper"),
            f"{context}.pressure-driven ambient upper",
        )
        tolerance = max(1.0e-9, 1.0e-12 * max(1.0, expected_upper))
        if (
            actual_lower > actual_nominal + tolerance
            or actual_nominal > actual_upper + tolerance
            or actual_lower < expected_lower - tolerance
            or actual_upper > expected_upper + tolerance
        ):
            raise ValueError(
                f"{context} pressure-driven ambient corner escapes strict input bounds"
            )
    expected_uncertainty = source.uncertainty_fields().get(
        "pressure_driven_ambient_pressure_pa"
    )
    report_uncertainty = _mapping(
        report_source.get("uncertainty"),
        f"{context}.scenario.source.uncertainty",
    )
    actual_uncertainty = _mapping(
        report_uncertainty.get("pressure_driven_ambient_pressure_pa"),
        f"{context}.scenario.source.uncertainty.pressure_driven_ambient_pressure_pa",
    )
    expected_uncertainty_record = (
        None if expected_uncertainty is None else expected_uncertainty.as_dict()
    )
    if expected_uncertainty_record is None:
        raise ValueError(
            f"{context} pressure-driven ambient uncertainty is missing from strict input"
        )
    if not allow_corner and actual_uncertainty != expected_uncertainty_record:
        raise ValueError(
            f"{context} pressure-driven ambient uncertainty does not match strict input"
        )
    if allow_corner:
        expected_lower = float(expected_uncertainty_record["lower"])
        expected_upper = float(expected_uncertainty_record["upper"])
        actual_lower = _report_number(
            actual_uncertainty.get("lower"),
            f"{context}.pressure-driven uncertainty lower",
        )
        actual_nominal = _report_number(
            actual_uncertainty.get("nominal"),
            f"{context}.pressure-driven uncertainty nominal",
        )
        actual_upper = _report_number(
            actual_uncertainty.get("upper"),
            f"{context}.pressure-driven uncertainty upper",
        )
        if (
            actual_uncertainty.get("unit") != expected_uncertainty_record.get("unit")
            or actual_uncertainty.get("source") != expected_uncertainty_record.get("source")
            or actual_lower < expected_lower
            or actual_nominal < actual_lower
            or actual_nominal > actual_upper
            or actual_upper > expected_upper
        ):
            raise ValueError(
                f"{context} pressure-driven ambient uncertainty escapes strict input bounds"
            )
    expected_metadata = _mapping(
        expected_source.get("metadata"), f"{context}.scenario.source.metadata"
    )
    actual_metadata = _mapping(
        report_source.get("metadata"), f"{context}.scenario.source.metadata"
    )
    if not allow_corner:
        for key, expected_value in expected_metadata.items():
            if isinstance(key, str) and key.startswith("mass_flow_derivation_"):
                if actual_metadata.get(key) != expected_value:
                    raise ValueError(
                        f"{context} pressure-driven derivation metadata does not match strict input"
                    )
    else:
        if (
            actual_metadata.get("mass_flow_derivation")
            != expected_metadata.get("mass_flow_derivation")
            or actual_metadata.get("mass_flow_derivation_source_id")
            != expected_boundary.get("source_id")
            or actual_metadata.get("mass_flow_derivation_allow_supercritical_gas")
            != str(expected_boundary.get("allow_supercritical_gas")).lower()
        ):
            raise ValueError(
                f"{context} pressure-driven derivation identity does not match strict input"
            )
        try:
            metadata_nominal = float(
                actual_metadata["mass_flow_derivation_ambient_pressure_pa"]
            )
            metadata_bounds = tuple(
                float(item)
                for item in str(
                    actual_metadata["mass_flow_derivation_ambient_pressure_bounds_pa"]
                ).split(",")
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f"{context} pressure-driven derivation ambient metadata is invalid"
            ) from error
        if len(metadata_bounds) != 2 or not all(math.isfinite(item) for item in metadata_bounds):
            raise ValueError(
                f"{context} pressure-driven derivation ambient metadata is invalid"
            )
        actual_ambient = _mapping(
            actual_boundary["ambient_pressure_pa"],
            f"{context}.pressure-driven ambient pressure",
        )
        boundary_nominal = _report_number(
            actual_ambient["nominal"], f"{context}.pressure-driven ambient nominal"
        )
        boundary_lower = _report_number(
            actual_ambient["lower"], f"{context}.pressure-driven ambient lower"
        )
        boundary_upper = _report_number(
            actual_ambient["upper"], f"{context}.pressure-driven ambient upper"
        )
        expected_lower = float(expected_boundary["ambient_pressure_pa"]["lower"])
        expected_upper = float(expected_boundary["ambient_pressure_pa"]["upper"])
        if (
            not math.isclose(metadata_nominal, boundary_nominal, rel_tol=1.0e-12, abs_tol=1.0e-9)
            or metadata_bounds[0] > metadata_bounds[1]
            or metadata_bounds[0] < expected_lower
            or metadata_bounds[1] > expected_upper
            or metadata_bounds[0] > boundary_lower
            or metadata_bounds[1] < boundary_upper
        ):
            raise ValueError(
                f"{context} pressure-driven derivation ambient metadata escapes strict input bounds"
            )


def _validate_pressure_driven_source_reports(
    value: object, *, scenario: object, context: str,
) -> None:
    """Validate typed pressure-derived source records in nested reports."""
    if isinstance(value, Mapping):
        if isinstance(value.get("scenario"), Mapping):
            _validate_pressure_driven_source_report(
                value, scenario=scenario, context=context, allow_corner=True,
            )
        for key, child in value.items():
            _validate_pressure_driven_source_reports(
                child, scenario=scenario, context=f"{context}.{key}",
            )
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _validate_pressure_driven_source_reports(
                child, scenario=scenario, context=f"{context}[{index}]",
            )


def _screen_option(data: Mapping[str, Any], key: str, expected: type) -> Any:
    value = data.get(key)
    if expected is bool:
        if not isinstance(value, bool):
            raise ValueError(f"execution input.{key} must be boolean")
    elif expected is int:
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"execution input.{key} must be a positive integer")
    elif expected is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"execution input.{key} must be finite")
        value = float(value)
    return value


def _recompute_supported_screening(
    wrapper: Mapping[str, Any], *, input_path: Path,
) -> tuple[dict[str, object], dict[str, object]] | None:
    """Recompute only screening modes whose numerical options are recorded."""
    input_data = _mapping(wrapper["input"], "execution input")
    required_options = {
        "refinement_requested", "uncertainty_envelope_requested",
        "sensor_array_envelope_requested", "joint_source_sensor_envelope_requested",
        "allow_conditional", "refinement_factors", "relative_tolerance",
        "max_cell_steps", "max_uncertainty_cases",
    }
    if not required_options.issubset(input_data):
        return None
    parsed = read_field_screening_case_json(input_path)
    refinement_requested = _screen_option(input_data, "refinement_requested", bool)
    uncertainty_requested = _screen_option(
        input_data, "uncertainty_envelope_requested", bool,
    )
    sensor_requested = _screen_option(input_data, "sensor_array_envelope_requested", bool)
    joint_requested = _screen_option(
        input_data, "joint_source_sensor_envelope_requested", bool,
    )
    allow_conditional = _screen_option(input_data, "allow_conditional", bool)
    kind = input_data.get("uncertainty_envelope_kind")
    has_phase_input = parsed.phase_routing_transport is not None
    has_history_input = parsed.imported_history is not None
    if has_phase_input and (
        not uncertainty_requested or kind != "phase_routing_transport"
    ):
        raise ValueError(
            "phase-routing case execution must declare its complete phase envelope"
        )
    if kind == "phase_routing_transport" and not has_phase_input:
        raise ValueError(
            "execution declares phase-routing transport but the case has none"
        )
    history_kinds = {"joint_measured_history", "joint_pressure_driven_history"}
    if has_history_input and uncertainty_requested and kind not in history_kinds:
        raise ValueError(
            "measured-history envelope execution must retain its joint-history kind"
        )
    if kind in history_kinds and not has_history_input:
        raise ValueError(
            "execution declares joint measured history but the case has none"
        )
    if (
        kind == "joint_pressure_driven_history"
        and has_history_input
        and not isinstance(parsed.imported_history, ImportedPressureDrivenMeasuredHistory)
    ):
        raise ValueError(
            "execution declares joint pressure-driven history but the case history is measured-flow"
        )
    if (
        kind == "joint_measured_history"
        and isinstance(parsed.imported_history, ImportedPressureDrivenMeasuredHistory)
    ):
        raise ValueError(
            "execution declares joint measured history but the case history is pressure-driven"
        )
    source_case = input_data.get("atmospheric_source_schedule_case")
    if sensor_requested and joint_requested:
        raise ValueError(
            "execution cannot request both sensor-array and source-sensor envelopes"
        )
    if (sensor_requested or joint_requested) and uncertainty_requested:
        raise ValueError(
            "sensor-array and source-sensor envelope executions cannot retain "
            "the standalone uncertainty flag"
        )
    if sensor_requested and kind != "sensor_calibration":
        raise ValueError(
            "sensor-array envelope execution must retain sensor_calibration kind"
        )
    if kind == "sensor_calibration" and not sensor_requested:
        raise ValueError(
            "execution declares sensor calibration without requesting its envelope"
        )
    if joint_requested and kind != "source_sensor":
        raise ValueError(
            "source-sensor envelope execution must retain source_sensor kind"
        )
    if kind == "source_sensor" and not joint_requested:
        raise ValueError(
            "execution declares source-sensor uncertainty without requesting its envelope"
        )
    if sensor_requested and source_case is not None:
        raise ValueError(
            "sensor-array envelope cannot retain an atmospheric source schedule"
        )
    if joint_requested and has_phase_input:
        raise ValueError(
            "source-sensor envelope cannot retain phase-routing input"
        )
    if joint_requested and has_history_input:
        raise ValueError(
            "source-sensor envelope cannot retain measured-history input"
        )
    refinement_factors = input_data.get("refinement_factors")
    if (
        not isinstance(refinement_factors, list)
        or not refinement_factors
        or any(isinstance(value, bool) or not isinstance(value, int) or value <= 0 for value in refinement_factors)
    ):
        raise ValueError("execution input.refinement_factors must be a non-empty positive-integer array")
    relative_tolerance = _screen_option(input_data, "relative_tolerance", float)
    max_cell_steps = _screen_option(input_data, "max_cell_steps", int)
    max_uncertainty_cases = _screen_option(input_data, "max_uncertainty_cases", int)
    if sensor_requested:
        result = run_field_operational_sensor_array_envelope(
            parsed.request,
            max_cases=max_uncertainty_cases,
            include_refinement=refinement_requested,
            refinement_factors=tuple(refinement_factors),
            relative_tolerance=relative_tolerance,
            max_cell_steps=max_cell_steps,
            allow_conditional=allow_conditional,
        )
        return (
            field_operational_sensor_array_envelope_report(result),
            field_operational_screening_decision_record(result.operational_decision),
        )
    if joint_requested:
        atmospheric_source_schedule = None
        if source_case is not None:
            source_data = _mapping(source_case, "atmospheric source case")
            source_path = _resolve_digest_file(
                {"path": source_data.get("path"), "sha256": source_data.get("sha256")},
                name="atmospheric source case",
                base_directory=input_path.parent,
            )
            atmospheric_source_schedule = read_field_source_schedule_json(source_path)
        result = run_field_operational_source_sensor_envelope(
            parsed.request,
            atmospheric_source_schedule=atmospheric_source_schedule,
            max_cases=max_uncertainty_cases,
            table_nodes=161,
            include_refinement=refinement_requested,
            refinement_factors=tuple(refinement_factors),
            relative_tolerance=relative_tolerance,
            max_cell_steps=max_cell_steps,
            allow_conditional=allow_conditional,
        )
        return (
            field_operational_source_sensor_envelope_report(result),
            field_operational_screening_decision_record(result.operational_decision),
        )
    if (
        uncertainty_requested
        and kind == "phase_routing_transport"
        and parsed.phase_routing_transport is not None
    ):
        phase_input = parsed.phase_routing_transport
        result = run_field_operational_phase_routing_transport_envelope(
            parsed.request,
            phase_input.config,
            phase_input.pool_launch,
            pool_vertical_sigma_m=phase_input.pool_vertical_sigma_m,
            pool_vertical_sigma_uncertainty=phase_input.pool_vertical_sigma_uncertainty,
            phase_uncertainty=phase_input.phase_uncertainty,
            max_cases=max_uncertainty_cases,
            table_nodes=phase_input.table_nodes,
            include_refinement=refinement_requested,
            refinement_factors=tuple(refinement_factors),
            relative_tolerance=relative_tolerance,
            max_cell_steps=max_cell_steps,
            allow_conditional=allow_conditional,
        )
        return (
            field_operational_phase_routing_transport_envelope_report(result),
            field_operational_screening_decision_record(result.operational_decision),
        )
    if uncertainty_requested and kind in history_kinds and parsed.imported_history is not None:
        schedule = parsed.measured_schedule
        source_request = parsed.source_request
        if (
            schedule is None
            or schedule.quality_assessment is None
            or source_request is None
        ):
            raise ValueError(
                "measured-history case is missing its approved source-envelope evidence"
            )
        operational_kwargs = {
            "quality_criteria": schedule.quality_assessment.criteria,
            "history_provenance": parsed.request.measured_history_provenance,
            "max_cases": max_uncertainty_cases,
            "include_refinement": refinement_requested,
            "refinement_factors": tuple(refinement_factors),
            "relative_tolerance": relative_tolerance,
            "max_cell_steps": max_cell_steps,
            "allow_conditional": allow_conditional,
        }
        if (
            kind == "joint_pressure_driven_history"
            and isinstance(parsed.imported_history, ImportedPressureDrivenMeasuredHistory)
        ):
            result = run_field_operational_joint_pressure_driven_history_envelope(
                source_request, parsed.imported_history.history, **operational_kwargs,
            )
        else:
            result = run_field_operational_joint_measured_history_envelope(
                source_request, parsed.imported_history.history, **operational_kwargs,
            )
        return (
            field_operational_joint_measured_history_envelope_report(result),
            field_operational_screening_decision_record(result.operational_decision),
        )
    if uncertainty_requested and kind == "nominal_field" and source_case is None:
        result = run_field_operational_uncertainty_envelope(
            parsed.request,
            max_cases=max_uncertainty_cases,
            include_refinement=refinement_requested,
            refinement_factors=tuple(refinement_factors),
            relative_tolerance=relative_tolerance,
            max_cell_steps=max_cell_steps,
            allow_conditional=allow_conditional,
        )
        return (
            field_operational_uncertainty_envelope_report(result),
            field_operational_screening_decision_record(result.operational_decision),
        )
    if not uncertainty_requested and kind is None and source_case is None:
        result = (
            run_field_semi_fv_refinement_study(
                parsed.request,
                refinement_factors=tuple(refinement_factors),
                relative_tolerance=relative_tolerance,
                max_cell_steps=max_cell_steps,
            )
            if refinement_requested
            else run_field_semi_fv_screening(parsed.request)
        )
        decision = evaluate_field_operational_screening(
            result, require_refinement=True, allow_conditional=allow_conditional,
        )
        report = (
            field_refinement_report(result)
            if refinement_requested else field_screening_report(result)
        )
        return report, field_operational_screening_decision_record(decision)
    if (
        uncertainty_requested
        and kind == "atmospheric_source"
        and source_case is not None
    ):
        source_data = _mapping(source_case, "atmospheric source case")
        source_path = _resolve_digest_file(
            {"path": source_data.get("path"), "sha256": source_data.get("sha256")},
            name="atmospheric source case",
            base_directory=input_path.parent,
        )
        schedule = read_field_source_schedule_json(source_path)
        result = run_field_operational_uncertainty_envelope(
            parsed.request,
            atmospheric_source_schedule=schedule,
            max_cases=max_uncertainty_cases,
            include_refinement=refinement_requested,
            refinement_factors=tuple(refinement_factors),
            relative_tolerance=relative_tolerance,
            max_cell_steps=max_cell_steps,
            allow_conditional=allow_conditional,
        )
        return (
            field_operational_uncertainty_envelope_report(result),
            field_operational_screening_decision_record(result.operational_decision),
        )
    return None


def _validate_measured_history_execution_provenance(
    report: Mapping[str, Any],
    parsed_case: object,
    input_data: Mapping[str, Any],
) -> None:
    """Bind measured-history reports to the CSV currently imported.

    Even when a report is replayable, the historian CSV remains an independent
    source boundary.  A changed CSV must invalidate the saved report rather
    than merely produce a fresh in-memory provenance.
    """
    imported = getattr(parsed_case, "imported_history", None)
    declared = input_data.get("measured_history_imported")
    if declared is not None and not isinstance(declared, bool):
        raise ValueError("execution input.measured_history_imported must be boolean")
    if imported is None:
        if declared is True:
            raise ValueError(
                "execution input claims measured-history import but the case has none"
            )
        return
    if declared is not True:
        raise ValueError(
            "execution input must explicitly retain measured_history_imported=true"
        )
    expected = dict(getattr(parsed_case.request, "measured_history_provenance", ()))
    if not expected:
        raise ValueError("measured-history case has no serialized CSV provenance")
    expected_event = expected.get("event_id")
    if input_data.get("measured_history_event_id") != expected_event:
        raise ValueError(
            "execution measured_history_event_id does not match the current case"
        )

    matches = 0

    def visit(value: object, context: str) -> None:
        nonlocal matches
        if isinstance(value, Mapping):
            if "measured_history_provenance" in value:
                provenance = value["measured_history_provenance"]
                if not isinstance(provenance, Mapping):
                    raise ValueError(
                        f"{context}.measured_history_provenance must be an object"
                    )
                observed = dict(provenance)
                # The nested, source-history-only envelope predates the
                # operational wrapper and intentionally has no field-request
                # provenance.  Empty records are acceptable there, but every
                # populated field report must match the current CSV exactly.
                if observed not in ({}, expected):
                    raise ValueError(
                        f"{context}.measured_history_provenance does not match the current CSV"
                    )
                if observed == expected:
                    matches += 1
            for key, child in value.items():
                visit(child, f"{context}.{key}")
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, f"{context}[{index}]")

    visit(report, "field execution field_report")
    if matches == 0:
        raise ValueError(
            "measured-history execution report is missing CSV provenance"
        )


def verify_field_execution_artifact(path: str | Path) -> dict[str, object]:
    """Verify one supported field execution record without promoting it."""
    artifact = Path(path)
    if not artifact.is_file():
        raise FileNotFoundError(artifact)
    try:
        payload = strict_json_loads(artifact.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid field execution JSON at {artifact}: {error.msg}") from error
    wrapper = _mapping(payload, "field execution artifact")
    schema = wrapper.get("schema")
    if schema not in _SUPPORTED_EXECUTION_SCHEMAS:
        allowed = ", ".join(sorted(_SUPPORTED_EXECUTION_SCHEMAS))
        raise ValueError(f"unsupported field execution schema {schema!r}; expected one of {allowed}")
    input_path = _input_file(
        wrapper,
        artifact=artifact,
        allow_extra=schema in {
            "degali.field-screening-execution.v1",
            "degali.field-batch-execution.v1",
        },
    )
    if schema == "degali.field-atmospheric-schedule-execution.v1":
        schedule = read_field_source_schedule_json(input_path)
        expected = field_source_schedule_record(schedule)
        actual = _source_schedule_record_with_defaults(
            _require_report(wrapper, "atmospheric_source_schedule"),
            "atmospheric_source_schedule",
        )
        report_recomputed = True
    elif schema == "degali.field-model-comparison-execution.v1":
        case = read_field_model_comparison_case_json(input_path)
        comparison = compare_field_model_sensor_sets(
            case.left,
            case.right,
            threshold_mole_fraction=case.threshold_mole_fraction,
            position_tolerance_m=case.position_tolerance_m,
        )
        expected = field_model_comparison_report(
            comparison, comparison_evidence=case.comparison_evidence,
        )
        actual = _require_report(wrapper, "field_model_comparison")
        report_recomputed = True
    elif schema == "degali.field-validation-execution.v1":
        case = read_field_validation_case_json(input_path)
        score = score_field_model_against_validation(
            case.model,
            case.dataset,
            threshold_mole_fraction=case.threshold_mole_fraction,
            position_tolerance_m=case.position_tolerance_m,
        )
        expected = field_validation_score_record(score)
        actual = _require_report(wrapper, "field_validation_score")
        report_recomputed = True
    elif schema == "degali.field-screening-execution.v1":
        parsed_case = read_field_screening_case_json(input_path)
        input_data = _mapping(wrapper["input"], "execution input")
        _validate_measured_history_execution_provenance(
            wrapper.get("field_report"), parsed_case, input_data,
        )
        source_case = input_data.get("atmospheric_source_schedule_case")
        if source_case is not None:
            source_path = _resolve_digest_file(
                {
                    "path": _mapping(source_case, "atmospheric source case").get("path"),
                    "sha256": _mapping(source_case, "atmospheric source case").get("sha256"),
                },
                name="atmospheric source case",
                base_directory=artifact.parent,
            )
            source = read_field_source_schedule_json(source_path)
            stored_source = _mapping(source_case, "atmospheric source case").get("schedule")
            if _source_schedule_record_with_defaults(
                stored_source, "atmospheric source case.schedule"
            ) != field_source_schedule_record(source):
                raise ValueError("field screening source schedule does not match its case")
        recomputed = _recompute_supported_screening(
            wrapper, input_path=input_path,
        )
        if recomputed is None:
            actual_report = _require_typed_report(
                wrapper, "field_report", _SUPPORTED_FIELD_REPORT_SCHEMAS,
            )
            decision = _require_typed_report(
                wrapper,
                "operational_screening",
                frozenset({FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA}),
            )
            _validate_operational_decision_record(
                decision,
                context="field execution operational_screening",
            )
            _validate_report_decision_applicability(
                actual_report,
                decision,
                context="field execution field_report",
            )
            _validate_report_transport_semantics(
                actual_report,
                context="field execution field_report",
                decision=decision,
            )
            _validate_pressure_driven_source_reports(
                actual_report,
                scenario=parsed_case.request.scenario,
                context="field execution field_report",
            )
            actual = expected = {}
            report_recomputed = False
        else:
            expected_report, expected_decision = recomputed
            actual_report = _require_typed_report(
                wrapper, "field_report", _SUPPORTED_FIELD_REPORT_SCHEMAS,
            )
            actual_decision = _require_typed_report(
                wrapper,
                "operational_screening",
                frozenset({FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA}),
            )
            _validate_operational_decision_record(
                actual_decision,
                context="field execution operational_screening",
            )
            _validate_report_decision_applicability(
                actual_report,
                actual_decision,
                context="field execution field_report",
            )
            actual = {"field_report": actual_report, "operational_screening": actual_decision}
            expected = {
                "field_report": expected_report,
                "operational_screening": expected_decision,
            }
            _validate_report_transport_semantics(
                actual_report,
                context="field execution field_report",
                decision=actual_decision,
            )
            report_recomputed = True
    else:
        read_field_batch_input_json(input_path)
        _verify_batch_artifacts(wrapper, artifact=artifact)
        actual = expected = {}
        report_recomputed = False
    if report_recomputed and dict(actual) != expected:
        raise ValueError(
            "field execution report does not match a fresh deterministic recomputation"
        )
    return {
        "schema": FIELD_EXECUTION_VERIFICATION_SCHEMA,
        "execution_schema": schema,
        "artifact_path": str(artifact.resolve()),
        "artifact_sha256": _sha256_file(artifact),
        "input_path": str(input_path),
        "input_sha256_verified": True,
        "referenced_artifacts_verified": True,
        "report_recomputed": report_recomputed,
        "promotion_allowed": False,
    }


__all__ = [
    "FIELD_EXECUTION_VERIFICATION_SCHEMA",
    "verify_field_execution_artifact",
]
