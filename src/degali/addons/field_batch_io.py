"""Strict JSON boundary for named field-screening batches.

The Python batch API intentionally accepts already-parsed requests.  This
module provides the deployment-facing file boundary around it: each named
case points to a strict field case JSON, and an optional atmospheric source
case points to the separately fingerprinted source schedule JSON.  Historian
and phase-routing cases are rejected here because they require their own
operational envelope adapters and must not be silently mixed into a nominal
batch.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

from .field_case_io import FieldScreeningCase, read_field_screening_case_json
from .field_json import strict_json_loads
from .field_decision import FieldConditionalReviewAuthorization
from .field_source_io import (
    FieldAtmosphericSourceSchedule,
    read_field_source_schedule_json,
)
from .field_workflow import FieldSemiFVRequest


FIELD_BATCH_INPUT_SCHEMA = "degali.field-batch-input.v1"
_SAFE_CASE_LABEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")


def _object(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _number(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite") from error
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _integer(value: object, name: str, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer of at least {minimum}")
    return value


def _boolean(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be boolean")
    return value


def _keys(
    value: Mapping[str, Any], name: str, *, required: set[str], optional: set[str],
) -> None:
    missing = sorted(required - set(value))
    unknown = sorted(set(value) - required - optional)
    if missing or unknown:
        details = []
        if missing:
            details.append("missing=" + ", ".join(missing))
        if unknown:
            details.append("unknown=" + ", ".join(unknown))
        raise ValueError(f"{name} keys are invalid: " + "; ".join(details))


@dataclass(frozen=True)
class FieldBatchExecutionOptions:
    """Execution controls retained in the batch input record."""

    include_refinement: bool = True
    refinement_factors: tuple[int, ...] = (1, 2)
    relative_tolerance: float = 0.05
    max_cell_steps: int = 20_000_000
    allow_conditional_operational_screening: bool = False
    require_refinement_for_operational_screening: bool = True
    require_in_plane_sensor_for_operational_screening: bool = True
    require_resolved_uncertainty_for_operational_screening: bool = True
    max_uncertainty_cases: int = 64
    table_nodes: int = 161

    def __post_init__(self) -> None:
        booleans = (
            self.include_refinement,
            self.allow_conditional_operational_screening,
            self.require_refinement_for_operational_screening,
            self.require_in_plane_sensor_for_operational_screening,
            self.require_resolved_uncertainty_for_operational_screening,
        )
        if not all(isinstance(value, bool) for value in booleans):
            raise TypeError("field batch execution boolean options must be boolean")
        if (
            not isinstance(self.refinement_factors, tuple)
            or not self.refinement_factors
            or any(
                isinstance(value, bool) or not isinstance(value, int) or value <= 0
                for value in self.refinement_factors
            )
        ):
            raise ValueError("refinement_factors must be a non-empty tuple of positive integers")
        if not math.isfinite(float(self.relative_tolerance)) or self.relative_tolerance < 0.0:
            raise ValueError("relative_tolerance must be finite and non-negative")
        if isinstance(self.max_cell_steps, bool) or not isinstance(self.max_cell_steps, int) or self.max_cell_steps < 1:
            raise ValueError("max_cell_steps must be a positive integer")
        if isinstance(self.max_uncertainty_cases, bool) or not isinstance(self.max_uncertainty_cases, int) or self.max_uncertainty_cases < 1:
            raise ValueError("max_uncertainty_cases must be a positive integer")
        if isinstance(self.table_nodes, bool) or not isinstance(self.table_nodes, int) or self.table_nodes < 2:
            raise ValueError("table_nodes must be an integer of at least two")

    def as_record(self) -> dict[str, object]:
        return {
            "include_refinement": self.include_refinement,
            "refinement_factors": list(self.refinement_factors),
            "relative_tolerance": self.relative_tolerance,
            "max_cell_steps": self.max_cell_steps,
            "allow_conditional_operational_screening": self.allow_conditional_operational_screening,
            "require_refinement_for_operational_screening": self.require_refinement_for_operational_screening,
            "require_in_plane_sensor_for_operational_screening": self.require_in_plane_sensor_for_operational_screening,
            "require_resolved_uncertainty_for_operational_screening": self.require_resolved_uncertainty_for_operational_screening,
            "max_uncertainty_cases": self.max_uncertainty_cases,
            "table_nodes": self.table_nodes,
        }

    def as_kwargs(self) -> dict[str, object]:
        return {
            "include_refinement": self.include_refinement,
            "refinement_factors": self.refinement_factors,
            "relative_tolerance": self.relative_tolerance,
            "max_cell_steps": self.max_cell_steps,
            "allow_conditional_operational_screening": self.allow_conditional_operational_screening,
            "require_refinement_for_operational_screening": self.require_refinement_for_operational_screening,
            "require_in_plane_sensor_for_operational_screening": self.require_in_plane_sensor_for_operational_screening,
            "require_resolved_uncertainty_for_operational_screening": self.require_resolved_uncertainty_for_operational_screening,
            "max_uncertainty_cases": self.max_uncertainty_cases,
            "table_nodes": self.table_nodes,
        }


@dataclass(frozen=True)
class FieldBatchInputCase:
    """One parsed named case and its input-artifact provenance."""

    label: str
    case_path: Path
    case_sha256: str
    request: FieldSemiFVRequest
    conditional_review: FieldConditionalReviewAuthorization | None = None
    source_case_path: Path | None = None
    source_case_sha256: str | None = None
    atmospheric_source_schedule: FieldAtmosphericSourceSchedule | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or _SAFE_CASE_LABEL.fullmatch(self.label) is None:
            raise ValueError(
                "batch case label must match the safe report-name pattern"
            )
        if not isinstance(self.case_path, Path) or not self.case_path.is_file():
            raise ValueError("batch case_path must point to an existing file")
        if len(self.case_sha256) != 64 or any(char not in "0123456789abcdef" for char in self.case_sha256):
            raise ValueError("batch case_sha256 must be a lowercase SHA-256 digest")
        if hashlib.sha256(self.case_path.read_bytes()).hexdigest() != self.case_sha256:
            raise ValueError("batch case_sha256 does not match case_path")
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("batch case request must be FieldSemiFVRequest")
        if self.conditional_review is not None and not isinstance(
            self.conditional_review, FieldConditionalReviewAuthorization
        ):
            raise TypeError("conditional_review must be FieldConditionalReviewAuthorization or None")
        if (self.source_case_path is None) != (self.source_case_sha256 is None):
            raise ValueError("source case path and SHA-256 must be supplied together")
        if self.source_case_path is not None:
            if not self.source_case_path.is_file():
                raise ValueError("source_case_path must point to an existing file")
            if len(self.source_case_sha256) != 64 or any(char not in "0123456789abcdef" for char in self.source_case_sha256):
                raise ValueError("source_case_sha256 must be a lowercase SHA-256 digest")
            if hashlib.sha256(self.source_case_path.read_bytes()).hexdigest() != self.source_case_sha256:
                raise ValueError("source_case_sha256 does not match source_case_path")
            if not isinstance(self.atmospheric_source_schedule, FieldAtmosphericSourceSchedule):
                raise TypeError("source case must produce FieldAtmosphericSourceSchedule")
        elif self.atmospheric_source_schedule is not None:
            raise ValueError("atmospheric source schedule requires source case provenance")

    def as_record(self) -> dict[str, object]:
        source = None
        if self.source_case_path is not None:
            source = {
                "path": str(self.source_case_path.resolve()),
                "sha256": self.source_case_sha256,
                "schedule": self.atmospheric_source_schedule.as_record(),
            }
        return {
            "label": self.label,
            "case_path": str(self.case_path.resolve()),
            "case_sha256": self.case_sha256,
            "conditional_review": (
                None if self.conditional_review is None else self.conditional_review.as_record()
            ),
            "atmospheric_source_case": source,
        }


@dataclass(frozen=True)
class FieldBatchInput:
    """Parsed strict batch input ready for ``export_field_screening_batch``."""

    cases: tuple[FieldBatchInputCase, ...]
    options: FieldBatchExecutionOptions = FieldBatchExecutionOptions()

    def __post_init__(self) -> None:
        if not isinstance(self.cases, tuple):
            raise TypeError("field batch input cases must be a tuple")
        if not self.cases:
            raise ValueError("field batch input requires at least one case")
        labels = [case.label for case in self.cases]
        if len(labels) != len(set(labels)):
            raise ValueError("field batch case labels must be unique")
        if not all(isinstance(case, FieldBatchInputCase) for case in self.cases):
            raise TypeError("field batch cases must be FieldBatchInputCase values")
        if not isinstance(self.options, FieldBatchExecutionOptions):
            raise TypeError("field batch options must be FieldBatchExecutionOptions")

    @property
    def requests(self) -> dict[str, FieldSemiFVRequest]:
        return {case.label: case.request for case in self.cases}

    @property
    def atmospheric_source_schedules(self) -> dict[str, FieldAtmosphericSourceSchedule]:
        return {
            case.label: case.atmospheric_source_schedule
            for case in self.cases
            if case.atmospheric_source_schedule is not None
        }

    @property
    def conditional_review_authorizations(self) -> dict[str, FieldConditionalReviewAuthorization]:
        return {
            case.label: case.conditional_review
            for case in self.cases
            if case.conditional_review is not None
        }

    def as_record(self) -> dict[str, object]:
        return {
            "schema": FIELD_BATCH_INPUT_SCHEMA,
            "options": self.options.as_record(),
            "cases": [case.as_record() for case in self.cases],
        }


def _parse_options(value: object) -> FieldBatchExecutionOptions:
    if value is None:
        return FieldBatchExecutionOptions()
    data = _object(value, "options")
    allowed = {
        "include_refinement", "refinement_factors", "relative_tolerance", "max_cell_steps",
        "allow_conditional_operational_screening",
        "require_refinement_for_operational_screening",
        "require_in_plane_sensor_for_operational_screening",
        "require_resolved_uncertainty_for_operational_screening",
        "max_uncertainty_cases", "table_nodes",
    }
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError("options has unknown keys: " + ", ".join(unknown))
    factors_value = data.get("refinement_factors", [1, 2])
    if not isinstance(factors_value, list) or not factors_value:
        raise ValueError("options.refinement_factors must be a non-empty array")
    factors = tuple(_integer(item, "options.refinement_factors", minimum=1) for item in factors_value)
    return FieldBatchExecutionOptions(
        include_refinement=_boolean(data.get("include_refinement", True), "options.include_refinement"),
        refinement_factors=factors,
        relative_tolerance=_number(data.get("relative_tolerance", 0.05), "options.relative_tolerance"),
        max_cell_steps=_integer(data.get("max_cell_steps", 20_000_000), "options.max_cell_steps", minimum=1),
        allow_conditional_operational_screening=_boolean(
            data.get("allow_conditional_operational_screening", False),
            "options.allow_conditional_operational_screening",
        ),
        require_refinement_for_operational_screening=_boolean(
            data.get("require_refinement_for_operational_screening", True),
            "options.require_refinement_for_operational_screening",
        ),
        require_in_plane_sensor_for_operational_screening=_boolean(
            data.get("require_in_plane_sensor_for_operational_screening", True),
            "options.require_in_plane_sensor_for_operational_screening",
        ),
        require_resolved_uncertainty_for_operational_screening=_boolean(
            data.get("require_resolved_uncertainty_for_operational_screening", True),
            "options.require_resolved_uncertainty_for_operational_screening",
        ),
        max_uncertainty_cases=_integer(
            data.get("max_uncertainty_cases", 64), "options.max_uncertainty_cases", minimum=1,
        ),
        table_nodes=_integer(data.get("table_nodes", 161), "options.table_nodes", minimum=2),
    )


def field_batch_input_from_mapping(
    value: object, *, base_directory: str | Path | None = None,
) -> FieldBatchInput:
    """Parse a strict batch mapping and load every referenced artifact."""
    data = _object(value, "field batch input")
    _keys(data, "field batch input", required={"schema", "cases"}, optional={"options"})
    if data["schema"] != FIELD_BATCH_INPUT_SCHEMA:
        raise ValueError(f"schema must be {FIELD_BATCH_INPUT_SCHEMA!r}")
    if base_directory is None:
        raise ValueError("field batch case paths require a batch file path")
    root = Path(base_directory).resolve()
    cases_value = data["cases"]
    if not isinstance(cases_value, list) or not cases_value:
        raise ValueError("field batch input cases must be a non-empty array")
    parsed: list[FieldBatchInputCase] = []
    labels: set[str] = set()
    for index, raw_case in enumerate(cases_value):
        item = _object(raw_case, f"cases[{index}]")
        _keys(
            item, f"cases[{index}]", required={"label", "case_path"},
            optional={"atmospheric_source_case_path"},
        )
        label = _string(item["label"], f"cases[{index}].label")
        if _SAFE_CASE_LABEL.fullmatch(label) is None:
            raise ValueError(
                f"cases[{index}].label must match the safe report-name pattern"
            )
        if label in labels:
            raise ValueError(f"duplicate field batch case label: {label}")
        labels.add(label)
        case_path = Path(_string(item["case_path"], f"cases[{index}].case_path"))
        if not case_path.is_absolute():
            case_path = root / case_path
        case_path = case_path.resolve()
        parsed_case: FieldScreeningCase = read_field_screening_case_json(case_path)
        if parsed_case.imported_history is not None or parsed_case.phase_routing_transport is not None:
            raise ValueError(
                f"cases[{index}] uses historian or phase-routing input; use its dedicated "
                "operational envelope instead of mixing it into field-batch"
            )
        source_path = None
        source_sha = None
        schedule = None
        if "atmospheric_source_case_path" in item and item["atmospheric_source_case_path"] is not None:
            source_path = Path(_string(
                item["atmospheric_source_case_path"],
                f"cases[{index}].atmospheric_source_case_path",
            ))
            if not source_path.is_absolute():
                source_path = root / source_path
            source_path = source_path.resolve()
            schedule = read_field_source_schedule_json(source_path)
            source_sha = hashlib.sha256(source_path.read_bytes()).hexdigest()
        parsed.append(FieldBatchInputCase(
            label=label,
            case_path=case_path,
            case_sha256=hashlib.sha256(case_path.read_bytes()).hexdigest(),
            request=parsed_case.request,
            conditional_review=parsed_case.conditional_review,
            source_case_path=source_path,
            source_case_sha256=source_sha,
            atmospheric_source_schedule=schedule,
        ))
    return FieldBatchInput(tuple(parsed), _parse_options(data.get("options")))


def read_field_batch_input_json(path: str | Path) -> FieldBatchInput:
    """Read a UTF-8 strict batch input and resolve paths relative to it."""
    source = Path(path)
    try:
        value = strict_json_loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid field batch input JSON: {error}") from error
    return field_batch_input_from_mapping(value, base_directory=source.parent)


__all__ = [
    "FIELD_BATCH_INPUT_SCHEMA", "FieldBatchExecutionOptions", "FieldBatchInputCase",
    "FieldBatchInput", "field_batch_input_from_mapping", "read_field_batch_input_json",
]
