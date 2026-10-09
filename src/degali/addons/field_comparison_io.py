"""Strict input contract for imported fixed-sensor model comparisons.

This module is intentionally separate from the field release-case parser.
It accepts already-generated prediction CSV files, but records enough of the
comparison operator that a result cannot silently become a DEGALI-versus-SLABx
or steady-versus-transient ranking when those runs used different conditions.
"""

from __future__ import annotations

import json
import hashlib
import math
import os
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping

from .field_comparison import (
    FieldComparisonBasis,
    FieldComparisonEvidence,
    FieldModelSensorSet,
    field_model_sensor_set_from_csv,
)
from .field_json import strict_json_loads


FIELD_MODEL_COMPARISON_INPUT_SCHEMA = "degali.field-model-comparison-input.v1"


@dataclass(frozen=True)
class FieldModelComparisonCase:
    """Two declared prediction sets and the numerical comparison tolerances."""

    left: FieldModelSensorSet
    right: FieldModelSensorSet
    threshold_mole_fraction: float = 0.04
    position_tolerance_m: float = 1.0e-9
    comparison_evidence: FieldComparisonEvidence | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.left, FieldModelSensorSet) or not isinstance(
            self.right, FieldModelSensorSet
        ):
            raise TypeError("comparison case models must be FieldModelSensorSet values")
        for name, value in {
            "threshold_mole_fraction": self.threshold_mole_fraction,
            "position_tolerance_m": self.position_tolerance_m,
        }.items():
            if isinstance(value, bool):
                raise ValueError(f"{name} must be finite and non-negative")
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
            if name == "threshold_mole_fraction" and float(value) > 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")
        if self.comparison_evidence is not None and not isinstance(
            self.comparison_evidence, FieldComparisonEvidence
        ):
            raise TypeError("comparison_evidence must be a FieldComparisonEvidence or None")


def _object(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _keys(
    value: Mapping[str, Any], name: str, *, required: set[str], optional: set[str] | None = None,
) -> None:
    optional = set() if optional is None else optional
    missing = sorted(required - set(value))
    unknown = sorted(set(value) - required - optional)
    if missing or unknown:
        parts = []
        if missing:
            parts.append("missing=" + ", ".join(missing))
        if unknown:
            parts.append("unknown=" + ", ".join(unknown))
        raise ValueError(f"{name} keys are invalid: " + "; ".join(parts))


def _string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _number(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _optional_string(value: Mapping[str, Any], key: str, name: str) -> str | None:
    if key not in value or value[key] is None:
        return None
    return _string(value[key], f"{name}.{key}")


def _basis(value: object, name: str) -> FieldComparisonBasis:
    data = _object(value, name)
    _keys(
        data, name,
        required={
            "source_boundary_id", "weather_id", "sensor_geometry_id", "temporal_operator_id",
        },
        optional={"averaging_time_s", "obstacle_representation_id"},
    )
    return FieldComparisonBasis(
        source_boundary_id=_string(data["source_boundary_id"], f"{name}.source_boundary_id"),
        weather_id=_string(data["weather_id"], f"{name}.weather_id"),
        sensor_geometry_id=_string(data["sensor_geometry_id"], f"{name}.sensor_geometry_id"),
        temporal_operator_id=_string(data["temporal_operator_id"], f"{name}.temporal_operator_id"),
        averaging_time_s=(
            None if "averaging_time_s" not in data or data["averaging_time_s"] is None
            else _number(data["averaging_time_s"], f"{name}.averaging_time_s")
        ),
        obstacle_representation_id=(
            None if "obstacle_representation_id" not in data
            or data["obstacle_representation_id"] is None
            else _string(
                data["obstacle_representation_id"],
                f"{name}.obstacle_representation_id",
            )
        ),
    )


def _model_set(
    value: object, name: str, *, base_directory: Path,
) -> FieldModelSensorSet:
    data = _object(value, name)
    _keys(
        data, name,
        required={
            "model_id", "temporal_mode", "csv_path", "prediction_column",
            "concentration_unit", "comparison_basis",
        },
        optional={
            "runtime_s", "sensor_id_column", "x_column", "y_column", "z_column",
            "averaging_time_column", "csv_sha256", "csv_row_count",
        },
    )
    temporal_mode = _string(data["temporal_mode"], f"{name}.temporal_mode")
    if temporal_mode not in {"steady", "transient"}:
        raise ValueError(f"{name}.temporal_mode must be steady or transient")
    concentration_unit = _string(data["concentration_unit"], f"{name}.concentration_unit")
    if concentration_unit not in {"mole_fraction", "volume_percent"}:
        raise ValueError(
            f"{name}.concentration_unit must be mole_fraction or volume_percent"
        )
    csv_path = Path(_string(data["csv_path"], f"{name}.csv_path"))
    if not csv_path.is_absolute():
        csv_path = base_directory / csv_path
    runtime_s = (
        None if "runtime_s" not in data or data["runtime_s"] is None
        else _number(data["runtime_s"], f"{name}.runtime_s")
    )
    expected_sha256 = _optional_string(data, "csv_sha256", name)
    expected_row_count = (
        None if "csv_row_count" not in data
        else _positive_integer(data["csv_row_count"], f"{name}.csv_row_count")
    )
    if (expected_sha256 is None) != (expected_row_count is None):
        raise ValueError(
            f"{name}.csv_sha256 and {name}.csv_row_count must be supplied together"
        )
    model = field_model_sensor_set_from_csv(
        csv_path,
        model_id=_string(data["model_id"], f"{name}.model_id"),
        temporal_mode=temporal_mode,
        basis=_basis(data["comparison_basis"], f"{name}.comparison_basis"),
        prediction_column=_string(data["prediction_column"], f"{name}.prediction_column"),
        concentration_unit=concentration_unit,
        sensor_id_column=_optional_string(data, "sensor_id_column", name) or "sensor",
        x_column=_optional_string(data, "x_column", name) or "x_downwind_m",
        y_column=_optional_string(data, "y_column", name) or "y_crosswind_m",
        z_column=_optional_string(data, "z_column", name) or "height_m",
        averaging_time_column=_optional_string(data, "averaging_time_column", name)
        if "averaging_time_column" in data else "averaging_time_s",
        runtime_s=runtime_s,
    )
    if expected_sha256 is not None:
        normalized_sha256 = expected_sha256.lower()
        if (
            expected_sha256 != normalized_sha256
            or len(normalized_sha256) != 64
            or any(character not in "0123456789abcdef" for character in normalized_sha256)
        ):
            raise ValueError(f"{name}.csv_sha256 must be a lowercase SHA-256 digest")
        if model.csv_provenance is None:  # defensive boundary if the importer changes
            raise ValueError(f"{name} CSV provenance is missing")
        if model.csv_provenance.sha256 != normalized_sha256:
            raise ValueError(f"{name} CSV SHA-256 does not match the supplied artifact")
        if model.csv_provenance.row_count != expected_row_count:
            raise ValueError(f"{name} CSV row count does not match the supplied artifact")
        model = replace(
            model,
            csv_provenance=replace(model.csv_provenance, integrity_pinned=True),
        )
    return model


def _comparison_evidence(
    value: object,
    *,
    base_directory: Path,
) -> FieldComparisonEvidence:
    data = _object(value, "comparison_evidence")
    _keys(
        data, "comparison_evidence",
        required={"manifest_path", "sha256", "qualification"},
    )
    manifest = Path(_string(data["manifest_path"], "comparison_evidence.manifest_path"))
    if not manifest.is_absolute():
        manifest = base_directory / manifest
    if not manifest.is_file():
        raise FileNotFoundError(manifest)
    expected = _string(data["sha256"], "comparison_evidence.sha256").lower()
    observed = hashlib.sha256(manifest.read_bytes()).hexdigest()
    if expected != observed:
        raise ValueError("comparison_evidence manifest SHA-256 does not match the supplied file")
    return FieldComparisonEvidence(
        path=str(manifest.resolve()), sha256=observed,
        qualification=_string(data["qualification"], "comparison_evidence.qualification"),
    )


def field_model_comparison_case_from_mapping(
    value: object, *, base_directory: str | Path | None = None,
) -> FieldModelComparisonCase:
    """Parse a fully declared comparison case without filename-based inference."""
    data = _object(value, "field model comparison case")
    _keys(
        data, "field model comparison case",
        required={"schema", "left", "right"},
        optional={"threshold_mole_fraction", "position_tolerance_m", "comparison_evidence"},
    )
    if data["schema"] != FIELD_MODEL_COMPARISON_INPUT_SCHEMA:
        raise ValueError(
            "schema must be " + repr(FIELD_MODEL_COMPARISON_INPUT_SCHEMA)
        )
    if base_directory is None:
        raise ValueError("field model comparison CSV paths require a case file path")
    directory = Path(base_directory).resolve()
    return FieldModelComparisonCase(
        left=_model_set(data["left"], "left", base_directory=directory),
        right=_model_set(data["right"], "right", base_directory=directory),
        threshold_mole_fraction=_number(
            data.get("threshold_mole_fraction", 0.04), "threshold_mole_fraction",
        ),
        position_tolerance_m=_number(
            data.get("position_tolerance_m", 1.0e-9), "position_tolerance_m",
        ),
        comparison_evidence=(
            None if "comparison_evidence" not in data
            else _comparison_evidence(data["comparison_evidence"], base_directory=directory)
        ),
    )


def read_field_model_comparison_case_json(path: str | Path) -> FieldModelComparisonCase:
    """Read a UTF-8 comparison case and resolve its CSVs relative to it."""
    source = Path(path)
    try:
        value = strict_json_loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid field model comparison JSON: {error}") from error
    return field_model_comparison_case_from_mapping(value, base_directory=source.parent)


def _comparison_basis_record(basis: FieldComparisonBasis) -> dict[str, object]:
    return {
        "source_boundary_id": basis.source_boundary_id,
        "weather_id": basis.weather_id,
        "sensor_geometry_id": basis.sensor_geometry_id,
        "temporal_operator_id": basis.temporal_operator_id,
        "averaging_time_s": basis.averaging_time_s,
        "obstacle_representation_id": basis.obstacle_representation_id,
    }


def _comparison_model_record(
    model: FieldModelSensorSet,
    *,
    name: str,
    base_directory: Path,
) -> dict[str, object]:
    provenance = model.csv_provenance
    if provenance is None:
        raise ValueError(f"{name} model has no CSV provenance to export")
    csv_path = Path(provenance.path).resolve()
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)
    observed_digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if observed_digest != provenance.sha256:
        raise ValueError(f"{name} CSV changed since its provenance was recorded")
    verified = field_model_sensor_set_from_csv(
        csv_path,
        model_id=model.model_id,
        temporal_mode=model.temporal_mode,
        basis=model.basis,
        prediction_column=provenance.prediction_column,
        concentration_unit=provenance.concentration_unit,
        sensor_id_column=provenance.sensor_id_column,
        x_column=provenance.x_column,
        y_column=provenance.y_column,
        z_column=provenance.z_column,
        averaging_time_column=provenance.averaging_time_column,
        runtime_s=model.runtime_s,
    )
    if verified.predictions != model.predictions:
        raise ValueError(f"{name} CSV predictions differ from the typed model set")
    relative_path = Path(os.path.relpath(csv_path, base_directory)).as_posix()
    return {
        "model_id": model.model_id,
        "temporal_mode": model.temporal_mode,
        "runtime_s": model.runtime_s,
        "csv_path": relative_path,
        "prediction_column": provenance.prediction_column,
        "concentration_unit": provenance.concentration_unit,
        "sensor_id_column": provenance.sensor_id_column,
        "x_column": provenance.x_column,
        "y_column": provenance.y_column,
        "z_column": provenance.z_column,
        "averaging_time_column": provenance.averaging_time_column,
        "csv_sha256": observed_digest,
        "csv_row_count": provenance.row_count,
        "comparison_basis": _comparison_basis_record(model.basis),
    }


def field_model_comparison_case_record(
    case: FieldModelComparisonCase,
    *,
    base_directory: str | Path,
) -> dict[str, object]:
    """Serialize a typed comparison case with fresh, pinned CSV fingerprints."""
    if not isinstance(case, FieldModelComparisonCase):
        raise TypeError("case must be a FieldModelComparisonCase value")
    directory = Path(base_directory).resolve()
    evidence = None
    if case.comparison_evidence is not None:
        evidence_path = Path(case.comparison_evidence.path).resolve()
        if not evidence_path.is_file():
            raise FileNotFoundError(evidence_path)
        observed_digest = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
        if observed_digest != case.comparison_evidence.sha256:
            raise ValueError("comparison evidence changed since its provenance was recorded")
        evidence = {
            "manifest_path": Path(os.path.relpath(evidence_path, directory)).as_posix(),
            "sha256": observed_digest,
            "qualification": case.comparison_evidence.qualification,
        }
    record: dict[str, object] = {
        "schema": FIELD_MODEL_COMPARISON_INPUT_SCHEMA,
        "threshold_mole_fraction": case.threshold_mole_fraction,
        "position_tolerance_m": case.position_tolerance_m,
        "left": _comparison_model_record(case.left, name="left", base_directory=directory),
        "right": _comparison_model_record(case.right, name="right", base_directory=directory),
    }
    if evidence is not None:
        record["comparison_evidence"] = evidence
    return record


def write_field_model_comparison_case_json(
    case: FieldModelComparisonCase,
    path: str | Path,
) -> FieldModelComparisonCase:
    """Write and re-read a strict comparison case with pinned CSV artefacts."""
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(
            f"refusing to overwrite existing field comparison case: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    record = field_model_comparison_case_record(
        case, base_directory=destination.parent,
    )
    payload = json.dumps(
        record, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False,
    ) + "\n"
    try:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing field comparison case: {destination}"
        ) from error
    return read_field_model_comparison_case_json(destination)


__all__ = [
    "FIELD_MODEL_COMPARISON_INPUT_SCHEMA", "FieldModelComparisonCase",
    "field_model_comparison_case_from_mapping", "read_field_model_comparison_case_json",
    "field_model_comparison_case_record", "write_field_model_comparison_case_json",
]
