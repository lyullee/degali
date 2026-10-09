"""Strict JSON case files for matched field validation scores."""

from __future__ import annotations

import json
import math
import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from .field_comparison import FieldComparisonBasis, FieldModelSensorSet
from .field_comparison_io import (
    _comparison_basis_record,
    _comparison_model_record,
    _model_set,
)
from .field_contracts import FieldValidationEvidence
from .field_evidence_manifest import (
    FieldEvidenceManifest,
    read_field_evidence_manifest_json,
)
from .field_validation import (
    FieldValidationDataset,
    field_validation_dataset_from_csv,
)
from .field_json import strict_json_loads


FIELD_VALIDATION_INPUT_SCHEMA = "degali.field-validation-input.v1"


@dataclass(frozen=True)
class FieldValidationManifestReference:
    """Hash-pinned manifest linked to a strict validation case."""

    path: str
    sha256: str
    manifest: FieldEvidenceManifest

    def __post_init__(self) -> None:
        if not isinstance(self.path, str) or not self.path.strip():
            raise ValueError("validation manifest path must be non-empty")
        if not Path(self.path).is_absolute():
            raise ValueError("validation manifest path must be absolute")
        digest = self.sha256.lower() if isinstance(self.sha256, str) else ""
        if (
            digest != self.sha256
            or len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
        ):
            raise ValueError("validation manifest sha256 must be lowercase SHA-256")
        if not isinstance(self.manifest, FieldEvidenceManifest):
            raise TypeError("validation manifest must be a FieldEvidenceManifest")

    def verify(self) -> None:
        path = Path(self.path).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed != self.sha256:
            raise ValueError(
                "validation manifest changed since its provenance was recorded"
            )
        self.manifest.verify_files()


@dataclass(frozen=True)
class FieldValidationCase:
    """One strict observed-data/model validation request."""

    model: FieldModelSensorSet
    dataset: FieldValidationDataset
    threshold_mole_fraction: float = 0.04
    position_tolerance_m: float = 1.0e-9
    manifest_reference: FieldValidationManifestReference | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.model, FieldModelSensorSet):
            raise TypeError("validation case model must be a FieldModelSensorSet")
        if not isinstance(self.dataset, FieldValidationDataset):
            raise TypeError("validation case dataset must be a FieldValidationDataset")
        if self.manifest_reference is not None:
            if not isinstance(self.manifest_reference, FieldValidationManifestReference):
                raise TypeError(
                    "manifest_reference must be a FieldValidationManifestReference or None"
                )
            self.manifest_reference.verify()
            if self.manifest_reference.manifest.as_validation_evidence() != self.dataset.evidence:
                raise ValueError(
                    "validation manifest evidence does not match validation_evidence"
                )
        if isinstance(self.threshold_mole_fraction, bool):
            raise ValueError("threshold_mole_fraction must be finite and non-negative")
        if not math.isfinite(float(self.threshold_mole_fraction)) or self.threshold_mole_fraction < 0.0:
            raise ValueError("threshold_mole_fraction must be finite and non-negative")
        if isinstance(self.position_tolerance_m, bool):
            raise ValueError("position_tolerance_m must be finite and non-negative")
        if not math.isfinite(float(self.position_tolerance_m)) or self.position_tolerance_m < 0.0:
            raise ValueError("position_tolerance_m must be finite and non-negative")


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


def _integer(value: object, name: str, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer of at least {minimum}")
    return value


def _manifest_reference(
    value: object, *, base_directory: Path,
) -> FieldValidationManifestReference:
    data = _object(value, "validation_manifest")
    _keys(data, "validation_manifest", required={"path", "sha256"})
    manifest_path = Path(_string(data["path"], "validation_manifest.path"))
    if not manifest_path.is_absolute():
        manifest_path = base_directory / manifest_path
    manifest_path = manifest_path.resolve()
    sha256 = _string(data["sha256"], "validation_manifest.sha256").lower()
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    observed = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if observed != sha256:
        raise ValueError("validation_manifest SHA-256 does not match its file")
    manifest = read_field_evidence_manifest_json(manifest_path)
    reference = FieldValidationManifestReference(
        path=str(manifest_path), sha256=observed, manifest=manifest,
    )
    reference.verify()
    return reference


def _evidence(value: object, *, base_directory: Path) -> FieldValidationEvidence:
    data = _object(value, "validation_evidence")
    _keys(
        data,
        "validation_evidence",
        required={
            "dataset_id", "path", "sha256", "row_count", "source_boundary_id",
            "weather_id", "obstacle_geometry_id", "receptor_geometry_id",
            "temporal_operator_id", "common_clock_id", "scope",
        },
    )
    evidence_path = Path(_string(data["path"], "validation_evidence.path"))
    if not evidence_path.is_absolute():
        evidence_path = base_directory / evidence_path
    return FieldValidationEvidence(
        dataset_id=_string(data["dataset_id"], "validation_evidence.dataset_id"),
        path=str(evidence_path.resolve()),
        sha256=_string(data["sha256"], "validation_evidence.sha256").lower(),
        row_count=_integer(data["row_count"], "validation_evidence.row_count", minimum=1),
        source_boundary_id=_string(data["source_boundary_id"], "validation_evidence.source_boundary_id"),
        weather_id=_string(data["weather_id"], "validation_evidence.weather_id"),
        obstacle_geometry_id=_string(data["obstacle_geometry_id"], "validation_evidence.obstacle_geometry_id"),
        receptor_geometry_id=_string(data["receptor_geometry_id"], "validation_evidence.receptor_geometry_id"),
        temporal_operator_id=_string(data["temporal_operator_id"], "validation_evidence.temporal_operator_id"),
        common_clock_id=_string(data["common_clock_id"], "validation_evidence.common_clock_id"),
        scope=_string(data["scope"], "validation_evidence.scope"),
    )


def _dataset(
    value: object, *, base_directory: Path, evidence: FieldValidationEvidence,
) -> FieldValidationDataset:
    data = _object(value, "dataset")
    _keys(
        data,
        "dataset",
        required={"csv_path", "temporal_mode", "comparison_basis"},
        optional={
            "concentration_column", "concentration_unit", "sensor_id_column",
            "x_column", "y_column", "z_column", "time_column",
            "averaging_time_column", "common_clock_column", "obstacle_geometry_column",
            "observation_kind_column",
        },
    )
    csv_path = Path(_string(data["csv_path"], "dataset.csv_path"))
    if not csv_path.is_absolute():
        csv_path = base_directory / csv_path
    temporal_mode = _string(data["temporal_mode"], "dataset.temporal_mode")
    if temporal_mode not in {"steady", "transient"}:
        raise ValueError("dataset.temporal_mode must be steady or transient")
    # Reuse the comparison input's fully checked basis/model vocabulary.
    from .field_comparison_io import _basis

    basis = _basis(data["comparison_basis"], "dataset.comparison_basis")
    concentration_unit = _string(
        data.get("concentration_unit", "mole_fraction"), "dataset.concentration_unit",
    )
    if concentration_unit not in {"mole_fraction", "volume_percent"}:
        raise ValueError("dataset.concentration_unit is unsupported")
    return field_validation_dataset_from_csv(
        csv_path,
        evidence=evidence,
        basis=basis,
        temporal_mode=temporal_mode,
        concentration_column=_string(
            data.get("concentration_column", "observed_mole_fraction"),
            "dataset.concentration_column",
        ),
        concentration_unit=concentration_unit,
        sensor_id_column=_string(data.get("sensor_id_column", "sensor"), "dataset.sensor_id_column"),
        x_column=_string(data.get("x_column", "x_downwind_m"), "dataset.x_column"),
        y_column=_string(data.get("y_column", "y_crosswind_m"), "dataset.y_column"),
        z_column=_string(data.get("z_column", "height_m"), "dataset.z_column"),
        time_column=_string(data.get("time_column", "observation_time_s"), "dataset.time_column"),
        averaging_time_column=_string(
            data.get("averaging_time_column", "averaging_time_s"),
            "dataset.averaging_time_column",
        ),
        common_clock_column=_string(
            data.get("common_clock_column", "common_clock_id"),
            "dataset.common_clock_column",
        ),
        obstacle_geometry_column=_string(
            data.get("obstacle_geometry_column", "obstacle_geometry_id"),
            "dataset.obstacle_geometry_column",
        ),
        observation_kind_column=(
            None if "observation_kind_column" not in data
            else _string(data["observation_kind_column"], "dataset.observation_kind_column")
        ),
    )


def field_validation_case_from_mapping(
    value: object, *, base_directory: str | Path | None = None,
) -> FieldValidationCase:
    """Parse a strict validation case and verify both CSV fingerprints."""
    data = _object(value, "field validation case")
    _keys(
        data,
        "field validation case",
        required={"schema", "model", "dataset", "validation_evidence"},
        optional={
            "threshold_mole_fraction", "position_tolerance_m", "validation_manifest",
        },
    )
    if data["schema"] != FIELD_VALIDATION_INPUT_SCHEMA:
        raise ValueError("schema must be " + repr(FIELD_VALIDATION_INPUT_SCHEMA))
    if base_directory is None:
        raise ValueError("field validation CSV paths require a case file path")
    directory = Path(base_directory).resolve()
    evidence = _evidence(data["validation_evidence"], base_directory=directory)
    manifest_reference = (
        _manifest_reference(data["validation_manifest"], base_directory=directory)
        if "validation_manifest" in data else None
    )
    model = _model_set(data["model"], "model", base_directory=directory)
    dataset = _dataset(data["dataset"], base_directory=directory, evidence=evidence)
    return FieldValidationCase(
        model=model,
        dataset=dataset,
        threshold_mole_fraction=_number(
            data.get("threshold_mole_fraction", 0.04), "threshold_mole_fraction",
        ),
        position_tolerance_m=_number(
            data.get("position_tolerance_m", 1.0e-9), "position_tolerance_m",
        ),
        manifest_reference=manifest_reference,
    )


def read_field_validation_case_json(path: str | Path) -> FieldValidationCase:
    """Read a UTF-8 strict validation case and resolve its CSV paths."""
    source = Path(path)
    try:
        value = strict_json_loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid field validation JSON: {error}") from error
    return field_validation_case_from_mapping(value, base_directory=source.parent)


def _validation_evidence_record(
    evidence: FieldValidationEvidence,
    *,
    base_directory: Path,
) -> dict[str, object]:
    evidence_path = Path(evidence.path).resolve()
    if not evidence_path.is_file():
        raise FileNotFoundError(evidence_path)
    observed_digest = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
    if observed_digest != evidence.sha256:
        raise ValueError("validation evidence CSV changed since its provenance was recorded")
    return {
        "dataset_id": evidence.dataset_id,
        "path": Path(os.path.relpath(evidence_path, base_directory)).as_posix(),
        "sha256": observed_digest,
        "row_count": evidence.row_count,
        "source_boundary_id": evidence.source_boundary_id,
        "weather_id": evidence.weather_id,
        "obstacle_geometry_id": evidence.obstacle_geometry_id,
        "receptor_geometry_id": evidence.receptor_geometry_id,
        "temporal_operator_id": evidence.temporal_operator_id,
        "common_clock_id": evidence.common_clock_id,
        "scope": evidence.scope,
    }


def _validation_dataset_record(
    dataset: FieldValidationDataset,
    *,
    evidence: FieldValidationEvidence,
    base_directory: Path,
) -> dict[str, object]:
    provenance = dataset.csv_provenance
    csv_path = Path(provenance.path).resolve()
    evidence_path = Path(evidence.path).resolve()
    if csv_path != evidence_path:
        raise ValueError("validation dataset CSV and evidence paths do not match")
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)
    observed_digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if observed_digest != provenance.sha256 or observed_digest != evidence.sha256:
        raise ValueError("validation dataset CSV changed since its provenance was recorded")
    averaging_column = provenance.averaging_time_column or "averaging_time_s"
    verified = field_validation_dataset_from_csv(
        csv_path,
        evidence=evidence,
        basis=dataset.basis,
        temporal_mode=dataset.temporal_mode,
        concentration_column=provenance.prediction_column,
        concentration_unit=provenance.concentration_unit,
        sensor_id_column=provenance.sensor_id_column,
        x_column=provenance.x_column,
        y_column=provenance.y_column,
        z_column=provenance.z_column,
        time_column=provenance.time_column,
        averaging_time_column=averaging_column,
        common_clock_column=provenance.common_clock_column,
        obstacle_geometry_column=provenance.obstacle_geometry_column,
        observation_kind_column=provenance.observation_kind_column,
    )
    if verified.observations != dataset.observations:
        raise ValueError("validation CSV observations differ from the typed dataset")
    record: dict[str, object] = {
        "csv_path": Path(os.path.relpath(csv_path, base_directory)).as_posix(),
        "temporal_mode": dataset.temporal_mode,
        "concentration_column": provenance.prediction_column,
        "concentration_unit": provenance.concentration_unit,
        "sensor_id_column": provenance.sensor_id_column,
        "x_column": provenance.x_column,
        "y_column": provenance.y_column,
        "z_column": provenance.z_column,
        "time_column": provenance.time_column,
        "averaging_time_column": averaging_column,
        "common_clock_column": provenance.common_clock_column,
        "obstacle_geometry_column": provenance.obstacle_geometry_column,
        "comparison_basis": _comparison_basis_record(dataset.basis),
    }
    if provenance.observation_kind_column is not None:
        record["observation_kind_column"] = provenance.observation_kind_column
    return record


def field_validation_case_record(
    case: FieldValidationCase,
    *,
    base_directory: str | Path,
) -> dict[str, object]:
    """Serialize a typed validation case after rechecking both CSV artefacts."""
    if not isinstance(case, FieldValidationCase):
        raise TypeError("case must be a FieldValidationCase value")
    directory = Path(base_directory).resolve()
    record: dict[str, object] = {
        "schema": FIELD_VALIDATION_INPUT_SCHEMA,
        "threshold_mole_fraction": case.threshold_mole_fraction,
        "position_tolerance_m": case.position_tolerance_m,
        "model": _comparison_model_record(
            case.model, name="model", base_directory=directory,
        ),
        "dataset": _validation_dataset_record(
            case.dataset,
            evidence=case.dataset.evidence,
            base_directory=directory,
        ),
        "validation_evidence": _validation_evidence_record(
            case.dataset.evidence, base_directory=directory,
        ),
    }
    if case.manifest_reference is not None:
        case.manifest_reference.verify()
        manifest_path = Path(case.manifest_reference.path).resolve()
        record["validation_manifest"] = {
            "path": Path(os.path.relpath(manifest_path, directory)).as_posix(),
            "sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        }
    return record


def write_field_validation_case_json(
    case: FieldValidationCase,
    path: str | Path,
) -> FieldValidationCase:
    """Write and re-read a strict validation case with fresh CSV fingerprints."""
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(
            f"refusing to overwrite existing field validation case: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    record = field_validation_case_record(case, base_directory=destination.parent)
    payload = json.dumps(
        record, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False,
    ) + "\n"
    try:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing field validation case: {destination}"
        ) from error
    return read_field_validation_case_json(destination)


__all__ = [
    "FIELD_VALIDATION_INPUT_SCHEMA", "FieldValidationManifestReference",
    "FieldValidationCase",
    "field_validation_case_from_mapping", "read_field_validation_case_json",
    "field_validation_case_record", "write_field_validation_case_json",
]
