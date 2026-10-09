"""Matched observed-sensor validation for the field model.

This module is deliberately narrower than the model-form comparison helper.
It imports one already-selected, common-clock observation per sensor and
scores a supplied model set against those observations.  It does not infer a
source, reconstruct a time history, interpolate a hazard distance, or promote
an obstacle result to design basis.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import math
from collections import Counter
from pathlib import Path
from typing import Literal

from .field_comparison import (
    FIELD_MODEL_COMPARISON_GATE_CODES,
    FieldComparisonBasis,
    FieldModelComparison,
    FieldModelSensorSet,
    FieldSensorCsvProvenance,
    FieldSensorPrediction,
    compare_field_model_sensor_sets,
)
from .field_contracts import FieldValidationEvidence


FIELD_VALIDATION_SCORE_SCHEMA = "degali.field-validation-score.v1"
ObservationKind = Literal["exact", "lower_bound"]
FIELD_VALIDATION_GATE_CODES = FIELD_MODEL_COMPARISON_GATE_CODES | frozenset({
    "comparison_blocked",
    "comparison_conditional",
    "obstacle_validation_unsupported",
    "lower_bound_present",
    "lower_bound_violated",
    "lower_bound_unresolved",
    "validation_qualified",
})


def _validate_validation_gate_codes(values: tuple[str, ...]) -> None:
    if not isinstance(values, tuple):
        raise TypeError("validation score gate_codes must be a tuple")
    if any(
        not isinstance(value, str)
        or not value.strip()
        or value not in FIELD_VALIDATION_GATE_CODES
        for value in values
    ):
        raise ValueError("validation score gate_codes contain an unsupported code")
    if len(set(values)) != len(values):
        raise ValueError("validation score gate_codes must be unique")


@dataclass(frozen=True)
class FieldValidationObservation:
    """One measured concentration summary at one fixed receptor."""

    sensor_id: str
    position_m: tuple[float, float, float]
    mole_fraction: float
    observation_time_s: float
    observation_kind: ObservationKind = "exact"

    def __post_init__(self) -> None:
        if not isinstance(self.sensor_id, str) or not self.sensor_id.strip():
            raise ValueError("validation sensor_id must be a non-empty string")
        if not isinstance(self.position_m, tuple):
            raise TypeError("validation sensor position must be a tuple")
        if any(isinstance(value, bool) for value in self.position_m):
            raise TypeError("validation sensor coordinates must be numeric, not boolean")
        if len(self.position_m) != 3 or not all(
            math.isfinite(float(value)) for value in self.position_m
        ):
            raise ValueError("validation sensor position must contain three finite values")
        if isinstance(self.mole_fraction, bool):
            raise TypeError("validation mole_fraction must be numeric, not boolean")
        if not math.isfinite(float(self.mole_fraction)) or not 0.0 <= self.mole_fraction <= 1.0:
            raise ValueError("validation mole_fraction must lie in [0, 1]")
        if isinstance(self.observation_time_s, bool):
            raise TypeError("validation observation_time_s must be numeric, not boolean")
        if not math.isfinite(float(self.observation_time_s)):
            raise ValueError("validation observation_time_s must be finite")
        if not isinstance(self.observation_kind, str) or self.observation_kind not in {
            "exact", "lower_bound"
        }:
            raise ValueError(
                "validation observation_kind must be exact or lower_bound"
            )

    @property
    def is_lower_bound(self) -> bool:
        """Whether ``mole_fraction`` is a one-sided lower-bound observation."""
        return self.observation_kind == "lower_bound"


@dataclass(frozen=True)
class FieldValidationDataset:
    """A fingerprinted, common-clock fixed-sensor observation set."""

    evidence: FieldValidationEvidence
    basis: FieldComparisonBasis
    temporal_mode: Literal["steady", "transient"]
    observations: tuple[FieldValidationObservation, ...]
    csv_provenance: FieldSensorCsvProvenance

    def __post_init__(self) -> None:
        if not isinstance(self.evidence, FieldValidationEvidence):
            raise TypeError("validation evidence must be a FieldValidationEvidence")
        if not isinstance(self.basis, FieldComparisonBasis):
            raise TypeError("validation basis must be a FieldComparisonBasis")
        if not isinstance(self.temporal_mode, str) or self.temporal_mode not in {
            "steady", "transient"
        }:
            raise ValueError("validation temporal_mode must be steady or transient")
        if not isinstance(self.observations, tuple):
            raise TypeError("validation observations must be a tuple")
        if not self.observations:
            raise ValueError("validation dataset must contain at least one observation")
        if any(not isinstance(item, FieldValidationObservation) for item in self.observations):
            raise TypeError("validation observations must be FieldValidationObservation values")
        if len({item.sensor_id for item in self.observations}) != len(self.observations):
            raise ValueError(
                "validation dataset requires exactly one summary observation per sensor"
            )
        if not isinstance(self.csv_provenance, FieldSensorCsvProvenance):
            raise TypeError("validation csv_provenance must be a FieldSensorCsvProvenance")
        if self.csv_provenance.row_count != self.evidence.row_count:
            raise ValueError("validation evidence row_count must match the CSV row count")
        if len(self.observations) != self.evidence.row_count:
            raise ValueError(
                "validation observation count must match the evidence row_count"
            )
        if self.csv_provenance.sha256 != self.evidence.sha256:
            raise ValueError("validation evidence sha256 must match the CSV fingerprint")
        if Path(self.csv_provenance.path).resolve() != Path(self.evidence.path).resolve():
            raise ValueError(
                "validation evidence path must match the CSV provenance path"
            )
        # The evidence record is the provenance boundary for the observed
        # dataset, while ``basis`` is the operator actually used for scoring.
        # Allowing these identifiers to disagree would let a correctly
        # fingerprinted CSV be scored as if it belonged to a different source,
        # weather record, receptor layout, or temporal operator.  Obstacle IDs
        # are intentionally not compared here: the physical geometry ID and a
        # model's mask/wake representation ID are different namespaces and are
        # checked separately by ``score_field_model_against_validation``.
        basis_evidence_pairs = {
            "source boundary": (
                self.basis.source_boundary_id,
                self.evidence.source_boundary_id,
            ),
            "weather": (self.basis.weather_id, self.evidence.weather_id),
            "sensor geometry": (
                self.basis.sensor_geometry_id,
                self.evidence.receptor_geometry_id,
            ),
            "temporal operator": (
                self.basis.temporal_operator_id,
                self.evidence.temporal_operator_id,
            ),
        }
        mismatches = tuple(
            f"validation basis {name} {basis!r} does not match evidence {evidence!r}"
            for name, (basis, evidence) in basis_evidence_pairs.items()
            if basis != evidence
        )
        if mismatches:
            raise ValueError("; ".join(mismatches))


@dataclass(frozen=True)
class FieldValidationScore:
    """A matched observed-versus-model score with an explicit disposition."""

    dataset: FieldValidationDataset
    model: FieldModelSensorSet
    comparison: FieldModelComparison
    status: Literal["qualified", "conditional", "withheld"]
    reasons: tuple[str, ...]
    mean_absolute_error_mole_fraction: float | None
    root_mean_square_error_mole_fraction: float | None
    mean_bias_model_minus_observed: float | None
    maximum_absolute_error_mole_fraction: float | None
    lower_bound_observation_count: int = 0
    lower_bound_satisfied_count: int | None = None
    lower_bound_satisfaction_fraction: float | None = None
    mean_lower_bound_deficit_mole_fraction: float | None = None
    maximum_lower_bound_deficit_mole_fraction: float | None = None
    lower_bound_constraint_status: Literal[
        "not_applicable", "satisfied", "violated", "withheld"
    ] = "not_applicable"
    gate_codes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.dataset, FieldValidationDataset):
            raise TypeError("validation score dataset must be a FieldValidationDataset")
        if not isinstance(self.model, FieldModelSensorSet):
            raise TypeError("validation score model must be a FieldModelSensorSet")
        if not isinstance(self.comparison, FieldModelComparison):
            raise TypeError("validation score comparison must be a FieldModelComparison")
        if self.status not in {"qualified", "conditional", "withheld"}:
            raise ValueError("validation score status is unsupported")
        if not isinstance(self.reasons, tuple):
            raise TypeError("validation score reasons must be a tuple")
        if any(not isinstance(value, str) or not value.strip() for value in self.reasons):
            raise ValueError("validation score reasons must be non-empty strings")
        _validate_validation_gate_codes(self.gate_codes)
        if self.status == "conditional" and not self.reasons:
            raise ValueError("conditional validation score requires a reason")
        if self.status == "withheld" and not self.reasons:
            raise ValueError("withheld validation score requires a reason")
        if self.model.model_id != self.comparison.left.model_id:
            raise ValueError("validation score model must match the comparison left model")
        if self.model != self.comparison.left:
            raise ValueError(
                "validation score model must match the complete comparison left model"
            )
        expected_observed = FieldModelSensorSet(
            model_id=f"observed:{self.dataset.evidence.dataset_id}",
            temporal_mode=self.dataset.temporal_mode,
            basis=self.dataset.basis,
            predictions=tuple(
                FieldSensorPrediction(
                    item.sensor_id, item.position_m, item.mole_fraction,
                )
                for item in self.dataset.observations
            ),
            csv_provenance=self.dataset.csv_provenance,
        )
        if self.comparison.right != expected_observed:
            raise ValueError(
                "validation score comparison right model must match the observed dataset"
            )
        if self.comparison.applicability.status == "blocked" and self.comparison.rows:
            raise ValueError("blocked validation comparison must not contain matched rows")
        metric_values = {
            "mean_absolute_error_mole_fraction": self.mean_absolute_error_mole_fraction,
            "root_mean_square_error_mole_fraction": self.root_mean_square_error_mole_fraction,
            "mean_bias_model_minus_observed": self.mean_bias_model_minus_observed,
            "maximum_absolute_error_mole_fraction": self.maximum_absolute_error_mole_fraction,
        }
        if (
            isinstance(self.lower_bound_observation_count, bool)
            or not isinstance(self.lower_bound_observation_count, int)
            or self.lower_bound_observation_count < 0
        ):
            raise ValueError("lower_bound_observation_count must be a non-negative integer")
        if self.lower_bound_satisfied_count is not None and (
            isinstance(self.lower_bound_satisfied_count, bool)
            or not isinstance(self.lower_bound_satisfied_count, int)
            or self.lower_bound_satisfied_count < 0
            or self.lower_bound_satisfied_count > self.lower_bound_observation_count
        ):
            raise ValueError(
                "lower_bound_satisfied_count must lie between zero and the lower-bound count"
            )
        if self.lower_bound_constraint_status not in {
            "not_applicable", "satisfied", "violated", "withheld",
        }:
            raise ValueError("lower_bound_constraint_status is unsupported")
        lower_bound_metrics = {
            "lower_bound_satisfaction_fraction": self.lower_bound_satisfaction_fraction,
            "mean_lower_bound_deficit_mole_fraction": self.mean_lower_bound_deficit_mole_fraction,
            "maximum_lower_bound_deficit_mole_fraction": self.maximum_lower_bound_deficit_mole_fraction,
        }
        if any(isinstance(value, bool) for value in lower_bound_metrics.values() if value is not None):
            raise TypeError("lower-bound metrics must be numeric, not boolean")
        if any(
            value is not None and not math.isfinite(float(value))
            for value in lower_bound_metrics.values()
        ):
            raise ValueError("lower-bound metrics must be finite when supplied")
        for name, value in lower_bound_metrics.items():
            if name != "lower_bound_satisfaction_fraction" and value is not None and float(value) < 0.0:
                raise ValueError(f"{name} must be non-negative")
        if self.lower_bound_satisfaction_fraction is not None and not 0.0 <= float(
            self.lower_bound_satisfaction_fraction
        ) <= 1.0:
            raise ValueError("lower_bound_satisfaction_fraction must lie in [0, 1]")
        if self.lower_bound_observation_count == 0:
            if self.lower_bound_constraint_status != "not_applicable":
                raise ValueError("lower-bound status must be not_applicable when no lower bounds exist")
            if self.lower_bound_satisfied_count is not None or any(
                value is not None for value in lower_bound_metrics.values()
            ):
                raise ValueError("lower-bound metrics must be None when no lower bounds exist")
        elif self.lower_bound_satisfied_count is None:
            if self.lower_bound_constraint_status != "withheld":
                raise ValueError("unscored lower bounds require a withheld status")
        expected_lower_bound_count = sum(
            item.is_lower_bound for item in self.dataset.observations
        )
        if self.lower_bound_observation_count != expected_lower_bound_count:
            raise ValueError("lower_bound_observation_count does not match the dataset")
        if any(isinstance(value, bool) for value in metric_values.values() if value is not None):
            raise TypeError("validation score metrics must be numeric, not boolean")
        if any(
            value is not None and not math.isfinite(float(value))
            for value in metric_values.values()
        ):
            raise ValueError("validation score metrics must be finite when supplied")
        for name in (
            "mean_absolute_error_mole_fraction",
            "root_mean_square_error_mole_fraction",
            "maximum_absolute_error_mole_fraction",
        ):
            value = metric_values[name]
            if value is not None and float(value) < 0.0:
                raise ValueError(f"{name} must be non-negative")
        if self.comparison.rows:
            lower_bound_ids = {
                item.sensor_id for item in self.dataset.observations if item.is_lower_bound
            }
            errors = tuple(
                row.left_mole_fraction - row.right_mole_fraction
                for row in self.comparison.rows
                if row.sensor_id not in lower_bound_ids
            )
            if not errors:
                if any(value is not None for value in metric_values.values()):
                    raise ValueError("symmetric validation metrics must be None when no exact rows exist")
                errors = ()
            if errors:
                expected = {
                    "mean_absolute_error_mole_fraction": math.fsum(abs(value) for value in errors) / len(errors),
                    "root_mean_square_error_mole_fraction": math.sqrt(
                        math.fsum(value * value for value in errors) / len(errors)
                    ),
                    "mean_bias_model_minus_observed": math.fsum(errors) / len(errors),
                    "maximum_absolute_error_mole_fraction": max(abs(value) for value in errors),
                }
                for name, expected_value in expected.items():
                    actual = metric_values[name]
                    if actual is None or not math.isclose(
                        float(actual), expected_value, rel_tol=1.0e-9, abs_tol=1.0e-12,
                    ):
                        raise ValueError(f"{name} does not match the validation comparison rows")
            if lower_bound_ids:
                lower_rows = tuple(
                    row for row in self.comparison.rows if row.sensor_id in lower_bound_ids
                )
                margins = tuple(
                    row.left_mole_fraction - row.right_mole_fraction for row in lower_rows
                )
                if len(lower_rows) != len(lower_bound_ids):
                    raise ValueError("comparison rows do not cover every lower-bound observation")
                expected_satisfied = sum(margin >= -1.0e-12 for margin in margins)
                expected_fraction = expected_satisfied / len(margins)
                expected_deficits = tuple(max(-margin, 0.0) for margin in margins)
                expected_status = (
                    "satisfied" if expected_satisfied == len(margins) else "violated"
                )
                lower_expected = {
                    "lower_bound_satisfied_count": expected_satisfied,
                    "lower_bound_satisfaction_fraction": expected_fraction,
                    "mean_lower_bound_deficit_mole_fraction": math.fsum(expected_deficits) / len(margins),
                    "maximum_lower_bound_deficit_mole_fraction": max(expected_deficits),
                }
                for name, expected_value in lower_expected.items():
                    actual = getattr(self, name)
                    if actual is None or (
                        isinstance(expected_value, float)
                        and not math.isclose(
                            float(actual), expected_value,
                            rel_tol=1.0e-9, abs_tol=1.0e-12,
                        )
                    ) or (
                        not isinstance(expected_value, float) and actual != expected_value
                    ):
                        raise ValueError(f"{name} does not match the lower-bound comparison rows")
                if self.lower_bound_constraint_status != expected_status:
                    raise ValueError("lower_bound_constraint_status does not match the comparison rows")
        elif any(value is not None for value in metric_values.values()):
            raise ValueError("validation score metrics must be None when no matched rows exist")
        if self.status == "qualified" and (
            self.reasons or self.comparison.applicability.status != "accepted"
            or self.lower_bound_observation_count
        ):
            raise ValueError("qualified validation score requires an accepted comparison without reasons")


def _string(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _number(value: str, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def field_validation_dataset_from_csv(
    path: str | Path,
    *,
    evidence: FieldValidationEvidence,
    basis: FieldComparisonBasis,
    temporal_mode: Literal["steady", "transient"],
    concentration_column: str = "observed_mole_fraction",
    concentration_unit: Literal["mole_fraction", "volume_percent"] = "mole_fraction",
    sensor_id_column: str = "sensor",
    x_column: str = "x_downwind_m",
    y_column: str = "y_crosswind_m",
    z_column: str = "height_m",
    time_column: str = "observation_time_s",
    averaging_time_column: str = "averaging_time_s",
    common_clock_column: str = "common_clock_id",
    obstacle_geometry_column: str = "obstacle_geometry_id",
    observation_kind_column: str | None = None,
) -> FieldValidationDataset:
    """Read one observed summary per sensor under an explicit evidence record.

    Every row must repeat the evidence common-clock and obstacle-geometry IDs.
    This prevents a CSV assembled from different clocks or obstacle cases from
    being scored as one validation experiment.  The file digest and row count
    must also match :class:`FieldValidationEvidence` exactly.
    """
    if not isinstance(evidence, FieldValidationEvidence):
        raise TypeError("evidence must be a FieldValidationEvidence")
    if not isinstance(basis, FieldComparisonBasis):
        raise TypeError("basis must be a FieldComparisonBasis")
    if temporal_mode not in {"steady", "transient"}:
        raise ValueError("temporal_mode must be steady or transient")
    if concentration_unit not in {"mole_fraction", "volume_percent"}:
        raise ValueError("concentration_unit must be mole_fraction or volume_percent")
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    # The evidence path is part of the external-artifact provenance boundary.
    # A matching digest alone is not enough here: accepting a copied file under
    # a different path would let a report name an artefact that was never the
    # file actually parsed for the observations.
    evidence_path = Path(evidence.path)
    if evidence_path.resolve() != source.resolve():
        raise ValueError(
            "validation evidence path does not match the observed CSV path: "
            f"{evidence_path.resolve()} != {source.resolve()}"
        )
    column_roles = [
        sensor_id_column, x_column, y_column, z_column, time_column,
        concentration_column, averaging_time_column, common_clock_column,
        obstacle_geometry_column,
    ]
    if observation_kind_column is not None:
        column_roles.append(observation_kind_column)
    if any(
        not isinstance(name, str) or not name.strip()
        or any(character in name for character in "\r\n")
        for name in column_roles
    ):
        raise ValueError("validation CSV column names must be non-empty and single-line")
    role_duplicates = sorted(
        name for name, count in Counter(column_roles).items() if count > 1
    )
    if role_duplicates:
        raise ValueError(
            "validation CSV column roles must be distinct: "
            + ", ".join(role_duplicates)
        )
    required = {
        sensor_id_column, x_column, y_column, z_column, time_column,
        concentration_column, common_clock_column, obstacle_geometry_column,
    }
    if basis.averaging_time_s is not None:
        required.add(averaging_time_column)
    if observation_kind_column is not None:
        required.add(observation_kind_column)
    with source.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        raw_fieldnames = next(reader, None)
        fieldnames = tuple(raw_fieldnames or ())
        if not fieldnames:
            raise ValueError("validation CSV must contain a header row")
        if any(not name.strip() for name in fieldnames):
            raise ValueError("validation CSV header contains an empty column name")
        normalized = tuple(name.strip() for name in fieldnames)
        header_duplicates = sorted(
            name for name, count in Counter(normalized).items() if count > 1
        )
        if header_duplicates:
            raise ValueError(
                "validation CSV header contains duplicate columns: "
                + ", ".join(header_duplicates)
            )
        available = set(fieldnames)
        missing = sorted(required - available)
        if missing:
            raise ValueError(
                "validation CSV is missing required columns: " + ", ".join(missing)
            )
        rows: list[dict[str, str]] = []
        for row_number, values in enumerate(reader, start=2):
            if len(values) != len(fieldnames):
                detail = (
                    "more values than header columns"
                    if len(values) > len(fieldnames)
                    else "fewer values than header columns"
                )
                raise ValueError(
                    f"invalid validation CSV row {row_number}: {detail}; "
                    f"expected {len(fieldnames)} values, got {len(values)}"
                )
            rows.append(dict(zip(fieldnames, values)))
    if not rows:
        raise ValueError("validation CSV contains no data rows")
    if len(rows) != evidence.row_count:
        raise ValueError(
            f"validation CSV row_count {len(rows)} does not match evidence {evidence.row_count}"
        )
    observed_digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if observed_digest != evidence.sha256:
        raise ValueError("validation CSV SHA-256 does not match validation evidence")

    observations: list[FieldValidationObservation] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        try:
            sensor_id = _string(row[sensor_id_column], f"row {row_number} sensor_id")
            if sensor_id in seen:
                raise ValueError("validation CSV must contain one row per sensor")
            seen.add(sensor_id)
            clock_id = _string(row[common_clock_column], f"row {row_number} common_clock_id")
            if clock_id != evidence.common_clock_id:
                raise ValueError(
                    f"common_clock_id {clock_id!r} does not match evidence {evidence.common_clock_id!r}"
                )
            obstacle_id = _string(
                row[obstacle_geometry_column], f"row {row_number} obstacle_geometry_id"
            )
            if obstacle_id != evidence.obstacle_geometry_id:
                raise ValueError(
                    f"obstacle_geometry_id {obstacle_id!r} does not match evidence "
                    f"{evidence.obstacle_geometry_id!r}"
                )
            if basis.averaging_time_s is not None:
                averaging = _number(row[averaging_time_column], f"row {row_number} averaging_time_s")
                if not math.isclose(
                    averaging, basis.averaging_time_s, rel_tol=1.0e-9, abs_tol=1.0e-9,
                ):
                    raise ValueError(
                        f"averaging time {averaging:g} s does not match declared "
                        f"basis {basis.averaging_time_s:g} s"
                    )
            value = _number(row[concentration_column], f"row {row_number} concentration")
            if concentration_unit == "volume_percent":
                if not 0.0 <= value <= 100.0:
                    raise ValueError("volume-percent concentration must lie in [0, 100]")
                value /= 100.0
            elif not 0.0 <= value <= 1.0:
                raise ValueError("mole-fraction concentration must lie in [0, 1]")
            observation_kind: ObservationKind = "exact"
            if observation_kind_column is not None:
                observation_kind = _string(
                    row[observation_kind_column],
                    f"row {row_number} observation_kind",
                )
                if observation_kind not in {"exact", "lower_bound"}:
                    raise ValueError(
                        "observation_kind must be exact or lower_bound"
                    )
            observations.append(FieldValidationObservation(
                sensor_id,
                tuple(_number(row[column], f"row {row_number} {column}")
                      for column in (x_column, y_column, z_column)),
                value,
                _number(row[time_column], f"row {row_number} observation_time_s"),
                observation_kind,
            ))
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"invalid validation CSV row {row_number}: {error}") from error

    provenance = FieldSensorCsvProvenance(
        path=str(source.resolve()), sha256=observed_digest, row_count=len(rows),
        prediction_column=concentration_column, concentration_unit=concentration_unit,
        averaging_time_column=(
            averaging_time_column if averaging_time_column in available else None
        ),
        observation_kind_column=(
            observation_kind_column
            if observation_kind_column is not None and observation_kind_column in available
            else None
        ),
        sensor_id_column=sensor_id_column,
        x_column=x_column,
        y_column=y_column,
        z_column=z_column,
        time_column=time_column,
        common_clock_column=common_clock_column,
        obstacle_geometry_column=obstacle_geometry_column,
    )
    return FieldValidationDataset(
        evidence=evidence, basis=basis, temporal_mode=temporal_mode,
        observations=tuple(observations), csv_provenance=provenance,
    )


def score_field_model_against_validation(
    model: FieldModelSensorSet,
    dataset: FieldValidationDataset,
    *,
    threshold_mole_fraction: float = 0.04,
    position_tolerance_m: float = 1.0e-9,
) -> FieldValidationScore:
    """Score model predictions against matched observations without extrapolation."""
    if not isinstance(model, FieldModelSensorSet):
        raise TypeError("model must be a FieldModelSensorSet")
    if not isinstance(dataset, FieldValidationDataset):
        raise TypeError("dataset must be a FieldValidationDataset")
    if isinstance(threshold_mole_fraction, bool) or not math.isfinite(float(threshold_mole_fraction)) or threshold_mole_fraction < 0.0:
        raise ValueError("threshold_mole_fraction must be finite and non-negative")
    if isinstance(position_tolerance_m, bool) or not math.isfinite(float(position_tolerance_m)) or position_tolerance_m < 0.0:
        raise ValueError("position_tolerance_m must be finite and non-negative")
    observed = FieldModelSensorSet(
        model_id=f"observed:{dataset.evidence.dataset_id}",
        temporal_mode=dataset.temporal_mode,
        basis=dataset.basis,
        predictions=tuple(
            FieldSensorPrediction(item.sensor_id, item.position_m, item.mole_fraction)
            for item in dataset.observations
        ),
        csv_provenance=dataset.csv_provenance,
    )
    comparison = compare_field_model_sensor_sets(
        model, observed,
        threshold_mole_fraction=threshold_mole_fraction,
        position_tolerance_m=position_tolerance_m,
    )
    reasons: list[str] = list(comparison.applicability.reasons)
    obstacle_validation_unsupported = False
    if comparison.applicability.status == "conditional":
        reasons.extend(comparison.applicability.warnings)
    if model.basis.obstacle_representation_id is not None and not dataset.evidence.supports_obstacle_transport:
        reasons.append("validation evidence does not cover obstacle transport")
        obstacle_validation_unsupported = True
    if dataset.evidence.supports_obstacle_transport and model.basis.obstacle_representation_id is None:
        reasons.append("obstacle validation requires a declared model obstacle_representation_id")
        obstacle_validation_unsupported = True
    reasons = list(dict.fromkeys(reasons))
    if reasons and comparison.applicability.status == "blocked":
        status: Literal["qualified", "conditional", "withheld"] = "withheld"
    elif reasons:
        status = "withheld" if any("obstacle" in item for item in reasons) else "conditional"
    else:
        status = "qualified"
    lower_bound_ids = {
        item.sensor_id for item in dataset.observations if item.is_lower_bound
    }
    exact_rows = [row for row in comparison.rows if row.sensor_id not in lower_bound_ids]
    lower_bound_rows = [row for row in comparison.rows if row.sensor_id in lower_bound_ids]
    if exact_rows:
        errors = [
            row.left_mole_fraction - row.right_mole_fraction for row in exact_rows
        ]
        absolute = [abs(value) for value in errors]
        mae = sum(absolute) / len(absolute)
        rmse = math.sqrt(sum(value * value for value in errors) / len(errors))
        bias = sum(errors) / len(errors)
        maximum = max(absolute)
    else:
        mae = rmse = bias = maximum = None
    lower_bound_count = len(lower_bound_ids)
    if lower_bound_count and len(lower_bound_rows) == lower_bound_count:
        margins = [row.left_mole_fraction - row.right_mole_fraction for row in lower_bound_rows]
        satisfied_count = sum(margin >= -1.0e-12 for margin in margins)
        deficits = [max(-margin, 0.0) for margin in margins]
        fraction = satisfied_count / lower_bound_count
        mean_deficit = math.fsum(deficits) / lower_bound_count
        max_deficit = max(deficits)
        bound_status: Literal["satisfied", "violated"] = (
            "satisfied" if satisfied_count == lower_bound_count else "violated"
        )
    elif lower_bound_count:
        satisfied_count = None
        fraction = mean_deficit = max_deficit = None
        bound_status = "withheld"
    else:
        satisfied_count = None
        fraction = mean_deficit = max_deficit = None
        bound_status = "not_applicable"
    if lower_bound_count and bound_status == "violated":
        violated = [
            row.sensor_id for row in lower_bound_rows
            if row.left_mole_fraction - row.right_mole_fraction < -1.0e-12
        ]
        reasons.append(
            "model prediction violates lower-bound observation at sensor(s): "
            + ", ".join(violated)
        )
    elif lower_bound_count and bound_status == "withheld":
        reasons.append(
            "lower-bound observations could not be evaluated because the matched comparison is incomplete"
        )
    elif lower_bound_count:
        reasons.append(
            "dataset includes lower-bound observations; symmetric accuracy metrics exclude censored rows"
        )
    reasons = list(dict.fromkeys(reasons))
    if bound_status == "violated" or bound_status == "withheld":
        status = "withheld"
    elif reasons and comparison.applicability.status == "blocked":
        status = "withheld"
    elif reasons:
        status = "withheld" if any("obstacle" in item for item in reasons) else "conditional"
    else:
        status = "qualified"
    if status == "withheld" and not reasons:
        reasons = ["validation score is withheld until the matched comparison basis is resolved"]
    gate_codes = list(comparison.gate_codes)
    if comparison.applicability.status == "blocked":
        gate_codes.append("comparison_blocked")
    elif comparison.applicability.status == "conditional":
        gate_codes.append("comparison_conditional")
    if obstacle_validation_unsupported:
        gate_codes.append("obstacle_validation_unsupported")
    if lower_bound_count:
        gate_codes.append("lower_bound_present")
        if bound_status == "violated":
            gate_codes.append("lower_bound_violated")
        elif bound_status == "withheld":
            gate_codes.append("lower_bound_unresolved")
    if status == "qualified":
        gate_codes.append("validation_qualified")
    return FieldValidationScore(
        dataset=dataset, model=model, comparison=comparison, status=status,
        reasons=tuple(reasons), mean_absolute_error_mole_fraction=mae,
        root_mean_square_error_mole_fraction=rmse,
        mean_bias_model_minus_observed=bias,
        maximum_absolute_error_mole_fraction=maximum,
        lower_bound_observation_count=lower_bound_count,
        lower_bound_satisfied_count=satisfied_count,
        lower_bound_satisfaction_fraction=fraction,
        mean_lower_bound_deficit_mole_fraction=mean_deficit,
        maximum_lower_bound_deficit_mole_fraction=max_deficit,
        lower_bound_constraint_status=bound_status,
        gate_codes=tuple(dict.fromkeys(gate_codes)),
    )


def field_validation_score_record(score: FieldValidationScore) -> dict[str, object]:
    """Return a JSON-native validation score with its non-extrapolation scope."""
    if not isinstance(score, FieldValidationScore):
        raise TypeError("score must be a FieldValidationScore")
    lower_bound_ids = {
        item.sensor_id for item in score.dataset.observations if item.is_lower_bound
    }
    exact_classification_disagreements = sum(
        row.left_above_threshold != row.right_above_threshold
        for row in score.comparison.rows
        if row.sensor_id not in lower_bound_ids
    )
    return {
        "schema": FIELD_VALIDATION_SCORE_SCHEMA,
        "status": score.status,
        "reasons": list(score.reasons),
        "gate_codes": list(score.gate_codes),
        "comparison_gate_codes": list(score.comparison.gate_codes),
        "dataset_id": score.dataset.evidence.dataset_id,
        "dataset_scope": score.dataset.evidence.scope,
        "csv_provenance": {
            "path": score.dataset.csv_provenance.path,
            "sha256": score.dataset.csv_provenance.sha256,
            "row_count": score.dataset.csv_provenance.row_count,
            "averaging_time_column": score.dataset.csv_provenance.averaging_time_column,
            "observation_kind_column": score.dataset.csv_provenance.observation_kind_column,
        },
        "model_id": score.model.model_id,
        "mean_absolute_error_mole_fraction": score.mean_absolute_error_mole_fraction,
        "root_mean_square_error_mole_fraction": score.root_mean_square_error_mole_fraction,
        "mean_bias_model_minus_observed": score.mean_bias_model_minus_observed,
        "maximum_absolute_error_mole_fraction": score.maximum_absolute_error_mole_fraction,
        "observation_counts": {
            "exact": len(score.dataset.observations) - score.lower_bound_observation_count,
            "lower_bound": score.lower_bound_observation_count,
        },
        "lower_bound_constraint": {
            "status": score.lower_bound_constraint_status,
            "satisfied_count": score.lower_bound_satisfied_count,
            "satisfaction_fraction": score.lower_bound_satisfaction_fraction,
            "mean_deficit_mole_fraction": score.mean_lower_bound_deficit_mole_fraction,
            "maximum_deficit_mole_fraction": score.maximum_lower_bound_deficit_mole_fraction,
        },
        "observation_kind_column": score.dataset.csv_provenance.observation_kind_column,
        "classification_disagreement_count": score.comparison.classification_disagreement_count,
        "exact_classification_disagreement_count": exact_classification_disagreements,
        "qualification": (
            "matched fixed-sensor observations only; no spatial interpolation, "
            "hazard-distance claim, design basis or approval decision"
        ),
    }


__all__ = [
    "FIELD_VALIDATION_SCORE_SCHEMA", "FIELD_VALIDATION_GATE_CODES", "ObservationKind", "FieldValidationObservation",
    "FieldValidationDataset", "FieldValidationScore",
    "field_validation_dataset_from_csv", "score_field_model_against_validation",
    "field_validation_score_record",
]
