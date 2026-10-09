"""Auditable model-form comparisons at matched field sensors.

This module compares supplied predictions; it does not pretend to execute
SLABx or transform a steady value into a transient measurement.  The explicit
comparison basis makes source, meteorology, geometry and temporal-operator
mismatches machine-visible in every result.
"""

from __future__ import annotations

import csv
import hashlib
import math
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Sequence

import numpy as np

from .field_contracts import FieldApplicability

if TYPE_CHECKING:
    from .field_superposition import FieldSensorSuperposition
    from .field_workflow import FieldSensorDeploymentResult, FieldSemiFVScreeningResult


# Stable machine-readable reasons for comparison applicability and the
# decision-facing impact.  Human-readable reasons remain in
# ``FieldApplicability``/``FieldModelDecisionImpact``; these codes are the
# contract consumed by batch gates and downstream audit tooling.
FIELD_MODEL_COMPARISON_GATE_CODES = frozenset({
    "sensor_set_mismatch",
    "sensor_geometry_mismatch",
    "basis_mismatch",
    "averaging_operator_unverified",
    "execution_blocked",
    "execution_conditional",
    "execution_uncertainty_unresolved",
    "temporal_mode_mismatch",
    "matched_basis",
    "model_selection_withheld",
    "classification_disagreement",
    "classification_aligned",
})


def _validate_comparison_gate_codes(
    values: tuple[str, ...],
    *,
    name: str,
) -> None:
    if not isinstance(values, tuple):
        raise TypeError(f"{name} must be a tuple")
    if any(
        not isinstance(value, str)
        or not value.strip()
        or value not in FIELD_MODEL_COMPARISON_GATE_CODES
        for value in values
    ):
        raise ValueError(f"{name} contains an unsupported code")
    if len(set(values)) != len(values):
        raise ValueError(f"{name} must contain unique codes")


@dataclass(frozen=True)
class FieldComparisonBasis:
    """Identifiers for conditions which must match before a direct score."""

    source_boundary_id: str
    weather_id: str
    sensor_geometry_id: str
    temporal_operator_id: str
    averaging_time_s: float | None = None
    obstacle_representation_id: str | None = None

    def __post_init__(self) -> None:
        for name, value in {
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "sensor_geometry_id": self.sensor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
        }.items():
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
        if isinstance(self.averaging_time_s, bool):
            raise TypeError("averaging_time_s must be numeric, not boolean")
        if self.averaging_time_s is not None and (
            not math.isfinite(float(self.averaging_time_s))
            or self.averaging_time_s <= 0.0
        ):
            raise ValueError("averaging_time_s must be positive and finite when supplied")
        if self.obstacle_representation_id is not None and (
            not isinstance(self.obstacle_representation_id, str)
            or not self.obstacle_representation_id.strip()
        ):
            raise ValueError(
                "obstacle_representation_id must be non-empty when supplied"
            )


@dataclass(frozen=True)
class FieldSensorPrediction:
    """One non-negative mole-fraction prediction at a fixed field sensor."""

    sensor_id: str
    position_m: tuple[float, float, float]
    mole_fraction: float

    def __post_init__(self) -> None:
        if not isinstance(self.sensor_id, str) or not self.sensor_id.strip():
            raise ValueError("sensor_id must be a non-empty string")
        if not isinstance(self.position_m, tuple):
            raise TypeError("sensor position must be a tuple")
        if any(isinstance(value, bool) for value in self.position_m):
            raise TypeError("sensor position coordinates must be numeric, not boolean")
        if len(self.position_m) != 3 or not all(math.isfinite(float(v)) for v in self.position_m):
            raise ValueError("sensor position must contain three finite values")
        if isinstance(self.mole_fraction, bool):
            raise TypeError("mole_fraction must be numeric, not boolean")
        if (
            not math.isfinite(float(self.mole_fraction))
            or not 0.0 <= self.mole_fraction <= 1.0
        ):
            raise ValueError("mole_fraction must be finite and lie in [0, 1]")


@dataclass(frozen=True)
class FieldSensorCsvProvenance:
    """Fingerprint of a sensor-prediction CSV read into a model set.

    The fingerprint describes an input artefact, not a claim that the external
    model or its physical source mapping has been independently validated.
    """

    path: str
    sha256: str
    row_count: int
    prediction_column: str
    concentration_unit: Literal["mole_fraction", "volume_percent"]
    averaging_time_column: str | None = None
    observation_kind_column: str | None = None
    integrity_pinned: bool = False
    sensor_id_column: str = "sensor"
    x_column: str = "x_downwind_m"
    y_column: str = "y_crosswind_m"
    z_column: str = "height_m"
    time_column: str = "observation_time_s"
    common_clock_column: str = "common_clock_id"
    obstacle_geometry_column: str = "obstacle_geometry_id"

    def __post_init__(self) -> None:
        if not isinstance(self.path, str) or not self.path.strip():
            raise ValueError("CSV provenance path must be a non-empty string")
        if not isinstance(self.sha256, str) or len(self.sha256) != 64 or any(
            char not in "0123456789abcdef" for char in self.sha256
        ):
            raise ValueError("CSV provenance sha256 must be a lowercase SHA-256 digest")
        if isinstance(self.row_count, bool) or not isinstance(self.row_count, int) or self.row_count <= 0:
            raise ValueError("CSV provenance row_count must be a positive integer")
        if not isinstance(self.prediction_column, str) or not self.prediction_column.strip():
            raise ValueError("CSV provenance prediction_column must be non-empty")
        if self.concentration_unit not in {"mole_fraction", "volume_percent"}:
            raise ValueError("CSV provenance concentration_unit is unsupported")
        if self.averaging_time_column is not None and (
            not isinstance(self.averaging_time_column, str)
            or not self.averaging_time_column.strip()
        ):
            raise ValueError(
                "CSV provenance averaging_time_column must be non-empty when supplied"
            )
        if self.observation_kind_column is not None and (
            not isinstance(self.observation_kind_column, str)
            or not self.observation_kind_column.strip()
        ):
            raise ValueError(
                "CSV provenance observation_kind_column must be non-empty when supplied"
            )
        if not isinstance(self.integrity_pinned, bool):
            raise TypeError("CSV provenance integrity_pinned must be boolean")
        for name, value in {
            "sensor_id_column": self.sensor_id_column,
            "x_column": self.x_column,
            "y_column": self.y_column,
            "z_column": self.z_column,
            "time_column": self.time_column,
            "common_clock_column": self.common_clock_column,
            "obstacle_geometry_column": self.obstacle_geometry_column,
        }.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or any(character in value for character in "\r\n")
            ):
                raise ValueError(f"CSV provenance {name} must be a non-empty single-line string")


_FIELD_SCREENING_PREDICTION_OPERATORS = frozenset({
    "peak_indicated",
    "peak_true",
    "time_average_indicated",
    "time_average_true",
    "final_indicated",
    "final_true",
})


@dataclass(frozen=True)
class FieldModelExecutionProvenance:
    """Provenance for predictions extracted from a DEGALI field execution.

    A comparison basis identifies what two model paths claim to share.  This
    record separately preserves whether the originating field execution was
    accepted or conditional, and which explicit trace operator produced each
    fixed-sensor prediction.  A conditional execution therefore cannot look
    like an accepted model path merely because its basis IDs happen to match.
    """

    source_applicability: FieldApplicability
    prediction_operator: str

    def __post_init__(self) -> None:
        if not isinstance(self.source_applicability, FieldApplicability):
            raise TypeError(
                "source_applicability must be a FieldApplicability"
            )
        if not isinstance(self.prediction_operator, str):
            raise TypeError("prediction_operator must be a string")
        if self.prediction_operator not in _FIELD_SCREENING_PREDICTION_OPERATORS:
            raise ValueError(
                "prediction_operator must be one of "
                + ", ".join(sorted(_FIELD_SCREENING_PREDICTION_OPERATORS))
            )


@dataclass(frozen=True)
class FieldComparisonEvidence:
    """Fingerprint of a non-tabular comparison qualification artefact.

    Sensor results often omit information such as a common averaging window or
    source mapping.  This object can preserve a checked manifest without
    claiming that its text makes mismatched model operators equivalent.
    """

    path: str
    sha256: str
    qualification: str

    def __post_init__(self) -> None:
        if not isinstance(self.path, str) or not self.path.strip():
            raise ValueError("comparison evidence path must be a non-empty string")
        if not isinstance(self.sha256, str) or len(self.sha256) != 64 or any(
            char not in "0123456789abcdef" for char in self.sha256
        ):
            raise ValueError("comparison evidence sha256 must be a lowercase SHA-256 digest")
        if not isinstance(self.qualification, str) or not self.qualification.strip():
            raise ValueError("comparison evidence qualification must be non-empty")


@dataclass(frozen=True)
class FieldModelSensorSet:
    """One model's predictions under an explicitly named comparison basis."""

    model_id: str
    temporal_mode: str
    basis: FieldComparisonBasis
    predictions: tuple[FieldSensorPrediction, ...]
    runtime_s: float | None = None
    csv_provenance: FieldSensorCsvProvenance | None = None
    execution_provenance: FieldModelExecutionProvenance | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.model_id, str) or not self.model_id.strip():
            raise ValueError("model_id must be a non-empty string")
        if not isinstance(self.temporal_mode, str) or self.temporal_mode not in {
            "steady", "transient"
        }:
            raise ValueError("temporal_mode must be steady or transient")
        if not isinstance(self.basis, FieldComparisonBasis):
            raise TypeError("basis must be a FieldComparisonBasis")
        if self.csv_provenance is not None and not isinstance(
            self.csv_provenance, FieldSensorCsvProvenance
        ):
            raise TypeError(
                "csv_provenance must be FieldSensorCsvProvenance or None"
            )
        if self.execution_provenance is not None and not isinstance(
            self.execution_provenance, FieldModelExecutionProvenance
        ):
            raise TypeError(
                "execution_provenance must be FieldModelExecutionProvenance or None"
            )
        if not isinstance(self.predictions, tuple):
            raise TypeError("predictions must be a tuple")
        if not self.predictions:
            raise ValueError("at least one sensor prediction is required")
        if any(not isinstance(item, FieldSensorPrediction) for item in self.predictions):
            raise TypeError("predictions must be FieldSensorPrediction values")
        labels = [item.sensor_id for item in self.predictions]
        if len(labels) != len(set(labels)):
            raise ValueError("sensor prediction identifiers must be unique")
        if self.csv_provenance is not None and (
            self.csv_provenance.row_count != len(self.predictions)
        ):
            raise ValueError(
                "CSV provenance row_count must match the number of model predictions"
            )
        if isinstance(self.runtime_s, bool):
            raise TypeError("runtime_s must be numeric, not boolean")
        if self.runtime_s is not None and (
            not math.isfinite(float(self.runtime_s)) or self.runtime_s < 0.0
        ):
            raise ValueError("runtime_s must be finite and non-negative when supplied")


def field_model_sensor_set_from_csv(
    path: str | Path,
    *,
    model_id: str,
    temporal_mode: Literal["steady", "transient"],
    basis: FieldComparisonBasis,
    prediction_column: str,
    concentration_unit: Literal["mole_fraction", "volume_percent"],
    sensor_id_column: str = "sensor",
    x_column: str = "x_downwind_m",
    y_column: str = "y_crosswind_m",
    z_column: str = "height_m",
    averaging_time_column: str | None = "averaging_time_s",
    runtime_s: float | None = None,
) -> FieldModelSensorSet:
    """Read a fixed-sensor model result without losing its comparison contract.

    ``volume_percent`` is converted exactly to mole fraction by dividing by
    100.  If an averaging-time column and basis value are both supplied, every
    row must match that declared operator.  The function intentionally does
    not infer a source boundary or temporal mode from an external filename.
    """
    if concentration_unit not in {"mole_fraction", "volume_percent"}:
        raise ValueError("concentration_unit must be mole_fraction or volume_percent")
    column_roles = [sensor_id_column, x_column, y_column, z_column, prediction_column]
    if averaging_time_column is not None:
        column_roles.append(averaging_time_column)
    if any(
        not isinstance(name, str) or not name.strip() or any(char in name for char in "\r\n")
        for name in column_roles
    ):
        raise ValueError("sensor-prediction CSV column names must be non-empty and single-line")
    role_duplicates = sorted(
        name for name, count in Counter(column_roles).items() if count > 1
    )
    if role_duplicates:
        raise ValueError(
            "sensor-prediction CSV column roles must be distinct: "
            + ", ".join(role_duplicates)
        )
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    required_columns = {
        sensor_id_column, x_column, y_column, z_column, prediction_column,
    }
    if averaging_time_column is not None:
        required_columns.add(averaging_time_column)
    with source.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        raw_fieldnames = next(reader, None)
        fieldnames = tuple(raw_fieldnames or ())
        if not fieldnames:
            raise ValueError("sensor-prediction CSV must contain a header row")
        if any(not name.strip() for name in fieldnames):
            raise ValueError("sensor-prediction CSV header contains an empty column name")
        normalized = tuple(name.strip() for name in fieldnames)
        duplicates = sorted(
            name for name, count in Counter(normalized).items() if count > 1
        )
        if duplicates:
            raise ValueError(
                "sensor-prediction CSV header contains duplicate columns: "
                + ", ".join(duplicates)
            )
        available = set(fieldnames)
        missing = sorted(required_columns - available)
        if missing:
            raise ValueError(
                f"sensor-prediction CSV is missing required columns: {', '.join(missing)}"
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
                    f"invalid sensor-prediction CSV row {row_number}: {detail}; "
                    f"expected {len(fieldnames)} values, got {len(values)}"
                )
            rows.append(dict(zip(fieldnames, values)))
    if not rows:
        raise ValueError("sensor-prediction CSV contains no data rows")

    predictions = []
    for row_number, row in enumerate(rows, start=2):
        try:
            sensor_id = row[sensor_id_column].strip()
            position = tuple(float(row[name]) for name in (x_column, y_column, z_column))
            value = float(row[prediction_column])
            if concentration_unit == "volume_percent":
                if not 0.0 <= value <= 100.0:
                    raise ValueError("volume-percent concentration must lie in [0, 100]")
                value /= 100.0
            if averaging_time_column is not None and basis.averaging_time_s is not None:
                observed_average = float(row[averaging_time_column])
                if not math.isclose(
                    observed_average, basis.averaging_time_s,
                    rel_tol=1.0e-9, abs_tol=1.0e-9,
                ):
                    raise ValueError(
                        f"averaging time {observed_average:g} s does not match declared "
                        f"basis {basis.averaging_time_s:g} s"
                    )
            predictions.append(FieldSensorPrediction(sensor_id, position, value))
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid sensor-prediction CSV row {row_number}: {exc}") from exc

    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    provenance = FieldSensorCsvProvenance(
        path=str(source.resolve()), sha256=digest, row_count=len(rows),
        prediction_column=prediction_column, concentration_unit=concentration_unit,
        averaging_time_column=(
            averaging_time_column
            if averaging_time_column is not None and averaging_time_column in available
            else None
        ),
        sensor_id_column=sensor_id_column,
        x_column=x_column,
        y_column=y_column,
        z_column=z_column,
    )
    return FieldModelSensorSet(
        model_id=model_id, temporal_mode=temporal_mode, basis=basis,
        predictions=tuple(predictions), runtime_s=runtime_s,
        csv_provenance=provenance,
    )


def write_field_model_sensor_set_csv(
    model: FieldModelSensorSet,
    path: str | Path,
    *,
    prediction_column: str = "mole_fraction",
) -> FieldModelSensorSet:
    """Persist a typed model sensor set as a fingerprinted comparison CSV.

    The file contains fixed-sensor predictions and geometry only; comparison
    basis identifiers remain explicit on the returned typed model and in the
    comparison case.  Existing files are never overwritten.
    """
    if not isinstance(model, FieldModelSensorSet):
        raise TypeError("model must be a FieldModelSensorSet value")
    if not isinstance(prediction_column, str) or not prediction_column.strip():
        raise ValueError("prediction_column must be a non-empty string")
    if any(character in prediction_column for character in "\r\n"):
        raise ValueError("prediction_column must not contain newlines")
    if prediction_column in {
        "sensor", "x_downwind_m", "y_crosswind_m", "height_m", "averaging_time_s",
    }:
        raise ValueError(
            "prediction_column must be distinct from reserved geometry/operator columns"
        )
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(
            f"refusing to overwrite existing sensor-prediction CSV: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    averaging_time_column = (
        "averaging_time_s" if model.basis.averaging_time_s is not None else None
    )
    fieldnames = [
        "sensor", "x_downwind_m", "y_crosswind_m", "height_m", prediction_column,
    ]
    if averaging_time_column is not None:
        fieldnames.append(averaging_time_column)
    try:
        handle = destination.open("x", newline="", encoding="utf-8")
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing sensor-prediction CSV: {destination}"
        ) from error
    with handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for prediction in model.predictions:
            row = {
                "sensor": prediction.sensor_id,
                "x_downwind_m": prediction.position_m[0],
                "y_crosswind_m": prediction.position_m[1],
                "height_m": prediction.position_m[2],
                prediction_column: prediction.mole_fraction,
            }
            if averaging_time_column is not None:
                row[averaging_time_column] = model.basis.averaging_time_s
            writer.writerow(row)
    observed_digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    provenance = FieldSensorCsvProvenance(
        path=str(destination.resolve()),
        sha256=observed_digest,
        row_count=len(model.predictions),
        prediction_column=prediction_column,
        concentration_unit="mole_fraction",
        averaging_time_column=averaging_time_column,
        sensor_id_column="sensor",
        x_column="x_downwind_m",
        y_column="y_crosswind_m",
        z_column="height_m",
    )
    return replace(model, csv_provenance=provenance)


def _screening_trace_value(trace: object, prediction_operator: str) -> float:
    """Apply one explicit, auditable fixed-sensor trace operator."""
    indicated = prediction_operator.endswith("_indicated")
    signal_name = "indicated_mole_fraction" if indicated else "true_mole_fraction"
    values = np.asarray(getattr(trace, signal_name), dtype=float)
    time = np.asarray(getattr(trace, "time_s"), dtype=float)
    if values.ndim != 1 or values.size == 0 or time.shape != values.shape:
        raise ValueError("field screening trace arrays must be non-empty and aligned")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(time)):
        raise ValueError("field screening trace arrays must be finite")
    if prediction_operator.startswith("peak_"):
        value = float(np.max(values))
    elif prediction_operator.startswith("final_"):
        value = float(values[-1])
    elif prediction_operator.startswith("time_average_"):
        duration = float(time[-1] - time[0])
        if duration <= 0.0:
            value = float(values[-1])
        else:
            value = float(
                np.sum(0.5 * (values[1:] + values[:-1]) * np.diff(time))
                / duration
            )
    else:  # pragma: no cover - constructor validation owns this branch
        raise ValueError(f"unsupported field screening prediction operator {prediction_operator!r}")
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(
            f"field screening {prediction_operator} prediction must lie in [0, 1]"
        )
    return value


def field_model_sensor_set_from_screening(
    screening: "FieldSemiFVScreeningResult",
    *,
    model_id: str,
    basis: FieldComparisonBasis,
    prediction_operator: Literal[
        "peak_indicated", "peak_true", "time_average_indicated",
        "time_average_true", "final_indicated", "final_true",
    ] = "peak_indicated",
    sensor_labels: Sequence[str] | None = None,
    runtime_s: float | None = None,
) -> FieldModelSensorSet:
    """Convert completed field traces into matched fixed-sensor predictions.

    The trace operator is mandatory and never inferred from the model name.
    Withheld sensors and blocked executions fail closed; indicated values are
    not clipped if a declared gain/bias pushes them outside the physical
    mole-fraction range.
    """
    from .field_workflow import FieldSemiFVScreeningResult

    if not isinstance(screening, FieldSemiFVScreeningResult):
        raise TypeError("screening must be a FieldSemiFVScreeningResult")
    if not isinstance(basis, FieldComparisonBasis):
        raise TypeError("basis must be a FieldComparisonBasis")
    if prediction_operator not in _FIELD_SCREENING_PREDICTION_OPERATORS:
        raise ValueError(
            "prediction_operator must be one of "
            + ", ".join(sorted(_FIELD_SCREENING_PREDICTION_OPERATORS))
        )
    if not screening.completed:
        raise ValueError(
            "field screening must be completed and non-blocked before model comparison"
        )

    results = {item.label: item for item in screening.sensor_results}
    if not results and screening.sensor_trace is not None and screening.request.scenario.sensor is not None:
        results = {
            "field_sensor": FieldSensorDeploymentResult(
                "field_sensor", screening.request.scenario.sensor,
                screening.sensor_trace,
            ),
        }
    if sensor_labels is None:
        labels = tuple(results)
    else:
        if isinstance(sensor_labels, str):
            raise TypeError("sensor_labels must be a sequence of sensor identifiers")
        labels = tuple(sensor_labels)
        if not labels:
            raise ValueError("sensor_labels must contain at least one sensor")
        if any(not isinstance(label, str) or not label.strip() for label in labels):
            raise ValueError("sensor_labels must contain non-empty strings")
        if len(set(labels)) != len(labels):
            raise ValueError("sensor_labels must be unique")
    if not labels:
        raise ValueError("completed field screening contains no sensor traces")

    missing = [label for label in labels if label not in results]
    if missing:
        raise ValueError("requested field screening sensors are unavailable: " + ", ".join(missing))
    withheld = [label for label in labels if results[label].trace is None]
    if withheld:
        raise ValueError(
            "field screening sensors are withheld and cannot enter model comparison: "
            + ", ".join(withheld)
        )
    if basis.averaging_time_s is not None:
        for label in labels:
            actual = results[label].sensor.averaging_time_s
            if not math.isclose(
                actual, basis.averaging_time_s, rel_tol=1.0e-9, abs_tol=1.0e-9,
            ):
                raise ValueError(
                    f"sensor {label!r} averaging time {actual:g} s does not match "
                    f"comparison basis {basis.averaging_time_s:g} s"
                )

    predictions = tuple(
        FieldSensorPrediction(
            label,
            tuple(float(value) for value in results[label].sensor.position_m),
            _screening_trace_value(results[label].trace, prediction_operator),
        )
        for label in labels
    )
    return FieldModelSensorSet(
        model_id=model_id,
        temporal_mode=screening.request.scenario.temporal_mode,
        basis=basis,
        predictions=predictions,
        runtime_s=runtime_s,
        execution_provenance=FieldModelExecutionProvenance(
            screening.applicability, prediction_operator,
        ),
    )


def field_model_sensor_set_from_superposition(
    superposition: "FieldSensorSuperposition",
    *,
    model_id: str,
    basis: FieldComparisonBasis,
    prediction_operator: Literal[
        "peak_indicated", "peak_true", "time_average_indicated",
        "time_average_true", "final_indicated", "final_true",
    ] = "peak_indicated",
    sensor_id: str = "field_sensor",
    runtime_s: float | None = None,
) -> FieldModelSensorSet:
    """Export one common-sensor source superposition for model comparison.

    Superposition deliberately remains conditional because source-source
    interactions and three-dimensional cold-cloud mixing are unresolved.  The
    provenance is retained on the returned model set so matched basis IDs
    cannot promote the result to an accepted model-selection conclusion.
    """
    from .field_superposition import FieldSensorSuperposition

    if not isinstance(superposition, FieldSensorSuperposition):
        raise TypeError("superposition must be a FieldSensorSuperposition")
    if not isinstance(basis, FieldComparisonBasis):
        raise TypeError("basis must be a FieldComparisonBasis")
    if prediction_operator not in _FIELD_SCREENING_PREDICTION_OPERATORS:
        raise ValueError(
            "prediction_operator must be one of "
            + ", ".join(sorted(_FIELD_SCREENING_PREDICTION_OPERATORS))
        )
    if not superposition.completed:
        raise ValueError(
            "field sensor superposition must be completed and non-blocked before model comparison"
        )
    if not isinstance(sensor_id, str) or not sensor_id.strip():
        raise ValueError("sensor_id must be a non-empty string")
    if basis.averaging_time_s is not None and not math.isclose(
        superposition.sensor_trace.sensor.averaging_time_s,
        basis.averaging_time_s,
        rel_tol=1.0e-9,
        abs_tol=1.0e-9,
    ):
        actual = superposition.sensor_trace.sensor.averaging_time_s
        raise ValueError(
            f"sensor {sensor_id!r} averaging time {actual:g} s does not match "
            f"comparison basis {basis.averaging_time_s:g} s"
        )
    if not superposition.branches:
        raise ValueError("field sensor superposition contains no source branches")
    temporal_modes = {
        branch.result.request.scenario.temporal_mode
        for branch in superposition.branches
    }
    if len(temporal_modes) != 1:
        raise ValueError("field sensor superposition branches must share one temporal mode")
    trace = superposition.sensor_trace
    return FieldModelSensorSet(
        model_id=model_id,
        temporal_mode=next(iter(temporal_modes)),
        basis=basis,
        predictions=(FieldSensorPrediction(
            sensor_id,
            tuple(float(value) for value in trace.sensor.position_m),
            _screening_trace_value(trace, prediction_operator),
        ),),
        runtime_s=runtime_s,
        execution_provenance=FieldModelExecutionProvenance(
            superposition.applicability, prediction_operator,
        ),
    )


@dataclass(frozen=True)
class FieldModelComparisonRow:
    sensor_id: str
    position_m: tuple[float, float, float]
    left_mole_fraction: float
    right_mole_fraction: float
    difference_mole_fraction: float
    ratio_right_to_left: float | None
    left_above_threshold: bool
    right_above_threshold: bool

    def __post_init__(self) -> None:
        if not isinstance(self.sensor_id, str) or not self.sensor_id.strip():
            raise ValueError("comparison row sensor_id must be a non-empty string")
        if not isinstance(self.position_m, tuple):
            raise TypeError("comparison row position must be a tuple")
        if any(isinstance(value, bool) for value in self.position_m):
            raise TypeError("comparison row coordinates must be numeric, not boolean")
        if len(self.position_m) != 3 or not all(
            math.isfinite(float(value)) for value in self.position_m
        ):
            raise ValueError("comparison row position must contain three finite values")
        numeric = {
            "left_mole_fraction": self.left_mole_fraction,
            "right_mole_fraction": self.right_mole_fraction,
            "difference_mole_fraction": self.difference_mole_fraction,
        }
        if any(isinstance(value, bool) for value in numeric.values()):
            raise TypeError("comparison row concentrations must be numeric, not boolean")
        if any(
            not math.isfinite(float(value))
            for value in numeric.values()
        ):
            raise ValueError("comparison row concentrations must be finite")
        if not all(0.0 <= float(numeric[name]) <= 1.0 for name in (
            "left_mole_fraction", "right_mole_fraction",
        )):
            raise ValueError("comparison row mole fractions must lie in [0, 1]")
        expected_difference = float(self.right_mole_fraction) - float(self.left_mole_fraction)
        if not math.isclose(
            float(self.difference_mole_fraction), expected_difference,
            rel_tol=1.0e-9, abs_tol=1.0e-12,
        ):
            raise ValueError("comparison row difference must equal right minus left")
        if self.ratio_right_to_left is not None:
            if isinstance(self.ratio_right_to_left, bool):
                raise TypeError("comparison row ratio must be numeric, not boolean")
            if not math.isfinite(float(self.ratio_right_to_left)) or self.ratio_right_to_left < 0.0:
                raise ValueError("comparison row ratio must be finite and non-negative")
            if self.left_mole_fraction == 0.0:
                raise ValueError("comparison row ratio must be None when left prediction is zero")
            expected_ratio = float(self.right_mole_fraction) / float(self.left_mole_fraction)
            if not math.isclose(
                float(self.ratio_right_to_left), expected_ratio,
                rel_tol=1.0e-9, abs_tol=1.0e-12,
            ):
                raise ValueError("comparison row ratio does not match the two predictions")
        elif self.left_mole_fraction != 0.0:
            raise ValueError("comparison row ratio is required when left prediction is non-zero")
        if not isinstance(self.left_above_threshold, bool) or not isinstance(
            self.right_above_threshold, bool
        ):
            raise TypeError("comparison row threshold flags must be boolean")


@dataclass(frozen=True)
class FieldModelComparison:
    left: FieldModelSensorSet
    right: FieldModelSensorSet
    applicability: FieldApplicability
    rows: tuple[FieldModelComparisonRow, ...]
    threshold_mole_fraction: float
    mean_absolute_difference_mole_fraction: float | None
    classification_disagreement_count: int
    position_tolerance_m: float = 1.0e-9
    gate_codes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.left, FieldModelSensorSet) or not isinstance(
            self.right, FieldModelSensorSet
        ):
            raise TypeError("comparison models must be FieldModelSensorSet values")
        if self.left.model_id == self.right.model_id:
            raise ValueError("comparison requires two distinct model identifiers")
        if not isinstance(self.applicability, FieldApplicability):
            raise TypeError("comparison applicability must be FieldApplicability")
        if not isinstance(self.rows, tuple):
            raise TypeError("comparison rows must be a tuple")
        if any(not isinstance(row, FieldModelComparisonRow) for row in self.rows):
            raise TypeError("comparison rows must be FieldModelComparisonRow values")
        if len({row.sensor_id for row in self.rows}) != len(self.rows):
            raise ValueError("comparison row sensor identifiers must be unique")
        if isinstance(self.threshold_mole_fraction, bool) or not math.isfinite(
            float(self.threshold_mole_fraction)
        ) or not 0.0 <= float(self.threshold_mole_fraction) <= 1.0:
            raise ValueError("comparison threshold must be finite and lie in [0, 1]")
        if isinstance(self.position_tolerance_m, bool) or not math.isfinite(
            float(self.position_tolerance_m)
        ) or self.position_tolerance_m < 0.0:
            raise ValueError("comparison position tolerance must be finite and non-negative")
        if isinstance(self.classification_disagreement_count, bool) or not isinstance(
            self.classification_disagreement_count, int
        ) or self.classification_disagreement_count < 0:
            raise ValueError("classification disagreement count must be a non-negative integer")
        _validate_comparison_gate_codes(self.gate_codes, name="comparison gate_codes")

        if self.applicability.status != "blocked":
            left_by_id = {item.sensor_id: item for item in self.left.predictions}
            right_by_id = {item.sensor_id: item for item in self.right.predictions}
            row_ids = {row.sensor_id for row in self.rows}
            if row_ids != set(left_by_id) or row_ids != set(right_by_id):
                raise ValueError("non-blocked comparison rows must cover every matched sensor")
            for row in self.rows:
                left = left_by_id[row.sensor_id]
                right = right_by_id[row.sensor_id]
                if math.dist(left.position_m, right.position_m) > self.position_tolerance_m:
                    raise ValueError("comparison row contains a sensor geometry mismatch")
                if any(
                    not math.isclose(
                        float(actual), float(expected),
                        rel_tol=1.0e-9, abs_tol=1.0e-12,
                    )
                    for actual, expected in zip(row.position_m, left.position_m)
                ):
                    raise ValueError("comparison row position does not match the model sensor geometry")
                if not math.isclose(row.left_mole_fraction, left.mole_fraction, rel_tol=1.0e-9, abs_tol=1.0e-12):
                    raise ValueError("comparison row left prediction does not match the model sensor set")
                if not math.isclose(row.right_mole_fraction, right.mole_fraction, rel_tol=1.0e-9, abs_tol=1.0e-12):
                    raise ValueError("comparison row right prediction does not match the model sensor set")
            expected_disagreements = sum(
                row.left_above_threshold != row.right_above_threshold for row in self.rows
            )
            if self.classification_disagreement_count != expected_disagreements:
                raise ValueError("classification disagreement count does not match comparison rows")
            expected_mean = math.fsum(
                abs(row.difference_mole_fraction) for row in self.rows
            ) / len(self.rows)
            if self.mean_absolute_difference_mole_fraction is None or not math.isclose(
                float(self.mean_absolute_difference_mole_fraction), expected_mean,
                rel_tol=1.0e-9, abs_tol=1.0e-12,
            ):
                raise ValueError("mean absolute difference does not match comparison rows")
        else:
            if self.rows:
                raise ValueError("blocked comparison must not contain comparison rows")
            if self.classification_disagreement_count != 0:
                raise ValueError("blocked comparison cannot claim classification disagreements")
            if self.mean_absolute_difference_mole_fraction is not None:
                raise ValueError("blocked comparison mean difference must be None")


@dataclass(frozen=True)
class FieldModelDecisionImpact:
    """Decision-facing summary of a matched model comparison.

    The reported distance is the farthest *monitored* threshold exceedance on
    the supplied sensor line, never an interpolated hazard distance. Any
    model-form or operator mismatch withholds a model-selection conclusion.
    """

    comparison: FieldModelComparison
    status: Literal["aligned", "different", "withheld"]
    reasons: tuple[str, ...]
    left_alert_sensor_ids: tuple[str, ...]
    right_alert_sensor_ids: tuple[str, ...]
    left_farthest_monitored_exceedance_x_m: float | None
    right_farthest_monitored_exceedance_x_m: float | None
    runtime_ratio_right_to_left: float | None
    gate_codes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.comparison, FieldModelComparison):
            raise TypeError("model decision impact comparison must be a FieldModelComparison")
        if self.status not in {"aligned", "different", "withheld"}:
            raise ValueError("model decision impact status is unsupported")
        if not isinstance(self.reasons, tuple):
            raise TypeError("model decision impact reasons must be a tuple")
        if any(not isinstance(value, str) or not value.strip() for value in self.reasons):
            raise ValueError("model decision impact reasons must be non-empty strings")
        _validate_comparison_gate_codes(
            self.gate_codes, name="model decision impact gate_codes",
        )
        if self.status == "withheld" and not self.reasons:
            raise ValueError("withheld model decision impact requires a reason")
        for name, values in (
            ("left_alert_sensor_ids", self.left_alert_sensor_ids),
            ("right_alert_sensor_ids", self.right_alert_sensor_ids),
        ):
            if not isinstance(values, tuple):
                raise TypeError(f"{name} must be a tuple")
            if any(not isinstance(value, str) or not value.strip() for value in values):
                raise ValueError(f"{name} must contain non-empty strings")
            if len(set(values)) != len(values):
                raise ValueError(f"{name} must contain unique sensor identifiers")
        for name, value in (
            ("left_farthest_monitored_exceedance_x_m", self.left_farthest_monitored_exceedance_x_m),
            ("right_farthest_monitored_exceedance_x_m", self.right_farthest_monitored_exceedance_x_m),
        ):
            if value is not None and (
                isinstance(value, bool) or not math.isfinite(float(value))
            ):
                raise ValueError(f"{name} must be finite when supplied")
        if self.runtime_ratio_right_to_left is not None and (
            isinstance(self.runtime_ratio_right_to_left, bool)
            or not math.isfinite(float(self.runtime_ratio_right_to_left))
            or self.runtime_ratio_right_to_left < 0.0
        ):
            raise ValueError("runtime_ratio_right_to_left must be finite and non-negative")
        row_ids = {row.sensor_id for row in self.comparison.rows}
        if not set(self.left_alert_sensor_ids).issubset(row_ids) or not set(
            self.right_alert_sensor_ids
        ).issubset(row_ids):
            raise ValueError("model decision impact alert identifiers must match comparison rows")
        if self.comparison.applicability.status != "accepted":
            if (
                self.left_alert_sensor_ids or self.right_alert_sensor_ids
                or self.left_farthest_monitored_exceedance_x_m is not None
                or self.right_farthest_monitored_exceedance_x_m is not None
            ):
                raise ValueError(
                    "withheld model decision impact cannot carry alert or distance claims"
                )
        else:
            expected_left = tuple(
                row.sensor_id for row in self.comparison.rows if row.left_above_threshold
            )
            expected_right = tuple(
                row.sensor_id for row in self.comparison.rows if row.right_above_threshold
            )
            if self.status in {"aligned", "different"} and (
                set(self.left_alert_sensor_ids) != set(expected_left)
                or set(self.right_alert_sensor_ids) != set(expected_right)
            ):
                raise ValueError(
                    "model decision impact alert identifiers do not match the comparison classifications"
                )
            for name, values, distance in (
                (
                    "left",
                    self.left_alert_sensor_ids,
                    self.left_farthest_monitored_exceedance_x_m,
                ),
                (
                    "right",
                    self.right_alert_sensor_ids,
                    self.right_farthest_monitored_exceedance_x_m,
                ),
            ):
                expected_distance = max(
                    (
                        row.position_m[0]
                        for row in self.comparison.rows
                        if row.sensor_id in values
                    ),
                    default=None,
                )
                if expected_distance is None:
                    if distance is not None:
                        raise ValueError(
                            f"{name} monitored distance requires an alert sensor"
                        )
                elif distance is None or not math.isclose(
                    float(distance), expected_distance,
                    rel_tol=1.0e-9, abs_tol=1.0e-12,
                ):
                    raise ValueError(
                        f"{name} monitored distance does not match alert sensors"
                    )
        if self.comparison.applicability.status != "accepted":
            if self.status != "withheld":
                raise ValueError("non-accepted model comparison impact must be withheld")
        elif self.status == "withheld":
            raise ValueError("accepted model comparison cannot have a withheld impact")
        elif self.status == "different" and self.comparison.classification_disagreement_count == 0:
            raise ValueError("different model impact requires a classification disagreement")
        elif self.status == "aligned" and self.comparison.classification_disagreement_count != 0:
            raise ValueError("aligned model impact requires no classification disagreement")


def field_model_decision_impact(
    comparison: FieldModelComparison,
) -> FieldModelDecisionImpact:
    """Translate matched sensor classifications into a fail-safe decision summary.

    ``accepted`` comparison basis is required before a disagreement is called
    a model-selection difference. Conditional comparisons, including a
    DEGALI–SLABx steady/transient operator mismatch, retain their numerical
    rows but return ``withheld`` rather than suggesting an operational model
    preference.
    """
    if not isinstance(comparison, FieldModelComparison):
        raise TypeError("comparison must be a FieldModelComparison")
    left_runtime = comparison.left.runtime_s
    right_runtime = comparison.right.runtime_s
    runtime_ratio = (
        None if left_runtime is None or right_runtime is None or left_runtime == 0.0
        else right_runtime / left_runtime
    )
    if comparison.applicability.status != "accepted":
        reasons = comparison.applicability.reasons or comparison.applicability.warnings
        gate_codes = tuple(dict.fromkeys(
            comparison.gate_codes + ("model_selection_withheld",)
        ))
        return FieldModelDecisionImpact(
            comparison, "withheld",
            tuple(reasons) + (
                "model-selection impact is withheld until source, weather, sensor geometry and temporal operator are matched",
            ),
            (), (), None, None, runtime_ratio,
            gate_codes,
        )
    left_alerts = tuple(row.sensor_id for row in comparison.rows if row.left_above_threshold)
    right_alerts = tuple(row.sensor_id for row in comparison.rows if row.right_above_threshold)
    left_distance = max(
        (row.position_m[0] for row in comparison.rows if row.left_above_threshold),
        default=None,
    )
    right_distance = max(
        (row.position_m[0] for row in comparison.rows if row.right_above_threshold),
        default=None,
    )
    different = comparison.classification_disagreement_count > 0
    return FieldModelDecisionImpact(
        comparison,
        "different" if different else "aligned",
        (
            ("matched sensor threshold classifications differ between model paths",)
            if different else
            ("matched sensor threshold classifications agree; this does not establish accuracy outside monitored locations",)
        ),
        left_alerts, right_alerts, left_distance, right_distance, runtime_ratio,
        ("classification_disagreement",) if different else ("classification_aligned",),
    )


def field_model_decision_impact_record(
    impact: FieldModelDecisionImpact,
) -> dict[str, object]:
    """Return a JSON-native operational summary without hiding its holdback."""
    if not isinstance(impact, FieldModelDecisionImpact):
        raise TypeError("impact must be a FieldModelDecisionImpact")
    return {
        "status": impact.status,
        "reasons": list(impact.reasons),
        "threshold_mole_fraction": impact.comparison.threshold_mole_fraction,
        "left_model_id": impact.comparison.left.model_id,
        "right_model_id": impact.comparison.right.model_id,
        "left_alert_sensor_ids": list(impact.left_alert_sensor_ids),
        "right_alert_sensor_ids": list(impact.right_alert_sensor_ids),
        "left_farthest_monitored_exceedance_x_m": impact.left_farthest_monitored_exceedance_x_m,
        "right_farthest_monitored_exceedance_x_m": impact.right_farthest_monitored_exceedance_x_m,
        "runtime_ratio_right_to_left": impact.runtime_ratio_right_to_left,
        "gate_codes": list(impact.gate_codes),
        "qualification": (
            "farthest monitored threshold exceedance only; no unmeasured spatial interpolation or hazard-distance claim"
        ),
    }


FIELD_MODEL_COMPARISON_REPORT_SCHEMA = "degali.field-model-comparison-report.v1"


def _comparison_basis_record(basis: FieldComparisonBasis) -> dict[str, object]:
    """Return the declared operator rather than inferring it from a model name."""
    return {
        "source_boundary_id": basis.source_boundary_id,
        "weather_id": basis.weather_id,
        "sensor_geometry_id": basis.sensor_geometry_id,
        "temporal_operator_id": basis.temporal_operator_id,
        "averaging_time_s": basis.averaging_time_s,
        "obstacle_representation_id": basis.obstacle_representation_id,
    }


def _model_sensor_set_record(model: FieldModelSensorSet) -> dict[str, object]:
    provenance = model.csv_provenance
    execution = model.execution_provenance
    return {
        "model_id": model.model_id,
        "temporal_mode": model.temporal_mode,
        "runtime_s": model.runtime_s,
        "comparison_basis": _comparison_basis_record(model.basis),
        "csv_provenance": None if provenance is None else {
            "path": provenance.path,
            "sha256": provenance.sha256,
            "row_count": provenance.row_count,
            "prediction_column": provenance.prediction_column,
            "concentration_unit": provenance.concentration_unit,
            "averaging_time_column": provenance.averaging_time_column,
            "observation_kind_column": provenance.observation_kind_column,
            "integrity_pinned": provenance.integrity_pinned,
            "sensor_id_column": provenance.sensor_id_column,
            "x_column": provenance.x_column,
            "y_column": provenance.y_column,
            "z_column": provenance.z_column,
            "time_column": provenance.time_column,
            "common_clock_column": provenance.common_clock_column,
            "obstacle_geometry_column": provenance.obstacle_geometry_column,
        },
        "execution_provenance": None if execution is None else {
            "prediction_operator": execution.prediction_operator,
            "source_applicability": {
                "status": execution.source_applicability.status,
                "reasons": list(execution.source_applicability.reasons),
                "warnings": list(execution.source_applicability.warnings),
                "uncertainty_complete": execution.source_applicability.uncertainty_complete,
            },
        },
    }


def field_model_comparison_report(
    comparison: FieldModelComparison,
    *,
    comparison_evidence: FieldComparisonEvidence | None = None,
) -> dict[str, object]:
    """Serialise a supplied-model comparison with all decision holdbacks.

    The report is deliberately a comparison of declared sensor predictions.
    It does not claim that an imported CSV represents a validated SLABx,
    CFD, PHAST or DEGALI execution; that provenance must be supplied outside
    this narrow input boundary.
    """
    if not isinstance(comparison, FieldModelComparison):
        raise TypeError("comparison must be a FieldModelComparison")
    if comparison_evidence is not None and not isinstance(
        comparison_evidence, FieldComparisonEvidence
    ):
        raise TypeError("comparison_evidence must be a FieldComparisonEvidence or None")
    impact = field_model_decision_impact(comparison)
    return {
        "schema": FIELD_MODEL_COMPARISON_REPORT_SCHEMA,
        "left": _model_sensor_set_record(comparison.left),
        "right": _model_sensor_set_record(comparison.right),
        "comparison_evidence": None if comparison_evidence is None else {
            "path": comparison_evidence.path,
            "sha256": comparison_evidence.sha256,
            "qualification": comparison_evidence.qualification,
        },
        "comparison": {
            "gate_codes": list(comparison.gate_codes),
            "applicability": {
                "status": comparison.applicability.status,
                "reasons": list(comparison.applicability.reasons),
                "warnings": list(comparison.applicability.warnings),
                "uncertainty_complete": comparison.applicability.uncertainty_complete,
                "gate_codes": list(comparison.gate_codes),
            },
            "threshold_mole_fraction": comparison.threshold_mole_fraction,
            "position_tolerance_m": comparison.position_tolerance_m,
            "mean_absolute_difference_mole_fraction": (
                comparison.mean_absolute_difference_mole_fraction
            ),
            "classification_disagreement_count": comparison.classification_disagreement_count,
            "rows": [
                {
                    "sensor_id": row.sensor_id,
                    "position_m": list(row.position_m),
                    "left_mole_fraction": row.left_mole_fraction,
                    "right_mole_fraction": row.right_mole_fraction,
                    "difference_mole_fraction": row.difference_mole_fraction,
                    "ratio_right_to_left": row.ratio_right_to_left,
                    "left_above_threshold": row.left_above_threshold,
                    "right_above_threshold": row.right_above_threshold,
                }
                for row in comparison.rows
            ],
        },
        "model_selection_impact": field_model_decision_impact_record(impact),
        "scope": (
            "supplied fixed-sensor predictions only; this is not an independent "
            "accuracy validation, a spatial hazard-distance interpolation, or an "
            "approval/design-basis decision"
        ),
    }


def _basis_mismatches(
    left: FieldComparisonBasis, right: FieldComparisonBasis,
) -> tuple[str, ...]:
    labels = {
        "source boundary": (left.source_boundary_id, right.source_boundary_id),
        "weather": (left.weather_id, right.weather_id),
        "sensor geometry": (left.sensor_geometry_id, right.sensor_geometry_id),
        "temporal operator": (left.temporal_operator_id, right.temporal_operator_id),
        "averaging time": (left.averaging_time_s, right.averaging_time_s),
        "obstacle representation": (
            left.obstacle_representation_id, right.obstacle_representation_id,
        ),
    }
    return tuple(
        f"{name} differs between model paths ({first!r} versus {second!r})"
        for name, (first, second) in labels.items() if first != second
    )


def _averaging_provenance_warnings(
    left: FieldModelSensorSet, right: FieldModelSensorSet,
) -> tuple[str, ...]:
    """Keep CSV averaging evidence aligned with the declared comparison basis."""
    warnings: list[str] = []
    for model in (left, right):
        provenance = model.csv_provenance
        if provenance is None:
            continue
        label = model.model_id
        if model.basis.averaging_time_s is not None and provenance.averaging_time_column is None:
            warnings.append(
                f"model {label!r} CSV has no averaging-time column to verify the declared "
                f"{model.basis.averaging_time_s:g} s operator"
            )
        elif model.basis.averaging_time_s is None and provenance.averaging_time_column is not None:
            warnings.append(
                f"model {label!r} CSV declares an averaging-time column but the comparison "
                "basis omits averaging_time_s"
            )
    return tuple(warnings)


def _execution_provenance_state(
    left: FieldModelSensorSet, right: FieldModelSensorSet,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return blocking reasons and conditional warnings from model executions."""
    reasons: list[str] = []
    warnings: list[str] = []
    for model in (left, right):
        provenance = model.execution_provenance
        if provenance is None:
            continue
        applicability = provenance.source_applicability
        label = f"model {model.model_id!r}"
        if applicability.status == "blocked":
            reasons.extend(
                f"{label} field execution was blocked: {reason}"
                for reason in applicability.reasons
            )
        elif applicability.status == "conditional" or not applicability.uncertainty_complete:
            warnings.append(
                f"{label} field execution is conditional or unresolved; model-selection impact remains withheld"
            )
            warnings.extend(
                f"{label} execution warning: {warning}"
                for warning in applicability.warnings
            )
    return tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(warnings))


def compare_field_model_sensor_sets(
    left: FieldModelSensorSet,
    right: FieldModelSensorSet,
    *,
    threshold_mole_fraction: float = 0.04,
    position_tolerance_m: float = 1.0e-9,
) -> FieldModelComparison:
    """Compare two sets only after checking the actual comparison operator."""
    if not isinstance(left, FieldModelSensorSet) or not isinstance(
        right, FieldModelSensorSet
    ):
        raise TypeError("left and right must be FieldModelSensorSet values")
    if isinstance(threshold_mole_fraction, bool) or not math.isfinite(
        float(threshold_mole_fraction)
    ) or not 0.0 <= float(threshold_mole_fraction) <= 1.0:
        raise ValueError("threshold_mole_fraction must be finite and lie in [0, 1]")
    if isinstance(position_tolerance_m, bool) or not math.isfinite(
        float(position_tolerance_m)
    ) or position_tolerance_m < 0.0:
        raise ValueError("position_tolerance_m must be finite and non-negative")
    if left.model_id == right.model_id:
        raise ValueError("comparison requires two distinct model identifiers")
    left_by_id = {item.sensor_id: item for item in left.predictions}
    right_by_id = {item.sensor_id: item for item in right.predictions}
    if set(left_by_id) != set(right_by_id):
        reasons = ("model paths do not contain the same sensor identifiers",)
        return FieldModelComparison(
            left, right, FieldApplicability("blocked", reasons), (),
            threshold_mole_fraction, None, 0, position_tolerance_m,
            ("sensor_set_mismatch",),
        )
    position_errors = []
    for sensor_id in left_by_id:
        distance = math.dist(left_by_id[sensor_id].position_m, right_by_id[sensor_id].position_m)
        if distance > position_tolerance_m:
            position_errors.append(
                f"sensor {sensor_id!r} differs by {distance:g} m between model paths"
            )
    if position_errors:
        return FieldModelComparison(
            left, right, FieldApplicability("blocked", tuple(position_errors)), (),
            threshold_mole_fraction, None, 0, position_tolerance_m,
            ("sensor_geometry_mismatch",),
        )
    mismatches = _basis_mismatches(left.basis, right.basis)
    warnings = list(mismatches)
    averaging_warnings = _averaging_provenance_warnings(left, right)
    warnings.extend(averaging_warnings)
    execution_reasons, execution_warnings = _execution_provenance_state(left, right)
    if execution_reasons:
        return FieldModelComparison(
            left, right, FieldApplicability("blocked", execution_reasons), (),
            threshold_mole_fraction, None, 0, position_tolerance_m,
            ("execution_blocked",),
        )
    warnings.extend(execution_warnings)
    gate_codes: list[str] = []
    if mismatches:
        gate_codes.append("basis_mismatch")
    if averaging_warnings:
        gate_codes.append("averaging_operator_unverified")
    execution_provenances = tuple(
        provenance
        for model in (left, right)
        if (provenance := model.execution_provenance) is not None
    )
    if execution_warnings:
        gate_codes.append("execution_conditional")
    if any(
        not provenance.source_applicability.uncertainty_complete
        for provenance in execution_provenances
    ):
        gate_codes.append("execution_uncertainty_unresolved")
    if left.temporal_mode != right.temporal_mode:
        warnings.append(
            "steady/transient model-form comparison; values are not equivalent unless the temporal operator matches"
        )
        gate_codes.append("temporal_mode_mismatch")
    if not warnings:
        gate_codes.append("matched_basis")
    applicability = FieldApplicability(
        "conditional" if warnings else "accepted",
        warnings=tuple(warnings),
        uncertainty_complete=not warnings,
    )
    rows = []
    for sensor_id in sorted(left_by_id):
        first, second = left_by_id[sensor_id], right_by_id[sensor_id]
        ratio = None if first.mole_fraction == 0.0 else second.mole_fraction / first.mole_fraction
        rows.append(FieldModelComparisonRow(
            sensor_id=sensor_id,
            position_m=first.position_m,
            left_mole_fraction=first.mole_fraction,
            right_mole_fraction=second.mole_fraction,
            difference_mole_fraction=second.mole_fraction - first.mole_fraction,
            ratio_right_to_left=ratio,
            left_above_threshold=first.mole_fraction >= threshold_mole_fraction,
            right_above_threshold=second.mole_fraction >= threshold_mole_fraction,
        ))
    differences = np.asarray([abs(item.difference_mole_fraction) for item in rows])
    disagreements = sum(
        item.left_above_threshold != item.right_above_threshold for item in rows
    )
    return FieldModelComparison(
        left, right, applicability, tuple(rows), threshold_mole_fraction,
        float(np.mean(differences)), disagreements, position_tolerance_m,
        tuple(dict.fromkeys(gate_codes)),
    )


__all__ = [
    "FieldComparisonBasis", "FieldSensorPrediction", "FieldSensorCsvProvenance",
    "FieldComparisonEvidence", "FieldModelExecutionProvenance",
    "FieldModelSensorSet", "field_model_sensor_set_from_csv",
    "write_field_model_sensor_set_csv",
    "field_model_sensor_set_from_screening",
    "field_model_sensor_set_from_superposition",
    "FieldModelComparisonRow", "FieldModelComparison",
    "FieldModelDecisionImpact", "compare_field_model_sensor_sets",
    "field_model_decision_impact", "field_model_decision_impact_record",
    "FIELD_MODEL_COMPARISON_GATE_CODES", "FIELD_MODEL_COMPARISON_REPORT_SCHEMA",
    "field_model_comparison_report",
]
