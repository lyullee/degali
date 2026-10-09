"""Deterministic FFI source-state sensitivity and residual maps.

This module keeps the validated free-field ``hydrogen_jet`` path unchanged.
It adds an explicit audit boundary around the public FFI/Spadeadam outdoor
horizontal-release records: source, wind and ambient bounds are crossed as
deterministic corners, the model is evaluated at the reported sensor
locations, and the resulting residual map retains the fact that arc/sensor
peaks are lower bounds because the array can miss the plume between sensors.

The result is a sensitivity envelope, not a confidence interval or an
operational validation.  No coefficient is fitted and no missing sensor value
is replaced by zero.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
import hashlib
import math
from typing import Literal
from pathlib import Path

from ..addons.field_contracts import BoundedValue, CircularBoundedValue
from .nearfield import REACH, STEP, Trajectory, hydrogen_jet
from .spadeadam import RELEASE_HEIGHT, Trial, Reading


FFI_SOURCE_STATE_ENVELOPE_SCHEMA = "degali.ffi-source-state-envelope.v1"


def ffi_reference_provenance(root: str | Path) -> dict[str, object]:
    """Fingerprint the public FFI tables used by the source-state audit.

    Only metadata is returned.  The extracted tables remain at ``root`` and
    are never copied into an execution report.
    """
    directory = Path(root).resolve()
    files = {}
    for name in ("conditions.csv", "sensors.csv"):
        path = directory / name
        if not path.is_file():
            raise FileNotFoundError(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        files[name] = {
            "path": str(path),
            "sha256": digest,
            "size_bytes": path.stat().st_size,
        }
    return {"root": str(directory), "files": files}


def _declared_bound(
    value: object,
    name: str,
    *,
    circular: bool = False,
) -> BoundedValue | CircularBoundedValue:
    expected = CircularBoundedValue if circular else BoundedValue
    if not isinstance(value, expected):
        label = "CircularBoundedValue" if circular else "BoundedValue"
        raise TypeError(f"{name} must be a {label}")
    source = value.source
    if not isinstance(source, str) or not source.strip() or source.strip().lower() == "unspecified":
        raise ValueError(f"{name} requires an explicit evidence source")
    return value


def _positive(value: BoundedValue | CircularBoundedValue, name: str, *, allow_zero: bool = False) -> None:
    lower = float(value.lower)
    if not math.isfinite(lower) or (lower < 0.0 if allow_zero else lower <= 0.0):
        relation = "non-negative" if allow_zero else "positive"
        raise ValueError(f"{name} lower bound must be finite and {relation}")


def _exact(value: BoundedValue, selected: float) -> BoundedValue:
    return BoundedValue(
        selected, selected, selected, unit=value.unit, source=value.source,
    )


def _exact_direction(value: CircularBoundedValue, selected: float) -> CircularBoundedValue:
    return CircularBoundedValue(
        selected, selected, selected, unit=value.unit, source=value.source,
    )


def _choices(value: BoundedValue | CircularBoundedValue) -> tuple[float, ...]:
    if value.is_exact:
        return (float(value.nominal),)
    return (float(value.lower), float(value.upper))


def _validate_selection(selection: object) -> None:
    if not isinstance(selection, tuple):
        raise TypeError("FFI source-state selection must be a tuple")
    seen: set[str] = set()
    for item in selection:
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError("FFI source-state selection entries must be pairs")
        name, value = item
        if not isinstance(name, str) or not name.strip() or name in seen:
            raise ValueError("FFI source-state selection names must be unique strings")
        if isinstance(value, bool) or not math.isfinite(float(value)):
            raise ValueError("FFI source-state selection values must be finite numbers")
        seen.add(name)


@dataclass(frozen=True)
class FfiSourceState:
    """Source, meteorology and ambient bounds for one FFI outdoor release.

    The bounds are deterministic engineering endpoints.  Each bound carries
    its own evidence string, and the model never assigns weights to corners.
    ``storage_pressure_barg`` is retained because it is part of the free-field
    source construction even though the public FFI record reports it as gauge
    pressure.
    """

    rate_kg_s: BoundedValue
    orifice_m: BoundedValue
    release_height_m: BoundedValue
    storage_pressure_barg: BoundedValue
    wind_m_s: BoundedValue
    wind_direction_from_deg: CircularBoundedValue
    ambient_temperature_k: BoundedValue
    relative_humidity_pct: BoundedValue
    ambient_pressure_pa: BoundedValue
    wind_reference_height_m: BoundedValue

    def __post_init__(self) -> None:
        fields = {
            "rate_kg_s": self.rate_kg_s,
            "orifice_m": self.orifice_m,
            "release_height_m": self.release_height_m,
            "storage_pressure_barg": self.storage_pressure_barg,
            "wind_m_s": self.wind_m_s,
            "wind_direction_from_deg": self.wind_direction_from_deg,
            "ambient_temperature_k": self.ambient_temperature_k,
            "relative_humidity_pct": self.relative_humidity_pct,
            "ambient_pressure_pa": self.ambient_pressure_pa,
            "wind_reference_height_m": self.wind_reference_height_m,
        }
        for name, value in fields.items():
            _declared_bound(
                value, name, circular=(name == "wind_direction_from_deg")
            )
        _positive(self.rate_kg_s, "rate_kg_s")
        _positive(self.orifice_m, "orifice_m")
        _positive(self.release_height_m, "release_height_m")
        _positive(self.storage_pressure_barg, "storage_pressure_barg", allow_zero=True)
        _positive(self.wind_m_s, "wind_m_s")
        _positive(self.ambient_temperature_k, "ambient_temperature_k")
        _positive(self.ambient_pressure_pa, "ambient_pressure_pa")
        _positive(self.wind_reference_height_m, "wind_reference_height_m")
        if self.relative_humidity_pct.lower < 0.0 or self.relative_humidity_pct.upper > 100.0:
            raise ValueError("relative_humidity_pct must lie in [0, 100]")

    @classmethod
    def from_trial(
        cls,
        trial: Trial,
        *,
        wind: float | None = None,
        ambient_temperature_k: float = 277.15,
        relative_humidity_pct: float = 90.0,
        ambient_pressure_pa: float = 101325.0,
        wind_reference_height_m: float = 1.5,
    ) -> "FfiSourceState":
        """Build an exact, explicitly labelled boundary from an FFI trial."""
        if not isinstance(trial, Trial):
            raise TypeError("trial must be a spadeadam Trial")
        selected_wind = trial.wind_low if wind is None else float(wind)
        return cls(
            rate_kg_s=BoundedValue(
                trial.rate, unit="kg/s", source="reference/spadeadam/conditions.csv",
            ),
            orifice_m=BoundedValue(
                trial.orifice, unit="m", source="reference/spadeadam/conditions.csv",
            ),
            release_height_m=BoundedValue(
                RELEASE_HEIGHT, unit="m", source="FFI report horizontal-release geometry",
            ),
            storage_pressure_barg=BoundedValue(
                trial.line_pressure, unit="barg", source="reference/spadeadam/conditions.csv",
            ),
            wind_m_s=BoundedValue(
                selected_wind, unit="m/s", source="reference/spadeadam/conditions.csv:low mast",
            ),
            wind_direction_from_deg=CircularBoundedValue(
                trial.wind_direction, unit="deg", source="reference/spadeadam/conditions.csv",
            ),
            ambient_temperature_k=BoundedValue(
                ambient_temperature_k, unit="K", source="FFI winter ambient assumption",
            ),
            relative_humidity_pct=BoundedValue(
                relative_humidity_pct, unit="%", source="FFI winter ambient assumption",
            ),
            ambient_pressure_pa=BoundedValue(
                ambient_pressure_pa, unit="Pa", source="standard ambient pressure assumption",
            ),
            wind_reference_height_m=BoundedValue(
                wind_reference_height_m, unit="m", source="FFI mast reference height assumption",
            ),
        )

    def uncertainty_fields(self) -> dict[str, BoundedValue | CircularBoundedValue]:
        return {
            "rate_kg_s": self.rate_kg_s,
            "orifice_m": self.orifice_m,
            "release_height_m": self.release_height_m,
            "storage_pressure_barg": self.storage_pressure_barg,
            "wind_m_s": self.wind_m_s,
            "wind_direction_from_deg": self.wind_direction_from_deg,
            "ambient_temperature_k": self.ambient_temperature_k,
            "relative_humidity_pct": self.relative_humidity_pct,
            "ambient_pressure_pa": self.ambient_pressure_pa,
            "wind_reference_height_m": self.wind_reference_height_m,
        }

    def corner_cases(self, *, max_cases: int = 64) -> tuple[tuple[tuple[str, float], "FfiSourceState"], ...]:
        if isinstance(max_cases, bool) or not isinstance(max_cases, int) or max_cases <= 0:
            raise ValueError("max_cases must be a positive integer")
        fields = self.uncertainty_fields()
        names = tuple(fields)
        choices = tuple(_choices(fields[name]) for name in names)
        count = math.prod(len(values) for values in choices)
        if count > max_cases:
            raise ValueError(
                f"{count} FFI source-state corners exceed max_cases={max_cases}"
            )
        out = []
        for combination in product(*choices):
            selected = dict(zip(names, combination))
            exact = replace(
                self,
                rate_kg_s=_exact(self.rate_kg_s, selected["rate_kg_s"]),
                orifice_m=_exact(self.orifice_m, selected["orifice_m"]),
                release_height_m=_exact(self.release_height_m, selected["release_height_m"]),
                storage_pressure_barg=_exact(self.storage_pressure_barg, selected["storage_pressure_barg"]),
                wind_m_s=_exact(self.wind_m_s, selected["wind_m_s"]),
                wind_direction_from_deg=_exact_direction(
                    self.wind_direction_from_deg, selected["wind_direction_from_deg"]
                ),
                ambient_temperature_k=_exact(self.ambient_temperature_k, selected["ambient_temperature_k"]),
                relative_humidity_pct=_exact(self.relative_humidity_pct, selected["relative_humidity_pct"]),
                ambient_pressure_pa=_exact(self.ambient_pressure_pa, selected["ambient_pressure_pa"]),
                wind_reference_height_m=_exact(self.wind_reference_height_m, selected["wind_reference_height_m"]),
            )
            out.append((tuple((name, float(selected[name])) for name in names), exact))
        return tuple(out)

    def as_record(self) -> dict[str, object]:
        return {
            name: value.as_dict()
            for name, value in self.uncertainty_fields().items()
        }


@dataclass(frozen=True)
class FfiResidualRow:
    """One reported sensor and its model value in one source-state corner."""

    sensor: str
    radius_m: float
    bearing_deg: float
    height_m: float
    alongwind_m: float
    crosswind_m: float
    observed_peak_vol_pct: float
    observed_is_lower_bound: bool
    predicted_vol_pct: float | None
    residual_vol_pct: float | None
    withheld_reason: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.sensor, str) or not self.sensor.strip():
            raise ValueError("FFI residual sensor must be a non-empty string")
        for name in (
            "radius_m", "bearing_deg", "height_m", "alongwind_m", "crosswind_m",
            "observed_peak_vol_pct",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"FFI residual {name} must be finite")
        if self.radius_m < 0.0 or self.height_m < 0.0 or self.observed_peak_vol_pct < 0.0:
            raise ValueError("FFI residual geometry and observation must be non-negative")
        if not isinstance(self.observed_is_lower_bound, bool):
            raise TypeError("observed_is_lower_bound must be boolean")
        if self.predicted_vol_pct is None:
            if self.residual_vol_pct is not None:
                raise ValueError("withheld FFI residual cannot have a residual value")
            if not isinstance(self.withheld_reason, str) or not self.withheld_reason.strip():
                raise ValueError("withheld FFI residual requires a reason")
        else:
            if not math.isfinite(float(self.predicted_vol_pct)) or self.predicted_vol_pct < 0.0:
                raise ValueError("predicted_vol_pct must be finite and non-negative")
            if self.residual_vol_pct is None or not math.isfinite(float(self.residual_vol_pct)):
                raise ValueError("available FFI residual requires a finite residual")
            expected = float(self.predicted_vol_pct) - float(self.observed_peak_vol_pct)
            if not math.isclose(float(self.residual_vol_pct), expected, rel_tol=1.0e-10, abs_tol=1.0e-12):
                raise ValueError("FFI residual_vol_pct is inconsistent with prediction and observation")
            if self.withheld_reason is not None:
                raise ValueError("available FFI residual cannot carry a withheld reason")

    def as_record(self) -> dict[str, object]:
        return {
            "sensor": self.sensor,
            "radius_m": self.radius_m,
            "bearing_deg": self.bearing_deg,
            "height_m": self.height_m,
            "alongwind_m": self.alongwind_m,
            "crosswind_m": self.crosswind_m,
            "observed_peak_vol_pct": self.observed_peak_vol_pct,
            "observed_is_lower_bound": self.observed_is_lower_bound,
            "predicted_vol_pct": self.predicted_vol_pct,
            "residual_vol_pct": self.residual_vol_pct,
            "withheld_reason": self.withheld_reason,
        }


def _case_metrics(rows: tuple[FfiResidualRow, ...]) -> dict[str, object]:
    available = [row for row in rows if row.predicted_vol_pct is not None]
    if not available:
        return {
            "available_rows": 0,
            "withheld_rows": len(rows),
            "mae_vol_pct_point": None,
            "rmse_vol_pct_point": None,
            "bias_vol_pct_point": None,
            "max_predicted_vol_pct": None,
            "max_observed_peak_vol_pct": (
                None if not rows else max(row.observed_peak_vol_pct for row in rows)
            ),
            "lower_bound_rows": 0,
            "lower_bound_satisfied_rows": 0,
            "lower_bound_satisfaction_fraction": None,
            "mean_lower_bound_deficit_vol_pct_point": None,
            "max_lower_bound_deficit_vol_pct_point": None,
            "lower_bound_constraint_status": "withheld",
        }
    residuals = [float(row.residual_vol_pct) for row in available]
    lower_bound_margins = [
        float(row.predicted_vol_pct) - float(row.observed_peak_vol_pct)
        for row in available if row.observed_is_lower_bound
    ]
    lower_bound_satisfied = sum(
        margin >= -1.0e-12 for margin in lower_bound_margins
    )
    deficits = [max(-margin, 0.0) for margin in lower_bound_margins]
    return {
        "available_rows": len(available),
        "withheld_rows": len(rows) - len(available),
        "mae_vol_pct_point": sum(abs(value) for value in residuals) / len(residuals),
        "rmse_vol_pct_point": math.sqrt(sum(value * value for value in residuals) / len(residuals)),
        "bias_vol_pct_point": sum(residuals) / len(residuals),
        "max_predicted_vol_pct": max(row.predicted_vol_pct for row in available),
        "max_observed_peak_vol_pct": max(row.observed_peak_vol_pct for row in rows),
        "lower_bound_rows": len(lower_bound_margins),
        "lower_bound_satisfied_rows": lower_bound_satisfied,
        "lower_bound_satisfaction_fraction": (
            None if not lower_bound_margins
            else lower_bound_satisfied / len(lower_bound_margins)
        ),
        "mean_lower_bound_deficit_vol_pct_point": (
            None if not deficits else sum(deficits) / len(deficits)
        ),
        "max_lower_bound_deficit_vol_pct_point": (
            None if not deficits else max(deficits)
        ),
        "lower_bound_constraint_status": (
            "not_applicable" if not lower_bound_margins
            else "satisfied" if max(deficits) <= 1.0e-12
            else "violated"
        ),
    }


def _operator_rows_for_case(
    case: FfiSourceStateEnvelopeCase,
    *,
    group: Literal["radius", "height"],
) -> list[dict[str, object]]:
    """Aggregate one corner through the finite FFI observation operator.

    Radius groups reproduce the reported arc maximum. Height groups retain the
    same exact-coordinate predictions while maximising only over bearings at a
    declared sensor height. Withheld predictions are counted explicitly and
    never replaced by zero; a group with no available prediction therefore has
    ``operator_status='withheld'``.
    """
    if group not in {"radius", "height"}:
        raise ValueError("operator group must be 'radius' or 'height'")
    buckets: dict[float, dict[str, object]] = {}
    for row in case.rows:
        key = round(float(row.radius_m if group == "radius" else row.height_m), 6)
        bucket = buckets.setdefault(
            key,
            {
                "observed": [],
                "available": [],
                "withheld": 0,
                "lower_bound": [],
                "margins": [],
            },
        )
        bucket["observed"].append(float(row.observed_peak_vol_pct))
        bucket["lower_bound"].append(bool(row.observed_is_lower_bound))
        if row.predicted_vol_pct is None:
            bucket["withheld"] += 1
            continue
        predicted = float(row.predicted_vol_pct)
        bucket["available"].append(predicted)
        if row.observed_is_lower_bound:
            bucket["margins"].append(predicted - float(row.observed_peak_vol_pct))

    output: list[dict[str, object]] = []
    for key in sorted(buckets):
        bucket = buckets[key]
        observed = max(bucket["observed"])
        available = bucket["available"]
        projected = max(available) if available else None
        margins = bucket["margins"]
        satisfied = sum(margin >= -1.0e-12 for margin in margins)
        withheld = int(bucket["withheld"])
        lower_bound_constraint_status = (
            "not_applicable" if not margins
            else "violated" if any(margin < -1.0e-12 for margin in margins)
            else "satisfied"
        )
        if not available:
            operator_status = "withheld"
        elif withheld:
            operator_status = "partial"
        else:
            operator_status = lower_bound_constraint_status
        output.append({
            ("radius_m" if group == "radius" else "height_m"): key,
            "sensor_count": len(bucket["observed"]),
            "available_sensor_count": len(available),
            "withheld_sensor_count": withheld,
            "observed_arc_max_vol_pct": observed,
            "projected_sensor_max_vol_pct": projected,
            "observed_over_projected": (
                None if projected is None or projected <= 0.0
                else observed / projected
            ),
            "observed_is_lower_bound": all(bucket["lower_bound"]),
            "lower_bound_satisfaction_fraction": (
                None if not margins else satisfied / len(margins)
            ),
            "worst_lower_bound_deficit_vol_pct_point": (
                None if not margins else max(-min(margins), 0.0)
            ),
            "lower_bound_constraint_status": lower_bound_constraint_status,
            "operator_status": operator_status,
        })
    return output


@dataclass(frozen=True)
class FfiSourceStateEnvelopeCase:
    """One exact source-state corner and its sensor residual map."""

    selection: tuple[tuple[str, float], ...]
    state: FfiSourceState
    rows: tuple[FfiResidualRow, ...]
    status: Literal["complete", "partial", "blocked"]
    reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _validate_selection(self.selection)
        if not isinstance(self.state, FfiSourceState):
            raise TypeError("FFI envelope case state must be FfiSourceState")
        expected_fields = self.state.uncertainty_fields()
        selected = dict(self.selection)
        if set(selected) != set(expected_fields):
            raise ValueError("FFI envelope selection must cover every source-state field")
        for name, value in selected.items():
            if not math.isclose(
                float(value), float(expected_fields[name].nominal),
                rel_tol=1.0e-10, abs_tol=1.0e-12,
            ):
                raise ValueError("FFI envelope selection does not match its exact state")
        if not self.rows:
            raise ValueError("FFI envelope case requires at least one residual row")
        if any(not isinstance(row, FfiResidualRow) for row in self.rows):
            raise TypeError("FFI envelope case rows must be FfiResidualRow values")
        sensors = [row.sensor for row in self.rows]
        if len(sensors) != len(set(sensors)):
            raise ValueError("FFI envelope case sensor identifiers must be unique")
        if self.status not in {"complete", "partial", "blocked"}:
            raise ValueError("invalid FFI envelope case status")
        if self.status == "blocked" and not self.reasons:
            raise ValueError("blocked FFI envelope case requires a reason")
        withheld = any(row.predicted_vol_pct is None for row in self.rows)
        if self.status == "complete" and withheld:
            raise ValueError("complete FFI envelope case cannot contain withheld rows")
        if self.status == "partial" and not withheld:
            raise ValueError("partial FFI envelope case requires a withheld row")

    @property
    def metrics(self) -> dict[str, object]:
        return _case_metrics(self.rows)

    def as_record(self) -> dict[str, object]:
        return {
            "selection": dict(self.selection),
            "status": self.status,
            "reasons": list(self.reasons),
            "metrics": self.metrics,
            "rows": [row.as_record() for row in self.rows],
        }


@dataclass(frozen=True)
class FfiSourceStateEnvelope:
    """Complete deterministic source-state envelope for one outdoor trial."""

    trial_test: int
    source_state: FfiSourceState
    cases: tuple[FfiSourceStateEnvelopeCase, ...]
    corrections: bool
    excluded_over_range_sensor_count: int
    status: Literal["complete", "partial", "blocked"]
    reasons: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if isinstance(self.trial_test, bool) or not isinstance(self.trial_test, int) or self.trial_test <= 0:
            raise ValueError("trial_test must be a positive integer")
        if not isinstance(self.source_state, FfiSourceState):
            raise TypeError("source_state must be FfiSourceState")
        if not self.cases:
            raise ValueError("FFI source-state envelope requires at least one case")
        if any(not isinstance(case, FfiSourceStateEnvelopeCase) for case in self.cases):
            raise TypeError("FFI envelope cases must be FfiSourceStateEnvelopeCase values")
        if not isinstance(self.corrections, bool):
            raise TypeError("corrections must be boolean")
        if (
            isinstance(self.excluded_over_range_sensor_count, bool)
            or not isinstance(self.excluded_over_range_sensor_count, int)
            or self.excluded_over_range_sensor_count < 0
        ):
            raise ValueError("excluded_over_range_sensor_count must be non-negative")
        if self.status not in {"complete", "partial", "blocked"}:
            raise ValueError("invalid FFI source-state envelope status")
        if self.status == "blocked" and not self.reasons:
            raise ValueError("blocked FFI source-state envelope requires a reason")
        expected_status = (
            "complete" if all(case.status == "complete" for case in self.cases)
            else "partial" if any(case.status != "blocked" for case in self.cases)
            else "blocked"
        )
        if self.status != expected_status:
            raise ValueError("FFI source-state envelope status is inconsistent with its cases")

    def sensor_envelope(self) -> dict[str, dict[str, object]]:
        values: dict[str, list[float]] = {}
        residuals: dict[str, list[tuple[float, tuple[tuple[str, float], ...]]]] = {}
        observed: dict[str, float] = {}
        lower_bound: dict[str, bool] = {}
        for case in self.cases:
            for row in case.rows:
                if row.sensor in observed and not math.isclose(
                    observed[row.sensor], row.observed_peak_vol_pct,
                    rel_tol=1.0e-12, abs_tol=1.0e-12,
                ):
                    raise ValueError(
                        f"FFI sensor {row.sensor!r} has inconsistent observed values"
                    )
                if row.sensor in lower_bound and lower_bound[row.sensor] != row.observed_is_lower_bound:
                    raise ValueError(
                        f"FFI sensor {row.sensor!r} has inconsistent lower-bound flags"
                    )
                observed[row.sensor] = row.observed_peak_vol_pct
                lower_bound[row.sensor] = row.observed_is_lower_bound
                if row.predicted_vol_pct is not None:
                    values.setdefault(row.sensor, []).append(float(row.predicted_vol_pct))
                    residuals.setdefault(row.sensor, []).append(
                        (float(row.residual_vol_pct), case.selection)
                    )
        output = {}
        for sensor in sorted(observed):
            sensor_residuals = residuals.get(sensor, [])
            lower_bound_margins = [
                (
                    residual,
                    selection,
                )
                for residual, selection in sensor_residuals
                if lower_bound[sensor]
            ]
            worst_lower_bound = min(
                lower_bound_margins,
                key=lambda item: (item[0], repr(item[1])),
                default=None,
            )
            satisfied_count = sum(
                margin >= -1.0e-12 for margin, _selection in lower_bound_margins
            )
            worst = max(
                sensor_residuals,
                key=lambda item: (abs(item[0]), item[0], repr(item[1])),
                default=None,
            )
            output[sensor] = {
                "observed_peak_vol_pct": observed[sensor],
                "observed_is_lower_bound": lower_bound[sensor],
                "predicted_lower_vol_pct": min(values[sensor]) if values.get(sensor) else None,
                "predicted_upper_vol_pct": max(values[sensor]) if values.get(sensor) else None,
                "residual_lower_vol_pct_point": (
                    min(item[0] for item in sensor_residuals)
                    if sensor_residuals else None
                ),
                "residual_upper_vol_pct_point": (
                    max(item[0] for item in sensor_residuals)
                    if sensor_residuals else None
                ),
                "worst_abs_residual_vol_pct_point": (
                    abs(worst[0]) if worst is not None else None
                ),
                "worst_residual_selection": (
                    dict(worst[1]) if worst is not None else None
                ),
                "available_case_count": len(values.get(sensor, [])),
                "lower_bound_satisfied_case_count": (
                    satisfied_count if lower_bound_margins else 0
                ),
                "lower_bound_satisfaction_fraction": (
                    None if not lower_bound_margins
                    else satisfied_count / len(lower_bound_margins)
                ),
                "worst_lower_bound_deficit_vol_pct_point": (
                    None if worst_lower_bound is None
                    else max(-worst_lower_bound[0], 0.0)
                ),
                "worst_lower_bound_selection": (
                    None if worst_lower_bound is None
                    else dict(worst_lower_bound[1])
                ),
                "lower_bound_constraint_status": (
                    "not_applicable" if not lower_bound_margins
                    else "satisfied"
                    if max(-item[0] for item in lower_bound_margins) <= 1.0e-12
                    else "violated"
                ),
            }
        return output

    def observation_operator_tables(self) -> dict[str, list[dict[str, object]]]:
        """Return corner-wise arc and sensor-height aggregate tables.

        The tables are deliberately additive to the sensor residual map. Each
        row retains the exact source-state selection and corner status, while
        withheld sensors remain visible in the counts and status fields.
        """
        tables = {
            "arc_max_table": [],
            "sensor_height_table": [],
        }
        for case_index, case in enumerate(self.cases):
            for table_name, group in (
                ("arc_max_table", "radius"),
                ("sensor_height_table", "height"),
            ):
                for row in _operator_rows_for_case(case, group=group):
                    row["case_index"] = case_index
                    row["case_status"] = case.status
                    row["selection"] = dict(case.selection)
                    tables[table_name].append(row)
        return tables

    def as_record(self) -> dict[str, object]:
        return {
            "schema": FFI_SOURCE_STATE_ENVELOPE_SCHEMA,
            "scope": (
                "deterministic source/weather/ambient sensitivity at reported FFI "
                "sensors; not a confidence interval or operational validation"
            ),
            "trial_test": self.trial_test,
            "status": self.status,
            "disposition": {
                "status": "research_only",
                "operational_screening_allowed": False,
                "validation_qualified": False,
                "reason": (
                    "source-state corners are deterministic sensitivity results; "
                    "reported FFI sensor peaks remain lower-bound observations"
                ),
            },
            "reasons": list(self.reasons),
            "warnings": list(self.warnings),
            "corrections": self.corrections,
            "excluded_over_range_sensor_count": self.excluded_over_range_sensor_count,
            "source_state": self.source_state.as_record(),
            "sensor_envelope": self.sensor_envelope(),
            "observation_operator": self.observation_operator_tables(),
            "cases": [case.as_record() for case in self.cases],
            "limitations": [
                "FFI arc/sensor peaks are lower bounds because the plume may pass between sensors",
                "corners have no probability weights",
                "no source or model coefficient was fitted",
                "the free-field steady jet has no plume storage or transient source history",
            ],
        }


def _coordinates(reading: Reading, wind_from_deg: float) -> tuple[float, float]:
    bearing = math.radians(float(reading.bearing))
    east = float(reading.radius) * math.sin(bearing)
    north = float(reading.radius) * math.cos(bearing)
    direction = math.radians(float(wind_from_deg))
    down_east, down_north = -math.sin(direction), -math.cos(direction)
    right_east, right_north = down_north, -down_east
    return (
        east * down_east + north * down_north,
        east * right_east + north * right_north,
    )


def _rows_for_case(
    trial: Trial,
    state: FfiSourceState,
    readings: tuple[Reading, ...],
    *,
    corrections: bool,
    max_distance: float,
) -> tuple[tuple[FfiResidualRow, ...], tuple[str, ...]]:
    coordinates = {
        reading.sensor: _coordinates(reading, state.wind_direction_from_deg.nominal)
        for reading in readings
    }
    try:
        setup, initial = hydrogen_jet(
            rate=state.rate_kg_s.nominal,
            diameter=state.orifice_m.nominal,
            wind=state.wind_m_s.nominal,
            height=state.release_height_m.nominal,
            ambient_temperature=state.ambient_temperature_k.nominal,
            relative_humidity=state.relative_humidity_pct.nominal,
            ambient_pressure=state.ambient_pressure_pa.nominal,
            storage_pressure_barg=state.storage_pressure_barg.nominal,
            wind_reference_height=state.wind_reference_height_m.nominal,
            corrections=corrections,
            ground_effect=False,
        )
        result = setup.run(
            initial, distmx=STEP, smax=max(REACH, 2.5 * max_distance)
        )
        trajectory = Trajectory(setup.th.table, result.rows)
        if not trajectory.ok:
            raise RuntimeError("free-field FFI trajectory did not complete")
    except (RuntimeError, ValueError, OverflowError, FloatingPointError) as error:
        reason = f"transport failed for this source-state corner: {error}"
        return tuple(
            FfiResidualRow(
                sensor=reading.sensor,
                radius_m=reading.radius,
                bearing_deg=reading.bearing,
                height_m=reading.height,
                alongwind_m=coordinates[reading.sensor][0],
                crosswind_m=coordinates[reading.sensor][1],
                observed_peak_vol_pct=reading.peak,
                observed_is_lower_bound=True,
                predicted_vol_pct=None,
                residual_vol_pct=None,
                withheld_reason=reason,
            )
            for reading in readings
        ), (reason,)

    rows = []
    reasons = []
    for reading in readings:
        alongwind, crosswind = coordinates[reading.sensor]
        reason = None
        predicted = None
        if alongwind <= 0.0:
            reason = "sensor is upwind of the source for this wind-direction corner"
        elif trajectory.at(alongwind) is None:
            reason = "trajectory does not reach the sensor alongwind coordinate"
        else:
            predicted = max(
                0.0,
                float(trajectory.concentration_at(alongwind, crosswind, reading.height)),
            )
            if not math.isfinite(predicted):
                reason = "trajectory returned a non-finite sensor prediction"
                predicted = None
        if reason is not None:
            reasons.append(f"{reading.sensor}: {reason}")
        rows.append(FfiResidualRow(
            sensor=reading.sensor,
            radius_m=reading.radius,
            bearing_deg=reading.bearing,
            height_m=reading.height,
            alongwind_m=alongwind,
            crosswind_m=crosswind,
            observed_peak_vol_pct=reading.peak,
            observed_is_lower_bound=True,
            predicted_vol_pct=predicted,
            residual_vol_pct=(None if predicted is None else predicted - reading.peak),
            withheld_reason=reason,
        ))
    return tuple(rows), tuple(reasons)


def run_ffi_source_state_envelope(
    trial: Trial,
    source_state: FfiSourceState,
    *,
    corrections: bool = True,
    readings: tuple[Reading, ...] | list[Reading] | None = None,
    max_cases: int = 64,
) -> FfiSourceStateEnvelope:
    """Run every declared source-state corner against one FFI outdoor trial.

    ``readings`` defaults to all non-over-range readings.  A custom subset is
    useful for a declared arc or height audit, but every retained prediction is
    still evaluated at the sensor coordinates rather than at the plume
    centreline.  A sensor that is upwind or outside the computed trajectory is
    kept as a withheld row.
    """
    if not isinstance(trial, Trial):
        raise TypeError("trial must be a spadeadam Trial")
    if not trial.outdoor or not trial.horizontal:
        raise ValueError(
            f"FFI source-state envelope requires an outdoor horizontal trial, got test {trial.test}"
        )
    if not isinstance(source_state, FfiSourceState):
        raise TypeError("source_state must be FfiSourceState")
    if not isinstance(corrections, bool):
        raise TypeError("corrections must be boolean")
    selected = tuple(trial.readings if readings is None else readings)
    if not selected:
        raise ValueError("FFI source-state envelope requires at least one sensor reading")
    if any(not isinstance(item, Reading) for item in selected):
        raise TypeError("readings must contain Reading values")
    sensors = [item.sensor for item in selected]
    if len(sensors) != len(set(sensors)):
        raise ValueError("FFI source-state envelope sensor identifiers must be unique")
    over_range = sum(1 for item in selected if item.over_range)
    usable = tuple(item for item in selected if not item.over_range)
    if not usable:
        raise ValueError("FFI source-state envelope has no non-over-range sensors")
    max_distance = max(float(item.radius) for item in usable)
    cases = []
    for selection, corner in source_state.corner_cases(max_cases=max_cases):
        rows, reasons = _rows_for_case(
            trial, corner, usable, corrections=corrections, max_distance=max_distance,
        )
        available = sum(row.predicted_vol_pct is not None for row in rows)
        if available == 0:
            status: Literal["complete", "partial", "blocked"] = "blocked"
        elif reasons:
            status = "partial"
        else:
            status = "complete"
        cases.append(FfiSourceStateEnvelopeCase(
            selection=selection, state=corner, rows=rows, status=status,
            reasons=reasons,
        ))
    complete_cases = tuple(cases)
    status: Literal["complete", "partial", "blocked"]
    if all(case.status == "complete" for case in complete_cases):
        status = "complete"
    elif any(case.status != "blocked" for case in complete_cases):
        status = "partial"
    else:
        status = "blocked"
    reasons = tuple(
        f"corner {dict(case.selection)}: {reason}"
        for case in complete_cases if case.status == "blocked"
        for reason in case.reasons
    )
    warnings = (
        "reported FFI peaks are lower bounds, not exact centreline observations",
        "source, wind and ambient corners are deterministic sensitivities without probability weights",
    )
    return FfiSourceStateEnvelope(
        trial_test=trial.test,
        source_state=source_state,
        cases=complete_cases,
        corrections=corrections,
        excluded_over_range_sensor_count=over_range,
        status=status,
        reasons=reasons,
        warnings=warnings,
    )


def ffi_source_state_envelope_report(result: FfiSourceStateEnvelope) -> dict[str, object]:
    """Return a JSON-safe audit record for a source-state envelope."""
    if not isinstance(result, FfiSourceStateEnvelope):
        raise TypeError("result must be FfiSourceStateEnvelope")
    return result.as_record()


__all__ = [
    "FFI_SOURCE_STATE_ENVELOPE_SCHEMA",
    "ffi_reference_provenance",
    "FfiSourceState",
    "FfiResidualRow",
    "FfiSourceStateEnvelopeCase",
    "FfiSourceStateEnvelope",
    "run_ffi_source_state_envelope",
    "ffi_source_state_envelope_report",
]
