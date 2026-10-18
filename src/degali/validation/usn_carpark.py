"""Read and score the public USN mechanically ventilated H2 experiments.

The data record at DOI ``10.23642/USN.14405903`` contains MATLAB v7.3
files.  Those files are HDF5 containers with a common local experiment clock
for the mass-flow meter, tank, ventilation boundary and concentration
sensors.  This module preserves the published values and exposes the dynamic
4 vol.% observation operator used in the accompanying publication.

The experiments use gaseous hydrogen in a confined ISO container.  They are
therefore component-level evidence for source/sensor timing and ventilation,
not end-to-end liquid-hydrogen field validation for DEGALI.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re
from typing import Literal

import numpy as np


DATASET_DOI = "10.23642/USN.14405903"
CONTAINER_VOLUME_M3 = 60.8
AMBIENT_HYDROGEN_DENSITY_KG_M3 = 0.08375
_CANBUS_ID = re.compile(r"^D(?P<number>\d{3})$")
_WIFI_ID = re.compile(r"^DC\d+$")


def _h5py():
    try:
        import h5py
    except ImportError as error:  # pragma: no cover - environment boundary
        raise RuntimeError(
            "USN car-park MATLAB v7.3 intake requires h5py; install DEGALI's "
            "test/research dependencies"
        ) from error
    return h5py


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matlab_value(handle, dataset):
    """Resolve one MATLAB v7.3 dataset, including char and cell values."""

    values = np.asarray(dataset)
    if dataset.attrs.get("MATLAB_class", b"") == b"char":
        return "".join(
            chr(int(value)) for value in values.ravel(order="F") if value
        )
    if values.dtype == object:
        return tuple(
            _matlab_value(handle, handle[reference])
            for reference in values.ravel(order="F")
        )
    return values


def _referenced(handle, group: str, field: str, index: int):
    reference = handle[f"{group}/{field}"][index, 0]
    return _matlab_value(handle, handle[reference])


def _numeric_vector(value, *, name: str) -> np.ndarray:
    if isinstance(value, str):
        raise ValueError(f"USN {name} is textual rather than numerical")
    result = np.asarray(value, dtype=float).ravel()
    if result.size == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"USN {name} must contain finite numerical values")
    return result


@dataclass(frozen=True)
class UsnCarparkSensor:
    """One published concentration sensor in native clock/coordinate units."""

    sensor_id: str
    position_raw_mm: tuple[float, float, float] | None
    time_s: np.ndarray
    concentration_vol_percent: np.ndarray
    family: Literal["canbus", "wifi", "unknown"]

    def __post_init__(self) -> None:
        if not self.sensor_id.strip():
            raise ValueError("USN sensor_id must be non-empty")
        if self.time_s.ndim != 1 or self.concentration_vol_percent.shape != self.time_s.shape:
            raise ValueError("USN sensor time and concentration shapes must match")
        if self.time_s.size < 2 or not np.all(np.isfinite(self.time_s)):
            raise ValueError("USN sensor needs at least two finite timestamps")
        if np.any(np.diff(self.time_s) < 0.0):
            raise ValueError("USN sensor timestamps must be nondecreasing")
        if self.position_raw_mm is not None and (
            len(self.position_raw_mm) != 3
            or not all(math.isfinite(value) for value in self.position_raw_mm)
        ):
            raise ValueError("USN sensor position must contain three finite values")

    @property
    def is_published_room_cloud_sensor(self) -> bool:
        """Whether this sensor belongs to the D139-D162 room-cloud array.

        D163-D168 are the under-table group.  The paper's Table 3 threshold
        operator uses the room array, which reproduces its first-arrival and
        duration values without including the high local under-table peaks.
        """

        match = _CANBUS_ID.fullmatch(self.sensor_id)
        return match is not None and 139 <= int(match.group("number")) <= 162


@dataclass(frozen=True)
class UsnCarparkRun:
    """One local USN experiment with its published common local clock."""

    source_path: Path
    experiment_id: str
    source_time_s: np.ndarray
    mass_flow_g_s: np.ndarray
    tank_pressure_bar: np.ndarray
    tank_temperature_c: np.ndarray
    ventilation_flow_m3_h: np.ndarray
    source_t0_s: float
    sensors: tuple[UsnCarparkSensor, ...]
    raw_sha256: str
    dataset_doi: str = DATASET_DOI

    def __post_init__(self) -> None:
        arrays = (
            self.source_time_s,
            self.mass_flow_g_s,
            self.tank_pressure_bar,
            self.tank_temperature_c,
            self.ventilation_flow_m3_h,
        )
        if any(array.ndim != 1 for array in arrays):
            raise ValueError("USN source arrays must be one-dimensional")
        if any(array.shape != self.source_time_s.shape for array in arrays[1:]):
            raise ValueError("USN source arrays must share one common clock")
        if self.source_time_s.size < 2 or np.any(np.diff(self.source_time_s) <= 0.0):
            raise ValueError("USN source clock must be strictly increasing")
        if not all(np.all(np.isfinite(array)) for array in arrays):
            raise ValueError("USN source arrays must be finite")
        if not math.isfinite(self.source_t0_s):
            raise ValueError("USN source t0 must be finite")
        if len(self.raw_sha256) != 64:
            raise ValueError("USN raw SHA-256 is malformed")
        if not self.sensors:
            raise ValueError("USN run contains no concentration sensors")

    @property
    def room_cloud_sensors(self) -> tuple[UsnCarparkSensor, ...]:
        return tuple(
            sensor for sensor in self.sensors
            if sensor.is_published_room_cloud_sensor
        )

    @property
    def canbus_sensors(self) -> tuple[UsnCarparkSensor, ...]:
        return tuple(sensor for sensor in self.sensors if sensor.family == "canbus")

    @property
    def wifi_sensors(self) -> tuple[UsnCarparkSensor, ...]:
        return tuple(sensor for sensor in self.sensors if sensor.family == "wifi")

    @property
    def validation_scope(self) -> str:
        return "gh2_confined_dynamic_transport_component_only"

    @property
    def quantitative_lh2_field_validation_allowed(self) -> bool:
        return False


@dataclass(frozen=True)
class UsnThresholdMetrics:
    experiment_id: str
    threshold_vol_percent: float
    sensor_scope: str
    sensor_count: int
    observed_peak_vol_percent: float
    first_arrival_after_source_s: float | None
    final_exceedance_after_source_s: float | None
    exceedance_duration_s: float


@dataclass(frozen=True)
class UsnWellMixedComparison:
    """No-fit 0-D ventilation baseline compared with the room sensor array."""

    experiment_id: str
    threshold_vol_percent: float
    observed_peak_vol_percent: float
    model_peak_vol_percent: float
    peak_ratio_model_to_observed: float
    observed_arrival_after_source_s: float | None
    model_arrival_after_source_s: float | None
    observed_exceedance_duration_s: float
    model_exceedance_duration_s: float
    mae_against_spatial_mean_vol_percent: float
    rmse_against_spatial_mean_vol_percent: float


def read_usn_carpark_mat(path: str | Path) -> UsnCarparkRun:
    """Read one unmodified public MATLAB v7.3 experiment file."""

    h5py = _h5py()
    source = Path(path)
    match = re.fullmatch(r"Exp(?P<number>\d+)\.mat", source.name, re.IGNORECASE)
    if match is None:
        raise ValueError("USN experiment filename must match ExpNN.mat")
    with h5py.File(source, "r") as handle:
        required = {"MFM", "S", "Tank", "Vent"}
        missing = required - set(handle)
        if missing:
            raise ValueError("USN MATLAB file is missing groups: " + ", ".join(sorted(missing)))
        source_time = np.asarray(handle["MFM/time"], dtype=float).ravel()
        ventilation_time = np.asarray(handle["Vent/time"], dtype=float).ravel()
        tank_time = np.asarray(handle["Tank/time"], dtype=float).ravel()
        if not np.array_equal(source_time, ventilation_time) or not np.array_equal(
            source_time, tank_time
        ):
            raise ValueError("USN MFM, ventilation and tank clocks differ")

        sensors: list[UsnCarparkSensor] = []
        for index in range(handle["S/ID"].shape[0]):
            sensor_id = str(_referenced(handle, "S", "ID", index)).strip()
            raw_position = _referenced(handle, "S", "pos", index)
            position = None
            if not isinstance(raw_position, str):
                position_values = np.asarray(raw_position, dtype=float).ravel()
                if position_values.size == 3 and np.all(np.isfinite(position_values)):
                    position = tuple(float(value) for value in position_values)
            time_s = _numeric_vector(
                _referenced(handle, "S", "time", index),
                name=f"{sensor_id} time",
            )
            concentration = _numeric_vector(
                _referenced(handle, "S", "conc", index),
                name=f"{sensor_id} concentration",
            )
            family: Literal["canbus", "wifi", "unknown"] = "unknown"
            if _CANBUS_ID.fullmatch(sensor_id):
                family = "canbus"
            elif _WIFI_ID.fullmatch(sensor_id):
                family = "wifi"
            sensors.append(UsnCarparkSensor(
                sensor_id=sensor_id,
                position_raw_mm=position,
                time_s=time_s,
                concentration_vol_percent=concentration,
                family=family,
            ))

        return UsnCarparkRun(
            source_path=source.resolve(),
            experiment_id=f"USN-Exp{int(match.group('number'))}",
            source_time_s=source_time,
            mass_flow_g_s=np.asarray(handle["MFM/mfr"], dtype=float).ravel(),
            tank_pressure_bar=np.asarray(handle["Tank/p"], dtype=float).ravel(),
            tank_temperature_c=np.asarray(handle["Tank/T"], dtype=float).ravel(),
            ventilation_flow_m3_h=np.asarray(handle["Vent/vfr"], dtype=float).ravel(),
            source_t0_s=float(np.asarray(handle["MFM/t0"], dtype=float).ravel()[0]),
            sensors=tuple(sensors),
            raw_sha256=_sha256(source),
        )


def _scope_sensors(
    run: UsnCarparkRun, scope: Literal["room_cloud", "canbus", "all"],
) -> tuple[UsnCarparkSensor, ...]:
    if scope == "room_cloud":
        sensors = run.room_cloud_sensors
    elif scope == "canbus":
        sensors = run.canbus_sensors
    elif scope == "all":
        sensors = tuple(
            sensor for sensor in run.sensors if sensor.position_raw_mm is not None
        )
    else:  # pragma: no cover - Literal guard for direct Python callers
        raise ValueError(f"unsupported USN sensor scope {scope!r}")
    if not sensors:
        raise ValueError(f"USN {scope} sensor scope is empty")
    reference = sensors[0].time_s
    if any(
        sensor.time_s.shape != reference.shape
        or not np.array_equal(sensor.time_s, reference)
        for sensor in sensors[1:]
    ):
        raise ValueError(f"USN {scope} sensors do not share one exact clock")
    return sensors


def threshold_metrics(
    run: UsnCarparkRun,
    *,
    threshold_vol_percent: float = 4.0,
    sensor_scope: Literal["room_cloud", "canbus", "all"] = "room_cloud",
) -> UsnThresholdMetrics:
    """Apply a declared any-sensor exceedance operator on one common clock."""

    if not math.isfinite(threshold_vol_percent) or threshold_vol_percent <= 0.0:
        raise ValueError("USN threshold must be finite and positive")
    sensors = _scope_sensors(run, sensor_scope)
    clock = sensors[0].time_s
    envelope = np.nanmax(
        np.vstack([sensor.concentration_vol_percent for sensor in sensors]), axis=0,
    )
    above = np.isfinite(envelope) & (envelope >= threshold_vol_percent)
    indices = np.flatnonzero(above)
    duration = float(np.sum(np.diff(clock) * above[:-1]))
    arrival = None
    departure = None
    if indices.size:
        arrival = float(clock[indices[0]] - run.source_t0_s)
        departure = float(clock[indices[-1]] - run.source_t0_s)
    return UsnThresholdMetrics(
        experiment_id=run.experiment_id,
        threshold_vol_percent=float(threshold_vol_percent),
        sensor_scope=sensor_scope,
        sensor_count=len(sensors),
        observed_peak_vol_percent=float(np.nanmax(envelope)),
        first_arrival_after_source_s=arrival,
        final_exceedance_after_source_s=departure,
        exceedance_duration_s=duration,
    )


def integrated_positive_source_mass_g(
    run: UsnCarparkRun, *, duration_s: float = 900.0,
) -> float:
    """Integrate the non-negative recorded source rate after published t0."""

    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise ValueError("USN source duration must be finite and positive")
    keep = (
        (run.source_time_s >= run.source_t0_s)
        & (run.source_time_s <= run.source_t0_s + duration_s)
    )
    time_s = run.source_time_s[keep]
    rate = np.maximum(run.mass_flow_g_s[keep], 0.0)
    if time_s.size < 2:
        raise ValueError("USN source window contains fewer than two samples")
    return float(np.sum(0.5 * (rate[:-1] + rate[1:]) * np.diff(time_s)))


def well_mixed_concentration_history(
    run: UsnCarparkRun,
    *,
    enclosure_volume_m3: float = CONTAINER_VOLUME_M3,
    hydrogen_density_kg_m3: float = AMBIENT_HYDROGEN_DENSITY_KG_M3,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the paper's no-fit, ideal well-mixed ventilation baseline.

    The explicit Euler update is performed on the 3 Hz CAN-bus room clock.
    It is a diagnostic mass-balance comparator, not the DEGALI free-field
    transport solver and not a CFD surrogate.
    """

    if enclosure_volume_m3 <= 0.0 or not math.isfinite(enclosure_volume_m3):
        raise ValueError("USN enclosure volume must be finite and positive")
    if hydrogen_density_kg_m3 <= 0.0 or not math.isfinite(hydrogen_density_kg_m3):
        raise ValueError("USN hydrogen density must be finite and positive")
    sensors = _scope_sensors(run, "room_cloud")
    clock = sensors[0].time_s
    source_rate_kg_s = np.interp(
        clock, run.source_time_s, np.maximum(run.mass_flow_g_s, 0.0) / 1000.0,
    )
    ventilation_m3_s = np.interp(
        clock, run.source_time_s,
        np.maximum(run.ventilation_flow_m3_h, 0.0) / 3600.0,
    )
    source_rate_kg_s = source_rate_kg_s.copy()
    source_rate_kg_s[clock < run.source_t0_s] = 0.0
    concentration = np.zeros_like(clock)
    for index in range(1, clock.size):
        dt = float(clock[index] - clock[index - 1])
        derivative = (
            source_rate_kg_s[index - 1]
            / (hydrogen_density_kg_m3 * enclosure_volume_m3)
            - ventilation_m3_s[index - 1]
            / enclosure_volume_m3
            * concentration[index - 1]
        )
        concentration[index] = max(0.0, concentration[index - 1] + dt * derivative)
    return clock.copy(), 100.0 * concentration


def compare_well_mixed_baseline(
    run: UsnCarparkRun, *, threshold_vol_percent: float = 4.0,
) -> UsnWellMixedComparison:
    """Compare the no-fit 0-D baseline to the exact room-array operator."""

    observed = threshold_metrics(
        run, threshold_vol_percent=threshold_vol_percent,
        sensor_scope="room_cloud",
    )
    sensors = _scope_sensors(run, "room_cloud")
    clock, model = well_mixed_concentration_history(run)
    spatial_mean = np.nanmean(
        np.vstack([sensor.concentration_vol_percent for sensor in sensors]), axis=0,
    )
    above = model >= threshold_vol_percent
    indices = np.flatnonzero(above)
    model_arrival = None if not indices.size else float(
        clock[indices[0]] - run.source_t0_s
    )
    model_duration = float(np.sum(np.diff(clock) * above[:-1]))
    model_peak = float(np.max(model))
    residual = model - spatial_mean
    return UsnWellMixedComparison(
        experiment_id=run.experiment_id,
        threshold_vol_percent=float(threshold_vol_percent),
        observed_peak_vol_percent=observed.observed_peak_vol_percent,
        model_peak_vol_percent=model_peak,
        peak_ratio_model_to_observed=(
            model_peak / observed.observed_peak_vol_percent
        ),
        observed_arrival_after_source_s=observed.first_arrival_after_source_s,
        model_arrival_after_source_s=model_arrival,
        observed_exceedance_duration_s=observed.exceedance_duration_s,
        model_exceedance_duration_s=model_duration,
        mae_against_spatial_mean_vol_percent=float(np.mean(np.abs(residual))),
        rmse_against_spatial_mean_vol_percent=float(np.sqrt(np.mean(residual ** 2))),
    )


__all__ = [
    "AMBIENT_HYDROGEN_DENSITY_KG_M3", "CONTAINER_VOLUME_M3", "DATASET_DOI",
    "UsnCarparkRun", "UsnCarparkSensor", "UsnThresholdMetrics",
    "UsnWellMixedComparison", "compare_well_mixed_baseline",
    "integrated_positive_source_mass_g", "read_usn_carpark_mat",
    "threshold_metrics", "well_mixed_concentration_history",
]
