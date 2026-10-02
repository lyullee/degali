"""Observation intake for the public FFI open-channel hydrogen experiments.

The record at DOI 10.23642/usn.26117989.v2 provides a useful, independently
timed hydrogen source and 29 concentration signals.  It is *not* an outdoor
LH2-pool experiment.  This module consequently reads a local public CSV and
reduces declared threshold windows, while making it impossible to describe
that reduction as a calibrated atmospheric-LH2 plume validation.

No measurement file, derived series, geometry factor, or fitted parameter is
distributed with DEGALI.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import math
from pathlib import Path

import numpy as np


_FLOW_TIME = "flow time [s]"
_FLOW_RATE = "mass flow meter 1 [g/s]"
_SENSOR_TIME = "h2 sensor time"
_SENSOR_SUFFIX = " h2 concentration [%]"


@dataclass(frozen=True)
class ThresholdWindow:
    """A threshold-defined concentration timing observation.

    The threshold is caller-declared in published sensor units.  No baseline
    estimate, smoothing window, or connected-component selection is hidden in
    this object: arrival and departure are the first and final recorded values
    at or above that threshold.
    """

    sensor: str
    threshold_percent: float
    arrival_s: float
    peak_time_s: float
    departure_s: float
    peak_percent: float

    @property
    def duration_s(self) -> float:
        return self.departure_s - self.arrival_s


@dataclass(frozen=True)
class OpenChannelHydrogenRun:
    """One local FFI open-ended-channel CSV in its published units."""

    source_path: Path
    flow_time_s: np.ndarray
    mass_flow_g_s: np.ndarray
    sensor_time_s: np.ndarray
    sensor_percent: dict[str, np.ndarray]
    dataset_doi: str = "10.23642/usn.26117989.v2"
    source_geometry: str = "open_ended_rectangular_channel"

    @property
    def quantitative_lh2_pool_validation_allowed(self) -> bool:
        """The experiment is hydrogen, but not an LH2-pool plume analogue."""

        return False

    @property
    def validation_scope(self) -> str:
        return "hydrogen_source_and_sensor_timing_in_channel_only"

    @property
    def sensor_names(self) -> tuple[str, ...]:
        return tuple(sorted(self.sensor_percent))


def _float(value: str | None) -> float:
    if value is None or not value.strip():
        return math.nan
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(f"invalid numerical value {value!r} in FFI H2 CSV") from error


def _finite_pair(time: np.ndarray, value: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    keep = np.isfinite(time) & np.isfinite(value)
    return time[keep], value[keep]


def read_open_channel_h2_csv(path: str | Path) -> OpenChannelHydrogenRun:
    """Read one user-supplied CSV from the FFI channel data record.

    The reader intentionally requires the separately logged flow and sensor
    clocks rather than assuming rows imply a common time base.  Extra metadata
    columns and the dataset's descriptive early rows are preserved as ignored
    CSV rows; only explicitly named physical columns are used.
    """

    source = Path(path)
    with source.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = tuple(reader.fieldnames or ())
        required = {_FLOW_TIME, _FLOW_RATE, _SENSOR_TIME}
        missing = required - set(fields)
        if missing:
            raise ValueError(
                "FFI H2 CSV is missing required columns: " + ", ".join(sorted(missing))
            )
        sensor_columns = tuple(field for field in fields if field.endswith(_SENSOR_SUFFIX))
        if not sensor_columns:
            raise ValueError("FFI H2 CSV contains no hydrogen concentration columns")
        rows = list(reader)

    flow_time = np.asarray([_float(row.get(_FLOW_TIME)) for row in rows], dtype=float)
    flow_rate = np.asarray([_float(row.get(_FLOW_RATE)) for row in rows], dtype=float)
    sensor_time = np.asarray([_float(row.get(_SENSOR_TIME)) for row in rows], dtype=float)
    flow_time, flow_rate = _finite_pair(flow_time, flow_rate)
    sensor_time_finite = np.isfinite(sensor_time)
    sensor_time = sensor_time[sensor_time_finite]
    if flow_time.size < 2 or sensor_time.size < 2:
        raise ValueError("FFI H2 CSV needs at least two finite source and sensor timestamps")
    if np.any(np.diff(flow_time) < 0.0) or np.any(np.diff(sensor_time) < 0.0):
        raise ValueError("FFI H2 source and sensor clocks must be nondecreasing")

    sensors: dict[str, np.ndarray] = {}
    for column in sensor_columns:
        values = np.asarray([_float(row.get(column)) for row in rows], dtype=float)
        values = values[sensor_time_finite]
        if not np.any(np.isfinite(values)):
            continue
        sensors[column.removesuffix(_SENSOR_SUFFIX)] = values
    if not sensors:
        raise ValueError("FFI H2 CSV has no finite sensor measurements")
    return OpenChannelHydrogenRun(
        source_path=source, flow_time_s=flow_time, mass_flow_g_s=flow_rate,
        sensor_time_s=sensor_time, sensor_percent=sensors,
    )


def threshold_window(
    run: OpenChannelHydrogenRun, sensor: str, *, threshold_percent: float,
) -> ThresholdWindow | None:
    """Return one declared threshold window, or ``None`` when it is never met."""

    if not math.isfinite(threshold_percent) or threshold_percent <= 0.0:
        raise ValueError("threshold_percent must be a finite positive value")
    try:
        concentrations = run.sensor_percent[sensor]
    except KeyError as error:
        raise ValueError(f"unknown FFI H2 sensor {sensor!r}") from error
    if concentrations.shape != run.sensor_time_s.shape:
        raise ValueError("FFI H2 sensor and clock arrays have inconsistent lengths")
    above = np.isfinite(concentrations) & (concentrations >= threshold_percent)
    if not np.any(above):
        return None
    indices = np.flatnonzero(above)
    peak_index = int(np.nanargmax(np.where(above, concentrations, np.nan)))
    return ThresholdWindow(
        sensor=sensor, threshold_percent=float(threshold_percent),
        arrival_s=float(run.sensor_time_s[indices[0]]),
        peak_time_s=float(run.sensor_time_s[peak_index]),
        departure_s=float(run.sensor_time_s[indices[-1]]),
        peak_percent=float(concentrations[peak_index]),
    )


__all__ = [
    "OpenChannelHydrogenRun", "ThresholdWindow", "read_open_channel_h2_csv",
    "threshold_window",
]
