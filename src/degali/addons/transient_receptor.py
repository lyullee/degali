"""Research-only transient observation operator for steady plume libraries.

The dispersion models in :mod:`degali` are spatial, steady integral models.
This module does not turn them into a transient CFD model.  It provides the
smaller operation needed for field-test comparison: rotate fixed receptors
into the instantaneous wind frame, sample a caller-supplied library of steady
solutions, and apply a declared first-order instrument response.

No source multiplier, meander width, or response time is inferred from an
observation.  Every time series and every steady response must be supplied by
the caller.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np


SteadySampler = Callable[[float, float, float], float]


def _one_dimensional_finite(values: Sequence[float], name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class WindHistory:
    """Measured wind record using the meteorological ``from`` convention."""

    time_s: Sequence[float]
    speed_m_s: Sequence[float]
    direction_from_deg: Sequence[float]

    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        time = _one_dimensional_finite(self.time_s, "time_s")
        speed = _one_dimensional_finite(self.speed_m_s, "speed_m_s")
        direction = _one_dimensional_finite(
            self.direction_from_deg, "direction_from_deg"
        )
        if speed.size != time.size or direction.size != time.size:
            raise ValueError("wind-history arrays must have equal lengths")
        if time.size > 1 and np.any(np.diff(time) <= 0.0):
            raise ValueError("wind-history time must be strictly increasing")
        if np.any(speed <= 0.0):
            raise ValueError("wind speed must be positive")
        return time, speed, np.mod(direction, 360.0)


@dataclass(frozen=True)
class FixedReceptor:
    """A fixed sensor in east/north/up coordinates relative to the source."""

    name: str
    east_m: float
    north_m: float
    height_m: float
    response_t90_s: float | None = None

    def validate(self) -> None:
        coordinates = (self.east_m, self.north_m, self.height_m)
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError(f"receptor {self.name!r} coordinates must be finite")
        if self.height_m < 0.0:
            raise ValueError(f"receptor {self.name!r} height must be non-negative")
        if self.response_t90_s is not None and (
            not math.isfinite(self.response_t90_s)
            or self.response_t90_s <= 0.0
        ):
            raise ValueError(
                f"receptor {self.name!r} response_t90_s must be positive"
            )


@dataclass(frozen=True)
class ReceptorTrace:
    """True and instrument-indicated histories for one fixed receptor."""

    receptor: FixedReceptor
    time_s: np.ndarray
    wind_speed_m_s: np.ndarray
    wind_direction_from_deg: np.ndarray
    alongwind_m: np.ndarray
    crosswind_m: np.ndarray
    true_mole_fraction: np.ndarray
    indicated_mole_fraction: np.ndarray

    def statistics(self, start_s: float, end_s: float) -> "WindowStatistics":
        """Return inclusive-window statistics without resampling the trace."""
        if not math.isfinite(start_s) or not math.isfinite(end_s):
            raise ValueError("statistics window must be finite")
        if end_s < start_s:
            raise ValueError("statistics window end must not precede its start")
        selected = (self.time_s >= start_s) & (self.time_s <= end_s)
        if not np.any(selected):
            raise ValueError("statistics window contains no samples")
        true = self.true_mole_fraction[selected]
        indicated = self.indicated_mole_fraction[selected]
        return WindowStatistics(
            start_s=float(start_s),
            end_s=float(end_s),
            sample_count=int(np.count_nonzero(selected)),
            true_mean=float(np.mean(true)),
            true_maximum=float(np.max(true)),
            indicated_mean=float(np.mean(indicated)),
            indicated_maximum=float(np.max(indicated)),
        )


@dataclass(frozen=True)
class WindowStatistics:
    """Mean and maximum before and after instrument response."""

    start_s: float
    end_s: float
    sample_count: int
    true_mean: float
    true_maximum: float
    indicated_mean: float
    indicated_maximum: float


class SteadyPlumeTable:
    """Linearly interpolate a small library of steady wind-speed solutions.

    Each sampler receives ``(alongwind_m, crosswind_m, height_m)`` and returns
    mole fraction.  Values outside the supplied wind range are rejected rather
    than silently extrapolated.  The expensive plume solve is therefore done
    once per table wind, not once per wind-history sample.
    """

    def __init__(
        self,
        wind_speeds_m_s: Sequence[float],
        samplers: Sequence[SteadySampler],
    ) -> None:
        speeds = _one_dimensional_finite(wind_speeds_m_s, "wind_speeds_m_s")
        if len(samplers) != speeds.size:
            raise ValueError("one steady sampler is required for each wind speed")
        if np.any(speeds <= 0.0) or (
            speeds.size > 1 and np.any(np.diff(speeds) <= 0.0)
        ):
            raise ValueError("table wind speeds must be positive and increasing")
        if not all(callable(sampler) for sampler in samplers):
            raise TypeError("every steady plume sampler must be callable")
        self.wind_speeds_m_s = speeds.copy()
        self.samplers = tuple(samplers)

    def __call__(
        self,
        wind_speed_m_s: float,
        alongwind_m: float,
        crosswind_m: float,
        height_m: float,
    ) -> float:
        speed = float(wind_speed_m_s)
        if not math.isfinite(speed):
            raise ValueError("wind speed must be finite")
        lower_bound = float(self.wind_speeds_m_s[0])
        upper_bound = float(self.wind_speeds_m_s[-1])
        tolerance = 1.0e-12 * max(upper_bound, 1.0)
        if speed < lower_bound - tolerance or speed > upper_bound + tolerance:
            raise ValueError(
                f"wind speed {speed:g} m/s is outside the steady table "
                f"[{lower_bound:g}, {upper_bound:g}] m/s"
            )
        speed = min(max(speed, lower_bound), upper_bound)
        upper = int(np.searchsorted(self.wind_speeds_m_s, speed, side="right"))
        if upper == 0:
            return self._sample(0, alongwind_m, crosswind_m, height_m)
        if upper == len(self.samplers):
            return self._sample(-1, alongwind_m, crosswind_m, height_m)
        lower = upper - 1
        low_speed = float(self.wind_speeds_m_s[lower])
        high_speed = float(self.wind_speeds_m_s[upper])
        weight = (speed - low_speed) / (high_speed - low_speed)
        low = self._sample(lower, alongwind_m, crosswind_m, height_m)
        high = self._sample(upper, alongwind_m, crosswind_m, height_m)
        return (1.0 - weight) * low + weight * high

    def _sample(
        self, index: int, alongwind_m: float, crosswind_m: float, height_m: float
    ) -> float:
        value = float(self.samplers[index](alongwind_m, crosswind_m, height_m))
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("steady plume sampler returned an invalid mole fraction")
        return value


def meteorological_from_to_math_radians(direction_from_deg: float) -> float:
    """Convert meteorological wind-from bearing to an east/north math angle.

    The returned angle is measured counter-clockwise from east and points in
    the direction the wind blows toward.  Thus wind from 270 degrees (west)
    returns zero: transport is toward east.
    """
    if not math.isfinite(direction_from_deg):
        raise ValueError("wind direction must be finite")
    return math.radians((270.0 - direction_from_deg) % 360.0)


def first_order_sensor_response(
    time_s: Sequence[float],
    signal: Sequence[float],
    *,
    t90_s: float | None,
    initial_value: float = 0.0,
) -> np.ndarray:
    """Apply an exact discrete first-order response to a sampled signal.

    The input is held at its current sample value over each interval.  ``t90``
    is converted to the time constant by ``tau = t90 / ln(10)``.  ``None``
    bypasses the response model.
    """
    time = _one_dimensional_finite(time_s, "time_s")
    values = _one_dimensional_finite(signal, "signal")
    if values.size != time.size:
        raise ValueError("time and signal must have equal lengths")
    if time.size > 1 and np.any(np.diff(time) <= 0.0):
        raise ValueError("time must be strictly increasing")
    if np.any(values < 0.0):
        raise ValueError("sensor input signal must be non-negative")
    if t90_s is None:
        return values.copy()
    if not math.isfinite(t90_s) or t90_s <= 0.0:
        raise ValueError("t90_s must be finite and positive")
    if not math.isfinite(initial_value) or initial_value < 0.0:
        raise ValueError("initial sensor value must be finite and non-negative")
    tau = t90_s / math.log(10.0)
    response = np.empty_like(values)
    response[0] = float(initial_value)
    for index in range(1, values.size):
        decay = math.exp(-(time[index] - time[index - 1]) / tau)
        response[index] = values[index] + (
            response[index - 1] - values[index]
        ) * decay
    return response


def replay_fixed_receptors(
    history: WindHistory,
    receptors: Sequence[FixedReceptor],
    plume: Callable[[float, float, float, float], float],
    *,
    source_east_m: float = 0.0,
    source_north_m: float = 0.0,
    initial_indicated_mole_fraction: float = 0.0,
) -> tuple[ReceptorTrace, ...]:
    """Replay fixed receptors through an instantaneous steady wind frame.

    ``plume`` receives ``(wind_speed, alongwind, crosswind, height)``.  A
    receptor instantaneously upwind is assigned zero before instrument lag.
    This is an observation replay only: it contains no travel-time storage,
    gust response, source history, or dynamic plume deformation.
    """
    time, speed, direction = history.arrays()
    if not callable(plume):
        raise TypeError("plume must be callable")
    if not math.isfinite(source_east_m) or not math.isfinite(source_north_m):
        raise ValueError("source coordinates must be finite")
    rows = []
    for receptor in receptors:
        receptor.validate()
        east = receptor.east_m - source_east_m
        north = receptor.north_m - source_north_m
        angles = np.radians(np.mod(270.0 - direction, 360.0))
        down_east = np.cos(angles)
        down_north = np.sin(angles)
        alongwind = east * down_east + north * down_north
        crosswind = -east * down_north + north * down_east
        true = np.zeros_like(time)
        upwind_tolerance = 1.0e-12 * max(math.hypot(east, north), 1.0)
        for index in range(time.size):
            if alongwind[index] <= upwind_tolerance:
                continue
            value = float(plume(
                float(speed[index]),
                float(alongwind[index]),
                float(crosswind[index]),
                float(receptor.height_m),
            ))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError("plume returned an invalid mole fraction")
            true[index] = value
        indicated = first_order_sensor_response(
            time,
            true,
            t90_s=receptor.response_t90_s,
            initial_value=initial_indicated_mole_fraction,
        )
        rows.append(ReceptorTrace(
            receptor=receptor,
            time_s=time.copy(),
            wind_speed_m_s=speed.copy(),
            wind_direction_from_deg=direction.copy(),
            alongwind_m=alongwind,
            crosswind_m=crosswind,
            true_mole_fraction=true,
            indicated_mole_fraction=indicated,
        ))
    return tuple(rows)
