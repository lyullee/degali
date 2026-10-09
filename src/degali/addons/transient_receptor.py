"""Research-only transient observation operator for steady plume libraries.

The dispersion models in :mod:`degali` are spatial, steady integral models.
This module does not turn them into a transient CFD model.  It provides the
smaller operations needed for field-test comparison: rotate fixed receptors
into an instantaneous wind frame, or carry finite source-history packets by
the integrated wind vector; sample caller-supplied kernels; and apply a
declared first-order instrument response.

No source multiplier, meander width, or response time is inferred from an
observation.  Every time series and every steady response must be supplied by
the caller.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Literal, Sequence

import numpy as np


SteadySampler = Callable[[float, float, float], float]
PacketSampler = Callable[["SourcePacket", float, float, float, float], float]


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
class SourcePacket:
    """One finite source-history packet passed to a caller-owned transport kernel.

    ``release_time_s`` and ``duration_s`` use the same event clock as
    :class:`WindHistory`.  Pressure and temperature are provenance/state
    values; this module does not convert them to a discharge rate or perform
    a flash calculation.  ``source_direction_to_deg`` is retained separately
    from the ambient wind direction so a fixed nozzle axis is not silently
    rotated with the weather.
    """

    release_time_s: float
    duration_s: float
    mass_kg: float
    pressure_pa: float | None = None
    temperature_k: float | None = None
    source_direction_to_deg: float = 0.0

    def validate(self) -> None:
        for name, value in {
            "release_time_s": self.release_time_s,
            "duration_s": self.duration_s,
            "mass_kg": self.mass_kg,
            "source_direction_to_deg": self.source_direction_to_deg,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"packet {name} must be finite and numeric")
        if self.release_time_s < 0.0:
            raise ValueError("packet release_time_s must be non-negative")
        if self.duration_s <= 0.0 or self.mass_kg <= 0.0:
            raise ValueError("packet duration_s and mass_kg must be positive")
        for name, value in {
            "pressure_pa": self.pressure_pa,
            "temperature_k": self.temperature_k,
        }.items():
            if value is not None and (
                isinstance(value, bool) or not math.isfinite(float(value))
                or float(value) <= 0.0
            ):
                raise ValueError(f"packet {name} must be positive and finite when supplied")


@dataclass(frozen=True)
class SourceHistory:
    """Measured/declared source clock converted to finite packet quadrature.

    ``mass_rate_kg_s`` is already an atmospheric H2 mass-rate boundary.  The
    history must end with an explicit zero rate before it can be transported.
    Pressure and temperature, when present, are carried into each packet but
    are never used here to invent a flow or phase closure.  Packet midpoint
    quadrature is deterministic and can be refined with
    ``to_packets(subdivisions_per_interval=...)``.
    """

    time_s: Sequence[float]
    mass_rate_kg_s: Sequence[float]
    pressure_pa: Sequence[float] | None = None
    temperature_k: Sequence[float] | None = None
    source_direction_to_deg: Sequence[float] | None = None
    rate_operator: Literal["piecewise_constant", "linear"] = "piecewise_constant"

    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray]:
        time = _one_dimensional_finite(self.time_s, "source_time_s")
        rate = _one_dimensional_finite(self.mass_rate_kg_s, "mass_rate_kg_s")
        if time.size < 2 or rate.size != time.size:
            raise ValueError("source history needs equal time/rate arrays of at least two values")
        if abs(float(time[0])) > 1.0e-12:
            raise ValueError("source history must start at event time zero")
        if np.any(np.diff(time) <= 0.0):
            raise ValueError("source history time must be strictly increasing")
        if np.any(rate < 0.0):
            raise ValueError("source mass rate must be non-negative")
        if rate[-1] != 0.0:
            raise ValueError("source history must end with an explicit zero rate")
        if not np.any(rate[:-1] > 0.0):
            raise ValueError("source history must contain positive released mass")
        if self.rate_operator not in {"piecewise_constant", "linear"}:
            raise ValueError("source rate_operator must be piecewise_constant or linear")

        def optional_channel(values: Sequence[float] | None, name: str) -> np.ndarray | None:
            if values is None:
                return None
            array = _one_dimensional_finite(values, name)
            if array.size != time.size:
                raise ValueError(f"{name} must have the same length as source time")
            if np.any(array <= 0.0):
                raise ValueError(f"{name} must be positive")
            return array

        pressure = optional_channel(self.pressure_pa, "pressure_pa")
        temperature = optional_channel(self.temperature_k, "temperature_k")
        if (pressure is None) != (temperature is None):
            raise ValueError("pressure_pa and temperature_k must be supplied together")
        if self.source_direction_to_deg is None:
            direction = np.zeros_like(time)
        else:
            direction = _one_dimensional_finite(
                self.source_direction_to_deg, "source_direction_to_deg"
            )
            if direction.size != time.size:
                raise ValueError("source_direction_to_deg must have the same length as source time")
            direction = np.mod(direction, 360.0)
        return time, rate, pressure, temperature, direction

    def to_packets(self, *, subdivisions_per_interval: int = 1) -> tuple[SourcePacket, ...]:
        """Return positive-mass midpoint packets with exactly conserved mass."""
        if isinstance(subdivisions_per_interval, bool) or not isinstance(
            subdivisions_per_interval, int
        ) or subdivisions_per_interval <= 0:
            raise ValueError("subdivisions_per_interval must be a positive integer")
        time, rate, pressure, temperature, direction = self.arrays()
        packets: list[SourcePacket] = []

        def circular_midpoint(left: float, right: float, fraction: float) -> float:
            left_rad = math.radians(left)
            right_rad = left_rad + math.atan2(
                math.sin(math.radians(right) - left_rad),
                math.cos(math.radians(right) - left_rad),
            )
            return math.degrees(left_rad + fraction * (right_rad - left_rad)) % 360.0

        for index, (left_time, right_time) in enumerate(zip(time[:-1], time[1:])):
            span = float(right_time - left_time)
            for subdivision in range(subdivisions_per_interval):
                fraction_left = subdivision / subdivisions_per_interval
                fraction_right = (subdivision + 1) / subdivisions_per_interval
                start = float(left_time + span * fraction_left)
                end = float(left_time + span * fraction_right)
                interval = end - start
                if self.rate_operator == "piecewise_constant":
                    left_rate = right_rate = float(rate[index])
                else:
                    left_rate = float(rate[index] + (rate[index + 1] - rate[index]) * fraction_left)
                    right_rate = float(rate[index] + (rate[index + 1] - rate[index]) * fraction_right)
                mass = 0.5 * (left_rate + right_rate) * interval
                if mass <= 0.0:
                    continue
                midpoint_fraction = 0.5 * (fraction_left + fraction_right)
                packet_pressure = None if pressure is None else float(
                    pressure[index] + (pressure[index + 1] - pressure[index]) * midpoint_fraction
                )
                packet_temperature = None if temperature is None else float(
                    temperature[index] + (temperature[index + 1] - temperature[index]) * midpoint_fraction
                )
                packets.append(SourcePacket(
                    release_time_s=0.5 * (start + end),
                    duration_s=interval,
                    mass_kg=float(mass),
                    pressure_pa=packet_pressure,
                    temperature_k=packet_temperature,
                    source_direction_to_deg=circular_midpoint(
                        float(direction[index]), float(direction[index + 1]), midpoint_fraction
                    ),
                ))
        if not packets:
            raise ValueError("source history produced no positive-mass packets")
        for packet in packets:
            packet.validate()
        expected = self.released_mass_kg()
        actual = math.fsum(packet.mass_kg for packet in packets)
        if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15):
            raise RuntimeError("source packet quadrature did not conserve declared mass")
        return tuple(packets)

    def released_mass_kg(self) -> float:
        time, rate, _pressure, _temperature, _direction = self.arrays()
        if self.rate_operator == "linear":
            return float(math.fsum(
                0.5 * (float(left) + float(right)) * float(end - start)
                for start, end, left, right in zip(time[:-1], time[1:], rate[:-1], rate[1:])
            ))
        return float(math.fsum(
            float(value) * float(end - start)
            for start, end, value in zip(time[:-1], time[1:], rate[:-1])
        ))


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


def _wind_vector_at(
    time_s: float,
    time: np.ndarray,
    speed: np.ndarray,
    direction_from_deg: np.ndarray,
) -> np.ndarray:
    """Linearly sample the horizontal wind vector on the common clock."""
    value = float(time_s)
    tolerance = 1.0e-12 * max(1.0, abs(float(time[0])), abs(float(time[-1])))
    if value < time[0] - tolerance or value > time[-1] + tolerance:
        raise ValueError("wind history does not cover the requested packet time")
    value = min(max(value, float(time[0])), float(time[-1]))
    transport_angle = np.unwrap(np.radians((270.0 - direction_from_deg) % 360.0))
    angle = float(np.interp(value, time, transport_angle))
    magnitude = float(np.interp(value, time, speed))
    return np.array([magnitude * math.cos(angle), magnitude * math.sin(angle)])


def _wind_displacement(
    start_s: float,
    end_s: float,
    time: np.ndarray,
    speed: np.ndarray,
    direction_from_deg: np.ndarray,
) -> np.ndarray:
    """Integrate the declared piecewise-linear wind vector by trapezoids."""
    if end_s < start_s:
        raise ValueError("wind displacement interval must be ordered")
    if end_s == start_s:
        return np.zeros(2, dtype=float)
    interior = time[(time > start_s) & (time < end_s)]
    nodes = np.concatenate((np.asarray([start_s]), interior, np.asarray([end_s])))
    vectors = np.vstack([
        _wind_vector_at(node, time, speed, direction_from_deg) for node in nodes
    ])
    return np.sum(0.5 * (vectors[:-1] + vectors[1:]) * np.diff(nodes)[:, None], axis=0)


def _observation_wind(
    observation_time_s: np.ndarray,
    time: np.ndarray,
    speed: np.ndarray,
    direction_from_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    transport_angle = np.unwrap(np.radians((270.0 - direction_from_deg) % 360.0))
    sampled_angle = np.interp(observation_time_s, time, transport_angle)
    sampled_direction = np.mod(270.0 - np.degrees(sampled_angle), 360.0)
    return np.interp(observation_time_s, time, speed), sampled_direction


def replay_fixed_receptors_with_source_history(
    source_history: SourceHistory,
    wind_history: WindHistory,
    receptors: Sequence[FixedReceptor],
    packet_sampler: PacketSampler,
    *,
    observation_time_s: Sequence[float] | None = None,
    source_east_m: float = 0.0,
    source_north_m: float = 0.0,
    initial_indicated_mole_fraction: float = 0.0,
    packet_subdivisions_per_interval: int = 1,
) -> tuple[ReceptorTrace, ...]:
    """Replay finite source-history packets with causal wind travel memory.

    The packet centre is advected by the integrated measured wind from its
    release midpoint.  ``packet_sampler`` receives
    ``(packet, age_s, relative_east_m, relative_north_m, receptor_height_m)``
    and returns that packet's mole-fraction contribution. Contributions are
    summed, bounded to ``[0, 1]``, and then passed through the declared sensor
    response. This is an observation/transport operator, not a transient CFD
    replacement or a source/parameter-fitting routine.
    """
    if not isinstance(source_history, SourceHistory):
        raise TypeError("source_history must be a SourceHistory")
    if not isinstance(wind_history, WindHistory):
        raise TypeError("wind_history must be a WindHistory")
    if not callable(packet_sampler):
        raise TypeError("packet_sampler must be callable")
    if not math.isfinite(source_east_m) or not math.isfinite(source_north_m):
        raise ValueError("source coordinates must be finite")
    wind_time, wind_speed, wind_direction = wind_history.arrays()
    packets = source_history.to_packets(
        subdivisions_per_interval=packet_subdivisions_per_interval
    )
    if observation_time_s is None:
        observation = wind_time.copy()
    else:
        observation = _one_dimensional_finite(observation_time_s, "observation_time_s")
        if observation.size > 1 and np.any(np.diff(observation) <= 0.0):
            raise ValueError("observation_time_s must be strictly increasing")
    tolerance = 1.0e-12 * max(1.0, abs(float(wind_time[0])), abs(float(wind_time[-1])))
    if observation[0] < wind_time[0] - tolerance or observation[-1] > wind_time[-1] + tolerance:
        raise ValueError("wind history must cover every observation time")
    if any(
        packet.release_time_s < wind_time[0] - tolerance
        or packet.release_time_s > wind_time[-1] + tolerance
        for packet in packets
    ):
        raise ValueError("wind history must cover every source packet release time")
    observation = np.clip(observation, wind_time[0], wind_time[-1])
    observed_speed, observed_direction = _observation_wind(
        observation, wind_time, wind_speed, wind_direction
    )
    rows = []
    for receptor in receptors:
        receptor.validate()
        true = np.zeros_like(observation)
        receptor_position = np.asarray(
            [receptor.east_m - source_east_m, receptor.north_m - source_north_m],
            dtype=float,
        )
        for index, current_time in enumerate(observation):
            total = 0.0
            for packet in packets:
                if current_time < packet.release_time_s:
                    continue
                displacement = _wind_displacement(
                    packet.release_time_s, float(current_time),
                    wind_time, wind_speed, wind_direction,
                )
                relative = receptor_position - displacement
                contribution = float(packet_sampler(
                    packet,
                    float(current_time - packet.release_time_s),
                    float(relative[0]),
                    float(relative[1]),
                    float(receptor.height_m),
                ))
                if not math.isfinite(contribution) or contribution < 0.0:
                    raise ValueError("packet sampler returned an invalid mole fraction")
                total += contribution
            if total > 1.0 + 1.0e-12:
                raise ValueError("packet contributions exceed physical mole fraction one")
            true[index] = min(total, 1.0)
        indicated = first_order_sensor_response(
            observation,
            true,
            t90_s=receptor.response_t90_s,
            initial_value=initial_indicated_mole_fraction,
        )
        rows.append(ReceptorTrace(
            receptor=receptor,
            time_s=observation.copy(),
            wind_speed_m_s=observed_speed.copy(),
            wind_direction_from_deg=observed_direction.copy(),
            # A packet cloud can have multiple local transport directions;
            # these fields are intentionally unavailable rather than guessed.
            alongwind_m=np.full_like(observation, np.nan),
            crosswind_m=np.full_like(observation, np.nan),
            true_mole_fraction=true,
            indicated_mole_fraction=indicated,
        ))
    return tuple(rows)


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
