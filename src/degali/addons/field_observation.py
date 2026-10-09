"""Observation operator for reduced-order field concentration traces."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from .field_contracts import SensorModel
from .transient_receptor import first_order_sensor_response


HYDROGEN_MOLAR_MASS_KG_MOL = 2.01588e-3
DRY_AIR_MOLAR_MASS_KG_MOL = 28.9652e-3


@dataclass(frozen=True)
class FieldSensorTrace:
    """True and indicated H2 mole-fraction histories at one field sensor."""

    sensor: SensorModel
    time_s: np.ndarray
    true_mole_fraction: np.ndarray
    response_mole_fraction: np.ndarray
    indicated_mole_fraction: np.ndarray
    ambient_air_density_kg_m3: float

    def __post_init__(self) -> None:
        if not isinstance(self.sensor, SensorModel):
            raise TypeError("field sensor trace sensor must be a SensorModel")
        arrays = {
            "time_s": self.time_s,
            "true_mole_fraction": self.true_mole_fraction,
            "response_mole_fraction": self.response_mole_fraction,
            "indicated_mole_fraction": self.indicated_mole_fraction,
        }
        for name, value in arrays.items():
            if not isinstance(value, np.ndarray):
                raise TypeError(f"field sensor trace {name} must be a numpy array")
            array = np.asarray(value, dtype=float)
            if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
                raise ValueError(f"field sensor trace {name} must be a finite non-empty one-dimensional array")
            if array.size > 1 and name == "time_s" and np.any(np.diff(array) <= 0.0):
                raise ValueError("field sensor trace time_s must be strictly increasing")
        lengths = {np.asarray(value).size for value in arrays.values()}
        if len(lengths) != 1:
            raise ValueError("field sensor trace arrays must have equal lengths")
        for name in ("true_mole_fraction", "response_mole_fraction"):
            values = np.asarray(arrays[name], dtype=float)
            if np.any(values < 0.0) or np.any(values > 1.0):
                raise ValueError(f"field sensor trace {name} must lie in [0, 1]")
        if isinstance(self.ambient_air_density_kg_m3, bool) or not math.isfinite(
            float(self.ambient_air_density_kg_m3)
        ) or self.ambient_air_density_kg_m3 <= 0.0:
            raise ValueError("field sensor trace ambient air density must be positive and finite")


def h2_mole_fraction_from_mass_concentration(
    concentration_kg_m3: Sequence[float],
    *,
    ambient_air_density_kg_m3: float,
) -> np.ndarray:
    """Convert a scalar H2 mass concentration using declared ambient air density.

    The semi-FV operator transports H2 inventory rather than a full mixture
    equation of state. This conversion is therefore an explicit ideal-volume
    observation approximation, not an EOS claim for the cold cloud interior.
    """
    concentration = np.asarray(concentration_kg_m3, dtype=float)
    if concentration.ndim != 1 or not np.all(np.isfinite(concentration)):
        raise ValueError("concentration must be a finite one-dimensional sequence")
    if np.any(concentration < 0.0):
        raise ValueError("concentration cannot be negative")
    density = float(ambient_air_density_kg_m3)
    if not math.isfinite(density) or density <= 0.0:
        raise ValueError("ambient_air_density_kg_m3 must be positive and finite")
    h2_moles = concentration / HYDROGEN_MOLAR_MASS_KG_MOL
    air_moles = density / DRY_AIR_MOLAR_MASS_KG_MOL
    return h2_moles / (h2_moles + air_moles)


def trailing_time_average(
    time_s: Sequence[float], signal: Sequence[float], *, window_s: float,
) -> np.ndarray:
    """Return a trailing piecewise-linear time average without resampling."""
    time = np.asarray(time_s, dtype=float)
    values = np.asarray(signal, dtype=float)
    if time.ndim != 1 or values.ndim != 1 or time.size == 0 or time.size != values.size:
        raise ValueError("time and signal must be non-empty equal-length one-dimensional sequences")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(values)):
        raise ValueError("time and signal must be finite")
    if time.size > 1 and np.any(np.diff(time) <= 0.0):
        raise ValueError("time must be strictly increasing")
    if not math.isfinite(float(window_s)) or window_s < 0.0:
        raise ValueError("window_s must be finite and non-negative")
    if window_s == 0.0:
        return values.copy()
    output = np.empty_like(values)
    for index, end in enumerate(time):
        start = max(float(time[0]), float(end - window_s))
        if end == start:
            output[index] = values[index]
            continue
        interior = time[(time > start) & (time < end)]
        knots = np.concatenate(([start], interior, [end]))
        samples = np.interp(knots, time, values)
        integral = np.sum(0.5 * (samples[1:] + samples[:-1]) * np.diff(knots))
        output[index] = integral / (end - start)
    return output


def apply_sensor_model(
    time_s: Sequence[float],
    concentration_kg_m3: Sequence[float],
    sensor: SensorModel,
    *,
    ambient_air_density_kg_m3: float,
    initial_indicated_mole_fraction: float = 0.0,
) -> FieldSensorTrace:
    """Apply declared t90, gain, bias and averaging to one true H2 trace."""
    time = np.asarray(time_s, dtype=float)
    true = h2_mole_fraction_from_mass_concentration(
        concentration_kg_m3,
        ambient_air_density_kg_m3=ambient_air_density_kg_m3,
    )
    if time.shape != true.shape:
        raise ValueError("time and concentration must have equal lengths")
    response = first_order_sensor_response(
        time,
        true,
        t90_s=sensor.response_time_s.nominal,
        initial_value=initial_indicated_mole_fraction,
    )
    calibrated = (
        sensor.gain.nominal * response + sensor.bias_mole_fraction.nominal
    )
    indicated = trailing_time_average(
        time, calibrated, window_s=sensor.averaging_time_s,
    )
    return FieldSensorTrace(
        sensor=sensor,
        time_s=time.copy(),
        true_mole_fraction=true,
        response_mole_fraction=response,
        indicated_mole_fraction=indicated,
        ambient_air_density_kg_m3=float(ambient_air_density_kg_m3),
    )


__all__ = [
    "DRY_AIR_MOLAR_MASS_KG_MOL", "HYDROGEN_MOLAR_MASS_KG_MOL",
    "FieldSensorTrace", "apply_sensor_model", "h2_mole_fraction_from_mass_concentration",
    "trailing_time_average",
]
