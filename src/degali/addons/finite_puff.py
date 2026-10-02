"""Native conservative finite Gaussian puff continuation for DEGALI.

This module continues :mod:`degali.addons.finite_release` after source
cessation.  It is a full-inventory, three-dimensional model: marked cloud
mass grows by ambient entrainment, hydrogen mass is conserved, ambient
momentum relaxes the puff toward the wind, buoyancy changes vertical momentum,
and a ground-contacting dense puff spreads in both horizontal directions.

The closure coefficients are declared in :class:`GaussianPuffConfig`.  Their
defaults are the published integral-model values already used elsewhere in
DEGALI/SLABx; they are not fitted to an LH2 concentration dataset.  The model
does not resolve obstacles, terrain, stochastic meander, combustion, or a
time-varying wind field.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

import numpy as np

from .finite_release import FiniteReleasePuffHandoff
from .transient_receptor import (
    WindHistory,
    meteorological_from_to_math_radians,
)


GRAVITY_M_S2 = 9.80665


@dataclass(frozen=True)
class GaussianPuffConfig:
    """Declared closure and integration controls for one puff continuation."""

    duration_s: float
    time_step_s: float = 0.05
    wind_velocity_m_s: tuple[float, float, float] = (1.0, 0.0, 0.0)
    shear_entrainment_coefficient: float = 0.1
    ambient_turbulence_coefficient: float = 0.1
    gravity_spreading_coefficient: float = 0.75
    form_drag_coefficient: float = 0.2
    width_cutoff_sigma: float = 3.0
    ground_heat_input_w: float = 0.0
    wind_history: WindHistory | None = None
    wind_history_time_offset_s: float = 0.0

    def validate(self) -> None:
        positive = {
            "duration_s": self.duration_s,
            "time_step_s": self.time_step_s,
            "width_cutoff_sigma": self.width_cutoff_sigma,
        }
        for name, value in positive.items():
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        nonnegative = {
            "shear_entrainment_coefficient": self.shear_entrainment_coefficient,
            "ambient_turbulence_coefficient": self.ambient_turbulence_coefficient,
            "gravity_spreading_coefficient": self.gravity_spreading_coefficient,
            "form_drag_coefficient": self.form_drag_coefficient,
            "ground_heat_input_w": self.ground_heat_input_w,
        }
        for name, value in nonnegative.items():
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
        wind = np.asarray(self.wind_velocity_m_s, dtype=float)
        if wind.shape != (3,) or not np.all(np.isfinite(wind)):
            raise ValueError("wind_velocity_m_s must contain three finite values")
        if (
            not math.isfinite(self.wind_history_time_offset_s)
            or self.wind_history_time_offset_s < 0.0
        ):
            raise ValueError(
                "wind_history_time_offset_s must be finite and non-negative"
            )
        if self.wind_history is not None:
            time, _speed, _direction = self.wind_history.arrays()
            start = self.wind_history_time_offset_s
            end = start + self.duration_s
            tolerance = 1.0e-12 * max(abs(start), abs(end), 1.0)
            if start < time[0] - tolerance or end > time[-1] + tolerance:
                raise ValueError(
                    "wind history must cover the complete puff time window"
                )


@dataclass(frozen=True)
class GaussianPuffState:
    """One time-resolved full-inventory puff state."""

    elapsed_s: float
    centre_position_m: tuple[float, float, float]
    bulk_velocity_m_s: tuple[float, float, float]
    longitudinal_half_length_m: float
    lateral_half_width_m: float
    vertical_half_width_m: float
    total_mass_kg: float
    hydrogen_mass_kg: float
    momentum_kg_m_s: tuple[float, float, float]
    relative_total_energy_j: float
    bulk_temperature_k: float
    bulk_density_kg_m3: float
    bulk_hydrogen_mass_fraction: float
    cumulative_entrained_air_kg: float
    ground_contact: bool


@dataclass(frozen=True)
class PuffReceptorTrace:
    """Concentration and temperature history at one fixed global point."""

    point_m: tuple[float, float, float]
    time_s: np.ndarray
    mole_fraction: np.ndarray
    temperature_k: np.ndarray

    @property
    def peak_mole_fraction(self) -> float:
        return float(np.max(self.mole_fraction))

    @property
    def peak_time_s(self) -> float:
        return float(self.time_s[int(np.argmax(self.mole_fraction))])


@dataclass
class GaussianPuffResult:
    """Native puff trajectory with conservative receptor observation methods."""

    states: tuple[GaussianPuffState, ...]
    handoff: FiniteReleasePuffHandoff
    config: GaussianPuffConfig
    thermodynamics: Any
    shape_volume_factor: float
    maximum_relative_mass_residual: float
    maximum_relative_hydrogen_residual: float

    def state_at(self, elapsed_s: float) -> GaussianPuffState:
        """Linearly interpolate a state; extrapolation is rejected."""
        value = float(elapsed_s)
        times = np.array([state.elapsed_s for state in self.states])
        if not math.isfinite(value) or value < times[0] or value > times[-1]:
            raise ValueError("puff time lies outside the integrated trajectory")
        upper = int(np.searchsorted(times, value, side="left"))
        if upper == 0:
            return self.states[0]
        if upper == len(times):
            return self.states[-1]
        lower = upper - 1
        if times[upper] == value:
            return self.states[upper]
        fraction = (value - times[lower]) / (times[upper] - times[lower])
        return _blend_state(self.states[lower], self.states[upper], fraction)

    def observation_at(
        self, elapsed_s: float, point_m: Sequence[float]
    ) -> tuple[float, float]:
        """Return H2 mole fraction and temperature at a fixed global point."""
        point = np.asarray(point_m, dtype=float)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("point_m must contain three finite values")
        if point[2] < 0.0:
            raise ValueError("receptor height cannot be negative")
        state = self.state_at(elapsed_s)
        centre = np.asarray(state.centre_position_m)
        velocity = np.asarray(state.bulk_velocity_m_s)
        horizontal = velocity[:2]
        if np.linalg.norm(horizontal) <= 1.0e-12:
            horizontal = _wind_velocity_at(
                self.config, state.elapsed_s
            )[:2]
        if np.linalg.norm(horizontal) <= 1.0e-12:
            along = np.array([1.0, 0.0])
        else:
            along = horizontal / np.linalg.norm(horizontal)
        cross = np.array([-along[1], along[0]])
        offset = point - centre
        local_x = float(offset[:2] @ along)
        local_y = float(offset[:2] @ cross)
        cutoff = self.config.width_cutoff_sigma
        sigma_x = state.longitudinal_half_length_m / cutoff
        sigma_y = state.lateral_half_width_m / cutoff
        sigma_z = state.vertical_half_width_m / cutoff
        horizontal_shape = math.exp(-0.5 * (
            (local_x / sigma_x) ** 2 + (local_y / sigma_y) ** 2
        ))
        direct = math.exp(-0.5 * ((point[2] - centre[2]) / sigma_z) ** 2)
        image = math.exp(-0.5 * ((point[2] + centre[2]) / sigma_z) ** 2)
        shape = horizontal_shape * (direct + image)
        gaussian_volume = (2.0 * math.pi) ** 1.5 * sigma_x * sigma_y * sigma_z
        hydrogen_density = state.hydrogen_mass_kg / gaussian_volume * shape
        if hydrogen_density <= np.finfo(float).tiny:
            return 0.0, float(self.thermodynamics.ambient_temperature)
        thermal_shape = min(shape, 1.0)
        local_temperature = (
            self.thermodynamics.ambient_temperature
            + (state.bulk_temperature_k - self.thermodynamics.ambient_temperature)
            * thermal_shape
        )
        mass_fraction = _mass_fraction_from_hydrogen_density(
            self.thermodynamics, hydrogen_density, local_temperature
        )
        numerator = mass_fraction / self.thermodynamics.fuel_molecular_weight
        denominator = numerator + (
            (1.0 - mass_fraction)
            / self.thermodynamics._humid_ambient_molecular_weight
        )
        return float(numerator / denominator), float(local_temperature)

    def receptor_trace(self, point_m: Sequence[float]) -> PuffReceptorTrace:
        values = [
            self.observation_at(state.elapsed_s, point_m)
            for state in self.states
        ]
        return PuffReceptorTrace(
            point_m=tuple(float(value) for value in point_m),
            time_s=np.array([state.elapsed_s for state in self.states]),
            mole_fraction=np.array([value[0] for value in values]),
            temperature_k=np.array([value[1] for value in values]),
        )


def _thermodynamic_state(
    thermodynamics: Any,
    mass: float,
    hydrogen_mass: float,
    momentum: np.ndarray,
    total_energy: float,
) -> tuple[float, float, float]:
    fraction = hydrogen_mass / mass
    kinetic = float(momentum @ momentum) / (2.0 * mass)
    relative_enthalpy = (total_energy - kinetic) / mass
    enthalpy = thermodynamics._ambient_enthalpy + relative_enthalpy
    temperature = float(
        thermodynamics._temperature_from_enthalpy(enthalpy, fraction)
    )
    density = float(
        thermodynamics._density_from_temperature(temperature, fraction)
    )
    if not all(math.isfinite(value) and value > 0.0 for value in (temperature, density)):
        raise RuntimeError("puff thermodynamic reconstruction is non-physical")
    return temperature, density, fraction


def _ellipsoid_surface_area(a: float, b: float, c: float) -> float:
    # Knud Thomsen approximation; exact for a sphere and below 1.1% error for
    # all ellipsoids used by the integral geometry.
    p = 1.6075
    return 4.0 * math.pi * (
        ((a * b) ** p + (a * c) ** p + (b * c) ** p) / 3.0
    ) ** (1.0 / p)


def _state(
    *, elapsed: float, centre: np.ndarray, momentum: np.ndarray,
    half_x: float, half_y: float, half_z: float, mass: float,
    hydrogen: float, energy: float, temperature: float, density: float,
    entrained: float,
) -> GaussianPuffState:
    velocity = momentum / mass
    return GaussianPuffState(
        elapsed_s=float(elapsed),
        centre_position_m=tuple(float(value) for value in centre),
        bulk_velocity_m_s=tuple(float(value) for value in velocity),
        longitudinal_half_length_m=float(half_x),
        lateral_half_width_m=float(half_y),
        vertical_half_width_m=float(half_z),
        total_mass_kg=float(mass),
        hydrogen_mass_kg=float(hydrogen),
        momentum_kg_m_s=tuple(float(value) for value in momentum),
        relative_total_energy_j=float(energy),
        bulk_temperature_k=float(temperature),
        bulk_density_kg_m3=float(density),
        bulk_hydrogen_mass_fraction=float(hydrogen / mass),
        cumulative_entrained_air_kg=float(entrained),
        ground_contact=bool(centre[2] <= half_z * (1.0 + 1.0e-10)),
    )


def integrate_finite_gaussian_puff(
    handoff: FiniteReleasePuffHandoff,
    thermodynamics: Any,
    config: GaussianPuffConfig,
) -> GaussianPuffResult:
    """Continue an accepted finite-release handoff through a native 3-D puff."""

    config.validate()
    if handoff.status != "transition_ready":
        raise ValueError("a transition-ready finite-release handoff is required")
    required = (
        handoff.centre_position_m, handoff.transverse_widths_m,
        handoff.longitudinal_half_length_m, handoff.total_mass_kg,
        handoff.hydrogen_mass_kg, handoff.momentum_kg_m_s,
        handoff.relative_energy_j,
    )
    if any(value is None for value in required):
        raise ValueError("finite-release handoff is incomplete")

    centre = np.asarray(handoff.centre_position_m, dtype=float)
    momentum = np.asarray(handoff.momentum_kg_m_s, dtype=float)
    mass = float(handoff.total_mass_kg)
    initial_mass = mass
    hydrogen = float(handoff.hydrogen_mass_kg)
    energy = float(handoff.relative_energy_j)
    half_x = float(handoff.longitudinal_half_length_m)
    half_y = config.width_cutoff_sigma * float(handoff.transverse_widths_m[0])
    half_z = config.width_cutoff_sigma * float(handoff.transverse_widths_m[1])
    temperature, density, _fraction = _thermodynamic_state(
        thermodynamics, mass, hydrogen, momentum, energy
    )
    shape_factor = mass / (density * half_x * half_y * half_z)
    if not math.isfinite(shape_factor) or shape_factor <= 0.0:
        raise RuntimeError("invalid puff geometry-volume mapping")
    if centre[2] < half_z:
        centre[2] = half_z
    elapsed = 0.0
    entrained = 0.0
    states = [_state(
        elapsed=elapsed, centre=centre, momentum=momentum,
        half_x=half_x, half_y=half_y, half_z=half_z, mass=mass,
        hydrogen=hydrogen, energy=energy, temperature=temperature,
        density=density, entrained=entrained,
    )]
    mass_residuals = [0.0]
    hydrogen_residuals = [0.0]

    while elapsed < config.duration_s - 1.0e-12:
        dt = min(config.time_step_s, config.duration_s - elapsed)
        if config.wind_history is not None:
            history_time = config.wind_history.arrays()[0]
            absolute_time = config.wind_history_time_offset_s + elapsed
            future = history_time[history_time > absolute_time + 1.0e-12]
            if future.size:
                dt = min(dt, float(future[0] - absolute_time))
        wind = _wind_velocity_at(config, elapsed + 0.5 * dt)
        velocity = momentum / mass
        relative = wind - velocity
        relative_speed = float(np.linalg.norm(relative))
        wind_speed = float(np.linalg.norm(wind))
        entrainment_velocity = (
            config.shear_entrainment_coefficient * relative_speed
            + config.ambient_turbulence_coefficient * wind_speed
        )
        surface = _ellipsoid_surface_area(half_x, half_y, half_z)
        ground_contact = centre[2] <= half_z * (1.0 + 1.0e-10)
        exposed_surface = surface * (0.5 if ground_contact else 1.0)
        dm = (
            thermodynamics.ambient_density * exposed_surface
            * entrainment_velocity * dt
        )
        volume = mass / density
        buoyancy_force = (
            (thermodynamics.ambient_density - density)
            * GRAVITY_M_S2 * volume
        )
        projected_areas = math.pi * np.array([
            half_y * half_z, half_x * half_z, half_x * half_y,
        ])
        drag_force = (
            0.5 * config.form_drag_coefficient
            * thermodynamics.ambient_density
            * projected_areas * relative * np.abs(relative)
        )
        force = drag_force.copy()
        force[2] += buoyancy_force
        old_velocity = velocity.copy()
        momentum = momentum + dm * wind + force * dt
        mass += dm
        entrained += dm
        energy += (
            0.5 * float(wind @ wind) * dm
            + float(force @ old_velocity) * dt
            + (config.ground_heat_input_w * dt if ground_contact else 0.0)
        )
        reduced_gravity = (
            GRAVITY_M_S2
            * max(density - thermodynamics.ambient_density, 0.0)
            / thermodynamics.ambient_density
        )
        gravity_spread = (
            config.gravity_spreading_coefficient
            * math.sqrt(max(2.0 * reduced_gravity * half_z, 0.0))
            if ground_contact else 0.0
        )
        half_x += (entrainment_velocity + gravity_spread) * dt
        half_y += (entrainment_velocity + gravity_spread) * dt
        temperature, density, _fraction = _thermodynamic_state(
            thermodynamics, mass, hydrogen, momentum, energy
        )
        half_z = mass / (density * shape_factor * half_x * half_y)
        if not math.isfinite(half_z) or half_z <= 0.0:
            raise RuntimeError("puff vertical scale became non-physical")
        new_velocity = momentum / mass
        centre = centre + 0.5 * (old_velocity + new_velocity) * dt
        if centre[2] < half_z:
            centre[2] = half_z
            if momentum[2] < 0.0:
                momentum[2] = 0.0
                new_velocity = momentum / mass
        elapsed += dt
        states.append(_state(
            elapsed=elapsed, centre=centre, momentum=momentum,
            half_x=half_x, half_y=half_y, half_z=half_z, mass=mass,
            hydrogen=hydrogen, energy=energy, temperature=temperature,
            density=density, entrained=entrained,
        ))
        mass_residuals.append(mass - initial_mass - entrained)
        hydrogen_residuals.append(states[-1].hydrogen_mass_kg - hydrogen)

    mass_scale = max(initial_mass, 1.0e-30)
    hydrogen_scale = max(hydrogen, 1.0e-30)
    return GaussianPuffResult(
        states=tuple(states), handoff=handoff, config=config,
        thermodynamics=thermodynamics, shape_volume_factor=float(shape_factor),
        maximum_relative_mass_residual=float(
            max(abs(value) for value in mass_residuals) / mass_scale
        ),
        maximum_relative_hydrogen_residual=float(
            max(abs(value) for value in hydrogen_residuals) / hydrogen_scale
        ),
    )


def _mass_fraction_from_hydrogen_density(
    thermodynamics: Any, hydrogen_density: float, temperature: float
) -> float:
    def residual(fraction: float) -> float:
        density = float(
            thermodynamics._density_from_temperature(temperature, fraction)
        )
        return density * fraction - hydrogen_density

    upper = residual(1.0)
    if upper < -1.0e-12:
        raise RuntimeError(
            "Gaussian puff requests an H2 vapour density above the local "
            "pressure/temperature closure"
        )
    low, high = 0.0, 1.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if residual(middle) < 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _wind_velocity_at(
    config: GaussianPuffConfig, elapsed_s: float
) -> np.ndarray:
    if config.wind_history is None:
        return np.asarray(config.wind_velocity_m_s, dtype=float)
    time, speed, direction = config.wind_history.arrays()
    angles = np.array([
        meteorological_from_to_math_radians(float(value))
        for value in direction
    ])
    east = speed * np.cos(angles)
    north = speed * np.sin(angles)
    absolute_time = config.wind_history_time_offset_s + elapsed_s
    return np.array([
        np.interp(absolute_time, time, east),
        np.interp(absolute_time, time, north),
        0.0,
    ])


def _blend_state(
    first: GaussianPuffState, second: GaussianPuffState, fraction: float
) -> GaussianPuffState:
    def scalar(name: str) -> float:
        return float(
            getattr(first, name)
            + fraction * (getattr(second, name) - getattr(first, name))
        )

    def vector(name: str) -> tuple[float, float, float]:
        a = np.asarray(getattr(first, name), dtype=float)
        b = np.asarray(getattr(second, name), dtype=float)
        return tuple(float(value) for value in a + fraction * (b - a))

    return GaussianPuffState(
        elapsed_s=scalar("elapsed_s"),
        centre_position_m=vector("centre_position_m"),
        bulk_velocity_m_s=vector("bulk_velocity_m_s"),
        longitudinal_half_length_m=scalar("longitudinal_half_length_m"),
        lateral_half_width_m=scalar("lateral_half_width_m"),
        vertical_half_width_m=scalar("vertical_half_width_m"),
        total_mass_kg=scalar("total_mass_kg"),
        hydrogen_mass_kg=scalar("hydrogen_mass_kg"),
        momentum_kg_m_s=vector("momentum_kg_m_s"),
        relative_total_energy_j=scalar("relative_total_energy_j"),
        bulk_temperature_k=scalar("bulk_temperature_k"),
        bulk_density_kg_m3=scalar("bulk_density_kg_m3"),
        bulk_hydrogen_mass_fraction=scalar("bulk_hydrogen_mass_fraction"),
        cumulative_entrained_air_kg=scalar("cumulative_entrained_air_kg"),
        ground_contact=first.ground_contact if fraction < 0.5 else second.ground_contact,
    )


__all__ = [
    "GaussianPuffConfig", "GaussianPuffState", "GaussianPuffResult",
    "PuffReceptorTrace", "integrate_finite_gaussian_puff",
]
