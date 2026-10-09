"""Three-dimensional transient dense-gas transport on a Cartesian FV mesh.

This module is an explicit, conservative reduced-order solver.  It is not a
replacement for a variable-density LES/RANS code, but it does solve the
time-dependent source--wind coupling that the steady plume paths cannot
represent:

* a finite atmospheric H2 source schedule is integrated at every step;
* the measured meteorological wind history is linearly interpolated at the
  same event-clock time and used as the horizontal advective velocity;
* the transported H2 inventory produces a bounded mixture-density anomaly;
* the density anomaly produces a local vertical buoyancy velocity;
* H2 is advected, diffused, sampled at receptors, and audited for mass
  conservation on a 3-D grid.

The solver deliberately keeps phase inference outside this module.  The
source schedule must already be an atmospheric H2 mass-rate boundary, and a
``SourceStateLedger`` can be handed off through
``SourceStateLedger.to_source_rate_schedule``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import numpy as np

from .semi_fv_obstacle import SourceRateSchedule
from .transient_receptor import WindHistory


TRANSIENT_DENSE_GAS_3D_SCHEMA = "degali.transient-dense-gas-3d.v1"


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_tuple(values: Sequence[object], name: str, length: int) -> tuple[float, ...]:
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} values")
    result = tuple(_finite(value, f"{name}[{index}]") for index, value in enumerate(values))
    if any(value <= 0.0 for value in result):
        raise ValueError(f"{name} values must be positive")
    return result


def _nonnegative_tuple(values: Sequence[object], name: str, length: int) -> tuple[float, ...]:
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} values")
    result = tuple(_finite(value, f"{name}[{index}]") for index, value in enumerate(values))
    if any(value < 0.0 for value in result):
        raise ValueError(f"{name} values must be non-negative")
    return result


@dataclass(frozen=True)
class DenseGas3DReceptor:
    """One fixed Cartesian receptor sampled from the transient field."""

    label: str
    position_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("receptor label must be non-empty")
        if len(self.position_m) != 3:
            raise ValueError("receptor position must contain three coordinates")
        values = tuple(_finite(value, f"receptor position[{index}]") for index, value in enumerate(self.position_m))
        object.__setattr__(self, "position_m", values)
        if any(value < 0.0 for value in values):
            raise ValueError("receptor coordinates must be non-negative")


@dataclass(frozen=True)
class TransientDenseGas3DConfig:
    """Numerical, source, ambient and buoyancy inputs for one 3-D run."""

    domain_m: tuple[float, float, float]
    cells: tuple[int, int, int]
    duration_s: float
    time_step_s: float
    source_schedule: SourceRateSchedule
    wind_history: WindHistory
    source_position_m: tuple[float, float, float]
    source_sigma_m: float
    source_density_kg_m3: float
    source_h2_mass_fraction: float
    ambient_density_kg_m3: float
    diffusivity_m2_s: float
    buoyancy_length_m: float = 1.0
    gravity_m_s2: float = 9.81
    maximum_buoyant_speed_m_s: float = 20.0
    store_fields: bool = False
    prognostic_velocity: bool = True
    wind_relaxation_time_s: float = 0.0
    momentum_diffusivity_m2_s: float | None = None
    store_velocity_fields: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "domain_m", _positive_tuple(self.domain_m, "domain_m", 3))
        if len(self.cells) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 2
            for value in self.cells
        ):
            raise ValueError("cells must contain three integers of at least two")
        object.__setattr__(self, "duration_s", _finite(self.duration_s, "duration_s"))
        object.__setattr__(self, "time_step_s", _finite(self.time_step_s, "time_step_s"))
        if self.duration_s <= 0.0 or self.time_step_s <= 0.0:
            raise ValueError("duration_s and time_step_s must be positive")
        if not isinstance(self.source_schedule, SourceRateSchedule):
            raise TypeError("source_schedule must be a SourceRateSchedule")
        self.source_schedule.require_zero_endpoint("3-D source schedule")
        if self.source_schedule.duration_s > self.duration_s + 1.0e-12:
            raise ValueError("source schedule cannot outlast the 3-D run")
        if not isinstance(self.wind_history, WindHistory):
            raise TypeError("wind_history must be a WindHistory")
        source_position = _nonnegative_tuple(self.source_position_m, "source_position_m", 3)
        object.__setattr__(self, "source_position_m", source_position)
        if any(value > bound for value, bound in zip(source_position, self.domain_m)):
            raise ValueError("source_position_m must lie inside the domain")
        for name in (
            "source_sigma_m", "source_density_kg_m3", "ambient_density_kg_m3",
            "diffusivity_m2_s", "buoyancy_length_m", "gravity_m_s2",
            "maximum_buoyant_speed_m_s",
        ):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.source_sigma_m <= 0.0 or self.source_density_kg_m3 <= 0.0:
            raise ValueError("source_sigma_m and source_density_kg_m3 must be positive")
        if not 0.0 < self.source_h2_mass_fraction <= 1.0:
            raise ValueError("source_h2_mass_fraction must lie in (0, 1]")
        if self.ambient_density_kg_m3 <= 0.0 or self.diffusivity_m2_s < 0.0:
            raise ValueError("ambient density must be positive and diffusivity non-negative")
        if self.buoyancy_length_m <= 0.0 or self.gravity_m_s2 <= 0.0:
            raise ValueError("buoyancy length and gravity must be positive")
        if self.maximum_buoyant_speed_m_s <= 0.0:
            raise ValueError("maximum_buoyant_speed_m_s must be positive")
        if not isinstance(self.store_fields, bool):
            raise TypeError("store_fields must be boolean")
        if not isinstance(self.prognostic_velocity, bool):
            raise TypeError("prognostic_velocity must be boolean")
        object.__setattr__(
            self,
            "wind_relaxation_time_s",
            _finite(self.wind_relaxation_time_s, "wind_relaxation_time_s"),
        )
        if self.wind_relaxation_time_s < 0.0:
            raise ValueError("wind_relaxation_time_s must be non-negative")
        if self.momentum_diffusivity_m2_s is not None:
            object.__setattr__(
                self,
                "momentum_diffusivity_m2_s",
                _finite(self.momentum_diffusivity_m2_s, "momentum_diffusivity_m2_s"),
            )
            if self.momentum_diffusivity_m2_s < 0.0:
                raise ValueError("momentum_diffusivity_m2_s must be non-negative")
        if not isinstance(self.store_velocity_fields, bool):
            raise TypeError("store_velocity_fields must be boolean")


@dataclass(frozen=True)
class TransientDenseGas3DResult:
    """Time history and conservation evidence from one 3-D run."""

    schema: str
    time_s: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    total_mass_kg: np.ndarray
    outflow_mass_kg: np.ndarray
    injected_mass_kg: np.ndarray
    mass_residual_kg: np.ndarray
    wind_vector_m_s: np.ndarray
    maximum_buoyant_speed_m_s: np.ndarray
    receptor_traces_kg_m3: Mapping[str, np.ndarray]
    fields_kg_m3: np.ndarray | None
    diagnostics: Mapping[str, object]
    mean_velocity_m_s: np.ndarray | None = None
    velocity_fields_m_s: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.schema != TRANSIENT_DENSE_GAS_3D_SCHEMA:
            raise ValueError("unsupported transient dense-gas schema")
        time = np.asarray(self.time_s, dtype=float)
        if time.ndim != 1 or time.size < 2 or not np.all(np.isfinite(time)):
            raise ValueError("time_s must be a finite one-dimensional history")
        if np.any(np.diff(time) <= 0.0):
            raise ValueError("time_s must be strictly increasing")
        object.__setattr__(self, "time_s", time)
        for name in ("x_m", "y_m", "z_m"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim != 1 or value.size < 2 or not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must be a finite grid vector")
            if np.any(np.diff(value) <= 0.0):
                raise ValueError(f"{name} must be strictly increasing")
            object.__setattr__(self, name, value)
        for name in (
            "total_mass_kg", "outflow_mass_kg", "injected_mass_kg",
            "mass_residual_kg", "maximum_buoyant_speed_m_s",
        ):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim != 1 or value.size != time.size or not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must be a finite vector matching time_s")
            object.__setattr__(self, name, value)
        wind = np.asarray(self.wind_vector_m_s, dtype=float)
        if wind.shape != (time.size, 3) or not np.all(np.isfinite(wind)):
            raise ValueError("wind_vector_m_s must have shape (time, 3)")
        object.__setattr__(self, "wind_vector_m_s", wind)
        for label, trace in self.receptor_traces_kg_m3.items():
            if not isinstance(label, str) or not label.strip():
                raise ValueError("receptor trace labels must be non-empty")
            values = np.asarray(trace, dtype=float)
            if values.shape != time.shape or not np.all(np.isfinite(values)):
                raise ValueError("receptor traces must be finite vectors matching time_s")
        if self.fields_kg_m3 is not None:
            fields = np.asarray(self.fields_kg_m3, dtype=float)
            if fields.ndim != 4 or fields.shape[0] != time.size or not np.all(np.isfinite(fields)):
                raise ValueError("fields_kg_m3 must have shape (time, z, y, x)")
            object.__setattr__(self, "fields_kg_m3", fields)
        if self.mean_velocity_m_s is not None:
            mean_velocity = np.asarray(self.mean_velocity_m_s, dtype=float)
            if mean_velocity.shape != (time.size, 3) or not np.all(np.isfinite(mean_velocity)):
                raise ValueError("mean_velocity_m_s must have shape (time, 3)")
            object.__setattr__(self, "mean_velocity_m_s", mean_velocity)
        if self.velocity_fields_m_s is not None:
            velocity_fields = np.asarray(self.velocity_fields_m_s, dtype=float)
            if (
                velocity_fields.ndim != 5
                or velocity_fields.shape[0] != time.size
                or velocity_fields.shape[1] != 3
                or not np.all(np.isfinite(velocity_fields))
            ):
                raise ValueError("velocity_fields_m_s must have shape (time, 3, z, y, x)")
            object.__setattr__(self, "velocity_fields_m_s", velocity_fields)

    @property
    def maximum_mass_residual_kg(self) -> float:
        return float(np.max(np.abs(self.mass_residual_kg)))

    def as_record(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "time_s": self.time_s.tolist(),
            "total_mass_kg": self.total_mass_kg.tolist(),
            "outflow_mass_kg": self.outflow_mass_kg.tolist(),
            "injected_mass_kg": self.injected_mass_kg.tolist(),
            "mass_residual_kg": self.mass_residual_kg.tolist(),
            "wind_vector_m_s": self.wind_vector_m_s.tolist(),
            "maximum_buoyant_speed_m_s": self.maximum_buoyant_speed_m_s.tolist(),
            "receptor_labels": sorted(self.receptor_traces_kg_m3),
            "stored_fields": self.fields_kg_m3 is not None,
            "stored_velocity_fields": self.velocity_fields_m_s is not None,
            "diagnostics": dict(self.diagnostics),
        }


def _wind_vector_at(history: WindHistory, time_s: float) -> np.ndarray:
    times, speeds, directions_from = history.arrays()
    if time_s < times[0] - 1.0e-12 or time_s > times[-1] + 1.0e-12:
        raise ValueError("wind history does not cover the requested event time")
    speed = float(np.interp(time_s, times, speeds))
    direction_radians = np.radians(directions_from)
    # Meteorological bearings are clockwise from north: 270° means the wind
    # comes from west, so the transport vector points east (+x).
    from_unit = np.column_stack((np.sin(direction_radians), np.cos(direction_radians)))
    from_x = float(np.interp(time_s, times, from_unit[:, 0]))
    from_y = float(np.interp(time_s, times, from_unit[:, 1]))
    norm = math.hypot(from_x, from_y)
    if norm <= 1.0e-12:
        raise ValueError("wind direction interpolation became undefined")
    # Meteorological direction is where the wind comes from; transport uses
    # the opposite (to) vector.
    return np.array((-speed * from_x / norm, -speed * from_y / norm, 0.0))


def _source_weights(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    position: tuple[float, float, float],
    sigma: float,
) -> np.ndarray:
    xx, yy, zz = np.meshgrid(x, y, z, indexing="xy")
    weights = np.exp(-0.5 * (
        ((xx - position[0]) / sigma) ** 2
        + ((yy - position[1]) / sigma) ** 2
        + ((zz - position[2]) / sigma) ** 2
    ))
    # meshgrid(indexing='xy') produces (ny,nx,nz); transpose to (nz,ny,nx).
    weights = np.transpose(weights, (2, 0, 1))
    total = float(np.sum(weights))
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("source Gaussian has no positive grid weight")
    return weights / total


def _advect_uniform_axis(
    field: np.ndarray,
    velocity: float,
    dt: float,
    spacing: float,
    axis: int,
) -> tuple[np.ndarray, float]:
    """Conservative first-order upwind advection with vacuum boundaries."""
    if abs(velocity) <= 0.0:
        return field, 0.0
    shape = list(field.shape)
    shape[axis] += 1
    flux = np.zeros(shape, dtype=float)
    source_slice = [slice(None)] * field.ndim
    face_slice = [slice(None)] * field.ndim
    if velocity > 0.0:
        source_slice[axis] = slice(None, -1)
        face_slice[axis] = slice(1, -1)
        flux[tuple(face_slice)] = velocity * field[tuple(source_slice)]
        outlet = [slice(None)] * field.ndim
        outlet[axis] = -1
        interior = [slice(None)] * field.ndim
        interior[axis] = -1
        flux[tuple(outlet)] = velocity * field[tuple(interior)]
    else:
        source_slice[axis] = slice(1, None)
        face_slice[axis] = slice(1, -1)
        flux[tuple(face_slice)] = velocity * field[tuple(source_slice)]
        inlet = [slice(None)] * field.ndim
        inlet[axis] = 0
        interior = [slice(None)] * field.ndim
        interior[axis] = 0
        flux[tuple(inlet)] = velocity * field[tuple(interior)]
    increment = -(dt / spacing) * np.diff(flux, axis=axis)
    outgoing = float(np.sum(np.maximum(flux.take(indices=-1, axis=axis), 0.0)))
    if velocity < 0.0:
        outgoing += float(np.sum(np.minimum(flux.take(indices=0, axis=axis), 0.0) * -1.0))
    # The first/last-face flux units are kg m^-2 s^-1; area factors are
    # applied by the caller through the cell-volume ratio.  Return only the
    # directional flux integral proxy for diagnostics.
    return field + increment, outgoing * dt


def _advect_variable_axis(
    field: np.ndarray,
    velocity: np.ndarray,
    dt: float,
    spacing: float,
    axis: int,
) -> tuple[np.ndarray, float]:
    """Conservative upwind advection for a cell-centred velocity field."""
    if velocity.shape != field.shape:
        raise ValueError("variable advection velocity must match the scalar field shape")
    face_shape = list(field.shape)
    face_shape[axis] += 1
    flux = np.zeros(face_shape, dtype=float)
    lower = [slice(None)] * field.ndim
    upper = [slice(None)] * field.ndim
    interior = [slice(None)] * field.ndim
    lower[axis] = slice(None, -1)
    upper[axis] = slice(1, None)
    interior[axis] = slice(1, -1)
    face_velocity = 0.5 * (velocity[tuple(lower)] + velocity[tuple(upper)])
    flux[tuple(interior)] = np.where(
        face_velocity >= 0.0,
        face_velocity * field[tuple(lower)],
        face_velocity * field[tuple(upper)],
    )
    left_face = [slice(None)] * field.ndim
    right_face = [slice(None)] * field.ndim
    left_face[axis] = 0
    right_face[axis] = -1
    left_cell = [slice(None)] * field.ndim
    right_cell = [slice(None)] * field.ndim
    left_cell[axis] = 0
    right_cell[axis] = -1
    # Vacuum boundaries: only outward flow carries scalar out; inward flow
    # sees zero concentration outside the domain.
    flux[tuple(left_face)] = np.minimum(velocity[tuple(left_cell)], 0.0) * field[tuple(left_cell)]
    flux[tuple(right_face)] = np.maximum(velocity[tuple(right_cell)], 0.0) * field[tuple(right_cell)]
    result = field - (dt / spacing) * np.diff(flux, axis=axis)
    outflow = (
        np.sum(np.maximum(flux[tuple(right_face)], 0.0))
        - np.sum(np.minimum(flux[tuple(left_face)], 0.0))
    ) * dt
    return result, float(outflow)


def _advect_vertical(
    field: np.ndarray,
    vertical_velocity: np.ndarray,
    dt: float,
    spacing: float,
) -> tuple[np.ndarray, float]:
    """Conservative upwind advection for a spatially varying vertical speed."""
    nz, ny, nx = field.shape
    flux = np.zeros((nz + 1, ny, nx), dtype=float)
    w_face = 0.5 * (vertical_velocity[:-1] + vertical_velocity[1:])
    flux[1:-1] = np.where(w_face >= 0.0, w_face * field[:-1], w_face * field[1:])
    flux[0] = np.minimum(vertical_velocity[0], 0.0) * field[0]
    flux[-1] = np.maximum(vertical_velocity[-1], 0.0) * field[-1]
    result = field - (dt / spacing) * np.diff(flux, axis=0)
    outflow = float(
        np.sum(np.maximum(flux[-1], 0.0))
        - np.sum(np.minimum(flux[0], 0.0))
    ) * dt
    return result, outflow


def _diffuse_neumann(field: np.ndarray, diffusivity: float, dt: float, spacings: tuple[float, float, float]) -> np.ndarray:
    if diffusivity <= 0.0:
        return field
    dx, dy, dz = spacings
    lap = np.zeros_like(field)
    for axis, spacing in zip((2, 1, 0), (dx, dy, dz)):
        plus = np.take(field, indices=range(1, field.shape[axis]), axis=axis)
        minus = np.take(field, indices=range(field.shape[axis] - 1), axis=axis)
        centre = field
        plus_full = np.concatenate((np.take(field, [0], axis=axis), plus), axis=axis)
        minus_full = np.concatenate((minus, np.take(field, [-1], axis=axis)), axis=axis)
        lap += (plus_full - 2.0 * centre + minus_full) / spacing**2
    return field + diffusivity * dt * lap


def _sample_trilinear(
    field: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    position: tuple[float, float, float],
) -> float:
    indices = []
    fractions = []
    for coordinates, value in zip((x, y, z), position):
        if value < coordinates[0] or value > coordinates[-1]:
            return 0.0
        right = int(np.clip(np.searchsorted(coordinates, value), 1, len(coordinates) - 1))
        left = right - 1
        span = coordinates[right] - coordinates[left]
        indices.append((left, right))
        fractions.append(0.0 if span <= 0.0 else (value - coordinates[left]) / span)
    value = 0.0
    for iz in (0, 1):
        for iy in (0, 1):
            for ix in (0, 1):
                weight = (
                    (fractions[2] if iz else 1.0 - fractions[2])
                    * (fractions[1] if iy else 1.0 - fractions[1])
                    * (fractions[0] if ix else 1.0 - fractions[0])
                )
                value += weight * field[
                    indices[2][iz], indices[1][iy], indices[0][ix]
                ]
    return float(value)


def solve_transient_dense_gas_3d(
    config: TransientDenseGas3DConfig,
    *,
    receptors: Sequence[DenseGas3DReceptor] = (),
) -> TransientDenseGas3DResult:
    """Run a source-history/wind-history coupled 3-D transient calculation."""
    if not isinstance(config, TransientDenseGas3DConfig):
        raise TypeError("config must be a TransientDenseGas3DConfig")
    if any(not isinstance(item, DenseGas3DReceptor) for item in receptors):
        raise TypeError("receptors must contain DenseGas3DReceptor values")
    labels = [item.label for item in receptors]
    if len(set(labels)) != len(labels):
        raise ValueError("receptor labels must be unique")
    wind_times, _speeds, _directions = config.wind_history.arrays()
    if wind_times[0] > 0.0 or wind_times[-1] < config.duration_s:
        raise ValueError("wind history must cover the complete solver duration")

    nx, ny, nz = config.cells
    lx, ly, lz = config.domain_m
    dx, dy, dz = lx / nx, ly / ny, lz / nz
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    z = (np.arange(nz) + 0.5) * dz
    cell_volume = dx * dy * dz
    weights = _source_weights(x, y, z, config.source_position_m, config.source_sigma_m)
    steps = int(math.ceil(config.duration_s / config.time_step_s))
    dt = config.duration_s / steps
    max_wind = float(np.max(_speeds))
    max_density_fraction = 1.0
    max_buoyancy = min(
        config.maximum_buoyant_speed_m_s,
        math.sqrt(
            config.gravity_m_s2 * config.buoyancy_length_m
            * abs(config.source_density_kg_m3 - config.ambient_density_kg_m3)
            * max_density_fraction / config.ambient_density_kg_m3
        ),
    )
    cfl = (max_wind * dt / dx) + (max_wind * dt / dy) + (max_buoyancy * dt / dz)
    diffusion_number = 2.0 * config.diffusivity_m2_s * dt * (
        1.0 / dx**2 + 1.0 / dy**2 + 1.0 / dz**2
    )
    if cfl > 1.0 + 1.0e-12:
        raise ValueError(f"3-D advection CFL exceeds one: {cfl:.6g}")
    if diffusion_number > 1.0 + 1.0e-12:
        raise ValueError(f"3-D diffusion stability number exceeds one: {diffusion_number:.6g}")

    field = np.zeros((nz, ny, nx), dtype=float)
    volume = cell_volume
    times = [0.0]
    total_mass = [0.0]
    outflow_mass = [0.0]
    injected_mass = [0.0]
    residuals = [0.0]
    wind_vectors = [_wind_vector_at(config.wind_history, 0.0)]
    max_buoyant_history = [0.0]
    initial_wind = _wind_vector_at(config.wind_history, 0.0)
    velocity_u = np.full_like(field, initial_wind[0], dtype=float)
    velocity_v = np.full_like(field, initial_wind[1], dtype=float)
    velocity_w = np.zeros_like(field)
    mean_velocity_history = [
        np.array((float(np.mean(velocity_u)), float(np.mean(velocity_v)), 0.0))
    ]
    stored_velocity_fields = (
        [np.stack((velocity_u.copy(), velocity_v.copy(), velocity_w.copy()))]
        if config.store_velocity_fields else None
    )
    receptor_traces = {item.label: [0.0] for item in receptors}
    stored_fields = [field.copy()] if config.store_fields else None
    cumulative_outflow = 0.0
    cumulative_injected = 0.0

    for index in range(steps):
        start = index * dt
        end = (index + 1) * dt
        wind = _wind_vector_at(config.wind_history, 0.5 * (start + end))
        source_mass = config.source_schedule.mass_between(start, end)
        if source_mass > 0.0:
            field += source_mass * weights / volume
        cumulative_injected += source_mass

        source_fraction = np.clip(
            field / (config.source_density_kg_m3 * config.source_h2_mass_fraction),
            0.0,
            1.0,
        )
        density = config.ambient_density_kg_m3 + (
            config.source_density_kg_m3 - config.ambient_density_kg_m3
        ) * source_fraction
        buoyancy_acceleration = (
            config.gravity_m_s2
            * (config.ambient_density_kg_m3 - density)
            / config.ambient_density_kg_m3
        )
        if config.prognostic_velocity:
            if config.wind_relaxation_time_s <= 0.0:
                velocity_u.fill(float(wind[0]))
                velocity_v.fill(float(wind[1]))
            else:
                relaxation = min(1.0, dt / config.wind_relaxation_time_s)
                velocity_u += relaxation * (float(wind[0]) - velocity_u)
                velocity_v += relaxation * (float(wind[1]) - velocity_v)
            velocity_w += dt * buoyancy_acceleration
            velocity_w = np.clip(
                velocity_w,
                -config.maximum_buoyant_speed_m_s,
                config.maximum_buoyant_speed_m_s,
            )
            momentum_diffusivity = (
                config.diffusivity_m2_s
                if config.momentum_diffusivity_m2_s is None
                else config.momentum_diffusivity_m2_s
            )
            if momentum_diffusivity > 0.0:
                velocity_u = _diffuse_neumann(velocity_u, momentum_diffusivity, dt, (dx, dy, dz))
                velocity_v = _diffuse_neumann(velocity_v, momentum_diffusivity, dt, (dx, dy, dz))
                velocity_w = _diffuse_neumann(velocity_w, momentum_diffusivity, dt, (dx, dy, dz))
            field, out_x = _advect_variable_axis(field, velocity_u, dt, dx, axis=2)
            field, out_y = _advect_variable_axis(field, velocity_v, dt, dy, axis=1)
            field, out_z = _advect_variable_axis(field, velocity_w, dt, dz, axis=0)
        else:
            velocity_u.fill(float(wind[0]))
            velocity_v.fill(float(wind[1]))
            velocity_w = np.sign(config.ambient_density_kg_m3 - density) * np.sqrt(
                config.gravity_m_s2 * config.buoyancy_length_m
                * np.abs(density - config.ambient_density_kg_m3)
                / config.ambient_density_kg_m3
            )
            velocity_w = np.clip(
                velocity_w,
                -config.maximum_buoyant_speed_m_s,
                config.maximum_buoyant_speed_m_s,
            )
            field, out_x = _advect_uniform_axis(field, float(wind[0]), dt, dx, axis=2)
            field, out_y = _advect_uniform_axis(field, float(wind[1]), dt, dy, axis=1)
            field, out_z = _advect_vertical(field, velocity_w, dt, dz)
        field = _diffuse_neumann(field, config.diffusivity_m2_s, dt, (dx, dy, dz))
        minimum = float(np.min(field))
        if minimum < -1.0e-10:
            raise RuntimeError(f"3-D transient transport became negative: {minimum:.6g}")
        field = np.maximum(field, 0.0)
        step_outflow = (out_x * dy * dz) + (out_y * dx * dz) + (out_z * dx * dy)
        cumulative_outflow += step_outflow
        now = end
        mass = float(np.sum(field) * volume)
        residual = mass + cumulative_outflow - cumulative_injected
        times.append(now)
        total_mass.append(mass)
        outflow_mass.append(cumulative_outflow)
        injected_mass.append(cumulative_injected)
        residuals.append(residual)
        wind_vectors.append(_wind_vector_at(config.wind_history, now))
        max_buoyant_history.append(float(np.max(np.abs(velocity_w))))
        mean_velocity_history.append(np.array((
            float(np.mean(velocity_u)),
            float(np.mean(velocity_v)),
            float(np.mean(velocity_w)),
        )))
        for receptor in receptors:
            receptor_traces[receptor.label].append(
                _sample_trilinear(field, x, y, z, receptor.position_m)
            )
        if stored_fields is not None:
            stored_fields.append(field.copy())
        if stored_velocity_fields is not None:
            stored_velocity_fields.append(
                np.stack((velocity_u.copy(), velocity_v.copy(), velocity_w.copy()))
            )

    return TransientDenseGas3DResult(
        schema=TRANSIENT_DENSE_GAS_3D_SCHEMA,
        time_s=np.asarray(times),
        x_m=x,
        y_m=y,
        z_m=z,
        total_mass_kg=np.asarray(total_mass),
        outflow_mass_kg=np.asarray(outflow_mass),
        injected_mass_kg=np.asarray(injected_mass),
        mass_residual_kg=np.asarray(residuals),
        wind_vector_m_s=np.asarray(wind_vectors),
        maximum_buoyant_speed_m_s=np.asarray(max_buoyant_history),
        receptor_traces_kg_m3={label: np.asarray(values) for label, values in receptor_traces.items()},
        fields_kg_m3=None if stored_fields is None else np.asarray(stored_fields),
        diagnostics={
            "cells": [nx, ny, nz],
            "cell_size_m": [dx, dy, dz],
            "time_step_s": dt,
            "cfl": cfl,
            "diffusion_number": diffusion_number,
            "source_id": config.source_schedule.source_id,
            "wind_history_start_s": float(wind_times[0]),
            "wind_history_end_s": float(wind_times[-1]),
            "maximum_mass_residual_kg": float(np.max(np.abs(residuals))),
            "prognostic_velocity": config.prognostic_velocity,
            "wind_relaxation_time_s": config.wind_relaxation_time_s,
            "solver_scope": "3-D explicit reduced-order finite-volume buoyant H2 transport with prognostic velocity closure; not LES/RANS",
        },
        mean_velocity_m_s=np.asarray(mean_velocity_history),
        velocity_fields_m_s=(
            None if stored_velocity_fields is None else np.asarray(stored_velocity_fields)
        ),
    )


__all__ = [
    "TRANSIENT_DENSE_GAS_3D_SCHEMA",
    "DenseGas3DReceptor",
    "TransientDenseGas3DConfig",
    "TransientDenseGas3DResult",
    "solve_transient_dense_gas_3d",
]
