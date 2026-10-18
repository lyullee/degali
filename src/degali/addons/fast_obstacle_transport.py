"""Fast local obstacle transport without non-local scalar routing.

This module is the first tier of the obstacle-physics rebuild.  A sparse
Laplace solve produces a divergence-free, wind-plane potential-flow field
around one ground-mounted rectangular obstacle.  Hydrogen is then transported
with conservative local face fluxes.  Solid faces have exactly zero advective
and diffusive scalar flux.

The model captures windward deceleration and continuous over-roof diversion.
Potential flow is inviscid and slip-walled, so it does *not* resolve a
separated roof shear layer, a recirculating lee wake, wall friction, lateral
three-dimensional bypass, or variable-density momentum feedback.  Those
belong to the pressure/Brinkman and three-dimensional tiers.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import time
from typing import Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .semi_fv_obstacle import RectangularObstacle2D, SourceRateSchedule


FAST_OBSTACLE_TRANSPORT_SCHEMA = "degali.fast-obstacle-transport.v1"
R_GAS_J_MOL_K = 8.314462618
HYDROGEN_MOLAR_MASS_KG_MOL = 2.01588e-3


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class FastObstacleReceptor:
    """One receptor in the wind-aligned x-z plane."""

    label: str
    x_m: float
    z_m: float

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("receptor label must be non-empty")
        object.__setattr__(self, "x_m", _finite(self.x_m, "receptor x_m"))
        object.__setattr__(self, "z_m", _finite(self.z_m, "receptor z_m"))


@dataclass(frozen=True)
class FastObstacleTransportConfig:
    """Numerical and physical inputs for the Tier-1 wind-plane screen.

    ``source_schedule`` is atmospheric hydrogen mass rate per metre of span.
    The two-dimensional inventory consequently has units kg/m.
    """

    length_m: float
    height_m: float
    nx: int
    nz: int
    duration_s: float
    time_step_s: float
    wind_speed_m_s: float
    diffusivity_m2_s: float
    source_schedule: SourceRateSchedule
    source_x_m: float
    source_height_m: float
    source_sigma_m: float
    obstacle: RectangularObstacle2D | None = None
    ambient_temperature_k: float = 293.15
    ambient_pressure_pa: float = 101_325.0
    store_fields: bool = False

    def __post_init__(self) -> None:
        for name in (
            "length_m", "height_m", "duration_s", "time_step_s",
            "wind_speed_m_s", "source_sigma_m", "ambient_temperature_k",
            "ambient_pressure_pa",
        ):
            value = _finite(getattr(self, name), name)
            object.__setattr__(self, name, value)
            if value <= 0.0:
                raise ValueError(f"{name} must be positive")
        object.__setattr__(
            self, "diffusivity_m2_s",
            _finite(self.diffusivity_m2_s, "diffusivity_m2_s"),
        )
        if self.diffusivity_m2_s < 0.0:
            raise ValueError("diffusivity_m2_s must be non-negative")
        for name in ("source_x_m", "source_height_m"):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if not (0.0 <= self.source_x_m <= self.length_m):
            raise ValueError("source_x_m must lie inside the domain")
        if not (0.0 <= self.source_height_m <= self.height_m):
            raise ValueError("source_height_m must lie inside the domain")
        if (
            isinstance(self.nx, bool) or isinstance(self.nz, bool)
            or not isinstance(self.nx, int) or not isinstance(self.nz, int)
            or self.nx < 8 or self.nz < 8
        ):
            raise ValueError("nx and nz must be integers of at least eight")
        if not isinstance(self.source_schedule, SourceRateSchedule):
            raise TypeError("source_schedule must be a SourceRateSchedule")
        self.source_schedule.require_zero_endpoint("fast obstacle source schedule")
        if self.source_schedule.duration_s > self.duration_s + 1.0e-12:
            raise ValueError("source schedule cannot outlast the transport run")
        if self.obstacle is not None:
            if not isinstance(self.obstacle, RectangularObstacle2D):
                raise TypeError("obstacle must be RectangularObstacle2D or None")
            if abs(self.obstacle.base_height_m) > 1.0e-12:
                raise ValueError(
                    "Tier-1 potential flow currently requires a ground-mounted obstacle"
                )
            if self.obstacle.x_min_m <= 0.0 or self.obstacle.x_max_m >= self.length_m:
                raise ValueError("obstacle must lie strictly inside the x domain")
            if self.obstacle.z_max_m >= self.height_m:
                raise ValueError("obstacle must leave a fluid path below the top boundary")
            if (
                self.obstacle.x_min_m <= self.source_x_m <= self.obstacle.x_max_m
                and self.source_height_m <= self.obstacle.z_max_m
            ):
                raise ValueError("source may not lie inside the obstacle")
        if not isinstance(self.store_fields, bool):
            raise TypeError("store_fields must be boolean")


@dataclass(frozen=True)
class FastObstacleTransportResult:
    """Velocity, pressure proxy, concentration and audit histories."""

    schema: str
    x_m: np.ndarray
    z_m: np.ndarray
    streamfunction_m2_s: np.ndarray
    u_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    velocity_u_m_s: np.ndarray
    velocity_w_m_s: np.ndarray
    pressure_coefficient: np.ndarray
    solid_mask: np.ndarray
    time_s: np.ndarray
    concentration_kg_m3: np.ndarray
    concentration_fields_kg_m3: np.ndarray | None
    receptor_traces_kg_m3: dict[str, np.ndarray]
    injected_mass_per_span_kg_m: np.ndarray
    domain_mass_per_span_kg_m: np.ndarray
    outflow_mass_per_span_kg_m: np.ndarray
    mass_residual_per_span_kg_m: np.ndarray
    diagnostics: dict[str, object]
    ambient_temperature_k: float
    ambient_pressure_pa: float

    @property
    def concentration_vol_percent(self) -> np.ndarray:
        molar_density = self.ambient_pressure_pa / (
            R_GAS_J_MOL_K * self.ambient_temperature_k
        )
        result = (
            self.concentration_kg_m3
            / HYDROGEN_MOLAR_MASS_KG_MOL
            / molar_density
            * 100.0
        )
        return np.clip(result, 0.0, 100.0)

    @property
    def maximum_mass_residual_per_span_kg_m(self) -> float:
        return float(np.max(np.abs(self.mass_residual_per_span_kg_m)))


def _masks(
    config: FastObstacleTransportConfig,
    x_nodes: np.ndarray,
    z_nodes: np.ndarray,
    x_cells: np.ndarray,
    z_cells: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if config.obstacle is None:
        return (
            np.zeros((len(z_nodes), len(x_nodes)), dtype=bool),
            np.zeros((len(z_cells), len(x_cells)), dtype=bool),
        )
    obstacle = config.obstacle
    cell_solid = (
        (x_cells[None, :] >= obstacle.x_min_m)
        & (x_cells[None, :] <= obstacle.x_max_m)
        & (z_cells[:, None] <= obstacle.z_max_m)
    )
    # Snap the declared rectangle to complete Cartesian cells, then mark all
    # nodes belonging to those cells.  A constant streamfunction on that
    # complete node boundary makes the normal velocity exactly zero and keeps
    # the discrete curl divergence-free in adjacent fluid cells.
    node_solid = np.zeros((len(z_nodes), len(x_nodes)), dtype=bool)
    solid_j, solid_i = np.where(cell_solid)
    if solid_i.size:
        node_solid[
            int(solid_j.min()):int(solid_j.max()) + 2,
            int(solid_i.min()):int(solid_i.max()) + 2,
        ] = True
    return node_solid, cell_solid


def _solve_streamfunction(
    config: FastObstacleTransportConfig,
) -> tuple[np.ndarray, ...]:
    """Solve Laplace's equation for a ground-obstacle streamfunction."""
    x_nodes = np.linspace(0.0, config.length_m, config.nx + 1)
    z_nodes = np.linspace(0.0, config.height_m, config.nz + 1)
    x_cells = 0.5 * (x_nodes[:-1] + x_nodes[1:])
    z_cells = 0.5 * (z_nodes[:-1] + z_nodes[1:])
    dx = config.length_m / config.nx
    dz = config.height_m / config.nz
    node_solid, cell_solid = _masks(
        config, x_nodes, z_nodes, x_cells, z_cells,
    )

    fixed = np.zeros_like(node_solid)
    fixed[0, :] = True
    fixed[-1, :] = True
    fixed[:, 0] = True
    fixed[:, -1] = True
    fixed |= node_solid
    psi = config.wind_speed_m_s * z_nodes[:, None] * np.ones_like(x_nodes)[None, :]
    psi[node_solid] = 0.0

    unknown_positions = np.argwhere(~fixed)
    lookup = -np.ones_like(node_solid, dtype=int)
    for index, (j, i) in enumerate(unknown_positions):
        lookup[j, i] = index
    count = len(unknown_positions)
    matrix = sparse.lil_matrix((count, count), dtype=float)
    rhs = np.zeros(count, dtype=float)
    coefficients = (
        (-1, 0, 1.0 / dz**2),
        (1, 0, 1.0 / dz**2),
        (0, -1, 1.0 / dx**2),
        (0, 1, 1.0 / dx**2),
    )
    diagonal = 2.0 / dx**2 + 2.0 / dz**2
    for row, (j, i) in enumerate(unknown_positions):
        matrix[row, row] = diagonal
        for dj, di, coefficient in coefficients:
            neighbour_j, neighbour_i = j + dj, i + di
            neighbour_index = lookup[neighbour_j, neighbour_i]
            if neighbour_index >= 0:
                matrix[row, neighbour_index] = -coefficient
            else:
                rhs[row] += coefficient * psi[neighbour_j, neighbour_i]
    if count:
        psi[~fixed] = spsolve(matrix.tocsr(), rhs)

    # The streamfunction gives face-normal velocities with an exactly
    # divergence-free Cartesian curl in every fluid cell.
    u_face = (psi[1:, :] - psi[:-1, :]) / dz
    w_face = -(psi[:, 1:] - psi[:, :-1]) / dx

    # Enforce impermeability exactly even when an obstacle edge falls between
    # nominal coordinates and its mask is snapped to the Cartesian grid.
    if config.nx > 1:
        blocked_x = cell_solid[:, :-1] | cell_solid[:, 1:]
        u_face[:, 1:-1][blocked_x] = 0.0
    u_face[:, 0][cell_solid[:, 0]] = 0.0
    u_face[:, -1][cell_solid[:, -1]] = 0.0
    if config.nz > 1:
        blocked_z = cell_solid[:-1, :] | cell_solid[1:, :]
        w_face[1:-1, :][blocked_z] = 0.0
    w_face[0, cell_solid[0, :]] = 0.0
    w_face[-1, cell_solid[-1, :]] = 0.0

    u_cell = 0.5 * (u_face[:, :-1] + u_face[:, 1:])
    w_cell = 0.5 * (w_face[:-1, :] + w_face[1:, :])
    u_cell[cell_solid] = 0.0
    w_cell[cell_solid] = 0.0
    speed_squared = u_cell**2 + w_cell**2
    pressure_coefficient = 1.0 - speed_squared / config.wind_speed_m_s**2
    pressure_coefficient[cell_solid] = np.nan
    return (
        x_cells, z_cells, psi, u_face, w_face, u_cell, w_cell,
        pressure_coefficient, cell_solid,
    )


def _source_weights(
    config: FastObstacleTransportConfig,
    x: np.ndarray,
    z: np.ndarray,
    solid: np.ndarray,
) -> np.ndarray:
    xx, zz = np.meshgrid(x, z)
    weights = np.exp(-0.5 * (
        ((xx - config.source_x_m) / config.source_sigma_m) ** 2
        + ((zz - config.source_height_m) / config.source_sigma_m) ** 2
    ))
    weights[solid] = 0.0
    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError("source Gaussian has no fluid cell")
    return weights / total


def _transport_step(
    concentration: np.ndarray,
    *,
    u_face: np.ndarray,
    w_face: np.ndarray,
    solid: np.ndarray,
    diffusivity: float,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, float]:
    """One conservative local advection--diffusion step."""
    nz, nx = concentration.shape
    flux_x = np.zeros((nz, nx + 1), dtype=float)
    interior_u = u_face[:, 1:-1]
    flux_x[:, 1:-1] = np.where(
        interior_u >= 0.0,
        interior_u * concentration[:, :-1],
        interior_u * concentration[:, 1:],
    )
    # Zero-concentration exterior at the inlet; convective outflow elsewhere.
    flux_x[:, 0] = np.minimum(u_face[:, 0], 0.0) * concentration[:, 0]
    flux_x[:, -1] = np.maximum(u_face[:, -1], 0.0) * concentration[:, -1]
    blocked_x = solid[:, :-1] | solid[:, 1:]
    flux_x[:, 1:-1][blocked_x] = 0.0

    flux_z = np.zeros((nz + 1, nx), dtype=float)
    interior_w = w_face[1:-1, :]
    flux_z[1:-1, :] = np.where(
        interior_w >= 0.0,
        interior_w * concentration[:-1, :],
        interior_w * concentration[1:, :],
    )
    # Ground and top are streamfunction boundaries and therefore impermeable
    # in this wind-plane screen.
    blocked_z = solid[:-1, :] | solid[1:, :]
    flux_z[1:-1, :][blocked_z] = 0.0

    result = concentration.copy()
    result -= (dt / dx) * (flux_x[:, 1:] - flux_x[:, :-1])
    result -= (dt / dz) * (flux_z[1:, :] - flux_z[:-1, :])

    if diffusivity > 0.0:
        diff_x = np.zeros_like(flux_x)
        diff_x[:, 1:-1] = -diffusivity * (
            concentration[:, 1:] - concentration[:, :-1]
        ) / dx
        diff_x[:, 1:-1][blocked_x] = 0.0
        diff_z = np.zeros_like(flux_z)
        diff_z[1:-1, :] = -diffusivity * (
            concentration[1:, :] - concentration[:-1, :]
        ) / dz
        diff_z[1:-1, :][blocked_z] = 0.0
        result -= (dt / dx) * (diff_x[:, 1:] - diff_x[:, :-1])
        result -= (dt / dz) * (diff_z[1:, :] - diff_z[:-1, :])

    result[solid] = 0.0
    outflow_rate_per_span = (
        float(np.sum(np.maximum(flux_x[:, -1], 0.0))) * dz
        - float(np.sum(np.minimum(flux_x[:, 0], 0.0))) * dz
    )
    return result, outflow_rate_per_span * dt


def solve_fast_obstacle_transport(
    config: FastObstacleTransportConfig,
    *,
    receptors: Sequence[FastObstacleReceptor] = (),
) -> FastObstacleTransportResult:
    """Run the Tier-1 potential-flow obstacle screen."""
    if not isinstance(config, FastObstacleTransportConfig):
        raise TypeError("config must be FastObstacleTransportConfig")
    checked_receptors = tuple(receptors)
    if any(not isinstance(item, FastObstacleReceptor) for item in checked_receptors):
        raise TypeError("receptors must contain FastObstacleReceptor values")
    labels = tuple(item.label for item in checked_receptors)
    if len(set(labels)) != len(labels):
        raise ValueError("receptor labels must be unique")

    wall_started = time.perf_counter()
    flow_started = time.perf_counter()
    (
        x, z, psi, u_face, w_face, u_cell, w_cell, pressure_coefficient,
        solid,
    ) = _solve_streamfunction(config)
    flow_runtime = time.perf_counter() - flow_started
    dx, dz = config.length_m / config.nx, config.height_m / config.nz
    cell_area = dx * dz

    receptor_indices: list[tuple[int, int]] = []
    for receptor in checked_receptors:
        if not (0.0 <= receptor.x_m <= config.length_m):
            raise ValueError(f"receptor {receptor.label!r} x lies outside the domain")
        if not (0.0 <= receptor.z_m <= config.height_m):
            raise ValueError(f"receptor {receptor.label!r} z lies outside the domain")
        i = int(np.clip(np.searchsorted(x, receptor.x_m), 0, config.nx - 1))
        j = int(np.clip(np.searchsorted(z, receptor.z_m), 0, config.nz - 1))
        if solid[j, i]:
            raise ValueError(f"receptor {receptor.label!r} lies inside the obstacle")
        receptor_indices.append((j, i))

    max_advective_number = config.time_step_s * (
        float(np.max(np.abs(u_face))) / dx
        + float(np.max(np.abs(w_face))) / dz
    )
    diffusion_number = 2.0 * config.diffusivity_m2_s * config.time_step_s * (
        1.0 / dx**2 + 1.0 / dz**2
    )
    if max_advective_number + diffusion_number > 1.0 + 1.0e-12:
        raise ValueError(
            "fast obstacle explicit stability limit exceeded: "
            f"advection={max_advective_number:.6g}, diffusion={diffusion_number:.6g}"
        )

    weights = _source_weights(config, x, z, solid)
    concentration = np.zeros((config.nz, config.nx), dtype=float)
    steps = int(math.ceil(config.duration_s / config.time_step_s))
    dt = config.duration_s / steps
    times = [0.0]
    injected = [0.0]
    domain_mass = [0.0]
    outflow = [0.0]
    residual = [0.0]
    traces = {item.label: [0.0] for item in checked_receptors}
    stored_fields = [concentration.copy()] if config.store_fields else None
    cumulative_injected = 0.0
    cumulative_outflow = 0.0

    transport_started = time.perf_counter()
    for step in range(steps):
        start = step * dt
        end = (step + 1) * dt
        source_mass_per_span = (
            config.source_schedule.mass_between(
                start, min(end, config.source_schedule.duration_s),
            )
            if start < config.source_schedule.duration_s
            else 0.0
        )
        if source_mass_per_span > 0.0:
            concentration += source_mass_per_span * weights / cell_area
        cumulative_injected += source_mass_per_span
        concentration, step_outflow = _transport_step(
            concentration,
            u_face=u_face,
            w_face=w_face,
            solid=solid,
            diffusivity=config.diffusivity_m2_s,
            dt=dt,
            dx=dx,
            dz=dz,
        )
        minimum = float(np.min(concentration))
        if minimum < -1.0e-12:
            raise RuntimeError(
                f"fast obstacle concentration became negative: {minimum:.6g}"
            )
        concentration = np.maximum(concentration, 0.0)
        cumulative_outflow += step_outflow
        mass = float(np.sum(concentration) * cell_area)
        now = end
        times.append(now)
        injected.append(cumulative_injected)
        domain_mass.append(mass)
        outflow.append(cumulative_outflow)
        residual.append(mass + cumulative_outflow - cumulative_injected)
        for receptor, (j, i) in zip(checked_receptors, receptor_indices):
            traces[receptor.label].append(float(concentration[j, i]))
        if stored_fields is not None:
            stored_fields.append(concentration.copy())
    transport_runtime = time.perf_counter() - transport_started

    divergence = (
        (u_face[:, 1:] - u_face[:, :-1]) / dx
        + (w_face[1:, :] - w_face[:-1, :]) / dz
    )
    divergence[solid] = 0.0
    total_runtime = time.perf_counter() - wall_started
    result_concentration = concentration.copy()
    result_concentration[solid] = np.nan
    return FastObstacleTransportResult(
        schema=FAST_OBSTACLE_TRANSPORT_SCHEMA,
        x_m=x,
        z_m=z,
        streamfunction_m2_s=psi,
        u_face_m_s=u_face,
        w_face_m_s=w_face,
        velocity_u_m_s=u_cell,
        velocity_w_m_s=w_cell,
        pressure_coefficient=pressure_coefficient,
        solid_mask=solid,
        time_s=np.asarray(times),
        concentration_kg_m3=result_concentration,
        concentration_fields_kg_m3=(
            None if stored_fields is None else np.asarray(stored_fields)
        ),
        receptor_traces_kg_m3={
            label: np.asarray(values) for label, values in traces.items()
        },
        injected_mass_per_span_kg_m=np.asarray(injected),
        domain_mass_per_span_kg_m=np.asarray(domain_mass),
        outflow_mass_per_span_kg_m=np.asarray(outflow),
        mass_residual_per_span_kg_m=np.asarray(residual),
        diagnostics={
            "fidelity_tier": "tier_1_fast_potential_flow",
            "applicability": "conditional",
            "quantitative_prediction_allowed": False,
            "nonlocal_scalar_routing": False,
            "flow_model": "2-D incompressible irrotational potential flow",
            "scalar_model": "local conservative finite-volume advection-diffusion",
            "solid_velocity_boundary": "no penetration; tangential slip",
            "solid_scalar_boundary": "zero normal advective and diffusive flux",
            "cell_size_m": [dx, dz],
            "time_step_s": dt,
            "advective_number": max_advective_number,
            "diffusion_number": diffusion_number,
            "maximum_absolute_divergence_s_1": float(np.max(np.abs(divergence))),
            "maximum_mass_residual_per_span_kg_m": float(
                np.max(np.abs(residual))
            ),
            "flow_solve_runtime_s": flow_runtime,
            "transport_runtime_s": transport_runtime,
            "total_runtime_s": total_runtime,
            "warnings": [
                "potential flow has no no-slip separation or recirculating wake",
                "wind-plane model has no finite-width lateral bypass",
                "temperature/density feedback is not active in Tier 1",
                "not validated for quantitative LH2 facility prediction",
            ],
        },
        ambient_temperature_k=config.ambient_temperature_k,
        ambient_pressure_pa=config.ambient_pressure_pa,
    )


__all__ = [
    "FAST_OBSTACLE_TRANSPORT_SCHEMA",
    "FastObstacleReceptor",
    "FastObstacleTransportConfig",
    "FastObstacleTransportResult",
    "solve_fast_obstacle_transport",
]
