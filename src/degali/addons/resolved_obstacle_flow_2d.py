"""Resolved two-dimensional viscous obstacle flow and hydrogen transport.

This is the second fidelity tier of the obstacle-physics rebuild.  The flow
solver is an isothermal, low-Mach D2Q9 BGK lattice-Boltzmann discretisation of
the incompressible Navier--Stokes limit.  A ground-mounted Cartesian obstacle
uses link bounce-back, so velocity is not prescribed around its roof or lee
side: stagnation, separation and a viscous recirculating wake emerge from the
momentum calculation.  Hydrogen then advances on the same grid with strictly
local conservative finite-volume face fluxes and zero scalar flux through
the solid.

The implementation is deliberately a validation-stage model.  Its effective
kinematic viscosity must be declared, the flow is two-dimensional and
isothermal, and hydrogen does not yet feed density or buoyancy back into the
momentum equation.  Consequently the result is not authorized for
quantitative facility decisions until matched obstacle validation is added.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import time
from typing import Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .fast_obstacle_transport import (
    HYDROGEN_MOLAR_MASS_KG_MOL,
    R_GAS_J_MOL_K,
    FastObstacleReceptor,
    _transport_step,
)
from .semi_fv_obstacle import RectangularObstacle2D, SourceRateSchedule


RESOLVED_OBSTACLE_FLOW_2D_SCHEMA = "degali.resolved-obstacle-flow-2d.v1"

_CX = np.asarray((0, 1, 0, -1, 0, 1, -1, -1, 1), dtype=int)
_CZ = np.asarray((0, 0, 1, 0, -1, 1, 1, -1, -1), dtype=int)
_WEIGHTS = np.asarray((4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9,
                       1 / 36, 1 / 36, 1 / 36, 1 / 36), dtype=float)
_OPPOSITE = np.asarray((0, 3, 4, 1, 2, 7, 8, 5, 6), dtype=int)
_CS2 = 1.0 / 3.0


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class ResolvedObstacleFlow2DConfig:
    """Physical and numerical inputs for the Tier-2 wind-plane model.

    ``kinematic_viscosity_m2_s`` is an effective molecular-plus-unresolved
    viscosity.  ``source_schedule`` is hydrogen mass rate per metre of span,
    so all two-dimensional mass ledgers are in kg/m.
    """

    length_m: float
    height_m: float
    nx: int
    nz: int
    wind_speed_m_s: float
    kinematic_viscosity_m2_s: float
    flow_steps: int
    lattice_inlet_speed: float
    duration_s: float
    scalar_time_step_s: float
    scalar_diffusivity_m2_s: float
    source_schedule: SourceRateSchedule
    source_x_m: float
    source_height_m: float
    source_sigma_m: float
    obstacle: RectangularObstacle2D
    inlet_profile_z_m: tuple[float, ...] | None = None
    inlet_profile_u_m_s: tuple[float, ...] | None = None
    flow_convergence_tolerance: float = 1.0e-7
    minimum_flow_steps: int = 500
    convergence_check_interval: int = 100
    ambient_temperature_k: float = 293.15
    ambient_pressure_pa: float = 101_325.0
    ambient_air_density_kg_m3: float = 1.204
    store_scalar_fields: bool = False

    def __post_init__(self) -> None:
        positive = (
            "length_m", "height_m", "wind_speed_m_s",
            "kinematic_viscosity_m2_s", "lattice_inlet_speed", "duration_s",
            "scalar_time_step_s", "source_sigma_m", "ambient_temperature_k",
            "ambient_pressure_pa", "ambient_air_density_kg_m3",
        )
        for name in positive:
            value = _finite(getattr(self, name), name)
            object.__setattr__(self, name, value)
            if value <= 0.0:
                raise ValueError(f"{name} must be positive")
        object.__setattr__(
            self, "scalar_diffusivity_m2_s",
            _finite(self.scalar_diffusivity_m2_s, "scalar_diffusivity_m2_s"),
        )
        if self.scalar_diffusivity_m2_s < 0.0:
            raise ValueError("scalar_diffusivity_m2_s must be non-negative")
        object.__setattr__(
            self, "flow_convergence_tolerance",
            _finite(self.flow_convergence_tolerance, "flow_convergence_tolerance"),
        )
        if self.flow_convergence_tolerance <= 0.0:
            raise ValueError("flow_convergence_tolerance must be positive")
        for name, minimum in (
            ("nx", 16), ("nz", 12), ("flow_steps", 1),
            ("minimum_flow_steps", 1), ("convergence_check_interval", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer of at least {minimum}")
        if self.minimum_flow_steps > self.flow_steps:
            raise ValueError("minimum_flow_steps cannot exceed flow_steps")
        dx = self.length_m / self.nx
        dz = self.height_m / self.nz
        if not math.isclose(dx, dz, rel_tol=1.0e-9, abs_tol=1.0e-12):
            raise ValueError("D2Q9 flow grid requires equal x and z cell sizes")
        if self.lattice_inlet_speed >= 0.1:
            raise ValueError("lattice_inlet_speed must remain below 0.1 for low-Mach flow")
        if not isinstance(self.source_schedule, SourceRateSchedule):
            raise TypeError("source_schedule must be a SourceRateSchedule")
        self.source_schedule.require_zero_endpoint("resolved obstacle source schedule")
        if self.source_schedule.duration_s > self.duration_s + 1.0e-12:
            raise ValueError("source schedule cannot outlast the scalar transport run")
        if not isinstance(self.obstacle, RectangularObstacle2D):
            raise TypeError("obstacle must be RectangularObstacle2D")
        if abs(self.obstacle.base_height_m) > 1.0e-12:
            raise ValueError("Tier-2 solver currently requires a ground-mounted obstacle")
        if self.obstacle.x_min_m <= 2.0 * dx or self.obstacle.x_max_m >= self.length_m - 3.0 * dx:
            raise ValueError("obstacle must leave resolved inlet and outlet fetch")
        if self.obstacle.z_max_m >= self.height_m - 2.0 * dz:
            raise ValueError("obstacle must leave a resolved flow path below the top boundary")
        for name in ("source_x_m", "source_height_m"):
            value = _finite(getattr(self, name), name)
            object.__setattr__(self, name, value)
        if not (0.0 <= self.source_x_m <= self.length_m):
            raise ValueError("source_x_m must lie inside the domain")
        if not (0.0 <= self.source_height_m <= self.height_m):
            raise ValueError("source_height_m must lie inside the domain")
        if (
            self.obstacle.x_min_m <= self.source_x_m <= self.obstacle.x_max_m
            and self.source_height_m <= self.obstacle.z_max_m
        ):
            raise ValueError("source may not lie inside the obstacle")
        if (self.inlet_profile_z_m is None) != (self.inlet_profile_u_m_s is None):
            raise ValueError(
                "inlet_profile_z_m and inlet_profile_u_m_s must be supplied together"
            )
        if self.inlet_profile_z_m is not None:
            if not isinstance(self.inlet_profile_z_m, tuple) or not isinstance(
                self.inlet_profile_u_m_s, tuple
            ):
                raise TypeError("inlet profile coordinates and velocities must be tuples")
            if len(self.inlet_profile_z_m) < 2 or len(self.inlet_profile_z_m) != len(
                self.inlet_profile_u_m_s
            ):
                raise ValueError(
                    "inlet profile requires equally sized tuples of at least two points"
                )
            profile_z = tuple(
                _finite(value, "inlet_profile_z_m") for value in self.inlet_profile_z_m
            )
            profile_u = tuple(
                _finite(value, "inlet_profile_u_m_s")
                for value in self.inlet_profile_u_m_s
            )
            if any(right <= left for left, right in zip(profile_z, profile_z[1:])):
                raise ValueError("inlet_profile_z_m must be strictly increasing")
            if profile_z[0] < 0.0 or any(value < 0.0 for value in profile_u):
                raise ValueError("inlet profile heights and velocities must be non-negative")
            object.__setattr__(self, "inlet_profile_z_m", profile_z)
            object.__setattr__(self, "inlet_profile_u_m_s", profile_u)
        if not isinstance(self.store_scalar_fields, bool):
            raise TypeError("store_scalar_fields must be boolean")


@dataclass(frozen=True)
class ResolvedObstacleFlow2DResult:
    """Resolved flow, scalar field, histories and validation diagnostics."""

    schema: str
    x_m: np.ndarray
    z_m: np.ndarray
    solid_mask: np.ndarray
    velocity_u_m_s: np.ndarray
    velocity_w_m_s: np.ndarray
    u_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    pressure_pa_gauge: np.ndarray
    pressure_coefficient: np.ndarray
    vorticity_s_1: np.ndarray
    time_s: np.ndarray
    concentration_kg_m3: np.ndarray
    maximum_concentration_kg_m3: np.ndarray
    concentration_fields_kg_m3: np.ndarray | None
    receptor_traces_kg_m3: dict[str, np.ndarray]
    injected_mass_per_span_kg_m: np.ndarray
    domain_mass_per_span_kg_m: np.ndarray
    outflow_mass_per_span_kg_m: np.ndarray
    mass_residual_per_span_kg_m: np.ndarray
    flow_residual_history: np.ndarray
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
    def maximum_concentration_vol_percent(self) -> np.ndarray:
        molar_density = self.ambient_pressure_pa / (
            R_GAS_J_MOL_K * self.ambient_temperature_k
        )
        result = (
            self.maximum_concentration_kg_m3
            / HYDROGEN_MOLAR_MASS_KG_MOL
            / molar_density
            * 100.0
        )
        return np.clip(result, 0.0, 100.0)

    @property
    def maximum_mass_residual_per_span_kg_m(self) -> float:
        return float(np.max(np.abs(self.mass_residual_per_span_kg_m)))


def _equilibrium(rho: np.ndarray, ux: np.ndarray, uz: np.ndarray) -> np.ndarray:
    cu = (
        _CX[:, None, None] * ux[None, :, :]
        + _CZ[:, None, None] * uz[None, :, :]
    )
    speed2 = ux**2 + uz**2
    return _WEIGHTS[:, None, None] * rho[None, :, :] * (
        1.0 + 3.0 * cu + 4.5 * cu**2 - 1.5 * speed2[None, :, :]
    )


def _moments(distribution: np.ndarray, solid: np.ndarray) -> tuple[np.ndarray, ...]:
    rho = np.sum(distribution, axis=0)
    rho = np.maximum(rho, 1.0e-12)
    ux = np.tensordot(_CX, distribution, axes=(0, 0)) / rho
    uz = np.tensordot(_CZ, distribution, axes=(0, 0)) / rho
    ux[solid] = 0.0
    uz[solid] = 0.0
    return rho, ux, uz


def _obstacle_mask(config: ResolvedObstacleFlow2DConfig) -> tuple[np.ndarray, ...]:
    dx = config.length_m / config.nx
    dz = config.height_m / config.nz
    x = (np.arange(config.nx) + 0.5) * dx
    z = (np.arange(config.nz) + 0.5) * dz
    obstacle = config.obstacle
    solid = (
        (x[None, :] >= obstacle.x_min_m)
        & (x[None, :] <= obstacle.x_max_m)
        & (z[:, None] >= obstacle.base_height_m)
        & (z[:, None] <= obstacle.z_max_m)
    )
    return x, z, solid


def _apply_outer_boundaries(
    distribution: np.ndarray,
    *,
    inlet_speed: np.ndarray,
    solid: np.ndarray,
) -> None:
    """Apply no-slip ground, slip top, Zou--He inlet and pressure outlet."""
    # Bottom no-slip bounce-back; assignments use copies to avoid aliasing.
    distribution[2, 0, :] = distribution[4, 0, :].copy()
    distribution[5, 0, :] = distribution[7, 0, :].copy()
    distribution[6, 0, :] = distribution[8, 0, :].copy()
    # Top specular reflection (zero normal velocity and zero tangential shear).
    distribution[4, -1, :] = distribution[2, -1, :].copy()
    distribution[7, -1, :] = distribution[6, -1, :].copy()
    distribution[8, -1, :] = distribution[5, -1, :].copy()

    fluid_inlet = ~solid[:, 0]
    f = distribution[:, :, 0]
    if inlet_speed.shape != (solid.shape[0],):
        raise ValueError("inlet_speed profile must have one value per z cell")
    rho = (
        f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])
    ) / (1.0 - inlet_speed)
    f[1, fluid_inlet] = (
        f[3, fluid_inlet]
        + (2.0 / 3.0) * rho[fluid_inlet] * inlet_speed[fluid_inlet]
    )
    f[5, fluid_inlet] = (
        f[7, fluid_inlet]
        + 0.5 * (f[4, fluid_inlet] - f[2, fluid_inlet])
        + (1.0 / 6.0) * rho[fluid_inlet] * inlet_speed[fluid_inlet]
    )
    f[8, fluid_inlet] = (
        f[6, fluid_inlet]
        + 0.5 * (f[2, fluid_inlet] - f[4, fluid_inlet])
        + (1.0 / 6.0) * rho[fluid_inlet] * inlet_speed[fluid_inlet]
    )

    # Unit-density pressure outlet with reconstructed streamwise velocity.
    fluid_outlet = ~solid[:, -1]
    f = distribution[:, :, -1]
    outlet_rho = 1.0
    ux = -1.0 + (
        f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])
    ) / outlet_rho
    f[3, fluid_outlet] = (
        f[1, fluid_outlet] - (2.0 / 3.0) * outlet_rho * ux[fluid_outlet]
    )
    f[7, fluid_outlet] = (
        f[5, fluid_outlet]
        + 0.5 * (f[2, fluid_outlet] - f[4, fluid_outlet])
        - (1.0 / 6.0) * outlet_rho * ux[fluid_outlet]
    )
    f[6, fluid_outlet] = (
        f[8, fluid_outlet]
        + 0.5 * (f[4, fluid_outlet] - f[2, fluid_outlet])
        - (1.0 / 6.0) * outlet_rho * ux[fluid_outlet]
    )


def _solve_flow(config: ResolvedObstacleFlow2DConfig) -> tuple[np.ndarray, ...]:
    x, z, solid = _obstacle_mask(config)
    dx = config.length_m / config.nx
    physical_dt = config.lattice_inlet_speed * dx / config.wind_speed_m_s
    if config.inlet_profile_z_m is None:
        inlet_u_physical = np.full(config.nz, config.wind_speed_m_s, dtype=float)
    else:
        assert config.inlet_profile_u_m_s is not None
        inlet_u_physical = np.interp(
            z,
            np.asarray(config.inlet_profile_z_m),
            np.asarray(config.inlet_profile_u_m_s),
        )
    inlet_u_lattice = (
        inlet_u_physical * config.lattice_inlet_speed / config.wind_speed_m_s
    )
    if float(np.max(inlet_u_lattice)) >= 0.1:
        raise ValueError(
            "mapped inlet profile exceeds the low-Mach lattice speed limit 0.1"
        )
    lattice_viscosity = config.kinematic_viscosity_m2_s * physical_dt / dx**2
    relaxation_time = 0.5 + 3.0 * lattice_viscosity
    if relaxation_time <= 0.505:
        raise ValueError(
            "D2Q9 relaxation time is too close to the stability limit; increase "
            "effective viscosity, lattice speed, or cell resolution"
        )
    rho = np.ones((config.nz, config.nx), dtype=float)
    ux = np.repeat(inlet_u_lattice[:, None], config.nx, axis=1)
    uz = np.zeros_like(rho)
    ux[solid] = 0.0
    distribution = _equilibrium(rho, ux, uz)
    residual_history: list[tuple[float, float]] = []
    previous_ux = ux.copy()
    previous_uz = uz.copy()
    actual_steps = config.flow_steps

    for step in range(1, config.flow_steps + 1):
        rho, ux, uz = _moments(distribution, solid)
        equilibrium = _equilibrium(rho, ux, uz)
        post_collision = distribution - (
            distribution - equilibrium
        ) / relaxation_time
        streamed = np.empty_like(distribution)
        for direction in range(9):
            streamed[direction] = np.roll(
                post_collision[direction],
                shift=(_CZ[direction], _CX[direction]),
                axis=(0, 1),
            )
        # Link-wise bounce-back at every fluid-to-solid link.
        for direction in range(1, 9):
            destination_solid = np.roll(
                solid,
                shift=(-_CZ[direction], -_CX[direction]),
                axis=(0, 1),
            )
            bounce = (~solid) & destination_solid
            streamed[_OPPOSITE[direction], bounce] = post_collision[direction, bounce]
        _apply_outer_boundaries(
            streamed, inlet_speed=inlet_u_lattice, solid=solid,
        )
        distribution = streamed

        if step % config.convergence_check_interval == 0 or step == config.flow_steps:
            if not np.all(np.isfinite(distribution)):
                raise RuntimeError(
                    "D2Q9 flow became non-finite; increase effective viscosity or "
                    "reduce the mapped lattice inlet speed"
                )
            _, current_ux, current_uz = _moments(distribution, solid)
            change = np.sqrt(
                (current_ux - previous_ux) ** 2
                + (current_uz - previous_uz) ** 2
            )
            reference = max(
                float(np.max(np.sqrt(current_ux**2 + current_uz**2))),
                float(np.max(np.abs(inlet_u_lattice))),
            )
            residual = float(np.max(change) / reference)
            residual_history.append((float(step), residual))
            previous_ux, previous_uz = current_ux.copy(), current_uz.copy()
            if step >= config.minimum_flow_steps and residual <= config.flow_convergence_tolerance:
                actual_steps = step
                break

    rho, ux, uz = _moments(distribution, solid)
    velocity_scale = dx / physical_dt
    u = ux * velocity_scale
    w = uz * velocity_scale
    u[solid] = 0.0
    w[solid] = 0.0

    # Isothermal LBM pressure; subtract the fluid outlet reference.
    pressure_scale = config.ambient_air_density_kg_m3 * velocity_scale**2
    outlet_reference = float(np.mean(rho[~solid[:, -1], -1]))
    pressure = _CS2 * (rho - outlet_reference) * pressure_scale
    dynamic_pressure = 0.5 * config.ambient_air_density_kg_m3 * config.wind_speed_m_s**2
    pressure_coefficient = pressure / dynamic_pressure
    pressure[solid] = np.nan
    pressure_coefficient[solid] = np.nan

    u_face = np.empty((config.nz, config.nx + 1), dtype=float)
    w_face = np.empty((config.nz + 1, config.nx), dtype=float)
    u_face[:, 0] = inlet_u_physical
    u_face[:, -1] = u[:, -1]
    u_face[:, 1:-1] = 0.5 * (u[:, :-1] + u[:, 1:])
    w_face[0, :] = 0.0
    w_face[-1, :] = 0.0
    w_face[1:-1, :] = 0.5 * (w[:-1, :] + w[1:, :])
    blocked_x = solid[:, :-1] | solid[:, 1:]
    u_face[:, 1:-1][blocked_x] = 0.0
    u_face[:, 0][solid[:, 0]] = 0.0
    u_face[:, -1][solid[:, -1]] = 0.0
    blocked_z = solid[:-1, :] | solid[1:, :]
    w_face[1:-1, :][blocked_z] = 0.0

    (
        u_face, w_face, divergence_before_projection,
        divergence_after_projection,
    ) = _mass_consistent_face_projection(u_face, w_face, solid, dx)

    # Preserve the D2Q9 momentum solution as the reported validation field.
    # The separately projected faces are a transport operator: rebuilding the
    # cell field from them can distort the pressure-outlet neighbourhood and
    # must not silently replace the flow quantity being validated.
    edge_order = 2 if min(config.nx, config.nz) >= 3 else 1
    du_dz = np.gradient(u, dx, axis=0, edge_order=edge_order)
    dw_dx = np.gradient(w, dx, axis=1, edge_order=edge_order)
    vorticity = dw_dx - du_dz
    vorticity[solid] = np.nan
    return (
        x, z, solid, u, w, u_face, w_face, pressure,
        pressure_coefficient, vorticity, np.asarray(residual_history),
        rho, relaxation_time, physical_dt, actual_steps,
        divergence_before_projection, divergence_after_projection,
    )


def _source_weights(
    config: ResolvedObstacleFlow2DConfig,
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


def _receptor_stencil(
    x: np.ndarray,
    z: np.ndarray,
    solid: np.ndarray,
    receptor: FastObstacleReceptor,
) -> tuple[tuple[int, int, float], ...]:
    """Return a bilinear, grid-continuous observation stencil in fluid cells."""
    if not (0.0 <= receptor.x_m <= float(x[-1] + 0.5 * (x[1] - x[0]))):
        raise ValueError(f"receptor {receptor.label!r} x lies outside the domain")
    if not (0.0 <= receptor.z_m <= float(z[-1] + 0.5 * (z[1] - z[0]))):
        raise ValueError(f"receptor {receptor.label!r} z lies outside the domain")

    def bracket(coordinates: np.ndarray, value: float) -> tuple[int, int, float]:
        if value <= coordinates[0]:
            return 0, 0, 0.0
        if value >= coordinates[-1]:
            last = len(coordinates) - 1
            return last, last, 0.0
        upper = int(np.searchsorted(coordinates, value, side="right"))
        lower = upper - 1
        fraction = float(
            (value - coordinates[lower])
            / (coordinates[upper] - coordinates[lower])
        )
        return lower, upper, fraction

    i0, i1, fx = bracket(x, receptor.x_m)
    j0, j1, fz = bracket(z, receptor.z_m)
    raw = (
        (j0, i0, (1.0 - fx) * (1.0 - fz)),
        (j0, i1, fx * (1.0 - fz)),
        (j1, i0, (1.0 - fx) * fz),
        (j1, i1, fx * fz),
    )
    combined: dict[tuple[int, int], float] = {}
    for j, i, weight in raw:
        if weight <= 1.0e-15:
            continue
        combined[(j, i)] = combined.get((j, i), 0.0) + weight
    if any(solid[j, i] for j, i in combined):
        raise ValueError(
            f"receptor {receptor.label!r} interpolation stencil intersects the obstacle"
        )
    return tuple((j, i, weight) for (j, i), weight in combined.items())


def _mass_consistent_face_projection(
    u_face: np.ndarray,
    w_face: np.ndarray,
    solid: np.ndarray,
    cell_size: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Project LBM-derived face velocity onto the local FV continuity space.

    LBM is weakly compressible and its cell moments are not exactly compatible
    with the separate finite-volume scalar stencil.  This Helmholtz projection
    keeps inlet, outlet and impermeable faces fixed, applies zero-normal
    correction at walls, and removes divergence from all non-outlet fluid
    cells to linear-solver precision.
    """
    u_projected = u_face.copy()
    w_projected = w_face.copy()
    nz, nx = solid.shape
    divergence_before = (
        (u_projected[:, 1:] - u_projected[:, :-1]) / cell_size
        + (w_projected[1:, :] - w_projected[:-1, :]) / cell_size
    )
    fixed = np.zeros_like(solid)
    fixed[:, -1] = ~solid[:, -1]
    unknown = (~solid) & (~fixed)
    positions = np.argwhere(unknown)
    lookup = -np.ones_like(solid, dtype=int)
    for index, (j, i) in enumerate(positions):
        lookup[j, i] = index
    matrix = sparse.lil_matrix((len(positions), len(positions)), dtype=float)
    rhs = np.zeros(len(positions), dtype=float)
    coefficient = 1.0 / cell_size**2
    for row, (j, i) in enumerate(positions):
        diagonal = 0.0
        for dj, di in ((0, -1), (0, 1), (-1, 0), (1, 0)):
            nj, ni = j + dj, i + di
            if nj < 0 or nj >= nz or ni < 0 or ni >= nx or solid[nj, ni]:
                continue
            diagonal += coefficient
            neighbour = lookup[nj, ni]
            if neighbour >= 0:
                matrix[row, neighbour] = -coefficient
        matrix[row, row] = diagonal
        rhs[row] = -divergence_before[j, i]
    potential = np.zeros_like(solid, dtype=float)
    if len(positions):
        potential[unknown] = spsolve(matrix.tocsr(), rhs)

    fluid_x = (~solid[:, :-1]) & (~solid[:, 1:])
    correction_x = (potential[:, 1:] - potential[:, :-1]) / cell_size
    u_projected[:, 1:-1][fluid_x] -= correction_x[fluid_x]
    fluid_z = (~solid[:-1, :]) & (~solid[1:, :])
    correction_z = (potential[1:, :] - potential[:-1, :]) / cell_size
    w_projected[1:-1, :][fluid_z] -= correction_z[fluid_z]

    divergence_after = (
        (u_projected[:, 1:] - u_projected[:, :-1]) / cell_size
        + (w_projected[1:, :] - w_projected[:-1, :]) / cell_size
    )
    core = (~solid).copy()
    core[:, -1] = False
    return (
        u_projected,
        w_projected,
        float(np.max(np.abs(divergence_before[core]))),
        float(np.max(np.abs(divergence_after[core]))),
    )


def solve_resolved_obstacle_flow_2d(
    config: ResolvedObstacleFlow2DConfig,
    *,
    receptors: Sequence[FastObstacleReceptor] = (),
) -> ResolvedObstacleFlow2DResult:
    """Solve viscous obstacle flow, then conservative hydrogen transport."""
    if not isinstance(config, ResolvedObstacleFlow2DConfig):
        raise TypeError("config must be ResolvedObstacleFlow2DConfig")
    checked_receptors = tuple(receptors)
    if any(not isinstance(item, FastObstacleReceptor) for item in checked_receptors):
        raise TypeError("receptors must contain FastObstacleReceptor values")
    labels = tuple(item.label for item in checked_receptors)
    if len(set(labels)) != len(labels):
        raise ValueError("receptor labels must be unique")

    wall_started = time.perf_counter()
    flow_started = time.perf_counter()
    (
        x, z, solid, u, w, u_face, w_face, pressure, pressure_coefficient,
        vorticity, residual_history, lattice_density, relaxation_time,
        flow_dt, actual_flow_steps,
        divergence_before_projection, divergence_after_projection,
    ) = _solve_flow(config)
    flow_runtime = time.perf_counter() - flow_started
    dx = config.length_m / config.nx
    cell_area = dx**2

    receptor_stencils: list[tuple[tuple[int, int, float], ...]] = []
    for receptor in checked_receptors:
        receptor_stencils.append(_receptor_stencil(x, z, solid, receptor))

    advective_number = config.scalar_time_step_s * (
        float(np.max(np.abs(u_face))) / dx
        + float(np.max(np.abs(w_face))) / dx
    )
    diffusion_number = (
        4.0 * config.scalar_diffusivity_m2_s
        * config.scalar_time_step_s / dx**2
    )
    if advective_number + diffusion_number > 1.0 + 1.0e-12:
        raise ValueError(
            "resolved obstacle scalar stability limit exceeded: "
            f"advection={advective_number:.6g}, diffusion={diffusion_number:.6g}"
        )

    weights = _source_weights(config, x, z, solid)
    concentration = np.zeros((config.nz, config.nx), dtype=float)
    maximum_concentration = concentration.copy()
    steps = int(math.ceil(config.duration_s / config.scalar_time_step_s))
    scalar_dt = config.duration_s / steps
    times = [0.0]
    injected = [0.0]
    domain_mass = [0.0]
    outflow = [0.0]
    residual = [0.0]
    traces = {item.label: [0.0] for item in checked_receptors}
    fields = [concentration.copy()] if config.store_scalar_fields else None
    cumulative_injected = 0.0
    cumulative_outflow = 0.0

    transport_started = time.perf_counter()
    for step in range(steps):
        start = step * scalar_dt
        end = (step + 1) * scalar_dt
        source_mass = (
            config.source_schedule.mass_between(
                start, min(end, config.source_schedule.duration_s),
            )
            if start < config.source_schedule.duration_s else 0.0
        )
        if source_mass > 0.0:
            concentration += source_mass * weights / cell_area
        cumulative_injected += source_mass
        concentration, step_outflow = _transport_step(
            concentration,
            u_face=u_face,
            w_face=w_face,
            solid=solid,
            diffusivity=config.scalar_diffusivity_m2_s,
            dt=scalar_dt,
            dx=dx,
            dz=dx,
        )
        minimum = float(np.min(concentration))
        if minimum < -1.0e-12:
            raise RuntimeError(
                f"resolved obstacle concentration became negative: {minimum:.6g}"
            )
        concentration = np.maximum(concentration, 0.0)
        np.maximum(maximum_concentration, concentration, out=maximum_concentration)
        cumulative_outflow += step_outflow
        mass = float(np.sum(concentration) * cell_area)
        times.append(end)
        injected.append(cumulative_injected)
        domain_mass.append(mass)
        outflow.append(cumulative_outflow)
        residual.append(mass + cumulative_outflow - cumulative_injected)
        for receptor, stencil in zip(checked_receptors, receptor_stencils):
            traces[receptor.label].append(float(sum(
                weight * concentration[j, i] for j, i, weight in stencil
            )))
        if fields is not None:
            fields.append(concentration.copy())
    transport_runtime = time.perf_counter() - transport_started

    divergence = (
        (u_face[:, 1:] - u_face[:, :-1]) / dx
        + (w_face[1:, :] - w_face[:-1, :]) / dx
    )
    divergence[solid] = 0.0
    obstacle = config.obstacle
    lee_region = (
        (x[None, :] > obstacle.x_max_m)
        & (x[None, :] <= obstacle.x_max_m + 4.0 * obstacle.height_m)
        & (z[:, None] <= 1.5 * obstacle.height_m)
        & ~solid
    )
    recirculating = lee_region & (u < 0.0)
    total_runtime = time.perf_counter() - wall_started
    result_concentration = concentration.copy()
    result_concentration[solid] = np.nan
    result_maximum_concentration = maximum_concentration.copy()
    result_maximum_concentration[solid] = np.nan
    density_variation = float(
        (np.max(lattice_density[~solid]) - np.min(lattice_density[~solid]))
        / np.mean(lattice_density[~solid])
    )
    return ResolvedObstacleFlow2DResult(
        schema=RESOLVED_OBSTACLE_FLOW_2D_SCHEMA,
        x_m=x,
        z_m=z,
        solid_mask=solid,
        velocity_u_m_s=u,
        velocity_w_m_s=w,
        u_face_m_s=u_face,
        w_face_m_s=w_face,
        pressure_pa_gauge=pressure,
        pressure_coefficient=pressure_coefficient,
        vorticity_s_1=vorticity,
        time_s=np.asarray(times),
        concentration_kg_m3=result_concentration,
        maximum_concentration_kg_m3=result_maximum_concentration,
        concentration_fields_kg_m3=(None if fields is None else np.asarray(fields)),
        receptor_traces_kg_m3={
            label: np.asarray(values) for label, values in traces.items()
        },
        injected_mass_per_span_kg_m=np.asarray(injected),
        domain_mass_per_span_kg_m=np.asarray(domain_mass),
        outflow_mass_per_span_kg_m=np.asarray(outflow),
        mass_residual_per_span_kg_m=np.asarray(residual),
        flow_residual_history=residual_history,
        diagnostics={
            "fidelity_tier": "tier_2_resolved_viscous_2d",
            "applicability": "conditional",
            "quantitative_prediction_allowed": False,
            "nonlocal_scalar_routing": False,
            "flow_model": "isothermal low-Mach D2Q9 BGK lattice-Boltzmann",
            "governing_limit": "two-dimensional incompressible Navier-Stokes",
            "pressure_model": "isothermal lattice equation of state p=rho*c_s^2",
            "scalar_model": "local conservative finite-volume advection-diffusion",
            "reported_cell_velocity": "unprojected D2Q9 momentum moments",
            "scalar_face_velocity": "mass-consistent Helmholtz projection",
            "receptor_operator": "bilinear interpolation over fluid cells only",
            "solid_velocity_boundary": "link bounce-back no-slip",
            "solid_scalar_boundary": "zero normal advective and diffusive flux",
            "ground_boundary": "no-slip bounce-back",
            "top_boundary": "impermeable free-slip",
            "inlet_boundary": (
                "Zou-He measured/interpolated velocity profile"
                if config.inlet_profile_z_m is not None
                else "Zou-He prescribed uniform velocity"
            ),
            "outlet_boundary": "Zou-He unit-density pressure",
            "cell_size_m": dx,
            "flow_physical_time_step_s": flow_dt,
            "scalar_time_step_s": scalar_dt,
            "flow_steps_requested": config.flow_steps,
            "flow_steps_completed": actual_flow_steps,
            "flow_converged": bool(
                residual_history.size
                and residual_history[-1, 1] <= config.flow_convergence_tolerance
            ),
            "final_flow_residual": (
                None if not residual_history.size else float(residual_history[-1, 1])
            ),
            "lattice_relaxation_time": relaxation_time,
            "maximum_inlet_speed_m_s": float(np.max(u_face[:, 0])),
            "minimum_inlet_speed_m_s": float(np.min(u_face[:, 0])),
            "maximum_lattice_mach_number": (
                float(np.max(np.abs(u_face[:, 0])))
                * config.lattice_inlet_speed
                / config.wind_speed_m_s
                / math.sqrt(_CS2)
            ),
            # Backward-compatible alias: this is the maximum when a profile is used.
            "lattice_mach_number": (
                float(np.max(np.abs(u_face[:, 0])))
                * config.lattice_inlet_speed
                / config.wind_speed_m_s
                / math.sqrt(_CS2)
            ),
            "maximum_relative_lattice_density_variation": density_variation,
            "reynolds_number_obstacle_height": (
                config.wind_speed_m_s * obstacle.height_m
                / config.kinematic_viscosity_m2_s
            ),
            "effective_kinematic_viscosity_m2_s": config.kinematic_viscosity_m2_s,
            "advective_number": advective_number,
            "diffusion_number": diffusion_number,
            "maximum_absolute_divergence_s_1": divergence_after_projection,
            "maximum_absolute_divergence_including_pressure_outlet_s_1": float(
                np.max(np.abs(divergence))
            ),
            "maximum_absolute_divergence_before_projection_s_1": (
                divergence_before_projection
            ),
            "maximum_absolute_core_divergence_s_1": divergence_after_projection,
            "recirculating_lee_cell_count": int(np.count_nonzero(recirculating)),
            "minimum_lee_streamwise_velocity_m_s": (
                None if not np.any(lee_region) else float(np.min(u[lee_region]))
            ),
            "maximum_mass_residual_per_span_kg_m": float(np.max(np.abs(residual))),
            "flow_solve_runtime_s": flow_runtime,
            "transport_runtime_s": transport_runtime,
            "total_runtime_s": total_runtime,
            "warnings": [
                "two-dimensional wind plane cannot resolve finite-width lateral bypass",
                "effective viscosity is prescribed and is not a calibrated turbulence closure",
                "isothermal flow omits cryogenic heat transfer and density/buoyancy feedback",
                "quantitative facility use is withheld pending matched obstacle validation",
            ],
        },
        ambient_temperature_k=config.ambient_temperature_k,
        ambient_pressure_pa=config.ambient_pressure_pa,
    )


__all__ = [
    "RESOLVED_OBSTACLE_FLOW_2D_SCHEMA",
    "ResolvedObstacleFlow2DConfig",
    "ResolvedObstacleFlow2DResult",
    "solve_resolved_obstacle_flow_2d",
]
