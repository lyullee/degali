"""Finite-width three-dimensional potential flow and local H2 transport.

This module is a fast geometry-resolved baseline, not a wake model.  A
Helmholtz projection diverts a prescribed stream through the fluid faces
around one or more finite axis-aligned cuboids.  The resulting face field is
mass-consistent, includes both roof and lateral bypass, and is then used by a
strictly local conservative finite-volume scalar operator.  Solid faces carry
neither advective nor diffusive hydrogen flux.

The potential-flow closure is inviscid and slip-walled.  It cannot generate
separation, recirculation, turbulent mixing, heat transfer, or buoyancy
feedback, so quantitative facility use remains prohibited.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import time
from typing import Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg
from scipy.special import erf

from .fast_obstacle_transport import (
    HYDROGEN_MOLAR_MASS_KG_MOL,
    R_GAS_J_MOL_K,
)
from .semi_fv_obstacle import SourceRateSchedule
from .site_geometry import (
    AxisAlignedCuboid,
    CircularArcWall,
    CrosswindCylinder,
    SolidObstacle3D,
    circular_arc_surface_face_masks,
    solid_obstacle_bounds_m,
    solid_obstacle_contains_point,
    solid_obstacle_voxel_mask,
)


FAST_OBSTACLE_TRANSPORT_3D_SCHEMA = "degali.fast-obstacle-transport-3d.v2"


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class FastObstacleReceptor3D:
    label: str
    position_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("receptor label must be non-empty")
        if not isinstance(self.position_m, tuple) or len(self.position_m) != 3:
            raise TypeError("receptor position_m must be a three-value tuple")
        values = tuple(
            _finite(value, f"receptor position_m[{index}]")
            for index, value in enumerate(self.position_m)
        )
        if any(value < 0.0 for value in values):
            raise ValueError("receptor coordinates must be non-negative")
        object.__setattr__(self, "position_m", values)


@dataclass(frozen=True)
class FastObstacleTransport3DConfig:
    domain_m: tuple[float, float, float]
    cells: tuple[int, int, int]
    duration_s: float
    time_step_s: float
    wind_speed_m_s: float
    diffusivity_m2_s: float
    source_schedule: SourceRateSchedule
    source_position_m: tuple[float, float, float]
    source_sigma_m: float
    source_sigma_xyz_m: tuple[float, float, float] | None = None
    source_kernel: str = "gaussian"
    obstacles: tuple[SolidObstacle3D, ...] = ()
    open_boundary_faces: tuple[str, ...] = (
        "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
    )
    inlet_profile_z_m: tuple[float, ...] | None = None
    inlet_profile_u_m_s: tuple[float, ...] | None = None
    ambient_temperature_k: float = 293.15
    ambient_pressure_pa: float = 101_325.0
    store_fields: bool = False
    projection_relative_tolerance: float = 1.0e-10
    projection_maximum_iterations: int = 20_000

    def __post_init__(self) -> None:
        if not isinstance(self.domain_m, tuple) or len(self.domain_m) != 3:
            raise TypeError("domain_m must be a three-value tuple")
        domain = tuple(
            _finite(value, f"domain_m[{index}]")
            for index, value in enumerate(self.domain_m)
        )
        if any(value <= 0.0 for value in domain):
            raise ValueError("domain dimensions must be positive")
        object.__setattr__(self, "domain_m", domain)
        if (
            not isinstance(self.cells, tuple) or len(self.cells) != 3
            or any(isinstance(value, bool) or not isinstance(value, int) or value < 8
                   for value in self.cells)
        ):
            raise ValueError("cells must contain three integers of at least eight")
        for name in (
            "duration_s", "time_step_s", "wind_speed_m_s", "source_sigma_m",
            "ambient_temperature_k", "ambient_pressure_pa",
            "projection_relative_tolerance",
        ):
            value = _finite(getattr(self, name), name)
            object.__setattr__(self, name, value)
            if value <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.source_sigma_xyz_m is not None:
            if not isinstance(self.source_sigma_xyz_m, tuple) or len(self.source_sigma_xyz_m) != 3:
                raise TypeError("source_sigma_xyz_m must be a three-value tuple or None")
            sigma_xyz = tuple(
                _finite(value, f"source_sigma_xyz_m[{index}]")
                for index, value in enumerate(self.source_sigma_xyz_m)
            )
            if any(value <= 0.0 for value in sigma_xyz):
                raise ValueError("source_sigma_xyz_m values must be positive")
            object.__setattr__(self, "source_sigma_xyz_m", sigma_xyz)
        if self.source_kernel not in {
            "gaussian", "variance_matched_box", "cosine_tapered_box",
        }:
            raise ValueError(
                "source_kernel must be 'gaussian', 'variance_matched_box' "
                "or 'cosine_tapered_box'"
            )
        diffusivity = _finite(self.diffusivity_m2_s, "diffusivity_m2_s")
        if diffusivity < 0.0:
            raise ValueError("diffusivity_m2_s must be non-negative")
        object.__setattr__(self, "diffusivity_m2_s", diffusivity)
        if not isinstance(self.source_schedule, SourceRateSchedule):
            raise TypeError("source_schedule must be a SourceRateSchedule")
        self.source_schedule.require_zero_endpoint("fast 3-D obstacle source schedule")
        if self.source_schedule.duration_s > self.duration_s + 1.0e-12:
            raise ValueError("source schedule cannot outlast the transport run")
        if not isinstance(self.source_position_m, tuple) or len(self.source_position_m) != 3:
            raise TypeError("source_position_m must be a three-value tuple")
        source = tuple(
            _finite(value, f"source_position_m[{index}]")
            for index, value in enumerate(self.source_position_m)
        )
        if any(value < 0.0 or value > bound for value, bound in zip(source, domain)):
            raise ValueError("source_position_m must lie inside the domain")
        object.__setattr__(self, "source_position_m", source)
        if not isinstance(self.obstacles, tuple):
            raise TypeError("obstacles must be a tuple")
        if any(
            not isinstance(item, (AxisAlignedCuboid, CircularArcWall, CrosswindCylinder))
            for item in self.obstacles
        ):
            raise TypeError(
                "obstacles must contain AxisAlignedCuboid or CircularArcWall values"
            )
        for obstacle in self.obstacles:
            bounds = solid_obstacle_bounds_m(obstacle)
            if bounds[4] < 0.0:
                raise ValueError("obstacle may not extend below the domain")
            if not (
                0.0 < bounds[0] < bounds[1] < domain[0]
                and 0.0 < bounds[2] < bounds[3] < domain[1]
                and 0.0 <= bounds[4] < bounds[5] < domain[2]
            ):
                raise ValueError("obstacles must lie inside the domain with x/y fetch")
            if solid_obstacle_contains_point(obstacle, source):
                raise ValueError("source may not lie inside an obstacle")
        valid_boundary_faces = {
            "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
        }
        if (
            not isinstance(self.open_boundary_faces, tuple)
            or len(set(self.open_boundary_faces)) != len(self.open_boundary_faces)
            or any(face not in valid_boundary_faces for face in self.open_boundary_faces)
        ):
            raise ValueError(
                "open_boundary_faces must be a tuple of unique Cartesian face names"
            )
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
        if (
            isinstance(self.projection_maximum_iterations, bool)
            or not isinstance(self.projection_maximum_iterations, int)
            or self.projection_maximum_iterations < 1
        ):
            raise ValueError("projection_maximum_iterations must be a positive integer")
        if not isinstance(self.store_fields, bool):
            raise TypeError("store_fields must be boolean")


@dataclass(frozen=True)
class FastObstacleTransport3DResult:
    schema: str
    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    solid_mask: np.ndarray
    velocity_u_m_s: np.ndarray
    velocity_v_m_s: np.ndarray
    velocity_w_m_s: np.ndarray
    u_face_m_s: np.ndarray
    v_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    pressure_coefficient: np.ndarray
    time_s: np.ndarray
    concentration_kg_m3: np.ndarray
    maximum_concentration_kg_m3: np.ndarray
    concentration_fields_kg_m3: np.ndarray | None
    receptor_traces_kg_m3: dict[str, np.ndarray]
    injected_mass_kg: np.ndarray
    domain_mass_kg: np.ndarray
    outflow_mass_kg: np.ndarray
    mass_residual_kg: np.ndarray
    diagnostics: dict[str, object]
    ambient_temperature_k: float
    ambient_pressure_pa: float

    @property
    def concentration_vol_percent(self) -> np.ndarray:
        molar_density = self.ambient_pressure_pa / (
            R_GAS_J_MOL_K * self.ambient_temperature_k
        )
        return np.clip(
            self.concentration_kg_m3
            / HYDROGEN_MOLAR_MASS_KG_MOL
            / molar_density
            * 100.0,
            0.0, 100.0,
        )

    @property
    def maximum_concentration_vol_percent(self) -> np.ndarray:
        molar_density = self.ambient_pressure_pa / (
            R_GAS_J_MOL_K * self.ambient_temperature_k
        )
        return np.clip(
            self.maximum_concentration_kg_m3
            / HYDROGEN_MOLAR_MASS_KG_MOL
            / molar_density
            * 100.0,
            0.0, 100.0,
        )

    @property
    def maximum_mass_residual_kg(self) -> float:
        return float(np.max(np.abs(self.mass_residual_kg)))


def _grid_and_mask(config: FastObstacleTransport3DConfig) -> tuple[np.ndarray, ...]:
    nx, ny, nz = config.cells
    lx, ly, lz = config.domain_m
    dx, dy, dz = lx / nx, ly / ny, lz / nz
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    z = (np.arange(nz) + 0.5) * dz
    solid = np.zeros((nz, ny, nx), dtype=bool)
    for obstacle in config.obstacles:
        obstacle_mask = solid_obstacle_voxel_mask(obstacle, x, y, z)
        # A zero-thickness CircularArcWall is represented by the embedded
        # blocked-face operator below, not by solid cell centres. Every other
        # declared solid must occupy at least one cell or it would silently
        # disappear from this screening calculation.
        if not np.any(obstacle_mask) and not (
            isinstance(obstacle, CircularArcWall) and obstacle.thickness_m == 0.0
        ):
            raise ValueError(
                f"obstacle {obstacle.label!r} has no voxelized solid cells on "
                "the configured grid; refine/realign the grid or revise the "
                "obstacle geometry before solving"
            )
        solid |= obstacle_mask
    return x, y, z, solid, np.asarray((dx, dy, dz))


def _circular_arc_blocked_faces(
    config: FastObstacleTransport3DConfig,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Union exact mid-surface face crossings for circular-arc obstacles."""
    blocked_x = np.zeros((z.size, y.size, x.size + 1), dtype=bool)
    blocked_y = np.zeros((z.size, y.size + 1, x.size), dtype=bool)
    blocked_z = np.zeros((z.size + 1, y.size, x.size), dtype=bool)
    for obstacle in config.obstacles:
        if not isinstance(obstacle, CircularArcWall):
            continue
        faces = circular_arc_surface_face_masks(
            x, y, z,
            center_xy_m=obstacle.center_m,
            radius_m=obstacle.radius_m,
            z_min_m=obstacle.z_min_m,
            z_max_m=obstacle.z_max_m,
            arc_start_deg=obstacle.arc_start_deg,
            arc_end_deg=obstacle.arc_end_deg,
        )
        blocked_x |= faces[0]
        blocked_y |= faces[1]
    return blocked_x, blocked_y, blocked_z


def _mass_consistent_face_projection_3d(
    u_face: np.ndarray,
    v_face: np.ndarray,
    w_face: np.ndarray,
    solid: np.ndarray,
    spacing: np.ndarray,
    *,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    relative_tolerance: float = 1.0e-10,
    maximum_iterations: int = 20_000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]:
    """Project arbitrary 3-D face fluxes onto the scalar continuity space.

    External faces and impermeable solid faces remain fixed.  The last fluid
    x-cell is the zero-potential outlet reference; all other fluid-cell
    divergence is removed by a Jacobi-preconditioned conjugate-gradient solve.
    """
    u_projected = u_face.copy()
    v_projected = v_face.copy()
    w_projected = w_face.copy()
    nz, ny, nx = solid.shape
    dx, dy, dz = (float(value) for value in spacing)
    if blocked_faces is None:
        blocked_x, blocked_y, blocked_z = (
            np.zeros_like(u_projected, dtype=bool),
            np.zeros_like(v_projected, dtype=bool),
            np.zeros_like(w_projected, dtype=bool),
        )
    else:
        if not isinstance(blocked_faces, tuple) or len(blocked_faces) != 3:
            raise ValueError("blocked_faces must be a tuple of x, y, and z face masks")
        masks = tuple(np.asarray(mask, dtype=bool) for mask in blocked_faces)
        expected_shapes = (u_projected.shape, v_projected.shape, w_projected.shape)
        if any(mask.shape != shape for mask, shape in zip(masks, expected_shapes)):
            raise ValueError("blocked face masks must match their staggered face arrays")
        blocked_x, blocked_y, blocked_z = masks
    u_projected[blocked_x] = 0.0
    v_projected[blocked_y] = 0.0
    w_projected[blocked_z] = 0.0
    divergence_before = (
        (u_projected[:, :, 1:] - u_projected[:, :, :-1]) / dx
        + (v_projected[:, 1:, :] - v_projected[:, :-1, :]) / dy
        + (w_projected[1:, :, :] - w_projected[:-1, :, :]) / dz
    )
    fixed = np.zeros_like(solid)
    fixed[:, :, -1] = ~solid[:, :, -1]
    unknown = (~solid) & (~fixed)
    lookup = -np.ones_like(solid, dtype=np.int64)
    unknown_count = int(np.count_nonzero(unknown))
    lookup[unknown] = np.arange(unknown_count, dtype=np.int64)
    padded_solid = np.pad(solid, 1, mode="constant", constant_values=True)
    padded_lookup = np.pad(lookup, 1, mode="constant", constant_values=-1)
    neighbours = (
        (padded_solid[1:-1, 1:-1, :-2], padded_lookup[1:-1, 1:-1, :-2], blocked_x[:, :, :nx], 1.0 / dx**2),
        (padded_solid[1:-1, 1:-1, 2:], padded_lookup[1:-1, 1:-1, 2:], blocked_x[:, :, 1:nx + 1], 1.0 / dx**2),
        (padded_solid[1:-1, :-2, 1:-1], padded_lookup[1:-1, :-2, 1:-1], blocked_y[:, :ny, :], 1.0 / dy**2),
        (padded_solid[1:-1, 2:, 1:-1], padded_lookup[1:-1, 2:, 1:-1], blocked_y[:, 1:ny + 1, :], 1.0 / dy**2),
        (padded_solid[:-2, 1:-1, 1:-1], padded_lookup[:-2, 1:-1, 1:-1], blocked_z[:nz, :, :], 1.0 / dz**2),
        (padded_solid[2:, 1:-1, 1:-1], padded_lookup[2:, 1:-1, 1:-1], blocked_z[1:nz + 1, :, :], 1.0 / dz**2),
    )
    diagonal_grid = np.zeros_like(solid, dtype=float)
    off_diagonal_rows: list[np.ndarray] = []
    off_diagonal_columns: list[np.ndarray] = []
    off_diagonal_values: list[np.ndarray] = []
    for neighbour_solid, neighbour_lookup, face_blocked, coefficient in neighbours:
        connected = unknown & ~neighbour_solid & ~face_blocked
        diagonal_grid[connected] += coefficient
        connected_rows = lookup[connected]
        connected_neighbours = neighbour_lookup[connected]
        has_unknown_neighbour = connected_neighbours >= 0
        off_diagonal_rows.append(connected_rows[has_unknown_neighbour])
        off_diagonal_columns.append(connected_neighbours[has_unknown_neighbour])
        off_diagonal_values.append(
            np.full(np.count_nonzero(has_unknown_neighbour), -coefficient)
        )
    diagonal_indices = np.arange(unknown_count, dtype=np.int64)
    rows = np.concatenate((*off_diagonal_rows, diagonal_indices))
    columns = np.concatenate((*off_diagonal_columns, diagonal_indices))
    values = np.concatenate((*off_diagonal_values, diagonal_grid[unknown]))
    rhs = -divergence_before[unknown]
    matrix = sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(unknown_count, unknown_count),
    )
    potential = np.zeros_like(solid, dtype=float)
    iterations = 0
    if unknown_count:
        diagonal = matrix.diagonal()
        if np.any(diagonal <= 0.0):
            raise RuntimeError("3-D face projection contains an isolated fluid cell")
        preconditioner = LinearOperator(
            matrix.shape, matvec=lambda vector: vector / diagonal,
        )

        def count_iteration(_vector: np.ndarray) -> None:
            nonlocal iterations
            iterations += 1

        solution, info = cg(
            matrix, rhs,
            rtol=relative_tolerance,
            atol=1.0e-12,
            maxiter=maximum_iterations,
            M=preconditioner,
            callback=count_iteration,
        )
        if info != 0 or not np.all(np.isfinite(solution)):
            raise RuntimeError(f"3-D face projection failed to converge (cg info={info})")
        potential[unknown] = solution
    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:]) & ~blocked_x[:, :, 1:-1]
    correction_x = (potential[:, :, 1:] - potential[:, :, :-1]) / dx
    u_projected[:, :, 1:-1][fluid_x] -= correction_x[fluid_x]
    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :]) & ~blocked_y[:, 1:-1, :]
    correction_y = (potential[:, 1:, :] - potential[:, :-1, :]) / dy
    v_projected[:, 1:-1, :][fluid_y] -= correction_y[fluid_y]
    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :]) & ~blocked_z[1:-1, :, :]
    correction_z = (potential[1:, :, :] - potential[:-1, :, :]) / dz
    w_projected[1:-1, :, :][fluid_z] -= correction_z[fluid_z]
    divergence_after = (
        (u_projected[:, :, 1:] - u_projected[:, :, :-1]) / dx
        + (v_projected[:, 1:, :] - v_projected[:, :-1, :]) / dy
        + (w_projected[1:, :, :] - w_projected[:-1, :, :]) / dz
    )
    return (
        u_projected, v_projected, w_projected,
        float(np.max(np.abs(divergence_before[unknown]))) if np.any(unknown) else 0.0,
        float(np.max(np.abs(divergence_after[unknown]))) if np.any(unknown) else 0.0,
        iterations,
    )


def _potential_face_flow(
    config: FastObstacleTransport3DConfig,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    solid: np.ndarray,
    spacing: np.ndarray,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, ...]:
    nz, ny, nx = solid.shape
    dx, dy, dz = (float(value) for value in spacing)
    if blocked_faces is None:
        blocked_x, blocked_y, blocked_z = (
            np.zeros((nz, ny, nx + 1), dtype=bool),
            np.zeros((nz, ny + 1, nx), dtype=bool),
            np.zeros((nz + 1, ny, nx), dtype=bool),
        )
    else:
        blocked_x, blocked_y, blocked_z = blocked_faces
    if config.inlet_profile_z_m is None:
        inlet_u = np.full(nz, config.wind_speed_m_s, dtype=float)
    else:
        assert config.inlet_profile_u_m_s is not None
        inlet_u = np.interp(
            z,
            np.asarray(config.inlet_profile_z_m),
            np.asarray(config.inlet_profile_u_m_s),
        )
    u_face = np.broadcast_to(
        inlet_u[:, None, None], (nz, ny, nx + 1),
    ).copy()
    v_face = np.zeros((nz, ny + 1, nx), dtype=float)
    w_face = np.zeros((nz + 1, ny, nx), dtype=float)
    u_face[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    u_face[:, :, 0][solid[:, :, 0]] = 0.0
    u_face[:, :, -1][solid[:, :, -1]] = 0.0
    v_face[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    w_face[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    u_face[blocked_x] = 0.0
    v_face[blocked_y] = 0.0
    w_face[blocked_z] = 0.0

    divergence_before = (
        (u_face[:, :, 1:] - u_face[:, :, :-1]) / dx
        + (v_face[:, 1:, :] - v_face[:, :-1, :]) / dy
        + (w_face[1:, :, :] - w_face[:-1, :, :]) / dz
    )
    fixed = np.zeros_like(solid)
    fixed[:, :, -1] = ~solid[:, :, -1]
    unknown = (~solid) & (~fixed)
    positions = np.argwhere(unknown)
    lookup = -np.ones_like(solid, dtype=np.int64)
    lookup[unknown] = np.arange(len(positions), dtype=np.int64)
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    rhs = np.empty(len(positions), dtype=float)
    neighbours = (
        (0, 0, -1, 1.0 / dx**2, 0), (0, 0, 1, 1.0 / dx**2, 1),
        (0, -1, 0, 1.0 / dy**2, 2), (0, 1, 0, 1.0 / dy**2, 3),
        (-1, 0, 0, 1.0 / dz**2, 4), (1, 0, 0, 1.0 / dz**2, 5),
    )
    for row, (k, j, i) in enumerate(positions):
        diagonal = 0.0
        for dk, dj, di, coefficient, face_direction in neighbours:
            nk, nj, ni = k + dk, j + dj, i + di
            if (
                nk < 0 or nk >= nz or nj < 0 or nj >= ny
                or ni < 0 or ni >= nx or solid[nk, nj, ni]
            ):
                continue
            if face_direction == 0 and blocked_x[k, j, i]:
                continue
            if face_direction == 1 and blocked_x[k, j, i + 1]:
                continue
            if face_direction == 2 and blocked_y[k, j, i]:
                continue
            if face_direction == 3 and blocked_y[k, j + 1, i]:
                continue
            if face_direction == 4 and blocked_z[k, j, i]:
                continue
            if face_direction == 5 and blocked_z[k + 1, j, i]:
                continue
            diagonal += coefficient
            neighbour = int(lookup[nk, nj, ni])
            if neighbour >= 0:
                rows.append(row)
                columns.append(neighbour)
                values.append(-coefficient)
        rows.append(row)
        columns.append(row)
        values.append(diagonal)
        rhs[row] = -divergence_before[k, j, i]
    matrix = sparse.csr_matrix(
        (np.asarray(values), (np.asarray(rows), np.asarray(columns))),
        shape=(len(positions), len(positions)),
    )
    potential = np.zeros_like(solid, dtype=float)
    iterations = 0
    if len(positions):
        diagonal = matrix.diagonal()
        if np.any(diagonal <= 0.0):
            raise RuntimeError("3-D potential projection contains an isolated fluid cell")
        preconditioner = LinearOperator(
            matrix.shape, matvec=lambda vector: vector / diagonal,
        )

        def count_iteration(_vector: np.ndarray) -> None:
            nonlocal iterations
            iterations += 1

        solution, info = cg(
            matrix, rhs,
            rtol=config.projection_relative_tolerance,
            atol=1.0e-12,
            maxiter=config.projection_maximum_iterations,
            M=preconditioner,
            callback=count_iteration,
        )
        if info != 0 or not np.all(np.isfinite(solution)):
            raise RuntimeError(
                f"3-D potential projection failed to converge (cg info={info})"
            )
        potential[unknown] = solution

    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:]) & ~blocked_x[:, :, 1:-1]
    correction_x = (potential[:, :, 1:] - potential[:, :, :-1]) / dx
    u_face[:, :, 1:-1][fluid_x] -= correction_x[fluid_x]
    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :]) & ~blocked_y[:, 1:-1, :]
    correction_y = (potential[:, 1:, :] - potential[:, :-1, :]) / dy
    v_face[:, 1:-1, :][fluid_y] -= correction_y[fluid_y]
    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :]) & ~blocked_z[1:-1, :, :]
    correction_z = (potential[1:, :, :] - potential[:-1, :, :]) / dz
    w_face[1:-1, :, :][fluid_z] -= correction_z[fluid_z]

    divergence_after = (
        (u_face[:, :, 1:] - u_face[:, :, :-1]) / dx
        + (v_face[:, 1:, :] - v_face[:, :-1, :]) / dy
        + (w_face[1:, :, :] - w_face[:-1, :, :]) / dz
    )
    core = unknown
    u = 0.5 * (u_face[:, :, :-1] + u_face[:, :, 1:])
    v = 0.5 * (v_face[:, :-1, :] + v_face[:, 1:, :])
    w = 0.5 * (w_face[:-1, :, :] + w_face[1:, :, :])
    u[solid] = v[solid] = w[solid] = 0.0
    speed2 = u**2 + v**2 + w**2
    pressure_coefficient = 1.0 - speed2 / config.wind_speed_m_s**2
    pressure_coefficient[solid] = np.nan
    return (
        u_face, v_face, w_face, u, v, w, pressure_coefficient,
        float(np.max(np.abs(divergence_before[core]))) if np.any(core) else 0.0,
        float(np.max(np.abs(divergence_after[core]))) if np.any(core) else 0.0,
        iterations,
    )


def _source_weights(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, solid: np.ndarray,
    position: tuple[float, float, float], sigma: float,
) -> np.ndarray:
    return _source_weights_xyz(
        x, y, z, solid, position, (sigma, sigma, sigma),
    )


def _source_weights_xyz(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, solid: np.ndarray,
    position: tuple[float, float, float],
    sigma_xyz: tuple[float, float, float],
) -> np.ndarray:
    """Return cell-integrated, exactly normalised anisotropic source weights.

    The separable Gaussian is integrated over every Cartesian cell before the
    solid mask is applied.  Renormalisation therefore preserves the declared
    source mass exactly even when a narrow source is only marginally resolved;
    resolution adequacy remains a separate model-applicability question.
    """
    if len(sigma_xyz) != 3 or any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in sigma_xyz
    ):
        raise ValueError("source Gaussian sigma_xyz values must be positive and finite")

    def integrated_axis_weights(
        coordinates: np.ndarray, centre: float, sigma: float,
    ) -> np.ndarray:
        if coordinates.ndim != 1 or coordinates.size < 2:
            raise ValueError("source Gaussian requires at least two cells per axis")
        spacing = float(coordinates[1] - coordinates[0])
        if not np.allclose(
            np.diff(coordinates), spacing, rtol=1.0e-12, atol=1.0e-14,
        ):
            raise ValueError("source Gaussian requires uniform Cartesian coordinates")
        scale = math.sqrt(2.0) * sigma
        lower = coordinates - 0.5 * spacing
        upper = coordinates + 0.5 * spacing
        # The common sigma*sqrt(pi/2) factor cancels during normalisation.
        return erf((upper - centre) / scale) - erf((lower - centre) / scale)

    wx = integrated_axis_weights(x, position[0], float(sigma_xyz[0]))
    wy = integrated_axis_weights(y, position[1], float(sigma_xyz[1]))
    wz = integrated_axis_weights(z, position[2], float(sigma_xyz[2]))
    weights = wz[:, None, None] * wy[None, :, None] * wx[None, None, :]
    weights = np.maximum(weights, 0.0)
    weights[solid] = 0.0
    total = float(np.sum(weights))
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("source Gaussian has no fluid-cell support")
    return weights / total


def _source_weights_variance_matched_box(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, solid: np.ndarray,
    position: tuple[float, float, float],
    sigma_xyz: tuple[float, float, float],
) -> np.ndarray:
    """Return normalized compact-support weights with Gaussian-matched variance.

    A uniform box with half-width ``sqrt(3)*sigma`` has the same one-dimensional
    variance as a Gaussian with standard deviation ``sigma``.  Cell overlap,
    rather than cell-centre membership, keeps the source mass conservative.
    This is an explicit source-shape sensitivity, never the default kernel.
    """
    if len(sigma_xyz) != 3 or any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in sigma_xyz
    ):
        raise ValueError("source box sigma_xyz values must be positive and finite")

    def overlaps(coordinates: np.ndarray, centre: float, sigma: float) -> np.ndarray:
        if coordinates.ndim != 1 or coordinates.size < 2:
            raise ValueError("source box requires at least two cells per axis")
        spacing = float(coordinates[1] - coordinates[0])
        if not np.allclose(np.diff(coordinates), spacing, rtol=1.0e-12, atol=1.0e-14):
            raise ValueError("source box requires uniform Cartesian coordinates")
        half_width = math.sqrt(3.0) * float(sigma)
        lower = coordinates - 0.5 * spacing
        upper = coordinates + 0.5 * spacing
        return np.maximum(0.0, np.minimum(upper, centre + half_width)
                          - np.maximum(lower, centre - half_width))

    wx = overlaps(x, position[0], float(sigma_xyz[0]))
    wy = overlaps(y, position[1], float(sigma_xyz[1]))
    wz = overlaps(z, position[2], float(sigma_xyz[2]))
    weights = wz[:, None, None] * wy[None, :, None] * wx[None, None, :]
    weights[solid] = 0.0
    total = float(np.sum(weights))
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("source box has no fluid-cell support")
    return weights / total


def _source_weights_cosine_tapered_box(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, solid: np.ndarray,
    position: tuple[float, float, float],
    sigma_xyz: tuple[float, float, float],
) -> np.ndarray:
    """Return a compact, variance-scaled source with a smooth cosine edge.

    The support half-width is ``sqrt(3)*sigma`` as in the uniform box, but the
    cell-integrated one-dimensional profile is proportional to
    ``cos(pi*r/(2*h))`` inside ``|r| <= h``.  This removes the discontinuous
    box edge without adding a fitted tail or changing the declared variance
    scale.  It is a diagnostic source-shape candidate, not the default.
    """
    if len(sigma_xyz) != 3 or any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in sigma_xyz
    ):
        raise ValueError(
            "source cosine-tapered-box sigma_xyz values must be positive and finite"
        )

    def integrated_axis_weights(
        coordinates: np.ndarray, centre: float, sigma: float,
    ) -> np.ndarray:
        if coordinates.ndim != 1 or coordinates.size < 2:
            raise ValueError("source cosine-tapered box requires at least two cells per axis")
        spacing = float(coordinates[1] - coordinates[0])
        if not np.allclose(
            np.diff(coordinates), spacing, rtol=1.0e-12, atol=1.0e-14,
        ):
            raise ValueError(
                "source cosine-tapered box requires uniform Cartesian coordinates"
            )
        half_width = math.sqrt(3.0) * float(sigma)
        lower = np.maximum(coordinates - 0.5 * spacing, centre - half_width)
        upper = np.minimum(coordinates + 0.5 * spacing, centre + half_width)
        overlap = upper > lower
        phase_lower = math.pi * (lower - centre) / (2.0 * half_width)
        phase_upper = math.pi * (upper - centre) / (2.0 * half_width)
        weights = np.zeros_like(coordinates, dtype=float)
        weights[overlap] = (2.0 * half_width / math.pi) * (
            np.sin(phase_upper[overlap]) - np.sin(phase_lower[overlap])
        )
        return np.maximum(weights, 0.0)

    wx = integrated_axis_weights(x, position[0], float(sigma_xyz[0]))
    wy = integrated_axis_weights(y, position[1], float(sigma_xyz[1]))
    wz = integrated_axis_weights(z, position[2], float(sigma_xyz[2]))
    weights = wz[:, None, None] * wy[None, :, None] * wx[None, None, :]
    weights[solid] = 0.0
    total = float(np.sum(weights))
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("source cosine-tapered box has no fluid-cell support")
    return weights / total


def _receptor_stencil(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, solid: np.ndarray,
    receptor: FastObstacleReceptor3D,
    *,
    vertical_sampling: str = "cell_center",
) -> tuple[tuple[int, int, int, float], ...]:
    if vertical_sampling not in {"cell_center", "ground_linear"}:
        raise ValueError(
            "vertical_sampling must be 'cell_center' or 'ground_linear'"
        )
    def bracket(coordinates: np.ndarray, value: float) -> tuple[int, int, float]:
        half = 0.5 * (coordinates[1] - coordinates[0])
        lower_boundary = coordinates[0] - half
        upper_boundary = coordinates[-1] + half
        tolerance = 1.0e-12 * max(1.0, abs(float(value)), abs(float(upper_boundary)))
        if value < lower_boundary - tolerance or value > upper_boundary + tolerance:
            raise ValueError(f"receptor {receptor.label!r} lies outside the domain")
        if value <= coordinates[0]:
            return 0, 0, 0.0
        if value >= coordinates[-1]:
            last = len(coordinates) - 1
            return last, last, 0.0
        upper = int(np.searchsorted(coordinates, value, side="right"))
        lower = upper - 1
        fraction = float(
            (value - coordinates[lower]) / (coordinates[upper] - coordinates[lower])
        )
        return lower, upper, fraction

    ix = bracket(x, receptor.position_m[0])
    iy = bracket(y, receptor.position_m[1])
    iz = bracket(z, receptor.position_m[2])
    if (
        vertical_sampling == "ground_linear"
        and receptor.position_m[2] < float(z[0])
        and len(z) >= 2
    ):
        # Reconstruct the sub-cell value between the ground and the first
        # cell centre using the first-centre/second-centre slope.  This is an
        # observation operator only; transport remains cell-centred and
        # non-negative values are clipped at sampling time by the caller.
        fraction = float(
            (receptor.position_m[2] - z[0]) / (z[1] - z[0])
        )
        iz = (0, 1, fraction)
    combined: dict[tuple[int, int, int], float] = {}
    for kz, wz in ((iz[0], 1.0 - iz[2]), (iz[1], iz[2])):
        for jy, wy in ((iy[0], 1.0 - iy[2]), (iy[1], iy[2])):
            for ix_value, wx in ((ix[0], 1.0 - ix[2]), (ix[1], ix[2])):
                weight = wx * wy * wz
                if abs(weight) > 1.0e-15:
                    key = (kz, jy, ix_value)
                    combined[key] = combined.get(key, 0.0) + weight
    if any(solid[index] for index in combined):
        raise ValueError(
            f"receptor {receptor.label!r} interpolation stencil intersects an obstacle"
        )
    return tuple((*index, weight) for index, weight in combined.items())


def _transport_step(
    concentration: np.ndarray,
    *,
    u_face: np.ndarray,
    v_face: np.ndarray,
    w_face: np.ndarray,
    solid: np.ndarray,
    diffusivity: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    dt: float,
    spacing: np.ndarray,
    signed_boundary_flux: bool = False,
    advection_scheme: str = "upwind",
    open_boundary_faces: tuple[str, ...] = ("x_min", "x_max"),
    prescribed_x_min_inflow_state: np.ndarray | None = None,
    prescribed_x_min_inflow_flux: np.ndarray | None = None,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    face_flux_callback=None,
) -> tuple[np.ndarray, float]:
    if face_flux_callback is not None and not callable(face_flux_callback):
        raise TypeError("face_flux_callback must be callable")
    if advection_scheme not in {
        "upwind", "muscl_van_leer", "muscl_ssprk2",
    }:
        raise ValueError(
            "advection_scheme must be 'upwind', 'muscl_van_leer', or "
            "'muscl_ssprk2'"
        )
    valid_boundary_faces = {
        "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
    }
    if (
        not isinstance(open_boundary_faces, tuple)
        or len(set(open_boundary_faces)) != len(open_boundary_faces)
        or any(face not in valid_boundary_faces for face in open_boundary_faces)
    ):
        raise ValueError(
            "open_boundary_faces must be a tuple of unique Cartesian face names"
        )
    open_faces = frozenset(open_boundary_faces)
    inlet_state = None
    if prescribed_x_min_inflow_state is not None:
        inlet_state = np.asarray(prescribed_x_min_inflow_state, dtype=float)
        if inlet_state.shape != concentration.shape[:2] or not np.all(np.isfinite(inlet_state)):
            raise ValueError(
                "prescribed_x_min_inflow_state must be finite with shape (nz, ny)"
            )
    inlet_flux = None
    if prescribed_x_min_inflow_flux is not None:
        inlet_flux = np.asarray(prescribed_x_min_inflow_flux, dtype=float)
        if inlet_flux.shape != concentration.shape[:2] or not np.all(np.isfinite(inlet_flux)):
            raise ValueError(
                "prescribed_x_min_inflow_flux must be finite with shape (nz, ny)"
            )
        if np.any(inlet_flux < 0.0):
            raise ValueError("prescribed_x_min_inflow_flux must be non-negative")
    if blocked_faces is None:
        blocked_x, blocked_y, blocked_z = (
            np.zeros_like(u_face, dtype=bool),
            np.zeros_like(v_face, dtype=bool),
            np.zeros_like(w_face, dtype=bool),
        )
    else:
        if not isinstance(blocked_faces, tuple) or len(blocked_faces) != 3:
            raise ValueError("blocked_faces must be a tuple of x, y, and z face masks")
        expected_shapes = (u_face.shape, v_face.shape, w_face.shape)
        masks = tuple(np.asarray(mask, dtype=bool) for mask in blocked_faces)
        if any(mask.shape != shape for mask, shape in zip(masks, expected_shapes)):
            raise ValueError("blocked face masks must match their staggered face arrays")
        blocked_x, blocked_y, blocked_z = masks
    supplied_diffusivities = (
        diffusivity if isinstance(diffusivity, tuple)
        else (diffusivity, diffusivity, diffusivity)
    )
    if len(supplied_diffusivities) != 3:
        raise ValueError("directional diffusivity must contain (Dx, Dy, Dz)")
    directional_fields: list[np.ndarray] = []
    for supplied in supplied_diffusivities:
        if np.ndim(supplied) == 0:
            scalar_diffusivity = float(supplied)
            if not math.isfinite(scalar_diffusivity) or scalar_diffusivity < 0.0:
                raise ValueError("diffusivity must be finite and non-negative")
            field = np.full_like(concentration, scalar_diffusivity, dtype=float)
        else:
            field = np.asarray(supplied, dtype=float)
            if field.shape != concentration.shape:
                raise ValueError(
                    "each diffusivity field must have the same shape as concentration"
                )
            if not np.all(np.isfinite(field)) or np.any(field < 0.0):
                raise ValueError("diffusivity field must be finite and non-negative")
        directional_fields.append(field)
    diffusivity_x_field, diffusivity_y_field, diffusivity_z_field = directional_fields

    def limited_slope(state: np.ndarray, axis: int) -> np.ndarray:
        """Return a van-Leer slope, zeroed at boundaries and solid neighbours."""
        slope = np.zeros_like(state)
        center = [slice(None)] * 3
        lower = [slice(None)] * 3
        upper = [slice(None)] * 3
        center[axis] = slice(1, -1)
        lower[axis] = slice(None, -2)
        upper[axis] = slice(2, None)
        center_key = tuple(center)
        lower_key = tuple(lower)
        upper_key = tuple(upper)
        backward = state[center_key] - state[lower_key]
        forward = state[upper_key] - state[center_key]
        same_sign = backward * forward > 0.0
        denominator = backward + forward
        local = np.zeros_like(backward)
        valid_denominator = same_sign & (np.abs(denominator) > 1.0e-300)
        local[valid_denominator] = (
            2.0 * backward[valid_denominator] * forward[valid_denominator]
            / denominator[valid_denominator]
        )
        valid_stencil = (
            (~solid[center_key])
            & (~solid[lower_key])
            & (~solid[upper_key])
        )
        local[~valid_stencil] = 0.0
        slope[center_key] = local
        return slope

    dx, dy, dz = (float(value) for value in spacing)

    def tendency_and_boundary_rate(
        state: np.ndarray, *, hancock: bool,
    ) -> tuple[np.ndarray, float]:
        if advection_scheme == "upwind":
            slope_x = slope_y = slope_z = None
        else:
            slope_x = limited_slope(state, 2)
            slope_y = limited_slope(state, 1)
            slope_z = limited_slope(state, 0)

        flux_x = np.zeros_like(u_face)
        fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:])
        interior_u = u_face[:, :, 1:-1]
        if slope_x is None:
            left_x = state[:, :, :-1]
            right_x = state[:, :, 1:]
        elif hancock:
            courant_x = interior_u * dt / dx
            left_x = state[:, :, :-1] + (
                0.5 * (1.0 - courant_x) * slope_x[:, :, :-1]
            )
            right_x = state[:, :, 1:] - (
                0.5 * (1.0 + courant_x) * slope_x[:, :, 1:]
            )
        else:
            left_x = state[:, :, :-1] + 0.5 * slope_x[:, :, :-1]
            right_x = state[:, :, 1:] - 0.5 * slope_x[:, :, 1:]
        advective = np.where(
            interior_u >= 0.0,
            interior_u * left_x,
            interior_u * right_x,
        )
        diffusivity_x = 0.5 * (
            diffusivity_x_field[:, :, :-1]
            + diffusivity_x_field[:, :, 1:]
        )
        diffusive = -diffusivity_x * (
            state[:, :, 1:] - state[:, :, :-1]
        ) / dx
        flux_x[:, :, 1:-1][fluid_x] = (advective + diffusive)[fluid_x]
        if "x_min" in open_faces:
            inlet_velocity = u_face[:, :, 0]
            flux_x[:, :, 0] = np.minimum(inlet_velocity, 0.0) * state[:, :, 0]
            if inlet_flux is not None:
                flux_x[:, :, 0] += inlet_flux
            elif inlet_state is not None:
                flux_x[:, :, 0] += np.maximum(inlet_velocity, 0.0) * inlet_state
        if "x_max" in open_faces:
            flux_x[:, :, -1] = (
                np.maximum(u_face[:, :, -1], 0.0) * state[:, :, -1]
            )

        flux_y = np.zeros_like(v_face)
        fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :])
        interior_v = v_face[:, 1:-1, :]
        if slope_y is None:
            left_y = state[:, :-1, :]
            right_y = state[:, 1:, :]
        elif hancock:
            courant_y = interior_v * dt / dy
            left_y = state[:, :-1, :] + (
                0.5 * (1.0 - courant_y) * slope_y[:, :-1, :]
            )
            right_y = state[:, 1:, :] - (
                0.5 * (1.0 + courant_y) * slope_y[:, 1:, :]
            )
        else:
            left_y = state[:, :-1, :] + 0.5 * slope_y[:, :-1, :]
            right_y = state[:, 1:, :] - 0.5 * slope_y[:, 1:, :]
        advective = np.where(
            interior_v >= 0.0,
            interior_v * left_y,
            interior_v * right_y,
        )
        diffusivity_y = 0.5 * (
            diffusivity_y_field[:, :-1, :]
            + diffusivity_y_field[:, 1:, :]
        )
        diffusive = -diffusivity_y * (
            state[:, 1:, :] - state[:, :-1, :]
        ) / dy
        flux_y[:, 1:-1, :][fluid_y] = (advective + diffusive)[fluid_y]
        if "y_min" in open_faces:
            flux_y[:, 0, :] = (
                np.minimum(v_face[:, 0, :], 0.0) * state[:, 0, :]
            )
        if "y_max" in open_faces:
            flux_y[:, -1, :] = (
                np.maximum(v_face[:, -1, :], 0.0) * state[:, -1, :]
            )

        flux_z = np.zeros_like(w_face)
        fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :])
        interior_w = w_face[1:-1, :, :]
        if slope_z is None:
            left_z = state[:-1, :, :]
            right_z = state[1:, :, :]
        elif hancock:
            courant_z = interior_w * dt / dz
            left_z = state[:-1, :, :] + (
                0.5 * (1.0 - courant_z) * slope_z[:-1, :, :]
            )
            right_z = state[1:, :, :] - (
                0.5 * (1.0 + courant_z) * slope_z[1:, :, :]
            )
        else:
            left_z = state[:-1, :, :] + 0.5 * slope_z[:-1, :, :]
            right_z = state[1:, :, :] - 0.5 * slope_z[1:, :, :]
        advective = np.where(
            interior_w >= 0.0,
            interior_w * left_z,
            interior_w * right_z,
        )
        diffusivity_z = 0.5 * (
            diffusivity_z_field[:-1, :, :]
            + diffusivity_z_field[1:, :, :]
        )
        diffusive = -diffusivity_z * (
            state[1:, :, :] - state[:-1, :, :]
        ) / dz
        flux_z[1:-1, :, :][fluid_z] = (advective + diffusive)[fluid_z]
        if "z_min" in open_faces:
            flux_z[0, :, :] = (
                np.minimum(w_face[0, :, :], 0.0) * state[0, :, :]
            )
        if "z_max" in open_faces:
            flux_z[-1, :, :] = (
                np.maximum(w_face[-1, :, :], 0.0) * state[-1, :, :]
            )

        # A geometric wall face carries neither advective nor diffusive scalar
        # flux. Apply this after boundary reconstruction so no later branch can
        # reopen a face masked by the embedded geometry.
        flux_x[blocked_x] = 0.0
        flux_y[blocked_y] = 0.0
        flux_z[blocked_z] = 0.0

        if face_flux_callback is not None:
            face_flux_callback(flux_x.copy(), flux_y.copy(), flux_z.copy())
        tendency = -(
            (flux_x[:, :, 1:] - flux_x[:, :, :-1]) / dx
            + (flux_y[:, 1:, :] - flux_y[:, :-1, :]) / dy
            + (flux_z[1:, :, :] - flux_z[:-1, :, :]) / dz
        )
        tendency[solid] = 0.0
        if signed_boundary_flux:
            boundary_rate = (
                np.sum(flux_x[:, :, -1]) - np.sum(flux_x[:, :, 0])
            ) * dy * dz + (
                np.sum(flux_y[:, -1, :]) - np.sum(flux_y[:, 0, :])
            ) * dx * dz + (
                np.sum(flux_z[-1, :, :]) - np.sum(flux_z[0, :, :])
            ) * dx * dy
        else:
            boundary_rate = (
                np.sum(np.maximum(flux_x[:, :, -1], 0.0))
                - np.sum(np.minimum(flux_x[:, :, 0], 0.0))
            ) * dy * dz + (
                np.sum(np.maximum(flux_y[:, -1, :], 0.0))
                - np.sum(np.minimum(flux_y[:, 0, :], 0.0))
            ) * dx * dz + (
                np.sum(np.maximum(flux_z[-1, :, :], 0.0))
                - np.sum(np.minimum(flux_z[0, :, :], 0.0))
            ) * dx * dy
        return tendency, float(boundary_rate)

    if advection_scheme == "muscl_ssprk2":
        first_tendency, first_boundary_rate = tendency_and_boundary_rate(
            concentration, hancock=False,
        )
        stage = concentration + dt * first_tendency
        stage[solid] = 0.0
        second_tendency, second_boundary_rate = tendency_and_boundary_rate(
            stage, hancock=False,
        )
        result = 0.5 * concentration + 0.5 * (
            stage + dt * second_tendency
        )
        boundary_rate = 0.5 * (
            first_boundary_rate + second_boundary_rate
        )
    else:
        tendency, boundary_rate = tendency_and_boundary_rate(
            concentration,
            hancock=(advection_scheme == "muscl_van_leer"),
        )
        result = concentration + dt * tendency
    result[solid] = 0.0
    return result, float(dt * boundary_rate)


def solve_fast_obstacle_transport_3d(
    config: FastObstacleTransport3DConfig,
    *,
    receptors: Sequence[FastObstacleReceptor3D] = (),
    face_flux_callback=None,
    interface_state_callback=None,
    interface_momentum_callback=None,
    interface_inlet_flux_callback=None,
    interface_state_mass_callback=None,
) -> FastObstacleTransport3DResult:
    if not isinstance(config, FastObstacleTransport3DConfig):
        raise TypeError("config must be FastObstacleTransport3DConfig")
    if face_flux_callback is not None and not callable(face_flux_callback):
        raise TypeError("face_flux_callback must be callable")
    if interface_state_callback is not None and not callable(interface_state_callback):
        raise TypeError("interface_state_callback must be callable")
    if interface_momentum_callback is not None and not callable(interface_momentum_callback):
        raise TypeError("interface_momentum_callback must be callable")
    if interface_inlet_flux_callback is not None and not callable(interface_inlet_flux_callback):
        raise TypeError("interface_inlet_flux_callback must be callable")
    if interface_state_mass_callback is not None and not callable(interface_state_mass_callback):
        raise TypeError("interface_state_mass_callback must be callable")
    checked_receptors = tuple(receptors)
    if any(not isinstance(item, FastObstacleReceptor3D) for item in checked_receptors):
        raise TypeError("receptors must contain FastObstacleReceptor3D values")
    if len({item.label for item in checked_receptors}) != len(checked_receptors):
        raise ValueError("receptor labels must be unique")

    wall_started = time.perf_counter()
    x, y, z, solid, spacing = _grid_and_mask(config)
    blocked_faces = _circular_arc_blocked_faces(config, x, y, z)
    flow_started = time.perf_counter()
    (
        u_face, v_face, w_face, u, v, w, pressure_coefficient,
        divergence_before, divergence_after, projection_iterations,
    ) = _potential_face_flow(
        config, x, y, z, solid, spacing, blocked_faces=blocked_faces,
    )
    flow_runtime = time.perf_counter() - flow_started
    dx, dy, dz = (float(value) for value in spacing)
    cell_volume = dx * dy * dz
    advective_number = config.time_step_s * (
        float(np.max(np.abs(u_face))) / dx
        + float(np.max(np.abs(v_face))) / dy
        + float(np.max(np.abs(w_face))) / dz
    )
    diffusion_number = 2.0 * config.diffusivity_m2_s * config.time_step_s * (
        1.0 / dx**2 + 1.0 / dy**2 + 1.0 / dz**2
    )
    if advective_number + diffusion_number > 1.0 + 1.0e-12:
        raise ValueError(
            "fast 3-D obstacle scalar stability limit exceeded: "
            f"advection={advective_number:.6g}, diffusion={diffusion_number:.6g}"
        )
    weights = (
        _source_weights_xyz(
            x, y, z, solid, config.source_position_m, config.source_sigma_xyz_m,
        )
        if config.source_sigma_xyz_m is not None
        else _source_weights(
            x, y, z, solid, config.source_position_m, config.source_sigma_m,
        )
    )
    stencils = [
        _receptor_stencil(x, y, z, solid, receptor)
        for receptor in checked_receptors
    ]
    steps = int(math.ceil(config.duration_s / config.time_step_s))
    dt = config.duration_s / steps
    concentration = np.zeros_like(solid, dtype=float)
    maximum_concentration = concentration.copy()
    times = [0.0]
    injected = [0.0]
    domain_mass = [0.0]
    outflow = [0.0]
    residual = [0.0]
    traces = {item.label: [0.0] for item in checked_receptors}
    fields = [concentration.copy()] if config.store_fields else None
    cumulative_injected = 0.0
    cumulative_outflow = 0.0
    transport_started = time.perf_counter()
    for step in range(steps):
        start = step * dt
        end = (step + 1) * dt
        source_mass = config.source_schedule.mass_between(start, end)
        if source_mass > 0.0:
            concentration += source_mass * weights / cell_volume
        cumulative_injected += source_mass
        if interface_state_callback is not None:
            updated = interface_state_callback(concentration.copy(), start, dt)
            updated = np.asarray(updated, dtype=float)
            if updated.shape != concentration.shape or not np.all(np.isfinite(updated)) or np.any(updated < 0.0):
                raise ValueError("interface_state_callback must return a finite non-negative field with solver shape")
            concentration = updated
        if interface_state_mass_callback is not None:
            callback_injected_mass = float(interface_state_mass_callback(start, dt))
            if not math.isfinite(callback_injected_mass) or callback_injected_mass < 0.0:
                raise ValueError("interface state mass callback must return finite non-negative mass")
            cumulative_injected += callback_injected_mass
        if interface_momentum_callback is not None:
            updated_momentum = interface_momentum_callback(
                u.copy(), v.copy(), w.copy(),
                u_face.copy(), v_face.copy(), w_face.copy(),
                start,
                dt,
            )
            if not isinstance(updated_momentum, tuple) or len(updated_momentum) != 6:
                raise TypeError(
                    "interface_momentum_callback must return six velocity arrays"
                )
            arrays = tuple(np.asarray(value, dtype=float) for value in updated_momentum)
            expected = (
                u.shape, v.shape, w.shape,
                u_face.shape, v_face.shape, w_face.shape,
            )
            if any(array.shape != shape for array, shape in zip(arrays, expected)):
                raise ValueError("interface momentum arrays must match solver shapes")
            if any(not np.all(np.isfinite(array)) for array in arrays):
                raise ValueError("interface momentum arrays must be finite")
            u, v, w, u_face, v_face, w_face = arrays
            u[solid] = v[solid] = w[solid] = 0.0
        inlet_flux = None
        if interface_inlet_flux_callback is not None:
            inlet_flux = interface_inlet_flux_callback(
                start, dt, u_face.copy(), concentration.copy(), solid.copy(),
            )
            inlet_flux = np.asarray(inlet_flux, dtype=float)
            if inlet_flux.shape != concentration.shape[:2] or not np.all(np.isfinite(inlet_flux)):
                raise ValueError(
                    "interface_inlet_flux_callback must return finite field with shape (nz, ny)"
                )
            if np.any(inlet_flux < 0.0):
                raise ValueError("interface inlet flux callback must return non-negative flux")
            inlet_flux[solid[:, :, 0]] = 0.0
            # The callback is a conservative boundary source, not a field
            # overwrite. Include its face-integrated mass in the same ledger
            # used for source injection and outflow residuals.
            cumulative_injected += float(np.sum(inlet_flux) * dy * dz * dt)
        concentration, step_outflow = _transport_step(
            concentration,
            u_face=u_face, v_face=v_face, w_face=w_face,
            solid=solid, diffusivity=config.diffusivity_m2_s,
            dt=dt, spacing=spacing,
            blocked_faces=blocked_faces,
            open_boundary_faces=config.open_boundary_faces,
            prescribed_x_min_inflow_flux=inlet_flux,
            face_flux_callback=face_flux_callback,
        )
        minimum = float(np.min(concentration))
        if minimum < -1.0e-12:
            raise RuntimeError(
                f"fast 3-D obstacle concentration became negative: {minimum:.6g}"
            )
        concentration = np.maximum(concentration, 0.0)
        np.maximum(maximum_concentration, concentration, out=maximum_concentration)
        cumulative_outflow += step_outflow
        mass = float(np.sum(concentration) * cell_volume)
        times.append(end)
        injected.append(cumulative_injected)
        domain_mass.append(mass)
        outflow.append(cumulative_outflow)
        residual.append(mass + cumulative_outflow - cumulative_injected)
        for receptor, stencil in zip(checked_receptors, stencils):
            traces[receptor.label].append(float(sum(
                weight * concentration[k, j, i]
                for k, j, i, weight in stencil
            )))
        if fields is not None:
            fields.append(concentration.copy())
    transport_runtime = time.perf_counter() - transport_started
    result_concentration = concentration.copy()
    result_concentration[solid] = np.nan
    result_maximum = maximum_concentration.copy()
    result_maximum[solid] = np.nan
    return FastObstacleTransport3DResult(
        schema=FAST_OBSTACLE_TRANSPORT_3D_SCHEMA,
        x_m=x, y_m=y, z_m=z, solid_mask=solid,
        velocity_u_m_s=u, velocity_v_m_s=v, velocity_w_m_s=w,
        u_face_m_s=u_face, v_face_m_s=v_face, w_face_m_s=w_face,
        pressure_coefficient=pressure_coefficient,
        time_s=np.asarray(times),
        concentration_kg_m3=result_concentration,
        maximum_concentration_kg_m3=result_maximum,
        concentration_fields_kg_m3=(None if fields is None else np.asarray(fields)),
        receptor_traces_kg_m3={
            label: np.asarray(values) for label, values in traces.items()
        },
        injected_mass_kg=np.asarray(injected),
        domain_mass_kg=np.asarray(domain_mass),
        outflow_mass_kg=np.asarray(outflow),
        mass_residual_kg=np.asarray(residual),
        diagnostics={
            "fidelity_tier": "tier_1_finite_width_3d_potential",
            "applicability": "conditional",
            "quantitative_prediction_allowed": False,
            "nonlocal_scalar_routing": False,
            "local_conservative_scalar_transport": True,
            "interface_momentum_callback_applied": interface_momentum_callback is not None,
            "interface_inlet_flux_callback_applied": interface_inlet_flux_callback is not None,
            "interface_state_mass_callback_applied": interface_state_mass_callback is not None,
            "flow_model": "incompressible irrotational Helmholtz projection",
            "scalar_model": "local conservative 3-D finite-volume advection-diffusion",
            "source_spatial_discretization": (
                "cell-integrated separable Gaussian over Cartesian cells"
            ),
            "source_sigma_m": config.source_sigma_m,
            "source_sigma_xyz_m": (
                None if config.source_sigma_xyz_m is None
                else list(config.source_sigma_xyz_m)
            ),
            "obstacle_count": len(config.obstacles),
            "embedded_surface_blocked_face_count": int(sum(
                np.count_nonzero(mask) for mask in blocked_faces
            )),
            "solid_velocity_boundary": "zero normal flux; inviscid slip tangential flow",
            "solid_scalar_boundary": "zero normal advective and diffusive flux",
            "inlet_boundary": (
                "measured/interpolated velocity profile"
                if config.inlet_profile_z_m is not None
                else "prescribed uniform velocity"
            ),
            "open_boundary_faces": list(config.open_boundary_faces),
            "minimum_inlet_speed_m_s": float(np.min(u_face[:, :, 0])),
            "maximum_inlet_speed_m_s": float(np.max(u_face[:, :, 0])),
            "finite_lateral_bypass": True,
            "recirculating_wake_resolved": False,
            "thermal_buoyancy_feedback": False,
            "cell_size_m": [dx, dy, dz],
            "time_step_s": dt,
            "advective_number": advective_number,
            "diffusion_number": diffusion_number,
            "projection_iterations": projection_iterations,
            "maximum_absolute_divergence_before_projection_s_1": divergence_before,
            "maximum_absolute_divergence_s_1": divergence_after,
            "maximum_mass_residual_kg": float(np.max(np.abs(residual))),
            "flow_solve_runtime_s": flow_runtime,
            "transport_runtime_s": transport_runtime,
            "total_runtime_s": time.perf_counter() - wall_started,
            "warnings": [
                "potential flow cannot resolve separation, recirculation or turbulent mixing",
                "slip walls omit the no-slip boundary layer",
                "isothermal transport omits cryogenic heat transfer and buoyancy feedback",
                "quantitative facility use is withheld pending resolved 3-D validation",
            ],
        },
        ambient_temperature_k=config.ambient_temperature_k,
        ambient_pressure_pa=config.ambient_pressure_pa,
    )


__all__ = [
    "FAST_OBSTACLE_TRANSPORT_3D_SCHEMA",
    "FastObstacleReceptor3D",
    "FastObstacleTransport3DConfig",
    "FastObstacleTransport3DResult",
    "solve_fast_obstacle_transport_3d",
]
