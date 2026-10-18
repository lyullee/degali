"""Discrete variable-density mass-flux projection on a Cartesian 3-D grid.

The operator enforces the finite-volume continuity equation

    (rho_new - rho_old) / dt + div(m_dot) = S_rho

in every fluid cell.  Prescribed inlet and impermeable external/solid faces
keep their declared flux.  One or more declared pressure-open faces are
zero-potential pressure references whose boundary mass flux is corrected, so
the cellwise equation and the global mass budget close together.

This is a numerical building block for a future low-Mach obstacle solver.  It
does not by itself advance momentum, species, enthalpy, or an equation of
state and therefore is not a dispersion model.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from threading import Lock

import numpy as np
import pyamg
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg
from scipy.sparse.linalg import splu
from scipy.sparse.csgraph import connected_components


VARIABLE_DENSITY_MASS_FLUX_PROJECTION_SCHEMA = (
    "degali.variable-density-mass-flux-projection-3d.v4"
)

_BOUNDARY_FACES = frozenset({
    "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
})


@dataclass(frozen=True)
class VariableDensityMassFluxProjection3DResult:
    schema: str
    mass_flux_x_kg_m2_s: np.ndarray
    mass_flux_y_kg_m2_s: np.ndarray
    mass_flux_z_kg_m2_s: np.ndarray
    correction_potential_kg_m_s: np.ndarray
    pressure_correction_pa: np.ndarray
    target_mass_flux_divergence_kg_m3_s: np.ndarray
    continuity_residual_before_kg_m3_s: np.ndarray
    continuity_residual_after_kg_m3_s: np.ndarray
    solver_method: str
    solver_iterations: int
    preconditioner: str
    preconditioner_level_count: int | None
    preconditioner_operator_complexity: float | None
    pyamg_version: str | None
    maximum_absolute_residual_before_kg_m3_s: float
    maximum_absolute_residual_after_kg_m3_s: float
    domain_mass_change_rate_kg_s: float
    integrated_source_rate_kg_s: float
    net_boundary_outflow_rate_kg_s: float
    global_mass_balance_residual_kg_s: float
    open_boundary_faces: tuple[str, ...]


@dataclass(frozen=True)
class _ProjectionGeometry:
    matrix: sparse.csr_matrix
    diagonal: np.ndarray


@dataclass(frozen=True)
class _DirectProjectionFactorization:
    factorization: object
    lock: Lock


@dataclass(frozen=True)
class _AmgProjectionPreconditioner:
    preconditioner: LinearOperator
    level_count: int
    operator_complexity: float


@lru_cache(maxsize=16)
def _projection_geometry(
    shape: tuple[int, int, int],
    solid_bytes: bytes,
    spacing: tuple[float, float, float],
    open_boundary_faces: tuple[str, ...],
    blocked_face_bytes: tuple[bytes, bytes, bytes],
) -> _ProjectionGeometry:
    """Build the immutable Poisson geometry once per grid and solid mask."""
    solid = np.frombuffer(solid_bytes, dtype=np.bool_).reshape(shape)
    nz, ny, nx = shape
    dx, dy, dz = spacing
    blocked_x = np.frombuffer(
        blocked_face_bytes[0], dtype=np.bool_,
    ).reshape((nz, ny, nx + 1))
    blocked_y = np.frombuffer(
        blocked_face_bytes[1], dtype=np.bool_,
    ).reshape((nz, ny + 1, nx))
    blocked_z = np.frombuffer(
        blocked_face_bytes[2], dtype=np.bool_,
    ).reshape((nz + 1, ny, nx))
    fluid = ~solid
    positions = np.argwhere(fluid)
    lookup = -np.ones_like(solid, dtype=np.int64)
    lookup[fluid] = np.arange(len(positions), dtype=np.int64)
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    neighbours = (
        (0, 0, -1, 1.0 / dx**2, 0), (0, 0, 1, 1.0 / dx**2, 1),
        (0, -1, 0, 1.0 / dy**2, 2), (0, 1, 0, 1.0 / dy**2, 3),
        (-1, 0, 0, 1.0 / dz**2, 4), (1, 0, 0, 1.0 / dz**2, 5),
    )
    for row, (k, j, i) in enumerate(positions):
        diagonal = 0.0
        for dk, dj, di, coefficient, face_direction in neighbours:
            nk, nj, ni = k + dk, j + dj, i + di
            if 0 <= nk < nz and 0 <= nj < ny and 0 <= ni < nx:
                if solid[nk, nj, ni]:
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
                rows.append(row)
                columns.append(int(lookup[nk, nj, ni]))
                values.append(-coefficient)
            else:
                face = None
                distance = None
                if di == -1 and i == 0:
                    face, distance = "x_min", dx
                elif di == 1 and i == nx - 1:
                    face, distance = "x_max", dx
                elif dj == -1 and j == 0:
                    face, distance = "y_min", dy
                elif dj == 1 and j == ny - 1:
                    face, distance = "y_max", dy
                elif dk == -1 and k == 0:
                    face, distance = "z_min", dz
                elif dk == 1 and k == nz - 1:
                    face, distance = "z_max", dz
                face_blocked = (
                    (face == "x_min" and blocked_x[k, j, 0])
                    or (face == "x_max" and blocked_x[k, j, nx])
                    or (face == "y_min" and blocked_y[k, 0, i])
                    or (face == "y_max" and blocked_y[k, ny, i])
                    or (face == "z_min" and blocked_z[0, j, i])
                    or (face == "z_max" and blocked_z[nz, j, i])
                )
                if face in open_boundary_faces and not face_blocked:
                    # Cell centre to the zero-potential pressure face is d/2.
                    assert distance is not None
                    diagonal += 2.0 / distance**2
        rows.append(row)
        columns.append(row)
        values.append(diagonal)
    matrix = sparse.csr_matrix(
        (np.asarray(values), (np.asarray(rows), np.asarray(columns))),
        shape=(len(positions), len(positions)),
    )
    diagonal = matrix.diagonal()
    if np.any(diagonal <= 0.0):
        raise ValueError(
            "every connected fluid component must reach a declared pressure-open "
            "boundary through unblocked faces"
        )
    diagonal.setflags(write=False)
    return _ProjectionGeometry(matrix=matrix, diagonal=diagonal)


@lru_cache(maxsize=4)
def _projection_component_labels(
    shape: tuple[int, int, int],
    solid_bytes: bytes,
    spacing: tuple[float, float, float],
    open_boundary_faces: tuple[str, ...],
    blocked_face_bytes: tuple[bytes, bytes, bytes],
) -> tuple[int, np.ndarray]:
    """Cache the static fluid connectivity for a fixed projection geometry."""
    matrix = _projection_geometry(
        shape, solid_bytes, spacing, open_boundary_faces, blocked_face_bytes,
    ).matrix
    adjacency = matrix.copy()
    adjacency.setdiag(0.0)
    adjacency.eliminate_zeros()
    component_count, labels = connected_components(
        adjacency, directed=False,
    )
    solid = np.frombuffer(solid_bytes, dtype=np.bool_).reshape(shape)
    grid_labels = np.full(shape, -1, dtype=np.int64)
    grid_labels[~solid] = labels
    grid_labels.setflags(write=False)
    return component_count, grid_labels


@lru_cache(maxsize=2)
def _direct_projection_factorization(
    shape: tuple[int, int, int],
    solid_bytes: bytes,
    spacing: tuple[float, float, float],
    open_boundary_faces: tuple[str, ...],
    blocked_face_bytes: tuple[bytes, bytes, bytes],
) -> _DirectProjectionFactorization:
    geometry = _projection_geometry(
        shape, solid_bytes, spacing, open_boundary_faces, blocked_face_bytes,
    )
    return _DirectProjectionFactorization(
        factorization=splu(geometry.matrix.tocsc()),
        lock=Lock(),
    )


@lru_cache(maxsize=4)
def _amg_projection_preconditioner(
    shape: tuple[int, int, int],
    solid_bytes: bytes,
    spacing: tuple[float, float, float],
    open_boundary_faces: tuple[str, ...],
    blocked_face_bytes: tuple[bytes, bytes, bytes],
) -> _AmgProjectionPreconditioner:
    """Cache a symmetric smoothed-aggregation hierarchy by grid geometry."""
    geometry = _projection_geometry(
        shape, solid_bytes, spacing, open_boundary_faces, blocked_face_bytes,
    )
    hierarchy = pyamg.smoothed_aggregation_solver(
        geometry.matrix,
        symmetry="symmetric",
        max_coarse=50,
    )
    return _AmgProjectionPreconditioner(
        preconditioner=hierarchy.aspreconditioner(cycle="V"),
        level_count=len(hierarchy.levels),
        operator_complexity=float(hierarchy.operator_complexity()),
    )


def _as_finite_array(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


def _boundary_outflow_rate(
    flux_x: np.ndarray,
    flux_y: np.ndarray,
    flux_z: np.ndarray,
    spacing: tuple[float, float, float],
) -> float:
    dx, dy, dz = spacing
    return float(
        (np.sum(flux_x[:, :, -1]) - np.sum(flux_x[:, :, 0])) * dy * dz
        + (np.sum(flux_y[:, -1, :]) - np.sum(flux_y[:, 0, :])) * dx * dz
        + (np.sum(flux_z[-1, :, :]) - np.sum(flux_z[0, :, :])) * dx * dy
    )


def project_variable_density_mass_flux_3d(
    mass_flux_x_kg_m2_s: np.ndarray,
    mass_flux_y_kg_m2_s: np.ndarray,
    mass_flux_z_kg_m2_s: np.ndarray,
    *,
    density_old_kg_m3: np.ndarray,
    density_new_kg_m3: np.ndarray,
    mass_source_rate_kg_m3_s: np.ndarray | float,
    solid_mask: np.ndarray,
    spacing_m: tuple[float, float, float] | np.ndarray,
    time_step_s: float,
    relative_tolerance: float = 1.0e-11,
    absolute_tolerance: float = 1.0e-12,
    maximum_iterations: int = 50_000,
    linear_solver: str = "auto",
    open_boundary_faces: tuple[str, ...] = ("x_max",),
    initial_correction_potential_kg_m_s: np.ndarray | None = None,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> VariableDensityMassFluxProjection3DResult:
    """Project provisional face mass fluxes onto variable-density continuity."""
    solid = np.asarray(solid_mask)
    if solid.dtype != np.bool_ or solid.ndim != 3:
        raise ValueError("solid_mask must be a three-dimensional boolean array")
    nz, ny, nx = solid.shape
    if nx < 2 or ny < 2 or nz < 2:
        raise ValueError("projection grid requires at least two cells per direction")
    if not isinstance(open_boundary_faces, tuple) or not open_boundary_faces:
        raise ValueError("open_boundary_faces must be a non-empty tuple")
    if (
        len(set(open_boundary_faces)) != len(open_boundary_faces)
        or any(face not in _BOUNDARY_FACES for face in open_boundary_faces)
    ):
        raise ValueError(
            "open_boundary_faces must contain unique values selected from "
            "x_min, x_max, y_min, y_max, z_min, and z_max"
        )
    open_faces = tuple(sorted(open_boundary_faces))
    expected_face_shapes = (
        (nz, ny, nx + 1), (nz, ny + 1, nx), (nz + 1, ny, nx),
    )
    if blocked_faces is None:
        blocked_x_face, blocked_y_face, blocked_z_face = tuple(
            np.zeros(shape, dtype=bool) for shape in expected_face_shapes
        )
    else:
        if not isinstance(blocked_faces, tuple) or len(blocked_faces) != 3:
            raise ValueError("blocked_faces must be a tuple of x, y, and z face masks")
        masks = tuple(np.asarray(mask) for mask in blocked_faces)
        if any(
            mask.dtype != np.bool_ or mask.shape != expected
            for mask, expected in zip(masks, expected_face_shapes)
        ):
            raise ValueError("blocked face masks must be boolean staggered face arrays")
        blocked_x_face, blocked_y_face, blocked_z_face = masks
    old = _as_finite_array(density_old_kg_m3, solid.shape, "density_old_kg_m3")
    new = _as_finite_array(density_new_kg_m3, solid.shape, "density_new_kg_m3")
    if np.any(old[~solid] <= 0.0) or np.any(new[~solid] <= 0.0):
        raise ValueError("fluid-cell densities must be positive")
    source = np.broadcast_to(
        np.asarray(mass_source_rate_kg_m3_s, dtype=float), solid.shape,
    ).copy()
    if not np.all(np.isfinite(source)):
        raise ValueError("mass_source_rate_kg_m3_s must be finite and broadcastable")
    source[solid] = 0.0

    spacing_array = np.asarray(spacing_m, dtype=float)
    if spacing_array.shape != (3,) or np.any(~np.isfinite(spacing_array)) or np.any(
        spacing_array <= 0.0
    ):
        raise ValueError("spacing_m must contain three positive finite values")
    dx, dy, dz = (float(value) for value in spacing_array)
    spacing = (dx, dy, dz)
    dt = float(time_step_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_s must be positive and finite")
    if not math.isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise ValueError("relative_tolerance must be positive and finite")
    if not math.isfinite(absolute_tolerance) or absolute_tolerance < 0.0:
        raise ValueError("absolute_tolerance must be finite and non-negative")
    if (
        isinstance(maximum_iterations, bool)
        or not isinstance(maximum_iterations, int)
        or maximum_iterations < 1
    ):
        raise ValueError("maximum_iterations must be a positive integer")
    if linear_solver not in {"auto", "cached_direct", "cg", "amg_cg"}:
        raise ValueError(
            "linear_solver must be 'auto', 'cached_direct', 'cg', or 'amg_cg'"
        )

    flux_x = _as_finite_array(
        mass_flux_x_kg_m2_s, (nz, ny, nx + 1), "mass_flux_x_kg_m2_s",
    ).copy()
    flux_y = _as_finite_array(
        mass_flux_y_kg_m2_s, (nz, ny + 1, nx), "mass_flux_y_kg_m2_s",
    ).copy()
    flux_z = _as_finite_array(
        mass_flux_z_kg_m2_s, (nz + 1, ny, nx), "mass_flux_z_kg_m2_s",
    ).copy()

    blocked_x = (solid[:, :, :-1] | solid[:, :, 1:]) | blocked_x_face[:, :, 1:-1]
    blocked_y = (solid[:, :-1, :] | solid[:, 1:, :]) | blocked_y_face[:, 1:-1, :]
    blocked_z = (solid[:-1, :, :] | solid[1:, :, :]) | blocked_z_face[1:-1, :, :]
    blocked_values = (
        np.max(np.abs(flux_x[:, :, 1:-1][blocked_x]), initial=0.0),
        np.max(np.abs(flux_y[:, 1:-1, :][blocked_y]), initial=0.0),
        np.max(np.abs(flux_z[1:-1, :, :][blocked_z]), initial=0.0),
        np.max(np.abs(flux_x[blocked_x_face]), initial=0.0),
        np.max(np.abs(flux_y[blocked_y_face]), initial=0.0),
        np.max(np.abs(flux_z[blocked_z_face]), initial=0.0),
        np.max(np.abs(flux_x[:, :, 0][solid[:, :, 0]]), initial=0.0),
        np.max(np.abs(flux_x[:, :, -1][solid[:, :, -1]]), initial=0.0),
        np.max(np.abs(flux_y[:, 0, :][solid[:, 0, :]]), initial=0.0),
        np.max(np.abs(flux_y[:, -1, :][solid[:, -1, :]]), initial=0.0),
        np.max(np.abs(flux_z[0, :, :][solid[0, :, :]]), initial=0.0),
        np.max(np.abs(flux_z[-1, :, :][solid[-1, :, :]]), initial=0.0),
    )
    if max(blocked_values) > absolute_tolerance:
        raise ValueError("provisional mass flux penetrates an impermeable solid face")
    flux_x[:, :, 1:-1][blocked_x] = 0.0
    flux_y[:, 1:-1, :][blocked_y] = 0.0
    flux_z[1:-1, :, :][blocked_z] = 0.0
    flux_x[blocked_x_face] = 0.0
    flux_y[blocked_y_face] = 0.0
    flux_z[blocked_z_face] = 0.0
    flux_x[:, :, 0][solid[:, :, 0]] = 0.0
    flux_x[:, :, -1][solid[:, :, -1]] = 0.0
    flux_y[:, 0, :][solid[:, 0, :]] = 0.0
    flux_y[:, -1, :][solid[:, -1, :]] = 0.0
    flux_z[0, :, :][solid[0, :, :]] = 0.0
    flux_z[-1, :, :][solid[-1, :, :]] = 0.0

    density_rate = (new - old) / dt
    target_divergence = source - density_rate
    divergence_before = (
        (flux_x[:, :, 1:] - flux_x[:, :, :-1]) / dx
        + (flux_y[:, 1:, :] - flux_y[:, :-1, :]) / dy
        + (flux_z[1:, :, :] - flux_z[:-1, :, :]) / dz
    )
    residual_before = density_rate + divergence_before - source
    residual_before[solid] = 0.0

    fluid = ~solid
    solid_bytes = np.ascontiguousarray(solid).tobytes()
    blocked_face_bytes = tuple(
        np.ascontiguousarray(mask).tobytes()
        for mask in (blocked_x_face, blocked_y_face, blocked_z_face)
    )
    geometry = _projection_geometry(
        solid.shape, solid_bytes, spacing, open_faces, blocked_face_bytes,
    )
    matrix = geometry.matrix
    diagonal = geometry.diagonal
    component_count, component_by_cell = _projection_component_labels(
        solid.shape, solid_bytes, spacing, open_faces, blocked_face_bytes,
    )
    boundary_labels: set[int] = set()
    component_faces = {
        "x_min": component_by_cell[:, :, 0],
        "x_max": component_by_cell[:, :, -1],
        "y_min": component_by_cell[:, 0, :],
        "y_max": component_by_cell[:, -1, :],
        "z_min": component_by_cell[0, :, :],
        "z_max": component_by_cell[-1, :, :],
    }
    for face in open_faces:
        boundary_labels.update(
            int(label) for label in np.unique(component_faces[face])
            if label >= 0
        )
    if len(boundary_labels) != component_count:
        raise ValueError(
            "every connected fluid component must reach a declared pressure-open "
            "boundary through unblocked faces"
        )
    rhs = (target_divergence - divergence_before)[fluid]
    initial_solution = None
    if initial_correction_potential_kg_m_s is not None:
        initial_potential = _as_finite_array(
            initial_correction_potential_kg_m_s,
            solid.shape,
            "initial_correction_potential_kg_m_s",
        )
        initial_solution = initial_potential[fluid]
    selected_solver = linear_solver
    if selected_solver == "auto":
        selected_solver = (
            "cached_direct" if matrix.shape[0] <= 25_000 else "amg_cg"
        )
    if selected_solver == "cached_direct":
        preconditioner_name = "superlu_cached_direct"
        preconditioner_level_count = None
        preconditioner_operator_complexity = None
        pyamg_version = None
        direct = _direct_projection_factorization(
            solid.shape, solid_bytes, spacing, open_faces, blocked_face_bytes,
        )
        # SuperLU solve workspaces are protected because the immutable cached
        # factorization can be shared by simultaneous scenario evaluations.
        with direct.lock:
            solution = direct.factorization.solve(rhs)
        iterations = 1
        info = 0
    else:
        if selected_solver == "amg_cg":
            amg = _amg_projection_preconditioner(
                solid.shape, solid_bytes, spacing, open_faces, blocked_face_bytes,
            )
            preconditioner = amg.preconditioner
            preconditioner_name = "symmetric_smoothed_aggregation_v_cycle"
            preconditioner_level_count = amg.level_count
            preconditioner_operator_complexity = amg.operator_complexity
            pyamg_version = pyamg.__version__
        else:
            preconditioner = LinearOperator(
                matrix.shape, matvec=lambda vector: vector / diagonal,
            )
            preconditioner_name = "jacobi_diagonal"
            preconditioner_level_count = None
            preconditioner_operator_complexity = None
            pyamg_version = None
        iterations = 0

        def count_iteration(_vector: np.ndarray) -> None:
            nonlocal iterations
            iterations += 1

        solution, info = cg(
            matrix, rhs,
            x0=initial_solution,
            rtol=relative_tolerance,
            atol=absolute_tolerance,
            maxiter=maximum_iterations,
            M=preconditioner,
            callback=count_iteration,
        )
    if info != 0 or not np.all(np.isfinite(solution)):
        raise RuntimeError(
            "variable-density mass-flux projection failed to converge "
            f"({selected_solver} info={info})"
        )
    potential = np.zeros_like(old)
    potential[fluid] = solution

    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:]) & ~blocked_x_face[:, :, 1:-1]
    correction_x = (potential[:, :, 1:] - potential[:, :, :-1]) / dx
    flux_x[:, :, 1:-1][fluid_x] -= correction_x[fluid_x]
    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :]) & ~blocked_y_face[:, 1:-1, :]
    correction_y = (potential[:, 1:, :] - potential[:, :-1, :]) / dy
    flux_y[:, 1:-1, :][fluid_y] -= correction_y[fluid_y]
    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :]) & ~blocked_z_face[1:-1, :, :]
    correction_z = (potential[1:, :, :] - potential[:-1, :, :]) / dz
    flux_z[1:-1, :, :][fluid_z] -= correction_z[fluid_z]
    if "x_min" in open_faces:
        boundary_fluid = ~solid[:, :, 0] & ~blocked_x_face[:, :, 0]
        flux_x[:, :, 0][boundary_fluid] -= (
            2.0 * potential[:, :, 0][boundary_fluid] / dx
        )
    if "x_max" in open_faces:
        boundary_fluid = ~solid[:, :, -1] & ~blocked_x_face[:, :, -1]
        flux_x[:, :, -1][boundary_fluid] += (
            2.0 * potential[:, :, -1][boundary_fluid] / dx
        )
    if "y_min" in open_faces:
        boundary_fluid = ~solid[:, 0, :] & ~blocked_y_face[:, 0, :]
        flux_y[:, 0, :][boundary_fluid] -= (
            2.0 * potential[:, 0, :][boundary_fluid] / dy
        )
    if "y_max" in open_faces:
        boundary_fluid = ~solid[:, -1, :] & ~blocked_y_face[:, -1, :]
        flux_y[:, -1, :][boundary_fluid] += (
            2.0 * potential[:, -1, :][boundary_fluid] / dy
        )
    if "z_min" in open_faces:
        boundary_fluid = ~solid[0, :, :] & ~blocked_z_face[0, :, :]
        flux_z[0, :, :][boundary_fluid] -= (
            2.0 * potential[0, :, :][boundary_fluid] / dz
        )
    if "z_max" in open_faces:
        boundary_fluid = ~solid[-1, :, :] & ~blocked_z_face[-1, :, :]
        flux_z[-1, :, :][boundary_fluid] += (
            2.0 * potential[-1, :, :][boundary_fluid] / dz
        )

    divergence_after = (
        (flux_x[:, :, 1:] - flux_x[:, :, :-1]) / dx
        + (flux_y[:, 1:, :] - flux_y[:, :-1, :]) / dy
        + (flux_z[1:, :, :] - flux_z[:-1, :, :]) / dz
    )
    residual_after = density_rate + divergence_after - source
    residual_after[solid] = 0.0
    cell_volume = dx * dy * dz
    domain_mass_change_rate = float(np.sum(density_rate[fluid]) * cell_volume)
    integrated_source_rate = float(np.sum(source[fluid]) * cell_volume)
    net_boundary_outflow = _boundary_outflow_rate(
        flux_x, flux_y, flux_z, spacing,
    )
    global_residual = (
        domain_mass_change_rate + net_boundary_outflow - integrated_source_rate
    )
    return VariableDensityMassFluxProjection3DResult(
        schema=VARIABLE_DENSITY_MASS_FLUX_PROJECTION_SCHEMA,
        mass_flux_x_kg_m2_s=flux_x,
        mass_flux_y_kg_m2_s=flux_y,
        mass_flux_z_kg_m2_s=flux_z,
        correction_potential_kg_m_s=potential,
        pressure_correction_pa=potential / dt,
        target_mass_flux_divergence_kg_m3_s=target_divergence,
        continuity_residual_before_kg_m3_s=residual_before,
        continuity_residual_after_kg_m3_s=residual_after,
        solver_method=selected_solver,
        solver_iterations=iterations,
        preconditioner=preconditioner_name,
        preconditioner_level_count=preconditioner_level_count,
        preconditioner_operator_complexity=(
            preconditioner_operator_complexity
        ),
        pyamg_version=pyamg_version,
        maximum_absolute_residual_before_kg_m3_s=float(
            np.max(np.abs(residual_before[fluid]))
        ),
        maximum_absolute_residual_after_kg_m3_s=float(
            np.max(np.abs(residual_after[fluid]))
        ),
        domain_mass_change_rate_kg_s=domain_mass_change_rate,
        integrated_source_rate_kg_s=integrated_source_rate,
        net_boundary_outflow_rate_kg_s=net_boundary_outflow,
        global_mass_balance_residual_kg_s=global_residual,
        open_boundary_faces=open_faces,
    )


__all__ = [
    "VARIABLE_DENSITY_MASS_FLUX_PROJECTION_SCHEMA",
    "VariableDensityMassFluxProjection3DResult",
    "project_variable_density_mass_flux_3d",
]
