"""Coupled finite-volume outer-domain step for a conservative handoff.

The resolved-inner handoff supplies normal mass, species and sensible-
enthalpy fluxes together with the two tangential momentum fluxes on an x-min
plane.  This module applies all of them to the same finite-volume step.  The
momentum predictor receives the reconstructed three-component inlet velocity;
the binary-mixture corrector then holds the prescribed normal mass flux fixed
while projecting the remaining open faces to variable-density continuity.

This is an execution component.  It does not claim that an analytic handoff
test validates a field experiment; that requires a separate overlap and
observation audit.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .binary_mixture_corrector_3d import (
    BinaryMixtureCorrector3DResult,
    solve_binary_mixture_initial_cloud_step_3d,
)
from .binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
    BinaryIdealGasMixtureState,
)
from .conservative_nested_handoff_3d import (
    ConservativeNestedHandoff3DResult,
    NestedOuterInletFrame3D,
    nested_outer_inlet_frame_at,
)
from .variable_density_momentum_3d import (
    VariableDensityMomentum3DResult,
    advance_variable_density_momentum_3d,
)


NESTED_OUTER_BINARY_FLOW_3D_SCHEMA = "degali.nested-outer-binary-flow-3d.v1"


@dataclass(frozen=True)
class NestedOuterBinaryFlow3DState:
    """Conserved and primitive state carried between outer-domain steps."""

    schema: str
    time_s: float
    species_mass_concentration_kg_m3: np.ndarray
    sensible_enthalpy_density_j_m3: np.ndarray
    thermodynamic_state: BinaryIdealGasMixtureState
    u_m_s: np.ndarray
    v_m_s: np.ndarray
    w_m_s: np.ndarray
    u_face_m_s: np.ndarray
    v_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    correction_potential_kg_m_s: np.ndarray | None
    cumulative_species_boundary_outflow_mass_kg: float
    cumulative_sensible_enthalpy_boundary_outflow_j: float
    cumulative_species_mass_residual_kg: float
    cumulative_sensible_enthalpy_residual_j: float
    cumulative_prescribed_inlet_species_mass_kg: float
    cumulative_prescribed_inlet_sensible_enthalpy_j: float


@dataclass(frozen=True)
class NestedOuterBinaryFlow3DStepResult:
    """One fully coupled outer-domain update and its conservation records."""

    schema: str
    state: NestedOuterBinaryFlow3DState
    inlet_frame: NestedOuterInletFrame3D
    momentum: VariableDensityMomentum3DResult
    binary_corrector: BinaryMixtureCorrector3DResult
    diagnostics: dict[str, object]


def _finite_field(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    field = np.asarray(value, dtype=float)
    if field.shape != shape or np.any(~np.isfinite(field)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return field


def _plane(
    value: float | np.ndarray, nz: int, ny: int, name: str,
) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        result = np.full((nz, ny), float(array), dtype=float)
    elif array.shape == (nz,):
        result = np.broadcast_to(array[:, None], (nz, ny)).copy()
    elif array.shape == (nz, ny):
        result = array.copy()
    else:
        raise ValueError(
            f"{name} must be scalar, shape {(nz,)}, or shape {(nz, ny)}"
        )
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _face_density(
    density: np.ndarray, solid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nz, ny, nx = solid.shape
    x = np.zeros((nz, ny, nx + 1), dtype=float)
    y = np.zeros((nz, ny + 1, nx), dtype=float)
    z = np.zeros((nz + 1, ny, nx), dtype=float)
    x[:, :, 0], x[:, :, -1] = density[:, :, 0], density[:, :, -1]
    y[:, 0, :], y[:, -1, :] = density[:, 0, :], density[:, -1, :]
    z[0, :, :], z[-1, :, :] = density[0, :, :], density[-1, :, :]
    x[:, :, 1:-1] = 0.5 * (density[:, :, :-1] + density[:, :, 1:])
    y[:, 1:-1, :] = 0.5 * (density[:, :-1, :] + density[:, 1:, :])
    z[1:-1, :, :] = 0.5 * (density[:-1, :, :] + density[1:, :, :])
    x[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    y[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    z[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    x[:, :, 0][solid[:, :, 0]] = x[:, :, -1][solid[:, :, -1]] = 0.0
    y[:, 0, :][solid[:, 0, :]] = y[:, -1, :][solid[:, -1, :]] = 0.0
    z[0, :, :][solid[0, :, :]] = z[-1, :, :][solid[-1, :, :]] = 0.0
    return x, y, z


def _increment_faces_from_cell_delta(
    *,
    previous_u_face: np.ndarray,
    previous_v_face: np.ndarray,
    previous_w_face: np.ndarray,
    delta_u: np.ndarray,
    delta_v: np.ndarray,
    delta_w: np.ndarray,
    solid: np.ndarray,
    inlet_u: np.ndarray,
    nonstreamwise_wall_boundary: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Add the predictor delta without repeatedly filtering the face field."""

    u_face = previous_u_face.copy()
    v_face = previous_v_face.copy()
    w_face = previous_w_face.copy()
    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:])
    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :])
    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :])
    increment = 0.5 * (delta_u[:, :, :-1] + delta_u[:, :, 1:])
    u_face[:, :, 1:-1][fluid_x] += increment[fluid_x]
    increment = 0.5 * (delta_v[:, :-1, :] + delta_v[:, 1:, :])
    v_face[:, 1:-1, :][fluid_y] += increment[fluid_y]
    increment = 0.5 * (delta_w[:-1, :, :] + delta_w[1:, :, :])
    w_face[1:-1, :, :][fluid_z] += increment[fluid_z]
    outlet_fluid = ~solid[:, :, -1]
    u_face[:, :, -1][outlet_fluid] += delta_u[:, :, -1][outlet_fluid]
    if nonstreamwise_wall_boundary == "all_pressure_open":
        v_face[:, 0, :][~solid[:, 0, :]] += delta_v[:, 0, :][~solid[:, 0, :]]
        v_face[:, -1, :][~solid[:, -1, :]] += delta_v[:, -1, :][~solid[:, -1, :]]
        w_face[0, :, :][~solid[0, :, :]] += delta_w[0, :, :][~solid[0, :, :]]
        w_face[-1, :, :][~solid[-1, :, :]] += delta_w[-1, :, :][~solid[-1, :, :]]
    u_face[:, :, 1:-1][~fluid_x] = 0.0
    v_face[:, 1:-1, :][~fluid_y] = 0.0
    w_face[1:-1, :, :][~fluid_z] = 0.0
    u_face[:, :, 0] = inlet_u
    u_face[:, :, 0][solid[:, :, 0]] = 0.0
    u_face[:, :, -1][solid[:, :, -1]] = 0.0
    if nonstreamwise_wall_boundary != "all_pressure_open":
        v_face[:, 0, :] = 0.0
        v_face[:, -1, :] = 0.0
        w_face[0, :, :] = 0.0
        w_face[-1, :, :] = 0.0
    else:
        v_face[:, 0, :][solid[:, 0, :]] = 0.0
        v_face[:, -1, :][solid[:, -1, :]] = 0.0
        w_face[0, :, :][solid[0, :, :]] = 0.0
        w_face[-1, :, :][solid[-1, :, :]] = 0.0
    return u_face, v_face, w_face


def _cell_velocity_from_faces(
    u_face: np.ndarray,
    v_face: np.ndarray,
    w_face: np.ndarray,
    solid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = 0.5 * (u_face[:, :, :-1] + u_face[:, :, 1:])
    v = 0.5 * (v_face[:, :-1, :] + v_face[:, 1:, :])
    w = 0.5 * (w_face[:-1, :, :] + w_face[1:, :, :])
    u[solid] = v[solid] = w[solid] = 0.0
    return u, v, w


def initialise_nested_outer_binary_flow_3d(
    *,
    cells: tuple[int, int, int],
    closure: BinaryIdealGasMixtureClosure,
    ambient_u_m_s: float | np.ndarray,
    ambient_v_m_s: float | np.ndarray = 0.0,
    ambient_w_m_s: float | np.ndarray = 0.0,
    solid_mask: np.ndarray | None = None,
) -> NestedOuterBinaryFlow3DState:
    """Create an ambient, divergence-free outer state.

    ``cells`` follows the public solver convention ``(nx, ny, nz)`` while
    arrays are stored as ``(nz, ny, nx)``.
    """

    if (
        not isinstance(cells, tuple) or len(cells) != 3
        or any(isinstance(value, bool) or int(value) != value or value < 1 for value in cells)
    ):
        raise ValueError("cells must contain three positive integers")
    if not isinstance(closure, BinaryIdealGasMixtureClosure):
        raise TypeError("closure must be BinaryIdealGasMixtureClosure")
    nx, ny, nz = (int(value) for value in cells)
    shape = (nz, ny, nx)
    solid = (
        np.zeros(shape, dtype=bool)
        if solid_mask is None else np.asarray(solid_mask)
    )
    if solid.dtype != np.bool_ or solid.shape != shape:
        raise ValueError(f"solid_mask must be boolean with shape {shape}")
    inlet_u = _plane(ambient_u_m_s, nz, ny, "ambient_u_m_s")
    inlet_v = _plane(ambient_v_m_s, nz, ny, "ambient_v_m_s")
    inlet_w = _plane(ambient_w_m_s, nz, ny, "ambient_w_m_s")
    u = np.broadcast_to(inlet_u[:, :, None], shape).copy()
    v = np.broadcast_to(inlet_v[:, :, None], shape).copy()
    w = np.broadcast_to(inlet_w[:, :, None], shape).copy()
    u[solid] = v[solid] = w[solid] = 0.0
    u_face = np.broadcast_to(inlet_u[:, :, None], (nz, ny, nx + 1)).copy()
    v_face = np.zeros((nz, ny + 1, nx), dtype=float)
    w_face = np.zeros((nz + 1, ny, nx), dtype=float)
    v_face[:, 1:-1, :] = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
    w_face[1:-1, :, :] = 0.5 * (w[:-1, :, :] + w[1:, :, :])
    u_face[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    v_face[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    w_face[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    u_face[:, :, 0][solid[:, :, 0]] = 0.0
    u_face[:, :, -1][solid[:, :, -1]] = 0.0
    species = np.zeros(shape, dtype=float)
    enthalpy = np.zeros(shape, dtype=float)
    state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        species, enthalpy,
    )
    return NestedOuterBinaryFlow3DState(
        schema=NESTED_OUTER_BINARY_FLOW_3D_SCHEMA,
        time_s=0.0,
        species_mass_concentration_kg_m3=species,
        sensible_enthalpy_density_j_m3=enthalpy,
        thermodynamic_state=state,
        u_m_s=u,
        v_m_s=v,
        w_m_s=w,
        u_face_m_s=u_face,
        v_face_m_s=v_face,
        w_face_m_s=w_face,
        correction_potential_kg_m_s=None,
        cumulative_species_boundary_outflow_mass_kg=0.0,
        cumulative_sensible_enthalpy_boundary_outflow_j=0.0,
        cumulative_species_mass_residual_kg=0.0,
        cumulative_sensible_enthalpy_residual_j=0.0,
        cumulative_prescribed_inlet_species_mass_kg=0.0,
        cumulative_prescribed_inlet_sensible_enthalpy_j=0.0,
    )


def advance_nested_outer_binary_flow_step_3d(
    state: NestedOuterBinaryFlow3DState,
    *,
    handoff: ConservativeNestedHandoff3DResult | object,
    ambient_total_mass_flux_kg_m2_s: float | np.ndarray,
    closure: BinaryIdealGasMixtureClosure,
    solid_mask: np.ndarray,
    spacing_m: tuple[float, float, float] | np.ndarray,
    time_step_s: float,
    kinematic_viscosity_m2_s: float | np.ndarray,
    species_diffusivity_m2_s: float | np.ndarray,
    thermal_diffusivity_m2_s: float | np.ndarray,
    advection_scheme: str = "muscl_van_leer",
    include_buoyancy: bool = True,
    nonstreamwise_wall_boundary: str = "ground_no_slip_lateral_top_free_slip",
    open_boundary_faces: tuple[str, ...] = ("x_max",),
    projection_linear_solver: str = "auto",
    maximum_corrector_iterations: int = 32,
    minimum_corrector_iterations: int = 2,
    relative_density_tolerance: float = 1.0e-6,
    relative_mass_flux_tolerance: float = 1.0e-6,
) -> NestedOuterBinaryFlow3DStepResult:
    """Advance momentum, mixture transport and continuity by one common step."""

    if not isinstance(state, NestedOuterBinaryFlow3DState):
        raise TypeError("state must be NestedOuterBinaryFlow3DState")
    if not isinstance(closure, BinaryIdealGasMixtureClosure):
        raise TypeError("closure must be BinaryIdealGasMixtureClosure")
    solid = np.asarray(solid_mask)
    shape = state.species_mass_concentration_kg_m3.shape
    if solid.dtype != np.bool_ or solid.shape != shape:
        raise ValueError(f"solid_mask must be boolean with shape {shape}")
    nz, ny, nx = shape
    _finite_field(state.u_m_s, shape, "state.u_m_s")
    _finite_field(state.v_m_s, shape, "state.v_m_s")
    _finite_field(state.w_m_s, shape, "state.w_m_s")
    _finite_field(state.u_face_m_s, (nz, ny, nx + 1), "state.u_face_m_s")
    _finite_field(state.v_face_m_s, (nz, ny + 1, nx), "state.v_face_m_s")
    _finite_field(state.w_face_m_s, (nz + 1, ny, nx), "state.w_face_m_s")
    spacing_array = np.asarray(spacing_m, dtype=float)
    if (
        spacing_array.shape != (3,) or np.any(~np.isfinite(spacing_array))
        or np.any(spacing_array <= 0.0)
    ):
        raise ValueError("spacing_m must contain three positive finite values")
    dt = float(time_step_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_s must be positive and finite")
    query_time = state.time_s + 0.5 * dt
    if isinstance(handoff, ConservativeNestedHandoff3DResult):
        frame = nested_outer_inlet_frame_at(
            handoff,
            time_s=query_time,
            ambient_total_mass_flux_kg_m2_s=ambient_total_mass_flux_kg_m2_s,
            closure=closure,
        )
    elif callable(getattr(handoff, "frame_at", None)):
        frame = handoff.frame_at(
            time_s=query_time,
            ambient_total_mass_flux_kg_m2_s=ambient_total_mass_flux_kg_m2_s,
            closure=closure,
        )
        if not isinstance(frame, NestedOuterInletFrame3D):
            raise TypeError("handoff frame_at() must return NestedOuterInletFrame3D")
    else:
        raise TypeError(
            "handoff must be a ConservativeNestedHandoff3DResult or lazy provider"
        )
    if frame.total_mass_flux_kg_m2_s.shape != (nz, ny):
        raise ValueError(
            "handoff target plane must match the outer domain z-y inlet plane"
        )
    density = state.thermodynamic_state.mixture_density_kg_m3.copy()
    density[solid] = closure.ambient_density_kg_m3
    rho_x, rho_y, rho_z = _face_density(density, solid)
    mass_flux_x = rho_x * state.u_face_m_s
    mass_flux_y = rho_y * state.v_face_m_s
    mass_flux_z = rho_z * state.w_face_m_s
    mass_flux_x[:, :, 0] = frame.total_mass_flux_kg_m2_s
    mass_flux_x[:, :, 0][solid[:, :, 0]] = 0.0
    acceleration_z: float | np.ndarray = 0.0
    if include_buoyancy:
        acceleration_z = state.thermodynamic_state.buoyancy_acceleration_m_s2.copy()
        acceleration_z[solid] = 0.0
    momentum = advance_variable_density_momentum_3d(
        density,
        state.u_m_s,
        state.v_m_s,
        state.w_m_s,
        mass_flux_x_kg_m2_s=mass_flux_x,
        mass_flux_y_kg_m2_s=mass_flux_y,
        mass_flux_z_kg_m2_s=mass_flux_z,
        solid_mask=solid,
        spacing_m=tuple(float(value) for value in spacing_array),
        time_step_s=dt,
        kinematic_viscosity_m2_s=kinematic_viscosity_m2_s,
        acceleration_z_m_s2=acceleration_z,
        inlet_u_m_s=frame.velocity_u_m_s,
        inlet_v_m_s=frame.velocity_v_m_s,
        inlet_w_m_s=frame.velocity_w_m_s,
        nonstreamwise_wall_boundary=nonstreamwise_wall_boundary,
    )
    provisional_faces = _increment_faces_from_cell_delta(
        previous_u_face=state.u_face_m_s,
        previous_v_face=state.v_face_m_s,
        previous_w_face=state.w_face_m_s,
        delta_u=momentum.provisional_u_m_s - state.u_m_s,
        delta_v=momentum.provisional_v_m_s - state.v_m_s,
        delta_w=momentum.provisional_w_m_s - state.w_m_s,
        solid=solid,
        inlet_u=frame.velocity_u_m_s,
        nonstreamwise_wall_boundary=nonstreamwise_wall_boundary,
    )
    corrector = solve_binary_mixture_initial_cloud_step_3d(
        state.species_mass_concentration_kg_m3,
        state.sensible_enthalpy_density_j_m3,
        provisional_u_face_m_s=provisional_faces[0],
        provisional_v_face_m_s=provisional_faces[1],
        provisional_w_face_m_s=provisional_faces[2],
        solid_mask=solid,
        spacing_m=tuple(float(value) for value in spacing_array),
        time_step_s=dt,
        species_diffusivity_m2_s=species_diffusivity_m2_s,
        thermal_diffusivity_m2_s=thermal_diffusivity_m2_s,
        closure=closure,
        advection_scheme=advection_scheme,
        maximum_iterations=maximum_corrector_iterations,
        minimum_iterations=minimum_corrector_iterations,
        relative_density_tolerance=relative_density_tolerance,
        relative_mass_flux_tolerance=relative_mass_flux_tolerance,
        open_boundary_faces=open_boundary_faces,
        projection_linear_solver=projection_linear_solver,
        initial_correction_potential_kg_m_s=state.correction_potential_kg_m_s,
        inlet_total_mass_flux_kg_m2_s=frame.total_mass_flux_kg_m2_s,
        inlet_species_mass_flux_kg_m2_s=frame.species_mass_flux_kg_m2_s,
        inlet_sensible_enthalpy_flux_w_m2=frame.sensible_enthalpy_flux_w_m2,
    )
    u, v, w = _cell_velocity_from_faces(
        corrector.u_face_m_s,
        corrector.v_face_m_s,
        corrector.w_face_m_s,
        solid,
    )
    next_state = NestedOuterBinaryFlow3DState(
        schema=NESTED_OUTER_BINARY_FLOW_3D_SCHEMA,
        time_s=state.time_s + dt,
        species_mass_concentration_kg_m3=(
            corrector.species_mass_concentration_kg_m3
        ),
        sensible_enthalpy_density_j_m3=(
            corrector.sensible_enthalpy_density_j_m3
        ),
        thermodynamic_state=corrector.thermodynamic_state,
        u_m_s=u,
        v_m_s=v,
        w_m_s=w,
        u_face_m_s=corrector.u_face_m_s,
        v_face_m_s=corrector.v_face_m_s,
        w_face_m_s=corrector.w_face_m_s,
        correction_potential_kg_m_s=(
            corrector.mass_projection.correction_potential_kg_m_s.copy()
        ),
        cumulative_species_boundary_outflow_mass_kg=(
            state.cumulative_species_boundary_outflow_mass_kg
            + corrector.species_boundary_outflow_mass_kg
        ),
        cumulative_sensible_enthalpy_boundary_outflow_j=(
            state.cumulative_sensible_enthalpy_boundary_outflow_j
            + corrector.sensible_enthalpy_boundary_outflow_j
        ),
        cumulative_species_mass_residual_kg=(
            state.cumulative_species_mass_residual_kg
            + corrector.species_mass_residual_kg
        ),
        cumulative_sensible_enthalpy_residual_j=(
            state.cumulative_sensible_enthalpy_residual_j
            + corrector.sensible_enthalpy_residual_j
        ),
        cumulative_prescribed_inlet_species_mass_kg=(
            state.cumulative_prescribed_inlet_species_mass_kg
            + corrector.prescribed_inlet_species_mass_kg
        ),
        cumulative_prescribed_inlet_sensible_enthalpy_j=(
            state.cumulative_prescribed_inlet_sensible_enthalpy_j
            + corrector.prescribed_inlet_sensible_enthalpy_j
        ),
    )
    area = spacing_array[1] * spacing_array[2]
    inlet_momentum_flux = np.asarray((
        np.sum(frame.total_mass_flux_kg_m2_s * frame.velocity_u_m_s) * area,
        np.sum(frame.total_mass_flux_kg_m2_s * frame.velocity_v_m_s) * area,
        np.sum(frame.total_mass_flux_kg_m2_s * frame.velocity_w_m_s) * area,
    ))
    return NestedOuterBinaryFlow3DStepResult(
        schema=NESTED_OUTER_BINARY_FLOW_3D_SCHEMA,
        state=next_state,
        inlet_frame=frame,
        momentum=momentum,
        binary_corrector=corrector,
        diagnostics={
            "inlet_evaluation_time_s": query_time,
            "inlet_momentum_flux_n": inlet_momentum_flux.tolist(),
            "maximum_absolute_momentum_balance_residual_kg_m_s": float(
                np.max(np.abs(momentum.momentum_balance_residual_kg_m_s))
            ),
            "maximum_continuity_residual_kg_m3_s": (
                corrector.mass_projection
                .maximum_absolute_residual_after_kg_m3_s
            ),
            "corrector_converged": corrector.converged,
            "pressure_projection_fed_back_to_cell_velocity": True,
            "three_component_inlet_momentum_applied": True,
        },
    )


__all__ = [
    "NESTED_OUTER_BINARY_FLOW_3D_SCHEMA",
    "NestedOuterBinaryFlow3DState",
    "NestedOuterBinaryFlow3DStepResult",
    "initialise_nested_outer_binary_flow_3d",
    "advance_nested_outer_binary_flow_step_3d",
]
