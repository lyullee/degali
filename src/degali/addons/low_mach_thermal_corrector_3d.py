"""Species/enthalpy/EOS/mass-flux fixed-point corrector for one 3-D step.

This operator couples two conservative scalar updates to the variable-density
continuity projection.  It accepts a provisional momentum velocity, predicts
hydrogen and sensible enthalpy, reconstructs EOS density, projects face mass
flux, and repeats scalar transport with the corrected volume velocity.

The current thermodynamic closure is the explicitly bounded active-scalar
H2/air approximation.  The algorithm and diagnostics are structured so that a
full low-Mach mixture closure can replace it without changing the discrete
mass constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .fast_obstacle_transport_3d import _transport_step
from .hydrogen_air_thermodynamics import (
    HydrogenAirActiveScalarClosure,
    HydrogenAirThermodynamicState,
)
from .variable_density_projection_3d import (
    VariableDensityMassFluxProjection3DResult,
    project_variable_density_mass_flux_3d,
)


LOW_MACH_THERMAL_CORRECTOR_3D_SCHEMA = (
    "degali.low-mach-thermal-corrector-3d.v4"
)


@dataclass(frozen=True)
class LowMachThermalCorrector3DResult:
    schema: str
    hydrogen_mass_concentration_kg_m3: np.ndarray
    sensible_enthalpy_density_j_m3: np.ndarray
    thermodynamic_state: HydrogenAirThermodynamicState
    u_face_m_s: np.ndarray
    v_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    pressure_correction_pa: np.ndarray
    mass_projection: VariableDensityMassFluxProjection3DResult
    iteration_history: np.ndarray
    projection_iterations_history: np.ndarray
    iterations: int
    converged: bool
    hydrogen_outflow_mass_kg: float
    sensible_enthalpy_outflow_j: float
    hydrogen_mass_residual_kg: float
    sensible_enthalpy_residual_j: float


def _finite_array(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


def _face_density(
    density: np.ndarray, solid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nz, ny, nx = solid.shape
    rho_x = np.zeros((nz, ny, nx + 1), dtype=float)
    rho_y = np.zeros((nz, ny + 1, nx), dtype=float)
    rho_z = np.zeros((nz + 1, ny, nx), dtype=float)
    rho_x[:, :, 0] = density[:, :, 0]
    rho_x[:, :, -1] = density[:, :, -1]
    rho_x[:, :, 1:-1] = 0.5 * (density[:, :, :-1] + density[:, :, 1:])
    rho_y[:, 0, :] = density[:, 0, :]
    rho_y[:, -1, :] = density[:, -1, :]
    rho_y[:, 1:-1, :] = 0.5 * (density[:, :-1, :] + density[:, 1:, :])
    rho_z[0, :, :] = density[0, :, :]
    rho_z[-1, :, :] = density[-1, :, :]
    rho_z[1:-1, :, :] = 0.5 * (density[:-1, :, :] + density[1:, :, :])
    rho_x[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    rho_y[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    rho_z[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    rho_x[:, :, 0][solid[:, :, 0]] = 0.0
    rho_x[:, :, -1][solid[:, :, -1]] = 0.0
    rho_y[:, 0, :][solid[:, 0, :]] = 0.0
    rho_y[:, -1, :][solid[:, -1, :]] = 0.0
    rho_z[0, :, :][solid[0, :, :]] = 0.0
    rho_z[-1, :, :][solid[-1, :, :]] = 0.0
    return rho_x, rho_y, rho_z


def _velocity_from_mass_flux(
    mass_flux: np.ndarray, face_density: np.ndarray,
) -> np.ndarray:
    velocity = np.zeros_like(mass_flux)
    positive = face_density > 0.0
    velocity[positive] = mass_flux[positive] / face_density[positive]
    return velocity


def solve_low_mach_thermal_corrector_3d(
    hydrogen_mass_concentration_old_kg_m3: np.ndarray,
    sensible_enthalpy_density_old_j_m3: np.ndarray,
    *,
    hydrogen_source_increment_kg_m3: np.ndarray,
    sensible_enthalpy_source_increment_j_m3: np.ndarray,
    provisional_u_face_m_s: np.ndarray,
    provisional_v_face_m_s: np.ndarray,
    provisional_w_face_m_s: np.ndarray,
    solid_mask: np.ndarray,
    spacing_m: tuple[float, float, float] | np.ndarray,
    time_step_s: float,
    species_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    thermal_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    closure: HydrogenAirActiveScalarClosure,
    advection_scheme: str = "upwind",
    source_time_integration: str = "prestep",
    maximum_iterations: int = 8,
    minimum_iterations: int = 2,
    relative_density_tolerance: float = 1.0e-6,
    relative_velocity_tolerance: float = 1.0e-6,
    inlet_sensible_enthalpy_density_j_m3: np.ndarray | None = None,
    open_boundary_faces: tuple[str, ...] = ("x_max",),
    initial_correction_potential_kg_m_s: np.ndarray | None = None,
    projection_linear_solver: str = "auto",
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    face_flux_callback=None,
    enthalpy_flux_callback=None,
) -> LowMachThermalCorrector3DResult:
    """Solve one scalar/EOS/continuity fixed point without advancing momentum."""
    if not isinstance(closure, HydrogenAirActiveScalarClosure):
        raise TypeError("closure must be HydrogenAirActiveScalarClosure")
    if advection_scheme not in {
        "upwind", "muscl_van_leer", "muscl_ssprk2",
    }:
        raise ValueError(
            "advection_scheme must be 'upwind', 'muscl_van_leer', or "
            "'muscl_ssprk2'"
        )
    if source_time_integration not in {"prestep", "strang"}:
        raise ValueError("source_time_integration must be 'prestep' or 'strang'")
    if face_flux_callback is not None and not callable(face_flux_callback):
        raise TypeError("face_flux_callback must be callable")
    if enthalpy_flux_callback is not None and not callable(enthalpy_flux_callback):
        raise TypeError("enthalpy_flux_callback must be callable")
    solid = np.asarray(solid_mask)
    if solid.dtype != np.bool_ or solid.ndim != 3:
        raise ValueError("solid_mask must be a three-dimensional boolean array")
    shape = solid.shape
    nz, ny, nx = shape
    old_h2 = _finite_array(
        hydrogen_mass_concentration_old_kg_m3, shape,
        "hydrogen_mass_concentration_old_kg_m3",
    )
    old_energy = _finite_array(
        sensible_enthalpy_density_old_j_m3, shape,
        "sensible_enthalpy_density_old_j_m3",
    )
    source_h2 = _finite_array(
        hydrogen_source_increment_kg_m3, shape,
        "hydrogen_source_increment_kg_m3",
    )
    source_energy = _finite_array(
        sensible_enthalpy_source_increment_j_m3, shape,
        "sensible_enthalpy_source_increment_j_m3",
    )
    inlet_energy = None
    if inlet_sensible_enthalpy_density_j_m3 is not None:
        inlet_energy = _finite_array(
            inlet_sensible_enthalpy_density_j_m3, (nz, ny),
            "inlet_sensible_enthalpy_density_j_m3",
        )
    if np.any(old_h2[~solid] < 0.0) or np.any(source_h2[~solid] < 0.0):
        raise ValueError("hydrogen state and source increment must be non-negative")
    for array, name, expected in (
        (provisional_u_face_m_s, "provisional_u_face_m_s", (nz, ny, nx + 1)),
        (provisional_v_face_m_s, "provisional_v_face_m_s", (nz, ny + 1, nx)),
        (provisional_w_face_m_s, "provisional_w_face_m_s", (nz + 1, ny, nx)),
    ):
        _finite_array(array, expected, name)
    u_base = np.asarray(provisional_u_face_m_s, dtype=float)
    v_base = np.asarray(provisional_v_face_m_s, dtype=float)
    w_base = np.asarray(provisional_w_face_m_s, dtype=float)
    spacing = np.asarray(spacing_m, dtype=float)
    if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("spacing_m must contain three positive finite values")
    dt = float(time_step_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_s must be positive and finite")
    diffusivities: dict[str, object] = {}
    for value, name in (
        (species_diffusivity_m2_s, "species_diffusivity_m2_s"),
        (thermal_diffusivity_m2_s, "thermal_diffusivity_m2_s"),
    ):
        if isinstance(value, tuple):
            if len(value) != 3:
                raise ValueError(f"{name} directional value must contain Dx, Dy, and Dz")
            fields = []
            for directional in value:
                if np.ndim(directional) == 0:
                    scalar = float(directional)
                    if not math.isfinite(scalar) or scalar < 0.0:
                        raise ValueError(f"{name} must be finite and non-negative")
                    fields.append(np.full(shape, scalar, dtype=float))
                else:
                    field = np.asarray(directional, dtype=float)
                    if field.shape != shape or not np.all(np.isfinite(field)) or np.any(field < 0.0):
                        raise ValueError(f"{name} directional fields must be finite, non-negative, and have shape {shape}")
                    fields.append(field)
            diffusivities[name] = tuple(fields)
            continue
        if np.ndim(value) == 0:
            scalar = float(value)
            if not math.isfinite(scalar) or scalar < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
            diffusivities[name] = scalar
        else:
            field = np.asarray(value, dtype=float)
            if (
                field.shape != shape
                or not np.all(np.isfinite(field))
                or np.any(field < 0.0)
            ):
                raise ValueError(
                    f"{name} must be finite and non-negative with shape {shape}"
                )
            diffusivities[name] = field
    for value, name in (
        (relative_density_tolerance, "relative_density_tolerance"),
        (relative_velocity_tolerance, "relative_velocity_tolerance"),
    ):
        if not math.isfinite(float(value)) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    for value, name in (
        (maximum_iterations, "maximum_iterations"),
        (minimum_iterations, "minimum_iterations"),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if minimum_iterations > maximum_iterations:
        raise ValueError("minimum_iterations cannot exceed maximum_iterations")

    old_h2 = old_h2.copy()
    old_energy = old_energy.copy()
    source_h2 = source_h2.copy()
    source_energy = source_energy.copy()
    old_h2[solid] = old_energy[solid] = 0.0
    source_h2[solid] = source_energy[solid] = 0.0
    source_pre_fraction = 0.5 if source_time_integration == "strang" else 1.0
    source_post_fraction = 1.0 - source_pre_fraction
    source_state_h2 = old_h2 + source_pre_fraction * source_h2
    source_state_energy = old_energy + source_pre_fraction * source_energy
    old_state = closure.state_from_conserved(old_h2, old_energy)
    density_old = old_state.mixture_density_kg_m3.copy()
    density_old[solid] = closure.ambient_air_density_kg_m3

    u_guess, v_guess, w_guess = u_base.copy(), v_base.copy(), w_base.copy()
    density_guess = density_old.copy()
    velocity_reference = max(
        float(np.max(np.abs(u_base))), float(np.max(np.abs(v_base))),
        float(np.max(np.abs(w_base))), 1.0e-8,
    )
    history: list[tuple[float, float, float]] = []
    projection_iterations_history: list[int] = []
    converged = False
    projection = None
    new_h2 = source_state_h2.copy()
    new_energy = source_state_energy.copy()
    state = closure.state_from_conserved(new_h2, new_energy)
    h2_outflow = energy_outflow = 0.0
    projection_guess = None
    if initial_correction_potential_kg_m_s is not None:
        projection_guess = _finite_array(
            initial_correction_potential_kg_m_s,
            shape,
            "initial_correction_potential_kg_m_s",
        ).copy()

    for iteration in range(1, maximum_iterations + 1):
        new_h2, h2_outflow = _transport_step(
            source_state_h2,
            u_face=u_guess, v_face=v_guess, w_face=w_guess,
            solid=solid, diffusivity=diffusivities["species_diffusivity_m2_s"],
            dt=dt, spacing=spacing,
            advection_scheme=advection_scheme,
            open_boundary_faces=open_boundary_faces,
            blocked_faces=blocked_faces,
            face_flux_callback=face_flux_callback,
        )
        if source_post_fraction:
            new_h2 += source_post_fraction * source_h2
        minimum_h2 = float(np.min(new_h2[~solid]))
        if minimum_h2 < -1.0e-12:
            raise RuntimeError(
                f"low-Mach corrector produced negative H2 concentration {minimum_h2:.6g}"
            )
        new_h2 = np.maximum(new_h2, 0.0)
        new_energy, energy_outflow = _transport_step(
            source_state_energy,
            u_face=u_guess, v_face=v_guess, w_face=w_guess,
            solid=solid, diffusivity=diffusivities["thermal_diffusivity_m2_s"],
            dt=dt, spacing=spacing,
            signed_boundary_flux=True,
            advection_scheme=advection_scheme,
            open_boundary_faces=open_boundary_faces,
            prescribed_x_min_inflow_state=inlet_energy,
            blocked_faces=blocked_faces,
            face_flux_callback=enthalpy_flux_callback,
        )
        if source_post_fraction:
            new_energy += source_post_fraction * source_energy
        state = closure.state_from_conserved(new_h2, new_energy)
        density_new = state.mixture_density_kg_m3.copy()
        density_new[solid] = closure.ambient_air_density_kg_m3
        rho_x, rho_y, rho_z = _face_density(density_new, solid)
        projection = project_variable_density_mass_flux_3d(
            rho_x * u_base, rho_y * v_base, rho_z * w_base,
            density_old_kg_m3=density_old,
            density_new_kg_m3=density_new,
            mass_source_rate_kg_m3_s=source_h2 / dt,
            solid_mask=solid,
            spacing_m=spacing,
            time_step_s=dt,
            open_boundary_faces=open_boundary_faces,
            initial_correction_potential_kg_m_s=projection_guess,
            linear_solver=projection_linear_solver,
            blocked_faces=blocked_faces,
        )
        projection_iterations_history.append(projection.solver_iterations)
        projection_guess = projection.correction_potential_kg_m_s
        u_new = _velocity_from_mass_flux(
            projection.mass_flux_x_kg_m2_s, rho_x,
        )
        v_new = _velocity_from_mass_flux(
            projection.mass_flux_y_kg_m2_s, rho_y,
        )
        w_new = _velocity_from_mass_flux(
            projection.mass_flux_z_kg_m2_s, rho_z,
        )
        density_change = float(np.max(
            np.abs(density_new[~solid] - density_guess[~solid])
            / np.maximum(np.abs(density_new[~solid]), 1.0e-12)
        ))
        velocity_change = max(
            float(np.max(np.abs(u_new - u_guess))),
            float(np.max(np.abs(v_new - v_guess))),
            float(np.max(np.abs(w_new - w_guess))),
        ) / velocity_reference
        continuity_residual = (
            projection.maximum_absolute_residual_after_kg_m3_s
        )
        history.append((density_change, velocity_change, continuity_residual))
        u_guess, v_guess, w_guess = u_new, v_new, w_new
        density_guess = density_new
        if (
            iteration >= minimum_iterations
            and density_change <= relative_density_tolerance
            and velocity_change <= relative_velocity_tolerance
        ):
            converged = True
            break

    assert projection is not None
    dx, dy, dz = (float(value) for value in spacing)
    cell_volume = dx * dy * dz
    fluid = ~solid
    hydrogen_residual = float(
        np.sum((new_h2 - old_h2)[fluid]) * cell_volume
        + h2_outflow
        - np.sum(source_h2[fluid]) * cell_volume
    )
    energy_residual = float(
        np.sum((new_energy - old_energy)[fluid]) * cell_volume
        + energy_outflow
        - np.sum(source_energy[fluid]) * cell_volume
    )
    return LowMachThermalCorrector3DResult(
        schema=LOW_MACH_THERMAL_CORRECTOR_3D_SCHEMA,
        hydrogen_mass_concentration_kg_m3=new_h2,
        sensible_enthalpy_density_j_m3=new_energy,
        thermodynamic_state=state,
        u_face_m_s=u_guess,
        v_face_m_s=v_guess,
        w_face_m_s=w_guess,
        pressure_correction_pa=projection.pressure_correction_pa,
        mass_projection=projection,
        iteration_history=np.asarray(history),
        projection_iterations_history=np.asarray(
            projection_iterations_history, dtype=int,
        ),
        iterations=len(history),
        converged=converged,
        hydrogen_outflow_mass_kg=h2_outflow,
        sensible_enthalpy_outflow_j=energy_outflow,
        hydrogen_mass_residual_kg=hydrogen_residual,
        sensible_enthalpy_residual_j=energy_residual,
    )


__all__ = [
    "LOW_MACH_THERMAL_CORRECTOR_3D_SCHEMA",
    "LowMachThermalCorrector3DResult",
    "solve_low_mach_thermal_corrector_3d",
]
