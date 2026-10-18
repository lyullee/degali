"""Conservative binary-mixture/EOS/mass-flux corrector.

The transported species variable is ``rho*Y`` rather than a tracer added to a
fixed air inventory.  Advection uses the projected *total mass flux* and binary
diffusion uses ``-rho D grad(Y)``.  Sensible enthalpy is advected as a specific
quantity and conducted as ``-rho cp alpha grad(T)``.  The variable-density
projection then enforces total continuity at the EOS density.

An optional volumetric species/enthalpy source is coupled into the same-step
EOS density and total-mass projection. This supports both finite initial
clouds and source-driven binary-mixture releases without treating injected
species as a fixed-density passive dye.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
    BinaryIdealGasMixtureState,
)
from .variable_density_projection_3d import (
    VariableDensityMassFluxProjection3DResult,
    project_variable_density_mass_flux_3d,
)


BINARY_MIXTURE_CORRECTOR_3D_SCHEMA = "degali.binary-mixture-corrector-3d.v4"


@dataclass(frozen=True)
class BinaryMixtureCorrector3DResult:
    schema: str
    species_mass_concentration_kg_m3: np.ndarray
    sensible_enthalpy_density_j_m3: np.ndarray
    thermodynamic_state: BinaryIdealGasMixtureState
    u_face_m_s: np.ndarray
    v_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    mass_projection: VariableDensityMassFluxProjection3DResult
    iteration_history: np.ndarray
    projection_iterations_history: np.ndarray
    iterations: int
    converged: bool
    species_boundary_outflow_mass_kg: float
    sensible_enthalpy_boundary_outflow_j: float
    species_mass_residual_kg: float
    sensible_enthalpy_residual_j: float
    prescribed_inlet_species_mass_kg: float
    prescribed_inlet_sensible_enthalpy_j: float


def _finite_field(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    field = np.asarray(value, dtype=float)
    if field.shape != shape or not np.all(np.isfinite(field)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return field


def _diffusivity_field(
    value: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    shape: tuple[int, ...], name: str,
) -> np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray]:
    if isinstance(value, tuple):
        if len(value) != 3:
            raise ValueError(f"{name} directional value must contain Dx, Dy, and Dz")
        return tuple(_diffusivity_field(item, shape, name) for item in value)  # type: ignore[return-value]
    if np.ndim(value) == 0:
        scalar = float(value)
        if not math.isfinite(scalar) or scalar < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
        return np.full(shape, scalar, dtype=float)
    field = _finite_field(value, shape, name)
    if np.any(field < 0.0):
        raise ValueError(f"{name} must be non-negative")
    return field


def _face_average(
    cell: np.ndarray, solid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nz, ny, nx = solid.shape
    x = np.zeros((nz, ny, nx + 1), dtype=float)
    y = np.zeros((nz, ny + 1, nx), dtype=float)
    z = np.zeros((nz + 1, ny, nx), dtype=float)
    x[:, :, 0], x[:, :, -1] = cell[:, :, 0], cell[:, :, -1]
    y[:, 0, :], y[:, -1, :] = cell[:, 0, :], cell[:, -1, :]
    z[0, :, :], z[-1, :, :] = cell[0, :, :], cell[-1, :, :]
    x[:, :, 1:-1] = 0.5 * (cell[:, :, :-1] + cell[:, :, 1:])
    y[:, 1:-1, :] = 0.5 * (cell[:, :-1, :] + cell[:, 1:, :])
    z[1:-1, :, :] = 0.5 * (cell[:-1, :, :] + cell[1:, :, :])
    x[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    y[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    z[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    x[:, :, 0][solid[:, :, 0]] = x[:, :, -1][solid[:, :, -1]] = 0.0
    y[:, 0, :][solid[:, 0, :]] = y[:, -1, :][solid[:, -1, :]] = 0.0
    z[0, :, :][solid[0, :, :]] = z[-1, :, :][solid[-1, :, :]] = 0.0
    return x, y, z


def _velocity_from_mass_flux(
    mass_flux: np.ndarray, face_density: np.ndarray,
) -> np.ndarray:
    velocity = np.zeros_like(mass_flux)
    fluid = face_density > 0.0
    velocity[fluid] = mass_flux[fluid] / face_density[fluid]
    return velocity


def _limited_slope(
    state: np.ndarray, solid: np.ndarray, axis: int,
) -> np.ndarray:
    """Return a van-Leer slope with solid and domain-edge stencils disabled."""

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
    valid = same_sign & (np.abs(denominator) > 1.0e-300)
    local[valid] = (
        2.0 * backward[valid] * forward[valid] / denominator[valid]
    )
    valid_stencil = (
        (~solid[center_key])
        & (~solid[lower_key])
        & (~solid[upper_key])
    )
    local[~valid_stencil] = 0.0
    slope[center_key] = local
    return slope


def _reconstructed_states(
    state: np.ndarray,
    solid: np.ndarray,
    *,
    axis: int,
    face_velocity: np.ndarray,
    dt: float,
    spacing: float,
    advection_scheme: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return left/right interior-face states for upwind or MUSCL-Hancock."""

    lower = [slice(None)] * 3
    upper = [slice(None)] * 3
    lower[axis] = slice(None, -1)
    upper[axis] = slice(1, None)
    lower_key = tuple(lower)
    upper_key = tuple(upper)
    if advection_scheme == "upwind":
        return state[lower_key], state[upper_key]
    slope = _limited_slope(state, solid, axis)
    courant = face_velocity * dt / spacing
    left = state[lower_key] + 0.5 * (1.0 - courant) * slope[lower_key]
    right = state[upper_key] - 0.5 * (1.0 + courant) * slope[upper_key]
    return left, right


def _transport_binary_conserved_fields(
    species_mass_concentration: np.ndarray,
    sensible_enthalpy_density: np.ndarray,
    state: BinaryIdealGasMixtureState,
    *,
    mass_flux_x: np.ndarray,
    mass_flux_y: np.ndarray,
    mass_flux_z: np.ndarray,
    solid: np.ndarray,
    species_diffusivity: np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray],
    thermal_diffusivity: np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray],
    spacing: tuple[float, float, float],
    dt: float,
    open_boundary_faces: tuple[str, ...],
    advection_scheme: str,
    prescribed_x_min_species_flux_kg_m2_s: np.ndarray | None = None,
    prescribed_x_min_enthalpy_flux_w_m2: np.ndarray | None = None,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """One first-order bounded FV update driven by total face mass flux."""
    dx, dy, dz = spacing
    open_faces = frozenset(open_boundary_faces)
    rho = state.mixture_density_kg_m3
    species_y = state.species_mass_fraction
    specific_enthalpy = sensible_enthalpy_density / rho
    temperature = state.temperature_k
    heat_capacity = state.mixture_heat_capacity_j_kg_k
    rho_x, rho_y, rho_z = _face_average(rho, solid)
    def directional_faces(value):
        if isinstance(value, tuple):
            return tuple(_face_average(field, solid)[axis] for axis, field in enumerate(value))
        return _face_average(value, solid)

    d_x, d_y, d_z = directional_faces(species_diffusivity)
    alpha_x, alpha_y, alpha_z = directional_faces(thermal_diffusivity)
    cp_x, cp_y, cp_z = _face_average(heat_capacity, solid)

    species_flux_x = np.zeros_like(mass_flux_x)
    energy_flux_x = np.zeros_like(mass_flux_x)
    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:])
    mdot = mass_flux_x[:, :, 1:-1]
    velocity = _velocity_from_mass_flux(mdot, rho_x[:, :, 1:-1])
    left_species, right_species = _reconstructed_states(
        species_y, solid, axis=2, face_velocity=velocity,
        dt=dt, spacing=dx, advection_scheme=advection_scheme,
    )
    left_enthalpy, right_enthalpy = _reconstructed_states(
        specific_enthalpy, solid, axis=2, face_velocity=velocity,
        dt=dt, spacing=dx, advection_scheme=advection_scheme,
    )
    adv_species = mdot * np.where(
        mdot >= 0.0, left_species, right_species,
    )
    diff_species = -rho_x[:, :, 1:-1] * d_x[:, :, 1:-1] * (
        species_y[:, :, 1:] - species_y[:, :, :-1]
    ) / dx
    adv_energy = mdot * np.where(
        mdot >= 0.0, left_enthalpy, right_enthalpy,
    )
    conductive = -(
        rho_x[:, :, 1:-1] * cp_x[:, :, 1:-1] * alpha_x[:, :, 1:-1]
        * (temperature[:, :, 1:] - temperature[:, :, :-1]) / dx
    )
    species_flux_x[:, :, 1:-1][fluid_x] = (
        adv_species + diff_species
    )[fluid_x]
    energy_flux_x[:, :, 1:-1][fluid_x] = (adv_energy + conductive)[fluid_x]
    if prescribed_x_min_species_flux_kg_m2_s is not None:
        species_flux_x[:, :, 0] = prescribed_x_min_species_flux_kg_m2_s
        energy_flux_x[:, :, 0] = prescribed_x_min_enthalpy_flux_w_m2
    elif "x_min" in open_faces:
        species_flux_x[:, :, 0] = np.minimum(mass_flux_x[:, :, 0], 0.0) * species_y[:, :, 0]
        energy_flux_x[:, :, 0] = np.minimum(mass_flux_x[:, :, 0], 0.0) * specific_enthalpy[:, :, 0]
    if "x_max" in open_faces:
        species_flux_x[:, :, -1] = np.maximum(mass_flux_x[:, :, -1], 0.0) * species_y[:, :, -1]
        energy_flux_x[:, :, -1] = np.maximum(mass_flux_x[:, :, -1], 0.0) * specific_enthalpy[:, :, -1]

    species_flux_y = np.zeros_like(mass_flux_y)
    energy_flux_y = np.zeros_like(mass_flux_y)
    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :])
    mdot = mass_flux_y[:, 1:-1, :]
    velocity = _velocity_from_mass_flux(mdot, rho_y[:, 1:-1, :])
    left_species, right_species = _reconstructed_states(
        species_y, solid, axis=1, face_velocity=velocity,
        dt=dt, spacing=dy, advection_scheme=advection_scheme,
    )
    left_enthalpy, right_enthalpy = _reconstructed_states(
        specific_enthalpy, solid, axis=1, face_velocity=velocity,
        dt=dt, spacing=dy, advection_scheme=advection_scheme,
    )
    adv_species = mdot * np.where(
        mdot >= 0.0, left_species, right_species,
    )
    diff_species = -rho_y[:, 1:-1, :] * d_y[:, 1:-1, :] * (
        species_y[:, 1:, :] - species_y[:, :-1, :]
    ) / dy
    adv_energy = mdot * np.where(
        mdot >= 0.0, left_enthalpy, right_enthalpy,
    )
    conductive = -(
        rho_y[:, 1:-1, :] * cp_y[:, 1:-1, :] * alpha_y[:, 1:-1, :]
        * (temperature[:, 1:, :] - temperature[:, :-1, :]) / dy
    )
    species_flux_y[:, 1:-1, :][fluid_y] = (
        adv_species + diff_species
    )[fluid_y]
    energy_flux_y[:, 1:-1, :][fluid_y] = (adv_energy + conductive)[fluid_y]
    if "y_min" in open_faces:
        species_flux_y[:, 0, :] = np.minimum(mass_flux_y[:, 0, :], 0.0) * species_y[:, 0, :]
        energy_flux_y[:, 0, :] = np.minimum(mass_flux_y[:, 0, :], 0.0) * specific_enthalpy[:, 0, :]
    if "y_max" in open_faces:
        species_flux_y[:, -1, :] = np.maximum(mass_flux_y[:, -1, :], 0.0) * species_y[:, -1, :]
        energy_flux_y[:, -1, :] = np.maximum(mass_flux_y[:, -1, :], 0.0) * specific_enthalpy[:, -1, :]

    species_flux_z = np.zeros_like(mass_flux_z)
    energy_flux_z = np.zeros_like(mass_flux_z)
    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :])
    mdot = mass_flux_z[1:-1, :, :]
    velocity = _velocity_from_mass_flux(mdot, rho_z[1:-1, :, :])
    left_species, right_species = _reconstructed_states(
        species_y, solid, axis=0, face_velocity=velocity,
        dt=dt, spacing=dz, advection_scheme=advection_scheme,
    )
    left_enthalpy, right_enthalpy = _reconstructed_states(
        specific_enthalpy, solid, axis=0, face_velocity=velocity,
        dt=dt, spacing=dz, advection_scheme=advection_scheme,
    )
    adv_species = mdot * np.where(
        mdot >= 0.0, left_species, right_species,
    )
    diff_species = -rho_z[1:-1, :, :] * d_z[1:-1, :, :] * (
        species_y[1:, :, :] - species_y[:-1, :, :]
    ) / dz
    adv_energy = mdot * np.where(
        mdot >= 0.0, left_enthalpy, right_enthalpy,
    )
    conductive = -(
        rho_z[1:-1, :, :] * cp_z[1:-1, :, :] * alpha_z[1:-1, :, :]
        * (temperature[1:, :, :] - temperature[:-1, :, :]) / dz
    )
    species_flux_z[1:-1, :, :][fluid_z] = (
        adv_species + diff_species
    )[fluid_z]
    energy_flux_z[1:-1, :, :][fluid_z] = (adv_energy + conductive)[fluid_z]
    if "z_min" in open_faces:
        species_flux_z[0, :, :] = np.minimum(mass_flux_z[0, :, :], 0.0) * species_y[0, :, :]
        energy_flux_z[0, :, :] = np.minimum(mass_flux_z[0, :, :], 0.0) * specific_enthalpy[0, :, :]
    if "z_max" in open_faces:
        species_flux_z[-1, :, :] = np.maximum(mass_flux_z[-1, :, :], 0.0) * species_y[-1, :, :]
        energy_flux_z[-1, :, :] = np.maximum(mass_flux_z[-1, :, :], 0.0) * specific_enthalpy[-1, :, :]

    if blocked_faces is not None:
        blocked_x, blocked_y, blocked_z = blocked_faces
        species_flux_x[blocked_x] = energy_flux_x[blocked_x] = 0.0
        species_flux_y[blocked_y] = energy_flux_y[blocked_y] = 0.0
        species_flux_z[blocked_z] = energy_flux_z[blocked_z] = 0.0

    species_tendency = -(
        (species_flux_x[:, :, 1:] - species_flux_x[:, :, :-1]) / dx
        + (species_flux_y[:, 1:, :] - species_flux_y[:, :-1, :]) / dy
        + (species_flux_z[1:, :, :] - species_flux_z[:-1, :, :]) / dz
    )
    energy_tendency = -(
        (energy_flux_x[:, :, 1:] - energy_flux_x[:, :, :-1]) / dx
        + (energy_flux_y[:, 1:, :] - energy_flux_y[:, :-1, :]) / dy
        + (energy_flux_z[1:, :, :] - energy_flux_z[:-1, :, :]) / dz
    )
    species_tendency[solid] = energy_tendency[solid] = 0.0
    new_species = species_mass_concentration + dt * species_tendency
    new_energy = sensible_enthalpy_density + dt * energy_tendency
    new_species[solid] = new_energy[solid] = 0.0

    species_outflow_rate = (
        (np.sum(species_flux_x[:, :, -1]) - np.sum(species_flux_x[:, :, 0])) * dy * dz
        + (np.sum(species_flux_y[:, -1, :]) - np.sum(species_flux_y[:, 0, :])) * dx * dz
        + (np.sum(species_flux_z[-1, :, :]) - np.sum(species_flux_z[0, :, :])) * dx * dy
    )
    energy_outflow_rate = (
        (np.sum(energy_flux_x[:, :, -1]) - np.sum(energy_flux_x[:, :, 0])) * dy * dz
        + (np.sum(energy_flux_y[:, -1, :]) - np.sum(energy_flux_y[:, 0, :])) * dx * dz
        + (np.sum(energy_flux_z[-1, :, :]) - np.sum(energy_flux_z[0, :, :])) * dx * dy
    )
    return (
        new_species,
        new_energy,
        float(dt * species_outflow_rate),
        float(dt * energy_outflow_rate),
    )


def solve_binary_mixture_initial_cloud_step_3d(
    species_mass_concentration_old_kg_m3: np.ndarray,
    sensible_enthalpy_density_old_j_m3: np.ndarray,
    *,
    provisional_u_face_m_s: np.ndarray,
    provisional_v_face_m_s: np.ndarray,
    provisional_w_face_m_s: np.ndarray,
    solid_mask: np.ndarray,
    spacing_m: tuple[float, float, float] | np.ndarray,
    time_step_s: float,
    species_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    thermal_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    closure: BinaryIdealGasMixtureClosure,
    advection_scheme: str = "upwind",
    source_time_integration: str = "prestep",
    maximum_iterations: int = 32,
    minimum_iterations: int = 2,
    relative_density_tolerance: float = 1.0e-6,
    relative_mass_flux_tolerance: float = 1.0e-6,
    open_boundary_faces: tuple[str, ...] = ("x_max",),
    projection_linear_solver: str = "auto",
    initial_correction_potential_kg_m_s: np.ndarray | None = None,
    species_source_increment_kg_m3: np.ndarray | float = 0.0,
    total_mass_source_increment_kg_m3: np.ndarray | float | None = None,
    sensible_enthalpy_source_increment_j_m3: np.ndarray | float = 0.0,
    inlet_total_mass_flux_kg_m2_s: np.ndarray | None = None,
    inlet_species_mass_flux_kg_m2_s: np.ndarray | None = None,
    inlet_sensible_enthalpy_flux_w_m2: np.ndarray | None = None,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> BinaryMixtureCorrector3DResult:
    """Advance one dense-mixture step and enforce total continuity.

    The optional x-min flux triplet is the conservative outer-domain boundary
    contract.  All three fields must be supplied together.  The prescribed
    total mass flux is held fixed while the pressure projection corrects the
    downstream/open faces.
    """
    if not isinstance(closure, BinaryIdealGasMixtureClosure):
        raise TypeError("closure must be BinaryIdealGasMixtureClosure")
    if advection_scheme not in {
        "upwind", "muscl_van_leer", "muscl_ssprk2",
    }:
        raise ValueError(
            "advection_scheme must be 'upwind', 'muscl_van_leer' or "
            "'muscl_ssprk2' for the full binary-mixture corrector"
        )
    if source_time_integration not in {"prestep", "strang"}:
        raise ValueError("source_time_integration must be 'prestep' or 'strang'")
    solid = np.asarray(solid_mask)
    if solid.dtype != np.bool_ or solid.ndim != 3:
        raise ValueError("solid_mask must be a three-dimensional boolean array")
    shape = solid.shape
    nz, ny, nx = shape
    if blocked_faces is not None:
        expected = ((nz, ny, nx + 1), (nz, ny + 1, nx), (nz + 1, ny, nx))
        if not isinstance(blocked_faces, tuple) or len(blocked_faces) != 3:
            raise ValueError("blocked_faces must be a tuple of x, y, and z face masks")
        blocked_faces = tuple(np.asarray(mask) for mask in blocked_faces)
        if any(
            mask.dtype != np.bool_ or mask.shape != target
            for mask, target in zip(blocked_faces, expected)
        ):
            raise ValueError("blocked face masks must be boolean staggered face arrays")
    old_species = _finite_field(
        species_mass_concentration_old_kg_m3, shape,
        "species_mass_concentration_old_kg_m3",
    ).copy()
    old_energy = _finite_field(
        sensible_enthalpy_density_old_j_m3, shape,
        "sensible_enthalpy_density_old_j_m3",
    ).copy()
    if np.any(old_species[~solid] < 0.0):
        raise ValueError("species mass concentration must be non-negative")
    old_species[solid] = old_energy[solid] = 0.0
    species_source = np.broadcast_to(
        np.asarray(species_source_increment_kg_m3, dtype=float), shape,
    ).copy()
    energy_source = np.broadcast_to(
        np.asarray(sensible_enthalpy_source_increment_j_m3, dtype=float), shape,
    ).copy()
    total_mass_source = np.broadcast_to(
        np.asarray(
            species_source_increment_kg_m3
            if total_mass_source_increment_kg_m3 is None
            else total_mass_source_increment_kg_m3,
            dtype=float,
        ),
        shape,
    ).copy()
    if not np.all(np.isfinite(species_source)) or not np.all(np.isfinite(energy_source)):
        raise ValueError("binary-mixture source increments must be finite")
    if not np.all(np.isfinite(total_mass_source)):
        raise ValueError("total mass source increment must be finite")
    if np.any(species_source < 0.0):
        raise ValueError("species source increment must be non-negative")
    if np.any(total_mass_source < species_source):
        raise ValueError("total mass source must be at least the tracked species source")
    if np.any(np.abs(species_source[solid]) > 0.0) or np.any(
        np.abs(energy_source[solid]) > 0.0
    ) or np.any(
        np.abs(total_mass_source[solid]) > 0.0
    ):
        raise ValueError("binary-mixture source cannot be applied inside solid cells")
    u_base = _finite_field(
        provisional_u_face_m_s, (nz, ny, nx + 1), "provisional_u_face_m_s",
    )
    v_base = _finite_field(
        provisional_v_face_m_s, (nz, ny + 1, nx), "provisional_v_face_m_s",
    )
    w_base = _finite_field(
        provisional_w_face_m_s, (nz + 1, ny, nx), "provisional_w_face_m_s",
    )
    spacing_array = np.asarray(spacing_m, dtype=float)
    if spacing_array.shape != (3,) or np.any(~np.isfinite(spacing_array)) or np.any(spacing_array <= 0.0):
        raise ValueError("spacing_m must contain three positive finite values")
    spacing = tuple(float(value) for value in spacing_array)
    cell_volume = math.prod(spacing)
    injected_species_mass = float(np.sum(species_source[~solid]) * cell_volume)
    injected_energy = float(np.sum(energy_source[~solid]) * cell_volume)
    dt = float(time_step_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_s must be positive and finite")
    for value, name in (
        (relative_density_tolerance, "relative_density_tolerance"),
        (relative_mass_flux_tolerance, "relative_mass_flux_tolerance"),
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
    species_diffusivity = _diffusivity_field(
        species_diffusivity_m2_s, shape, "species_diffusivity_m2_s",
    )
    thermal_diffusivity = _diffusivity_field(
        thermal_diffusivity_m2_s, shape, "thermal_diffusivity_m2_s",
    )
    inlet_values = (
        inlet_total_mass_flux_kg_m2_s,
        inlet_species_mass_flux_kg_m2_s,
        inlet_sensible_enthalpy_flux_w_m2,
    )
    if any(value is not None for value in inlet_values) and not all(
        value is not None for value in inlet_values
    ):
        raise ValueError("all three binary-mixture inlet flux fields are required")
    inlet_total = inlet_species = inlet_enthalpy = None
    if all(value is not None for value in inlet_values):
        if "x_min" in open_boundary_faces:
            raise ValueError(
                "prescribed x-min mass flux requires pressure correction on other faces"
            )
        inlet_shape = (nz, ny)
        inlet_total = _finite_field(
            inlet_total_mass_flux_kg_m2_s, inlet_shape,
            "inlet_total_mass_flux_kg_m2_s",
        ).copy()
        inlet_species = _finite_field(
            inlet_species_mass_flux_kg_m2_s, inlet_shape,
            "inlet_species_mass_flux_kg_m2_s",
        ).copy()
        inlet_enthalpy = _finite_field(
            inlet_sensible_enthalpy_flux_w_m2, inlet_shape,
            "inlet_sensible_enthalpy_flux_w_m2",
        ).copy()
        inlet_fluid = ~solid[:, :, 0]
        if np.any(inlet_total[inlet_fluid] <= 0.0):
            raise ValueError("fluid x-min total mass flux must be positive")
        if np.any(inlet_species[inlet_fluid] < 0.0) or np.any(
            inlet_species[inlet_fluid] > inlet_total[inlet_fluid]
        ):
            raise ValueError(
                "x-min species flux must lie between zero and total mass flux"
            )
        if any(
            np.any(np.abs(field[~inlet_fluid]) > 0.0)
            for field in (inlet_total, inlet_species, inlet_enthalpy)
        ):
            raise ValueError("prescribed inlet flux cannot penetrate a solid cell")

    old_state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        old_species, old_energy,
    )
    source_pre_fraction = 0.5 if source_time_integration == "strang" else 1.0
    transport_species = old_species + source_pre_fraction * species_source
    transport_energy = old_energy + source_pre_fraction * energy_source
    transport_species[solid] = transport_energy[solid] = 0.0
    transport_state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        transport_species, transport_energy,
    )
    old_density = old_state.mixture_density_kg_m3.copy()
    old_density[solid] = closure.ambient_density_kg_m3
    rho_x, rho_y, rho_z = _face_average(old_density, solid)
    mass_flux_x = rho_x * u_base
    mass_flux_y = rho_y * v_base
    mass_flux_z = rho_z * w_base
    if blocked_faces is not None:
        mass_flux_x[blocked_faces[0]] = 0.0
        mass_flux_y[blocked_faces[1]] = 0.0
        mass_flux_z[blocked_faces[2]] = 0.0
    if inlet_total is not None:
        mass_flux_x[:, :, 0] = inlet_total
    density_guess = old_density.copy()
    flux_scale = max(
        float(np.max(np.abs(mass_flux_x))),
        float(np.max(np.abs(mass_flux_y))),
        float(np.max(np.abs(mass_flux_z))),
        # A zero provisional velocity can still acquire a small EOS-expansion
        # flux.  Normalise convergence by one ambient cell mass per step, not
        # by an arbitrary near-zero velocity scale.
        closure.ambient_density_kg_m3 * min(spacing) / dt,
    )
    history: list[tuple[float, float, float]] = []
    projection_iterations: list[int] = []
    projection_guess = None
    if initial_correction_potential_kg_m_s is not None:
        projection_guess = _finite_field(
            initial_correction_potential_kg_m_s,
            shape,
            "initial_correction_potential_kg_m_s",
        ).copy()
    converged = False
    projection = None
    new_species = old_species.copy()
    new_energy = old_energy.copy()
    state = old_state
    species_outflow = energy_outflow = 0.0
    for iteration in range(1, maximum_iterations + 1):
        new_species, new_energy, species_outflow, energy_outflow = (
            _transport_binary_conserved_fields(
                transport_species, transport_energy, transport_state,
                mass_flux_x=mass_flux_x,
                mass_flux_y=mass_flux_y,
                mass_flux_z=mass_flux_z,
                solid=solid,
                species_diffusivity=species_diffusivity,
                thermal_diffusivity=thermal_diffusivity,
                spacing=spacing,
                dt=dt,
                open_boundary_faces=open_boundary_faces,
                advection_scheme=advection_scheme,
                prescribed_x_min_species_flux_kg_m2_s=inlet_species,
                prescribed_x_min_enthalpy_flux_w_m2=inlet_enthalpy,
                blocked_faces=blocked_faces,
            )
        )
        minimum_species = float(np.min(new_species[~solid]))
        if minimum_species < -1.0e-12:
            raise RuntimeError(
                "binary-mixture step produced negative species concentration: "
                f"{minimum_species:.6g} kg/m3"
            )
        source_post_fraction = 1.0 - source_pre_fraction
        new_species = np.maximum(
            new_species + source_post_fraction * species_source, 0.0,
        )
        new_energy = new_energy + source_post_fraction * energy_source
        state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
            new_species, new_energy,
        )
        density_new = state.mixture_density_kg_m3.copy()
        density_new[solid] = closure.ambient_density_kg_m3
        new_rho_x, new_rho_y, new_rho_z = _face_average(density_new, solid)
        provisional_mass_flux_x = new_rho_x * u_base
        if inlet_total is not None:
            provisional_mass_flux_x[:, :, 0] = inlet_total
        projection = project_variable_density_mass_flux_3d(
            provisional_mass_flux_x,
            new_rho_y * v_base,
            new_rho_z * w_base,
            density_old_kg_m3=old_density,
            density_new_kg_m3=density_new,
            mass_source_rate_kg_m3_s=total_mass_source / dt,
            solid_mask=solid,
            spacing_m=spacing,
            time_step_s=dt,
            open_boundary_faces=open_boundary_faces,
            linear_solver=projection_linear_solver,
            initial_correction_potential_kg_m_s=projection_guess,
            blocked_faces=blocked_faces,
        )
        projection_iterations.append(projection.solver_iterations)
        projection_guess = projection.correction_potential_kg_m_s.copy()
        next_x = projection.mass_flux_x_kg_m2_s
        next_y = projection.mass_flux_y_kg_m2_s
        next_z = projection.mass_flux_z_kg_m2_s
        density_change = float(np.max(
            np.abs(density_new[~solid] - density_guess[~solid])
            / np.maximum(np.abs(density_new[~solid]), 1.0e-12)
        ))
        flux_change = max(
            float(np.max(np.abs(next_x - mass_flux_x))),
            float(np.max(np.abs(next_y - mass_flux_y))),
            float(np.max(np.abs(next_z - mass_flux_z))),
        ) / flux_scale
        history.append((
            density_change,
            flux_change,
            projection.maximum_absolute_residual_after_kg_m3_s,
        ))
        density_guess = density_new
        mass_flux_x, mass_flux_y, mass_flux_z = next_x, next_y, next_z
        if (
            iteration >= minimum_iterations
            and density_change <= relative_density_tolerance
            and flux_change <= relative_mass_flux_tolerance
        ):
            converged = True
            break

    assert projection is not None
    final_rho_x, final_rho_y, final_rho_z = _face_average(
        state.mixture_density_kg_m3, solid,
    )
    u_face = _velocity_from_mass_flux(mass_flux_x, final_rho_x)
    v_face = _velocity_from_mass_flux(mass_flux_y, final_rho_y)
    w_face = _velocity_from_mass_flux(mass_flux_z, final_rho_z)
    fluid = ~solid
    species_residual = float(
        np.sum((new_species - old_species)[fluid]) * cell_volume
        + species_outflow - injected_species_mass
    )
    energy_residual = float(
        np.sum((new_energy - old_energy)[fluid]) * cell_volume
        + energy_outflow - injected_energy
    )
    return BinaryMixtureCorrector3DResult(
        schema=BINARY_MIXTURE_CORRECTOR_3D_SCHEMA,
        species_mass_concentration_kg_m3=new_species,
        sensible_enthalpy_density_j_m3=new_energy,
        thermodynamic_state=state,
        u_face_m_s=u_face,
        v_face_m_s=v_face,
        w_face_m_s=w_face,
        mass_projection=projection,
        iteration_history=np.asarray(history),
        projection_iterations_history=np.asarray(projection_iterations, dtype=int),
        iterations=len(history),
        converged=converged,
        species_boundary_outflow_mass_kg=species_outflow,
        sensible_enthalpy_boundary_outflow_j=energy_outflow,
        species_mass_residual_kg=species_residual,
        sensible_enthalpy_residual_j=energy_residual,
        prescribed_inlet_species_mass_kg=(
            0.0 if inlet_species is None
            else float(np.sum(inlet_species) * spacing[1] * spacing[2] * dt)
        ),
        prescribed_inlet_sensible_enthalpy_j=(
            0.0 if inlet_enthalpy is None
            else float(np.sum(inlet_enthalpy) * spacing[1] * spacing[2] * dt)
        ),
    )


__all__ = [
    "BINARY_MIXTURE_CORRECTOR_3D_SCHEMA",
    "BinaryMixtureCorrector3DResult",
    "solve_binary_mixture_initial_cloud_step_3d",
]
