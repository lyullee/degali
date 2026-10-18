"""Extract a conservative binary-mixture handoff frame from an x face.

The plane operator uses the same mass-fraction and specific-enthalpy variables
as the binary-mixture corrector.  It records flux densities, not cell values,
so the result can be area-remapped without reconstructing a concentration from
an arbitrary velocity.  Interior species and heat diffusion are retained;
the two tangential momentum fluxes use the momentum predictor's upwind state.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .binary_mixture_thermodynamics import BinaryIdealGasMixtureState


BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA = (
    "degali.binary-mixture-handoff-plane-3d.v2"
)


@dataclass(frozen=True)
class BinaryMixtureHandoffPlane3DFrame:
    schema: str
    time_s: float
    x_face_index: int
    total_mass_flux_anomaly_kg_m2_s: np.ndarray
    advective_species_mass_flux_kg_m2_s: np.ndarray
    species_mass_flux_kg_m2_s: np.ndarray
    advective_sensible_enthalpy_flux_w_m2: np.ndarray
    sensible_enthalpy_flux_w_m2: np.ndarray
    lateral_momentum_flux_anomaly_pa: np.ndarray
    vertical_momentum_flux_anomaly_pa: np.ndarray
    diagnostics: dict[str, object]
    dense_layer_packet: dict[str, np.ndarray] | None = None


@dataclass(frozen=True)
class BinaryMixtureHandoffPlane3DRequest:
    """Grid-independent sampling request for a resolved-inner run."""

    x_face_m: float
    sample_interval_s: float
    include_initial: bool = True

    def __post_init__(self) -> None:
        for name in ("x_face_m", "sample_interval_s"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
            object.__setattr__(self, name, value)
        if not isinstance(self.include_initial, bool):
            raise TypeError("include_initial must be boolean")


def _field(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.shape != shape or np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return result


def _plane(value: float | np.ndarray, shape: tuple[int, int], name: str) -> np.ndarray:
    result = np.broadcast_to(np.asarray(value, dtype=float), shape).copy()
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be finite and broadcastable to {shape}")
    return result


def _limited_x_slope(state: np.ndarray, solid: np.ndarray) -> np.ndarray:
    slope = np.zeros_like(state)
    backward = state[:, :, 1:-1] - state[:, :, :-2]
    forward = state[:, :, 2:] - state[:, :, 1:-1]
    denominator = backward + forward
    valid = (
        backward * forward > 0.0
    ) & (np.abs(denominator) > 1.0e-300)
    valid &= (
        (~solid[:, :, :-2])
        & (~solid[:, :, 1:-1])
        & (~solid[:, :, 2:])
    )
    local = np.zeros_like(backward)
    local[valid] = (
        2.0 * backward[valid] * forward[valid] / denominator[valid]
    )
    slope[:, :, 1:-1] = local
    return slope


def extract_binary_mixture_handoff_plane_3d(
    state: BinaryIdealGasMixtureState,
    sensible_enthalpy_density_j_m3: np.ndarray,
    *,
    u_face_m_s: np.ndarray,
    cell_v_m_s: np.ndarray,
    cell_w_m_s: np.ndarray,
    solid_mask: np.ndarray,
    x_face_index: int,
    spacing_m: tuple[float, float, float] | np.ndarray,
    time_s: float,
    time_step_s: float,
    species_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    thermal_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    ambient_total_mass_flux_kg_m2_s: float | np.ndarray,
    ambient_v_m_s: float | np.ndarray = 0.0,
    ambient_w_m_s: float | np.ndarray = 0.0,
    advection_scheme: str = "muscl_van_leer",
) -> BinaryMixtureHandoffPlane3DFrame:
    """Return conserved flux densities on one solver x face.

    Boundary faces use the adjacent cell as a zero-gradient exterior state.
    Interior faces use the declared upwind or bounded MUSCL--Hancock scalar
    reconstruction.  ``muscl_ssprk2`` uses the same bounded MUSCL spatial
    reconstruction at the exported face; SSPRK2 is the time integrator of
    the transient transport step, not a separate face reconstruction.
    Tangential momentum follows the first-order upwind momentum predictor.
    """

    if not isinstance(state, BinaryIdealGasMixtureState):
        raise TypeError("state must be BinaryIdealGasMixtureState")
    solid = np.asarray(solid_mask)
    if solid.dtype != np.bool_ or solid.ndim != 3:
        raise ValueError("solid_mask must be a three-dimensional boolean array")
    shape = solid.shape
    nz, ny, nx = shape
    if (
        isinstance(x_face_index, bool)
        or not isinstance(x_face_index, int)
        or not 0 <= x_face_index <= nx
    ):
        raise ValueError(f"x_face_index must lie in [0, {nx}]")
    if advection_scheme not in {
        "upwind", "muscl_van_leer", "muscl_ssprk2",
    }:
        raise ValueError(
            "advection_scheme must be 'upwind', 'muscl_van_leer' or "
            "'muscl_ssprk2'"
        )
    density = _field(state.mixture_density_kg_m3, shape, "mixture density")
    mass_fraction = _field(state.species_mass_fraction, shape, "mass fraction")
    heat_capacity = _field(
        state.mixture_heat_capacity_j_kg_k, shape, "mixture heat capacity",
    )
    temperature = _field(state.temperature_k, shape, "temperature")
    enthalpy_density = _field(
        sensible_enthalpy_density_j_m3, shape, "sensible enthalpy density",
    )
    cell_v = _field(cell_v_m_s, shape, "cell_v_m_s")
    cell_w = _field(cell_w_m_s, shape, "cell_w_m_s")
    u_face = _field(u_face_m_s, (nz, ny, nx + 1), "u_face_m_s")
    spacing = np.asarray(spacing_m, dtype=float)
    if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("spacing_m must contain three positive finite values")
    dt = float(time_step_s)
    query_time = float(time_s)
    if not math.isfinite(dt) or dt <= 0.0 or not math.isfinite(query_time):
        raise ValueError("time and time_step_s must be finite; step must be positive")
    # Handoff plane is normal to x; its reported diffusive flux uses Dx.
    species_x_diffusivity = (
        species_diffusivity_m2_s[0]
        if isinstance(species_diffusivity_m2_s, tuple)
        else species_diffusivity_m2_s
    )
    thermal_x_diffusivity = (
        thermal_diffusivity_m2_s[0]
        if isinstance(thermal_diffusivity_m2_s, tuple)
        else thermal_diffusivity_m2_s
    )
    diffusivity = np.broadcast_to(
        np.asarray(species_x_diffusivity, dtype=float), shape,
    ).copy()
    thermal_diffusivity = np.broadcast_to(
        np.asarray(thermal_x_diffusivity, dtype=float), shape,
    ).copy()
    if (
        np.any(~np.isfinite(diffusivity)) or np.any(diffusivity < 0.0)
        or np.any(~np.isfinite(thermal_diffusivity))
        or np.any(thermal_diffusivity < 0.0)
    ):
        raise ValueError("diffusivities must be finite and non-negative")
    plane_shape = (nz, ny)
    ambient_total = _plane(
        ambient_total_mass_flux_kg_m2_s, plane_shape,
        "ambient_total_mass_flux_kg_m2_s",
    )
    ambient_v = _plane(ambient_v_m_s, plane_shape, "ambient_v_m_s")
    ambient_w = _plane(ambient_w_m_s, plane_shape, "ambient_w_m_s")
    specific_enthalpy = enthalpy_density / density

    if x_face_index == 0:
        rho_face = density[:, :, 0]
        mdot = rho_face * u_face[:, :, 0]
        left_species = right_species = mass_fraction[:, :, 0]
        left_enthalpy = right_enthalpy = specific_enthalpy[:, :, 0]
        tangential_v = cell_v[:, :, 0]
        tangential_w = cell_w[:, :, 0]
        fluid_face = ~solid[:, :, 0]
        species_diffusion = np.zeros(plane_shape)
        heat_conduction = np.zeros(plane_shape)
    elif x_face_index == nx:
        rho_face = density[:, :, -1]
        mdot = rho_face * u_face[:, :, -1]
        left_species = right_species = mass_fraction[:, :, -1]
        left_enthalpy = right_enthalpy = specific_enthalpy[:, :, -1]
        tangential_v = cell_v[:, :, -1]
        tangential_w = cell_w[:, :, -1]
        fluid_face = ~solid[:, :, -1]
        species_diffusion = np.zeros(plane_shape)
        heat_conduction = np.zeros(plane_shape)
    else:
        left = x_face_index - 1
        right = x_face_index
        rho_face = 0.5 * (density[:, :, left] + density[:, :, right])
        mdot = rho_face * u_face[:, :, x_face_index]
        fluid_face = (~solid[:, :, left]) & (~solid[:, :, right])
        velocity = np.zeros(plane_shape)
        positive_density = rho_face > 0.0
        velocity[positive_density] = mdot[positive_density] / rho_face[positive_density]
        if advection_scheme == "upwind":
            left_species = mass_fraction[:, :, left]
            right_species = mass_fraction[:, :, right]
            left_enthalpy = specific_enthalpy[:, :, left]
            right_enthalpy = specific_enthalpy[:, :, right]
        else:
            species_slope = _limited_x_slope(mass_fraction, solid)
            enthalpy_slope = _limited_x_slope(specific_enthalpy, solid)
            courant = velocity * dt / spacing[0]
            left_species = mass_fraction[:, :, left] + (
                0.5 * (1.0 - courant) * species_slope[:, :, left]
            )
            right_species = mass_fraction[:, :, right] - (
                0.5 * (1.0 + courant) * species_slope[:, :, right]
            )
            left_enthalpy = specific_enthalpy[:, :, left] + (
                0.5 * (1.0 - courant) * enthalpy_slope[:, :, left]
            )
            right_enthalpy = specific_enthalpy[:, :, right] - (
                0.5 * (1.0 + courant) * enthalpy_slope[:, :, right]
            )
        tangential_v = np.where(
            mdot >= 0.0, cell_v[:, :, left], cell_v[:, :, right],
        )
        tangential_w = np.where(
            mdot >= 0.0, cell_w[:, :, left], cell_w[:, :, right],
        )
        face_diffusivity = 0.5 * (
            diffusivity[:, :, left] + diffusivity[:, :, right]
        )
        species_diffusion = -rho_face * face_diffusivity * (
            mass_fraction[:, :, right] - mass_fraction[:, :, left]
        ) / spacing[0]
        face_thermal = 0.5 * (
            thermal_diffusivity[:, :, left]
            + thermal_diffusivity[:, :, right]
        )
        face_cp = 0.5 * (
            heat_capacity[:, :, left] + heat_capacity[:, :, right]
        )
        heat_conduction = -rho_face * face_cp * face_thermal * (
            temperature[:, :, right] - temperature[:, :, left]
        ) / spacing[0]

    advective_species_flux = mdot * np.where(
        mdot >= 0.0, left_species, right_species,
    )
    species_flux = advective_species_flux + species_diffusion
    advective_enthalpy_flux = mdot * np.where(
        mdot >= 0.0, left_enthalpy, right_enthalpy,
    )
    enthalpy_flux = advective_enthalpy_flux + heat_conduction
    total_anomaly = mdot - ambient_total
    lateral_anomaly = mdot * tangential_v - ambient_total * ambient_v
    vertical_anomaly = mdot * tangential_w - ambient_total * ambient_w
    # Explicit marked-branch closure for the 2D receiver. Its mass is distinct
    # from the signed ambient anomaly. All momentum components use this SAME
    # branch mass and the solved donor velocity, never a wind-speed guess.
    branch_mass = np.maximum(total_anomaly, np.maximum(species_flux, 0.0))
    normal_velocity = u_face[:, :, x_face_index]
    branch = {
        'mixture_mass_flux_kg_m2_s': branch_mass,
        'species_mass_flux_kg_m2_s': species_flux.copy(),
        'momentum_x_flux_pa': branch_mass * normal_velocity,
        'momentum_y_flux_pa': branch_mass * tangential_v,
        'momentum_z_flux_pa': branch_mass * tangential_w,
    }
    for field in branch.values():
        field[~fluid_face] = 0.0
    for field in (
        total_anomaly, advective_species_flux, species_flux,
        advective_enthalpy_flux, enthalpy_flux,
        lateral_anomaly, vertical_anomaly,
    ):
        field[~fluid_face] = 0.0
    area = spacing[1] * spacing[2]
    return BinaryMixtureHandoffPlane3DFrame(
        schema=BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
        time_s=query_time,
        x_face_index=x_face_index,
        total_mass_flux_anomaly_kg_m2_s=total_anomaly,
        advective_species_mass_flux_kg_m2_s=advective_species_flux,
        species_mass_flux_kg_m2_s=species_flux,
        advective_sensible_enthalpy_flux_w_m2=advective_enthalpy_flux,
        sensible_enthalpy_flux_w_m2=enthalpy_flux,
        lateral_momentum_flux_anomaly_pa=lateral_anomaly,
        vertical_momentum_flux_anomaly_pa=vertical_anomaly,
        diagnostics={
            "advection_scheme": advection_scheme,
            "species_flux_decomposition": "advective plus molecular/SGS diffusive",
            "sensible_enthalpy_flux_decomposition": "advective plus conductive",
            "boundary_or_interior": (
                "boundary" if x_face_index in {0, nx} else "interior"
            ),
            "fluid_face_count": int(np.count_nonzero(fluid_face)),
            "backflow_face_count": int(np.count_nonzero(fluid_face & (mdot < 0.0))),
            "integrated_total_mass_flux_anomaly_kg_s": float(
                np.sum(total_anomaly) * area
            ),
            "integrated_species_mass_flux_kg_s": float(
                np.sum(species_flux) * area
            ),
            "integrated_sensible_enthalpy_flux_w": float(
                np.sum(enthalpy_flux) * area
            ),
            "integrated_lateral_momentum_flux_n": float(
                np.sum(lateral_anomaly) * area
            ),
            "integrated_vertical_momentum_flux_n": float(
                np.sum(vertical_anomaly) * area
            ),
            "clipping_or_posthoc_renormalisation": False,
        },
        dense_layer_packet=branch,
    )


__all__ = [
    "BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA",
    "BinaryMixtureHandoffPlane3DFrame",
    "BinaryMixtureHandoffPlane3DRequest",
    "extract_binary_mixture_handoff_plane_3d",
]
