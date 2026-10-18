"""Conservative finite-volume momentum predictor for low-Mach 3-D flow.

The predictor advances cell momentum with the same face mass flux that carries
species and enthalpy.  Pressure is intentionally absent: a single subsequent
variable-density projection supplies the pressure correction.  This avoids
stacking an isothermal lattice pressure and a thermodynamic pressure solve.

The present operator is first-order upwind in momentum advection and explicit
in viscosity.  It is a verified building block, not yet a qualified obstacle
dispersion solver.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


VARIABLE_DENSITY_MOMENTUM_3D_SCHEMA = (
    "degali.variable-density-momentum-predictor-3d.v2"
)


@dataclass(frozen=True)
class VariableDensityMomentum3DResult:
    schema: str
    provisional_u_m_s: np.ndarray
    provisional_v_m_s: np.ndarray
    provisional_w_m_s: np.ndarray
    provisional_momentum_x_kg_m2_s: np.ndarray
    provisional_momentum_y_kg_m2_s: np.ndarray
    provisional_momentum_z_kg_m2_s: np.ndarray
    advective_divergence_x_n_m3: np.ndarray
    advective_divergence_y_n_m3: np.ndarray
    advective_divergence_z_n_m3: np.ndarray
    viscous_divergence_x_n_m3: np.ndarray
    viscous_divergence_y_n_m3: np.ndarray
    viscous_divergence_z_n_m3: np.ndarray
    momentum_initial_kg_m_s: np.ndarray
    momentum_provisional_kg_m_s: np.ndarray
    net_advective_outflow_force_n: np.ndarray
    viscous_boundary_force_n: np.ndarray
    body_force_n: np.ndarray
    source_momentum_force_n: np.ndarray
    momentum_balance_residual_kg_m_s: np.ndarray
    advective_courant_number: float
    viscous_stability_number: float


def _finite_array(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


def smooth_wall_log_law_force_density_3d(
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    w_m_s: np.ndarray,
    density_kg_m3: np.ndarray,
    *,
    surface_area_m2: np.ndarray,
    area_normal_x_m2: np.ndarray,
    area_normal_y_m2: np.ndarray,
    wall_distance_m: np.ndarray,
    spacing_m: tuple[float, float, float],
    molecular_kinematic_viscosity_m2_s: float,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], dict[str, float | int]]:
    """Return a smooth-wall immersed-boundary shear force per fluid volume.

    A Reichardt smooth-wall law maps the cell-centred tangential velocity and
    wall distance to friction velocity.  The resulting stress is applied over
    the exact embedded surface area supplied by the geometry layer.  This is
    an equilibrium smooth-wall model, not a resolved viscous sublayer; rough
    wall sensitivity and separated-flow validation remain necessary.
    """
    velocity = tuple(np.asarray(item, dtype=float) for item in (
        u_m_s, v_m_s, w_m_s,
    ))
    if velocity[0].ndim != 3:
        raise ValueError("velocity fields must be three-dimensional arrays")
    shape = velocity[0].shape
    for name, field in zip(("u_m_s", "v_m_s", "w_m_s"), velocity):
        if field.shape != shape:
            raise ValueError(f"{name} must have shape {shape}")
    density = np.asarray(density_kg_m3, dtype=float)
    area = np.asarray(surface_area_m2, dtype=float)
    normal_x_area = np.asarray(area_normal_x_m2, dtype=float)
    normal_y_area = np.asarray(area_normal_y_m2, dtype=float)
    for name, field in (
        ("density_kg_m3", density), ("surface_area_m2", area),
        ("area_normal_x_m2", normal_x_area),
        ("area_normal_y_m2", normal_y_area),
    ):
        if field.shape != shape:
            raise ValueError(f"{name} must have shape {shape}")
    distance = np.asarray(wall_distance_m, dtype=float)
    if distance.shape != shape or np.any(np.isnan(distance)):
        raise ValueError(f"wall_distance_m must have shape {shape} and no NaNs")
    area_flat = area.ravel()
    if np.any(~np.isfinite(area_flat)) or np.any(area_flat < 0.0):
        raise ValueError("surface_area_m2 must be non-negative")
    if len(spacing_m) != 3 or any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in spacing_m
    ):
        raise ValueError("spacing_m must contain three positive finite values")
    viscosity = float(molecular_kinematic_viscosity_m2_s)
    if not math.isfinite(viscosity) or viscosity <= 0.0:
        raise ValueError("molecular_kinematic_viscosity_m2_s must be positive")
    surface_indices = np.flatnonzero(area_flat > 0.0)
    force = tuple(np.zeros(shape, dtype=float) for _ in range(3))
    if surface_indices.size == 0:
        return force, {
            "surface_cell_count": 0,
            "maximum_tangential_speed_m_s": 0.0,
            "maximum_friction_velocity_m_s": 0.0,
            "maximum_wall_shear_stress_pa": 0.0,
            "net_wall_force_x_n": 0.0,
            "net_wall_force_y_n": 0.0,
            "net_wall_force_z_n": 0.0,
        }
    flat = tuple(item.ravel() for item in velocity)
    density_at_surface = density.ravel()[surface_indices]
    area_at_surface = area_flat[surface_indices]
    normal_x_at_surface = normal_x_area.ravel()[surface_indices]
    normal_y_at_surface = normal_y_area.ravel()[surface_indices]
    distance_at_surface = distance.ravel()[surface_indices]
    velocity_at_surface = tuple(item[surface_indices] for item in flat)
    if any(not np.all(np.isfinite(item)) for item in (
        density_at_surface, area_at_surface, normal_x_at_surface,
        normal_y_at_surface, distance_at_surface, *velocity_at_surface,
    )):
        raise ValueError("wall-cell state and geometry must be finite")
    if np.any(density_at_surface <= 0.0):
        raise ValueError("density_kg_m3 must be positive at wall cells")
    if np.any(distance_at_surface < 0.0):
        raise ValueError("wall distances must be finite and non-negative on surface cells")
    normal_magnitude = np.hypot(normal_x_at_surface, normal_y_at_surface)
    has_normal = normal_magnitude > 0.0
    if not np.any(has_normal):
        return force, {
            "surface_cell_count": 0,
            "maximum_tangential_speed_m_s": 0.0,
            "maximum_friction_velocity_m_s": 0.0,
            "maximum_wall_shear_stress_pa": 0.0,
            "net_wall_force_x_n": 0.0,
            "net_wall_force_y_n": 0.0,
            "net_wall_force_z_n": 0.0,
        }

    surface_indices = surface_indices[has_normal]
    density_at_surface = density_at_surface[has_normal]
    area_at_surface = area_at_surface[has_normal]
    distance_at_surface = distance_at_surface[has_normal]
    nx = normal_x_at_surface[has_normal] / normal_magnitude[has_normal]
    ny = normal_y_at_surface[has_normal] / normal_magnitude[has_normal]
    velocity_at_surface = tuple(item[has_normal] for item in velocity_at_surface)
    normal_velocity = velocity_at_surface[0] * nx + velocity_at_surface[1] * ny
    tangential = (
        velocity_at_surface[0] - normal_velocity * nx,
        velocity_at_surface[1] - normal_velocity * ny,
        velocity_at_surface[2],
    )
    speed = np.sqrt(sum(component * component for component in tangential))
    moving = speed > 0.0
    if not np.any(moving):
        return force, {
            "surface_cell_count": 0,
            "maximum_tangential_speed_m_s": 0.0,
            "maximum_friction_velocity_m_s": 0.0,
            "maximum_wall_shear_stress_pa": 0.0,
            "net_wall_force_x_n": 0.0,
            "net_wall_force_y_n": 0.0,
            "net_wall_force_z_n": 0.0,
        }

    # Keep the first-cell distance finite when a grid point lies exactly on
    # the mathematical surface; this is an explicit unresolved-wall floor.
    surface_indices = surface_indices[moving]
    density_active = density_at_surface[moving]
    area_active = area_at_surface[moving]
    distance_eff = np.maximum(
        distance_at_surface[moving],
        0.25 * min(float(value) for value in spacing_m),
    )
    tangential_active = tuple(item[moving] for item in tangential)
    speed_active = speed[moving]
    low = np.zeros_like(speed_active)
    high = speed_active.copy()
    kappa = 0.41
    for _ in range(40):
        friction = 0.5 * (low + high)
        y_plus = friction * distance_eff / viscosity
        scaled_y = y_plus / 11.0
        u_plus = (
            np.log1p(kappa * y_plus) / kappa
            + 7.8 * (
                1.0 - np.exp(-scaled_y)
                - scaled_y * np.exp(-y_plus / 3.0)
            )
        )
        predicted_speed = friction * u_plus
        below = predicted_speed < speed_active
        low = np.where(below, friction, low)
        high = np.where(below, high, friction)
    friction_velocity = 0.5 * (low + high)
    shear_stress = density_active * friction_velocity**2
    cell_volume = math.prod(float(value) for value in spacing_m)
    force_scale = shear_stress * area_active / (cell_volume * speed_active)
    net_force = []
    for index in range(3):
        force_values = -force_scale * tangential_active[index]
        force[index].ravel()[surface_indices] = force_values
        net_force.append(float(np.sum(force_values) * cell_volume))
    return force, {
        "surface_cell_count": int(surface_indices.size),
        "maximum_tangential_speed_m_s": float(np.max(speed_active)),
        "maximum_friction_velocity_m_s": float(np.max(friction_velocity)),
        "maximum_wall_shear_stress_pa": float(np.max(shear_stress)),
        "net_wall_force_x_n": net_force[0],
        "net_wall_force_y_n": net_force[1],
        "net_wall_force_z_n": net_force[2],
    }


def _inlet_plane(
    value: float | np.ndarray, nz: int, ny: int, name: str,
) -> np.ndarray:
    """Broadcast a scalar, z profile, or complete z-y inlet plane."""
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        plane = np.full((nz, ny), float(array), dtype=float)
    elif array.shape == (nz,):
        plane = np.broadcast_to(array[:, None], (nz, ny)).copy()
    elif array.shape == (nz, ny):
        plane = array.copy()
    else:
        raise ValueError(
            f"{name} must be scalar, shape {(nz,)}, or shape {(nz, ny)}"
        )
    if not np.all(np.isfinite(plane)):
        raise ValueError(f"{name} must be finite")
    return plane


def _upwind_momentum_fluxes(
    scalar: np.ndarray,
    mass_flux_x: np.ndarray,
    mass_flux_y: np.ndarray,
    mass_flux_z: np.ndarray,
    solid: np.ndarray,
    *,
    component_axis: int,
    inlet_value: np.ndarray,
    nonstreamwise_wall_boundary: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return face fluxes of one velocity component, kg/(m s2)."""
    nz, ny, nx = solid.shape
    flux_x = np.zeros((nz, ny, nx + 1), dtype=float)
    flux_y = np.zeros((nz, ny + 1, nx), dtype=float)
    flux_z = np.zeros((nz + 1, ny, nx), dtype=float)

    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:])
    interior = mass_flux_x[:, :, 1:-1]
    face_value = np.where(interior >= 0.0, scalar[:, :, :-1], scalar[:, :, 1:])
    flux_x[:, :, 1:-1][fluid_x] = (interior * face_value)[fluid_x]
    inlet_state = inlet_value
    low_x_value = np.where(
        mass_flux_x[:, :, 0] >= 0.0, inlet_state, scalar[:, :, 0],
    )
    flux_x[:, :, 0] = mass_flux_x[:, :, 0] * low_x_value
    # The high-x pressure outlet exchanges flow with quiescent ambient air.
    # Use the interior velocity for outflow, but ambient-at-rest velocity for
    # backflow; otherwise negative mass flux spuriously imports the interior
    # momentum state instead of the exterior boundary state.
    high_x_mass_flux = mass_flux_x[:, :, -1]
    high_x_value = np.where(high_x_mass_flux >= 0.0, scalar[:, :, -1], 0.0)
    flux_x[:, :, -1] = high_x_mass_flux * high_x_value
    flux_x[:, :, 0][solid[:, :, 0]] = 0.0
    flux_x[:, :, -1][solid[:, :, -1]] = 0.0

    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :])
    interior = mass_flux_y[:, 1:-1, :]
    face_value = np.where(interior >= 0.0, scalar[:, :-1, :], scalar[:, 1:, :])
    flux_y[:, 1:-1, :][fluid_y] = (interior * face_value)[fluid_y]
    pressure_open = nonstreamwise_wall_boundary == "all_pressure_open"
    low_y_boundary = (
        scalar[:, 0, :]
        if pressure_open or component_axis != 1 else np.zeros((nz, nx))
    )
    high_y_boundary = (
        scalar[:, -1, :]
        if pressure_open or component_axis != 1 else np.zeros((nz, nx))
    )
    low_y_value = np.where(
        mass_flux_y[:, 0, :] >= 0.0, low_y_boundary, scalar[:, 0, :],
    )
    high_y_value = np.where(
        mass_flux_y[:, -1, :] >= 0.0, scalar[:, -1, :], high_y_boundary,
    )
    flux_y[:, 0, :] = mass_flux_y[:, 0, :] * low_y_value
    flux_y[:, -1, :] = mass_flux_y[:, -1, :] * high_y_value
    flux_y[:, 0, :][solid[:, 0, :]] = 0.0
    flux_y[:, -1, :][solid[:, -1, :]] = 0.0

    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :])
    interior = mass_flux_z[1:-1, :, :]
    face_value = np.where(interior >= 0.0, scalar[:-1, :, :], scalar[1:, :, :])
    flux_z[1:-1, :, :][fluid_z] = (interior * face_value)[fluid_z]
    # Pressure-open faces use a zero-gradient exterior velocity state.  Wall
    # and symmetry modes retain their declared normal-velocity constraints.
    low_z_boundary = (
        scalar[0, :, :] if pressure_open else np.zeros((ny, nx))
    )
    high_z_boundary = (
        scalar[-1, :, :]
        if pressure_open or component_axis != 2 else np.zeros((ny, nx))
    )
    low_z_value = np.where(
        mass_flux_z[0, :, :] >= 0.0, low_z_boundary, scalar[0, :, :],
    )
    high_z_value = np.where(
        mass_flux_z[-1, :, :] >= 0.0, scalar[-1, :, :], high_z_boundary,
    )
    flux_z[0, :, :] = mass_flux_z[0, :, :] * low_z_value
    flux_z[-1, :, :] = mass_flux_z[-1, :, :] * high_z_value
    flux_z[0, :, :][solid[0, :, :]] = 0.0
    flux_z[-1, :, :][solid[-1, :, :]] = 0.0
    return flux_x, flux_y, flux_z


def _viscous_fluxes(
    scalar: np.ndarray,
    dynamic_viscosity: np.ndarray,
    solid: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    component_axis: int,
    inlet_value: np.ndarray,
    nonstreamwise_wall_boundary: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return mu grad(u_i) on faces with declared wall conditions."""
    nz, ny, nx = solid.shape
    dx, dy, dz = spacing
    flux_x = np.zeros((nz, ny, nx + 1), dtype=float)
    flux_y = np.zeros((nz, ny + 1, nx), dtype=float)
    flux_z = np.zeros((nz + 1, ny, nx), dtype=float)

    left_fluid = ~solid[:, :, :-1]
    right_fluid = ~solid[:, :, 1:]
    both = left_fluid & right_fluid
    mu_face = 0.5 * (
        dynamic_viscosity[:, :, :-1] + dynamic_viscosity[:, :, 1:]
    )
    gradient = (scalar[:, :, 1:] - scalar[:, :, :-1]) / dx
    flux_x[:, :, 1:-1][both] = (mu_face * gradient)[both]
    left_only = left_fluid & ~right_fluid
    right_only = ~left_fluid & right_fluid
    flux_x[:, :, 1:-1][left_only] = (
        dynamic_viscosity[:, :, :-1][left_only]
        * (-2.0 * scalar[:, :, :-1][left_only] / dx)
    )
    flux_x[:, :, 1:-1][right_only] = (
        dynamic_viscosity[:, :, 1:][right_only]
        * (2.0 * scalar[:, :, 1:][right_only] / dx)
    )
    inlet = inlet_value
    flux_x[:, :, 0] = (
        dynamic_viscosity[:, :, 0] * 2.0
        * (scalar[:, :, 0] - inlet) / dx
    )
    flux_x[:, :, 0][solid[:, :, 0]] = 0.0
    # Zero-gradient outflow.

    low_fluid = ~solid[:, :-1, :]
    high_fluid = ~solid[:, 1:, :]
    both = low_fluid & high_fluid
    mu_face = 0.5 * (
        dynamic_viscosity[:, :-1, :] + dynamic_viscosity[:, 1:, :]
    )
    gradient = (scalar[:, 1:, :] - scalar[:, :-1, :]) / dy
    flux_y[:, 1:-1, :][both] = (mu_face * gradient)[both]
    low_only = low_fluid & ~high_fluid
    high_only = ~low_fluid & high_fluid
    flux_y[:, 1:-1, :][low_only] = (
        dynamic_viscosity[:, :-1, :][low_only]
        * (-2.0 * scalar[:, :-1, :][low_only] / dy)
    )
    flux_y[:, 1:-1, :][high_only] = (
        dynamic_viscosity[:, 1:, :][high_only]
        * (2.0 * scalar[:, 1:, :][high_only] / dy)
    )
    if nonstreamwise_wall_boundary != "all_pressure_open" and (
        component_axis == 1
        or nonstreamwise_wall_boundary == "all_no_slip_tunnel"
    ):
        flux_y[:, 0, :] = (
            dynamic_viscosity[:, 0, :] * 2.0 * scalar[:, 0, :] / dy
        )
        flux_y[:, -1, :] = (
            dynamic_viscosity[:, -1, :] * -2.0 * scalar[:, -1, :] / dy
        )
    flux_y[:, 0, :][solid[:, 0, :]] = 0.0
    flux_y[:, -1, :][solid[:, -1, :]] = 0.0

    low_fluid = ~solid[:-1, :, :]
    high_fluid = ~solid[1:, :, :]
    both = low_fluid & high_fluid
    mu_face = 0.5 * (
        dynamic_viscosity[:-1, :, :] + dynamic_viscosity[1:, :, :]
    )
    gradient = (scalar[1:, :, :] - scalar[:-1, :, :]) / dz
    flux_z[1:-1, :, :][both] = (mu_face * gradient)[both]
    low_only = low_fluid & ~high_fluid
    high_only = ~low_fluid & high_fluid
    flux_z[1:-1, :, :][low_only] = (
        dynamic_viscosity[:-1, :, :][low_only]
        * (-2.0 * scalar[:-1, :, :][low_only] / dz)
    )
    flux_z[1:-1, :, :][high_only] = (
        dynamic_viscosity[1:, :, :][high_only]
        * (2.0 * scalar[1:, :, :][high_only] / dz)
    )
    # The facility default has a no-slip ground.  A fully symmetric/free-slip
    # surrogate retains zero tangential gradient at the lower face; its normal
    # velocity is still zero through the face-flux boundary.
    if nonstreamwise_wall_boundary != "all_pressure_open" and (
        component_axis == 2
        or nonstreamwise_wall_boundary != "all_symmetry_free_slip"
    ):
        flux_z[0, :, :] = (
            dynamic_viscosity[0, :, :] * 2.0 * scalar[0, :, :] / dz
        )
    # The open-atmosphere mode keeps a free-slip top.  A tunnel roof is a
    # viscous no-slip wall for every velocity component.
    if nonstreamwise_wall_boundary != "all_pressure_open" and (
        component_axis == 2
        or nonstreamwise_wall_boundary == "all_no_slip_tunnel"
    ):
        flux_z[-1, :, :] = (
            dynamic_viscosity[-1, :, :] * -2.0 * scalar[-1, :, :] / dz
        )
    flux_z[0, :, :][solid[0, :, :]] = 0.0
    flux_z[-1, :, :][solid[-1, :, :]] = 0.0
    return flux_x, flux_y, flux_z


def _divergence(
    flux_x: np.ndarray,
    flux_y: np.ndarray,
    flux_z: np.ndarray,
    spacing: tuple[float, float, float],
) -> np.ndarray:
    dx, dy, dz = spacing
    return (
        (flux_x[:, :, 1:] - flux_x[:, :, :-1]) / dx
        + (flux_y[:, 1:, :] - flux_y[:, :-1, :]) / dy
        + (flux_z[1:, :, :] - flux_z[:-1, :, :]) / dz
    )


def _net_boundary_flux(
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


def advance_variable_density_momentum_3d(
    density_kg_m3: np.ndarray,
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    w_m_s: np.ndarray,
    *,
    mass_flux_x_kg_m2_s: np.ndarray,
    mass_flux_y_kg_m2_s: np.ndarray,
    mass_flux_z_kg_m2_s: np.ndarray,
    solid_mask: np.ndarray,
    spacing_m: tuple[float, float, float] | np.ndarray,
    time_step_s: float,
    kinematic_viscosity_m2_s: float | np.ndarray,
    acceleration_x_m_s2: float | np.ndarray = 0.0,
    acceleration_y_m_s2: float | np.ndarray = 0.0,
    acceleration_z_m_s2: float | np.ndarray = 0.0,
    mass_source_rate_kg_m3_s: float | np.ndarray = 0.0,
    source_velocity_m_s: tuple[float, float, float] = (0.0, 0.0, 0.0),
    inlet_u_m_s: float | np.ndarray = 0.0,
    inlet_v_m_s: float | np.ndarray = 0.0,
    inlet_w_m_s: float | np.ndarray = 0.0,
    nonstreamwise_wall_boundary: str = "ground_no_slip_lateral_top_free_slip",
) -> VariableDensityMomentum3DResult:
    """Advance conservative momentum without pressure, ready for projection."""
    solid = np.asarray(solid_mask)
    if solid.dtype != np.bool_ or solid.ndim != 3:
        raise ValueError("solid_mask must be a three-dimensional boolean array")
    shape = solid.shape
    nz, ny, nx = shape
    density = _finite_array(density_kg_m3, shape, "density_kg_m3").copy()
    if np.any(density[~solid] <= 0.0):
        raise ValueError("fluid density must be positive")
    velocities = [
        _finite_array(value, shape, name).copy()
        for value, name in (
            (u_m_s, "u_m_s"), (v_m_s, "v_m_s"), (w_m_s, "w_m_s"),
        )
    ]
    for velocity in velocities:
        velocity[solid] = 0.0
    mass_fluxes = (
        _finite_array(
            mass_flux_x_kg_m2_s, (nz, ny, nx + 1),
            "mass_flux_x_kg_m2_s",
        ),
        _finite_array(
            mass_flux_y_kg_m2_s, (nz, ny + 1, nx),
            "mass_flux_y_kg_m2_s",
        ),
        _finite_array(
            mass_flux_z_kg_m2_s, (nz + 1, ny, nx),
            "mass_flux_z_kg_m2_s",
        ),
    )
    spacing_array = np.asarray(spacing_m, dtype=float)
    if spacing_array.shape != (3,) or np.any(~np.isfinite(spacing_array)) or np.any(
        spacing_array <= 0.0
    ):
        raise ValueError("spacing_m must contain three positive finite values")
    spacing = tuple(float(value) for value in spacing_array)
    dt = float(time_step_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_s must be positive and finite")
    viscosity = np.broadcast_to(
        np.asarray(kinematic_viscosity_m2_s, dtype=float), shape,
    ).copy()
    if np.any(~np.isfinite(viscosity)) or np.any(viscosity < 0.0):
        raise ValueError("kinematic_viscosity_m2_s must be finite and non-negative")
    viscosity[solid] = 0.0
    dynamic_viscosity = density * viscosity
    accelerations = [
        np.broadcast_to(np.asarray(value, dtype=float), shape).copy()
        for value in (
            acceleration_x_m_s2, acceleration_y_m_s2, acceleration_z_m_s2,
        )
    ]
    if any(np.any(~np.isfinite(value)) for value in accelerations):
        raise ValueError("accelerations must be finite and broadcastable")
    source = np.broadcast_to(
        np.asarray(mass_source_rate_kg_m3_s, dtype=float), shape,
    ).copy()
    if np.any(~np.isfinite(source)):
        raise ValueError("mass_source_rate_kg_m3_s must be finite and broadcastable")
    source[solid] = 0.0
    if (
        not isinstance(source_velocity_m_s, tuple)
        or len(source_velocity_m_s) != 3
        or any(not math.isfinite(float(value)) for value in source_velocity_m_s)
    ):
        raise ValueError("source_velocity_m_s must contain three finite values")
    inlet_values = [
        _inlet_plane(value, nz, ny, name)
        for value, name in (
            (inlet_u_m_s, "inlet_u_m_s"),
            (inlet_v_m_s, "inlet_v_m_s"),
            (inlet_w_m_s, "inlet_w_m_s"),
        )
    ]
    if nonstreamwise_wall_boundary not in {
        "ground_no_slip_lateral_top_free_slip",
        "all_no_slip_tunnel",
        "all_symmetry_free_slip",
        "all_pressure_open",
    }:
        raise ValueError(
            "nonstreamwise_wall_boundary must be "
            "'ground_no_slip_lateral_top_free_slip', 'all_no_slip_tunnel', "
            "'all_symmetry_free_slip', or 'all_pressure_open'"
        )

    fluid = ~solid
    cell_volume = spacing[0] * spacing[1] * spacing[2]
    advective_divergences = []
    viscous_divergences = []
    advective_outflows = []
    viscous_forces = []
    new_momenta = []
    initial_integrals = []
    final_integrals = []
    body_forces = []
    source_forces = []
    residuals = []
    for axis, (velocity, acceleration, inlet, source_velocity) in enumerate(zip(
        velocities, accelerations, inlet_values, source_velocity_m_s,
    )):
        advective_fluxes = _upwind_momentum_fluxes(
            velocity, *mass_fluxes, solid,
            component_axis=axis, inlet_value=inlet,
            nonstreamwise_wall_boundary=nonstreamwise_wall_boundary,
        )
        viscous_fluxes = _viscous_fluxes(
            velocity, dynamic_viscosity, solid, spacing,
            component_axis=axis, inlet_value=inlet,
            nonstreamwise_wall_boundary=nonstreamwise_wall_boundary,
        )
        advective_divergence = _divergence(*advective_fluxes, spacing)
        viscous_divergence = _divergence(*viscous_fluxes, spacing)
        advective_divergence[solid] = 0.0
        viscous_divergence[solid] = 0.0
        momentum_old = density * velocity
        body = density * acceleration
        source_momentum = source * float(source_velocity)
        momentum_new = momentum_old + dt * (
            -advective_divergence + viscous_divergence
            + body + source_momentum
        )
        momentum_new[solid] = 0.0
        initial = float(np.sum(momentum_old[fluid]) * cell_volume)
        final = float(np.sum(momentum_new[fluid]) * cell_volume)
        advective_outflow = _net_boundary_flux(*advective_fluxes, spacing)
        viscous_force = float(np.sum(viscous_divergence[fluid]) * cell_volume)
        body_force = float(np.sum(body[fluid]) * cell_volume)
        source_force = float(np.sum(source_momentum[fluid]) * cell_volume)
        residual = final - initial - dt * (
            -advective_outflow + viscous_force + body_force + source_force
        )
        advective_divergences.append(advective_divergence)
        viscous_divergences.append(viscous_divergence)
        advective_outflows.append(advective_outflow)
        viscous_forces.append(viscous_force)
        new_momenta.append(momentum_new)
        initial_integrals.append(initial)
        final_integrals.append(final)
        body_forces.append(body_force)
        source_forces.append(source_force)
        residuals.append(residual)

    provisional_velocities = []
    for momentum in new_momenta:
        velocity = np.zeros_like(momentum)
        velocity[fluid] = momentum[fluid] / density[fluid]
        provisional_velocities.append(velocity)

    minimum_density = float(np.min(density[fluid]))
    dx, dy, dz = spacing
    advective_courant = dt * (
        float(np.max(np.abs(mass_fluxes[0]))) / (minimum_density * dx)
        + float(np.max(np.abs(mass_fluxes[1]))) / (minimum_density * dy)
        + float(np.max(np.abs(mass_fluxes[2]))) / (minimum_density * dz)
    )
    viscous_stability = 2.0 * dt * float(np.max(viscosity[fluid])) * (
        1.0 / dx**2 + 1.0 / dy**2 + 1.0 / dz**2
    )
    if advective_courant > 1.0 + 1.0e-12:
        raise ValueError(
            "momentum advection exceeds the explicit Courant limit: "
            f"{advective_courant:.6g}"
        )
    if viscous_stability > 1.0 + 1.0e-12:
        raise ValueError(
            "momentum diffusion exceeds the explicit stability limit: "
            f"{viscous_stability:.6g}"
        )
    return VariableDensityMomentum3DResult(
        schema=VARIABLE_DENSITY_MOMENTUM_3D_SCHEMA,
        provisional_u_m_s=provisional_velocities[0],
        provisional_v_m_s=provisional_velocities[1],
        provisional_w_m_s=provisional_velocities[2],
        provisional_momentum_x_kg_m2_s=new_momenta[0],
        provisional_momentum_y_kg_m2_s=new_momenta[1],
        provisional_momentum_z_kg_m2_s=new_momenta[2],
        advective_divergence_x_n_m3=advective_divergences[0],
        advective_divergence_y_n_m3=advective_divergences[1],
        advective_divergence_z_n_m3=advective_divergences[2],
        viscous_divergence_x_n_m3=viscous_divergences[0],
        viscous_divergence_y_n_m3=viscous_divergences[1],
        viscous_divergence_z_n_m3=viscous_divergences[2],
        momentum_initial_kg_m_s=np.asarray(initial_integrals),
        momentum_provisional_kg_m_s=np.asarray(final_integrals),
        net_advective_outflow_force_n=np.asarray(advective_outflows),
        viscous_boundary_force_n=np.asarray(viscous_forces),
        body_force_n=np.asarray(body_forces),
        source_momentum_force_n=np.asarray(source_forces),
        momentum_balance_residual_kg_m_s=np.asarray(residuals),
        advective_courant_number=advective_courant,
        viscous_stability_number=viscous_stability,
    )


__all__ = [
    "VARIABLE_DENSITY_MOMENTUM_3D_SCHEMA",
    "VariableDensityMomentum3DResult",
    "advance_variable_density_momentum_3d",
]
