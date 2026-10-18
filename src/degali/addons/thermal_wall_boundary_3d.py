"""Conservative first-cell heat exchange for isothermal 3-D walls."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ThermalWallHeatSource3D:
    sensible_enthalpy_source_rate_w_m3: np.ndarray
    wall_heat_rate_w: float
    obstacle_wall_heat_rate_w: float
    ground_wall_heat_rate_w: float
    obstacle_wall_face_count: int
    ground_wall_cell_count: int


def isothermal_wall_heat_source_3d(
    temperature_k: np.ndarray,
    density_kg_m3: np.ndarray,
    heat_capacity_j_kg_k: np.ndarray,
    thermal_diffusivity_m2_s: float | np.ndarray | tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray],
    solid_mask: np.ndarray,
    *,
    spacing_m: tuple[float, float, float],
    obstacle_wall_temperature_k: float | None = None,
    ground_temperature_k: float | None = None,
) -> ThermalWallHeatSource3D:
    """Return wall conduction as a conservative heat-source density.

    The voxel solid/fluid center separation is one cell, so wall-to-first-fluid
    center distance is approximated as half a cell. Every exposed solid face is
    counted once. A ground wall uses the analogous half-cell distance at z_min.
    This is a first-order wall boundary approximation, not a resolved thermal
    boundary layer or conjugate heat-transfer model.
    """
    solid = np.asarray(solid_mask)
    temperature = np.asarray(temperature_k, dtype=float)
    density = np.asarray(density_kg_m3, dtype=float)
    heat_capacity = np.asarray(heat_capacity_j_kg_k, dtype=float)
    if solid.dtype != np.bool_ or solid.ndim != 3:
        raise ValueError("solid_mask must be a three-dimensional boolean array")
    shape = solid.shape
    for value, name in (
        (temperature, "temperature_k"), (density, "density_kg_m3"),
        (heat_capacity, "heat_capacity_j_kg_k"),
    ):
        if value.shape != shape or not np.all(np.isfinite(value)):
            raise ValueError(f"{name} must be finite with the solid-mask shape")
    if np.any(temperature <= 0.0) or np.any(density <= 0.0) or np.any(heat_capacity <= 0.0):
        raise ValueError("temperature, density, and heat capacity must be positive")
    if isinstance(thermal_diffusivity_m2_s, tuple):
        if len(thermal_diffusivity_m2_s) != 3:
            raise ValueError("directional thermal diffusivity must contain Dx, Dy, and Dz")
        supplied = thermal_diffusivity_m2_s
    else:
        supplied = (
            thermal_diffusivity_m2_s,
            thermal_diffusivity_m2_s,
            thermal_diffusivity_m2_s,
        )
    diffusivities = []
    for value in supplied:
        if np.ndim(value) == 0:
            field = np.full(shape, float(value), dtype=float)
        else:
            field = np.asarray(value, dtype=float)
        if field.shape != shape or not np.all(np.isfinite(field)) or np.any(field < 0.0):
            raise ValueError("each thermal diffusivity must be finite, non-negative, and match the grid")
        diffusivities.append(field)
    spacing = tuple(float(item) for item in spacing_m)
    if len(spacing) != 3 or any(not np.isfinite(item) or item <= 0.0 for item in spacing):
        raise ValueError("spacing_m must contain three positive finite values")
    for value, name in (
        (obstacle_wall_temperature_k, "obstacle_wall_temperature_k"),
        (ground_temperature_k, "ground_temperature_k"),
    ):
        if value is not None and (not np.isfinite(value) or value <= 0.0):
            raise ValueError(f"{name} must be positive and finite")

    source = np.zeros(shape, dtype=float)
    # Diffusivity tuples are in xyz order; array axes are zyx.
    axis_conductivity = tuple(
        density * heat_capacity * diffusivities[2 - axis]
        for axis in range(3)
    )
    obstacle_faces = 0
    obstacle_power = 0.0
    if obstacle_wall_temperature_k is not None:
        for axis in range(3):
            distance = spacing[2 - axis]
            distance = spacing[2 - axis]
            low = [slice(None)] * 3
            high = [slice(None)] * 3
            low[axis] = slice(None, -1)
            high[axis] = slice(1, None)
            low_key, high_key = tuple(low), tuple(high)
            low_is_fluid_high_is_solid = (~solid[low_key]) & solid[high_key]
            high_is_fluid_low_is_solid = solid[low_key] & (~solid[high_key])
            for fluid_key, mask in (
                (low_key, low_is_fluid_high_is_solid),
                (high_key, high_is_fluid_low_is_solid),
            ):
                if not np.any(mask):
                    continue
                delta_t = obstacle_wall_temperature_k - temperature[fluid_key]
                cell_rate = (
                    2.0 * axis_conductivity[axis][fluid_key]
                    * delta_t / distance**2
                )
                local_source = source[fluid_key]
                local_source[mask] += cell_rate[mask]
                source[fluid_key] = local_source
                obstacle_power += float(np.sum(cell_rate[mask]) * np.prod(spacing))
                obstacle_faces += int(np.count_nonzero(mask))

    ground_cells = 0
    ground_power = 0.0
    if ground_temperature_k is not None:
        fluid_bottom = ~solid[0]
        if np.any(fluid_bottom):
            dz = spacing[2]
            ground_rate = 2.0 * axis_conductivity[0][0] * (
                ground_temperature_k - temperature[0]
            ) / dz**2
            source[0][fluid_bottom] += ground_rate[fluid_bottom]
            ground_power = float(np.sum(ground_rate[fluid_bottom]) * np.prod(spacing))
            ground_cells = int(np.count_nonzero(fluid_bottom))
    source[solid] = 0.0
    return ThermalWallHeatSource3D(
        sensible_enthalpy_source_rate_w_m3=source,
        wall_heat_rate_w=obstacle_power + ground_power,
        obstacle_wall_heat_rate_w=obstacle_power,
        ground_wall_heat_rate_w=ground_power,
        obstacle_wall_face_count=obstacle_faces,
        ground_wall_cell_count=ground_cells,
    )
