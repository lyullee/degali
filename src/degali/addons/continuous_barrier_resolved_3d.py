"""Fast resolved-physics configuration for a continuous H2 wall release.

This module is the authoritative path for the presentation barrier scenario.
It deliberately uses the variable-density finite-volume momentum predictor,
pressure projection and thermal/species corrector from
``resolved_obstacle_flow_3d``.  The older animation-only path blocked scalar
flux at the wall while relaxing the carrier velocity back to the ambient wind;
that path cannot generate stagnation, separation or a wake and must not be
used for obstacle-physics claims.

The grid is intentionally coarse (5 m cubic cells) so a 120 s demonstration
can be calculated interactively.  It resolves the 15 m wall by three cells and
keeps lateral and top bypass open.  This is a physics-bearing reduced-order
calculation, not an externally validated LES/RANS prediction.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    ResolvedObstacleFlow3DResult,
)
from .semi_fv_obstacle import SourceRateSchedule
from .site_geometry import AxisAlignedCuboid


def _finite_positive(value: float, name: str) -> float:
    checked = float(value)
    if not math.isfinite(checked) or checked <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return checked


def _release_schedule(
    *, duration_s: float, release_duration_s: float, release_rate_kg_s: float,
    layout: str,
) -> SourceRateSchedule:
    duration = _finite_positive(duration_s, "duration_s")
    release_duration = _finite_positive(
        release_duration_s, "release_duration_s",
    )
    rate = _finite_positive(release_rate_kg_s, "release_rate_kg_s")
    if release_duration > duration:
        raise ValueError("release_duration_s cannot exceed duration_s")
    taper = min(release_duration - 1.0e-3, duration - 1.0e-3)
    if taper <= 0.0:
        taper = 0.5 * release_duration
    if math.isclose(release_duration, duration):
        times = (0.0, taper, duration)
        rates = (rate, rate, 0.0)
    else:
        times = (0.0, taper, release_duration, duration)
        rates = (rate, rate, 0.0, 0.0)
    return SourceRateSchedule(
        times,
        rates,
        source_id=f"continuous-h2-{layout}-resolved-wall-v1",
    )


def build_continuous_barrier_resolved_config(
    *,
    duration_s: float = 120.0,
    release_duration_s: float = 120.0,
    release_rate_kg_s: float = 0.10,
    wind_speed_m_s: float = 2.2,
    barrier_layout: str = "blocking-wall",
    store_scalar_fields: bool = True,
) -> tuple[ResolvedObstacleFlow3DConfig, dict[str, Any]]:
    """Build the interactive, pressure-projected barrier calculation.

    ``blocking-wall`` means a wide wall, not a sealed cross-section.  A sealed
    wall would force all inflow through the roof opening and turn the chosen
    open-atmosphere problem into a duct contraction.  Fifteen-metre side fetch
    is retained so pressure can relax and lateral bypass remains physical.
    """
    if barrier_layout not in {
        "single", "split-gap", "blocking-wall", "building-cluster",
    }:
        raise ValueError("unsupported barrier_layout")
    duration = _finite_positive(duration_s, "duration_s")
    release_duration = _finite_positive(
        release_duration_s, "release_duration_s",
    )
    release_rate = _finite_positive(
        release_rate_kg_s, "release_rate_kg_s",
    )
    wind_speed = _finite_positive(wind_speed_m_s, "wind_speed_m_s")

    domain = (140.0, 80.0, 40.0)
    cells = (28, 16, 8)  # exactly 5 m cubic cells
    source = (20.0, 40.0, 2.5)
    barrier: dict[str, Any] = {
        "layout": barrier_layout,
        "x_m": 75.0,
        "thickness_m": 5.0,
        "height_m": 15.0,
        "y_min_m": 15.0,
        "y_max_m": 65.0,
    }
    if barrier_layout == "single":
        barrier.update({"y_min_m": 20.0, "y_max_m": 60.0})
    elif barrier_layout == "split-gap":
        barrier.update({
            "segments": (
                {"y_min_m": 15.0, "y_max_m": 30.0},
                {"y_min_m": 50.0, "y_max_m": 65.0},
            ),
            "gap_y_min_m": 30.0,
            "gap_y_max_m": 50.0,
            "gap_width_m": 20.0,
        })
    elif barrier_layout == "building-cluster":
        # Staggered urban block: the front-centre building divides the plume,
        # the two middle buildings form offset street canyons and the taller
        # rear building forces a second impact/overtopping event.
        domain = (160.0, 120.0, 50.0)
        cells = (32, 24, 10)  # retain exactly 5 m cubic cells
        source = (20.0, 60.0, 2.5)
        buildings = (
            {
                "id": "front-centre", "x_m": 65.0,
                "thickness_m": 15.0, "y_min_m": 45.0,
                "y_max_m": 75.0, "height_m": 20.0,
            },
            {
                "id": "middle-south", "x_m": 90.0,
                "thickness_m": 20.0, "y_min_m": 15.0,
                "y_max_m": 40.0, "height_m": 15.0,
            },
            {
                "id": "middle-north", "x_m": 90.0,
                "thickness_m": 20.0, "y_min_m": 80.0,
                "y_max_m": 105.0, "height_m": 25.0,
            },
            {
                "id": "rear-centre", "x_m": 120.0,
                "thickness_m": 20.0, "y_min_m": 45.0,
                "y_max_m": 75.0, "height_m": 30.0,
            },
        )
        barrier = {
            "layout": barrier_layout,
            # The leading building remains the audit reference obstacle.
            **buildings[0],
            "buildings": buildings,
        }

    if barrier_layout == "building-cluster":
        obstacle_blocks = tuple(barrier["buildings"])
    else:
        segment_data = barrier.get("segments")
        segments = (
            tuple(segment_data)
            if segment_data is not None
            else ({
                "y_min_m": barrier["y_min_m"],
                "y_max_m": barrier["y_max_m"],
            },)
        )
        obstacle_blocks = tuple({
            "id": f"protective-wall-{index + 1}",
            "x_m": barrier["x_m"],
            "thickness_m": barrier["thickness_m"],
            "y_min_m": segment["y_min_m"],
            "y_max_m": segment["y_max_m"],
            "height_m": barrier["height_m"],
        } for index, segment in enumerate(segments))
    obstacles = tuple(
        AxisAlignedCuboid(
            float(block["x_m"]),
            float(block["x_m"] + block["thickness_m"]),
            float(block["y_min_m"]),
            float(block["y_max_m"]),
            0.0,
            float(block["height_m"]),
            str(block["id"]),
        )
        for block in obstacle_blocks
    )

    config = ResolvedObstacleFlow3DConfig(
        domain_m=domain,
        cells=cells,
        wind_speed_m_s=wind_speed,
        velocity_reference_scale_m_s=wind_speed,
        kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=1,
        minimum_flow_steps=1,
        lattice_inlet_speed=0.03,
        duration_s=duration,
        scalar_time_step_s=0.20,
        scalar_diffusivity_m2_s=6.1e-5,
        source_schedule=_release_schedule(
            duration_s=duration,
            release_duration_s=release_duration,
            release_rate_kg_s=release_rate,
            layout=barrier_layout,
        ),
        source_position_m=source,
        source_sigma_m=5.0,
        source_sigma_xyz_m=(5.0, 5.0, 5.0),
        source_kernel="gaussian",
        source_momentum_velocity_m_s=(0.0, 0.0, 0.0),
        obstacles=obstacles,
        scalar_advection_scheme="muscl_van_leer",
        source_time_integration="strang",
        subgrid_model="smagorinsky",
        smagorinsky_constant=0.12,
        nonstreamwise_wall_boundary=(
            "ground_no_slip_lateral_top_free_slip"
        ),
        thermal_buoyancy_coupling=True,
        source_temperature_k=293.15,
        thermal_diffusivity_m2_s=2.2e-5,
        effective_schmidt_number=0.70,
        effective_prandtl_number=0.85,
        scalar_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        thermal_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        maximum_active_scalar_density_departure=0.10,
        low_mach_mass_flux_correction=True,
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        active_provisional_face_projection="variable_density_only",
        variable_density_projection_linear_solver="auto",
        thermodynamic_corrector_maximum_iterations=4,
        thermodynamic_corrector_minimum_iterations=2,
        thermodynamic_corrector_density_tolerance=1.0e-6,
        thermodynamic_corrector_velocity_tolerance=1.0e-6,
        require_thermodynamic_corrector_convergence=True,
        store_scalar_fields=bool(store_scalar_fields),
    )
    return config, barrier


def audit_continuous_barrier_physics(
    result: ResolvedObstacleFlow3DResult,
    barrier: dict[str, Any],
) -> dict[str, Any]:
    """Audit flow mechanisms that the old animation solver could not model."""
    if not isinstance(result, ResolvedObstacleFlow3DResult):
        raise TypeError("result must be ResolvedObstacleFlow3DResult")
    x, y, z = result.x_m, result.y_m, result.z_m
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    solid = np.asarray(result.solid_mask, dtype=bool)
    fluid = ~solid
    dx = float(x[1] - x[0])
    x0 = float(barrier["x_m"])
    x1 = x0 + float(barrier["thickness_m"])
    y0 = float(barrier["y_min_m"])
    y1 = float(barrier["y_max_m"])
    height = float(barrier["height_m"])
    lateral_margin = dx

    reference = (
        (xx >= x0 - 7.0 * dx) & (xx < x0 - 4.0 * dx)
        & (yy >= y0 + lateral_margin) & (yy < y1 - lateral_margin)
        & (zz < height) & fluid
    )
    windward = (
        (xx >= x0 - 2.0 * dx) & (xx < x0)
        & (yy >= y0 + lateral_margin) & (yy < y1 - lateral_margin)
        & (zz < height) & fluid
    )
    roof = (
        (xx >= x0 - dx) & (xx < x1 + 2.0 * dx)
        & (yy >= y0 + lateral_margin) & (yy < y1 - lateral_margin)
        & (zz >= height) & (zz < height + 2.0 * dx) & fluid
    )
    lee = (
        (xx >= x1) & (xx < x1 + 6.0 * dx)
        & (yy >= y0 + lateral_margin) & (yy < y1 - lateral_margin)
        & (zz < height) & fluid
    )
    for name, mask in (
        ("reference", reference), ("windward", windward),
        ("roof", roof), ("lee", lee),
    ):
        if not np.any(mask):
            raise ValueError(f"{name} audit region contains no fluid cells")

    u = np.asarray(result.velocity_u_m_s, dtype=float)
    w = np.asarray(result.velocity_w_m_s, dtype=float)
    pressure = result.pressure_correction_pa
    if pressure is None:
        pressure = result.pressure_pa_gauge
    pressure = np.asarray(pressure, dtype=float)
    reference_u = float(np.mean(u[reference]))
    windward_u = float(np.mean(u[windward]))
    slowdown = 1.0 - windward_u / max(reference_u, 1.0e-12)
    windward_pressure = float(np.mean(pressure[windward]))
    lee_pressure = float(np.mean(pressure[lee]))
    roof_u_max = float(np.max(u[roof]))
    roof_vertical_shear_max = float(np.max(np.abs(
        np.gradient(u, dx, axis=0)[roof]
    )))
    lee_u_min = float(np.min(u[lee]))
    lee_w_min = float(np.min(w[lee]))
    lee_reverse_cells = int(np.count_nonzero(u[lee] < 0.0))
    lee_downwash_cells = int(np.count_nonzero(w[lee] < 0.0))

    blocked_x = np.pad(
        solid[:, :, :-1] | solid[:, :, 1:],
        ((0, 0), (0, 0), (1, 1)),
        constant_values=False,
    )
    blocked_y = np.pad(
        solid[:, :-1, :] | solid[:, 1:, :],
        ((0, 0), (1, 1), (0, 0)),
        constant_values=False,
    )
    blocked_z = np.pad(
        solid[:-1, :, :] | solid[1:, :, :],
        ((1, 1), (0, 0), (0, 0)),
        constant_values=False,
    )
    maximum_wall_normal_velocity = max(
        float(np.max(np.abs(result.u_face_m_s[blocked_x]))),
        float(np.max(np.abs(result.v_face_m_s[blocked_y]))),
        float(np.max(np.abs(result.w_face_m_s[blocked_z]))),
    )
    solid_values = np.asarray(result.concentration_kg_m3[solid], dtype=float)
    finite_solid_values = solid_values[np.isfinite(solid_values)]
    solid_concentration = (
        0.0 if finite_solid_values.size == 0
        else float(np.max(np.abs(finite_solid_values)))
    )
    injected = max(float(result.injected_mass_kg[-1]), 1.0e-12)
    relative_mass_residual = result.maximum_mass_residual_kg / injected
    continuity_residual = float(result.diagnostics.get(
        "maximum_variable_density_continuity_residual_kg_m3_s", math.inf,
    ))

    # Resolve the downstream scalar transport at the first x face behind the
    # wall into physically distinct paths.  A zero blocked-face velocity only
    # proves solid impermeability; it does not explain hydrogen found in the
    # lee.  The route audit therefore distinguishes direct wall-projection
    # transport from roof overtopping and lateral/open-gap bypass.
    x_face_index = int(round(x1 / dx))
    if not (1 <= x_face_index < x.size):
        raise ValueError("barrier downstream face is outside the x grid")
    dy = float(y[1] - y[0])
    dz = float(z[1] - z[0])
    concentration = np.asarray(result.concentration_kg_m3, dtype=float)
    left_concentration = concentration[:, :, x_face_index - 1]
    right_concentration = concentration[:, :, x_face_index]
    face_velocity = np.asarray(
        result.u_face_m_s[:, :, x_face_index], dtype=float,
    )
    donor_concentration = np.where(
        face_velocity >= 0.0, left_concentration, right_concentration,
    )
    effective_diffusivity = (
        6.1e-5 + np.asarray(result.eddy_viscosity_m2_s, dtype=float) / 0.70
    )
    face_diffusivity = 0.5 * (
        effective_diffusivity[:, :, x_face_index - 1]
        + effective_diffusivity[:, :, x_face_index]
    )
    scalar_flux_density = (
        face_velocity * donor_concentration
        - face_diffusivity
        * (right_concentration - left_concentration)
        / dx
    )
    blocked_face = (
        solid[:, :, x_face_index - 1] | solid[:, :, x_face_index]
    )
    scalar_flux_density[blocked_face] = 0.0
    downstream_flux = np.maximum(scalar_flux_density, 0.0) * dy * dz

    face_z, face_y = np.meshgrid(z, y, indexing="ij")
    projected_wall = np.zeros_like(face_y, dtype=bool)
    segments = barrier.get("segments") or ({
        "y_min_m": y0, "y_max_m": y1,
    },)
    for segment in segments:
        projected_wall |= (
            (face_y >= float(segment["y_min_m"]))
            & (face_y < float(segment["y_max_m"]))
        )
    direct_path = projected_wall & (face_z < height)
    overtop_path = projected_wall & (face_z >= height)
    lateral_or_gap_path = ~projected_wall
    route_flux = {
        "direct_through_wall_projection_kg_s": float(
            np.sum(downstream_flux[direct_path])
        ),
        "over_wall_top_kg_s": float(np.sum(downstream_flux[overtop_path])),
        "around_sides_or_open_gap_kg_s": float(
            np.sum(downstream_flux[lateral_or_gap_path])
        ),
    }
    route_flux["total_downstream_kg_s"] = float(sum(route_flux.values()))
    lee_mask = (xx >= x1) & fluid
    lee_inventory = float(np.sum(concentration[lee_mask]) * dx * dy * dz)

    checks = {
        "wall_impermeability": maximum_wall_normal_velocity <= 1.0e-10,
        "solid_scalar_exclusion": solid_concentration <= 1.0e-14,
        "mass_conservation": relative_mass_residual <= 1.0e-9,
        "pressure_projection": continuity_residual <= 1.0e-8,
        "windward_stagnation": slowdown >= 0.15,
        "windward_pressure_rise": windward_pressure > lee_pressure,
        "roof_acceleration_and_shear": (
            roof_u_max > 1.05 * reference_u
            and roof_vertical_shear_max > 0.02
        ),
        "lee_recirculation": lee_u_min < -0.02 and lee_reverse_cells > 0,
        "lee_downwash": lee_w_min < -0.02 and lee_downwash_cells > 0,
        "no_direct_scalar_transport_through_wall_projection": (
            route_flux["direct_through_wall_projection_kg_s"] <= 1.0e-14
        ),
    }
    return {
        "schema": "degali.continuous-barrier-physics-audit.v1",
        "all_internal_physics_checks_passed": bool(all(checks.values())),
        "checks": checks,
        "metrics": {
            "reference_streamwise_velocity_m_s": reference_u,
            "windward_streamwise_velocity_m_s": windward_u,
            "windward_slowdown_fraction": slowdown,
            "windward_pressure_correction_pa": windward_pressure,
            "lee_pressure_correction_pa": lee_pressure,
            "windward_to_lee_pressure_gain_pa": (
                windward_pressure - lee_pressure
            ),
            "roof_maximum_streamwise_velocity_m_s": roof_u_max,
            "roof_maximum_vertical_streamwise_shear_s_inv": (
                roof_vertical_shear_max
            ),
            "lee_minimum_streamwise_velocity_m_s": lee_u_min,
            "lee_reverse_flow_cell_count": lee_reverse_cells,
            "lee_minimum_vertical_velocity_m_s": lee_w_min,
            "lee_downwash_cell_count": lee_downwash_cells,
            "maximum_wall_normal_velocity_m_s": maximum_wall_normal_velocity,
            "maximum_solid_concentration_kg_m3": solid_concentration,
            "maximum_relative_mass_residual": relative_mass_residual,
            "maximum_continuity_residual_kg_m3_s": continuity_residual,
            "leeward_hydrogen_inventory_kg": lee_inventory,
        },
        "final_downstream_transport_routes": {
            "plane_x_m": x1,
            **route_flux,
            "interpretation": (
                "Hydrogen behind the wall is acceptable only when it is "
                "accounted for by the over-top or lateral/open-gap routes; "
                "the direct wall-projection route must remain zero."
            ),
        },
        "scope": (
            "internal mechanism verification on a coarse resolved grid; "
            "external concentration validation is still required for "
            "quantitative safety or digital-twin promotion"
        ),
    }


__all__ = [
    "audit_continuous_barrier_physics",
    "build_continuous_barrier_resolved_config",
]
