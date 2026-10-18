from dataclasses import replace

import numpy as np
import pytest

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.resolved_obstacle_flow_3d import (
    RESOLVED_OBSTACLE_FLOW_3D_SCHEMA,
    ResolvedInletVelocityHistory3D,
    ResolvedSourceMomentumHistory3D,
    ResolvedObstacleFlow3DConfig,
    _C3,
    _active_inlet_plane,
    _active_inlet_velocity_planes,
    _apply_outer_boundaries_d3q19,
    _boundary_acoustic_fraction,
    _circular_arc_blocked_faces,
    _convective_boundary_distribution,
    _equilibrium_d3q19,
    _effective_transport_diffusivity,
    _increment_projected_faces_from_cell_delta,
    _limited_boundary_extrapolation,
    _multidimensional_characteristic_state,
    _nscbc_xmax_boundary_distribution,
    _obstacle_voxelization_diagnostics,
    _restart_configuration_sha256,
    _source_kernel_fwhm_cells,
    _solid_mask,
    _stream_d3q19,
    _circular_arc_blocked_links,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import (
    AxisAlignedCuboid, CircularArcWall, CrosswindCylinder,
)


def _config(**changes):
    values = {
        "domain_m": (4.0, 2.4, 2.0),
        "cells": (40, 24, 20),
        "wind_speed_m_s": 2.0,
        "kinematic_viscosity_m2_s": 0.5,
        "flow_steps": 500,
        "minimum_flow_steps": 500,
        "lattice_inlet_speed": 0.03,
        "duration_s": 0.1,
        "scalar_time_step_s": 0.0005,
        "scalar_diffusivity_m2_s": 0.03,
        "source_schedule": SourceRateSchedule(
            (0.0, 0.05, 0.1), (0.02, 0.0, 0.0),
        ),
        "source_position_m": (0.5, 1.2, 0.4),
        "source_sigma_m": 0.15,
        "obstacles": (
            AxisAlignedCuboid(1.5, 2.0, 0.8, 1.6, 0.0, 1.0, "finite block"),
        ),
        "smagorinsky_constant": 0.10,
    }
    values.update(changes)
    return ResolvedObstacleFlow3DConfig(**values)


@pytest.fixture(scope="module")
def resolved_3d_result():
    return solve_resolved_obstacle_flow_3d(
        _config(),
        receptors=(
            FastObstacleReceptor3D("windward", (1.2, 1.2, 0.5)),
            FastObstacleReceptor3D("side", (1.8, 0.6, 0.5)),
            FastObstacleReceptor3D("roof", (1.8, 1.2, 1.2)),
            FastObstacleReceptor3D("lee", (2.4, 1.2, 0.5)),
        ),
    )


def test_resolved_3d_generates_finite_width_viscous_wake(resolved_3d_result):
    result = resolved_3d_result
    diagnostics = result.diagnostics

    assert result.schema == RESOLVED_OBSTACLE_FLOW_3D_SCHEMA
    assert np.all(result.velocity_u_m_s[result.solid_mask] == 0.0)
    assert np.all(result.velocity_v_m_s[result.solid_mask] == 0.0)
    assert np.all(result.velocity_w_m_s[result.solid_mask] == 0.0)
    assert diagnostics["finite_lateral_bypass"] is True
    assert diagnostics["recirculating_lee_cell_count"] > 0
    assert diagnostics["minimum_lee_streamwise_velocity_m_s"] < 0.0
    assert diagnostics["maximum_relative_lattice_density_variation"] < 0.05
    assert diagnostics["maximum_absolute_divergence_s_1"] < 1.0e-7
    assert diagnostics["maximum_absolute_raw_face_divergence_s_1"] > 0.0
    assert diagnostics["obstacle_grid_alignment_exact"] is True
    assert diagnostics["maximum_obstacle_face_error_cells"] < 1.0e-10


def test_source_kernel_resolution_diagnostic_distinguishes_underresolved_width():
    narrow = _source_kernel_fwhm_cells((0.001, 0.001, 0.001), 0.025)
    grid_scaled = _source_kernel_fwhm_cells((0.025, 0.025, 0.025), 0.025)

    assert narrow == pytest.approx((0.09419, 0.09419, 0.09419), rel=2.0e-4)
    assert min(narrow) < 2.0
    assert grid_scaled == pytest.approx((2.35482, 2.35482, 2.35482), rel=2.0e-4)
    assert min(grid_scaled) >= 2.0


def test_resolved_solver_warns_when_source_kernel_is_subgrid():
    result = solve_resolved_obstacle_flow_3d(_config(
        flow_steps=12,
        minimum_flow_steps=12,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        source_sigma_m=0.001,
    ))

    assert result.diagnostics["source_kernel_underresolved"] is True
    assert result.diagnostics["source_kernel_fwhm_cell_count_xyz"] == pytest.approx(
        [0.02355, 0.02355, 0.02355], rel=3.0e-4,
    )
    assert any("source Gaussian FWHM is below two cells" in warning
               for warning in result.diagnostics["warnings"])


def test_resolved_solver_initializes_and_inflows_measured_temperature_profile():
    config = _config(
        flow_steps=12,
        minimum_flow_steps=12,
        duration_s=0.003,
        scalar_time_step_s=0.0015,
        source_schedule=SourceRateSchedule((0.0, 0.003), (0.0, 0.0)),
        thermal_buoyancy_coupling=True,
        source_temperature_k=291.15,
        inlet_temperature_profile_z_m=(0.0, 1.0, 2.0),
        inlet_temperature_profile_k=(294.0, 291.15, 288.0),
        obstacle_wall_temperature_k=297.6,
        ground_temperature_k=318.95,
    )

    result = solve_resolved_obstacle_flow_3d(config)

    assert result.temperature_k is not None
    assert np.nanmin(result.temperature_k) < 291.15
    assert np.nanmax(result.temperature_k) > 291.15
    assert result.diagnostics["energy_equation"] is True
    assert result.diagnostics["maximum_relative_energy_residual"] < 1.0e-10
    assert result.diagnostics["initial_temperature_profile_applied"] is True
    assert result.diagnostics["inlet_temperature_boundary_model"].startswith(
        "conservative_x_min_sensible_enthalpy_inflow"
    )
    assert result.diagnostics["isothermal_wall_heat_boundary_model"] == (
        "first_order_half_cell_conduction_to_fluid_cells"
    )
    assert result.diagnostics["cumulative_wall_heat_added_j"] > 0.0


def test_resolved_solver_applies_conserved_state_and_momentum_callbacks_before_transport():
    state_calls = []
    momentum_calls = []

    def state_callback(concentration, sensible_enthalpy, start, dt):
        state_calls.append((start, dt, sensible_enthalpy is not None))
        return concentration, sensible_enthalpy

    def momentum_callback(u, v, w, u_face, v_face, w_face, start, dt):
        momentum_calls.append((start, dt))
        return u, v, w, u_face, v_face, w_face

    result = solve_resolved_obstacle_flow_3d(
        _config(
            flow_steps=12,
            minimum_flow_steps=12,
            duration_s=0.003,
            scalar_time_step_s=0.0015,
            source_schedule=SourceRateSchedule((0.0, 0.003), (0.0, 0.0)),
            thermal_buoyancy_coupling=True,
            source_temperature_k=291.15,
        ),
        interface_conserved_state_callback=state_callback,
        interface_momentum_callback=momentum_callback,
    )

    assert state_calls
    assert len(state_calls) == len(momentum_calls)
    assert all(item[2] for item in state_calls)
    assert result.diagnostics["interface_conserved_state_callback_applied"] is True
    assert result.diagnostics["interface_momentum_callback_applied"] is True


def test_resolved_solver_applies_opt_in_stratification_diffusion_candidate():
    config = _config(
        flow_steps=12,
        minimum_flow_steps=12,
        duration_s=0.003,
        scalar_time_step_s=0.0015,
        source_schedule=SourceRateSchedule((0.0, 0.003), (0.0, 0.0)),
        thermal_buoyancy_coupling=True,
        source_temperature_k=291.15,
        inlet_temperature_profile_z_m=(0.0, 1.0, 2.0),
        inlet_temperature_profile_k=(291.14, 291.15, 291.16),
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.7,
        stratification_damping_alpha=5.0,
        stratification_shear_regularization_s_inv=0.041,
    )

    result = solve_resolved_obstacle_flow_3d(config)

    assert result.diagnostics["stratification_damping_alpha"] == 5.0
    assert result.diagnostics["stratification_shear_regularization_s_inv"] == 0.041
    assert result.diagnostics["scalar_diffusivity_anisotropy"].startswith("horizontal")
    assert result.diagnostics["maximum_gradient_richardson_number"] > 0.0
    assert result.diagnostics["maximum_relative_energy_residual"] < 1.0e-10
    assert any("cross-domain sensitivity candidate" in warning for warning in result.diagnostics["warnings"])


def test_stratification_candidate_parameters_are_restart_hash_pinned():
    config = _config(
        thermal_buoyancy_coupling=True,
        source_temperature_k=291.15,
        effective_schmidt_number=0.7,
        stratification_damping_alpha=5.0,
        stratification_shear_regularization_s_inv=0.041,
    )
    changed = replace(config, stratification_damping_alpha=6.0)

    assert _restart_configuration_sha256(config) != _restart_configuration_sha256(changed)


def test_resolved_3d_voxelization_is_half_open_and_reports_geometry_error():
    config = _config(obstacles=(
        AxisAlignedCuboid(1.55, 2.05, 0.8, 1.6, 0.0, 1.0, "centre-cut block"),
    ))
    x, y, z, solid = _solid_mask(config)
    x_indices = np.flatnonzero(np.any(solid, axis=(0, 1)))
    assert x_indices.size == 5
    records, face_error, size_error = _obstacle_voxelization_diagnostics(
        config, x, y, z
    )
    assert records[0]["effective_voxel_bounds_m"][:2] == pytest.approx((1.5, 2.0))
    assert records[0]["solid_cell_counts_xyz"][0] == 5
    assert face_error == pytest.approx(0.05)
    assert size_error < 1.0e-12


def test_resolved_3d_voxelizes_curved_wall_without_false_exact_alignment():
    config = _config(obstacles=(CircularArcWall(
        center_m=(2.0, 1.2), radius_m=0.5, thickness_m=0.2,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
        arc_end_deg=90.0, label="semicircular fence",
    ),))
    x, y, z, solid = _solid_mask(config)
    assert np.any(solid)
    records, face_error, _ = _obstacle_voxelization_diagnostics(
        config, x, y, z
    )
    assert records[0]["geometry_type"] == "CircularArcWall"
    assert records[0]["curved_surface_voxelization"] is True
    assert records[0]["solid_cell_count"] == int(np.count_nonzero(solid))
    assert face_error >= 0.05


def test_resolved_3d_accepts_round_crosswind_rods_for_obstacle_arrays():
    rod = CrosswindCylinder(
        center_x_m=2.0, center_z_m=0.5, radius_m=0.1,
        y_min_m=0.4, y_max_m=2.0, label="transverse-rod-1",
    )
    config = _config(obstacles=(rod,))
    x, y, z, solid = _solid_mask(config)
    records, _, _ = _obstacle_voxelization_diagnostics(config, x, y, z)

    assert np.any(solid)
    assert records[0]["geometry_type"] == "CrosswindCylinder"
    assert records[0]["solid_cell_count"] == int(np.count_nonzero(solid))


def test_d3q19_stream_bounces_back_on_exactly_blocked_arc_link():
    populations = np.arange(19 * 8 * 8 * 8, dtype=float).reshape(19, 8, 8, 8)
    solid = np.zeros((8, 8, 8), dtype=bool)
    fractions = np.full(populations.shape, np.nan)
    fractions[1, 3, 3, 3] = 0.25  # +x link; wall is one quarter along it
    streamed = _stream_d3q19(populations, solid, fractions)

    assert streamed[2, 3, 3, 3] == pytest.approx(
        0.5 * populations[1, 3, 3, 3]
        + 0.5 * populations[1, 3, 3, 2]
    )
    assert streamed[1, 3, 3, 4] == populations[1, 3, 3, 3]


def test_d3q19_curved_link_interpolation_uses_far_formula_and_safe_fallback():
    populations = np.arange(19 * 8 * 8 * 8, dtype=float).reshape(19, 8, 8, 8)
    solid = np.zeros((8, 8, 8), dtype=bool)
    fractions = np.full(populations.shape, np.nan)
    fractions[1, 3, 3, 3] = 0.75
    fractions[1, 4, 3, 3] = 0.25
    solid[4, 3, 2] = True  # invalid upstream interpolation sample for second link
    streamed = _stream_d3q19(populations, solid, fractions)

    q = 0.75
    expected_far = (
        populations[1, 3, 3, 3] / (2.0 * q)
        + ((2.0 * q - 1.0) / (2.0 * q)) * populations[2, 3, 3, 3]
    )
    assert streamed[2, 3, 3, 3] == pytest.approx(expected_far)
    # Upstream solid forces q=1/2, recovering ordinary halfway bounce-back.
    assert streamed[2, 4, 3, 3] == populations[1, 4, 3, 3]


def test_zero_thickness_arc_representation_uses_link_and_face_geometry_not_voxels():
    wall = CircularArcWall(
        center_m=(2.0, 1.2), radius_m=0.5, thickness_m=0.0,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
        arc_end_deg=90.0, label="idealized zero-thickness fence",
    )
    config = _config(
        obstacles=(wall,), circular_arc_representation="zero_thickness_surface",
    )
    x, y, z, solid = _solid_mask(config)
    links = _circular_arc_blocked_links(config, x, y, z)
    records, _, _ = _obstacle_voxelization_diagnostics(config, x, y, z)

    assert not np.any(solid)
    assert links is not None and np.count_nonzero(np.isfinite(links)) > 0
    assert records[0]["geometry_representation"] == "zero_thickness_surface"
    assert records[0]["declared_radial_thickness_m"] == pytest.approx(0.0)
    assert records[0]["solid_cell_count"] == 0


def test_zero_thickness_arc_cannot_be_silently_used_as_a_voxel_obstacle():
    wall = CircularArcWall(
        center_m=(2.0, 1.2), radius_m=0.5, thickness_m=0.0,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
        arc_end_deg=90.0,
    )
    with pytest.raises(ValueError, match="finite_thickness_voxel requires positive"):
        _config(obstacles=(wall,))


def test_resolved_d3q19_run_applies_selected_zero_thickness_arc_boundary():
    wall = CircularArcWall(
        center_m=(2.0, 1.2), radius_m=0.5, thickness_m=0.0,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
        arc_end_deg=90.0, label="idealized zero-thickness fence",
    )
    result = solve_resolved_obstacle_flow_3d(_config(
        obstacles=(wall,),
        circular_arc_representation="zero_thickness_surface",
        flow_steps=1, minimum_flow_steps=1,
        duration_s=0.0015, scalar_time_step_s=0.0015,
        source_schedule=SourceRateSchedule((0.0, 0.0015), (0.0, 0.0)),
    ))

    assert result.diagnostics["d3q19_zero_thickness_arc_link_bounce_back_applied"]
    assert result.diagnostics["d3q19_zero_thickness_arc_link_count"] > 0
    assert result.diagnostics["circular_arc_representation"] == "zero_thickness_surface"
    assert result.diagnostics["obstacle_grid_alignment_exact"] is None
    assert result.diagnostics["maximum_mass_residual_kg"] < 1.0e-12


def test_resolved_fv_arc_faces_remain_closed_in_flow_and_scalar_state():
    config = _config(
        obstacles=(CircularArcWall(
            center_m=(2.0, 1.2), radius_m=0.5, thickness_m=0.2,
            z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
            arc_end_deg=90.0, label="semicircular fence",
        ),),
        flow_steps=1,
        minimum_flow_steps=1,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        maximum_active_scalar_density_departure=1.0,
        active_momentum_solver="finite_volume_projection",
    )
    result = solve_resolved_obstacle_flow_3d(config)
    blocked = _circular_arc_blocked_faces(
        config, result.x_m, result.y_m, result.z_m,
    )

    assert result.diagnostics["embedded_arc_blocked_face_count"] > 0
    assert result.diagnostics["embedded_arc_face_masks_applied_to_scalar_fluxes"]
    assert not result.diagnostics["embedded_arc_fractional_cut_cell_areas"]
    assert np.all(result.u_face_m_s[blocked[0]] == 0.0)
    assert np.all(result.v_face_m_s[blocked[1]] == 0.0)
    assert np.all(result.w_face_m_s[blocked[2]] == 0.0)
    assert result.diagnostics[
        "maximum_variable_density_continuity_residual_kg_m3_s"
    ] < 1.0e-8
    assert result.diagnostics["maximum_mass_residual_kg"] < 1.0e-12
    assert result.diagnostics["maximum_energy_residual_j"] < 1.0e-12


def test_zero_thickness_arc_fv_projection_starts_without_voxel_solid_leakage():
    wall = CircularArcWall(
        center_m=(2.0, 1.2), radius_m=0.5, thickness_m=0.0,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
        arc_end_deg=90.0, label="idealized zero-thickness fence",
    )
    config = _config(
        obstacles=(wall,),
        circular_arc_representation="zero_thickness_surface",
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        thermal_buoyancy_coupling=True,
        source_temperature_k=293.15,
        flow_steps=1, minimum_flow_steps=1,
        duration_s=0.001, scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
    )
    result = solve_resolved_obstacle_flow_3d(config)
    blocked = _circular_arc_blocked_faces(config, result.x_m, result.y_m, result.z_m)

    assert not np.any(result.solid_mask)
    assert np.count_nonzero(blocked[0]) + np.count_nonzero(blocked[1]) > 0
    assert np.all(result.u_face_m_s[blocked[0]] == 0.0)
    assert np.all(result.v_face_m_s[blocked[1]] == 0.0)
    assert not result.diagnostics["d3q19_zero_thickness_arc_link_bounce_back_applied"]
    assert result.diagnostics["embedded_arc_tangential_no_slip_stress_applied"] is True
    assert "smooth-wall Reichardt" in result.diagnostics[
        "embedded_arc_momentum_boundary_model"
    ]
    assert result.diagnostics["embedded_arc_surface_area_m2"] == pytest.approx(
        np.pi * 0.5 * 1.0, rel=1.0e-12,
    )
    assert result.diagnostics["embedded_arc_domain_area_coverage_fraction"] == pytest.approx(
        1.0, abs=1.0e-12,
    )
    assert result.diagnostics["embedded_arc_surface_cell_count"] > 0
    assert np.isfinite(result.diagnostics["embedded_arc_minimum_wall_distance_m"])
    assert result.diagnostics["maximum_embedded_arc_wall_shear_stress_pa"] > 0.0
    assert result.diagnostics["maximum_embedded_arc_wall_friction_velocity_m_s"] > 0.0
    assert result.diagnostics["maximum_momentum_balance_residual_kg_m_s"] < 1.0e-12
    assert result.diagnostics["maximum_mass_residual_kg"] < 1.0e-12


def test_zero_ambient_wind_fv_mode_keeps_closed_inlet_and_generates_buoyant_h2():
    result = solve_resolved_obstacle_flow_3d(_config(
        wind_speed_m_s=0.0,
        velocity_reference_scale_m_s=2.0,
        inlet_profile_z_m=(0.0, 2.0),
        inlet_profile_u_m_s=(0.0, 0.0),
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        thermal_buoyancy_coupling=True,
        source_temperature_k=293.15,
        maximum_active_scalar_density_departure=1.0,
        flow_steps=1, minimum_flow_steps=1,
        duration_s=0.001, scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.01, 0.0)),
    ))

    assert result.diagnostics["ambient_wind_speed_m_s"] == 0.0
    assert result.diagnostics["velocity_reference_scale_m_s"] == 2.0
    assert result.diagnostics["reference_dynamic_pressure_pa"] == pytest.approx(2.408)
    assert result.diagnostics["zero_ambient_wind_mode"] is True
    assert result.diagnostics["active_inlet_plane_minimum_m_s"] == 0.0
    assert result.diagnostics["active_inlet_plane_maximum_m_s"] == 0.0
    assert result.u_face_m_s[:, :, 0] == pytest.approx(0.0)
    assert np.all(np.isfinite(result.pressure_coefficient[~result.solid_mask]))
    assert result.injected_mass_kg[-1] > 0.0
    assert np.max(
        result.buoyancy_acceleration_m_s2[~result.solid_mask]
    ) > 0.0
    assert result.diagnostics[
        "maximum_variable_density_continuity_residual_kg_m3_s"
    ] < 1.0e-8
    assert result.diagnostics["maximum_mass_residual_kg"] < 1.0e-12
    assert result.diagnostics["maximum_energy_residual_j"] < 1.0e-12


def test_zero_ambient_wind_requires_reference_scale_and_projected_fv():
    with pytest.raises(ValueError, match="require velocity_reference_scale_m_s"):
        _config(wind_speed_m_s=0.0)
    with pytest.raises(ValueError, match="requires projected_inlet initialization"):
        _config(
            wind_speed_m_s=0.0,
            velocity_reference_scale_m_s=2.0,
        )


def test_resolved_3d_supports_an_open_domain_reference_case():
    result = solve_resolved_obstacle_flow_3d(_config(
        obstacles=(),
        circular_arc_representation="zero_thickness_surface",
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        effective_schmidt_number=0.8,
        effective_prandtl_number=0.9,
    ))

    assert not np.any(result.solid_mask)
    assert result.diagnostics["embedded_arc_tangential_no_slip_stress_applied"] is True
    assert result.diagnostics["embedded_arc_momentum_boundary_model"].startswith(
        "not applicable"
    )
    assert result.diagnostics["obstacle_count"] == 0
    assert result.diagnostics["finite_lateral_bypass"] is False
    assert result.diagnostics["obstacle_voxelization"] == []
    assert result.diagnostics["recirculating_lee_cell_count"] == 0
    assert result.diagnostics["minimum_lee_streamwise_velocity_m_s"] is None
    assert result.diagnostics["applied_scalar_diffusivity_m2_s"] == pytest.approx(
        0.5 / 0.8
    )
    assert result.diagnostics["applied_thermal_diffusivity_m2_s"] == pytest.approx(
        0.5 / 0.9
    )
    assert "effective Schmidt" in result.diagnostics["warnings"][-2]


def test_incremental_face_map_is_exact_identity_for_zero_momentum_delta():
    shape = (4, 5, 6)
    rng = np.random.default_rng(731)
    solid = np.zeros(shape, dtype=bool)
    solid[1:3, 2:4, 3:5] = True
    u_face = rng.normal(size=(4, 5, 7))
    v_face = rng.normal(size=(4, 6, 6))
    w_face = rng.normal(size=(5, 5, 6))
    u_face[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    v_face[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    w_face[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    v_face[:, (0, -1), :] = 0.0
    w_face[(0, -1), :, :] = 0.0
    zero = np.zeros(shape)

    result = _increment_projected_faces_from_cell_delta(
        previous_u_face=u_face,
        previous_v_face=v_face,
        previous_w_face=w_face,
        delta_u=zero,
        delta_v=zero,
        delta_w=zero,
        solid=solid,
        spacing_m=0.2,
    )

    assert result[0] == pytest.approx(u_face)
    assert result[1] == pytest.approx(v_face)
    assert result[2] == pytest.approx(w_face)


def test_incremental_face_map_updates_pressure_open_lateral_faces():
    shape = (4, 5, 6)
    solid = np.zeros(shape, dtype=bool)
    zero = np.zeros(shape)
    delta_v = np.zeros(shape)
    delta_w = np.zeros(shape)
    delta_v[:, 0, :] = -0.03
    delta_v[:, -1, :] = 0.04
    delta_w[0, :, :] = -0.05
    delta_w[-1, :, :] = 0.06

    result = _increment_projected_faces_from_cell_delta(
        previous_u_face=np.zeros((4, 5, 7)),
        previous_v_face=np.zeros((4, 6, 6)),
        previous_w_face=np.zeros((5, 5, 6)),
        delta_u=zero,
        delta_v=delta_v,
        delta_w=delta_w,
        solid=solid,
        spacing_m=0.2,
        nonstreamwise_wall_boundary="all_pressure_open",
    )

    assert result[1][:, 0, :] == pytest.approx(-0.03)
    assert result[1][:, -1, :] == pytest.approx(0.04)
    assert result[2][0, :, :] == pytest.approx(-0.05)
    assert result[2][-1, :, :] == pytest.approx(0.06)


def test_resolved_3d_scalar_is_impermeable_and_mass_conservative(
    resolved_3d_result,
):
    result = resolved_3d_result
    solid = result.solid_mask
    blocked_x = solid[:, :, :-1] | solid[:, :, 1:]
    blocked_y = solid[:, :-1, :] | solid[:, 1:, :]
    blocked_z = solid[:-1, :, :] | solid[1:, :, :]

    assert np.all(result.u_face_m_s[:, :, 1:-1][blocked_x] == 0.0)
    assert np.all(result.v_face_m_s[:, 1:-1, :][blocked_y] == 0.0)
    assert np.all(result.w_face_m_s[1:-1, :, :][blocked_z] == 0.0)
    assert result.injected_mass_kg[-1] == pytest.approx(0.001)
    assert result.maximum_mass_residual_kg < 1.0e-12
    assert np.nanmax(result.maximum_concentration_vol_percent) > 0.0
    assert set(result.receptor_traces_kg_m3) == {
        "windward", "side", "roof", "lee",
    }


def test_resolved_3d_reports_unvalidated_closures_and_runtime(resolved_3d_result):
    diagnostics = resolved_3d_result.diagnostics

    assert diagnostics["quantitative_prediction_allowed"] is False
    assert diagnostics["nonlocal_scalar_routing"] is False
    assert diagnostics["local_conservative_species_transport"] is True
    assert diagnostics["turbulence_model"].startswith("explicit Smagorinsky")
    assert diagnostics["smagorinsky_constant"] == pytest.approx(0.10)
    assert diagnostics["wall_damping"] is False
    assert diagnostics["energy_equation"] is False
    assert diagnostics["density_buoyancy_feedback"] is False
    assert diagnostics["scalar_advection_scheme"] == "muscl_van_leer"
    assert diagnostics["maximum_eddy_viscosity_m2_s"] > 0.0
    assert diagnostics["total_runtime_s"] >= diagnostics["flow_solve_runtime_s"]


def test_resolved_3d_requires_cubic_cells_and_complete_inlet_profile():
    with pytest.raises(ValueError, match="cubic cells"):
        _config(cells=(40, 20, 20))
    with pytest.raises(ValueError, match="supplied together"):
        _config(inlet_profile_z_m=(0.0, 2.0))
    with pytest.raises(ValueError, match="cannot exceed"):
        _config(maximum_lattice_mach_number=0.31)
    with pytest.raises(ValueError, match="scalar_advection_scheme"):
        _config(scalar_advection_scheme="unbounded")
    with pytest.raises(ValueError, match="subgrid_model"):
        _config(subgrid_model="invented")
    with pytest.raises(ValueError, match="wale_constant"):
        _config(wale_constant=1.1)
    with pytest.raises(ValueError, match="flow_boundary_scheme"):
        _config(flow_boundary_scheme="unverified")
    assert _config(
        flow_boundary_scheme="characteristic_multidimensional_extrapolated",
    ).flow_boundary_scheme == "characteristic_multidimensional_extrapolated"
    assert _config(
        flow_boundary_scheme="characteristic_multidimensional_vortical_limited",
    ).flow_boundary_scheme == "characteristic_multidimensional_vortical_limited"
    assert _config(
        flow_boundary_scheme="characteristic_nscbc_convective",
    ).flow_boundary_scheme == "characteristic_nscbc_convective"
    with pytest.raises(ValueError, match="nonstreamwise_wall_boundary"):
        _config(nonstreamwise_wall_boundary="frictionless_tunnel")
    with pytest.raises(ValueError, match="effective_schmidt_number"):
        _config(effective_schmidt_number=0.0)
    with pytest.raises(ValueError, match="effective_prandtl_number"):
        _config(effective_prandtl_number=-0.1)
    with pytest.raises(ValueError, match="scalar_turbulent_diffusivity_model"):
        _config(scalar_turbulent_diffusivity_model="fitted_constant")
    with pytest.raises(ValueError, match="thermal_turbulent_diffusivity_model"):
        _config(thermal_turbulent_diffusivity_model="fitted_constant")
    with pytest.raises(ValueError, match="source_time_integration"):
        _config(source_time_integration="future")
    with pytest.raises(ValueError, match="active_provisional_face_projection"):
        _config(active_provisional_face_projection="double_projection")
    with pytest.raises(ValueError, match="active_momentum_solver"):
        _config(active_momentum_solver="invented")
    with pytest.raises(ValueError, match="requires thermal coupling"):
        _config(active_momentum_solver="finite_volume_projection")
    with pytest.raises(TypeError, match="three-value tuple"):
        _config(source_momentum_velocity_m_s=(1.0, 0.0))
    with pytest.raises(ValueError, match="requires the finite_volume_projection"):
        _config(source_momentum_velocity_m_s=(1.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="supplied together"):
        _config(active_inlet_profile_y_m=(0.0, 2.4))
    with pytest.raises(ValueError, match="strictly increasing"):
        _config(
            active_inlet_profile_y_m=(0.0, 1.0, 0.5),
            active_inlet_profile_u_m_s=(2.0, 2.1, 2.2),
            thermal_buoyancy_coupling=True,
            source_temperature_k=80.0,
            low_mach_mass_flux_correction=True,
            active_momentum_solver="finite_volume_projection",
        )
    with pytest.raises(ValueError, match="active y-dependent inlet profiles"):
        _config(
            active_inlet_profile_y_m=(0.0, 2.4),
            active_inlet_profile_u_m_s=(1.8, 2.2),
        )
    with pytest.raises(ValueError, match="pressure_correction_momentum_coupling"):
        _config(pressure_correction_momentum_coupling="unbounded_impulse")
    with pytest.raises(ValueError, match="variable_density_projection_linear_solver"):
        _config(variable_density_projection_linear_solver="multigrid")
    assert _config(gravity_m_s2=0.0).gravity_m_s2 == 0.0
    with pytest.raises(ValueError, match="gravity_m_s2"):
        _config(gravity_m_s2=-1.0)
    with pytest.raises(ValueError, match="outlet_lattice_density"):
        _config(outlet_lattice_density=0.0)
    with pytest.raises(ValueError, match="zero or at least two"):
        _config(flow_average_last_steps=1)
    with pytest.raises(ValueError, match="less than flow_steps"):
        _config(flow_average_last_steps=500)
    with pytest.raises(ValueError, match="minimum_flow_steps to equal"):
        _config(flow_average_last_steps=100, minimum_flow_steps=400)


def test_sgs_diffusivity_separates_molecular_and_local_eddy_terms():
    eddy = np.asarray([[[0.0, 0.02], [0.04, 0.08]]])
    physical = _effective_transport_diffusivity(
        molecular_diffusivity_m2_s=1.0e-4,
        effective_dimensionless_number=0.5,
        base_kinematic_viscosity_m2_s=0.5,
        subgrid_eddy_viscosity_m2_s=eddy,
        model="molecular_plus_subgrid_eddy_viscosity",
    )
    legacy = _effective_transport_diffusivity(
        molecular_diffusivity_m2_s=1.0e-4,
        effective_dimensionless_number=0.5,
        base_kinematic_viscosity_m2_s=0.5,
        subgrid_eddy_viscosity_m2_s=eddy,
        model="legacy_total_kinematic_viscosity_ratio",
    )

    assert physical == pytest.approx(1.0e-4 + eddy / 0.5)
    assert legacy == pytest.approx(1.0)


def test_active_inlet_plane_combines_declared_y_and_z_mean_profiles():
    config = _config(
        wind_speed_m_s=2.0,
        active_inlet_profile_y_m=(0.0, 1.0, 2.0),
        active_inlet_profile_u_m_s=(1.0, 2.0, 3.0),
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        low_mach_mass_flux_correction=True,
        active_momentum_solver="finite_volume_projection",
    )
    y = np.asarray((0.0, 0.5, 1.0, 1.5, 2.0))
    z_profile = np.asarray((1.0, 2.0, 3.0))

    plane = _active_inlet_plane(config, y, z_profile)

    expected_y = np.asarray((1.0, 1.5, 2.0, 2.5, 3.0))
    assert plane == pytest.approx(
        z_profile[:, None] * expected_y[None, :] / 2.0
    )


def test_vector_inlet_history_interpolates_and_maps_declared_spatial_shape():
    history = ResolvedInletVelocityHistory3D(
        time_s=(0.0, 1.0),
        u_m_s=(1.0, 3.0),
        v_m_s=(-0.2, 0.2),
        w_m_s=(0.0, 0.4),
        history_id="measured-vector-v1",
    )
    assert history.velocity_at(0.5) == pytest.approx((2.0, 0.0, 0.2))
    config = _config(
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        active_momentum_solver="finite_volume_projection",
        active_inlet_velocity_history=history,
    )
    reference = np.asarray([[1.0, 2.0], [0.5, 1.5]])
    u_plane, v_plane, w_plane = _active_inlet_velocity_planes(
        config, reference, 0.5,
    )
    assert u_plane == pytest.approx(reference)
    assert v_plane == pytest.approx(np.zeros_like(reference))
    assert w_plane == pytest.approx(0.1 * reference)

    block_history = ResolvedInletVelocityHistory3D(
        time_s=(0.0, 1.0, 2.0),
        u_m_s=(1.0, 2.0, 2.0),
        v_m_s=(0.1, 0.2, 0.2),
        w_m_s=(0.0, -0.1, -0.1),
        history_id="block-mean-v1",
        interpolation="piecewise_constant",
    )
    assert block_history.velocity_at(0.999) == pytest.approx((1.0, 0.1, 0.0))
    assert block_history.velocity_at(1.0) == pytest.approx((2.0, 0.2, -0.1))


def test_vector_inlet_history_requires_supported_active_solver_and_coverage():
    short = ResolvedInletVelocityHistory3D(
        time_s=(0.0, 0.05),
        u_m_s=(2.0, 2.0),
        v_m_s=(0.0, 0.0),
        w_m_s=(0.0, 0.0),
        history_id="short",
    )
    with pytest.raises(ValueError, match="cover duration_s"):
        _config(active_inlet_velocity_history=short)
    complete = ResolvedInletVelocityHistory3D(
        time_s=(0.0, 0.1),
        u_m_s=(2.0, 2.0),
        v_m_s=(0.0, 0.0),
        w_m_s=(0.0, 0.0),
        history_id="complete",
    )
    with pytest.raises(ValueError, match="finite_volume_projection"):
        _config(active_inlet_velocity_history=complete)


def test_source_momentum_history_interpolates_and_requires_projection_solver():
    history = ResolvedSourceMomentumHistory3D(
        time_s=(0.0, 0.1),
        u_m_s=(0.0, 0.0),
        v_m_s=(0.0, 0.0),
        w_m_s=(0.0, -20.0),
        history_id="source-jet-flow-scaled-v1",
    )
    assert history.velocity_at(0.05) == pytest.approx((0.0, 0.0, -10.0))
    with pytest.raises(ValueError, match="finite_volume_projection"):
        _config(source_momentum_velocity_history=history)
    config = _config(
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        source_momentum_velocity_history=history,
    )
    assert config.source_momentum_velocity_history is history


def test_nonequilibrium_extrapolation_enforces_open_boundary_moments():
    nz, ny, nx = 5, 6, 7
    zz, yy, xx = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij"
    )
    rho = 1.0 + 2.0e-4 * xx + 1.0e-4 * yy
    ux = 0.012 + 4.0e-4 * xx
    uy = 8.0e-4 * np.sin(yy)
    uz = 6.0e-4 * np.cos(zz)
    distribution = _equilibrium_d3q19(rho, ux, uy, uz)
    # Add a zero-mass, zero-momentum non-equilibrium stress perturbation.
    distribution[0] += 2.0e-5
    distribution[1] -= 1.0e-5
    distribution[2] -= 1.0e-5
    outlet_neighbor = distribution[:, :, :, -2].copy()
    outlet_neighbor_rho = np.sum(outlet_neighbor, axis=0)
    outlet_neighbor_u = np.stack([
        np.tensordot(_C3[:, component], outlet_neighbor, axes=(0, 0))
        / outlet_neighbor_rho
        for component in range(3)
    ])
    inlet_profile = np.linspace(0.0, 0.035, nz)

    _apply_outer_boundaries_d3q19(
        distribution,
        inlet_profile,
        boundary_scheme="nonequilibrium_extrapolation",
        outlet_lattice_density=1.003,
    )

    interior = np.s_[1:-1, 1:-1]
    inlet = distribution[:, :, :, 0]
    inlet_rho = np.sum(inlet, axis=0)
    inlet_velocity = np.stack([
        np.tensordot(_C3[:, component], inlet, axes=(0, 0)) / inlet_rho
        for component in range(3)
    ])
    assert inlet_velocity[0][interior] == pytest.approx(
        np.broadcast_to(inlet_profile[:, None], (nz, ny))[interior]
    )
    assert inlet_velocity[1][interior] == pytest.approx(0.0, abs=1.0e-14)
    assert inlet_velocity[2][interior] == pytest.approx(0.0, abs=1.0e-14)

    outlet = distribution[:, :, :, -1]
    outlet_rho = np.sum(outlet, axis=0)
    outlet_velocity = np.stack([
        np.tensordot(_C3[:, component], outlet, axes=(0, 0)) / outlet_rho
        for component in range(3)
    ])
    assert outlet_rho[interior] == pytest.approx(1.003, abs=1.0e-14)
    assert outlet_velocity[:, 1:-1, 1:-1] == pytest.approx(
        outlet_neighbor_u[:, 1:-1, 1:-1], abs=1.0e-14
    )


def test_tunnel_wall_boundary_uses_full_not_specular_reflection():
    nz, ny, nx = 5, 6, 7
    rng = np.random.default_rng(1842)
    original = rng.uniform(0.001, 0.05, size=(19, nz, ny, nx))
    tunnel = original.copy()
    free = original.copy()
    inlet = np.full(nz, 0.02)
    _apply_outer_boundaries_d3q19(
        tunnel, inlet,
        boundary_scheme="equilibrium_copy",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary="all_no_slip_tunnel",
    )
    _apply_outer_boundaries_d3q19(
        free, inlet,
        boundary_scheme="equilibrium_copy",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary=(
            "ground_no_slip_lateral_top_free_slip"
        ),
    )
    direction = next(
        index for index, vector in enumerate(_C3)
        if tuple(vector) == (1, 0, -1)
    )
    opposite = next(
        index for index, vector in enumerate(_C3)
        if tuple(vector) == (-1, 0, 1)
    )
    specular = next(
        index for index, vector in enumerate(_C3)
        if tuple(vector) == (1, 0, 1)
    )
    point = (nz - 1, ny // 2, nx // 2)
    assert tunnel[(direction, *point)] == pytest.approx(
        original[(opposite, *point)]
    )
    assert free[(direction, *point)] == pytest.approx(
        original[(specular, *point)]
    )


def test_all_symmetry_boundary_removes_the_no_slip_lower_wall():
    nz, ny, nx = 5, 6, 7
    rng = np.random.default_rng(2841)
    original = rng.uniform(0.001, 0.05, size=(19, nz, ny, nx))
    symmetry = original.copy()
    ground = original.copy()
    inlet = np.full(nz, 0.02)
    _apply_outer_boundaries_d3q19(
        symmetry, inlet,
        boundary_scheme="equilibrium_copy",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary="all_symmetry_free_slip",
    )
    _apply_outer_boundaries_d3q19(
        ground, inlet,
        boundary_scheme="equilibrium_copy",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary=(
            "ground_no_slip_lateral_top_free_slip"
        ),
    )
    direction = next(
        index for index, vector in enumerate(_C3)
        if tuple(vector) == (1, 0, 1)
    )
    opposite = next(
        index for index, vector in enumerate(_C3)
        if tuple(vector) == (-1, 0, -1)
    )
    specular = next(
        index for index, vector in enumerate(_C3)
        if tuple(vector) == (1, 0, -1)
    )
    point = (0, ny // 2, nx // 2)
    assert symmetry[(direction, *point)] == pytest.approx(
        original[(specular, *point)]
    )
    assert ground[(direction, *point)] == pytest.approx(
        original[(opposite, *point)]
    )


@pytest.mark.parametrize(
    "boundary_scheme",
    [
        "nonequilibrium_extrapolation",
        "characteristic_nonreflecting",
        "characteristic_multidimensional",
        "characteristic_multidimensional_extrapolated",
        "characteristic_multidimensional_vortical_limited",
        "characteristic_nscbc_convective",
    ],
)
def test_pressure_open_boundary_preserves_uniform_equilibrium_state(
    boundary_scheme,
):
    nz, ny, nx = 5, 6, 7
    rho = np.full((nz, ny, nx), 1.003)
    ux = np.full((nz, ny, nx), 0.02)
    zero = np.zeros_like(rho)
    distribution = _equilibrium_d3q19(rho, ux, zero, zero)

    _apply_outer_boundaries_d3q19(
        distribution,
        np.full(nz, 0.02),
        boundary_scheme=boundary_scheme,
        outlet_lattice_density=1.003,
        nonstreamwise_wall_boundary="all_pressure_open",
    )

    result_rho = np.sum(distribution, axis=0)
    result_velocity = np.stack([
        np.tensordot(_C3[:, component], distribution, axes=(0, 0))
        / result_rho
        for component in range(3)
    ])
    assert result_rho == pytest.approx(rho, abs=1.0e-14)
    assert result_velocity[0] == pytest.approx(ux, abs=1.0e-14)
    assert result_velocity[1:] == pytest.approx(0.0, abs=1.0e-14)


def test_pressure_open_boundary_sets_lateral_density_from_ambient_reference():
    nz, ny, nx = 5, 6, 7
    zz, yy, xx = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij",
    )
    rho = 1.0 + 1.0e-3 * yy + 2.0e-4 * zz
    ux = np.full_like(rho, 0.02)
    zero = np.zeros_like(rho)
    distribution = _equilibrium_d3q19(rho, ux, zero, zero)

    _apply_outer_boundaries_d3q19(
        distribution,
        np.full(nz, 0.02),
        boundary_scheme="nonequilibrium_extrapolation",
        outlet_lattice_density=1.001,
        nonstreamwise_wall_boundary="all_pressure_open",
    )

    result_rho = np.sum(distribution, axis=0)
    # x faces have streamwise precedence; inspect the lateral-face interiors.
    assert result_rho[:, 0, 1:-1] == pytest.approx(1.001, abs=1.0e-14)
    assert result_rho[:, -1, 1:-1] == pytest.approx(1.001, abs=1.0e-14)
    assert result_rho[0, 1:-1, 1:-1] == pytest.approx(1.001, abs=1.0e-14)
    assert result_rho[-1, 1:-1, 1:-1] == pytest.approx(1.001, abs=1.0e-14)


def test_characteristic_outlet_sets_incoming_isothermal_invariant():
    nz, ny, nx = 5, 6, 7
    zz, yy, xx = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij",
    )
    rho = 1.0 + 2.0e-4 * xx
    sound_speed = np.sqrt(1.0 / 3.0)
    ux = sound_speed * np.log(rho)
    zero = np.zeros_like(rho)
    distribution = _equilibrium_d3q19(rho, ux, zero, zero)

    _apply_outer_boundaries_d3q19(
        distribution,
        np.zeros(nz),
        boundary_scheme="characteristic_nonreflecting",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary="all_symmetry_free_slip",
    )

    outlet = distribution[:, :, :, -1]
    outlet_rho = np.sum(outlet, axis=0)
    outlet_ux = np.tensordot(_C3[:, 0], outlet, axes=(0, 0)) / outlet_rho
    incoming_invariant = outlet_ux - sound_speed * np.log(outlet_rho)
    assert incoming_invariant[1:-1, 1:-1] == pytest.approx(
        0.0, abs=1.0e-14,
    )


def test_multidimensional_characteristic_preserves_oblique_outgoing_impedance():
    sound_speed = np.sqrt(1.0 / 3.0)
    angle = np.radians(45.0)
    log_density = np.asarray([[2.0e-4, -3.0e-4]])
    density = np.exp(log_density)
    outward = sound_speed * np.cos(angle) * log_density
    tangential = sound_speed * np.sin(angle) * log_density

    target_density, target_outward, cosine = (
        _multidimensional_characteristic_state(
            density,
            outward,
            tangential,
            np.zeros_like(outward),
            reference_density=1.0,
            reference_outward_velocity=0.0,
        )
    )

    assert cosine == pytest.approx(np.cos(angle), abs=1.0e-13)
    assert target_density == pytest.approx(density, abs=1.0e-13)
    assert target_outward == pytest.approx(outward, abs=1.0e-13)


def test_multidimensional_characteristic_preserves_uniform_carrier():
    density = np.ones((2, 3))
    carrier = np.full_like(density, 0.04)
    target_density, target_outward, cosine = (
        _multidimensional_characteristic_state(
            density,
            carrier,
            np.zeros_like(density),
            np.zeros_like(density),
            reference_density=1.0,
            reference_outward_velocity=carrier,
        )
    )
    assert target_density == pytest.approx(density, abs=1.0e-14)
    assert target_outward == pytest.approx(carrier, abs=1.0e-14)
    assert cosine == pytest.approx(1.0, abs=1.0e-14)


def test_zero_acoustic_fraction_preserves_pressure_free_rotation():
    density = np.ones((2, 3))
    normal = np.full_like(density, 0.012)
    tangent = np.full_like(density, -0.018)
    target_density, target_normal, _ = _multidimensional_characteristic_state(
        density,
        normal,
        tangent,
        np.zeros_like(density),
        reference_density=1.0,
        reference_outward_velocity=0.0,
        acoustic_fraction=np.zeros_like(density),
    )
    assert target_density == pytest.approx(density, abs=1.0e-14)
    assert target_normal == pytest.approx(normal, abs=1.0e-14)


def test_unit_acoustic_fraction_retains_pure_acoustic_correction():
    sound_speed = np.sqrt(1.0 / 3.0)
    log_density = np.asarray([[2.0e-4, -3.0e-4]])
    density = np.exp(log_density)
    normal = sound_speed * log_density
    base = _multidimensional_characteristic_state(
        density,
        normal,
        np.zeros_like(normal),
        np.zeros_like(normal),
        reference_density=1.0,
        reference_outward_velocity=0.0,
    )
    weighted = _multidimensional_characteristic_state(
        density,
        normal,
        np.zeros_like(normal),
        np.zeros_like(normal),
        reference_density=1.0,
        reference_outward_velocity=0.0,
        acoustic_fraction=np.ones_like(density),
    )
    assert weighted[0] == pytest.approx(base[0], abs=1.0e-14)
    assert weighted[1] == pytest.approx(base[1], abs=1.0e-14)


def test_boundary_mode_fraction_separates_curl_and_acoustic_gradient():
    shape = (3, 5, 1)
    density = np.ones(shape)
    y = np.arange(shape[1], dtype=float)[None, :, None]
    ux_vortex = np.broadcast_to(0.001 * y, shape).copy()
    zero = np.zeros(shape)
    vortex = _equilibrium_d3q19(density, ux_vortex, zero, zero)
    vortical_fraction = _boundary_acoustic_fraction(
        vortex, vortex,
        normal_spatial_axis=2,
        normal_coordinate_sign=1.0,
    )
    assert vortical_fraction == pytest.approx(0.0, abs=5.0e-13)

    sound_speed = np.sqrt(1.0 / 3.0)
    log_density = 2.0e-4
    acoustic_first = _equilibrium_d3q19(
        np.full(shape, np.exp(log_density)),
        np.full(shape, sound_speed * log_density), zero, zero,
    )
    acoustic_second = _equilibrium_d3q19(density, zero, zero, zero)
    acoustic_fraction = _boundary_acoustic_fraction(
        acoustic_first, acoustic_second,
        normal_spatial_axis=2,
        normal_coordinate_sign=1.0,
    )
    assert acoustic_fraction == pytest.approx(1.0, abs=1.0e-12)


def test_convective_boundary_uses_local_outward_courant_and_stops_backflow():
    shape = (2, 3, 1)
    zero = np.zeros(shape)
    boundary = _equilibrium_d3q19(
        np.ones(shape), np.full(shape, 0.08), zero, zero,
    )
    interior = _equilibrium_d3q19(
        np.full(shape, 0.99), np.full(shape, 0.07), zero, zero,
    )
    advanced, courant = _convective_boundary_distribution(
        boundary, interior,
        normal_velocity_component=0,
        outward_sign=1.0,
    )
    assert courant == pytest.approx(0.08, abs=1.0e-14)
    assert advanced == pytest.approx(
        boundary - 0.08 * (boundary - interior), abs=1.0e-14,
    )

    backflow = _equilibrium_d3q19(
        np.ones(shape), np.full(shape, -0.02), zero, zero,
    )
    stopped, backflow_courant = _convective_boundary_distribution(
        backflow, interior,
        normal_velocity_component=0,
        outward_sign=1.0,
    )
    assert backflow_courant == pytest.approx(0.0, abs=1.0e-14)
    assert stopped == pytest.approx(backflow, abs=1.0e-14)


def test_nscbc_xmax_preserves_uniform_carrier_exactly():
    shape = (3, 5, 1)
    density = np.ones(shape)
    ux = np.full(shape, 0.08)
    zero = np.zeros(shape)
    equilibrium = _equilibrium_d3q19(density, ux, zero, zero)
    result = _nscbc_xmax_boundary_distribution(
        equilibrium,
        equilibrium,
        reference_density=1.0,
        reference_velocity_x=ux,
    )
    assert result == pytest.approx(equilibrium, abs=1.0e-14)


def test_nscbc_xmax_advects_outgoing_invariant_without_incoming_normal_wave():
    shape = (3, 5, 1)
    sound_speed = np.sqrt(1.0 / 3.0)
    zero = np.zeros(shape)
    boundary_log_density = 2.0e-4
    interior_log_density = 3.0e-4
    boundary_density = np.full(shape, np.exp(boundary_log_density))
    interior_density = np.full(shape, np.exp(interior_log_density))
    boundary_ux = np.full(shape, 0.08 + sound_speed * boundary_log_density)
    interior_ux = np.full(shape, 0.08 + sound_speed * interior_log_density)
    boundary = _equilibrium_d3q19(
        boundary_density, boundary_ux, zero, zero,
    )
    interior = _equilibrium_d3q19(
        interior_density, interior_ux, zero, zero,
    )
    result = _nscbc_xmax_boundary_distribution(
        boundary,
        interior,
        reference_density=1.0,
        reference_velocity_x=np.full(shape, 0.08),
    )
    result_rho = np.sum(result, axis=0)
    result_ux = np.tensordot(_C3[:, 0], result, axes=(0, 0)) / result_rho
    incoming_before = boundary_ux - sound_speed * boundary_log_density
    incoming_after = result_ux - sound_speed * np.log(result_rho)
    assert incoming_after == pytest.approx(incoming_before, abs=1.0e-13)


def test_limited_boundary_extrapolation_preserves_linear_slope_and_extrema():
    linear = _limited_boundary_extrapolation(
        np.asarray([3.0, -3.0]),
        np.asarray([2.0, -2.0]),
        np.asarray([1.0, -1.0]),
    )
    assert linear == pytest.approx([4.0, -4.0])
    extremum = _limited_boundary_extrapolation(
        np.asarray([3.0, -3.0]),
        np.asarray([2.0, -2.0]),
        np.asarray([2.5, -2.5]),
    )
    assert extremum == pytest.approx([3.0, -3.0])
    limited = _limited_boundary_extrapolation(
        np.asarray([3.0]), np.asarray([2.0]), np.asarray([1.5]),
    )
    assert limited == pytest.approx([3.5])


def test_selected_lateral_characteristic_face_can_be_isolated_for_audit():
    nz, ny, nx = 5, 7, 6
    zz, yy, xx = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij",
    )
    rho = 1.0 + 2.0e-4 * yy
    sound_speed = np.sqrt(1.0 / 3.0)
    uy = sound_speed * np.log(rho)
    zero = np.zeros_like(rho)
    distribution = _equilibrium_d3q19(rho, zero, uy, zero)
    original_x_min_interior = distribution[:, 1:-1, 1:-1, 0].copy()

    _apply_outer_boundaries_d3q19(
        distribution,
        np.zeros(nz),
        boundary_scheme="characteristic_nonreflecting",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary="all_pressure_open",
        pressure_open_nonstream_faces=("y_max",),
        apply_streamwise_boundaries=False,
    )

    lateral = distribution[:, :, -1, :]
    lateral_rho = np.sum(lateral, axis=0)
    lateral_uy = np.tensordot(_C3[:, 1], lateral, axes=(0, 0)) / lateral_rho
    incoming_invariant = lateral_uy - sound_speed * np.log(lateral_rho)
    assert incoming_invariant[1:-1, :] == pytest.approx(0.0, abs=1.0e-14)
    assert distribution[:, 1:-1, 1:-1, 0] == pytest.approx(
        original_x_min_interior,
    )


def test_pressure_open_mode_runs_and_reports_unqualified_candidate_status():
    result = solve_resolved_obstacle_flow_3d(_config(
        obstacles=(),
        flow_steps=2,
        minimum_flow_steps=2,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        flow_boundary_scheme="characteristic_nonreflecting",
        nonstreamwise_wall_boundary="all_pressure_open",
    ))

    assert result.diagnostics["pressure_projection_open_boundary_faces"] == [
        "x_max", "y_min", "y_max", "z_min", "z_max",
    ]
    assert result.diagnostics["pressure_open_boundary_qualification"] == (
        "normal_and_30deg_linear_acoustic_gates_passed_"
        "angle_sweep_failed_corner_and_domain_pending"
    )
    assert result.maximum_mass_residual_kg < 1.0e-12


def test_resolved_3d_wale_is_opt_in_and_reported():
    result = solve_resolved_obstacle_flow_3d(_config(
        flow_steps=2,
        minimum_flow_steps=2,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        subgrid_model="wale",
        smagorinsky_constant=0.0,
        wale_constant=0.5,
    ))

    assert result.diagnostics["subgrid_model"] == "wale"
    assert result.diagnostics["turbulence_model"].startswith("wall-adapting")
    assert result.diagnostics["wall_adapting_subgrid_model"] is True
    assert result.diagnostics["wall_damping"] is False
    assert np.nanmax(result.eddy_viscosity_m2_s) >= 0.0


def test_resolved_3d_nonequilibrium_boundary_is_opt_in_and_reported():
    result = solve_resolved_obstacle_flow_3d(_config(
        flow_steps=2,
        minimum_flow_steps=2,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        flow_boundary_scheme="nonequilibrium_extrapolation",
        outlet_lattice_density=1.0,
    ))

    assert result.diagnostics["flow_boundary_scheme"] == (
        "nonequilibrium_extrapolation"
    )
    assert "non-equilibrium extrapolation" in result.diagnostics["inlet_boundary"]
    assert "fixed lattice density" in result.diagnostics["outlet_boundary"]
    assert result.diagnostics["outlet_lattice_density"] == pytest.approx(1.0)


def test_resolved_3d_time_averages_flow_and_uses_statistical_gate():
    result = solve_resolved_obstacle_flow_3d(_config(
        flow_steps=6,
        minimum_flow_steps=6,
        flow_average_last_steps=4,
        flow_statistical_convergence_tolerance=10.0,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
    ))

    diagnostics = result.diagnostics
    assert diagnostics["flow_field_statistic"] == "time_mean"
    assert diagnostics["flow_average_last_steps"] == 4
    assert diagnostics["flow_average_sample_count"] == 4
    assert diagnostics["flow_statistical_residual"] is not None
    assert diagnostics["flow_statistical_residual"] >= 0.0
    assert diagnostics["flow_convergence_criterion"].startswith("first-half")
    assert diagnostics["flow_converged"] is True
    assert diagnostics["final_flow_residual"] == pytest.approx(
        diagnostics["flow_statistical_residual"]
    )


def test_resolved_3d_rejects_receptor_inside_obstacle():
    with pytest.raises(ValueError, match="interpolation stencil intersects"):
        solve_resolved_obstacle_flow_3d(
            _config(
                flow_steps=1,
                minimum_flow_steps=1,
                duration_s=0.001,
                scalar_time_step_s=0.001,
                source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
            ),
            receptors=(FastObstacleReceptor3D("solid", (1.75, 1.2, 0.5)),),
        )


def test_resolved_3d_active_thermal_scalar_conserves_energy_and_forces_momentum():
    dt = 0.03 * 0.1 / 2.0
    result = solve_resolved_obstacle_flow_3d(
        _config(
            duration_s=10.0 * dt,
            scalar_time_step_s=dt,
            source_schedule=SourceRateSchedule(
                (0.0, 5.0 * dt, 10.0 * dt), (0.004, 0.0, 0.0),
            ),
            thermal_buoyancy_coupling=True,
            source_temperature_k=80.0,
            thermal_diffusivity_m2_s=0.03,
            maximum_active_scalar_density_departure=1.0e-6,
            store_scalar_fields=True,
        ),
        receptors=(FastObstacleReceptor3D("source", (0.5, 1.2, 0.4)),),
    )

    assert result.temperature_k is not None
    assert result.hydrogen_mole_fraction is not None
    assert result.maximum_hydrogen_mole_fraction is not None
    assert result.mixture_density_kg_m3 is not None
    assert result.buoyancy_acceleration_m_s2 is not None
    assert result.energy_residual_j is not None
    assert result.temperature_fields_k is not None
    assert np.nanmin(result.temperature_k) < result.ambient_temperature_k
    assert np.nanmax(result.maximum_concentration_vol_percent) > 0.0
    assert result.maximum_mass_residual_kg < 1.0e-12
    assert result.maximum_energy_residual_j is not None
    assert result.maximum_energy_residual_j < 1.0e-8
    assert result.diagnostics["energy_equation"] is True
    assert result.diagnostics["density_buoyancy_feedback"] is True
    assert result.diagnostics["variable_density_continuity_equation"] is True
    assert result.diagnostics["pressure_correction_to_momentum"] is True
    assert result.diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
    assert result.diagnostics[
        "thermodynamic_corrector_maximum_iterations_used"
    ] <= 3
    assert result.diagnostics[
        "maximum_variable_density_continuity_residual_kg_m3_s"
    ] < 1.0e-8
    assert result.pressure_correction_pa is not None
    assert result.continuity_residual_kg_m3_s is not None
    assert np.nanmax(np.abs(result.pressure_correction_pa)) > 0.0
    assert np.nanmax(np.abs(result.continuity_residual_kg_m3_s)) < 1.0e-8
    assert result.diagnostics["maximum_absolute_buoyancy_acceleration_m_s2"] > 0.0
    assert result.diagnostics["coupled_transient_steps"] == 10
    assert result.diagnostics["active_scalar_applicability_passed"] is False
    assert result.diagnostics["active_scalar_out_of_range_cell_steps"] > 0
    assert result.diagnostics["maximum_relative_mass_residual"] < 1.0e-12
    assert result.diagnostics["maximum_relative_energy_residual"] < 1.0e-12
    assert result.diagnostics["maximum_lattice_mach_number"] < 0.15
    assert result.diagnostics["quantitative_prediction_allowed"] is False
    assert result.receptor_traces_vol_percent["source"].shape == result.time_s.shape
    assert result.receptor_temperature_traces_k["source"].shape == result.time_s.shape


def test_resolved_3d_fv_active_path_applies_declared_source_momentum():
    result = solve_resolved_obstacle_flow_3d(_config(
        flow_steps=1,
        minimum_flow_steps=1,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule(
            (0.0, 0.001), (0.01, 0.0),
        ),
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        maximum_active_scalar_density_departure=1.0,
        active_momentum_solver="finite_volume_projection",
        source_momentum_velocity_m_s=(1.5, 0.0, 0.2),
        active_inlet_profile_y_m=(0.0, 1.2, 2.4),
        active_inlet_profile_u_m_s=(1.0, 2.0, 3.0),
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.85,
        scalar_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        thermal_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
    ))

    assert result.diagnostics["active_momentum_solver"] == (
        "finite_volume_projection"
    )
    assert result.diagnostics["active_inlet_profile_model"] == (
        "declared_y_profile_times_z_profile"
    )
    assert result.diagnostics["active_inlet_plane_minimum_m_s"] == pytest.approx(
        1.0 + 1.0 / 24.0
    )
    assert result.diagnostics["active_inlet_plane_maximum_m_s"] == pytest.approx(
        3.0 - 1.0 / 24.0
    )
    assert result.diagnostics["active_inlet_plane_mean_m_s"] == pytest.approx(2.0)
    assert result.u_face_m_s[:, :, 0] == pytest.approx(
        np.broadcast_to(
            np.linspace(1.0 + 1.0 / 24.0, 3.0 - 1.0 / 24.0, 24)[None, :],
            (20, 24),
        )
    )
    assert result.diagnostics["source_momentum_velocity_m_s"] == [1.5, 0.0, 0.2]
    assert result.diagnostics["maximum_momentum_balance_residual_kg_m_s"] < 1.0e-12
    assert result.diagnostics["applied_scalar_diffusivity_m2_s"] is None
    assert result.diagnostics[
        "applied_scalar_diffusivity_maximum_m2_s"
    ] > result.diagnostics["applied_scalar_diffusivity_minimum_m2_s"]
    assert result.diagnostics[
        "maximum_eddy_viscosity_over_transient_m2_s"
    ] > 0.0


def test_resolved_3d_fv_active_path_applies_time_dependent_vector_inlet():
    history = ResolvedInletVelocityHistory3D(
        time_s=(0.0, 0.001),
        u_m_s=(0.8, 1.2),
        v_m_s=(0.2, 0.2),
        w_m_s=(0.1, 0.1),
        history_id="synthetic-10hz-vector-v1",
    )
    result = solve_resolved_obstacle_flow_3d(_config(
        flow_steps=1,
        minimum_flow_steps=1,
        duration_s=0.001,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule(
            (0.0, 0.001), (0.0, 0.0),
        ),
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        maximum_active_scalar_density_departure=1.0,
        active_momentum_solver="finite_volume_projection",
        active_inlet_profile_y_m=(0.0, 1.2, 2.4),
        active_inlet_profile_u_m_s=(1.0, 2.0, 3.0),
        active_inlet_velocity_history=history,
    ))

    expected_y = np.linspace(
        0.5 + 0.5 / 24.0, 1.5 - 0.5 / 24.0, 24,
    )
    assert result.u_face_m_s[:, :, 0] == pytest.approx(
        np.broadcast_to(expected_y[None, :], (20, 24))
    )
    diagnostics = result.diagnostics
    assert diagnostics["active_inlet_velocity_history_id"] == (
        "synthetic-10hz-vector-v1"
    )
    assert diagnostics["active_inlet_velocity_history_sample_count"] == 2
    assert diagnostics["active_inlet_temporal_mapping"] == (
        "linear_interpolation_single_vector_times_declared_spatial_shape"
    )
    assert diagnostics["minimum_applied_transient_inlet_u_m_s"] == pytest.approx(
        0.8 * (0.5 + 0.5 / 24.0)
    )
    assert diagnostics["maximum_applied_transient_inlet_u_m_s"] == pytest.approx(
        1.5 - 0.5 / 24.0
    )
    assert any(
        "does not reconstruct spatial turbulence" in warning
        for warning in diagnostics["warnings"]
    )


def test_resolved_3d_active_thermal_scalar_requires_lbm_time_alignment():
    with pytest.raises(ValueError, match="scalar_time_step_s to equal"):
        solve_resolved_obstacle_flow_3d(
            _config(
                duration_s=0.015,
                scalar_time_step_s=0.001,
                source_schedule=SourceRateSchedule(
                    (0.0, 0.0075, 0.015), (0.0, 0.0, 0.0),
                ),
                thermal_buoyancy_coupling=True,
                source_temperature_k=80.0,
            )
        )
