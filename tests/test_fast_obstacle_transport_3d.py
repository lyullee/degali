import numpy as np
import pytest

from degali.addons.fast_obstacle_transport_3d import (
    FAST_OBSTACLE_TRANSPORT_3D_SCHEMA,
    FastObstacleReceptor3D,
    FastObstacleTransport3DConfig,
    _circular_arc_blocked_faces,
    _grid_and_mask,
    _source_weights,
    _source_weights_xyz,
    _source_weights_cosine_tapered_box,
    _source_weights_variance_matched_box,
    _transport_step,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import (
    AxisAlignedCuboid, CircularArcWall, CrosswindCylinder,
)


def _config(**changes):
    values = {
        "domain_m": (4.0, 2.4, 2.0),
        "cells": (40, 24, 20),
        "duration_s": 0.2,
        "time_step_s": 0.001,
        "wind_speed_m_s": 2.0,
        "diffusivity_m2_s": 0.03,
        "source_schedule": SourceRateSchedule(
            (0.0, 0.05, 0.2), (0.02, 0.0, 0.0),
        ),
        "source_position_m": (0.5, 1.2, 0.4),
        "source_sigma_m": 0.15,
        "obstacles": (
            AxisAlignedCuboid(1.5, 2.0, 0.8, 1.6, 0.0, 1.0, "finite block"),
        ),
    }
    values.update(changes)
    return FastObstacleTransport3DConfig(**values)


@pytest.fixture(scope="module")
def finite_width_result():
    return solve_fast_obstacle_transport_3d(
        _config(),
        receptors=(
            FastObstacleReceptor3D("windward", (1.2, 1.2, 0.5)),
            FastObstacleReceptor3D("side", (1.8, 0.6, 0.5)),
            FastObstacleReceptor3D("roof", (1.8, 1.2, 1.2)),
            FastObstacleReceptor3D("lee", (2.4, 1.2, 0.5)),
        ),
    )


def test_open_domain_reproduces_uniform_flow():
    result = solve_fast_obstacle_transport_3d(_config(
        obstacles=(), duration_s=0.001, time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
    ))

    assert result.velocity_u_m_s == pytest.approx(2.0)
    assert result.velocity_v_m_s == pytest.approx(0.0)
    assert result.velocity_w_m_s == pytest.approx(0.0)
    assert result.diagnostics["maximum_absolute_divergence_s_1"] < 1.0e-12


def test_fast_solver_accepts_conserved_state_and_momentum_interface_callbacks():
    calls = []

    def state_callback(concentration, start, dt):
        calls.append(("state", start, dt))
        return concentration

    def momentum_callback(u, v, w, u_face, v_face, w_face, start, dt):
        calls.append(("momentum", start, dt))
        return u, v, w, u_face, v_face, w_face

    result = solve_fast_obstacle_transport_3d(
        _config(
            duration_s=0.002,
            time_step_s=0.001,
            source_schedule=SourceRateSchedule((0.0, 0.002), (0.0, 0.0)),
        ),
        interface_state_callback=state_callback,
        interface_momentum_callback=momentum_callback,
    )

    assert len(calls) == 4
    assert result.diagnostics["interface_momentum_callback_applied"] is True


def test_fast_3d_tier_rejects_quiescent_ffi_channel_condition():
    with pytest.raises(ValueError, match="wind_speed_m_s must be positive"):
        _config(wind_speed_m_s=0.0)


def test_fast_3d_serializes_declared_open_boundary_faces():
    result = solve_fast_obstacle_transport_3d(
        _config(
            duration_s=0.001,
            time_step_s=0.001,
            open_boundary_faces=("x_max",),
            source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        )
    )

    assert result.diagnostics["open_boundary_faces"] == ["x_max"]


def test_receptor_stencil_accepts_exact_domain_wall_coordinate():
    result = solve_fast_obstacle_transport_3d(
        _config(
            duration_s=0.001,
            time_step_s=0.001,
            source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        ),
        receptors=(FastObstacleReceptor3D("wall", (0.0, 0.0, 2.0)),),
    )

    assert result.receptor_traces_kg_m3["wall"].shape == (2,)


def test_transport_step_accepts_anisotropic_diffusivity_and_conserves_mass():
    concentration = np.zeros((3, 2, 2), dtype=float)
    concentration[0, 0, 0] = 1.0
    solid = np.zeros_like(concentration, dtype=bool)
    u_face = np.zeros((3, 2, 3), dtype=float)
    v_face = np.zeros((3, 3, 2), dtype=float)
    w_face = np.zeros((4, 2, 2), dtype=float)
    horizontal = np.zeros_like(concentration)
    vertical = np.full_like(concentration, 0.05)

    updated, boundary_rate = _transport_step(
        concentration,
        u_face=u_face,
        v_face=v_face,
        w_face=w_face,
        solid=solid,
        diffusivity=(horizontal, horizontal, vertical),
        dt=0.1,
        spacing=np.ones(3),
        open_boundary_faces=(),
    )

    assert updated[1, 0, 0] > 0.0
    assert updated[0, 0, 1] == pytest.approx(concentration[0, 0, 1])
    assert updated[0, 1, 0] == pytest.approx(concentration[0, 1, 0])
    assert np.sum(updated) == pytest.approx(np.sum(concentration))
    assert boundary_rate == pytest.approx(0.0)


def test_source_weights_xyz_are_normalised_and_resolved_by_cell_integral():
    config = _config(
        domain_m=(1.0, 1.0, 1.0), cells=(10, 10, 10),
        source_position_m=(0.5, 0.5, 0.05), source_sigma_m=0.05,
        source_sigma_xyz_m=(0.2, 0.2, 0.1), obstacles=(),
    )
    x, y, z, solid, _ = _grid_and_mask(config)
    weights = _source_weights_xyz(
        x, y, z, solid, config.source_position_m, config.source_sigma_xyz_m,
    )

    assert np.sum(weights) == pytest.approx(1.0)
    assert np.all(weights >= 0.0)
    assert np.sum(weights[0, :, :]) > np.sum(weights[-1, :, :])


def test_variance_matched_box_source_is_compact_and_conservative():
    x = (np.arange(20, dtype=float) + 0.5) / 20.0
    y = (np.arange(20, dtype=float) + 0.5) / 20.0
    z = (np.arange(20, dtype=float) + 0.5) / 20.0
    solid = np.zeros((20, 20, 20), dtype=bool)
    weights = _source_weights_variance_matched_box(
        x, y, z, solid, (0.5, 0.5, 0.5), (0.08, 0.08, 0.05),
    )
    assert np.sum(weights) == pytest.approx(1.0)
    assert np.all(weights >= 0.0)
    assert np.count_nonzero(weights) < weights.size


def test_cosine_tapered_box_source_is_compact_smooth_and_conservative():
    x = (np.arange(20, dtype=float) + 0.5) / 20.0
    y = (np.arange(20, dtype=float) + 0.5) / 20.0
    z = (np.arange(20, dtype=float) + 0.5) / 20.0
    solid = np.zeros((20, 20, 20), dtype=bool)
    weights = _source_weights_cosine_tapered_box(
        x, y, z, solid, (0.5, 0.5, 0.5), (0.08, 0.08, 0.05),
    )
    assert np.sum(weights) == pytest.approx(1.0)
    assert np.all(weights >= 0.0)
    support = np.flatnonzero(np.sum(weights, axis=(0, 1)) > 0.0)
    assert support.size < x.size
    assert weights[10, 10, 10] > weights[10, 10, support[0]]


def test_open_domain_preserves_a_declared_vertical_inlet_profile():
    result = solve_fast_obstacle_transport_3d(_config(
        obstacles=(), duration_s=0.001, time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        inlet_profile_z_m=(0.0, 1.0, 2.0),
        inlet_profile_u_m_s=(0.0, 1.0, 2.0),
    ))

    expected = np.interp(result.z_m, (0.0, 1.0, 2.0), (0.0, 1.0, 2.0))
    assert result.u_face_m_s[:, 0, 0] == pytest.approx(expected)
    assert result.velocity_u_m_s[:, 0, 0] == pytest.approx(expected)
    assert "measured/interpolated" in result.diagnostics["inlet_boundary"]


def test_cuboid_voxelization_uses_half_open_cell_membership():
    config = _config(
        obstacles=(AxisAlignedCuboid(
            1.55, 2.05, 0.85, 1.65, 0.05, 1.05, "half-open block",
        ),),
    )

    x, y, z, solid, _ = _grid_and_mask(config)

    expected = (
        (x[None, None, :] >= 1.55) & (x[None, None, :] < 2.05)
        & (y[None, :, None] >= 0.85) & (y[None, :, None] < 1.65)
        & (z[:, None, None] >= 0.05) & (z[:, None, None] < 1.05)
    )
    assert np.array_equal(solid, expected)
    assert np.count_nonzero(np.any(solid, axis=(0, 1))) == 5


def test_cell_integrated_gaussian_source_is_nested_grid_consistent():
    domain = (0.8, 1.2, 1.0)

    def coordinates(count: int, length: float) -> np.ndarray:
        spacing = length / count
        return (np.arange(count) + 0.5) * spacing

    coarse_shape = (5, 6, 4)
    fine_shape = tuple(2 * count for count in coarse_shape)
    coarse_z, coarse_y, coarse_x = (
        coordinates(count, length)
        for count, length in zip(coarse_shape, domain[::-1])
    )
    fine_z, fine_y, fine_x = (
        coordinates(count, length)
        for count, length in zip(fine_shape, domain[::-1])
    )
    position = (0.31, 0.47, 0.29)
    sigma = 0.14
    coarse = _source_weights(
        coarse_x, coarse_y, coarse_z, np.zeros(coarse_shape, dtype=bool),
        position, sigma,
    )
    fine = _source_weights(
        fine_x, fine_y, fine_z, np.zeros(fine_shape, dtype=bool),
        position, sigma,
    )
    aggregated = fine.reshape(5, 2, 6, 2, 4, 2).sum(axis=(1, 3, 5))

    assert np.sum(coarse) == pytest.approx(1.0)
    assert np.all(coarse >= 0.0)
    assert aggregated == pytest.approx(coarse, rel=2.0e-14, abs=2.0e-15)


def test_finite_width_flow_has_local_roof_and_lateral_bypass(finite_width_result):
    result = finite_width_result
    solid = result.solid_mask

    assert result.schema == FAST_OBSTACLE_TRANSPORT_3D_SCHEMA
    assert np.max(np.abs(result.velocity_v_m_s[~solid])) > 0.05
    assert np.max(result.velocity_w_m_s[~solid]) > 0.05
    assert np.min(result.velocity_v_m_s[~solid]) < -0.05
    assert np.max(result.velocity_v_m_s[~solid]) > 0.05
    assert result.diagnostics["finite_lateral_bypass"] is True
    assert result.diagnostics["recirculating_wake_resolved"] is False
    assert result.diagnostics["maximum_absolute_divergence_s_1"] < 1.0e-8


def test_finite_width_scalar_flux_is_impermeable_and_mass_conservative(
    finite_width_result,
):
    result = finite_width_result
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


def test_fast_3d_result_is_explicitly_nonquantitative(finite_width_result):
    diagnostics = finite_width_result.diagnostics

    assert diagnostics["quantitative_prediction_allowed"] is False
    assert diagnostics["nonlocal_scalar_routing"] is False
    assert diagnostics["local_conservative_scalar_transport"] is True
    assert diagnostics["applicability"] == "conditional"
    assert diagnostics["source_spatial_discretization"].startswith(
        "cell-integrated"
    )
    assert diagnostics["flow_solve_runtime_s"] > 0.0
    assert diagnostics["total_runtime_s"] >= diagnostics["flow_solve_runtime_s"]
    assert "separation" in diagnostics["warnings"][0]


def test_fast_3d_rejects_a_receptor_inside_the_obstacle():
    with pytest.raises(ValueError, match="interpolation stencil intersects"):
        solve_fast_obstacle_transport_3d(
            _config(duration_s=0.001, time_step_s=0.001,
                    source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0))),
            receptors=(FastObstacleReceptor3D("solid", (1.75, 1.2, 0.5)),),
        )


def test_fast_3d_accepts_a_curved_wall_as_a_local_solid_mask():
    config = _config(
        obstacles=(CircularArcWall(
            center_m=(1.2, 1.2), radius_m=0.8, thickness_m=0.2,
            z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
            arc_end_deg=90.0, label="curved fence",
        ),),
        duration_s=0.001,
        time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
    )
    result = solve_fast_obstacle_transport_3d(config)
    assert np.any(result.solid_mask)
    assert np.all(np.isnan(result.concentration_kg_m3[result.solid_mask]))
    assert np.all(result.velocity_u_m_s[result.solid_mask] == 0.0)
    assert np.all(result.velocity_v_m_s[result.solid_mask] == 0.0)
    assert np.all(result.velocity_w_m_s[result.solid_mask] == 0.0)
    assert result.diagnostics["obstacle_count"] == 1
    assert result.diagnostics["embedded_surface_blocked_face_count"] > 0
    blocked = _circular_arc_blocked_faces(
        config, result.x_m, result.y_m, result.z_m,
    )
    assert np.all(result.u_face_m_s[blocked[0]] == 0.0)
    assert np.all(result.v_face_m_s[blocked[1]] == 0.0)
    assert np.all(result.w_face_m_s[blocked[2]] == 0.0)


def test_fast_3d_accepts_finite_crosswind_cylinder_obstacles():
    rod = CrosswindCylinder(
        center_x_m=1.55, center_z_m=0.45, radius_m=0.05,
        y_min_m=0.2, y_max_m=2.2, label="horizontal-array-rod",
    )
    config = _config(
        obstacles=(rod,), duration_s=0.001, time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
    )
    result = solve_fast_obstacle_transport_3d(config)

    assert np.any(result.solid_mask)
    assert result.diagnostics["obstacle_count"] == 1
    assert np.all(np.isnan(result.concentration_kg_m3[result.solid_mask]))


def test_fast_3d_rejects_obstacle_that_disappears_between_cell_centres():
    rod = CrosswindCylinder(
        center_x_m=1.50, center_z_m=0.40, radius_m=0.05,
        y_min_m=0.2, y_max_m=2.2, label="subgrid-rod",
    )
    config = _config(
        obstacles=(rod,), duration_s=0.001, time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
    )

    with pytest.raises(ValueError, match="subgrid-rod.*no voxelized solid cells"):
        _grid_and_mask(config)

def test_transport_step_exposes_face_flux_for_coupling():
    seen = []
    state = np.ones((2, 2, 2), dtype=float)
    zeros_u = np.zeros((2, 2, 3), dtype=float)
    zeros_v = np.zeros((2, 3, 2), dtype=float)
    zeros_w = np.zeros((3, 2, 2), dtype=float)
    _transport_step(state, u_face=zeros_u, v_face=zeros_v, w_face=zeros_w,
                    solid=np.zeros_like(state, dtype=bool), diffusivity=0.01,
                    dt=0.01, spacing=np.array((1., 1., 1.)),
                    face_flux_callback=lambda fx, fy, fz: seen.append((fx, fy, fz)))
    assert len(seen) == 1
    assert seen[0][0].shape == zeros_u.shape

def test_solver_forwards_face_flux_callback():
    seen = []
    solve_fast_obstacle_transport_3d(
        _config(duration_s=0.001, time_step_s=0.001,
                source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0))),
        face_flux_callback=lambda fx, fy, fz: seen.append((fx, fy, fz)))
    assert len(seen) == 1

def test_solver_accepts_interface_state_callback():
    seen = []
    def callback(field, t, dt):
        seen.append((t, dt))
        return field
    solve_fast_obstacle_transport_3d(
        _config(duration_s=0.001, time_step_s=0.001,
                source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0))),
        interface_state_callback=callback)
    assert seen == [(0.0, 0.001)]

def test_interface_state_callback_precedes_face_flux_callback():
    events = []
    def state_callback(field, t, dt):
        events.append("state")
        return field
    def flux_callback(fx, fy, fz):
        events.append("flux")
    solve_fast_obstacle_transport_3d(
        _config(duration_s=0.001, time_step_s=0.001,
                source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0))),
        interface_state_callback=state_callback, face_flux_callback=flux_callback)
    assert events == ["state", "flux"]

def test_solver_counts_conservative_interface_inlet_flux_in_mass_ledger():
    calls = []
    def inlet_flux(t, dt, u_face, concentration, solid):
        calls.append((t, dt))
        return np.full(concentration.shape[:2], 1.0e-4)
    result = solve_fast_obstacle_transport_3d(
        _config(duration_s=0.004, time_step_s=0.001,
                source_schedule=SourceRateSchedule((0.0, 0.004), (0.0, 0.0))),
        interface_inlet_flux_callback=inlet_flux,
    )
    assert len(calls) == 4
    assert result.diagnostics["interface_inlet_flux_callback_applied"] is True
    assert result.maximum_mass_residual_kg < 1.0e-12
