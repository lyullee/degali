import numpy as np
import pytest

from degali.addons.fast_obstacle_transport import FastObstacleReceptor
from degali.addons.resolved_obstacle_flow_2d import (
    RESOLVED_OBSTACLE_FLOW_2D_SCHEMA,
    ResolvedObstacleFlow2DConfig,
    solve_resolved_obstacle_flow_2d,
)
from degali.addons.semi_fv_obstacle import (
    RectangularObstacle2D,
    SourceRateSchedule,
)


def _config(**changes):
    values = {
        "length_m": 8.0,
        "height_m": 3.0,
        "nx": 80,
        "nz": 30,
        "wind_speed_m_s": 2.0,
        "kinematic_viscosity_m2_s": 0.3,
        "flow_steps": 800,
        "lattice_inlet_speed": 0.03,
        "duration_s": 1.5,
        "scalar_time_step_s": 0.005,
        "scalar_diffusivity_m2_s": 0.03,
        "source_schedule": SourceRateSchedule(
            (0.0, 0.5, 1.5), (0.02, 0.0, 0.0),
        ),
        "source_x_m": 0.6,
        "source_height_m": 0.5,
        "source_sigma_m": 0.2,
        "obstacle": RectangularObstacle2D(2.75, 0.7, 1.5, label="block"),
        "minimum_flow_steps": 500,
    }
    values.update(changes)
    return ResolvedObstacleFlow2DConfig(**values)


@pytest.fixture(scope="module")
def resolved_result():
    return solve_resolved_obstacle_flow_2d(
        _config(),
        receptors=(
            FastObstacleReceptor("windward", 2.3, 0.8),
            FastObstacleReceptor("roof", 2.75, 1.7),
            FastObstacleReceptor("lee", 3.4, 0.8),
        ),
    )


def test_resolved_flow_has_no_slip_obstacle_and_emergent_lee_recirculation(
    resolved_result,
):
    result = resolved_result
    diagnostics = result.diagnostics

    assert result.schema == RESOLVED_OBSTACLE_FLOW_2D_SCHEMA
    assert np.all(result.velocity_u_m_s[result.solid_mask] == 0.0)
    assert np.all(result.velocity_w_m_s[result.solid_mask] == 0.0)
    assert diagnostics["recirculating_lee_cell_count"] > 0
    assert diagnostics["minimum_lee_streamwise_velocity_m_s"] < 0.0
    assert diagnostics["maximum_absolute_divergence_s_1"] < 1.0e-9
    assert diagnostics["nonlocal_scalar_routing"] is False


def test_resolved_scalar_transport_is_local_impermeable_and_mass_conservative(
    resolved_result,
):
    result = resolved_result
    solid = result.solid_mask

    blocked_x = solid[:, :-1] | solid[:, 1:]
    blocked_z = solid[:-1, :] | solid[1:, :]
    assert np.all(result.u_face_m_s[:, 1:-1][blocked_x] == 0.0)
    assert np.all(result.w_face_m_s[1:-1, :][blocked_z] == 0.0)
    assert result.maximum_mass_residual_per_span_kg_m < 1.0e-12
    assert result.injected_mass_per_span_kg_m[-1] == pytest.approx(0.01)
    assert np.nanmax(result.concentration_vol_percent) > 0.0
    assert np.all(
        result.maximum_concentration_kg_m3[~solid]
        >= result.concentration_kg_m3[~solid]
    )
    assert np.nanmax(result.maximum_concentration_vol_percent) >= np.nanmax(
        result.concentration_vol_percent
    )
    assert set(result.receptor_traces_kg_m3) == {"windward", "roof", "lee"}


def test_resolved_result_exposes_pressure_vorticity_runtime_and_scope(resolved_result):
    result = resolved_result
    diagnostics = result.diagnostics

    assert np.isfinite(result.pressure_coefficient[~result.solid_mask]).all()
    assert np.isfinite(result.vorticity_s_1[~result.solid_mask]).all()
    assert diagnostics["fidelity_tier"] == "tier_2_resolved_viscous_2d"
    assert diagnostics["quantitative_prediction_allowed"] is False
    assert diagnostics["lattice_mach_number"] < 0.1
    assert diagnostics["maximum_relative_lattice_density_variation"] < 0.1
    assert diagnostics["flow_solve_runtime_s"] > 0.0
    assert diagnostics["total_runtime_s"] >= diagnostics["flow_solve_runtime_s"]


def test_resolved_config_enforces_d2q9_and_low_mach_stability_limits():
    with pytest.raises(ValueError, match="equal x and z"):
        _config(nz=29)
    with pytest.raises(ValueError, match="below 0.1"):
        _config(lattice_inlet_speed=0.1)
    with pytest.raises(ValueError, match="relaxation time"):
        solve_resolved_obstacle_flow_2d(
            _config(kinematic_viscosity_m2_s=0.01, flow_steps=1, minimum_flow_steps=1)
        )


def test_resolved_solver_preserves_a_declared_measured_inlet_profile():
    result = solve_resolved_obstacle_flow_2d(_config(
        flow_steps=10,
        minimum_flow_steps=10,
        duration_s=0.01,
        scalar_time_step_s=0.001,
        source_schedule=SourceRateSchedule((0.0, 0.01), (0.0, 0.0)),
        inlet_profile_z_m=(0.0, 1.5, 3.0),
        inlet_profile_u_m_s=(0.0, 1.5, 3.0),
    ))

    expected = np.interp(
        result.z_m, (0.0, 1.5, 3.0), (0.0, 1.5, 3.0),
    )
    assert result.u_face_m_s[:, 0] == pytest.approx(expected)
    assert "measured/interpolated" in result.diagnostics["inlet_boundary"]


def test_resolved_config_rejects_incomplete_or_nonmonotone_inlet_profiles():
    with pytest.raises(ValueError, match="supplied together"):
        _config(inlet_profile_z_m=(0.0, 3.0))
    with pytest.raises(ValueError, match="strictly increasing"):
        _config(
            inlet_profile_z_m=(0.0, 2.0, 1.0),
            inlet_profile_u_m_s=(0.0, 2.0, 1.0),
        )


def test_resolved_solver_rejects_a_receptor_inside_the_obstacle():
    with pytest.raises(ValueError, match="stencil intersects the obstacle"):
        solve_resolved_obstacle_flow_2d(
            _config(flow_steps=1, minimum_flow_steps=1),
            receptors=(FastObstacleReceptor("solid", 2.75, 0.5),),
        )
