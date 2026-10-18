import numpy as np

from degali.addons.continuous_barrier_resolved_3d import (
    audit_continuous_barrier_physics,
    build_continuous_barrier_resolved_config,
)
from degali.addons.resolved_obstacle_flow_3d import (
    solve_resolved_obstacle_flow_3d,
)


def test_continuous_barrier_uses_coupled_projected_momentum_path():
    config, barrier = build_continuous_barrier_resolved_config(
        duration_s=30.0,
        release_duration_s=30.0,
        barrier_layout="blocking-wall",
        store_scalar_fields=False,
    )
    assert config.active_momentum_solver == "finite_volume_projection"
    assert config.carrier_flow_initialization == "projected_inlet"
    assert config.thermal_buoyancy_coupling is True
    assert config.low_mach_mass_flux_correction is True
    assert config.cells == (28, 16, 8)
    assert np.allclose(
        np.asarray(config.domain_m) / np.asarray(config.cells), 5.0
    )
    assert barrier["height_m"] / 5.0 == 3.0

    result = solve_resolved_obstacle_flow_3d(config)
    audit = audit_continuous_barrier_physics(result, barrier)

    assert audit["all_internal_physics_checks_passed"] is True
    assert all(audit["checks"].values())
    metrics = audit["metrics"]
    assert metrics["windward_slowdown_fraction"] > 0.15
    assert metrics["windward_to_lee_pressure_gain_pa"] > 0.0
    assert metrics["roof_maximum_streamwise_velocity_m_s"] > (
        metrics["reference_streamwise_velocity_m_s"]
    )
    assert metrics["roof_maximum_vertical_streamwise_shear_s_inv"] > 0.02
    assert metrics["lee_minimum_streamwise_velocity_m_s"] < 0.0
    assert metrics["lee_minimum_vertical_velocity_m_s"] < 0.0
    assert metrics["maximum_wall_normal_velocity_m_s"] < 1.0e-10
    assert metrics["maximum_relative_mass_residual"] < 1.0e-9
    assert metrics["maximum_continuity_residual_kg_m3_s"] < 1.0e-8
    routes = audit["final_downstream_transport_routes"]
    assert routes["plane_x_m"] == 80.0
    assert routes["direct_through_wall_projection_kg_s"] == 0.0
    assert routes["over_wall_top_kg_s"] >= 0.0
    assert routes["around_sides_or_open_gap_kg_s"] >= 0.0
    assert routes["total_downstream_kg_s"] == (
        routes["over_wall_top_kg_s"]
        + routes["around_sides_or_open_gap_kg_s"]
    )


def test_building_cluster_uses_four_grounded_cuboids_on_cubic_grid():
    config, geometry = build_continuous_barrier_resolved_config(
        duration_s=1.0,
        release_duration_s=1.0,
        barrier_layout="building-cluster",
        store_scalar_fields=False,
    )
    assert config.domain_m == (160.0, 120.0, 50.0)
    assert config.cells == (32, 24, 10)
    assert np.allclose(
        np.asarray(config.domain_m) / np.asarray(config.cells), 5.0
    )
    assert len(config.obstacles) == 4
    assert len(geometry["buildings"]) == 4
    assert min(block["height_m"] for block in geometry["buildings"]) == 15.0
    assert max(block["height_m"] for block in geometry["buildings"]) == 30.0

