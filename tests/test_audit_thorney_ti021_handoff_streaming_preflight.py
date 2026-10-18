import json
import math
from pathlib import Path

from tools.audit_thorney_ti021_handoff_streaming_preflight import (
    SCHEMA,
    build_configuration,
)


def test_streaming_preflight_freezes_actual_inner_and_handoff_coordinates():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, closure, coordinates = build_configuration(properties)

    assert SCHEMA.endswith(".v1")
    assert config.domain_m == (200.0, 200.0, 40.0)
    assert config.cells == (40, 40, 8)
    assert config.duration_s == 0.25
    assert config.scalar_time_step_s == 0.125
    assert config.scalar_advection_scheme == "muscl_van_leer"
    assert config.thermodynamic_corrector_density_tolerance == 1.0e-6
    assert config.initial_binary_cloud_momentum == "quiescent"
    assert coordinates["handoff_inner_solver_x_m"] == 160.0
    assert coordinates["source_solver_xyz_m"] == [60.0, 100.0, 1.0]
    assert coordinates[
        "source_transverse_q_equals_minus_site_x_bounds_m"
    ] == [-100.0, 100.0]
    assert coordinates[
        "nested_outer_transverse_q_equals_minus_site_x_bounds_m"
    ] == [-320.0, 320.0]
    assert closure.species_name == "R12"
    assert config.background_surface_layer_mixing is not None


def test_thorney_fence_can_use_grid_independent_zero_thickness_surface():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    coarse, _, _ = build_configuration(
        properties, cell_size_m=5.0, circular_arc_representation="zero_thickness_surface",
    )
    fine, _, _ = build_configuration(
        properties, cell_size_m=1.25, circular_arc_representation="zero_thickness_surface",
    )

    coarse_wall, = coarse.obstacles
    fine_wall, = fine.obstacles
    assert coarse.circular_arc_representation == "zero_thickness_surface"
    assert fine.circular_arc_representation == "zero_thickness_surface"
    assert coarse_wall.thickness_m == fine_wall.thickness_m == 0.0
    assert coarse_wall.center_m == fine_wall.center_m == (60.0, 100.0)
    assert coarse_wall.radius_m == fine_wall.radius_m == 50.0
    assert coarse_wall.z_min_m == fine_wall.z_min_m == 0.0
    assert coarse_wall.z_max_m == fine_wall.z_max_m == 5.0


def test_configuration_supports_sensor_free_handoff_with_shared_downwind_buffer():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, coordinates = build_configuration(
        properties, cell_size_m=5.0, time_step_s=0.125, active_steps=2,
        inner_domain_m=(300.0, 200.0, 40.0),
        handoff_site_downwind_y_m=125.0,
    )
    assert config.domain_m == (300.0, 200.0, 40.0)
    assert coordinates["handoff_site_downwind_y_m"] == 125.0
    assert coordinates["handoff_inner_solver_x_m"] == 185.0


def test_finite_source_mass_override_is_explicit_and_changes_schedule_rate():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=80,
        finite_source_release_duration_s=10.0,
        finite_source_mass_kg=4962.0,
    )
    assert config.source_schedule is not None
    assert config.source_schedule.rate_kg_s[0] == 496.2


def test_finite_source_start_delay_preserves_mass_and_declares_zero_prefix():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=80,
        finite_source_release_duration_s=5.0,
        finite_source_start_delay_s=5.0,
        finite_source_mass_kg=4962.0,
    )
    schedule = config.source_schedule
    assert schedule is not None
    assert schedule.time_s == (0.0, 5.0, 10.0)
    assert schedule.rate_kg_s == (0.0, 992.4, 0.0)
    assert schedule.mass_between(0.0, 5.0) == 0.0
    assert schedule.released_mass_kg == 4962.0


def test_finite_source_smooth_ramp_preserves_mass_and_zero_endpoints():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=80,
        finite_source_release_duration_s=5.0,
        finite_source_start_delay_s=5.0,
        finite_source_mass_kg=4962.0,
        finite_source_release_shape="smooth_ramp",
    )
    schedule = config.source_schedule
    assert schedule is not None
    assert schedule.rate_operator == "linear"
    assert schedule.time_s == (0.0, 5.0, 6.25, 7.5, 8.75, 10.0)
    assert schedule.rate_kg_s[0] == 0.0
    assert schedule.rate_kg_s[-1] == 0.0
    assert schedule.mass_between(0.0, 5.0) == 0.0
    assert schedule.released_mass_kg == 4962.0


def test_finite_source_momentum_velocity_is_explicitly_recorded():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=80,
        finite_source_release_duration_s=10.0,
        finite_source_mass_kg=4962.0,
        finite_source_momentum_velocity_m_s=(1.0, 0.0, 0.0),
    )
    assert config.source_momentum_velocity_m_s == (1.0, 0.0, 0.0)


def test_finite_source_mixture_ratio_is_explicitly_recorded():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=80,
        finite_source_release_duration_s=10.0,
        finite_source_mass_kg=4962.0,
        finite_source_total_mass_to_species_mass_ratio=1.5090569459032526,
    )
    assert config.source_total_mass_to_species_mass_ratio == 1.5090569459032526


def test_total_mass_basis_is_converted_to_tracked_species_schedule_mass():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, coordinates = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=80,
        finite_source_release_duration_s=10.0,
        finite_source_mass_kg=4962.0,
        finite_source_mass_basis="total",
        finite_source_total_mass_to_species_mass_ratio=1.5090569459032526,
    )
    schedule = config.source_schedule
    assert schedule is not None
    assert math.isclose(schedule.released_mass_kg, 4962.0 / 1.5090569459032526)
    assert coordinates["finite_source_mass_basis"] == "total"
    assert math.isclose(
        coordinates["finite_source_schedule_total_mass_kg"], 4962.0,
    )


def test_wind_prior_is_explicit_and_not_silently_replaced():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        time_step_s=0.125,
        active_steps=2,
        wind_speed_m_s=3.4,
        wind_heading_relative_deg=15.8,
    )
    assert config.wind_speed_m_s == 3.4
    assert config.active_inlet_velocity_history.history_id == (
        "TI021-explicit-wind-3.4mps-15.8deg-relative"
    )


def test_widened_inner_domain_keeps_source_obstacle_and_handoff_coordinates_centered():
    property_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
    )
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _, coordinates = build_configuration(
        properties,
        cell_size_m=5.0,
        time_step_s=0.125,
        active_steps=2,
        inner_domain_m=(450.0, 640.0, 40.0),
    )
    assert config.domain_m == (450.0, 640.0, 40.0)
    assert config.initial_binary_cloud.center_xy_m == (60.0, 320.0)
    assert config.obstacles[0].center_m == (60.0, 320.0)
    assert config.source_position_m == (60.0, 320.0, 1.0)
    assert coordinates["inner_crosswind_origin_m"] == 320.0
    assert coordinates["source_transverse_q_equals_minus_site_x_bounds_m"] == [-320.0, 320.0]
