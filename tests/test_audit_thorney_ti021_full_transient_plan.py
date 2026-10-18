from pathlib import Path

import pytest

from degali.validation.thorney_ti021 import (
    ThorneyTi021Dataset,
    Ti021ArcObservation,
    Ti021DoseSensor,
)
from tools.audit_thorney_ti021_full_transient_plan import (
    SCHEMA,
    build_plan,
    site_local_to_nested_outer,
    site_local_to_solver,
)


def _dataset() -> ThorneyTi021Dataset:
    arc_sensors = [
        Ti021DoseSensor(
            400.0, 200.0 + distance, 0.0, distance, height,
            0.0, None, None,
        )
        for distance in (75.0, 100.0)
        for height in (0.4, 2.4, 4.4)
    ]
    remaining = [
        Ti021DoseSensor(
            400.0 + (-300.0 + 600.0 * index / 138.0), 700.0,
            -300.0 + 600.0 * index / 138.0, 500.0,
            0.4 if index % 2 == 0 else 2.4, 0.0, None, None,
        )
        for index in range(139)
    ]
    arcs = tuple(
        Ti021ArcObservation(
            distance, height, 0.0, 0.0, 0.0, 1.0, 0.0, None,
        )
        for distance in (75.0, 100.0)
        for height in (0.4, 2.4, 4.4)
    )
    return ThorneyTi021Dataset(
        workbook_path=Path("ti21.xls"), workbook_sha256="a" * 64,
        release_origin_site_m=(400.0, 200.0, 0.0), substance="R12",
        release_type="puff", initial_r12_mole_fraction=0.32,
        effective_molecular_weight_kg_mol=0.0585,
        initial_temperature_k=291.15, initial_pressure_pa=100_500.0,
        workbook_mass_proxy_kg=4_962.0, workbook_mass_proxy_duration_s=1.0,
        workbook_wind_speed_m_s=3.4, workbook_wind_reference_height_m=10.0,
        workbook_wind_from_north_deg=194.2, ideal_wind_from_north_deg=210.0,
        site_y_axis_from_north_deg=30.0, friction_velocity_m_s=0.205,
        roughness_length_m=0.005, stability_class="D/E",
        workbook_wall_radius_m=50.0,
        workbook_wall_center_site_m=(400.0, 200.0),
        workbook_wall_height_m=6.0, workbook_porosity_field=1.0,
        arc_observations=arcs,
        dose_sensors=tuple(arc_sensors + remaining),
    )


def _benchmark(dx: float, cells: int, steps: int, runtime: float):
    return {
        "configuration": {
            "cell_size_m": dx, "cell_count": cells, "active_steps": steps,
        },
        "runtime": {"total_s": runtime},
    }


def test_site_to_solver_rotation_matches_frozen_wind_frame():
    assert site_local_to_solver(
        -14.0, 14.0, 0.4, source_solver_xyz_m=(60.0, 320.0, 0.0),
    ) == pytest.approx((74.0, 334.0, 0.4))
    assert site_local_to_nested_outer(
        -14.0, 100.0, 0.4,
        handoff_downwind_y_m=100.0,
        outer_crosswind_origin_m=320.0,
    ) == pytest.approx((0.0, 334.0, 0.4))


def test_full_transient_plan_exposes_preflight_coverage_and_freezes_nested_domain():
    record = build_plan(
        _dataset(),
        dx2_benchmark=_benchmark(2.0, 36_400, 40, 60.0),
        dx1_benchmark=_benchmark(1.0, 291_200, 1, 34.0),
    )

    assert record["schema"] == SCHEMA
    assert record["status"] == (
        "accepted_coordinate_consistent_nested_execution_plan"
    )
    assert record["all_gates_pass"]
    assert record["coverage"]["current_preflight"]["covered_sensor_count"] < 145
    assert record["coverage"]["monolithic_sensor_domain"][
        "covered_sensor_count"
    ] == 145
    assert record["coverage"]["nested_union"][
        "combined_unique_sensor_count"
    ] == 145
    assert record["coverage"]["nested_union"]["overlap_sensor_count"] > 0
    assert record["coverage"]["nested_union"][
        "handoff_inner_solver_x_m"
    ] == 160.0
    assert record["coverage"]["nested_outer_sensor_domain"][
        "local_site_bounds_m"
    ]["downwind_y"] == [100.0, 540.0]
    assert record["coverage"]["resolved_inner"][
        "covered_arc_centerline_count"
    ] == 6
    inner = record["runtime_model"]["resolved_inner_three_grid"]
    assert [item["wall_height_cell_count"] for item in inner] == [1.0, 2.0, 4.0]
    assert record["runtime_model"]["monolithic_all_sensor_three_grid"][-1][
        "runtime_estimate_hours"
    ]["lower"] > 1_000.0
    assert record["runtime_model"]["nested_outer_three_grid"][-1][
        "runtime_estimate_hours"
    ]["lower"] < record["runtime_model"]["monolithic_all_sensor_three_grid"][-1][
        "runtime_estimate_hours"
    ]["lower"]
    assert record["dataset"]["observation_value_fields_accessed_for_planning"] == []
