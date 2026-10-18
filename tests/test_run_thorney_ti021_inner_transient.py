from pathlib import Path

from degali.validation.thorney_ti021 import (
    ThorneyTi021Dataset,
    Ti021DoseSensor,
)
from tools.run_thorney_ti021_inner_transient import (
    SCHEMA, _handoff_crosswind_edges, inner_receptors,
)


def test_inner_receptors_use_coordinates_only_and_freeze_solver_rotation():
    inside = Ti021DoseSensor(
        400.0, 300.0, 0.0, 100.0, 2.4, 999.0, 12.0, 50.0,
    )
    outside = Ti021DoseSensor(
        400.0, 500.1, 0.0, 300.1, 2.4, 1.0, 1.0, 2.0,
    )
    dataset = ThorneyTi021Dataset(
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
        arc_observations=(), dose_sensors=(inside, outside),
    )

    receptors, metadata = inner_receptors(dataset)

    assert SCHEMA.endswith(".v1")
    assert len(receptors) == 1
    assert receptors[0].position_m == (160.0, 100.0, 2.4)
    assert metadata[0]["site_local_xyz_m"] == [0.0, 100.0, 2.4]
    assert "dose" not in metadata[0]
    assert "arrival" not in metadata[0]


def test_inner_receptors_include_shared_downwind_sensor_rows_when_domain_extends():
    sensors = tuple(
        Ti021DoseSensor(
            400.0, float(300 + index * 50), float(index * 50), float(100 + index * 50),
            2.4, 999.0, 12.0, 50.0,
        )
        for index in range(3)
    )
    dataset = ThorneyTi021Dataset(
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
        arc_observations=(), dose_sensors=sensors,
    )
    default, _ = inner_receptors(dataset)
    extended, metadata = inner_receptors(dataset, domain_m=(300.0, 200.0, 40.0))
    wide, wide_metadata = inner_receptors(
        dataset, domain_m=(300.0, 640.0, 40.0),
    )
    assert len(default) == 1
    assert len(extended) == 3
    assert len(wide) == 3
    assert [item["site_local_xyz_m"][1] for item in metadata] == [100.0, 150.0, 200.0]
    assert [item.position_m[1] for item in extended] == [100.0, 50.0, 0.0]
    assert [item.position_m[1] for item in wide] == [320.0, 270.0, 220.0]
    assert [item["solver_xyz_m"][1] for item in wide_metadata] == [320.0, 270.0, 220.0]


def test_inner_handoff_crosswind_edges_follow_domain_width():
    narrow = _handoff_crosswind_edges(200.0, 40)
    wide = _handoff_crosswind_edges(640.0, 128)
    assert narrow[[0, -1]].tolist() == [-100.0, 100.0]
    assert wide[[0, -1]].tolist() == [-320.0, 320.0]
    assert len(wide) == 129
    assert (wide[1] - wide[0]) == 5.0
