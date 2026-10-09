import json
import math
from dataclasses import replace

import pytest

from degali.addons.field_contracts import (
    BoundedValue,
    FieldCoordinateReference,
    FieldScenario,
    FieldValidationEvidence,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_meteorology import StabilityScalarMixingClosure
from degali.addons.field_meteorology import FieldWindHistory
from degali.addons.field_distributed_source import FieldDistributedVapourSource
from degali.addons.transient_receptor import WindHistory
from degali.addons.field_report import (
    FIELD_SCREENING_REPORT_SCHEMA,
    FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA,
    field_refinement_report,
    field_sensor_array_uncertainty_envelope_report,
    field_semi_fv_envelope_report,
    field_screening_report,
    write_field_refinement_report,
    write_field_screening_report,
)
from degali.addons.field_workflow import (
    FieldSensorDeployment,
    FieldSemiFVRequest,
    run_field_semi_fv_refinement_study,
    run_field_semi_fv_envelope,
    run_field_sensor_array_uncertainty_envelope,
    run_field_semi_fv_screening,
)
from degali.addons.semi_fv_obstacle import SemiFVConfig
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid, OrientedCuboid


def _request(*, model_family="degali"):
    source = ReleaseSource(
        fluid="lh2",
        location_m=(0.0, 0.0, 0.5),
        upstream_pressure=BoundedValue(0.4e6),
        upstream_temperature=BoundedValue(26.084),
        mass_flow_kg_s=BoundedValue(0.265),
        opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8),
        discharge_coefficient=BoundedValue(0.8),
        liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium",
        duration_s=1.0,
    )
    scenario = FieldScenario(
        source=source,
        weather=WeatherState(
            speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
        ),
        sensor=SensorModel((1.0, 0.0, 0.5), averaging_time_s=0.1),
        model_family=model_family,
        temporal_mode="transient",
    )
    transport = SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        time_step_s=0.005, duration_s=1.0, source_sigma_m=0.2,
    )
    return FieldSemiFVRequest(
        scenario,
        coordinate_reference=FieldCoordinateReference(
            "plant-grid-rev-A", "transfer-skid-origin", 0.0, "site-grade", "layout-drawing-A",
        ),
        transport=transport,
    )


def test_field_report_exposes_inputs_numerics_and_conditional_scope(tmp_path):
    result = run_field_semi_fv_screening(_request())
    report = field_screening_report(result)

    assert report["schema"] == FIELD_SCREENING_REPORT_SCHEMA
    assert report["applicability"]["status"] == "conditional"
    assert report["scenario"]["sensor"]["position_m"] == [1.0, 0.0, 0.5]
    assert report["transport_result"]["diagnostics"]["x_advection_courant"] <= 1.0
    assert report["transport_result"]["diagnostics"]["source_mass_ledger_residual_kg"] == 0.0
    assert report["source_preparation"]["flash"]["conservative"]
    assert report["sensor_result"]["maximum_indicated_mole_fraction"] >= 0.0
    assert report["transport_input"]["post_release_duration_s"] == 0.0
    assert report["coordinate_reference"]["coordinate_system_id"] == "plant-grid-rev-A"
    assert report["wind_frame"]["wind_to_local_math_radians"] == report["wind_frame"]["wind_to_earth_math_radians"]

    output = write_field_screening_report(result, tmp_path / "field-report.json")
    saved = json.loads(output.read_text(encoding="utf-8"))
    assert saved["applicability"] == report["applicability"]
    with pytest.raises(FileExistsError):
        write_field_screening_report(result, output)


def test_field_report_retains_ambient_boundary_uncertainty():
    request = replace(
        _request(),
        ambient_temperature_uncertainty_k=BoundedValue(
            295.0, 290.0, 300.0, "K", "ambient-mast-A",
        ),
    )
    report = field_screening_report(run_field_semi_fv_screening(request))

    assert report["transport_input"]["ambient_boundary"]["temperature_k"] == pytest.approx(295.0)
    assert report["transport_input"]["ambient_boundary"]["uncertainty"][
        "ambient_temperature_k"
    ]["source"] == "ambient-mast-A"


def test_blocked_field_report_omits_fabricated_transport_and_sensor_values():
    result = run_field_semi_fv_screening(_request(model_family="slabx"))
    report = field_screening_report(result)

    assert report["applicability"]["status"] == "blocked"
    assert report["transport_result"] is None
    assert report["sensor_result"] is None


def test_field_report_retains_validation_evidence_fingerprint():
    evidence = FieldValidationEvidence(
        dataset_id="lh2-obstacle-trial-A",
        path="evidence/trial-A.csv",
        sha256="0" * 64,
        row_count=12,
        source_boundary_id="source-A",
        weather_id="met-A",
        obstacle_geometry_id="obstacle-A",
        receptor_geometry_id="detectors-A",
        temporal_operator_id="t90-average-A",
        common_clock_id="clock-A",
    )
    result = run_field_semi_fv_screening(
        replace(_request(), lh2_validation_available=True, validation_evidence=evidence)
    )
    report = field_screening_report(result)

    assert report["validation_evidence"] == evidence.as_record()


def test_refinement_report_attaches_resolution_evidence_to_the_same_case(tmp_path):
    result = run_field_semi_fv_refinement_study(
        _request(), refinement_factors=(1, 2), relative_tolerance=1.0,
    )
    report = field_refinement_report(result)

    assert report["numerical_refinement"] is not None
    assert report["numerical_refinement"]["estimated_cell_steps"] == 680_000
    assert len(report["numerical_refinement"]["cases"]) == 2
    assert report["refinement_applicability"]["status"] == "conditional"
    output = write_field_refinement_report(result, tmp_path / "refinement-report.json")
    assert json.loads(output.read_text(encoding="utf-8"))["numerical_refinement"] == report[
        "numerical_refinement"
    ]
    with pytest.raises(FileExistsError):
        write_field_refinement_report(result, output)


def test_field_report_records_selected_stability_mixing_provenance():
    request = _request()
    scenario = replace(
        request.scenario,
        weather=replace(request.scenario.weather, stability="stable"),
    )
    result = run_field_semi_fv_screening(
        replace(
            request,
            scenario=scenario,
            stability_mixing_closure=StabilityScalarMixingClosure(
                {"stable": 0.125}, evidence_id="site-mast-closure-2026-10-05",
            ),
        )
    )
    report = field_screening_report(result)

    closure = report["transport_input"]["stability_mixing_closure"]
    assert closure == {
        "model_id": "declared_stability_scalar_diffusivity",
        "evidence_id": "site-mast-closure-2026-10-05",
        "selected_stability": "stable",
        "selected_diffusivity_m2_s": 0.125,
    }
    assert report["transport_result"]["diagnostics"]["scalar_diffusivity_m2_s"] == 0.125


def test_field_envelope_report_retains_every_deterministic_corner_without_probability_claim():
    request = _request()
    scenario = replace(
        request.scenario,
        weather=replace(
            request.scenario.weather,
            speed_m_s=BoundedValue(2.0, 1.8, 2.2, unit="m/s"),
        ),
    )
    envelope = run_field_semi_fv_envelope(
        replace(request, scenario=scenario), max_cases=4,
    )
    report = field_semi_fv_envelope_report(envelope)

    assert report["schema"] == FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA
    assert report["completed_case_count"] == 4
    assert report["property_table_used"]
    assert len(report["cases"]) == 4
    assert "not a probability interval" in report["warnings"][-1]
    json.dumps(report, allow_nan=False)


def test_field_report_records_wind_history_as_a_static_gate_only():
    request = _request()
    result = run_field_semi_fv_screening(
        replace(
            request,
            wind_history=FieldWindHistory(
                WindHistory(
                    time_s=(0.0, 0.5, 1.0),
                    speed_m_s=(2.0, 2.1, 2.0),
                    direction_from_deg=(270.0, 271.0, 270.0),
                ),
                source_id="met-mast-01",
                evidence_id="met-mast-cal-1",
            ),
        )
    )
    report = field_screening_report(result)

    assert report["wind_history"]["source_id"] == "met-mast-01"
    assert report["wind_history"]["steady_wind_assessment"]["applicable"]
    assert report["wind_history"]["nominal_speed_relative_difference"] < 0.1
    assert report["wind_history"]["nominal_direction_difference_deg"] < 5.0
    assert "static-wind applicability gate only" in report["wind_history"]["use"]


def test_field_report_records_every_declared_obstacle_projection():
    request = _request()
    result = run_field_semi_fv_screening(
        replace(
            request,
            obstacles=(
                AxisAlignedCuboid(1.5, 2.0, -1.0, 1.0, 0.0, 1.0, "shed-a"),
                AxisAlignedCuboid(2.5, 3.0, -1.0, 1.0, 0.0, 1.0, "shed-b"),
            ),
        )
    )
    report = field_screening_report(result)

    assert report["obstacle_projection"] is None
    assert [item["local_obstacle"]["label"] for item in report["obstacle_projections"]] == [
        "shed-a", "shed-b",
    ]
    assert [item["geometry_type"] for item in report["declared_global_obstacles"]] == [
        "AxisAlignedCuboid", "AxisAlignedCuboid",
    ]


def test_field_report_preserves_oriented_obstacle_geometry_before_projection():
    request = _request()
    result = run_field_semi_fv_screening(
        replace(
            request,
            obstacles=(OrientedCuboid(
                center_m=(2.0, 0.0), length_m=2.0, width_m=1.0,
                z_min_m=0.0, z_max_m=1.0, long_axis_bearing_deg=45.0,
                label="angled-skid",
            ),),
        )
    )
    report = field_screening_report(result)

    declared = report["declared_global_obstacles"]
    assert declared[0]["geometry_type"] == "OrientedCuboid"
    assert declared[0]["geometry"]["long_axis_bearing_deg"] == 45.0
    assert report["obstacle_projections"][0]["local_obstacle"]["label"] == "angled-skid"


def test_field_report_records_declared_internal_vapour_source():
    request = _request()
    source = FieldDistributedVapourSource(
        label="pool-vapour-01", position_m=(1.0, 0.0, 0.1), vertical_sigma_m=0.15,
        schedule=SourceRateSchedule((0.0, 1.0), (0.01, 0.0), source_id="pool-ledger-01"),
        evidence_id="pool-ledger-audit-01", source_kind="pool_vapour",
    )
    report = field_screening_report(
        run_field_semi_fv_screening(replace(request, distributed_vapour_sources=(source,)))
    )

    assert report["transport_input"]["distributed_vapour_sources"] == [{
        "label": "pool-vapour-01",
        "source_kind": "pool_vapour",
        "position_m": [1.0, 0.0, 0.1],
            "vertical_sigma_m": 0.15,
            "schedule_source_id": "pool-ledger-01",
            "schedule_rate_operator": "piecewise_constant",
            "schedule_duration_s": 1.0,
        "released_mass_kg": 0.01,
        "evidence_id": "pool-ledger-audit-01",
        "warnings": [],
    }]


def test_field_report_retains_each_sensor_deployment_or_withheld_reason():
    request = replace(
        _request(),
        sensor_deployments=(
            FieldSensorDeployment("detector-02", SensorModel((2.0, 0.0, 0.5))),
            FieldSensorDeployment("detector-off-plane", SensorModel((2.0, 1.0, 0.5))),
        ),
    )
    report = field_screening_report(run_field_semi_fv_screening(request))

    assert [item["label"] for item in report["transport_input"]["additional_sensor_deployments"]] == [
        "detector-02", "detector-off-plane",
    ]
    sensors = {item["label"]: item for item in report["sensor_results"]}
    assert sensors["field_sensor"]["withheld"] is False
    assert sensors["detector-02"]["withheld"] is False
    assert sensors["detector-off-plane"]["withheld"] is True
    assert "outside the local wind plane" in sensors["detector-off-plane"]["warnings"][0]


def test_sensor_array_uncertainty_report_reuses_the_base_transport_audit():
    request = replace(
        _request(),
        sensor_deployments=(FieldSensorDeployment(
            "detector-uncertain",
            SensorModel(
                (2.0, 0.0, 0.5),
                response_time_s=BoundedValue(0.3, 0.1, 0.5),
            ),
        ),),
    )
    envelope = run_field_sensor_array_uncertainty_envelope(request, max_cases=4)
    report = field_sensor_array_uncertainty_envelope_report(envelope)

    assert report["completed_case_count"] == 2
    assert len(report["cases"]) == 2
    assert report["base_field_screening"]["transport_result"] is not None
    assert all(len(case["sensor_results"]) == 2 for case in report["cases"])
    json.dumps(report, allow_nan=False)
