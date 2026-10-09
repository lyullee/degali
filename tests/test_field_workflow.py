import math
from dataclasses import replace

import numpy as np
import pytest
import degali.addons.field_workflow as field_workflow

from degali.addons.field_contracts import (
    BoundedValue,
    FieldCoordinateReference,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
)
from degali.addons.field_meteorology import StabilityScalarMixingClosure
from degali.addons.field_meteorology import FieldWindHistory
from degali.addons.field_distributed_source import FieldDistributedVapourSource
from degali.addons.field_source_io import (
    FieldAtmosphericSourceSchedule,
    FieldSourceScheduleEvidence,
)
from degali.addons.field_geometry import FieldObstacleGeometryUncertainty
from degali.addons.transient_receptor import WindHistory
from degali.addons.field_workflow import (
    FieldSensorDeployment,
    FieldSensorDeploymentResult,
    FieldSensorArrayUncertaintyCase,
    FieldSensorArrayUncertaintyEnvelope,
    FieldSemiFVEnvelope,
    FieldSemiFVEnvelopeCase,
    FieldSemiFVRefinementResult,
    FieldSemiFVRequest,
    run_field_semi_fv_envelope,
    run_field_sensor_array_uncertainty_envelope,
    run_field_semi_fv_refinement_study,
    run_field_semi_fv_screening,
)
from degali.addons.field_report import field_semi_fv_envelope_report
from degali.addons.field_decision import evaluate_field_operational_screening
from degali.addons.semi_fv_obstacle import SemiFVConfig, SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


def _scenario(*, sensor=True, model_family="degali"):
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
    return FieldScenario(
        source=source,
        weather=WeatherState(
            speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
        ),
        sensor=(
            SensorModel((1.0, 0.0, 0.5), response_time_s=BoundedValue(0.2))
            if sensor else None
        ),
        model_family=model_family,
        temporal_mode="transient",
    )


def _transport():
    return SemiFVConfig(
        length_m=10.0, height_m=5.0, nx=80, nz=40,
        duration_s=3.0, time_step_s=0.005, source_sigma_m=0.25,
    )


def test_field_workflow_routes_direct_flash_vapour_to_transport_and_sensor():
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )
    assert result.completed
    assert result.source_preparation is not None
    assert result.transport is not None
    assert result.sensor_trace is not None
    assert np.max(result.sensor_trace.true_mole_fraction) > 0.0
    assert result.applicability.status == "conditional"
    assert any("post-flash liquid" in item for item in result.applicability.warnings)


def test_field_runtime_gate_blocks_final_inventory_residual():
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )
    assert result.transport is not None
    tampered_diagnostics = replace(
        result.transport.diagnostics,
        final_mass_residual_kg=0.1,
    )
    tampered_transport = replace(
        result.transport,
        diagnostics=tampered_diagnostics,
    )
    applicability = field_workflow._transport_conservation_applicability(
        tampered_transport
    )
    assert applicability.status == "blocked"
    assert any("final mass-conservation residual" in reason for reason in applicability.reasons)


def test_field_request_rejects_untyped_scenario_and_transport():
    with pytest.raises(TypeError, match="scenario must be"):
        FieldSemiFVRequest("not-a-scenario")
    with pytest.raises(TypeError, match="transport must be"):
        FieldSemiFVRequest(_scenario(), transport="not-a-transport")


def test_field_screening_result_contract_rejects_malformed_direct_records():
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )

    with pytest.raises(TypeError, match="applicability"):
        replace(result, applicability="conditional")
    with pytest.raises(ValueError, match="sensor traces require"):
        replace(result, transport=None)
    with pytest.raises(ValueError, match="sensor labels must be unique"):
        replace(result, sensor_results=(result.sensor_results[0], result.sensor_results[0]))
    with pytest.raises(ValueError, match="withheld field sensor result"):
        FieldSensorDeploymentResult(
            "detector-withheld", SensorModel((2.0, 0.0, 0.5)), None,
        )


def test_field_envelope_and_refinement_contracts_reject_malformed_records():
    screening = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )
    with pytest.raises(TypeError, match="cannot be boolean"):
        FieldSemiFVEnvelopeCase((("wind", True),), screening)
    with pytest.raises(TypeError, match="corner values must be a tuple"):
        FieldSemiFVEnvelopeCase([("wind", 1.0)], screening)
    with pytest.raises(TypeError, match="cases must be a tuple"):
        FieldSemiFVEnvelope(
            screening.request,
            [FieldSemiFVEnvelopeCase((), screening)],
            False,
        )
    with pytest.raises(ValueError, match="at least one case"):
        FieldSemiFVEnvelope(screening.request, (), False)
    duplicate_case = FieldSemiFVEnvelopeCase((('corner', 1.0),), screening)
    with pytest.raises(ValueError, match="corner selections must be unique"):
        FieldSemiFVEnvelope(
            screening.request, (duplicate_case, duplicate_case), False,
        )
    with pytest.raises(TypeError, match="refinement screening"):
        FieldSemiFVRefinementResult("not-a-screening", None, screening.applicability)
    with pytest.raises(ValueError, match="needs at least one sensor result"):
        FieldSensorArrayUncertaintyCase((("gain", 1.0),), ())
    with pytest.raises(TypeError, match="screening must be"):
        FieldSensorArrayUncertaintyEnvelope("not-a-screening", (), ("withheld",))


def test_field_workflow_withholds_nonconservative_transport_diagnostics(monkeypatch):
    from degali.addons import field_workflow

    original_solver = field_workflow.solve_semi_fv_obstacle

    def nonconservative_solver(*args, **kwargs):
        result = original_solver(*args, **kwargs)
        diagnostics = replace(
            result.diagnostics,
            maximum_mass_residual_kg=1.0,
        )
        return replace(result, diagnostics=diagnostics)

    monkeypatch.setattr(field_workflow, "solve_semi_fv_obstacle", nonconservative_solver)
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )

    assert result.transport is not None
    assert not result.completed
    assert result.applicability.status == "blocked"
    assert any("mass-conservation residual" in reason
               for reason in result.applicability.reasons)
    from degali.addons.field_decision import evaluate_field_operational_screening
    decision = evaluate_field_operational_screening(
        result, require_refinement=False,
    )
    assert "transport_mass_residual_exceeded" in decision.gate_codes


def test_field_workflow_withholds_an_unclosed_source_mass_ledger(monkeypatch):
    from degali.addons import field_workflow

    original_solver = field_workflow.solve_semi_fv_obstacle

    def unclosed_ledger_solver(*args, **kwargs):
        result = original_solver(*args, **kwargs)
        diagnostics = replace(
            result.diagnostics,
            source_mass_ledger_residual_kg=1.0,
        )
        return replace(result, diagnostics=diagnostics)

    monkeypatch.setattr(field_workflow, "solve_semi_fv_obstacle", unclosed_ledger_solver)
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )

    assert result.transport is not None
    assert not result.completed
    assert result.applicability.status == "blocked"
    assert any("source-mass ledger residual" in reason
               for reason in result.applicability.reasons)
    from degali.addons.field_decision import evaluate_field_operational_screening
    decision = evaluate_field_operational_screening(
        result, require_refinement=False,
    )
    assert "source_ledger_mismatch" in decision.gate_codes


def test_field_workflow_withholds_a_source_schedule_mass_mismatch(monkeypatch):
    from degali.addons import field_workflow

    original_solver = field_workflow.solve_semi_fv_obstacle

    def mismatched_schedule_solver(*args, **kwargs):
        result = original_solver(*args, **kwargs)
        diagnostics = replace(
            result.diagnostics,
            maximum_source_mass_schedule_residual_kg=1.0,
        )
        return replace(result, diagnostics=diagnostics)

    monkeypatch.setattr(field_workflow, "solve_semi_fv_obstacle", mismatched_schedule_solver)
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(), transport=_transport())
    )

    assert result.transport is not None
    assert not result.completed
    assert result.applicability.status == "blocked"
    assert any("source-mass schedule residual" in reason
               for reason in result.applicability.reasons)
    decision = evaluate_field_operational_screening(
        result, require_refinement=False,
    )
    assert "source_schedule_mass_mismatch" in decision.gate_codes


def test_field_workflow_withholds_unbacked_validation_availability_flag():
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(), transport=_transport(), lh2_validation_available=True,
        )
    )

    assert not result.completed
    assert result.transport is None
    assert any("requires an explicit FieldValidationEvidence" in reason
               for reason in result.applicability.reasons)


def test_field_workflow_rotates_meteorological_wind_into_a_declared_site_grid():
    reference = FieldCoordinateReference(
        "plant-grid-rev-A", "transfer-skid-origin", 90.0, "site-grade", "layout-drawing-A",
    )
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(sensor=False), coordinate_reference=reference, transport=_transport(),
        )
    )

    assert result.completed
    assert result.wind_frame is not None
    assert result.wind_frame.direction_rad == pytest.approx(1.5 * math.pi)


def test_field_workflow_evaluates_named_sensor_array_and_withholds_off_plane_detector():
    request = FieldSemiFVRequest(
        _scenario(),
        transport=_transport(),
        sensor_deployments=(
            FieldSensorDeployment("detector-downwind", SensorModel((2.0, 0.0, 0.5))),
            FieldSensorDeployment("detector-off-plane", SensorModel((2.0, 1.0, 0.5))),
        ),
    )
    result = run_field_semi_fv_screening(request)

    assert result.completed
    assert result.transport is not None
    assert result.transport.receptor_traces[0].receptor.label == "field_sensor"
    assert result.transport.diagnostics is not None
    assert len(result.transport.receptor_traces) == 2
    by_label = {item.label: item for item in result.sensor_results}
    assert result.sensor_trace is by_label["field_sensor"].trace
    assert by_label["detector-downwind"].trace is not None
    assert by_label["detector-off-plane"].trace is None
    assert "outside the local wind plane" in by_label["detector-off-plane"].warnings[0]
    with pytest.raises(ValueError, match="reserved"):
        FieldSemiFVRequest(
            _scenario(),
            sensor_deployments=(FieldSensorDeployment("field_sensor", SensorModel((2.0, 0.0, 0.5))),),
        )
    with pytest.raises(ValueError, match="sensor-array uncertainty study"):
        run_field_semi_fv_envelope(FieldSemiFVRequest(
            _scenario(),
            sensor_deployments=(FieldSensorDeployment(
                "uncertain-detector",
                SensorModel((2.0, 0.0, 0.5), response_time_s=BoundedValue(1.0, 0.5, 1.5)),
            ),),
        ))


def test_sensor_array_uncertainty_reuses_one_transport_and_varies_only_calibration():
    request = FieldSemiFVRequest(
        _scenario(),
        transport=_transport(),
        sensor_deployments=(FieldSensorDeployment(
            "detector-uncertain",
            SensorModel(
                (2.0, 0.0, 0.5),
                response_time_s=BoundedValue(0.3, 0.1, 0.5),
            ),
        ),),
    )
    envelope = run_field_sensor_array_uncertainty_envelope(request, max_cases=4)

    assert envelope.screening.completed
    assert len(envelope.cases) == 2
    assert {
        case.corner["detector-uncertain.sensor_response_time_s"]
        for case in envelope.cases
    } == {0.1, 0.5}
    first, second = envelope.cases
    first_by_label = {item.label: item for item in first.sensor_results}
    second_by_label = {item.label: item for item in second.sensor_results}
    assert first_by_label["detector-uncertain"].trace is not None
    assert second_by_label["detector-uncertain"].trace is not None
    assert np.array_equal(
        first_by_label["detector-uncertain"].trace.true_mole_fraction,
        second_by_label["detector-uncertain"].trace.true_mole_fraction,
    )
    assert not np.array_equal(
        first_by_label["detector-uncertain"].trace.indicated_mole_fraction,
        second_by_label["detector-uncertain"].trace.indicated_mole_fraction,
    )
    with pytest.raises(ValueError, match="withheld detectors"):
        run_field_sensor_array_uncertainty_envelope(FieldSemiFVRequest(
            _scenario(),
            transport=_transport(),
            sensor_deployments=(FieldSensorDeployment(
                "off-plane-uncertain",
                SensorModel(
                    (2.0, 1.0, 0.5),
                    response_time_s=BoundedValue(0.3, 0.1, 0.5),
                ),
            ),),
        ))


def test_sensor_array_uncertainty_envelope_rejects_empty_or_duplicate_completed_cases():
    request = FieldSemiFVRequest(
        _scenario(),
        transport=_transport(),
        sensor_deployments=(FieldSensorDeployment(
            "detector-uncertain",
            SensorModel(
                (2.0, 0.0, 0.5),
                response_time_s=BoundedValue(0.3, 0.1, 0.5),
            ),
        ),),
    )
    envelope = run_field_sensor_array_uncertainty_envelope(request, max_cases=4)

    with pytest.raises(ValueError, match="requires at least one case"):
        replace(envelope, cases=())
    with pytest.raises(ValueError, match="corner selections must be unique"):
        replace(envelope, cases=(envelope.cases[0], envelope.cases[0]))


def test_field_workflow_projects_directional_obstacle_and_withholds_off_plane_sensor():
    scenario = _scenario(sensor=True)
    moved_sensor = SensorModel((1.0, 2.0, 0.5))
    scenario = FieldScenario(
        source=scenario.source, weather=scenario.weather, surface=scenario.surface,
        sensor=moved_sensor, temporal_mode="transient",
    )
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            scenario,
            transport=_transport(),
            obstacle=AxisAlignedCuboid(3.0, 4.0, -1.0, 1.0, 0.0, 2.0),
            lateral_sensor_tolerance_m=0.01,
        )
    )
    assert result.completed
    assert result.obstacle_projection is not None
    assert result.obstacle_projection.obstacle is not None
    assert result.sensor_trace is None
    assert any("outside the local wind plane" in item for item in result.applicability.warnings)


def test_field_workflow_projects_multiple_facility_obstacles_into_one_conservative_mask():
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(sensor=False),
            transport=_transport(),
            obstacles=(
                AxisAlignedCuboid(2.0, 3.0, -1.0, 1.0, 0.0, 1.0, "shed-a"),
                AxisAlignedCuboid(5.0, 6.0, -1.0, 1.0, 0.0, 1.5, "shed-b"),
            ),
        )
    )

    assert result.completed
    assert result.obstacle_projection is None  # plural path has no fabricated primary obstacle
    assert len(result.obstacle_projections) == 2
    assert result.transport is not None
    assert result.transport.diagnostics.obstacle_count == 2
    assert result.transport.diagnostics.maximum_mass_residual_kg < 1.0e-9
    from degali.addons.field_decision import evaluate_field_operational_screening
    decision = evaluate_field_operational_screening(
        result,
        require_refinement=False,
        require_in_plane_sensor=False,
        allow_conditional=True,
    )
    assert "obstacle_transport_conditional" in decision.gate_codes

    with pytest.raises(ValueError, match="either obstacle or obstacles"):
        FieldSemiFVRequest(
            _scenario(sensor=False),
            obstacle=AxisAlignedCuboid(2.0, 3.0, -1.0, 1.0, 0.0, 1.0, "shed-a"),
            obstacles=(AxisAlignedCuboid(5.0, 6.0, -1.0, 1.0, 0.0, 1.0, "shed-b"),),
        )


def test_field_envelope_propagates_bounded_obstacle_geometry_corners():
    obstacle = AxisAlignedCuboid(2.0, 3.0, -1.0, 1.0, 0.0, 1.0, "shed-a")
    geometry = FieldObstacleGeometryUncertainty(
        obstacle,
        (
            ("x_min_m", BoundedValue(2.0, 1.8, 2.2, "m", "layout-A")),
            ("x_max_m", BoundedValue(3.0, unit="m", source="layout-A")),
            ("y_min_m", BoundedValue(-1.0, unit="m", source="layout-A")),
            ("y_max_m", BoundedValue(1.0, unit="m", source="layout-A")),
            ("z_min_m", BoundedValue(0.0, unit="m", source="layout-A")),
            ("z_max_m", BoundedValue(1.0, unit="m", source="layout-A")),
        ),
        "layout-revision-A",
    )
    request = FieldSemiFVRequest(
        _scenario(sensor=False),
        transport=_transport(),
        obstacle=obstacle,
        obstacle_geometry_uncertainty=(geometry,),
    )

    envelope = run_field_semi_fv_envelope(request, max_cases=4)

    assert len(envelope.cases) == 4  # two obstacle corners x two default surface corners
    assert {
        case.corner["obstacle.shed-a.x_min_m"]
        for case in envelope.cases
    } == {1.8, 2.2}
    assert all(case.result.request.obstacle_geometry_uncertainty == () for case in envelope.cases)
    report = field_semi_fv_envelope_report(envelope)
    assert report["declared_obstacle_geometry_uncertainty"][0]["evidence_id"] == "layout-revision-A"


def test_field_envelope_propagates_ambient_boundary_corners_and_selection():
    request = FieldSemiFVRequest(
        _scenario(),
        transport=_transport(),
        ambient_temperature_uncertainty_k=BoundedValue(
            295.0, 290.0, 300.0, "K", "ambient-mast-A",
        ),
    )

    envelope = run_field_semi_fv_envelope(request, max_cases=8)

    assert len(envelope.cases) == 4  # ambient temperature x default surface corners
    assert {
        case.corner["ambient_temperature_k"] for case in envelope.cases
    } == {290.0, 300.0}
    assert all(
        case.result.request.ambient_temperature_uncertainty_k is None
        for case in envelope.cases
    )
    report = field_semi_fv_envelope_report(envelope)
    assert report["cases"][0]["selection"]["ambient_temperature_k"] in {290.0, 300.0}


@pytest.mark.parametrize(
    ("field_name", "unit", "message"),
    [
        ("ambient_temperature_uncertainty_k", "Pa", "ambient_temperature_uncertainty_k"),
        ("ambient_pressure_uncertainty_pa", "K", "ambient_pressure_uncertainty_pa"),
        ("ambient_air_density_uncertainty_kg_m3", "m", "ambient_air_density_uncertainty_kg_m3"),
    ],
)
def test_field_request_rejects_ambient_uncertainty_unit_mismatch(
    field_name, unit, message,
):
    with pytest.raises(ValueError, match=message):
        FieldSemiFVRequest(
            _scenario(),
            transport=_transport(),
            **{field_name: BoundedValue(
                {
                    "ambient_temperature_uncertainty_k": 295.0,
                    "ambient_pressure_uncertainty_pa": 101325.0,
                    "ambient_air_density_uncertainty_kg_m3": 1.2,
                }[field_name],
                unit=unit,
                source="ambient-unit-test",
            )},
        )


def test_field_envelope_builds_table_over_ambient_pressure_bounds(monkeypatch):
    from degali.addons.field_lh2 import build_lh2_saturation_table_for_release as build_table

    captured = {}

    def wrapped(release, **kwargs):
        captured.update(kwargs)
        return build_table(release, **kwargs)

    monkeypatch.setattr(field_workflow, "build_lh2_saturation_table_for_release", wrapped)
    request = FieldSemiFVRequest(
        _scenario(),
        transport=_transport(),
        ambient_pressure_uncertainty_pa=BoundedValue(
            101325.0, 80_000.0, 120_000.0, "Pa", "ambient-mast-A",
        ),
    )

    envelope = run_field_semi_fv_envelope(request, max_cases=4)

    assert len(envelope.cases) == 4
    assert captured["ambient_pressure_bounds_pa"] == (80_000.0, 120_000.0)


def test_operational_gate_withholds_nominal_result_with_unpropagated_ambient_bounds():
    request = FieldSemiFVRequest(
        _scenario(),
        transport=_transport(),
        ambient_air_density_uncertainty_kg_m3=BoundedValue(
            1.2, 1.1, 1.3, "kg/m3", "ambient-mast-A",
        ),
    )
    screening = run_field_semi_fv_screening(request)

    decision = evaluate_field_operational_screening(
        screening,
        require_refinement=False,
        require_resolved_uncertainty=True,
    )

    assert not decision.screening_allowed
    assert any("ambient_air_density_kg_m3" in reason for reason in decision.reasons)


def test_field_workflow_does_not_present_itself_as_slabx():
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(_scenario(model_family="slabx"), transport=_transport())
    )
    assert not result.completed
    assert result.applicability.status == "blocked"
    assert "not a SLABx solver" in result.applicability.reasons[0]


def test_field_refinement_reuses_the_declared_field_source_and_sensor_plane():
    result = run_field_semi_fv_refinement_study(
        FieldSemiFVRequest(
            _scenario(),
            transport=SemiFVConfig(
                length_m=4.0, height_m=2.0, nx=10, nz=10,
                duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
            ),
        ),
        refinement_factors=(1, 2), relative_tolerance=1.0,
    )

    assert result.screening.completed
    assert result.study is not None
    assert result.study.estimated_cell_steps == 340_000
    assert result.study.receptor_changes[0].receptor_label == "field_sensor"
    # Physical-validation warnings from the field screen remain visible even
    # when the deliberately loose numerical check converges.
    assert result.applicability.status == "conditional"


def test_field_refinement_blocks_without_a_declared_sensor():
    result = run_field_semi_fv_refinement_study(
        FieldSemiFVRequest(_scenario(sensor=False), transport=_transport())
    )
    assert not result.completed
    assert result.applicability.status == "blocked"
    assert "requires a declared sensor" in result.applicability.reasons[-1]


def test_field_workflow_accepts_only_an_explicit_post_flash_vapour_schedule():
    schedule = SourceRateSchedule(
        time_s=(0.0, 0.5, 1.0), rate_kg_s=(0.01, 0.02, 0.0),
        source_id="plant-tag:direct-vapour-h2",
    )
    request = FieldSemiFVRequest(
        _scenario(), transport=_transport(), direct_vapour_schedule=schedule,
    )
    result = run_field_semi_fv_screening(request)

    assert result.completed
    assert result.transport is not None
    assert result.transport.diagnostics.source_mode == "scheduled"
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(0.015)
    assert any("was not re-flashed" in item for item in result.applicability.warnings)
    with pytest.raises(ValueError, match="time-aligned source uncertainty"):
        run_field_semi_fv_envelope(request)


def test_field_workflow_rejects_a_nonzero_direct_schedule_endpoint():
    with pytest.raises(ValueError, match="final endpoint rate must be zero"):
        FieldSemiFVRequest(
            _scenario(),
            transport=_transport(),
            direct_vapour_schedule=SourceRateSchedule(
                time_s=(0.0, 0.5, 1.0),
                rate_kg_s=(0.01, 0.02, 0.001),
                source_id="plant-tag:nonzero-endpoint",
            ),
        )


def test_field_envelope_propagates_fingerprinted_atmospheric_source_corners():
    source_id = "source-schedule-a"
    evidence = FieldSourceScheduleEvidence(
        dataset_id="schedule-a", path="schedule.csv", sha256="0" * 64,
        row_count=3, source_boundary_id=source_id, common_clock_id="clock-a",
        source_kind="post_flash_atmospheric_vapour",
    )
    nominal = SourceRateSchedule((0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id)
    lower = SourceRateSchedule((0.0, 0.5, 1.0), (0.005, 0.01, 0.0), source_id=source_id)
    upper = SourceRateSchedule((0.0, 0.5, 1.0), (0.02, 0.04, 0.0), source_id=source_id)
    atmospheric = FieldAtmosphericSourceSchedule(
        nominal, evidence, "time_s", "rate_kg_s", lower, upper,
        "rate_lower_kg_s", "rate_upper_kg_s",
    )
    envelope = run_field_semi_fv_envelope(
        FieldSemiFVRequest(_scenario(sensor=False), transport=_transport()),
        atmospheric_source_schedule=atmospheric,
        max_cases=8,
        table_nodes=17,
    )

    assert len(envelope.cases) == 6  # two surface corners x three source-rate corners
    assert {case.corner["source_schedule"] for case in envelope.cases} == {
        "lower", "nominal", "upper",
    }
    assert envelope.atmospheric_source_schedule is atmospheric
    report = field_semi_fv_envelope_report(envelope)
    assert report["atmospheric_source_schedule"]["evidence"]["dataset_id"] == "schedule-a"


def test_field_envelope_accepts_a_declared_atmospheric_source_boundary():
    source_id = "declared-source-schedule"
    atmospheric = FieldAtmosphericSourceSchedule(
        SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id,
        ),
        FieldSourceScheduleEvidence(
            dataset_id="declared-schedule-a",
            path="schedule.csv",
            sha256="3" * 64,
            row_count=3,
            source_boundary_id=source_id,
            common_clock_id="clock-declared-a",
            source_kind="declared_atmospheric_vapour",
        ),
        "time_s",
        "rate_kg_s",
    )

    envelope = run_field_semi_fv_envelope(
        FieldSemiFVRequest(_scenario(sensor=False), transport=_transport()),
        atmospheric_source_schedule=atmospheric,
        max_cases=2,
        table_nodes=17,
    )

    assert envelope.atmospheric_source_schedule is atmospheric
    assert all(
        case.result.request.direct_vapour_schedule is not None
        for case in envelope.cases
    )


def test_field_atmospheric_schedule_preserves_source_geometry_corners():
    base = _scenario(sensor=False)
    scenario = replace(
        base,
        source=replace(
            base.source,
            location_m=(0.0, 0.0, 0.5),
            location_uncertainty_m=(
                BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-schedule"),
                BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-schedule"),
                BoundedValue(0.5, 0.5, 0.5, "m", "layout-review-schedule"),
            ),
        ),
    )
    request = FieldSemiFVRequest(scenario, transport=_transport())
    source_id = "schedule-geometry-a"
    atmospheric = FieldAtmosphericSourceSchedule(
        SourceRateSchedule((0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id),
        FieldSourceScheduleEvidence(
            dataset_id="schedule-geometry-a", path="schedule.csv", sha256="1" * 64,
            row_count=3, source_boundary_id=source_id, common_clock_id="clock-a",
            source_kind="post_flash_atmospheric_vapour",
        ),
        "time_s", "rate_kg_s",
    )

    envelope = run_field_semi_fv_envelope(
        request, atmospheric_source_schedule=atmospheric, max_cases=8, table_nodes=17,
    )

    assert len(envelope.cases) == 4  # two source-location x two surface corners
    assert {dict(case.values)["source_location_x_m"] for case in envelope.cases} == {
        -0.1, 0.1,
    }
    assert all(case.result.request.scenario.source.location_uncertainty_m is None for case in envelope.cases)


def test_field_atmospheric_schedule_rejects_unrepresented_duration_uncertainty():
    base = _scenario(sensor=False)
    source = replace(
        base.source,
        duration_uncertainty=BoundedValue(1.0, 0.8, 1.2, "s", "timing-review-schedule"),
    )
    request = FieldSemiFVRequest(
        replace(base, source=source), transport=replace(_transport(), duration_s=1.2),
    )
    source_id = "schedule-duration-a"
    atmospheric = FieldAtmosphericSourceSchedule(
        SourceRateSchedule((0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id),
        FieldSourceScheduleEvidence(
            dataset_id="schedule-duration-a", path="schedule.csv", sha256="2" * 64,
            row_count=3, source_boundary_id=source_id, common_clock_id="clock-a",
            source_kind="post_flash_atmospheric_vapour",
        ),
        "time_s", "rate_kg_s",
    )

    with pytest.raises(ValueError, match="cannot represent non-exact release-duration"):
        run_field_semi_fv_envelope(
            request, atmospheric_source_schedule=atmospheric, max_cases=8, table_nodes=17,
        )


def test_field_envelope_does_not_place_pool_or_droplet_schedule_at_the_nozzle():
    source_id = "source-schedule-a"
    evidence = FieldSourceScheduleEvidence(
        dataset_id="schedule-pool-a", path="schedule.csv", sha256="0" * 64,
        row_count=3, source_boundary_id=source_id, common_clock_id="clock-a",
        source_kind="pool_vapour",
    )
    atmospheric = FieldAtmosphericSourceSchedule(
        SourceRateSchedule((0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id),
        evidence, "time_s", "rate_kg_s",
    )

    with pytest.raises(ValueError, match="pool_vapour and droplet_evaporation"):
        run_field_semi_fv_envelope(
            FieldSemiFVRequest(_scenario(sensor=False), transport=_transport()),
            atmospheric_source_schedule=atmospheric,
        )


def test_field_workflow_rejects_an_anonymous_post_flash_schedule():
    schedule = SourceRateSchedule(
        time_s=(0.0, 0.5, 1.0), rate_kg_s=(0.01, 0.02, 0.0),
    )
    with pytest.raises(ValueError, match="explicit source_id provenance"):
        FieldSemiFVRequest(
            _scenario(), transport=_transport(), direct_vapour_schedule=schedule,
        )


def test_field_workflow_propagates_source_location_corners_into_transport_cases():
    base = _scenario()
    scenario = replace(
        base,
        source=replace(
            base.source,
            location_m=(0.0, 0.0, 0.5),
            location_uncertainty_m=(
                BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-01"),
                BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-01"),
                BoundedValue(0.5, 0.5, 0.5, "m", "elevation-review-01"),
            ),
        ),
        surface=SurfaceBoundary(heat_transfer_w_m2_k=BoundedValue(10.0)),
    )
    request = FieldSemiFVRequest(scenario, transport=_transport())

    envelope = run_field_semi_fv_envelope(request, max_cases=2)

    assert len(envelope.cases) == 2
    assert {dict(case.values)["source_location_x_m"] for case in envelope.cases} == {
        -0.1, 0.1,
    }
    assert all(case.result.completed for case in envelope.cases)


def test_field_workflow_adds_an_evidenced_internal_vapour_source_conservatively():
    source = FieldDistributedVapourSource(
        label="pool-vapour-01",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 0.2, 1.0), (0.0, 0.02, 0.0),
            source_id="pool-ledger-01",
        ),
        evidence_id="pool-ledger-audit-01",
        source_kind="pool_vapour",
    )
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(),
            distributed_vapour_sources=(source,),
        )
    )

    assert result.completed
    assert result.transport is not None
    assert result.source_preparation is not None
    assert result.source_preparation.flash_result is not None
    assert result.transport.diagnostics.source_mode == "primary_plus_distributed"
    assert result.transport.diagnostics.distributed_source_count == 1
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(
        result.source_preparation.flash_result.flash.vapour_mass_flow + 0.016
    )
    assert dict(result.transport.diagnostics.source_mass_injected_kg) == {
        "primary": pytest.approx(
            result.source_preparation.flash_result.flash.vapour_mass_flow
        ),
        "distributed:pool-vapour-01": pytest.approx(0.016),
    }
    assert any("pool-vapour-01" in item for item in result.applicability.warnings)


def test_field_distributed_source_rejects_an_unknown_source_kind():
    with pytest.raises(ValueError, match="source_kind must be one of"):
        FieldDistributedVapourSource(
            label="unknown-kind",
            position_m=(1.0, 0.0, 0.1),
            vertical_sigma_m=0.15,
            schedule=SourceRateSchedule(
                (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="unknown-kind-ledger",
            ),
            evidence_id="unknown-kind-evidence",
            source_kind="unvalidated_custom_source",
        )


def test_field_request_rejects_mutable_corner_collections():
    source = FieldDistributedVapourSource(
        label="mutable-source",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="mutable-ledger",
        ),
        evidence_id="mutable-evidence",
    )
    with pytest.raises(TypeError, match="distributed_vapour_sources must be a tuple"):
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(),
            distributed_vapour_sources=[source],
        )
    with pytest.raises(TypeError, match="obstacles must be a tuple"):
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(), obstacles=[],
        )
    with pytest.raises(TypeError, match="direct_vapour_warnings must be a tuple"):
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(),
            direct_vapour_warnings=["manual warning"],
        )
    with pytest.raises(TypeError, match="position_m must be a tuple"):
        FieldDistributedVapourSource(
            label="mutable-position",
            position_m=[1.0, 0.0, 0.1],
            vertical_sigma_m=0.15,
            schedule=source.schedule,
            evidence_id="mutable-position-evidence",
        )
    with pytest.raises(TypeError, match="warnings must be a tuple"):
        FieldDistributedVapourSource(
            label="mutable-warnings",
            position_m=(1.0, 0.0, 0.1),
            vertical_sigma_m=0.15,
            schedule=source.schedule,
            evidence_id="mutable-warnings-evidence",
            warnings=["manual warning"],
        )


def test_field_workflow_blocks_an_off_plane_distributed_source():
    source = FieldDistributedVapourSource(
        label="off-plane-vapour",
        position_m=(1.0, 1.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 1.0), (0.01, 0.0), source_id="off-plane-ledger",
        ),
        evidence_id="off-plane-audit-01",
    )
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(),
            lateral_source_tolerance_m=0.01,
            distributed_vapour_sources=(source,),
        )
    )

    assert not result.completed
    assert result.source_preparation is None
    assert "outside the local wind plane" in result.applicability.reasons[0]


def test_field_envelope_propagates_distributed_source_geometry_and_width_bounds():
    source = FieldDistributedVapourSource(
        label="distributed-release",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="distributed-ledger-A",
        ),
        evidence_id="distributed-evidence-A",
        position_uncertainty_m=(
            BoundedValue(1.0, 0.9, 1.1, "m", "distributed-layout-A"),
            BoundedValue(0.0, 0.0, 0.0, "m", "distributed-layout-A"),
            BoundedValue(0.1, 0.1, 0.1, "m", "distributed-layout-A"),
        ),
        vertical_sigma_uncertainty=BoundedValue(
            0.15, 0.1, 0.2, "m", "distributed-width-A",
        ),
    )
    request = FieldSemiFVRequest(
        _scenario(sensor=False), transport=_transport(),
        distributed_vapour_sources=(source,),
    )

    envelope = run_field_semi_fv_envelope(request, max_cases=8, table_nodes=17)
    assert len(envelope.cases) == 8  # two distributed-width x two surface corners
    assert {
        dict(case.values)["distributed_source.distributed-release.position_x_m"]
        for case in envelope.cases
    } == {0.9, 1.1}
    assert {
        dict(case.values)["distributed_source.distributed-release.vertical_sigma_m"]
        for case in envelope.cases
    } == {0.1, 0.2}
    assert {
        field_source.position_m[0]
        for case in envelope.cases
        for field_source in case.result.request.distributed_vapour_sources
    } == {0.9, 1.1}
    assert {
        field_source.vertical_sigma_m
        for case in envelope.cases
        for field_source in case.result.request.distributed_vapour_sources
    } == {0.1, 0.2}
    report = field_semi_fv_envelope_report(envelope)
    declared = report["declared_distributed_vapour_sources"][0]
    assert declared["position_uncertainty_m"][0]["source"] == "distributed-layout-A"
    assert declared["vertical_sigma_uncertainty"]["source"] == "distributed-width-A"


def test_field_envelope_propagates_distributed_source_rate_schedule_corners():
    nominal = SourceRateSchedule(
        (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="distributed-ledger-rate-A",
    )
    source = FieldDistributedVapourSource(
        label="rate-bounded-release",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=nominal,
        lower_schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.01, 0.005, 0.0), source_id="distributed-ledger-rate-A",
        ),
        upper_schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.03, 0.015, 0.0), source_id="distributed-ledger-rate-A",
        ),
        evidence_id="distributed-rate-evidence-A",
    )
    request = FieldSemiFVRequest(
        _scenario(sensor=False), transport=_transport(),
        distributed_vapour_sources=(source,),
    )

    envelope = run_field_semi_fv_envelope(request, max_cases=8, table_nodes=17)

    assert len(envelope.cases) == 6  # three rate corners x two surface corners
    labels = {
        dict(case.values)["distributed_source.rate-bounded-release.schedule"]
        for case in envelope.cases
    }
    assert labels == {"lower", "nominal", "upper"}
    masses = {
        dict(case.values)["distributed_source.rate-bounded-release.schedule"]: {
            field_source.schedule.released_mass_kg
            for field_source in case.result.request.distributed_vapour_sources
        }
        for case in envelope.cases
    }
    assert next(iter(masses["lower"])) == pytest.approx(0.0075)
    assert next(iter(masses["nominal"])) == pytest.approx(0.015)
    assert next(iter(masses["upper"])) == pytest.approx(0.0225)
    report = field_semi_fv_envelope_report(envelope)
    uncertainty = report["declared_distributed_vapour_sources"][0]["schedule_uncertainty"]
    assert uncertainty["lower_rate_kg_s"] == [0.01, 0.005, 0.0]


def test_field_envelope_rejects_a_distributed_schedule_longer_than_a_duration_corner():
    base = _scenario(sensor=False)
    scenario = replace(
        base,
        source=replace(
            base.source,
            duration_uncertainty=BoundedValue(
                1.0, 0.8, 1.2, "s", "timing-review-A",
            ),
        ),
    )
    distributed = FieldDistributedVapourSource(
        label="late-pool-vapour",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="late-pool-ledger-A",
        ),
        evidence_id="late-pool-evidence-A",
        source_kind="pool_vapour",
    )
    request = FieldSemiFVRequest(
        scenario, transport=_transport(), distributed_vapour_sources=(distributed,),
    )
    with pytest.raises(ValueError, match="duration uncertainty corner"):
        run_field_semi_fv_envelope(request, max_cases=2)


def test_field_workflow_transports_a_finite_release_after_shutoff_without_extra_mass():
    base = FieldSemiFVRequest(_scenario(), transport=_transport())
    result = run_field_semi_fv_screening(
        replace(base, post_release_duration_s=1.0)
    )

    assert result.completed
    assert result.transport is not None
    assert result.source_preparation is not None
    assert result.source_preparation.flash_result is not None
    assert result.transport.diagnostics.source_mode == "scheduled"
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(
        result.source_preparation.flash_result.flash.vapour_mass_flow
    )
    assert result.sensor_trace is not None
    assert result.sensor_trace.time_s[-1] == pytest.approx(2.0)
    assert any("after source shutoff" in item for item in result.applicability.warnings)


def test_field_workflow_runs_untruncated_declared_uncertainty_envelope():
    request = FieldSemiFVRequest(_scenario(sensor=False), transport=_transport())
    envelope = run_field_semi_fv_envelope(request, max_cases=4, table_nodes=49)
    assert len(envelope.cases) == 2  # default surface heat-transfer lower/upper bounds
    assert envelope.property_table_used
    assert envelope.completed_case_count == 2
    assert envelope.blocked_case_count == 0
    assert {case.corner["surface_heat_transfer_w_m2_k"] for case in envelope.cases} == {5.0, 20.0}


def test_field_workflow_uses_only_declared_stability_scalar_mixing_value():
    base = _scenario(sensor=False)
    scenario = replace(
        base,
        weather=replace(base.weather, stability="stable"),
    )
    closure = StabilityScalarMixingClosure(
        {"stable": 0.125}, evidence_id="site-mast-closure-2026-10-05",
    )

    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            scenario, transport=_transport(), stability_mixing_closure=closure,
        )
    )

    assert result.completed
    assert result.transport is not None
    assert result.transport.diagnostics.scalar_diffusivity_m2_s == pytest.approx(0.125)
    assert any("caller-declared scalar mixing closure" in item for item in result.applicability.warnings)


def test_field_workflow_rejects_mixing_closure_without_the_scenario_class():
    base = _scenario(sensor=False)
    scenario = replace(base, weather=replace(base.weather, stability="stable"))

    with pytest.raises(ValueError, match="no diffusivity for scenario stability"):
        FieldSemiFVRequest(
            scenario,
            stability_mixing_closure=StabilityScalarMixingClosure(
                {"neutral": 0.25}, evidence_id="site-mast-closure-2026-10-05",
            ),
        )


def test_field_workflow_blocks_variable_measured_wind_before_flash_calculation():
    wind_history = FieldWindHistory(
        WindHistory(
            time_s=(0.0, 0.5, 1.0),
            speed_m_s=(2.0, 2.0, 2.0),
            direction_from_deg=(270.0, 300.0, 310.0),
        ),
        source_id="met-mast-01",
        evidence_id="met-mast-cal-1",
        maximum_direction_span_deg=20.0,
    )
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(), wind_history=wind_history,
        )
    )

    assert not result.completed
    assert result.source_preparation is None
    assert result.wind_history_assessment is not None
    assert any("wind_direction_span_exceeds_steady_limit" in item for item in result.applicability.reasons)


def test_field_workflow_retains_nominal_wind_after_a_passed_measured_wind_gate():
    wind_history = FieldWindHistory(
        WindHistory(
            time_s=(0.0, 0.5, 1.0),
            speed_m_s=(2.0, 2.1, 2.0),
            direction_from_deg=(270.0, 271.0, 270.0),
        ),
        source_id="met-mast-01",
        evidence_id="met-mast-cal-1",
    )
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            _scenario(sensor=False), transport=_transport(), wind_history=wind_history,
        )
    )

    assert result.completed
    assert result.wind_history_assessment is not None
    assert result.wind_history_assessment.applicable
    assert any("nominal wind consistent" in item for item in result.applicability.warnings)


def test_field_workflow_blocks_weather_that_disagrees_with_the_measured_wind_record():
    wind_history = FieldWindHistory(
        WindHistory(
            time_s=(0.0, 0.5, 1.0),
            speed_m_s=(2.0, 2.0, 2.0),
            direction_from_deg=(270.0, 270.0, 270.0),
        ),
        source_id="met-mast-01",
        evidence_id="met-mast-cal-1",
        maximum_nominal_speed_relative_difference=0.1,
    )
    base = _scenario(sensor=False)
    scenario = replace(base, weather=replace(base.weather, speed_m_s=BoundedValue(3.0)))
    result = run_field_semi_fv_screening(
        FieldSemiFVRequest(
            scenario, transport=_transport(), wind_history=wind_history,
        )
    )

    assert not result.completed
    assert result.source_preparation is None
    assert any("wind_speed_nominal_difference" in item for item in result.applicability.reasons)
