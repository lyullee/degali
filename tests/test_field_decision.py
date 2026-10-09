import json
import math
from dataclasses import replace

import pytest

from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_decision import (
    FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA,
    evaluate_field_operational_screening,
    field_operational_screening_decision_record,
)
from degali.addons.field_distributed_source import FieldDistributedVapourSource
from degali.addons.field_geometry import FieldObstacleGeometryUncertainty
from degali.addons.field_workflow import (
    FieldSensorDeployment,
    FieldSemiFVRequest,
    run_field_semi_fv_refinement_study,
    run_field_semi_fv_screening,
)
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.addons.semi_fv_obstacle import SemiFVConfig, SourceRateSchedule


def _request(*, extra_sensor=None):
    source = ReleaseSource(
        fluid="lh2", location_m=(0.0, 0.0, 0.5),
        upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.265, unit="kg/s"),
        opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"),
        discharge_coefficient=BoundedValue(0.8), liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium", duration_s=1.0,
    )
    scenario = FieldScenario(
        source=source,
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        sensor=SensorModel((0.8, 0.0, 0.5)), temporal_mode="transient",
    )
    return FieldSemiFVRequest(
        scenario,
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        sensor_deployments=() if extra_sensor is None else (
            FieldSensorDeployment("detector-02", extra_sensor),
        ),
    )


def test_operational_gate_withholds_a_raw_screen_without_refinement_evidence():
    decision = evaluate_field_operational_screening(run_field_semi_fv_screening(_request()))

    assert decision.status == "withheld"
    assert not decision.screening_allowed
    assert not decision.design_basis_allowed
    assert any("refinement evidence" in reason for reason in decision.reasons)
    assert "refinement_missing" in decision.gate_codes
    with pytest.raises(ValueError, match="operational screening is withheld"):
        decision.require_screening()


def test_operational_decision_contract_rejects_malformed_status_and_holdback():
    common = dict(
        screening_allowed=False,
        design_basis_allowed=False,
        reasons=("test holdback",),
        required_actions=(),
        refinement_required=True,
        refinement_available=False,
        uncertainty_resolution_required=True,
        uncertainty_resolved=False,
    )
    with pytest.raises(ValueError, match="status is unsupported"):
        from degali.addons.field_decision import FieldOperationalScreeningDecision
        FieldOperationalScreeningDecision("unknown", **common)
    with pytest.raises(ValueError, match="requires at least one reason"):
        from degali.addons.field_decision import FieldOperationalScreeningDecision
        FieldOperationalScreeningDecision("withheld", reasons=(), **{
            key: value for key, value in common.items() if key != "reasons"
        })


def test_allowed_operational_decision_cannot_bypass_required_gates():
    from degali.addons.field_decision import FieldOperationalScreeningDecision

    with pytest.raises(ValueError, match="required refinement evidence"):
        FieldOperationalScreeningDecision(
            "screening_allowed", True, False, (), (),
            refinement_required=True, refinement_available=False,
            uncertainty_resolution_required=False, uncertainty_resolved=True,
        )
    with pytest.raises(ValueError, match="uncertainty remains unresolved"):
        FieldOperationalScreeningDecision(
            "conditional_allowed", True, False, (), (),
            refinement_required=False, refinement_available=False,
            uncertainty_resolution_required=True, uncertainty_resolved=False,
        )


def test_operational_decision_direct_boundary_requires_immutable_explanations():
    from degali.addons.field_decision import FieldOperationalScreeningDecision

    common = dict(
        status="screening_allowed", screening_allowed=True,
        design_basis_allowed=False, reasons=("accepted basis",),
        required_actions=(), refinement_required=False,
        refinement_available=False, uncertainty_resolution_required=False,
        uncertainty_resolved=True,
    )
    with pytest.raises(TypeError, match="reasons must be a tuple"):
        FieldOperationalScreeningDecision(**{**common, "reasons": ["accepted basis"]})
    with pytest.raises(TypeError, match="actions must be a tuple"):
        FieldOperationalScreeningDecision(**{**common, "required_actions": ["review"]})
    with pytest.raises(ValueError, match="conditional_allowed screening requires"):
        FieldOperationalScreeningDecision(
            status="conditional_allowed", screening_allowed=True,
            design_basis_allowed=False, reasons=(), required_actions=(),
            refinement_required=False, refinement_available=False,
            uncertainty_resolution_required=False, uncertainty_resolved=True,
        )
    with pytest.raises(ValueError, match="recorded review action"):
        FieldOperationalScreeningDecision(
            status="conditional_allowed", screening_allowed=True,
            design_basis_allowed=False,
            reasons=("conditional scope",), required_actions=(),
            refinement_required=False, refinement_available=False,
            uncertainty_resolution_required=False, uncertainty_resolved=True,
        )
    with pytest.raises(ValueError, match="unsupported code"):
        FieldOperationalScreeningDecision(
            **common, gate_codes=("not-a-real-gate",),
        )


def test_operational_gate_requires_explicit_conditional_review_even_after_refinement():
    refined = run_field_semi_fv_refinement_study(
        _request(), refinement_factors=(1, 2), relative_tolerance=1.0,
    )
    withheld = evaluate_field_operational_screening(refined)
    allowed = evaluate_field_operational_screening(refined, allow_conditional=True)
    record = field_operational_screening_decision_record(allowed)

    assert refined.completed
    assert withheld.status == "withheld"
    assert allowed.status == "conditional_allowed"
    assert allowed.screening_allowed
    assert not allowed.design_basis_allowed
    assert not allowed.approval_allowed
    assert record["schema"] == FIELD_OPERATIONAL_SCREENING_DECISION_SCHEMA
    assert record["refinement_available"] is True
    assert "physical_applicability_conditional" in record["gate_codes"]
    json.dumps(record, allow_nan=False)


def test_operational_gate_withholds_an_off_plane_declared_detector():
    result = run_field_semi_fv_screening(_request(extra_sensor=SensorModel((1.0, 1.0, 0.5))))
    decision = evaluate_field_operational_screening(
        result, require_refinement=False, allow_conditional=True,
    )

    assert result.completed
    assert decision.status == "withheld"
    assert any("detector-02" in reason for reason in decision.reasons)
    assert "sensor_trace_missing" in decision.gate_codes


def test_operational_gate_withholds_unpropagated_source_weather_or_sensor_bounds():
    request = _request()
    uncertain_source = replace(
        request.scenario.source,
        mass_flow_kg_s=BoundedValue(0.265, 0.24, 0.29, unit="kg/s"),
    )
    uncertain_weather = replace(
        request.scenario.weather,
        speed_m_s=BoundedValue(2.0, 1.8, 2.2, unit="m/s"),
    )
    uncertain_request = replace(
        request,
        scenario=replace(
            request.scenario, source=uncertain_source, weather=uncertain_weather,
        ),
    )
    refined = run_field_semi_fv_refinement_study(
        uncertain_request, refinement_factors=(1, 2), relative_tolerance=1.0,
    )
    decision = evaluate_field_operational_screening(refined, allow_conditional=True)

    assert refined.completed
    assert decision.status == "withheld"
    assert decision.uncertainty_resolution_required
    assert not decision.uncertainty_resolved
    assert "mass_flow_kg_s" in " ".join(decision.reasons)
    assert "wind_speed_m_s" in " ".join(decision.reasons)
    assert "uncertainty_unresolved" in decision.gate_codes


def test_direct_schedule_does_not_hide_source_geometry_or_timing_bounds():
    request = _request()
    source = replace(
        request.scenario.source,
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, unit="m", source="layout-direct-schedule"),
            BoundedValue(0.0, 0.0, 0.0, unit="m", source="layout-direct-schedule"),
            BoundedValue(0.5, 0.5, 0.5, unit="m", source="layout-direct-schedule"),
        ),
        direction_uncertainty_m=(
            BoundedValue(1.0, 1.0, 1.0, unit="1", source="direction-direct-schedule"),
            BoundedValue(0.0, -0.1, 0.1, unit="1", source="direction-direct-schedule"),
            BoundedValue(0.0, 0.0, 0.0, unit="1", source="direction-direct-schedule"),
        ),
    )
    request = replace(
        request,
        scenario=replace(request.scenario, source=source),
        direct_vapour_schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.01, 0.02, 0.0),
            source_id="history-direct-schedule",
        ),
    )

    result = run_field_semi_fv_screening(request)
    decision = evaluate_field_operational_screening(
        result, require_refinement=False, allow_conditional=True,
    )

    assert result.completed
    assert decision.status == "withheld"
    assert not decision.uncertainty_resolved
    reason_text = " ".join(decision.reasons)
    assert "source_location_x_m" in reason_text
    assert "source_direction_y" in reason_text


def test_operational_gate_withholds_unpropagated_distributed_source_geometry_bounds():
    source = FieldDistributedVapourSource(
        label="distributed-release",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="distributed-ledger-A",
        ),
        evidence_id="distributed-evidence-A",
        position_uncertainty_m=(
            BoundedValue(1.0, 0.9, 1.1, unit="m", source="distributed-layout-A"),
            BoundedValue(0.0, 0.0, 0.0, unit="m", source="distributed-layout-A"),
            BoundedValue(0.1, 0.1, 0.1, unit="m", source="distributed-layout-A"),
        ),
        vertical_sigma_uncertainty=BoundedValue(
            0.15, 0.1, 0.2, unit="m", source="distributed-width-A",
        ),
    )
    request = replace(_request(), distributed_vapour_sources=(source,))
    result = run_field_semi_fv_screening(request)
    decision = evaluate_field_operational_screening(
        result, require_refinement=False, allow_conditional=True,
    )

    assert result.completed
    assert decision.status == "withheld"
    assert not decision.uncertainty_resolved
    reason_text = " ".join(decision.reasons)
    assert "distributed_source.distributed-release.position_x_m" in reason_text
    assert "distributed_source.distributed-release.vertical_sigma_m" in reason_text


def test_operational_gate_withholds_unpropagated_obstacle_geometry_bounds():
    obstacle = AxisAlignedCuboid(2.0, 3.0, -1.0, 1.0, 0.0, 1.0, "shed-a")
    uncertainty = FieldObstacleGeometryUncertainty(
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
    request = replace(
        _request(), obstacle=obstacle,
        obstacle_geometry_uncertainty=(uncertainty,),
    )
    result = run_field_semi_fv_screening(request)
    decision = evaluate_field_operational_screening(
        result, require_refinement=False,
    )

    assert result.completed
    assert decision.status == "withheld"
    assert not decision.uncertainty_resolved
    assert "obstacle.shed-a.x_min_m" in " ".join(decision.reasons)


def test_operational_gate_withholds_a_declared_cuboid_outside_the_local_plane():
    request = replace(
        _request(),
        obstacles=(
            AxisAlignedCuboid(1.5, 2.0, 2.0, 3.0, 0.0, 1.0, "off-plane-skid"),
        ),
    )
    refined = run_field_semi_fv_refinement_study(
        request, refinement_factors=(1, 2), relative_tolerance=1.0,
    )
    decision = evaluate_field_operational_screening(refined, allow_conditional=True)

    assert refined.completed
    assert decision.status == "withheld"
    assert "off-plane-skid" in " ".join(decision.reasons)
    assert "three-dimensional obstacle/wake" in " ".join(decision.required_actions)


def test_operational_gate_withholds_when_declared_sensor_result_is_missing():
    screening = run_field_semi_fv_screening(_request())
    malformed = replace(screening, sensor_results=())
    decision = evaluate_field_operational_screening(
        malformed, require_refinement=False, allow_conditional=True,
    )

    assert screening.completed
    assert decision.status == "withheld"
    assert "field_sensor" in " ".join(decision.reasons)
    assert "three-dimensional model" in " ".join(decision.required_actions)


def test_operational_gate_does_not_drop_obstacles_on_projection_cardinality_mismatch():
    request = replace(
        _request(),
        obstacles=(AxisAlignedCuboid(1.5, 2.0, -0.2, 0.2, 0.0, 1.0, "skid-A"),),
    )
    screening = run_field_semi_fv_screening(request)
    malformed = replace(
        screening, obstacle_projection=None, obstacle_projections=(),
    )
    decision = evaluate_field_operational_screening(
        malformed, require_refinement=False, allow_conditional=True,
    )

    assert screening.completed
    assert decision.status == "withheld"
    assert "skid-A" in " ".join(decision.reasons)
