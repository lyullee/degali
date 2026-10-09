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
from degali.addons.field_decision import FieldOperationalScreeningDecision
from degali.addons.field_operational_sensor_array import (
    FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA,
    field_operational_sensor_array_envelope_report,
    run_field_operational_sensor_array_envelope,
)
from degali.addons.field_workflow import (
    FieldSensorDeployment,
    FieldSemiFVRequest,
)
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _request() -> FieldSemiFVRequest:
    source = ReleaseSource(
        fluid="lh2", location_m=(0.0, 0.0, 0.5),
        upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.265, unit="kg/s"),
        opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"),
        discharge_coefficient=BoundedValue(0.8), liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium", duration_s=1.0,
    )
    return FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.8, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        sensor_deployments=(FieldSensorDeployment(
            "detector-uncertain",
            SensorModel(
                (0.8, 0.0, 0.5),
                response_time_s=BoundedValue(0.3, 0.1, 0.5),
            ),
        ),),
    )


def test_operational_sensor_array_refines_and_aggregates_every_calibration_corner():
    result = run_field_operational_sensor_array_envelope(
        _request(), max_cases=4, refinement_factors=(1, 2),
        relative_tolerance=1.0, allow_conditional=True,
    )
    report = field_operational_sensor_array_envelope_report(result)

    assert len(result.cases) == 2
    assert all(case.refinement is not None for case in result.cases)
    assert result.operational_decision.status == "conditional_allowed"
    assert result.operational_decision.uncertainty_resolved
    assert "physical_applicability_conditional" in result.operational_decision.gate_codes
    assert "conditional_review_required" in result.operational_decision.gate_codes
    assert report["schema"] == FIELD_OPERATIONAL_SENSOR_ARRAY_ENVELOPE_SCHEMA
    assert report["deterministic_sensor_envelope"]["status"] == "complete"
    assert report["deterministic_sensor_envelope"]["sensor_count"] == 2
    assert all(
        item["available_case_count"] == 2
        for item in report["deterministic_sensor_envelope"]["sensors"]
    )
    assert len(report["cases"]) == 2
    json.dumps(report, allow_nan=False)
    with pytest.raises(ValueError, match="missing refinement"):
        replace(
            result.cases[0],
            refinement=None,
            decision=FieldOperationalScreeningDecision(
                "screening_allowed", True, False, (), (), True, True, True, True,
            ),
        )


def test_operational_sensor_array_withholds_without_refinement_but_keeps_calibration_resolved():
    result = run_field_operational_sensor_array_envelope(
        _request(), max_cases=4, include_refinement=False,
    )

    assert result.operational_decision.status == "withheld"
    assert result.operational_decision.uncertainty_resolved
    assert "refinement_missing" in result.operational_decision.gate_codes
    assert all(case.refinement is None for case in result.cases)


def test_operational_sensor_array_does_not_waive_unpropagated_source_uncertainty():
    request = _request()
    source = replace(
        request.scenario.source,
        mass_flow_kg_s=BoundedValue(0.265, 0.24, 0.29, unit="kg/s"),
    )
    request = replace(
        request,
        scenario=replace(request.scenario, source=source),
    )
    result = run_field_operational_sensor_array_envelope(
        request, max_cases=4, include_refinement=False,
    )

    assert result.operational_decision.status == "withheld"
    assert not result.operational_decision.uncertainty_resolved
    assert "uncertainty" in " ".join(result.operational_decision.reasons)
