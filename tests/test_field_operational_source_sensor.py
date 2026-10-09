import json
import math

from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_operational_source_sensor import (
    FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA,
    field_operational_source_sensor_envelope_report,
    run_field_operational_source_sensor_envelope,
)
from degali.addons.field_workflow import FieldSensorDeployment, FieldSemiFVRequest
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _request() -> FieldSemiFVRequest:
    source = ReleaseSource(
        fluid="lh2", location_m=(0.0, 0.0, 0.5),
        upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.265, 0.24, 0.29, unit="kg/s"),
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


def test_joint_source_sensor_envelope_crosses_physical_and_calibration_corners():
    result = run_field_operational_source_sensor_envelope(
        _request(), max_cases=8, refinement_factors=(1, 2),
        relative_tolerance=1.0, allow_conditional=True,
    )
    report = field_operational_source_sensor_envelope_report(result)

    assert len(result.physical_envelope.cases) == 4
    assert len(result.sensor_selections) == 2
    assert len(result.cases) == 8
    assert result.operational_decision.status == "conditional_allowed"
    assert result.operational_decision.uncertainty_resolved
    assert "physical_applicability_conditional" in result.operational_decision.gate_codes
    assert "conditional_review_required" in result.operational_decision.gate_codes
    assert report["schema"] == FIELD_OPERATIONAL_SOURCE_SENSOR_ENVELOPE_SCHEMA
    assert report["deterministic_sensor_envelope"]["status"] == "complete"
    assert report["deterministic_sensor_envelope"]["sensor_count"] == 2
    assert len(report["cases"]) == 8
    json.dumps(report, allow_nan=False)
