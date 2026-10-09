import math
from dataclasses import replace

import numpy as np
import pytest

from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_comparison import (
    FieldComparisonBasis,
    field_model_sensor_set_from_superposition,
)
from degali.addons.field_observation import apply_sensor_model
from degali.addons.field_report import (
    FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA,
    field_sensor_superposition_report,
)
from degali.addons.field_superposition import (
    FieldSensorBranch,
    superpose_field_sensor_branches,
)
from degali.addons.field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _scenario():
    return FieldScenario(
        source=ReleaseSource(
            fluid="lh2", location_m=(0.0, 0.0, 0.5),
            upstream_pressure=BoundedValue(0.4e6),
            upstream_temperature=BoundedValue(26.084),
            mass_flow_kg_s=BoundedValue(0.265),
            opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8),
            discharge_coefficient=BoundedValue(0.8),
            liquid_fraction=BoundedValue(0.922),
            flash_model="homogeneous_equilibrium", duration_s=1.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        sensor=SensorModel(
            (1.6, 0.0, 0.5), response_time_s=BoundedValue(0.2),
            gain=BoundedValue(1.1), bias_mole_fraction=BoundedValue(0.001),
            averaging_time_s=0.1,
        ),
        temporal_mode="transient",
    )


def _request(scenario):
    return FieldSemiFVRequest(
        scenario,
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        post_release_duration_s=0.5,
    )


def test_common_sensor_field_branches_are_summed_before_sensor_calibration():
    scenario = _scenario()
    direct = run_field_semi_fv_screening(_request(scenario))
    offset_source = replace(scenario.source, location_m=(0.2, 0.0, 0.5))
    pool_like = run_field_semi_fv_screening(_request(replace(scenario, source=offset_source)))
    combined = superpose_field_sensor_branches((
        FieldSensorBranch("direct_flash_vapour", direct),
        FieldSensorBranch("pool_vapour", pool_like),
    ))

    assert combined.completed
    assert combined.applicability.status == "conditional"
    assert direct.transport is not None
    assert pool_like.transport is not None
    expected_concentration = (
        np.interp(
            combined.time_s,
            direct.transport.receptor_traces[0].time_s,
            direct.transport.receptor_traces[0].concentration_kg_m3,
        )
        + np.interp(
            combined.time_s,
            pool_like.transport.receptor_traces[0].time_s,
            pool_like.transport.receptor_traces[0].concentration_kg_m3,
        )
    )
    assert combined.concentration_kg_m3 == pytest.approx(expected_concentration)
    expected_sensor = apply_sensor_model(
        combined.time_s,
        expected_concentration,
        scenario.sensor,
        ambient_air_density_kg_m3=1.2,
    )
    assert combined.sensor_trace.indicated_mole_fraction == pytest.approx(
        expected_sensor.indicated_mole_fraction
    )
    report = field_sensor_superposition_report(combined)
    assert report["schema"] == FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA
    assert [item["label"] for item in report["branches"]] == [
        "direct_flash_vapour", "pool_vapour",
    ]


def test_superposition_exports_conditional_model_comparison_provenance():
    scenario = _scenario()
    direct = run_field_semi_fv_screening(_request(scenario))
    offset_source = replace(scenario.source, location_m=(0.2, 0.0, 0.5))
    second = run_field_semi_fv_screening(_request(replace(scenario, source=offset_source)))
    combined = superpose_field_sensor_branches((
        FieldSensorBranch("direct", direct),
        FieldSensorBranch("second", second),
    ))

    model = field_model_sensor_set_from_superposition(
        combined,
        model_id="DEGALI-superposed",
        basis=FieldComparisonBasis(
            "phase-superposition", "weather-a", "sensor-a", "peak-indicated",
        ),
        prediction_operator="peak_indicated",
    )
    assert model.predictions[0].sensor_id == "field_sensor"
    assert model.predictions[0].mole_fraction == pytest.approx(
        np.max(combined.sensor_trace.indicated_mole_fraction)
    )
    assert model.execution_provenance is not None
    assert model.execution_provenance.source_applicability.status == "conditional"


def test_sensor_superposition_marks_unresolved_branch_bounds_incomplete():
    scenario = _scenario()
    uncertain_source = replace(
        scenario.source,
        upstream_pressure=BoundedValue(0.4e6, 0.35e6, 0.45e6, "Pa", "pressure-review"),
    )
    first = run_field_semi_fv_screening(
        _request(replace(scenario, source=uncertain_source))
    )
    second = run_field_semi_fv_screening(
        _request(replace(scenario, source=replace(uncertain_source, location_m=(0.2, 0.0, 0.5))))
    )
    combined = superpose_field_sensor_branches((
        FieldSensorBranch("first", first),
        FieldSensorBranch("second", second),
    ))

    assert not combined.applicability.uncertainty_complete
    assert any("upstream_pressure" in warning for warning in combined.warnings)


def test_sensor_superposition_refuses_missing_common_trace_coverage():
    scenario = _scenario()
    direct = run_field_semi_fv_screening(_request(scenario))
    short = run_field_semi_fv_screening(
        replace(_request(scenario), post_release_duration_s=0.0)
    )

    with pytest.raises(ValueError, match="trace end time"):
        superpose_field_sensor_branches((
            FieldSensorBranch("direct", direct),
            FieldSensorBranch("short", short),
        ))


def test_sensor_superposition_refuses_different_internal_time_grids():
    scenario = _scenario()
    direct = run_field_semi_fv_screening(_request(scenario))
    different_grid_request = _request(scenario)
    different_grid = run_field_semi_fv_screening(
        replace(
            different_grid_request,
            transport=replace(different_grid_request.transport, time_step_s=0.01),
        )
    )

    with pytest.raises(ValueError, match="identical sensor time grid"):
        superpose_field_sensor_branches((
            FieldSensorBranch("direct", direct),
            FieldSensorBranch("different-grid", different_grid),
        ))


def test_sensor_superposition_refuses_steady_transient_temporal_mismatch():
    scenario = _scenario()
    transient = run_field_semi_fv_screening(_request(scenario))
    steady_request = _request(scenario)
    steady_request = replace(
        steady_request,
        scenario=replace(scenario, temporal_mode="steady"),
        post_release_duration_s=0.0,
        transport=replace(steady_request.transport, duration_s=1.5),
    )
    steady = run_field_semi_fv_screening(
        steady_request
    )

    with pytest.raises(ValueError, match="temporal mode"):
        superpose_field_sensor_branches((
            FieldSensorBranch("transient", transient),
            FieldSensorBranch("steady", steady),
        ))


def test_sensor_trace_and_superposition_contracts_reject_malformed_arrays():
    scenario = _scenario()
    direct = run_field_semi_fv_screening(_request(scenario))
    offset_source = replace(scenario.source, location_m=(0.2, 0.0, 0.5))
    second = run_field_semi_fv_screening(_request(replace(scenario, source=offset_source)))
    combined = superpose_field_sensor_branches((
        FieldSensorBranch("direct", direct),
        FieldSensorBranch("second", second),
    ))

    with pytest.raises(ValueError, match="does not match the true sensor trace"):
        replace(combined, concentration_kg_m3=combined.concentration_kg_m3 * 2.0)
    assert direct.sensor_trace is not None
    with pytest.raises(TypeError, match="numpy array"):
        replace(direct.sensor_trace, time_s=list(direct.sensor_trace.time_s))
