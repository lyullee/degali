import json

import pytest

from degali.addons.droplet_rainout import (
    DropletClass,
    DropletTransportInput,
    transport_droplet_population,
)
from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_droplet_launch import (
    DropletVapourLaunchBoundary,
    distributed_sources_from_droplet_evaporation,
    run_field_droplet_handoff_refinement_study,
)
from degali.addons.field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening
from degali.addons.field_report import field_droplet_handoff_refinement_report
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _boundary():
    return DropletTransportInput(
        liquid_mass_flow_kg_s=0.2,
        classes=(
            DropletClass(diameter_m=1.0e-4, mass_fraction=0.25),
            DropletClass(diameter_m=1.0e-3, mass_fraction=0.75),
        ),
        liquid_density_kg_m3=70.0,
        initial_velocity_m_s=(0.0, 0.0, 0.0),
        release_position_m=(0.0, 0.0, 0.1),
        air_density_kg_m3=1.2,
        air_kinematic_viscosity_m2_s=1.5e-5,
        evaporation_coefficient_m2_s=1.0e-7,
        schmidt_number=0.7,
        wind_velocity_m_s=(5.0, 0.0, 0.0),
    )


def test_droplet_trajectory_handoff_closes_vapour_mass_without_nozzle_collapse():
    droplets = transport_droplet_population(
        _boundary(), max_time_s=2.0, trajectory_segment_duration_s=0.05,
    )
    sources = distributed_sources_from_droplet_evaporation(
        droplets,
        release_duration_s=3.0,
        launch=DropletVapourLaunchBoundary(
            vertical_sigma_m=0.2, evidence_id="droplet-flight-closure-01",
        ),
    )

    assert sources
    assert sum(item.schedule.released_mass_kg for item in sources) == pytest.approx(
        droplets.airborne_vapour_mass_flow_kg_s * 3.0,
        abs=1.0e-12,
    )
    assert all(item.source_kind == "in_flight_droplet_evaporation" for item in sources)
    assert all(item.position_m[2] >= 0.0 for item in sources)
    assert any(item.position_m[0] > 0.0 for item in sources)
    assert all(item.schedule.time_s[0] == 0.0 for item in sources)


def test_droplet_handoff_refuses_endpoint_only_and_source_count_overflow():
    endpoint_only = transport_droplet_population(_boundary(), max_time_s=2.0)
    launch = DropletVapourLaunchBoundary(
        vertical_sigma_m=0.2, evidence_id="droplet-flight-closure-01",
    )
    with pytest.raises(ValueError, match="trajectory_segment_duration_s"):
        distributed_sources_from_droplet_evaporation(
            endpoint_only, release_duration_s=3.0, launch=launch,
        )


def test_trajectory_sources_enter_the_field_screen_as_separate_delayed_sources():
    droplets = transport_droplet_population(
        _boundary(), max_time_s=2.0, trajectory_segment_duration_s=0.05,
    )
    sources = distributed_sources_from_droplet_evaporation(
        droplets,
        release_duration_s=3.0,
        launch=DropletVapourLaunchBoundary(
            vertical_sigma_m=0.2, evidence_id="droplet-flight-closure-01",
        ),
    )
    scenario = FieldScenario(
        source=ReleaseSource(
            fluid="lh2", location_m=(0.0, 0.0, 0.1),
            upstream_pressure=BoundedValue(0.4e6),
            upstream_temperature=BoundedValue(26.084),
            mass_flow_kg_s=BoundedValue(0.265),
            opening_area_m2=BoundedValue(3.534e-4),
            discharge_coefficient=BoundedValue(0.8),
            liquid_fraction=BoundedValue(0.922),
            flash_model="homogeneous_equilibrium", duration_s=3.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        temporal_mode="transient",
    )
    result = run_field_semi_fv_screening(FieldSemiFVRequest(
        scenario,
        transport=SemiFVConfig(
            length_m=10.0, height_m=3.0, nx=40, nz=20,
            duration_s=5.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        post_release_duration_s=2.0,
        distributed_vapour_sources=sources,
    ))

    assert result.completed
    assert result.transport is not None
    assert result.source_preparation is not None
    assert result.source_preparation.flash_result is not None
    assert result.transport.diagnostics.distributed_source_count == len(sources)
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(
        3.0 * (
            result.source_preparation.flash_result.flash.vapour_mass_flow
            + droplets.airborne_vapour_mass_flow_kg_s
        ),
        abs=1.0e-10,
    )


def test_droplet_handoff_refinement_holds_field_inputs_fixed_and_reports_sensor_change():
    scenario = FieldScenario(
        source=ReleaseSource(
            fluid="lh2", location_m=(0.0, 0.0, 0.1),
            upstream_pressure=BoundedValue(0.4e6),
            upstream_temperature=BoundedValue(26.084),
            mass_flow_kg_s=BoundedValue(0.265),
            opening_area_m2=BoundedValue(3.534e-4),
            discharge_coefficient=BoundedValue(0.8),
            liquid_fraction=BoundedValue(0.922),
            flash_model="homogeneous_equilibrium", duration_s=3.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        sensor=SensorModel((1.0, 0.0, 0.1)),
        temporal_mode="transient",
    )
    request = FieldSemiFVRequest(
        scenario,
        transport=SemiFVConfig(
            length_m=10.0, height_m=3.0, nx=40, nz=20,
            duration_s=5.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        post_release_duration_s=2.0,
    )
    study = run_field_droplet_handoff_refinement_study(
        request,
        _boundary(),
        release_duration_s=3.0,
        launch=DropletVapourLaunchBoundary(
            vertical_sigma_m=0.2, evidence_id="droplet-flight-closure-01",
        ),
        trajectory_segment_durations_s=(0.2, 0.1, 0.05),
        maximum_droplet_time_s=2.0,
        relative_tolerance=0.05,
    )

    assert len(study.cases) == 3
    assert study.reference_trajectory_segment_duration_s == pytest.approx(0.05)
    assert study.cases[0].source_count <= study.cases[-1].source_count
    assert study.cases[0].in_flight_vapour_mass_kg == pytest.approx(
        study.cases[-1].in_flight_vapour_mass_kg, abs=1.0e-12,
    )
    assert study.peak_relative_changes[-1][1] == pytest.approx(0.0)
    assert study.dose_relative_changes[-1][1] == pytest.approx(0.0)
    assert study.converged
    report = field_droplet_handoff_refinement_report(study)
    assert report["converged"] is True
    assert report["reference_trajectory_segment_duration_s"] == pytest.approx(0.05)
    assert len(report["cases"]) == 3
    assert report["cases"][-1]["field_screening"]["transport_result"] is not None
    json.dumps(report, allow_nan=False)

    recorded = transport_droplet_population(
        _boundary(), max_time_s=2.0, trajectory_segment_duration_s=0.05,
    )
    with pytest.raises(ValueError, match="maximum_sources"):
        distributed_sources_from_droplet_evaporation(
            recorded,
            release_duration_s=3.0,
            launch=DropletVapourLaunchBoundary(
                vertical_sigma_m=0.2, evidence_id="droplet-flight-closure-01",
                maximum_sources=1,
            ),
        )
