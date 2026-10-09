import math
from dataclasses import replace
from types import SimpleNamespace

import pytest

from degali.addons.field_contracts import (
    BoundedValue,
    FieldApplicability,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
)
from degali.addons.field_phase_routing import FieldPhaseRoutingConfig, FieldPhaseRoutingResult
from degali.addons.field_pool_launch import (
    FieldPoolVapourSchedule,
    PoolVapourLaunchBoundary,
    distributed_source_from_pool_vapour_schedule,
    pool_vapour_schedule_from_phase_routing,
    request_with_pool_vapour_schedule,
)
from degali.addons.field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening
from degali.addons.semi_fv_obstacle import SemiFVConfig, SourceRateSchedule


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
        surface=SurfaceBoundary(
            heat_transfer_w_m2_k=BoundedValue(10.0),
            surface_temperature_k=BoundedValue(293.15),
            substrate="concrete", evidence_id="surface-log-01",
        ),
        sensor=SensorModel((1.6, 0.0, 0.5)),
        temporal_mode="transient",
    )


def _phase_result():
    pool = SimpleNamespace(steps=(
        SimpleNamespace(elapsed_s=0.5, evaporated_mass_kg=0.1, wet_area_m2=0.2),
        SimpleNamespace(elapsed_s=1.0, evaporated_mass_kg=0.3, wet_area_m2=0.4),
    ))
    coupling = SimpleNamespace(pool=pool, pool_vapour_mass_kg=0.3)
    phase = SimpleNamespace(
        pool_coupling=coupling,
        droplets=SimpleNamespace(impact_centroid_m=(1.0, 0.0, 0.0)),
    )
    coupled = SimpleNamespace(conservative=True, phase_routing=phase)
    return FieldPhaseRoutingResult(
        _scenario(),
        FieldPhaseRoutingConfig(
            post_release_duration_s=1.0, puff_duration_s=1.0, pool_area_m2=0.5,
            pool_time_step_s=0.1, evaporation_coefficient_m2_s=1.0e-7,
        ),
        FieldApplicability("conditional"), None, 100.0, coupled,
    )


def test_pool_vapour_schedule_requires_typed_zero_terminated_source():
    with pytest.raises(ValueError, match="final endpoint rate"):
        FieldPoolVapourSchedule(
            schedule=SourceRateSchedule((0.0, 1.0), (0.2, 0.1), source_id="pool-source"),
            position_m=(1.0, 0.0, 0.0),
            maximum_wet_area_m2=1.0,
            evaporated_mass_kg=0.2,
            warnings=(),
        )
    with pytest.raises(TypeError, match="SourceRateSchedule"):
        FieldPoolVapourSchedule(
            schedule="not-a-schedule",
            position_m=(1.0, 0.0, 0.0),
            maximum_wet_area_m2=1.0,
            evaporated_mass_kg=0.2,
            warnings=(),
        )


def test_pool_vapour_schedule_preserves_dynamic_pool_mass_and_timing():
    schedule = pool_vapour_schedule_from_phase_routing(_phase_result())

    assert schedule.schedule.time_s == pytest.approx((0.0, 0.5, 1.0))
    assert schedule.schedule.rate_kg_s == pytest.approx((0.2, 0.4, 0.0))
    assert schedule.schedule.released_mass_kg == pytest.approx(0.3)
    assert schedule.position_m == (1.0, 0.0, 0.0)
    assert schedule.maximum_wet_area_m2 == pytest.approx(0.4)


def test_pool_vapour_launch_requires_evidence_and_stays_a_separate_scalar_screen():
    schedule = pool_vapour_schedule_from_phase_routing(_phase_result())
    with pytest.raises(ValueError, match="evidence_id"):
        PoolVapourLaunchBoundary()
    request = FieldSemiFVRequest(
        _scenario(),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    pool_request = request_with_pool_vapour_schedule(
        request, schedule, PoolVapourLaunchBoundary(evidence_id="pool-launch-closure-01"),
    )
    result = run_field_semi_fv_screening(pool_request)

    assert pool_request.scenario.source.location_m == (1.0, 0.0, 0.0)
    assert pool_request.direct_vapour_effective_area_m2 == pytest.approx(0.4)
    assert result.completed
    assert result.transport is not None
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(0.3)
    assert any("vertical momentum" in warning for warning in result.applicability.warnings)


def test_pool_vapour_schedule_resolves_an_original_release_duration_bound():
    schedule = pool_vapour_schedule_from_phase_routing(_phase_result())
    source = _scenario().source
    bounded_scenario = FieldScenario(
        source=replace(
            source,
            duration_uncertainty=BoundedValue(
                1.0, 0.8, 1.2, "s", "timing-review-01",
            ),
        ),
        weather=_scenario().weather,
        surface=_scenario().surface,
        sensor=_scenario().sensor,
        temporal_mode="transient",
    )
    request = FieldSemiFVRequest(
        bounded_scenario,
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )

    pool_request = request_with_pool_vapour_schedule(
        request, schedule, PoolVapourLaunchBoundary(evidence_id="pool-launch-closure-01"),
    )

    assert pool_request.scenario.source.duration_s == pytest.approx(1.0)
    assert pool_request.scenario.source.duration_uncertainty is None


def test_pool_vapour_can_be_an_internal_source_without_replacing_direct_flash():
    schedule = pool_vapour_schedule_from_phase_routing(_phase_result())
    pool_source = distributed_source_from_pool_vapour_schedule(
        schedule,
        PoolVapourLaunchBoundary(evidence_id="pool-launch-closure-01"),
        vertical_sigma_m=0.2,
        label="pool-internal-01",
    )
    request = FieldSemiFVRequest(
        _scenario(),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        distributed_vapour_sources=(pool_source,),
    )
    result = run_field_semi_fv_screening(request)

    assert pool_source.position_m == (1.0, 0.0, 0.0)
    assert pool_source.source_kind == "pool_vapour"
    assert pool_source.schedule.released_mass_kg == pytest.approx(0.3)
    assert result.completed
    assert result.source_preparation is not None
    assert result.source_preparation.flash_result is not None
    assert result.transport is not None
    assert result.transport.diagnostics.source_mode == "primary_plus_distributed"
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(
        result.source_preparation.flash_result.flash.vapour_mass_flow + 0.3
    )
