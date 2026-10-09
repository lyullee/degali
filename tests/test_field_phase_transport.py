import json
import math
from types import SimpleNamespace

import pytest
from dataclasses import replace

from degali.addons.field_contracts import (
    BoundedValue,
    FieldApplicability,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
)
from degali.addons.field_phase_routing import (
    FieldPhaseRoutingConfig,
    FieldPhaseRoutingEnvelope,
    FieldPhaseRoutingEnvelopeCase,
    FieldPhaseRoutingResult,
)
from degali.addons.field_decision import FieldOperationalScreeningDecision
from degali.addons.field_phase_transport import (
    FIELD_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA,
    FieldPhaseRoutingTransportEnvelopeCase,
    field_phase_routing_transport_envelope_report,
    run_field_phase_routing_transport_envelope,
)
from degali.addons.field_operational_phase_transport import (
    FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA,
    FieldOperationalPhaseRoutingTransportEnvelopeCase,
    field_operational_phase_routing_transport_envelope_report,
    run_field_operational_phase_routing_transport_envelope,
)
from degali.addons.field_pool_launch import PoolVapourLaunchBoundary
from degali.addons.field_geometry import FieldObstacleGeometryUncertainty
from degali.addons.field_workflow import FieldSemiFVRequest
from degali.addons.semi_fv_obstacle import SemiFVConfig
from degali.addons.site_geometry import AxisAlignedCuboid
import degali.addons.field_phase_transport as field_phase_transport
import degali.addons.field_operational_phase_transport as field_operational_phase_transport


def test_phase_transport_withheld_corners_keep_specific_gate_codes():
    phase_failed = field_operational_phase_transport._withheld_decision(
        "phase-routing ledger did not complete; field transport withheld",
    )
    pool_failed = field_operational_phase_transport._withheld_decision(
        "pool-vapour field transport withheld: source is off-plane",
    )

    assert phase_failed.gate_codes == (
        "transport_incomplete", "phase_routing_incomplete",
    )
    assert pool_failed.gate_codes == (
        "transport_incomplete", "pool_transport_withheld",
    )


def _scenario():
    return FieldScenario(
        source=ReleaseSource(
            fluid="lh2", location_m=(0.0, 0.0, 0.5),
            upstream_pressure=BoundedValue(0.4e6),
            upstream_temperature=BoundedValue(26.084),
            mass_flow_kg_s=BoundedValue(0.265),
            opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8),
            discharge_coefficient=BoundedValue(0.8), liquid_fraction=BoundedValue(0.922),
            flash_model="homogeneous_equilibrium", duration_s=1.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        surface=SurfaceBoundary(
            heat_transfer_w_m2_k=BoundedValue(10.0),
            surface_temperature_k=BoundedValue(293.15),
            substrate="concrete", evidence_id="surface-log-01",
        ),
        sensor=SensorModel((1.6, 0.0, 0.5)), temporal_mode="transient",
    )


def _config():
    return FieldPhaseRoutingConfig(
        post_release_duration_s=1.0, puff_duration_s=1.0, pool_area_m2=0.5,
        pool_time_step_s=0.1, evaporation_coefficient_m2_s=1.0e-7,
    )


def _phase_result(scenario, *, timed=True, pool_duration_s=1.0):
    midpoint = pool_duration_s / 2.0
    steps = (
        SimpleNamespace(elapsed_s=midpoint, evaporated_mass_kg=0.1, wet_area_m2=0.2),
        SimpleNamespace(elapsed_s=pool_duration_s, evaporated_mass_kg=0.3, wet_area_m2=0.4),
    ) if timed else ()
    pool = SimpleNamespace(steps=steps)
    phase = SimpleNamespace(
        pool_coupling=SimpleNamespace(pool=pool, pool_vapour_mass_kg=0.3),
        droplets=SimpleNamespace(impact_centroid_m=(1.0, 0.0, 0.0)),
    )
    coupled = SimpleNamespace(conservative=True, phase_routing=phase)
    return FieldPhaseRoutingResult(
        scenario, _config(), FieldApplicability("conditional"),
        None, 100.0, coupled,
    )


def _request():
    return FieldSemiFVRequest(
        _scenario(),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )


def _envelope(request, phase_result):
    return FieldPhaseRoutingEnvelope(
        request.scenario, _config(),
        (FieldPhaseRoutingEnvelopeCase(
            (("surface_heat_transfer_w_m2_k", 10.0),), phase_result,
        ),),
        property_table_used=False,
    )


def test_phase_routing_corner_contracts_reject_malformed_or_duplicate_records():
    scenario = _scenario()
    phase = _phase_result(scenario)
    with pytest.raises(TypeError, match="corner values cannot be boolean"):
        FieldPhaseRoutingEnvelopeCase((("corner", True),), phase)
    with pytest.raises(TypeError, match="envelope result"):
        FieldPhaseRoutingEnvelopeCase((("corner", 1.0),), object())
    duplicate = FieldPhaseRoutingEnvelopeCase((("corner", 1.0),), phase)
    with pytest.raises(ValueError, match="corner selections must be unique"):
        FieldPhaseRoutingEnvelope(
            scenario, _config(), (duplicate, duplicate), property_table_used=False,
        )
    with pytest.raises(TypeError, match="property_table_used"):
        FieldPhaseRoutingEnvelope(
            scenario, _config(), (duplicate,), property_table_used=1,
        )


def test_phase_transport_case_contract_rejects_malformed_records():
    scenario = _scenario()
    phase = _phase_result(scenario)
    with pytest.raises(TypeError, match="selection values cannot be boolean"):
        FieldPhaseRoutingTransportEnvelopeCase((("corner", True),), phase, None, None)
    with pytest.raises(TypeError, match="phase_routing"):
        FieldPhaseRoutingTransportEnvelopeCase((("corner", 1.0),), object(), None, None)
    with pytest.raises(TypeError, match="pool_schedule"):
        FieldPhaseRoutingTransportEnvelopeCase((("corner", 1.0),), phase, object(), None)
    with pytest.raises(TypeError, match="screening"):
        FieldPhaseRoutingTransportEnvelopeCase((("corner", 1.0),), phase, None, object())
    with pytest.raises(TypeError, match="warnings"):
        FieldPhaseRoutingTransportEnvelopeCase(
            (("corner", 1.0),), phase, None, None, warnings=("ok", 1),
        )


def test_phase_transport_adds_conservative_pool_vapour_without_replacing_direct_flash(monkeypatch):
    request = _request()
    phase = _phase_result(request.scenario)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )

    result = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2,
    )
    report = field_phase_routing_transport_envelope_report(result)
    case = result.cases[0]

    assert case.completed
    assert case.pool_schedule is not None
    assert case.pool_schedule.evaporated_mass_kg == 0.3
    assert case.screening is not None and case.screening.transport is not None
    assert case.screening.transport.diagnostics.source_mode == "primary_plus_distributed"
    assert report["schema"] == FIELD_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA
    assert report["completed_case_count"] == 1
    assert report["cases"][0]["pool_vapour"]["evaporated_mass_kg"] == 0.3
    assert "in-flight vapour" in report["scope"]
    json.dumps(report, allow_nan=False)


def test_phase_transport_propagates_bounded_obstacle_geometry_corners(monkeypatch):
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
    phase = _phase_result(request.scenario)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )

    result = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2, max_cases=2,
    )

    assert len(result.cases) == 2
    assert {
        dict(case.selection)["obstacle.shed-a.x_min_m"]
        for case in result.cases
    } == {1.8, 2.2}


def test_phase_transport_withholds_a_nonzero_pool_without_timed_conservative_source(monkeypatch):
    request = _request()
    phase = _phase_result(request.scenario, timed=False)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )

    result = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2,
    )

    assert result.completed_case_count == 0
    assert result.withheld_case_count == 1
    assert result.cases[0].screening is None
    assert "time-resolved dynamic-pool ledger" in result.cases[0].warnings[-1]


def test_phase_transport_rejects_a_duration_that_would_clip_pool_evaporation():
    request = replace(
        _request(),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    with pytest.raises(ValueError, match="must cover source.duration_s"):
        run_field_phase_routing_transport_envelope(
            request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
            pool_vertical_sigma_m=0.2,
        )


def test_phase_transport_aligns_pool_ledger_with_phase_post_release_window(monkeypatch):
    request = _request()
    phase = _phase_result(request.scenario, pool_duration_s=2.0)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )

    result = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2,
    )
    case = result.cases[0]

    assert case.completed
    assert case.pool_schedule is not None
    assert case.pool_schedule.schedule.duration_s == 2.0
    assert case.screening is not None
    assert case.screening.request.post_release_duration_s == 1.0


def test_phase_transport_propagates_declared_pool_vertical_width_bounds(monkeypatch):
    request = _request()
    phase = _phase_result(request.scenario)
    captured = []
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_semi_fv_screening",
        lambda field_request: captured.append(field_request) or SimpleNamespace(completed=True),
    )

    result = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2,
        pool_vertical_sigma_uncertainty=BoundedValue(
            0.2, 0.1, 0.3, unit="m", source="pool-width-review-A",
        ),
        max_cases=2,
    )

    assert len(result.cases) == 2
    assert {case.selection[-1] for case in result.cases} == {
        ("pool_vertical_sigma_m", 0.1),
        ("pool_vertical_sigma_m", 0.3),
    }
    assert {
        source.vertical_sigma_m
        for field_request in captured
        for source in field_request.distributed_vapour_sources
    } == {0.1, 0.3}
    assert result.pool_vertical_sigma_uncertainty is not None
    assert result.pool_vertical_sigma_uncertainty.source == "pool-width-review-A"


def test_phase_transport_propagates_bounded_ambient_boundaries(monkeypatch):
    request = replace(
        _request(),
        ambient_temperature_uncertainty_k=BoundedValue(
            295.0, 290.0, 300.0, unit="K", source="met-boundary-A",
        ),
    )
    phase = _phase_result(request.scenario)
    phase_calls = []
    field_requests = []

    def phase_envelope(*_args, **kwargs):
        phase_calls.append(kwargs["ambient_temperature_k"])
        return _envelope(request, phase)

    monkeypatch.setattr(field_phase_transport, "run_field_phase_routing_envelope", phase_envelope)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_semi_fv_screening",
        lambda field_request: field_requests.append(field_request)
        or SimpleNamespace(completed=True, request=field_request),
    )

    result = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2, max_cases=2,
    )

    assert len(result.cases) == 2
    assert {dict(case.selection)["ambient_temperature_k"] for case in result.cases} == {
        290.0, 300.0,
    }
    assert set(phase_calls) == {290.0, 300.0}
    assert {
        field_request.ambient_temperature_k for field_request in field_requests
    } == {290.0, 300.0}
    assert all(
        field_request.ambient_temperature_uncertainty_k is None
        for field_request in field_requests
    )


def test_phase_transport_rejects_a_pool_width_bound_with_the_wrong_nominal():
    with pytest.raises(ValueError, match="nominal must match"):
        run_field_phase_routing_transport_envelope(
            _request(), _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
            pool_vertical_sigma_m=0.2,
            pool_vertical_sigma_uncertainty=BoundedValue(
                0.25, 0.1, 0.3, unit="m", source="pool-width-review-A",
            ),
        )


def test_operational_phase_transport_refines_every_pool_and_sensor_corner(monkeypatch):
    scenario = replace(
        _scenario(),
        sensor=SensorModel(
            (1.6, 0.0, 0.5), gain=BoundedValue(1.0, 0.9, 1.1),
        ),
    )
    request = replace(_request(), scenario=scenario)
    phase = _phase_result(scenario)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )
    physical = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2,
    )
    monkeypatch.setattr(
        field_operational_phase_transport,
        "run_field_phase_routing_transport_envelope",
        lambda *_args, **_kwargs: physical,
    )

    result = run_field_operational_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2, max_cases=2,
        refinement_factors=(1, 2), relative_tolerance=1.0, allow_conditional=True,
    )
    report = field_operational_phase_routing_transport_envelope_report(result)

    assert len(result.cases) == 2
    assert all(case.refinement is not None for case in result.cases)
    assert all(case.decision.screening_allowed for case in result.cases)
    assert all(
        case.screening is not None
        and case.screening.request.phase_routing_uncertainty_resolved
        for case in result.cases
    )
    assert result.operational_decision.status in {"screening_allowed", "conditional_allowed"}
    expected_gate_codes = {
        code for case in result.cases for code in case.decision.gate_codes
    }
    assert set(result.operational_decision.gate_codes) == expected_gate_codes
    if result.operational_decision.status == "conditional_allowed":
        assert "conditional_review_required" in result.operational_decision.gate_codes
    assert report["schema"] == FIELD_OPERATIONAL_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA
    assert report["operational_screening"]["uncertainty_resolved"]
    assert report["deterministic_sensor_envelope"]["status"] == "complete"
    assert report["deterministic_sensor_envelope"]["sensor_count"] == 1
    assert all(
        item["field_result"]["transport_input"]["phase_routing_uncertainty_resolved"]
        for item in report["cases"]
    )
    with pytest.raises(TypeError, match="sensor_selections must be a tuple"):
        replace(result, sensor_selections=list(result.sensor_selections))
    with pytest.raises(TypeError, match="cases must be a tuple"):
        replace(result, cases=list(result.cases))


def test_operational_phase_aggregate_keeps_uncertainty_flag_when_only_refinement_is_missing(monkeypatch):
    request = _request()
    phase = _phase_result(request.scenario)
    monkeypatch.setattr(
        field_phase_transport,
        "run_field_phase_routing_envelope",
        lambda *_args, **_kwargs: _envelope(request, phase),
    )
    physical = run_field_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2,
    )
    monkeypatch.setattr(
        field_operational_phase_transport,
        "run_field_phase_routing_transport_envelope",
        lambda *_args, **_kwargs: physical,
    )

    result = run_field_operational_phase_routing_transport_envelope(
        request, _config(), PoolVapourLaunchBoundary(evidence_id="pool-launch-A"),
        pool_vertical_sigma_m=0.2, max_cases=1, include_refinement=False,
    )

    assert result.operational_decision.status == "withheld"
    assert result.operational_decision.uncertainty_resolved
    assert "refinement_missing" in result.operational_decision.gate_codes


def test_operational_phase_case_contract_rejects_malformed_corner_records():
    decision = FieldOperationalScreeningDecision(
        "withheld", False, False, ("missing",), (), True, False, True, False,
    )
    with pytest.raises(TypeError, match="cannot be boolean"):
        FieldOperationalPhaseRoutingTransportEnvelopeCase(
            (("phase", True),), (), None, None, decision,
        )
    with pytest.raises(TypeError, match="finite numbers"):
        FieldOperationalPhaseRoutingTransportEnvelopeCase(
            (("phase", "nominal"),), (("sensor", "nominal"),), None, None, decision,
        )
    with pytest.raises(TypeError, match="refinement must"):
        FieldOperationalPhaseRoutingTransportEnvelopeCase(
            (("phase", "nominal"),), (), None,
            object(), decision,
        )
    with pytest.raises(ValueError, match="cannot allow"):
        FieldOperationalPhaseRoutingTransportEnvelopeCase(
            (("phase", "nominal"),), (), None, None,
            FieldOperationalScreeningDecision(
                "screening_allowed", True, False, ("allowed",), (),
                True, True, True, True,
            ),
        )


def test_operational_phase_aggregate_preserves_uncertainty_for_conditional_holdback():
    conditional = FieldOperationalScreeningDecision(
        "conditional_allowed", True, False,
        ("conditional obstacle/launch scope requires review",),
        ("review conditional obstacle/launch scope before use",),
        True, True, True, True,
    )
    result = field_operational_phase_transport._aggregate(
        (SimpleNamespace(label="conditional", decision=conditional),),
        allow_conditional=False,
    )

    assert result.status == "withheld"
    assert result.uncertainty_resolved
    assert result.gate_codes == ("conditional_review_required",)
