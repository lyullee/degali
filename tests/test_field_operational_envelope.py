import json
import math
from dataclasses import replace

import pytest

from degali.addons.field_contracts import BoundedValue, FieldScenario, ReleaseSource, SensorModel, SurfaceBoundary, WeatherState
from degali.addons.field_distributed_source import FieldDistributedVapourSource
from degali.addons.field_geometry import FieldObstacleGeometryUncertainty
from degali.addons.field_decision import (
    FieldOperationalScreeningDecision,
    evaluate_field_operational_screening,
)
from degali.addons.field_meteorology import FieldStabilityAlternatives, StabilityScalarMixingClosure
from degali.addons.field_source_io import FieldAtmosphericSourceSchedule, FieldSourceScheduleEvidence
from degali.addons.field_operational_envelope import (
    FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA,
    FieldOperationalUncertaintyEnvelopeCase,
    _aggregate_decision,
    deterministic_sensor_envelope_record,
    field_operational_uncertainty_envelope_report,
    run_field_operational_uncertainty_envelope,
)
from degali.addons.field_workflow import (
    FieldSemiFVRequest,
    run_field_semi_fv_envelope,
    run_field_semi_fv_refinement_study,
)
from degali.addons.semi_fv_obstacle import SemiFVConfig
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


def _request():
    source = ReleaseSource(
        fluid="lh2", location_m=(0.0, 0.0, 0.5),
        upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.265, 0.24, 0.29, unit="kg/s"),
        opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"),
        discharge_coefficient=BoundedValue(0.8),
        liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium", duration_s=1.0,
    )
    scenario = FieldScenario(
        source=source,
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        surface=SurfaceBoundary(
            heat_transfer_w_m2_k=BoundedValue(10.0, unit="W/m2/K"),
            surface_temperature_k=BoundedValue(293.15, unit="K"),
            substrate="concrete", evidence_id="surface-review-A",
        ),
        sensor=SensorModel((0.8, 0.0, 0.5)), temporal_mode="transient",
    )
    return FieldSemiFVRequest(
        scenario,
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )


def test_operational_uncertainty_envelope_refines_and_decides_every_corner():
    result = run_field_operational_uncertainty_envelope(
        _request(), max_cases=2, refinement_factors=(1, 2),
        relative_tolerance=1.0, allow_conditional=True,
    )
    report = field_operational_uncertainty_envelope_report(result)

    assert len(result.cases) == 2
    assert result.operational_decision.status == "conditional_allowed"
    assert result.operational_decision.uncertainty_resolved
    assert "physical_applicability_conditional" in result.operational_decision.gate_codes
    assert "conditional_review_required" in result.operational_decision.gate_codes
    assert all(case.refinement is not None and case.decision.screening_allowed for case in result.cases)
    assert report["schema"] == FIELD_OPERATIONAL_UNCERTAINTY_ENVELOPE_SCHEMA
    assert report["operational_screening"]["status"] == "conditional_allowed"
    assert len(report["cases"]) == 2
    sensor_envelope = report["deterministic_sensor_envelope"]
    assert sensor_envelope["status"] == "complete"
    assert sensor_envelope["case_count"] == 2
    assert sensor_envelope["sensor_count"] == 1
    sensor_record = sensor_envelope["sensors"][0]
    assert sensor_record["label"] == "field_sensor"
    assert sensor_record["available_case_count"] == 2
    assert sensor_record["withheld_case_count"] == 0
    assert sensor_record["bounds"]["peak_true_mole_fraction"]["minimum"] <= sensor_record["bounds"]["peak_true_mole_fraction"]["maximum"]
    json.dumps(report, allow_nan=False)


def test_operational_uncertainty_envelope_rejects_a_fabricated_aggregate_status():
    result = run_field_operational_uncertainty_envelope(
        _request(), max_cases=2, refinement_factors=(1, 2),
        relative_tolerance=1.0, allow_conditional=True,
    )
    fabricated = FieldOperationalScreeningDecision(
        "screening_allowed", True, False,
        ("fabricated aggregate",), ("review aggregate",),
        refinement_required=True, refinement_available=True,
        uncertainty_resolution_required=True, uncertainty_resolved=True,
    )
    with pytest.raises(ValueError, match="status does not match"):
        replace(result, operational_decision=fabricated)
    incomplete_screening = replace(
        result.cases[0].screening,
        transport=None, sensor_trace=None, sensor_results=(),
    )
    allowed_corner = FieldOperationalScreeningDecision(
        "screening_allowed", True, False, (), (),
        refinement_required=True, refinement_available=True,
        uncertainty_resolution_required=True, uncertainty_resolved=True,
    )
    with pytest.raises(ValueError, match="completed screening result"):
        replace(
            result.cases[0],
            screening=incomplete_screening,
            decision=allowed_corner,
        )
    with pytest.raises(ValueError, match="accepted physical applicability"):
        replace(result.cases[0], decision=allowed_corner)


def test_operational_sensor_envelope_withholds_missing_corner_trace():
    request = replace(
        _request(),
        scenario=replace(
            _request().scenario,
            sensor=SensorModel((0.8, 0.25, 0.5)),
        ),
    )
    result = run_field_operational_uncertainty_envelope(
        request, max_cases=2, include_refinement=False,
    )
    report = field_operational_uncertainty_envelope_report(result)

    sensor_envelope = report["deterministic_sensor_envelope"]
    assert sensor_envelope["status"] == "withheld"
    sensor_record = sensor_envelope["sensors"][0]
    assert sensor_record["available_case_count"] == 0
    assert sensor_record["withheld_case_count"] == 2
    assert sensor_record["bounds"] == {}
    assert sensor_record["withheld_reasons"]


def test_operational_sensor_envelope_keeps_legacy_primary_trace_fallback():
    result = run_field_operational_uncertainty_envelope(
        _request(), max_cases=2, include_refinement=False,
    )
    legacy_screening = replace(result.cases[0].screening, sensor_results=())
    legacy_case = replace(result.cases[0], screening=legacy_screening)
    record = deterministic_sensor_envelope_record((legacy_case,))

    assert record["status"] == "complete"
    assert record["sensors"][0]["label"] == "field_sensor"
    assert record["sensors"][0]["available_case_count"] == 1


def test_operational_sensor_envelope_keeps_a_missing_screening_corner_withheld():
    record = deterministic_sensor_envelope_record((object(),), request=_request())

    assert record["status"] == "withheld"
    sensor_record = record["sensors"][0]
    assert sensor_record["available_case_count"] == 0
    assert sensor_record["withheld_case_count"] == 1
    assert sensor_record["bounds"] == {}
    assert sensor_record["withheld_reasons"]


def test_operational_uncertainty_envelope_keeps_diagnostic_corners_withheld_without_refinement():
    result = run_field_operational_uncertainty_envelope(
        _request(), max_cases=2, include_refinement=False,
    )

    assert result.operational_decision.status == "withheld"
    assert "refinement_missing" in result.operational_decision.gate_codes
    assert result.operational_decision.uncertainty_resolved
    assert all(case.refinement is None for case in result.cases)
    assert all(not case.decision.screening_allowed for case in result.cases)


def test_operational_uncertainty_aggregate_preserves_an_unresolved_corner_flag():
    result = run_field_operational_uncertainty_envelope(
        _request(), max_cases=2, include_refinement=False,
    )
    unresolved = replace(
        result.cases[0].decision,
        uncertainty_resolved=False,
        reasons=("unresolved synthetic corner",),
    )
    aggregate = _aggregate_decision(
        (replace(result.cases[0], decision=unresolved),),
        allow_conditional=False,
    )

    assert aggregate.status == "withheld"
    assert not aggregate.uncertainty_resolved


def test_operational_uncertainty_envelope_refines_fingerprinted_source_corners():
    request = _request()
    request = replace(
        request,
        scenario=replace(
            request.scenario,
            surface=replace(
                request.scenario.surface,
                heat_transfer_w_m2_k=BoundedValue(
                    10.0, 5.0, 20.0, unit="W/m2/K", source="surface-review-A",
                ),
            ),
        ),
    )
    source_id = "source-schedule-operational"
    evidence = FieldSourceScheduleEvidence(
        dataset_id="operational-schedule-A",
        path="schedule.csv",
        sha256="a" * 64,
        row_count=3,
        source_boundary_id=source_id,
        common_clock_id="clock-operational-A",
        source_kind="post_flash_atmospheric_vapour",
    )
    atmospheric = FieldAtmosphericSourceSchedule(
        SourceRateSchedule((0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id),
        evidence,
        "time_s",
        "rate_kg_s",
        SourceRateSchedule((0.0, 0.5, 1.0), (0.005, 0.01, 0.0), source_id=source_id),
        SourceRateSchedule((0.0, 0.5, 1.0), (0.02, 0.04, 0.0), source_id=source_id),
        "rate_lower_kg_s",
        "rate_upper_kg_s",
    )

    result = run_field_operational_uncertainty_envelope(
        request,
        atmospheric_source_schedule=atmospheric,
        max_cases=8,
        refinement_factors=(1, 2),
        relative_tolerance=1.0,
        allow_conditional=True,
    )

    assert len(result.cases) == 6  # two surface corners x three source corners
    assert {dict(case.selection)["source_schedule"] for case in result.cases} == {
        "lower", "nominal", "upper",
    }
    assert result.envelope.atmospheric_source_schedule is atmospheric
    assert result.operational_decision.status == "conditional_allowed"
    report = field_operational_uncertainty_envelope_report(result)
    assert report["field_uncertainty_envelope"]["atmospheric_source_schedule"]["evidence"][
        "dataset_id"
    ] == "operational-schedule-A"


def test_operational_source_schedule_keeps_source_geometry_corners():
    request = _request()
    source = replace(
        request.scenario.source,
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, unit="m", source="layout-review-A"),
            BoundedValue(0.0, 0.0, 0.0, unit="m", source="layout-review-A"),
            BoundedValue(0.5, 0.5, 0.5, unit="m", source="layout-review-A"),
        ),
    )
    request = replace(request, scenario=replace(request.scenario, source=source))
    source_id = "source-schedule-geometry-operational"
    atmospheric = FieldAtmosphericSourceSchedule(
        SourceRateSchedule((0.0, 0.5, 1.0), (0.01, 0.02, 0.0), source_id=source_id),
        FieldSourceScheduleEvidence(
            dataset_id="operational-schedule-geometry-A",
            path="schedule.csv",
            sha256="b" * 64,
            row_count=3,
            source_boundary_id=source_id,
            common_clock_id="clock-operational-geometry-A",
            source_kind="post_flash_atmospheric_vapour",
        ),
        "time_s",
        "rate_kg_s",
    )

    result = run_field_operational_uncertainty_envelope(
        request,
        atmospheric_source_schedule=atmospheric,
        max_cases=8,
        include_refinement=False,
    )

    assert len(result.cases) == 2  # two source-location corners
    assert {dict(case.selection)["source_location_x_m"] for case in result.cases} == {
        -0.1, 0.1,
    }
    assert all(
        case.screening.request.scenario.source.location_uncertainty_m is None
        for case in result.cases
    )


def test_operational_envelope_resolves_distributed_source_geometry_corners():
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
    request = replace(
        _request(),
        scenario=replace(_request().scenario, sensor=None),
        distributed_vapour_sources=(source,),
    )

    result = run_field_operational_uncertainty_envelope(
        request, max_cases=8, table_nodes=17, include_refinement=False,
    )

    assert len(result.cases) == 8
    assert all(case.decision.uncertainty_resolved for case in result.cases)
    assert all(
        "distributed_source.distributed-release.position_x_m"
        not in " ".join(case.decision.reasons)
        for case in result.cases
    )
    assert result.operational_decision.status == "withheld"


def test_operational_envelope_resolves_distributed_source_rate_corners():
    base = _request()
    source = FieldDistributedVapourSource(
        label="rate-bounded-release",
        position_m=(1.0, 0.0, 0.1),
        vertical_sigma_m=0.15,
        schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.02, 0.01, 0.0), source_id="distributed-rate-A",
        ),
        lower_schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.01, 0.005, 0.0), source_id="distributed-rate-A",
        ),
        upper_schedule=SourceRateSchedule(
            (0.0, 0.5, 1.0), (0.03, 0.015, 0.0), source_id="distributed-rate-A",
        ),
        evidence_id="distributed-rate-evidence-A",
    )
    scenario = replace(
        base.scenario,
        source=replace(base.scenario.source, mass_flow_kg_s=BoundedValue(0.265, unit="kg/s")),
        sensor=None,
    )
    result = run_field_operational_uncertainty_envelope(
        replace(base, scenario=scenario, distributed_vapour_sources=(source,)),
        max_cases=4, table_nodes=17, include_refinement=False,
    )

    assert len(result.cases) == 3
    assert {
        dict(case.selection)["distributed_source.rate-bounded-release.schedule"]
        for case in result.cases
    } == {"lower", "nominal", "upper"}
    assert all(case.decision.uncertainty_resolved for case in result.cases)
    assert all(
        "atmospheric rate schedule bounds" not in " ".join(case.decision.reasons)
        for case in result.cases
    )
    assert result.operational_decision.status == "withheld"


def test_operational_envelope_crosses_independent_distributed_source_rate_corners_and_ledgers():
    base = _request()
    scenario = replace(
        base.scenario,
        source=replace(
            base.scenario.source,
            mass_flow_kg_s=BoundedValue(0.265, unit="kg/s"),
        ),
        sensor=None,
    )

    def source(label: str, scale: float) -> FieldDistributedVapourSource:
        source_id = f"{label}-ledger"
        return FieldDistributedVapourSource(
            label=label,
            position_m=(1.0, 0.0, 0.1),
            vertical_sigma_m=0.15,
            schedule=SourceRateSchedule(
                (0.0, 0.5, 1.0),
                (0.02 * scale, 0.01 * scale, 0.0),
                source_id=source_id,
            ),
            lower_schedule=SourceRateSchedule(
                (0.0, 0.5, 1.0),
                (0.01 * scale, 0.005 * scale, 0.0),
                source_id=source_id,
            ),
            upper_schedule=SourceRateSchedule(
                (0.0, 0.5, 1.0),
                (0.03 * scale, 0.015 * scale, 0.0),
                source_id=source_id,
            ),
            evidence_id=f"{label}-evidence",
        )

    result = run_field_operational_uncertainty_envelope(
        replace(
            base,
            scenario=scenario,
            distributed_vapour_sources=(source("pool-a", 1.0), source("pool-b", 2.0)),
        ),
        max_cases=9,
        table_nodes=17,
        include_refinement=False,
    )

    assert len(result.cases) == 9
    assert {
        (
            dict(case.selection)["distributed_source.pool-a.schedule"],
            dict(case.selection)["distributed_source.pool-b.schedule"],
        )
        for case in result.cases
    } == {
        (left, right)
        for left in ("lower", "nominal", "upper")
        for right in ("lower", "nominal", "upper")
    }
    for case in result.cases:
        selection = dict(case.selection)
        expected = {
            "lower": 0.0075,
            "nominal": 0.015,
            "upper": 0.0225,
        }
        expected_b = {
            "lower": 0.015,
            "nominal": 0.03,
            "upper": 0.045,
        }
        ledger = dict(case.screening.transport.diagnostics.source_mass_injected_kg)
        assert ledger["distributed:pool-a"] == pytest.approx(
            expected[selection["distributed_source.pool-a.schedule"]]
        )
        assert ledger["distributed:pool-b"] == pytest.approx(
            expected_b[selection["distributed_source.pool-b.schedule"]]
        )
        assert case.screening.transport.diagnostics.source_mass_ledger_residual_kg == pytest.approx(
            0.0, abs=1.0e-12,
        )
        assert case.decision.uncertainty_resolved
    assert result.operational_decision.status == "withheld"


def test_operational_envelope_resolves_obstacle_geometry_corners():
    request = _request()
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
        request, obstacle=obstacle,
        obstacle_geometry_uncertainty=(uncertainty,),
    )

    result = run_field_operational_uncertainty_envelope(
        request, max_cases=4, include_refinement=False,
    )

    assert len(result.cases) == 4  # two obstacle corners x two source-flow corners
    assert all(case.decision.uncertainty_resolved for case in result.cases)
    assert all(
        "obstacle.shed-a.x_min_m" not in " ".join(case.decision.reasons)
        for case in result.cases
    )
    assert result.operational_decision.status == "withheld"


def test_operational_uncertainty_envelope_propagates_release_duration_corners():
    request = _request()
    source = replace(
        request.scenario.source,
        mass_flow_kg_s=BoundedValue(0.265, unit="kg/s"),
        duration_uncertainty=BoundedValue(1.0, 0.8, 1.2, unit="s", source="timing-review-A"),
    )
    scenario = replace(request.scenario, source=source)
    request = replace(
        request,
        scenario=scenario,
        transport=replace(request.transport, duration_s=1.5),
    )

    result = run_field_semi_fv_envelope(request, max_cases=2)

    assert {dict(case.values)["duration_s"] for case in result.cases} == {0.8, 1.2}
    assert {case.result.request.scenario.source.duration_s for case in result.cases} == {
        0.8, 1.2,
    }
    assert all(case.result.request.scenario.source.duration_uncertainty is None for case in result.cases)


def test_stability_alternatives_require_declared_closure_and_are_propagated():
    nominal = _request()
    alternatives = FieldStabilityAlternatives(
        ("stable",), "met-stability-classification-A",
    )
    with_coverage = replace(
        nominal,
        stability_mixing_closure=StabilityScalarMixingClosure(
            {"neutral": 0.5, "stable": 0.2}, evidence_id="met-mixing-A",
        ),
        stability_alternatives=alternatives,
    )
    envelope = run_field_semi_fv_envelope(with_coverage, max_cases=4)

    assert len(envelope.cases) == 4
    assert {dict(case.values)["weather_stability"] for case in envelope.cases} == {
        "neutral", "stable",
    }
    refined_nominal = run_field_semi_fv_refinement_study(
        with_coverage, refinement_factors=(1, 2), relative_tolerance=1.0,
    )
    withheld = evaluate_field_operational_screening(
        refined_nominal, allow_conditional=True,
    )
    assert withheld.status == "withheld"
    assert "stability alternatives" in " ".join(withheld.reasons)

    operational = run_field_operational_uncertainty_envelope(
        with_coverage, max_cases=4, refinement_factors=(1, 2),
        relative_tolerance=1.0, allow_conditional=True,
    )
    assert operational.operational_decision.status == "conditional_allowed"
    assert all(case.decision.uncertainty_resolved for case in operational.cases)
    assert all(
        case.screening.request.stability_uncertainty_resolved
        for case in operational.cases
    )
    report = field_operational_uncertainty_envelope_report(operational)
    assert all(
        item["field_result"]["transport_input"]["stability_alternatives"]
        ["uncertainty_resolved_for_this_case"]
        for item in report["cases"]
    )


def test_operational_uncertainty_case_contract_rejects_malformed_selection():
    decision = FieldOperationalScreeningDecision(
        "withheld", False, False, ("missing",), (), True, False, True, False,
    )
    with pytest.raises(TypeError, match="cannot be boolean"):
        FieldOperationalUncertaintyEnvelopeCase(
            (("source", True),), None, None, decision,
        )
    with pytest.raises(ValueError, match="keys must be unique"):
        FieldOperationalUncertaintyEnvelopeCase(
            (("source", 1.0), ("source", 2.0)), None, None, decision,
        )
