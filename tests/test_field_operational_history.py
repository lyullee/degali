import json
import math
from dataclasses import replace
from types import SimpleNamespace

import pytest

from degali.addons.field_contracts import BoundedValue, FieldScenario, ReleaseSource, SensorModel, WeatherState
from degali.addons.field_history import (
    MeasuredHistoryQualityCriteria,
    MeasuredReleaseHistory,
    MeasuredTimeSeries,
    PressureDrivenMeasuredHistory,
)
from degali.addons.field_meteorology import FieldStabilityAlternatives, StabilityScalarMixingClosure
from degali.addons.field_operational_history import (
    FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA,
    FieldOperationalMeasuredHistoryEnvelopeCase,
    _aggregate,
    field_operational_joint_measured_history_envelope_report,
    run_field_operational_joint_measured_history_envelope,
    run_field_operational_joint_pressure_driven_history_envelope,
)
from degali.addons.field_decision import FieldOperationalScreeningDecision
from degali.addons.field_workflow import FieldSemiFVRequest
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _request():
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
            weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
            sensor=SensorModel((0.8, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )


def _history():
    return MeasuredReleaseHistory(
        MeasuredTimeSeries((0.0, 1.0), (0.4e6, 0.401e6), unit="Pa", source_id="PT-01"),
        MeasuredTimeSeries((0.0, 1.0), (26.084, 26.2), unit="K", source_id="TT-01"),
        MeasuredTimeSeries(
            (0.0, 1.0), (0.265, 0.0), lower=(0.24, 0.0), upper=(0.29, 0.0),
            unit="kg/s", source_id="FT-01",
        ),
        MeasuredTimeSeries((0.0, 1.0), (0.922, 0.922), unit="1", source_id="LT-01"),
    )


def _quality():
    return MeasuredHistoryQualityCriteria(
        maximum_sample_interval_s=1.0, maximum_response_time_s=1.0,
        maximum_absolute_time_offset_s=0.1, maximum_relative_half_width=0.2,
        evidence_id="historian-quality-procedure-A",
    )


def _pressure_history():
    return PressureDrivenMeasuredHistory(
        MeasuredTimeSeries(
            (0.0, 1.0), (0.4e6, 0.39e6), unit="Pa", source_id="PT-ORIFICE-01",
        ),
        MeasuredTimeSeries(
            (0.0, 1.0), (26.076, 25.949), unit="K", source_id="TT-ORIFICE-01",
        ),
        event_id="pressure-event-A", phase_evidence_id="phase-evidence-A",
    )


def test_operational_joint_measured_history_envelope_refines_every_source_corner():
    base = _request()
    source = replace(
        base.scenario.source,
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-measured"),
            BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-measured"),
            BoundedValue(0.5, 0.5, 0.5, "m", "layout-review-measured"),
        ),
    )
    request = replace(base, scenario=replace(base.scenario, source=source))
    result = run_field_operational_joint_measured_history_envelope(
        request, _history(), quality_criteria=_quality(),
        history_provenance=(("event_id", "vent-event-A"), ("source_sha256", "a" * 64)),
        max_cases=4, refinement_factors=(1, 2), relative_tolerance=1.0,
        allow_conditional=True,
    )
    report = field_operational_joint_measured_history_envelope_report(result)

    assert len(result.cases) == 4
    assert result.operational_decision.status == "conditional_allowed"
    assert result.operational_decision.uncertainty_resolved
    assert "physical_applicability_conditional" in result.operational_decision.gate_codes
    assert "conditional_review_required" in result.operational_decision.gate_codes
    assert all(case.refinement is not None and case.decision.screening_allowed for case in result.cases)
    assert report["schema"] == FIELD_OPERATIONAL_MEASURED_HISTORY_ENVELOPE_SCHEMA
    assert report["operational_screening"]["status"] == "conditional_allowed"
    assert report["deterministic_sensor_envelope"]["status"] == "complete"
    assert report["deterministic_sensor_envelope"]["sensors"][0]["available_case_count"] == 4
    assert {
        dict(case.field_selection)["source_location_x_m"]
        for case in result.cases
    } == {-0.1, 0.1}
    provenance = report["cases"][0]["field_result"]["transport_input"]["measured_history_provenance"]
    assert provenance["event_id"] == "vent-event-A"
    assert provenance["source_sha256"] == "a" * 64
    json.dumps(report, allow_nan=False)


def test_operational_joint_pressure_driven_history_recomputes_and_gates_every_corner():
    base = _request()
    area = base.scenario.source.opening_area_m2.nominal
    source = replace(
        base.scenario.source,
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        opening_area_m2=BoundedValue(area, 0.9 * area, 1.1 * area, "m2", "area-review"),
        discharge_coefficient=BoundedValue(0.8, 0.72, 0.88, "1", "cd-review"),
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-pressure"),
            BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-pressure"),
            BoundedValue(0.5, 0.5, 0.5, "m", "layout-review-pressure"),
        ),
    )
    request = replace(base, scenario=replace(base.scenario, source=source))
    result = run_field_operational_joint_pressure_driven_history_envelope(
        request,
        _pressure_history(),
        quality_criteria=_quality(),
        history_provenance=(
            ("history_kind", "pressure_driven_orifice"),
            ("event_id", "pressure-event-A"),
        ),
        max_cases=8,
        refinement_factors=(1, 2),
        relative_tolerance=1.0,
        allow_conditional=True,
    )
    report = field_operational_joint_measured_history_envelope_report(result)

    assert len(result.cases) == 8
    assert result.cases[0].refinement is not None
    assert result.operational_decision.status == "conditional_allowed"
    assert result.operational_decision.uncertainty_resolved
    assert {dict(case.source_selection)["opening_area_m2"] for case in result.cases} == {
        "lower", "upper",
    }
    assert {dict(case.source_selection)["discharge_coefficient"] for case in result.cases} == {
        "lower", "upper",
    }
    assert {
        dict(case.field_selection)["source_location_x_m"]
        for case in result.cases
    } == {-0.1, 0.1}
    provenance = report["cases"][0]["field_result"]["transport_input"][
        "measured_history_provenance"
    ]
    assert provenance["history_kind"] == "pressure_driven_orifice"
    assert report["joint_measured_history_envelope"]["warnings"]
    json.dumps(report, allow_nan=False)


def test_operational_joint_measured_history_envelope_withholds_when_refinement_is_diagnostic_only():
    result = run_field_operational_joint_measured_history_envelope(
        _request(), _history(), quality_criteria=_quality(), max_cases=2,
        include_refinement=False,
    )

    assert result.operational_decision.status == "withheld"
    assert "refinement_missing" in result.operational_decision.gate_codes
    assert all(case.refinement is None for case in result.cases)


def test_operational_history_propagates_ambient_boundary_before_operational_gate():
    request = replace(
        _request(),
        ambient_air_density_uncertainty_kg_m3=BoundedValue(
            1.2, 1.0, 1.4, unit="kg/m3", source="met-density-boundary-A",
        ),
    )
    result = run_field_operational_joint_measured_history_envelope(
        request, _history(), quality_criteria=_quality(), max_cases=4,
        include_refinement=False,
    )

    assert result.operational_decision.status == "withheld"
    assert "refinement_missing" in result.operational_decision.gate_codes
    assert all(case.decision.uncertainty_resolved for case in result.cases)
    assert all(
        case.screening.request.ambient_air_density_uncertainty_kg_m3 is None
        for case in result.cases
    )


def test_joint_measured_history_envelope_includes_declared_stability_cases():
    request = replace(
        _request(),
        stability_mixing_closure=StabilityScalarMixingClosure(
            {"neutral": 0.5, "stable": 0.2}, evidence_id="met-mixing-A",
        ),
        stability_alternatives=FieldStabilityAlternatives(
            ("stable",), "met-stability-classification-A",
        ),
    )
    result = run_field_operational_joint_measured_history_envelope(
        request, _history(), quality_criteria=_quality(), max_cases=4,
        refinement_factors=(1, 2), relative_tolerance=1.0, allow_conditional=True,
    )

    assert len(result.cases) == 4
    assert {dict(case.field_selection)["weather_stability"] for case in result.cases} == {
        "neutral", "stable",
    }
    assert result.operational_decision.status == "conditional_allowed"
    assert all(case.decision.uncertainty_resolved for case in result.cases)


def test_operational_history_aggregate_preserves_unresolved_corner_flag():
    allowed = FieldOperationalScreeningDecision(
        "screening_allowed", True, False, (), (), True, True, True, True,
    )
    unresolved = FieldOperationalScreeningDecision(
        "withheld", False, False, ("unresolved",), (), True, True, True, False,
    )
    result = _aggregate(
        (
            SimpleNamespace(label="resolved", decision=allowed),
            SimpleNamespace(label="unresolved", decision=unresolved),
        ),
        allow_conditional=False,
    )

    assert result.status == "withheld"
    assert not result.uncertainty_resolved
    assert result.gate_codes == ("corner_withheld",)


def test_operational_history_case_contract_rejects_malformed_selections():
    decision = FieldOperationalScreeningDecision(
        "withheld", False, False, ("missing",), (), True, False, True, False,
    )
    with pytest.raises(TypeError, match="source_selection values"):
        FieldOperationalMeasuredHistoryEnvelopeCase(
            (("source", 1.0),), (), None, None, decision,
        )
    with pytest.raises(ValueError, match="field_selection keys must be unique"):
        FieldOperationalMeasuredHistoryEnvelopeCase(
            (("source", "lower"),), (("weather", 1.0), ("weather", 2.0)),
            None, None, decision,
        )
