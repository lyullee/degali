import json
import math
from dataclasses import replace
from types import SimpleNamespace

import pytest

from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_history import (
    FieldMeasuredFlashSchedule,
    MeasuredHistoryQualityCriteria,
    MeasuredReleaseHistory,
    MeasuredTimeSeries,
    PressureDrivenMeasuredHistory,
    assess_measured_history_quality,
    direct_vapour_schedule_envelope_from_measured_history,
    direct_vapour_schedule_from_measured_history,
    direct_vapour_schedule_envelope_from_pressure_driven_history,
    direct_vapour_schedule_from_pressure_driven_history,
    request_with_measured_flash_schedule,
    run_field_joint_measured_history_envelope,
    run_field_joint_pressure_driven_history_envelope,
    run_field_measured_history_envelope,
    run_field_pressure_driven_history_envelope,
)
from degali.addons.field_workflow import (
    FieldSensorDeployment,
    FieldSemiFVRequest,
    run_field_semi_fv_screening,
    run_field_semi_fv_refinement_study,
)
import degali.addons.field_history as field_history
from degali.addons.field_decision import evaluate_field_operational_screening
from degali.addons.field_report import (
    field_measured_history_envelope_report,
    field_joint_measured_history_envelope_report,
    field_screening_report,
)
from degali.addons.semi_fv_obstacle import SemiFVConfig
from degali.addons.semi_fv_obstacle import SourceRateSchedule


def _release(duration_s=2.0):
    return ReleaseSource(
        fluid="lh2",
        upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.265, unit="kg/s"),
        opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"),
        discharge_coefficient=BoundedValue(0.8),
        liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium",
        duration_s=duration_s,
    )


def _history(*, pressure_end=2.0):
    pressure = MeasuredTimeSeries(
        time_s=(0.0, 0.5 * pressure_end, pressure_end),
        nominal=(0.4e6, 0.41e6, 0.42e6),
        unit="Pa", source_id="PT-101",
    )
    # The temperature channel was recorded 0.1 s before the release clock;
    # its explicit offset aligns it at the same physical nodes as flow.
    temperature = MeasuredTimeSeries(
        time_s=(-0.1, 0.9, 1.9),
        nominal=(26.084, 26.2, 26.3),
        unit="K", source_id="TT-101", time_offset_s=0.1, response_time_s=0.25,
    )
    flow = MeasuredTimeSeries(
        time_s=(0.0, 1.0, 2.0),
        nominal=(0.265, 0.20, 0.0),
        lower=(0.25, 0.18, 0.0),
        upper=(0.28, 0.22, 0.0),
        unit="kg/s", source_id="FT-101",
    )
    return MeasuredReleaseHistory(pressure, temperature, flow)


def _quality(*, maximum_response_time_s=0.5):
    return MeasuredHistoryQualityCriteria(
        maximum_sample_interval_s=1.0,
        maximum_response_time_s=maximum_response_time_s,
        maximum_absolute_time_offset_s=0.2,
        maximum_relative_half_width=0.2,
        evidence_id="historian-quality-procedure-01",
    )


def test_measured_flash_schedule_requires_closed_typed_source_and_matching_mass():
    schedule = SourceRateSchedule((0.0, 1.0), (1.0, 0.0), source_id="history-source")
    with pytest.raises(ValueError, match="schedule mass"):
        FieldMeasuredFlashSchedule(
            schedule=schedule,
            history_duration_s=1.0,
            selection=(),
            total_measured_mass_kg=1.0,
            direct_vapour_mass_kg=0.5,
            unrouted_postflash_liquid_mass_kg=0.5,
            alignment_method="flow-clock",
            warnings=(),
        )
    with pytest.raises(ValueError, match="final endpoint rate"):
        FieldMeasuredFlashSchedule(
            schedule=SourceRateSchedule((0.0, 1.0), (1.0, 0.1), source_id="history-source"),
            history_duration_s=1.0,
            selection=(),
            total_measured_mass_kg=1.0,
            direct_vapour_mass_kg=1.0,
            unrouted_postflash_liquid_mass_kg=0.0,
            alignment_method="flow-clock",
            warnings=(),
        )


def _pressure_driven_history(*, bounded=False):
    pressure = MeasuredTimeSeries(
        time_s=(0.0, 1.0, 2.0),
        nominal=(0.4e6, 0.39e6, 0.38e6),
        lower=((0.39e6, 0.38e6, 0.37e6) if bounded else None),
        upper=((0.41e6, 0.40e6, 0.39e6) if bounded else None),
        unit="Pa", source_id="PT-ORIFICE-01",
    )
    temperature = MeasuredTimeSeries(
        time_s=(0.0, 1.0, 2.0),
        nominal=(26.076, 25.949, 25.819),
        lower=((26.03, 25.90, 25.77) if bounded else None),
        upper=((26.13, 26.00, 25.87) if bounded else None),
        unit="K", source_id="TT-ORIFICE-01",
    )
    return PressureDrivenMeasuredHistory(
        pressure, temperature,
        event_id="event-pressure-driven-01",
        phase_evidence_id="phase-evidence-01",
    )


def test_pressure_driven_history_recomputes_mass_flow_and_flash_coherently():
    release = replace(_release(), mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"))
    schedule = direct_vapour_schedule_from_pressure_driven_history(
        release,
        _pressure_driven_history(),
        source_id="orifice-history-01",
        quality_criteria=_quality(),
    )

    assert schedule.schedule.source_id == "orifice-history-01"
    assert schedule.total_measured_mass_kg > schedule.direct_vapour_mass_kg > 0.0
    assert schedule.unrouted_postflash_liquid_mass_kg > 0.0
    assert schedule.quality_assessment is not None
    assert schedule.quality_assessment.approved
    assert schedule.source_uncertainty_resolved
    assert "pressure trend alone was not treated as a leak rate" in schedule.warnings[0]
    assert dict(schedule.provenance)["history_kind"] == "pressure_driven_orifice"


def test_pressure_driven_history_envelope_keeps_pressure_and_temperature_corners_coupled():
    release = replace(_release(), mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"))
    envelope = direct_vapour_schedule_envelope_from_pressure_driven_history(
        release,
        _pressure_driven_history(bounded=True),
        source_id="orifice-history-envelope-01",
        quality_criteria=_quality(),
        max_cases=4,
    )

    assert len(envelope.cases) == 4
    assert {dict(case.selection)["pressure_pa"] for case in envelope.cases} == {
        "lower", "upper",
    }
    assert {dict(case.selection)["temperature_k"] for case in envelope.cases} == {
        "lower", "upper",
    }
    masses = {case.schedule.total_measured_mass_kg for case in envelope.cases}
    assert len(masses) == 4
    assert all(case.schedule.quality_assessment.approved for case in envelope.cases)


def test_pressure_driven_history_envelope_crosses_orifice_source_bounds_coherently():
    base = _release()
    area = base.opening_area_m2.nominal
    release = replace(
        base,
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        opening_area_m2=BoundedValue(area, 0.9 * area, 1.1 * area, "m2", "area-review"),
        discharge_coefficient=BoundedValue(0.8, 0.72, 0.88, "1", "cd-review"),
    )
    envelope = direct_vapour_schedule_envelope_from_pressure_driven_history(
        release,
        _pressure_driven_history(bounded=True),
        source_id="orifice-source-envelope-01",
        quality_criteria=_quality(),
        max_cases=16,
    )

    assert len(envelope.cases) == 16
    assert {dict(case.selection)["opening_area_m2"] for case in envelope.cases} == {
        "lower", "upper",
    }
    assert {dict(case.selection)["discharge_coefficient"] for case in envelope.cases} == {
        "lower", "upper",
    }
    assert len({case.schedule.total_measured_mass_kg for case in envelope.cases}) > 4


def test_pressure_driven_history_refuses_unpropagated_or_measured_rate_inputs():
    history = _pressure_driven_history()
    with pytest.raises(ValueError, match="exact opening_area_m2"):
        direct_vapour_schedule_from_pressure_driven_history(
            replace(
                replace(_release(), mass_flow_kg_s=BoundedValue(0.0, unit="kg/s")),
                opening_area_m2=BoundedValue(1.0e-4, 0.9e-4, 1.1e-4, "m2", "area-review"),
            ),
            history,
        )
    with pytest.raises(ValueError, match="will not overwrite"):
        direct_vapour_schedule_from_pressure_driven_history(_release(), history)


def test_measured_history_aligns_channels_and_builds_direct_vapour_schedule():
    history = _history()
    aligned = history.temperature_k.values_at(history.mass_flow_kg_s.physical_time_s)
    result = direct_vapour_schedule_from_measured_history(_release(), history)

    assert aligned == pytest.approx((26.084, 26.2, 26.3))
    assert result.schedule.time_s == pytest.approx((0.0, 1.0, 2.0))
    assert result.schedule.rate_kg_s[-1] == 0.0
    assert result.total_measured_mass_kg == pytest.approx(0.465)
    assert result.direct_vapour_mass_kg > 0.0
    assert result.unrouted_postflash_liquid_mass_kg > 0.0
    assert result.direct_vapour_mass_kg + result.unrouted_postflash_liquid_mass_kg == pytest.approx(
        result.total_measured_mass_kg, rel=1.0e-8
    )
    assert any("not deconvolved" in warning for warning in result.warnings)


def test_measured_history_blocks_an_unclosed_flash_interval(monkeypatch):
    nominal = field_history.lh2_flash_source_from_release(_release())
    invalid = replace(
        nominal,
        flash=replace(nominal.flash, energy_residual=1.0e-3),
    )
    monkeypatch.setattr(
        field_history,
        "lh2_flash_source_from_release",
        lambda *args, **kwargs: invalid,
    )

    with pytest.raises(ValueError, match="not conservatively closed"):
        direct_vapour_schedule_from_measured_history(_release(), _history())


def test_measured_history_blocks_a_nonclosing_phase_ledger(monkeypatch):
    nominal = field_history.lh2_flash_source_from_release(_release())
    malformed = SimpleNamespace(
        conservative=True,
        flash=SimpleNamespace(
            vapour_mass_flow=0.5 * nominal.flash.vapour_mass_flow,
            liquid_mass_flow=0.0,
        ),
    )
    monkeypatch.setattr(
        field_history,
        "lh2_flash_source_from_release",
        lambda *args, **kwargs: malformed,
    )

    with pytest.raises(ValueError, match="measured-mass ledger"):
        direct_vapour_schedule_from_measured_history(_release(), _history())


def test_measured_history_refuses_interpolation_outside_declared_channel_coverage():
    with pytest.raises(ValueError, match="does not cover"):
        direct_vapour_schedule_from_measured_history(
            _release(), _history(pressure_end=1.0),
        )


def test_measured_history_envelope_uses_global_measurement_bound_corners():
    envelope = direct_vapour_schedule_envelope_from_measured_history(
        _release(), _history(), max_cases=4,
    )

    assert len(envelope.cases) == 2
    assert {dict(case.selection)["mass_flow_kg_s"] for case in envelope.cases} == {
        "lower", "upper",
    }
    assert all(case.schedule.direct_vapour_mass_kg > 0.0 for case in envelope.cases)
    assert "not a probability interval" in envelope.warnings[0]


def test_measured_history_boundaries_reject_invalid_and_duplicate_selections():
    schedule = direct_vapour_schedule_from_measured_history(_release(), _history())
    with pytest.raises(ValueError, match="values must be nominal"):
        replace(schedule, selection=(("mass_flow_kg_s", "fabricated"),))
    with pytest.raises(ValueError, match="keys must be unique"):
        replace(
            schedule,
            selection=(
                ("mass_flow_kg_s", "lower"),
                ("mass_flow_kg_s", "upper"),
            ),
        )
    envelope = direct_vapour_schedule_envelope_from_measured_history(
        _release(), _history(), max_cases=4,
    )
    with pytest.raises(ValueError, match="selections must be unique"):
        replace(envelope, cases=(envelope.cases[0], envelope.cases[0]))


def test_measured_history_quality_gate_requires_declared_timing_and_uncertainty_limits():
    approved = assess_measured_history_quality(_history(), _quality())
    blocked = assess_measured_history_quality(
        _history(), _quality(maximum_response_time_s=0.1),
    )

    assert approved.approved
    assert not blocked.approved
    assert "temperature_k:response_time_exceeds_declared_limit" in blocked.reasons


def test_measured_flash_schedule_enters_field_transport_with_its_alignment_warnings():
    history = _history()
    schedule = direct_vapour_schedule_from_measured_history(
        _release(), history, quality_criteria=_quality(),
    )
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    result = run_field_semi_fv_screening(
        request_with_measured_flash_schedule(request, schedule)
    )

    assert result.completed
    assert result.transport is not None
    assert result.transport.diagnostics.source_mode == "scheduled"
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(
        schedule.direct_vapour_mass_kg
    )
    assert any("not deconvolved" in warning for warning in result.applicability.warnings)
    assert result.request.measured_history_quality is not None
    assert result.request.measured_history_quality.approved
    report = field_screening_report(result)
    quality = report["transport_input"]["measured_history_quality"]
    assert quality is not None
    assert quality["approved"] is True
    assert quality["criteria"]["evidence_id"] == "historian-quality-procedure-01"
    assert {item["channel"] for item in quality["channel_metrics"]} == {
        "pressure_pa", "temperature_k", "mass_flow_kg_s",
    }


def test_measured_history_rejects_unpropagated_release_duration_uncertainty():
    release = replace(
        _release(),
        duration_uncertainty=BoundedValue(
            _release().duration_s, _release().duration_s - 0.1,
            _release().duration_s + 0.1, "s", "timing-review-A",
        ),
    )
    with pytest.raises(ValueError, match="time-aligned source envelope"):
        direct_vapour_schedule_from_measured_history(
            release, _history(), quality_criteria=_quality(),
        )


def test_measured_history_schedule_cannot_enter_field_transport_without_approved_quality():
    schedule = direct_vapour_schedule_from_measured_history(_release(), _history())
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )

    with pytest.raises(ValueError, match="approved quality assessment"):
        request_with_measured_flash_schedule(request, schedule)


def test_operational_gate_withholds_a_measured_nominal_schedule_until_its_source_bounds_are_propagated():
    schedule = direct_vapour_schedule_from_measured_history(
        _release(), _history(), quality_criteria=_quality(),
    )
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    attached = request_with_measured_flash_schedule(request, schedule)
    refined = run_field_semi_fv_refinement_study(
        attached, refinement_factors=(1, 2), relative_tolerance=1.0,
    )
    decision = evaluate_field_operational_screening(refined, allow_conditional=True)

    assert refined.completed
    assert not attached.measured_history_source_uncertainty_resolved
    assert decision.status == "withheld"
    assert "measured-history source" in " ".join(decision.reasons)


def test_measured_history_envelope_propagates_source_bounds_to_field_transport():
    source = replace(
        _release(),
        location_m=(0.0, 0.0, 0.5),
        location_uncertainty_m=(
            BoundedValue(0.0, lower=-0.1, upper=0.1, unit="m", source="survey"),
            BoundedValue(0.0, unit="m", source="survey"),
            BoundedValue(0.5, unit="m", source="survey"),
        ),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    envelope = run_field_measured_history_envelope(
        request, _history(), quality_criteria=_quality(), max_cases=4,
    )

    assert envelope.quality_assessment.approved
    assert len(envelope.cases) == 4
    assert envelope.completed_case_count == 4
    assert envelope.property_table_used
    assert {dict(case.selection)["mass_flow_kg_s"] for case in envelope.cases} == {
        "lower", "upper",
    }
    assert {dict(case.selection)["source_location_x_m"] for case in envelope.cases} == {
        -0.1, 0.1,
    }
    for case in envelope.cases:
        assert case.result.transport is not None
        assert case.result.transport.diagnostics.mass_injected_kg == pytest.approx(
            case.schedule.direct_vapour_mass_kg,
        )
    assert "not a probability interval" in envelope.warnings[0]
    assert "source location/direction corners" in " ".join(envelope.warnings)
    report = field_measured_history_envelope_report(envelope)
    assert report["completed_case_count"] == 4
    assert len(report["cases"]) == 4
    assert report["quality_assessment"]["approved"] is True
    assert report["property_table_used"] is True
    json.dumps(report, allow_nan=False)


def test_measured_history_field_envelope_refuses_an_unapproved_historian():
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )

    with pytest.raises(ValueError, match="approved quality criteria"):
        run_field_measured_history_envelope(
            request, _history(), quality_criteria=_quality(maximum_response_time_s=0.1),
        )


def test_pressure_driven_history_envelope_recomputes_source_before_field_transport():
    source = replace(
        _release(),
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        location_m=(0.0, 0.0, 0.5),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    envelope = run_field_pressure_driven_history_envelope(
        request,
        _pressure_driven_history(bounded=True),
        quality_criteria=_quality(),
        max_cases=4,
    )

    assert envelope.quality_assessment.approved
    assert len(envelope.cases) == 4
    assert envelope.completed_case_count == 4
    assert envelope.property_table_used
    assert all(case.schedule.total_measured_mass_kg > 0.0 for case in envelope.cases)
    assert all(
        case.result.transport is not None
        and case.result.transport.diagnostics.mass_injected_kg == pytest.approx(
            case.schedule.direct_vapour_mass_kg,
        )
        for case in envelope.cases
    )
    assert all(
        any("history was re-flashed through the explicit pressure-driven throat closure" in warning
            for warning in case.result.applicability.warnings)
        for case in envelope.cases
    )
    assert any("recomputes throat flow" in warning for warning in envelope.warnings)
    report = field_measured_history_envelope_report(envelope)
    assert all(
        case["field_screening"]["transport_input"]["measured_history_provenance"]["history_kind"]
        == "pressure_driven_orifice"
        for case in report["cases"]
    )
    json.dumps(report, allow_nan=False)


def test_pressure_driven_history_field_envelope_propagates_source_geometry_corners():
    source = replace(
        _release(),
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        location_m=(0.0, 0.0, 0.5),
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-pressure"),
            BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-pressure"),
            BoundedValue(0.5, 0.5, 0.5, "m", "layout-review-pressure"),
        ),
        direction_m=(1.0, 0.0, 0.0),
        direction_uncertainty_m=(
            BoundedValue(1.0, 1.0, 1.0, "1", "direction-review-pressure"),
            BoundedValue(0.0, -0.1, 0.1, "1", "direction-review-pressure"),
            BoundedValue(0.0, 0.0, 0.0, "1", "direction-review-pressure"),
        ),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    envelope = run_field_pressure_driven_history_envelope(
        request,
        _pressure_driven_history(),
        quality_criteria=_quality(),
        max_cases=4,
    )

    assert len(envelope.cases) == 4
    assert {
        dict(case.selection)["source_location_x_m"]
        for case in envelope.cases
    } == {-0.1, 0.1}
    assert all(
        case.result.request.scenario.source.location_uncertainty_m is None
        and case.result.request.scenario.source.direction_uncertainty_m is None
        for case in envelope.cases
    )
    assert {
        case.result.request.scenario.source.location_m[0]
        for case in envelope.cases
    } == {-0.1, 0.1}
    assert {
        case.result.request.scenario.source.direction_m[1]
        for case in envelope.cases
    } == {-0.1, 0.1}
    assert any("source location and direction corners" in warning for warning in envelope.warnings)


def test_joint_pressure_driven_history_envelope_keeps_source_geometry_in_field_selection():
    area = _release().opening_area_m2.nominal
    source = replace(
        _release(),
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        opening_area_m2=BoundedValue(area, 0.9 * area, 1.1 * area, "m2", "area-review"),
        discharge_coefficient=BoundedValue(0.8, 0.72, 0.88, "1", "cd-review"),
        location_m=(0.0, 0.0, 0.5),
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-pressure"),
            BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-pressure"),
            BoundedValue(0.5, 0.5, 0.5, "m", "layout-review-pressure"),
        ),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    envelope = run_field_joint_pressure_driven_history_envelope(
        request,
        _pressure_driven_history(),
        quality_criteria=_quality(),
        max_cases=8,
    )

    assert len(envelope.cases) == 8
    assert {
        dict(case.field_values)["source_location_x_m"]
        for case in envelope.cases
    } == {-0.1, 0.1}
    assert all(
        case.result.request.scenario.source.location_uncertainty_m is None
        for case in envelope.cases
    )
    assert {
        item["opening_area_m2"] for item in
        (dict(case.source_selection) for case in envelope.cases)
    } == {"lower", "upper"}


def test_pressure_driven_field_envelope_rejects_unresolved_ambient_pressure():
    source = replace(
        _release(),
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        location_m=(0.0, 0.0, 0.5),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        ambient_pressure_uncertainty_pa=BoundedValue(
            101325.0, 100000.0, 102000.0, "Pa", "ambient-pressure-review",
        ),
    )

    with pytest.raises(ValueError, match="one fixed ambient pressure"):
        run_field_pressure_driven_history_envelope(
            request, _pressure_driven_history(), quality_criteria=_quality(),
        )


def test_joint_pressure_driven_history_envelope_propagates_ambient_source_corners():
    area = _release().opening_area_m2.nominal
    source = replace(
        _release(),
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        location_m=(0.0, 0.0, 0.5),
        opening_area_m2=BoundedValue(area, 0.9 * area, 1.1 * area, "m2", "area-review"),
        discharge_coefficient=BoundedValue(0.8, 0.72, 0.88, "1", "cd-review"),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        ambient_temperature_uncertainty_k=BoundedValue(
            295.0, 294.0, 296.0, "K", "ambient-temperature-review",
        ),
        ambient_pressure_uncertainty_pa=BoundedValue(
            101325.0, 100000.0, 102000.0, "Pa", "ambient-pressure-review",
        ),
    )
    envelope = run_field_joint_pressure_driven_history_envelope(
        request,
        _pressure_driven_history(bounded=True),
        quality_criteria=_quality(),
        max_cases=64,
    )

    assert envelope.quality_assessment.approved
    assert len(envelope.cases) == 64
    assert envelope.completed_case_count == 64
    ambient_pressures = {
        dict(case.field_values)["ambient_pressure_pa"]
        for case in envelope.cases
    }
    ambient_temperatures = {
        dict(case.field_values)["ambient_temperature_k"]
        for case in envelope.cases
    }
    assert ambient_pressures == {100000.0, 102000.0}
    assert ambient_temperatures == {294.0, 296.0}
    source_selections = [dict(case.source_selection) for case in envelope.cases]
    assert {item["opening_area_m2"] for item in source_selections} == {"lower", "upper"}
    assert {item["discharge_coefficient"] for item in source_selections} == {
        "lower", "upper",
    }
    assert any("throat recomputation" in warning for warning in envelope.warnings)


def test_joint_measured_history_envelope_combines_only_independent_weather_and_sensor_bounds():
    source = replace(
        _release(),
        location_m=(0.0, 0.0, 0.5),
        location_uncertainty_m=(
            BoundedValue(0.0, -0.1, 0.1, "m", "layout-review-measured"),
            BoundedValue(0.0, 0.0, 0.0, "m", "layout-review-measured"),
            BoundedValue(0.5, 0.5, 0.5, "m", "layout-review-measured"),
        ),
    )
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0, 1.5, 2.5, "m/s", "met-mast-01"),
                direction_deg=BoundedValue(270.0, unit="deg"),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )
    envelope = run_field_joint_measured_history_envelope(
        request, _history(), quality_criteria=_quality(), max_cases=8,
    )

    assert envelope.quality_assessment.approved
    assert len(envelope.cases) == 8
    assert envelope.completed_case_count == 8
    assert envelope.property_table_used
    assert {dict(case.source_selection)["mass_flow_kg_s"] for case in envelope.cases} == {
        "lower", "upper",
    }
    assert {dict(case.field_values)["wind_speed_m_s"] for case in envelope.cases} == {
        1.5, 2.5,
    }
    assert {
        dict(case.field_values)["source_location_x_m"]
        for case in envelope.cases
    } == {-0.1, 0.1}
    assert all(
        "opening_area_m2" not in dict(case.field_values)
        and "mass_flow_kg_s" not in dict(case.field_values)
        for case in envelope.cases
    )
    assert any("opening area and discharge coefficient" in warning for warning in envelope.warnings)
    assert any("source location/direction" in warning for warning in envelope.warnings)
    report = field_joint_measured_history_envelope_report(envelope)
    assert report["completed_case_count"] == 8
    assert len(report["cases"]) == 8
    assert report["schema"] == "degali.field-joint-measured-history-envelope.v1"
    assert report["property_table_used"] is True
    json.dumps(report, allow_nan=False)


def test_joint_measured_history_envelope_propagates_ambient_flash_and_sensor_bounds():
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0),
                direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        ambient_temperature_uncertainty_k=BoundedValue(
            295.0, 294.0, 296.0, "K", "ambient-mast-01",
        ),
        ambient_pressure_uncertainty_pa=BoundedValue(
            101325.0, 90_000.0, 110_000.0, "Pa", "ambient-mast-01",
        ),
        ambient_air_density_uncertainty_kg_m3=BoundedValue(
            1.2, 1.1, 1.3, "kg/m3", "ambient-mast-01",
        ),
    )
    envelope = run_field_joint_measured_history_envelope(
        request, _history(), quality_criteria=_quality(), max_cases=16,
    )

    assert len(envelope.cases) == 16  # two source bounds x 2^3 ambient corners
    assert envelope.completed_case_count == 16
    assert envelope.property_table_used
    selections = [dict(case.field_values) for case in envelope.cases]
    assert {selection["ambient_temperature_k"] for selection in selections} == {294.0, 296.0}
    assert {selection["ambient_pressure_pa"] for selection in selections} == {90_000.0, 110_000.0}
    assert {selection["ambient_air_density_kg_m3"] for selection in selections} == {1.1, 1.3}
    assert all(
        case.result.request.ambient_temperature_uncertainty_k is None
        and case.result.request.ambient_pressure_uncertainty_pa is None
        and case.result.request.ambient_air_density_uncertainty_kg_m3 is None
        for case in envelope.cases
    )
    assert any("ambient temperature, pressure and air-density" in warning
               for warning in envelope.warnings)


def test_joint_measured_history_envelope_fails_before_running_too_many_cases():
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0, 1.5, 2.5, "m/s", "met-mast-01"),
                direction_deg=BoundedValue(270.0, unit="deg"),
            ), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )

    with pytest.raises(ValueError, match="joint measured-history/field corners exceed max_cases=2"):
        run_field_joint_measured_history_envelope(
            request, _history(), quality_criteria=_quality(), max_cases=2,
        )


def test_joint_measured_history_envelope_reuses_transport_for_named_sensor_calibration():
    source = replace(_release(), location_m=(0.0, 0.0, 0.5))
    request = FieldSemiFVRequest(
        FieldScenario(
            source=source,
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=2.0, time_step_s=0.005, source_sigma_m=0.2,
        ),
        sensor_deployments=(
            FieldSensorDeployment(
                "detector-02",
                SensorModel(
                    (0.8, 0.0, 0.5),
                    gain=BoundedValue(1.0, 0.8, 1.2, "1", "detector-02-cal"),
                ),
            ),
        ),
    )
    envelope = run_field_joint_measured_history_envelope(
        request, _history(), quality_criteria=_quality(), max_cases=4,
    )

    assert len(envelope.cases) == 4  # two flow-history × two detector gains
    assert {dict(case.field_values)["detector-02.sensor_gain"] for case in envelope.cases} == {
        0.8, 1.2,
    }
    for source_bound in ("lower", "upper"):
        pair = [
            case for case in envelope.cases
            if dict(case.source_selection)["mass_flow_kg_s"] == source_bound
        ]
        assert len(pair) == 2
        assert pair[0].result.transport is pair[1].result.transport
        gain_to_peak = {
            dict(case.field_values)["detector-02.sensor_gain"]: next(
                result for result in case.result.sensor_results
                if result.label == "detector-02"
            ).trace.indicated_mole_fraction.max()
            for case in pair
        }
        assert gain_to_peak[1.2] > gain_to_peak[0.8]
        assert all(
            "uses nominal calibration values in this single screen"
            not in warning for case in pair for warning in case.result.applicability.warnings
        )
