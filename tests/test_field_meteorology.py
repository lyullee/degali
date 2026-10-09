import pytest

from degali.addons.field_contracts import BoundedValue, WeatherState
from degali.addons.field_meteorology import (
    FieldStabilityAlternatives,
    FieldWindHistory,
    StabilityScalarMixingClosure,
)
from degali.addons.transient_receptor import WindHistory


def test_stability_scalar_mixing_closure_rejects_implicit_or_invalid_values():
    with pytest.raises(ValueError, match="evidence_id"):
        StabilityScalarMixingClosure({"neutral": 0.2}, evidence_id="unspecified")
    with pytest.raises(ValueError, match="unsupported stability"):
        StabilityScalarMixingClosure({"class-f": 0.2}, evidence_id="evidence")
    with pytest.raises(ValueError, match="positive and finite"):
        StabilityScalarMixingClosure({"neutral": 0.0}, evidence_id="evidence")


def test_stability_scalar_mixing_closure_never_falls_back_between_classes():
    values = {"neutral": 0.2}
    closure = StabilityScalarMixingClosure(values, evidence_id="evidence")
    values["neutral"] = 9.0

    assert closure.selected_diffusivity_m2_s("neutral") == pytest.approx(0.2)
    with pytest.raises(TypeError):
        closure.diffusivity_m2_s["stable"] = 0.1
    with pytest.raises(ValueError, match="no diffusivity"):
        closure.selected_diffusivity_m2_s("stable")


def test_stability_alternatives_require_evidence_and_retain_nominal_once():
    with pytest.raises(ValueError, match="at least one"):
        FieldStabilityAlternatives((), "met-review-A")
    with pytest.raises(ValueError, match="unsupported"):
        FieldStabilityAlternatives(("class-f",), "met-review-A")

    alternatives = FieldStabilityAlternatives(
        ("neutral", "stable"), "met-stability-review-A",
    )
    assert alternatives.classes_for("neutral") == ("neutral", "stable")


def test_field_wind_history_requires_provenance_and_uses_declared_steady_limits():
    history = WindHistory(
        time_s=(0.0, 0.5, 1.0),
        speed_m_s=(2.0, 2.1, 2.0),
        direction_from_deg=(270.0, 271.0, 270.0),
    )
    with pytest.raises(ValueError, match="source_id"):
        FieldWindHistory(history, source_id="unspecified", evidence_id="met-mast-cal-1")

    declared = FieldWindHistory(
        history,
        source_id="met-mast-01",
        evidence_id="met-mast-cal-1",
        maximum_direction_span_deg=5.0,
        maximum_speed_range_fraction=0.1,
    )
    assessment = declared.assess(1.0)

    assert assessment.applicable
    assert assessment.direction_span_deg == pytest.approx(1.0)
    assert assessment.speed_range_fraction == pytest.approx(0.1 / 2.033333333333333)
    record_speed, record_direction, speed_difference, direction_difference = (
        declared.nominal_weather_difference(
            WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            )
        )
    )
    assert record_speed > 0.0
    assert min(record_direction, 360.0 - record_direction) < 1.0
    assert speed_difference < 0.1
    assert direction_difference < 5.0
