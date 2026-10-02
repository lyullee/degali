import pytest

from degali.validation.model_comparison import (
    ComparisonCase,
    ModelPrediction,
    compare_models,
)


def _case(**changes):
    values = dict(
        source_mode="horizontal_jet",
        source_rate_kg_s=0.12,
        wind_speed_m_s=2.0,
        wind_direction_rad=0.0,
        release_height_m=1.0,
        receptor_operator="arc_max_1s_peak",
        averaging_time_s=1.0,
        phase_closure="measured_throat",
        geometry="open_horizontal",
    )
    values.update(changes)
    return ComparisonCase(**values)


def test_matching_uncalibrated_cases_can_be_ranked_without_raw_data_in_report():
    case = _case()
    report = compare_models(
        [1.0, 2.0, 4.0],
        [
            ModelPrediction("degali", [1.0, 2.0, 4.0], case, provenance="local"),
            ModelPrediction("oracle", [1.1, 2.2, 3.8], case, provenance="controlled"),
        ],
        case=case,
    )
    assert report.ranking_allowed
    assert report.best_fac2_model == "degali"
    assert report.eligible_models == ("degali", "oracle")
    payload = report.as_dict()
    assert "values" not in payload["results"]["degali"]
    assert payload["results"]["degali"]["n"] == 3


def test_different_receptor_or_wind_cannot_support_superiority_claim():
    case = _case()
    mismatch = _case(wind_speed_m_s=0.5, receptor_operator="centreline_60s_mean")
    report = compare_models(
        [1.0, 2.0],
        [
            ModelPrediction("degali", [1.0, 2.0], case),
            ModelPrediction("other", [1.1, 2.1], mismatch),
        ],
        case=case,
    )
    assert not report.ranking_allowed
    assert report.best_fac2_model is None
    assert "wind_speed_m_s" in report.excluded_models["other"]
    assert "receptor_operator" in report.excluded_models["other"]


def test_calibrated_prediction_is_reported_but_excluded_from_ranking():
    case = _case()
    report = compare_models(
        [1.0, 2.0],
        [
            ModelPrediction("fitted", [1.0, 2.0], case, calibrated_to_observations=True),
            ModelPrediction("independent", [1.2, 1.8], case),
        ],
        case=case,
    )
    assert not report.ranking_allowed
    assert "calibrated" in report.excluded_models["fitted"]
    assert "independent" in report.eligible_models


def test_invalid_case_and_duplicate_model_names_are_rejected():
    with pytest.raises(ValueError, match="positive"):
        _case(source_rate_kg_s=0.0)
    case = _case()
    with pytest.raises(ValueError, match="unique"):
        compare_models(
            [1.0],
            [ModelPrediction("same", [1.0], case), ModelPrediction("same", [1.0], case)],
            case=case,
        )
