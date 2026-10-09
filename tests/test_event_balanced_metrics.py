import pytest

from degali.addons.event_balanced_metrics import (
    EventMetricObservation,
    score_event_balanced_metrics,
)


def test_event_balanced_score_reports_pooled_and_macro_event_metrics():
    score = score_event_balanced_metrics(
        (
            EventMetricObservation("large", "a", 0.0, 0.1),
            EventMetricObservation("large", "b", 0.0, 0.1),
            EventMetricObservation("small", "a", 0.0, 0.0),
        ),
        threshold_mole_fraction=0.04,
    )

    assert score.event_count == 2
    assert score.pooled_mean_absolute_error_mole_fraction == pytest.approx(0.0666666666667)
    assert score.macro_mean_absolute_error_mole_fraction == pytest.approx(0.05)
    assert score.threshold_counts == {"tp": 0, "fn": 0, "fp": 2, "tn": 1, "scored_exact_rows": 3}
    assert score.events[0].event_id == "large"
    assert score.as_dict()["schema"] == "degali.event-balanced-metrics.v1"


def test_lower_bound_rows_are_not_used_as_symmetric_error():
    score = score_event_balanced_metrics(
        (
            EventMetricObservation("event", "exact", 0.01, 0.02),
            EventMetricObservation("event", "bound", 0.20, 0.15, "lower_bound"),
        ),
    )

    assert score.exact_row_count == 1
    assert score.lower_bound_row_count == 1
    assert score.pooled_mean_absolute_error_mole_fraction == pytest.approx(0.01)
    assert score.lower_bound_satisfied_count == 0
    assert score.lower_bound_satisfaction_fraction == pytest.approx(0.0)
    assert score.threshold_counts["scored_exact_rows"] == 1


def test_duplicate_event_sensor_pairs_are_rejected():
    rows = (
        EventMetricObservation("event", "sensor", 0.01, 0.02),
        EventMetricObservation("event", "sensor", 0.01, 0.03),
    )
    with pytest.raises(ValueError, match="pairs must be unique"):
        score_event_balanced_metrics(rows)


def test_invalid_threshold_and_concentration_are_rejected():
    with pytest.raises(ValueError, match=r"lie in \[0, 1\]"):
        EventMetricObservation("e", "s", 1.1, 0.0)
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        score_event_balanced_metrics(
            (EventMetricObservation("e", "s", 0.0, 0.0),),
            threshold_mole_fraction=1.1,
        )


def test_event_without_exact_rows_has_no_symmetric_macro_metric():
    score = score_event_balanced_metrics(
        (EventMetricObservation("event", "sensor", 0.2, 0.1, "lower_bound"),),
    )
    assert score.pooled_mean_absolute_error_mole_fraction is None
    assert score.macro_mean_absolute_error_mole_fraction is None
    assert score.lower_bound_satisfaction_fraction == pytest.approx(0.0)
