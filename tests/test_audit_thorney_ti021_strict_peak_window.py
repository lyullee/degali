import numpy as np
import pytest

from tools.audit_thorney_ti021_strict_peak_window import select_model_peak_segment


def test_model_peak_selects_the_committed_segment_without_observations():
    times = np.asarray((0.0, 5.0, 10.0, 15.0, 20.0))
    traces = np.asarray(((0.0, 1.0, 2.0, 5.0, 4.0), (0.0, 0.5, 1.0, 2.0, 3.0)))
    segments = [
        {"segment_index": 1, "start_time_s": 0.0, "stop_time_s": 10.0},
        {"segment_index": 2, "start_time_s": 10.0, "stop_time_s": 20.0},
    ]
    result = select_model_peak_segment(times, traces, segments)
    assert result["segment_index"] == 2
    assert result["model_peak_time_s"] == pytest.approx(15.0)
    assert result["model_peak_receptor_index"] == 0
    assert result["model_peak_vol_percent"] == pytest.approx(5.0)


def test_model_peak_assigns_shared_boundary_to_later_segment():
    times = np.asarray((0.0, 10.0, 20.0))
    traces = np.asarray(((0.0, 3.0, 1.0),))
    segments = [
        {"segment_index": 1, "start_time_s": 0.0, "stop_time_s": 10.0},
        {"segment_index": 2, "start_time_s": 10.0, "stop_time_s": 20.0},
    ]
    assert select_model_peak_segment(times, traces, segments)["segment_index"] == 2


def test_model_peak_rejects_nonfinite_trace():
    with pytest.raises(ValueError, match="finite traces"):
        select_model_peak_segment(
            np.asarray((0.0, 1.0)), np.asarray(((0.0, np.nan),)),
            [{"segment_index": 1, "start_time_s": 0.0, "stop_time_s": 1.0}],
        )
