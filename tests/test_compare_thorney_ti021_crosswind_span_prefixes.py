from __future__ import annotations

import numpy as np

from tools.compare_thorney_ti021_crosswind_span_prefixes import _relative_summary


def test_relative_summary_uses_declared_peak_and_dose_floors():
    pair = {
        "inner_index": {"S1": 0, "S2": 1},
        "outer_index": {"S1": 0, "S2": 1},
        "inner_metrics": {
            "peak": np.asarray([0.0, 1.0]),
            "dose": np.asarray([0.0, 2.0]),
        },
        "outer_metrics": {
            "peak": np.asarray([0.05, 1.2]),
            "dose": np.asarray([0.05, 2.5]),
        },
    }
    result = _relative_summary(pair, ["S1", "S2"], peak_floor=0.1, dose_floor=0.1)
    assert result["sensor_count"] == 2
    assert result["peak_change"]["maximum_percent"] == 50.0
    assert result["dose_change"]["maximum_percent"] == 50.0
    assert result["peak_change"]["median_percent"] == 35.0
    assert result["dose_change"]["median_percent"] == 37.5


def test_relative_summary_handles_empty_selection():
    pair = {
        "inner_index": {}, "outer_index": {},
        "inner_metrics": {"peak": np.asarray([]), "dose": np.asarray([])},
        "outer_metrics": {"peak": np.asarray([]), "dose": np.asarray([])},
    }
    result = _relative_summary(pair, [], peak_floor=0.1, dose_floor=0.1)
    assert result["sensor_count"] == 0
    assert result["peak_change"]["maximum_percent"] is None
    assert result["dose_change"]["p95_percent"] is None
