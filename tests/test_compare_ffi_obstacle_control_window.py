from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from tools.compare_ffi_obstacle_control_window import _window_mean, build_report


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"
DATA = ROOT / "reference/hydrogen/open_channel_usn/extracted"


def test_window_mean_uses_time_weighting_and_interpolated_edges():
    # y=t on [0, 10], so the 2--8 s time mean is exactly 5.
    assert _window_mean(np.array([0.0, 10.0]), np.array([0.0, 10.0]), 2.0, 8.0) == pytest.approx(5.0)


def test_real_ffi_60_120_baseline_is_hash_pinned_but_remains_withheld():
    report = build_report(REGISTRY, DATA)
    assert report["window"]["start_s"] == 60.0
    assert report["window"]["end_s"] == 120.0
    assert report["promotion_allowed"] is False
    assert report["degali_score_allowed"] is False
    assert len(report["comparisons"]) == 2
    for comparison in report["comparisons"]:
        assert len(comparison["repeats"]) == 2
        for repeat in comparison["repeats"]:
            assert len(repeat["sensor_means_vol_pct"]) == 29
            assert len(repeat["control_sha256"]) == 64
            assert len(repeat["obstacle_sha256"]) == 64
            assert repeat["source_flow"]["control"]["sample_count"] > 0
            assert repeat["source_flow"]["obstacle"]["sample_count"] > 0
            groups = repeat["aggregate"]["by_sensor_placement"]
            assert sum(group["sensor_count"] for group in groups.values()) == 29
            assert {"top", "mid-high", "mid-low", "bottom"} <= set(groups)
            assert all(
                np.isfinite(group["mean_obstacle_minus_control_vol_pct"])
                for group in groups.values()
            )


def test_real_ffi_baseline_rejects_sensor_window_beyond_available_source_history():
    with pytest.raises(ValueError, match="source-flow record does not cover"):
        build_report(REGISTRY, DATA, 190.0, 210.0)
