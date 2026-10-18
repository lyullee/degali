from __future__ import annotations

from pathlib import Path
import runpy

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = runpy.run_path(ROOT / "tools/audit_hytunnel_source_mixing_sensitivity.py")


def test_sampling_recovers_bilinear_plane_at_boundary() -> None:
    x = np.asarray([0.25, 0.75])
    y = np.asarray([0.25, 0.75])
    xx, yy = np.meshgrid(x, y)
    field = 2.0 * xx + 3.0 * yy
    assert TOOL["_sample"](field, x, y, 1.0, 0.5, (1.0, 1.0)) == pytest.approx(3.5)


def test_metrics_keep_false_positives_visible() -> None:
    observed = np.asarray([5.0, 1.0, 0.2, 6.0])
    predicted = np.asarray([5.5, 4.5, 0.3, 3.5])
    result = TOOL["_metrics"](observed, predicted)
    assert result["four_vol_confusion"] == {
        "true_positive": 1,
        "true_negative": 1,
        "false_positive": 1,
        "false_negative": 1,
        "sensitivity": 0.5,
        "specificity": 0.5,
    }


def test_frozen_source_mixing_sensitivity_is_diagnostic_not_calibration(
    tmp_path: Path,
) -> None:
    record, arrays = TOOL["audit_source_mixing_sensitivity"](ROOT)
    assert record["promotion_allowed"] is False
    assert record["calibration_or_best_variant_selection_allowed"] is False
    assert record["common_point_accounting"]["common_supported_point_count"] == 45
    baseline = record["cases"]["handoff0p06_sc0p7"]["metrics_on_common_points"]
    lower_sc = record["cases"]["handoff0p06_sc0p5"]["metrics_on_common_points"]
    higher_sc = record["cases"]["handoff0p06_sc1p0"]["metrics_on_common_points"]
    conservative_handoff = record["cases"]["handoff0p05_sc0p7"]
    assert baseline["mae_vol_pct_point"] == pytest.approx(2.1130664702)
    assert lower_sc["mae_vol_pct_point"] == pytest.approx(1.6697072118)
    assert higher_sc["mae_vol_pct_point"] == pytest.approx(2.5290350563)
    assert lower_sc["four_vol_confusion"]["false_positive"] == 13
    assert baseline["four_vol_confusion"]["false_positive"] == 20
    assert conservative_handoff["derived_handoff"]["distance_m"] > 1.30
    destination = tmp_path / "sensitivity.png"
    TOOL["plot_sensitivity"](record, arrays, destination)
    assert destination.stat().st_size > 10_000
