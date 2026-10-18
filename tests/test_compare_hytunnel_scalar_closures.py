from pathlib import Path
import runpy

import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = runpy.run_path(ROOT / "tools/compare_hytunnel_scalar_closures.py")


def test_frozen_scalar_closure_comparison_rejects_local_sgs_candidate(tmp_path):
    record = TOOL["compare_closures"]()
    delta = record["candidate_minus_legacy"]

    assert record["promotion_allowed"] is False
    assert record["fixed_basis"]["common_supported_point_count"] == 45
    assert delta["mae_vol_pct_point"] == pytest.approx(
        3.4248066947066316 - 2.1854529230530475
    )
    assert delta["four_vol_specificity"] == pytest.approx(0.25 - 0.45)
    assert delta["pearson_r"] < 0.0
    fan_delta = record["measured_fan_minus_uniform_local_sgs"]
    assert fan_delta["mae_vol_pct_point"] == pytest.approx(
        3.432068099515342 - 3.4248066947066316
    )
    assert fan_delta["rmse_vol_pct_point"] > 0.0
    assert fan_delta["pearson_r"] < 0.0
    assert fan_delta["four_vol_specificity"] == pytest.approx(0.0)
    assert record["cases"][
        "measured_mean_fan_plus_local_sgs_over_sc"
    ]["fan_inlet_profile"] == "measured_first_station"
    destination = tmp_path / "comparison.png"
    TOOL["plot_comparison"](record, destination)
    assert destination.stat().st_size > 10_000
