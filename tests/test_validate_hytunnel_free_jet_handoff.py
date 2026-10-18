from pathlib import Path
import runpy

import pytest


ROOT = Path("reference/datasets/hytunnel-ps-4668554/raw")


@pytest.mark.skipif(
    not (ROOT / "HTE2440PS002MIXED_300321.xlsx").is_file(),
    reason="requires the locally controlled public HyTunnel data",
)
def test_no_fit_hytunnel_free_jet_route_is_rejected_without_hiding_correlation(
    tmp_path,
):
    tool = runpy.run_path("tools/validate_hytunnel_free_jet_handoff.py")
    record, arrays = tool["compare_hytunnel_free_jet_handoff"](ROOT)

    assert record["status"] == "diagnostic_component_route_rejected"
    assert record["source"]["source_model_match_to_official_cfd"] is False
    assert record["source"]["handoff_maximum_conservation_residual"] < 1.0e-12
    assert record["source"]["first_model_coordinate_m"] == pytest.approx(
        record["source"]["handoff_distance_m"]
    )
    metrics = record["paired_average_concentrations"]["all_downstream_points"]
    assert metrics["pair_count"] == 108
    assert metrics["mae_vol_pct_point"] == pytest.approx(2.68784, rel=2.0e-5)
    assert metrics["rmse_vol_pct_point"] == pytest.approx(3.15966, rel=2.0e-5)
    assert metrics["mean_bias_predicted_minus_observed_vol_pct_point"] > 0.0
    assert metrics["pearson_r"] > 0.93
    hanna = record["paired_average_concentrations"]["hanna_statistics"]
    assert hanna["fac2"] == pytest.approx(0.518868, rel=2.0e-5)
    assert hanna["acceptable"] is False
    threshold = record["threshold_classification"]["four_vol_pct"]
    assert threshold["sensitivity"] == 1.0
    assert threshold["false_positive"] == 39
    assert threshold["precision"] == pytest.approx(27.0 / 66.0)
    lengths = record["centreline_threshold_lengths"]
    assert lengths["degali_free_jet_4_vol_pct_m"] == pytest.approx(3.41307, rel=2e-5)
    assert lengths["degali_relative_error_vs_experiment"] > 1.0
    assert not record["diagnostic_gates"]["fac2_at_floor_0p10_at_least_0p70"]
    assert arrays["observed_vol_pct"].shape == arrays["predicted_vol_pct"].shape

    plot = tmp_path / "comparison.png"
    tool["plot_comparison"](arrays, plot)
    assert plot.is_file() and plot.stat().st_size > 10_000

