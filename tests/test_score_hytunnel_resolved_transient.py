from __future__ import annotations

from pathlib import Path
import runpy

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def _tool():
    return runpy.run_path(ROOT / "tools/score_hytunnel_resolved_transient.py")


def test_bilinear_cell_centre_sample_interpolates_and_reaches_face() -> None:
    sample = _tool()["bilinear_cell_centre_sample"]
    x = np.asarray([0.25, 0.75])
    y = np.asarray([0.25, 0.75])
    xx, yy = np.meshgrid(x, y)
    field = 2.0 * xx + 3.0 * yy
    assert sample(field, x, y, 0.5, 0.5, domain_m=(1.0, 1.0)) == pytest.approx(2.5)
    assert sample(field, x, y, 1.0, 0.5, domain_m=(1.0, 1.0)) == pytest.approx(3.5)
    with pytest.raises(ValueError, match="outside"):
        sample(field, x, y, 1.01, 0.5, domain_m=(1.0, 1.0))


def test_frozen_resolved_score_has_honest_overlap_and_withholds_promotion() -> None:
    record, arrays = _tool()["score_hytunnel_resolved_transient"](
        ROOT,
        ROOT / "outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json",
    )
    accounting = record["point_accounting"]
    assert accounting["published_average_point_count"] == 153
    assert accounting["supported_downstream_point_count"] == 54
    assert accounting["excluded_upstream_of_handoff_count"] == 81
    assert accounting["excluded_beyond_common_domain_count"] == 18
    assert accounting["included_x_stations_m"] == [1.25, 1.5, 1.75, 2.0, 2.25, 2.5]
    metrics = record["case_scores"]["published_box"]["all_supported_points"]
    assert metrics["mae_vol_pct_point"] == pytest.approx(2.9656679729)
    assert metrics["rmse_vol_pct_point"] == pytest.approx(3.1545223931)
    assert metrics["pearson_r"] == pytest.approx(0.7604521220)
    four = record["case_scores"]["published_box"]["threshold_classification"]["four_vol_pct"]
    assert four["true_positive"] == 8
    assert four["false_negative"] == 0
    assert four["false_positive"] == 32
    assert four["true_negative"] == 14
    assert record["exact_operator_match"] is False
    assert record["qualified_steady_endpoint_proxy"] is False
    assert record["promotion_allowed"] is False
    assert record["status"] == "resolved_transient_validation_withheld"
    assert arrays["observed_vol_pct"].shape == (54,)


def test_resolved_score_plot_is_written(tmp_path: Path) -> None:
    tool = _tool()
    _, arrays = tool["score_hytunnel_resolved_transient"](
        ROOT,
        ROOT / "outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json",
    )
    destination = tmp_path / "score.png"
    tool["plot_score"](arrays, destination)
    assert destination.stat().st_size > 10_000


def test_stationary_endpoint_replay_is_reproduced_and_scored() -> None:
    tool = _tool()
    replay = (
        ROOT
        / "outputs/hytunnel-resolved-published-box-0p707-final-2026-10-11-v1.json"
    )
    record, _ = tool["score_hytunnel_resolved_transient"](
        ROOT,
        ROOT / "outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json",
        replay,
    )
    assert record["resolved_endpoint_replay"]["maximum_slice_reproduced"] is True
    assert record["qualified_steady_endpoint_proxy"] is True
    assert record["published_vs_expanded_prediction_comparable"] is False
    metrics = record["case_scores"]["published_box"]["all_supported_points"]
    assert metrics["mae_vol_pct_point"] == pytest.approx(2.1792203023)
    assert metrics["rmse_vol_pct_point"] == pytest.approx(2.4905014328)
    subsets = record["case_scores"]["published_box"]["spatial_subsets"]
    assert subsets["centreline_y_zero"]["pair_count"] == 6
    assert subsets["centreline_y_zero"]["mae_vol_pct_point"] == pytest.approx(
        3.1914146621
    )
    by_x = record["case_scores"]["published_box"]["by_x_station"]
    assert by_x["1.75"]["mean_bias_predicted_minus_observed_vol_pct_point"] == pytest.approx(
        2.8909732477
    )
    four = record["case_scores"]["published_box"]["threshold_classification"]["four_vol_pct"]
    assert four["sensitivity"] == 1.0
    assert four["specificity"] == pytest.approx(24.0 / 46.0)
    assert record["promotion_allowed"] is False


def test_changed_physics_candidate_is_scored_without_claiming_frozen_replay() -> None:
    tool = _tool()
    candidate = (
        ROOT / "outputs/hytunnel-resolved-sgs-diffusivity-handoff0p05-"
        "official-0p707-2026-10-11-v1.json"
    )
    record, _ = tool["score_hytunnel_resolved_transient"](
        ROOT,
        ROOT / "outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json",
        candidate,
    )

    evidence = record["resolved_endpoint_replay"]
    assert evidence["frozen_replay_required"] is False
    assert evidence["maximum_slice_reproduced"] is None
    assert evidence["comparison_status"] == (
        "candidate_configuration_not_a_frozen_replay"
    )
    assert record["qualified_steady_endpoint_proxy"] is True
    assert record["published_vs_expanded_prediction_comparable"] is False
    assert set(record["case_scores"]) == {"published_box"}


def test_measured_fan_v2_candidate_is_not_mislabelled_as_frozen_replay() -> None:
    tool = _tool()
    candidate = (
        ROOT / "outputs/hytunnel-measured-fan-sgs-handoff0p05-"
        "official-0p707-2026-10-11-v1.json"
    )
    record, _ = tool["score_hytunnel_resolved_transient"](
        ROOT,
        ROOT / "outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json",
        candidate,
    )

    evidence = record["resolved_endpoint_replay"]
    assert evidence["candidate_fan_inlet_profile"] == "measured_first_station"
    assert evidence["frozen_replay_required"] is False
    assert evidence["maximum_slice_reproduced"] is None
    assert record["case_scores"]["published_box"]["all_supported_points"][
        "pair_count"
    ] == 45
