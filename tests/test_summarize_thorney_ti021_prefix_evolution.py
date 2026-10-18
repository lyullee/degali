import json

import numpy as np

from tools.summarize_thorney_ti021_prefix_evolution import summarize


def _write_prefix(root, stop_time, scale):
    root.mkdir()
    solid = np.zeros((2, 2, 2), dtype=bool)
    solid[0, 0, 0] = True
    peak = np.full((2, 2, 2), 2.0 * scale)
    final = np.full((2, 2, 2), 3.0 * scale)
    peak[solid] = np.nan
    final[solid] = np.nan
    np.savez_compressed(
        root / "ti021-inner-final-field.npz",
        solid_mask=solid,
        maximum_concentration_vol_percent=peak,
        concentration_vol_percent=final,
    )
    labels = np.asarray(["S1", "S2"])
    traces = np.full((2, int(stop_time * 10) + 1), scale)
    np.savez_compressed(
        root / "ti021-inner-receptor-traces.npz",
        time_s=np.linspace(0.0, stop_time, traces.shape[1]),
        labels=labels,
        concentration_vol_percent=traces,
    )
    (root / "ti021-inner-transient.json").write_text(
        json.dumps({
            "stop_time_s": stop_time,
            "numerics": {"corrector_nonconverged_steps": 0},
        }),
        encoding="utf-8",
    )


def test_prefix_evolution_is_monotonic_and_hash_pinned(tmp_path):
    earlier = tmp_path / "earlier"
    later = tmp_path / "later"
    _write_prefix(earlier, 1.0, 1.0)
    _write_prefix(later, 2.0, 2.0)
    record = summarize(earlier, later)
    assert record["status"] == "accepted_computed_window_evolution"
    assert record["delta_later_minus_earlier"]["stop_time_s"] == 1.0
    assert record["later"]["receptor_global_peak_vol_percent"] == 2.0
    assert record["earlier"]["record_sha256"] != record["later"]["record_sha256"]
    assert record["promotion_allowed"] is False
