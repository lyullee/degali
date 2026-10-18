from pathlib import Path
import runpy

import numpy as np
import pytest


ROOT = Path("reference/datasets/usn-carpark-14405903/raw")


@pytest.mark.slow
@pytest.mark.skipif(
    not (ROOT / "Exp19.mat").is_file(),
    reason="requires the locally controlled public USN car-park data",
)
def test_usn_enclosure_parameter_holdout_is_bounded_and_reproducible():
    tool = runpy.run_path("tools/validate_usn_enclosure_dynamic.py")
    cases = {
        number: tool["_case_data"](ROOT / f"Exp{number}.mat")
        for number in tool["ALL_CASES"]
    }
    parameters, diagnostics = tool["_fit"](cases[tool["TRAINING_CASE"]])

    assert parameters == pytest.approx((0.547503, 0.586732, 5.1), abs=2.0e-4)
    assert diagnostics["delay_grid_count"] == 101
    arrival_errors = []
    duration_errors = []
    for number in tool["HOLDOUT_CASES"]:
        result = tool["_predict"](cases[number], parameters)
        observed = cases[number]["observed_vol_percent"]
        score = tool["_score"](result.concentration_vol_percent, observed)
        assert score["mae_vol_percent_points"] < 0.35
        assert score["pearson_r"] > 0.97
        assert np.all(result.concentration_vol_percent >= 0.0)
        assert np.all(result.concentration_vol_percent <= 100.0)
        assert result.maximum_mass_residual_kg < 1.0e-12
        measured_threshold = tool["_threshold"](observed, cases[number]["time_s"])
        predicted_threshold = tool["_threshold"](
            result.concentration_vol_percent, cases[number]["time_s"]
        )
        arrival_errors.append(
            predicted_threshold["arrival_s"] - measured_threshold["arrival_s"]
        )
        duration_errors.append(
            predicted_threshold["duration_s"] - measured_threshold["duration_s"]
        )
    assert np.mean(np.abs(arrival_errors)) == pytest.approx(3.7, abs=0.01)
    assert np.mean(np.abs(duration_errors)) == pytest.approx(27.1, abs=0.01)
