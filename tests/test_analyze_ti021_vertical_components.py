from __future__ import annotations

import numpy as np
import pytest

from tools.analyze_ti021_vertical_components import _summary


def test_summary_attributes_component_dose_without_sensor_fitting() -> None:
    observed = np.array([10.0, 20.0, 0.0])
    predicted = np.array([12.0, 16.0, 99.0])
    components = np.array([[7.0, 10.0, 50.0], [5.0, 6.0, 49.0]])
    result = _summary(
        np.ones(3, dtype=bool), observed, predicted, components, ["core", "wake"]
    )
    assert result["positive_sensor_count"] == 2
    assert result["mape_percent"] == pytest.approx(20.0)
    assert result["component_dose_sum_vol_percent_s"] == {
        "core": 17.0,
        "wake": 11.0,
    }
    assert result["component_share_of_predicted_dose"] == {
        "core": 17.0 / 28.0,
        "wake": 11.0 / 28.0,
    }
