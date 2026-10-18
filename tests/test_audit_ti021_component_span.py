import numpy as np
import pytest

from tools.audit_ti021_component_span import minimum_global_component_mape


def test_global_component_span_uses_one_multiplier_for_every_sensor() -> None:
    component_dose = np.array([
        [0.0, 0.0],
        [1.0, 2.0],
        [2.0, 1.0],
    ])
    observed = np.array([3.0, 3.0])

    mape, multipliers = minimum_global_component_mape(
        component_dose, observed, upper_multiplier=1.0
    )

    assert mape == pytest.approx(0.0, abs=1.0e-10)
    assert multipliers == pytest.approx([1.0, 1.0, 1.0])


def test_global_component_span_reports_unreachable_sensor_shape() -> None:
    component_dose = np.array([
        [0.0, 0.0],
        [1.0, 1.0],
    ])
    observed = np.array([1.0, 2.0])

    mape, multipliers = minimum_global_component_mape(
        component_dose, observed, upper_multiplier=5.0
    )

    assert mape == pytest.approx(25.0)
    assert 1.0 <= multipliers[1] <= 2.0
