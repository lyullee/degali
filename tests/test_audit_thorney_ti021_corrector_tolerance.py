import numpy as np
import pytest

from tools.audit_thorney_ti021_corrector_tolerance import _difference


def test_tolerance_difference_retains_absolute_and_relative_scales():
    reference = np.asarray((1.0, 2.0, np.nan, 3.0))
    candidate = np.asarray((1.1, 1.8, np.nan, 3.0))

    result = _difference(reference, candidate)

    assert result["maximum_absolute_difference"] == pytest.approx(0.2)
    assert result["root_mean_square_difference"] == pytest.approx(
        np.sqrt((0.1**2 + 0.2**2) / 3.0)
    )
    assert result["relative_l1_difference"] == pytest.approx(0.3 / 6.0)
