import numpy as np
import pytest

from tools.plot_thorney_ti021_final_field import (
    _projection_max_ignoring_nonfinite,
)


def test_projection_max_ignores_solid_nans_and_preserves_all_solid_cells():
    field = np.asarray(
        [
            [[np.nan, 1.0], [2.0, np.nan]],
            [[np.nan, 3.0], [4.0, np.nan]],
        ]
    )

    top = _projection_max_ignoring_nonfinite(field, axis=0)
    downwind_vertical = _projection_max_ignoring_nonfinite(field, axis=1)
    crosswind_vertical = _projection_max_ignoring_nonfinite(field, axis=2)

    np.testing.assert_allclose(top, [[np.nan, 3.0], [4.0, np.nan]], equal_nan=True)
    np.testing.assert_allclose(
        downwind_vertical, [[2.0, 1.0], [4.0, 3.0]], equal_nan=True
    )
    np.testing.assert_allclose(
        crosswind_vertical, [[1.0, 2.0], [3.0, 4.0]], equal_nan=True
    )
    with pytest.raises(ValueError, match=r"\[z, y, x\]"):
        _projection_max_ignoring_nonfinite(np.ones((2, 2)), axis=0)
