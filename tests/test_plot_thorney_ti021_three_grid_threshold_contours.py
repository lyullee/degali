import numpy as np
import pytest

from tools.plot_thorney_ti021_three_grid_threshold_contours import (
    _interpolate_cell_center_slice,
)


def test_shared_plane_interpolation_is_linear_and_does_not_cross_solids():
    values = np.asarray((
        ((0.0, 2.0), (4.0, 6.0)),
        ((10.0, 12.0), (14.0, 16.0)),
        ((20.0, 22.0), (24.0, 26.0)),
    ))
    solid = np.zeros_like(values, dtype=bool)
    solid[0, 0, 1] = True
    solid[1, 1, 0] = True
    solid[0, 1, 1] = True
    solid[1, 1, 1] = True

    plane = _interpolate_cell_center_slice(
        values, solid, np.asarray((0.5, 1.5, 2.5)), axis=0, target=1.0,
    )

    assert plane[0, 0] == pytest.approx(5.0)
    assert plane[0, 1] == pytest.approx(12.0)  # use adjacent fluid cell only
    assert plane[1, 0] == pytest.approx(4.0)   # use adjacent fluid cell only
    assert np.ma.is_masked(plane[1, 1])


def test_shared_plane_rejects_out_of_domain_targets():
    values = np.zeros((2, 2, 2))
    solid = np.zeros_like(values, dtype=bool)
    with pytest.raises(ValueError, match="outside cell-centre range"):
        _interpolate_cell_center_slice(
            values, solid, np.asarray((0.5, 1.5)), axis=0, target=2.0,
        )
