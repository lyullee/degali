import numpy as np
import pytest

from degali.addons.integral_obstacle_patch import (
    IntegralObstaclePatchResult,
    all_boundary_trace,
    downstream_trace,
)


def _result():
    steps, nz, ny, nx = 4, 2, 3, 4
    rates = {
        'x_min': np.zeros((steps, nz, ny)),
        'x_max': np.ones((steps, nz, ny)),
        'y_min': np.zeros((steps, nz, nx)),
        'y_max': np.zeros((steps, nz, nx)),
        'z_min': np.zeros((steps, ny, nx)),
        'z_max': np.zeros((steps, ny, nx)),
    }
    return IntegralObstaclePatchResult(
        np.arange(steps + 1, dtype=float), np.arange(nx) + .5,
        np.arange(ny) + .5, np.arange(nz) + .5,
        np.zeros((nz, ny, nx), dtype=bool),
        np.zeros((nz, ny, nx)), np.zeros((nz, ny, nx)),
        np.zeros(steps + 1), np.zeros(steps + 1), np.zeros(steps + 1),
        np.zeros(steps + 1), np.ones((steps, nz, ny)), {},
        {'cell_size_m': 1.0}, rates,
    )


def test_all_face_operator_matches_xmax_route_for_x_wind():
    result = _result()
    point = (5.0, 1.5, .5)
    old = downstream_trace(result, point_m=point, wind_speed_m_s=1.,
                           diffusivity_m2_s=.4 * .205 * .5 / .7)
    new = all_boundary_trace(
        result, point_m=point, wind_xy_m_s=(1., 0.),
        horizontal_diffusivity_m2_s=.4 * .205 * .5 / .7,
    )
    assert new == pytest.approx(old, rel=1.e-10, abs=1.e-12)
