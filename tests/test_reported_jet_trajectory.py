import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.validation.reported_jet_trajectory import ReportedJetTrajectory


class LinearTable:
    def from_concentration(self, c):
        return SimpleNamespace(yc=c/100, temp=300-c)


def reported():
    rows = np.ones((3, 12))
    rows[:, 0] = [1, 2, 3]
    rows[:, 1] = [0.5, 0.8, 1.2]
    rows[:, 3:5] = [[0.6, 0.9], [0.7, 1.1], [0.8, 1.4]]
    rows[:, 2] = np.array([2, 4, 8])*(1+np.exp(-2*(rows[:, 1]/rows[:, 4])**2))
    return rows


def test_reported_centreline_recovered_and_input_unchanged():
    rows = reported()
    saved = rows.copy()
    trajectory = ReportedJetTrajectory(LinearTable(), rows)
    for row in rows:
        assert trajectory.concentration_at(row[0], 0, row[1]) == pytest.approx(row[2])
        assert trajectory.temperature_at(row[0], 0, row[1]) == pytest.approx(300-row[2])
    np.testing.assert_array_equal(rows, saved)
    np.testing.assert_allclose(trajectory.rows[:, 2], [2, 4, 8])


def test_receptors_and_interpolation_use_bare_node_values():
    model = ReportedJetTrajectory(LinearTable(), reported())
    x, y, z = 1.5, 0.4, 0.3
    cc, height, sy, sz = 3, 0.65, 0.65, 1.0
    expected = cc*math.exp(-0.5*(y/sy)**2)*(math.exp(-0.5*((z-height)/sz)**2)+math.exp(-0.5*((z+height)/sz)**2))
    assert model.concentration_at(x, y, z) == pytest.approx(expected)
    assert model.concentration_at(x, -y, z) == pytest.approx(expected)


def test_ground_symmetry_has_zero_normal_slope():
    model = ReportedJetTrajectory(LinearTable(), reported())
    eps = 1e-5
    slope = (model.concentration_at(2, 0, eps)-model.concentration_at(2, 0, 0))/eps
    assert abs(slope) < 2e-5


@pytest.mark.parametrize('change', ['touchdown', 'negative_width', 'nan', 'duplicate_x', 'negative_c'])
def test_invalid_reported_rows_rejected(change):
    rows = reported()
    column, value = {'touchdown': (1, 0), 'negative_width': (4, -1), 'nan': (2, np.nan),
                     'duplicate_x': (0, 2), 'negative_c': (2, -1)}[change]
    rows[0, column] = value
    with pytest.raises(ValueError):
        ReportedJetTrajectory(LinearTable(), rows)


@pytest.mark.parametrize('q', [(0, 0, 1), (4, 0, 1), (2, 0, -1), (2, np.nan, 1)])
def test_invalid_receptors_never_use_ambient_fallback(q):
    model = ReportedJetTrajectory(LinearTable(), reported())
    with pytest.raises(ValueError):
        model.concentration_at(*q)
    with pytest.raises(ValueError):
        model.temperature_at(*q)


def test_wrong_row_schema_rejected():
    with pytest.raises(ValueError):
        ReportedJetTrajectory(LinearTable(), reported()[:, :5])
