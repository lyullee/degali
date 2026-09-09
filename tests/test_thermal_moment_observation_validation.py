"""Decision and refinement guards for the thermal-moment observation audit."""

import copy

import pytest

from tools.audit_thermal_moment_observation_validation import (
    directionally_better,
    sensor_refinement_difference,
)


def _aggregate(amplitude, variance):
    return {
        "thermal_deficit": {
            "centre_amplitude": {
                "complete": True, "median_absolute_log_ratio": amplitude,
            },
            "variance": {
                "complete": True, "median_absolute_log_ratio": variance,
            },
        }
    }


def test_directional_decision_requires_no_worse_amplitude_and_variance():
    control = _aggregate(.2, .3)
    assert directionally_better(_aggregate(.2, .29), control)
    assert not directionally_better(_aggregate(.19, .31), control)
    assert not directionally_better(_aggregate(.2, .3), control)


def _rows(delta=0.):
    rows = []
    for trial in (10, 23):
        for station in (1.78, 4.):
            rows.append({
                "trial": trial,
                "station_m": station,
                "projection": {
                    "thermal_deficit_K": [1., 2., 4., 2., 1.],
                    "hydrogen_vol_pct": [1., 2., 4.+delta, 2., 1.],
                },
            })
    return rows


def test_sensor_refinement_uses_frozen_floors_and_population():
    coarse = _rows()
    refined = _rows(.004)
    assert sensor_refinement_difference(coarse, refined) == pytest.approx(.004/4.004)
    missing = copy.deepcopy(refined[:-1])
    with pytest.raises(ValueError, match="do not match"):
        sensor_refinement_difference(coarse, missing)
