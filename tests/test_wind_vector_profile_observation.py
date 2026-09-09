from types import SimpleNamespace

import pytest

from audit_wind_vector_profile_observation import (
    YawedPercentTrajectory,
    displacement_gate,
    group_errors,
    step_convergence,
)


class FakeYawedTrajectory:
    def __init__(self):
        thermodynamics = SimpleNamespace(ambient_temperature=289.0)
        self.model = SimpleNamespace(base=SimpleNamespace(thermodynamics=thermodynamics))

    def section(self, station, lateral):
        if station < 0.0:
            raise ValueError("outside trajectory")
        return "state", lateral

    def temperature_at(self, station, lateral, height):
        return station + lateral + height

    def concentration_at(self, station, lateral, height):
        return 0.123


def test_yawed_adapter_exposes_common_percent_api():
    trajectory = YawedPercentTrajectory(FakeYawedTrajectory())
    assert trajectory.model.thermodynamics.ambient_temperature == 289.0
    assert trajectory.state_at(1.0) == "state"
    assert trajectory.state_at(-1.0) is None
    assert trajectory.temperature_at(1.0, 0.2, 0.5) == pytest.approx(1.7)
    assert trajectory.concentration_at(1.0, 0.2, 0.5) == pytest.approx(12.3)


def comparison_row(value=0.2, within=False):
    comparison = {}
    for scalar in ("thermal_deficit", "hydrogen"):
        comparison[scalar] = {}
        for field in ("centre_amplitude", "zeroth_moment", "variance"):
            comparison[scalar][field] = {
                "absolute_log_ratio": value,
                "within_observed_iqr": within,
            }
    return {"comparison": comparison}


def test_group_errors_keeps_scalar_diagnostics_separate():
    first = comparison_row(0.1, True)
    second = comparison_row(0.3, False)
    summary = group_errors([first, second])
    assert summary["thermal_deficit"]["centre_amplitude"] == {
        "profiles": 2,
        "median_absolute_log_ratio": pytest.approx(0.2),
        "maximum_absolute_log_ratio": pytest.approx(0.3),
        "within_observed_iqr": 1,
    }
    assert summary["hydrogen"]["zeroth_moment"]["profiles"] == 2


def summary(value):
    return {
        scalar: {
            field: {"median_absolute_log_ratio": value}
            for field in ("centre_amplitude", "zeroth_moment", "variance")
        }
        for scalar in ("thermal_deficit", "hydrogen")
    }


def test_displacement_gate_requires_each_amplitude_and_variance_condition():
    control = summary(0.3)
    candidate = summary(0.2)
    assert displacement_gate(control, candidate)["passed"]

    candidate["hydrogen"]["variance"]["median_absolute_log_ratio"] = 0.31
    failed = displacement_gate(control, candidate)
    assert not failed["passed"]
    assert failed["amplitude_passed"]
    assert not failed["variance_passed"]


def projection_row(station, multiplier=1.0):
    moment = {
        "zeroth": 2.0 * multiplier,
        "centre": 0.1 * multiplier,
        "variance": 0.4 * multiplier,
    }
    return {
        "station_m": station,
        "projection": {
            "thermal_deficit_K": [1.0 * multiplier, 2.0 * multiplier],
            "hydrogen_vol_pct": [3.0 * multiplier, 4.0 * multiplier],
            "thermal_moment": moment,
            "hydrogen_moment": moment,
        },
    }


def test_step_convergence_applies_frozen_relative_limit():
    reference = [projection_row(1.78), projection_row(4.0)]
    close = [projection_row(1.78, 1.004), projection_row(4.0, 1.004)]
    far = [projection_row(1.78, 1.006), projection_row(4.0, 1.006)]
    assert step_convergence(reference, close)["passed"]
    assert not step_convergence(reference, far)["passed"]
    with pytest.raises(ValueError, match="matching stations"):
        step_convergence(reference, [projection_row(1.78)])
