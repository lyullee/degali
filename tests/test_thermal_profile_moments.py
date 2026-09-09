import numpy as np
import pytest

from degali.validation.thermal_profile_moments import (
    best_integer_lag,
    project_model_profile,
    summarize_paired_profile_moments,
    truncated_profile_moment,
)


def test_truncated_profile_moment_is_symmetric_and_scale_invariant():
    z = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    value = np.array([1.0, 3.0, 5.0, 3.0, 1.0])
    one = truncated_profile_moment(z, value)
    two = truncated_profile_moment(z, 17.0 * value)
    assert one.centre == pytest.approx(0.0, abs=1e-15)
    assert one.variance == pytest.approx(two.variance, abs=1e-15)
    assert two.zeroth == pytest.approx(17.0 * one.zeroth)


@pytest.mark.parametrize(
    "coordinate,value,message",
    [
        ([0.0, 1.0], [1.0, 1.0], "at least 3"),
        ([0.0, 1.0, 0.5], [1.0, 1.0, 1.0], "increasing"),
        ([0.0, 1.0, 2.0], [1.0, -1.0, 1.0], "nonnegative"),
        ([0.0, 1.0, 2.0], [0.0, 0.0, 0.0], "positive"),
    ],
)
def test_truncated_profile_moment_rejects_invalid_inputs(coordinate, value, message):
    with pytest.raises(ValueError, match=message):
        truncated_profile_moment(coordinate, value)


def test_best_integer_lag_recovers_delayed_signal():
    x = np.sin(np.linspace(0.0, 4.0 * np.pi, 80)) + np.linspace(0.0, 0.4, 80)
    delayed = np.r_[np.zeros(3), x[:-3]]
    result = best_integer_lag(x, delayed, maximum_seconds=5)
    assert result.seconds == 3
    assert result.correlation > 0.99


def test_summary_keeps_negative_growth_and_requires_common_delay():
    time = np.arange(20.0)
    z = np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    pulse = 20.0 + 5.0 * np.sin(time / 3.0)
    narrow = np.exp(-0.5 * (z / 0.16) ** 2)
    broad = np.exp(-0.5 * (z / 0.25) ** 2)
    thermal_up = pulse[:, None] * narrow
    thermal_down = np.r_[np.zeros((2, 5)), (pulse[:-2, None] * broad)]
    hydrogen_up = 0.2 * pulse[:, None] * broad
    hydrogen_down = np.r_[np.zeros((2, 5)), (0.2 * pulse[:-2, None] * narrow)]
    data = dict(
        coordinate=z,
        time_seconds=time,
        profiles={
            "1.78": {"thermal_deficit": thermal_up, "hydrogen": hydrogen_up},
            "4.0": {"thermal_deficit": thermal_down, "hydrogen": hydrogen_down},
        },
    )
    result = summarize_paired_profile_moments(data)
    assert result["common_lag_seconds"] == 2
    assert result["common_delay_passed"]
    assert result["paired_profiles"] >= 10
    assert result["growth_summary"]["thermal_positive"] > 0
    assert result["growth_summary"]["species_nonpositive"] > 0
    assert result["growth_summary"]["median_ratio"] < 0.0


class _FakeTrajectory:
    class model:
        class thermodynamics:
            ambient_temperature = 300.0

    def __init__(self, *, covered=True, invalid_hydrogen=False):
        self.covered = covered
        self.invalid_hydrogen = invalid_hydrogen
        self.sampled_heights = []

    def state_at(self, station):
        return np.ones(7) if self.covered and station == 4.0 else None

    def temperature_at(self, station, lateral, height):
        self.sampled_heights.append(height)
        offset = height - 1.5
        return 300.0 - 20.0 * np.exp(-0.5 * (offset / 0.25) ** 2)

    def concentration_at(self, station, lateral, height):
        if self.invalid_hydrogen:
            return 101.0
        offset = height - 1.5
        return 10.0 * np.exp(-0.5 * (offset / 0.20) ** 2)


def test_project_model_profile_uses_absolute_sensor_heights_and_direct_moments():
    trajectory = _FakeTrajectory()
    result = project_model_profile(trajectory, 4.0, 1.5)
    assert result.absolute_height == pytest.approx([1.0, 1.25, 1.5, 1.75, 2.0])
    assert trajectory.sampled_heights == pytest.approx(result.absolute_height)
    assert result.thermal_deficit[2] == pytest.approx(20.0)
    assert result.hydrogen[2] == pytest.approx(10.0)
    assert result.thermal_moment.variance > result.hydrogen_moment.variance


def test_project_model_profile_rejects_missing_or_nonphysical_receptor():
    with pytest.raises(ValueError, match="does not cover"):
        project_model_profile(_FakeTrajectory(covered=False), 4.0, 1.5)
    with pytest.raises(ValueError, match="physical range"):
        project_model_profile(_FakeTrajectory(invalid_hydrogen=True), 4.0, 1.5)
    with pytest.raises(ValueError, match="below ground"):
        project_model_profile(_FakeTrajectory(), 4.0, 0.25)
