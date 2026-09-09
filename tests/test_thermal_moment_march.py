import math

import numpy as np
import pytest

from degali.addons.thermal_moment_march import GuardedThermalMomentMarch


def evaluator(parameters):
    rate = np.zeros(8)
    rate[0] = parameters[0]
    return {
        "rates": rate,
        "ledger": {"sources": np.r_[parameters[0], np.zeros(4)]},
        "moment_rate": 2.0 * parameters[0],
        "weak_budget_scaled_error": 1.0e-10,
        "minimum_chi_species": 0.2,
        "minimum_chi_momentum": 0.3,
        "edge_gradient_defects": {"heat": 0.1},
        "valid": True,
    }


def test_adaptive_midpoint_marches_state_and_cumulative_budgets():
    march = GuardedThermalMomentMarch(evaluator).march(
        np.r_[1.0, np.zeros(7)],
        0.1,
        initial_step=0.02,
        maximum_step=0.02,
        minimum_step=1.0e-6,
        tolerance=1.0e-5,
    )
    assert march.stop_reason == "arc-length ceiling reached"
    assert march.states[-1, 0] == pytest.approx(math.exp(0.1), rel=2.0e-5)
    assert march.states[-1, 8] == pytest.approx(march.states[-1, 0] - 1.0, rel=2.0e-5)
    assert march.states[-1, 13] == pytest.approx(2.0 * (march.states[-1, 0] - 1.0), rel=2.0e-5)
    assert march.minimum_sampled_diffusivity == pytest.approx(0.2)
    assert march.maximum_weak_residual == pytest.approx(1.0e-10)


def test_target_predicate_stops_after_first_accepted_crossing():
    result = GuardedThermalMomentMarch(evaluator).march(
        np.r_[1.0, np.zeros(7)],
        1.0,
        initial_step=0.01,
        maximum_step=0.01,
        target=lambda parameters: parameters[0] >= math.exp(0.05),
    )
    assert result.reached_target
    assert result.stop_reason == "target reached"
    assert result.arc_length[-1] == pytest.approx(0.05, abs=0.011)


def test_invalid_local_closure_is_not_clipped_or_extrapolated():
    def invalid(parameters):
        if parameters[0] > 1.005:
            raise ValueError("counter-gradient mixing")
        return evaluator(parameters)

    result = GuardedThermalMomentMarch(invalid).march(
        np.r_[1.0, np.zeros(7)],
        0.1,
        initial_step=0.01,
        maximum_step=0.01,
        minimum_step=0.001,
    )
    assert not result.reached_target
    assert "minimum step rejected" in result.stop_reason
    assert "counter-gradient mixing" in result.stop_reason
    assert result.arc_length[-1] < 0.01
    assert result.rejected_steps > 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"length": 0.0},
        {"length": 1.0, "initial_step": 0.1, "maximum_step": 0.01},
        {"length": 1.0, "maximum_steps": 0},
    ],
)
def test_invalid_march_options_are_rejected(kwargs):
    with pytest.raises(ValueError):
        GuardedThermalMomentMarch(evaluator).march(np.ones(8), **kwargs)


def test_incomplete_evaluator_result_is_rejected_at_minimum_step():
    result = GuardedThermalMomentMarch(lambda parameters: {"valid": True}).march(
        np.ones(8),
        0.01,
        initial_step=0.001,
        maximum_step=0.001,
        minimum_step=0.001,
    )
    assert "incomplete result" in result.stop_reason
