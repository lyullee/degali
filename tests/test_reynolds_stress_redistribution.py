import numpy as np
import pytest

from degali.addons.reynolds_stress_redistribution import (
    maximum_psd_step, slow_pressure_strain,
)


def test_slow_pressure_strain_is_trace_free_and_returns_to_isotropy():
    covariance = np.diag([4., 1., 1.])
    out = slow_pressure_strain(covariance, 3., 2., 1.8)
    assert out['trace'] == pytest.approx(0., abs=1e-14)
    assert out['tke_change_rate'] == 0.
    assert out['pressure_strain'][0, 0] < 0.
    assert out['pressure_strain'][1, 1] > 0.
    assert out['closure_inputs_explicit'] and not out['adopted']


def test_redistribution_rejects_trace_inconsistent_tke():
    with pytest.raises(ValueError, match='trace'):
        slow_pressure_strain(np.eye(3), 2., 1., 1.8)


def test_psd_step_limit_is_exact_for_a_diagonal_decay_rate():
    out = maximum_psd_step(np.diag([4., 2., 1.]), np.diag([-2., .5, 0.]))
    assert out['maximum_step'] == pytest.approx(2.)
    assert out['limiting_normalized_eigenvalue'] == pytest.approx(-.5)
    assert out['covariance_projected'] is False


def test_psd_step_rejects_semidefinite_initial_state_and_nonsymmetric_rate():
    with pytest.raises(ValueError, match='positive definite'):
        maximum_psd_step(np.diag([1., 1., 0.]), np.eye(3))
    with pytest.raises(ValueError, match='symmetric'):
        maximum_psd_step(np.eye(3), [[0., 1., 0.], [0., 0., 0.], [0., 0., 0.]])
