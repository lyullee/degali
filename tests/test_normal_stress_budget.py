import numpy as np
import pytest
from degali.addons.normal_stress_budget import (normal_work_interval,
    minimum_tke_for_normal_work_cap, covariance_completion, axial_normal_fluxes)


def test_sharp_normal_interval_and_large_budget_stability():
    bound = normal_work_interval(5., 13.)
    assert bound['minimum'] == pytest.approx(1., rel=1e-15)
    assert bound['maximum'] == pytest.approx(25., rel=1e-15)
    assert normal_work_interval(5., 5.)['minimum'] == 5.
    assert normal_work_interval(0., 0.)['maximum'] == 0.
    assert not normal_work_interval(5., 4.)['feasible']
    assert normal_work_interval(1., 1e12)['minimum'] == pytest.approx(5e-13, rel=1e-15)


def test_normal_work_cap_and_inverse_budget():
    assert minimum_tke_for_normal_work_cap(5., 1.) == 13.
    assert minimum_tke_for_normal_work_cap(5., 5.) == 5.
    assert minimum_tke_for_normal_work_cap(5., 20.) == 5.
    assert minimum_tke_for_normal_work_cap(0., .2) == 0.
    for cap in (.01, .03, .08):
        required = minimum_tke_for_normal_work_cap(.11, cap)
        assert normal_work_interval(.11, required)['minimum'] == pytest.approx(cap)


def test_covariance_completion_is_psd_and_has_exact_trace_and_shear():
    r = np.array([[3., 4.], [3., 4.], [0., 0.], [.2, -.4]])
    k, a = np.array([5., 13., 2., 1.]), np.array([5., 1., 0., .8])
    result = covariance_completion(r, k, a)
    cov = result['covariance']
    assert np.min(np.linalg.eigvalsh(cov)) >= -1e-14
    assert np.trace(cov, axis1=1, axis2=2) == pytest.approx(2*k, abs=1e-14)
    np.testing.assert_array_equal(cov[:, 0, 1:], r)
    assert not result['physical_closure']
    with pytest.raises(ValueError, match='no clipping'):
        covariance_completion(np.array([[3., 4.]]), np.array([4.]), np.array([4.]))
    with pytest.raises(ValueError, match='zero axial'):
        covariance_completion(np.array([[3., 4.]]), np.array([5.]), np.array([0.]))


def test_weighted_schur_budget_is_satisfied_by_random_covariances():
    rng = np.random.default_rng(92006)
    factors = rng.normal(size=(100, 3, 3))
    cov = factors @ factors.transpose(0, 2, 1)
    weight = rng.uniform(.01, 3., 100)
    k = .5*np.trace(cov, axis1=1, axis2=2)
    shear = np.linalg.norm(cov[:, 0, 1:], axis=1)
    b, total, normal = weight @ shear, weight @ k, weight @ cov[:, 0, 0]
    bounds = normal_work_interval(b, total)
    assert bounds['minimum'] <= normal <= bounds['maximum']
    assert total >= .5*(normal+b*b/normal)


def test_weighted_budget_can_be_attained_without_fitting():
    r = np.array([[3., 4.], [6., 8.], [.3, .4]])
    s = np.linalg.norm(r, axis=1)
    a = .2*s
    k = .5*(a+s*s/a)
    cov = covariance_completion(r, k, a)['covariance']
    weight = np.array([.5, 2., 4.])
    b, total, normal = weight @ s, weight @ k, weight @ cov[:, 0, 0]
    assert normal_work_interval(b, total)['minimum'] == pytest.approx(normal, rel=1e-14)


def test_axial_momentum_and_energy_work_are_both_retained():
    w, rho, u, a = np.array([.3, .7]), np.array([1., 2.]), np.array([10., 20.]), np.array([.4, .8])
    out = axial_normal_fluxes(.1, w, rho, u, a)
    assert out['axial_momentum_correction'] == pytest.approx(.124)
    assert out['axial_stress_work'] == pytest.approx(2.36)
    assert out['positive_advective_weights'] and not out['coupled_pressure_production_closed']
    negative = axial_normal_fluxes(.1, w, rho, -u, a)
    assert negative['axial_stress_work'] == pytest.approx(-2.36)
    assert not negative['positive_advective_weights']


@pytest.mark.parametrize('b,k', [(-1., 2.), (1., np.nan), (np.inf, 2.)])
def test_bad_budgets_rejected(b, k):
    with pytest.raises(ValueError):
        normal_work_interval(b, k)
    with pytest.raises(ValueError):
        minimum_tke_for_normal_work_cap(1., 0.)
