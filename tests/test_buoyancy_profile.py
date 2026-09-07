"""Independent thermal width has a sixth constraint, never a fake transport law."""

import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection


@pytest.fixture
def section(candidate):
    return BuoyancyConstrainedEnthalpySection(candidate.jetplume, candidate.thermodynamics,
                                             quadrature_points=256)


STATE = np.array([1.13, .15, .02, .01, 25., 1., 1.5])


def test_unit_thermal_width_reproduces_previous_candidate(section, candidate):
    shape = np.array([1., .5, .01, 1e-6, 0.])
    for old, new in zip(candidate.thermodynamic_profile(STATE, shape), section.thermodynamic_profile(STATE, shape)):
        assert new == pytest.approx(old, rel=1e-12, abs=1e-10)
    old = candidate._as_array(candidate.integral_fluxes(STATE, quadrature_points=256))
    assert section.moments(STATE)[:5] == pytest.approx(old, rel=1e-12)


def test_distinct_gaussians_obey_independent_analytic_enthalpy_integral(section):
    section.thermal_width_ratio = 1.08
    moment = section.moments(STATE)
    u, rho, y, _, h = section.profiles(STATE)
    q, w = section._quadrature(section.quadrature_points)
    hc, _ = section.phase_inverse.enthalpy_and_slope(STATE[0], STATE[0]*STATE[1])
    exponent = 1./1.08**2
    assert h == pytest.approx(hc*np.exp(-exponent*q), rel=1e-12)
    analytic = hc*STATE[2]*(2.*np.cos(STATE[3])*section.k.profile_integral(exponent)
               + STATE[4]*section.k.profile_integral(exponent+1.16**2))
    kinetic = .5*STATE[2]*np.sum(w*rho*u**3)
    assert moment[4]-kinetic == pytest.approx(analytic, rel=1e-10)
    assert moment[5] == pytest.approx(section.buoyancy_force(STATE), rel=1e-12)


def test_synthetic_six_moment_boundary_closes(section):
    known = STATE.copy()
    known[:2] = 1.25, .1537  # Away from interpolation knots and neutral force.
    section.thermal_width_ratio = 1.08
    target = section.moments(known, quadrature_points=1024)
    section.thermal_width_ratio = 1.
    result = section.project_buoyancy(target, known)
    if not result.success:
        pytest.fail(f"{result.message}: flux={result.relative_residuals.tolist()}, "
                    f"quadrature={result.quadrature_residuals.tolist()}, beta={result.thermal_width_ratio}")
    assert max(result.relative_residuals) < 1e-8
    assert max(result.quadrature_residuals) < 1e-5
    assert result.thermal_width_ratio == pytest.approx(1.08, rel=1e-4)
    assert result.state == pytest.approx(known, rel=1e-4, abs=1e-6)


def test_near_neutral_force_is_not_approved_from_flux_accuracy_alone(section):
    # This case reconstructs the state, but force cancellation makes its
    # 1024/2048 quadrature difference fail the separately frozen force gate.
    section.thermal_width_ratio = 1.08
    target = section.moments(STATE, quadrature_points=1024)
    result = section.project_buoyancy(target, STATE)
    assert max(result.relative_residuals) < 1e-8
    assert result.quadrature_residuals[5] > 1e-5
    assert not result.success


def test_phase_partial_derivative_at_fixed_density(section):
    rho, c = np.array([1.18, 1.4, 2.]), np.array([.03, .4, .793])
    _, partial = section._phase_partials(rho, c)
    step = 1e-7
    forward = section.thermodynamics._condensed_air_state(rho, (c+step)/rho)[1]
    backward = section.thermodynamics._condensed_air_state(rho, (c-step)/rho)[1]
    assert partial == pytest.approx((forward-backward)/(2*step), rel=2e-7)


def test_reduced_moment_jacobian_matches_fixed_species_difference(section):
    # A centred difference straddles the table slope jump at Y=.15; compare
    # differentiable states, leaving the knot case to the solve/gate test.
    state = STATE.copy()
    state[:2] = 1.1312, .1537
    section.thermal_width_ratio = 1.04
    analytic = section.reduced_moment_jacobian(state)
    target_h2 = section.moments(state)[1]
    params = np.log([state[0], state[2], state[4], 1.04])
    def evaluate(p):
        rho, area, uc, beta = np.exp(p)
        section.thermal_width_ratio = beta
        trial = state.copy()
        trial[[0, 2, 4]] = rho, area, uc
        denominator = area*rho*(section._wind(trial)*np.cos(trial[3])*section.k.profile_integral(1.)
                      + uc*section.k.profile_integral(1.+section.velocity_shape_exponent))
        trial[1] = target_h2/denominator
        moments = section.moments(trial)
        return np.array([moments[0], np.hypot(moments[2], moments[3]), moments[4], moments[5]])
    step = 1e-7
    numerical = np.column_stack([(evaluate(params+np.eye(4)[i]*step)-evaluate(params-np.eye(4)[i]*step))/(2*step)
                                 for i in range(4)])
    assert analytic == pytest.approx(numerical, rel=5e-5, abs=1e-6)


@pytest.mark.parametrize("value", [0., -.1, np.nan, np.inf])
def test_invalid_thermal_width_is_rejected(section, value):
    with pytest.raises(ValueError, match="thermal width"):
        section.thermal_width_ratio = value


@pytest.mark.parametrize("method", ["solve", "derivatives", "source_terms", "_match_flux_array", "_receptor", "project"])
def test_unclosed_transport_and_five_constraint_shortcut_are_blocked(section, method):
    with pytest.raises(NotImplementedError):
        getattr(section, method)(STATE)


def test_invalid_six_moment_targets_are_rejected(section):
    with pytest.raises(ValueError, match="six finite moments"):
        section.project_buoyancy(np.ones(5), STATE)
    with pytest.raises(ValueError, match="two orders"):
        section.project_buoyancy(section.moments(STATE), STATE, maximum_quadrature_points=256)
