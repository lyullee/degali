"""Conservation identities for the Gaussian volumetric-enthalpy candidate."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("CoolProp")
from CoolProp.CoolProp import PropsSI

from degali.addons.axisymmetric_jet import AxisymmetricJetSource, ConservedGaussianJet, phase_ambient_from_rh
from degali.addons.enthalpy_profile import GaussianEnthalpyCrosswind, PhaseMassEnthalpyInverter
from degali.core.jetplume import JetCoefficients


@pytest.fixture(scope="module")
def phase():
    mw, humidity, rho = phase_ambient_from_rh(288.65, 101325., 53.6666667)
    return ConservedGaussianJet(
        AxisymmetricJetSource(diameter=.01, velocity=100., density=1., temperature=100.),
        ambient_temperature=288.65, ambient_pressure=101325., ambient_density=rho,
        fuel_molecular_weight=PropsSI("M", "Hydrogen"), ambient_molecular_weight=mw,
        fuel_heat_capacity=14300., ambient_heat_capacity=1006.,
        ambient_absolute_humidity=humidity, consistent_phase_ambient=True,
        equilibrium_air_condensation=True, temperature_dependent_phase_enthalpy=True,
        conservative_establishment="entrained_mass", radial_points=41,
    )


@pytest.fixture
def candidate(phase):
    jp = SimpleNamespace(
        th=None, k=JetCoefficients(sc=1.16**2), deltay=0., deltaz=0., betay=1., betaz=1., gammaz=0.,
        _split=lambda w, *_: (math.sqrt(w), math.sqrt(w)),
        _wind=lambda *_: 2., _wind_profile=lambda *_: (2., 0.),
    )
    return GaussianEnthalpyCrosswind(jp, phase, quadrature_points=64)


def test_local_phase_inverse_roundtrip_across_cold_and_warm_table(phase):
    inverse = PhaseMassEnthalpyInverter(phase)
    ti, y = np.meshgrid([20., 40., 80., 150., 250., 288.65, 300.], [.001, .02, .2, .5, 1.])
    rho = inverse.a / ti / (1.+inverse.k*y)
    t, h = phase._condensed_air_state(rho, y)
    actual_rho, actual_y, actual_t = inverse.state(rho*y, h)
    assert actual_rho == pytest.approx(rho, rel=2e-8, abs=2e-10)
    assert actual_y == pytest.approx(y, rel=2e-8)
    assert actual_t == pytest.approx(t, rel=2e-8, abs=2e-6)


def test_bilinear_slope_matches_independent_finite_difference(phase):
    inverse = PhaseMassEnthalpyInverter(phase)
    # Off table knots: the piecewise bilinear derivative is one-sided at a
    # knot, whereas a centred finite difference would average two slopes.
    rho, c = np.array([1.18, 1.4, 2.]), np.array([.03, .4, .793])
    _, slope = inverse.enthalpy_and_slope(rho, c)
    step = 1e-6
    forward = phase._condensed_air_state(rho+step, c/(rho+step))[1]
    backward = phase._condensed_air_state(rho-step, c/(rho-step))[1]
    assert slope == pytest.approx((forward-backward)/(2*step), rel=2e-7)
    assert np.all(slope < 0.)


def test_random_table_states_have_no_clipped_inverse_solutions(phase):
    inverse = PhaseMassEnthalpyInverter(phase)
    rng = np.random.default_rng(20260905)
    ti = rng.uniform(14.5, 300., 1024)
    y = rng.uniform(0., 1., 1024)
    rho = inverse.a / ti / (1.+inverse.k*y)
    expected_t, h = phase._condensed_air_state(rho, y)
    found_rho, found_y, found_t = inverse.state(rho*y, h)
    assert found_rho == pytest.approx(rho, rel=2e-8, abs=2e-10)
    assert found_y == pytest.approx(y, rel=2e-8)
    assert found_t == pytest.approx(expected_t, rel=2e-8, abs=2e-6)


@pytest.mark.parametrize("c,h", [(-1., 0.), (math.nan, 0.), (.1, math.inf), (10., -1e6), (.1, 1e9)])
def test_inverse_rejects_invalid_or_unbracketed_states(phase, c, h):
    with pytest.raises(ValueError):
        PhaseMassEnthalpyInverter(phase).state(c, h)


def test_ambient_fixed_point_and_finite_dilution_limit(candidate):
    state = np.array([1.13, .15, .02, .01, 25., 1., 1.5])
    g = np.array([1., .5, .01, 1e-4, 1e-8, 0.])
    rho, y, t, h = candidate.thermodynamic_profile(state, g)
    assert rho[-1] == pytest.approx(candidate.rhoa, abs=1e-11)
    assert t[-1] == pytest.approx(288.65, abs=1e-8)
    assert y[-1] == h[-1] == 0.
    assert rho[0] == pytest.approx(state[0], rel=1e-10)
    assert y[0] == pytest.approx(state[1], rel=1e-10)
    assert h[:-1]/(rho[:-1]*y[:-1]) == pytest.approx(np.full(5, h[0]/(state[0]*state[1])), rel=1e-10)


def test_profile_energy_and_species_integrals_obey_same_scalar_shape(candidate):
    state = np.array([1.13, .15, .02, .01, 25., 1., 1.5])
    u, rho, y, t, h = candidate.profiles(state)
    exponent, weights = candidate._quadrature(candidate.quadrature_points)
    cc = state[0]*state[1]
    hc, _ = candidate.phase_inverse.enthalpy_and_slope(state[0], cc)
    assert rho*y == pytest.approx(cc*np.exp(-exponent), rel=1e-12)
    assert h == pytest.approx(hc*np.exp(-exponent), rel=1e-12)
    flux = candidate.integral_fluxes(state)
    kinetic = .5*state[2]*np.sum(rho*u**3*weights)
    assert flux.energy-kinetic == pytest.approx(hc/cc*flux.contaminant_mass, rel=1e-12)
    true_h = candidate.thermodynamics._condensed_air_state(rho, y)[1]
    assert true_h == pytest.approx(h, rel=1e-8, abs=1e-8)
    assert candidate.buoyancy_force(state) == pytest.approx(9.81*state[2]*np.sum((candidate.rhoa-rho)*weights))


def test_five_flux_projection_closes_without_density_gaussian_shortcut(candidate):
    state = np.array([1.13, .15, .02, .01, 25., 1., 1.5])
    target = candidate.integral_fluxes(state)
    guess = state.copy()
    guess[:5] *= [1.002, .998, 1.002, 1., .998]
    result = candidate.project(target, guess)
    assert result.success
    assert max(result.relative_residuals.values()) < 1e-8
    assert result.quadrature_residual < 1e-5


def test_dedicated_flux_inverse_closes_all_five_balances(candidate):
    state = np.array([1.13, .15, .02, .01, 25., 1., 1.5])
    target = candidate._as_array(candidate.integral_fluxes(state))
    guess = state.copy()
    guess[:5] *= [1.002, .998, 1.002, 1., .998]
    found = candidate._match_flux_array(target, guess)
    actual = candidate._as_array(candidate.integral_fluxes(found))
    assert actual == pytest.approx(target, rel=1e-9, abs=1e-9)


def test_short_march_is_conservative_and_step_refines(candidate):
    initial = np.array([1.13, .15, .02, .01, 25., 1., 1.5])
    coarse = candidate.solve(initial, maximum_distance=.02, maximum_step=.01)
    fine = candidate.solve(initial, maximum_distance=.02, maximum_step=.005)
    assert coarse.maximum_relative_balance_residual < 1e-8
    assert fine.maximum_relative_balance_residual < 1e-8
    assert fine.states[-1] == pytest.approx(coarse.states[-1], rel=1e-5, abs=1e-7)
    assert fine.states[-1, 1] < initial[1]
    assert np.all(fine.states[:, [0, 1, 2, 4]] > 0.)


def test_receptor_uses_the_same_thermal_profile_as_fluxes(candidate):
    state = np.array([1.13, .15, .02, .01, 25., 1., 1.5])
    y, z = .05, 1.6
    sy, sz = candidate.section_widths(state)
    shape = math.exp(-.5*(y/sy)**2) * (
        math.exp(-.5*((z-state[-1])/sz)**2) + math.exp(-.5*((z+state[-1])/sz)**2)
    )
    rho, fraction, t, h = candidate.thermodynamic_profile(state, np.array([shape]))
    assert candidate.point_temperature(state, y, z) == t[0]
    mw_h, mw_air = candidate.thermodynamics.fuel_molecular_weight, candidate.thermodynamics._humid_ambient_molecular_weight
    expected_x = (fraction[0]/mw_h) / (fraction[0]/mw_h + (1.-fraction[0])/mw_air)
    assert candidate.point_mole_fraction(state, y, z) == pytest.approx(expected_x)


@pytest.mark.parametrize("order", [16, 32, 64, 128])
def test_polar_square_quadrature_preserves_domain_and_gaussian_integrals(candidate, order):
    q, weight = candidate._quadrature(order)
    assert len(q) == 2*order
    assert np.all(weight > 0.)
    assert sum(weight) == pytest.approx(math.pi*candidate.k.delta**2, rel=2e-14)
    for power in (0., .1, 1., 1.16**2, 2., 3.*1.16**2):
        # Use full-precision pi, as the original tensor quadrature does;
        # legacy profile_integral retains a rounded Fortran pi constant.
        analytic = (math.pi*candidate.k.delta**2 if power == 0. else
                    2.*math.pi/power*math.erf(candidate.k.delta*math.sqrt(math.pi*power/8.))**2)
        assert np.sum(weight*np.exp(-power*q)) == pytest.approx(analytic, rel=3e-14)


def test_quadrature_gate_catches_mass_error_even_when_energy_matches(candidate):
    from degali.core.jetplume import JetIntegralFluxes
    fine = JetIntegralFluxes(total_mass=2., contaminant_mass=.1, momentum_x=3., momentum_z=.1, energy=-1e6)
    coarse = JetIntegralFluxes(total_mass=2.001, contaminant_mass=.1, momentum_x=3., momentum_z=.1, energy=-1e6)
    assert candidate.quadrature_error(coarse, fine) == pytest.approx(.0005)
