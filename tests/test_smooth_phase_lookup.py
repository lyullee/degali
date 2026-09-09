"""C1 phase lookup identities and inverse consistency."""

import numpy as np
import pytest

from test_enthalpy_profile import candidate, phase
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.enthalpy_profile import PhaseMassEnthalpyInverter
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.smooth_phase_lookup import (
    HermitePhaseSurface,
    SmoothPhaseMassEnthalpyInverter,
)


def test_hermite_surface_reproduces_nodes_and_has_continuous_first_derivatives():
    x = np.linspace(1., 4., 7)
    y = np.linspace(.1, .9, 9)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    values = np.log(xx)+yy*yy+.2*xx*yy
    surface = HermitePhaseSurface(x, y, values)
    actual = surface.evaluate(xx, yy)[0]
    assert actual == pytest.approx(values, abs=2e-15)
    epsilon = 1e-9
    left = surface.evaluate(np.full(5, x[3]-epsilon), np.linspace(.2, .8, 5))
    right = surface.evaluate(np.full(5, x[3]+epsilon), np.linspace(.2, .8, 5))
    assert left[0] == pytest.approx(right[0], abs=2e-8)
    assert left[1] == pytest.approx(right[1], abs=2e-8)
    assert left[2] == pytest.approx(right[2], abs=2e-8)
    below = surface.evaluate(np.linspace(1.2, 3.8, 5), np.full(5, y[4]-epsilon))
    above = surface.evaluate(np.linspace(1.2, 3.8, 5), np.full(5, y[4]+epsilon))
    assert below[0] == pytest.approx(above[0], abs=2e-8)
    assert below[1] == pytest.approx(above[1], abs=2e-8)
    assert below[2] == pytest.approx(above[2], abs=2e-8)


def test_smooth_phase_inverse_roundtrips_its_own_enthalpy_surface(phase):
    inverse = SmoothPhaseMassEnthalpyInverter(phase)
    rng = np.random.default_rng(240912)
    ideal = rng.uniform(inverse.t_grid[2], inverse.t_grid[-3], 256)
    fraction = rng.uniform(inverse.y_grid[2], inverse.y_grid[-3], 256)
    density = inverse.a/ideal/(1.+inverse.k*fraction)
    enthalpy = inverse.enthalpy_surface.evaluate(ideal, fraction)[0]
    found_density, found_fraction, found_temperature = inverse.state(
        density*fraction, enthalpy,
    )
    expected_temperature = inverse.temperature_surface.evaluate(ideal, fraction)[0]
    assert found_density == pytest.approx(density, rel=2e-9, abs=2e-11)
    assert found_fraction == pytest.approx(fraction, rel=2e-9, abs=2e-11)
    assert found_temperature == pytest.approx(expected_temperature, rel=2e-9, abs=2e-9)


def test_smooth_enthalpy_partials_match_independent_differences(phase):
    inverse = SmoothPhaseMassEnthalpyInverter(phase)
    density = np.array([1.18, 1.4, 2.])
    fuel_density = np.array([.03, .4, .793])
    h_rho, h_c = inverse.enthalpy_partials(density, fuel_density)
    step = 1e-6
    rho_plus = inverse.enthalpy_and_slope(density+step, fuel_density)[0]
    rho_minus = inverse.enthalpy_and_slope(density-step, fuel_density)[0]
    c_plus = inverse.enthalpy_and_slope(density, fuel_density+step)[0]
    c_minus = inverse.enthalpy_and_slope(density, fuel_density-step)[0]
    assert h_rho == pytest.approx((rho_plus-rho_minus)/(2*step), rel=2e-7)
    assert h_c == pytest.approx((c_plus-c_minus)/(2*step), rel=2e-7)
    assert np.all(h_rho < 0.)


def test_existing_phase_inverse_remains_the_default(phase):
    assert type(PhaseMassEnthalpyInverter(phase)) is PhaseMassEnthalpyInverter


def _smooth_section(candidate):
    return BuoyancyConstrainedEnthalpySection(
        candidate.jetplume, candidate.thermodynamics, quadrature_points=256,
        thermal_width_ratio=1.0235, phase_interpolation="c1_hermite",
    )


def test_smooth_reduced_moment_jacobian_matches_actual_flux_differences(candidate):
    section = _smooth_section(candidate)
    state = np.array([1.20255, .309375, .002295, -6.32e-5, 51.72, .414, 1.5])
    target_hydrogen = section.moments(state)[1]
    parameters = np.log([state[0], state[2], state[4], section.thermal_width_ratio])
    analytic = section.reduced_moment_jacobian(state)

    def evaluate(reduced):
        rho, area, velocity, beta = np.exp(reduced)
        section.thermal_width_ratio = beta
        trial = state.copy()
        trial[[0, 2, 4]] = rho, area, velocity
        denominator = area*rho*(
            section._wind(trial)*np.cos(trial[3])*section.k.profile_integral(1.)
            +velocity*section.k.profile_integral(1.+section.velocity_shape_exponent)
        )
        trial[1] = target_hydrogen/denominator
        moments = section.moments(trial)
        return np.array([
            moments[0], np.hypot(moments[2], moments[3]), moments[4], moments[5],
        ])

    step = 1e-7
    numerical = np.column_stack([
        (evaluate(parameters+step*unit)-evaluate(parameters-step*unit))/(2.*step)
        for unit in np.eye(4)
    ])
    assert analytic == pytest.approx(numerical, rel=7e-5, abs=2e-5)


def test_smooth_five_flux_jacobian_matches_fixed_partition_differences(candidate):
    section = _smooth_section(candidate)
    state = np.array([1.20255, .309375, .002295, -6.32e-5, 51.72, .414, 1.5])
    mixing = ConservativeTransverseMixing(section, state)
    analytic = mixing.flux_jacobian(order=16)
    q, weights = mixing.partition.nodes(16, square=True)
    q, weights = q.ravel(), weights.ravel()
    parameters = np.array([
        np.log(state[0]), np.log(state[0]*state[1]), np.log(state[2]), state[3],
        np.log(state[4]), state[5], state[6], np.log(section.thermal_width_ratio),
    ])

    def evaluate(encoded):
        rho, c, area, velocity, beta = np.exp(encoded[[0, 1, 2, 4, 7]])
        section.thermal_width_ratio = beta
        trial = np.array([rho, c/rho, area, encoded[3], velocity, encoded[5], encoded[6]])
        rr, yy, _, hh = section.thermodynamic_profile(trial, np.exp(-q))
        u = section._wind(trial)*np.cos(trial[3])+velocity*np.exp(-section.velocity_shape_exponent*q)
        return area*np.array([
            np.sum(weights*rr*u), np.sum(weights*rr*yy*u),
            np.sum(weights*rr*u*u)*np.cos(trial[3]),
            np.sum(weights*rr*u*u)*np.sin(trial[3]),
            np.sum(weights*(hh*u+.5*rr*u**3)),
        ])

    step = 1e-7
    numerical = np.column_stack([
        (evaluate(parameters+step*unit)-evaluate(parameters-step*unit))/(2.*step)
        for unit in np.eye(8)
    ])
    assert analytic == pytest.approx(numerical, rel=1e-4, abs=3e-4)
