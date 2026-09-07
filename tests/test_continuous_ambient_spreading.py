import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.core.jetplume import JetCoefficients
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.continuous_ambient_spreading import with_continuous_ambient_spreading


def model(
    *,
    spread=.1,
    wind=2.,
    energy='total',
    thermal_relaxation_rate=0.0,
    phase_transition_lag_rate=0.0,
    turbulence_heat_exchange_rate=0.0,
    equilibrium_air_condensation=False,
    condensed_enthalpy_specific=0.0,
    condensed_air_temperature=80.0,
):
    def temperature_from_density(density, fraction):
        density = np.asarray(density, dtype=float)
        fraction = np.asarray(fraction, dtype=float)
        mixture_molecular_weight = 1.0 / (
            fraction / 0.002016 + (1.0 - fraction) / 0.028965
        )
        return 101325.0 * mixture_molecular_weight / (8.314462618 * density)

    jp = SimpleNamespace(th=None, k=JetCoefficients(sc=1.16**2),
        deltay=spread, deltaz=.8*spread, betay=.92, betaz=.87, gammaz=-.012,
        _split=lambda area, *_: (math.sqrt(area), math.sqrt(area)),
        _wind=lambda *_: wind, _wind_profile=lambda *_: (wind, 0.))
    th = SimpleNamespace(
        ambient_density=1.2, spreading_ratio=1.16, ambient_pressure=101325.0,
        ambient_temperature=295.0, equilibrium_air_condensation=False,
        fuel_molecular_weight=0.002016, _humid_ambient_molecular_weight=0.028965,
        _ambient_enthalpy=0.0, _buoyancy_coefficient=.01,
        plume_entrainment_limit=.12,
        _temperature_from_density=temperature_from_density,
        _mixture_enthalpy=lambda temperature, fraction: 1000.0 * (np.asarray(temperature)-295.0),
        _condensed_air_state=lambda density, fraction: (
            np.full_like(np.asarray(density, dtype=float), condensed_air_temperature, dtype=float),
            np.full_like(np.asarray(density, dtype=float), condensed_enthalpy_specific, dtype=float),
        ),
    )
    th.equilibrium_air_condensation = bool(equilibrium_air_condensation)
    return IndependentEnergyCrosswind(
        jp,
        th,
        energy_transport=energy,
        thermal_relaxation_rate=thermal_relaxation_rate,
        phase_transition_lag_rate=phase_transition_lag_rate,
        turbulence_heat_exchange_rate=turbulence_heat_exchange_rate,
    )


def state(x):
    return np.array([1.1, .02, .03, .1, 15., x, 1.5])


def product(jp, x):
    return jp.deltay*x**jp.betay * jp.deltaz*x**jp.betaz * math.exp(jp.gammaz*math.log(x)**2)


@pytest.mark.parametrize('x', [.2, .5, 1., 2.])
def test_background_area_growth_matches_independent_finite_difference(x):
    old = model()
    new = with_continuous_ambient_spreading(old)
    q = state(x)
    g = new._geometry(q)
    represented = g[8]*g[9] + g[7]*g[10]
    for step in (2e-5, 1e-5):
        dx = step*math.cos(q[3])
        finite = (product(old.jetplume, x+dx)-product(old.jetplume, x-dx))/(2*step)
        assert represented == pytest.approx(finite, rel=5e-8)


@pytest.mark.parametrize('energy', ['total', 'enthalpy'])
@pytest.mark.parametrize('x', [.25, .75, 1.])
def test_increment_carries_ambient_mass_momentum_and_energy_once(x, energy):
    old = model(energy=energy)
    new = with_continuous_ambient_spreading(old)
    q = state(x)
    before = old._geometry(q)
    actual = new.source_terms(q)-old.source_terms(q)
    h = 1e-5
    dx = h*math.cos(q[3])
    darea = (product(old.jetplume, x+dx)-product(old.jetplume, x-dx))/(2*h)
    ua = before[2]
    dm = old.rhoa*old.k.profile_integral(0.)*ua*darea
    expected = [dm, 0., ua*dm, 0., .5*ua*ua*dm if energy == 'total' else 0.]
    assert actual == pytest.approx(expected, rel=5e-8, abs=1e-12)
    assert old._geometry(q) == before


@pytest.mark.parametrize('x', [1.0001, 2., 6.])
def test_above_one_metre_is_bitwise_unchanged(x):
    old = model()
    new = with_continuous_ambient_spreading(old)
    assert np.array_equal(old.source_terms(state(x)), new.source_terms(state(x)))
    assert old._geometry(state(x)) == new._geometry(state(x))


@pytest.mark.parametrize('spread,wind', [(0., 2.), (.1, 0.)])
def test_no_extra_entrainment_without_background_spread_or_wind(spread, wind):
    old = model(spread=spread, wind=wind)
    new = with_continuous_ambient_spreading(old)
    assert np.array_equal(old.source_terms(state(.5)), new.source_terms(state(.5)))


def test_no_one_metre_source_jump_in_new_option():
    old = model()
    new = with_continuous_ambient_spreading(old)
    left, right = state(1.-1e-8), state(1.+1e-8)
    assert np.max(abs(new.source_terms(left)-new.source_terms(right))) < 1e-7
    assert old.source_terms(right)[0]-old.source_terms(left)[0] > 1e-3


@pytest.mark.parametrize('x', [0., -.1])
def test_nonpositive_coordinate_is_refused(x):
    with pytest.raises(ValueError):
        with_continuous_ambient_spreading(model()).source_terms(state(x))


def test_factory_refuses_to_drop_other_model_subclass_behavior():
    class Unknown(IndependentEnergyCrosswind):
        pass
    old = model()
    other = Unknown(old.jetplume, old.thermodynamics)
    with pytest.raises(ValueError):
        with_continuous_ambient_spreading(other)


def test_total_energy_work_option_adds_buoyancy_projection_term():
    release = model()
    q = state(1.2)
    base = release.source_terms(q)
    with_work = release.source_terms(q, include_work=True)
    _, unit_weights = release._quadrature(release.quadrature_points)
    # Reconstruct the same profile basis used by the model for this source term
    exponent, _ = release._quadrature(release.quadrature_points)
    scalar_shape = np.exp(-exponent)
    velocity_shape = np.exp(-release.velocity_shape_exponent * exponent)
    rho_c, _fraction_c, _sysz, theta, uc, _x, _z = q
    density = release.rhoa + (rho_c - release.rhoa) * scalar_shape
    velocity = release._wind(q) * math.cos(theta) + uc * velocity_shape
    expected = (
        9.81 * math.sin(q[3])
        * np.sum((release.rhoa - density) * velocity * (q[2] * unit_weights))
    )
    assert with_work[4] == pytest.approx(base[4] + expected)
    assert release.source_terms(q, include_work=True)[4] > base[4]


def test_enthalpy_transport_keeps_work_term_off():
    release = model(energy="enthalpy")
    q = state(1.2)
    assert release.source_terms(q, include_work=True)[4] == pytest.approx(
        release.source_terms(q, include_work=False)[4]
    )


def test_thermal_relaxation_default_is_noop():
    release = model()
    q = state(1.2)
    assert release.source_terms(
        q, include_thermal_relaxation=True
    ) == pytest.approx(release.source_terms(q))


def test_thermal_relaxation_adds_ambient_enthalpy_closure():
    release = model(thermal_relaxation_rate=0.75)
    q = state(1.2)
    base = release.source_terms(q, include_thermal_relaxation=False)
    relaxed = release.source_terms(q, include_thermal_relaxation=True)
    velocity, density, _fraction, _temperature, rho_h = release.profiles(q)
    _, unit_weights = release._quadrature(release.quadrature_points)
    weights = q[2] * unit_weights
    mass_flux = float(np.sum(density * velocity * weights))
    mean_relative_enthalpy = float(np.sum(
        (rho_h - release.thermodynamics._ambient_enthalpy * density)
        * velocity * weights
    ) / mass_flux)
    expected = base[4] - release.thermal_relaxation_rate * base[0] * mean_relative_enthalpy
    assert relaxed[4] == pytest.approx(expected)
    if mean_relative_enthalpy < 0.0:
        assert relaxed[4] > base[4]


def test_work_and_thermal_relaxation_terms_are_additive():
    release = model(thermal_relaxation_rate=0.75)
    q = state(1.2)
    base = release.source_terms(q)
    with_work = release.source_terms(q, include_work=True, include_thermal_relaxation=False)
    with_thermal = release.source_terms(
        q, include_work=False, include_thermal_relaxation=True,
    )
    with_both = release.source_terms(
        q, include_work=True, include_thermal_relaxation=True,
    )
    assert with_both[4] == pytest.approx(
        base[4] + (with_work[4] - base[4]) + (with_thermal[4] - base[4])
    )


def test_phase_transition_relaxation_default_is_noop():
    release = model()
    q = state(1.2)
    assert release.source_terms(
        q, include_phase_transition=True
    ) == pytest.approx(release.source_terms(q))


def test_phase_transition_relaxation_adds_frozen_to_equilibrium_correction():
    release = model(
        equilibrium_air_condensation=True,
        phase_transition_lag_rate=0.75,
        condensed_enthalpy_specific=-5_000.0,
    )
    q = state(1.2)
    base = release.source_terms(q, include_phase_transition=False)
    corrected = release.source_terms(q, include_phase_transition=True)
    mean_frozen, *_ = release._mean_section_quantities(
        q, use_condensed_air=False
    )
    mean_equilibrium, equilibrium_temperature, *_ = (
        release._mean_section_quantities(
            q, use_condensed_air=True
        )
    )
    phase_gap = mean_frozen - mean_equilibrium
    subcooling = max(
        0.0, release.thermodynamics.ambient_temperature - equilibrium_temperature
    )
    activity = min(subcooling / (2.0 * release.thermodynamics.ambient_temperature), 1.0)
    expected = base[4] + 0.75 * base[0] * phase_gap * activity
    assert corrected[4] == pytest.approx(expected)
    if mean_frozen > mean_equilibrium and activity > 0.0:
        assert corrected[4] > base[4]


def test_turbulent_heat_exchange_default_is_noop():
    release = model()
    q = state(1.2)
    assert release.source_terms(
        q, include_turbulent_heat_exchange=True
    ) == pytest.approx(release.source_terms(q))


def test_turbulent_heat_exchange_reduces_cold_energy_source():
    release = model(turbulence_heat_exchange_rate=0.05)
    # Force a warm relative specific enthalpy so the sign of exchange is clear.
    release.thermodynamics._mixture_enthalpy = lambda temperature, fraction: (
        1000.0 * (np.asarray(temperature)-273.15)
    )
    q = state(1.2)
    base = release.source_terms(q)
    warmed = release.source_terms(
        q, include_turbulent_heat_exchange=True
    )
    mean_relative_enthalpy, mean_temperature, mean_fraction, mean_cp = (
        release._mean_section_quantities(
            q, use_condensed_air=release.thermodynamics.equilibrium_air_condensation
        )
    )
    thermal_deficit = mean_cp * (
        mean_temperature - release.thermodynamics.ambient_temperature
    )
    _, sysz, _ua, perimeter, *_rest = release._geometry(q)
    transfer_scale = perimeter / (4.0 * math.sqrt(sysz))
    expected = base[4] - release.turbulence_heat_exchange_rate * base[0] * (
        thermal_deficit * transfer_scale
    )
    assert warmed[4] == pytest.approx(expected)
    if thermal_deficit > 0.0:
        assert warmed[4] < base[4]


def test_mean_section_quantities_tracks_existing_helpers():
    release = model()
    q = state(1.2)
    mean_h, mean_t, mean_fraction, mean_cp = release._mean_section_quantities(
        q, use_condensed_air=False
    )
    assert mean_h == pytest.approx(release._mean_relative_specific_enthalpy(q))
    assert mean_t == pytest.approx(release._mean_temperature(q))
    assert 0.0 < mean_fraction <= 1.0
    assert mean_cp > 0.0


def test_phase_transition_and_turbulent_terms_are_additive():
    release = model(
        equilibrium_air_condensation=True,
        phase_transition_lag_rate=0.4,
        turbulence_heat_exchange_rate=0.05,
        condensed_enthalpy_specific=-5_000.0,
    )
    release.thermodynamics._mixture_enthalpy = lambda temperature, fraction: (
        1000.0 * (np.asarray(temperature)-273.15)
    )
    q = state(1.2)
    base = release.source_terms(q)
    with_both = release.source_terms(
        q, include_phase_transition=True, include_turbulent_heat_exchange=True,
    )
    with_phase = release.source_terms(
        q, include_phase_transition=True, include_turbulent_heat_exchange=False,
    )
    with_turb = release.source_terms(
        q, include_phase_transition=False, include_turbulent_heat_exchange=True,
    )
    assert with_both[4] == pytest.approx(
        base[4]
        + (with_phase[4] - base[4])
        + (with_turb[4] - base[4])
    )
