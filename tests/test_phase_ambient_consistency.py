"""No-release fixed point and dilute limit, independent of field observations."""

from functools import lru_cache
import math

import numpy as np
import pytest

pytest.importorskip("CoolProp")
from CoolProp.CoolProp import PropsSI

from degali.addons.axisymmetric_jet import (
    AxisymmetricJetSource, ConservedGaussianJet, R_UNIVERSAL,
    phase_ambient_from_rh, phase_dry_air_composition,
)


@lru_cache(maxsize=None)
def _model(consistent=True, argon=False, rh=53.6666667):
    ta, pa = 288.65, 101325.0
    mw = PropsSI("M", "Air")
    mw_w = PropsSI("M", "Water")
    pv = rh / 100.0 * PropsSI("P", "T", ta, "Q", 0, "Water")
    humidity = mw_w / mw * pv / (pa - pv)
    mw_humid = (1 + humidity) / (1 / mw + humidity / mw_w)
    density = PropsSI("D", "T", ta, "P", pa, "Air") * mw_humid / mw
    if consistent:
        mw, humidity, density = phase_ambient_from_rh(ta, pa, rh, include_argon=argon)
    return ConservedGaussianJet(
        AxisymmetricJetSource(diameter=0.01, velocity=100., density=1., temperature=100.),
        ambient_temperature=ta, ambient_pressure=pa, ambient_density=density,
        fuel_molecular_weight=PropsSI("M", "Hydrogen"), ambient_molecular_weight=mw,
        fuel_heat_capacity=PropsSI("C", "T", ta, "P", pa, "Hydrogen"),
        ambient_heat_capacity=PropsSI("C", "T", ta, "P", pa, "Air"),
        equilibrium_air_condensation=True, temperature_dependent_phase_enthalpy=True,
        equilibrium_argon_condensation=argon, ambient_absolute_humidity=humidity,
        consistent_phase_ambient=consistent, radial_points=41,
        conservative_establishment="entrained_mass",
    )


@pytest.mark.parametrize("argon", [False, True])
@pytest.mark.parametrize("rh", [0., 53.6666667, 100.])
def test_exact_and_interpolated_ambient_fixed_point(argon, rh):
    model = _model(True, argon, rh)
    rho, y = np.array([model.ambient_density]), np.array([0.])
    for inverse in (model._condensed_air_state_exact, model._condensed_air_state):
        t, rho_h = inverse(rho, y)
        assert t[0] == pytest.approx(model.ambient_temperature, abs=1e-7)
        assert rho_h[0] == pytest.approx(0., abs=1e-4)
    mw, fractions = phase_dry_air_composition(argon)
    assert sum(fractions) == pytest.approx(1.)
    assert model.ambient_molecular_weight == mw
    assert model.ambient_density == pytest.approx(
        model.ambient_pressure * model._humid_ambient_molecular_weight
        / (R_UNIVERSAL * model.ambient_temperature), rel=1e-14,
    )


@pytest.mark.parametrize("argon,offset", [(False, -1.196714), (True, -0.173999)])
def test_legacy_ambient_defect_is_reproduced_not_silently_changed(argon, offset):
    model = _model(False, argon)
    t, rho_h = model._condensed_air_state_exact(np.array([model.ambient_density]), np.array([0.]))
    assert t[0] - model.ambient_temperature == pytest.approx(offset, abs=2e-6)
    assert rho_h[0] < -200.


@pytest.mark.parametrize("argon", [False, True])
def test_dilute_profile_has_no_background_enthalpy_offset(argon):
    model = _model(True, argon)
    g = np.array([1e-4, 1e-6, 1e-8, 1e-10])
    rho = model.ambient_density + (1.11 - model.ambient_density) * g
    mass_density = 0.1 * g
    t, rho_h = model._condensed_air_state(rho, mass_density / rho)
    # A finite, bounded enthalpy/fuel ratio as Y->0 rules out an ambient
    # energy pedestal. Its value is not fitted to a measured temperature.
    assert np.all(np.abs(rho_h / mass_density) < 1e7)
    assert abs(t[-1] - model.ambient_temperature) < 1e-6
    assert abs(rho_h[-1]) < abs(rho_h[0]) * 2e-6


def test_exact_warm_gas_and_pure_hydrogen_limits():
    model = _model()
    y = np.array([0., 0.01, 0.2, 1.])
    temperature = 295.0
    mw = 1 / (y / model.fuel_molecular_weight + (1-y) / model._humid_ambient_molecular_weight)
    rho = model.ambient_pressure * mw / (R_UNIVERSAL * temperature)
    t, rho_h = model._condensed_air_state_exact(rho, y)
    assert t == pytest.approx(np.full(4, temperature), abs=1e-7)
    assert rho_h / rho == pytest.approx(model._mixture_enthalpy(temperature, y), abs=1e-3)


@pytest.mark.parametrize("ta,pa,rh", [
    (301., 101325., 50.), (math.nan, 101325., 50.),
    (295., 0., 50.), (295., 101325., 101.),
])
def test_phase_ambient_rejects_invalid_inputs(ta, pa, rh):
    with pytest.raises(ValueError):
        phase_ambient_from_rh(ta, pa, rh)
