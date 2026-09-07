"""Independent ledgers and limiting states for the bounded liquid screen."""
import math

import numpy as np
import pytest
from scipy.optimize import root

from degali.addons.mixed_air_liquid import ideal_flash, liquid_properties, nbs_activity_coefficients, R


@pytest.mark.parametrize('inert,n,o,kn,ko', [
    (3., 2., 1., .2, .03), (1e-10, .79, .21, .6, .1),
    (2., 3., 4., 1.5, .08), (1., .5, .3, .1, .05),
])
@pytest.mark.parametrize('q', [0., .01])
def test_amounts_and_independent_nonlinear_equilibrium(inert, n, o, kn, ko, q):
    f = ideal_flash(inert, n, o, kn, ko, reservoir_vapor_fraction=q)
    assert f.gas_nitrogen_mol+f.liquid_nitrogen_mol == pytest.approx(n, abs=1e-13)
    assert f.gas_oxygen_mol+f.liquid_oxygen_mol == pytest.approx(o, abs=1e-13)
    assert f.reservoir_vapor_mol/f.gas_mol == pytest.approx(q, abs=1e-13)
    if f.liquid_mol > 0:
        # Solve the TWO component balances at prescribed gas/liq composition,
        # independently of the Rachford-Rice implementation.
        def equations(values):
            v, x = values
            liq = n+o-v*(kn*x+ko*(1-x))
            return [inert-v*(1-q-kn*x-ko*(1-x)), n-v*kn*x-liq*x]
        solution = root(equations, [f.gas_mol*1.01, f.liquid_nitrogen_mol/f.liquid_mol*.99])
        assert solution.success
        assert solution.x[0] == pytest.approx(f.gas_mol, rel=1e-9)
        assert solution.x[1] == pytest.approx(f.liquid_nitrogen_mol/f.liquid_mol, rel=1e-9)


def test_shared_liquid_can_condense_before_either_pure_species():
    # yN/KN=.6, yO/KO=.6, each unsaturated against pure condensate,
    # yet a shared liquid is stable (sum1.2).
    f = ideal_flash(.52, .36, .12, .6, .2)
    assert .36 < .6 and .12 < .2
    assert f.dew_index == pytest.approx(1.2)
    assert f.liquid_mol > 0


@pytest.mark.parametrize('gas', [(1.,0.), (0.,1.)])
def test_pure_condensable_reduces_to_pure_saturation(gas):
    f = ideal_flash(1., *gas, .1, .1)
    assert f.gas_nitrogen_mol+f.gas_oxygen_mol == pytest.approx(1/9.)
    assert f.liquid_mol == pytest.approx(8/9.)


def test_inert_only_and_all_gas():
    f = ideal_flash(2., 0., 0., .1, .01)
    assert f.gas_mol == 2 and f.liquid_mol == 0
    g = ideal_flash(1000., 1., .1, .1, .01)
    assert g.liquid_mol == 0 and g.gas_mol == 1001.1


@pytest.mark.parametrize('scale', [1e-15, 1., 1e15])
def test_extensivity(scale):
    a = ideal_flash(2., 1., .5, .2, .03)
    b = ideal_flash(2*scale, scale, .5*scale, .2, .03)
    assert b.gas_mol/scale == pytest.approx(a.gas_mol, rel=1e-12)
    assert b.liquid_mol/scale == pytest.approx(a.liquid_mol, rel=1e-12)


def test_zero_inert_liquid_and_vapor_limits():
    assert ideal_flash(0., .8, .2, .1, .01).vapor_fraction_without_reservoir == 0
    assert ideal_flash(0., .8, .2, 2., 3.).vapor_fraction_without_reservoir == 1
    mixed = ideal_flash(0., .8, .2, 1.5, .1)
    assert 0 < mixed.vapor_fraction_without_reservoir < 1


@pytest.mark.parametrize('temperature', [20., 60., 63., 121., float('nan')])
def test_no_silent_solid_or_supercritical_extrapolation(temperature):
    with pytest.raises(ValueError):
        liquid_properties(temperature)


def test_activity_reference_and_gibbs_duhem():
    gn, go = nbs_activity_coefficients(65., .5)
    phi = 32.58/(32.58+25.32)
    assert math.log(gn) == pytest.approx(32.58*1.47*(1-phi)**2/(R/4.184*65.))
    assert math.log(go) == pytest.approx(25.32*1.47*phi**2/(R/4.184*65.))
    x, eps = .37, 1e-5
    plus = np.log(nbs_activity_coefficients(70., x+eps))
    minus = np.log(nbs_activity_coefficients(70., x-eps))
    derivative = (plus-minus)/(2*eps)
    assert abs(x*derivative[0]+(1-x)*derivative[1]) < 1e-9


@pytest.mark.parametrize('args', [(0.,0.,0.,1.,1.), (-1.,1.,1.,1.,1.), (1.,1.,1.,0.,1.),
                                   (1.,1.,1.,float('nan'),1.)])
def test_invalid_flash_inputs(args):
    with pytest.raises(ValueError):
        ideal_flash(*args)
