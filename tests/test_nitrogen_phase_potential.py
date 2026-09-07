"""Independent thermodynamic tests, not dispersion validation."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from degali.addons.nitrogen_phase_potential import (
    MIN_T,MAX_T,MAX_P,FUSION_J_MOL,VOLUME_M3_MOL,CP_COEFFICIENTS,
    solid,ideal_gas,triple_reference,equilibrium_partial_pressure,
    pure_vapor_pressure,ideal_gas_to_solid_enthalpy,
)


@pytest.mark.parametrize('t',[35.61,63.152,0,float('nan'),float('inf')])
def test_domain_t(t):
    for fn in [lambda:solid(t,100000),lambda:ideal_gas(t,10),lambda:pure_vapor_pressure(t)]:
        with pytest.raises(ValueError): fn()


@pytest.mark.parametrize('p',[-1,0,200001,float('nan'),float('inf')])
def test_domain_p(p):
    with pytest.raises(ValueError): solid(60,p)
    with pytest.raises(ValueError): equilibrium_partial_pressure(60,p)


def test_reference_and_fusion():
    r = triple_reference()
    s = solid(MAX_T,r['pressure_Pa'])
    assert s.gibbs_J_mol == pytest.approx(r['liquid_g_J_mol'],abs=1e-10)
    assert r['liquid_h_J_mol']-s.enthalpy_J_mol == pytest.approx(FUSION_J_MOL,abs=1e-10)
    assert r['liquid_s_J_mol_K']-s.entropy_J_mol_K == pytest.approx(FUSION_J_MOL/MAX_T,abs=1e-12)
    assert 0.98 < pure_vapor_pressure(MAX_T)/r['pressure_Pa'] < 1
    r['solid_h_J_mol'] = 0
    assert triple_reference()['solid_h_J_mol'] != 0


@pytest.mark.parametrize('t',[36.,45.,50.,59.,60.,61.,62.,63.])
def test_gibbs_derivatives(t):
    p,dt,dp = 101325.,.001,100.
    s = solid(t,p)
    left,right = solid(t-dt,p),solid(t+dt,p)
    down,up = solid(t,p-dp),solid(t,p+dp)
    assert s.enthalpy_J_mol == pytest.approx(s.gibbs_J_mol+t*s.entropy_J_mol_K,abs=2e-12)
    assert -(right.gibbs_J_mol-left.gibbs_J_mol)/(2*dt) == pytest.approx(s.entropy_J_mol_K,rel=2e-9)
    assert (right.enthalpy_J_mol-left.enthalpy_J_mol)/(2*dt) == pytest.approx(s.heat_capacity_J_mol_K,rel=2e-9)
    assert (up.gibbs_J_mol-down.gibbs_J_mol)/(2*dp) == pytest.approx(VOLUME_M3_MOL,rel=1e-8)
    assert up.entropy_J_mol_K == down.entropy_J_mol_K


def test_caloric_integral_against_quad():
    r = triple_reference()
    def cp(t): return sum(a*t**j for j,a in enumerate(CP_COEFFICIENTS))
    for t in [MIN_T,40,50,60,MAX_T-1e-7,MAX_T]:
        s = solid(t,r['pressure_Pa'])
        assert s.enthalpy_J_mol == pytest.approx(r['solid_h_J_mol']+quad(cp,MAX_T,t)[0],abs=5e-10)
        assert s.entropy_J_mol_K == pytest.approx(r['solid_s_J_mol_K']+quad(lambda x:cp(x)/x,MAX_T,t)[0],abs=1e-10)


@pytest.mark.parametrize('t',[36.,45.,50.,59.,60.,61.,62.,63.])
def test_equilibrium_and_clapeyron(t):
    p = pure_vapor_pressure(t)
    s,g = solid(t,p),ideal_gas(t,p)
    assert s.gibbs_J_mol == pytest.approx(g.gibbs_J_mol,abs=2e-9)
    dt = .001
    dpdt = (pure_vapor_pressure(t+dt)-pure_vapor_pressure(t-dt))/(2*dt)
    implied = t*(g.volume_m3_mol-s.volume_m3_mol)*dpdt
    assert implied == pytest.approx(g.enthalpy_J_mol-s.enthalpy_J_mol,rel=1e-7)
    assert ideal_gas_to_solid_enthalpy(t,p) == pytest.approx(implied,rel=1e-7)


def test_poynting_total_partial_distinction():
    t = 60.
    ratio = equilibrium_partial_pressure(t,150000)/equilibrium_partial_pressure(t,50000)
    r = triple_reference()['gas_constant']
    assert ratio == pytest.approx(math.exp(VOLUME_M3_MOL*100000/(r*t)),rel=3e-14)


def test_ideal_gas_reference():
    t,p,dt = 60.,100.,.001
    s = ideal_gas(t,p)
    left,right = ideal_gas(t-dt,p),ideal_gas(t+dt,p)
    assert -(right.gibbs_J_mol-left.gibbs_J_mol)/(2*dt) == pytest.approx(s.entropy_J_mol_K,rel=1e-9)
    assert (right.enthalpy_J_mol-left.enthalpy_J_mol)/(2*dt) == pytest.approx(s.heat_capacity_J_mol_K,rel=1e-9)


def test_full_domain_positive_monotonic():
    previous = 0.
    for t in np.linspace(MIN_T,MAX_T,401):
        for p in [1e-8,12519.,101325.,MAX_P]:
            s = solid(float(t),p)
            assert s.heat_capacity_J_mol_K > 0
            assert s.volume_m3_mol > 0
        p = pure_vapor_pressure(float(t))
        assert p > previous
        previous = p
