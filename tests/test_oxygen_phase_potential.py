"""Independent thermodynamic identities for bounded gamma-O2 candidate."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from degali.addons.oxygen_phase_potential import (
    MIN_T, MAX_T, MAX_P, FUSION_J_MOL, CP_COEFFICIENTS,
    solid, ideal_gas, triple_reference, volume_and_derivatives,
    equilibrium_partial_pressure, pure_vapor_pressure,
    ideal_gas_to_solid_enthalpy, nbs1977_vapor_pressure,
)


@pytest.mark.parametrize('t', [43.0,54.362,float('nan'),float('inf'),0.0])
def test_temperature_domain(t):
    for operation in (lambda:solid(t,100000), lambda:ideal_gas(t,10),
                      lambda:pure_vapor_pressure(t)):
        with pytest.raises(ValueError):
            operation()


@pytest.mark.parametrize('p', [-1,0,200001,float('nan'),float('inf')])
def test_pressure_domain(p):
    with pytest.raises(ValueError):
        solid(50,p)
    with pytest.raises(ValueError):
        equilibrium_partial_pressure(50,p)


def test_triple_reference_and_fusion():
    r = triple_reference()
    s = solid(MAX_T,r['pressure_Pa'])
    assert s.gibbs_J_mol == pytest.approx(r['liquid_g_J_mol'],abs=1e-10)
    assert r['liquid_h_J_mol']-s.enthalpy_J_mol == pytest.approx(FUSION_J_MOL,abs=1e-10)
    assert r['liquid_s_J_mol_K']-s.entropy_J_mol_K == pytest.approx(FUSION_J_MOL/MAX_T,abs=1e-12)
    # Fugacity is not mislabeled as a real-gas saturation pressure.
    peq = pure_vapor_pressure(MAX_T)
    assert 0.995 < peq/r['pressure_Pa'] < 1


@pytest.mark.parametrize('t', [44.0,47.0,50.0,54.0])
def test_gibbs_derivatives_and_maxwell(t):
    p, dt, dp = 101325., .001, 100.
    c = solid(t,p)
    tm,tp = solid(t-dt,p),solid(t+dt,p)
    pm,pp = solid(t,p-dp),solid(t,p+dp)
    assert c.enthalpy_J_mol == pytest.approx(c.gibbs_J_mol+t*c.entropy_J_mol_K,abs=2e-12)
    assert -(tp.gibbs_J_mol-tm.gibbs_J_mol)/(2*dt) == pytest.approx(c.entropy_J_mol_K,rel=1e-9)
    assert (pp.gibbs_J_mol-pm.gibbs_J_mol)/(2*dp) == pytest.approx(c.volume_m3_mol,rel=1e-8)
    assert (tp.enthalpy_J_mol-tm.enthalpy_J_mol)/(2*dt) == pytest.approx(c.heat_capacity_J_mol_K,rel=1e-9)
    assert (pp.entropy_J_mol_K-pm.entropy_J_mol_K)/(2*dp) == pytest.approx(-volume_and_derivatives(t)[1],rel=1e-8)


def test_caloric_integrals_independent_quad():
    r = triple_reference()
    a,b,c = CP_COEFFICIENTS
    for t in [MIN_T,47,50,MAX_T-1e-7,MAX_T]:
        s = solid(t,r['pressure_Pa'])
        assert s.enthalpy_J_mol == pytest.approx(r['solid_h_J_mol']+quad(lambda x:a+b*x+c*x*x,MAX_T,t)[0],abs=1e-10)
        assert s.entropy_J_mol_K == pytest.approx(r['solid_s_J_mol_K']+quad(lambda x:(a+b*x+c*x*x)/x,MAX_T,t)[0],abs=1e-10)


@pytest.mark.parametrize('t', [44.,47.,50.,54.])
def test_equilibrium_and_clapeyron(t):
    p = pure_vapor_pressure(t)
    s,g = solid(t,p),ideal_gas(t,p)
    assert s.gibbs_J_mol == pytest.approx(g.gibbs_J_mol,abs=2e-9)
    dt = 1e-3
    derivative = (pure_vapor_pressure(t+dt)-pure_vapor_pressure(t-dt))/(2*dt)
    latent = g.enthalpy_J_mol-s.enthalpy_J_mol
    clapeyron = t*(g.volume_m3_mol-s.volume_m3_mol)*derivative
    assert clapeyron == pytest.approx(latent,rel=1e-7)
    assert ideal_gas_to_solid_enthalpy(t,p) == pytest.approx(latent,abs=1e-10)


def test_total_pressure_is_not_partial_pressure():
    t,pm,pp = 50.,50000.,150000.
    v = volume_and_derivatives(t)[0]
    r = triple_reference()['gas_constant']
    ratio = equilibrium_partial_pressure(t,pp)/equilibrium_partial_pressure(t,pm)
    assert ratio == pytest.approx(math.exp(v*(pp-pm)/(r*t)),rel=2e-14)
    peq = equilibrium_partial_pressure(t,pp)
    assert ideal_gas(t,peq).gibbs_J_mol == pytest.approx(solid(t,pp).gibbs_J_mol,abs=1e-9)


def test_ideal_gas_identities():
    t,p,dt = 50.,100.,1e-3
    c = ideal_gas(t,p)
    left,right = ideal_gas(t-dt,p),ideal_gas(t+dt,p)
    assert -(right.gibbs_J_mol-left.gibbs_J_mol)/(2*dt) == pytest.approx(c.entropy_J_mol_K,rel=1e-9)
    assert (right.enthalpy_J_mol-left.enthalpy_J_mol)/(2*dt) == pytest.approx(c.heat_capacity_J_mol_K,rel=1e-9)


def test_full_domain_positive_and_pressure_monotonicity():
    previous = 0
    for t in np.linspace(MIN_T,MAX_T,301):
        for p in [1e-8,146.3,101325.,MAX_P]:
            state = solid(float(t),p)
            assert state.volume_m3_mol > 0
            assert state.heat_capacity_J_mol_K > 0
        current = pure_vapor_pressure(float(t))
        assert current > previous
        previous = current


def test_published_comparator_transcription():
    assert nbs1977_vapor_pressure(44)/101325 == pytest.approx(.000019200,rel=1e-7)
    assert nbs1977_vapor_pressure(54.359)/101325 == pytest.approx(.0014451,rel=1e-7)
    with pytest.raises(ValueError):
        nbs1977_vapor_pressure(MAX_T)


def test_reference_cannot_be_mutated_by_caller():
    reference = triple_reference()
    original = reference['solid_h_J_mol']
    reference['solid_h_J_mol'] = 0
    assert triple_reference()['solid_h_J_mol'] == original
