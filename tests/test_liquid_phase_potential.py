"""Independent identities/minimization for the common liquid reference."""
import math

import CoolProp as CP
import numpy as np
import pytest
from scipy.optimize import minimize

from degali.addons import liquid_phase_potential as lp
from degali.addons import nitrogen_phase_potential as ns
from degali.addons import oxygen_phase_potential as os


@pytest.mark.parametrize('sp,t', [('Nitrogen',63),('Oxygen',54),('Argon',70),('Nitrogen',101),('Nitrogen',float('nan'))])
def test_liquid_domain(sp,t):
    with pytest.raises(ValueError):
        lp.liquid(sp,t,101325)


@pytest.mark.parametrize('p', [0,-1,200001,float('nan'),float('inf')])
def test_pressure_domain(p):
    with pytest.raises(ValueError):
        lp.liquid('Nitrogen',70,p)
    with pytest.raises(ValueError):
        lp.flash(70,p,.3,.5,.2)


@pytest.mark.parametrize('sp,t', [('Nitrogen',64),('Nitrogen',75),('Oxygen',55),('Oxygen',70),('Oxygen',90)])
def test_exact_saturation_reference(sp,t):
    st = CP.AbstractState('HEOS',sp)
    st.update(CP.QT_INPUTS,0.,t)
    c = lp.liquid(sp,t,st.p())
    assert c.gibbs_J_mol == pytest.approx(st.gibbsmolar(),abs=1e-12)
    assert c.enthalpy_J_mol == pytest.approx(st.hmolar(),abs=1e-12)
    assert c.entropy_J_mol_K == pytest.approx(st.smolar(),abs=1e-12)
    assert c.volume_m3_mol == pytest.approx(1/st.rhomolar(),rel=1e-14)


@pytest.mark.parametrize('sp,t', [('Nitrogen',64),('Nitrogen',70),('Nitrogen',95),('Oxygen',55),('Oxygen',80),('Oxygen',99)])
def test_potential_derivatives(sp,t):
    p,dt,dp = 101325.,1e-3,100.
    c = lp.liquid(sp,t,p)
    tm,tp = lp.liquid(sp,t-dt,p),lp.liquid(sp,t+dt,p)
    pm,pp = lp.liquid(sp,t,p-dp),lp.liquid(sp,t,p+dp)
    ref = lp.saturated_reference(sp,t)
    assert c.gibbs_J_mol+t*c.entropy_J_mol_K == pytest.approx(c.enthalpy_J_mol,abs=2e-10)
    assert -(tp.gibbs_J_mol-tm.gibbs_J_mol)/(2*dt) == pytest.approx(c.entropy_J_mol_K,rel=2e-8)
    assert (tp.enthalpy_J_mol-tm.enthalpy_J_mol)/(2*dt) == pytest.approx(c.heat_capacity_J_mol_K,rel=2e-8)
    assert (pp.gibbs_J_mol-pm.gibbs_J_mol)/(2*dp) == pytest.approx(c.volume_m3_mol,rel=1e-8)
    assert (pp.entropy_J_mol_K-pm.entropy_J_mol_K)/(2*dp) == pytest.approx(-ref.volume_T_m3_mol_K,rel=1e-8)
    left,right = lp.saturated_reference(sp,t-dt),lp.saturated_reference(sp,t+dt)
    assert (right.volume_T_m3_mol_K-left.volume_T_m3_mol_K)/(2*dt) == pytest.approx(ref.volume_TT_m3_mol_K2,rel=2e-7)


@pytest.mark.parametrize('sp,solid', [('Nitrogen',ns),('Oxygen',os)])
def test_solid_liquid_common_triple_reference(sp,solid):
    ref = solid.triple_reference()
    t,p = ref['temperature_K'],ref['pressure_Pa']
    l,s = lp.liquid(sp,t,p),solid.solid(t,p)
    assert l.gibbs_J_mol == pytest.approx(s.gibbs_J_mol,abs=1e-9)
    assert l.enthalpy_J_mol-s.enthalpy_J_mol == pytest.approx(solid.FUSION_J_MOL,abs=1e-9)
    assert l.entropy_J_mol_K-s.entropy_J_mol_K == pytest.approx(solid.FUSION_J_MOL/t,abs=1e-10)


@pytest.mark.parametrize('sp', lp.GAS_SPECIES)
def test_common_ideal_gas(sp):
    t,p,dt,dp = 70.,10000.,1e-3,1.
    c = lp.ideal_gas(sp,t,p)
    tm,tp = lp.ideal_gas(sp,t-dt,p),lp.ideal_gas(sp,t+dt,p)
    pm,pp = lp.ideal_gas(sp,t,p-dp),lp.ideal_gas(sp,t,p+dp)
    assert -(tp.gibbs_J_mol-tm.gibbs_J_mol)/(2*dt) == pytest.approx(c.entropy_J_mol_K,rel=1e-9)
    assert (tp.enthalpy_J_mol-tm.enthalpy_J_mol)/(2*dt) == pytest.approx(c.heat_capacity_J_mol_K,rel=1e-9)
    assert (pp.gibbs_J_mol-pm.gibbs_J_mol)/(2*dp) == pytest.approx(lp.R*t/p,rel=1e-8)


@pytest.mark.parametrize('x', [0.,.1,.5,.9,1.])
def test_ideal_liquid_mixing(x):
    t,p,dt = 70.,101325.,1e-3
    c = lp.ideal_liquid_mixture(t,p,x)
    tm,tp = lp.ideal_liquid_mixture(t-dt,p,x),lp.ideal_liquid_mixture(t+dt,p,x)
    assert -(tp.gibbs_J_mol-tm.gibbs_J_mol)/(2*dt) == pytest.approx(c.entropy_J_mol_K,rel=2e-8)
    assert (tp.enthalpy_J_mol-tm.enthalpy_J_mol)/(2*dt) == pytest.approx(c.heat_capacity_J_mol_K,rel=2e-8)
    mu = lp.liquid_chemical_potentials(t,p,x)
    assert sum(z*m for z,m in zip((x,1-x),mu) if z>0) == pytest.approx(c.gibbs_J_mol,abs=2e-10)
    if x in (0.,1.):
        pure = lp.liquid('Nitrogen' if x==1 else 'Oxygen',t,p)
        assert c == pure


def test_composition_derivatives_and_gibbs_duhem():
    t,p,n,delta = 70.,101325.,np.array([.6,.3]),1e-5
    def total_g(moles):
        return sum(moles)*lp.ideal_liquid_mixture(t,p,moles[0]/sum(moles)).gibbs_J_mol
    mu = lp.liquid_chemical_potentials(t,p,n[0]/sum(n))
    for i in range(2):
        dn = np.zeros(2)
        dn[i] = delta
        assert (total_g(n+dn)-total_g(n-dn))/(2*delta) == pytest.approx(mu[i],rel=1e-10)
    x,dx = .6,1e-5
    mm,mp = np.array(lp.liquid_chemical_potentials(t,p,x-dx)),np.array(lp.liquid_chemical_potentials(t,p,x+dx))
    assert np.dot([x,1-x],(mp-mm)/(2*dx)) == pytest.approx(0.,abs=2e-7)


@pytest.mark.parametrize('t', [64.,70.,77.5])
def test_flash_mu_equality_and_independent_minimum(t):
    p,hi,totals = 101325.,.2,np.array([.65,.15])
    f = lp.flash(t,p,hi,*totals)
    assert f.liquid_mol > 0 and f.gas_mol > 0
    expected = np.array([f.liquid_nitrogen_mol,f.liquid_oxygen_mol])
    mu_l = lp.liquid_chemical_potentials(t,p,expected[0]/sum(expected))
    mu_g = [lp.ideal_gas(sp,t,p*n/f.gas_mol).gibbs_J_mol
            for sp,n in zip(('Nitrogen','Oxygen'),(f.gas_nitrogen_mol,f.gas_oxygen_mol))]
    assert np.array(mu_l) == pytest.approx(mu_g,abs=1e-8)
    def objective(amounts):
        gas = totals-amounts
        vg = hi+sum(gas)
        gg = sum(n*lp.ideal_gas(sp,t,p*n/vg).gibbs_J_mol
                 for sp,n in zip(lp.GAS_SPECIES,[hi,*gas]))
        gl = sum(amounts)*lp.ideal_liquid_mixture(t,p,amounts[0]/sum(amounts)).gibbs_J_mol
        return (gg+gl)/(lp.R*t)
    def gradient(amounts):
        gas = totals-amounts
        vg = hi+sum(gas)
        return (np.array(lp.liquid_chemical_potentials(t,p,amounts[0]/sum(amounts)))
                -np.array([lp.ideal_gas(sp,t,p*n/vg).gibbs_J_mol
                           for sp,n in zip(('Nitrogen','Oxygen'),gas)]))/(lp.R*t)
    result = minimize(objective,totals*.5,jac=gradient,method='SLSQP',
        bounds=[(1e-10,z-1e-10) for z in totals],options={'ftol':1e-13,'maxiter':300})
    assert result.success, result.message
    assert result.x == pytest.approx(expected,abs=1e-7)
    assert f.maximum_inventory_error < 1e-12


def test_flash_scale_and_no_air_limits():
    a,b = lp.flash(70,101325,.2,.65,.15),lp.flash(70,101325,2,6.5,1.5)
    assert b.liquid_mol == pytest.approx(10*a.liquid_mol,rel=1e-12)
    f = lp.flash(70,101325,1,0,0)
    assert f.gas_mol == 1 and f.liquid_mol == 0
    assert lp.flash(100,101325,.2,.65,.15).liquid_mol == 0


def test_full_declared_domain_positive():
    for sp,tmin in lp.TRIPLE_T.items():
        for t in np.linspace(tmin,lp.MAX_T,101):
            for p in (1.,101325.,lp.MAX_P):
                c = lp.liquid(sp,float(t),p)
                assert c.volume_m3_mol > 0 and c.heat_capacity_J_mol_K > 0


def test_latent_uses_ideal_reference_not_real_vapor():
    sp,t,p = 'Nitrogen',70.,101325.
    c = lp.liquid(sp,t,p)
    h0 = lp.ideal_standard(sp,t)[0]
    assert lp.ideal_gas_to_liquid_enthalpy(sp,t,p) == pytest.approx(h0-c.enthalpy_J_mol,abs=1e-12)
    peq = lp.equilibrium_partial_pressure(sp,t,p)
    assert lp.ideal_gas(sp,t,peq).gibbs_J_mol == pytest.approx(c.gibbs_J_mol,abs=1e-9)
