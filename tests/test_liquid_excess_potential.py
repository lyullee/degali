"""Independent potential identities, phase stability and total-G minimization."""
import math

import numpy as np
import pytest
from scipy.optimize import minimize

from degali.addons import liquid_excess_potential as ep
from degali.addons import liquid_phase_potential as lp
from degali.addons.mixed_air_liquid import nbs_activity_coefficients


@pytest.mark.parametrize('t',[64.999,77.501,float('nan'),float('inf')])
def test_no_temperature_extrapolation(t):
    with pytest.raises(ValueError): ep.excess(t,.5,interpolation='quadratic')


@pytest.mark.parametrize('x',[-.01,1.01,float('nan')])
def test_composition_guard(x):
    with pytest.raises(ValueError): ep.excess(70,x,interpolation='quadratic')


def test_interpolation_is_explicit():
    with pytest.raises(TypeError): ep.excess(70,.5)
    with pytest.raises(ValueError): ep.excess(70,.5,interpolation='extrapolated')


@pytest.mark.parametrize('method',['quadratic','pchip'])
@pytest.mark.parametrize('t',[65.,70.,77.5])
def test_exact_original_isotherm_activities(method,t):
    for x in (0.,.01,.25,.5,.75,.99,1.):
        assert ep.activity_coefficients(t,x,interpolation=method) == pytest.approx(nbs_activity_coefficients(t,x),rel=1e-14)


@pytest.mark.parametrize('method',['quadratic','pchip'])
@pytest.mark.parametrize('t',[65.2,68.,72.,77.2])
def test_temperature_and_pressure_derivatives(method,t):
    p,x,dt,dp = 101325.,.63,1e-3,100.
    e = ep.excess(t,x,interpolation=method)
    em,ee = [ep.excess(t+delta,x,interpolation=method) for delta in (-dt,dt)]
    assert e.gibbs_J_mol+t*e.entropy_J_mol_K == pytest.approx(e.enthalpy_J_mol,abs=1e-13)
    assert -(ee.gibbs_J_mol-em.gibbs_J_mol)/(2*dt) == pytest.approx(e.entropy_J_mol_K,rel=3e-7)
    assert (ee.enthalpy_J_mol-em.enthalpy_J_mol)/(2*dt) == pytest.approx(e.heat_capacity_J_mol_K,rel=3e-7)
    c = ep.liquid_mixture(t,p,x,interpolation=method)
    cm,cp = [ep.liquid_mixture(t+delta,p,x,interpolation=method) for delta in (-dt,dt)]
    pm,pp = [ep.liquid_mixture(t,p+delta,x,interpolation=method) for delta in (-dp,dp)]
    assert c.gibbs_J_mol+t*c.entropy_J_mol_K == pytest.approx(c.enthalpy_J_mol,abs=1e-9)
    assert -(cp.gibbs_J_mol-cm.gibbs_J_mol)/(2*dt) == pytest.approx(c.entropy_J_mol_K,rel=3e-8)
    assert (cp.enthalpy_J_mol-cm.enthalpy_J_mol)/(2*dt) == pytest.approx(c.heat_capacity_J_mol_K,rel=3e-8)
    assert (pp.gibbs_J_mol-pm.gibbs_J_mol)/(2*dp) == pytest.approx(c.volume_m3_mol,rel=3e-8)


@pytest.mark.parametrize('method',['quadratic','pchip'])
def test_composition_euler_gibbs_duhem_and_curvature(method):
    t,p,n,dn = 69.,101325.,np.array([.7,.2]),1e-5
    x=n[0]/sum(n)
    mu=ep.chemical_potentials(t,p,x,interpolation=method)
    g=ep.liquid_mixture(t,p,x,interpolation=method).gibbs_J_mol
    assert np.dot([x,1-x],mu) == pytest.approx(g,abs=2e-10)
    def total_g(amounts):
        return sum(amounts)*ep.liquid_mixture(t,p,amounts[0]/sum(amounts),interpolation=method).gibbs_J_mol
    for i in range(2):
        delta=np.zeros(2);delta[i]=dn
        assert (total_g(n+delta)-total_g(n-delta))/(2*dn) == pytest.approx(mu[i],rel=2e-9)
    mm,mp=[np.array(ep.chemical_potentials(t,p,x+dx,interpolation=method)) for dx in (-dn,dn)]
    deriv=(mp-mm)/(2*dn)
    assert np.dot([x,1-x],deriv) == pytest.approx(0.,abs=2e-6)
    assert deriv[0]-deriv[1] == pytest.approx(ep.composition_curvature(t,x,interpolation=method),rel=2e-8)


@pytest.mark.parametrize('x',[0.,1.])
def test_pure_limit(x):
    e=ep.excess(70,x,interpolation='quadratic')
    assert (e.gibbs_J_mol,e.enthalpy_J_mol,e.entropy_J_mol_K,e.heat_capacity_J_mol_K)==(0.,0.,0.,0.)
    assert ep.liquid_mixture(70,101325,x,interpolation='quadratic') == lp.liquid('Nitrogen' if x else 'Oxygen',70,101325)


def test_domain_convexity_and_heat_capacity():
    for method in ('quadratic','pchip'):
        for t in np.linspace(65,77.5,26):
            for x in np.linspace(0,1,31):
                assert ep.composition_curvature(t,x,interpolation=method)>0
                for p in (1.,200000.):
                    assert ep.liquid_mixture(t,p,x,interpolation=method).heat_capacity_J_mol_K>0


def test_interpolation_sensitivity_not_hidden():
    # PCHIP is C1, NOT a calorimetrically established C2 function at70K.
    left=ep.excess(70-1e-7,.5,interpolation='pchip')
    right=ep.excess(70+1e-7,.5,interpolation='pchip')
    assert abs(left.enthalpy_J_mol-right.enthalpy_J_mol)<1e-6
    assert abs(left.heat_capacity_J_mol_K-right.heat_capacity_J_mol_K)>1e-4
    assert ep.excess(70,.5,interpolation='quadratic').enthalpy_J_mol != ep.excess(70,.5,interpolation='pchip').enthalpy_J_mol


@pytest.mark.parametrize('t',[65.,70.,77.5])
@pytest.mark.parametrize('method',['quadratic','pchip'])
def test_flash_mu_equality_independent_total_g_minimum(t,method):
    p,hi,totals=101325.,.2,np.array([.65,.15])
    f=ep.flash(t,p,hi,*totals,interpolation=method)
    a=f.amounts
    expected=np.array([a.liquid_nitrogen_mol,a.liquid_oxygen_mol])
    assert sum(expected)>0
    assert f.minimum_liquid_tangent_distance_RT<0
    assert f.chemical_equilibrium_error_RT<1e-10
    assert a.maximum_inventory_error<1e-12
    def objective(amounts):
        gas=totals-amounts;vg=hi+sum(gas)
        gg=sum(n*lp.ideal_gas(sp,t,p*n/vg).gibbs_J_mol for sp,n in zip(lp.GAS_SPECIES,[hi,*gas]))
        gl=sum(amounts)*ep.liquid_mixture(t,p,amounts[0]/sum(amounts),interpolation=method).gibbs_J_mol
        return (gg+gl)/(lp.R*t)
    def gradient(amounts):
        gas=totals-amounts;vg=hi+sum(gas)
        return (np.array(ep.chemical_potentials(t,p,amounts[0]/sum(amounts),interpolation=method))
            -np.array([lp.ideal_gas(sp,t,p*n/vg).gibbs_J_mol for sp,n in zip(('Nitrogen','Oxygen'),gas)]))/(lp.R*t)
    result=minimize(objective,totals*.5,jac=gradient,method='SLSQP',bounds=[(1e-10,z-1e-10) for z in totals],options={'ftol':1e-13,'maxiter':300})
    assert result.success,result.message
    assert result.x==pytest.approx(expected,abs=1e-7)
    assert objective(expected)<=objective(result.x)+1e-12


def test_flash_scale_gas_and_pure_air_limits():
    f=ep.flash(70,101325,.2,.65,.15,interpolation='quadratic').amounts
    ff=ep.flash(70,101325,2,6.5,1.5,interpolation='quadratic').amounts
    assert ff.liquid_mol==pytest.approx(f.liquid_mol*10,rel=1e-12)
    g=ep.flash(77.5,101325,1000,.65,.15,interpolation='quadratic')
    assert g.minimum_liquid_tangent_distance_RT>0 and g.amounts.liquid_mol==0
    for amounts in ((1,0,0),(.2,.8,0),(.2,0,.8)):
        actual=ep.flash(70,101325,*amounts,interpolation='quadratic').amounts
        old=lp.flash(70,101325,*amounts)
        assert actual==old
    with pytest.raises(ValueError): ep.flash(70,101325,0,.8,.2,interpolation='quadratic')


def test_reservoir_and_material_balance():
    f=ep.flash(70,101325,.3,.5,.2,interpolation='quadratic',reservoir_vapor_fraction=.01)
    a=f.amounts
    assert a.reservoir_vapor_mol/a.gas_mol==pytest.approx(.01,rel=1e-12)
    assert a.gas_nitrogen_mol+a.liquid_nitrogen_mol==pytest.approx(.5,abs=1e-13)
    assert a.gas_oxygen_mol+a.liquid_oxygen_mol==pytest.approx(.2,abs=1e-13)
    assert f.chemical_equilibrium_error_RT<1e-10


@pytest.mark.parametrize('args',[(70,0,.2,.6,.2),(70,200001,.2,.6,.2),(70,101325,-1,.6,.2),(70,101325,.2,-.1,.2)])
def test_flash_invalid_inputs(args):
    with pytest.raises(ValueError): ep.flash(*args,interpolation='quadratic')
