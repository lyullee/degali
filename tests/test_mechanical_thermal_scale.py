import math
from types import SimpleNamespace

import numpy as np
import pytest

from audit_mechanical_thermal_scale import velocity_moment,legacy_forward
from run_preslhy_ambient_profile_audit import thermodynamics


class ManufacturedGaussian:
    def __init__(self,ambient_density,wind,lam,cutoff):
        self.rhoa=ambient_density
        self.wind=wind
        self.velocity_shape_exponent=lam
        self.cutoff=cutoff
        self.k=SimpleNamespace(profile_integral=self.integral)
    def integral(self,power):
        if power==0: return 4*self.cutoff**2
        return math.pi/power*math.erf(self.cutoff*math.sqrt(power))**2
    def _physical(self,state): return state
    def _wind(self,state): return self.wind


@pytest.mark.parametrize('power',[1,2,3])
@pytest.mark.parametrize('wind',[0.,2.])
@pytest.mark.parametrize('rho',[.2,2.])
def test_analytic_velocity_moment_against_tensor_quadrature(power,wind,rho):
    model=ManufacturedGaussian(1.1,wind,1.16**2,2.1)
    state=np.array([rho,.1,.08,.2,30.,1.,.5])
    node,weight=np.polynomial.legendre.leggauss(96)
    xx,yy=np.meshgrid(node*model.cutoff,node*model.cutoff,indexing='ij')
    q=np.exp(-(xx*xx+yy*yy))
    density=model.rhoa+(rho-model.rhoa)*q
    velocity=wind*math.cos(state[3])+state[4]*q**model.velocity_shape_exponent
    expected=state[2]*model.cutoff**2*np.sum(density*velocity**power*weight[:,None]*weight[None,:])
    assert velocity_moment(model,state,power)==pytest.approx(expected,rel=2e-13)


@pytest.fixture(scope='module')
def th():
    return thermodynamics(dict(T_C=18.,RH_pct=60.,release_height_m=.5),consistent=True)


@pytest.mark.parametrize('t,y',[(30.,.6),(50.,.3),(60.,.28),(63.151,.27),(70.,.2),(100.,.1),(200.,.05),(291.,.001),(300.,.1),(100.,1.)])
def test_forward_old_eos_exact_inverse(th,t,y):
    state=legacy_forward(th,t,y)
    rho=state['density_kg_m3']
    recovered,hv=th._condensed_air_state_exact(np.array([rho]),np.array([y]))
    assert recovered[0]==pytest.approx(t,abs=8e-8)
    # Original exact inverse has32bisections, giving a finite caloric roundoff.
    assert hv[0]/rho==pytest.approx(state['specific_enthalpy_J_kg'],abs=.002)


@pytest.mark.parametrize('t,y',[(14.,.1),(301.,.1),(100.,0.),(100.,1.1),(float('nan'),.1)])
def test_outside_diagnostic_domain(th,t,y):
    with pytest.raises(ValueError): legacy_forward(th,t,y)
