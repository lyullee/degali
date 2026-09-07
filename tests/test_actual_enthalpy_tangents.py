import numpy as np
import pytest

from audit_actual_enthalpy_tangents import caloric_tangent,profile_q
from run_preslhy_ambient_profile_audit import thermodynamics


@pytest.fixture(scope='module')
def th():
    return thermodynamics(dict(T_C=18.,RH_pct=60.,release_height_m=.5),consistent=True)


@pytest.mark.parametrize('t,y',[(285.,.01),(285.,.2),(291.,.01),(291.,.6),(295.,.5)])
def test_gaseous_enthalpy_composition_derivative(th,t,y):
    ht,hy=caloric_tangent(th,t,y)
    expected=th._mixture_enthalpy(t,1.)-th._mixture_enthalpy(t,0.)
    assert hy==pytest.approx(expected,rel=1e-7,abs=1e-7)
    assert ht>0
    fine=caloric_tangent(th,t,y,.5)
    assert ht==pytest.approx(fine[0],rel=1e-3)
    assert hy==pytest.approx(fine[1],rel=1e-6,abs=1e-6)


def test_q_profile_and_species_gradient(th):
    class Model:
        rhoa=th.ambient_density
        thermodynamics=th
    model=Model()
    state=np.array([1.1,.15,.02,0.,25.,1.2,.5])
    for q in (0.,.2,1.,4.):
        rho,c,y,t,h=profile_q(model,state,q)
        dq=1e-5
        _,_,yp,_,_=profile_q(model,state,np.array([q-dq,q+dq]))
        expected=-c[0]*model.rhoa/rho[0]**2
        assert (yp[1]-yp[0])/(2*dq)==pytest.approx(expected,rel=1e-8)
        assert np.all(np.isfinite(t)) and np.all(np.isfinite(h))
