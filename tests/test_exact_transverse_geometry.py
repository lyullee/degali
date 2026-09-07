import math
from types import SimpleNamespace
import numpy as np
import pytest
from degali.addons.exact_transverse_geometry import exact_split,exact_geometry_tangent,ExactConstraintGeometry
from degali.addons.stable_enthalpy_difference import mp_wind


def model(rml=0.,ustar=.2,spread_floor=False):
    return SimpleNamespace(ustar=ustar,zr=.1,rml=rml,k=SimpleNamespace(delta=2.15),
        deltay=.02,betay=.9,deltaz=.01,betaz=.8,gammaz=.01,spread_floor=spread_floor)


@pytest.mark.parametrize('sya,sza',[(.1,.2),(.2,.1),(0.,0.),(100.,.001)])
def test_exact_split_recovers_area_and_excess_spread(sya,sza):
    sy,sn=exact_split(.02,sya,sza)
    assert (sy*sn).real==pytest.approx(.02,rel=1e-14)
    assert (sy*sy-sn*sn).real==pytest.approx(sya*sya-sza*sza,rel=1e-13,abs=1e-14)


@pytest.mark.parametrize('rml',[-20.,0.,20.])
def test_complex_step_geometry_matches_independent_high_precision(rml):
    mp=pytest.importorskip('mpmath').mp
    jp=model(rml); state=np.array([1.,.1,.02,.15,20.,1.,.5])
    out=exact_geometry_tangent(jp,state)
    with mp.workdps(60):
        point=[mp.mpf(math.log(.02)),mp.mpf(.15),mp.mpf(1.),mp.mpf(.5)]
        def value(logarea,theta,x,z):
            wind,sy,sn=mp_wind(SimpleNamespace(jetplume=jp),[0.,-2.,logarea,theta,3.,x,z],mp)
            return sy,sn,wind
        expected=np.array([float(v) for v in value(*point)])
        jac=np.array([[float(mp.diff(lambda *v:value(*v)[row],point,tuple(int(j==col) for j in range(4))))
                       for col in range(4)] for row in range(3)])
    assert out['values']==pytest.approx(expected,rel=2e-14,abs=1e-14)
    assert out['jacobian']==pytest.approx(jac,rel=1e-10,abs=1e-12)


def test_quiescent_geometry_still_preserves_nonisotropic_spreads():
    jp=model(ustar=0.)
    out=exact_geometry_tangent(jp,[1.,.1,.02,.15,20.,1.,.5])
    sy,sn,wind=out['values']
    assert sy!=sn and sy*sn==pytest.approx(.02)
    assert wind==0. and np.all(out['jacobian'][2]==0.)


def test_wrapper_does_not_mutate_original_geometry():
    jp=model(spread_floor=True); before=vars(jp).copy()
    wrapper=ExactConstraintGeometry(jp)
    sy,sn=wrapper._split(.001,.2,.05)
    assert sy==.2 and sy*sn==pytest.approx(.001)
    assert vars(jp)==before and not hasattr(jp,'_split')
