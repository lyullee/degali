import numpy as np
import pytest
import math
from types import SimpleNamespace
pytest.importorskip('mpmath')
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.stable_enthalpy_difference import scalar_pair,stable_enthalpy_difference,mp_wind
from degali.addons.phase_radial_quadrature import gauss_rule


@pytest.fixture
def synthetic_geometry_only(monkeypatch,transport):
    # The inherited transport fixture deliberately uses a toy wind interface.
    # Test center/enthalpy algebra independently of the real-wind tests below.
    m=transport.projection.mixing
    monkeypatch.setattr('degali.addons.stable_enthalpy_difference.mp_wind',
        lambda section,encoded,mp:(mp.mpf(m.wind),mp.mpf(m.sy),mp.mpf(m.sn)))


def test_high_precision_scalar_values_keep_existing_center(transport,synthetic_geometry_only):
    p=transport.projection
    pair=scalar_pair(p,transport.parameters,np.zeros(transport.count),1e-6)
    assert pair['h']['mean']==pytest.approx(p.hc,rel=2e-14)
    assert pair['h']['difference']==0.
    assert pair['wind']['mean']==pytest.approx(p.mixing.wind,rel=1e-14)


def test_smooth_enthalpy_actual_difference_matches_integrated_derivative(transport,synthetic_geometry_only):
    rng=np.random.default_rng(55); rates=rng.normal(0.,.01,transport.count)
    rates[3]=rates[5]=rates[6]=0.; rates[2]=0. # no geometry policy difference in this unit test
    x,w=gauss_rule(48); a,b=np.meshgrid(x,x,indexing='ij')
    d=transport.local(a.ravel(),b.ravel())
    expected=(8*transport.projection.q0*np.outer(w,w)).ravel()@d['bh']@rates
    out=stable_enthalpy_difference(transport.projection,transport.parameters,rates,1e-6,order=48)
    assert out['derivative']==pytest.approx(expected,rel=1e-9,abs=1e-6)


@pytest.mark.parametrize('rml',[-20.,0.,20.])
def test_high_precision_wind_matches_float_formula_at_exact_width(rml):
    from mpmath import mp
    from degali.core.jetplume import JetPlume
    jp=SimpleNamespace(ustar=.2,zr=.1,rml=rml,k=SimpleNamespace(delta=2.15),
        deltay=.02,betay=.9,deltaz=.01,betaz=.8,gammaz=.01,spread_floor=False)
    encoded=[0.,-2.,math.log(.02),.15,math.log(20.),1.,.5]
    with mp.workdps(50):
        wind,sy,sn=mp_wind(SimpleNamespace(jetplume=jp),list(map(mp.mpf,encoded)),mp)
    expected=JetPlume._wind_profile(jp,.5,float(sn),math.cos(.15))[0]
    assert float(wind)==pytest.approx(expected,rel=2e-14)
    assert float(sy*sn)==pytest.approx(.02,rel=1e-14)
