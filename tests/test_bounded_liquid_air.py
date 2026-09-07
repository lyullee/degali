"""Independent checks before the in-domain actual handoff experiment."""
import numpy as np
import pytest

from degali.addons.bounded_liquid_air import BoundedLiquidAir,binary_flash,MIN_T
from degali.addons.mixed_air_liquid import ideal_flash
from degali.addons import liquid_phase_potential as lp
from run_preslhy_ambient_profile_audit import thermodynamics


@pytest.fixture(scope='module')
def adapter():
    return BoundedLiquidAir(thermodynamics(dict(T_C=20.,RH_pct=40.,release_height_m=1.),consistent=True))


@pytest.mark.parametrize('case',[(.2,.65,.15,.1,.02,0),(.2,.65,.15,1,.02,0),
    (0,.65,.35,2,.2,0),(0,.65,.35,.2,.02,0),(.2,.65,.15,2,3,0),
    (1,0,0,1e-20,1e-20,0),(.2,.65,.15,.1,.02,.1),(.2,0,.8,3,.1,0),
    (1e-12,.65,.35,.1,.01,0)])
def test_quadratic_flash_against_independent_bracket(case):
    out = binary_flash(*case)
    if case[1]+case[2]==0:
        assert all(float(a)==0 for a in out)
        return
    f = ideal_flash(*case[:5],reservoir_vapor_fraction=case[5])
    expect = (f.gas_nitrogen_mol,f.gas_oxygen_mol,f.liquid_nitrogen_mol,
              f.liquid_oxygen_mol,f.reservoir_vapor_mol)
    assert np.array(out) == pytest.approx(expect,rel=1e-9,abs=1e-13)


def test_random_flash_vector_and_scale():
    rng = np.random.default_rng(606)
    moles = rng.uniform(.001,3.,(3,300))
    ks = np.exp(rng.uniform(-7,3,(2,300)))
    q = rng.uniform(0,.05,300)
    actual = binary_flash(*moles,*ks,q)
    for i in range(300):
        expected = ideal_flash(*moles[:,i],*ks[:,i],reservoir_vapor_fraction=q[i])
        vals = (expected.gas_nitrogen_mol,expected.gas_oxygen_mol,expected.liquid_nitrogen_mol,
                expected.liquid_oxygen_mol,expected.reservoir_vapor_mol)
        assert np.array(actual)[:,i] == pytest.approx(vals,rel=1e-9,abs=1e-12)
    assert np.array(binary_flash(*(moles*1e6),*ks,q)) == pytest.approx(np.array(actual)*1e6,rel=1e-10,abs=1e-7)


@pytest.mark.parametrize('t',[64.017,66.173,70.231,77.501,90.037,99.981])
def test_potential_interpolation_vs_exact(adapter,t):
    k,h,v = adapter.air_properties(t)
    for i,sp in enumerate(('Nitrogen','Oxygen')):
        assert k[i] == pytest.approx(lp.equilibrium_partial_pressure(sp,t,adapter.ambient_pressure)/adapter.ambient_pressure,rel=1e-10)
        assert h[i] == pytest.approx(lp.ideal_gas_to_liquid_enthalpy(sp,t,adapter.ambient_pressure),rel=2e-9)
        assert v[i] == pytest.approx(lp.liquid(sp,t,adapter.ambient_pressure).volume_m3_mol,rel=1e-10)


def test_density_monotonic_forward_inverse(adapter):
    t = np.linspace(MIN_T,300,801)
    for y in (0.,.01,.1,.3,.8,1.):
        out = adapter.forward(t,y)
        assert np.all(np.diff(out['density'])<0)
        assert np.all(np.diff(out['enthalpy'])>0)
        assert np.all(out['liquid_N2']>=0) and np.all(out['liquid_O2']>=0)
        assert np.all(out['ice_water']>=0)
        tn,rhoh = adapter._condensed_air_state(out['density'],np.full_like(t,y))
        assert tn == pytest.approx(t,abs=1e-8)
        assert rhoh/out['density'] == pytest.approx(out['enthalpy'],abs=2e-5)


def test_legacy_warm_gas_water_parity(adapter):
    t = np.array([101.,120.,160.,220.,280.,293.15])
    for y in (0.,.01,.3,.8):
        out = adapter.forward(t,y)
        ot,oh = adapter.original._condensed_air_state_exact(out['density'],np.full_like(t,y))
        assert ot == pytest.approx(t,abs=1e-7)
        assert oh/out['density'] == pytest.approx(out['enthalpy'],abs=1e-3)


def test_100K_branch_continuity(adapter):
    eps=1e-6
    for y in (0.,.01,.3,1.):
        a,b = adapter.forward(100-eps,y),adapter.forward(100+eps,y)
        assert a['liquid_N2']==a['liquid_O2']==b['liquid_N2']==b['liquid_O2']==0
        assert a['density']==pytest.approx(b['density'],rel=3e-8)
        assert a['enthalpy']==pytest.approx(b['enthalpy'],abs=.04)


def test_domain_guard_is_not_a_phase_splice(adapter):
    with pytest.raises(ValueError): adapter.forward(63.999,.3)
    rho = adapter.forward(64.,.3)['density']
    old_count = adapter.domain_failures
    with pytest.raises(ValueError,match='outside64--300K'):
        adapter._condensed_air_state(rho*1.01,.3)
    assert adapter.domain_failures==old_count+1
    assert adapter.last_domain_failure['fraction']==.3
    with pytest.raises(ValueError): adapter.air_properties(100.001)


def test_ambient_and_original_are_unchanged(adapter):
    assert adapter.original.equilibrium_air_condensation is True
    t,rhoh = adapter._condensed_air_state(adapter.ambient_density,0.)
    assert t==adapter.ambient_temperature and rhoh==0
    assert adapter.original.__class__.__name__=='ConservedGaussianJet'


def test_original_local_caloric_comparison(adapter):
    from audit_liquid_phase_potential import local_state
    for t,y in ((64.2,.3),(66.7,.2),(72.3,.01),(90.1,.3)):
        direct = local_state(adapter.original,t,y)
        actual = adapter.forward(t,y)
        assert actual['density']==pytest.approx(direct['density_kg_m3'],rel=2e-10)
        assert actual['enthalpy']==pytest.approx(direct['specific_enthalpy_J_kg'],abs=2e-5)
