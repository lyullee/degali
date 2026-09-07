import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.axisymmetric_jet import AxisymmetricJetSource, ConservedGaussianJet
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.yawed_crosswind import YawedCrosswind, YawedTrajectory, direction, normal_drag
from degali.core.jetplume import JetCoefficients


def base_model(wind=2., energy='total', ground='free'):
    jp = SimpleNamespace(th=None,k=JetCoefficients(sc=1.16**2),
        deltay=.1,deltaz=.08,betay=.92,betaz=.87,gammaz=-.012,ustar=.2,
        _split=lambda area,*_: (math.sqrt(area),math.sqrt(area)),
        _wind=lambda *_:wind,_wind_profile=lambda *_:(wind,0.))
    source = AxisymmetricJetSource(diameter=.01,velocity=100.,density=.6,temperature=50.)
    th = ConservedGaussianJet(source,ambient_temperature=295.,ambient_pressure=101325.,
        ambient_density=1.197,fuel_molecular_weight=.00201588,ambient_molecular_weight=.02896546,
        fuel_heat_capacity=14294.8,ambient_heat_capacity=1006.2,radial_points=41)
    return IndependentEnergyCrosswind(jp,th,quadrature_points=32,
        energy_transport=energy,ground_interaction=ground)


def old_state(pitch=.1):
    return np.array([1.1,.02,.03,pitch,15.,2.,1.5])


def lift_flux(f):
    return np.array([f[0],f[1],f[2],0.,f[3],f[4]])


@pytest.mark.parametrize('pitch',[-.25,0.,.25])
@pytest.mark.parametrize('energy',['total','enthalpy'])
def test_coplanar_flux_and_sources_recover_original(pitch,energy):
    base = base_model(energy=energy)
    candidate = YawedCrosswind(base,0.)
    q = old_state(pitch)
    state = candidate.lift(q)
    assert candidate.fluxes(state) == pytest.approx(lift_flux(base._as_array(base.integral_fluxes(q))),rel=2e-14,abs=1e-12)
    assert candidate.sources(state) == pytest.approx(lift_flux(base.source_terms(q)),rel=2e-14,abs=1e-12)


@pytest.mark.parametrize('ground',['geometry','surface_layer'])
def test_ground_handling_coplanar_limit(ground):
    base = base_model(ground=ground)
    candidate = YawedCrosswind(base,0.)
    old = old_state(.15)
    old[6] = .15
    assert candidate.sources(candidate.lift(old)) == pytest.approx(lift_flux(base.source_terms(old)),rel=2e-14,abs=1e-12)


@pytest.mark.parametrize('angle',[-.8,.43,1.7])
def test_horizontal_rotation_covariance(angle):
    base = base_model()
    initial = YawedCrosswind.lift(old_state(),yaw=.2)
    original = YawedCrosswind(base,.4)
    rotated = YawedCrosswind(base,.4+angle)
    state = YawedCrosswind.lift(old_state(),yaw=.2+angle)
    c,s = math.cos(angle),math.sin(angle)
    rot = np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
    for method in ('fluxes','sources'):
        a = getattr(original,method)(initial)
        b = getattr(rotated,method)(state)
        expected = a.copy(); expected[2:5] = rot@a[2:5]
        assert b == pytest.approx(expected,rel=2e-13,abs=1e-11)


def test_mirror_wind_mirrors_transverse_momentum():
    base = base_model()
    a,b = YawedCrosswind(base,.35),YawedCrosswind(base,-.35)
    sa,sb = a.lift(old_state(),.1),b.lift(old_state(),-.1)
    for method in ('sources','fluxes'):
        expected = getattr(a,method)(sa).copy(); expected[3] *= -1
        assert getattr(b,method)(sb) == pytest.approx(expected,rel=2e-14,abs=1e-11)


def test_drag_is_normal_and_rotates():
    n = direction(.3,.2)
    wind = np.array([1.,2.,0.])
    force = normal_drag(wind,n,.7)
    assert float(force@n) == pytest.approx(0.,abs=1e-14)
    assert normal_drag(np.zeros(3),n,.7) == pytest.approx(np.zeros(3))
    with pytest.raises(ValueError): normal_drag(wind,n*2,.7)


def test_six_flux_inversion_and_base_not_mutated():
    base = base_model()
    before = base._as_array(base.integral_fluxes(old_state()))
    candidate = YawedCrosswind(base,.35)
    state = candidate.lift(old_state(),.12)
    target = candidate.fluxes(state)
    guess = state.copy();guess[2] *= 1.01;guess[5] *= .99
    fitted = candidate.match(target,guess)
    assert candidate.fluxes(fitted) == pytest.approx(target,rel=1e-9,abs=1e-9)
    assert fitted == pytest.approx(state,rel=2e-7,abs=1e-8)
    assert np.array_equal(base._as_array(base.integral_fluxes(old_state())),before)
    assert '_wind' not in vars(base)


def test_short_coplanar_rk4_replays_original():
    base = base_model()
    old = old_state()
    expected = base.solve(old,maximum_distance=.04,maximum_step=.02)
    candidate = YawedCrosswind(base,0.)
    result = candidate.solve(candidate.lift(old),distance=.04,step=.02)
    actual = np.array([candidate.proxy(s) for s in result['states']])
    assert actual == pytest.approx(expected.states,rel=5e-9,abs=1e-9)
    assert result['states'][:,4] == pytest.approx(0.,abs=1e-14)
    assert result['states'][:,8] == pytest.approx(0.,abs=1e-14)
    assert result['maximum_relative_balance_residual'] < 1e-7


def test_receptor_zero_yaw_limit_and_no_extrapolation():
    base = base_model()
    model = YawedCrosswind(base,0.)
    left,right = old_state(),old_state()
    left[5],right[5] = 1.,3.
    tr = YawedTrajectory(model,[model.lift(left),model.lift(right)])
    for y in (-1.,0.,1.):
        assert tr.temperature_at(2.,y,1.2) == pytest.approx(base.point_temperature(old_state(),y,1.2),abs=1e-11)
        assert tr.concentration_at(2.,y,1.2) == pytest.approx(base.point_mole_fraction(old_state(),y,1.2),abs=1e-14)
    with pytest.raises(ValueError): tr.section(4.,0.)


def test_receptor_horizontal_rotation_covariance():
    base = base_model();angle=.45
    a,b = YawedCrosswind(base,0.),YawedCrosswind(base,angle)
    left,right = old_state(),old_state();left[5],right[5] = 1.,3.
    ta = YawedTrajectory(a,[a.lift(left),a.lift(right)])
    tb = YawedTrajectory(b,[b.lift(left,angle),b.lift(right,angle)])
    x,y = 2.,.2
    xp,yp = x*math.cos(angle)-y*math.sin(angle),x*math.sin(angle)+y*math.cos(angle)
    assert tb.temperature_at(xp,yp,1.2) == pytest.approx(ta.temperature_at(x,y,1.2),abs=1e-10)


@pytest.mark.parametrize('angle',[math.pi,math.radians(105),-math.pi])
def test_reverse_ambient_is_rejected(angle):
    model = YawedCrosswind(base_model(),angle)
    with pytest.raises(ValueError,match='reverse ambient'):
        model.fluxes(model.lift(old_state()))


@pytest.mark.parametrize('angle',[.5,math.pi])
def test_zero_wind_flux_and_sources_independent_of_wind_bearing(angle):
    base = base_model(wind=0.)
    a,b = YawedCrosswind(base,0.),YawedCrosswind(base,angle)
    state = a.lift(old_state())
    assert a.fluxes(state) == pytest.approx(b.fluxes(state),rel=1e-14)
    assert a.sources(state) == pytest.approx(b.sources(state),rel=1e-14)
