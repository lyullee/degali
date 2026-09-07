"""Gauge, phase-aware cumulative transport and differential identities."""

import numpy as np
import pytest
from test_edge_enrichment import phase, candidate, section, mixing, reservoir, projection
from degali.addons.enriched_transport import (PrescribedRadialMixing, EnrichedModalTransport,
    panel_cumulative,advective_moments,shifted_advective_moments)
from degali.addons.phase_radial_quadrature import gauss_rule


@pytest.fixture
def transport(projection):
    supplied = PrescribedRadialMixing.from_weak_baseline(projection, order=8)
    return EnrichedModalTransport(projection, np.zeros(projection.count), scalar_mixing=supplied,
        thermal_species_ratio=1.,mechanical_work="reduced_buoyancy_work_immediate_shear_heat")


@pytest.mark.parametrize("order",[4,8,16])
def test_panel_antiderivative_recovers_polynomials(order):
    knots=np.array([0.,.04,.3,1.])
    x,_=gauss_rule(order)
    q=(knots[:-1,None]+np.diff(knots)[:,None]*x).ravel()
    values=np.column_stack([q**j for j in range(order)])
    actual,end=panel_cumulative(values,knots,order)
    expected=np.column_stack([q**(j+1)/(j+1) for j in range(order)])
    assert actual == pytest.approx(expected,abs=2e-13)
    assert end == pytest.approx(1./np.arange(1,order+1),abs=2e-13)


def test_thermal_width_change_is_a_shape_gauge_not_a_new_state(transport):
    p=transport.projection
    a,b=np.array([0.,.2,.6,1.]),np.array([0.,.4,.1,.7])
    d=transport.local(a,b)
    assert d["psi"]@transport.q_mode == pytest.approx(d["q"],abs=1e-12)
    old=np.zeros(8)
    old[7]=.17
    rates=transport.old_rate_embedding(old)
    assert d["dh"]@rates == pytest.approx(2*d["q"]*d["h"]*.17/p.section.thermal_width_ratio**2,rel=1e-12,abs=1e-8)
    assert transport.count == 23


def test_new_scalar_partials_match_independent_shape_change(transport):
    p=transport.projection
    a,b=np.array([.15,.44,.71]),np.array([.24,.55,.84])
    d=transport.local(a,b)
    direction=np.random.default_rng(3).normal(0.,.1,p.count)
    step=1e-6
    plus=p.fields(transport.parameters+step*direction,p.prepare(a,b))
    minus=p.fields(transport.parameters-step*direction,p.prepare(a,b))
    for key,derivative in [("rho","drho"),("c","dc"),("h","dh")]:
        assert d[derivative][:,7:]@direction == pytest.approx((plus[key]-minus[key])/(2*step),rel=2e-5,abs=1e-5)


def test_mass_momentum_kinetic_product_identity(transport):
    d=transport.local(np.array([.1,.3,.7]),np.array([.2,.4,.8]))
    assert d["bk"] == pytest.approx(d["u"][:,None]*d["bp"]-.5*d["u"][:,None]**2*d["bm"],rel=1e-12,abs=1e-8)


def test_scalar_gradient_has_angular_component(transport):
    transport.parameters[2]=1e-3
    a,b=np.array([.3,.7]),np.array([.6,.2])
    d=transport.local(a,b)
    angular=-b*d["grad_y"][:,0]+a*d["grad_y"][:,1]
    assert max(abs(angular))>1e-7
    step=1e-6
    plus=transport.local(a+step,b)["y"]
    minus=transport.local(a-step,b)["y"]
    assert d["grad_y"][:,0] == pytest.approx((plus-minus)/(2*step),rel=1e-6,abs=1e-8)


def test_explicit_mixing_preserves_phase_panel_jumps():
    knots=np.array([0.,.2,1.])
    values=np.zeros((2,8,2))
    values[0,:,0],values[1,:,0]=1.,2.
    values[:,:,1]=-1.
    model=PrescribedRadialMixing(knots,values)
    assert model.evaluate(np.array([.2-1e-10,.2,.2+1e-10]))[:,0] == pytest.approx([1.,2.,2.])
    values[0,:,0]=-1.
    with pytest.raises(ValueError,match="positive"):
        PrescribedRadialMixing(knots,values)


def test_omitted_closure_or_energy_assumption_is_rejected(projection,transport):
    with pytest.raises(TypeError):
        EnrichedModalTransport(projection,np.zeros(projection.count),scalar_mixing=None,
            thermal_species_ratio=1.,mechanical_work="reduced_buoyancy_work_immediate_shear_heat")
    with pytest.raises(ValueError):
        EnrichedModalTransport(projection,np.zeros(projection.count),scalar_mixing=transport.mixing,
            thermal_species_ratio=1.,mechanical_work="unselected")


def test_angle_dependent_radial_flux_still_has_conservative_divergence():
    def flux(a,b):
        q=a*a+b*b
        phi=np.arctan2(b,a)
        f=-.5*(1.+.5*q*np.cos(4*phi)+q*q/3*np.cos(8*phi))
        return np.array([a*f,b*f])
    a,b,step=.4,.7,1e-6
    divergence=(flux(a+step,b)[0]-flux(a-step,b)[0]+flux(a,b+step)[1]-flux(a,b-step)[1])/(2*step)
    q,phi=a*a+b*b,np.arctan2(b,a)
    assert divergence == pytest.approx(-1.-q*np.cos(4*phi)-q*q*np.cos(8*phi),rel=1e-9)


def test_physical_velocity_adds_deforming_coordinate_motion(transport):
    a,b=np.array([.2,.7]),np.array([.1,.6])
    rates=np.zeros(transport.count)
    rates[2]=.2
    d=transport.local(a,b)
    # Isotropic fixture: both widths grow at half the area rate.
    mass=-transport.area*d["rho"]*d["u"]*.1
    vy,vn=transport.physical_velocity(a,b,rates,mass)
    assert vy == pytest.approx(np.zeros(2),abs=1e-10)
    assert vn == pytest.approx(np.zeros(2),abs=1e-10)


def test_coupled_weak_rates_and_actual_moment_direction(transport):
    out=transport.assemble(order=4,angular_order=4)
    assert out["linear_scaled_error"] < 1e-8
    assert out["weak_heat_scaled_error"] < 1e-8
    assert abs(out["mass_boundary_error"]) < 1e-8
    assert max(abs(out["rates"][7:])) > 1e-4
    direction=np.random.default_rng(23).normal(0.,.01,transport.count)
    ds=1e-5
    plus=shifted_advective_moments(transport.projection,transport.parameters,direction,ds,order=4,angular_order=4)
    minus=shifted_advective_moments(transport.projection,transport.parameters,direction,-ds,order=4,angular_order=4)
    scale=np.maximum(abs(out["fluxes"]),1.)
    assert max(abs((plus-minus)/(2*ds)-out["moment_jacobian"]@direction)/scale) < 2e-5
