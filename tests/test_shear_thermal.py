"""Shear production, thermal compatibility, and independent moving-face tests."""

import math
from types import SimpleNamespace
import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_transverse_mixing import mixing
from degali.addons.transverse_mixing import ConservativeTransverseMixing, TangentFamily
from degali.addons.phase_radial_quadrature import PhaseRadialQuadrature
from degali.addons.shear_thermal import ReducedShearThermal, affine_root
from degali.addons.thermal_moments import PlanarMovingSection


@pytest.fixture
def shear(mixing):
    return ReducedShearThermal(mixing, axial_force_density=lambda q, d: 9.81*(mixing.section.rhoa-d["rho"])*math.sin(mixing.state[3]))


@pytest.fixture
def family(mixing):
    return mixing.tangent_family(np.array([.1, 0., .2, .03, .4]))


def test_axial_kinetic_product_rule_is_not_a_heat_assumption(shear, family):
    a = shear.axial_fields(np.array([0., .02, .41, 1.3, 2.8]), family)
    u = a["local"]["u"][:, None]
    assert a["kinetic"] == pytest.approx(u*a["momentum"]-.5*u*u*a["mass"], rel=1e-10, abs=1e-9)


def test_reconstructed_momentum_and_kinetic_local_identities(shear, family):
    q, step = np.array([.013, .27, .75, 1.32, 2.12, 3.31]), 2e-6
    def sample(qq):
        return shear.radial_fields(qq, family, thermal_species_ratio=1.)
    d, plus, minus = sample(q), sample(q+step), sample(q-step)
    div = lambda key: 2*d[key]+2*q[:, None]*(plus[key]-minus[key])/(2*step)
    force = np.column_stack([d["force"], np.zeros(len(q))])
    assert d["momentum"]+div("momentum_flux") == pytest.approx(force, rel=2e-4, abs=3e-4)
    assert d["kinetic"]+div("kinetic_flux") == pytest.approx(d["work"]-d["production"], rel=2e-4, abs=2e-3)


def test_mass_and_species_are_identical_to_the_previous_operator(shear, family):
    q = np.array([0., .04, .37, 1.2, 3.4])
    old = shear.mixing.mixing_family(q, family)
    new = shear.radial_fields(q, family, thermal_species_ratio=1.)
    for previous, current in (("mass", "mass_flux"), ("species", "species_flux"), ("chi", "chi_species")):
        assert new[current] == pytest.approx(old[previous], rel=1e-10, abs=1e-10)


def test_kinetic_work_budget_and_energy_discrepancy_are_independent(shear, family):
    b = shear.equilibrium_budgets(family, thermal_species_ratio=1.)
    assert b["kinetic_identity"] == pytest.approx(np.zeros(2), abs=1e-6)
    assert b["total_energy_discrepancy"] == pytest.approx(b["zeroth_residual"], rel=1e-9, abs=1e-6)
    gamma = affine_root(b["second_residual"])
    assert b["second_residual"]@np.array([1., gamma]) == pytest.approx(0., abs=1e-9)


def manufactured_jet():
    """Constant-density self-similar neutral jet; passive heat is NOT Q=P."""
    m = object.__new__(ConservativeTransverseMixing)
    rho, sigma, uc, alpha = 1.3, .4, 5., .12
    m.area, m.sy, m.sn = sigma*sigma, sigma, sigma
    m.state = np.array([rho, .02, m.area, 0., uc, 1., 1.])
    m.section = SimpleNamespace(velocity_shape_exponent=1.)
    m.log_sy_partials = np.array([0., 0., .5, 0., 0., 0., 0., 0.])
    m.log_sn_partials = m.log_sy_partials.copy()
    p = object.__new__(PhaseRadialQuadrature)
    p.check = lambda: None
    p.q0, p.qmax = math.pi*2.15**2/8., math.pi*2.15**2/4.
    p.knots, p._rules = np.array([0., p.q0, p.qmax]), {}
    m.partition = p
    def local(q):
        u, y, h = uc*np.exp(-q), .02*np.exp(-q), -1200*np.exp(-q)
        drho, dc, dh, du = [np.zeros((len(q), 8)) for _ in range(4)]
        drho[:, 0], dc[:, 1], dh[:, 1], du[:, 4] = rho, rho*y, h, u
        return dict(rho=np.full_like(q, rho), u=u, y=y, c=rho*y, h=h, yq=-y, hq=-h/rho,
                    drho=drho, dc=dc, dh=dh, du=du)
    m.local = local
    rates = np.array([0., -alpha, 2*alpha, 0., -alpha, 1., 0., 0.])
    family = TangentFamily(rates, np.zeros(8), np.zeros((5, 8)), np.zeros(5), 0.)
    return ReducedShearThermal(m, axial_force_density=lambda q, d: 0.), family, alpha


def test_manufactured_neutral_jet_recovers_positive_shear_and_production():
    op, family, alpha = manufactured_jet()
    q = np.array([0., .04, .3, 1.2, 3.])
    d = op.radial_fields(q, family, thermal_species_ratio=1.)
    m = op.mixing
    mean = np.ones_like(q)
    mean[1:] = -np.expm1(-q[1:])/q[1:]
    fm = -.5*alpha*m.area*m.state[0]*m.state[4]*mean
    assert d["mass_flux"][:, 0] == pytest.approx(fm, rel=1e-10)
    assert d["momentum_flux"] == pytest.approx(np.zeros((len(q), 2)), abs=1e-13)
    assert d["chi_momentum"][:, 0] == pytest.approx(-fm/(m.area*m.state[0]), rel=1e-10)
    assert np.all(d["chi_momentum"][:, 0] > 0.)
    assert d["production"][0, 0] == 0.
    assert np.all(d["production"][1:, 0] > 0.)
    assert d["enthalpy_flux"] == pytest.approx(np.zeros((len(q), 2)), abs=1e-12)


def test_passive_heat_does_not_silently_gain_all_shear_energy():
    op, family, _ = manufactured_jet()
    b = op.equilibrium_budgets(family, thermal_species_ratio=1., order=32)
    assert b["production"][0] > 0.
    assert b["kinetic_identity"] == pytest.approx(np.zeros(2), abs=1e-10)
    assert b["zeroth_residual"] == pytest.approx(-b["production"], abs=1e-10)
    assert b["second_residual"] == pytest.approx(-b["weighted_production"], abs=1e-10)


def test_negative_shear_is_retained_not_converted_to_positive_heating():
    op, family, _ = manufactured_jet()
    reverse = TangentFamily(-family.origin, family.response, family.jacobian, family.sources, 0.)
    d = op.radial_fields(np.array([.1, .7, 2.]), reverse, thermal_species_ratio=1.)
    assert np.all(d["chi_momentum"][:, 0] < 0.)
    assert np.all(d["production"][:, 0] < 0.)


def test_weighted_radial_budget_matches_independent_curved_moving_four_faces():
    op, family, _ = manufactured_jet()
    m = op.mixing
    m.sy *= 2.
    m.sn *= .5
    rates = family.origin.copy()
    rates[3] = .3
    family = TangentFamily(rates, family.response, family.jacobian, family.sources, 0.)
    # This verifies the weak geometry identity, not the full curved momentum model.
    ratio = 1.17
    def field(y, n):
        q = .5*((y/m.sy)**2+(n/m.sn)**2)
        d = op.radial_fields(q.ravel(), family, thermal_species_ratio=ratio, order=32)
        local = d["local"]
        shape = y.shape
        metric = 1-rates[3]*n
        ay, an = m.log_sy_partials@rates, m.log_sn_partials@rates
        hu = (local["h"]*local["u"]).reshape(shape)
        fh = d["enthalpy_flux"][:, 0].reshape(shape)
        source = d["production"][:, 0].reshape(shape)/(m.area*metric)
        return hu, y/metric*(hu*ay+fh/m.area), n/metric*(hu*an+fh/m.area), source
    geometry = PlanarMovingSection(m.sy, m.sn, 2.15, rates[3],
        m.log_sy_partials@rates, m.log_sn_partials@rates)
    physical = geometry.budget(field, points=40)
    b = op.equilibrium_budgets(family, thermal_species_ratio=ratio, order=32)
    spread_rate = 2*m.sy*m.sy*(m.log_sy_partials@rates)+2*m.sn*m.sn*(m.log_sn_partials@rates)
    iq = m.partition.integrate(lambda q: q*m.area*m.local(q)["h"]*m.local(q)["u"], 32, square=True)
    actual_m2_rate = b["weighted_enthalpy_axial"][0]+spread_rate*iq
    assert actual_m2_rate-physical.derivative == pytest.approx(b["second_residual"][0], rel=1e-9, abs=1e-9)


@pytest.mark.parametrize("pair", [[1., 0.], [0., 0.], [math.nan, 1.], [1., math.inf], [1.]])
def test_undetermined_or_invalid_width_equation_is_rejected(pair):
    with pytest.raises(ValueError):
        affine_root(pair)


def test_force_and_thermal_ratio_must_be_explicit(mixing, shear, family):
    with pytest.raises(TypeError):
        ReducedShearThermal(mixing)
    with pytest.raises(TypeError):
        shear.equilibrium_budgets(family)
    with pytest.raises(ValueError):
        shear.radial_fields(np.array([.1]), family, thermal_species_ratio=-1.)
    with pytest.raises(TypeError):
        ReducedShearThermal(mixing, axial_force_density=0.)
    invalid = ReducedShearThermal(mixing, axial_force_density=lambda q, d: math.nan)
    with pytest.raises(ValueError):
        invalid.axial_fields(np.array([.1]), family)
    assert not hasattr(shear, "solve")
