"""Conservative transverse reconstruction, not a fitted thermal-width ODE."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_thermal_moments import STATE
from degali.addons.phase_radial_quadrature import PhaseRadialQuadrature
from degali.addons.transverse_mixing import (
    ConservativeTransverseMixing, TangentFamily, nonnegative_affine_interval,
)


@pytest.fixture
def mixing(section):
    section.thermal_width_ratio = 1.04
    return ConservativeTransverseMixing(section, STATE)


def test_partition_retains_every_monotone_phase_cell(mixing):
    p = mixing.partition
    assert len(p.knots) > 20
    q, _ = p.nodes(8)
    ti, y = p.coordinates(q.ravel())
    inv = mixing.section.phase_inverse
    ti_cells = np.searchsorted(inv.t_grid, ti).reshape(q.shape)
    y_cells = np.searchsorted(inv.y_grid, y).reshape(q.shape)
    assert np.all(ti_cells == ti_cells[:, :1])
    assert np.all(y_cells == y_cells[:, :1])
    other = PhaseRadialQuadrature(mixing.section, STATE, probes=1025)
    assert other.knots == pytest.approx(p.knots, abs=1e-11)


def test_split_radial_square_and_edge_measures_are_exact(mixing):
    p = mixing.partition
    power = .83
    expected = 2*math.pi/power*math.erf(math.sqrt(power*p.q0))**2
    for order in (8, 16):
        actual = p.integrate(lambda q: np.exp(-power*q), order, square=True)
        assert actual == pytest.approx(expected, rel=1e-10)
        assert p.integrate(lambda q: np.ones_like(q), order, square=True) == pytest.approx(8*p.q0, rel=1e-12)
        assert p.edge_integrate(lambda q: q, order) == pytest.approx(4*p.q0/3., rel=1e-12)
    q = np.array([0., .01, .3, 1.7, p.qmax])
    actual = p.cumulative(lambda q: np.column_stack([np.exp(-power*q), q*q]), q)
    expected = np.column_stack([-np.expm1(-power*q)/power, q**3/3.])
    assert actual == pytest.approx(expected, rel=1e-10, abs=1e-12)


def test_full_flux_jacobian_matches_primitive_state_difference(mixing):
    model, state = mixing.section, mixing.state
    analytic = mixing.flux_jacobian()
    # Use the common fixed phase-split integration nodes for finite differences.
    q, weights = mixing.partition.nodes(16, square=True)
    q, weights = q.ravel(), weights.ravel()
    params = np.array([math.log(state[0]), math.log(state[0]*state[1]), math.log(state[2]), state[3],
                       math.log(state[4]), state[5], state[6], math.log(model.thermal_width_ratio)])
    def evaluate(par):
        rho, c, area, theta, uc, x, z, beta = par.copy()
        rho, c, area, uc, beta = np.exp([rho, c, area, uc, beta])
        model.thermal_width_ratio = beta
        trial = np.array([rho, c/rho, area, theta, uc, x, z])
        rr, yy, _, hh = model.thermodynamic_profile(trial, np.exp(-q))
        u = model._wind(trial)*math.cos(theta)+uc*np.exp(-model.velocity_shape_exponent*q)
        return area*np.array([np.sum(weights*rr*u), np.sum(weights*rr*yy*u),
            np.sum(weights*rr*u*u)*math.cos(theta), np.sum(weights*rr*u*u)*math.sin(theta),
            np.sum(weights*(hh*u+.5*rr*u**3))])
    step = 1e-7
    numerical = np.column_stack([(evaluate(params+step*e)-evaluate(params-step*e))/(2*step) for e in np.eye(8)])
    model.thermal_width_ratio = math.exp(params[-1])
    assert analytic == pytest.approx(numerical, rel=1e-4, abs=2e-4)


def test_tangent_retains_free_width_rate_and_five_fluxes(mixing):
    source = np.array([.1, 0., .2, .03, .4])
    family = mixing.tangent_family(source)
    assert family.maximum_scaled_residual < 1e-8
    for gamma in (-.2, 0., .3):
        rates = family.at(gamma)
        assert family.jacobian@rates == pytest.approx(source, abs=1e-8)
        assert rates[7] == gamma
        assert rates[5:7] == pytest.approx([math.cos(STATE[3]), math.sin(STATE[3])])
    assert not hasattr(mixing, "solve")


def test_radial_reconstruction_satisfies_local_mass_and_species_divergence(mixing):
    q, step = np.array([.013, .27, .75, 1.32, 2.12, 3.31]), 1e-5
    f = mixing.transverse_basis(q)
    df = (mixing.transverse_basis(q+step)-mixing.transverse_basis(q-step))/(2*step)
    divergence = 2*(f+q[:, None, None]*df)
    expected = -mixing.axial_rate_density(q)
    assert divergence == pytest.approx(expected, rel=2e-4, abs=2e-5)
    centre = mixing.transverse_basis(np.array([0.]))
    assert centre == pytest.approx(-.5*mixing.axial_rate_density(np.array([0.])), rel=1e-12)


def manufactured_diffusion():
    """Constant-density exact Gaussian heat/species diffusion, in moving axes."""
    rho, u, diffusion, sigma = 1.3, 4., .07, .4
    area, alpha = sigma*sigma, diffusion/(u*sigma*sigma)
    obj = object.__new__(ConservativeTransverseMixing)
    obj.area, obj.sy, obj.sn = area, sigma, sigma
    p = object.__new__(PhaseRadialQuadrature)
    p.check = lambda: None
    p.q0, p.qmax = math.pi*2.15**2/8., math.pi*2.15**2/4.
    p.knots, p._rules = np.array([0., p.q0, p.qmax]), {}
    obj.partition = p
    obj.log_sy_partials = np.array([0., 0., .5, 0., 0., 0., 0., 0.])
    obj.log_sn_partials = obj.log_sy_partials.copy()
    def local(q):
        y = .02*np.exp(-q)
        h = -1200.*np.exp(-q)
        return dict(rho=np.full_like(q, rho), u=np.full_like(q, u), y=y, c=rho*y,
                    h=h, yq=-y, hq=-h/rho)
    obj.local = local
    def axial(q):
        result = np.zeros((len(q), 2, 8))
        result[:, 0, 2] = area*rho*u
        result[:, 1, 1] = area*rho*.02*np.exp(-q)*u
        result[:, 1, 2] = result[:, 1, 1]
        return result
    obj.axial_rate_density = axial
    rates = np.array([0., -2*alpha, 2*alpha, 0., 0., 1., 0., 0.])
    family = TangentFamily(rates, np.zeros(8), np.zeros((5, 8)), np.zeros(5), 0.)
    return obj, family, diffusion


def test_exact_gaussian_diffusion_recovers_velocity_diffusivity_and_heat_flux():
    obj, family, diffusion = manufactured_diffusion()
    y, n = np.array([0., .05, .2, .4]), np.array([0., -.04, .1, -.1])
    result = obj.physical_transport(y, n, family, gamma=0., thermal_species_ratio=1.)
    q = .5*((y/obj.sy)**2+(n/obj.sn)**2)
    data = obj.local(q)
    assert result["velocity_y"] == pytest.approx(np.zeros(4), abs=1e-12)
    assert result["velocity_n"] == pytest.approx(np.zeros(4), abs=1e-12)
    assert result["diffusivity_y"] == pytest.approx(diffusion, rel=1e-12)
    assert result["diffusivity_n"] == pytest.approx(diffusion, rel=1e-12)
    assert result["species_flux_y"] == pytest.approx(diffusion*y/obj.sy**2*data["c"], rel=1e-12)
    assert result["enthalpy_flux_n"] == pytest.approx(diffusion*n/obj.sn**2*data["h"], rel=1e-12)


def test_negative_inferred_diffusion_is_not_clipped():
    obj, family, _ = manufactured_diffusion()
    negative = TangentFamily(-family.origin, family.response, family.jacobian, family.sources, 0.)
    with pytest.raises(ValueError, match="counter-gradient"):
        obj.physical_transport(np.array([.1]), np.array([0.]), negative, gamma=0., thermal_species_ratio=1.)


def test_curved_elliptic_physical_flux_obeys_local_species_balance():
    obj, family, diffusion = manufactured_diffusion()
    obj.sy *= 2.
    obj.sn *= .5
    rates = family.origin.copy()
    rates[3] = .3
    family = TangentFamily(rates, family.response, family.jacobian, family.sources, 0.)
    y, n, step = np.array([.04, .2, -.3]), np.array([.03, -.1, .08]), 1e-6
    def flux(yy, nn, component):
        result = obj.physical_transport(yy, nn, family, gamma=0., thermal_species_ratio=1.)
        return (1.-rates[3]*nn)*result[component]
    divergence = (flux(y+step, n, "species_flux_y")-flux(y-step, n, "species_flux_y"))/(2*step)
    divergence += (flux(y, n+step, "species_flux_n")-flux(y, n-step, "species_flux_n"))/(2*step)
    q = .5*((y/obj.sy)**2+(n/obj.sn)**2)
    d = obj.local(q)
    alpha = .5*rates[2]
    axial_rate_fixed_physical_coordinates = 2.*alpha*d["c"]*d["u"]*(q-1.)
    assert divergence == pytest.approx(-axial_rate_fixed_physical_coordinates, rel=1e-7, abs=1e-10)
    result = obj.physical_transport(y, n, family, gamma=0., thermal_species_ratio=1.)
    assert result["diffusivity_y"]/result["diffusivity_n"] == pytest.approx(16.)


def test_finite_moving_faces_conserve_species_in_known_diffusion():
    obj, family, _ = manufactured_diffusion()
    from scipy.special import roots_legendre
    nodes, weights = roots_legendre(24)
    length = math.sqrt(2*obj.partition.q0)
    ly, ln = length*obj.sy, length*obj.sn
    alpha = .5*family.origin[2]
    outward = 0.
    for axis in (0, 1):
        for sign in (-1., 1.):
            y = np.full_like(nodes, sign*ly) if axis == 0 else ly*nodes
            n = ln*nodes if axis == 0 else np.full_like(nodes, sign*ln)
            result = obj.physical_transport(y, n, family, gamma=0., thermal_species_ratio=1.)
            q = .5*((y/obj.sy)**2+(n/obj.sn)**2)
            data = obj.local(q)
            normal = sign*result["species_flux_y" if axis == 0 else "species_flux_n"]
            speed = alpha*(ly if axis == 0 else ln)
            outward += np.sum(weights*(ln if axis == 0 else ly)*(normal-data["c"]*data["u"]*speed))
    assert outward == pytest.approx(0., abs=1e-12)


@pytest.mark.parametrize("ratio", [0., -1., math.nan, math.inf])
def test_unphysical_thermal_species_ratio_is_rejected(ratio):
    obj, family, _ = manufactured_diffusion()
    with pytest.raises(ValueError, match="ratio"):
        obj.physical_transport(np.array([.1]), np.array([0.]), family, gamma=0., thermal_species_ratio=ratio)


def test_missing_thermal_ratio_is_not_replaced_by_a_default():
    obj, family, _ = manufactured_diffusion()
    with pytest.raises(TypeError):
        obj.physical_transport(np.array([.1]), np.array([0.]), family, gamma=0.)


def test_positivity_interval_and_impossible_constraints():
    assert nonnegative_affine_interval([2., 3.], [1., -1.]) == (-2., 3.)
    lo, hi = nonnegative_affine_interval([-2., -3.], [1., -1.])
    assert lo > hi
    assert nonnegative_affine_interval([-1.], [0.]) == (math.inf, -math.inf)


def test_partition_cannot_silently_follow_changed_width(mixing):
    mixing.section.thermal_width_ratio *= 1.01
    with pytest.raises(RuntimeError, match="cannot be reused"):
        mixing.partition.nodes()
