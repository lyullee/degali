"""Manufactured and analytic tests; no fitted thermal transport closure."""

import math
import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from degali.addons.thermal_moments import (
    EnthalpyMomentOperators, PlanarMovingSection, specific_enthalpy_slope, eddy_enthalpy_flux,
)


STATE = np.array([1.1312, .1537, .02, .01, 25., 1., 1.5])


def analytic_q_integral(p, delta):
    z = math.sqrt(p*math.pi*delta**2/8.)
    return (2.*math.pi*math.erf(z)**2 - 4.*math.sqrt(math.pi)*math.erf(z)*z*math.exp(-z*z))/p**2


@pytest.mark.parametrize("beta", [.7, 1., 1.08, 1.7])
def test_advective_moment_matches_finite_square_analytic_integral(section, beta):
    section.thermal_width_ratio = beta
    # Anisotropy must survive the radial quadrature reduction.
    section.jetplume._split = lambda a, *_: (2.*math.sqrt(a), .5*math.sqrt(a))
    op = EnthalpyMomentOperators(section)
    hc = section.phase_inverse.enthalpy_and_slope(STATE[0], STATE[0]*STATE[1])[0]
    p = 1./beta**2
    sy, sn = section.section_widths(STATE)
    expected = STATE[2]*(sy*sy+sn*sn)*hc*(
        section._wind(STATE)*math.cos(STATE[3])*analytic_q_integral(p, section.k.delta)
        + STATE[4]*analytic_q_integral(p+section.velocity_shape_exponent, section.k.delta))
    assert op.enthalpy_second_moment(STATE) == pytest.approx(expected, rel=1e-10)
    def axial(y, n):
        q = .5*((y/sy)**2+(n/sn)**2)
        return hc*np.exp(-p*q)*(2.*math.cos(STATE[3])+STATE[4]*np.exp(-section.velocity_shape_exponent*q))
    assert op.geometry(STATE).second_moment(axial) == pytest.approx(expected, rel=1e-10)


def test_reduced_thermal_moment_derivative_at_fixed_species_flux(section):
    section.thermal_width_ratio = 1.04
    op = EnthalpyMomentOperators(section)
    expected = op.reduced_moment_jacobian(STATE)
    target_species = section.moments(STATE)[1]
    params = np.log([STATE[0], STATE[2], STATE[4], 1.04])
    def evaluate(p):
        state = STATE.copy()
        rho, area, uc, beta = np.exp(p)
        state[[0, 2, 4]] = rho, area, uc
        section.thermal_width_ratio = beta
        den = area*rho*(section._wind(state)*np.cos(state[3])*section.k.profile_integral(1.)
                       + uc*section.k.profile_integral(1.+section.velocity_shape_exponent))
        state[1] = target_species/den
        return op.enthalpy_second_moment(state)
    step = 1e-7
    actual = [(evaluate(params+step*e)-evaluate(params-step*e))/(2*step) for e in np.eye(4)]
    assert actual == pytest.approx(expected, rel=5e-5, abs=1e-6)


def gaussian_heat_fields(s):
    """Exact U*d_s h = D*laplacian(h), with negative ambient-relative h."""
    rho, velocity, diffusion = 1.3, 4., .7
    vy, vn = .16+2*diffusion*s/velocity, .25+2*diffusion*s/velocity
    def fields(y, n):
        h = -120.*math.sqrt(.16*.25/(vy*vn))*np.exp(-.5*(y*y/vy+n*n/vn))
        return rho*velocity*h, rho*diffusion*y/vy*h, rho*diffusion*n/vn*h, 0.
    return fields


def heat_geometry(s, rates):
    return PlanarMovingSection(.41*math.exp(rates[0]*s), .53*math.exp(rates[1]*s),
                               2.15, 0., *rates)


@pytest.mark.parametrize("rates", [(0., 0.), (.12, .07), (-.12, -.07), (.1, -.2)])
def test_known_diffusion_balance_keeps_finite_edges_and_deformation(rates):
    s, step = .3, 1e-5
    geom = heat_geometry(s, rates)
    budget = geom.budget(gaussian_heat_fields(s))
    def moment(at):
        return heat_geometry(at, rates).second_moment(lambda y, n: gaussian_heat_fields(at)(y, n)[0])
    expected = (moment(s+step)-moment(s-step))/(2*step)
    assert budget.derivative == pytest.approx(expected, rel=1e-7, abs=1e-7)
    assert abs(budget.outward_boundary_flux) > 1.
    if rates != (0., 0.):
        assert abs(budget.moving_boundary) > .1
    # Dropping the edge does NOT reproduce the exact solution.
    assert abs(budget.transverse_volume+budget.moving_boundary-expected) > 1.


def test_known_diffusion_short_march_has_fourth_order_refinement():
    # Only the known heat equation is marched here, not the unclosed plume.
    rates = (.12, -.07)
    def exact(s):
        return heat_geometry(s, rates).second_moment(lambda y, n: gaussian_heat_fields(s)(y, n)[0])
    def rhs(s):
        return heat_geometry(s, rates).budget(gaussian_heat_fields(s)).derivative
    errors = []
    for count in (2, 4, 8):
        ds, total = .4/count, exact(.2)
        for i in range(count):
            s = .2+i*ds
            total += ds/6.*(rhs(s)+4.*rhs(s+.5*ds)+rhs(s+ds))
        errors.append(abs(total-exact(.6)))
    assert errors[0]/errors[1] > 12.
    assert errors[1]/errors[2] > 12.
    assert errors[2]/abs(exact(.6)) < 1e-7


@pytest.mark.parametrize("curvature", [-.7, 0., .7])
def test_curved_moving_section_manufactured_conservation(curvature):
    rates, a, b, c = (.13, -.09), .4, 1.3, -.8
    geom = PlanarMovingSection(.4, .6, 2.15, curvature, *rates)
    def fields(y, n):
        metric = 1.-curvature*n
        js = 1.+y*y+.5*n*n
        jy, jn = b*y*(1.+.2*n), c*n
        source = (a*js+metric*b*(1.+.2*n)+c*(1.-2*curvature*n))/metric
        return js, jy, jn, source
    actual = geom.budget(fields)
    def moment(s):
        domain = PlanarMovingSection(.4*math.exp(rates[0]*s), .6*math.exp(rates[1]*s),
                                     2.15, curvature, *rates)
        return domain.second_moment(lambda y, n: math.exp(a*s)*(1.+y*y+.5*n*n))
    step = 1e-5
    expected = (moment(step)-moment(-step))/(2*step)
    assert actual.derivative == pytest.approx(expected, rel=1e-7, abs=1e-7)


def test_density_gradient_alone_is_not_a_specific_enthalpy_gradient():
    rho, rho_q, h = np.array([1., 2., 3.]), np.array([.2, -.4, .7]), -500.
    correct = specific_enthalpy_slope(rho, rho*h, rho_q, rho_q*h)
    assert correct == pytest.approx(np.zeros(3), abs=1e-12)
    assert np.all(np.abs(rho_q*h) > 1.)  # The incorrect volumetric gradient is nonzero.
    assert eddy_enthalpy_flux(rho, correct, .7) == pytest.approx(np.zeros(3), abs=1e-12)


def test_phase_enthalpy_gradient_matches_independent_difference(section):
    section.thermal_width_ratio = 1.04
    op = EnthalpyMomentOperators(section)
    q, step = np.array([.013, .237, .713, 1.173, 2.319, 3.417]), 1e-7
    def specific(q):
        rho, _, _, h = section.thermodynamic_profile(STATE, np.exp(-q))
        return h/rho
    rho = section.thermodynamic_profile(STATE, np.exp(-q))[0]
    numerical = -rho*(specific(q+step)-specific(q-step))/(2*step)
    assert op.radial_diffusion_coefficient(STATE, q) == pytest.approx(numerical, rel=1e-4)


def test_unit_response_matches_tensor_and_four_face_integrals(section):
    # A smooth manufactured H/rho field isolates geometry from table knots.
    op = EnthalpyMomentOperators(section)
    op.radial_diffusion_coefficient = lambda state, q: -3.*np.exp(-.8*q)
    section.jetplume._split = lambda a, *_: (2.*math.sqrt(a), .5*math.sqrt(a))
    sy, sn = section.section_widths(STATE)
    def fields(y, n):
        q = .5*((y/sy)**2+(n/sn)**2)
        k = op.radial_diffusion_coefficient(STATE, q)
        return 0., k*y/sy**2, k*n/sn**2, 0.
    tensor = op.geometry(STATE).budget(fields)
    radial = op.unit_diffusion_response(STATE)
    assert radial.transverse_volume == pytest.approx(tensor.transverse_volume, rel=1e-10)
    assert radial.outward_boundary_flux == pytest.approx(tensor.outward_boundary_flux, rel=1e-10)
    assert radial.derivative == pytest.approx(tensor.derivative, rel=1e-10)


@pytest.mark.parametrize("diffusion", [-.1, math.nan, math.inf])
def test_invalid_diffusivity_rejected(diffusion):
    with pytest.raises(ValueError, match="diffusivity"):
        eddy_enthalpy_flux(1., .2, diffusion)


@pytest.mark.parametrize("rho", [0., -1., math.nan])
def test_invalid_density_rejected(rho):
    with pytest.raises(ValueError, match="density"):
        specific_enthalpy_slope(rho, 1., 1., 1.)


def test_invalid_geometry_and_unspecified_fields_are_rejected():
    with pytest.raises(ValueError, match="positive"):
        PlanarMovingSection(1., 1., 2.15, 1., 0., 0.)
    geom = PlanarMovingSection(1., 1., 2.15, 0., 0., 0.)
    with pytest.raises(ValueError, match="source"):
        geom.budget(lambda y, n: (1., 0., 0.))
    with pytest.raises(ValueError, match="integer"):
        geom.nodes(8.5)


def test_moment_operators_do_not_enable_unclosed_transport_or_ground(section):
    op = EnthalpyMomentOperators(section)
    assert not hasattr(op, "solve")
    with pytest.raises(NotImplementedError, match="boundary-only"):
        section.solve(STATE)
    section.ground_interaction = "geometry"
    with pytest.raises(NotImplementedError, match="ground"):
        EnthalpyMomentOperators(section)
