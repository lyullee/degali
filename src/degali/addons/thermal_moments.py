"""Thermal-shape moment operators, NOT a closed downstream plume model.

The mean advective enthalpy moment, a supplied-flux moving-section identity,
and the response per unit turbulent enthalpy diffusivity are kept separate.
No diffusivity, transverse mean velocity, mechanical heating distribution,
or zero derivative of beta_H is silently selected.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from .buoyancy_profile import BuoyancyConstrainedEnthalpySection


@lru_cache(maxsize=24)
def _legendre(points):
    if isinstance(points, bool) or int(points) != points or points < 8:
        raise ValueError("moment quadrature requires an integer of at least eight")
    nodes, weights = np.polynomial.legendre.leggauss(int(points))
    nodes.flags.writeable = weights.flags.writeable = False
    return nodes, weights


def specific_enthalpy_slope(rho, h_density, rho_slope, h_density_slope):
    """d(H/rho)/dq; q is the dimensionless Gaussian radius, not distance."""
    rho, h_density, rho_slope, h_density_slope = np.broadcast_arrays(
        *map(lambda x: np.asarray(x, float), (rho, h_density, rho_slope, h_density_slope)))
    if (not all(np.all(np.isfinite(x)) for x in (rho, h_density, rho_slope, h_density_slope))
            or np.any(rho <= 0.)):
        raise ValueError("enthalpy slope requires finite values and positive density")
    return (h_density_slope - h_density/rho*rho_slope)/rho


def eddy_enthalpy_flux(density, specific_enthalpy_gradient, diffusivity):
    """One flux component for an EXPLICIT supplied diffusivity [m²/s]."""
    rho, gradient, diffusion = np.broadcast_arrays(
        *map(lambda x: np.asarray(x, float), (density, specific_enthalpy_gradient, diffusivity)))
    if (not all(np.all(np.isfinite(x)) for x in (rho, gradient, diffusion))
            or np.any(rho <= 0.) or np.any(diffusion < 0.)):
        raise ValueError("finite positive density and nonnegative diffusivity are required")
    return -rho*diffusion*gradient


@dataclass(frozen=True)
class SecondMomentBudget:
    """All signed contributions to dM2/ds, in W m for enthalpy fluxes."""

    transverse_volume: float
    source_volume: float
    outward_boundary_flux: float
    moving_boundary: float

    @property
    def derivative(self):
        return (self.transverse_volume + self.source_volume
                - self.outward_boundary_flux + self.moving_boundary)


@dataclass(frozen=True)
class PlanarMovingSection:
    """Finite rectangle in a planar curved centreline's local normal plane.

    Flux components are physical components in the orthonormal (s,y,n)
    directions. sigma_n is a NORMAL width, not a projected ground height.
    This geometry has no ground cut, torsion or rotating transverse axes.
    """

    sigma_y: float
    sigma_n: float
    delta: float
    curvature: float
    log_sigma_y_rate: float
    log_sigma_n_rate: float

    def __post_init__(self):
        values = (self.sigma_y, self.sigma_n, self.delta, self.curvature,
                  self.log_sigma_y_rate, self.log_sigma_n_rate)
        if not all(math.isfinite(x) for x in values) or min(values[:3]) <= 0.:
            raise ValueError("section geometry must be finite with positive widths and cutoff")
        if abs(self.curvature)*self.half_widths[1] >= 1.:
            raise ValueError("curved section metric must remain positive everywhere")

    @property
    def half_widths(self):
        factor = math.sqrt(math.pi)/2.*self.delta
        return factor*self.sigma_y, factor*self.sigma_n

    def nodes(self, points):
        nodes, weights = _legendre(points)
        ly, ln = self.half_widths
        yy, nn = np.meshgrid(ly*nodes, ln*nodes, indexing="ij")
        return yy, nn, ly*ln*np.outer(weights, weights)

    def second_moment(self, axial_flux, *, points=48):
        y, n, weights = self.nodes(points)
        js = np.broadcast_to(np.asarray(axial_flux(y, n), float), y.shape)
        if not np.all(np.isfinite(js)):
            raise ValueError("axial flux must be finite")
        return float(np.sum(weights*(y*y+n*n)*js))

    def budget(self, flux_and_source, *, points=48):
        """Evaluate the exact weak identity for EXPLICITLY supplied fields.

        The callable must return Js, Jy, Jn, S at every (y,n), including
        the four faces. It must satisfy the governing balance if the result
        is to equal a physical M2 derivative; this operator cannot close it.
        """
        def sample(y, n):
            fields = flux_and_source(y, n)
            if len(fields) != 4:
                raise ValueError("supply axial/transverse fluxes AND the local source")
            arrays = [np.broadcast_to(np.asarray(f, float), np.shape(y)) for f in fields]
            if not all(np.all(np.isfinite(f)) for f in arrays):
                raise ValueError("fluxes and source must be finite")
            return arrays

        y, n, weights = self.nodes(points)
        _, jy, jn, source = sample(y, n)
        metric = 1.-self.curvature*n
        volume = float(np.sum(weights*2.*(y*jy+n*jn)*metric))
        source_term = float(np.sum(weights*(y*y+n*n)*metric*source))
        nodes, weights = _legendre(points)
        ly, ln = self.half_widths
        outward = motion = 0.
        for axis in (0, 1):
            for sign in (-1., 1.):
                y = np.full_like(nodes, sign*ly) if axis == 0 else ly*nodes
                n = ln*nodes if axis == 0 else np.full_like(nodes, sign*ln)
                js, jy, jn, _ = sample(y, n)
                line_weight = weights*(ln if axis == 0 else ly)
                normal_flux = sign*(jy if axis == 0 else jn)
                boundary_speed = (self.log_sigma_y_rate*ly if axis == 0
                                  else self.log_sigma_n_rate*ln)
                radius2 = y*y+n*n
                outward += float(np.sum(line_weight*radius2*(1.-self.curvature*n)*normal_flux))
                motion += float(np.sum(line_weight*radius2*js*boundary_speed))
        return SecondMomentBudget(volume, source_term, outward, motion)


class EnthalpyMomentOperators:
    """Diagnostics on a six-moment boundary; no downstream solve method."""

    def __init__(self, section):
        if not isinstance(section, BuoyancyConstrainedEnthalpySection):
            raise TypeError("thermal moments require an independent-enthalpy-width section")
        if section.ground_interaction != "free":
            raise NotImplementedError("thermal moment geometry has no ground-cut boundary closure")
        self.section = section

    def geometry(self, state):
        sy, sn = self.section.section_widths(state)
        return PlanarMovingSection(sy, sn, self.section.k.delta, 0., 0., 0.)

    def enthalpy_second_moment(self, state, *, quadrature_points=None):
        """Mean advective integral of (y²+n²)*H*u [W m²].

        Square symmetry makes the angular mean of y²+n² equal to
        (sigma_y²+sigma_n²)*q, also in the four corner sectors.
        """
        model = self.section
        _, _, area, theta, uc, *_ = model._physical(state)
        points = model.quadrature_points if quadrature_points is None else quadrature_points
        q, weights = model._quadrature(points)
        hc, _ = model.phase_inverse.enthalpy_and_slope(state[0], state[0]*state[1])
        h = hc*np.exp(-q/model.thermal_width_ratio**2)
        u = model._wind(state)*math.cos(theta)+uc*np.exp(-model.velocity_shape_exponent*q)
        sy, sn = model.section_widths(state)
        return float(area*(sy*sy+sn*sn)*np.sum(weights*q*h*u))

    def reduced_moment_jacobian(self, state, *, quadrature_points=None):
        """dM2/dlog(rho_c,A,uc,beta_H) at fixed H2 flux, theta,x,z.

        Only the established wind/width geometry is differenced, as in the
        boundary solver. Phase derivatives are analytic within each cell.
        """
        model = self.section
        rc, yc, area, theta, uc, *_ = model._physical(state)
        points = model.quadrature_points if quadrature_points is None else quadrature_points
        q, weights = model._quadrature(points)
        sy, sn = model.section_widths(state)
        spread = sy*sy+sn*sn
        v = model._wind(state)*math.cos(theta)
        plus, minus = np.array(state), np.array(state)
        step = 1e-5
        plus[2] *= math.exp(step)
        minus[2] *= math.exp(-step)
        dv = (model._wind(plus)-model._wind(minus))*math.cos(theta)/(2*step)
        dspread = (sum(w*w for w in model.section_widths(plus))
                   - sum(w*w for w in model.section_widths(minus)))/(2*step)
        i1 = model.k.profile_integral(1.)
        iu = model.k.profile_integral(1.+model.velocity_shape_exponent)
        cc = rc*yc
        den = v*i1+uc*iu
        dcc = np.array([0., -cc*(1.+dv*i1/den), -cc*uc*iu/den, 0.])
        hc, _ = model.phase_inverse.enthalpy_and_slope(rc, cc)
        hr, hc_partial = model._phase_partials(np.asarray(rc), np.asarray(cc))
        dhc = hc_partial*dcc
        dhc[0] = hr*rc
        p = 1./model.thermal_width_ratio**2
        h = hc*np.exp(-p*q)
        dh = np.exp(-p*q)[:, None]*dhc
        dh[:, 3] += 2.*p*q*h
        excess = uc*np.exp(-model.velocity_shape_exponent*q)
        u = v+excess
        du = np.zeros_like(dh)
        du[:, 1], du[:, 2] = dv, excess
        weight = area*spread*q*weights
        derivative = np.sum(weight[:, None]*(u[:, None]*dh+h[:, None]*du), axis=0)
        derivative[1] += np.sum(weight*h*u)*(1.+dspread/spread)
        return derivative

    def radial_diffusion_coefficient(self, state, q):
        """Return K=-rho*d(h)/dq, so Q_perp/D_h = K*grad(q).

        This is an eddy-enthalpy-gradient hypothesis, not a molecular Fourier
        law. No value of D_h, Pr_t, or Sc_t is selected by this method.
        """
        model = self.section
        q = np.asarray(q, float)
        if not np.all(np.isfinite(q)) or np.any(q < 0.):
            raise ValueError("Gaussian radius must be finite and nonnegative")
        rho, y, _, h = model.thermodynamic_profile(state, np.exp(-q))
        c = rho*y
        h_q = -h/model.thermal_width_ratio**2
        hr, hc = model._phase_partials(rho, c)
        if np.any(hr >= 0.) or not np.all(np.isfinite(hr)):
            raise ValueError("phase enthalpy inverse must have a finite negative density slope")
        rho_q = (h_q+hc*c)/hr
        return -rho*specific_enthalpy_slope(rho, h, rho_q, h_q)

    def unit_diffusion_response(self, state, *, quadrature_points=None):
        """Only the diffusion contribution to dM2/ds, DIVIDED BY D_h.

        Fixed free rectangle; no advection, source or moving-edge estimate
        is supplied. These are separated components, never a plume RHS.
        """
        model = self.section
        _, _, area, *_ = model._physical(state)
        points = model.quadrature_points if quadrature_points is None else quadrature_points
        q, weights = model._quadrature(points)
        volume = area*np.sum(weights*4.*q*self.radial_diffusion_coefficient(state, q))
        sy, sn = model.section_widths(state)
        nodes, weights = _legendre(points)
        t, weights = .5*(nodes+1.), .5*weights
        length = math.sqrt(math.pi)/2.*model.k.delta
        edge_q = .5*length**2*(1.+t*t)
        coefficient = self.radial_diffusion_coefficient(state, edge_q)
        outward = 4.*length**4*area*np.sum(
            weights*(2.+((sn/sy)**2+(sy/sn)**2)*t*t)*coefficient)
        return SecondMomentBudget(float(volume), 0., float(outward), 0.)
