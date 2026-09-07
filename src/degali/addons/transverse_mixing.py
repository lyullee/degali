"""Conditional transverse mass/species reconstruction for the LH2 section.

This supplies a conservative radial-in-deforming-coordinates shape ansatz.
It retains gamma=dlog(beta_H)/ds as an unknown, and does not choose thermal
diffusivity or a mechanical energy-exchange distribution. No plume solver.
"""

from dataclasses import dataclass
import math

import numpy as np

from .phase_radial_quadrature import PhaseRadialQuadrature


@dataclass(frozen=True)
class TangentFamily:
    """Rates of log(rho_c),log(Cc),log(A),theta,log(uc),x,z,log(beta)."""
    origin: np.ndarray
    response: np.ndarray
    jacobian: np.ndarray
    sources: np.ndarray
    maximum_scaled_residual: float

    def at(self, gamma):
        if not math.isfinite(gamma):
            raise ValueError("thermal width rate must be finite")
        return self.origin+gamma*self.response


def nonnegative_affine_interval(origin, response):
    """Sampled necessary condition a+b*gamma>=0, with no clipping."""
    a, b = np.broadcast_arrays(np.asarray(origin, float), np.asarray(response, float))
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("positivity constraints must be finite")
    if np.any((b == 0.) & (a < 0.)):
        return math.inf, -math.inf
    lower = float(np.max(-a[b > 0.]/b[b > 0.], initial=-math.inf))
    upper = float(np.min(-a[b < 0.]/b[b < 0.], initial=math.inf))
    return lower, upper


class ConservativeTransverseMixing:
    def __init__(self, section, state, *, probes=513):
        self.partition = PhaseRadialQuadrature(section, state, probes=probes)
        self.section, self.state = section, np.array(state, float)
        self.area = self.state[2]
        self.sy, self.sn = section.section_widths(self.state)
        self.wind = section._wind(self.state)
        self.wind_partials = np.zeros(8)
        self.log_sy_partials = np.zeros(8)
        for col, state_col, logarithmic in ((2, 2, True), (3, 3, False), (5, 5, False), (6, 6, False)):
            plus, minus = self.state.copy(), self.state.copy()
            step = 1e-5 if logarithmic else 1e-5*max(abs(self.state[state_col]), 1.)
            if logarithmic:
                plus[state_col] *= math.exp(step)
                minus[state_col] *= math.exp(-step)
            else:
                plus[state_col] += step
                minus[state_col] -= step
            self.wind_partials[col] = (section._wind(plus)-section._wind(minus))/(2*step)
            self.log_sy_partials[col] = math.log(section.section_widths(plus)[0]/section.section_widths(minus)[0])/(2*step)
        self.log_sn_partials = -self.log_sy_partials.copy()
        self.log_sn_partials[2] += 1.  # Enforce sigma_y*sigma_n=A exactly.

    def local(self, q):
        self.partition.check()
        q = np.asarray(q, float)
        if q.ndim != 1 or np.any(q < 0.) or not np.all(np.isfinite(q)):
            raise ValueError("radial samples must be a finite nonnegative vector")
        model, state = self.section, self.state
        rc, yc, _, theta, uc, *_ = state
        rho, y, _, h = model.thermodynamic_profile(state, np.exp(-q))
        c = rho*y
        hr, hc = model._phase_partials(rho, c)
        if np.any(hr >= 0.):
            raise ValueError("phase density slope must be negative")
        p = 1./model.thermal_width_ratio**2
        hq = -p*h
        rhoq = (hq+hc*c)/hr
        yq = (-c-y*rhoq)/rho
        specific_hq = (hq-h/rho*rhoq)/rho
        center_hr, center_hc = model._phase_partials(np.asarray(rc), np.asarray(rc*yc))
        dh = np.zeros((len(q), 8))
        dh[:, 0] = center_hr*rc*np.exp(-p*q)
        dh[:, 1] = center_hc*rc*yc*np.exp(-p*q)
        dh[:, 7] = 2.*p*q*h
        dc = np.zeros_like(dh)
        dc[:, 1] = c
        drho = (dh-hc[:, None]*dc)/hr[:, None]
        excess = uc*np.exp(-model.velocity_shape_exponent*q)
        u = self.wind*math.cos(theta)+excess
        du = np.broadcast_to(self.wind_partials*math.cos(theta), dh.shape).copy()
        du[:, 3] -= self.wind*math.sin(theta)
        du[:, 4] += excess
        return dict(rho=rho, y=y, c=c, h=h, u=u, yq=yq, hq=specific_hq,
                    drho=drho, dc=dc, dh=dh, du=du)

    def axial_rate_density(self, q):
        """B_M and B_C derivative matrices, including moving area A."""
        d = self.local(q)
        u, rho, c = d["u"], d["rho"], d["c"]
        mass = self.area*(d["drho"]*u[:, None]+rho[:, None]*d["du"])
        species = self.area*(d["dc"]*u[:, None]+c[:, None]*d["du"])
        mass[:, 2] += self.area*rho*u
        species[:, 2] += self.area*c*u
        return np.stack([mass, species], axis=1)

    def flux_rate_density(self, q):
        d = self.local(q)
        u, rho, h, c = d["u"], d["rho"], d["h"], d["c"]
        dm = self.area*(d["drho"]*u[:, None]+rho[:, None]*d["du"])
        dc = self.area*(d["dc"]*u[:, None]+c[:, None]*d["du"])
        dm[:, 2] += self.area*rho*u
        dc[:, 2] += self.area*c*u
        momentum = self.area*rho*u*u
        dp = self.area*(d["drho"]*(u*u)[:, None]+(2.*rho*u)[:, None]*d["du"])
        dp[:, 2] += momentum
        theta = self.state[3]
        dx, dz = dp*math.cos(theta), dp*math.sin(theta)
        dx[:, 3] -= momentum*math.sin(theta)
        dz[:, 3] += momentum*math.cos(theta)
        de = self.area*(d["dh"]*u[:, None]+h[:, None]*d["du"]
            +.5*d["drho"]*(u**3)[:, None]+(1.5*rho*u*u)[:, None]*d["du"])
        de[:, 2] += self.area*(h*u+.5*rho*u**3)
        return np.stack([dm, dc, dx, dz, de], axis=1)

    def flux_jacobian(self, order=16):
        return self.partition.integrate(self.flux_rate_density, order, square=True)

    def tangent_family(self, sources, *, order=16):
        sources = np.asarray(sources, float)
        if sources.shape != (5,) or not np.all(np.isfinite(sources)):
            raise ValueError("supply all five finite conserved-flux source terms")
        jac = self.flux_jacobian(order)
        origin, response = np.zeros(8), np.zeros(8)
        origin[5:7] = math.cos(self.state[3]), math.sin(self.state[3])
        response[7] = 1.
        scale = np.maximum(abs(sources), [1., 1., 1., 1., 1.])
        matrix = jac[:, :5]/scale[:, None]
        origin[:5] = np.linalg.solve(matrix, (sources-jac@origin)/scale)
        response[:5] = np.linalg.solve(matrix, -(jac@response)/scale)
        residual = max(np.max(abs(jac@origin-sources)/scale), np.max(abs(jac@response)/scale))
        return TangentFamily(origin, response, jac, sources, float(residual))

    def transverse_basis(self, q, *, order=16):
        """f_M/f_C matrices. No physical gamma is selected."""
        q = np.asarray(q, float)
        integral = self.partition.cumulative(self.axial_rate_density, q, order)
        result = np.empty_like(integral)
        nonzero = q > 0.
        result[nonzero] = -integral[nonzero]/(2.*q[nonzero, None, None])
        if np.any(~nonzero):
            result[~nonzero] = -.5*self.axial_rate_density(np.array([0.]))[0]
        return result

    def mixing_family(self, q, family, *, order=16):
        q = np.asarray(q, float)
        d = self.local(q)
        f = self.transverse_basis(q, order=order)@np.column_stack([family.origin, family.response])
        residual = f[:, 1]-d["y"][:, None]*f[:, 0]
        denominator = self.area*d["rho"]*d["yq"]
        if np.any(abs(denominator) < 1e-14):
            raise ValueError("species gradient does not identify a finite mixing diffusivity")
        chi = -residual/denominator[:, None]
        return dict(mass=f[:, 0], species=f[:, 1], chi=chi, y=d["y"], yq=d["yq"])

    def physical_transport(self, y, n, family, *, gamma, thermal_species_ratio, order=16):
        """Conditional physical velocity and heat/species fluxes, not an RHS.

        A ratio and gamma are REQUIRED. Negative diffusion and singular
        curved coordinates are rejected, not clipped. No heat source is
        supplied, so a total-energy/thermal moment budget remains unclosed.
        """
        if not math.isfinite(thermal_species_ratio) or thermal_species_ratio <= 0.:
            raise ValueError("supply a positive finite thermal/species diffusivity ratio")
        y, n = np.broadcast_arrays(np.asarray(y, float), np.asarray(n, float))
        shape = y.shape
        y, n = y.ravel(), n.ravel()
        q = .5*((y/self.sy)**2+(n/self.sn)**2)
        rates = family.at(gamma)
        metric = 1.-rates[3]*n
        if np.any(metric <= 0.):
            raise ValueError("curved section metric is not positive")
        data = self.local(q)
        mix = self.mixing_family(q, family, order=order)
        fm, fc, chi = [mix[key]@np.array([1., gamma]) for key in ("mass", "species", "chi")]
        if np.any(chi < 0.):
            raise ValueError("implied mixing is counter-gradient; no positive diffusivity closure")
        h, rho, u = data["h"]/data["rho"], data["rho"], data["u"]
        fh = h*fm-self.area*rho*thermal_species_ratio*chi*data["hq"]
        ay, an = self.log_sy_partials@rates, self.log_sn_partials@rates
        result = {
            "velocity_y": y/metric*(u*ay+fm/(self.area*rho)),
            "velocity_n": n/metric*(u*an+fm/(self.area*rho)),
            "species_flux_y": y/metric*(data["c"]*u*ay+fc/self.area),
            "species_flux_n": n/metric*(data["c"]*u*an+fc/self.area),
            "enthalpy_flux_y": y/metric*(data["h"]*u*ay+fh/self.area),
            "enthalpy_flux_n": n/metric*(data["h"]*u*an+fh/self.area),
            "diffusivity_y": self.sy**2*chi/metric,
            "diffusivity_n": self.sn**2*chi/metric,
            "axial_enthalpy_flux": data["h"]*u,
        }
        return {key: value.reshape(shape) for key, value in result.items()}
