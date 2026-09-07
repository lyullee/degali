"""Boundary-only independent enthalpy width constrained by six moments.

Buoyancy force per unit length is a boundary moment, NOT a conserved flux.
No downstream thermal-width transport equation is supplied or assumed.
"""

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import least_squares

from .enthalpy_profile import GaussianEnthalpyCrosswind


@dataclass(frozen=True)
class BuoyancyWidthProjection:
    state: np.ndarray
    thermal_width_ratio: float
    moments: np.ndarray
    relative_residuals: np.ndarray
    quadrature_residuals: np.ndarray
    quadrature_points: int
    at_width_bound: bool
    success: bool
    message: str
    solver_message: str
    solver_evaluations: int
    solver_optimality: float


class BuoyancyConstrainedEnthalpySection(GaussianEnthalpyCrosswind):
    """C=Cc*G and H=Hc*G**(1/beta_H**2), at a single section only."""

    def __init__(self, *args, thermal_width_ratio=1., **kwargs):
        super().__init__(*args, **kwargs)
        self.thermal_width_ratio = thermal_width_ratio

    @property
    def thermal_width_ratio(self):
        return self._thermal_width_ratio

    @thermal_width_ratio.setter
    def thermal_width_ratio(self, value):
        value = float(value)
        if not math.isfinite(value) or value <= 0.:
            raise ValueError("thermal width ratio must be finite and positive")
        self._thermal_width_ratio = value

    def thermodynamic_profile(self, state, shape):
        rho_c, yc, *_ = self._physical(state)
        shape = np.asarray(shape, float)
        if not np.all(np.isfinite(shape)) or np.any(shape < 0.) or np.any(shape > 1.):
            raise ValueError("boundary Gaussian shape must lie in [0, 1]")
        cc = rho_c*yc
        hc, _ = self.phase_inverse.enthalpy_and_slope(rho_c, cc)
        c = cc*shape
        h = hc*shape**(1./self.thermal_width_ratio**2)
        guess = self.rhoa + (rho_c-self.rhoa)*shape
        rho, y, t = self.phase_inverse.state(c, h, density_guess=guess)
        return rho, y, t, h

    def moments(self, state, *, quadrature_points=None):
        """Mass, H2, Px, Pz, total energy, buoyancy force per metre."""
        if self.energy_transport != "total":
            raise ValueError("six-moment boundary requires total-energy transport")
        _, _, area, theta, *_ = self._physical(state)
        points = self.quadrature_points if quadrature_points is None else quadrature_points
        weight = area*self._quadrature(points)[1]
        u, rho, y, _, h = self.profiles(state, quadrature_points=points)
        momentum = np.sum(weight*rho*u*u)
        return np.array([
            np.sum(weight*rho*u), np.sum(weight*rho*y*u),
            momentum*math.cos(theta), momentum*math.sin(theta),
            np.sum(weight*(h*u + .5*rho*u**3)),
            9.81*np.sum(weight*(self.rhoa-rho)),
        ])

    @staticmethod
    def moment_scales(values):
        return np.maximum(np.abs(values), [1e-12, 1e-12, 1., 1., 1., 1e-3])

    def _phase_partials(self, rho, c):
        """dH/drho at fixed C, dH/dC at fixed rho on the same bilinear table."""
        inv = self.phase_inverse
        _, h_rho = inv.enthalpy_and_slope(rho, c)
        y, ti = c/rho, inv.a/(rho+inv.k*c)
        i = np.clip(np.searchsorted(inv.t_grid, ti, side="right")-1, 0, len(inv.t_grid)-2)
        j = np.clip(np.searchsorted(inv.y_grid, y, side="right")-1, 0, len(inv.y_grid)-2)
        dt, dy = inv.t_grid[i+1]-inv.t_grid[i], inv.y_grid[j+1]-inv.y_grid[j]
        f, g = (ti-inv.t_grid[i])/dt, (y-inv.y_grid[j])/dy
        h00, h01 = inv.h_values[i, j], inv.h_values[i, j+1]
        h10, h11 = inv.h_values[i+1, j], inv.h_values[i+1, j+1]
        h_t = ((1-g)*(h10-h00)+g*(h11-h01))/dt
        h_y = ((1-f)*(h01-h00)+f*(h11-h10))/dy
        return h_rho, -h_t*ti*inv.k/(rho+inv.k*c)+h_y/rho

    def reduced_moment_jacobian(self, state, *, quadrature_points=None):
        """Four-moment derivatives in log(rho_c,A,uc,beta), at fixed H2 flux.

        Analytic implicit phase derivatives avoid differencing a tight phase
        root through table knots. Only the pre-existing wind geometry uses a
        scalar finite difference. Theta and the receptor coordinates are fixed.
        """
        rc, yc, area, theta, uc, *_ = self._physical(state)
        points = self.quadrature_points if quadrature_points is None else quadrature_points
        q, weights = self._quadrature(points)
        weight = area*weights
        u, rho, y, _, h = self.profiles(state, quadrature_points=points)
        cc, c = rc*yc, rho*y
        p = 1./self.thermal_width_ratio**2
        v = self._wind(state)*math.cos(theta)
        plus, minus = np.array(state), np.array(state)
        step = 1e-5
        plus[2] *= math.exp(step)
        minus[2] *= math.exp(-step)
        dv = (self._wind(plus)-self._wind(minus))*math.cos(theta)/(2*step)
        i1, iu = self.k.profile_integral(1.), self.k.profile_integral(1.+self.velocity_shape_exponent)
        denominator = v*i1+uc*iu
        dcc = np.array([0., -cc*(1.+dv*i1/denominator), -cc*uc*iu/denominator, 0.])
        hc_rho, hc_c = self._phase_partials(np.asarray(rc), np.asarray(cc))
        dhc = hc_c*dcc
        dhc[0] = hc_rho*rc
        dc = np.exp(-q)[:, None]*dcc
        dh = np.exp(-p*q)[:, None]*dhc
        dh[:, 3] += 2.*p*q*h
        h_rho, h_c = self._phase_partials(rho, c)
        drho = (dh-h_c[:, None]*dc)/h_rho[:, None]
        du = np.zeros_like(drho)
        du[:, 1], du[:, 2] = dv, uc*np.exp(-self.velocity_shape_exponent*q)
        integrands = np.array([rho*u, rho*u*u, h*u+.5*rho*u**3, 9.81*(self.rhoa-rho)])
        derivatives = np.array([
            drho*u[:, None]+rho[:, None]*du,
            drho*(u*u)[:, None]+(2*rho*u)[:, None]*du,
            dh*u[:, None]+h[:, None]*du+.5*drho*(u**3)[:, None]+(1.5*rho*u*u)[:, None]*du,
            -9.81*drho,
        ])
        jacobian = np.sum(derivatives*weight[None, :, None], axis=1)
        jacobian[:, 1] += np.sum(integrands*weight, axis=1)
        return jacobian

    def project_buoyancy(self, target, initial_state, *, maximum_quadrature_points=2048):
        """Solve four unknowns after eliminating H2 fraction and direction."""
        target = np.asarray(target, float)
        initial = np.asarray(initial_state, float)
        self._physical(initial)
        if (target.shape != (6,) or not np.all(np.isfinite(target))
                or min(target[0], target[1], target[2]) <= 0.):
            raise ValueError("boundary target must have six finite moments and positive mass/H2/Px")
        points = self.quadrature_points
        if maximum_quadrature_points < 2*points:
            raise ValueError("quadrature screen requires at least two orders")
        theta = math.atan2(target[3], target[2])
        axial = math.hypot(target[2], target[3])
        reduced_target = target[[0, 2, 4, 5]].copy()
        reduced_target[1] = axial
        reduced_scales = self.moment_scales(target)[[0, 2, 4, 5]]
        reduced_scales[1] = max(axial, 1.)
        bounds = (np.log([1e-3, 1e-14, 1e-8, .5]), np.log([100., 100., 1e4, 2.]))
        i1, iu = self.k.profile_integral(1.), self.k.profile_integral(1.+self.velocity_shape_exponent)

        def decode(parameters):
            rho, area, uc, beta = np.exp(parameters)
            self.thermal_width_ratio = beta
            state = np.array([rho, .1, area, theta, uc, initial[5], initial[6]])
            wind = self._wind(state)*math.cos(theta)
            state[1] = target[1]/(area*rho*(wind*i1+uc*iu))
            return state

        def residual(parameters):
            try:
                moment = self.moments(decode(parameters), quadrature_points=points)
                reduced = moment[[0, 2, 4, 5]].copy()
                reduced[1] = math.hypot(moment[2], moment[3])
                return (reduced-reduced_target)/reduced_scales
            except (ValueError, RuntimeError, FloatingPointError):
                return np.full(4, 1e6)

        def fit(start):
            def jacobian(parameters):
                try:
                    return self.reduced_moment_jacobian(decode(parameters), quadrature_points=points)/reduced_scales[:, None]
                except (ValueError, RuntimeError, FloatingPointError):
                    return np.zeros((4, 4))
            result = least_squares(residual, np.clip(start, bounds[0]+1e-10, bounds[1]-1e-10),
                jac=jacobian, bounds=bounds, x_scale="jac",
                xtol=1e-11, ftol=1e-11, gtol=1e-11, max_nfev=1000)
            return result, float(np.max(np.abs(residual(result.x))))

        def feasible_start(beta):
            # Close the five fluxes at a fixed width before asking the fourth
            # degree of freedom to match force. This is a numerical warm start,
            # not a second physical closure or downstream transport equation.
            self.thermal_width_ratio = beta
            self.quadrature_points = points
            try:
                start = GaussianEnthalpyCrosswind._match_flux_array(self, target[:5], initial)
            except (ValueError, RuntimeError):
                start = initial
            return np.log([start[0], start[2], start[4], beta])

        parameters = feasible_start(1.)
        best = fit(parameters)
        if not best[0].success or best[1] > 1e-8:
            for beta in (.8, 1.2):
                start = feasible_start(beta)
                attempt = fit(start)
                if attempt[1] < best[1]:
                    best = attempt
        parameters = best[0].x
        quadrature_errors = np.full(6, math.inf)
        message = "six-moment fit failed"
        while best[0].success and best[1] <= 1e-8:
            state = decode(parameters)
            coarse = self.moments(state, quadrature_points=points)
            fine_points = min(2*points, maximum_quadrature_points)
            fine = self.moments(state, quadrature_points=fine_points)
            quadrature_errors = np.abs(fine-coarse)/self.moment_scales(fine)
            old_points = points
            points = fine_points
            best = fit(parameters)
            parameters = best[0].x
            state = decode(parameters)
            # Recheck at the REFITTED state; no stale pre-refit convergence claim.
            coarse = self.moments(state, quadrature_points=old_points)
            fine = self.moments(state, quadrature_points=points)
            quadrature_errors = np.abs(fine-coarse)/self.moment_scales(fine)
            if np.max(quadrature_errors) <= 1e-5:
                message = "six moments and quadrature converged"
                break
            if points >= maximum_quadrature_points:
                message = "six-moment quadrature did not converge"
                break
        state = decode(parameters)
        actual = self.moments(state, quadrature_points=points)
        errors = np.abs(actual-target)/self.moment_scales(target)
        beta = self.thermal_width_ratio
        at_bound = min(abs(beta-.5), abs(beta-2.)) < 1e-5
        success = bool(best[0].success and max(errors) <= 1e-8
                       and max(quadrature_errors) <= 1e-5 and not at_bound)
        if not best[0].success or max(errors) > 1e-8:
            message = "six-moment fit failed"
        self.quadrature_points = points
        return BuoyancyWidthProjection(state, beta, actual, errors, quadrature_errors,
                                       points, at_bound, success, message, str(best[0].message),
                                       int(best[0].nfev), float(best[0].optimality))

    def _no_transport(self, *args, **kwargs):
        raise NotImplementedError("boundary-only profile: independent thermal-width transport is not closed")

    solve = _no_transport
    derivatives = _no_transport
    source_terms = _no_transport
    _match_flux_array = _no_transport
    _receptor = _no_transport

    def project(self, *args, **kwargs):
        raise NotImplementedError("use project_buoyancy with all six boundary moments")
