"""Opt-in Gaussian volumetric-enthalpy cross-section screening.

The candidate shares one Gaussian for H2 mass density and ambient-relative
enthalpy density, not for bulk density. It is a structural hypothesis, not
an experimentally validated turbulent transport law. No fitted coefficient
or new phase/EOS table is introduced here.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import least_squares

from .axisymmetric_jet import ConservedGaussianJet, R_UNIVERSAL
from .energy_crosswind import IndependentEnergyCrosswind


class PhaseMassEnthalpyInverter:
    """Invert local H2 mass density C and enthalpy density H at fixed pressure.

    The existing bilinear (ideal-temperature, mass-fraction) phase table is
    the property oracle. A safeguarded Newton solve uses its analytic slope
    along Y=C/rho. Coordinates outside the table are rejected, never clamped
    into an apparently physical solution.
    """

    def __init__(self, thermodynamics: ConservedGaussianJet):
        if not (thermodynamics.consistent_phase_ambient
                and thermodynamics.temperature_dependent_phase_enthalpy
                and thermodynamics.equilibrium_air_condensation):
            raise ValueError("enthalpy profile requires consistent component phase thermodynamics")
        if thermodynamics._ambient_enthalpy != 0.0:
            raise ValueError("enthalpy profile requires zero ambient-reference enthalpy")
        self.th = thermodynamics
        self.t_grid, self.y_grid, _, h_table = thermodynamics._condensed_lookup
        self.h_values = h_table.values
        ma, mh = thermodynamics._humid_ambient_molecular_weight, thermodynamics.fuel_molecular_weight
        self.a = thermodynamics.ambient_pressure * ma / R_UNIVERSAL
        self.k = ma / mh - 1.0

    def enthalpy_and_slope(self, density, fuel_density):
        """H(rho,C) and dH/drho at fixed C using the exact bilinear table."""
        rho, c = np.broadcast_arrays(np.asarray(density, float), np.asarray(fuel_density, float))
        if not (np.all(np.isfinite(rho)) and np.all(np.isfinite(c))):
            raise ValueError("local mass/enthalpy inputs must be finite")
        if np.any(rho <= 0.) or np.any(c < 0.):
            raise ValueError("local density must be positive and H2 mass density nonnegative")
        y = c / rho
        ti = self.a / (rho + self.k * c)
        eps = 32 * np.finfo(float).eps
        if (np.any(y > 1. + eps) or np.any(ti < self.t_grid[0] * (1.-eps))
                or np.any(ti > self.t_grid[-1] * (1.+eps))):
            raise ValueError("local C/rho lies outside the phase interpolation domain")
        # Only arithmetic roundoff at a valid endpoint is rounded inward.
        ti = np.clip(ti, self.t_grid[0], self.t_grid[-1])
        y = np.minimum(y, 1.)
        i = np.clip(np.searchsorted(self.t_grid, ti, side="right")-1, 0, len(self.t_grid)-2)
        j = np.clip(np.searchsorted(self.y_grid, y, side="right")-1, 0, len(self.y_grid)-2)
        dt, dy = self.t_grid[i+1]-self.t_grid[i], self.y_grid[j+1]-self.y_grid[j]
        f, g = (ti-self.t_grid[i])/dt, (y-self.y_grid[j])/dy
        h00, h01 = self.h_values[i, j], self.h_values[i, j+1]
        h10, h11 = self.h_values[i+1, j], self.h_values[i+1, j+1]
        low, high = h00 + g*(h01-h00), h10 + g*(h11-h10)
        h = low + f*(high-low)
        dh_dt = (high-low)/dt
        dh_dy = ((1.-f)*(h01-h00) + f*(h11-h10))/dy
        slope = -dh_dt * ti/(rho+self.k*c) - dh_dy*y/rho
        return h, slope

    def state(self, fuel_density, enthalpy_density, *, density_guess=None):
        """Return rho, Y, T; reject unbracketed or inaccurate roots."""
        c, h = np.broadcast_arrays(np.asarray(fuel_density, float), np.asarray(enthalpy_density, float))
        shape = c.shape
        c, h = c.ravel(), h.ravel()
        if not (np.all(np.isfinite(c)) and np.all(np.isfinite(h))) or np.any(c < 0.):
            raise ValueError("local C and H must be finite, with C nonnegative")
        lo = np.maximum(c, self.a/self.t_grid[-1] - self.k*c)
        hi = self.a/self.t_grid[0] - self.k*c
        if np.any(hi < lo) or np.any(lo <= 0.):
            raise ValueError("H2 mass density exceeds the phase interpolation domain")
        h_lo, _ = self.enthalpy_and_slope(lo, c)
        h_hi, _ = self.enthalpy_and_slope(hi, c)
        tolerance = 1e-10 * np.maximum(np.abs(h), 1.)
        if np.any(h > h_lo+tolerance) or np.any(h < h_hi-tolerance):
            raise ValueError("local C,H has no bracketed phase state")
        guess = self.th.ambient_density if density_guess is None else density_guess
        rho = np.clip(np.broadcast_to(np.asarray(guess, float), shape).ravel(), lo, hi)
        if not np.all(np.isfinite(rho)):
            raise ValueError("density guess must be finite")
        for _ in range(48):
            value, slope = self.enthalpy_and_slope(rho, c)
            residual = value-h
            done = np.abs(residual) <= tolerance
            if np.all(done):
                break
            lo = np.where((residual > 0.) & ~done, rho, lo)
            hi = np.where((residual < 0.) & ~done, rho, hi)
            with np.errstate(divide="ignore", invalid="ignore"):
                newton = rho-residual/slope
            safe = (slope < 0.) & np.isfinite(newton) & (newton > lo) & (newton < hi)
            rho = np.where(done, rho, np.where(safe, newton, .5*(lo+hi)))
        else:
            raise RuntimeError("local mass/enthalpy phase inversion did not converge")
        y = c/rho
        t, represented_h = self.th._condensed_air_state(rho, y)
        if np.any(np.abs(represented_h-h) > 1.01*tolerance):
            raise RuntimeError("phase interpolation and enthalpy root disagree")
        return rho.reshape(shape), y.reshape(shape), t.reshape(shape)


class GaussianEnthalpyCrosswind(IndependentEnergyCrosswind):
    """Five-state section with Gaussian C and H, and non-Gaussian density.

    The downstream inverse eliminates H2 fraction using the Gaussian species
    integral, then solves density, area and velocity against mass, momentum
    and energy. It never uses the density-Gaussian analytic mass shortcut.
    """

    def __init__(self, jetplume, thermodynamics, **kwargs):
        super().__init__(jetplume, thermodynamics, **kwargs)
        self.phase_inverse = PhaseMassEnthalpyInverter(thermodynamics)

    def _quadrature(self, points):
        """Exact polar reduction of JETPLU's finite square Gaussian domain.

        Every flux integrand depends only on q=pi/8*(a*a+b*b). The square's
        angular measure is 2*pi below q0=pi*delta**2/8, and
        2*pi-8*atan(t) above it, with q=q0*(1+t*t), 0<t<1.
        This preserves the domain and its corners; it is not a disk cutoff.
        """
        if points not in self._quadrature_cache:
            nodes, weights = np.polynomial.legendre.leggauss(points)
            s, w = .5*(nodes+1.), .5*weights
            q0 = math.pi*self.k.delta**2/8.
            self._quadrature_cache[points] = (
                np.r_[q0*s, q0*(1.+s*s)],
                np.r_[2.*math.pi*q0*w, (2.*math.pi-8.*np.arctan(s))*2.*q0*s*w],
            )
        return self._quadrature_cache[points]

    def quadrature_error(self, coarse, fine):
        """Non-Gaussian density requires checking every conserved flux."""
        c, f = self._as_array(coarse), self._as_array(fine)
        return float(np.max(np.abs(c-f)/np.maximum(np.abs(f), [1e-12, 1e-12, 1., 1., 1.])))

    def thermodynamic_profile(self, state, shape):
        rho_c, yc, *_ = self._physical(state)
        shape = np.asarray(shape, float)
        if not np.all(np.isfinite(shape)) or np.any(shape < 0.):
            raise ValueError("Gaussian shape must be finite and nonnegative")
        cc = rho_c*yc
        hc, _ = self.phase_inverse.enthalpy_and_slope(rho_c, cc)
        c, h = cc*shape, hc*shape
        guess = self.rhoa + (rho_c-self.rhoa)*shape
        rho, y, t = self.phase_inverse.state(c, h, density_guess=guess)
        return rho, y, t, h

    def profiles(self, state, *, quadrature_points=None):
        _, _, _, theta, uc, *_ = self._physical(state)
        points = self.quadrature_points if quadrature_points is None else int(quadrature_points)
        exponent, _ = self._quadrature(points)
        rho, y, t, h = self.thermodynamic_profile(state, np.exp(-exponent))
        velocity = self._wind(state)*math.cos(theta) + uc*np.exp(-self.velocity_shape_exponent*exponent)
        if np.any(velocity <= 0.):
            raise ValueError("enthalpy Gaussian profile has reverse axial flow")
        return velocity, rho, y, t, h

    def centre_temperature(self, state):
        # Strict domain validation, including centre trial states in projection.
        return float(self.thermodynamic_profile(state, np.array([1.]))[2][0])

    def buoyancy_force(self, state):
        _, weights = self._quadrature(self.quadrature_points)
        density = self.profiles(state)[1]
        return float(9.81 * state[2] * np.sum((self.rhoa-density)*weights))

    def _receptor(self, state, lateral, height):
        *_, centre_z = self._physical(state)
        sy, sz = self.section_widths(state)
        shape = math.exp(-.5*(float(lateral)/sy)**2) * (
            math.exp(-.5*((float(height)-centre_z)/sz)**2)
            + math.exp(-.5*((float(height)+centre_z)/sz)**2)
        )
        return self.thermodynamic_profile(state, np.array([shape]))

    def point_temperature(self, state, lateral, height):
        return float(self._receptor(state, lateral, height)[2][0])

    def point_mole_fraction(self, state, lateral, height):
        y = float(self._receptor(state, lateral, height)[1][0])
        numerator = y/self.thermodynamics.fuel_molecular_weight
        return numerator/(numerator + (1.-y)/self.thermodynamics._humid_ambient_molecular_weight)

    def _match_flux_array(self, target_values, initial_state, *, relative_tolerance=1e-9):
        """Three-unknown flux inverse valid for the non-Gaussian density."""
        initial_state = np.asarray(initial_state, float)
        rho0, _, area0, _, uc0, x, z = self._physical(initial_state)
        target = np.asarray(target_values, float)
        if target.shape != (5,) or not np.all(np.isfinite(target)):
            raise ValueError("enthalpy-profile targets must be five finite fluxes")
        mass, fuel, px, pz, energy = target
        momentum = math.hypot(px, pz)
        theta = math.atan2(pz, px)
        if min(mass, fuel, momentum) <= 0.:
            raise ValueError("target mass, hydrogen and momentum must be positive")
        i1 = self.k.profile_integral(1.)
        i1u = self.k.profile_integral(1.+self.velocity_shape_exponent)
        scales = np.maximum(np.abs(target), [1e-12, 1e-12, 1., 1., 1.])
        residual_scales = np.array([scales[0], max(momentum, 1.), scales[4]])
        lower, upper = np.log([1e-3, 1e-14, 1e-8]), np.log([100., 100., 1e4])

        def decode(parameters):
            rho, area, uc = np.exp(parameters)
            trial = np.array([rho, .1, area, theta, uc, x, z])
            ambient_axial = self._wind(trial)*math.cos(theta)
            trial[1] = fuel/(area*rho*(ambient_axial*i1 + uc*i1u))
            return trial

        def residual(parameters):
            try:
                actual = self._as_array(self.integral_fluxes(decode(parameters)))
                return (np.array([actual[0], math.hypot(actual[2], actual[3]), actual[4]])
                        - [mass, momentum, energy])/residual_scales
            except (RuntimeError, ValueError, FloatingPointError):
                return np.full(3, 1e6)

        def fit_from(start):
            fit = least_squares(
                residual, np.clip(start, lower+1e-12, upper-1e-12), bounds=(lower, upper),
                xtol=1e-11, ftol=1e-11, gtol=1e-11, diff_step=2e-6, max_nfev=200,
            )
            try:
                state = decode(fit.x)
                actual = self._as_array(self.integral_fluxes(state))
                score = float(np.max(np.abs(actual-target)/scales))
            except (RuntimeError, ValueError, FloatingPointError):
                state, score = np.full(7, np.nan), math.inf
            return fit, state, score

        initial = np.log([rho0, area0, uc0])
        best = fit_from(initial)
        if not best[0].success or best[2] > relative_tolerance:
            offset = math.log(2.)
            for da in (-offset, 0., offset):
                for du in (-offset, 0., offset):
                    if da == du == 0.:
                        continue
                    trial = fit_from(initial+[0., da, du])
                    if trial[2] < best[2]:
                        best = trial
                    if best[0].success and best[2] <= relative_tolerance:
                        return best[1]
        if not best[0].success or best[2] > relative_tolerance:
            raise RuntimeError(f"enthalpy-profile flux inverse failed; maximum residual {best[2]:.3e}")
        return best[1]
