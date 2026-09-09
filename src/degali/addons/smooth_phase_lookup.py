"""C1 phase-property lookup retaining the established equilibrium-table nodes."""

from __future__ import annotations

import weakref

import numpy as np
from scipy.interpolate import PchipInterpolator

from .enthalpy_profile import PhaseMassEnthalpyInverter


class HermitePhaseSurface:
    """Vectorized C1 bicubic Hermite surface with precomputed nodal slopes."""

    def __init__(self, x, y, values):
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
        values = np.asarray(values, float)
        if (self.x.ndim != 1 or self.y.ndim != 1
                or values.shape != (len(self.x), len(self.y))
                or len(self.x) < 4 or len(self.y) < 4
                or not np.all(np.diff(self.x) > 0.)
                or not np.all(np.diff(self.y) > 0.)
                or not np.all(np.isfinite(values))):
            raise ValueError("Hermite phase surface requires a finite rectangular grid")
        fx = PchipInterpolator(self.x, values, axis=0).derivative()(self.x)
        fy = PchipInterpolator(self.y, values, axis=1).derivative()(self.y)
        fxy_x = PchipInterpolator(self.y, fx, axis=1).derivative()(self.y)
        fxy_y = PchipInterpolator(self.x, fy, axis=0).derivative()(self.x)
        fxy = .5*(fxy_x+fxy_y)
        dx, dy = np.diff(self.x)[:, None], np.diff(self.y)[None, :]
        p = np.empty((len(self.x)-1, len(self.y)-1, 4, 4))
        p[:, :, 0, 0], p[:, :, 0, 1] = values[:-1, :-1], values[:-1, 1:]
        p[:, :, 1, 0], p[:, :, 1, 1] = values[1:, :-1], values[1:, 1:]
        p[:, :, 0, 2], p[:, :, 0, 3] = fy[:-1, :-1]*dy, fy[:-1, 1:]*dy
        p[:, :, 1, 2], p[:, :, 1, 3] = fy[1:, :-1]*dy, fy[1:, 1:]*dy
        p[:, :, 2, 0], p[:, :, 2, 1] = fx[:-1, :-1]*dx, fx[:-1, 1:]*dx
        p[:, :, 3, 0], p[:, :, 3, 1] = fx[1:, :-1]*dx, fx[1:, 1:]*dx
        scale = dx*dy
        p[:, :, 2, 2], p[:, :, 2, 3] = fxy[:-1, :-1]*scale, fxy[:-1, 1:]*scale
        p[:, :, 3, 2], p[:, :, 3, 3] = fxy[1:, :-1]*scale, fxy[1:, 1:]*scale
        matrix = np.array([
            [1., 0., 0., 0.], [0., 0., 1., 0.],
            [-3., 3., -2., -1.], [2., -2., 1., 1.],
        ])
        self.coefficients = np.einsum("ab,ijbc,cd->ijad", matrix, p, matrix.T)

    def evaluate(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        if (not np.all(np.isfinite(x)) or not np.all(np.isfinite(y))
                or np.any(x < self.x[0]) or np.any(x > self.x[-1])
                or np.any(y < self.y[0]) or np.any(y > self.y[-1])):
            raise ValueError("phase surface coordinates lie outside the lookup domain")
        shape = x.shape
        x, y = x.ravel(), y.ravel()
        i = np.clip(np.searchsorted(self.x, x, side="right")-1, 0, len(self.x)-2)
        j = np.clip(np.searchsorted(self.y, y, side="right")-1, 0, len(self.y)-2)
        dx, dy = self.x[i+1]-self.x[i], self.y[j+1]-self.y[j]
        tx, ty = (x-self.x[i])/dx, (y-self.y[j])/dy
        bx = np.column_stack([np.ones_like(tx), tx, tx*tx, tx*tx*tx])
        by = np.column_stack([np.ones_like(ty), ty, ty*ty, ty*ty*ty])
        dbx = np.column_stack([
            np.zeros_like(tx), np.ones_like(tx), 2.*tx, 3.*tx*tx,
        ])/dx[:, None]
        dby = np.column_stack([
            np.zeros_like(ty), np.ones_like(ty), 2.*ty, 3.*ty*ty,
        ])/dy[:, None]
        coefficients = self.coefficients[i, j]
        value = np.einsum("ni,nij,nj->n", bx, coefficients, by)
        derivative_x = np.einsum("ni,nij,nj->n", dbx, coefficients, by)
        derivative_y = np.einsum("ni,nij,nj->n", bx, coefficients, dby)
        return (
            value.reshape(shape), derivative_x.reshape(shape),
            derivative_y.reshape(shape),
        )


_SURFACE_CACHE = weakref.WeakKeyDictionary()


class SmoothPhaseMassEnthalpyInverter(PhaseMassEnthalpyInverter):
    """Invert the same phase table through a nodally exact C1 surface."""

    def __init__(self, thermodynamics):
        super().__init__(thermodynamics)
        if thermodynamics not in _SURFACE_CACHE:
            lookup = thermodynamics._condensed_lookup
            _SURFACE_CACHE[thermodynamics] = (
                HermitePhaseSurface(self.t_grid, self.y_grid, lookup[2].values),
                HermitePhaseSurface(self.t_grid, self.y_grid, lookup[3].values),
            )
        self.temperature_surface, self.enthalpy_surface = _SURFACE_CACHE[thermodynamics]

    def _coordinates(self, density, fuel_density):
        rho, c = np.broadcast_arrays(
            np.asarray(density, float), np.asarray(fuel_density, float),
        )
        if (not np.all(np.isfinite(rho)) or not np.all(np.isfinite(c))
                or np.any(rho <= 0.) or np.any(c < 0.)):
            raise ValueError("local mass/enthalpy inputs must be finite and physical")
        y = c/rho
        ti = self.a/(rho+self.k*c)
        eps = 32.*np.finfo(float).eps
        if (np.any(y > 1.+eps) or np.any(ti < self.t_grid[0]*(1.-eps))
                or np.any(ti > self.t_grid[-1]*(1.+eps))):
            raise ValueError("local C/rho lies outside the smooth phase domain")
        return rho, c, np.clip(ti, self.t_grid[0], self.t_grid[-1]), np.minimum(y, 1.)

    def enthalpy_partials(self, density, fuel_density):
        rho, c, ti, y = self._coordinates(density, fuel_density)
        _, h_t, h_y = self.enthalpy_surface.evaluate(ti, y)
        h_rho = -h_t*ti/(rho+self.k*c)-h_y*y/rho
        h_c = -h_t*ti*self.k/(rho+self.k*c)+h_y/rho
        return h_rho, h_c

    def enthalpy_and_slope(self, density, fuel_density):
        rho, c, ti, y = self._coordinates(density, fuel_density)
        h, h_t, h_y = self.enthalpy_surface.evaluate(ti, y)
        slope = -h_t*ti/(rho+self.k*c)-h_y*y/rho
        return h, slope

    def state(self, fuel_density, enthalpy_density, *, density_guess=None):
        c, h = np.broadcast_arrays(
            np.asarray(fuel_density, float), np.asarray(enthalpy_density, float),
        )
        shape = c.shape
        c, h = c.ravel(), h.ravel()
        if not np.all(np.isfinite(c)) or not np.all(np.isfinite(h)) or np.any(c < 0.):
            raise ValueError("local C and H must be finite, with C nonnegative")
        lo = np.maximum(c, self.a/self.t_grid[-1]-self.k*c)
        hi = self.a/self.t_grid[0]-self.k*c
        if np.any(hi < lo) or np.any(lo <= 0.):
            raise ValueError("H2 mass density exceeds the smooth phase domain")
        h_lo, _ = self.enthalpy_and_slope(lo, c)
        h_hi, _ = self.enthalpy_and_slope(hi, c)
        tolerance = 1e-10*np.maximum(abs(h), 1.)
        if np.any(h > h_lo+tolerance) or np.any(h < h_hi-tolerance):
            raise ValueError("local C,H has no bracketed smooth phase state")
        guess = self.th.ambient_density if density_guess is None else density_guess
        rho = np.clip(np.broadcast_to(np.asarray(guess, float), shape).ravel(), lo, hi)
        if not np.all(np.isfinite(rho)):
            raise ValueError("density guess must be finite")
        for _ in range(48):
            value, slope = self.enthalpy_and_slope(rho, c)
            residual = value-h
            done = abs(residual) <= tolerance
            if np.all(done):
                break
            lo = np.where((residual > 0.) & ~done, rho, lo)
            hi = np.where((residual < 0.) & ~done, rho, hi)
            with np.errstate(divide="ignore", invalid="ignore"):
                newton = rho-residual/slope
            safe = (slope < 0.) & np.isfinite(newton) & (newton > lo) & (newton < hi)
            rho = np.where(done, rho, np.where(safe, newton, .5*(lo+hi)))
        else:
            raise RuntimeError("smooth phase inversion did not converge")
        y = c/rho
        ti = self.a/(rho+self.k*c)
        temperature = self.temperature_surface.evaluate(ti, y)[0]
        represented = self.enthalpy_surface.evaluate(ti, y)[0]
        if np.any(abs(represented-h) > 1.01*tolerance):
            raise RuntimeError("smooth phase interpolation and enthalpy root disagree")
        return rho.reshape(shape), y.reshape(shape), temperature.reshape(shape)


__all__ = ["HermitePhaseSurface", "SmoothPhaseMassEnthalpyInverter"]
