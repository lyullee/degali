"""Necessary joint normal-stress/TKE budgets and explicit PSD completions.

These are algebraic covariance constraints, NOT a turbulence closure. Flux
tradeoffs require nonnegative weights A*rho*u. They bound flux values only,
not derivatives or net pressure-adjusted forces/work.
"""

import math
import numpy as np


def minimum_tke_from_axial_rms(shear_covariance, axial_rms):
    """Sharp pointwise TKE lower bound from one measured normal stress.

    ``axial_rms`` supplies sqrt(R_ss), while ``shear_covariance`` supplies
    the two R_s-perp components. Positive-semidefiniteness requires

        k >= 0.5 * (R_ss + ||R_s-perp||**2 / R_ss).

    This is a realizability constraint, not an isotropy assumption or a
    dissipation closure. When only axial RMS is available, pass zero shear;
    the resulting ``0.5*axial_rms**2`` is the weakest defensible bound.
    """
    r = np.asarray(shear_covariance, float)
    rms = np.asarray(axial_rms, float)
    if r.ndim != 2 or r.shape[1] != 2 or rms.shape != (len(r),):
        raise ValueError('matching two-component shear and axial-RMS arrays required')
    if (not np.all(np.isfinite(r)) or not np.all(np.isfinite(rms))
            or np.any(rms < 0.)):
        raise ValueError('finite covariance components and nonnegative axial RMS required')
    variance = rms*rms
    shear_squared = np.sum(r*r, axis=1)
    positive = variance > 0.
    if np.any(shear_squared[~positive] > 0.):
        raise ValueError('zero axial RMS cannot support nonzero shear covariance')
    minimum = .5*variance
    minimum[positive] += .5*shear_squared[positive]/variance[positive]
    if not np.all(np.isfinite(minimum)):
        raise ValueError('TKE lower bound overflows the supplied numerical range')
    return dict(minimum_tke=minimum, axial_variance=variance,
                uses_isotropy=False, physical_closure=False)


def minimum_tke_from_two_normal_rms(axial_rms, transverse_rms):
    """Sharp TKE lower bound from two observed normal RMS components.

    If axial and one transverse velocity fluctuations are observed,

        k >= 0.5 * (R_ss + R_perp_perp).

    The unobserved second transverse normal stress and every covariance are
    not completed: their smallest realizable values are zero.  This is only a
    lower bound, not an isotropy assumption, Reynolds-stress model, or
    dissipation closure.  It is useful when a publication supplies plotted
    axial and radial RMS profiles but no covariance tensor.
    """
    axial = np.asarray(axial_rms, float)
    transverse = np.asarray(transverse_rms, float)
    if (axial.ndim != 1 or transverse.shape != axial.shape
            or not np.all(np.isfinite(axial))
            or not np.all(np.isfinite(transverse))
            or np.any(axial < 0.) or np.any(transverse < 0.)):
        raise ValueError('matching finite nonnegative axial and transverse RMS arrays required')
    axial_variance = axial*axial
    transverse_variance = transverse*transverse
    minimum = .5*(axial_variance + transverse_variance)
    if not np.all(np.isfinite(minimum)):
        raise ValueError('TKE lower bound overflows the supplied numerical range')
    return dict(minimum_tke=minimum, axial_variance=axial_variance,
                transverse_variance=transverse_variance,
                unmeasured_normal_variance_lower_bound=np.zeros_like(minimum),
                uses_isotropy=False, physical_closure=False)


def normal_work_interval(shear_bound_flux, tke_flux):
    """Sharp admissible N interval given B and K; stable at large K/B.

    B=int(A*rho*u*|R_st|), K=int(A*rho*u*k), N=int(A*rho*u*R_ss).
    The same algebra also gives the pointwise axial covariance interval.
    """
    b, k = float(shear_bound_flux), float(tke_flux)
    if not all(math.isfinite(v) and v >= 0. for v in (b, k)):
        raise ValueError('finite nonnegative shear and TKE fluxes required')
    if k < b:
        return dict(feasible=False, minimum=None, maximum=None, physical_closure=False)
    if k == 0.:
        return dict(feasible=True, minimum=0., maximum=0., physical_closure=False)
    ratio = b/k
    root = math.sqrt(1.-ratio)*math.sqrt(1.+ratio)
    upper = k*(1.+root)
    lower = b*(b/upper)
    if not math.isfinite(upper):
        raise ValueError('normal-stress interval overflows the supplied numerical range')
    return dict(feasible=True, minimum=lower, maximum=upper, physical_closure=False)


def minimum_tke_for_normal_work_cap(shear_bound_flux, normal_work_cap):
    b, cap = float(shear_bound_flux), float(normal_work_cap)
    if not math.isfinite(b) or b < 0. or not math.isfinite(cap) or cap <= 0.:
        raise ValueError('finite nonnegative shear flux and positive normal-work cap required')
    lower = b if cap >= b else .5*cap+.5*b*(b/cap)
    if not math.isfinite(lower):
        raise ValueError('required TKE flux overflows the supplied numerical range')
    return lower


def covariance_completion(shear_covariance, tke, axial_variance):
    """One PSD covariance at specified k and R_ss; no negative clipping.

    The remaining transverse variance is split equally as a mathematical
    construction. This choice must not be interpreted as measured isotropy.
    """
    r = np.asarray(shear_covariance, float)
    k, a = np.asarray(tke, float), np.asarray(axial_variance, float)
    if r.ndim != 2 or r.shape[1] != 2 or k.shape != (len(r),) or a.shape != k.shape:
        raise ValueError('matching two-component shear, k and axial variance arrays required')
    if (not all(np.all(np.isfinite(v)) for v in (r, k, a))
            or np.any(k < 0.) or np.any(a < 0.)):
        raise ValueError('finite covariance components and nonnegative variances required')
    out = np.zeros((len(r), 3, 3))
    out[:, 0, 0] = a
    out[:, 0, 1:], out[:, 1:, 0] = r, r
    positive = a > 0.
    if np.any(np.linalg.norm(r[~positive], axis=1) > 0.):
        raise ValueError('zero axial variance cannot support nonzero shear covariance')
    residual = 2*k-a
    residual[positive] -= np.sum(r[positive]**2, axis=1)/a[positive]
    if np.any(residual < 0.):
        raise ValueError('specified k and axial variance cannot realize the shear; no clipping')
    out[positive, 1:, 1:] = np.einsum('ni,nj->nij', r[positive], r[positive])/a[positive, None, None]
    out[:, 1, 1] += .5*residual
    out[:, 2, 2] += .5*residual
    if not np.all(np.isfinite(out)):
        raise ValueError('covariance completion overflow')
    return dict(covariance=out, transverse_remainder=residual, physical_closure=False,
                completion_is_mathematical_only=True)


def axial_normal_fluxes(area, weights, density, velocity, axial_variance):
    """Evaluate BOTH missing mean axial momentum and stress-work terms.

    weights integrate normalized physical cross-section coordinates, WITHOUT
    area. The caller must also supply pressure and matching production in a
    coupled transport equation; these two values alone do not close it.
    """
    w, rho, u, a = [np.asarray(v, float) for v in (weights, density, velocity, axial_variance)]
    if w.ndim != 1 or any(v.shape != w.shape for v in (rho, u, a)):
        raise ValueError('matching one-dimensional quadrature fields required')
    if (not math.isfinite(area) or area <= 0. or not all(np.all(np.isfinite(v)) for v in (w, rho, u, a))
            or np.any(w < 0.) or np.any(rho <= 0.) or np.any(a < 0.)):
        raise ValueError('finite nonnegative quadrature/variance and positive area/density required')
    momentum = float(area*(w @ (rho*a)))
    work = float(area*(w @ (rho*u*a)))
    if not math.isfinite(momentum) or not math.isfinite(work):
        raise ValueError('normal-stress flux overflow')
    return dict(axial_momentum_correction=momentum, axial_stress_work=work,
                positive_advective_weights=bool(np.all(u >= 0.)),
                coupled_pressure_production_closed=False, adopted=False)
