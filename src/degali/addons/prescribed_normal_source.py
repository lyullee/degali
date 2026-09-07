"""Boundary-only Q plus axial Reynolds-normal-stress source reallocation.

Explicit variance ratio and ambient-pressure assumption; no real-trial input
is supplied or inferred. Both normal momentum and stress work are accounted.
"""

import math
import numpy as np
from .phase_radial_quadrature import gauss_rule
from .prescribed_tke_source import prescribed_q_flux, retract_source_with_prescribed_tke


def prescribed_normal_fluxes(projection, tke_parameters, *, axial_variance_fraction, order=64):
    ratio = float(axial_variance_fraction)
    if not math.isfinite(ratio) or not 0. <= ratio <= 2.:
        raise ValueError('explicit axial variance / k in [0,2] required')
    q_flux = prescribed_q_flux(projection, tke_parameters, order=order)
    par = np.asarray(tke_parameters, float)
    x, w = gauss_rule(order)
    a, b = np.meshgrid(x, x, indexing='ij')
    prep = projection.prepare(a.ravel(), b.ravel())
    q = np.exp(par[0]-projection.section.velocity_shape_exponent*prep['q']+prep['psi'] @ par[1:])
    weights = (8*projection.q0*projection.mixing.area*np.outer(w, w)).ravel()
    momentum = float(ratio*(weights @ q))
    work = ratio*q_flux
    if not math.isfinite(momentum) or not math.isfinite(work):
        raise ValueError('normal source flux overflows')
    return dict(tke_flux=q_flux, axial_normal_momentum=momentum, axial_normal_work=work,
                axial_variance_fraction=ratio, shear_covariance_not_checked=True)


def retract_source_with_normal_stress(projection, parameters, tke_parameters, *, total_moment_target,
                                     axial_variance_fraction, pressure_assumption,
                                     callback=None, maximum_iterations=8):
    if pressure_assumption != 'ambient_pressure_no_compensation':
        raise ValueError('explicit ambient-pressure/no-compensation assumption required; no pressure field is solved')
    target = np.asarray(total_moment_target, float)
    if target.shape != (6,) or not np.all(np.isfinite(target)) or np.any(target[:3] <= 0.):
        raise ValueError('six finite enriched-projection target moments required')
    scales = np.maximum(abs(target), [1e-12, 1e-12, 1., 1., 1e-3, 1.])
    coarse, fine = [prescribed_normal_fluxes(projection, tke_parameters,
        axial_variance_fraction=axial_variance_fraction, order=order) for order in (32, 64)]
    refinement = {key: abs(fine[key]-coarse[key])/max(abs(fine[key]), 1.)
                  for key in ('tke_flux', 'axial_normal_momentum', 'axial_normal_work')}
    if max(refinement.values()) > 1e-10:
        raise ValueError('prescribed normal-stress source integrals did not refine')
    reduced = target.copy()
    reduced[2] -= fine['axial_normal_momentum']
    reduced[3] -= fine['axial_normal_work']
    if reduced[2] <= 0.:
        raise ValueError('prescribed normal stress exhausts the forward momentum budget')
    inner = retract_source_with_prescribed_tke(projection, parameters, tke_parameters,
        total_moment_target=reduced, callback=callback, maximum_iterations=maximum_iterations)
    actual = inner['actual_total_moments'].copy()
    actual[2] += fine['axial_normal_momentum']
    actual[3] += fine['axial_normal_work']
    error = abs(actual-target)/scales
    return dict(parameters=inner['parameters'], tke_parameters=inner['tke_parameters'],
                total_moment_target=target.copy(), actual_total_moments=actual,
                original_moment_errors=error, source_fluxes=fine, smooth_refinement=refinement,
                q_only_inner_retraction=inner, pressure_assumption=pressure_assumption,
                numerical_passed=bool(inner['numerical_passed'] and max(error) <= 1e-8),
                physical_initialization_passed=False, normal_transport_closed=False,
                adopted=False, field_scored=False, variance_ratio_was_fitted=False)
