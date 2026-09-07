"""Boundary-only source energy matching for an EXPLICIT positive Q profile.

The six targets use the enriched projection convention (one axial momentum
and a thermal second moment), not the distinct two-component source array.
No initial turbulence intensity, time scale or mixing law is selected here.
"""

import copy
import numpy as np
from .phase_radial_quadrature import gauss_rule
from .edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit


def prescribed_q_flux(projection, tke_parameters, *, order=64):
    parameters = np.asarray(tke_parameters, float)
    if parameters.shape != (projection.basis.size+1,) or not np.all(np.isfinite(parameters)):
        raise ValueError('explicit finite log-Q center and square modes required')
    x, w = gauss_rule(order)
    a, b = np.meshgrid(x, x, indexing='ij')
    prep = projection.prepare(a.ravel(), b.ravel())
    correction = prep['psi'] @ parameters[1:]
    if np.any(abs(correction) > .1):
        raise ValueError('prescribed Q shape leaves the log-shape trust region')
    with np.errstate(over='ignore', under='ignore'):
        q = np.exp(parameters[0]-projection.section.velocity_shape_exponent*prep['q']+correction)
    if np.any(q <= 0.) or not np.all(np.isfinite(q)) or np.any(prep['u'] <= 0.):
        raise ValueError('positive finite Q and forward source velocity required')
    weights = (8*projection.q0*projection.mixing.area*np.outer(w, w)).ravel()
    flux = float(weights @ (q*prep['u']))
    if not np.isfinite(flux) or flux <= 0.:
        raise ValueError('positive finite prescribed turbulent energy flux required')
    return flux


def retract_source_with_prescribed_tke(projection, parameters, tke_parameters, *,
                                      total_moment_target, callback=None, maximum_iterations=8):
    """Reallocate an unchanged total source energy to H+meanKE+Q.

    Fixed centers, geometry and velocity; only C/H shapes may change within
    the existing bound. Q is not optimized. Legacy edge evaluation inside
    the reused retraction checks domains only, NOT the new TKE boundary.
    """
    target = np.asarray(total_moment_target, float)
    if target.shape != (6,) or not np.all(np.isfinite(target)) or np.any(target[:3] <= 0.):
        raise ValueError('six finite enriched-projection moments and positive mass/H2/axial momentum required')
    scales = np.maximum(abs(target), [1e-12, 1e-12, 1., 1., 1e-3, 1.])
    coarse_q, fine_q = [prescribed_q_flux(projection, tke_parameters, order=order) for order in (32, 64)]
    q_refinement = abs(fine_q-coarse_q)/max(abs(fine_q), 1.)
    if q_refinement > 1e-10:
        raise ValueError('prescribed Q source flux did not refine')
    view = copy.copy(projection)
    view.target, view.scales = target.copy(), scales.copy()
    view.target[3] -= fine_q
    view._grids, view._edges = {}, {}
    retracted, history = ConservativeEdgeRefit(view).retract(parameters, callback=callback,
                                                            maximum_iterations=maximum_iterations)
    quadrature = FaceSplitSquareMoments(view)
    coarse = quadrature.moments(retracted, order=8, angular_order=8)
    fine = quadrature.moments(retracted, order=8, angular_order=16)
    total = fine.copy()
    total[3] += fine_q
    error = abs(total-target)/scales
    refinement = abs(fine-coarse)/scales
    return dict(parameters=retracted, tke_parameters=np.array(tke_parameters, float),
                total_moment_target=target.copy(), mean_thermal_target=view.target.copy(),
                actual_total_moments=total, actual_mean_thermal_moments=fine, tke_flux=fine_q,
                moment_errors=error, moment_refinement=refinement, tke_flux_refinement=q_refinement,
                history=history, numerical_passed=bool(max(error) <= 1e-8 and max(refinement) <= 1e-5),
                physical_initialization_passed=False, adopted=False, field_scored=False,
                q_was_fitted=False, centers_geometry_velocity_fixed=True)
