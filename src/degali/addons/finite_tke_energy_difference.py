"""Actual smooth Q-flux differences on independent rectangular quadrature.

Uses actual +/- scalar values with cancellation-safe product differences,
not the Q Jacobian or any EOS derivative. Optional mpmath stays in the
existing scalar-pair verification helper; no solver dependency is added.
"""

import math
import numpy as np
from .stable_enthalpy_difference import scalar_pair
from .phase_radial_quadrature import gauss_rule


def _product_pair(left, right):
    lm, ld = left
    rm, rd = right
    return lm*rm+.25*ld*rd, lm*rd+rm*ld


def actual_tke_moment_difference(model, rates, step, *, order=64, precision=70):
    rates = np.asarray(rates, float)
    if (rates.shape != (model.count,) or not np.all(np.isfinite(rates))
            or not np.isfinite(step) or step <= 0.):
        raise ValueError('matching finite direction and positive difference distance required')
    b = model.base
    pair = scalar_pair(b.projection, b.parameters, rates[:model.mean_count], step,
                       precision=precision, wind_policy='exact_constraint')
    q0 = math.exp(model.parameters[0])
    center_shift = step*rates[model.mean_count]
    q_pair = q0*math.cosh(center_shift), 2*q0*math.sinh(center_shift)
    aq = _product_pair((pair['area']['mean'], pair['area']['difference']), q_pair)
    aqw = _product_pair(aq, (pair['parallel']['mean'], pair['parallel']['difference']))
    aqu = _product_pair(aq, (pair['uc']['mean'], pair['uc']['difference']))
    x, w = gauss_rule(order)
    a, c = np.meshgrid(x, x, indexing='ij')
    prep = b.projection.prepare(a.ravel(), c.ravel())
    psi, q = prep['psi'], prep['q']
    correction = psi @ model.parameters[1:]
    shift = step*(psi @ rates[model.mean_count+1:])
    if np.any(abs(correction)+abs(shift) > .1):
        raise ValueError('actual perturbed Q shape leaves the trust region')
    shape = np.exp(-b.projection.section.velocity_shape_exponent*q+correction)
    velocity_shape = np.exp(-b.projection.section.velocity_shape_exponent*q)
    cs, sn = np.cosh(shift), np.sinh(shift)
    delta = shape*(aqw[1]*cs+2*aqw[0]*sn+velocity_shape*(aqu[1]*cs+2*aqu[0]*sn))
    theta = np.column_stack([np.ones(len(q)), psi])
    weight = (8*b.projection.q0*np.outer(w, w)).ravel()
    derivatives = np.array([math.fsum(weight*delta*col)/(2*step) for col in theta.T])
    return dict(derivatives=derivatives, derivative=float(derivatives[0]), step_m=step, order=order,
                actual_values_not_EOS_derivatives=True, all_fields_arbitrary_precision=False,
                scalar_pair=pair)
