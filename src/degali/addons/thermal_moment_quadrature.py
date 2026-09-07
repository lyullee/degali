"""Error-controlled quadrature of phase-table thermal-gradient components.

Separate from thermal_moments so the original fixed-order audit is retained.
An accurate advective moment does not imply accurate phase derivatives.
Requires the pinned scientific research environment (SciPy cubature).
"""

from dataclasses import dataclass
import math

import numpy as np

from .thermal_moments import EnthalpyMomentOperators, SecondMomentBudget


@dataclass(frozen=True)
class AdaptiveDiffusionResponse:
    response: SecondMomentBudget
    estimated_absolute_errors: tuple[float, float]
    estimated_scaled_errors: tuple[float, float]
    component_status: tuple[str, str]
    evaluations: int
    tolerance: float
    converged: bool


def adaptive_unit_diffusion_response(operators: EnthalpyMomentOperators, state, *,
                                     tolerance=1e-7, max_subdivisions=10000):
    """Adapt volume and edge independently, with no chosen physical D_h.

    The error request applies to EACH signed component, scaled by
    max(abs(component),1); their near-cancelling net is not used as scale.
    Returned error estimates are numerical estimates, not rigorous bounds.
    """
    from scipy.integrate import cubature

    if not math.isfinite(tolerance) or not 0. < tolerance < 1.:
        raise ValueError("quadrature tolerance must be finite and lie in (0,1)")
    if isinstance(max_subdivisions, bool) or int(max_subdivisions) != max_subdivisions or max_subdivisions < 2:
        raise ValueError("at least two integer subdivisions are required")
    model = operators.section
    _, _, area, *_ = model._physical(state)
    sy, sn = model.section_widths(state)
    length = math.sqrt(math.pi)/2.*model.k.delta
    q0 = .5*length*length
    pilot = operators.unit_diffusion_response(state, quadrature_points=64)
    scales = np.maximum(np.abs([pilot.transverse_volume, pilot.outward_boundary_flux]), 1.)
    evaluations = [0]

    def volume_integrand(points):
        q = points[:, 0]
        evaluations[0] += len(q)
        # Direct q coordinate, distinct from the square-mapped GL rule.
        angular = 2.*math.pi-8.*np.arccos(np.sqrt(np.minimum(q0/q, 1.)))
        k = operators.radial_diffusion_coefficient(state, q)
        return (area*angular*4.*q*k)[:, None]/scales[0]

    def edge_integrand(points):
        # Physical normal distance on the y face, rescaled only by sigma_n.
        t = points[:, 0]
        q = .5*(length*length+t*t)
        evaluations[0] += len(q)
        k = operators.radial_diffusion_coefficient(state, q)
        shape = 2.*length*length+((sn/sy)**2+(sy/sn)**2)*t*t
        return (4.*length*area*shape*k)[:, None]/scales[1]

    # A margin handles a pilot scale differing from the final magnitude.
    options = dict(atol=tolerance*.25, rtol=0., rule="gk21", max_subdivisions=max_subdivisions)
    volume = cubature(volume_integrand, [0.], [2.*q0], points=[[q0]], **options)
    edge = cubature(edge_integrand, [0.], [length], **options)
    values = np.array([volume.estimate[0], edge.estimate[0]])*scales
    absolute = np.array([volume.error[0], edge.error[0]])*scales
    scaled = absolute/np.maximum(np.abs(values), 1.)
    statuses = (volume.status, edge.status)
    result = SecondMomentBudget(float(values[0]), 0., float(values[1]), 0.)
    return AdaptiveDiffusionResponse(result, tuple(map(float, absolute)), tuple(map(float, scaled)),
        statuses, evaluations[0], tolerance,
        bool(all(s == "converged" for s in statuses) and np.max(scaled) <= tolerance))
