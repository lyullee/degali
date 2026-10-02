"""Opt-in pressure--strain redistribution with explicit realizability checks.

Pressure--strain is trace-free in the incompressible Reynolds-stress budget:
it redistributes TKE among normal/shear components and must not create total
TKE.  This module supplies that identity and a no-clipping PSD time-step
guard.  It is deliberately not wired into the default DEGALI transport.
"""

import numpy as np


def _symmetric_tensor(value, name, *, positive_definite=False,
                      require_positive_semidefinite=True):
    tensor = np.asarray(value, float)
    if tensor.shape != (3, 3) or not np.all(np.isfinite(tensor)):
        raise ValueError(f'{name} must be a finite 3 by 3 tensor')
    if not np.allclose(tensor, tensor.T, rtol=2e-13, atol=2e-13):
        raise ValueError(f'{name} must be symmetric')
    eigenvalues = np.linalg.eigvalsh(tensor)
    if positive_definite:
        if np.any(eigenvalues <= 0.):
            raise ValueError(f'{name} must be positive definite for a finite step bound')
    elif (require_positive_semidefinite
          and np.any(eigenvalues < -1e-12*max(1., np.max(abs(eigenvalues))))):
        raise ValueError(f'{name} must be positive semidefinite')
    return tensor


def slow_pressure_strain(covariance, tke, dissipation_rate,
                         return_to_isotropy_coefficient):
    """Return the explicit Rotta-type slow redistribution tensor.

    ``Pi=-C_phi*(epsilon/k)*(R-2*k*I/3)`` has zero trace.  The user must
    provide every input including ``C_phi``; this is a comparison closure,
    not a parameter selected for cryogenic hydrogen by DEGALI.
    """
    r = _symmetric_tensor(covariance, 'covariance')
    values = [float(v) for v in (tke, dissipation_rate,
                                  return_to_isotropy_coefficient)]
    if not all(np.isfinite(v) and v > 0. for v in values):
        raise ValueError('TKE, dissipation rate and redistribution coefficient must be positive')
    k, epsilon, coefficient = values
    if not np.isclose(np.trace(r), 2.*k, rtol=2e-12, atol=2e-12*max(1., k)):
        raise ValueError('covariance trace must equal 2*TKE')
    pi = -coefficient*epsilon/k*(r - 2.*k*np.eye(3)/3.)
    return dict(pressure_strain=pi, trace=float(np.trace(pi)),
                tke_change_rate=0., closure_inputs_explicit=True,
                adopted=False, physical_closure_validated=False)


def maximum_psd_step(covariance, rate):
    """Largest forward step keeping ``R + dt*dR/dt`` positive semidefinite.

    No tensor is projected, clipped, or silently repaired.  The result is an
    algebraic step bound for an explicitly supplied rate; it is not a time
    integrator or a pressure--strain closure.
    """
    r = _symmetric_tensor(covariance, 'covariance', positive_definite=True)
    derivative = _symmetric_tensor(rate, 'rate', require_positive_semidefinite=False)
    inv_half = np.linalg.inv(np.linalg.cholesky(r))
    normalized = inv_half @ derivative @ inv_half.T
    minimum = float(np.min(np.linalg.eigvalsh(normalized)))
    maximum = np.inf if minimum >= 0. else -1./minimum
    return dict(maximum_step=maximum, limiting_normalized_eigenvalue=minimum,
                covariance_projected=False, realizability_enforced_by_rejection=True,
                adopted=False, physical_closure_validated=False)
