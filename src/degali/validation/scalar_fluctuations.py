"""Scalar intermittency diagnostics; not a stochastic dispersion closure.

The public cryogenic Raman experiment reports a relation between relative
scalar RMS and plume occupancy.  These functions preserve the relation and
reject impossible inputs without converting a steady DEGALI mean field into a
fitted time series.
"""

import numpy as np


def intermittency_from_relative_rms(relative_rms, shape_factor):
    """Evaluate ``gamma=(K+1)/(1+(f_rms/f_mean)^2)`` without clipping.

    ``K`` is an explicitly supplied experiment-specific shape factor. A
    returned occupancy outside ``[0, 1]`` is reported as incompatible, rather
    than silently forced into a probability.
    """
    ratio, factor = np.broadcast_arrays(
        np.asarray(relative_rms, float), np.asarray(shape_factor, float)
    )
    if (not np.all(np.isfinite(ratio)) or not np.all(np.isfinite(factor))
            or np.any(ratio < 0.) or np.any(factor < 0.)):
        raise ValueError('relative RMS and shape factor must be finite and nonnegative')
    occupancy = (1.+factor)/(1.+ratio*ratio)
    physical = (occupancy >= 0.) & (occupancy <= 1.)
    return dict(occupancy=occupancy, occupancy_physical=physical,
                relative_rms=ratio, shape_factor=factor,
                clipped=False, adopted=False, stochastic_closure_validated=False)


def conditional_active_scalar(mean_scalar, occupancy):
    """Return active-state mean implied by a positive supplied occupancy.

    This is the identity ``mean=occupancy*conditional_mean``. It is useful for
    reporting a two-state envelope, not a claim that a cryogenic plume has a
    two-state PDF.
    """
    mean, gamma = np.broadcast_arrays(
        np.asarray(mean_scalar, float), np.asarray(occupancy, float)
    )
    if (not np.all(np.isfinite(mean)) or not np.all(np.isfinite(gamma))
            or np.any(mean < 0.) or np.any(gamma <= 0.) or np.any(gamma > 1.)):
        raise ValueError('mean must be nonnegative and occupancy must lie in (0, 1]')
    return dict(conditional_active_mean=mean/gamma, mean=mean, occupancy=gamma,
                two_state_interpretation_only=True, adopted=False)
