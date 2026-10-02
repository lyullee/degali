"""Transparent finite-TKE parameter conversions; none selects a closure.

The reduced transport stores ``Q=rho*k`` and needs a TKE diffusivity and a
dissipation time. These helpers express conventional k-epsilon relations,
but require every physical scale and coefficient from the caller.
"""

import numpy as np

from .normal_stress_budget import (
    minimum_tke_from_axial_rms,
    minimum_tke_from_two_normal_rms,
)
from .transport_evidence import TURBULENCE_BOUNDARY


def _matching_positive(*fields):
    """Return broadcast fields after rejecting nonphysical inputs."""
    values = [np.asarray(v, float) for v in fields]
    try:
        values = np.broadcast_arrays(*values)
    except ValueError as exc:
        raise ValueError('physical fields must be broadcast-compatible') from exc
    if not all(np.all(np.isfinite(v)) and np.all(v > 0.) for v in values):
        raise ValueError('all physical fields must be finite and strictly positive')
    return values


def dissipation_time_from_integral_scale(integral_length, tke,
                                         dissipation_coefficient):
    """Return ``tau=L/(C_e sqrt(k))`` from explicit scale inputs.

    This equals ``epsilon=C_e*k**1.5/L`` and ``tau=k/epsilon``. ``C_e`` is
    neither fitted nor defaulted. In conventional high-Re k-epsilon notation,
    it is often written ``C_mu**(3/4)``; transfer to a cryogenic jet remains
    a study hypothesis.
    """
    length, energy, coefficient = _matching_positive(
        integral_length, tke, dissipation_coefficient
    )
    tau = length/(coefficient*np.sqrt(energy))
    return dict(dissipation_time=tau, integral_length=length,
                tke=energy, dissipation_coefficient=coefficient,
                closure_inputs_explicit=True, adopted=False,
                **TURBULENCE_BOUNDARY.audit_fields())


def tke_diffusivity_from_eddy_viscosity(eddy_viscosity, tke_prandtl):
    """Return ``chi_k=nu_t/sigma_k`` with explicit transport-ratio input.

    ``sigma_k`` is sometimes called the TKE turbulent Prandtl number. This
    does not equate it to a scalar Schmidt number or provide a default.
    """
    viscosity, prandtl = _matching_positive(eddy_viscosity, tke_prandtl)
    return dict(tke_diffusivity=viscosity/prandtl,
                eddy_viscosity=viscosity, tke_prandtl=prandtl,
                closure_inputs_explicit=True, adopted=False,
                **TURBULENCE_BOUNDARY.audit_fields())


def k_epsilon_dissipation_coefficient(c_mu):
    """Algebraically convert a supplied k-epsilon ``C_mu`` to ``C_e``."""
    (coefficient,) = _matching_positive(c_mu)
    return coefficient**.75


def tke_realizability_margin(tke, axial_rms, shear_covariance):
    """Screen an explicit TKE field against measured axial-RMS information.

    This is deliberately only a pass/fail diagnostic.  The sharp
    positive-semidefinite lower bound is evaluated from the supplied axial
    RMS and axial/transverse shear covariances; a deficient ``k`` is reported
    as deficient and is never raised to the bound.  It supplies neither an
    isotropy completion nor an epsilon/length-scale closure.
    """
    energy = np.asarray(tke, float)
    if energy.ndim != 1 or not np.all(np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("TKE must be a one-dimensional finite nonnegative field")
    bound = minimum_tke_from_axial_rms(shear_covariance, axial_rms)
    required = bound["minimum_tke"]
    if energy.shape != required.shape:
        raise ValueError("TKE and axial-RMS realizability fields must have matching shape")
    margin = energy - required
    return dict(
        tke=energy,
        minimum_tke=required,
        realizability_margin=margin,
        realizable=bool(np.all(margin >= 0.0)),
        uses_isotropy=False,
        adopted=False,
        **TURBULENCE_BOUNDARY.audit_fields(),
    )


def tke_two_normal_rms_realizability_margin(tke, axial_rms, transverse_rms):
    """Screen TKE against two observed normal-RMS components without repair.

    This companion to :func:`tke_realizability_margin` is for a published
    axial/radial RMS profile with no shear covariance.  It applies only the
    exact PSD trace lower bound and never creates the missing azimuthal
    stress, covariance, length scale, or epsilon.
    """
    energy = np.asarray(tke, float)
    if energy.ndim != 1 or not np.all(np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("TKE must be a one-dimensional finite nonnegative field")
    bound = minimum_tke_from_two_normal_rms(axial_rms, transverse_rms)
    required = bound["minimum_tke"]
    if energy.shape != required.shape:
        raise ValueError("TKE and normal-RMS realizability fields must have matching shape")
    margin = energy - required
    return dict(
        tke=energy,
        minimum_tke=required,
        realizability_margin=margin,
        realizable=bool(np.all(margin >= 0.0)),
        uses_isotropy=False,
        adopted=False,
        **TURBULENCE_BOUNDARY.audit_fields(),
    )
