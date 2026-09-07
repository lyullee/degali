"""Reference-covariant binary enthalpy flux for unequal scalar diffusivities.

This is an opt-in constitutive operator, NOT a turbulence closure or EOS.
The caller must provide h_Y at fixed temperature/pressure in the SAME
enthalpy reference as grad(h). Multiphase slip and pressure gradients are
not included. No diffusivity or phase-partition assumption is selected here.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BinaryEnthalpyFlux:
    species_mass_flux: np.ndarray
    heat_only_flux: np.ndarray
    species_enthalpy_flux: np.ndarray
    total_enthalpy_flux: np.ndarray
    naive_enthalpy_flux: np.ndarray
    composition_correction_flux: np.ndarray


def binary_enthalpy_flux(*,density,species_diffusivity,thermal_diffusivity,
                         enthalpy_gradient,mass_fraction_gradient,enthalpy_composition_derivative):
    """Return one spatial flux component, using broadcast-compatible inputs.

    rho:kg/m³; D:m²/s; grad(h):J/(kg*m); grad(Y):1/m; h_Y:J/kg.
    Returned species flux is kg/(m²*s), all other fluxes are W/m².
    For a gas mixture h_Y=h_fuel-h_carrier. More general effective binary
    thermodynamics require an independently justified h_Y and phase transport.
    """
    rho,dy,dt,gh,gy,hy=np.broadcast_arrays(*[
        np.asarray(v,dtype=float) for v in (density,species_diffusivity,thermal_diffusivity,
            enthalpy_gradient,mass_fraction_gradient,enthalpy_composition_derivative)])
    if not all(np.all(np.isfinite(v)) for v in (rho,dy,dt,gh,gy,hy)):
        raise ValueError('all binary flux inputs must be finite')
    if np.any(rho<=0) or np.any(dy<0) or np.any(dt<0):
        raise ValueError('positive density and nonnegative diffusivities required')
    jy=-rho*dy*gy
    heat=-rho*dt*(gh-hy*gy)
    carried=hy*jy
    naive=-rho*dt*gh
    correction=-rho*(dy-dt)*hy*gy
    # Sum the physically separate terms. Equality to naive+correction is
    # independently tested; equal-diffusivity recovery is up to float roundoff.
    return BinaryEnthalpyFlux(jy,heat,carried,heat+carried,naive,correction)
