"""Bounded equimolar N2/H2 reference data and conditional binary gas fluxes.

Not a humid/multiphase diffusivity, a turbulent closure, or an energy model.
Van Heijningen 1967, Chapter I Table I, printed p11; p15 flags the 65 K datum.
"""
import numpy as np

TEMPERATURE_K = np.array([65.25,77.35,90.2,169.3,294.8])
DP_M2_S_ATM = np.array([4.70,6.71,9.00,28.94,76.64])*1e-6
DP_PRINTED_90_INTERVAL = np.array([.04,.02,.05,.15,.20])*1e-6
for _table in (TEMPERATURE_K,DP_M2_S_ATM,DP_PRINTED_90_INTERVAL):
    _table.setflags(write=False)


def equimolar_reference_diffusivity(temperature,pressure,*,allow_flagged_65K=False):
    """Log-linear table interpolant in m²/s; dilute-gas inverse-P assumption.

    No correction for composition, O2/water, condensed matter, or dense gas.
    Published +/- values are not interpolated into a fictitious uncertainty.
    """
    t,p=np.broadcast_arrays(np.asarray(temperature,float),np.asarray(pressure,float))
    low=0 if allow_flagged_65K else 1
    if np.any(~np.isfinite(t)) or np.any(~np.isfinite(p)) or np.any(p<=0):
        raise ValueError('finite temperature and positive absolute pressure required')
    if np.any(t<TEMPERATURE_K[low]) or np.any(t>TEMPERATURE_K[-1]):
        raise ValueError('outside explicitly selected reference temperature range')
    return np.exp(np.interp(np.log(t),np.log(TEMPERATURE_K[low:]),np.log(DP_M2_S_ATM[low:])))*101325./p


def binary_gas_species_flux(*,density,diffusivity,mass_fraction,
        mass_fraction_gradient,log_temperature_gradient,thermal_diffusion_factor):
    """Mass-frame binary Fick+Soret, with a REQUIRED signed alpha convention.

    j1 = -rho D [gradY + alpha Y(1-Y) grad(log T)]; j2=-j1.
    Coefficients are caller inputs. No thermal/Dufour or phase-exchange closure.
    """
    rho,d,y,gy,gt,a=np.broadcast_arrays(*[np.asarray(v,float) for v in (
        density,diffusivity,mass_fraction,mass_fraction_gradient,
        log_temperature_gradient,thermal_diffusion_factor)])
    if any(np.any(~np.isfinite(v)) for v in (rho,d,y,gy,gt,a)):
        raise ValueError('finite values required')
    if np.any(rho<=0) or np.any(d<0) or np.any(y<0) or np.any(y>1):
        raise ValueError('positive density, nonnegative D and bounded Y required')
    fick=-rho*d*gy
    soret=-rho*d*a*y*(1-y)*gt
    j1=fick+soret
    return dict(fick=fick,soret=soret,species1=j1,species2=-j1)


def constant_reference_spreading(time,width,diffusivity):
    """Conditional one-dimensional Fick variance scale, not a plume solver."""
    t,w,d=np.broadcast_arrays(*[np.asarray(v,float) for v in (time,width,diffusivity)])
    if any(np.any(~np.isfinite(v)) for v in (t,w,d)) or np.any(t<0) or np.any(w<=0) or np.any(d<0):
        raise ValueError('nonnegative time/D and positive width required')
    variance=2*d*t
    r=variance/w**2
    # Stable evaluation of sqrt(1+r)-1 without cancellation at small r.
    return dict(diffusion_length=np.sqrt(variance),variance_ratio=r,
        relative_width_increment=r/(np.sqrt(1+r)+1),
        diffusivity_for_one_percent_width=np.divide((1.01**2-1)*w**2,2*t,
            out=np.full_like(t,np.inf),where=t>0))
