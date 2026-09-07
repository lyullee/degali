"""Verification-only Newton polishing of the SAME bilinear phase oracle.

No new thermodynamic values, changed root branch, or default replacement.
"""

import numpy as np
from .enthalpy_profile import PhaseMassEnthalpyInverter


class PolishedPhaseMassEnthalpyInverter(PhaseMassEnthalpyInverter):
    def __init__(self,thermodynamics):
        super().__init__(thermodynamics)
        self.statistics=dict(calls=0,points=0,improved_points=0,
            maximum_relative_residual_before=0.,maximum_relative_residual_after=0.)

    def state(self,fuel_density,enthalpy_density,*,density_guess=None):
        rho,_,_=super().state(fuel_density,enthalpy_density,density_guess=density_guess)
        c,h=np.broadcast_arrays(np.asarray(fuel_density,float),np.asarray(enthalpy_density,float))
        shape=c.shape
        c,h,rho=c.ravel(),h.ravel(),rho.ravel().copy()
        lower=np.maximum(c,self.a/self.t_grid[-1]-self.k*c)
        upper=self.a/self.t_grid[0]-self.k*c
        value,slope=self.enthalpy_and_slope(rho,c)
        error=value-h
        before=abs(error).copy()
        for _ in range(3):
            candidate=rho-error/slope
            valid=(candidate>=lower)&(candidate<=upper)&np.isfinite(candidate)&(slope<0.)
            candidate=np.where(valid,candidate,rho)
            new_value,new_slope=self.enthalpy_and_slope(candidate,c)
            new_error=new_value-h
            better=abs(new_error)<abs(error)
            rho,error,slope=np.where(better,candidate,rho),np.where(better,new_error,error),np.where(better,new_slope,slope)
            # If Newton rounds back to the same float, test its adjacent float
            # in the sign-correct direction; accept only a strictly better root.
            neighbor=np.nextafter(rho,np.where(error>0.,np.inf,-np.inf))
            valid=(neighbor>=lower)&(neighbor<=upper)&(error!=0.)
            neighbor=np.where(valid,neighbor,rho)
            new_value,new_slope=self.enthalpy_and_slope(neighbor,c)
            new_error=new_value-h
            better=abs(new_error)<abs(error)
            rho,error,slope=np.where(better,neighbor,rho),np.where(better,new_error,error),np.where(better,new_slope,slope)
        if np.any(abs(error)>before):
            raise RuntimeError('verification polishing made a phase root worse')
        y=c/rho
        temperature,represented=self.th._condensed_air_state(rho,y)
        # Preserve the original physical branch/representation check.
        if np.any(abs(represented-h)>1.01e-10*np.maximum(abs(h),1.)):
            raise RuntimeError('polished phase root disagrees with the original property oracle')
        stat=self.statistics
        stat['calls']+=1; stat['points']+=len(c)
        stat['improved_points']+=int(np.count_nonzero(abs(error)<before))
        stat['maximum_relative_residual_before']=max(stat['maximum_relative_residual_before'],float(np.max(before/np.maximum(abs(h),1.),initial=0.)))
        stat['maximum_relative_residual_after']=max(stat['maximum_relative_residual_after'],float(np.max(abs(error)/np.maximum(abs(h),1.),initial=0.)))
        return rho.reshape(shape),y.reshape(shape),temperature.reshape(shape)
