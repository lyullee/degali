"""Conservative density-weighted ensembles; no presumed PDF or turbulence law.

Finite-support linear extrema are NOT bounds over unsampled physical states.
No observation is a constraint and no extremising distribution is a model.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog


def _states(temperature, fraction, enthalpy, specific_volume):
    arrays=tuple(np.asarray(v,dtype=float) for v in (temperature,fraction,enthalpy,specific_volume))
    if (any(a.ndim!=1 for a in arrays) or len({a.shape for a in arrays})!=1
            or len(arrays[0])==0 or any(not np.all(np.isfinite(a)) for a in arrays)):
        raise ValueError('matching nonempty finite one-dimensional state arrays required')
    t,y,h,v=arrays
    if np.any(t<=0) or np.any(v<=0) or np.any(y<0) or np.any(y>1):
        raise ValueError('positive temperature/volume and fraction[0,1] required')
    return arrays


def ensemble(masses, temperature, fraction, enthalpy, specific_volume):
    """Convert explicit nonnegative parcel masses to Favre/Reynolds averages."""
    t,y,h,v=_states(temperature,fraction,enthalpy,specific_volume)
    masses=np.asarray(masses,dtype=float)
    if masses.shape!=t.shape or not np.all(np.isfinite(masses)) or np.any(masses<0):
        raise ValueError('matching finite nonnegative parcel masses required')
    total=float(np.sum(masses))
    if not np.isfinite(total) or total<=0: raise ValueError('finite positive total mass required')
    w=masses/total
    mean_v=float(w@v)
    pv=w*v/mean_v
    ty=float(w@y);tf=float(w@t)
    return dict(mass_weights=w,volume_weights=pv,total_mass=total,total_volume=total*mean_v,
        mean_density=1/mean_v,favre_mass_fraction=ty,favre_enthalpy=float(w@h),
        favre_temperature=tf,reynolds_temperature=float(pv@t),
        favre_fraction_variance=float(w@((y-ty)**2)),
        favre_temperature_variance=float(w@((t-tf)**2)),
        favre_temperature_fraction_covariance=float(w@((t-tf)*(y-ty))))


@dataclass(frozen=True)
class FiniteSupportExtremum:
    objective: float
    weights: np.ndarray
    maximum_scaled_moment_error: float
    maximum_scaled_dual_violation: float
    scaled_duality_gap: float
    moment_residuals: np.ndarray
    objective_replay_error: float


def extreme(fraction, enthalpy, specific_volume, *, mean_fraction, mean_enthalpy,
            mean_specific_volume, values, maximize=False,
            additional_moments=None, additional_targets=None) -> FiniteSupportExtremum:
    """Extremise sum(w*values) at fixed mass, Y, h and volume on this support.

    For Reynolds-temperature use values=T*v/mean_v. For an Eulerian CDF
    use values=v/mean_v*indicator(T<=threshold). Returned PDFs are witnesses.
    """
    y,h,v=(np.asarray(a,dtype=float) for a in (fraction,enthalpy,specific_volume))
    _,y,h,v=_states(np.ones_like(y),y,h,v)
    c=np.asarray(values,dtype=float)
    target=np.array([1.,mean_fraction,mean_enthalpy,mean_specific_volume],dtype=float)
    if (c.shape!=y.shape or not np.all(np.isfinite(c)) or not np.all(np.isfinite(target))
            or not 0<=mean_fraction<=1 or mean_specific_volume<=0):
        raise ValueError('finite objective and physical mean targets required')
    matrix=np.vstack([np.ones_like(y),y,h,v])
    extra_scale=np.empty(0)
    if additional_moments is not None or additional_targets is not None:
        extra=np.asarray(additional_moments,dtype=float)
        targets=np.asarray(additional_targets,dtype=float)
        if (extra.ndim!=2 or extra.shape[1]!=len(y) or targets.shape!=(extra.shape[0],)
                or not np.all(np.isfinite(extra)) or not np.all(np.isfinite(targets))):
            raise ValueError('matching finite additional linear moments and targets required')
        matrix=np.vstack([matrix,extra]);target=np.r_[target,targets]
        extra_scale=np.maximum(abs(targets),1e-6)
    scale=np.maximum(np.max(abs(matrix),axis=1),1e-12)
    a=matrix/scale[:,None];b=target/scale
    cs=max(float(np.max(abs(c))),1.)
    signed=(-1. if maximize else 1.)*c/cs
    result=linprog(signed,A_eq=a,b_eq=b,bounds=(0,None),method='highs',
        options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    if not result.success: raise ValueError('finite-support moment problem failed: '+result.message)
    w=result.x
    raw=matrix@w-target
    check_scale=np.r_[1.,1.,max(abs(mean_enthalpy),1000.),mean_specific_volume,extra_scale]
    error=float(max(abs(raw)/check_scale))
    dual=result.eqlin.marginals
    dual_violation=float(max(0.,np.max(a.T@dual-signed)))
    gap=abs(float(signed@w-b@dual))
    replay=abs(float(signed@w-result.fun))
    if np.any(w<0) or max(error,dual_violation,gap,replay)>1e-8:
        raise RuntimeError('finite-support independent primal/dual check failed')
    return FiniteSupportExtremum(float(c@w),w,error,dual_violation,gap,raw,replay)
