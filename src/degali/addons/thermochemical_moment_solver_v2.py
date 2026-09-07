"""Conditioned finite-support moment solver; same original physical gates.

SVD changes only the equality-row basis, not support, objective or constraints.
Original-unit primal and original scaled dual certificates are still required.
"""
import numpy as np
from scipy.optimize import linprog

from .thermochemical_ensemble import _states,FiniteSupportExtremum


def whitened_extreme(fraction, enthalpy, specific_volume, *, mean_fraction,
                     mean_enthalpy,mean_specific_volume,values,maximize=False,
                     additional_moments=None,additional_targets=None):
    y,h,v=(np.asarray(a,dtype=float) for a in (fraction,enthalpy,specific_volume))
    _,y,h,v=_states(np.ones_like(y),y,h,v)
    c=np.asarray(values,dtype=float)
    target=np.array([1.,mean_fraction,mean_enthalpy,mean_specific_volume],dtype=float)
    if (c.shape!=y.shape or not np.all(np.isfinite(c)) or not np.all(np.isfinite(target))
            or not 0<=mean_fraction<=1 or mean_specific_volume<=0):
        raise ValueError('finite objective and physical mean targets required')
    matrix=np.vstack([np.ones_like(y),y,h,v]);extra_scale=np.empty(0)
    if additional_moments is not None or additional_targets is not None:
        extra=np.asarray(additional_moments,dtype=float);targets=np.asarray(additional_targets,dtype=float)
        if (extra.ndim!=2 or extra.shape[1]!=len(y) or targets.shape!=(extra.shape[0],)
                or not np.all(np.isfinite(extra)) or not np.all(np.isfinite(targets))):
            raise ValueError('matching finite additional moments/targets required')
        matrix=np.vstack([matrix,extra]);target=np.r_[target,targets]
        extra_scale=np.maximum(abs(targets),1e-6)
    scale=np.maximum(np.max(abs(matrix),axis=1),1e-12)
    a,b=matrix/scale[:,None],target/scale
    u,s,_=np.linalg.svd(a,full_matrices=False)
    keep=s>1e-12*s[0]
    transform=(u[:,keep]/s[keep]).T
    aw,bw=transform@a,transform@b
    # Dropping a mathematically redundant row cannot change its target.
    if np.max(abs(b-u[:,keep]@(u[:,keep].T@b)))>1e-10:
        raise ValueError('targets incompatible with dependent equality rows')
    cs=max(float(np.max(abs(c))),1.);signed=(-1. if maximize else 1.)*c/cs
    result=linprog(signed,A_eq=aw,b_eq=bw,bounds=(0,None),method='highs',
        options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
    if not result.success: raise ValueError('conditioned moment problem failed: '+result.message)
    w=result.x;raw=matrix@w-target
    check_scale=np.r_[1.,1.,max(abs(mean_enthalpy),1000.),mean_specific_volume,extra_scale]
    error=float(max(abs(raw)/check_scale))
    dual=transform.T@result.eqlin.marginals
    dual_violation=float(max(0.,np.max(a.T@dual-signed)))
    gap=abs(float(signed@w-b@dual));replay=abs(float(signed@w-result.fun))
    if np.any(w<0) or max(error,dual_violation,gap,replay)>1e-8:
        raise RuntimeError(f'conditioned original-gate failure: primal={error},dual={dual_violation},gap={gap},minw={min(w)}')
    return FiniteSupportExtremum(float(c@w),w,error,dual_violation,gap,raw,replay)
