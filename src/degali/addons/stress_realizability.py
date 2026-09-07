"""Necessary TKE bound from known axial/transverse Reynolds shear covariance.

This supplies no turbulence closure or fitted coefficient. For a positive
semidefinite velocity covariance R, k=trace(R)/2 >= norm((R_sy,R_sn)).
The bound is sharp without other specified covariance entries.
"""
import numpy as np


def shear_tke_lower_bound(normalized_stress,area,density,transverse_lengths):
    stress=np.asarray(normalized_stress,float)
    rho=np.asarray(density,float)
    lengths=np.asarray(transverse_lengths,float)
    if stress.ndim!=2 or stress.shape[1]!=2 or rho.shape!=(len(stress),) or lengths.shape!=(2,):
        raise ValueError('two stress components, matching density and two metric lengths required')
    if not np.isfinite(area) or area<=0. or np.any(rho<=0.) or np.any(lengths<=0.) or not all(np.all(np.isfinite(v)) for v in (stress,rho,lengths)):
        raise ValueError('finite stresses and positive physical area/density/lengths required')
    # FM=A*rho*v_relative/L; FP-u*FM=A*rho*R_si/L in this reduced system.
    shear=stress*lengths[None,:]/(area*rho[:,None])
    minimum=np.linalg.norm(shear,axis=1)
    factor=np.zeros((len(stress),3))
    positive=minimum>0.
    factor[positive,0]=np.sqrt(minimum[positive])
    factor[positive,1:]=shear[positive]/factor[positive,0,None]
    covariance=np.einsum('ni,nj->nij',factor,factor)
    return dict(minimum_tke=minimum,physical_shear_covariance=shear,
        minimum_total_fluctuation_rms=np.sqrt(2*minimum),attaining_covariance=covariance,
        covariance_is_physical_closure=False,normal_stresses_not_measured=True)
