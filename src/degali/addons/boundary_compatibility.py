"""Rate-free necessary compatibility of three reservoir boundary conditions.

Diagnostic only: optimal fluxes are NOT a conservative transverse solution.
"""

import math
import numpy as np
from scipy.optimize import linprog


def boundary_compatibility(*, y, specific_h, velocity, ambient_speed, theta,
                           diffusion_y, diffusion_h, scales):
    """Sharp unconstrained and inflow/nonnegative-stress minimax bounds.

    Residual order C,H,P. diffusion inputs are +A*rho*chi*normal_gradient/(2q0).
    Positive inferred momentum diffusivity needs fP-u*fM >= 0 when u_q < 0.
    """
    values=np.array([y,specific_h,velocity,ambient_speed,theta,diffusion_y,diffusion_h],float)
    scales=np.asarray(scales,float)
    if not np.all(np.isfinite(values)) or not 0.<y<1. or scales.shape!=(3,) or np.any(scales<=0.) or not np.all(np.isfinite(scales)):
        raise ValueError('finite physical fraction and positive three residual scales required')
    parallel=ambient_speed*math.cos(theta)
    effective_h=specific_h-.5*(ambient_speed**2+velocity**2)+velocity*parallel
    matrix=np.array([[y,0.],[specific_h-.5*(ambient_speed**2+velocity**2),velocity],[-parallel,1.]])
    target=np.array([diffusion_y,diffusion_h,0.])
    null=np.array([-effective_h/y,1.,-velocity])*scales
    compatibility=diffusion_h-effective_h/y*diffusion_y
    amplitude=-compatibility/np.sum(abs(null))
    optimal_residual=amplitude*np.sign(null)
    mass=(diffusion_y+scales[0]*optimal_residual[0])/y
    momentum=parallel*mass+scales[2]*optimal_residual[2]
    actual=(matrix@np.array([mass,momentum])-target)/scales
    b=matrix/scales[:,None]
    g=target/scales
    # t bounds all three normalized residuals. This local LP does not impose
    # interior conservation; even its constrained optimum is only a LOWER bound.
    aub=np.vstack([np.c_[b,-np.ones(3)],np.c_[-b,-np.ones(3)],np.array([[velocity,-1.,0.]])])
    bub=np.r_[g,-g,0.]
    out=linprog([0.,0.,1.],A_ub=aub,b_ub=bub,bounds=[(None,0.),(None,None),(0.,None)],method='highs')
    if not out.success:
        raise ValueError(f'local compatibility LP failed: {out.message}')
    return dict(effective_specific_h=effective_h,compatibility_numerator=compatibility,
        unavoidable_minimax=abs(amplitude),unconstrained_optimal_flux=np.array([mass,momentum]),
        unconstrained_residual=actual,unconstrained_has_inflow=bool(mass<0.),
        unconstrained_has_nonnegative_stress=bool(momentum-velocity*mass>=0.),
        physical_minimax=float(out.fun),physical_optimal_flux=out.x[:2],
        physical_residual=(matrix@out.x[:2]-target)/scales,
        scaled_left_null=null,scaled_matrix=b,scaled_target=g,
        sufficient_for_conservative_flow=False)
