"""Mass circulation with a conservative, direction-preserving momentum flux.

Only a kinematic construction. It does not determine circulation amplitudes
or prove the resulting scalar momentum diffusivity is nonnegative.
"""

import math
import numpy as np
from numpy.polynomial.legendre import leg2poly
from scipy.special import gammainc,gamma
from .phase_radial_quadrature import gauss_rule
from .solenoidal_transport import solenoidal_basis,stream_polynomial


def aligned_mode_vectors(a,b,*,q0,ambient_parallel,excess_velocity,velocity_exponent,flux_modes=4,order=32):
    a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float))
    scalars=np.array([q0,ambient_parallel,excess_velocity,velocity_exponent])
    if not np.all(np.isfinite(scalars)) or q0<=0. or excess_velocity<0. or velocity_exponent<=0.:
        raise ValueError('finite positive geometry/exponent and nonnegative excess speed required')
    vm=solenoidal_basis(a,b,flux_modes)
    q=q0*(a*a+b*b)
    x,w=gauss_rule(order)
    aq=(a[:,None]*np.sqrt(x)).ravel()
    bq=(b[:,None]*np.sqrt(x)).ravel()
    vq=solenoidal_basis(aq,bq,flux_modes)
    uq=-velocity_exponent*excess_velocity*np.exp(-velocity_exponent*(q[:,None]*x).ravel())
    gradient=(2*q0*uq)[:,None]*np.column_stack([aq,bq])
    integrand=np.einsum('ni,nik->nk',gradient,vq).reshape(len(a),order,flux_modes)
    # r=-integral_0^q(grad(u).V)dq/(2q). Mapping to [0,1] cancels q,
    # giving the regular value0 at the origin without a division by zero.
    radial=-.5*np.einsum('j,njk->nk',w,integrand)
    u=ambient_parallel+excess_velocity*np.exp(-velocity_exponent*q)
    vp=u[:,None,None]*vm+np.column_stack([a,b])[:,:,None]*radial[:,None,:]
    return dict(mass_modes=vm,momentum_modes=vp,radial_stress_modes=radial,
        circulation_amplitudes_determined=False,positive_diffusivity_established=False)


def aligned_mode_vectors_analytic(a,b,*,q0,ambient_parallel,excess_velocity,velocity_exponent,flux_modes=4):
    """Same radial completion via exact Gaussian polynomial integrals.

    A small-argument series avoids cancellation/0/0 in the lower gamma ratio.
    No EOS or momentum-diffusivity value is introduced by this calculation.
    """
    a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float))
    scalars=np.array([q0,ambient_parallel,excess_velocity,velocity_exponent])
    if not np.all(np.isfinite(scalars)) or q0<=0. or excess_velocity<0. or velocity_exponent<=0.:
        raise ValueError('finite physical geometry and velocity parameters required')
    vm=solenoidal_basis(a,b,flux_modes)
    q=q0*(a*a+b*b); z=velocity_exponent*q
    radial=np.zeros((len(a),flux_modes))
    for k in range(1,flux_modes+1):
        coefficients=leg2poly(stream_polynomial(k)[0])
        for j in range(1,k+1):
            shape=a*a*b**(2*j)+b*b*a**(2*j)-(a**(2*j+2)+b**(2*j+2))/(2*j+1)
            integral=np.empty_like(z)
            small=z<1e-3
            term=np.ones_like(z[small]); series=term/(j+2)
            for power in range(1,9):
                term*=(-z[small])/power
                series+=term/(j+2+power)
            integral[small]=series
            integral[~small]=gamma(j+2)*gammainc(j+2,z[~small])/z[~small]**(j+2)
            radial[:,k-1]+=q0*velocity_exponent*excess_velocity*coefficients[2*j]*shape*integral
    u=ambient_parallel+excess_velocity*np.exp(-z)
    vp=u[:,None,None]*vm+np.column_stack([a,b])[:,:,None]*radial[:,None,:]
    return dict(mass_modes=vm,momentum_modes=vp,radial_stress_modes=radial,
        circulation_amplitudes_determined=False,positive_diffusivity_established=False)
