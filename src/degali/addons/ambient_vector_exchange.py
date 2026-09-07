"""Conserved rates carried by entrained ambient air in arbitrary coordinates.

Not drag, pressure work, turbulent closure, or a complete 3D plume model.
Output ordering: mass, H2, Px, Py, Pz, total energy per path length.
"""
import numpy as np


def entrained_ambient_rates(mass_rate,velocity,*,specific_enthalpy):
    v=np.asarray(velocity,dtype=float)
    if v.ndim<1 or v.shape[-1]!=3:
        raise ValueError('ambient velocity needs a final dimension of three')
    m,h=np.broadcast_arrays(np.asarray(mass_rate,float),np.asarray(specific_enthalpy,float))
    shape=np.broadcast_shapes(m.shape,v.shape[:-1])
    m,h=np.broadcast_to(m,shape),np.broadcast_to(h,shape)
    v=np.broadcast_to(v,shape+(3,))
    if any(np.any(~np.isfinite(a)) for a in (m,h,v)) or np.any(m<0):
        raise ValueError('finite state and nonnegative entrained mass rate required')
    out=np.zeros(shape+(6,))
    out[...,0]=m
    out[...,2:5]=m[...,None]*v
    out[...,5]=m*(h+.5*np.sum(v*v,axis=-1))
    return out
