"""Actual-value enthalpy differences with high-precision scalar prefactors.

Optional mpmath is confined to this verification helper. The EOS table and
legacy geometry remain unchanged. Exact-root wind is a separately labelled
evaluation of the same width constraint, not a replacement of JETPLU.
"""

import math
import numpy as np
from .enriched_segments import encode_enriched,shifted_field_view
from .phase_radial_quadrature import gauss_rule
from ..core.constants import VKC


def mp_center_enthalpy(inverter,log_density,log_fuel,mp):
    rho,c=mp.exp(log_density),mp.exp(log_fuel)
    y=c/rho; ti=mp.mpf(inverter.a)/(rho+mp.mpf(inverter.k)*c)
    if not (mp.mpf(inverter.t_grid[0])<=ti<=mp.mpf(inverter.t_grid[-1]) and 0<=y<=1):
        raise ValueError('high-precision center outside unchanged phase table')
    i=int(np.clip(np.searchsorted(inverter.t_grid,float(ti),side='right')-1,0,len(inverter.t_grid)-2))
    j=int(np.clip(np.searchsorted(inverter.y_grid,float(y),side='right')-1,0,len(inverter.y_grid)-2))
    f=(ti-mp.mpf(inverter.t_grid[i]))/(mp.mpf(inverter.t_grid[i+1])-mp.mpf(inverter.t_grid[i]))
    g=(y-mp.mpf(inverter.y_grid[j]))/(mp.mpf(inverter.y_grid[j+1])-mp.mpf(inverter.y_grid[j]))
    h00,h01,h10,h11=[mp.mpf(inverter.h_values[ii,jj]) for ii,jj in ((i,j),(i,j+1),(i+1,j),(i+1,j+1))]
    return h00+g*(h01-h00)+f*(h10+g*(h11-h10)-h00-g*(h01-h00))


def mp_wind(section,encoded,mp):
    jp=section.jetplume
    area,theta,x,z=mp.exp(encoded[2]),encoded[3],encoded[5],encoded[6]
    if jp.ustar<=0.:
        return mp.mpf(0),mp.sqrt(area),mp.sqrt(area)
    sya=mp.mpf(jp.deltay)*max(x,mp.mpf(0))**mp.mpf(jp.betay)
    sza=(mp.mpf(jp.deltaz)*x**mp.mpf(jp.betaz)*mp.exp(mp.mpf(jp.gammaz)*mp.log(x)**2)
         if x>0 else mp.mpf(0))
    delta=sya*sya-sza*sza; norm=mp.sqrt(delta*delta+4*area*area)
    if delta>=0:
        sy=mp.sqrt((norm+delta)/2); sn=area/sy
    else:
        sn=mp.sqrt((norm-delta)/2); sy=area/sn
    if jp.spread_floor and sy<sya:
        sy=sya; sn=area/sy
    zr=mp.mpf(jp.zr); z=max(z,zr); ct=mp.cos(theta)
    top=z+mp.mpf(jp.k.delta)*sn*ct+mp.mpf(.01)
    bot=max(z-mp.mpf(jp.k.delta)*sn*ct,zr)
    rml=mp.mpf(jp.rml)
    def stability(height):
        if rml<0:
            aa=(1-15*height/rml)**mp.mpf(.25)
            # Preserve the original rounded math.pi constant, not a new one.
            return 2*mp.log((1+aa)/2)+mp.log((1+aa*aa)/2)-2*mp.atan(aa)+mp.mpf(math.pi)/2
        return -mp.mpf(4.7)*height/rml if rml>0 else mp.mpf(0)
    top_u=mp.mpf(jp.ustar)/mp.mpf(VKC)*(mp.log((top+zr)/zr)-stability(top))
    mid_u=mp.mpf(jp.ustar)/mp.mpf(VKC)*(mp.log((z+zr)/zr)-stability(z))
    alpha=mp.log(top_u/mid_u)/mp.log(top/z)
    wind=mid_u/((1+alpha)*(top-bot)*z**alpha)*(top**(1+alpha)-bot**(1+alpha))
    return wind,sy,sn


def scalar_pair(projection,parameters,rates,step,*,precision=50,wind_policy='binary_legacy'):
    from mpmath import mp
    if precision<35 or step<=0. or wind_policy not in ('binary_legacy','exact_constraint'):
        raise ValueError('35+ digits, positive step and named wind policy required')
    origin=encode_enriched(projection,parameters)
    rates=np.asarray(rates,float)
    if rates.shape!=origin.shape or not np.all(np.isfinite(rates)):
        raise ValueError('finite matching direction required')
    with mp.workdps(precision):
        states=[[mp.mpf(float(x))+sign*mp.mpf(step)*mp.mpf(float(r)) for x,r in zip(origin,rates)] for sign in (1,-1)]
        values=[]
        for encoded in states:
            h=mp_center_enthalpy(projection.section.phase_inverse,encoded[0],encoded[1],mp)
            area,uc,theta=mp.exp(encoded[2]),mp.exp(encoded[4]),encoded[3]
            exact,sy,sn=mp_wind(projection.section,encoded,mp)
            view,_=shifted_field_view(projection,np.array([float(v) for v in encoded]))
            wind=mp.mpf(view.mixing.wind) if wind_policy=='binary_legacy' else exact
            values.append(dict(ah=area*h,ahw=area*h*wind*mp.cos(theta),ahu=area*h*uc,
                h=h,area=area,wind=wind,parallel=wind*mp.cos(theta),uc=uc,
                legacy_wind=mp.mpf(view.mixing.wind),exact_wind=exact,sy=sy,sn=sn,
                legacy_sy=mp.mpf(view.mixing.sy),legacy_sn=mp.mpf(view.mixing.sn)))
        pairs={key:dict(mean=float((values[0][key]+values[1][key])/2),
            difference=float(values[0][key]-values[1][key])) for key in ('ah','ahw','ahu','area','parallel','uc','h','wind')}
        pairs.update(precision_digits=precision,wind_policy=wind_policy,
            maximum_wind_value_change=float(max(abs(v['exact_wind']-v['legacy_wind']) for v in values)),
            maximum_width_value_change=float(max(abs(v[key]-v['legacy_'+key]) for v in values for key in ('sy','sn'))))
    return pairs


def stable_enthalpy_difference(projection,parameters,rates,step,*,order=64,precision=50,wind_policy='binary_legacy'):
    pair=scalar_pair(projection,parameters,rates,step,precision=precision,wind_policy=wind_policy)
    p=projection; size=p.basis.size
    x,w=gauss_rule(order); a,b=np.meshgrid(x,x,indexing='ij')
    prep=p.prepare(a.ravel(),b.ravel()); q=prep['q']; psi=prep['psi']
    base=np.exp(-q/p.section.thermal_width_ratio**2+psi@parameters[size:])
    shift=step*(psi@rates[7+size:]); cs,sn=np.cosh(shift),np.sinh(shift)
    velocity_shape=np.exp(-p.section.velocity_shape_exponent*q)
    dh=base*(pair['ahw']['difference']*cs+2*pair['ahw']['mean']*sn
        +velocity_shape*(pair['ahu']['difference']*cs+2*pair['ahu']['mean']*sn))
    weight=(8*p.q0*np.outer(w,w)).ravel()
    derivative=math.fsum(weight*dh)/(2*step)
    return dict(derivative=derivative,order=order,scalar_pair=pair,actual_values_not_EOS_derivatives=True,
        density_or_kinetic_verified=False,all_fields_arbitrary_precision=False)
