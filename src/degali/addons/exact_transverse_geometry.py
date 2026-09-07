"""Opt-in exact width constraint and consistent complex-step geometry tangent.

No changes to the legacy JetPlume object or its Fortran-compatible root finder.
Branches are evaluated at the real state; derivatives at a branch kink are not
defined by this helper. This is numerical geometry, not a turbulence closure.
"""

import cmath
import copy
import math
import numpy as np
from ..core.constants import VKC
from .buoyancy_profile import BuoyancyConstrainedEnthalpySection
from .transverse_mixing import ConservativeTransverseMixing
from .enriched_segments import decode_enriched


def exact_split(area,sya,sza,*,spread_floor=False):
    area,sya,sza=complex(area),complex(sya),complex(sza)
    if area.real<=0. or min(sya.real,sza.real)<0.:
        raise ValueError('positive area and nonnegative ambient spreads required')
    delta=(sya-sza)*(sya+sza)
    norm=cmath.sqrt(delta*delta+4*area*area)
    if delta.real>=0.:
        sy=cmath.sqrt((norm+delta)/2); sn=area/sy
    else:
        sn=cmath.sqrt((norm-delta)/2); sy=area/sn
    if spread_floor and sy.real<sya.real:
        sy=sya; sn=area/sy
    return sy,sn


def stable_wind_profile(jp,z,sn,ct):
    if jp.ustar<=0.:
        return 0j,0j
    z,sn,ct=complex(z),complex(sn),complex(ct)
    if z.real<jp.zr:
        z=complex(jp.zr)
    top=z+jp.k.delta*sn*ct+.01
    bot=z-jp.k.delta*sn*ct
    if bot.real<jp.zr:
        bot=complex(jp.zr)
    def stability(height):
        if jp.rml<0.:
            aa=(1.-15.*height/jp.rml)**.25
            return 2*cmath.log((1+aa)/2)+cmath.log((1+aa*aa)/2)-2*cmath.atan(aa)+math.pi/2
        return -4.7*height/jp.rml if jp.rml>0. else 0j
    top_u=jp.ustar/VKC*(cmath.log((top+jp.zr)/jp.zr)-stability(top))
    mid_u=jp.ustar/VKC*(cmath.log((z+jp.zr)/jp.zr)-stability(z))
    alpha=cmath.log(top_u/mid_u)/cmath.log(top/z)
    power=1+alpha
    # A difference of powers is evaluated as expm1, including complex-step input.
    integral=bot**power*np.expm1(power*cmath.log(top/bot))
    wind=mid_u/(power*(top-bot)*z**alpha)*integral
    return complex(wind),alpha


def geometry_from_coordinates(jp,coordinates):
    log_area,theta,x,z=map(complex,coordinates)
    area=cmath.exp(log_area)
    if x.real>0.:
        sya=jp.deltay*x**jp.betay
        sza=jp.deltaz*x**jp.betaz*cmath.exp(jp.gammaz*cmath.log(x)**2)
    else:
        sya=sza=0j
    sy,sn=exact_split(area,sya,sza,spread_floor=jp.spread_floor)
    wind,_=stable_wind_profile(jp,z,sn,cmath.cos(theta))
    return np.array([sy,sn,wind],complex)


def exact_geometry_tangent(jp,state,*,step=1e-24):
    state=np.asarray(state,float)
    if state.shape!=(7,) or not np.all(np.isfinite(state)) or state[2]<=0. or not 0.<step<=1e-10:
        raise ValueError('finite physical state and small positive complex step required')
    point=np.array([math.log(state[2]),state[3],state[5],state[6]])
    values=geometry_from_coordinates(jp,point).real
    jacobian=np.column_stack([geometry_from_coordinates(jp,point+1j*step*e).imag/step for e in np.eye(4)])
    # Reject exact nonsmooth switching states rather than claiming a derivative.
    if state[5]==0. or (jp.ustar>0. and state[6]==jp.zr):
        raise ValueError('geometry derivative at branch switch is not defined')
    return dict(values=values,jacobian=jacobian,coordinate_names=['log_area','theta','x','z'],
        exact_width_constraint=True,legacy_modified=False,complex_step=step)


class ExactConstraintGeometry:
    """Read-only geometry override; unrelated physics delegates to original jp."""
    def __init__(self,original):
        self.original=original

    def __getattr__(self,name):
        return getattr(self.original,name)

    def _split(self,area,sya,sza):
        return tuple(v.real for v in exact_split(area,sya,sza,spread_floor=self.spread_floor))

    def _wind_profile(self,z,sn,ct):
        return tuple(v.real for v in stable_wind_profile(self,z,sn,ct))

    def _wind(self,z,sn,ct):
        return self._wind_profile(z,sn,ct)[0]


def exact_geometry_field_view(template,encoded):
    state,parameters=decode_enriched(encoded,template.count)
    old=template.section
    original_jp=old.jetplume.original if isinstance(old.jetplume,ExactConstraintGeometry) else old.jetplume
    section=BuoyancyConstrainedEnthalpySection(ExactConstraintGeometry(original_jp),old.thermodynamics,
        thermal_width_ratio=old.thermal_width_ratio,quadrature_points=old.quadrature_points)
    section._quadrature_cache=old._quadrature_cache
    mixing=ConservativeTransverseMixing(section,state)
    exact=exact_geometry_tangent(section.jetplume,state)
    sy,sn,wind=exact['values']; jac=exact['jacobian']
    if max(abs(np.array([mixing.sy,mixing.sn,mixing.wind])-exact['values']))>1e-12:
        raise ValueError('exact geometry value and tangent path disagree')
    columns=[2,3,5,6]
    mixing.wind_partials[:]=0.; mixing.wind_partials[columns]=jac[2]
    mixing.log_sy_partials[:]=0.; mixing.log_sy_partials[columns]=jac[0]/sy
    mixing.log_sn_partials[:]=0.; mixing.log_sn_partials[columns]=jac[1]/sn
    identity=mixing.log_sy_partials+mixing.log_sn_partials
    expected=np.zeros(8); expected[2]=1.
    if max(abs(identity-expected))>1e-12:
        raise ValueError('exact area-product derivative failed')
    view=copy.copy(template); view.section=section; view.mixing=mixing
    view.hc=section.phase_inverse.enthalpy_and_slope(state[0],state[0]*state[1])[0]
    view._grids,view._edges={},{}
    return view,parameters,exact
