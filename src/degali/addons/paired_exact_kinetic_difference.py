"""Actual density/kinetic values on a union phase mesh with exact geometry.

Density remains binary64 with same-table root polishing. Small area and speed
differences use high-precision scalar values, not analytic EOS derivatives.
"""
import math
import numpy as np
from .enriched_segments import encode_enriched
from .exact_transverse_geometry import exact_geometry_field_view
from .edge_conservative_refit import FaceSplitSquareMoments
from .phase_radial_quadrature import gauss_rule
from .polished_phase_inverse import PolishedPhaseMassEnthalpyInverter
from .stable_enthalpy_difference import scalar_pair


def _union(values,scale):
    out=np.sort(np.concatenate(values))
    return out[np.r_[True,np.diff(out)>scale*2e-13]]


def paired_exact_kinetic_difference(projection,parameters,rates,step,*,order=8,angular_order=16,batch_size=8,callback=None):
    origin=encode_enriched(projection,parameters)
    rates=np.asarray(rates,float)
    if rates.shape!=origin.shape or not np.all(np.isfinite(rates)) or not math.isfinite(step) or step<=0.:
        raise ValueError('finite matching rates and positive difference distance required')
    if not isinstance(batch_size,int) or isinstance(batch_size,bool) or batch_size<1:
        raise ValueError('positive batch size required')
    plus,pp,_=exact_geometry_field_view(projection,origin+step*rates)
    minus,pm,_=exact_geometry_field_view(projection,origin-step*rates)
    plus.section.phase_inverse=PolishedPhaseMassEnthalpyInverter(plus.section.thermodynamics)
    minus.section.phase_inverse=PolishedPhaseMassEnthalpyInverter(minus.section.thermodynamics)
    scalar=scalar_pair(projection,parameters,rates,step,precision=70,wind_policy='exact_constraint')
    qp,qm=FaceSplitSquareMoments(plus),FaceSplitSquareMoments(minus)
    ak=_union([qp.angular_knots(pp),qm.angular_knots(pm)],1.); ak[-1]=math.pi/4
    ax,aw0=gauss_rule(angular_order); x,w=gauss_rule(order)
    angles=(ak[:-1,None]+np.diff(ak)[:,None]*ax).ravel(); aw=(np.diff(ak)[:,None]*aw0).ravel()
    sums=[]; peaks=0; absolute=0.
    for start in range(0,len(angles),batch_size):
        phi=angles[start:start+batch_size]; kp,km=qp.partitions(pp,phi),qm.partitions(pm,phi)
        qs,phis,weights=[],[],[]
        for j,(a,b) in enumerate(zip(kp,km)):
            cuts=_union([a,b],a[-1]); cuts[-1]=a[-1]
            q=(cuts[:-1,None]+np.diff(cuts)[:,None]*x).ravel()
            qs.append(q); phis.append(np.full(len(q),phi[j]))
            weights.append((8*aw[start+j]*np.diff(cuts)[:,None]*w).ravel())
        q,phi,weight=np.concatenate(qs),np.concatenate(phis),np.concatenate(weights); peaks=max(peaks,len(q))
        prep_p,prep_m=qp.prepare_ray(q,phi),qm.prepare_ray(q,phi)
        rho_p=plus.fields(pp,prep_p)['rho']; rho_m=minus.fields(pm,prep_m)['rho']
        up,um=prep_p['u'],prep_m['u']; drho=rho_p-rho_m; rmean=.5*(rho_p+rho_m)
        amean,da=scalar['area']['mean'],scalar['area']['difference']
        d_arho=amean*drho+rmean*da; mean_arho=amean*rmean+.25*da*drho
        du=scalar['parallel']['difference']+np.exp(-plus.section.velocity_shape_exponent*q)*scalar['uc']['difference']
        value=.5*(d_arho*.5*(up**3+um**3)+mean_arho*du*(up*up+up*um+um*um))
        sums.append(math.fsum(weight*value))
        absolute+=math.fsum(weight*.25*(abs(plus.mixing.area*rho_p*up**3)+abs(minus.mixing.area*rho_m*um**3)))
        if callback:
            callback(dict(completed_angles=min(start+batch_size,len(angles)),total_angles=len(angles),peak_batch_nodes=peaks))
    return dict(derivative=math.fsum(sums)/(2*step),step_m=step,order=order,angular_order=angular_order,
        angles=len(angles),peak_batch_nodes=peaks,kinetic_absolute_flux=absolute,
        scalar_pair=scalar,density_binary64=True,actual_values_not_EOS_derivatives=True,
        root_statistics=[plus.section.phase_inverse.statistics,minus.section.phase_inverse.statistics])
