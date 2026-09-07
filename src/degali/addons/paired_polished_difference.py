"""Same common-grid value differences, with explicit same-table root polishing."""

import math
import numpy as np
from .enriched_segments import encode_enriched,shifted_field_view
from .edge_conservative_refit import FaceSplitSquareMoments
from .phase_radial_quadrature import gauss_rule
from .polished_phase_inverse import PolishedPhaseMassEnthalpyInverter


def _union(values,scale):
    out=np.sort(np.concatenate(values))
    out=out[np.r_[True,np.diff(out)>scale*2e-13]]
    return out


def _product_difference(ap,am,bp,bm):
    return (ap-am)*(.5*(bp+bm))+(.5*(ap+am))*(bp-bm)


def paired_polished_moment_difference(projection,parameters,rates,step,*,order=8,angular_order=8,batch_size=8):
    rates=np.asarray(rates,float)
    if rates.shape!=(7+projection.count,) or not np.all(np.isfinite(rates)) or not math.isfinite(step) or step<=0.:
        raise ValueError('finite matching direction and positive difference distance required')
    if not isinstance(batch_size,int) or isinstance(batch_size,bool) or batch_size<1:
        raise ValueError('positive batch size required')
    origin=encode_enriched(projection,parameters)
    plus,pp=shifted_field_view(projection,origin+step*rates)
    minus,pm=shifted_field_view(projection,origin-step*rates)
    plus.section.phase_inverse=PolishedPhaseMassEnthalpyInverter(plus.section.thermodynamics)
    minus.section.phase_inverse=PolishedPhaseMassEnthalpyInverter(minus.section.thermodynamics)
    if plus.q0!=minus.q0:
        raise ValueError('paired normalized quadrature requires identical q0')
    qp,qm=FaceSplitSquareMoments(plus),FaceSplitSquareMoments(minus)
    ak=_union([qp.angular_knots(pp),qm.angular_knots(pm)],1.)
    ak[-1]=math.pi/4
    ax,aw0=gauss_rule(angular_order); x,w=gauss_rule(order)
    angles=(ak[:-1,None]+np.diff(ak)[:,None]*ax).ravel()
    aw=(np.diff(ak)[:,None]*aw0).ravel()
    sums=[]; magnitude=[]; peak=0
    for start in range(0,len(angles),batch_size):
        phi=angles[start:start+batch_size]
        kp,km=qp.partitions(pp,phi),qm.partitions(pm,phi)
        qs,phis,weights=[],[],[]
        for j,(a,b) in enumerate(zip(kp,km)):
            cuts=_union([a,b],a[-1]); cuts[-1]=a[-1]
            q=(cuts[:-1,None]+np.diff(cuts)[:,None]*x).ravel()
            qs.append(q); phis.append(np.full(len(q),phi[j]))
            weights.append((8*aw[start+j]*np.diff(cuts)[:,None]*w).ravel())
        q,phi,weight=np.concatenate(qs),np.concatenate(phis),np.concatenate(weights)
        peak=max(peak,len(q))
        prep_p,prep_m=qp.prepare_ray(q,phi),qm.prepare_ray(q,phi)
        dp,dm=plus.fields(pp,prep_p),minus.fields(pm,prep_m)
        ap,am=plus.mixing.area,minus.mixing.area
        up,um=prep_p['u'],prep_m['u']
        rp,rm=ap*dp['rho'],am*dm['rho']
        cp,cm=ap*dp['c'],am*dm['c']
        hp,hm=ap*dp['h'],am*dm['h']
        mass=_product_difference(rp,rm,up,um)
        species=_product_difference(cp,cm,up,um)
        enthalpy=_product_difference(hp,hm,up,um)
        du=up-um
        u2p,u2m=up*up,um*um
        u3p,u3m=u2p*up,u2m*um
        momentum=(rp-rm)*(.5*(u2p+u2m))+(.5*(rp+rm))*du*(up+um)
        kinetic=(rp-rm)*(.5*(u3p+u3m))+(.5*(rp+rm))*du*(u2p+up*um+u2m)
        theta_p,theta_m=plus.mixing.state[3],minus.mixing.state[3]
        cp0,cm0,sp0,sm0=math.cos(theta_p),math.cos(theta_m),math.sin(theta_p),math.sin(theta_m)
        momentum_mean=.5*(rp*u2p+rm*u2m)
        px=momentum*.5*(cp0+cm0)+momentum_mean*(cp0-cm0)
        pz=momentum*.5*(sp0+sm0)+momentum_mean*(sp0-sm0)
        # The modal tests have fixed normalized coordinates and unchanged basis.
        if not np.array_equal(prep_p['psi'],prep_m['psi']):
            raise ValueError('paired modal test coordinates unexpectedly changed')
        values=np.column_stack([mass,species,px,pz,enthalpy+.5*kinetic,
            species[:,None]*prep_p['psi'],enthalpy[:,None]*prep_p['psi']])
        sums.append(np.array([math.fsum(column) for column in (weight[:,None]*values).T]))
        magnitude.append(math.fsum(weight*(.5*(abs(hp*up)+abs(hm*um))+.25*(abs(rp*u3p)+abs(rm*u3m)))))
    total=np.array([math.fsum(column) for column in np.array(sums).T])
    return dict(derivative=total/(2*step),difference=total,step_m=step,
        angles=len(angles),peak_batch_nodes=peak,energy_flux_absolute_terms=math.fsum(magnitude),
        fields_are_binary64=True,common_union_phase_grid=True,local_difference_before_sum=True,
        roots_polished_on_same_table=True,root_statistics=[plus.section.phase_inverse.statistics,minus.section.phase_inverse.statistics])
