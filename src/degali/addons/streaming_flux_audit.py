"""Bounded-memory retained weak operator; no hard boundary-rate solve."""

import math
import numpy as np
from .phase_radial_quadrature import gauss_rule
from .enriched_transport import panel_cumulative
from .solenoidal_transport import solenoidal_basis,redistributed_fluxes
from .reservoir_thermal import _PhaseForceView
from .energy_crosswind import IndependentEnergyCrosswind


def stream_retained_operator(model,*,flux_modes=4,order=8,angular_order=16,
                             split_angles=True,batch_size=8,witness=None,callback=None):
    if not isinstance(batch_size,int) or isinstance(batch_size,bool) or batch_size<1:
        raise ValueError('positive integer streaming batch size required')
    p,m,n,size=model.projection,model.projection.mixing,model.count,model.size
    k=flux_modes
    if k not in (4,6):
        raise ValueError('declared4 or6 flux modes required')
    if split_angles:
        ak=model.quadrature.angular_knots(model.parameters)
        crossed=model.mixing.knots[(model.mixing.knots>p.q0)&(model.mixing.knots<2*p.q0)]
        ak=np.sort(np.r_[ak,np.arccos(np.sqrt(p.q0/crossed))])
        ak=ak[np.r_[True,np.diff(ak)>2e-13]]
    else:
        ak=np.array([0.,math.pi/4])
    ax,aw0=gauss_rule(angular_order)
    angles=(ak[:-1,None]+np.diff(ak)[:,None]*ax).ravel()
    aws=(np.diff(ak)[:,None]*aw0).ravel()
    x,w=gauss_rule(order)
    jf,jm=np.zeros((5,n)),np.zeros((2*size,n))
    rhs,extra=np.zeros((2*size,n+1)),np.zeros((2*size,2*k))
    force,work=0.,0.
    heat_axial,heat_production,heat_boundary,mass_boundary,p_mean,p_scale=[0.]*6
    max_edge=np.zeros(3)
    min_prod,min_diff,max_mass,max_nonradial=np.inf,np.inf,-np.inf,0.
    peak_nodes=0
    if witness is not None:
        rates,ma,pa=[np.asarray(witness[key],float) for key in ('rates','mass_modes','momentum_modes')]
        if rates.shape!=(n,) or ma.shape!=(k,) or pa.shape!=(k,) or not all(np.all(np.isfinite(a)) for a in (rates,ma,pa)):
            raise ValueError('finite matching fixed witness required')
    for start in range(0,len(angles),batch_size):
        phi,aw=angles[start:start+batch_size],aws[start:start+batch_size]
        cuts=model.quadrature.partitions(model.parameters,phi)
        knots,qs,weights,phis,slices=[],[],[],[],[]
        offset=0
        for angle,weight,c in zip(phi,aw,cuts):
            crossings=model.mixing.knots[(model.mixing.knots>0.)&(model.mixing.knots<c[-1])]
            c=np.sort(np.r_[c,crossings]); c=c[np.r_[True,np.diff(c)>c[-1]*2e-13]]
            q=(c[:-1,None]+np.diff(c)[:,None]*x).ravel()
            knots.append(c); qs.append(q); phis.append(np.full(len(q),angle))
            weights.append((8*weight*np.diff(c)[:,None]*w).ravel())
            slices.append(slice(offset,offset+len(q))); offset+=len(q)
        q,weight,ph=np.concatenate(qs),np.concatenate(weights),np.concatenate(phis)
        peak_nodes=max(peak_nodes,len(q))
        radius=np.sqrt(q/p.q0)
        d=model.local(np.minimum(radius*np.cos(ph),1.),np.minimum(radius*np.sin(ph),1.))
        t=np.tan(phi)
        ed=model.local(np.ones_like(t),t)
        integrand=np.zeros((len(q),2,n+1))
        integrand[:,0,1:]=-d['bm']
        integrand[:,1,0],integrand[:,1,1:]=d['force'],-d['bp']
        f,fe=np.empty_like(integrand),np.empty((len(phi),2,n+1))
        for j,(s,c) in enumerate(zip(slices,knots)):
            cum,end=panel_cumulative(integrand[s],c,order)
            f[s],fe[j]=cum/(2*q[s,None,None]),end/(2*c[-1])
        fm,fp=f[:,0],f[:,1]
        rho,c,h,u=[d[key] for key in ('rho','c','h','u')]
        chi=model.mixing.evaluate(q)[:,0]
        prod=-(2*q*d['uq'])[:,None]*(fp-u[:,None]*fm)
        dot=np.einsum('nik,ni->nk',d['grad_psi'],d['coordinate'])
        dc=np.einsum('nik,ni->nk',d['grad_psi'],d['grad_y'])
        dh=np.einsum('nik,ni->nk',d['grad_psi'],d['grad_h'])
        rhs[:size]+=np.einsum('n,nk,nr->kr',weight*d['y'],dot,fm,optimize=True)
        rhs[size:]+=np.einsum('n,nk,nr->kr',weight*d['specific_h'],dot,fm,optimize=True)
        rhs[size:]+=np.einsum('n,nk,nr->kr',weight,d['psi'],prod,optimize=True)
        rhs[:size,0]-=(weight*model.area*rho*chi/(2*p.q0))@dc
        rhs[size:,0]-=(weight*model.area*rho*chi*model.ratio/(2*p.q0))@dh
        jm[:size]+=np.einsum('n,nk,nr->kr',weight,d['psi'],d['bc'],optimize=True)
        jm[size:]+=np.einsum('n,nk,nr->kr',weight,d['psi'],d['bh'],optimize=True)
        theta=m.state[3]
        px,pz=d['bp']*math.cos(theta),d['bp']*math.sin(theta)
        px[:,3]-=model.area*rho*u*u*math.sin(theta)
        pz[:,3]+=model.area*rho*u*u*math.cos(theta)
        jf+=np.stack([weight@v for v in (d['bm'],d['bc'],px,pz,d['bh']+d['bk'])])
        force+=float(weight@(9.81*model.area*(p.section.rhoa-rho)))
        work+=float(weight@(d['force']*u))
        ue=ed['u']; kinetic=.5*(m.wind**2+ue**2)
        fh=kinetic[:,None]*fe[:,0]-ue[:,None]*fe[:,1]
        ew=16*p.q0*aw/np.cos(phi)**2
        rhs[size:]-=np.einsum('n,nk,nr->kr',ew,ed['psi'],fh,optimize=True)
        v=solenoidal_basis(d['a'],d['b'],k)
        ev=solenoidal_basis(np.ones_like(t),t,k)[:,0]
        vd=np.einsum('nim,nik->nmk',d['grad_psi'],v,optimize=True)
        ug=(2*p.q0*d['uq'])[:,None]*d['coordinate']
        uv=np.einsum('ni,nik->nk',ug,v)
        extra[:size,:k]+=np.einsum('n,nmk->mk',weight*d['y'],vd,optimize=True)
        extra[size:,:k]+=np.einsum('n,nmk->mk',weight*d['specific_h'],vd,optimize=True)
        extra[size:,:k]+=np.einsum('n,nm,nk->mk',weight*u,d['psi'],uv,optimize=True)
        extra[size:,:k]-=np.einsum('n,nm,nk->mk',ew*kinetic,ed['psi'],ev,optimize=True)
        extra[size:,k:]-=np.einsum('n,nm,nk->mk',weight,d['psi'],uv,optimize=True)
        extra[size:,k:]+=np.einsum('n,nm,nk->mk',ew*ue,ed['psi'],ev,optimize=True)
        if witness is not None:
            vec=redistributed_fluxes(model,d,fm,fp,rates,ma,pa)
            edge=redistributed_fluxes(model,ed,fe[:,0],fe[:,1],rates,ma,pa)
            em,ep=edge['mass'][:,0],edge['momentum'][:,0]
            heat_axial+=float(weight@d['bh']@rates)
            heat_production+=float(weight@vec['production'])
            heat_boundary+=float(ew@(kinetic*em-ue*ep))
            mass_boundary+=float(ew@em)
            ce,oldfm=model.mixing.evaluate(p.q0*(1+t*t)).T
            old=m.local(p.q0*(1+t*t))
            edge_scales=np.column_stack([np.maximum(abs(old['y']*oldfm),1e-12),
                np.maximum(abs(old['h']/old['rho']*oldfm),1.),np.maximum(abs(old['u']*oldfm),1.)])
            ec=ed['y']*em-model.area*ed['rho']*ce*ed['grad_y'][:,0]/(2*p.q0)
            eh=ed['specific_h']*em-model.area*ed['rho']*ce*model.ratio*ed['grad_h'][:,0]/(2*p.q0)-kinetic*em+ue*ep
            epx=ep-m.wind*math.cos(theta)*em
            max_edge=np.maximum(max_edge,np.max(abs(np.column_stack([ec,eh,epx]))/edge_scales,axis=0))
            p_mean+=float((aw/np.cos(phi)**2)@epx)
            p_scale+=float((aw/np.cos(phi)**2)@edge_scales[:,2])
            for value in (vec,edge):
                min_prod=min(min_prod,float(min(value['production'])))
                min_diff=min(min_diff,float(min(value['radial_momentum_diffusivity'])))
                max_nonradial=max(max_nonradial,float(max(value['nonradial_stress_fraction'])))
            max_mass=max(max_mass,float(max(em)))
        if callback:
            callback(dict(completed_angles=min(start+batch_size,len(angles)),total_angles=len(angles),peak_batch_nodes=peak_nodes))
    source=IndependentEnergyCrosswind.source_terms(_PhaseForceView(m,force),m.state)
    source[4]+=work
    matrix=np.zeros((5+2*size,n+2*k))
    matrix[:,:n]=np.vstack([jf,jm-rhs[:,1:]])
    matrix[5:,n:]=-extra
    right=np.r_[source,rhs[:,0]]
    out=dict(matrix=matrix,right=right,source=source,advective_jacobian=np.vstack([jf,jm]),
        force=force,work=work,angles=len(angles),peak_batch_nodes=peak_nodes,batch_size=batch_size,
        order=order,angular_order=angular_order,phase_split_angles=split_angles)
    if witness is not None:
        residual=(matrix@np.r_[rates,ma,pa]-right)/np.maximum(abs(right),1.)
        out.update(retained_scaled_residual=residual,maximum_retained_error=float(max(abs(residual))),
            edge_defects=dict(zip(('hydrogen','heat','momentum'),max_edge)),
            minimum_shear_production=min_prod,minimum_radial_momentum_diffusivity=min_diff,
            maximum_outward_mass=max_mass,maximum_nonradial_stress_fraction=max_nonradial,
            mass_boundary_error=mass_boundary+source[0],momentum_mean_scaled_error=abs(p_mean)/max(p_scale,1.),
            weak_heat_scaled_error=abs(heat_axial-heat_production+heat_boundary)/max(1.,abs(heat_axial),abs(heat_production),abs(heat_boundary)),
            heat_terms=dict(axial=heat_axial,production=heat_production,boundary=heat_boundary),
            curvature_half_width=float(abs(rates[3])*math.sqrt(2*p.q0)*m.sn))
    return out
