"""Direction-preserving circulation feasibility, NOT an evolution closure."""

import math
import numpy as np
from .aligned_momentum_flux import aligned_mode_vectors_analytic
from .admissible_flux_witness import minimum_amplitude_witness
from .phase_radial_quadrature import gauss_rule
from .enriched_transport import panel_cumulative


def modes(model, data, count=4):
    p,m=model.projection,model.projection.mixing
    exponent=p.section.velocity_shape_exponent
    return aligned_mode_vectors_analytic(data['a'],data['b'],q0=p.q0,
        ambient_parallel=m.wind*math.cos(m.state[3]),excess_velocity=m.state[4],
        velocity_exponent=exponent,flux_modes=count)


def stream_aligned_columns(model,*,count=4,order=8,angular_order=16,batch_size=8,callback=None):
    """Integrate new columns only on the same independently split phase volume."""
    if count!=4 or batch_size<1:
        raise ValueError('registered four-mode pilot and positive batch size required')
    p,size=model.projection,model.size
    ak=model.quadrature.angular_knots(model.parameters)
    crossed=model.mixing.knots[(model.mixing.knots>p.q0)&(model.mixing.knots<2*p.q0)]
    ak=np.sort(np.r_[ak,np.arccos(np.sqrt(p.q0/crossed))])
    ak=ak[np.r_[True,np.diff(ak)>2e-13]]
    ax,aw0=gauss_rule(angular_order); x,w=gauss_rule(order)
    angles=(ak[:-1,None]+np.diff(ak)[:,None]*ax).ravel()
    aws=(np.diff(ak)[:,None]*aw0).ravel()
    extra=np.zeros((5+2*size,count))
    global_heat=np.zeros(count); mass_face=np.zeros(count); momentum_face=np.zeros(count)
    peak=0
    for start in range(0,len(angles),batch_size):
        phi,aw=angles[start:start+batch_size],aws[start:start+batch_size]
        qs,weights,phis=[],[],[]
        for angle,weight,c in zip(phi,aw,model.quadrature.partitions(model.parameters,phi)):
            crosses=model.mixing.knots[(model.mixing.knots>0.)&(model.mixing.knots<c[-1])]
            c=np.sort(np.r_[c,crosses]); c=c[np.r_[True,np.diff(c)>c[-1]*2e-13]]
            q=(c[:-1,None]+np.diff(c)[:,None]*x).ravel()
            qs.append(q); phis.append(np.full(len(q),angle))
            weights.append((8*weight*np.diff(c)[:,None]*w).ravel())
        q,ph,weight=np.concatenate(qs),np.concatenate(phis),np.concatenate(weights)
        peak=max(peak,len(q)); radius=np.sqrt(q/p.q0)
        d=model.local(np.minimum(radius*np.cos(ph),1.),np.minimum(radius*np.sin(ph),1.))
        t=np.tan(phi); ed=model.local(np.ones_like(t),t)
        v,e=modes(model,d,count),modes(model,ed,count)
        dot=np.einsum('nim,nik->nmk',d['grad_psi'],v['mass_modes'],optimize=True)
        production=-(2*q*d['uq'])[:,None]*v['radial_stress_modes']
        kinetic=.5*(model.projection.mixing.wind**2+ed['u']**2)
        boundary=kinetic[:,None]*e['mass_modes'][:,0]-ed['u'][:,None]*e['momentum_modes'][:,0]
        ew=16*p.q0*aw/np.cos(phi)**2
        extra[5:5+size]-=np.einsum('n,nmk->mk',weight*d['y'],dot,optimize=True)
        extra[5+size:]-=np.einsum('n,nmk->mk',weight*d['specific_h'],dot,optimize=True)
        extra[5+size:]-=np.einsum('n,nm,nk->mk',weight,d['psi'],production,optimize=True)
        extra[5+size:]+=np.einsum('n,nm,nk->mk',ew,ed['psi'],boundary,optimize=True)
        global_heat+=weight@production-ew@boundary
        mass_face+=ew@e['mass_modes'][:,0]
        momentum_face+=ew@e['momentum_modes'][:,0]
        if callback:
            callback(dict(completed_angles=min(start+batch_size,len(angles)),total_angles=len(angles),peak_batch_nodes=peak))
    return dict(columns=extra,global_heat_correction=global_heat,mass_face_correction=mass_face,
        momentum_face_correction=momentum_face,angles=len(angles),peak_batch_nodes=peak,
        order=order,angular_order=angular_order)


def affine_vectors(model,data,fm,fp,rates,response):
    k=response.shape[1]
    v=modes(model,data,k)
    at=np.r_[1.,rates]; jr=np.vstack([np.zeros((1,k)),response])
    coord=data['coordinate']
    mass=coord*(fm@at)[:,None]; momentum=coord*(fp@at)[:,None]
    mass_j=coord[:,:,None]*(fm@jr)[:,None,:]+v['mass_modes']
    momentum_j=coord[:,:,None]*(fp@jr)[:,None,:]+v['momentum_modes']
    grad=(2*model.projection.q0*data['uq'])[:,None]*coord
    prod=-np.sum(grad*(momentum-data['u'][:,None]*mass),axis=1)
    prod_j=-np.einsum('ni,nik->nk',grad,momentum_j-data['u'][:,None,None]*mass_j)
    return dict(mass=mass,momentum=momentum,mass_j=mass_j,momentum_j=momentum_j,
        production=prod,production_j=prod_j)


def edge_affine(model,data,v,scales):
    p,m=model.projection,model.projection.mixing
    em,ep,jm,jp=[v[key][:,0] for key in ('mass','momentum','mass_j','momentum_j')]
    y,h,u,rho=[data[key] for key in ('y','specific_h','u','rho')]
    kinetic=.5*(m.wind**2+u*u)
    chi=model.mixing.evaluate(data['q'])[:,0]
    residual=np.column_stack([y*em-model.area*rho*chi*data['grad_y'][:,0]/(2*p.q0),
        (h-kinetic)*em+u*ep-model.area*rho*chi*model.ratio*data['grad_h'][:,0]/(2*p.q0),
        ep-m.wind*math.cos(m.state[3])*em])/scales
    jac=np.stack([y[:,None]*jm,(h-kinetic)[:,None]*jm+u[:,None]*jp,
        jp-m.wind*math.cos(m.state[3])*jm],axis=1)/scales[:,:,None]
    return residual,jac


def aligned_witness(container,operator,columns):
    model,mesh=container.model,container.mesh
    p,m,n,k=model.projection,model.projection.mixing,model.count,4
    matrix=np.asarray(operator['matrix'])[:,:n]; right=np.asarray(operator['right'])
    extra=np.asarray(columns['columns'])
    scales=np.maximum(np.max(abs(matrix[:,model.active]),axis=1),1e-12)
    scaled=matrix[:,model.active]/scales[:,None]
    singular=np.linalg.svd(scaled,compute_uv=False)
    if singular[-1]<=1e-12*singular[0]:
        raise ValueError('retained aligned state matrix lost rank')
    r0=np.zeros(n); r0[5:7]=math.cos(m.state[3]),math.sin(m.state[3])
    r0[model.active]=np.linalg.solve(scaled,(right-matrix@r0)/scales)
    response=np.zeros((n,k)); response[model.active]=np.linalg.solve(scaled,-extra/scales[:,None])
    d=model.local(mesh.a,mesh.b); ed=model.local(np.ones_like(mesh.face_t),mesh.face_t)
    q=d['q']; integrand=np.zeros((len(q),2,n+1))
    integrand[:,0,1:]=-d['bm']; integrand[:,1,0]=d['force']; integrand[:,1,1:]=-d['bp']
    f,fe=np.empty_like(integrand),np.empty((len(mesh.angles),2,n+1))
    for j,(s,cuts) in enumerate(zip(mesh.slices,mesh.knots)):
        cum,end=panel_cumulative(integrand[s],cuts,mesh.order)
        f[s],fe[j]=cum/(2*q[s,None,None]),end/(2*cuts[-1])
    volume=affine_vectors(model,d,f[:,0],f[:,1],r0,response)
    edge=affine_vectors(model,ed,fe[:,0],fe[:,1],r0,response)
    residual,jac=edge_affine(model,ed,edge,mesh.edge_scales)
    constraints=[jac.reshape(-1,k),-jac.reshape(-1,k)]
    bounds=[(.04-residual).ravel(),(.04+residual).ravel()]
    for local,v in ((d,volume),(ed,edge)):
        chi=model.mixing.evaluate(local['q'])[:,0]
        expected=2*local['q']*model.area*local['rho']*chi*local['uq']**2
        ps=np.maximum.reduce([abs(v['production']),expected,np.full_like(expected,1e-12)])
        constraints.append(-v['production_j']/ps[:,None]); bounds.append(v['production']/ps)
    old_mass=np.maximum(abs(model.mixing.evaluate(ed['q'])[:,1]),1e-12)
    constraints.append(edge['mass_j'][:,0]/old_mass[:,None])
    bounds.append(-1e-6-edge['mass'][:,0]/old_mass)
    halfwidth=math.sqrt(2*p.q0)*m.sn
    constraints.extend([response[3:4]*halfwidth,-response[3:4]*halfwidth])
    bounds.extend([np.array([.1-r0[3]*halfwidth]),np.array([.1+r0[3]*halfwidth])])
    amplitude_scales=np.full(k,max(abs(operator['source'][0])/(16*p.q0),1e-12))
    out=minimum_amplitude_witness(np.vstack(constraints),np.concatenate(bounds),amplitude_scales)
    out.update(edge_target=.04,amplitude_scales=amplitude_scales,baseline_rates=r0,response=response,
        baseline_edge_defects=np.max(abs(residual),axis=0),stress_direction_preserved=True,
        amplitudes_are_turbulence_closure=False,adopted=False,field_scored=False)
    if out['feasible']:
        rates=r0+response@out['amplitudes']
        out.update(rates=rates,edge_defects=np.max(abs(residual+np.einsum('njk,k->nj',jac,out['amplitudes'])),axis=0),
            retained_error=float(max(abs(matrix@rates+extra@out['amplitudes']-right)/np.maximum(abs(right),1.))))
    return out


def independent_aligned_rays(model,witness):
    p,m=model.projection,model.projection.mixing
    angles=np.linspace(0.,math.pi/4,65)
    defect=np.zeros(3); minimum_prod,minimum_diff,max_mass,nonradial=np.inf,np.inf,-np.inf,0.
    signatures=[]
    for angle,cuts in zip(angles,model.quadrature.partitions(model.parameters,angles)):
        crosses=model.mixing.knots[(model.mixing.knots>0.)&(model.mixing.knots<cuts[-1])]
        cuts=np.sort(np.r_[cuts,crosses]); cuts=cuts[np.r_[True,np.diff(cuts)>cuts[-1]*2e-13]]
        d,fm,fp,fe,_=model._ray(angle,cuts,16)
        t=math.tan(angle); ed=model.local(np.array([1.]),np.array([t]))
        volume=affine_vectors(model,d,fm,fp,witness['rates'],np.zeros((model.count,4)))
        edge=affine_vectors(model,ed,fe[None,0],fe[None,1],witness['rates'],np.zeros((model.count,4)))
        oldfm=model.mixing.evaluate(ed['q'])[:,1]; old=m.local(ed['q'])
        scales=np.column_stack([np.maximum(abs(old['y']*oldfm),1e-12),
            np.maximum(abs(old['h']/old['rho']*oldfm),1.),np.maximum(abs(old['u']*oldfm),1.)])
        er,ej=edge_affine(model,ed,edge,scales)
        defect=np.maximum(defect,np.max(abs(er+np.einsum('njk,k->nj',ej,witness['amplitudes'])),axis=0))
        for local,v in ((d,volume),(ed,edge)):
            mass=v['mass']+np.einsum('nik,k->ni',v['mass_j'],witness['amplitudes'])
            momentum=v['momentum']+np.einsum('nik,k->ni',v['momentum_j'],witness['amplitudes'])
            stress=momentum-local['u'][:,None]*mass; coord=local['coordinate']
            norm2=np.sum(coord**2,axis=1); radial=np.sum(coord*stress,axis=1)/norm2
            fraction=np.linalg.norm(stress-coord*radial[:,None],axis=1)/np.maximum(np.linalg.norm(stress,axis=1),1e-30)
            production=v['production']+v['production_j']@witness['amplitudes']
            chi=-radial/(model.area*local['rho']*local['uq'])
            minimum_prod=min(minimum_prod,float(min(production))); minimum_diff=min(minimum_diff,float(min(chi)))
            nonradial=max(nonradial,float(max(fraction)))
            if local is ed:
                max_mass=max(max_mass,float(mass[0,0]))
            else:
                signatures.append(np.column_stack([mass,momentum]))
    return dict(edge_defects=defect,minimum_shear_production=minimum_prod,
        minimum_radial_momentum_diffusivity=minimum_diff,maximum_outward_mass=max_mass,
        maximum_nonradial_stress_fraction=nonradial,signature=np.vstack(signatures),angles=65,order=16)
