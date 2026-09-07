"""Feasibility witnesses, deliberately NOT a downstream constitutive law."""

import math
import numpy as np
from scipy.optimize import linprog
from .enriched_transport import panel_cumulative
from .solenoidal_transport import solenoidal_basis,redistributed_fluxes


def minimum_amplitude_witness(matrix,right,scales):
    """Minimize max(abs(alpha/scales)) subject to matrix@alpha<=right."""
    a,b,s=np.asarray(matrix,float),np.asarray(right,float),np.asarray(scales,float)
    if a.ndim!=2 or b.shape!=(a.shape[0],) or s.shape!=(a.shape[1],) or np.any(s<=0.) or not all(np.all(np.isfinite(x)) for x in (a,b,s)):
        raise ValueError('finite matching inequalities and positive amplitude scales required')
    n=len(s)
    aub=np.vstack([np.c_[a*s,np.zeros(len(b))],np.c_[np.eye(n),-np.ones(n)],np.c_[-np.eye(n),-np.ones(n)]])
    bub=np.r_[b,np.zeros(2*n)]
    out=linprog(np.r_[np.zeros(n),1.],A_ub=aub,b_ub=bub,bounds=[(None,None)]*n+[(0.,None)],method='highs')
    result=dict(feasible=bool(out.success),status=int(out.status),message=out.message)
    if out.success:
        alpha=out.x[:n]*s
        result.update(amplitudes=alpha,maximum_normalized_amplitude=float(out.fun),
            maximum_inequality_violation=float(np.max(a@alpha-b,initial=0.)))
    return result


def admissible_witness(container,assembled,*,edge_target=.04):
    if edge_target!=.04:
        raise ValueError('this registered diagnostic fixes the target at4 percent')
    model,mesh,k=container.model,container.mesh,container.flux_modes
    p,m,n=model.projection,model.projection.mixing,model.count
    a=assembled['matrix'][:n-2,:n]
    extra=assembled['matrix'][:n-2,n:]
    right=assembled['right'][:n-2]
    r0=assembled['original_rates']
    scales=np.maximum(np.max(abs(a[:,model.active]),axis=1),1e-12)
    response=np.zeros((n,2*k))
    response[model.active]=np.linalg.solve(a[:,model.active]/scales[:,None],-extra/scales[:,None])
    d=model.local(mesh.a,mesh.b)
    ed=model.local(np.ones_like(mesh.face_t),mesh.face_t)
    q=d['q']
    integrand=np.zeros((len(q),2,n+1))
    integrand[:,0,1:]=-d['bm']
    integrand[:,1,0],integrand[:,1,1:]=d['force'],-d['bp']
    f,fe=np.empty_like(integrand),np.empty((len(mesh.angles),2,n+1))
    for j,(s,cuts) in enumerate(zip(mesh.slices,mesh.knots)):
        cum,end=panel_cumulative(integrand[s],cuts,mesh.order)
        f[s],fe[j]=cum/(2*q[s,None,None]),end/(2*cuts[-1])
    def affine(data,fm,fp):
        v=solenoidal_basis(data['a'],data['b'],k)
        at0=np.r_[1.,r0]; atr=np.vstack([np.zeros((1,2*k)),response])
        mass0=data['coordinate']*(fm@at0)[:,None]
        mom0=data['coordinate']*(fp@at0)[:,None]
        massj=data['coordinate'][:,:,None]*(fm@atr)[:,None,:]
        momj=data['coordinate'][:,:,None]*(fp@atr)[:,None,:]
        massj[:,:,:k]+=v; momj[:,:,k:]+=v
        grad=(2*p.q0*data['uq'])[:,None]*data['coordinate']
        prod0=-np.sum(grad*(mom0-data['u'][:,None]*mass0),axis=1)
        prodj=-np.einsum('ni,nik->nk',grad,momj-data['u'][:,None,None]*massj)
        chi=model.mixing.evaluate(data['q'])[:,0]
        expected=2*data['q']*model.area*data['rho']*chi*data['uq']**2
        prodscale=np.maximum.reduce([abs(prod0),expected,np.full_like(prod0,1e-12)])
        return mass0,mom0,massj,momj,prod0,prodj,prodscale
    volume=affine(d,f[:,0],f[:,1])
    edge=affine(ed,fe[:,0],fe[:,1])
    em,ep,jm,jp=[x[:,0] for x in edge[:4]]
    u,rho,y,h=[ed[key] for key in ('u','rho','y','specific_h')]
    kinetic=.5*(m.wind**2+u*u)
    ec=y*em-model.area*rho*mesh.edge_chi*ed['grad_y'][:,0]/(2*p.q0)
    eh=(h-kinetic)*em+u*ep-model.area*rho*mesh.edge_chi*ed['grad_h'][:,0]/(2*p.q0)
    epx=ep-m.wind*math.cos(m.state[3])*em
    residual=np.column_stack([ec,eh,epx])/mesh.edge_scales
    jac=np.stack([y[:,None]*jm,(h-kinetic)[:,None]*jm+u[:,None]*jp,jp-m.wind*math.cos(m.state[3])*jm],axis=1)/mesh.edge_scales[:,:,None]
    constraints=[jac.reshape(-1,2*k),-jac.reshape(-1,2*k)]
    bounds=[(edge_target-residual).ravel(),(edge_target+residual).ravel()]
    for data in (volume,edge):
        prod0,prodj,scale=data[4:]
        constraints.append(-prodj/scale[:,None]); bounds.append(prod0/scale)
    old_mass=np.maximum(abs(model.mixing.evaluate(ed['q'])[:,1]),1e-12)
    constraints.append(jm/old_mass[:,None]); bounds.append(-1e-6-em/old_mass)
    halfwidth=math.sqrt(2*p.q0)*m.sn
    constraints.extend([response[3:4]*halfwidth,-response[3:4]*halfwidth])
    bounds.extend([np.array([.1-r0[3]*halfwidth]),np.array([.1+r0[3]*halfwidth])])
    matrix,limit=np.vstack(constraints),np.concatenate(bounds)
    massscale=max(abs(assembled['source'][0])/(16*p.q0),1e-12)
    momscale=massscale*max(m.state[4]+m.wind*math.cos(m.state[3]),m.wind,1.)
    amplitude_scales=np.r_[np.full(k,massscale),np.full(k,momscale)]
    result=minimum_amplitude_witness(matrix,limit,amplitude_scales)
    result.update(constraints=len(limit),edge_target=edge_target,amplitude_scales=amplitude_scales,
        feasible_flow_is_a_closure=False,adopted=False)
    if not result['feasible']:
        return result
    alpha=result['amplitudes']; rates=r0+response@alpha
    ma,pa=alpha[:k],alpha[k:]
    vec=redistributed_fluxes(model,d,f[:,0],f[:,1],rates,ma,pa)
    eve=redistributed_fluxes(model,ed,fe[:,0],fe[:,1],rates,ma,pa)
    actual=residual+np.einsum('njk,k->nj',jac,alpha)
    em,ep=eve['mass'][:,0],eve['momentum'][:,0]
    ew=16*p.q0*mesh.angular_weights/np.cos(mesh.angles)**2
    axial=mesh.weight@d['bh']@rates
    production=mesh.weight@vec['production']
    boundary=ew@(kinetic*em-u*ep)
    linear=float(max(abs(a@rates+extra@alpha-right)/np.maximum(abs(right),1.)))
    heat=float(abs(axial-production+boundary)/max(1.,abs(axial),abs(production),abs(boundary)))
    mass=float(ew@em+assembled['source'][0])
    result.update(rates=rates,mass_modes=ma,momentum_modes=pa,source=assembled['source'],
        linear_scaled_error=linear,weak_heat_scaled_error=heat,mass_boundary_error=mass,
        heat_terms=dict(axial=axial,production=production,boundary=boundary),
        edge_defects=dict(zip(('hydrogen','heat','momentum'),np.max(abs(actual),axis=0))),
        minimum_shear_production=float(min(np.min(vec['production']),np.min(eve['production']))),
        minimum_radial_momentum_diffusivity=float(min(np.min(vec['radial_momentum_diffusivity']),np.min(eve['radial_momentum_diffusivity']))),
        maximum_nonradial_stress_fraction=float(max(np.max(vec['nonradial_stress_fraction']),np.max(eve['nonradial_stress_fraction']))),
        maximum_outward_mass=float(np.max(em)),curvature_half_width=float(abs(rates[3])*halfwidth),
        local_numerics_passed=bool(max(linear,heat,abs(mass)/max(abs(assembled['source'][0]),1.),result['maximum_inequality_violation'])<=1e-8))
    return result
