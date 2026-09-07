"""Conditional conservative transverse redistribution, not a fitted eddy law.

Retains all global and scalar weak equations. New flux amplitudes may imply
non-radial momentum stress, explicitly beyond scalar isotropic viscosity.
"""

import math
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import legval,legint,legvander
from .enriched_transport import EnrichedModalTransport,panel_cumulative
from .coupled_shape_initialization import FixedMeshCoupledTransport


@lru_cache(maxsize=16)
def stream_polynomial(k):
    if isinstance(k,bool) or int(k)!=k or k<1:
        raise ValueError('positive integer solenoidal mode required')
    coefficients=np.zeros(2*k+1)
    coefficients[-1]=math.sqrt(4*k+1)
    primitive=legint(coefficients)
    primitive[0]-=legval(0.,primitive)
    return coefficients,primitive


def solenoidal_basis(a,b,count):
    """V=(psi_b,-psi_a), psi=a*phi(b)-b*phi(a), phi'=normalized P_2k.

    V is exactly divergence-free; face normal is P_2k(t), with zero mean.
    """
    a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float))
    if a.ndim!=1 or not np.all(np.isfinite(a+b)) or np.any(abs(a)>1.) or np.any(abs(b)>1.):
        raise ValueError('finite vector coordinates inside the square required')
    if isinstance(count,bool) or int(count)!=count or count<1:
        raise ValueError('positive integer solenoidal basis count required')
    vectors=[]
    for k in range(1,count+1):
        c,primitive=stream_polynomial(k)
        pa,pb=legval(a,primitive),legval(b,primitive)
        da,db=legval(a,c),legval(b,c)
        vectors.append(np.column_stack([a*db-pa,b*da-pb]))
    return np.stack(vectors,axis=2)


def redistributed_fluxes(model,d,fm,fp,rates,mass_modes,momentum_modes):
    """Actual normalized vector fluxes; no pointwise clipping or refitting."""
    at=np.r_[1.,rates]
    v=solenoidal_basis(d['a'],d['b'],len(mass_modes))
    mass=d['coordinate']*(fm@at)[:,None]+np.einsum('nik,k->ni',v,mass_modes)
    momentum=d['coordinate']*(fp@at)[:,None]+np.einsum('nik,k->ni',v,momentum_modes)
    stress=momentum-d['u'][:,None]*mass
    grad_u=(2*model.projection.q0*d['uq'])[:,None]*d['coordinate']
    production=-np.sum(grad_u*stress,axis=1)
    norm2=np.sum(d['coordinate']**2,axis=1)
    radial_dot=np.sum(d['coordinate']*stress,axis=1)
    radial_diffusivity=-radial_dot/(model.area*d['rho']*d['uq']*norm2)
    perpendicular=stress-d['coordinate']*(radial_dot/norm2)[:,None]
    transverse_fraction=np.linalg.norm(perpendicular,axis=1)/np.maximum(np.linalg.norm(stress,axis=1),1e-30)
    return dict(mass=mass,momentum=momentum,stress=stress,production=production,
        radial_momentum_diffusivity=radial_diffusivity,nonradial_stress_fraction=transverse_fraction)


class SolenoidalTransport:
    def __init__(self,projection,parameters,*,scalar_mixing,thermal_species_ratio,
                 mechanical_work,flux_modes=4,order=8,angular_order=48,split_angles=False):
        if flux_modes not in (4,6):
            raise ValueError('select the declared4 or6 solenoidal-mode screen')
        self.model=EnrichedModalTransport(projection,parameters,scalar_mixing=scalar_mixing,
            thermal_species_ratio=thermal_species_ratio,mechanical_work=mechanical_work)
        self.mesh=FixedMeshCoupledTransport(projection,parameters,scalar_mixing=scalar_mixing,
            thermal_species_ratio=thermal_species_ratio,mechanical_work=mechanical_work,
            order=order,angular_order=angular_order,face_samples=65,split_angles=split_angles)
        self.flux_modes=flux_modes

    def assemble(self):
        model,mesh,k=self.model,self.mesh,self.flux_modes
        p,m,size,n=model.projection,model.projection.mixing,model.size,model.count
        base=mesh.evaluate(model.parameters)
        d=model.local(mesh.a,mesh.b)
        ed=model.local(np.ones_like(mesh.face_t),mesh.face_t)
        v=solenoidal_basis(mesh.a,mesh.b,k)
        ev=solenoidal_basis(np.ones_like(mesh.face_t),mesh.face_t,k)[:,0,:]
        q,u=d['q'],d['u']
        integrand=np.zeros((len(q),2,n+1))
        integrand[:,0,1:]=-d['bm']
        integrand[:,1,0],integrand[:,1,1:]=d['force'],-d['bp']
        f,fe=np.empty_like(integrand),np.empty((len(mesh.angles),2,n+1))
        for j,(s,cuts) in enumerate(zip(mesh.slices,mesh.knots)):
            cum,end=panel_cumulative(integrand[s],cuts,mesh.order)
            f[s],fe[j]=cum/(2*q[s,None,None]),end/(2*cuts[-1])
        dot=np.einsum('nim,nik->nmk',d['grad_psi'],v,optimize=True)
        gradient=(2*p.q0*d['uq'])[:,None]*d['coordinate']
        u_dot_v=np.einsum('ni,nik->nk',gradient,v)
        weights=mesh.weight
        ew=16*p.q0*mesh.angular_weights/np.cos(mesh.angles)**2
        ue=ed['u']
        kinetic=.5*(m.wind**2+ue**2)
        rc=np.einsum('n,nmk->mk',weights*d['y'],dot,optimize=True)
        rhm=np.einsum('n,nmk->mk',weights*d['specific_h'],dot,optimize=True)
        rhm+=np.einsum('n,nm,nk->mk',weights*u,d['psi'],u_dot_v,optimize=True)
        rhm-=np.einsum('n,nm,nk->mk',ew*kinetic,ed['psi'],ev,optimize=True)
        rhp=-np.einsum('n,nm,nk->mk',weights,d['psi'],u_dot_v,optimize=True)
        rhp+=np.einsum('n,nm,nk->mk',ew*ue,ed['psi'],ev,optimize=True)
        count=n+2*k
        matrix=np.zeros((5+2*size,count))
        matrix[:,:n]=base['matrix']
        matrix[5:5+size,n:n+k]=-rc
        matrix[5+size:,n:n+k],matrix[5+size:,n+k:]=-rhm,-rhp
        # Boundary vectors are affine in [old state rates, new flux amplitudes].
        em,ep=np.zeros((len(mesh.angles),count+1)),np.zeros((len(mesh.angles),count+1))
        em[:,:n+1],ep[:,:n+1]=fe[:,0],fe[:,1]
        em[:,n+1:n+k+1],ep[:,n+k+1:]=ev,ev
        fh=kinetic[:,None]*em-ue[:,None]*ep
        ec=ed['y'][:,None]*em
        eh=ed['specific_h'][:,None]*em-fh
        ec[:,0]-=model.area*ed['rho']*mesh.edge_chi*ed['grad_y'][:,0]/(2*p.q0)
        eh[:,0]-=model.area*ed['rho']*mesh.edge_chi*model.ratio*ed['grad_h'][:,0]/(2*p.q0)
        epx=ep-m.wind*math.cos(m.state[3])*em
        scalar_count=k//2+1
        indices=np.arange(max(scalar_count,k-1))
        leg=legvander(mesh.face_t,2*indices[-1])[:,2*indices]*np.sqrt(4*indices+1)
        face_weight=mesh.angular_weights/np.cos(mesh.angles)**2
        projected=[np.einsum('n,nj,nr->jr',face_weight,leg,e,optimize=True) for e in (ec,eh,epx)]
        bc=np.vstack([projected[0][:scalar_count],projected[1][:scalar_count],projected[2][1:k-1]])
        matrix=np.vstack([matrix,bc[:,1:]])
        right=np.r_[base['right'],-bc[:,0]]
        solution=np.zeros(count)
        solution[5:7]=math.cos(m.state[3]),math.sin(m.state[3])
        active=np.r_[model.active,np.arange(n,count)]
        scales=np.maximum(np.max(abs(matrix[:,active]),axis=1),1e-12)
        scaled=matrix[:,active]/scales[:,None]
        left,singular,_=np.linalg.svd(scaled,full_matrices=False)
        self.linear_diagnostics=dict(matrix=matrix,right=right,scales=scales,singular_values=singular,
            left_smallest_vector=left[:,-1],active=active)
        if singular[-1]<=1e-12*singular[0]:
            raise ValueError('conservative redistribution lost rank; no pseudoinverse fallback')
        solution[active]=np.linalg.solve(scaled,(right-matrix@solution)/scales)
        rates,ma,pa=solution[:n],solution[n:n+k],solution[n+k:]
        at=np.r_[1.,solution]
        vectors=redistributed_fluxes(model,d,f[:,0],f[:,1],rates,ma,pa)
        # Edge stress uses the complete normal+tangential vector, not just fP_n.
        edge_vectors=redistributed_fluxes(model,ed,fe[:,0],fe[:,1],rates,ma,pa)
        edge=np.column_stack([ec@at,eh@at,epx@at])/mesh.edge_scales
        heat_axial=weights@d['bh']@rates
        heat_production=weights@vectors['production']
        heat_boundary=ew@fh@at
        linear_error=float(max(abs(matrix@solution-right)/np.maximum(abs(right),1.)))
        heat_error=float(abs(heat_axial-heat_production+heat_boundary)/max(1.,abs(heat_axial),abs(heat_production),abs(heat_boundary)))
        mass_error=float(ew@em@at+base['source'][0])
        mean_p=float(face_weight@epx@at)/max(float(face_weight@mesh.edge_scales[:,2]),1.)
        minimum_diff=float(min(np.min(vectors['radial_momentum_diffusivity']),np.min(edge_vectors['radial_momentum_diffusivity'])))
        minimum_production=float(min(np.min(vectors['production']),np.min(edge_vectors['production'])))
        maximum_nonradial=float(max(np.max(vectors['nonradial_stress_fraction']),np.max(edge_vectors['nonradial_stress_fraction'])))
        defects=dict(zip(('hydrogen','heat','momentum'),np.max(abs(edge),axis=0)))
        return dict(rates=rates,mass_modes=ma,momentum_modes=pa,solution=solution,matrix=matrix,right=right,
            source=base['source'],original_rates=base['rates'],original_edge_defects=base['edge_defects'],
            linear_scaled_error=linear_error,weak_heat_scaled_error=heat_error,mass_boundary_error=mass_error,
            momentum_mean_scaled_error=abs(mean_p),heat_terms=dict(axial=heat_axial,production=heat_production,boundary=heat_boundary),
            matrix_condition=float(singular[0]/singular[-1]),edge_defects=defects,edge_residual=edge,
            minimum_radial_momentum_diffusivity=minimum_diff,minimum_shear_production=minimum_production,
            maximum_nonradial_stress_fraction=maximum_nonradial,maximum_outward_mass=float(np.max(em@at)),
            curvature_half_width=float(abs(rates[3])*math.sqrt(2*p.q0)*m.sn),
            local_numerics_passed=bool(max(linear_error,heat_error,abs(mean_p),abs(mass_error)/max(abs(base['source'][0]),1.))<=1e-8),
            all_original_weak_equations_retained=True,scalar_isotropic_momentum_closure=False,
            reynolds_tensor_realizability_established=False,adopted=False)
