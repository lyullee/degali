"""Independent phase-angle weak-equation audit of one FROZEN LP witness."""

import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK
from audit_solenoidal_transport import independent_vectors
from degali.addons.enriched_transport import PrescribedRadialMixing,shifted_advective_moments,panel_cumulative
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.solenoidal_transport import SolenoidalTransport,redistributed_fluxes


def diagnostics(container,witness,operator):
    model,mesh=container.model,container.mesh
    p,m,n=model.projection,model.projection.mixing,model.count
    rates,ma,pa=[witness[k] for k in ('rates','mass_modes','momentum_modes')]
    solution=np.r_[rates,ma,pa]
    matrix,right=operator['matrix'][:n-2],operator['right'][:n-2]
    residual=(matrix@solution-right)/np.maximum(abs(right),1.)
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
    vec=redistributed_fluxes(model,d,f[:,0],f[:,1],rates,ma,pa)
    eve=redistributed_fluxes(model,ed,fe[:,0],fe[:,1],rates,ma,pa)
    em,ep=eve['mass'][:,0],eve['momentum'][:,0]
    u,rho,y,h=[ed[k] for k in ('u','rho','y','specific_h')]
    fh=(.5*m.wind**2+.5*u*u)*em-u*ep
    ec=y*em-model.area*rho*mesh.edge_chi*ed['grad_y'][:,0]/(2*p.q0)
    eh=h*em-model.area*rho*mesh.edge_chi*ed['grad_h'][:,0]/(2*p.q0)-fh
    epx=ep-m.wind*math.cos(m.state[3])*em
    edge=np.column_stack([ec,eh,epx])/mesh.edge_scales
    ew=16*p.q0*mesh.angular_weights/np.cos(mesh.angles)**2
    axial=mesh.weight@d['bh']@rates
    production=mesh.weight@vec['production']
    boundary=ew@fh
    jac=np.vstack([matrix[:5,:n],*(np.einsum('n,nk,nr->kr',mesh.weight,d['psi'],d[key],optimize=True) for key in ('bc','bh'))])
    return dict(matrix=matrix,right=right,retained_scaled_residual=residual,
        maximum_retained_error=float(max(abs(residual))),advective_jacobian=jac,
        edge_defects=dict(zip(('hydrogen','heat','momentum'),np.max(abs(edge),axis=0))),
        minimum_shear_production=float(min(np.min(vec['production']),np.min(eve['production']))),
        minimum_radial_momentum_diffusivity=float(min(np.min(vec['radial_momentum_diffusivity']),np.min(eve['radial_momentum_diffusivity']))),
        maximum_outward_mass=float(np.max(em)),curvature_half_width=float(abs(rates[3])*math.sqrt(2*p.q0)*m.sn),
        mass_boundary_error=float(ew@em+operator['source'][0]),
        weak_heat_scaled_error=float(abs(axial-production+boundary)/max(1.,abs(axial),abs(production),abs(boundary))),
        heat_terms=dict(axial=axial,production=production,boundary=boundary))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite independent witness verification')
    ref=ROOT/'reference/preslhy'
    origin=ref/'admissible_flux_witness_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed immutable witness screen required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),ROOT/'docs/prereg-admissible-flux-independent.md']
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes=hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    seed=next(r for r in read('edge_conservative_refit_combined_2026-09-05.json')['rows'] if r['trial']==10)
    initial=next(r for r in read('coupled_initialization_combined_2026-09-06.json')['rows'] if r['trial']==10)
    measured={r['trial']:r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    p,replay=construct_projection(seed,read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json'),
        read('e35_reduced.json')['trials'],measured,{})
    par=np.array(initial['parameters'])
    supplied=PrescribedRadialMixing.from_weak_baseline(p,order=32)
    driver=EnrichedShortSegment(p,par,scalar_mixing=supplied,thermal_species_ratio=1.,mechanical_work=WORK,
        mixing_update='equilibrium_velocity_width')
    failure=read('enriched_boundary_exit_2026-09-06.json')
    failed=np.array(next(r for r in failure['rows'] if r['divisions']==8)['endpoint'])
    p,par,supplied,_=driver.context(failed)
    old=next(r for r in upstream['rows'] if r['state']=='rejected_8step_endpoint')['evaluations'][-1]
    if not old['feasible']:
        raise ValueError('selected witness was not feasible')
    witness={k:np.array(old[k]) for k in ('rates','mass_modes','momentum_modes')}
    data=dict(phase='independent_phase_split_admissible_witness',completed=False,evaluations=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,passed=False,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        for order,angular in ((4,8),(8,16)):
            data['stage']=f'phase_split_{order}_{angular}'; checkpoint()
            model=SolenoidalTransport(p,par,scalar_mixing=supplied,thermal_species_ratio=1.,mechanical_work=WORK,
                flux_modes=4,order=order,angular_order=angular,split_angles=True)
            operator=model.assemble()
            out=diagnostics(model,witness,operator)
            data['evaluations'].append(out); checkpoint()
            print('independent phase split',order,angular,'weak',out['maximum_retained_error'],'edges',out['edge_defects'],flush=True)
        coarse,fine=data['evaluations']
        data['matrix_refinement']=float(np.max(abs(fine['matrix']-coarse['matrix'])/np.maximum(abs(fine['matrix']),1.)))
        data['right_refinement']=float(max(abs(fine['right']-coarse['right'])/np.maximum(abs(fine['right']),1.)))
        data['stage']='actual_moment_differences'; checkpoint()
        rates=witness['rates']
        t=np.linspace(0.,1.,129); a,b=np.meshgrid(t,t,indexing='ij')
        psi=p.prepare(a.ravel(),b.ravel())['psi']
        steps=[1e-6/max(max(abs(rates)),1.)]
        for v,speed in zip(np.split(par,2),np.split(rates[7:],2)):
            value,change=psi@v,psi@speed
            mask=abs(change)>1e-12
            steps.append(float(np.min(.1*(.1-abs(value[mask]))/abs(change[mask]),initial=np.inf)))
        ds=min(steps)
        if not np.isfinite(ds) or ds<=0.:
            raise ValueError('no finite two-sided verification step inside shape domain')
        expected=fine['advective_jacobian']@rates
        fds=[]
        for step in (ds,ds/2):
            plus=shifted_advective_moments(p,par,rates,step,order=8,angular_order=8)
            minus=shifted_advective_moments(p,par,rates,-step,order=8,angular_order=8)
            fds.append((plus-minus)/(2*step))
        data['finite_difference']=dict(step_m=ds,expected=expected,derivatives=fds,
            errors=[float(max(abs(fd-expected)/np.maximum(abs(expected),1.))) for fd in fds],
            refinement=float(max(abs(fds[-1]-fds[0])/np.maximum(abs(expected),1.))))
        data['independent_rays']=independent_vectors(model.model,witness)
        rays=data['independent_rays']; fd=data['finite_difference']
        data['passed']=bool(max(o['maximum_retained_error'] for o in data['evaluations'])<=1e-8
            and max(data['matrix_refinement'],data['right_refinement'],*fd['errors'],fd['refinement'])<=1e-5
            and max(fine['weak_heat_scaled_error'],abs(fine['mass_boundary_error'])/max(abs(old['source'][0]),1.))<=1e-8
            and max(*fine['edge_defects'].values(),*rays['edge_defects'].values())<=.05
            and min(fine['minimum_shear_production'],rays['minimum_shear_production'],fine['minimum_radial_momentum_diffusivity'],rays['minimum_radial_momentum_diffusivity'])>=0.
            and max(fine['maximum_outward_mass'],rays['maximum_outward_mass'])<0. and fine['curvature_half_width']<.1)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('independent witness inputs changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('independent witness finished',data['passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
