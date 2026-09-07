"""Coefficient-free necessary TKE from the fixed aligned trial10 shear field."""
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
from degali.addons.enriched_transport import PrescribedRadialMixing,panel_cumulative
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.solenoidal_transport import SolenoidalTransport
from degali.addons.aligned_flux_witness import affine_vectors
from degali.addons.stress_realizability import shear_tke_lower_bound


def evaluate(container,witness):
    model,mesh=container.model,container.mesh; p,m,n=model.projection,model.projection.mixing,model.count
    d=model.local(mesh.a,mesh.b); q=d['q']; f=np.empty((len(q),2,n+1))
    integrand=np.zeros_like(f); integrand[:,0,1:]=-d['bm']
    integrand[:,1,0]=d['force']; integrand[:,1,1:]=-d['bp']
    for s,cuts in zip(mesh.slices,mesh.knots):
        cum,_=panel_cumulative(integrand[s],cuts,mesh.order)
        f[s]=cum/(2*q[s,None,None])
    v=affine_vectors(model,d,f[:,0],f[:,1],np.array(witness['rates']),np.zeros((n,4)))
    alpha=np.array(witness['amplitudes'])
    mass=v['mass']+np.einsum('nik,k->ni',v['mass_j'],alpha)
    momentum=v['momentum']+np.einsum('nik,k->ni',v['momentum_j'],alpha)
    stress=momentum-d['u'][:,None]*mass
    lengths=math.sqrt(2*p.q0)*np.array([m.sy,m.sn])
    bound=shear_tke_lower_bound(stress,model.area,d['rho'],lengths)
    k=bound['minimum_tke']; rms=bound['minimum_total_fluctuation_rms']/d['u']
    minimum_flux=float(model.area*mesh.weight@(d['rho']*d['u']*k))
    mean_ke=float(.5*model.area*mesh.weight@(d['rho']*d['u']**3))
    idx=int(np.argmax(k)); jdx=int(np.argmax(rms))
    return dict(order=mesh.order,angular_order=len(mesh.angles),nodes=len(q),
        minimum_tke_axial_flux=minimum_flux,mean_axial_kinetic_flux=mean_ke,
        minimum_tke_flux_over_mean_kinetic=minimum_flux/mean_ke,
        maximum_pointwise_minimum_tke=float(k[idx]),maximum_total_fluctuation_rms_over_axial_speed=float(rms[jdx]),
        max_k_location=dict(a=float(d['a'][idx]),b=float(d['b'][idx]),local_axial_speed=float(d['u'][idx]),
            physical_transverse_lengths=lengths,physical_shear_covariance=bound['physical_shear_covariance'][idx],
            attaining_covariance=bound['attaining_covariance'][idx]),
        max_rms_location=dict(a=float(d['a'][jdx]),b=float(d['b'][jdx]),local_axial_speed=float(d['u'][jdx])),
        actual_k_or_dissipation_established=False,axial_derivative_lower_bound=False)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('output',type=Path); args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite TKE lower-bound evidence')
    ref=ROOT/'reference/preslhy'; origin=ref/'aligned_flux_witness_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or not upstream['necessary_screen_passed']:
        raise ValueError('completed aligned local witness required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/stress_realizability.py',ROOT/'tests/test_stress_realizability.py',
        ROOT/'docs/prereg-shear-tke-lower-bound.md']
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
    par=np.array(initial['parameters']); mix=PrescribedRadialMixing.from_weak_baseline(p,order=32)
    driver=EnrichedShortSegment(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
        mixing_update='equilibrium_velocity_width')
    failed=np.array(next(r for r in read('enriched_boundary_exit_2026-09-06.json')['rows'] if r['divisions']==8)['endpoint'])
    p,par,mix,_=driver.context(failed)
    witness=upstream['evaluations'][-1]['witness']
    data=dict(phase='necessary_shear_tke_storage',completed=False,evaluations=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,resolved=False,
        adopted=False,field_scored=False,turbulence_closure=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        for order,angular in ((8,48),(16,96)):
            container=SolenoidalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
                flux_modes=4,order=order,angular_order=angular)
            out=evaluate(container,witness); data['evaluations'].append(out); checkpoint()
            print('TKE lower bound',order,angular,out['minimum_tke_flux_over_mean_kinetic'],flush=True)
        a,b=data['evaluations']
        data['integrated_lower_flux_refinement']=abs(a['minimum_tke_axial_flux']-b['minimum_tke_axial_flux'])/max(abs(b['minimum_tke_axial_flux']),1.)
        data['resolved']=bool(data['integrated_lower_flux_refinement']<=1e-3)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('TKE lower-bound dependencies changed')
    data['completed']=True
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('TKE lower bound finished',data['resolved'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
