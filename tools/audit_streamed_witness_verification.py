"""Memory-bounded restart with the original fixed candidate and gates."""

import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK
from audit_solenoidal_transport import independent_vectors
from degali.addons.enriched_transport import PrescribedRadialMixing,EnrichedModalTransport,shifted_advective_moments
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.streaming_flux_audit import stream_retained_operator


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite streamed verification')
    ref=ROOT/'reference/preslhy'
    origin=ref/'admissible_flux_independent_2026-09-06.partial.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if upstream['completed'] or len(upstream['evaluations'])!=1:
        raise ValueError('expected preserved memory-stopped audit with one coarse evaluation')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/streaming_flux_audit.py',ROOT/'tests/test_streaming_flux_audit.py',
        ROOT/'docs/prereg-streamed-witness-verification.md']
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
    mix=PrescribedRadialMixing.from_weak_baseline(p,order=32)
    driver=EnrichedShortSegment(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
        mixing_update='equilibrium_velocity_width')
    failed=np.array(next(r for r in read('enriched_boundary_exit_2026-09-06.json')['rows'] if r['divisions']==8)['endpoint'])
    p,par,mix,_=driver.context(failed)
    old=next(r for r in read('admissible_flux_witness_2026-09-06.json')['rows'] if r['state']=='rejected_8step_endpoint')['evaluations'][-1]
    witness={k:np.array(old[k]) for k in ('rates','mass_modes','momentum_modes')}
    model=EnrichedModalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK)
    data=dict(phase='streamed_independent_fixed_witness',completed=False,evaluations=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,passed=False,
        original_memory_stopped=True,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        for order,angular in ((4,8),(8,16)):
            data['stage']=f'stream_{order}_{angular}'; checkpoint()
            def progress(item):
                if item['completed_angles']%128==0 or item['completed_angles']==item['total_angles']:
                    data['progress']=item; checkpoint()
                    print('stream',order,angular,item,flush=True)
            out=stream_retained_operator(model,flux_modes=4,order=order,angular_order=angular,batch_size=8,
                witness=witness,callback=progress)
            data['evaluations'].append(out); checkpoint()
            print('stream result',order,angular,'weak',out['maximum_retained_error'],'edges',out['edge_defects'],flush=True)
        coarse,fine=data['evaluations']
        dense=upstream['evaluations'][0]
        data['dense_coarse_parity']={key:float(np.max(abs(coarse[key]-np.array(dense[key]))/np.maximum(abs(np.array(dense[key])),1.)))
            for key in ('matrix','right','advective_jacobian')}
        data['matrix_refinement']=float(np.max(abs(fine['matrix']-coarse['matrix'])/np.maximum(abs(fine['matrix']),1.)))
        data['right_refinement']=float(max(abs(fine['right']-coarse['right'])/np.maximum(abs(fine['right']),1.)))
        data['stage']='actual_moment_differences'; checkpoint()
        rates=witness['rates']
        t=np.linspace(0.,1.,129); a,b=np.meshgrid(t,t,indexing='ij')
        psi=p.prepare(a.ravel(),b.ravel())['psi']
        steps=[1e-6/max(max(abs(rates)),1.)]
        for value,speed in zip(np.split(par,2),np.split(rates[7:],2)):
            value,change=psi@value,psi@speed
            mask=abs(change)>1e-12
            steps.append(float(np.min(.1*(.1-abs(value[mask]))/abs(change[mask]),initial=np.inf)))
        ds=min(steps)
        if not np.isfinite(ds) or ds<=0.:
            raise ValueError('no valid two-sided finite difference distance')
        expected=fine['advective_jacobian']@rates
        derivatives=[]
        for step in (ds,ds/2):
            plus=shifted_advective_moments(p,par,rates,step,order=8,angular_order=8)
            minus=shifted_advective_moments(p,par,rates,-step,order=8,angular_order=8)
            derivatives.append((plus-minus)/(2*step))
        fd=dict(step_m=ds,expected=expected,derivatives=derivatives,
            errors=[float(max(abs(x-expected)/np.maximum(abs(expected),1.))) for x in derivatives],
            refinement=float(max(abs(derivatives[1]-derivatives[0])/np.maximum(abs(expected),1.))))
        data['finite_difference']=fd
        data['independent_rays']=independent_vectors(model,witness)
        rays=data['independent_rays']
        data['passed']=bool(max(data['dense_coarse_parity'].values())<=1e-8
            and max(o['maximum_retained_error'] for o in data['evaluations'])<=1e-8
            and max(data['matrix_refinement'],data['right_refinement'],*fd['errors'],fd['refinement'])<=1e-5
            and max(fine['weak_heat_scaled_error'],abs(fine['mass_boundary_error'])/max(abs(old['source'][0]),1.))<=1e-8
            and max(*fine['edge_defects'].values(),*rays['edge_defects'].values())<=.05
            and min(fine['minimum_shear_production'],rays['minimum_shear_production'],fine['minimum_radial_momentum_diffusivity'],rays['minimum_radial_momentum_diffusivity'])>=0.
            and max(fine['maximum_outward_mass'],rays['maximum_outward_mass'])<0. and fine['curvature_half_width']<.1)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('streamed witness dependencies changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('streamed witness finished',data['passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
