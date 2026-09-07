"""New local candidates using the independently resolved volume equations."""

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
from degali.addons.enriched_transport import PrescribedRadialMixing,EnrichedModalTransport,shifted_advective_moments
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.solenoidal_transport import SolenoidalTransport
from degali.addons.admissible_flux_witness import admissible_witness
from degali.addons.streaming_flux_audit import stream_retained_operator


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite phase-consistent witness evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'streamed_witness_verification_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if (not upstream['completed'] or len(upstream['evaluations'])!=2
            or max(upstream['dense_coarse_parity'].values())>1e-8
            or max(upstream['matrix_refinement'],upstream['right_refinement'])>1e-5):
        raise ValueError('completed matched/refined independent volume operator required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),ROOT/'docs/prereg-phase-consistent-flux-witness.md']
    initial_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
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
    model=EnrichedModalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK)
    container=SolenoidalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
        flux_modes=4,order=16,angular_order=96)
    data=dict(phase='phase_consistent_local_flux_witness',completed=False,evaluations=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,passed=False,
        turbulence_closure=False,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        operators=[]
        for previous in upstream['evaluations']:
            saved={k:np.array(previous[k]) for k in ('matrix','right','source')}
            matrix,right=saved['matrix'][:,:model.count],saved['right']
            rates=np.zeros(model.count)
            rates[5:7]=math.cos(p.mixing.state[3]),math.sin(p.mixing.state[3])
            scales=np.maximum(np.max(abs(matrix[:,model.active]),axis=1),1e-12)
            scaled=matrix[:,model.active]/scales[:,None]
            singular=np.linalg.svd(scaled,compute_uv=False)
            if singular[-1]<=1e-12*singular[0]:
                raise ValueError('phase-resolved baseline state matrix lost rank')
            rates[model.active]=np.linalg.solve(scaled,(right-matrix@rates)/scales)
            saved['original_rates']=rates
            operators.append(saved)
            data['stage']=f"witness_on_volume_{previous['order']}_{previous['angular_order']}"; checkpoint()
            out=admissible_witness(container,saved)
            out['volume_order'],out['volume_angular_order']=previous['order'],previous['angular_order']
            data['evaluations'].append(out); checkpoint()
            print('phase-consistent',out['volume_order'],out['volume_angular_order'],'feasible',out['feasible'],
                'amplitude',out.get('maximum_normalized_amplitude'),'edges',out.get('edge_defects'),flush=True)
        coarse,fine=data['evaluations']
        if not all(x['feasible'] for x in (coarse,fine)):
            raise ValueError('both phase volume inputs must give feasible candidates')
        for key in ('rates','amplitudes'):
            data[key+'_refinement']=float(max(abs(fine[key]-coarse[key])/np.maximum(abs(fine[key]),1.)))
        solution=np.r_[fine['rates'],fine['mass_modes'],fine['momentum_modes']]
        data['cross_volume_errors']=[float(max(abs(o['matrix']@solution-o['right'])/np.maximum(abs(o['right']),1.))) for o in operators]
        data['stage']='fresh_streamed_fine_balance'; checkpoint()
        def progress(item):
            if item['completed_angles']%256==0 or item['completed_angles']==item['total_angles']:
                data['progress']=item; checkpoint()
                print('new-witness stream',item,flush=True)
        fixed={k:fine[k] for k in ('rates','mass_modes','momentum_modes')}
        resolved=stream_retained_operator(model,flux_modes=4,order=8,angular_order=16,batch_size=8,
            witness=fixed,callback=progress)
        data['independent_volume']=resolved
        data['independent_rays']=independent_vectors(model,fixed)
        data['stage']='actual_moment_differences'; checkpoint()
        rates=fine['rates']
        t=np.linspace(0.,1.,129); a,b=np.meshgrid(t,t,indexing='ij')
        psi=p.prepare(a.ravel(),b.ravel())['psi']
        steps=[1e-6/max(max(abs(rates)),1.)]
        for value,speed in zip(np.split(par,2),np.split(rates[7:],2)):
            value,change=psi@value,psi@speed
            mask=abs(change)>1e-12
            steps.append(float(np.min(.1*(.1-abs(value[mask]))/abs(change[mask]),initial=np.inf)))
        ds=min(steps)
        if not np.isfinite(ds) or ds<=0.:
            raise ValueError('no positive finite verification step')
        expected=resolved['advective_jacobian']@rates
        derivatives=[]
        for step in (ds,ds/2):
            plus=shifted_advective_moments(p,par,rates,step,order=8,angular_order=8)
            minus=shifted_advective_moments(p,par,rates,-step,order=8,angular_order=8)
            derivatives.append((plus-minus)/(2*step))
        fd=dict(step_m=ds,expected=expected,derivatives=derivatives,
            errors=[float(max(abs(x-expected)/np.maximum(abs(expected),1.))) for x in derivatives],
            refinement=float(max(abs(derivatives[1]-derivatives[0])/np.maximum(abs(expected),1.))))
        data['finite_difference']=fd
        rays=data['independent_rays']
        data['passed']=bool(max(data['rates_refinement'],data['amplitudes_refinement'],*fd['errors'],fd['refinement'])<=1e-5
            and max(*data['cross_volume_errors'],resolved['maximum_retained_error'],fine['maximum_inequality_violation'])<=1e-8
            and max(resolved['weak_heat_scaled_error'],abs(resolved['mass_boundary_error'])/max(abs(resolved['source'][0]),1.))<=1e-8
            and max(*resolved['edge_defects'].values(),*rays['edge_defects'].values())<=.05
            and min(resolved['minimum_shear_production'],rays['minimum_shear_production'],resolved['minimum_radial_momentum_diffusivity'],rays['minimum_radial_momentum_diffusivity'])>=0.
            and max(resolved['maximum_outward_mass'],rays['maximum_outward_mass'])<0. and resolved['curvature_half_width']<.1)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    final_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if final_hashes!=initial_hashes:
        raise RuntimeError('phase-consistent witness dependencies changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('phase-consistent witness finished',data['passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
