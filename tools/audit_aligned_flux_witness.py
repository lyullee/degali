"""Fixed failed trial10: can a direction-preserving circulation be admissible?"""

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
from degali.addons.enriched_transport import PrescribedRadialMixing
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.solenoidal_transport import SolenoidalTransport
from degali.addons.aligned_flux_witness import stream_aligned_columns,aligned_witness,independent_aligned_rays


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('output',type=Path); args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite aligned witness evidence')
    ref=ROOT/'reference/preslhy'; origin=ref/'streamed_witness_verification_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or max(upstream['dense_coarse_parity'].values())>1e-8:
        raise ValueError('completed matched phase volume required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/aligned_momentum_flux.py',ROOT/'tests/test_aligned_momentum_flux.py',
        ROOT/'src/degali/addons/aligned_flux_witness.py',ROOT/'tests/test_aligned_flux_witness.py',
        ROOT/'docs/prereg-aligned-flux-witness.md']
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
    container=SolenoidalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
        flux_modes=4,order=16,angular_order=96)
    model=container.model
    data=dict(phase='direction_preserving_circulation_witness',completed=False,evaluations=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,
        necessary_screen_passed=False,energy_finite_difference_blocker_remains=True,
        turbulence_closure=False,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        for old in upstream['evaluations']:
            order,angular=old['order'],old['angular_order']
            data['stage']=f'aligned_columns_{order}_{angular}'; checkpoint()
            def progress(item):
                if item['completed_angles']%256==0 or item['completed_angles']==item['total_angles']:
                    data['progress']=item; checkpoint(); print('aligned stream',order,angular,item,flush=True)
            columns=stream_aligned_columns(model,order=order,angular_order=angular,callback=progress)
            row=dict(columns=columns,order=order,angular_order=angular)
            data['evaluations'].append(row); checkpoint()
            conservation=max(float(max(abs(columns[key]))) for key in
                ('global_heat_correction','mass_face_correction','momentum_face_correction'))
            row['global_correction_absolute_error']=conservation
            if conservation>1e-7:
                raise ValueError('aligned unit-amplitude global conservation failed')
            row['witness']=aligned_witness(container,old,columns); checkpoint()
            print('aligned witness',order,angular,row['witness']['feasible'],row['witness'].get('edge_defects'),flush=True)
        coarse,fine=data['evaluations']
        data['columns_refinement']=float(np.max(abs(fine['columns']['columns']-coarse['columns']['columns'])/np.maximum(abs(fine['columns']['columns']),1.)))
        if all(row['witness']['feasible'] for row in data['evaluations']):
            f=fine['witness']; c=coarse['witness']
            for key in ('rates','amplitudes'):
                data[key+'_refinement']=float(max(abs(f[key]-c[key])/np.maximum(abs(f[key]),1.)))
            errors=[]
            for old,row in zip(upstream['evaluations'],data['evaluations']):
                matrix=np.array(old['matrix'])[:,:model.count]; right=np.array(old['right'])
                errors.append(float(max(abs(matrix@f['rates']+row['columns']['columns']@f['amplitudes']-right)/np.maximum(abs(right),1.))))
            data['cross_volume_errors']=errors
            data['independent_rays']=independent_aligned_rays(model,f)
            rays=data['independent_rays']
            data['necessary_screen_passed']=bool(max(data['columns_refinement'],data['rates_refinement'],data['amplitudes_refinement'])<=1e-5
                and max(*errors,f['maximum_inequality_violation'])<=1e-8 and max(rays['edge_defects'])<=.05
                and min(rays['minimum_shear_production'],rays['minimum_radial_momentum_diffusivity'])>=0.
                and rays['maximum_outward_mass']<0. and rays['maximum_nonradial_stress_fraction']<=1e-10)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('aligned witness dependencies changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('aligned witness completed',data['necessary_screen_passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
