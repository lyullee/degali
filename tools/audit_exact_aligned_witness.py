"""Conservative local aligned flow on exact geometry, with new-direction FD."""
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
from degali.addons.exact_transverse_geometry import exact_geometry_field_view
from degali.addons.solenoidal_transport import SolenoidalTransport
from degali.addons.aligned_flux_witness import stream_aligned_columns,aligned_witness,independent_aligned_rays
from degali.addons.stable_enthalpy_difference import stable_enthalpy_difference
from degali.addons.paired_exact_kinetic_difference import paired_exact_kinetic_difference


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('output',type=Path); args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite exact aligned witness evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'exact_geometry_energy_2026-09-06.json'; aligned_origin=ref/'aligned_flux_witness_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8')); history=json.loads(aligned_origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or not upstream['energy_difference_passed']:
        raise ValueError('coherent geometry full actual energy gate must pass first')
    for evidence in (upstream,history):
        for name,digest in evidence['sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in set(upstream['sha256'])|set(history['sha256'])]+[origin,aligned_origin,Path(__file__).resolve(),
        ROOT/'tests/test_paired_exact_kinetic_difference.py',ROOT/'docs/prereg-exact-aligned-witness.md']
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
    _,_,mix,_=driver.context(failed)
    p,par,geometry=exact_geometry_field_view(p,failed)
    container=SolenoidalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
        flux_modes=4,order=16,angular_order=96); model=container.model
    data=dict(phase='exact_geometry_conservative_aligned_witness',completed=False,evaluations=[],differences=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,
        geometry=geometry,necessary_screen_passed=False,energy_difference_passed=False,passed=False,
        adopted=False,field_scored=False,turbulence_closure=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        def progress(item):
            if item['completed_angles']%256==0 or item['completed_angles']==item['total_angles']:
                data['progress']=item; checkpoint(); print('exact aligned',data['stage'],item,flush=True)
        for operator in upstream['operators']:
            order,angular=operator['order'],operator['angular_order']
            data['stage']=f'columns_{order}_{angular}'; checkpoint()
            columns=stream_aligned_columns(model,order=order,angular_order=angular,callback=progress)
            conservation=max(float(max(abs(columns[key]))) for key in ('global_heat_correction','mass_face_correction','momentum_face_correction'))
            row=dict(columns=columns,order=order,angular_order=angular,global_correction_absolute_error=conservation)
            data['evaluations'].append(row); checkpoint()
            if conservation>1e-7:
                raise ValueError('unit-amplitude global conservation failed')
            row['witness']=aligned_witness(container,operator,columns); checkpoint()
        coarse,fine=data['evaluations']
        if not all(row['witness']['feasible'] for row in data['evaluations']):
            raise ValueError('aligned flow infeasible on coherent exact geometry')
        c,f=coarse['witness'],fine['witness']
        for key in ('rates','amplitudes'):
            data[key+'_refinement']=float(max(abs(f[key]-c[key])/np.maximum(abs(f[key]),1.)))
        data['columns_refinement']=float(np.max(abs(fine['columns']['columns']-coarse['columns']['columns'])/np.maximum(abs(fine['columns']['columns']),1.)))
        errors=[]
        for operator,row in zip(upstream['operators'],data['evaluations']):
            matrix=np.array(operator['matrix'])[:,:model.count]; right=np.array(operator['right'])
            errors.append(float(max(abs(matrix@f['rates']+row['columns']['columns']@f['amplitudes']-right)/np.maximum(abs(right),1.))))
        data['cross_volume_errors']=errors; data['independent_rays']=independent_aligned_rays(model,f)
        rays=data['independent_rays']
        data['necessary_screen_passed']=bool(max(data['rates_refinement'],data['amplitudes_refinement'],data['columns_refinement'])<=1e-5
            and max(*errors,f['maximum_inequality_violation'])<=1e-8 and max(rays['edge_defects'])<=.05
            and min(rays['minimum_shear_production'],rays['minimum_radial_momentum_diffusivity'])>=0.
            and rays['maximum_outward_mass']<0. and rays['maximum_nonradial_stress_fraction']<=1e-10)
        checkpoint()
        if not data['necessary_screen_passed']:
            raise ValueError('exact aligned necessary local screen failed')
        expected=np.array(upstream['operators'][-1]['advective_jacobian'])@f['rates']
        data['expected_moment_derivative']=expected
        data['source_energy_error']=float(abs(expected[4]-upstream['operators'][-1]['source'][4])/max(abs(upstream['operators'][-1]['source'][4]),1.))
        scale=max(abs(expected[4]),1.); step=upstream['original_step']
        for factor in (.5,1.):
            enthalpy=stable_enthalpy_difference(p,par,f['rates'],step*factor,order=64,precision=70,wind_policy='exact_constraint')
            for order,angular in ((4,8),(8,16)):
                data['stage']=f'actual_energy_{factor}_{order}_{angular}'; checkpoint()
                kinetic=paired_exact_kinetic_difference(p,par,f['rates'],step*factor,order=order,angular_order=angular,callback=progress)
                derivative=enthalpy['derivative']+kinetic['derivative']
                row=dict(factor=factor,enthalpy=enthalpy,kinetic=kinetic,total_derivative=derivative,
                    energy_relative_error=abs(derivative-expected[4])/scale)
                data['differences'].append(row); checkpoint()
                print('exact aligned new energy',factor,order,angular,row['energy_relative_error'],flush=True)
        values=[r['total_derivative'] for r in data['differences']]
        data['all_difference_refinement']=(max(values)-min(values))/scale
        data['energy_difference_passed']=bool(max(data['all_difference_refinement'],*(r['energy_relative_error'] for r in data['differences']))<=1e-5)
        data['passed']=bool(data['necessary_screen_passed'] and data['energy_difference_passed'] and data['source_energy_error']<=1e-8)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('exact aligned dependencies changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('exact aligned completed',data['passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
