"""Coherent exact geometry and actual full-energy directional differences."""
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
from degali.addons.enriched_transport import PrescribedRadialMixing,EnrichedModalTransport
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.exact_transverse_geometry import exact_geometry_field_view,exact_geometry_tangent
from degali.addons.streaming_flux_audit import stream_retained_operator
from degali.addons.stable_enthalpy_difference import stable_enthalpy_difference
from degali.addons.paired_exact_kinetic_difference import paired_exact_kinetic_difference


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('output',type=Path); args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite coherent geometry evidence')
    ref=ROOT/'reference/preslhy'; origin=ref/'smooth_enthalpy_conditioning_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed smooth component diagnosis required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/exact_transverse_geometry.py',ROOT/'tests/test_exact_transverse_geometry.py',
        ROOT/'src/degali/addons/paired_exact_kinetic_difference.py',ROOT/'docs/prereg-exact-geometry-energy.md']
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
    old,par,mix,_=driver.context(failed)
    p,par,geometry=exact_geometry_field_view(p,failed)
    model=EnrichedModalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK)
    witness=read('phase_consistent_flux_witness_2026-09-06.json')
    rates=np.array(witness['evaluations'][-1]['rates']); step=witness['finite_difference']['step_m']
    data=dict(phase='coherent_exact_geometry_energy',completed=False,operators=[],differences=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,
        geometry=geometry,original_step=step,old_expected=witness['finite_difference']['expected'],
        geometry_value_change=geometry['values']-np.array([old.mixing.sy,old.mixing.sn,old.mixing.wind]),
        energy_difference_passed=False,old_rates_resolved_for_new_geometry=False,
        adopted=False,field_scored=False,turbulence_closure=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        other=exact_geometry_tangent(p.section.jetplume,p.mixing.state,step=1e-28)
        data['geometry_tangent_step_change']=float(np.max(abs(other['jacobian']-geometry['jacobian'])/np.maximum(abs(geometry['jacobian']),1.)))
        def progress(item):
            if item['completed_angles']%256==0 or item['completed_angles']==item['total_angles']:
                data['progress']=item; checkpoint(); print('exact geometry',data['stage'],item,flush=True)
        for order,angular in ((4,8),(8,16)):
            data['stage']=f'operator_{order}_{angular}'; checkpoint()
            out=stream_retained_operator(model,order=order,angular_order=angular,callback=progress)
            out['directional_moment_derivative']=out['advective_jacobian']@rates
            data['operators'].append(out); checkpoint()
        expected=data['operators'][-1]['directional_moment_derivative'][4]; scale=max(abs(expected),1.)
        data['new_expected_energy']=expected
        data['old_direction_new_source_energy_residual']=expected-data['operators'][-1]['source'][4]
        data['expected_energy_refinement']=abs(data['operators'][0]['directional_moment_derivative'][4]-expected)/scale
        for factor in (.5,1.):
            enthalpy=stable_enthalpy_difference(p,par,rates,step*factor,order=64,precision=70,wind_policy='exact_constraint')
            for order,angular in ((4,8),(8,16)):
                data['stage']=f'actual_difference_{factor}_{order}_{angular}'; checkpoint()
                kinetic=paired_exact_kinetic_difference(p,par,rates,step*factor,order=order,angular_order=angular,callback=progress)
                derivative=enthalpy['derivative']+kinetic['derivative']
                row=dict(factor=factor,enthalpy=enthalpy,kinetic=kinetic,total_derivative=derivative,
                    energy_relative_error=abs(derivative-expected)/scale)
                data['differences'].append(row); checkpoint()
                print('coherent E',factor,order,angular,'error',row['energy_relative_error'],flush=True)
        values=[r['total_derivative'] for r in data['differences']]
        data['all_difference_refinement']=(max(values)-min(values))/scale
        data['energy_difference_passed']=bool(data['geometry_tangent_step_change']<=1e-10
            and max(data['expected_energy_refinement'],data['all_difference_refinement'],
                *(r['energy_relative_error'] for r in data['differences']))<=1e-5)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('exact geometry dependencies changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('coherent exact geometry energy finished',data['energy_difference_passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
