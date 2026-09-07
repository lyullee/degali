"""Independent common-grid actual-field differences for a frozen witness."""

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
from degali.addons.paired_polished_difference import paired_polished_moment_difference


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite paired energy evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'paired_energy_difference_2026-09-06.json'
    history=json.loads(origin.read_text(encoding='utf-8'))
    if not history['completed']:
        raise ValueError('completed unpolished paired diagnosis required')
    upstream=json.loads((ref/'phase_consistent_flux_witness_2026-09-06.json').read_text(encoding='utf-8'))
    if not upstream['completed'] or 'independent_volume' not in upstream:
        raise ValueError('completed phase-consistent candidate with independent volume required')
    for name,digest in history['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in history['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/paired_polished_difference.py',ROOT/'tests/test_paired_polished_difference.py',
        ROOT/'src/degali/addons/polished_phase_inverse.py',ROOT/'tests/test_polished_phase_inverse.py',
        ROOT/'docs/prereg-polished-energy-difference.md']
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
    p,par,_,_=driver.context(failed)
    rates=np.array(upstream['evaluations'][-1]['rates'])
    expected=np.array(upstream['finite_difference']['expected'])
    h=upstream['finite_difference']['step_m']
    data=dict(phase='paired_polished_actual_energy_difference',completed=False,rows=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,
        expected=expected,original_step=h,same_step_pair_passed=False,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        for factor in (.5,1.,2.,4.):
            data['stage']=f'paired_factor_{factor}'; checkpoint()
            out=paired_polished_moment_difference(p,par,rates,h*factor,order=8,angular_order=8,batch_size=8)
            error=abs(out['derivative']-expected)/np.maximum(abs(expected),1.)
            out.update(factor=factor,relative_errors=error,maximum_error=float(max(error)),
                energy_error=float(error[4]),nonenergy_maximum_error=float(max(np.delete(error,4))))
            data['rows'].append(out); checkpoint()
            print('polished paired factor',factor,'energy',out['energy_error'],'other',out['nonenergy_maximum_error'],flush=True)
        half,full=data['rows'][:2]
        data['original_pair_refinement']=float(max(abs(half['derivative']-full['derivative'])/np.maximum(abs(expected),1.)))
        data['same_step_pair_passed']=bool(max(half['maximum_error'],full['maximum_error'],data['original_pair_refinement'])<=1e-5)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    final_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if final_hashes!=initial_hashes:
        raise RuntimeError('paired energy dependencies changed')
    data['completed']=True; data['stage']='completed'
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('polished paired same-step result',data['same_step_pair_passed'],data.get('failure'),flush=True)


if __name__=='__main__':
    main()
