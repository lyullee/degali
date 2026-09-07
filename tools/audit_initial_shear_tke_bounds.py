"""Necessary shear/TKE bounds at all seven frozen accepted initial sections."""
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
from audit_shear_tke_lower_bound import evaluate
from degali.addons.enriched_transport import PrescribedRadialMixing
from degali.addons.solenoidal_transport import SolenoidalTransport


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('output',type=Path); args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite initial TKE evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'shear_tke_lower_bound_2026-09-06.json'; initial_origin=ref/'coupled_initialization_combined_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8')); initials=json.loads(initial_origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or not upstream['resolved'] or len(initials['rows'])!=7 or not all(r['passed'] for r in initials['rows']):
        raise ValueError('resolved bound pilot and seven accepted initial fields required')
    for evidence in (upstream,initials):
        for name,digest in evidence['sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in set(upstream['sha256'])|set(initials['sha256'])]+[origin,initial_origin,Path(__file__).resolve(),
        ROOT/'docs/prereg-initial-shear-tke-bounds.md']
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes=hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    seeds={r['trial']:r for r in read('edge_conservative_refit_combined_2026-09-05.json')['rows']}
    measured={r['trial']:r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    boundary=read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json'); reduced=read('e35_reduced.json')['trials']
    data=dict(phase='initial_seven_shear_tke_bounds',completed=False,rows=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),adopted=False,field_scored=False,turbulence_closure=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint(); caches={}
    for initial in initials['rows']:
        trial=initial['trial']; row=dict(trial=trial,resolved=False,evaluations=[]); data['rows'].append(row); checkpoint()
        try:
            p,row['source_reconstruction']=construct_projection(seeds[trial],boundary,reduced,measured,caches)
            parameters=np.array(initial['parameters']); rates=np.array(initial['evaluations'][-1]['rates'])
            if parameters.shape!=(p.count,) or rates.shape!=(p.count+7,):
                raise ValueError('reconstructed initial basis differs from frozen dimensions')
            mixing=PrescribedRadialMixing.from_weak_baseline(p,order=32)
            witness=dict(rates=rates,amplitudes=np.zeros(4))
            for order,angular in ((8,48),(16,96)):
                row['stage']=f'bound_{order}_{angular}'; checkpoint()
                container=SolenoidalTransport(p,parameters,scalar_mixing=mixing,thermal_species_ratio=1.,mechanical_work=WORK,
                    flux_modes=4,order=order,angular_order=angular)
                out=evaluate(container,witness)
                out['total_sampled_rays']=out.pop('angular_order'); out['requested_angular_order']=angular
                row['evaluations'].append(out); checkpoint()
                print('initial TKE',trial,order,angular,out['minimum_tke_flux_over_mean_kinetic'],flush=True)
            a,b=row['evaluations']
            row['integrated_lower_flux_refinement']=abs(a['minimum_tke_axial_flux']-b['minimum_tke_axial_flux'])/max(abs(b['minimum_tke_axial_flux']),1.)
            row['resolved']=bool(row['integrated_lower_flux_refinement']<=1e-3)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            row['failure']=str(exc)
        row['stage']='completed'; checkpoint()
    if hashes()!=initial_hashes:
        raise RuntimeError('initial TKE dependencies changed')
    data['completed']=True; data['resolved_trials']=[r['trial'] for r in data['rows'] if r['resolved']]
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('initial TKE bounds completed',data['resolved_trials'],flush=True)


if __name__=='__main__':
    main()
