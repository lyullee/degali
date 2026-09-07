"""Conservative admissibility witnesses, not a tolerance-tuned closure law."""

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
from degali.addons.enriched_transport import PrescribedRadialMixing
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.solenoidal_transport import SolenoidalTransport
from degali.addons.admissible_flux_witness import admissible_witness


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite admissible-witness evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'solenoidal_transport_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed immutable solenoidal pilot required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/admissible_flux_witness.py',ROOT/'tests/test_admissible_flux_witness.py',
        ROOT/'docs/prereg-admissible-flux-witness.md']
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
    data=dict(phase='conservative_admissible_flux_witnesses',completed=False,rows=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,
        turbulence_closure=False,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    for label,encoded in (('initial',driver.initial),('rejected_8step_endpoint',failed)):
        row=dict(trial=10,state=label,evaluations=[],verified_witness=False)
        data['rows'].append(row); checkpoint()
        try:
            current,parameters,mix,factor=driver.context(encoded)
            old=next(r for r in upstream['rows'] if r['state']==label)['levels'][0]['evaluations']
            for index,(order,angular) in enumerate(((8,48),(16,96))):
                model=SolenoidalTransport(current,parameters,scalar_mixing=mix,thermal_species_ratio=1.,
                    mechanical_work=WORK,flux_modes=4,order=order,angular_order=angular)
                saved={key:np.array(value) if key in ('matrix','right','original_rates','source') else value for key,value in old[index].items()}
                out=admissible_witness(model,saved)
                row['evaluations'].append(out); checkpoint()
                print('admissibility',label,order,angular,'feasible',out['feasible'],
                    'amplitude',out.get('maximum_normalized_amplitude'),'edge',out.get('edge_defects'),flush=True)
            coarse,fine=row['evaluations']
            if fine['feasible']:
                rays=independent_vectors(model.model,fine)
                row['independent']=rays
                row['verified_witness']=bool(fine['local_numerics_passed'] and max(rays['edge_defects'].values())<=.05
                    and min(rays['minimum_shear_production'],fine['minimum_shear_production'])>=0.
                    and min(rays['minimum_radial_momentum_diffusivity'],fine['minimum_radial_momentum_diffusivity'])>=0.
                    and max(rays['maximum_outward_mass'],fine['maximum_outward_mass'])<0.
                    and fine['curvature_half_width']<.1)
            if coarse['feasible'] and fine['feasible']:
                for key in ('rates','amplitudes'):
                    row[key+'_refinement']=float(max(abs(fine[key]-coarse[key])/np.maximum(abs(fine[key]),1.)))
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            row['failure']=str(exc)
        row['completed']=True; checkpoint()
        print('admissibility',label,'verified',row['verified_witness'],'failure',row.get('failure'),flush=True)
    if hashes()!=initial_hashes:
        raise RuntimeError('admissibility dependencies changed during screen')
    data['completed']=True
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
