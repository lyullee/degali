"""Post-screen local boundary-row prototype; retain every failed outcome."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK,independent_rays
from degali.addons.enriched_transport import PrescribedRadialMixing,EnrichedModalTransport
from degali.addons.boundary_moment_transport import BoundaryMomentTransport,lift_same_fields


def rate_fields(model,rates):
    t=np.linspace(0.,1.,17)
    a,b=np.meshgrid(t,t,indexing='ij')
    d=model.local(a.ravel(),b.ravel())
    return np.column_stack([d['dc']@rates/d['c'],d['dh']@rates/d['h'],
        d['drho']@rates/d['rho'],d['du']@rates/np.maximum(abs(d['u']),1.)])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite boundary-row evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'coupled_initialization_combined_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or not upstream['all_passed']:
        raise ValueError('completed unchanged initialization required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen upstream changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/boundary_moment_transport.py',ROOT/'tests/test_boundary_moment_transport.py',
        ROOT/'docs/prereg-boundary-moment-transport.md']
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes=hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    seeds=read('edge_conservative_refit_combined_2026-09-05.json')['rows']
    boundary=read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json')
    reduced=read('e35_reduced.json')['trials']
    measured={r['trial']:r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    results=dict(phase='postscreen_boundary_row_local_screen',completed=False,selected_trials=[10,23],rows=[],
        sha256=initial_hashes,started_utc=datetime.now(timezone.utc).isoformat(),adopted=False,field_scored=False)
    caches={}
    def checkpoint():
        partial.write_text(json.dumps(results,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    for n in results['selected_trials']:
        row=dict(trial=n,levels=[],local_screen_passed=False)
        results['rows'].append(row)
        checkpoint()
        try:
            p,row['source_reconstruction']=construct_projection(next(r for r in seeds if r['trial']==n),boundary,reduced,measured,caches)
            par=np.array(next(r for r in upstream['rows'] if r['trial']==n)['parameters'])
            supplied=PrescribedRadialMixing.from_weak_baseline(p,order=32)
            kwargs=dict(scalar_mixing=supplied,thermal_species_ratio=1.,mechanical_work=WORK)
            field_rates=[]
            for degree in (p.basis.degree,p.basis.degree+2):
                view,parameters=lift_same_fields(p,par,degree)
                level=dict(degree=degree,evaluations=[],passed=False)
                row['levels'].append(level)
                checkpoint()
                for order,angular in ((8,48),(16,96)):
                    model=BoundaryMomentTransport(view,parameters,order=order,angular_order=angular,split_angles=False,**kwargs)
                    try:
                        out=model.assemble()
                    except ValueError as exc:
                        level['linear_failure']=getattr(model,'linear_diagnostics',{})
                        raise ValueError(f'degree{degree} order{order}/{angular}: {exc}') from exc
                    level['evaluations'].append(out)
                    checkpoint()
                    print('boundary rows trial',n,'degree',degree,'order',order,angular,
                        'edges',out['edge_defects'],'chi',out['minimum_chi_momentum'],
                        'numerics',out['local_numerics_passed'],flush=True)
                coarse,fine=level['evaluations']
                physical=EnrichedModalTransport(view,parameters,**kwargs)
                rays=independent_rays(physical,fine['rates'])
                level['independent_rays']=rays
                level['rate_refinement']=float(max(abs(fine['rates']-coarse['rates'])/np.maximum(abs(fine['rates']),1.)))
                level['matrix_refinement']=float(np.max(abs(fine['matrix']-coarse['matrix'])/np.maximum(abs(fine['matrix']),1.)))
                level['passed']=bool(all(o['local_numerics_passed'] for o in (coarse,fine))
                    and max(level['rate_refinement'],level['matrix_refinement'])<=1e-5
                    and max(*fine['edge_defects'].values(),*rays['edge_defects'].values())<=.05
                    and min(fine['minimum_chi_momentum'],rays['minimum_chi_momentum'])>=0.
                    and max(fine['maximum_outward_mass'],rays['maximum_outward_mass'])<0.
                    and fine['curvature_half_width']<.1)
                level['rate_fields']=rate_fields(physical,fine['rates'])
                field_rates.append(level['rate_fields'])
                checkpoint()
                if not level['passed']:
                    break
            if len(row['levels'])==2 and all(l['passed'] for l in row['levels']):
                a,b=field_rates
                row['degree_rate_change']=float(np.max(abs(a-b)/np.maximum(abs(b),1.)))
                a,b=[l['evaluations'][-1]['rates'][:7] for l in row['levels']]
                row['degree_base_change']=float(max(abs(a-b)/np.maximum(abs(b),1.)))
                row['local_screen_passed']=bool(max(row['degree_rate_change'],row['degree_base_change'])<=1e-5)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            row['failure']=str(exc)
        row['completed']=True
        checkpoint()
        print('boundary rows trial',n,'finished',row.get('failure'),row['local_screen_passed'],flush=True)
    if hashes()!=initial_hashes:
        raise RuntimeError('boundary-row dependencies changed during screen')
    results['completed']=True
    args.output.write_text(json.dumps(results,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
