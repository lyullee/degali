"""Necessary local profile constraints, independent of the modal-rate solve."""

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
from degali.addons.enriched_transport import PrescribedRadialMixing,EnrichedModalTransport
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.boundary_compatibility import boundary_compatibility


def inspect(p,par,supplied):
    m=p.mixing
    model=EnrichedModalTransport(p,par,scalar_mixing=supplied,thermal_species_ratio=1.,mechanical_work=WORK)
    t=np.linspace(0.,1.,1025)
    d=model.local(np.ones_like(t),t)
    q=p.q0*(1+t*t)
    chi,oldfm=supplied.evaluate(q).T
    old=m.local(q)
    scales=np.column_stack([np.maximum(abs(old['y']*oldfm),1e-12),
        np.maximum(abs(old['h']/old['rho']*oldfm),1.),np.maximum(abs(old['u']*oldfm),1.)])
    dy=model.area*d['rho']*chi*d['grad_y'][:,0]/(2*p.q0)
    dh=model.area*d['rho']*chi*d['grad_h'][:,0]/(2*p.q0)
    rows=[]
    for j in range(len(t)):
        out=boundary_compatibility(y=d['y'][j],specific_h=d['specific_h'][j],velocity=d['u'][j],
            ambient_speed=m.wind,theta=m.state[3],diffusion_y=dy[j],diffusion_h=dh[j],scales=scales[j])
        row=dict(t=t[j],y=d['y'][j],specific_h=d['specific_h'][j],velocity=d['u'][j],
            temperature=d['temperature'][j],density=d['rho'][j],diffusion_y=dy[j],diffusion_h=dh[j],
            scales=scales[j],physical_y_m=math.sqrt(2*p.q0)*m.sy,physical_n_m=math.sqrt(2*p.q0)*m.sn*t[j],
            **out)
        row['lp_identity_error']=abs(max(abs(out['physical_residual']))-out['physical_minimax'])
        row['sharp_identity_error']=abs(max(abs(out['unconstrained_residual']))-out['unavoidable_minimax'])
        rows.append(row)
    worst=max(rows,key=lambda r:r['physical_minimax'])
    return dict(samples=1025,worst=worst,rows=rows,
        maximum_unconstrained_bound=max(r['unavoidable_minimax'] for r in rows),
        maximum_physical_bound=worst['physical_minimax'],
        samples_proving_incompatible_at_5percent=sum(r['physical_minimax']>.05 for r in rows),
        numerical_identity_error=max(max(r['lp_identity_error'],r['sharp_identity_error']) for r in rows),
        sufficient_for_conservative_flow=False)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite compatibility evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'enriched_segments_a_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed immutable groupA required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/boundary_compatibility.py',ROOT/'tests/test_boundary_compatibility.py',
        ROOT/'docs/prereg-boundary-compatibility.md']
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes=hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    seeds=read('edge_conservative_refit_combined_2026-09-05.json')['rows']
    initial=read('coupled_initialization_combined_2026-09-06.json')['rows']
    boundary=read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json')
    reduced=read('e35_reduced.json')['trials']
    measured={r['trial']:r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    data=dict(phase='postscreen_rate_free_boundary_compatibility',completed=False,rows=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),adopted=False,field_scored=False)
    caches={}
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    for n in (10,11,12,22,23,24,25):
        p,replay=construct_projection(next(r for r in seeds if r['trial']==n),boundary,reduced,measured,caches)
        par=np.array(next(r for r in initial if r['trial']==n)['parameters'])
        supplied=PrescribedRadialMixing.from_weak_baseline(p,order=32)
        row=dict(trial=n,state='initial',source_reconstruction=replay,**inspect(p,par,supplied))
        data['rows'].append(row)
        print('compatibility trial',n,'initial lower-bound',row['maximum_physical_bound'],flush=True)
        checkpoint()
        if n==10:
            driver=EnrichedShortSegment(p,par,scalar_mixing=supplied,thermal_species_ratio=1.,mechanical_work=WORK,
                mixing_update='equilibrium_velocity_width')
            original=next(r for r in upstream['rows'] if r['trial']==10)
            segment=next(s for s in original['segments'] if s['divisions']==8)
            if segment['completed'] or segment['failure_stage']!='endpoint':
                raise ValueError('expected immutable8-step failed endpoint')
            state=np.array(segment['parameters'])
            ds=original['length_m']/8
            first=driver.evaluate(state)
            middle=driver.evaluate(state+.5*ds*first['rates'],enforce=False)
            endpoint=state+ds*middle['rates']
            moved,new,newmix,factor=driver.context(endpoint)
            row=dict(trial=n,state='rejected_8step_endpoint',parameters=endpoint,
                reached_m=segment['reached_m']+ds,mixing_amplitude_ratio=factor,**inspect(moved,new,newmix))
            data['rows'].append(row)
            print('compatibility trial10 rejected endpoint lower-bound',row['maximum_physical_bound'],flush=True)
            checkpoint()
    if hashes()!=initial_hashes:
        raise RuntimeError('compatibility dependencies changed during evaluation')
    data['completed']=True
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
