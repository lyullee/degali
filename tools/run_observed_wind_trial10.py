"""Two prescribed wind-input cells; same six source fluxes and scoring keys."""
import argparse
from datetime import datetime,timezone
import math
from pathlib import Path

import numpy as np

from degali.addons.yawed_crosswind import YawedCrosswind,YawedTrajectory
from degali.addons.observed_wind_profile import with_observed_wind
from stage2_matched_source import ROOT,REF,read,digest,write_new
from stage1 import verify
from run_preslhy_ambient_profile_audit import replay
from run_bounded_liquid_handoff import original_target_check
from run_yawed_crosswind_trial10 import score


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--mode',choices=['magnitude_only','measured_vector'],required=True)
    parser.add_argument('--step',type=float,choices=[.01,.02],default=.02)
    args=parser.parse_args()
    directory=args.directory.resolve()
    if directory.exists(): raise FileExistsError(directory)
    verify()
    paths=[REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
        REF/'gaussian_enthalpy_profile_allflux_interface_2026-09-05.json',
        REF/'lateral_temperature_symmetry_2026-09-06.json',REF/'temperature_mean_audit_2026-09-06.json',
        REF/'stage1_single_image_2026-09-06.json',REF/'e35_reduced.json',
        ROOT/'src/degali/addons/yawed_crosswind.py',ROOT/'src/degali/addons/observed_wind_profile.py',
        ROOT/'src/degali/core/atmosphere.py',ROOT/'src/degali/core/jetplume.py',
        ROOT/'tests/test_observed_wind_profile.py',ROOT/'tests/test_yawed_crosswind.py',
        ROOT/'docs/prereg-observed-wind-magnitude.md',ROOT/'tools/run_yawed_crosswind_trial10.py',
        ROOT/'tools/run_preslhy_ambient_profile_audit.py',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    field,targets,wind,means,frozen=map(read,paths[:5])
    fast=wind['wind']['10']['fast']
    if fast['status']!='valid_report_no_fault' or fast['range']!='CT30:CU78' or fast['samples']!=49:
        raise ValueError('observed input window changed')
    angle=math.remainder(math.radians(fast['mean_direction_from_deg']+180-75),2*math.pi)
    if args.mode=='magnitude_only': angle=0.
    speed=fast['vector_mean_speed_m_s'];height=fast['height_m']
    directory.mkdir(parents=True,exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),
        mode=args.mode,step_m=args.step,sha256=hashes,wind_speed=speed,reference_height=height,
        relative_wind_angle_deg=math.degrees(angle),source_unchanged=True,candidate_promoted=False))
    report=dict(completed=False,accepted=False,candidate_promoted=False,sha256=hashes)
    try:
        trial=next(t for t in read(paths[5])['trials'] if t['trial']==10)
        original,check=replay(field,trial)
        old=original.result.states[0].copy()
        target5,target_check=original_target_check(original.model,old,field['interfaces']['10'],targets['interfaces']['10'])
        original_f=original.model._as_array(target5)
        f=original_f;target6=np.array([f[0],f[1],f[2],0.,f[3],f[4]])
        copied=with_observed_wind(original.model,speed=speed,reference_height=height)
        model=YawedCrosswind(copied,angle)
        initial=model.match(target6,model.lift(old))
        coarse,fine=model.fluxes(initial),model.fluxes(initial,quadrature_points=2*copied.quadrature_points)
        qe=abs(fine[-1]-coarse[-1])/max(abs(fine[-1]),1)
        fe=float(np.max(abs(coarse-target6)/np.maximum(abs(target6),1)))
        temperature=copied.centre_temperature(model.proxy(initial))
        dt=abs(temperature-targets['interfaces']['10']['target_temperature_K'])
        dw=abs(math.sqrt(2*math.log(2)*initial[2])/targets['interfaces']['10']['target_halfwidth_m']-1)
        handoff=dict(accepted=bool(qe<=1e-5 and fe<=1e-8 and dt<=2 and dw<=.05),initial_state=initial,
            target_fluxes=target6,temperature_K=temperature,temperature_residual_K=dt,halfwidth_residual=dw,
            six_flux_residual=fe,quadrature_residual=qe,original_replay=check,target_check=target_check,
            original_u0=original.model.jetplume.u0,original_z0=original.model.jetplume.z0,
            original_ustar=original.model.jetplume.ustar,new_ustar=copied.jetplume.ustar)
        write_new(directory/'handoff.json',handoff)
        print('Observed wind',args.mode,'handoff:',handoff,flush=True)
        if not handoff['accepted']: raise RuntimeError('original same-source handoff gates failed')
        report['accepted']=True
        def checkpoint(index,data,error):
            write_new(directory/f'checkpoint_{index:05d}{"_failed" if error else ""}.json',
                dict(saved_utc=datetime.now(timezone.utc).isoformat(),error=error,downstream=data,sha256=hashes))
            if index%100==0 or error: print(args.mode,'step',index,error,flush=True)
        arc=field['interfaces']['10']['downstream']['arc_length']
        result=model.solve(initial,distance=arc[-1]-arc[0],step=args.step,checkpoint=checkpoint)
        write_new(directory/'field.json',result)
        if result['maximum_relative_balance_residual']>1e-5: raise RuntimeError('six-flux global ledger failed')
        report.update(scores=score(original,YawedTrajectory(model,result['states']),trial,frozen,means),
            balance_residual=result['maximum_relative_balance_residual'],steps=len(result['states'])-1)
        # Recompute original flux to verify shared model state was not modified.
        after=original.model._as_array(original.model.integral_fluxes(old))
        if np.max(abs(after-original_f)/np.maximum(abs(original_f),1))>1e-8:
            raise ValueError('original model changed')
        report['completed']=True
    except Exception as error:
        report['failure']=repr(error)
        print('FAILED',repr(error),flush=True)
    for name,value in hashes.items():
        if digest(ROOT/name)!=value: raise ValueError('sealed input changed: '+name)
    verify()
    report['finished_utc']=datetime.now(timezone.utc).isoformat()
    write_new(directory/('complete.json' if report['completed'] else 'failure.json'),report)
    if report['completed']:
        for key in ('temperature_summaries','concentration_summaries','geometry_summaries'):
            print(key,report['scores'].get(key),flush=True)


if __name__=='__main__': main()
