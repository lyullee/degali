"""Checkpointed actual field for individually accepted bounded-liquid trial10."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from run_bounded_liquid_handoff import clone_with_liquid
from run_preslhy_ambient_profile_audit import replay
from stage2_matched_source import ROOT,REF,read,digest,write_new,fixed_scores
from stage1 import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    if args.directory.exists(): raise FileExistsError(args.directory)
    verify()
    boundary_path=REF/'bounded_liquid_handoff_2026-09-06.json'
    boundary=read(boundary_path)
    row=boundary['interfaces']['10']
    if not row['accepted'] or not row['projection']['success']:
        raise ValueError('trial10 does not have a passed unchanged interface gate')
    for path,expected in boundary['input_sha256'].items():
        if digest(ROOT/path)!=expected: raise ValueError(f'boundary input changed:{path}')
    field_path=REF/'phase_ambient_consistency_field_complete_2026-09-05.json'
    frozen=read(field_path)
    trial=next(t for t in read(REF/'e35_reduced.json')['trials'] if t['trial']==10)
    trajectory,replay_check=replay(frozen,trial)
    model=clone_with_liquid(trajectory.model)
    model.quadrature_points=row['projection']['quadrature_points']
    state=np.array(row['projection']['state'])
    target=np.array([row['target_fluxes'][k] for k in ('total_mass','contaminant_mass','momentum_x','momentum_z','energy')])
    represented=model._as_array(model.integral_fluxes(state))
    if np.max(abs(represented-target)/np.maximum(abs(target),1))>1e-8:
        raise ValueError('accepted boundary replay changed')
    duration=frozen['interfaces']['10']['downstream']['arc_length'][-1]
    steps=int(math.ceil(duration/.02))
    step=duration/steps
    paths=[boundary_path,field_path,Path(__file__).resolve(),ROOT/'tools/stage2_matched_source.py']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    args.directory.mkdir(parents=True,exist_ok=False)
    write_new(args.directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=hashes,boundary_input_sha256=boundary['input_sha256'],trial=10,
        maximum_step_m=.02,actual_step_m=step,planned_steps=steps,
        original_replay_check=replay_check,individual_trial_only=True,candidate_promoted=False))
    parts=[]
    done=0
    try:
        while done<steps:
            count=min(10,steps-done)
            result=model.solve(state,maximum_distance=count*step,maximum_step=step*(1+1e-12),relative_tolerance=1e-5)
            if len(result.states)!=count+1: raise RuntimeError('chunk step count changed')
            write_new(args.directory/f'chunk_{done:04d}_{done+count:04d}.json',dict(start_step=done,
                end_step=done+count,downstream=vars(result),finished_utc=datetime.now(timezone.utc).isoformat()))
            parts.append(result)
            done+=count
            state=result.states[-1].copy()
            print(f'Trial10 bounded-liquid saved {done}/{steps}, x={state[5]:.4f}m, T={result.temperatures[-1]:.4f}K',flush=True)
    except (ValueError,RuntimeError,FloatingPointError) as err:
        write_new(args.directory/'failure.json',dict(finished_utc=datetime.now(timezone.utc).isoformat(),
            completed_steps=done,last_saved_state=state,last_saved_x_m=state[5],error=str(err),
            rejected_solver_domain_probes=model.thermodynamics.domain_failures,
            last_domain_probe=model.thermodynamics.last_domain_failure,
            new_complete_field=False,limitation='Domain probes may be rejected nonlinear solver trials, not a proven physical trajectory crossing.'))
        print('Trial10 bounded-liquid field stopped:',err,flush=True)
        verify()
        return
    def join(name):
        return np.concatenate([getattr(p,name) if i==0 else getattr(p,name)[1:] for i,p in enumerate(parts)])
    arc=np.linspace(0,duration,steps+1)
    cumulative=[]
    offset=np.zeros(5)
    for i,p in enumerate(parts):
        cs=p.cumulative_sources+offset
        cumulative.extend(cs if i==0 else cs[1:])
        offset=cs[-1].copy()
    flux=join('fluxes')
    sources=np.array(cumulative)
    residuals=flux-flux[0]-sources
    maximum=float(np.max(abs(residuals)/np.maximum(abs(flux[0]),1)))
    combined=SimpleNamespace(arc_length=arc,states=join('states'),fluxes=flux,
        sources=join('sources'),cumulative_sources=sources,balance_residuals=residuals,
        temperatures=join('temperatures'),centre_mole_fractions=join('centre_mole_fractions'),
        maximum_relative_balance_residual=maximum)
    interface=SimpleNamespace(model=model,accepted=True,downstream_result=combined)
    scores=fixed_scores(SimpleNamespace(handoffs={10:interface},failures={}),[10])
    assert hashes=={str(p.relative_to(ROOT)):digest(p) for p in paths}
    for path,expected in boundary['input_sha256'].items():
        if digest(ROOT/path)!=expected: raise ValueError(f'boundary input changed during run:{path}')
    verify()
    write_new(args.directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=hashes,boundary_input_sha256=boundary['input_sha256'],downstream=vars(combined),
        fixed_scores=scores,energy_quadrature_points=model.quadrature_points,
        candidate_promoted=False,seven_case_claim=False,trial23_unsupported=True))
    print('Trial10 bounded-liquid field complete. Conservation:',maximum,flush=True)


if __name__=='__main__': main()
