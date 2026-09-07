"""Resume a sealed yaw checkpoint without overwriting the interrupted run."""
import argparse
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

from degali.addons.yawed_crosswind import YawedCrosswind,YawedTrajectory
from stage2_matched_source import ROOT,REF,read,digest,write_new
from stage1 import verify
from run_preslhy_ambient_profile_audit import replay
from run_yawed_crosswind_trial10 import score


def join_segments(before,after):
    """Keep the shared boundary once, adding the old cumulative six sources."""
    for key in ('states','fluxes','sources'):
        a,b = np.asarray(before[key]),np.asarray(after[key])
        if not np.allclose(a[-1],b[0],rtol=1e-10,atol=1e-10):
            raise ValueError('checkpoint boundary mismatch: '+key)
    if after['arc_length'][0] != 0:
        raise ValueError('resumed segment must start at zero relative arc')
    out = {key:np.concatenate([np.asarray(before[key]),np.asarray(after[key])[1:]])
           for key in ('states','fluxes','sources')}
    out['arc_length'] = np.r_[before['arc_length'],np.asarray(after['arc_length'])[1:]+before['arc_length'][-1]]
    out['cumulative_sources'] = np.concatenate([np.asarray(before['cumulative_sources']),
        np.asarray(after['cumulative_sources'])[1:]+np.asarray(before['cumulative_sources'])[-1]])
    f,c = out['fluxes'],out['cumulative_sources']
    scale = np.maximum.reduce([abs(f),abs(c),np.ones_like(f)])
    out['maximum_relative_balance_residual'] = float(np.max(abs(f-f[0]-c)/scale))
    if not np.all(np.diff(out['arc_length'])>0): raise ValueError('non-increasing joined arc')
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint',type=Path)
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    directory=args.directory.resolve()
    if directory.exists(): raise FileExistsError(directory)
    verify()
    config=read(args.checkpoint.parent/'inputs.json')
    saved=read(args.checkpoint)
    for name,value in config['sha256'].items():
        if digest(ROOT/name)!=value: raise ValueError('sealed original changed: '+name)
    if saved['sha256']!=config['sha256'] or saved['error'] is not None:
        raise ValueError('only an intact accepted checkpoint may resume')
    before=saved['downstream']
    field=read(REF/'phase_ambient_consistency_field_complete_2026-09-05.json')
    trial=next(t for t in read(REF/'e35_reduced.json')['trials'] if t['trial']==10)
    original,check=replay(field,trial)
    model=YawedCrosswind(original.model,config['wind_angle_rad'])
    first=np.asarray(before['states'][-1])
    if not np.allclose(model.fluxes(first),before['fluxes'][-1],rtol=1e-10,atol=1e-10):
        raise ValueError('last checkpoint flux does not replay')
    arc=field['interfaces']['10']['downstream']['arc_length']
    remaining=arc[-1]-arc[0]-before['arc_length'][-1]
    if remaining<=0: raise ValueError('checkpoint already at terminal distance')
    paths=[args.checkpoint,args.checkpoint.parent/'inputs.json',Path(__file__),ROOT/'tests/test_yawed_resume.py']
    hashes={str(p.resolve().relative_to(ROOT)):digest(p) for p in paths}
    directory.mkdir(parents=True,exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),
        original_config=config,resume_sha256=hashes,checkpoint_steps=len(before['states'])-1,
        remaining_arc_m=remaining,source_checkpoint=str(args.checkpoint.resolve()),original_replay=check))
    def checkpoint(index,data,error):
        write_new(directory/f'resumed_checkpoint_{index:05d}{"_failed" if error else ""}.json',
            dict(error=error,downstream=data,sha256=hashes))
        print('resumed segment step',index,flush=True)
    try:
        after=model.solve(first,distance=remaining,step=config['step_m'],checkpoint=checkpoint)
        result=join_segments(before,after)
        if result['maximum_relative_balance_residual']>1e-5: raise ValueError('joined six-flux balance failed')
        write_new(directory/'field.json',result)
        scores=score(original,YawedTrajectory(model,result['states']),trial,
            read(REF/'stage1_single_image_2026-09-06.json'),read(REF/'temperature_mean_audit_2026-09-06.json'))
        for name,value in config['sha256'].items():
            if digest(ROOT/name)!=value: raise ValueError('sealed original changed during resume: '+name)
        for name,value in hashes.items():
            if digest(ROOT/name)!=value: raise ValueError('resume input changed: '+name)
        verify()
        write_new(directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),
            scores=scores,balance_residual=result['maximum_relative_balance_residual'],
            steps=len(result['states'])-1,sha256=config['sha256'],resume_sha256=hashes,
            candidate_promoted=False,stage1_unchanged=True,interrupted_run_preserved=True))
        print(scores['temperature_summaries'],flush=True)
    except Exception as error:
        write_new(directory/'failure.json',dict(completed=False,failure=repr(error),sha256=hashes))
        raise


if __name__=='__main__': main()
