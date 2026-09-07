"""Actual same-source trial10 yaw pilot with immutable ten-step checkpoints."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import math
from pathlib import Path

import numpy as np

from degali.addons.yawed_crosswind import YawedCrosswind,YawedTrajectory
from stage2_matched_source import ROOT,REF,read,digest,write_new
from stage1 import verify
from run_preslhy_ambient_profile_audit import replay,statistics
from run_bounded_liquid_handoff import original_target_check
from audit_stage1_existing_evidence import concentration_metrics,geometry_metrics


def score(original,candidate,trial,frozen,means):
    temperatures,pairs,geometry,missing = [],[],[],[]
    for old in means['rows']:
        if old['trial'] != 10: continue
        row = dict(trial=10,channel=old['channel'],x=old['x'],y=old['y'],z=old['z'],
            **{f'observed_{key}_K':old['stats'][key+'_K'] for key in ('minimum','p05','median','mean')})
        row['control_K'] = original.temperature_at(row['x'],row['y'],row['z'])
        if abs(row['control_K']-old['predictions_K']['CONTROL_K']) > 1e-6:
            raise ValueError('original mean-audit temperature replay failed')
        try: row['candidate_K'] = candidate.temperature_at(row['x'],row['y'],row['z'])
        except ValueError as error:
            missing.append(dict(group='temperature',channel=row['channel'],error=str(error)))
            continue
        temperatures.append(row)
    for old in frozen['concentration_pairs']:
        if old['trial'] != 10: continue
        sensors = [s for s in trial['sensors'] if round(s['x'],3)==old['x']]
        try:
            c = [100*candidate.concentration_at(old['x'],s['y'],s['z']) for s in sensors]
            b = [original.concentration_at(old['x'],s['y'],s['z']) for s in sensors]
        except ValueError as error:
            missing.append(dict(group='concentration',x=old['x'],error=str(error)))
            continue
        pairs.append(dict(old,control=max(b),candidate=max(c),sensor_control=b,sensor_candidate=c))
    for old in frozen['vertical_profiles']:
        if old['trial'] != 10: continue
        baseline = original.state_at(old['x'])
        try: state,lateral = candidate.section(old['x'],0.)
        except ValueError as error:
            missing.append(dict(group='geometry',x=old['x'],error=str(error)))
            continue
        geometry.append(dict(old,control_centre=float(baseline[6]),
            control_sigma=original.model.section_widths(baseline)[1],
            candidate_centre=float(state[9]),candidate_sigma=original.model.section_widths(candidate.model.proxy(state))[1],
            candidate_path_y=float(state[8]),candidate_yaw_deg=math.degrees(state[4]),
            horizontal_normal_lateral=float(lateral)))
    if (len(temperatures),len(pairs),len(geometry)) != (19,5,3) and not missing:
        raise ValueError('unexpected frozen trial10 observation counts')
    result = dict(temperature_rows=temperatures,concentration_pairs=pairs,vertical_profiles=geometry,
        missing=missing,coverage_complete=not missing)
    if not missing:
        result['temperature_summaries'] = {key:{v:statistics(temperatures,v+'_K',f'observed_{key}_K')
            for v in ('control','candidate')} for key in ('minimum','p05','median','mean')}
        result['concentration_summaries'] = {v:concentration_metrics([dict(r,predicted=r[v]) for r in pairs])
                                            for v in ('control','candidate')}
        result['geometry_summaries'] = {v:geometry_metrics([dict(r,modelled_centre=r[v+'_centre'],
            modelled_sigma_z=r[v+'_sigma']) for r in geometry]) for v in ('control','candidate')}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--step',type=float,choices=[.01,.02],default=.02)
    args = parser.parse_args()
    directory = args.directory.resolve()
    if directory.exists(): raise FileExistsError(directory)
    verify()
    source = REF/'phase_ambient_consistency_field_complete_2026-09-05.json'
    targets_path = REF/'gaussian_enthalpy_profile_allflux_interface_2026-09-05.json'
    wind_path = REF/'lateral_temperature_symmetry_2026-09-06.json'
    mean_path = REF/'temperature_mean_audit_2026-09-06.json'
    frozen_path = REF/'stage1_single_image_2026-09-06.json'
    paths = [source,targets_path,wind_path,mean_path,frozen_path,REF/'e35_reduced.json',
        ROOT/'src/degali/addons/yawed_crosswind.py',ROOT/'tests/test_yawed_crosswind.py',
        ROOT/'src/degali/addons/energy_crosswind.py',ROOT/'docs/prereg-yawed-crosswind.md',
        ROOT/'tools/run_preslhy_ambient_profile_audit.py',ROOT/'tools/run_bounded_liquid_handoff.py',
        Path(__file__).resolve()]
    hashes = {str(p.relative_to(ROOT)):digest(p) for p in paths}
    directory.mkdir(parents=True,exist_ok=False)
    field,targets,wind,means,frozen = map(read,[source,targets_path,wind_path,mean_path,frozen_path])
    fast = wind['wind']['10']['fast']
    if fast['status'] != 'valid_report_no_fault' or fast['samples'] != 49 or fast['range'] != 'CT30:CU78':
        raise ValueError('actual wind provenance differs from preregistered valid window')
    angle = math.remainder(math.radians(fast['mean_direction_from_deg']+180-75),2*math.pi)
    trial = next(t for t in read(REF/'e35_reduced.json')['trials'] if t['trial']==10)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        wind_angle_rad=angle,wind_angle_deg=math.degrees(angle),wind_magnitude='unchanged CONTROL height profile',
        selected_trials=[10],step_m=args.step,candidate_promoted=False,trial23_faulty_wind_excluded=True))
    report = dict(completed=False,accepted=False,stage1_unchanged=True,candidate_promoted=False)
    try:
        print('Replaying actual CONTROL trial10',flush=True)
        original,check = replay(field,trial)
        old = original.result.states[0].copy()
        entry = field['interfaces']['10']
        target5,target_check = original_target_check(original.model,old,entry,targets['interfaces']['10'])
        f = original.model._as_array(target5)
        target6 = np.array([f[0],f[1],f[2],0.,f[3],f[4]])
        candidate = YawedCrosswind(original.model,angle)
        zero = YawedCrosswind(original.model,0.)
        checks = []
        for i in np.unique(np.linspace(0,len(original.result.states)-1,9,dtype=int)):
            q = original.result.states[i]
            of = original.model._as_array(original.model.integral_fluxes(q))
            os = original.model.source_terms(q)
            indices = [0,1,2,4,5]
            checks.append(dict(index=int(i),flux_error=float(np.max(abs(zero.fluxes(zero.lift(q))[indices]-of)/np.maximum(abs(of),1))),
                source_error=float(np.max(abs(zero.sources(zero.lift(q))[indices]-os)/np.maximum(abs(os),1)))))
        if max(max(r['flux_error'],r['source_error']) for r in checks)>1e-12:
            raise ValueError('actual coplanar limit failed')
        initial = candidate.match(target6,candidate.lift(old))
        coarse,fine = candidate.fluxes(initial),candidate.fluxes(initial,quadrature_points=2*original.model.quadrature_points)
        q_error = abs(fine[-1]-coarse[-1])/max(abs(fine[-1]),1)
        flux_error = np.max(abs(coarse-target6)/np.maximum(abs(target6),1))
        temp = original.model.centre_temperature(candidate.proxy(initial))
        target_record = targets['interfaces']['10']
        dt = abs(temp-target_record['target_temperature_K'])
        dw = abs(math.sqrt(2*math.log(2)*initial[2])/target_record['target_halfwidth_m']-1)
        handoff = dict(original_replay=check,target_check=target_check,actual_coplanar_checks=checks,
            target_fluxes=target6,initial_state=initial,original_initial_state=old,
            centre_temperature_K=temp,temperature_residual_K=dt,halfwidth_residual=dw,
            six_flux_residual=float(flux_error),quadrature_residual=q_error,
            quadrature_points=original.model.quadrature_points,accepted=bool(dt<=2 and dw<=.05 and flux_error<=1e-8 and q_error<=1e-5))
        write_new(directory/'handoff.json',handoff)
        print('Actual yaw handoff:',handoff['accepted'],'dT=',dt,'width=',dw,'q=',q_error,flush=True)
        if not handoff['accepted']: raise RuntimeError('unchanged original handoff gates failed')
        report['accepted'] = True
        def checkpoint(index,data,error):
            name = f'checkpoint_{index:05d}.json' if error is None else f'failed_checkpoint_{index:05d}.json'
            write_new(directory/name,dict(saved_utc=datetime.now(timezone.utc).isoformat(),
                error=error,downstream=data,sha256=hashes))
            print('Yaw trial10 step',index,'X/Y=',data['states'][-1,7:9], 'error=',error,flush=True)
        arc = entry['downstream']['arc_length']
        result = candidate.solve(initial,distance=arc[-1]-arc[0],step=args.step,checkpoint=checkpoint)
        write_new(directory/'field.json',result)
        if result['maximum_relative_balance_residual']>1e-5:
            raise RuntimeError('six-flux integrated conservation gate failed')
        report['scores'] = score(original,YawedTrajectory(candidate,result['states']),trial,frozen,means)
        report['balance_residual'] = result['maximum_relative_balance_residual']
        report['steps'] = len(result['states'])-1
        report['completed'] = True
    except Exception as error:
        report['failure'] = repr(error)
        print('FAILED:',repr(error),flush=True)
    changed = [name for name,value in hashes.items() if digest(ROOT/name)!=value]
    if changed:
        report['changed_inputs'] = changed
        report['completed'] = False
    verify()
    report.update(finished_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes)
    write_new(directory/('complete.json' if report['completed'] else 'failure.json'),report)
    if report['completed']:
        print(report['scores'].get('temperature_summaries'),flush=True)
        print(report['scores'].get('concentration_summaries'),flush=True)


if __name__ == '__main__': main()
