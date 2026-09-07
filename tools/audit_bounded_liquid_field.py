"""Verify and compare actual bounded-liquid trial10, not the full7-case set."""
from __future__ import annotations

import argparse
from copy import copy
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from run_bounded_liquid_handoff import clone_with_liquid
from run_preslhy_ambient_profile_audit import replay
from stage2_matched_source import ROOT,REF,read,digest,write_new,fixed_scores
from stage1 import verify


def key(row,kind):
    if kind=='temperature_rows': return (row['trial'],row['channel'],row['x'],row['y'],row['z'])
    return (row['trial'],row['x'])


def audit(directory,output):
    directory=directory.resolve()
    output=output.resolve()
    if output.exists(): raise FileExistsError(output)
    complete_path=directory/'complete.json'
    completed=read(complete_path)
    if not completed.get('completed') or not completed['fixed_scores']['coverage_complete']:
        raise ValueError('complete actual trial10 coverage required before this audit')
    if completed['seven_case_claim']: raise ValueError('single case was mislabeled seven-case')
    verify()
    for mapping in ('input_sha256','boundary_input_sha256'):
        for path,expected in completed[mapping].items():
            if digest(ROOT/path)!=expected: raise ValueError(f'candidate input changed:{path}')
    field_path=REF/'phase_ambient_consistency_field_complete_2026-09-05.json'
    trial=next(t for t in read(REF/'e35_reduced.json')['trials'] if t['trial']==10)
    old,oldcheck=replay(read(field_path),trial)
    newmodel=clone_with_liquid(old.model)
    newmodel.quadrature_points=completed['energy_quadrature_points']
    run=SimpleNamespace(**{k:np.asarray(v) if isinstance(v,list) else v for k,v in completed['downstream'].items()})
    ts,_=newmodel.thermodynamics._condensed_air_state(run.states[:,0],run.states[:,1])
    terror=float(np.max(abs(ts-run.temperatures)))
    indices=np.unique(np.linspace(0,len(run.states)-1,7,dtype=int))
    flux_errors=[]
    for i in indices:
        f=newmodel._as_array(newmodel.integral_fluxes(run.states[i]))
        flux_errors.append(float(np.max(abs(f-run.fluxes[i])/np.maximum(abs(run.fluxes[i]),1))))
    if terror>1e-6 or max(flux_errors)>1e-9: raise ValueError('new field replay changed')
    residual=run.fluxes-run.fluxes[0]-run.cumulative_sources
    balance=float(np.max(abs(residual)/np.maximum(abs(run.fluxes[0]),1)))
    if balance>1e-5: raise ValueError('global conservation failed')
    old_i=SimpleNamespace(model=old.model,accepted=True,downstream_result=old.result)
    # frozen replay's lightweight object has only states; fixed_scores needs
    # the old complete field's independently stored conservation maximum.
    old_i.downstream_result.maximum_relative_balance_residual=read(field_path)['interfaces']['10']['downstream']['maximum_relative_balance_residual']
    new_i=SimpleNamespace(model=newmodel,accepted=True,downstream_result=run)
    oldscore=fixed_scores(SimpleNamespace(handoffs={10:old_i},failures={}),[10])
    newscore=fixed_scores(SimpleNamespace(handoffs={10:new_i},failures={}),[10])
    # Isolate the legacy2D lookup's measurement effect on IDENTICAL old
    # states. This is NOT an exact-EOS reintegration or a new control field.
    exact_model=copy(old.model)
    exact_model.thermodynamics=copy(old.model.thermodynamics)
    exact_model.thermodynamics._condensed_air_state=old.model.thermodynamics._condensed_air_state_exact
    exact_i=SimpleNamespace(model=exact_model,accepted=True,downstream_result=old.result)
    exact_score=fixed_scores(SimpleNamespace(handoffs={10:exact_i},failures={}),[10])
    exact_t,_=old.model.thermodynamics._condensed_air_state_exact(old.result.states[:,0],old.result.states[:,1])
    stored_t=np.array(read(field_path)['interfaces']['10']['downstream']['temperatures'])
    lookup_diagnostic=dict(same_old_states_only=True,not_a_reintegrated_control=True,
        maximum_centre_temperature_lookup_difference_K=float(np.max(abs(exact_t-stored_t))),
        maximum_sensor_temperature_lookup_difference_K=max(abs(a['matched_K']-b['matched_K'])
            for a,b in zip(oldscore['temperature_rows'],exact_score['temperature_rows'])),
        exact_operator_temperature_summaries=exact_score['temperature_summaries'])
    if not oldscore['coverage_complete'] or not newscore['coverage_complete']:
        raise ValueError('matched comparison has missing trial10 observations')
    rows={}
    for kind in ('concentration_pairs','vertical_profiles','temperature_rows'):
        a,b=oldscore[kind],newscore[kind]
        stored=completed['fixed_scores'][kind]
        if [key(r,kind) for r in a]!=[key(r,kind) for r in b] or [key(r,kind) for r in b]!=[key(r,kind) for r in stored]:
            raise ValueError('comparison observation keys changed')
        paired=[]
        for ra,rb,rs in zip(a,b,stored):
            if kind=='temperature_rows':
                if abs(rb['matched_K']-rs['matched_K'])>1e-7: raise ValueError('stored sensorT does not replay')
                paired.append(dict(trial=10,channel=ra['channel'],x=ra['x'],y=ra['y'],z=ra['z'],
                    observed_minimum_K=ra['observed_minimum_K'],observed_p05_K=ra['observed_p05_K'],
                    observed_median_K=ra['observed_median_K'],original_K=ra['matched_K'],candidate_K=rb['matched_K'],
                    candidate_minus_original_K=rb['matched_K']-ra['matched_K']))
            elif kind=='concentration_pairs':
                if abs(rb['predicted']/rs['predicted']-1)>1e-10: raise ValueError('stored arcC does not replay')
                paired.append(dict(trial=10,x=ra['x'],observed=ra['observed'],original=ra['predicted'],candidate=rb['predicted']))
            else:
                if any(abs(rb[k]-rs[k])>1e-10 for k in ('modelled_centre','modelled_sigma_z')):
                    raise ValueError('stored geometry does not replay')
                paired.append(dict(original=ra,candidate=rb))
        rows[kind]=paired
    files=[complete_path,field_path,Path(__file__).resolve(),ROOT/'reference/preslhy/stage1_single_image_2026-09-06.json']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in files}
    result=dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),input_sha256=hashes,
        underlying_candidate_hashes_verified=True,original_replay=oldcheck,
        candidate_replay=dict(temperature_max_abs_K=terror,sampled_flux_max_scaled_error=max(flux_errors),
            global_conservation_max_scaled=balance,field_steps=len(run.states)-1,x_last_m=float(run.states[-1,5])),
        counts={k:len(v) for k,v in rows.items()},paired_rows=rows,
        original_metrics=oldscore['metrics'],candidate_metrics=newscore['metrics'],
        original_temperature_summaries=oldscore['temperature_summaries'],candidate_temperature_summaries=newscore['temperature_summaries'],
        maximum_sensor_temperature_change_K=max(abs(r['candidate_minus_original_K']) for r in rows['temperature_rows']),
        legacy_lookup_measurement_diagnostic=lookup_diagnostic,
        trial10_only=True,full_seven_case_validated=False,resolution_test_complete=False,candidate_promoted=False)
    verify()
    write_new(output,result)
    print({k:v for k,v in result.items() if k not in ('paired_rows','input_sha256')},flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('output',type=Path)
    a=parser.parse_args()
    audit(a.directory,a.output)


if __name__=='__main__': main()
