"""Finite-support, moment-preserving thermochemical variability screen.

Extremisers are mathematical witnesses, NOT inferred or fitted physical PDFs.
Each sensor's original prediction is preserved; no new dispersion score.
"""
from datetime import datetime,timezone
from pathlib import Path
import math

import numpy as np

from stage2_matched_source import ROOT,REF,read,digest,write_new
from run_preslhy_ambient_profile_audit import replay
from audit_downstream_cold_envelope import scalar_shape
from audit_mechanical_thermal_scale import legacy_forward
from degali.addons.thermochemical_ensemble import ensemble,extreme
from stage1 import verify


def reported_x(th,y):
    return (y/th.fuel_molecular_weight)/(y/th.fuel_molecular_weight+(1-y)/th._humid_ambient_molecular_weight)


def parcel(th,t,y,kind='grid'):
    if y==0:
        if t!=th.ambient_temperature: raise ValueError('zero-fuel parcel only at actual ambient')
        return dict(T_K=t,Y_H2=y,h_J_kg=0.,v_m3_kg=1/th.ambient_density,
                    reported_X=0.,gas_X=0.,condensed_kg_kg=[0.,0.,0.],kind=kind)
    f=legacy_forward(th,t,y)
    dry=(1-y)/(1+th.ambient_absolute_humidity)
    masses=dry*np.array([th._dry_nitrogen_mass_fraction,th._dry_oxygen_mass_fraction,th.ambient_absolute_humidity])
    mws=np.array([.0280134,.0319988,th.water_molecular_weight])
    hydrogen=y/th.fuel_molecular_weight
    gas=hydrogen+float(sum((masses-np.array(f['condensed_kg_kg']))/mws))
    return dict(T_K=t,Y_H2=y,h_J_kg=f['specific_enthalpy_J_kg'],v_m3_kg=1/f['density_kg_m3'],
                reported_X=reported_x(th,y),gas_X=hydrogen/gas,condensed_kg_kg=f['condensed_kg_kg'],kind=kind)


def support(th,floor,level):
    count=33 if level=='coarse' else 65
    temperatures=np.linspace(floor,th.ambient_temperature,count)
    fractions=np.unique(np.r_[np.linspace(0,1,count)[1:],np.geomspace(1e-6,.1,count)])
    rows=[parcel(th,float(t),float(y)) for t in temperatures for y in fractions]
    rows.append(parcel(th,th.ambient_temperature,0.,kind='ambient'))
    return rows,fractions


def point(trajectory,raw):
    model=trajectory.model;th=model.thermodynamics
    state=trajectory.state_at(raw['x'])
    sy,sz=model.section_widths(state)
    shape=scalar_shape(state,sy,sz,raw['y'],raw['z'])
    rho=model.rhoa+(state[0]-model.rhoa)*shape
    y=state[0]*state[1]*shape/rho
    t,hv=th._condensed_air_state_exact(np.array([rho]),np.array([y]))
    t,h=float(t[0]),float(hv[0]/rho)
    replayed=trajectory.temperature_at(raw['x'],raw['y'],raw['z'])
    if abs(replayed-raw['CONTROL_K'])>1e-6: raise ValueError('sensor lookup replay mismatch')
    p=parcel(th,t,y,kind='original_homogeneous_witness')
    errors=dict(forward_rho_relative=abs(1/p['v_m3_kg']/rho-1),forward_h_J_kg=abs(p['h_J_kg']-h),
        lookup_replay_K=abs(replayed-raw['CONTROL_K']),exact_minus_lookup_K=t-replayed)
    if errors['forward_rho_relative']>1e-7 or errors['forward_h_J_kg']>1e-5:
        raise ValueError('independent original EOS replay failed')
    # Keep exactly the original inverse-EOS state, whose finite inversion
    # tolerance is explicitly measured above; don't alter the target moments.
    p.update(h_J_kg=h,v_m3_kg=1/rho)
    original_x=model.point_mole_fraction(state,raw['y'],raw['z'])
    if abs(original_x-p['reported_X'])>1e-12: raise ValueError('original concentration operator differs')
    return p,errors


def evaluate(grid,original,raw,hold_reported):
    keys=('T_K','Y_H2','h_J_kg','v_m3_kg','reported_X','gas_X')
    t,y,h,v,x,gx=(np.array([r[k] for r in grid]) for k in keys)
    vm=original['v_m3_kg']
    targets=dict(mean_fraction=original['Y_H2'],mean_enthalpy=original['h_J_kg'],mean_specific_volume=vm)
    if hold_reported:
        targets.update(additional_moments=[v*x],additional_targets=[vm*original['reported_X']])
    objectives=dict(mean_T_min=(t*v/vm,False),mean_T_max=(t*v/vm,True))
    for label in ('p05','median'):
        objectives[f'max_CDF_at_observed_{label}']=(v/vm*(t<=raw[f'observed_{label}_K']),True)
    results={}
    for label,(cost,maximize) in objectives.items():
        r=extreme(y,h,v,values=cost,maximize=maximize,**targets)
        e=ensemble(r.weights,t,y,h,v)
        active=np.flatnonzero(r.weights>0)
        witness=[dict(**grid[int(i)],mass_weight=float(r.weights[i]),volume_weight=float(e['volume_weights'][i])) for i in active]
        results[label]=dict(value=r.objective,maximum_scaled_moment_error=r.maximum_scaled_moment_error,
            maximum_scaled_dual_violation=r.maximum_scaled_dual_violation,scaled_duality_gap=r.scaled_duality_gap,
            moment_residuals=r.moment_residuals.tolist(),support_size=len(grid),witness=witness,
            favre_temperature_K=e['favre_temperature'],reynolds_temperature_K=e['reynolds_temperature'],
            favre_Y_variance=e['favre_fraction_variance'],favre_T_variance=e['favre_temperature_variance'],
            favre_Y_T_covariance=e['favre_temperature_fraction_covariance'],
            reynolds_reported_X=float(e['volume_weights']@x),reynolds_gas_X=float(e['volume_weights']@gx))
    return results


def main():
    directory=REF/'thermochemical_ensemble_2026-09-06'
    paths=[REF/'phase_ambient_consistency_field_complete_2026-09-05.json',REF/'stage1_single_image_2026-09-06.json',
        REF/'e35_reduced.json',ROOT/'docs/prereg-thermochemical-ensemble.md',
        ROOT/'src/degali/addons/thermochemical_ensemble.py',ROOT/'tests/test_thermochemical_ensemble.py',
        ROOT/'tools/audit_mechanical_thermal_scale.py',ROOT/'tools/audit_downstream_cold_envelope.py',
        ROOT/'tools/run_preslhy_ambient_profile_audit.py',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    directory.mkdir(exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        floors_K=[50.,65.],nested_levels=[33,65],source_and_observations_unchanged=True))
    verify()
    field,frozen=read(paths[0]),read(paths[1]);trials={r['trial']:r for r in read(paths[2])['trials']}
    rows=[];checks={};failures=[]
    for n in (10,23):
        trajectory,checks[str(n)]=replay(field,trials[n]);th=trajectory.model.thermodynamics
        selected=[r for r in frozen['temperature_rows'] if r['trial']==n]
        for floor in (50.,65.):
            grids={}
            for level in ('coarse','fine'):
                print('Building support',n,floor,level,flush=True)
                grids[level]=support(th,floor,level)
                write_new(directory/f'support_{n}_{int(floor)}_{level}.json',dict(trial=n,floor_K=floor,level=level,
                    rows=grids[level][0],fractions=grids[level][1]))
            for index,raw in enumerate(selected):
                original,errors=point(trajectory,raw)
                row=dict(trial=n,channel=raw['channel'],x=raw['x'],y=raw['y'],z=raw['z'],floor_K=floor,
                    original=original,original_replay=errors,original_outside_nominal_support=not floor<=original['T_K']<=th.ambient_temperature,
                    observed_p05_K=raw['observed_p05_K'],observed_median_K=raw['observed_median_K'],levels={})
                for level,(base,fractions) in grids.items():
                    extra=[]
                    for threshold in sorted({raw['observed_p05_K'],raw['observed_median_K']}):
                        if floor<=threshold<=th.ambient_temperature:
                            extra.extend(parcel(th,float(threshold),float(y),kind='CDF_threshold') for y in fractions)
                    grid=base+extra+[original]
                    row['levels'][level]={}
                    for hold in (False,True):
                        mode='hold_reported_X' if hold else 'four_moments'
                        try:
                            row['levels'][level][mode]=evaluate(grid,original,raw,hold)
                        except (ValueError,RuntimeError,FloatingPointError) as error:
                            row['levels'][level][mode]=dict(failed=True,error=str(error))
                            failures.append(dict(trial=n,channel=raw['channel'],floor_K=floor,level=level,mode=mode,error=str(error)))
                refinement={}
                for mode in ('four_moments','hold_reported_X'):
                    coarse,fine=(row['levels'][level][mode] for level in ('coarse','fine'))
                    if coarse.get('failed') or fine.get('failed'):
                        refinement[mode]=dict(resolved=False,failed=True);continue
                    diffs={key:abs(fine[key]['value']-coarse[key]['value']) for key in coarse}
                    resolved=max(diffs['mean_T_min'],diffs['mean_T_max'])<=.5 and max(diffs[k] for k in diffs if k.startswith('max_CDF'))<=.01
                    # A nested larger feasible support cannot contract extrema.
                    monotone=(fine['mean_T_min']['value']<=coarse['mean_T_min']['value']+1e-6 and all(fine[k]['value']>=coarse[k]['value']-1e-6 for k in coarse if k!='mean_T_min'))
                    if not monotone: raise RuntimeError('nested-support monotonicity failed')
                    refinement[mode]=dict(resolved=resolved,differences=diffs,nested_extrema_monotone=monotone)
                row['refinement']=refinement
                rows.append(row)
                write_new(directory/f'point_{n}_{int(floor)}_{index:02d}.json',row)
                print('Point',n,floor,index+1,'/',len(selected),raw['channel'],{k:v['resolved'] for k,v in refinement.items()},flush=True)
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('ensemble input changed')
    verify()
    summary={}
    for floor in (50.,65.):
        for mode in ('four_moments','hold_reported_X'):
            valid=[r for r in rows if r['floor_K']==floor and not r['levels']['fine'][mode].get('failed')]
            summary[f'{int(floor)}_{mode}']=dict(n=len(valid),resolved=sum(r['refinement'][mode]['resolved'] for r in valid),
                maximum_mean_T_reduction_K=max((r['original']['T_K']-r['levels']['fine'][mode]['mean_T_min']['value'] for r in valid),default=None),
                maximum_mean_T_increase_K=max((r['levels']['fine'][mode]['mean_T_max']['value']-r['original']['T_K'] for r in valid),default=None),
                grid_witness_can_reach_median=sum(r['levels']['fine'][mode]['max_CDF_at_observed_median']['value']>=.5-1e-8 for r in valid),
                grid_witness_can_reach_p05=sum(r['levels']['fine'][mode]['max_CDF_at_observed_p05']['value']>=.05-1e-8 for r in valid))
    write_new(directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        rows=rows,summary=summary,failures=failures,replay_checks=checks,expected_points=82,actual_points=len(rows),
        phase_closure='original independent pure condensates; no new mixed-liquid or solid EOS',
        finite_support_only=True,source_distribution_not_identified=True,new_dispersion_field=False,candidate_promoted=False))
    print('Complete:',summary,'failures',len(failures),flush=True)


if __name__=='__main__': main()
