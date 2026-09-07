"""Same-moment row-conditioning and targeted nested-grid refinement."""
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

from stage2_matched_source import ROOT,REF,read,digest,write_new
from run_preslhy_ambient_profile_audit import thermodynamics
from audit_thermochemical_ensemble import parcel
from degali.addons.thermochemical_ensemble import ensemble
from degali.addons.thermochemical_moment_solver_v2 import whitened_extreme
from stage1 import verify


OBJECTIVES=('mean_T_min','mean_T_max','max_CDF_at_observed_p05','max_CDF_at_observed_median')


def solve(grid,row,mode,label):
    keys=('T_K','Y_H2','h_J_kg','v_m3_kg','reported_X','gas_X')
    t,y,h,v,x,gx=(np.array([r[k] for r in grid]) for k in keys)
    original=row['original'];vm=original['v_m3_kg']
    target=dict(mean_fraction=original['Y_H2'],mean_enthalpy=original['h_J_kg'],mean_specific_volume=vm)
    if mode=='hold_reported_X':
        target.update(additional_moments=[v*x],additional_targets=[vm*original['reported_X']])
    if label.startswith('mean_'): cost=t*v/vm
    else: cost=v/vm*(t<=row['observed_'+label.split('_')[-1]+'_K'])
    r=whitened_extreme(y,h,v,values=cost,maximize=label!='mean_T_min',**target)
    e=ensemble(r.weights,t,y,h,v)
    active=np.flatnonzero(r.weights>0)
    return dict(value=r.objective,maximum_scaled_moment_error=r.maximum_scaled_moment_error,
        maximum_scaled_dual_violation=r.maximum_scaled_dual_violation,scaled_duality_gap=r.scaled_duality_gap,
        moment_residuals=r.moment_residuals.tolist(),support_size=len(grid),
        witness=[dict(**grid[int(i)],mass_weight=float(r.weights[i]),volume_weight=float(e['volume_weights'][i])) for i in active],
        reynolds_temperature_K=e['reynolds_temperature'],favre_temperature_K=e['favre_temperature'],
        reynolds_reported_X=float(e['volume_weights']@x),reynolds_gas_X=float(e['volume_weights']@gx),
        favre_Y_variance=e['favre_fraction_variance'],favre_T_variance=e['favre_temperature_variance'])


def main():
    olddir=REF/'thermochemical_ensemble_2026-09-06'
    directory=REF/'thermochemical_ensemble_refinement_2026-09-06'
    paths=[olddir/'complete.json',olddir/'inputs.json',REF/'e35_reduced.json',
        ROOT/'docs/prereg-thermochemical-refinement.md',ROOT/'src/degali/addons/thermochemical_moment_solver_v2.py',
        ROOT/'tests/test_thermochemical_moment_solver_v2.py',ROOT/'tools/audit_thermochemical_ensemble.py',Path(__file__).resolve()]
    paths.extend(sorted(olddir.glob('support_*.json')))
    prior=read(paths[0]);hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    for path,expected in prior['sha256'].items():
        if digest(ROOT/path)!=expected: raise ValueError('original ensemble input changed')
    directory.mkdir(exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        original_failures_preserved=len(prior['failures']),same_physical_model_and_original_gates=True))
    verify()
    trials={r['trial']:r for r in read(paths[2])['trials']};thermos={};grids={};rows=[];failures=[]
    def base(n,floor,level):
        key=(n,floor,level)
        if key not in grids:
            if level!='extra_fine':
                d=read(olddir/f'support_{n}_{int(floor)}_{level}.json')
                grids[key]=(d['rows'],np.array(d['fractions']))
            else:
                th=thermos[n];count=129
                ts=np.linspace(floor,th.ambient_temperature,count)
                ys=np.unique(np.r_[np.linspace(0,1,count)[1:],np.geomspace(1e-6,.1,count)])
                print('Extra-fine support',n,floor,flush=True)
                g=[parcel(th,float(t),float(y)) for t in ts for y in ys]
                g.append(parcel(th,th.ambient_temperature,0.,kind='ambient'))
                write_new(directory/f'support_{n}_{int(floor)}_extra_fine.json',dict(rows=g,fractions=ys))
                grids[key]=(g,ys)
        return grids[key]
    for index,old in enumerate(prior['rows']):
        n,floor=old['trial'],old['floor_K']
        if n not in thermos: thermos[n]=thermodynamics(trials[n],consistent=True)
        th=thermos[n]
        row={k:old[k] for k in ('trial','channel','x','y','z','floor_K','original','original_replay','observed_p05_K','observed_median_K')}
        row.update(levels={},resolution={})
        for level in ('coarse','fine','extra_fine'):
            wanted={mode:OBJECTIVES if level!='extra_fine' else tuple(label for label,data in row['resolution'][mode].items() if not data['resolved']) for mode in ('four_moments','hold_reported_X')}
            if not any(wanted.values()): continue
            b,ys=base(n,floor,level);extra=[]
            for threshold in sorted({row['observed_p05_K'],row['observed_median_K']}):
                if floor<=threshold<=th.ambient_temperature:
                    extra.extend(parcel(th,float(threshold),float(y),kind='CDF_threshold') for y in ys)
            grid=b+extra+[row['original']];row['levels'][level]={}
            for mode,labels in wanted.items():
                row['levels'][level][mode]={}
                for label in labels:
                    try: row['levels'][level][mode][label]=solve(grid,row,mode,label)
                    except (ValueError,RuntimeError,FloatingPointError) as error:
                        row['levels'][level][mode][label]=dict(failed=True,error=str(error))
                        failures.append(dict(index=index,trial=n,channel=row['channel'],floor_K=floor,level=level,mode=mode,objective=label,error=str(error)))
                if level!='coarse':
                    previous='coarse' if level=='fine' else 'fine'
                    if level=='fine': row['resolution'][mode]={}
                    for label in labels:
                        a=row['levels'][previous][mode][label];z=row['levels'][level][mode][label]
                        limit=.5 if label.startswith('mean_') else .01
                        if a.get('failed') or z.get('failed'):
                            row['resolution'][mode][label]=dict(resolved=False,last_level=level,solver_failure=True);continue
                        change=abs(a['value']-z['value'])
                        monotone=z['value']<=a['value']+1e-6 if label=='mean_T_min' else z['value']>=a['value']-1e-6
                        if not monotone: raise RuntimeError('conditioned nested extrema are not monotone')
                        row['resolution'][mode][label]=dict(resolved=change<=limit,last_level=level,absolute_difference=change,gate=limit)
        rows.append(row);write_new(directory/f'point_{index:02d}.json',row)
        print('Refined',index+1,'/82',n,floor,row['channel'],{mode:sum(v['resolved'] for v in data.values()) for mode,data in row['resolution'].items()},flush=True)
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('refinement input changed')
    verify()
    summary={}
    for floor in (50.,65.):
        for mode in ('four_moments','hold_reported_X'):
            selected=[r for r in rows if r['floor_K']==floor]
            report={}
            for label in OBJECTIVES:
                good=[]
                for row in selected:
                    level=row['resolution'][mode][label]['last_level']
                    item=row['levels'][level][mode][label]
                    if not item.get('failed'): good.append((row,item))
                report[label]=dict(successful=len(good),resolved=sum(r['resolution'][mode][label]['resolved'] for r,_ in good),
                    range=[min(v['value'] for _,v in good),max(v['value'] for _,v in good)] if good else None)
                if label=='mean_T_min': report[label]['maximum_reduction_K']=max(r['original']['T_K']-v['value'] for r,v in good)
                if label.startswith('max_CDF'):
                    q=.05 if label.endswith('p05') else .5
                    report[label]['grid_witness_reaches_quantile']=sum(v['value']>=q-1e-8 for _,v in good)
            summary[f'{int(floor)}_{mode}']=report
    write_new(directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        rows=rows,summary=summary,failures=failures,original_records_preserved=True,
        finite_support_only=True,no_physical_parameters_added=True,new_dispersion_field=False,candidate_promoted=False))
    print('Done',summary,'solver failures',len(failures),flush=True)


if __name__=='__main__': main()
