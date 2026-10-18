"""Matched 600 s, 145-sensor audit of the named-interface iteration."""
import json
from pathlib import Path
import numpy as np
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

ROOT=Path(__file__).resolve().parents[1]


def main():
    dataset=read_thorney_ti021_workbook(ROOT/'ti21_input_ascii.xls')
    paths={
        'baseline':ROOT/'outputs/ti021-baseline-after-align-default-2026-10-16.json',
        'named_packet':ROOT/'outputs/ti021-named-packet-coupled600s-2026-10-17.json',
        'control_subtracted':ROOT/'outputs/ti021-control-subtracted-coupled600s-2026-10-17.json'}
    records={key:json.loads(path.read_text(encoding='utf-8')) for key,path in paths.items()}
    observed=np.array([s.dose_percent_s for s in dataset.dose_sensors])
    positive=observed>0
    labels=[r['label'] for r in records['baseline']['sensor_rows']]
    x=np.array([60+s.local_y_m for s in dataset.dose_sensors])
    h=np.array([s.height_m for s in dataset.dose_sensors])
    groups={'upstream_of_interface':x<130.,'downstream_of_interface':x>=130.}
    groups.update({f'height_{v:g}m':h==v for v in np.unique(h)})
    predictions={}
    results={}
    for key,record in records.items():
        if labels!=[r['label'] for r in record['sensor_rows']]:raise ValueError('sensor order differs')
        actual_obs=np.array([r['observed_dose_vol_percent_s'] for r in record['sensor_rows']])
        if not np.array_equal(observed,actual_obs):raise ValueError('observations differ')
        if record['duration_s']!=600.:raise ValueError('prefix cannot support this audit')
        predictions[key]=np.array([r['predicted_dose_vol_percent_s'] for r in record['sensor_rows']])
        with np.load(paths[key].with_suffix('.npz')) as traces:
            if traces['time_s'][0]!=0. or not np.isclose(traces['time_s'][-1],600.):
                raise ValueError('trace coverage differs')
            if not np.allclose(np.trapezoid(traces['concentration_vol_percent'],traces['time_s'],axis=1),
                               predictions[key],atol=1e-10,rtol=1e-12):raise ValueError('dose not reproduced by traces')
        scores={}
        for name,mask in groups.items():
            active=mask&positive
            scores[name]={'count':int(active.sum()),
                'dose_mape_percent':float(np.mean(np.abs(predictions[key][active]/observed[active]-1))*100)}
        results[key]={'source':str(paths[key].resolve()),'dose_mape_percent':record['dose_mape_percent'],
            'normalized_mean_bias':record['normalized_mean_bias'],'groups':scores,
            'receiver_runtime_s':record['runtime_including_output_s'],
            'species_residual_kg':record['diagnostics']['maximum_species_mass_residual_kg'],
            'horizontal_momentum_residual_kg_m_s':record['diagnostics']['maximum_momentum_budget_residual_kg_m_s'],
            'interface':record.get('interface_diagnostics')}
    for key in ('named_packet','control_subtracted'):
        results[key]['all_groups_nondegrading']=all(
            results[key]['groups'][name]['dose_mape_percent']<=results['baseline']['groups'][name]['dose_mape_percent']
            for name in groups)
        results[key]['accepted']=False
    result={'schema':'degali.ti021-interface-cycle-audit.v1',
        'positive_sensor_count':84,'total_sensor_count':145,
        'trace_dose_reproduction_gate':True,'models':results,
        'full_goal_achieved':False,'active_baseline':'baseline',
        'next_physics':['signed full-mixture anomaly and ambient carrier state rather than a positive-only marked branch',
            'vertical momentum and height-resolved donor state retained by the receiving model',
            'near-ground residence and source-to-wall transport in the 3D donor'],
        'control_diagnostics_recovered_from_histories_only':True}
    path=ROOT/'outputs/ti021-interface-cycle-audit-2026-10-17.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'path':str(path),'scores':{k:v['dose_mape_percent'] for k,v in results.items()},
        'full_goal_achieved':False},indent=2))


if __name__=='__main__':main()
