"""Evidence-based diagnosis figure and summary; no model parameter fitting."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/integral-local-patch-2026-10-14'


def main():
    names=['diagnostic-flow1000','diagnostic-flow5000','diagnostic-impermeable-ground-flow5000']
    diagnoses=[json.loads((OUT/(n+'.json')).read_text(encoding='utf-8')) for n in names]
    reports=[json.loads(Path(d['run_json']).read_text(encoding='utf-8')) for d in diagnoses]
    dataset=read_thorney_ti021_workbook(ROOT/'ti21_input_ascii.xls')
    old_path=ROOT/'outputs/thorney-ti021-wide-nested-145-sensor-score-600s-x450y640-centered-2026-10-13-v2'
    old=json.loads((old_path/'ti021-nested-145-sensor-score.json').read_text(encoding='utf-8'))
    config=json.loads(Path(old['provenance']['inner_record']['path']).read_text(encoding='utf-8'))['configuration']
    with np.load(old_path/'ti021-nested-145-sensor-traces.npz',allow_pickle=False) as f:
        index=list(f['labels']).index('TI021-S001')
        old_window=(float(f['predicted_arrival'][index]),float(f['predicted_departure'][index]))
        old_doses=[float(f['predicted_dose'][list(f['labels']).index(label)]) for label in ('TI021-S001','TI021-S002')]
        alias=float(np.max(np.abs(f['concentration_vol_percent'][index]
                        -f['concentration_vol_percent'][list(f['labels']).index('TI021-S002')])))
    with np.load(Path(diagnoses[-1]['run_json']).with_suffix('.npz'),allow_pickle=False) as f:
        t=f['time_s'];c=f['predicted_concentration_vol_percent'][0]
        hit=np.flatnonzero(c>=.1)
        new_window=(float(t[hit[0]]),float(t[hit[-1]]))
    observation=dataset.dose_sensors[0]
    summary={
        'confirmed_implementation_errors':[
            {'error':'air ground is pressure-open, species ground is closed',
             'ground_normal_speed_max_m_s':diagnoses[1]['actual_transport_faces']['maximum_abs_ground_velocity_m_s'],
             'species_ground_flux_kg_m2_s':0.},
            {'error':'roof/side exported packets not propagated into exterior concentration',
             'roof_species_export_fraction':diagnoses[1]['final_exit_fraction_of_initial']['z_max']},
            {'error':'0.4 m and 2.4 m observations collapse to identical predictions',
             'old_cell_size_m':config['cell_size_m'],'old_max_abs_trace_difference_vol_percent':alias,
             'observed_sensor1_to_sensor2_dose_ratio':dataset.dose_sensors[0].dose_percent_s/dataset.dose_sensors[1].dose_percent_s},
        ],
        'confirmed_missing_fast_path_physics':diagnoses[1]['missing_density_physics_scale'],
        'ground_only_diagnostic':{
            'intervention':'impermeable free-slip ground; no-slip diagnostic failed separately',
            'before_roof_export_percent':100*diagnoses[1]['final_exit_fraction_of_initial']['z_max'],
            'after_roof_export_percent':100*diagnoses[2]['final_exit_fraction_of_initial']['z_max'],
            'before_dose_mape_percent':diagnoses[1]['dose_score'],
            'after_dose_mape_percent':diagnoses[2]['dose_score'],
            'after_carrier_converged':diagnoses[2]['carrier']['flow_converged'],
        },
        'pre_obstacle_failure':{
            'sensor_local_xy_m':[observation.local_x_m,observation.local_y_m],
            'fence_radius_m':50.,'observed_arrival_departure_s':[observation.arrival_time_s,observation.departure_time_s],
            'old_predicted_arrival_departure_s':list(old_window),
            'sealed_ground_passive_predicted_arrival_departure_s':list(new_window),
            'observed_dose_vol_percent_s':observation.dose_percent_s,
            'old_predicted_dose_vol_percent_s':old_doses[0],
            'new_predicted_dose_vol_percent_s':reports[-1]['sensor_rows'][0]['predicted_dose_vol_percent_s'],
        },
        'interpretation':[
            'Initial inventory and R12 molar conversion are consistent; no factor-of-100 unit error identified.',
            'Ground boundary error and missing roof continuation explain the pathological 5000-step fast carrier.',
            'Sealing the ground does not restore measured exposure; source-region gravity spreading and residence are major remaining candidates.',
            'Density feedback is explicitly absent in the fast path, but present in the old 84.9% path; do not attribute both to an identical missing term.',
            'The common near-ground observation alias is confirmed in both old and new paths.',
            'This is R12 diagnosis, not direct hydrogen accuracy or calibrated performance.',
        ],
        'literature':[
            'https://www.sciencedirect.com/science/article/abs/pii/S0304389403000116',
            'https://www.mdpi.com/2311-5521/6/1/2',
            'https://admlc.com/thorney-island/',
        ],
    }
    (OUT/'failure-cause-summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    labels=['Open ground, 1000','Open ground, 5000','Impermeable ground, 5000']
    for d,label in zip(diagnoses,labels):
        snap=[s for s in d['snapshots'] if s['mass_centroid_height_m'] is not None]
        axes[0,0].plot([s['time_s'] for s in snap],[s['mass_centroid_height_m'] for s in snap],marker='o',label=label)
    axes[0,0].set(title='Where did the gas move?',xlabel='Time (s)',ylabel='Species mass centroid height (m)',xlim=(0,80))
    axes[0,0].legend(fontsize=8)
    axes[0,1].bar(['Open ground','Impermeable ground'],[100*diagnoses[1]['final_exit_fraction_of_initial']['z_max'],100*diagnoses[2]['final_exit_fraction_of_initial']['z_max']])
    axes[0,1].set(title='Uncontinued roof export',ylabel='Initial species mass exported at roof (%)')
    positions=np.arange(2);width=.25
    axes[1,0].bar(positions-width,[dataset.dose_sensors[i].dose_percent_s for i in (0,1)],width,label='Observed')
    axes[1,0].bar(positions,old_doses,width,label='Old model')
    axes[1,0].bar(positions+width,[reports[-1]['sensor_rows'][i]['predicted_dose_vol_percent_s'] for i in (0,1)],width,label='Fast, impermeable ground')
    axes[1,0].set(xticks=positions,xticklabels=['S001: 0.4 m','S002: 2.4 m'],yscale='log',ylabel='Dose (vol% s)',title='Lost near-ground concentration structure')
    axes[1,0].legend(fontsize=8)
    windows=[(observation.arrival_time_s,observation.departure_time_s),old_window,new_window]
    for row,(start,end) in enumerate(windows):
        axes[1,1].plot([start,end],[row,row],lw=8,solid_capstyle='butt')
    axes[1,1].set(yticks=range(3),yticklabels=['Observed S001','Old model','Fast, impermeable ground'],xlabel='Time (s)',title='Detection windows, not measured concentration traces',xlim=(0,320))
    for a in axes.flat:a.grid(alpha=.2)
    fig.suptitle('TI021 R12 failure diagnosis: boundary error plus missing dense-cloud dynamics')
    fig.savefig(OUT/'failure-cause-summary.png',dpi=170)
    plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
