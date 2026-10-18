"""Summarize actual run files; never promote a model from a favorable MAPE."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/integral-local-patch-2026-10-14'


def main():
    old=json.loads((ROOT/'outputs/thorney-ti021-wide-nested-145-sensor-score-600s-x450y640-centered-2026-10-13-v2/ti021-nested-145-sensor-score.json').read_text(encoding='utf-8'))
    passive=json.loads((OUT/'ti021-ground-fixed.json').read_text(encoding='utf-8'))
    names=['ti021-dense-layer-ground-contact-fixed-h2p5.json',
           'ti021-dense-layer-ground-contact-fixed-h2p5-dt0p125.json',
           'ti021-dense-layer-ground-contact-fixed-h2p5-source0p0125.json',
           'ti021-dense-layer-ground-contact-fixed-h1p25.json']
    runs=[json.loads((OUT/name).read_text(encoding='utf-8')) for name in names]
    current=runs[0]
    b=np.array([row['predicted_dose_vol_percent_s'] for row in current['sensor_rows']])
    obs=np.array([row['observed_dose_vol_percent_s'] for row in current['sensor_rows']])
    checks=[]
    for name,run in zip(names[1:],runs[1:]):
        p=np.array([row['predicted_dose_vol_percent_s'] for row in run['sensor_rows']])
        checks.append({'run_file':name,'dose_prediction_relative_l1_change':float(np.abs(p-b).sum()/b.sum()),
                       'dose_mape_change_pp':run['dose_mape_percent']-current['dose_mape_percent']})
    summary={'status':'experimental_not_qualified_not_goal_complete','gas':'R12',
        'duration_s':600.,'sensor_count':145,'positive_observation_count':84,
        'metric':'positive-observed-dose MAPE, not instantaneous H2 concentration MAPE',
        'baseline_mape_percent':100*old['dose_sensor_score']['mean_absolute_percentage_error'],
        'passive_ground_fixed_mape_percent':passive['new_dose_mape_percent'],
        'current_mape_percent':current['dose_mape_percent'],
        'current_end_to_end_runtime_s':current['runtime_including_output_s'],
        'current_report':str(OUT/names[0]),'checks':checks,
        'remaining_bias':current['normalized_mean_bias'],
        'first_height_pair':current['sensor_rows'][:2],
        'unmapped_source_tail_fraction':current['unmapped_species_fraction'],
        'entrained_air_mass_is_not_grid_converged':True}
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    ax=axes[0]
    vals=[summary['baseline_mape_percent'],summary['passive_ground_fixed_mape_percent'],summary['current_mape_percent']]
    bars=ax.bar(['Old\ndensity FV','Fixed ground\npassive patch','Dense layer\nhybrid'],vals,color=['#7f8c8d','#cc6677','#00897b'])
    for bar,value in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,value+1,f'{value:.1f}%',ha='center')
    ax.set_ylim(0,110);ax.set_ylabel('Dose MAPE (%)');ax.set_title('Same 84 positive observed doses')
    ax=axes[1]
    floor=.01
    ax.scatter(obs[obs>0],np.maximum(b[obs>0],floor),s=18,color='#00897b',alpha=.8)
    ax.plot([.01,3000],[.01,3000],color='#555',ls='--',label='1:1')
    ax.set(xscale='log',yscale='log',xlabel='Observed dose (vol% s)',ylabel='Predicted dose (vol% s)',title='Underprediction still substantial')
    ax.set_ylim(.007,3000)
    ax.text(.02,.02,f'{np.sum((obs>0)&(b<floor))} predictions <0.01 shown at floor\nScores use unmodified predictions',
            transform=ax.transAxes,fontsize=8)
    ax.legend()
    ax=axes[2]
    x=np.arange(2);width=.35
    ax.bar(x-width/2,obs[:2],width,label='Observed',color='#444')
    ax.bar(x+width/2,b[:2],width,label='Dense-layer hybrid',color='#00897b')
    ax.set_xticks(x,['S001: 0.4 m','S002: 2.4 m']);ax.set_ylabel('Dose (vol% s)')
    ax.set_title('Height alias removed; doses still too low');ax.legend()
    fig.suptitle('TI021 R12 research hybrid | 600 s event | not H2/CFD field qualification')
    fig.savefig(OUT/'implemented-improvements-summary.png',dpi=170)
    (OUT/'implemented-improvements-summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
