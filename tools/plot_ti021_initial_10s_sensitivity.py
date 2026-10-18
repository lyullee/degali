"""Plot saved short-window sensitivities, with censored observation scope."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/integral-local-patch-2026-10-14'


def main():
    report=json.loads((OUT/'initial10s-physics-sensitivity.json').read_text(encoding='utf-8'))
    records={r['case']:r for r in report['cases']}
    data=np.load(OUT/'initial10s-physics-sensitivity.npz')
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    keys=['hybrid_baseline','gravity_off','gravity_half','gravity_double','shear_entrainment_off',
          'ambient_entrainment_off','all_entrainment_off','form_drag_off','compact_observer_only']
    labels=['Baseline','Gravity off','Gravity 0.5x','Gravity 2x','Shear entrainment off',
            'Ambient entrainment off','All entrainment off','Drag off','Compact observer only*']
    ax=axes[0,0];y=np.arange(len(keys));width=.36
    raw=[records[k]['scoring']['0.1']['early_detected_sensor_count'] for k in keys]
    filtered=[records[k]['causal_0p6s_averaged_scoring']['early_detected_sensor_count'] for k in keys]
    ax.barh(y-width/2,raw,width,label='Raw prediction',color='#267eaa')
    ax.barh(y+width/2,filtered,width,label='0.6 s sensor-average proxy',color='#b4cbd5')
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Too-early detections within 0-10 s (count)')
    ax.set_title('Gravity / source shape dominate onset');ax.legend(fontsize=8)
    ax=axes[0,1]
    styles=[('hybrid_baseline','Baseline','#267eaa','-'),('gravity_off','Gravity off','#d98b23','-'),
            ('compact_observer_only','Compact observer*','#00897b','-'),
            ('compact_observer_gravity_off','Compact + gravity off*','#775199','--')]
    for key,label,color,ls in styles:ax.plot(data['time_s'],data[key][0],label=label,color=color,ls=ls)
    ax.axvspan(10,12.7,color='#eeeeee');ax.axvline(12,color='#333',ls=':',label='Observed arrival: 12 s')
    ax.axhline(.1,color='#777',ls=':',lw=.8)
    ax.set(xlim=(0,12.7),xlabel='Time since release (s)',ylabel='S001 concentration (vol%)',title='S001: prediction already detectable at 0.525 s')
    ax.text(10.1,ax.get_ylim()[1]*.65,'Not\ncomputed',fontsize=8)
    ax.legend(fontsize=8,loc='upper left')
    ax=axes[1,0];audit=report['source_profile_audit']
    values=[32.,audit['gaussian_initial_centre_concentration_vol_percent'],audit['compact_initial_centre_concentration_vol_percent']]
    bars=ax.bar(['Uniform source\n(as reported)','Gaussian\n(current mapping)','Compact cylinder*\n(moment matched)'],values,color=['#444','#cc6677','#00897b'])
    for b,v in zip(bars,values):ax.text(b.get_x()+b.get_width()/2,v+1,f'{v:.2f}',ha='center')
    ax.set(ylim=(0,100),ylabel='Initial centre concentration (vol%)',title='Inventory matched does NOT mean profile matched')
    ax=axes[1,1]
    keys2=['hybrid_baseline','gravity_off','gravity_half','all_entrainment_off','ri_entrainment_suppression']
    labels2=['Baseline','Gravity\noff','Gravity\n0.5x','Mixing\noff','Ri\nsuppression']
    values=[records[k]['scoring']['0.1']['prearrival_threshold_excess_vol_percent_s'] for k in keys2]
    bars=ax.bar(labels2,values,color=['#267eaa','#d98b23','#87adbb','#cc6677','#ad719b'])
    for b,v in zip(bars,values):ax.text(b.get_x()+b.get_width()/2,v+8,f'{v:.0f}',ha='center')
    ax.set_ylabel('Sum of pre-arrival excess (vol% s)')
    ax.set_title('Removing mixing can worsen early overprediction')
    fig.suptitle('TI021 R12 | 38 cases, 10 physical seconds each | NOT concentration MAPE',fontsize=14)
    fig.supxlabel('* Structural control, not a validated replacement. No measured 0-10 s concentration traces; timing is censored.',fontsize=9)
    fig.savefig(OUT/'initial10s-physics-sensitivity.png',dpi=170)
    print(OUT/'initial10s-physics-sensitivity.png')


if __name__=='__main__':main()
