from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/integral-local-patch-2026-10-14'


def main():
    report=json.loads((OUT/'initial10s-repair-iterations.json').read_text(encoding='utf-8'))
    records={r['case']:r for r in report['cases']}
    old=np.load(OUT/'initial10s-physics-sensitivity.npz')
    new=np.load(OUT/'initial10s-repair-iterations.npz')
    fig,axes=plt.subplots(1,3,figsize=(13,4.3),layout='constrained')
    ax=axes[0]
    ax.plot(old['time_s'],old['hybrid_baseline'][0],label='Old Gaussian',color='#cc6677')
    ax.plot(new['time_s'],new['compact_legacy_radial_handoff'][0],label='Uniform prism',color='#d99735')
    ax.plot(new['time_s'],new['inertial_vertical_and_radial'][0],label='Prism + initial inertia',color='#00897b')
    ax.axhline(.1,lw=.7,color='#555',ls=':')
    ax.set(xlim=(0,10),xlabel='Time (s)',ylabel='S001 prediction (vol%)',title='Observed arrival: 12 s (outside window)')
    ax.legend(fontsize=8)
    ax=axes[1]
    keys=['compact_legacy_radial_handoff','inertial_vertical_and_radial']
    vals=[report['baseline']['scoring']['0.1']['early_detected_sensor_count']]+[records[k]['scores']['0.1']['early_detected_sensor_count'] for k in keys]
    bars=ax.bar(['Old\nGaussian','Uniform\nprism','Prism +\ninertia'],vals,color=['#cc6677','#d99735','#00897b'])
    for b,v in zip(bars,vals):ax.text(b.get_x()+b.get_width()/2,v+.3,str(v),ha='center')
    ax.set(ylim=(0,20),ylabel='Too-early detection count / 82',title='Censored at 10 s, not full arrival MAE')
    ax=axes[2]
    vals=[report['baseline']['initial_centre_vol_percent'],32.,32.]
    bars=ax.bar(['Old\nGaussian','Uniform\nprism','Prism +\ninertia'],vals,color=['#cc6677','#d99735','#00897b'])
    ax.axhline(32,color='#333',ls='--',label='Required source: 32 vol%')
    for b,v in zip(bars,vals):ax.text(b.get_x()+b.get_width()/2,v+1,str(round(v,2)),ha='center')
    ax.set(ylim=(0,100),ylabel='Initial centre (vol%)',title='Initial profile fixed without clipping mass')
    ax.legend(fontsize=8)
    fig.suptitle('TI021 R12 | 10-second repair iterations | research candidate, not validated concentration accuracy')
    fig.savefig(OUT/'initial10s-repair-iterations.png',dpi=170)
    print(OUT/'initial10s-repair-iterations.png')


if __name__=='__main__':main()
