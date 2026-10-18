"""Compare local-grid levels on one 600 s scenario; emit figure and audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np


def compare(coarse: Path, fine: Path, output: Path):
    a,b = (json.loads(p.read_text(encoding='utf-8')) for p in (coarse,fine))
    for key in ('duration_s','patch_origin_m','patch_extent_m','source_mass_kg','obstacle_present'):
        if a[key]!=b[key]:
            raise ValueError(f'comparison scenarios differ: {key}')
    # Fixed lattice inlet speed means physical flow time scales with h.
    horizons = [r['patch_diagnostics']['carrier_flow_steps']*
                r['patch_diagnostics']['cell_size_m']*.03/2. for r in (a,b)]
    if not np.isclose(*horizons):
        raise ValueError('carrier physical spin-up horizons differ')
    report = {
        'schema':'degali.integral-local-patch-grid-comparison.v1',
        'coarse_report':str(coarse.resolve()),'fine_report':str(fine.resolve()),
        'runtime_including_output_s':{
            '1m':a['runtime_including_output_s'],'0p5m':b['runtime_including_output_s']},
        'reference':'0.5 m result of the same hybrid model, not CFD/experimental truth',
        'matched_carrier_spinup_physical_s':horizons[0],
        'sensor_peak_relative_change_percent':{
            k:100*abs(v-b['sensor_peak_vol_percent'][k])/max(abs(b['sensor_peak_vol_percent'][k]),1e-30)
            for k,v in a['sensor_peak_vol_percent'].items()},
        'sensor_dose_relative_change_percent':{
            k:100*abs(v-b['sensor_dose_vol_percent_s'][k])/max(abs(b['sensor_dose_vol_percent_s'][k]),1e-30)
            for k,v in a['sensor_dose_vol_percent_s'].items()},
    }
    report['peak_refinement_20_percent_gate_passed']=max(report['sensor_peak_relative_change_percent'].values())<=20
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with np.load(coarse.with_suffix('.npz'),allow_pickle=False) as ca, np.load(fine.with_suffix('.npz'),allow_pickle=False) as fi:
        labels=list(ca['sensor_labels'])
        fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
        for index,label in enumerate(labels):
            axis=axes.flat[index]
            axis.plot(ca['time_s'],ca['sensor_vol_percent'][index],label='1 m local grid')
            fine_index=list(fi['sensor_labels']).index(label)
            axis.plot(fi['time_s'],fi['sensor_vol_percent'][fine_index],label='0.5 m local grid')
            axis.set(title=str(label),xlabel='Time (s)',ylabel='H2 (vol%)',xlim=(0,220))
            axis.grid(alpha=.2)
        axes.flat[0].legend(fontsize=8)
        axis=axes.flat[-1]
        axis.bar(['1 m','0.5 m'],[a['runtime_including_output_s'],b['runtime_including_output_s']])
        axis.axhline(30,color='red',linestyle='--',label='30 s target')
        axis.set(ylabel='Wall time incl. output (s)',title='Cold-start 600 s simulation')
        axis.legend(fontsize=8)
        fig.suptitle('Integral exterior + local resolved obstacle patch\nSteady wind / passive isothermal H2; numerical comparison')
        output.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(output.with_suffix('.png'),dpi=160)
        plt.close(fig)
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--coarse',type=Path,required=True)
    parser.add_argument('--fine',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(compare(args.coarse,args.fine,args.output),indent=2))
