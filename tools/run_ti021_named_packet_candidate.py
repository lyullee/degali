"""Reproduce the layer baseline with one named donor interface change.

No sensor values select the parameters. Read only the recorded configuration
keys before execution; compare dose labels after the full 600 s run.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools.score_ti021_dense_layer_hybrid import run


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--donor',type=Path,required=True)
    parser.add_argument('--interface-x-m',type=float,default=130.)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,default=ROOT/'outputs/ti021-baseline-after-align-default-2026-10-16.json')
    args=parser.parse_args()
    baseline=json.loads(args.baseline.read_text(encoding='utf-8'))
    params={key:value for key,value in baseline.items() if key in inspect.signature(run).parameters}
    params.update(h=baseline['horizontal_cell_size_m'],dt=baseline['max_dt_s'],
        source_dt=baseline['source_time_step_s'],entrainment=baseline['dense_entrainment_model'],
        handoff_interface_path=args.donor.resolve(),interface_x_m=args.interface_x_m)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    config_path=args.output.with_suffix('.configuration.json')
    config_path.write_text(json.dumps({key:str(value) if isinstance(value,Path) else value
        for key,value in params.items()},indent=2),encoding='utf-8')
    started=time.perf_counter()
    result=run(args.output,**params)
    pred=np.array([r['predicted_dose_vol_percent_s'] for r in result['sensor_rows']])
    old=np.array([r['predicted_dose_vol_percent_s'] for r in baseline['sensor_rows']])
    observed=np.array([r['observed_dose_vol_percent_s'] for r in baseline['sensor_rows']])
    positive=observed>0
    same_labels=[r['label'] for r in result['sensor_rows']]==[r['label'] for r in baseline['sensor_rows']]
    interface=result['interface_diagnostics']
    result['paired_baseline']={
        'source':str(args.baseline.resolve()),'identical_sensor_order':same_labels,
        'baseline_all_positive_mape_percent':float(np.mean(np.abs(old[positive]/observed[positive]-1))*100),
        'candidate_all_positive_mape_percent':float(np.mean(np.abs(pred[positive]/observed[positive]-1))*100),
        'configuration_path':str(config_path.resolve())}
    result['donor_file_sha256']=hashlib.sha256(args.donor.read_bytes()).hexdigest()
    result['donor_packet_closure']='marked branch max(total_mass_anomaly, nonnegative species); advective momentum from solved donor velocity'
    result['full_3d_coverage']=False
    result['promotion_allowed']=False
    result['unresolved_coupling']=[
        'vertical impulse and sensible enthalpy do not enter the horizontal receiver',
        'lateral/top donor boundary exports are not coupled',
        'native coarse layer remains upstream; this is an interior flux-replacement diagnostic',
        'near-field sensor histories remain from the layer rather than substituted 3D values']
    result['runtime_with_audit_s']=time.perf_counter()-started
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'paired_baseline':result['paired_baseline'],
        'interface_diagnostics':interface,'promotion_allowed':False},indent=2))


if __name__=='__main__':
    main()
