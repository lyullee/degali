"""Audit the 0--10 s TI021 timing gap by receptor height.

This audit does not turn arrival summaries into concentration MAPE.  It only
groups the existing censored timing diagnostic by sensor height so a future
multi-layer state has a measurable target.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]


def _summary(rows: list[dict], height: float) -> dict:
    selected=[row for row in rows if abs(float(row["height_m"])-height)<1e-9]
    observed=np.array([row["observed_arrival_s"] for row in selected],float)
    predicted=np.array([row["predicted_arrival_s"] if row["predicted_arrival_s"] is not None else np.nan
                        for row in selected],float)
    known=np.isfinite(observed)
    observed_by_10=known&(observed<=10.)
    detected=np.isfinite(predicted)
    comparable=known&detected
    early=comparable&(predicted<observed)&(predicted<=10.)
    return {
        "height_m":float(height),
        "sensor_count":len(selected),
        "known_arrival_count":int(known.sum()),
        "model_detected_by_10s_count":int(detected.sum()),
        "model_early_count":int(early.sum()),
        "observed_arrival_median_s":float(np.nanmedian(observed[known])) if known.any() else None,
        "observed_arrival_by_10s_count":int(observed_by_10.sum()),
        "observed_arrival_by_10s_median_s":float(np.nanmedian(observed[observed_by_10])) if observed_by_10.any() else None,
        "model_arrival_median_s":float(np.nanmedian(predicted[detected])) if detected.any() else None,
        "median_model_minus_observed_s":float(np.nanmedian(predicted[comparable]-observed[comparable])) if comparable.any() else None,
        "model_non_detection_count":int((known&~detected).sum()),
    }


def audit(output: Path, *, sweep: Path, cases: tuple[str, ...]) -> dict:
    dataset=read_thorney_ti021_workbook(ROOT/'ti21_input_ascii.xls')
    raw=json.loads(sweep.read_text(encoding='utf-8'))
    heights=sorted({float(sensor.height_m) for sensor in dataset.dose_sensors})
    selected=[]
    for case in cases:
        record=next((item for item in raw["cases"] if item["case"]==case),None)
        if record is None:
            raise ValueError(f"case not found: {case}")
        rows=[]
        for item,sensor in zip(record["sensor_rows"],dataset.dose_sensors):
            rows.append({**item,"height_m":sensor.height_m})
        selected.append({
            "case":case,
            "operator":"absolute 0.1 vol percent diagnostic threshold",
            "height_groups":[_summary(rows,height) for height in heights],
        })
    result={
        "schema":"degali.ti021-vertical-arrival-gap-audit.v1",
        "status":"diagnostic_target_for_multilayer_model",
        "concentration_mape_available":False,
        "sweep":str(sweep),
        "cases":selected,
        "interpretation":(
            "Height groups are censored timing diagnostics. They identify which "
            "vertical receptors are missed or reached prematurely; they do not "
            "validate concentration amplitudes."
        ),
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--sweep',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--case',action='append',dest='cases')
    args=parser.parse_args()
    result=audit(args.output,sweep=args.sweep,cases=tuple(args.cases or (
        'inertial_vertical_and_radial','inertial_height_floor_8m_damping_2',
        'inertial_height_floor_8m_damping_3')))
    for case in result['cases']:
        print(case['case'],[(g['height_m'],g['model_arrival_median_s'],g['observed_arrival_median_s'],g['model_non_detection_count'])
            for g in case['height_groups']])
