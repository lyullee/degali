"""Diagnostic screen for a volume-preserving initial source footprint.

This is not a fitted production calibration.  The TI021 workbook has arrival
summaries but no first-10-second concentration traces, so the screen only
tests premature threshold detections.  For each horizontal footprint W the
height is changed as H=H0*(W0/W)**2 to keep the initial prism volume fixed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from iterate_ti021_initial10s import run


def screen(output: Path, *, flat_to_flat_m: tuple[float, ...],
           reference_flat_to_flat_m: float = 14.0,
           reference_height_m: float = 13.0) -> dict:
    if reference_flat_to_flat_m <= 0.0 or reference_height_m <= 0.0:
        raise ValueError("reference dimensions must be positive")
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for width in flat_to_flat_m:
        if width <= 0.0:
            raise ValueError("candidate widths must be positive")
        height = reference_height_m * (reference_flat_to_flat_m / width) ** 2
        candidate = output.parent / f".candidate-volume-flat-{width:g}.json"
        report = run(candidate, only_case="inertial_vertical_and_radial",
                     initial_flat_to_flat_m=width,
                     initial_height_m=height)
        case = report["cases"][0]
        score = case["scores"]["0.1"]
        rows.append({
            "initial_cloud_flat_to_flat_m": float(width),
            "initial_cloud_height_m": float(height),
            "early_detected_sensor_count": score["early_detected_sensor_count"],
            "mean_censored_earliness_s": score[
                "mean_censored_earliness_s_over_known_arrivals"],
            "prearrival_threshold_excess_vol_percent_s": score[
                "prearrival_threshold_excess_vol_percent_s"],
            "sensor_1_predicted_detection_s": case["sensor_rows"][0][
                "predicted_arrival_s"],
            "sensor_4_predicted_detection_s": case["sensor_rows"][3][
                "predicted_arrival_s"],
            "integral_to_layer_time_s": case["integral_to_layer_time_s"],
        })
        for sidecar in (candidate, candidate.with_suffix(".npz")):
            if sidecar.exists():
                sidecar.unlink()
    ranked = sorted(rows, key=lambda row: (
        row["early_detected_sensor_count"],
        row["mean_censored_earliness_s"],
        row["prearrival_threshold_excess_vol_percent_s"],
    ))
    result = {
        "schema": "degali.ti021-initial10s-volume-preserving-footprint.v1",
        "status": "diagnostic_only_no_production_promotion",
        "physical_horizon_s": 10.0,
        "concentration_mape_available": False,
        "objective": "screen premature threshold detections while preserving initial volume",
        "reference_flat_to_flat_m": float(reference_flat_to_flat_m),
        "reference_height_m": float(reference_height_m),
        "height_rule": "H=H0*(W0/W)^2",
        "candidates": rows,
        "diagnostic_rank": ranked,
        "limitations": [
            "arrival summaries are censored timing diagnostics, not concentration MAPE",
            "a volume-preserving width change changes source aspect ratio and is not evidence of the true release geometry",
            "no candidate is promoted without measured release geometry, mass flow and injection momentum",
        ],
    }
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--flat-to-flat-m", type=float, nargs="+",
                        default=[4., 6., 8., 10., 12., 14.])
    args = parser.parse_args()
    result = screen(args.output, flat_to_flat_m=tuple(args.flat_to_flat_m))
    for row in result["diagnostic_rank"]:
        print(row["initial_cloud_flat_to_flat_m"],
              row["initial_cloud_height_m"],
              row["early_detected_sensor_count"],
              row["mean_censored_earliness_s"])
