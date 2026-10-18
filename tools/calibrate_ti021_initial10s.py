"""Screen a short-window source-footprint correction for TI021.

This is deliberately a *diagnostic* calibration screen.  It uses only the
published arrival summaries to reject premature detections in the first ten
seconds; it does not fit concentration data (the workbook has no 0--10 s
concentration traces) and it never promotes the winning geometry.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from iterate_ti021_initial10s import run


def screen(output: Path, *, flat_to_flat_m: tuple[float, ...]) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for value in flat_to_flat_m:
        candidate_path = output.parent / f".candidate-flat-{value:g}.json"
        report = run(candidate_path, only_case="inertial_vertical_and_radial",
                     initial_flat_to_flat_m=value)
        case = report["cases"][0]
        score = case["scores"]["0.1"]
        rows.append({
            "initial_cloud_flat_to_flat_m": float(value),
            "initial_cloud_volume_m3": float(
                12 * (value / 2.0) ** 2 * math.tan(math.pi / 12)
                * 13.0
            ),
            "early_detected_sensor_count": score["early_detected_sensor_count"],
            "mean_censored_earliness_s": score[
                "mean_censored_earliness_s_over_known_arrivals"
            ],
            "prearrival_threshold_excess_vol_percent_s": score[
                "prearrival_threshold_excess_vol_percent_s"
            ],
            "sensor_1_predicted_detection_s": case["sensor_rows"][0][
                "predicted_arrival_s"
            ],
            "sensor_2_predicted_detection_s": case["sensor_rows"][1][
                "predicted_arrival_s"
            ],
            "sensor_3_predicted_detection_s": case["sensor_rows"][2][
                "predicted_arrival_s"
            ],
            "sensor_4_predicted_detection_s": case["sensor_rows"][3][
                "predicted_arrival_s"
            ],
            "integral_to_layer_time_s": case["integral_to_layer_time_s"],
        })
        for sidecar in (candidate_path, candidate_path.with_suffix(".npz")):
            if sidecar.exists():
                sidecar.unlink()
    ranked = sorted(rows, key=lambda row: (
        row["early_detected_sensor_count"],
        row["mean_censored_earliness_s"],
        row["prearrival_threshold_excess_vol_percent_s"],
    ))
    result = {
        "schema": "degali.ti021-initial10s-source-footprint-screen.v1",
        "status": "diagnostic_only_no_production_promotion",
        "physical_horizon_s": 10.0,
        "observed_operator": "first concentration reaching 5 percent of local peak",
        "model_operator": "first model trace crossing absolute 0.1 vol percent",
        "concentration_mape_available": False,
        "objective": "minimize premature absolute-threshold detections before 10 s",
        "interpretation": (
            "A smaller initial footprint suppresses the known early detections, "
            "but it also changes the released inventory.  It is therefore a "
            "source-state hypothesis, not an adopted fitted correction."
        ),
        "candidates": rows,
        "diagnostic_rank": ranked,
        "recommended_next_measurement": (
            "Provide the release-to-sensor clock and the first 10 s sensor traces "
            "before selecting a source footprint."
        ),
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
              row["early_detected_sensor_count"],
              row["mean_censored_earliness_s"])
