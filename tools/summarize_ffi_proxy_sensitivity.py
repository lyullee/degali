"""Summarize corrected FFI O1 conditional proxy-width replays.

This is a diagnostic utility only.  It deliberately preserves the replay
status and refuses to turn proxy-width sensitivity into a facility score.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        default=[
            ROOT / "outputs/ffi-o1-proxy-0p4-corrected.json",
            ROOT / "outputs/ffi-o1-proxy-0p6-corrected.json",
            ROOT / "outputs/ffi-o1-proxy-0p8-corrected.json",
        ],
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs/ffi-o1-proxy-sensitivity-corrected.json",
    )
    args = parser.parse_args()
    rows = []
    for path in args.inputs:
        report = json.loads(path.read_text(encoding="utf-8"))
        diagnostic = report["diagnostic_sensor_error_not_validation_score"]
        rows.append(
            {
                "input": str(path),
                "proxy_width_m": report["proxy_width_m"],
                "runtime_s": report["runtime_s"],
                "maximum_predicted_vol_percent": report[
                    "maximum_predicted_vol_percent"
                ],
                "mae_vol_percent": diagnostic["mae_vol_percent"],
                "rmse_vol_percent": diagnostic["rmse_vol_percent"],
                "bias_vol_percent": diagnostic["bias_vol_percent"],
                "mape_percent_floor_0p1_vol_percent": diagnostic[
                    "mape_percent_floor_0p1_vol_percent"
                ],
                "observed_ge_4_mae_vol_percent": diagnostic[
                    "observed_ge_4_vol_percent"
                ]["mae_vol_percent"],
            }
        )
    rows.sort(key=lambda item: item["proxy_width_m"])
    mae_values = [item["mae_vol_percent"] for item in rows]
    output = {
        "schema": "degali.ffi-o1-proxy-sensitivity.v1",
        "status": "conditional_geometry_sensitivity_not_validation_score",
        "replay_count": len(rows),
        "rows": rows,
        "mae_range_vol_percent": [min(mae_values), max(mae_values)],
        "interpretation": (
            "Proxy width alone does not remove the observed underprediction; "
            "geometry, source/clock alignment, boundary condition, and sensor "
            "calibration gates remain unresolved."
        ),
        "promotion_allowed": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
