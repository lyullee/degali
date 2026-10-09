"""Conditional FFI Test 4 fixed-sensor comparison with an external SLABx run.

Bulk source, weather, geometry and nominal averaging setting are matched to
the separately archived SLABx Test 4 275 s mean calculation. DEGALI remains
quasi-steady; its 275 s parameter is not a synchronized time-series operator.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from degali.validation import nearfield
from degali.addons.field_comparison import (
    FieldComparisonBasis,
    compare_field_model_sensor_sets,
    field_model_decision_impact,
    field_model_decision_impact_record,
    field_model_sensor_set_from_csv,
)


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comparison_record(comparison):
    """Serialise the gate result without relabelling a conditional result."""
    decision_impact = field_model_decision_impact(comparison)
    return {
        "scope": "model-form comparison at fixed sensors; not a common-source validation score",
        "applicability": {
            "status": comparison.applicability.status,
            "reasons": list(comparison.applicability.reasons),
            "warnings": list(comparison.applicability.warnings),
            "uncertainty_complete": comparison.applicability.uncertainty_complete,
        },
        "threshold": {
            "mole_fraction": comparison.threshold_mole_fraction,
            "volume_percent": comparison.threshold_mole_fraction * 100.0,
        },
        "mean_absolute_difference": {
            "mole_fraction": comparison.mean_absolute_difference_mole_fraction,
            "volume_percent_point": (
                None if comparison.mean_absolute_difference_mole_fraction is None
                else comparison.mean_absolute_difference_mole_fraction * 100.0
            ),
        },
        "classification_disagreement_count": comparison.classification_disagreement_count,
        "decision_impact": field_model_decision_impact_record(decision_impact),
        "models": {
            "left": {
                "model_id": comparison.left.model_id,
                "temporal_mode": comparison.left.temporal_mode,
                "basis": comparison.left.basis.__dict__,
                "csv_provenance": (
                    None if comparison.left.csv_provenance is None
                    else comparison.left.csv_provenance.__dict__
                ),
            },
            "right": {
                "model_id": comparison.right.model_id,
                "temporal_mode": comparison.right.temporal_mode,
                "basis": comparison.right.basis.__dict__,
                "csv_provenance": (
                    None if comparison.right.csv_provenance is None
                    else comparison.right.csv_provenance.__dict__
                ),
            },
        },
        "sensor_rows": [
            {
                "sensor_id": row.sensor_id,
                "position_m": row.position_m,
                "left_mole_fraction": row.left_mole_fraction,
                "right_mole_fraction": row.right_mole_fraction,
                "difference_mole_fraction": row.difference_mole_fraction,
                "ratio_right_to_left": row.ratio_right_to_left,
                "left_above_threshold": row.left_above_threshold,
                "right_above_threshold": row.right_above_threshold,
            }
            for row in comparison.rows
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slabx-sensor-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(args.output_dir)
    rows = list(csv.DictReader(args.slabx_sensor_csv.open(encoding="utf-8-sig")))
    if len(rows) != 30 or {int(row["test"]) for row in rows} != {4}:
        raise ValueError("expected the accepted 30-sensor FFI Test 4 SLABx file")
    if {float(row["averaging_time_s"]) for row in rows} != {275.0}:
        raise ValueError("SLABx receptor window must be 275 s")
    start = time.perf_counter()
    jet, initial = nearfield.hydrogen_jet(
        rate=0.828, diameter=0.0254, wind=6.7,
        height=0.5, ambient_temperature=273.15 + 3.3,
        relative_humidity=75.0, storage_pressure_barg=2.12,
        roughness=0.03, stability="D", averaging=275.0,
        wind_reference_height=10.0, corrections=True, ground_effect=False,
    )
    trajectory = nearfield.Trajectory(
        jet.th.table, jet.run(initial, distmx=nearfield.STEP, smax=250.0).rows
    )
    if not trajectory.ok:
        raise RuntimeError("DEGALI Test 4 integration failed")
    paired = []
    for row in rows:
        x, y, z = (float(row[name]) for name in
                   ("x_downwind_m", "y_crosswind_m", "height_m"))
        prediction = trajectory.concentration_at(x, y, z)
        paired.append({
            "sensor": row["sensor"], "radius_m": float(row["radius_m"]),
            "x_downwind_m": x, "y_crosswind_m": y, "height_m": z,
            "averaging_time_s": 275.0,
            "observed_275s_mean_vol_pct": float(row["observed_time_mean_vol_pct"]),
            "SLABx_275s_mean_vol_pct": float(row["SLABx_time_mean_vol_pct"]),
            "DEGALI_steady_275s_setting_vol_pct": prediction,
            "SLABx_LFL": int(float(row["SLABx_time_mean_vol_pct"]) >= 4.0),
            "DEGALI_LFL": int(prediction >= 4.0),
            "observed_LFL": int(float(row["observed_time_mean_vol_pct"]) >= 4.0),
        })
    runtime = time.perf_counter() - start
    summary = {}
    obs = np.array([row["observed_275s_mean_vol_pct"] for row in paired])
    for name, field in (("SLABx", "SLABx_275s_mean_vol_pct"),
                        ("DEGALI", "DEGALI_steady_275s_setting_vol_pct")):
        pred = np.array([row[field] for row in paired])
        observed_positive = obs >= 4.0
        predicted_positive = pred >= 4.0
        summary[name] = {
            "n_sensors": len(rows), "bias_vol_pct_point": float(np.mean(pred - obs)),
            "mae_vol_pct_point": float(np.mean(abs(pred - obs))),
            "rmse_vol_pct_point": float(math.sqrt(np.mean((pred - obs) ** 2))),
            "true_positive": int(np.sum(observed_positive & predicted_positive)),
            "false_negative": int(np.sum(observed_positive & ~predicted_positive)),
            "false_positive": int(np.sum(~observed_positive & predicted_positive)),
        }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "test4_slabx_degali_paired_sensors.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(paired[0]))
        writer.writeheader()
        writer.writerows(paired)
    common_weather = "ffi-test4:wind-6.7m-s-at-10m:neutral:roughness-0.03m"
    common_sensors = "ffi-test4:30-fixed-sensor-positions"
    degali_set = field_model_sensor_set_from_csv(
        out,
        model_id="DEGALI-corrected-quasi-steady",
        temporal_mode="steady",
        basis=FieldComparisonBasis(
            "ffi-test4:degali-native-two-phase-source",
            common_weather, common_sensors,
            "degali-steady-averaging-setting-275s", 275.0,
        ),
        prediction_column="DEGALI_steady_275s_setting_vol_pct",
        concentration_unit="volume_percent",
    )
    slabx_set = field_model_sensor_set_from_csv(
        out,
        model_id="SLABx-Test4-external",
        temporal_mode="transient",
        basis=FieldComparisonBasis(
            "ffi-test4:slabx-native-two-phase-source",
            common_weather, common_sensors,
            "slabx-time-mean-275s", 275.0,
        ),
        prediction_column="SLABx_275s_mean_vol_pct",
        concentration_unit="volume_percent",
    )
    comparison = compare_field_model_sensor_sets(degali_set, slabx_set)
    comparison_path = args.output_dir / "test4_field_model_comparison.json"
    comparison_path.write_text(
        json.dumps(comparison_record(comparison), indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "scope": "conditional matched bulk inputs and fixed receptor positions",
        "source_file": str(args.slabx_sensor_csv.resolve()),
        "source_sha256": sha(args.slabx_sensor_csv),
        "matched_inputs": {"rate_kg_s": 0.828, "orifice_m": 0.0254,
                           "wind_m_s_at_10m": 6.7, "ambient_C": 3.3, "RH_pct": 75,
                           "roughness_m": 0.03, "stability": "D", "release_height_m": 0.5,
                           "averaging_setting_s": 275},
        "limits": ["SLABx uses its own two-phase source mapping; DEGALI uses its own flash source",
                   "DEGALI 275 s parameter adjusts steady meteorological closure and is not a synchronized 275 s time-series mean",
                   "one trial only; neither branch was fitted here"],
        "DEGALI_one_run_wall_time_s": runtime,
        "scores": summary,
        "field_model_comparison": {
            "file": comparison_path.name,
            "status": comparison.applicability.status,
            "warnings": list(comparison.applicability.warnings),
            "decision_impact_status": field_model_decision_impact(comparison).status,
        },
    }
    (args.output_dir / "test4_external_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
