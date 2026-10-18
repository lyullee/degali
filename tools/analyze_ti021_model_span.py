"""Audit the best global convex combination already spanned by TI021 runs.

This is a diagnostic, not a promoted calibration.  One non-negative weight is
assigned to each whole physical-model run and the weights sum to one.  No
sensor-specific multiplier or offset is allowed.  The result answers whether
the existing model family contains enough distinct spatial structure to reach
the requested dose-MAPE target before another expensive solver branch is
implemented.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outputs", type=Path, default=ROOT / "outputs")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--group-nondegradation", action="store_true")
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    rows = baseline["sensor_rows"]
    labels = [row["label"] for row in rows]
    observed = np.asarray(
        [row["observed_dose_vol_percent_s"] for row in rows], dtype=float
    )
    positive = observed > 0.0
    observed_positive = observed[positive]

    candidates: list[dict] = []
    seen: set[str] = set()
    baseline_observation_mode = baseline.get(
        "observation_species_mode", "r12_component"
    )
    baseline_gas = baseline.get("gas")
    skipped = {
        "unreadable": 0, "schema": 0, "duration": 0,
        "observable": 0, "duplicate": 0,
    }
    for path in sorted(args.outputs.rglob("*.json")):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            skipped["unreadable"] += 1
            continue
        sensor_rows = record.get("sensor_rows")
        if (
            not isinstance(sensor_rows, list)
            or len(sensor_rows) != len(rows)
            or [item.get("label") for item in sensor_rows] != labels
        ):
            skipped["schema"] += 1
            continue
        if abs(float(record.get("duration_s", 600.0)) - 600.0) > 1.0e-9:
            skipped["duration"] += 1
            continue
        if (record.get("observation_species_mode", "r12_component")
                != baseline_observation_mode
                or (baseline_gas is not None and record.get("gas") != baseline_gas)):
            skipped["observable"] += 1
            continue
        candidate_observed = np.asarray(
            [item.get("observed_dose_vol_percent_s", np.nan) for item in sensor_rows]
        )
        prediction = np.asarray(
            [item.get("predicted_dose_vol_percent_s", np.nan) for item in sensor_rows]
        )
        if (
            not np.all(np.isfinite(prediction))
            or not np.allclose(candidate_observed, observed, rtol=0.0, atol=1.0e-9)
        ):
            skipped["schema"] += 1
            continue
        digest = hashlib.sha256(prediction.tobytes()).hexdigest()
        if digest in seen:
            skipped["duplicate"] += 1
            continue
        seen.add(digest)
        relative = prediction[positive] / observed_positive
        candidates.append(
            {
                "path": str(path.resolve()),
                "prediction": prediction[positive],
                "individual_mape_percent": float(np.mean(np.abs(relative - 1.0)) * 100.0),
                "reported_parameters": {
                    key: record.get(key)
                    for key in (
                        "wind_angle_deg",
                        "wind_speed_m_s",
                        "layer_wind_speed_m_s",
                        "horizontal_diffusivity_m2_s",
                        "crosswind_diffusivity_factor",
                        "vertical_tail_diffusivity_m2_s",
                        "cavity_shelter_strength",
                        "wake_deficit_strength",
                        "wall_overflow_capture_fraction",
                        "wall_overflow_return_rate_s_inv",
                    )
                    if key in record
                },
            }
        )
    if not candidates:
        raise RuntimeError("no compatible 600 s TI021 result files were found")

    # R[j, i] is the relative prediction at sensor j from complete run i.
    relative_matrix = np.column_stack(
        [item["prediction"] / observed_positive for item in candidates]
    )
    sensor_count, model_count = relative_matrix.shape
    # Variables are model weights followed by per-sensor absolute residuals.
    objective = np.r_[np.zeros(model_count), np.full(sensor_count, 1.0 / sensor_count)]
    upper = np.block(
        [
            [relative_matrix, -np.eye(sensor_count)],
            [-relative_matrix, -np.eye(sensor_count)],
        ]
    )
    upper_rhs = np.r_[np.ones(sensor_count), -np.ones(sensor_count)]
    constrained_groups: dict[str, list[int]] = {}
    if args.group_nondegradation:
        dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
        positive_indices = np.flatnonzero(positive)
        sensors = [dataset.dose_sensors[index] for index in positive_indices]
        constrained_groups = {
            "height_0.4_m": [
                i for i,sensor in enumerate(sensors)
                if abs(sensor.height_m-0.4)<1.0e-9
            ],
            "height_2.4_m": [
                i for i,sensor in enumerate(sensors)
                if abs(sensor.height_m-2.4)<1.0e-9
            ],
            "height_4.4_m": [
                i for i,sensor in enumerate(sensors)
                if abs(sensor.height_m-4.4)<1.0e-9
            ],
            "height_6.4_m": [
                i for i,sensor in enumerate(sensors)
                if abs(sensor.height_m-6.4)<1.0e-9
            ],
            "downwind_lt_150_m": [
                i for i,sensor in enumerate(sensors)
                if sensor.local_y_m<150.0
            ],
            "downwind_150_to_300_m": [
                i for i,sensor in enumerate(sensors)
                if 150.0<=sensor.local_y_m<300.0
            ],
            "downwind_ge_300_m": [
                i for i,sensor in enumerate(sensors)
                if sensor.local_y_m>=300.0
            ],
        }
        baseline_relative=(
            np.asarray([
                row["predicted_dose_vol_percent_s"] for row in rows
            ])[positive]/observed_positive)
        group_rows=[]
        group_rhs=[]
        for indices in constrained_groups.values():
            row=np.zeros(model_count+sensor_count)
            row[model_count+np.asarray(indices,dtype=int)]=1.0/len(indices)
            group_rows.append(row)
            group_rhs.append(float(np.mean(
                np.abs(baseline_relative[indices]-1.0))))
        upper=np.vstack((upper,np.asarray(group_rows)))
        upper_rhs=np.r_[upper_rhs,np.asarray(group_rhs)]
    equality = np.zeros((1, model_count + sensor_count))
    equality[0, :model_count] = 1.0
    solution = linprog(
        objective,
        A_ub=upper,
        b_ub=upper_rhs,
        A_eq=equality,
        b_eq=np.array([1.0]),
        bounds=[(0.0, None)] * (model_count + sensor_count),
        method="highs",
    )
    if not solution.success:
        raise RuntimeError(f"convex ensemble optimization failed: {solution.message}")
    weights = solution.x[:model_count]
    ensemble_relative = relative_matrix @ weights
    active = np.flatnonzero(weights > 1.0e-7)
    active_rows = sorted(
        (
            {
                "weight": float(weights[index]),
                "path": candidates[index]["path"],
                "individual_mape_percent": candidates[index]["individual_mape_percent"],
                "reported_parameters": candidates[index]["reported_parameters"],
            }
            for index in active
        ),
        key=lambda item: item["weight"],
        reverse=True,
    )
    best = min(candidates, key=lambda item: item["individual_mape_percent"])
    report = {
        "schema": "degali.ti021-model-span-audit.v1",
        "purpose": "global physical-model-family diagnostic; not a promoted fit",
        "constraints": {
            "non_negative_whole_run_weights": True,
            "weights_sum_to_one": True,
            "sensor_specific_multiplier_or_offset": False,
            "metric": "600 s positive-observed dose MAPE",
            "height_and_distance_group_nondegradation":bool(
                args.group_nondegradation),
        },
        "compatible_unique_run_count": model_count,
        "positive_observed_sensor_count": sensor_count,
        "observation_species_mode": baseline_observation_mode,
        "gas": baseline_gas,
        "skipped": skipped,
        "constrained_groups":list(constrained_groups),
        "best_individual": {
            key: value for key, value in best.items() if key != "prediction"
        },
        "convex_ensemble": {
            "mape_percent": float(np.mean(np.abs(ensemble_relative - 1.0)) * 100.0),
            "normalized_mean_bias": float(np.mean(ensemble_relative - 1.0)),
            "maximum_absolute_percentage_error": float(
                np.max(np.abs(ensemble_relative - 1.0)) * 100.0
            ),
            "active_member_count": int(active.size),
            "members": active_rows,
            "target_10_percent_reached": bool(
                np.mean(np.abs(ensemble_relative - 1.0)) <= 0.10
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
