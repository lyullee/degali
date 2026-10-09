"""Read-only FFI diagnostics to support a DEGALI development strategy.

No coefficient or source input is selected from observed concentrations.
Wind-history results are instantaneous steady-field replay, not transient
transport predictions. Generated numerical artifacts are research diagnostics.
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

from degali.addons.transient_receptor import (
    FixedReceptor, SteadyPlumeTable, WindHistory, replay_fixed_receptors,
)
from degali.addons.event_balanced_metrics import (
    EventMetricObservation,
    score_event_balanced_metrics,
)
from degali.validation import nearfield


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def fingerprint(path):
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def scores(observed, predicted):
    observed, predicted = np.asarray(observed), np.asarray(predicted)
    residual = predicted - observed
    result = {"n": len(observed), "mae_vol_pct_point": float(np.mean(abs(residual))),
              "rmse_vol_pct_point": float(np.sqrt(np.mean(residual ** 2))),
              "bias_vol_pct_point": float(np.mean(residual))}
    for threshold in (2.0, 4.0):
        o, p = observed >= threshold, predicted >= threshold
        result[f"threshold_{threshold:g}_vol_pct"] = {
            "tp": int(np.sum(o & p)), "fn": int(np.sum(o & ~p)),
            "fp": int(np.sum(~o & p)), "tn": int(np.sum(~o & ~p)),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slabx-root", type=Path, required=True)
    parser.add_argument("--paired-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--event-id", default="FFI-Test4",
                        help="release/event identifier used by event-balanced metrics")
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True)
    rows = read_csv(args.paired_csv)
    observed = np.asarray([float(row["observed_275s_mean_vol_pct"]) for row in rows])
    research_root = args.slabx_root / "논문화/PSEP_rebuild_2026-10-01/PSEP_ready_to_submit_2026-10-07/AUTHOR_REVIEW/PRIVATE_SUPPLEMENTARY_DATA"
    research_file = research_root / "Data_S5_Test4_wind_route_receptor_comparison.csv"
    research_by_sensor = {row["sensor"]: row for row in read_csv(research_file)}
    for row in rows:
        research = research_by_sensor[row["sensor"]]
        if not math.isclose(float(research["observed_mean_vol_pct"]), float(row["observed_275s_mean_vol_pct"]), abs_tol=1e-12):
            raise ValueError("observation mismatch")
    models = {
        "DEGALI_static_high_10m": np.asarray([float(row["DEGALI_steady_275s_setting_vol_pct"]) for row in rows]),
        "SLABx_native": np.asarray([float(row["SLABx_275s_mean_vol_pct"]) for row in rows]),
        "SLABx_LH2_retained_research": np.asarray([float(research_by_sensor[row["sensor"]]["v1165_model_vol_pct"]) for row in rows]),
    }
    strata = {}
    for name, mask in [
        *( (f"radius_{r:g}m", np.asarray([float(row["radius_m"]) == r for row in rows])) for r in (30.0, 50.0, 100.0)),
        ("axis", np.asarray([abs(float(row["y_crosswind_m"])) < 1e-6 for row in rows])),
        ("off_axis", np.asarray([abs(float(row["y_crosswind_m"])) >= 1e-6 for row in rows])),
        ("positive_observations", observed > 0.0), ("zero_observations", observed == 0.0),
    ]:
        strata[name] = {model: scores(observed[mask], values[mask]) for model, values in models.items()}
    static_absolute = abs(models["DEGALI_static_high_10m"] - observed)
    axial = np.asarray([abs(float(row["y_crosswind_m"])) < 1e-6 for row in rows])
    axial_error_share = float(sum(static_absolute[axial]) / sum(static_absolute))
    top_residuals = sorted([
        {"sensor": row["sensor"], "radius_m": float(row["radius_m"]), "height_m": float(row["height_m"]),
         "observed_vol_pct": float(observed[i]), "DEGALI_vol_pct": float(models["DEGALI_static_high_10m"][i]),
         "residual_vol_pct_point": float(models["DEGALI_static_high_10m"][i] - observed[i])}
        for i, row in enumerate(rows)
    ], key=lambda row: abs(row["residual_vol_pct_point"]), reverse=True)
    receptors = tuple(FixedReceptor(row["sensor"], float(row["x_downwind_m"]),
                                   -float(row["y_crosswind_m"]), float(row["height_m"])) for row in rows)
    trajectories = {}
    solve_log = []

    def trajectory(speed, ref_height):
        key = (speed, ref_height)
        if key not in trajectories:
            start = time.perf_counter()
            jet, initial = nearfield.hydrogen_jet(
                rate=0.828, diameter=0.0254, wind=speed, height=0.5,
                ambient_temperature=276.45, relative_humidity=75.0,
                storage_pressure_barg=2.12, roughness=0.03, stability="D",
                averaging=275.0, wind_reference_height=ref_height,
                corrections=True, ground_effect=False,
            )
            result = jet.run(initial, distmx=nearfield.STEP, smax=250.0)
            plume = nearfield.Trajectory(jet.th.table, result.rows)
            if not plume.ok:
                raise RuntimeError(f"trajectory failed at {key}")
            if any(plume.at(x) is None for x in (30.0, 50.0, 100.0)):
                raise RuntimeError(f"trajectory coverage failed at {key}")
            trajectories[key] = plume
            solve_log.append({"wind_speed_m_s": speed, "reference_height_m": ref_height,
                              "wall_seconds": time.perf_counter() - start})
        return trajectories[key]

    baseline_plume = trajectory(6.7, 10.0)
    baseline_reproduced = np.asarray([baseline_plume.concentration_at(float(row["x_downwind_m"]), float(row["y_crosswind_m"]), float(row["height_m"])) for row in rows])
    if not np.allclose(baseline_reproduced, models["DEGALI_static_high_10m"], rtol=1e-10, atol=1e-10):
        raise RuntimeError("static baseline did not reproduce")
    sections = {str(x): vars(baseline_plume.at(x)) for x in (30.0, 50.0, 100.0)}
    fingerprints = [fingerprint(args.paired_csv), fingerprint(research_file)]
    histories = {}
    replay_metadata = {}
    for channel, ref_height, mean_speed in (("high", 10.0, 6.7), ("low", 5.0, 5.0)):
        path = args.slabx_root / f"model-comparison/inputs/ffi_test4_{channel}_anemometer_digitized_5s.csv"
        bins = [row for row in read_csv(path) if float(row["start_s"]) >= 25.0 and float(row["end_s"]) <= 300.0]
        starts = np.asarray([float(row["start_s"]) for row in bins])
        ends = np.asarray([float(row["end_s"]) for row in bins])
        if starts[0] != 25.0 or ends[-1] != 300.0 or not np.allclose(starts[1:], ends[:-1]):
            raise ValueError("wind history does not cover contiguous 25-300 s window")
        durations = ends - starts
        if not np.allclose(durations, 5.0):
            raise ValueError("this diagnostic requires equal 5s digitization bins")
        times = (starts + ends) / 2.0
        speeds = np.asarray([float(row["speed_m_s"]) for row in bins])
        directions = np.asarray([float(row["from_bearing_deg"]) for row in bins])
        fingerprints.append(fingerprint(path))
        histories[channel] = (times, speeds, directions, ref_height)
        metadata = {"window_s": [25.0, 300.0], "bin_count": len(bins), "reference_height_m": ref_height,
                    "speed_min_max_m_s": [float(min(speeds)), float(max(speeds))],
                    "duration_weighted_speed_m_s": float(np.average(speeds, weights=durations)),
                    "source": "digitized meteorology, not digitized sensor concentrations"}
        # A prescribed 1 m/s numerical table, with a later 0.5 m/s check.
        lower = max(0.5, math.floor(float(min(speeds))))
        upper = math.ceil(float(max(speeds)))
        nodes = np.arange(lower, upper + 0.1, 1.0)
        samplers = []
        for speed in nodes:
            plume = trajectory(float(speed), ref_height)
            samplers.append(lambda x, y, z, plume=plume: plume.concentration_at(x, y, z) / 100.0)
        table = SteadyPlumeTable(nodes, samplers)
        static = trajectory(mean_speed, ref_height)
        models[f"DEGALI_static_{channel}_reported_height"] = np.asarray([
            static.concentration_at(float(row["x_downwind_m"]), float(row["y_crosswind_m"]), float(row["height_m"])) for row in rows])
        for variant in ("direction_only", "speed_only", "speed_direction"):
            used_speeds = np.full_like(speeds, mean_speed) if variant == "direction_only" else speeds
            used_directions = np.full_like(directions, 270.0) if variant == "speed_only" else directions
            if variant == "direction_only":
                sampler = lambda speed, x, y, z: static.concentration_at(x, y, z) / 100.0
            else:
                sampler = table
            traces = replay_fixed_receptors(WindHistory(times, used_speeds, used_directions), receptors, sampler)
            label = f"DEGALI_replay_{channel}_{variant}"
            models[label] = np.asarray([float(np.average(trace.true_mole_fraction, weights=durations)) * 100.0 for trace in traces])
            print(label, json.dumps(scores(observed, models[label])), flush=True)
        replay_metadata[channel] = metadata
    # Numerical table refinement for the high-channel combined replay.
    times, speeds, directions, ref_height = histories["high"]
    lower = max(0.5, math.floor(float(min(speeds))))
    upper = math.ceil(float(max(speeds)))
    refined_nodes = np.arange(lower, upper + 0.1, 0.5)
    refined_samplers = []
    for speed in refined_nodes:
        plume = trajectory(float(speed), ref_height)
        refined_samplers.append(lambda x, y, z, plume=plume: plume.concentration_at(x, y, z) / 100.0)
    refined_table = SteadyPlumeTable(refined_nodes, refined_samplers)
    traces = replay_fixed_receptors(WindHistory(times, speeds, directions), receptors, refined_table)
    refined = np.asarray([float(np.mean(trace.true_mole_fraction)) * 100.0 for trace in traces])
    coarse = models["DEGALI_replay_high_speed_direction"]
    models["DEGALI_replay_high_speed_direction_refined"] = refined
    refinement = {"coarse_speed_step_m_s": 1.0, "refined_speed_step_m_s": 0.5,
                  "max_sensor_delta_vol_pct_point": float(max(abs(coarse - refined))),
                  "mae_delta_vol_pct_point": scores(observed, refined)["mae_vol_pct_point"] - scores(observed, coarse)["mae_vol_pct_point"]}
    event_balanced_scores = {
        name: score_event_balanced_metrics(
            tuple(
                EventMetricObservation(
                    args.event_id, row["sensor"], float(observed[index]) / 100.0,
                    float(values[index]) / 100.0,
                )
                for index, row in enumerate(rows)
            ),
        ).as_dict()
        for name, values in models.items()
    }
    report = {
        "scope": f"{args.event_id} exploratory mechanism diagnostics; no independently held-out validation",
        "concentration_unit": "vol%; errors are volume-percentage points",
        "source_input": {"rate_kg_s": 0.828, "diameter_m": 0.0254, "release_height_m": 0.5,
                         "ambient_temperature_k": 276.45, "humidity_pct": 75.0, "P04_barg": 2.12,
                         "roughness_m": 0.03, "stability": "D", "averaging_setting_s": 275.0},
        "limitations": ["instantaneous steady replay has no travel-time memory or source transient",
                        "5s digitization cannot resolve 10Hz concentration peaks",
                        "sensor response not fitted; no lag assumed in true-field diagnostic",
                        "model source boundaries differ; scores do not qualify model selection",
                        "one release, 30 dependent receptors; no sensor-level confidence claim"],
        "fingerprints": fingerprints, "baseline_strata": strata,
        "axial_share_of_DEGALI_static_absolute_error": axial_error_share,
        "largest_static_residuals": top_residuals[:10], "static_sections": sections,
        "wind_replay_metadata": replay_metadata, "wind_table_refinement": refinement,
        "scores": {name: scores(observed, values) for name, values in models.items()},
        "event_balanced_scores": event_balanced_scores,
        "solve_log": solve_log,
        "development_priorities": [
            {"priority": 0, "work": "freeze native and research competitors; harmonize atmospheric flux boundary and sensor time operator"},
            {"priority": 1, "work": "unify causal wind/source history with conservative packet memory; isolate replay benefit first"},
            {"priority": 2, "work": "independent-energy and ground-contact state evolution, initialized before near-field handoff; no Test6-specific switch"},
            {"priority": 3, "work": "conserve flash/direct-vapour/droplet/pool mass and energy, validate source and substrate components separately"},
            {"priority": 4, "work": "replace prescribed obstacle mass routing with validated velocity and turbulence transport; validate tracer before cold H2"},
            {"priority": 5, "work": "cache thermodynamic tables and packet library after accuracy protocol passes"},
        ],
        "proposed_acceptance_targets": {
            "status": "planning targets, not demonstrated results or formal standards",
            "native_Test4_bridge_MAE_vol_pct_point": 0.45,
            "retained_research_Test4_MAE_vol_pct_point": 0.085,
            "retained_research_all210_MAE_vol_pct_point": 0.069,
            "retained_research_all210_RMSE_vol_pct_point": 0.106,
            "external_validation": "release-group split; baseline-relative metrics with release-level uncertainty, threshold FN no worse",
            "final_claim": "external transient detector and obstacle capability; FFI-only improvement insufficient",
        },
    }
    (args.output_dir / "analysis.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    paired_path = args.output_dir / "test4_replay_paired_predictions.csv"
    with paired_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["event_id", "sensor", "radius_m", "height_m",
                                                    "observed_vol_pct", *models])
        writer.writeheader()
        for i, row in enumerate(rows):
            writer.writerow({"event_id": args.event_id, "sensor": row["sensor"],
                             "radius_m": row["radius_m"], "height_m": row["height_m"],
                             "observed_vol_pct": observed[i], **{name: values[i] for name, values in models.items()}})
    print(json.dumps({"scores": report["scores"], "axial_error_share": axial_error_share,
                      "refinement": refinement, "output": str(args.output_dir.resolve())}, indent=2), flush=True)


if __name__ == "__main__":
    main()
