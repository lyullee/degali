"""Freeze a scoring policy, then score a completed PRESLHY E3.4 replay.

``--freeze-policy`` reads only the already sealed replay metadata.  ``--score``
is allowed only after a hash-bound full model output has completed and passed
its numerical gates; that action is the first step that opens the H2 values.
No time shift, source change, scalar-closure change, or parameter selection is
performed from the observations.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.validation.preslhy_e34 import (  # noqa: E402
    read_e34_dynamic_wind_hydrogen_history,
)


POLICY_SCHEMA = "degali.preslhy-e34-dynamic-component-scoring-policy.v1"
SCORE_SCHEMA = "degali.preslhy-e34-dynamic-component-score.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def _frozen_dependencies(directory: Path) -> tuple[Path, dict[str, Any], Path, dict[str, Any]]:
    freeze_path = directory / "model-configuration-freeze.json"
    freeze = _load_json(freeze_path)
    if freeze.get("schema") != (
        "degali.preslhy-e34-dynamic-component-configuration-freeze.v1"
    ):
        raise ValueError("unsupported model-configuration freeze schema")
    control = freeze.get("observation_control", {})
    if control.get("hydrogen_observation_values_opened") is not False:
        raise ValueError("configuration freeze is not observation blind")
    replay_path = Path(freeze["replay_input"]["path"]).resolve()
    if _sha256(replay_path) != freeze["replay_input"]["sha256"]:
        raise ValueError("dynamic replay input changed after model freeze")
    replay = _load_json(replay_path)
    return freeze_path, freeze, replay_path, replay


def freeze_scoring_policy(freeze_directory: str | Path) -> dict[str, Any]:
    """Write the observation-blind scoring policy before model scoring."""

    directory = Path(freeze_directory).resolve()
    policy_path = directory / "scoring-policy-freeze.json"
    if policy_path.exists():
        raise FileExistsError(policy_path)
    freeze_path, freeze, replay_path, replay = _frozen_dependencies(directory)
    sealed = replay["sealed_observation_operator"]
    policy = {
        "schema": POLICY_SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "configuration_freeze": {
            "path": str(freeze_path),
            "sha256": _sha256(freeze_path),
            "model_configuration_sha256": freeze["model_configuration_sha256"],
        },
        "replay_input": {"path": str(replay_path), "sha256": _sha256(replay_path)},
        "observation_control": {
            "hydrogen_observation_values_opened": False,
            "hydrogen_observation_values_used": False,
            "scoring_policy_frozen_before_observation_values_opened": True,
        },
        "policy": {
            "channels": sealed["usable_channels"],
            "excluded_channels": sealed["excluded_channels"],
            "observation_clock": sealed["physical_observation_clock"],
            "declared_sample_transport_delay_s": sealed[
                "declared_sample_transport_delay_s"
            ],
            "evaluation_interval_s": [0.0, float(replay["duration_s"])],
            "model_to_observation_interpolation": "linear",
            "concentration_unit": "vol_percent",
            "pointwise_metrics": [
                "mae", "rmse", "mean_bias", "pearson_r", "fac2",
            ],
            "log_metric_floor_vol_percent": 0.10,
            "peak_operator": "maximum_raw_10hz_sample",
            "dose_operator": "trapezoidal_integral_over_complete_common_clock",
            "arrival_threshold_vol_percent": 1.0,
            "arrival_normalization": (
                "absolute_error_divided_by_max_observed_arrival_1s"
            ),
            "flammability_threshold_vol_percent": 4.0,
            "four_vol_classification_primary": "pooled_channel_time_samples",
            "four_vol_classification_secondary": "sensor_event_detected_or_not",
            "no_fit_controls": [
                "no_additional_time_shift",
                "no_source_adjustment",
                "no_closure_selection",
                "no_sensor_response_fit",
                "no_best_window_selection",
            ],
        },
        "predeclared_accuracy_targets": {
            "sensor_inclusion": "observed_peak_at_least_1_vol_percent",
            "median_absolute_peak_error_vol_percent_point_max": 1.0,
            "median_absolute_relative_dose_error_max": 0.30,
            "median_normalized_arrival_time_error_max": 0.20,
            "pooled_four_vol_sensitivity_min": 0.90,
            "pooled_four_vol_specificity_min": 0.90,
        },
        "scope": {
            "dynamic_component_only": True,
            "obstacle_region_scoring_available": False,
            "obstacle_transport_model_qualification_allowed": False,
            "quantitative_facility_prediction_allowed": False,
        },
    }
    policy_path.write_text(
        json.dumps(policy, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return policy


def _arrival(time_s: np.ndarray, values: np.ndarray, threshold: float) -> float | None:
    indices = np.flatnonzero(values >= threshold)
    return None if not indices.size else float(time_s[indices[0]])


def _safe_pearson(observed: np.ndarray, model: np.ndarray) -> float | None:
    if observed.size < 2 or np.std(observed) == 0.0 or np.std(model) == 0.0:
        return None
    return float(np.corrcoef(observed, model)[0, 1])


def _binary_counts(observed: np.ndarray, model: np.ndarray, threshold: float) -> dict[str, Any]:
    obs = observed >= threshold
    pred = model >= threshold
    tp = int(np.count_nonzero(obs & pred))
    tn = int(np.count_nonzero(~obs & ~pred))
    fp = int(np.count_nonzero(~obs & pred))
    fn = int(np.count_nonzero(obs & ~pred))
    return {
        "true_positive": tp,
        "true_negative": tn,
        "false_positive": fp,
        "false_negative": fn,
        "sensitivity": None if tp + fn == 0 else tp / (tp + fn),
        "specificity": None if tn + fp == 0 else tn / (tn + fp),
    }


def _channel_score(
    time_s: np.ndarray,
    observed: np.ndarray,
    model: np.ndarray,
    *,
    log_floor: float,
    arrival_threshold: float,
    flammability_threshold: float,
) -> dict[str, Any]:
    delta = model - observed
    observed_peak = float(np.max(observed))
    model_peak = float(np.max(model))
    observed_dose = float(np.trapezoid(observed, time_s))
    model_dose = float(np.trapezoid(model, time_s))
    observed_arrival = _arrival(time_s, observed, arrival_threshold)
    model_arrival = _arrival(time_s, model, arrival_threshold)
    if observed_arrival is None or model_arrival is None:
        normalized_arrival_error = None
    else:
        normalized_arrival_error = abs(model_arrival - observed_arrival) / max(
            observed_arrival, 1.0,
        )
    observed_log = np.maximum(observed, log_floor)
    model_log = np.maximum(model, log_floor)
    ratio = model_log / observed_log
    return {
        "samples": int(time_s.size),
        "mae_vol_percent_point": float(np.mean(np.abs(delta))),
        "rmse_vol_percent_point": float(np.sqrt(np.mean(delta**2))),
        "mean_bias_vol_percent_point": float(np.mean(delta)),
        "pearson_r": _safe_pearson(observed, model),
        "fac2": float(np.mean((ratio >= 0.5) & (ratio <= 2.0))),
        "observed_peak_vol_percent": observed_peak,
        "model_peak_vol_percent": model_peak,
        "absolute_peak_error_vol_percent_point": abs(model_peak - observed_peak),
        "observed_dose_vol_percent_s": observed_dose,
        "model_dose_vol_percent_s": model_dose,
        "absolute_relative_dose_error": abs(model_dose - observed_dose) / max(
            abs(observed_dose), 1.0e-30,
        ),
        "observed_arrival_1vol_s": observed_arrival,
        "model_arrival_1vol_s": model_arrival,
        "normalized_arrival_time_error": normalized_arrival_error,
        "observed_four_vol_event": bool(np.any(observed >= flammability_threshold)),
        "model_four_vol_event": bool(np.any(model >= flammability_threshold)),
        "four_vol_sample_classification": _binary_counts(
            observed, model, flammability_threshold,
        ),
    }


def score_completed_component(freeze_directory: str | Path) -> dict[str, Any]:
    """Open H2 values only after the complete frozen model output exists."""

    directory = Path(freeze_directory).resolve()
    score_path = directory / "dynamic-component-score.json"
    arrays_path = directory / "dynamic-component-score-arrays.npz"
    plot_path = directory / "dynamic-component-score.png"
    if any(path.exists() for path in (score_path, arrays_path, plot_path)):
        raise FileExistsError("one or more score outputs already exist")
    freeze_path, freeze, replay_path, replay = _frozen_dependencies(directory)
    policy_path = directory / "scoring-policy-freeze.json"
    policy = _load_json(policy_path)
    if policy.get("schema") != POLICY_SCHEMA:
        raise ValueError("unsupported scoring-policy freeze schema")
    if policy["configuration_freeze"]["sha256"] != _sha256(freeze_path):
        raise ValueError("configuration freeze changed after scoring-policy freeze")
    if policy["replay_input"]["sha256"] != _sha256(replay_path):
        raise ValueError("replay input changed after scoring-policy freeze")

    full_path = directory / "full-run-result.json"
    full = _load_json(full_path)
    if full.get("schema") != "degali.preslhy-e34-dynamic-component-full-run.v1":
        raise ValueError("unsupported full-run schema")
    if not full.get("all_predeclared_numerical_gates_passed", False):
        raise ValueError("full model run failed one or more numerical gates")
    if not full.get("evidence_decision", {}).get(
        "sealed_observation_scoring_allowed", False,
    ):
        raise ValueError("full model run withholds sealed observation scoring")
    model_output_path = Path(full["model_output"]["path"]).resolve()
    if _sha256(model_output_path) != full["model_output"]["sha256"]:
        raise ValueError("full model output fingerprint changed")
    provenance_correction_path = (
        directory / "model-configuration-provenance-correction.json"
    )
    provenance_correction = None
    if provenance_correction_path.is_file():
        provenance_correction = _load_json(provenance_correction_path)
        if provenance_correction.get("configuration_freeze_sha256") != _sha256(
            freeze_path
        ):
            raise ValueError("provenance correction does not bind this freeze")
        if provenance_correction.get("numerical_configuration_changed") is not False:
            raise ValueError("provenance correction must not change model inputs")
        if provenance_correction.get(
            "hydrogen_observation_values_opened_before_correction"
        ) is not False:
            raise ValueError("provenance correction was not observation blind")

    workbook = Path(replay["workbook_path"]).resolve()
    if _sha256(workbook) != replay["workbook_sha256"]:
        raise ValueError("PRESLHY workbook fingerprint changed")
    dynamic = read_e34_dynamic_wind_hydrogen_history(workbook)
    start_s, end_s = replay["source_window_synchronized_clock_s"]
    selected = (
        (dynamic.hydrogen_sample_time_s >= start_s)
        & (dynamic.hydrogen_sample_time_s <= end_s)
    )
    observation_time = dynamic.hydrogen_sample_time_s[selected] - start_s
    if observation_time.size != replay["sealed_observation_operator"][
        "selected_row_count"
    ]:
        raise ValueError("selected H2 row count changed after replay freeze")
    with np.load(model_output_path) as model_data:
        model_time = np.asarray(model_data["time_s"], dtype=float)
        observed_arrays: dict[str, np.ndarray] = {}
        model_arrays: dict[str, np.ndarray] = {}
        channel_scores: dict[str, dict[str, Any]] = {}
        scoring = policy["policy"]
        for label in scoring["channels"]:
            observed = np.asarray(dynamic.hydrogen_vol_percent[label][selected], dtype=float)
            model_native = np.asarray(
                model_data[f"model_receptor_{label}_h2_vol_percent"], dtype=float,
            )
            model = np.interp(observation_time, model_time, model_native)
            observed_arrays[label] = observed
            model_arrays[label] = model
            channel_scores[label] = _channel_score(
                observation_time,
                observed,
                model,
                log_floor=float(scoring["log_metric_floor_vol_percent"]),
                arrival_threshold=float(scoring["arrival_threshold_vol_percent"]),
                flammability_threshold=float(
                    scoring["flammability_threshold_vol_percent"]
                ),
            )

    included = [
        label for label, metric in channel_scores.items()
        if metric["observed_peak_vol_percent"] >= 1.0
    ]
    peak_errors = [
        channel_scores[label]["absolute_peak_error_vol_percent_point"]
        for label in included
    ]
    dose_errors = [
        channel_scores[label]["absolute_relative_dose_error"] for label in included
    ]
    arrival_errors = [
        channel_scores[label]["normalized_arrival_time_error"] for label in included
    ]
    arrival_complete = bool(
        included and all(value is not None for value in arrival_errors)
    )
    pooled_observed = np.concatenate([observed_arrays[label] for label in scoring["channels"]])
    pooled_model = np.concatenate([model_arrays[label] for label in scoring["channels"]])
    pooled_four = _binary_counts(pooled_observed, pooled_model, 4.0)
    event_observed = np.asarray([
        channel_scores[label]["observed_four_vol_event"] for label in scoring["channels"]
    ])
    event_model = np.asarray([
        channel_scores[label]["model_four_vol_event"] for label in scoring["channels"]
    ])
    event_four = _binary_counts(event_observed.astype(float), event_model.astype(float), 0.5)
    targets = policy["predeclared_accuracy_targets"]
    aggregate = {
        "included_sensor_count": len(included),
        "included_sensors": included,
        "median_absolute_peak_error_vol_percent_point": (
            None if not peak_errors else float(np.median(peak_errors))
        ),
        "median_absolute_relative_dose_error": (
            None if not dose_errors else float(np.median(dose_errors))
        ),
        "median_normalized_arrival_time_error": (
            None if not arrival_complete else float(np.median(arrival_errors))
        ),
        "pooled_four_vol_sample_classification": pooled_four,
        "sensor_event_four_vol_classification": event_four,
    }
    target_gates = {
        "eligible_sensor_set_nonempty": bool(included),
        "median_absolute_peak_error": bool(
            peak_errors
            and aggregate["median_absolute_peak_error_vol_percent_point"]
            <= targets["median_absolute_peak_error_vol_percent_point_max"]
        ),
        "median_absolute_relative_dose_error": bool(
            dose_errors
            and aggregate["median_absolute_relative_dose_error"]
            <= targets["median_absolute_relative_dose_error_max"]
        ),
        "median_normalized_arrival_time_error": bool(
            arrival_complete
            and aggregate["median_normalized_arrival_time_error"]
            <= targets["median_normalized_arrival_time_error_max"]
        ),
        "pooled_four_vol_sensitivity": bool(
            pooled_four["sensitivity"] is not None
            and pooled_four["sensitivity"] >= targets["pooled_four_vol_sensitivity_min"]
        ),
        "pooled_four_vol_specificity": bool(
            pooled_four["specificity"] is not None
            and pooled_four["specificity"] >= targets["pooled_four_vol_specificity_min"]
        ),
    }
    np.savez_compressed(
        arrays_path,
        time_s=observation_time,
        **{
            f"observed_{label}_h2_vol_percent": values
            for label, values in observed_arrays.items()
        },
        **{
            f"model_{label}_h2_vol_percent": values
            for label, values in model_arrays.items()
        },
    )
    record = {
        "schema": SCORE_SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "configuration_freeze": {"path": str(freeze_path), "sha256": _sha256(freeze_path)},
        "scoring_policy_freeze": {"path": str(policy_path), "sha256": _sha256(policy_path)},
        "full_model_run": {"path": str(full_path), "sha256": _sha256(full_path)},
        "configuration_provenance_correction": (
            None
            if provenance_correction is None
            else {
                "path": str(provenance_correction_path),
                "sha256": _sha256(provenance_correction_path),
                "numerical_configuration_changed": False,
            }
        ),
        "score_arrays": {"path": str(arrays_path), "sha256": _sha256(arrays_path)},
        "observation_control": {
            "hydrogen_observation_values_opened_after_model_and_policy_freeze": True,
            "observation_used_for_parameter_selection": False,
        },
        "channel_scores": channel_scores,
        "aggregate_targets": aggregate,
        "target_gate_results": target_gates,
        "all_predeclared_dynamic_accuracy_targets_passed": bool(
            all(target_gates.values())
        ),
        "evidence_decision": {
            "status": (
                "accepted_dynamic_component_accuracy"
                if all(target_gates.values())
                else "withheld_dynamic_component_accuracy_gate_failed"
            ),
            "obstacle_transport_model_qualification_allowed": False,
            "complete_g3_dynamic_hydrogen_evidence": False,
            "quantitative_facility_prediction_allowed": False,
            "limitations": [
                "no_site_registered_obstacle_geometry",
                "liquid_mass_loss_is_a_conditional_vapour_source_proxy",
                "point_wind_block_mean_is_applied_over_the_complete_inlet",
                "event_specific_sensor_calibration_is_unavailable",
                "only_one_vertical_pool-centre_sensor_region_is_scored",
            ],
        },
        "review_plot": str(plot_path),
    }

    figure, axes = plt.subplots(len(scoring["channels"]), 1, figsize=(10, 6), sharex=True)
    axes = np.atleast_1d(axes)
    for axis, label in zip(axes, scoring["channels"]):
        axis.plot(observation_time, observed_arrays[label], label=f"{label} observed", lw=1.2)
        axis.plot(observation_time, model_arrays[label], label=f"{label} DEGALI", lw=1.2)
        axis.axhline(4.0, color="tab:red", ls="--", lw=0.9, label="4 vol%")
        axis.set_ylabel("H2 (vol%)")
        axis.grid(alpha=0.25)
        axis.legend(loc="upper right")
    axes[-1].set_xlabel("Event time (s)")
    figure.suptitle(
        "PRESLHY E3.4 Sand03 window 1 — frozen dynamic component score\n"
        "No obstacle qualification; no concentration-based fitting"
    )
    figure.tight_layout()
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    score_path.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze-directory", type=Path, required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--freeze-policy", action="store_true")
    action.add_argument("--score", action="store_true")
    arguments = parser.parse_args(argv)
    record = (
        freeze_scoring_policy(arguments.freeze_directory)
        if arguments.freeze_policy
        else score_completed_component(arguments.freeze_directory)
    )
    print(json.dumps(record, indent=2, allow_nan=False))
    if arguments.score:
        return 0 if record["all_predeclared_dynamic_accuracy_targets_passed"] else 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
