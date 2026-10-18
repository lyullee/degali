"""Independently score completed TI021 resolved-inner receptor histories."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-inner-transient-score.v1"


def _resolution_classification(cell_size_m: float) -> dict[str, object]:
    """Return resolution-specific labels without implying convergence."""

    dx = float(cell_size_m)
    if not math.isfinite(dx) or dx <= 0.0:
        raise ValueError("cell_size_m must be positive and finite")
    wall_height_cells = 5.0 / dx
    near_ground_alias = int(math.floor(0.4 / dx)) == int(
        math.floor(2.4 / dx)
    )
    if math.isclose(dx, 5.0):
        status = "scored_no_fit_inner_dx5_not_grid_converged"
        scope = (
            "coarse 5 m resolved-inner diagnostic score; not a promoted "
            "grid-converged TI021 validation result"
        )
        reason = (
            "the 5 m execution is conservative and stable but materially "
            "under-resolves the 5 m fence and near-ground sensor heights; "
            "the 2.5 m and 1.25 m comparisons remain required"
        )
    elif math.isclose(dx, 2.5):
        status = "scored_no_fit_inner_dx2p5_pending_finest_grid"
        scope = (
            "2.5 m resolved-inner no-fit diagnostic score; two wall-height "
            "cells improve geometry but the 1.25 m impact-window result is pending"
        )
        reason = (
            "the 2.5 m execution has two cells over the 5 m fence but still "
            "aliases the 0.4/2.4 m receptors; comparison with the frozen 1.25 m "
            "obstacle-impact member is required before numerical promotion"
        )
    else:
        token = str(dx).replace(".", "p")
        status = f"scored_no_fit_inner_dx{token}_pending_grid_audit"
        scope = (
            f"{dx:g} m resolved-inner no-fit diagnostic score; promotion is "
            "controlled by the separate grid/time audit"
        )
        reason = (
            "a completed no-fit score does not itself establish grid convergence "
            "or external model-form accuracy"
        )
    return {
        "status": status,
        "scope": scope,
        "reason": reason,
        "cell_size_m": dx,
        "wall_height_cells": wall_height_cells,
        "0p4m_and_2p4m_sensors_share_vertical_cell": near_ground_alias,
    }


def _run_completion_metadata(run_status: str) -> dict[str, bool]:
    """Separate terminal solver completion from experiment-horizon coverage."""

    accepted = {
        "accepted_completed_inner_transient",
        "accepted_inner_transient_qualification",
        "accepted_local_patch_transient_qualification",
    }
    if run_status not in accepted:
        raise ValueError("cannot score a run that has not completed")
    return {
        "scoring_opened_after_completed_run": True,
        "observation_horizon_complete": (
            run_status == "accepted_completed_inner_transient"
        ),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def trailing_boxcar_average(
    time_s: np.ndarray,
    traces: np.ndarray,
    *,
    window_s: float,
) -> np.ndarray:
    """Average piecewise-linear traces over a trailing physical-time window."""

    time = np.asarray(time_s, dtype=float)
    values = np.asarray(traces, dtype=float)
    if time.ndim != 1 or time.size < 2 or np.any(np.diff(time) <= 0.0):
        raise ValueError("time_s must be a strictly increasing vector")
    if values.shape[-1] != time.size or np.any(~np.isfinite(values)):
        raise ValueError("trace final dimension must match finite time_s")
    window = float(window_s)
    if not math.isfinite(window) or window <= 0.0:
        raise ValueError("window_s must be positive and finite")
    flat = values.reshape((-1, time.size))
    output = np.empty_like(flat)
    for row, trace in enumerate(flat):
        for index, end in enumerate(time):
            start = max(time[0], end - window)
            interior = time[(time > start) & (time < end)]
            sample_time = np.concatenate(([start], interior, [end]))
            sample_value = np.interp(sample_time, time, trace)
            span = end - start
            output[row, index] = (
                trace[0] if span <= 0.0
                else np.trapezoid(sample_value, sample_time) / span
            )
    return output.reshape(values.shape)


def _finite_correlation(x: np.ndarray, y: np.ndarray) -> dict[str, float | None]:
    valid = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(valid) < 3:
        return {"pearson_r": None, "spearman_rho": None}
    return {
        "pearson_r": float(pearsonr(x[valid], y[valid]).statistic),
        "spearman_rho": float(spearmanr(x[valid], y[valid]).statistic),
    }


def _ratio_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, object]:
    positive_observed = observed > 0.0
    paired_positive = positive_observed & (predicted > 0.0)
    ratio = predicted[positive_observed] / observed[positive_observed]
    positive_ratio = predicted[paired_positive] / observed[paired_positive]
    log_ratio = np.log(positive_ratio)
    return {
        "positive_observation_count": int(np.count_nonzero(positive_observed)),
        "paired_positive_count": int(np.count_nonzero(paired_positive)),
        "zero_prediction_for_positive_observation_count": int(np.count_nonzero(
            positive_observed & (predicted <= 0.0)
        )),
        "mean_absolute_error_percent_s": float(np.mean(np.abs(
            predicted[positive_observed] - observed[positive_observed]
        ))),
        "mean_absolute_percentage_error": float(np.mean(
            np.abs(ratio - 1.0)
        )),
        "normalized_mean_bias": float(
            np.sum(predicted[positive_observed] - observed[positive_observed])
            / np.sum(observed[positive_observed])
        ),
        "factor_of_two_fraction": float(np.mean(
            (ratio >= 0.5) & (ratio <= 2.0)
        )),
        "geometric_mean_ratio_paired_positive_only": (
            None if log_ratio.size == 0 else float(np.exp(np.mean(log_ratio)))
        ),
        "geometric_standard_deviation_paired_positive_only": (
            None if log_ratio.size < 2 else float(np.exp(np.std(log_ratio, ddof=1)))
        ),
        **_finite_correlation(
            observed[positive_observed], predicted[positive_observed],
        ),
    }


def score_run(
    workbook: Path,
    run_directory: Path,
    *,
    allow_horizon_limited: bool = False,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    """Open observations only after a completed no-fit run and score it."""

    workbook_path = workbook.resolve()
    run = run_directory.resolve()
    record_path = run / "ti021-inner-transient.json"
    trace_path = run / "ti021-inner-receptor-traces.npz"
    run_record = json.loads(record_path.read_text(encoding="utf-8"))
    run_status = run_record.get("status")
    accepted_statuses = {"accepted_completed_inner_transient"}
    if allow_horizon_limited:
        accepted_statuses.update({
            "accepted_inner_transient_qualification",
            "accepted_local_patch_transient_qualification",
        })
    if run_status not in accepted_statuses:
        raise ValueError("TI021 inner transient is not a completed accepted run")
    completion_metadata = _run_completion_metadata(run_status)
    if run_record.get("configuration", {}).get(
        "observation_values_used_to_configure_run"
    ) != []:
        raise ValueError("run is not observation-blind")
    if run_record.get("configuration", {}).get(
        "initial_binary_cloud_momentum"
    ) != "quiescent":
        raise ValueError("run did not use the report-derived stationary cloud")
    dataset = read_thorney_ti021_workbook(workbook_path)
    with np.load(trace_path) as bundle:
        time_s = np.asarray(bundle["time_s"], dtype=float)
        labels = np.asarray(bundle["labels"]).astype(str)
        vol_percent = np.asarray(
            bundle["concentration_vol_percent"], dtype=float,
        )
    metadata = run_record["configuration"]["receptors"]
    if labels.tolist() != [item["label"] for item in metadata]:
        raise ValueError("receptor labels do not match the frozen run record")
    workbook_indices = np.asarray([
        item["workbook_row_index_zero_based"] for item in metadata
    ], dtype=int)
    sensors = [dataset.dose_sensors[index] for index in workbook_indices]
    observed_dose = np.asarray([sensor.dose_percent_s for sensor in sensors])
    predicted_dose = np.trapezoid(vol_percent, time_s, axis=1)
    filtered = trailing_boxcar_average(time_s, vol_percent, window_s=0.6)
    predicted_peak = np.max(filtered, axis=1)
    predicted_peak_time = time_s[np.argmax(filtered, axis=1)]
    threshold = 0.1
    predicted_arrival = np.full(labels.size, np.nan)
    predicted_departure = np.full(labels.size, np.nan)
    for index, trace in enumerate(filtered):
        detected = np.flatnonzero(trace >= threshold)
        if detected.size:
            predicted_arrival[index] = time_s[detected[0]]
            predicted_departure[index] = time_s[detected[-1]]
    observed_arrival = np.asarray([
        np.nan if sensor.arrival_time_s is None else sensor.arrival_time_s
        for sensor in sensors
    ])
    observed_departure = np.asarray([
        np.nan if sensor.departure_time_s is None else sensor.departure_time_s
        for sensor in sensors
    ])

    arrival_available = np.isfinite(observed_arrival)
    arrival_detected = arrival_available & np.isfinite(predicted_arrival)
    departure_available = np.isfinite(observed_departure)
    departure_detected = departure_available & np.isfinite(predicted_departure)
    arc_indices: list[int] = []
    observed_arc_peak: list[float] = []
    observed_arc_peak_time: list[float] = []
    observed_arc_dose: list[float] = []
    arc_labels: list[str] = []
    for arc in dataset.arc_observations:
        match = [
            index for index, sensor in enumerate(sensors)
            if math.isclose(sensor.local_x_m, 0.0, abs_tol=1.0e-12)
            and math.isclose(sensor.local_y_m, arc.distance_m, abs_tol=1.0e-12)
            and math.isclose(sensor.height_m, arc.height_m, abs_tol=1.0e-12)
        ]
        if len(match) > 1:
            raise ValueError("could not uniquely map an arc observation to a sensor")
        if not match:
            # A local patch can legitimately omit the six far-field arc
            # receptors.  Keep the subset score usable, but never fabricate
            # an arc metric from a missing sensor.
            continue
        arc_indices.append(match[0])
        observed_arc_peak.append(arc.maximum_concentration_percent)
        observed_arc_peak_time.append(arc.peak_time_s)
        observed_arc_dose.append(arc.maximum_dose_percent_s)
        arc_labels.append(f"{arc.distance_m:g}m/{arc.height_m:g}m")
    arc_indices_array = np.asarray(arc_indices, dtype=int)
    observed_arc_peak_array = np.asarray(observed_arc_peak)
    observed_arc_peak_time_array = np.asarray(observed_arc_peak_time)
    observed_arc_dose_array = np.asarray(observed_arc_dose)
    predicted_arc_peak = predicted_peak[arc_indices_array]
    predicted_arc_peak_time = predicted_peak_time[arc_indices_array]
    predicted_arc_dose = predicted_dose[arc_indices_array]
    arc_peak_relative_error = (
        predicted_arc_peak / observed_arc_peak_array - 1.0
    )
    arc_dose_relative_error = (
        predicted_arc_dose / observed_arc_dose_array - 1.0
    )
    observed_four_percent = observed_arc_peak_array >= 4.0
    predicted_four_percent = predicted_arc_peak >= 4.0
    if len(arc_indices) == len(dataset.arc_observations):
        arc_coverage_status = "complete"
    elif arc_indices:
        arc_coverage_status = "partial_local_patch_coverage"
    else:
        arc_coverage_status = "unavailable_local_patch_coverage"

    def _mean_or_none(values: np.ndarray) -> float | None:
        return None if values.size == 0 else float(np.mean(values))

    dose_metrics = _ratio_metrics(observed_dose, predicted_dose)
    resolution = _resolution_classification(
        run_record["configuration"]["cell_size_m"]
    )
    horizon_limited = not completion_metadata["observation_horizon_complete"]
    score_status = resolution["status"]
    if horizon_limited:
        score_status = "scored_horizon_limited_" + score_status.removeprefix("scored_")
    score = {
        "schema": SCHEMA,
        "status": score_status,
        "run_status": run_status,
        "horizon_limited": horizon_limited,
        "scope": (
            resolution["scope"]
            if run_record["configuration"].get("coverage_scope")
            != "local_patch_sensor_subset"
            else resolution["scope"] + "; local sensor subset"
        ),
        "provenance": {
            "workbook_path": str(workbook_path),
            "workbook_sha256": _sha256(workbook_path),
            "run_record_path": str(record_path),
            "run_record_sha256": _sha256(record_path),
            "trace_path": str(trace_path),
            "trace_sha256": _sha256(trace_path),
        },
        "protocol": {
            "configuration_accessed_observation_values": [],
            **completion_metadata,
            "horizon_limited_warning": (
                "metrics are provisional because the run ends before the "
                "published observation horizon" if horizon_limited else None
            ),
            "concentration_unit": "vol%",
            "dose_unit": "vol% s",
            "sensor_response_average_s": 0.6,
            "arrival_and_departure_threshold_vol_percent": threshold,
            "dose_integration": (
                "trapezoidal over the available run horizon; not a complete "
                "0-600 s dose" if horizon_limited else
                "trapezoidal over the complete 0-600 s trace"
            ),
            "zero_observations_retained": True,
        },
        "dose_sensor_score": {
            "sensor_count": labels.size,
            "observed_positive_count": int(np.count_nonzero(observed_dose > 0.0)),
            "predicted_detected_count_at_0p1vol_percent": int(np.count_nonzero(
                predicted_peak >= threshold
            )),
            **dose_metrics,
        },
        "arrival_score": {
            "observed_available_count": int(np.count_nonzero(arrival_available)),
            "model_detected_among_available_count": int(np.count_nonzero(
                arrival_detected
            )),
            "missed_count": int(np.count_nonzero(
                arrival_available & ~np.isfinite(predicted_arrival)
            )),
            "mean_error_s": float(np.mean(
                predicted_arrival[arrival_detected]
                - observed_arrival[arrival_detected]
            )),
            "mean_absolute_error_s": float(np.mean(np.abs(
                predicted_arrival[arrival_detected]
                - observed_arrival[arrival_detected]
            ))),
        },
        "departure_score": {
            "observed_available_count": int(np.count_nonzero(departure_available)),
            "model_detected_among_available_count": int(np.count_nonzero(
                departure_detected
            )),
            "mean_error_s": float(np.mean(
                predicted_departure[departure_detected]
                - observed_departure[departure_detected]
            )),
            "mean_absolute_error_s": float(np.mean(np.abs(
                predicted_departure[departure_detected]
                - observed_departure[departure_detected]
            ))),
        },
        "six_arc_score": {
            "coverage_status": arc_coverage_status,
            "available_arc_count": len(arc_indices),
            "declared_arc_count": len(dataset.arc_observations),
            "labels": arc_labels,
            "observed_peak_vol_percent": observed_arc_peak_array.tolist(),
            "predicted_peak_0p6s_average_vol_percent": predicted_arc_peak.tolist(),
            "peak_relative_error": arc_peak_relative_error.tolist(),
            "peak_mean_absolute_percentage_error": _mean_or_none(
                np.abs(arc_peak_relative_error)
            ),
            "observed_peak_time_s": observed_arc_peak_time_array.tolist(),
            "predicted_peak_time_s": predicted_arc_peak_time.tolist(),
            "peak_time_mean_error_s": _mean_or_none(
                predicted_arc_peak_time - observed_arc_peak_time_array
            ),
            "peak_time_mean_absolute_error_s": _mean_or_none(np.abs(
                predicted_arc_peak_time - observed_arc_peak_time_array
            )),
            "observed_dose_percent_s": observed_arc_dose_array.tolist(),
            "predicted_dose_percent_s": predicted_arc_dose.tolist(),
            "dose_relative_error": arc_dose_relative_error.tolist(),
            "dose_mean_absolute_percentage_error": _mean_or_none(
                np.abs(arc_dose_relative_error)
            ),
            "four_vol_percent_confusion": {
                "true_positive": int(np.count_nonzero(
                    observed_four_percent & predicted_four_percent
                )),
                "true_negative": int(np.count_nonzero(
                    ~observed_four_percent & ~predicted_four_percent
                )),
                "false_positive": int(np.count_nonzero(
                    ~observed_four_percent & predicted_four_percent
                )),
                "false_negative": int(np.count_nonzero(
                    observed_four_percent & ~predicted_four_percent
                )),
                "note": (
                    "all six published TI021 arc peaks are below 4 vol%; this "
                    "tests true-negative behavior only, not 4 vol% sensitivity"
                ),
            },
        },
        "resolution_limit": {
            "cell_size_m": resolution["cell_size_m"],
            "wall_height_cells": resolution["wall_height_cells"],
            "0p4m_and_2p4m_sensors_share_vertical_cell": resolution[
                "0p4m_and_2p4m_sensors_share_vertical_cell"
            ],
            "grid_converged": False,
            "promotion_allowed": False,
        },
        "decision": {
            "result": "withheld_as_horizon_limited" if horizon_limited else "withheld_as_quantitative_validation",
            "reason": (
                "run horizon is shorter than the published TI021 observation "
                "window; do not interpret dose errors as final validation"
                if horizon_limited else resolution["reason"]
            ),
        },
    }
    arrays = {
        "time_s": time_s,
        "labels": labels,
        "observed_dose": observed_dose,
        "predicted_dose": predicted_dose,
        "observed_arrival": observed_arrival,
        "predicted_arrival": predicted_arrival,
        "filtered_vol_percent": filtered,
        "arc_indices": arc_indices_array,
        "observed_arc_peak": observed_arc_peak_array,
        "predicted_arc_peak": predicted_arc_peak,
        "observed_arc_peak_time": observed_arc_peak_time_array,
        "predicted_arc_peak_time": predicted_arc_peak_time,
        "observed_arc_dose": observed_arc_dose_array,
        "predicted_arc_dose": predicted_arc_dose,
        "arc_labels": np.asarray(arc_labels),
    }
    return score, arrays


def _plot(score: dict[str, object], arrays: dict[str, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), constrained_layout=True)
    observed_dose = arrays["observed_dose"]
    predicted_dose = arrays["predicted_dose"]
    positive = observed_dose > 0.0
    floor = 0.1
    x = np.maximum(observed_dose[positive], floor)
    y = np.maximum(predicted_dose[positive], floor)
    axes[0, 0].scatter(x, y, s=28, color="tab:blue", edgecolor="black", linewidth=0.3)
    bounds = [floor, max(np.max(x), np.max(y)) * 1.3]
    axes[0, 0].plot(bounds, bounds, color="black", label="1:1")
    axes[0, 0].plot(bounds, np.asarray(bounds) * 2.0, "--", color="0.5", label="factor 2")
    axes[0, 0].plot(bounds, np.asarray(bounds) * 0.5, "--", color="0.5")
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("observed dose (vol% s)")
    axes[0, 0].set_ylabel("predicted dose (vol% s)")
    axes[0, 0].set_title("60-sensor inner-domain dose")
    axes[0, 0].legend(fontsize=8)

    arc_count = len(arrays["arc_labels"])
    positions = np.arange(arc_count)
    width = 0.38
    if arc_count:
        axes[0, 1].bar(
            positions - width / 2, arrays["observed_arc_peak"], width,
            label="observed", color="black", alpha=0.75,
        )
        axes[0, 1].bar(
            positions + width / 2, arrays["predicted_arc_peak"], width,
            label=(
                "DEGALI dx="
                f"{score['resolution_limit']['cell_size_m']:g} m"
            ), color="tab:orange",
        )
        axes[0, 1].axhline(4.0, color="crimson", linestyle="--", label="4 vol%")
        axes[0, 1].set_xticks(positions, arrays["arc_labels"], rotation=25)
        axes[0, 1].set_ylabel("0.6 s averaged peak (vol%)")
        axes[0, 1].set_title("published centreline arc peaks")
        axes[0, 1].legend(fontsize=8)
    else:
        axes[0, 1].text(
            0.5, 0.5, "Arc sensors are outside this local patch",
            ha="center", va="center", transform=axes[0, 1].transAxes,
        )
        axes[0, 1].set_axis_off()

    for index, label in zip(arrays["arc_indices"], arrays["arc_labels"]):
        axes[1, 0].plot(
            arrays["time_s"], arrays["filtered_vol_percent"][index],
            label=label,
        )
    axes[1, 0].axhline(0.1, color="0.4", linestyle=":", label="0.1 vol%")
    axes[1, 0].axhline(4.0, color="crimson", linestyle="--", label="4 vol%")
    axes[1, 0].set_xlim(0.0, 180.0)
    axes[1, 0].set_xlabel("time (s)")
    axes[1, 0].set_ylabel("R12 concentration (vol%)")
    axes[1, 0].set_title(
        "model arc traces" if arc_count else "no arc traces in local patch"
    )
    if arc_count:
        axes[1, 0].legend(fontsize=7, ncol=2)

    available = np.isfinite(arrays["observed_arrival"])
    detected = available & np.isfinite(arrays["predicted_arrival"])
    axes[1, 1].scatter(
        arrays["observed_arrival"][detected],
        arrays["predicted_arrival"][detected],
        s=28, color="tab:green", edgecolor="black", linewidth=0.3,
    )
    limit = max(
        float(np.nanmax(arrays["observed_arrival"])),
        float(np.nanmax(arrays["predicted_arrival"][detected])),
    ) * 1.1
    axes[1, 1].plot((0.0, limit), (0.0, limit), color="black")
    axes[1, 1].set_xlim(0.0, limit)
    axes[1, 1].set_ylim(0.0, limit)
    axes[1, 1].set_xlabel("observed arrival at 0.1 vol% (s)")
    axes[1, 1].set_ylabel("predicted arrival (s)")
    axes[1, 1].set_title("arrival-time comparison")
    fig.suptitle(
        "TI021 no-fit inner score — "
        f"dx={score['resolution_limit']['cell_size_m']:g} m, "
        "quantitative promotion withheld"
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--run-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument(
        "--allow-horizon-limited",
        action="store_true",
        help=(
            "score a qualification run that ends before the published "
            "observation horizon; metrics remain provisional"
        ),
    )
    args = parser.parse_args()
    score, arrays = score_run(
        args.workbook,
        args.run_directory,
        allow_horizon_limited=args.allow_horizon_limited,
    )
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "ti021-inner-transient-score.json"
    plot_path = output / "ti021-inner-transient-score.png"
    json_path.write_text(
        json.dumps(score, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _plot(score, arrays, plot_path)
    print(json.dumps({
        "status": score["status"],
        "decision": score["decision"],
        "dose_sensor_score": score["dose_sensor_score"],
        "arrival_score": score["arrival_score"],
        "six_arc_score": score["six_arc_score"],
        "json": str(json_path),
        "plot": str(plot_path),
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
