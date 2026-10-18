"""Diagnose, but never calibrate, the frozen PRESLHY E3.4 component score."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "degali.preslhy-e34-dynamic-component-diagnosis.v1"


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


def _safe_correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    if left.size < 2 or np.std(left) == 0.0 or np.std(right) == 0.0:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def diagnose_score(freeze_directory: str | Path) -> dict[str, Any]:
    directory = Path(freeze_directory).resolve()
    output_path = directory / "dynamic-component-diagnosis.json"
    plot_path = directory / "dynamic-component-diagnosis.png"
    if output_path.exists() or plot_path.exists():
        raise FileExistsError("diagnostic output already exists")
    score_path = directory / "dynamic-component-score.json"
    score = _load_json(score_path)
    if score.get("schema") != "degali.preslhy-e34-dynamic-component-score.v1":
        raise ValueError("unsupported dynamic component score schema")
    arrays_path = Path(score["score_arrays"]["path"]).resolve()
    if _sha256(arrays_path) != score["score_arrays"]["sha256"]:
        raise ValueError("score-array fingerprint changed")
    freeze_path = Path(score["configuration_freeze"]["path"]).resolve()
    freeze = _load_json(freeze_path)
    replay_path = Path(freeze["replay_input"]["path"]).resolve()
    replay = _load_json(replay_path)

    history = replay["inlet_velocity_history"]
    inlet_time = np.asarray(history["time_s"], dtype=float)
    inlet_u = np.asarray(history["u_m_s"], dtype=float)
    with np.load(arrays_path) as arrays:
        time_s = np.asarray(arrays["time_s"], dtype=float)
        sample_interval_s = float(np.median(np.diff(time_s)))
        index = np.searchsorted(inlet_time, time_s, side="right") - 1
        index = np.clip(index, 0, inlet_u.size - 1)
        applied_u = inlet_u[index]
        channels: dict[str, dict[str, Any]] = {}
        observed_series: dict[str, np.ndarray] = {}
        model_series: dict[str, np.ndarray] = {}
        for label in replay["sealed_observation_operator"]["usable_channels"]:
            observed = np.asarray(
                arrays[f"observed_{label}_h2_vol_percent"], dtype=float,
            )
            model = np.asarray(
                arrays[f"model_{label}_h2_vol_percent"], dtype=float,
            )
            observed_series[label] = observed
            model_series[label] = model
            observed_peak_index = int(np.argmax(observed))
            model_peak_index = int(np.argmax(model))
            low_wind = applied_u <= np.quantile(applied_u, 0.25)
            channels[label] = {
                "observed_initial_vol_percent": float(observed[0]),
                "model_initial_vol_percent": float(model[0]),
                "observed_peak_time_s": float(time_s[observed_peak_index]),
                "model_peak_time_s": float(time_s[model_peak_index]),
                "inlet_u_at_model_peak_m_s": float(applied_u[model_peak_index]),
                "observed_mean_vol_percent": float(np.mean(observed)),
                "model_mean_vol_percent": float(np.mean(model)),
                "observed_four_vol_sample_duration_s": float(
                    np.count_nonzero(observed >= 4.0) * sample_interval_s
                ),
                "model_four_vol_sample_duration_s": float(
                    np.count_nonzero(model >= 4.0) * sample_interval_s
                ),
                "model_mean_lowest_wind_quartile_vol_percent": float(
                    np.mean(model[low_wind])
                ),
                "model_mean_other_wind_quartiles_vol_percent": float(
                    np.mean(model[~low_wind])
                ),
                "observed_mean_lowest_wind_quartile_vol_percent": float(
                    np.mean(observed[low_wind])
                ),
                "observed_mean_other_wind_quartiles_vol_percent": float(
                    np.mean(observed[~low_wind])
                ),
                "model_correlation_with_inverse_inlet_u": _safe_correlation(
                    model, 1.0 / applied_u,
                ),
                "observed_correlation_with_inverse_inlet_u": _safe_correlation(
                    observed, 1.0 / applied_u,
                ),
            }

    initial_condition_mismatch = any(
        values["observed_initial_vol_percent"] >= 1.0
        and values["model_initial_vol_percent"] < 0.01
        for values in channels.values()
    )
    model_peak_indices = [
        int(np.argmax(model_series[label])) for label in model_series
    ]
    model_peaks_in_lowest_wind_quartile = all(
        applied_u[index_value] <= np.quantile(applied_u, 0.25)
        for index_value in model_peak_indices
    )
    record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "score": {"path": str(score_path), "sha256": _sha256(score_path)},
        "score_arrays": {"path": str(arrays_path), "sha256": _sha256(arrays_path)},
        "diagnostic_controls": {
            "parameter_adjustment_performed": False,
            "alternate_model_selected": False,
            "time_shift_fitted": False,
            "primary_frozen_score_replaced": False,
        },
        "channels": channels,
        "cross_channel_findings": {
            "initial_condition_mismatch_observed": initial_condition_mismatch,
            "model_peaks_in_lowest_inlet_wind_quartile": (
                model_peaks_in_lowest_wind_quartile
            ),
            "source_proxy_mass_relative_difference": replay[
                "source_schedule_mass_relative_difference"
            ],
            "primary_median_absolute_relative_dose_error": score[
                "aggregate_targets"
            ]["median_absolute_relative_dose_error"],
            "initial_condition_interpretation": (
                "The post-fill decay window begins after release activity during "
                "filling. A zero-concentration numerical initial field therefore "
                "does not represent the observed nonzero plume at event time zero."
            ),
            "low_wind_interpretation": (
                "Both model sensor peaks occur during the lowest inlet-U quartile. "
                "Applying one 10 s point-wind block mean over the whole inlet can "
                "turn a local lull into domain-wide accumulation."
            ),
            "source_magnitude_interpretation": (
                "The frozen source-proxy mass discrepancy is much smaller than the "
                "observed median dose error; it cannot by itself establish the "
                "dominant cause, and no source rescaling is authorized."
            ),
        },
        "next_evidence_actions": [
            "replay a source-only prehistory beginning before the post-fill scoring window and retain this run as the failed primary zero-initial-field case",
            "obtain spatial wind or CFD inflow evidence before replacing the whole-inlet point-history assumption",
            "run source-flux uncertainty as a declared bound rather than fitting mass rate to H2 observations",
            "do not use E3.4 to qualify obstacle transport because site-registered obstacles are unavailable",
        ],
        "status": "diagnosed_failed_primary_without_calibration",
        "review_plot": str(plot_path),
    }

    figure, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for label in observed_series:
        axes[0].plot(
            time_s, model_series[label] - observed_series[label], label=label,
        )
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].set_ylabel("Model - obs.\n(vol%p)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].step(time_s, applied_u, where="post", color="tab:green")
    axes[1].set_ylabel("Applied inlet U\n(m/s)")
    axes[1].grid(alpha=0.25)
    for label in observed_series:
        axes[2].plot(time_s, observed_series[label], label=f"{label} observed")
        axes[2].plot(time_s, model_series[label], label=f"{label} model", alpha=0.8)
    axes[2].axhline(4.0, color="tab:red", ls="--", lw=0.9)
    axes[2].set_ylabel("H2 (vol%)")
    axes[2].set_xlabel("Event time (s)")
    axes[2].grid(alpha=0.25)
    axes[2].legend(ncol=2)
    figure.suptitle(
        "PRESLHY E3.4 frozen-score diagnosis — failed primary retained, no fitting"
    )
    figure.tight_layout()
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    output_path.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze-directory", type=Path, required=True)
    arguments = parser.parse_args(argv)
    record = diagnose_score(arguments.freeze_directory)
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
