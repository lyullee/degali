"""Compare AIJ Case H scores over matched short and experiment-length windows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SHORT = ROOT / (
    "outputs/aij-case-h-tracer-fv-momentum-sgs-12s-dx2p5cm-dt1ms-"
    "2026-10-12-v1"
)
DEFAULT_FULL = ROOT / (
    "outputs/aij-case-h-tracer-fv-source-sigma-25mm-120s-full-window-"
    "2026-10-13-v1"
)
_ALLOWED_CONFIGURATION_CHANGES = {
    "scalar_duration_s",
    "scalar_mean_window_s",
    "integrated_released_species_mass_kg",
}


def _load(directory: Path) -> tuple[dict, dict[int, tuple[float, float]]]:
    record = json.loads(
        (directory / "aij-case-h-tracer-score.json").read_text(encoding="utf-8")
    )
    with (directory / "aij-case-h-tracer-point-scores.csv").open(
        "r", encoding="utf-8-sig", newline="",
    ) as stream:
        points = {
            int(row["point_id"]): (
                float(row["measured_ethylene_ppmv"]),
                float(row["predicted_mean_ethylene_ppmv"]),
            )
            for row in csv.DictReader(stream)
        }
    return record, points


def _metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    error = predicted - measured
    return {
        "measured_mean_ppmv": float(np.mean(measured)),
        "predicted_mean_ppmv": float(np.mean(predicted)),
        "bias_ppmv": float(np.mean(error)),
        "mae_ppmv": float(np.mean(np.abs(error))),
        "rmse_ppmv": float(np.sqrt(np.mean(error**2))),
        "predicted_p95_ppmv": float(np.percentile(predicted, 95)),
        "predicted_max_ppmv": float(np.max(predicted)),
    }


def compare(short_dir: Path, full_dir: Path, output_dir: Path) -> dict:
    short, short_points = _load(short_dir)
    full, full_points = _load(full_dir)
    short_config = short["configuration"]
    full_config = full["configuration"]
    if set(short_config) != set(full_config):
        raise ValueError("time-window pair has different configuration fields")
    changed = {
        key for key in short_config if short_config[key] != full_config[key]
    }
    if not changed <= _ALLOWED_CONFIGURATION_CHANGES or not {
        "scalar_duration_s", "scalar_mean_window_s",
    } <= changed:
        raise ValueError(
            "time-window pair may change only duration, averaging window, and "
            f"integrated released mass; found {sorted(changed)}"
        )
    if short["data"]["result_sha256"] != full["data"]["result_sha256"]:
        raise ValueError("AIJ observation data differ between the runs")
    common = sorted(short_points.keys() & full_points.keys())
    if not common:
        raise ValueError("short and full-window runs have no common points")
    measured = np.asarray([short_points[key][0] for key in common])
    full_measured = np.asarray([full_points[key][0] for key in common])
    if not np.array_equal(measured, full_measured):
        raise ValueError("measured concentration values differ on common points")
    short_pred = np.asarray([short_points[key][1] for key in common])
    full_pred = np.asarray([full_points[key][1] for key in common])
    delta = full_pred - short_pred
    output = {
        "schema": "degali.aij-case-h-time-window-sensitivity.v1",
        "status": "experiment_length_window_sensitivity_not_validation",
        "calibration_used": False,
        "allowed_configuration_changes": sorted(_ALLOWED_CONFIGURATION_CHANGES),
        "common_measurement_point_count": len(common),
        "source_sigma_m": short_config["source_sigma_m"],
        "short_case": {
            "duration_s": short_config["scalar_duration_s"],
            "mean_window_s": short_config["scalar_mean_window_s"],
            "integrated_released_species_mass_kg": short_config[
                "integrated_released_species_mass_kg"
            ],
            "scores": _metrics(measured, short_pred),
            "wall_clock_s": short["runtime"]["wall_clock_s"],
            "species_mass_residual": short["numerics"][
                "maximum_relative_species_mass_residual"
            ],
        },
        "full_case": {
            "duration_s": full_config["scalar_duration_s"],
            "mean_window_s": full_config["scalar_mean_window_s"],
            "integrated_released_species_mass_kg": full_config[
                "integrated_released_species_mass_kg"
            ],
            "scores": _metrics(measured, full_pred),
            "wall_clock_s": full["runtime"]["wall_clock_s"],
            "species_mass_residual": full["numerics"][
                "maximum_relative_species_mass_residual"
            ],
        },
        "paired_prediction_change": {
            "full_minus_short_mean_ppmv": float(np.mean(delta)),
            "full_minus_short_rmse_ppmv": float(np.sqrt(np.mean(delta**2))),
            "points_increased": int(np.count_nonzero(delta > 0.0)),
            "points_decreased": int(np.count_nonzero(delta < 0.0)),
        },
        "interpretation": (
            "This matched-input comparison isolates extending the modeled source "
            "and averaging window toward the 120 s experiment record. It does "
            "not isolate stationarity from transient accumulation, and remains "
            "a passive-ethylene component sensitivity rather than H2 validation."
        ),
        "promotion_allowed": False,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "aij-case-h-time-window-sensitivity.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    figure, axis = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    axis.scatter(short_pred, full_pred, c=measured, cmap="viridis", s=18, alpha=0.7)
    high = float(max(np.max(short_pred), np.max(full_pred)))
    axis.plot([0.0, high], [0.0, high], "k--", lw=1)
    axis.set(xlabel="12 s model mean, 8–12 s window (ppmv)",
             ylabel="120 s model mean, 0.001–120 s window (ppmv)")
    axis.set_title("AIJ Case H — matched inputs, experiment-length window\ncolor = measured 120 s mean; not validation")
    axis.grid(alpha=0.2)
    figure.savefig(output_dir / "aij-case-h-time-window-sensitivity.png", dpi=180)
    plt.close(figure)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--short", type=Path, default=DEFAULT_SHORT)
    parser.add_argument("--full", type=Path, default=DEFAULT_FULL)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.short, args.full, args.output_directory), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
