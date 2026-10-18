"""Compare AIJ Case H scores while changing only Gaussian source sigma."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = ROOT / (
    "outputs/aij-case-h-tracer-fv-momentum-sgs-12s-dx2p5cm-dt1ms-"
    "2026-10-12-v1"
)
DEFAULT_NARROW = ROOT / (
    "outputs/aij-case-h-tracer-fv-source-sigma-1mm-matched-12s-"
    "2026-10-13-v1"
)


def _load(directory: Path) -> tuple[dict, dict[int, tuple[float, float]]]:
    record = json.loads(
        (directory / "aij-case-h-tracer-score.json").read_text(encoding="utf-8")
    )
    with (directory / "aij-case-h-tracer-point-scores.csv").open(
        "r", encoding="utf-8-sig", newline="",
    ) as stream:
        rows = {
            int(row["point_id"]): (
                float(row["measured_ethylene_ppmv"]),
                float(row["predicted_mean_ethylene_ppmv"]),
            )
            for row in csv.DictReader(stream)
        }
    return record, rows


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


def compare(baseline_dir: Path, narrow_dir: Path, output_dir: Path) -> dict:
    base, base_rows = _load(baseline_dir)
    narrow, narrow_rows = _load(narrow_dir)
    base_config = base["configuration"]
    narrow_config = narrow["configuration"]
    if set(base_config) != set(narrow_config):
        raise ValueError("source-width pair has different configuration fields")
    changed = [
        key for key in base_config
        if base_config[key] != narrow_config[key]
    ]
    if changed != ["source_sigma_m"]:
        raise ValueError(
            "matched source-regularization comparison requires source_sigma_m "
            f"to be the only configuration change; found {changed}"
        )
    if base["data"]["result_sha256"] != narrow["data"]["result_sha256"]:
        raise ValueError("AIJ experimental source files differ between runs")
    common = sorted(base_rows.keys() & narrow_rows.keys())
    if not common:
        raise ValueError("source-width runs have no common sensor points")
    measured = np.asarray([base_rows[key][0] for key in common])
    measured_narrow = np.asarray([narrow_rows[key][0] for key in common])
    if not np.array_equal(measured, measured_narrow):
        raise ValueError("measured values differ on common sensor IDs")
    baseline = np.asarray([base_rows[key][1] for key in common])
    narrow_prediction = np.asarray([narrow_rows[key][1] for key in common])
    delta = narrow_prediction - baseline
    base_metrics = _metrics(measured, baseline)
    narrow_metrics = _metrics(measured, narrow_prediction)
    fwhm_per_sigma = float(2.0 * np.sqrt(2.0 * np.log(2.0)))
    cell_size = float(base_config["cell_size_m"])
    source_resolution = {
        "criterion": "Gaussian FWHM at least 2 grid cells in every axis; diagnostic only",
        "baseline_fwhm_cell_count_xyz": [
            fwhm_per_sigma * base_config["source_sigma_m"] / cell_size
        ] * 3,
        "baseline_underresolved": bool(
            fwhm_per_sigma * base_config["source_sigma_m"] / cell_size < 2.0
        ),
        "narrow_fwhm_cell_count_xyz": [
            fwhm_per_sigma * narrow_config["source_sigma_m"] / cell_size
        ] * 3,
        "narrow_underresolved": bool(
            fwhm_per_sigma * narrow_config["source_sigma_m"] / cell_size < 2.0
        ),
    }
    output = {
        "schema": "degali.aij-case-h-source-regularization-sensitivity.v1",
        "status": "matched_source_width_sensitivity_not_validation_or_convergence",
        "calibration_used": False,
        "only_changed_configuration_field": "source_sigma_m",
        "comparison_window_s": base_config["scalar_mean_window_s"],
        "common_measurement_point_count": len(common),
        "source_sigma_m": {
            "baseline": base_config["source_sigma_m"],
            "narrow": narrow_config["source_sigma_m"],
            "ratio_narrow_to_baseline": (
                narrow_config["source_sigma_m"] / base_config["source_sigma_m"]
            ),
            "physical_scale_note": (
                "1 mm equals D/4 for a 4 mm uniform circular outlet's radial "
                "second moment, but remains sub-grid (25 mm cells) and is not "
                "a resolved nozzle boundary"
            ),
        },
        "source_shape_resolution": source_resolution,
        "baseline": {
            "metrics": base_metrics,
            "wall_clock_s": base["runtime"]["wall_clock_s"],
            "species_mass_residual": base["numerics"][
                "maximum_relative_species_mass_residual"
            ],
        },
        "narrow_source": {
            "metrics": narrow_metrics,
            "wall_clock_s": narrow["runtime"]["wall_clock_s"],
            "species_mass_residual": narrow["numerics"][
                "maximum_relative_species_mass_residual"
            ],
        },
        "paired_prediction_change": {
            "narrow_minus_baseline_mean_ppmv": float(np.mean(delta)),
            "narrow_minus_baseline_rmse_ppmv": float(np.sqrt(np.mean(delta**2))),
            "points_increased": int(np.count_nonzero(delta > 0.0)),
            "points_decreased": int(np.count_nonzero(delta < 0.0)),
        },
        "interpretation": (
            "Changing the sub-grid Gaussian width alone changes the concentration "
            "field. The narrower source is an under-resolved source-model "
            "sensitivity, not an observation-fitted result, resolved-orifice "
            "calculation, or physical validation."
        ),
        "promotion_allowed": False,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "aij-case-h-source-regularization-sensitivity.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    figure, axis = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    axis.scatter(baseline, narrow_prediction, c=measured, cmap="viridis", s=18, alpha=0.7)
    high = float(max(np.max(baseline), np.max(narrow_prediction)))
    axis.plot([0.0, high], [0.0, high], "k--", lw=1)
    axis.set(
        xlabel=f"Baseline prediction, σ={base_config['source_sigma_m'] * 1000:g} mm (ppmv)",
        ylabel=f"Narrow-source prediction, σ={narrow_config['source_sigma_m'] * 1000:g} mm (ppmv)",
    )
    axis.set_title("AIJ Case H — matched 8–12 s source-width sensitivity\ncolor = measured 120 s mean; not validation")
    axis.grid(alpha=0.2)
    figure.savefig(output_dir / "aij-case-h-source-regularization-sensitivity.png", dpi=180)
    plt.close(figure)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--narrow", type=Path, default=DEFAULT_NARROW)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.baseline, args.narrow, args.output_directory), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
