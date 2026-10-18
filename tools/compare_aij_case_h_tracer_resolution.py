"""Compare the matched 8--12 s AIJ tracer prefixes on common measurement IDs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COARSE = ROOT / "outputs/aij-case-h-resolved-tracer-matched-dx2p5cm-2026-10-12-v1"
DEFAULT_FINE = ROOT / "outputs/aij-case-h-resolved-tracer-fine-2026-10-12-v1"


def _read_points(path: Path) -> dict[int, dict[str, float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = csv.DictReader(stream)
        return {
            int(row["point_id"]): {
                key: float(row[key])
                for key in (
                    "measured_ethylene_ppmv", "predicted_mean_ethylene_ppmv",
                    "x_m", "y_m", "z_m",
                )
            }
            for row in rows
        }


def _score(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    error = predicted - measured
    correlation = (
        float(np.corrcoef(measured, predicted)[0, 1])
        if np.std(measured) > 0.0 and np.std(predicted) > 0.0 else 0.0
    )
    return {
        "bias_ppmv": float(np.mean(error)),
        "mae_ppmv": float(np.mean(np.abs(error))),
        "rmse_ppmv": float(np.sqrt(np.mean(error ** 2))),
        "pearson_r": correlation,
        "mean_measured_ppmv": float(np.mean(measured)),
        "mean_predicted_ppmv": float(np.mean(predicted)),
    }


def compare(coarse_directory: Path, fine_directory: Path, output_directory: Path) -> dict[str, object]:
    coarse_record = json.loads((coarse_directory / "aij-case-h-tracer-score.json").read_text(encoding="utf-8"))
    fine_record = json.loads((fine_directory / "aij-case-h-tracer-score.json").read_text(encoding="utf-8"))
    if coarse_record["configuration"]["scalar_mean_window_s"] != fine_record["configuration"]["scalar_mean_window_s"]:
        raise ValueError("matched comparison requires identical scalar-mean windows")
    coarse = _read_points(coarse_directory / "aij-case-h-tracer-point-scores.csv")
    fine = _read_points(fine_directory / "aij-case-h-tracer-point-scores.csv")
    common = sorted(coarse.keys() & fine.keys())
    if not common:
        raise ValueError("no common point IDs between the two tracer score files")
    measured = np.asarray([coarse[i]["measured_ethylene_ppmv"] for i in common])
    if not np.allclose(
        measured,
        np.asarray([fine[i]["measured_ethylene_ppmv"] for i in common]),
        rtol=0.0, atol=1.0e-12,
    ):
        raise ValueError("measured concentration values differ for common AIJ IDs")
    coarse_prediction = np.asarray([coarse[i]["predicted_mean_ethylene_ppmv"] for i in common])
    fine_prediction = np.asarray([fine[i]["predicted_mean_ethylene_ppmv"] for i in common])
    output = {
        "schema": "degali.aij-case-h-tracer-resolution-sensitivity.v1",
        "status": "joint_grid_timestep_source_regularization_sensitivity_not_convergence",
        "calibration_used": False,
        "comparison_window_s": coarse_record["configuration"]["scalar_mean_window_s"],
        "common_measurement_point_count": len(common),
        "common_point_ids": common,
        "coarse": {
            "cell_size_m": coarse_record["configuration"]["cell_size_m"],
            "scalar_time_step_s": coarse_record["configuration"]["scalar_time_step_s"],
            "source_sigma_m": coarse_record["configuration"]["source_sigma_m"],
            "metrics_on_common_points": _score(measured, coarse_prediction),
            "flow_converged": coarse_record["numerics"]["flow_converged"],
            "flow_statistical_residual": coarse_record["numerics"]["flow_statistical_residual"],
            "relative_species_mass_residual": coarse_record["numerics"]["maximum_relative_species_mass_residual"],
        },
        "fine": {
            "cell_size_m": fine_record["configuration"]["cell_size_m"],
            "scalar_time_step_s": fine_record["configuration"]["scalar_time_step_s"],
            "source_sigma_m": fine_record["configuration"]["source_sigma_m"],
            "metrics_on_common_points": _score(measured, fine_prediction),
            "flow_converged": fine_record["numerics"]["flow_converged"],
            "flow_statistical_residual": fine_record["numerics"]["flow_statistical_residual"],
            "relative_species_mass_residual": fine_record["numerics"]["maximum_relative_species_mass_residual"],
        },
        "paired_prediction_change": {
            "fine_minus_coarse_mean_ppmv": float(np.mean(fine_prediction - coarse_prediction)),
            "fine_vs_coarse_rmse_ppmv": float(np.sqrt(np.mean((fine_prediction - coarse_prediction) ** 2))),
            "higher_at_common_points": int(np.count_nonzero(fine_prediction > coarse_prediction)),
            "lower_at_common_points": int(np.count_nonzero(fine_prediction < coarse_prediction)),
            "unchanged_at_common_points": int(np.count_nonzero(fine_prediction == coarse_prediction)),
        },
        "interpretation": (
            "Both grids materially overpredict measured mean ethylene. Refinement does not improve the error; "
            "because the scalar time step and Gaussian source width also change with grid size, this is not a pure spatial-convergence test. "
            "The neutral tracer component is withheld pending a physically represented nozzle source, source momentum, and turbulent scalar transport."
        ),
        "promotion_allowed": False,
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    axis.scatter(coarse_prediction, fine_prediction, c=measured, cmap="viridis", s=20, alpha=0.7)
    maximum = float(max(np.max(coarse_prediction), np.max(fine_prediction)))
    axis.plot([0.0, maximum], [0.0, maximum], "k--", lw=1)
    axis.set(xlabel="Coarse prediction (ppmv)", ylabel="Fine prediction (ppmv)")
    axis.set_title("AIJ tracer common points — same 8–12 s mean window\ncolor = measured ethylene (ppmv)")
    axis.grid(alpha=0.2)
    figure.savefig(output_directory / "aij-case-h-tracer-resolution-sensitivity.png", dpi=180)
    plt.close(figure)
    (output_directory / "aij-case-h-tracer-resolution-sensitivity.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", type=Path, default=DEFAULT_COARSE)
    parser.add_argument("--fine", type=Path, default=DEFAULT_FINE)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.coarse, args.fine, args.output_directory), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
