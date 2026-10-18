"""Compare AIJ FV tracer predictions on identical point IDs for lateral-domain sensitivity."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _read(path: Path) -> dict[int, dict[str, float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return {
            int(row["point_id"]): {
                "measured": float(row["measured_ethylene_ppmv"]),
                "predicted": float(row["predicted_mean_ethylene_ppmv"]),
            }
            for row in csv.DictReader(stream)
        }


def _metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    error = predicted - measured
    return {
        "measured_mean_ppmv": float(np.mean(measured)),
        "predicted_mean_ppmv": float(np.mean(predicted)),
        "bias_ppmv": float(np.mean(error)),
        "mae_ppmv": float(np.mean(np.abs(error))),
        "rmse_ppmv": float(np.sqrt(np.mean(error**2))),
        "pearson_r": float(np.corrcoef(measured, predicted)[0, 1]),
    }


def compare(base_directory: Path, wide_directory: Path, output_directory: Path) -> dict[str, object]:
    base = _read(base_directory / "aij-case-h-tracer-point-scores.csv")
    wide = _read(wide_directory / "aij-case-h-tracer-point-scores.csv")
    ids = sorted(set(base) & set(wide))
    if not ids:
        raise ValueError("no common AIJ point IDs")
    measured_base = np.asarray([base[i]["measured"] for i in ids])
    measured_wide = np.asarray([wide[i]["measured"] for i in ids])
    if not np.array_equal(measured_base, measured_wide):
        raise ValueError("measured values differ for common point IDs")
    base_pred = np.asarray([base[i]["predicted"] for i in ids])
    wide_pred = np.asarray([wide[i]["predicted"] for i in ids])
    base_score = _metrics(measured_base, base_pred)
    wide_score = _metrics(measured_base, wide_pred)
    delta = wide_pred - base_pred
    output = {
        "schema": "degali.aij-case-h-lateral-domain-tracer-comparison.v1",
        "status": "matched_point_lateral_boundary_sensitivity_not_validation",
        "calibration_used": False,
        "common_point_count": len(ids),
        "common_point_ids": ids,
        "base_directory": str(base_directory),
        "wide_directory": str(wide_directory),
        "base_common_point_metrics": base_score,
        "wide_common_point_metrics": wide_score,
        "wide_minus_base_prediction_change": {
            "mean_ppmv": float(np.mean(delta)),
            "mae_ppmv": float(np.mean(np.abs(delta))),
            "rmse_ppmv": float(np.sqrt(np.mean(delta**2))),
            "maximum_abs_ppmv": float(np.max(np.abs(delta))),
        },
        "rmse_reduction_fraction": float((base_score["rmse_ppmv"] - wide_score["rmse_ppmv"]) / base_score["rmse_ppmv"]),
        "interpretation": (
            "The paired comparison changes only lateral domain and coordinate origin while keeping source, mesh spacing, FV/SGS and averaging window fixed. "
            "It diagnoses side-boundary sensitivity; it is not H2 or facility validation."
        ),
        "promotion_allowed": False,
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "aij-case-h-lateral-domain-tracer-comparison.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    figure, axis = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    axis.scatter(base_pred, wide_pred, c=measured_base, cmap="viridis", s=18, alpha=0.7)
    high = float(max(np.max(base_pred), np.max(wide_pred)))
    axis.plot([0.0, high], [0.0, high], "k--", lw=1)
    axis.set(xlabel="Base 0.30 m lateral prediction (ppmv)", ylabel="Wide 0.50 m lateral prediction (ppmv)")
    axis.set_title(f"AIJ Case H — paired lateral-domain sensitivity, n={len(ids)}")
    axis.grid(alpha=0.2)
    figure.savefig(output_directory / "aij-case-h-lateral-domain-tracer-comparison.png", dpi=180)
    plt.close(figure)
    print(json.dumps(output, indent=2))
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--wide", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    compare(args.base, args.wide, args.output_directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
