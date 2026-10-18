"""Compare no-SGS and SGS-mixed AIJ Case H tracer predictions on common IDs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = ROOT / "outputs/aij-case-h-resolved-tracer-matched-dx2p5cm-2026-10-12-v1"
DEFAULT_SGS = ROOT / "outputs/aij-case-h-tracer-sgs-sensitivity-dx2p5cm-2026-10-12-v1"


def _read(directory: Path) -> tuple[dict, dict[int, dict[str, float]]]:
    record = json.loads((directory / "aij-case-h-tracer-score.json").read_text(encoding="utf-8"))
    with (directory / "aij-case-h-tracer-point-scores.csv").open(
        "r", encoding="utf-8-sig", newline="",
    ) as stream:
        rows = {
            int(row["point_id"]): {
                key: float(row[key]) for key in (
                    "measured_ethylene_ppmv", "predicted_mean_ethylene_ppmv",
                )
            }
            for row in csv.DictReader(stream)
        }
    return record, rows


def _score(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    error = predicted - measured
    return {
        "bias_ppmv": float(np.mean(error)),
        "mae_ppmv": float(np.mean(np.abs(error))),
        "rmse_ppmv": float(np.sqrt(np.mean(error ** 2))),
        "measured_mean_ppmv": float(np.mean(measured)),
        "predicted_mean_ppmv": float(np.mean(predicted)),
        "pearson_r": float(np.corrcoef(measured, predicted)[0, 1]),
    }


def compare(baseline_directory: Path, sgs_directory: Path, output_directory: Path) -> dict:
    baseline_record, baseline_rows = _read(baseline_directory)
    sgs_record, sgs_rows = _read(sgs_directory)
    if baseline_record["configuration"]["scalar_mean_window_s"] != sgs_record["configuration"]["scalar_mean_window_s"]:
        raise ValueError("baseline and SGS runs must use the same averaging window")
    ids = sorted(baseline_rows.keys() & sgs_rows.keys())
    if not ids:
        raise ValueError("no common sensor observations")
    measured = np.asarray([baseline_rows[i]["measured_ethylene_ppmv"] for i in ids])
    if not np.allclose(measured, [sgs_rows[i]["measured_ethylene_ppmv"] for i in ids], rtol=0.0, atol=1.0e-12):
        raise ValueError("sensor values differ between runs")
    base = np.asarray([baseline_rows[i]["predicted_mean_ethylene_ppmv"] for i in ids])
    sgs = np.asarray([sgs_rows[i]["predicted_mean_ethylene_ppmv"] for i in ids])
    base_score, sgs_score = _score(measured, base), _score(measured, sgs)
    output = {
        "schema": "degali.aij-case-h-tracer-sgs-sensitivity.v1",
        "status": "adverse_component_score_with_limited_improvement_from_uncalibrated_sgs_sensitivity",
        "calibration_used": False,
        "comparison_window_s": baseline_record["configuration"]["scalar_mean_window_s"],
        "common_measurement_point_count": len(ids),
        "baseline_molecular_only": {
            "configuration": baseline_record["configuration"], "scores": base_score,
            "wall_clock_s": baseline_record["runtime"]["wall_clock_s"],
            "mass_residual": baseline_record["numerics"]["maximum_relative_species_mass_residual"],
        },
        "sgs_gradient_diffusion_sensitivity": {
            "configuration": sgs_record["configuration"], "scores": sgs_score,
            "wall_clock_s": sgs_record["runtime"]["wall_clock_s"],
            "mass_residual": sgs_record["numerics"]["maximum_relative_species_mass_residual"],
            "mean_applied_scalar_diffusivity_m2_s": sgs_record["numerics"]["applied_scalar_diffusivity_mean_m2_s"],
        },
        "rmse_reduction_fraction": float((base_score["rmse_ppmv"] - sgs_score["rmse_ppmv"]) / base_score["rmse_ppmv"]),
        "interpretation": (
            "SGS gradient diffusion lowers this same-window mean RMSE modestly, but both runs grossly overpredict; "
            "the SGS diffusivity is uncalibrated, source momentum/nozzle boundary is unresolved, and neither run validates H2 or facility performance. "
            "This is a closure sensitivity, not validation or promotion evidence."
        ),
        "promotion_allowed": False,
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    axis.scatter(base, sgs, c=measured, cmap="viridis", s=20, alpha=0.7)
    high = float(max(np.max(base), np.max(sgs)))
    axis.plot([0.0, high], [0.0, high], "k--", lw=1)
    axis.set(xlabel="Molecular-only prediction (ppmv)", ylabel="SGS-mixed prediction (ppmv)")
    axis.set_title(f"AIJ Case H — same 8–12 s window, n={len(ids)}\ncolor = measured ethylene; no calibration")
    axis.grid(alpha=0.2)
    figure.savefig(output_directory / "aij-case-h-tracer-sgs-sensitivity.png", dpi=180)
    plt.close(figure)
    (output_directory / "aij-case-h-tracer-sgs-sensitivity.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8",
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--sgs", type=Path, default=DEFAULT_SGS)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.baseline, args.sgs, args.output_directory), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
