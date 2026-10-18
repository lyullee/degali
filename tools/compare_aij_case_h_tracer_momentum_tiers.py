"""No-fit matched-window comparison of D3Q19 and binary-FV AIJ tracer runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_D3Q19 = ROOT / "outputs/aij-case-h-tracer-d3q19-sgs-matcheddt-12s-dx2p5cm-dt1ms-2026-10-12-v1"
DEFAULT_FV = ROOT / "outputs/aij-case-h-tracer-fv-momentum-sgs-12s-dx2p5cm-dt1ms-2026-10-12-v1"


def _read(directory: Path) -> tuple[dict, dict[int, dict[str, float]]]:
    record = json.loads((directory / "aij-case-h-tracer-score.json").read_text(encoding="utf-8"))
    with (directory / "aij-case-h-tracer-point-scores.csv").open(
        "r", encoding="utf-8-sig", newline="",
    ) as stream:
        rows = {
            int(row["point_id"]): {
                "measured": float(row["measured_ethylene_ppmv"]),
                "predicted": float(row["predicted_mean_ethylene_ppmv"]),
                "z_m": float(row["z_m"]),
            }
            for row in csv.DictReader(stream)
        }
    return record, rows


def _scores(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    error = predicted - measured
    return {
        "bias_ppmv": float(np.mean(error)),
        "mae_ppmv": float(np.mean(np.abs(error))),
        "rmse_ppmv": float(np.sqrt(np.mean(error**2))),
        "mean_measured_ppmv": float(np.mean(measured)),
        "mean_predicted_ppmv": float(np.mean(predicted)),
        "pearson_r": float(np.corrcoef(measured, predicted)[0, 1]),
    }


def compare(d3q_directory: Path, fv_directory: Path, output_directory: Path) -> dict:
    d3q, d3q_rows = _read(d3q_directory)
    fv, fv_rows = _read(fv_directory)
    dc, fc = d3q["configuration"], fv["configuration"]
    for key in (
        "domain_m", "cells", "cell_size_m", "scalar_duration_s",
        "scalar_time_step_s", "scalar_mean_window_s", "source_mass_rate_kg_s",
        "source_sigma_m", "smagorinsky_constant", "effective_schmidt_number",
    ):
        if dc.get(key) != fc.get(key):
            raise ValueError(f"matched comparison requires identical {key}: {dc.get(key)} != {fc.get(key)}")
    if dc["momentum_model"] != "d3q19_bgk" or fc["momentum_model"] != "finite_volume_projection":
        raise ValueError("expected D3Q19 baseline and finite-volume candidate")
    ids = sorted(d3q_rows.keys() & fv_rows.keys())
    if not ids:
        raise ValueError("the two runs share no measurement IDs")
    measured = np.asarray([d3q_rows[i]["measured"] for i in ids])
    if not np.allclose(measured, [fv_rows[i]["measured"] for i in ids], rtol=0.0, atol=1.0e-12):
        raise ValueError("sensor measurement values do not match")
    d3q_pred = np.asarray([d3q_rows[i]["predicted"] for i in ids])
    fv_pred = np.asarray([fv_rows[i]["predicted"] for i in ids])
    reduction = (
        _scores(measured, d3q_pred)["rmse_ppmv"]
        - _scores(measured, fv_pred)["rmse_ppmv"]
    ) / _scores(measured, d3q_pred)["rmse_ppmv"]
    output = {
        "schema": "degali.aij-case-h-tracer-momentum-tier-comparison.v1",
        "status": "matched_8to12s_component_candidate_improves_score_but_validation_withheld",
        "calibration_used": False,
        "matched_conditions": {
            "domain_m": dc["domain_m"], "cells": dc["cells"],
            "cell_size_m": dc["cell_size_m"],
            "duration_s": dc["scalar_duration_s"],
            "time_step_s": dc["scalar_time_step_s"],
            "mean_window_s": dc["scalar_mean_window_s"],
            "source_mass_rate_kg_s": dc["source_mass_rate_kg_s"],
            "source_sigma_m": dc["source_sigma_m"],
            "subgrid_model": dc["subgrid_model"],
            "smagorinsky_constant": dc["smagorinsky_constant"],
            "effective_schmidt_number": dc["effective_schmidt_number"],
            "comparison_ids": len(ids),
        },
        "D3Q19_molecular_air_carrier_with_SGS_scalar": {
            "scores": _scores(measured, d3q_pred),
            "runtime_wall_s": d3q["runtime"]["wall_clock_s"],
            "species_mass_residual": d3q["numerics"]["maximum_relative_species_mass_residual"],
        },
        "finite_volume_binary_ethylene_air_with_inferred_nozzle_momentum": {
            "scores": _scores(measured, fv_pred),
            "runtime_wall_s": fv["runtime"]["wall_clock_s"],
            "species_mass_residual": fv["numerics"]["maximum_relative_species_mass_residual"],
            "continuity_residual_kg_m3_s": fv["numerics"]["maximum_variable_density_continuity_residual_kg_m3_s"],
            "global_mass_residual_kg_s": fv["numerics"]["maximum_variable_density_global_residual_kg_s"],
            "momentum_balance_residual_kg_m_s": fv["numerics"]["maximum_momentum_balance_residual_kg_m_s"],
            "source_exit_velocity_m_s_inferred": fc["source_exit_speed_from_reported_flow_m_s"],
        },
        "rmse_reduction_fraction": float(reduction),
        "paired_field_change": {
            "mean_fv_minus_d3q19_ppmv": float(np.mean(fv_pred - d3q_pred)),
            "rmse_between_predictions_ppmv": float(np.sqrt(np.mean((fv_pred - d3q_pred)**2))),
            "maximum_absolute_change_ppmv": float(np.max(np.abs(fv_pred - d3q_pred))),
        },
        "interpretation": (
            "Same mesh, source mass, time step/window and SGS scalar parameters give a direct model-tier comparison. "
            "FV density/pressure coupling and inferred source momentum are both changed together, so the improvement cannot be attributed to nozzle momentum alone. "
            "The 12 s model mean is still shorter than the experiment's 120 s acquisition; this component result is not H2/facility validation."
        ),
        "promotion_allowed": False,
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.5), constrained_layout=True)
    height = np.asarray([d3q_rows[i]["z_m"] for i in ids])
    measured_max = float(np.max(measured))
    predicted_max = float(max(np.max(d3q_pred), np.max(fv_pred)))
    comparison_max = min(measured_max * 1.08, predicted_max)
    for axis, prediction, label, score in (
        (axes[0], d3q_pred, "D3Q19 carrier + SGS scalar", output["D3Q19_molecular_air_carrier_with_SGS_scalar"]["scores"]),
        (axes[1], fv_pred, "FV ethylene-air + nozzle momentum", output["finite_volume_binary_ethylene_air_with_inferred_nozzle_momentum"]["scores"]),
    ):
        points = axis.scatter(measured, prediction, c=height, cmap="viridis", s=17, alpha=0.65)
        axis.plot([0.0, comparison_max], [0.0, comparison_max], "k--", lw=1)
        axis.set_xlim(0.0, measured_max * 1.08)
        axis.set_ylim(0.0, predicted_max * 1.04)
        axis.set_title(f"{label}\nRMSE {score['rmse_ppmv']:.0f} ppmv, r={score['pearson_r']:.3f}")
        axis.set_xlabel("Measured mean ethylene (ppmv)")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Predicted mean ethylene (ppmv)")
    figure.colorbar(points, ax=axes, label="Sensor height (m)")
    figure.suptitle("AIJ Case H — matched 8–12 s window; no calibration")
    figure.savefig(output_directory / "aij-case-h-momentum-tier-comparison.png", dpi=180)
    plt.close(figure)
    (output_directory / "aij-case-h-momentum-tier-comparison.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8",
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--d3q19", type=Path, default=DEFAULT_D3Q19)
    parser.add_argument("--finite-volume", type=Path, default=DEFAULT_FV)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.d3q19, args.finite_volume, args.output_directory), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
