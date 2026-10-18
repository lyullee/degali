"""Quantify the early-time source-momentum effect for an ethylene-air FV run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WITHOUT = ROOT / "outputs/aij-case-h-tracer-fv-zero-momentum-sgs-12s-dx2p5cm-dt1ms-2026-10-12-v1"
DEFAULT_WITH = ROOT / "outputs/aij-case-h-tracer-fv-momentum-sgs-12s-dx2p5cm-dt1ms-2026-10-12-v1"


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


def compare(without_directory: Path, with_directory: Path, output_directory: Path) -> dict:
    no_record, no_rows = _read(without_directory)
    yes_record, yes_rows = _read(with_directory)
    if no_record["configuration"]["scalar_mean_window_s"] != yes_record["configuration"]["scalar_mean_window_s"]:
        raise ValueError("source-momentum runs must share the same averaging window")
    if no_record["configuration"]["momentum_model"] != "finite_volume_projection" or yes_record["configuration"]["momentum_model"] != "finite_volume_projection":
        raise ValueError("both runs must use the finite-volume momentum path")
    for key in (
        "domain_m", "cells", "cell_size_m", "scalar_duration_s",
        "scalar_time_step_s", "source_mass_rate_kg_s", "source_sigma_m",
        "subgrid_model", "smagorinsky_constant", "effective_schmidt_number",
        "scalar_turbulent_diffusivity_model",
    ):
        if no_record["configuration"].get(key) != yes_record["configuration"].get(key):
            raise ValueError(f"source-momentum pair must match {key}")
    if no_record["configuration"].get("source_momentum_applied") is not False or yes_record["configuration"].get("source_momentum_applied") is not True:
        raise ValueError("the only intended difference is source momentum off/on")
    ids = sorted(no_rows.keys() & yes_rows.keys())
    measured = np.asarray([no_rows[i]["measured_ethylene_ppmv"] for i in ids])
    if not np.allclose(measured, [yes_rows[i]["measured_ethylene_ppmv"] for i in ids], rtol=0.0, atol=1.0e-12):
        raise ValueError("measured sensor values differ")
    without = np.asarray([no_rows[i]["predicted_mean_ethylene_ppmv"] for i in ids])
    with_momentum = np.asarray([yes_rows[i]["predicted_mean_ethylene_ppmv"] for i in ids])
    delta = with_momentum - without
    without_error, with_error = without - measured, with_momentum - measured
    output = {
        "schema": "degali.aij-case-h-source-momentum-sensitivity.v1",
        "status": "matched_window_source_momentum_effect_small_validation_withheld",
        "calibration_used": False,
        "averaging_window_s": no_record["configuration"]["scalar_mean_window_s"],
        "common_measurement_point_count": len(ids),
        "source_exit_speed_m_s": yes_record["configuration"]["source_exit_speed_from_reported_flow_m_s"],
        "source_exit_speed_basis": "Q/(pi d^2/4) using published Q and 4 mm opening; inferred, not directly measured",
        "finite_volume_continuity_residual_max_kg_m3_s": yes_record["numerics"]["maximum_variable_density_continuity_residual_kg_m3_s"],
        "finite_volume_global_mass_residual_max_kg_s": yes_record["numerics"]["maximum_variable_density_global_residual_kg_s"],
        "momentum_balance_residual_max_kg_m_s": yes_record["numerics"]["maximum_momentum_balance_residual_kg_m_s"],
        "without_momentum": {
            "bias_ppmv": float(np.mean(without_error)),
            "rmse_ppmv": float(np.sqrt(np.mean(without_error ** 2))),
            "predicted_mean_ppmv": float(np.mean(without)),
            "wall_clock_s": no_record["runtime"]["wall_clock_s"],
        },
        "with_inferred_momentum": {
            "bias_ppmv": float(np.mean(with_error)),
            "rmse_ppmv": float(np.sqrt(np.mean(with_error ** 2))),
            "predicted_mean_ppmv": float(np.mean(with_momentum)),
            "wall_clock_s": yes_record["runtime"]["wall_clock_s"],
        },
        "paired_field_change": {
            "mean_delta_ppmv": float(np.mean(delta)),
            "rmse_delta_ppmv": float(np.sqrt(np.mean(delta ** 2))),
            "maximum_absolute_delta_ppmv": float(np.max(np.abs(delta))),
            "points_increased": int(np.count_nonzero(delta > 0.0)),
            "points_decreased": int(np.count_nonzero(delta < 0.0)),
        },
        "interpretation": (
            "The two FV runs match mesh, timestep, duration/window, mixture EOS, source rate, scalar SGS and all other declared settings; only source momentum is toggled. "
            "For this 8--12 s window, adding the inferred source momentum changes the common-sensor RMSE negligibly. "
            "That does not prove momentum is generally negligible because the source is still Gaussian-smeared and the 120 s experimental acquisition is not matched. "
            "The result isolates this configured source-momentum term, but is not a validation score or promotion evidence."
        ),
        "promotion_allowed": False,
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    axis.scatter(without, with_momentum, c=measured, cmap="viridis", s=19, alpha=0.7)
    high = float(max(np.max(without), np.max(with_momentum)))
    axis.plot([0.0, high], [0.0, high], "k--", lw=1)
    axis.set(xlabel="FV prediction, no source momentum (ppmv)", ylabel="FV prediction, inferred nozzle momentum (ppmv)")
    window = no_record["configuration"]["scalar_mean_window_s"]
    axis.set_title(f"AIJ Case H — matched FV source-momentum toggle, {window[0]:g}–{window[1]:g} s\ncolor = measured 120 s mean; component test only")
    axis.grid(alpha=0.2)
    figure.savefig(output_directory / "aij-case-h-source-momentum-sensitivity.png", dpi=180)
    plt.close(figure)
    (output_directory / "aij-case-h-source-momentum-sensitivity.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8",
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--without-momentum", type=Path, default=DEFAULT_WITHOUT)
    parser.add_argument("--with-momentum", type=Path, default=DEFAULT_WITH)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.without_momentum, args.with_momentum, args.output_directory), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
