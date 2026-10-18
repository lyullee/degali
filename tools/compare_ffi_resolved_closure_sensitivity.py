"""Compare conditional FFI resolved-physics closure witnesses.

This is deliberately a diagnostic comparison, not a validation scorer: the
published O1 beam CAD, source/sensor clock operator, and event calibration
certificate remain unresolved.  The tool makes it difficult to mistake a
small closure perturbation for a solution to those upstream gates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") != "conditional_resolved_physics_witness_not_validation_score":
        raise ValueError(f"unexpected report status in {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--variant", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    baseline = load(args.baseline)
    base_diag = baseline["diagnostic_sensor_error_not_validation_score"]
    base_rmse = float(base_diag["rmse_vol_percent"])
    base_mae = float(base_diag["mae_vol_percent"])
    base_mape = float(base_diag["mape_percent_floor_0p1_vol_percent"])
    rows: list[dict] = []
    for path in [args.baseline, *args.variant]:
        report = baseline if path == args.baseline else load(path)
        diag = report["diagnostic_sensor_error_not_validation_score"]
        physics = report["physics"]
        rows.append(
            {
                "artifact": str(path),
                "runtime_s": float(report["runtime_s"]),
                "mae_vol_percent": float(diag["mae_vol_percent"]),
                "rmse_vol_percent": float(diag["rmse_vol_percent"]),
                "mape_percent_floor_0p1_vol_percent": float(
                    diag["mape_percent_floor_0p1_vol_percent"]
                ),
                "bias_vol_percent": float(diag["bias_vol_percent"]),
                "mass_residual_kg": float(report["maximum_mass_residual_kg"]),
                "energy_residual_j": float(report["maximum_energy_residual_j"]),
                "closure": {
                    "subgrid_model": physics.get("subgrid_model"),
                    "smagorinsky_constant": physics.get("smagorinsky_constant"),
                    "wale_constant": physics.get("wale_constant"),
                    "effective_schmidt_number": physics.get("effective_schmidt_number"),
                    "effective_prandtl_number": physics.get("effective_prandtl_number"),
                    "scalar_turbulent_diffusivity_model": physics.get(
                        "scalar_turbulent_diffusivity_model"
                    ),
                    "thermal_turbulent_diffusivity_model": physics.get(
                        "thermal_turbulent_diffusivity_model"
                    ),
                    "source_sigma_m": physics.get("source_sigma_m"),
                    "source_momentum_velocity_m_s": physics.get(
                        "source_momentum_velocity_m_s"
                    ),
                },
            }
        )

    for row in rows:
        row["delta_vs_baseline"] = {
            "mae_vol_percent": row["mae_vol_percent"] - base_mae,
            "rmse_vol_percent": row["rmse_vol_percent"] - base_rmse,
            "mape_percentage_points": (
                row["mape_percent_floor_0p1_vol_percent"] - base_mape
            ),
            "runtime_s": row["runtime_s"] - float(rows[0]["runtime_s"]),
        }

    best = min(rows, key=lambda row: row["rmse_vol_percent"])
    report = {
        "schema": "degali.ffi-o1-resolved-closure-sensitivity.v1",
        "status": "conditional_physics_sensitivity_not_validation",
        "baseline": str(args.baseline),
        "rows": rows,
        "best_rmse_artifact": best["artifact"],
        "interpretation": {
            "closure_effect": (
                "The tested SGS/source-shape perturbations are compared on the "
                "same conditional surrogate geometry and release-relative prefix."
            ),
            "promotion": (
                "No closure variant is promoted: the machine-registered beam CAD, "
                "source/sensor clock operator, and event calibration remain missing."
            ),
            "decision_rule": (
                "A lower diagnostic RMSE is a physics-path signal only; it is not "
                "evidence that the facility prediction is validated."
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
