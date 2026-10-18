"""Compare two isothermal FFI hybrid prefixes for time-step sensitivity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"report must be an object: {path}")
    return value


def _field_path(report: dict[str, object], report_path: Path) -> Path:
    raw = report.get("field_artifact")
    if not isinstance(raw, str):
        raise ValueError("report is missing field_artifact")
    path = Path(raw)
    if path.is_absolute():
        return path
    return (ROOT / path) if path.parts and path.parts[0] == "outputs" else report_path.parent / path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    coarse = _read(args.coarse)
    fine = _read(args.fine)
    coarse_means = coarse.get("predicted_sensor_mean_vol_percent", {})
    fine_means = fine.get("predicted_sensor_mean_vol_percent", {})
    if not isinstance(coarse_means, dict) or not isinstance(fine_means, dict):
        raise ValueError("reports must contain predicted sensor means")
    labels = sorted(set(coarse_means) & set(fine_means))
    delta = np.asarray(
        [float(fine_means[label]) - float(coarse_means[label]) for label in labels],
        dtype=float,
    )
    coarse_field_path = _field_path(coarse, args.coarse)
    fine_field_path = _field_path(fine, args.fine)
    with np.load(coarse_field_path) as coarse_data, np.load(fine_field_path) as fine_data:
        coarse_field = np.asarray(coarse_data["maximum_concentration_vol_percent"], dtype=float)
        fine_field = np.asarray(fine_data["maximum_concentration_vol_percent"], dtype=float)
    if coarse_field.shape != fine_field.shape:
        raise ValueError(f"field shapes differ: {coarse_field.shape} vs {fine_field.shape}")
    field_delta = fine_field - coarse_field
    report = {
        "schema": "degali.ffi-hybrid-temporal-convergence.v1",
        "status": "isothermal_time_convergence_audit_not_facility_validation",
        "coarse_report": str(args.coarse),
        "fine_report": str(args.fine),
        "coarse_time_step_s": coarse.get("time_step_s"),
        "fine_time_step_s": fine.get("time_step_s"),
        "common_sensor_count": len(labels),
        "sensor_mean_delta_vol_percent": {
            "mean_abs": float(np.mean(np.abs(delta))),
            "rmse": float(np.sqrt(np.mean(delta**2))),
            "bias": float(np.mean(delta)),
            "maximum_abs": float(np.max(np.abs(delta))),
        },
        "field_delta_vol_percent": {
            "mean_abs": float(np.mean(np.abs(field_delta))),
            "rmse": float(np.sqrt(np.mean(field_delta**2))),
            "maximum_abs": float(np.max(np.abs(field_delta))),
        },
        "runtime_s": {
            "coarse": coarse.get("runtime_s"),
            "fine": fine.get("runtime_s"),
            "fine_over_coarse": float(float(fine["runtime_s"]) / float(coarse["runtime_s"])),
        },
        "mass_residual_kg": {
            "coarse": coarse.get("maximum_partitioned_mass_balance_residual_kg"),
            "fine": fine.get("maximum_partitioned_mass_balance_residual_kg"),
        },
        "warnings": [
            "both runs are intentionally isothermal because the active thermal solver requires the D3Q19 physical step",
            "this is a numerical time/interface audit, not an accuracy score or facility promotion",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
