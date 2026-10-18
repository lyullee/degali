"""Compare conditional fast and resolved FFI O1 replay reports.

This is a diagnostic comparison only.  The fast report uses a conditional
beam proxy and the resolved report uses a separately declared lower-blockage
sensitivity geometry, so the result must not be presented as a facility score.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"report must be an object: {path}")
    return payload


def _common_mean_error(
    fast: dict[str, object], resolved: dict[str, object],
) -> dict[str, float | int]:
    fast_mean = fast.get("predicted_sensor_mean_vol_percent", {})
    resolved_mean = resolved.get("predicted_sensor_mean_vol_percent", {})
    if not isinstance(fast_mean, dict) or not isinstance(resolved_mean, dict):
        raise ValueError("reports must contain predicted sensor means")
    labels = sorted(set(fast_mean) & set(resolved_mean))
    if not labels:
        raise ValueError("reports have no common sensor labels")
    delta = np.asarray(
        [float(resolved_mean[label]) - float(fast_mean[label]) for label in labels],
        dtype=float,
    )
    return {
        "common_sensor_count": len(labels),
        "mean_abs_resolved_minus_fast_vol_percent": float(np.mean(np.abs(delta))),
        "rmse_resolved_minus_fast_vol_percent": float(np.sqrt(np.mean(delta**2))),
        "bias_resolved_minus_fast_vol_percent": float(np.mean(delta)),
        "maximum_abs_resolved_minus_fast_vol_percent": float(np.max(np.abs(delta))),
    }


def build(fast_path: Path, resolved_path: Path) -> dict[str, object]:
    fast = _read(fast_path)
    resolved = _read(resolved_path)
    return {
        "schema": "degali.ffi-fast-resolved-comparison.v1",
        "status": "conditional_model_comparison_not_validation_score",
        "fast_report": str(fast_path),
        "resolved_report": str(resolved_path),
        "fast": {
            "geometry_mode": fast.get("geometry_mode"),
            "sensor_count": fast.get("sensor_count"),
            "runtime_s": fast.get("runtime_s"),
            "maximum_predicted_vol_percent": fast.get("maximum_predicted_vol_percent"),
            "diagnostic": fast.get("diagnostic_sensor_error_not_validation_score"),
            "mass_residual_kg": fast.get("maximum_mass_residual_kg"),
        },
        "resolved": {
            "geometry_mode": resolved.get("geometry_mode"),
            "sensor_count": resolved.get("sensor_count"),
            "runtime_s": resolved.get("runtime_s"),
            "maximum_predicted_vol_percent": resolved.get("maximum_predicted_vol_percent"),
            "diagnostic": resolved.get("diagnostic_sensor_error_not_validation_score"),
            "mass_residual_kg": resolved.get("maximum_mass_residual_kg"),
        },
        "common_sensor_mean_comparison": _common_mean_error(fast, resolved),
        "warnings": [
            "fast and resolved geometries are not machine-registered identical facility CAD",
            "resolved source momentum/Mach closure is a conditional sensitivity",
            "comparison does not authorize promotion or accuracy ranking",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fast", type=Path,
        default=ROOT / "outputs/ffi-o1-crossbeam-conditional-replay-2026-10-14.json",
    )
    parser.add_argument(
        "--resolved", type=Path,
        default=ROOT / "outputs/ffi-o1-resolved-buoyant-prefix-60s-dynamic-momz20-gate016-2026-10-14.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "outputs/ffi-o1-fast-resolved-comparison-2026-10-14.json",
    )
    args = parser.parse_args()
    report = build(args.fast, args.resolved)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
