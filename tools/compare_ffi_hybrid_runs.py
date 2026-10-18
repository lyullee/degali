"""Compare fast, resolved and non-overlapping hybrid FFI O1 diagnostics.

All values remain conditional model diagnostics: the geometry and source
closures are not a machine-registered facility case.
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


def _delta(left: dict[str, object], right: dict[str, object], label: str) -> dict[str, object]:
    a = left.get("predicted_sensor_mean_vol_percent", {})
    b = right.get("predicted_sensor_mean_vol_percent", {})
    if not isinstance(a, dict) or not isinstance(b, dict):
        raise ValueError(f"{label} reports must contain predicted sensor means")
    labels = sorted(set(a) & set(b))
    if not labels:
        raise ValueError(f"{label} reports have no common sensor labels")
    values = np.asarray([float(b[key]) - float(a[key]) for key in labels], dtype=float)
    return {
        "common_sensor_count": len(labels),
        "mean_abs_right_minus_left_vol_percent": float(np.mean(np.abs(values))),
        "rmse_right_minus_left_vol_percent": float(np.sqrt(np.mean(values**2))),
        "bias_right_minus_left_vol_percent": float(np.mean(values)),
        "maximum_abs_right_minus_left_vol_percent": float(np.max(np.abs(values))),
    }


def _summary(report: dict[str, object]) -> dict[str, object]:
    return {
        "runtime_s": report.get("runtime_s"),
        "sensor_count": report.get("sensor_count"),
        "maximum_predicted_vol_percent": report.get("maximum_predicted_vol_percent"),
        "diagnostic": report.get("diagnostic_sensor_error_not_validation_score"),
        "mass_residual_kg": report.get(
            "maximum_partitioned_mass_balance_residual_kg",
            report.get("maximum_mass_residual_kg"),
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast", type=Path, default=ROOT / "outputs/ffi-o1-crossbeam-conditional-replay-2026-10-14.json")
    parser.add_argument("--resolved", type=Path, default=ROOT / "outputs/ffi-o1-resolved-buoyant-prefix-60s-dynamic-momz20-gate016-2026-10-14.json")
    parser.add_argument("--hybrid", type=Path, default=ROOT / "outputs/ffi-o1-hybrid-patch-60s-2026-10-14.json")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/ffi-o1-fast-resolved-hybrid-comparison-2026-10-14.json")
    args = parser.parse_args()
    fast, resolved, hybrid = (_read(path) for path in (args.fast, args.resolved, args.hybrid))
    report = {
        "schema": "degali.ffi-fast-resolved-hybrid-comparison.v1",
        "status": "conditional_model_comparison_not_validation_score",
        "reports": {"fast": str(args.fast), "resolved": str(args.resolved), "hybrid": str(args.hybrid)},
        "fast": _summary(fast),
        "resolved": _summary(resolved),
        "hybrid": _summary(hybrid),
        "common_sensor_mean_comparison": {
            "resolved_minus_fast": _delta(fast, resolved, "fast/resolved"),
            "hybrid_minus_fast": _delta(fast, hybrid, "fast/hybrid"),
            "hybrid_minus_resolved": _delta(resolved, hybrid, "resolved/hybrid"),
        },
        "warnings": [
            "all three paths use conditional or non-identical geometry/source closures",
            "the hybrid patch applies characteristic packet velocity; its geometry and carrier closure remain conditional",
            "no metric here is a facility validation score or accuracy ranking",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
