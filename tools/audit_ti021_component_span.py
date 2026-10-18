"""Audit whether existing TI021 observation branches span the measured dose.

This is a diagnostic, not a calibration path.  It integrates the additive
component traces from a completed 600 s run and solves a global, bounded
linear-programming lower bound on dose MAPE.  Every component receives one
event-wide multiplier; no receptor-specific factor is permitted.  The result
therefore answers whether observation-amplitude tuning alone could reach the
target before another expensive physical sweep is attempted.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def minimum_global_component_mape(
    component_dose: np.ndarray,
    observed_dose: np.ndarray,
    *,
    fixed_component_index: int = 0,
    fixed_component_multiplier: float = 1.0,
    upper_multiplier: float = 1.0,
) -> tuple[float, np.ndarray]:
    """Return the exact bounded global-component MAPE lower bound.

    The fixed component normally represents the pre-handoff source.  All
    other components share one multiplier across every receptor.  Absolute
    percentage error is linearized with one non-negative auxiliary variable
    per positive-dose receptor.
    """

    dose = np.asarray(component_dose, dtype=float)
    observed = np.asarray(observed_dose, dtype=float)
    if dose.ndim != 2 or observed.shape != (dose.shape[1],):
        raise ValueError("component_dose must be (components, sensors)")
    if not 0 <= fixed_component_index < dose.shape[0]:
        raise ValueError("fixed_component_index is outside component_dose")
    if not np.isfinite(upper_multiplier) or upper_multiplier <= 0.0:
        raise ValueError("upper_multiplier must be positive and finite")
    positive = observed > 0.0
    if not np.any(positive):
        raise ValueError("at least one positive observed dose is required")
    if np.any(~np.isfinite(dose)) or np.any(dose < -1.0e-12):
        raise ValueError("component doses must be finite and non-negative")

    free = np.arange(dose.shape[0]) != fixed_component_index
    normalized = dose[free][:, positive].T / observed[positive, None]
    fixed = (
        float(fixed_component_multiplier)
        * dose[fixed_component_index, positive]
        / observed[positive]
    )
    target = 1.0 - fixed
    sensor_count = int(np.count_nonzero(positive))
    component_count = int(np.count_nonzero(free))

    # |A w - target| <= t
    objective = np.r_[np.zeros(component_count), np.full(sensor_count, 1.0 / sensor_count)]
    upper = np.hstack((normalized, -np.eye(sensor_count)))
    lower = np.hstack((-normalized, -np.eye(sensor_count)))
    result = linprog(
        objective,
        A_ub=np.vstack((upper, lower)),
        b_ub=np.r_[target, -target],
        bounds=[(0.0, float(upper_multiplier))] * component_count
        + [(0.0, None)] * sensor_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"component-span linear program failed: {result.message}")
    multipliers = np.empty(dose.shape[0], dtype=float)
    multipliers[fixed_component_index] = float(fixed_component_multiplier)
    multipliers[free] = result.x[:component_count]
    prediction = multipliers @ dose[:, positive]
    mape = float(np.mean(np.abs(prediction / observed[positive] - 1.0)) * 100.0)
    return mape, multipliers


def audit(input_path: Path, output_path: Path) -> dict:
    with np.load(input_path, allow_pickle=False) as archive:
        required = {
            "time_s",
            "observation_component_labels",
            "observation_component_concentration_vol_percent",
            "observed_dose_vol_percent_s",
            "predicted_dose_vol_percent_s",
        }
        missing = sorted(required.difference(archive.files))
        if missing:
            raise ValueError(f"component archive lacks {missing}")
        time_s = np.asarray(archive["time_s"], dtype=float)
        labels = np.asarray(archive["observation_component_labels"]).astype(str)
        traces = np.asarray(
            archive["observation_component_concentration_vol_percent"],
            dtype=float,
        )
        observed = np.asarray(archive["observed_dose_vol_percent_s"], dtype=float)
        reported = np.asarray(archive["predicted_dose_vol_percent_s"], dtype=float)
    if traces.shape != (labels.size, observed.size, time_s.size):
        raise ValueError("component trace dimensions do not match labels/sensors/time")
    component_dose = np.trapezoid(traces, time_s, axis=2)
    reconstructed = component_dose.sum(axis=0)
    closure = float(np.max(np.abs(reconstructed - reported)))
    if not np.allclose(reconstructed, reported, rtol=1.0e-12, atol=1.0e-10):
        raise RuntimeError(f"component dose attribution does not close: {closure:.17g}")
    positive = observed > 0.0
    component_total = component_dose[:, positive].sum(axis=1)
    total = float(component_total.sum())
    bounds = {}
    for upper_multiplier in (1.0, 2.0, 5.0):
        mape, multipliers = minimum_global_component_mape(
            component_dose,
            observed,
            fixed_component_index=0,
            upper_multiplier=upper_multiplier,
        )
        bounds[f"0_to_{upper_multiplier:g}"] = {
            "minimum_mape_percent": mape,
            "multipliers": {
                str(label): float(value)
                for label, value in zip(labels, multipliers, strict=True)
            },
        }
    baseline_mape = float(
        np.mean(np.abs(reported[positive] / observed[positive] - 1.0)) * 100.0
    )
    report = {
        "schema": "degali.ti021-component-span-audit.v1",
        "status": "diagnostic_only_not_a_model_calibration",
        "input": str(input_path.resolve()),
        "positive_observed_sensor_count": int(np.count_nonzero(positive)),
        "component_labels": labels.tolist(),
        "reported_mape_percent": baseline_mape,
        "maximum_component_dose_closure_vol_percent_s": closure,
        "component_positive_sensor_dose": {
            str(label): {
                "sum_vol_percent_s": float(value),
                "share": float(value / total) if total > 0.0 else 0.0,
            }
            for label, value in zip(labels, component_total, strict=True)
        },
        "global_component_multiplier_lower_bounds": bounds,
        "target_mape_percent": 10.0,
        "target_reachable_by_global_component_amplitudes": bool(
            min(item["minimum_mape_percent"] for item in bounds.values()) <= 10.0
        ),
        "interpretation": (
            "A false target flag means the existing component footprints cannot "
            "reach 10% MAPE by event-wide branch amplitudes alone; new transport "
            "structure is required. Multipliers are diagnostic and are not fed "
            "back to the solver."
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.input, args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
