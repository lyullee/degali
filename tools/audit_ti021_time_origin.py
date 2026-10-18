"""Audit TI021 event-clock and transport-delay mismatch.

The TI021 workbook contains arrival summaries, not measured concentration
traces over the first ten seconds.  This utility therefore reports two
separate quantities:

* a workbook-only advection sanity check using the release-aligned local
  downwind coordinate; and
* an optional constant time-origin sensitivity for a resolved run.

The latter is never promoted to a physical calibration.  It is an explicit
observation-operator diagnostic: a clock offset may be applied only when an
event-time definition or synchronized sensor log supports it.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


def _finite_arrival_rows(workbook: Path) -> tuple[object, list[dict[str, float]]]:
    dataset = read_thorney_ti021_workbook(workbook)
    rows: list[dict[str, float]] = []
    for index, sensor in enumerate(dataset.dose_sensors):
        if sensor.arrival_time_s is None or sensor.dose_percent_s <= 0.0:
            continue
        rows.append({
            "workbook_row_index_zero_based": float(index),
            "local_y_m": float(sensor.local_y_m),
            "height_m": float(sensor.height_m),
            "observed_arrival_s": float(sensor.arrival_time_s),
        })
    return dataset, rows


def _summary(residual_s: np.ndarray) -> dict[str, float | None]:
    if residual_s.size == 0:
        return {"count": 0, "mean_s": None, "median_s": None,
                "mae_s": None, "rmse_s": None, "p90_abs_s": None}
    return {
        "count": int(residual_s.size),
        "mean_s": float(np.mean(residual_s)),
        "median_s": float(np.median(residual_s)),
        "mae_s": float(np.mean(np.abs(residual_s))),
        "rmse_s": float(np.sqrt(np.mean(residual_s * residual_s))),
        "p90_abs_s": float(np.percentile(np.abs(residual_s), 90.0)),
    }


def _run_time_offset_sensitivity(
    workbook: Path,
    run: Path,
    *,
    threshold_vol_percent: float,
    max_offset_s: float,
    offset_step_s: float,
) -> dict[str, object]:
    record = json.loads((run / "ti021-inner-transient.json").read_text(encoding="utf-8"))
    with np.load(run / "ti021-inner-receptor-traces.npz") as bundle:
        time_s = np.asarray(bundle["time_s"], dtype=float)
        traces = np.asarray(bundle["concentration_vol_percent"], dtype=float)
    dataset = read_thorney_ti021_workbook(workbook)
    indices = [int(item["workbook_row_index_zero_based"])
               for item in record["configuration"]["receptors"]]
    observed = np.asarray([
        np.nan if dataset.dose_sensors[index].arrival_time_s is None
        else dataset.dose_sensors[index].arrival_time_s
        for index in indices
    ], dtype=float)
    predicted = np.asarray([
        time_s[hit[0]] if (hit := np.flatnonzero(trace >= threshold_vol_percent)).size
        else np.nan
        for trace in traces
    ], dtype=float)
    valid = np.isfinite(observed) & np.isfinite(predicted)
    if not np.any(valid):
        return {"run": str(run.resolve()), "valid_pair_count": 0,
                "status": "no_comparable_pairs", "diagnostic_only": True}

    offsets = np.arange(0.0, max_offset_s + 0.5 * offset_step_s, offset_step_s)
    maes = np.asarray([
        np.mean(np.abs(predicted[valid] + offset - observed[valid]))
        for offset in offsets
    ])
    best_index = int(np.argmin(maes))
    best_offset = float(offsets[best_index])
    return {
        "run": str(run.resolve()),
        "threshold_vol_percent": float(threshold_vol_percent),
        "valid_pair_count": int(np.count_nonzero(valid)),
        "raw_mae_s": float(np.mean(np.abs(predicted[valid] - observed[valid]))),
        "best_nonnegative_offset_s": best_offset,
        "best_offset_mae_s": float(maes[best_index]),
        "offset_scan": [
            {"offset_s": float(offset), "mae_s": float(mae)}
            for offset, mae in zip(offsets, maes)
        ],
        "interpretation": (
            "A constant offset is an observation-clock sensitivity only. It is "
            "not a solver correction and must not be adopted without synchronized "
            "release and sensor timestamps."
        ),
        "diagnostic_only": True,
    }


def audit(
    workbook: Path,
    *,
    run: Path | None = None,
    threshold_vol_percent: float = 0.1,
    max_offset_s: float = 40.0,
    offset_step_s: float = 0.25,
) -> dict[str, object]:
    dataset, rows = _finite_arrival_rows(workbook)
    local_y = np.asarray([row["local_y_m"] for row in rows], dtype=float)
    observed = np.asarray([row["observed_arrival_s"] for row in rows], dtype=float)
    by_wind_speed: dict[str, object] = {}
    for speed in (dataset.workbook_wind_speed_m_s, 3.9):
        advective = local_y / float(speed)
        residual = observed - advective
        by_wind_speed[str(float(speed))] = {
            "wind_speed_m_s": float(speed),
            "residual_observed_minus_local_y_over_wind_s": _summary(residual),
            "rows": [
                {**row, "advective_time_s": float(t), "residual_s": float(r)}
                for row, t, r in zip(rows, advective, residual)
            ],
        }

    # A transparent linear fit is a descriptive transport sanity check, not a
    # fitted solver parameter.  Keep it separate from the fixed wind-speed
    # residuals above.
    design = np.column_stack((np.ones(local_y.size), local_y))
    coefficients, *_ = np.linalg.lstsq(design, observed, rcond=None)
    fitted = design @ coefficients
    fit_residual = observed - fitted
    slope = float(coefficients[1])
    linear_fit = {
        "intercept_s": float(coefficients[0]),
        "slope_s_per_m": slope,
        "effective_speed_m_s": None if slope <= 0.0 else float(1.0 / slope),
        "residual": _summary(fit_residual),
        "description": "descriptive observed-arrival versus local downwind distance fit",
    }
    result: dict[str, object] = {
        "schema": "degali.ti021-time-origin-audit.v1",
        "status": "diagnostic_only",
        "workbook": str(workbook.resolve()),
        "workbook_wind_speed_m_s": float(dataset.workbook_wind_speed_m_s),
        "positive_dose_arrival_row_count": len(rows),
        "advection_sanity_by_speed": by_wind_speed,
        "descriptive_linear_fit": linear_fit,
        "interpretation": (
            "The workbook arrival field is an observation operator, not a 0--10 s "
            "concentration trace. A positive residual indicates unrepresented "
            "release/transport/sensor delay; it is not permission to shift the "
            "resolved concentration field."
        ),
        "diagnostic_only": True,
    }
    if run is not None:
        result["resolved_run_time_offset_sensitivity"] = _run_time_offset_sensitivity(
            workbook, run, threshold_vol_percent=threshold_vol_percent,
            max_offset_s=max_offset_s, offset_step_s=offset_step_s,
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--run-directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold-vol-percent", type=float, default=0.1)
    parser.add_argument("--max-offset-s", type=float, default=40.0)
    parser.add_argument("--offset-step-s", type=float, default=0.25)
    args = parser.parse_args()
    result = audit(
        args.workbook,
        run=args.run_directory,
        threshold_vol_percent=args.threshold_vol_percent,
        max_offset_s=args.max_offset_s,
        offset_step_s=args.offset_step_s,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
