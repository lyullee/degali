"""Build a hash-pinned FFI obstacle/control observation baseline on 60--120 s.

This extracts directly comparable measured sensor means and source-flow
statistics from paired FFI trials. It is an observation reference only: it
does not score DEGALI, infer clock drift, or resolve obstacle geometry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from degali.validation.open_channel_h2_obstacle import (
    load_open_channel_h2_obstacle_benchmark,
)


SCHEMA = "degali.ffi-obstacle-control-window-baseline.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _window_mean(time: np.ndarray, values: np.ndarray, start: float, end: float) -> float:
    if time.ndim != 1 or values.shape != time.shape or len(time) < 2:
        raise ValueError("sensor time and value arrays must be aligned one-dimensional series")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(values)):
        raise ValueError("sensor series contains non-finite values")
    if np.any(np.diff(time) <= 0.0):
        raise ValueError("sensor clock must be strictly increasing")
    if time[0] > start or time[-1] < end or not end > start:
        raise ValueError("sensor record does not cover the requested averaging window")
    inside = (time > start) & (time < end)
    t = np.concatenate(([start], time[inside], [end]))
    y = np.interp(t, time, values)
    return float(np.trapezoid(y, t) / (end - start))


def _flow_window(run, start: float, end: float) -> dict[str, float | int]:
    time = np.asarray(run.flow_time_s, dtype=float)
    rate = np.asarray(run.mass_flow_g_s, dtype=float)
    if time[0] > start or time[-1] < end:
        raise ValueError("source-flow record does not cover the requested window")
    inside = (time > start) & (time < end)
    t = np.concatenate(([start], time[inside], [end]))
    q = np.interp(t, time, rate)
    integral = float(np.trapezoid(q, t))
    return {
        "mean_g_s": integral / (end - start),
        "integrated_mass_g": integral,
        "sample_count": int(np.count_nonzero((time >= start) & (time <= end))),
    }


def build_report(registry: str | Path, data_root: str | Path,
                 start_s: float = 60.0, end_s: float = 120.0) -> dict[str, object]:
    if not math.isfinite(start_s) or not math.isfinite(end_s) or end_s <= start_s:
        raise ValueError("window must be finite and have end > start")
    benchmark = load_open_channel_h2_obstacle_benchmark(registry, data_root)
    cases = {case.article_case: case for case in benchmark.cases}
    comparisons = (("A", "O1"), ("E", "O2"))
    sensor_rows = {sensor.sensor_id: sensor for sensor in benchmark.sensor_receptors}
    pair_reports = []
    for control_id, obstacle_id in comparisons:
        control, obstacle = cases[control_id], cases[obstacle_id]
        if len(control.runs) != 2 or len(obstacle.runs) != 2:
            raise ValueError(f"expected two paired repeats for {control_id}/{obstacle_id}")
        repeat_reports = []
        for index, (control_run, obstacle_run) in enumerate(zip(control.runs, obstacle.runs), start=1):
            cflow = _flow_window(control_run, start_s, end_s)
            oflow = _flow_window(obstacle_run, start_s, end_s)
            cmeans = {
                sid: _window_mean(control_run.sensor_time_s, values, start_s, end_s)
                for sid, values in control_run.sensor_percent.items()
            }
            omeans = {
                sid: _window_mean(obstacle_run.sensor_time_s, values, start_s, end_s)
                for sid, values in obstacle_run.sensor_percent.items()
            }
            if set(cmeans) != set(omeans) or set(cmeans) != set(sensor_rows):
                raise ValueError("control/obstacle sensor identities do not match the registry")
            deltas = {sid: omeans[sid] - cmeans[sid] for sid in cmeans}
            # A transparent source-rate normalization is provided only as a
            # sensitivity, not asserted as a validated linear correction.
            rate_ratio = float(oflow["mean_g_s"]) / float(cflow["mean_g_s"])
            rate_adjusted = {
                sid: omeans[sid] / rate_ratio - cmeans[sid] for sid in cmeans
            }
            placement_groups: dict[str, list[str]] = {}
            for sid, receptor in sensor_rows.items():
                placement_groups.setdefault(receptor.placement, []).append(sid)
            by_placement = {
                placement: {
                    "sensor_count": len(ids),
                    "mean_control_vol_pct": float(np.mean([cmeans[sid] for sid in ids])),
                    "mean_obstacle_vol_pct": float(np.mean([omeans[sid] for sid in ids])),
                    "mean_obstacle_minus_control_vol_pct": float(np.mean([deltas[sid] for sid in ids])),
                    "mean_source_rate_adjusted_delta_vol_pct": float(np.mean([rate_adjusted[sid] for sid in ids])),
                }
                for placement, ids in sorted(placement_groups.items())
            }
            repeat_reports.append({
                "repeat_index": index,
                "control_trial": control.trials[index - 1],
                "obstacle_trial": obstacle.trials[index - 1],
                "control_csv": control_run.source_path.name,
                "obstacle_csv": obstacle_run.source_path.name,
                "control_sha256": _sha256(control_run.source_path),
                "obstacle_sha256": _sha256(obstacle_run.source_path),
                "source_flow": {"control": cflow, "obstacle": oflow,
                                "obstacle_to_control_mean_rate_ratio": rate_ratio},
                "sensor_means_vol_pct": {
                    sid: {"control": cmeans[sid], "obstacle": omeans[sid],
                          "obstacle_minus_control": deltas[sid],
                          "source_rate_adjusted_obstacle_minus_control": rate_adjusted[sid],
                          "position_m": list(sensor_rows[sid].position_m),
                          "placement": sensor_rows[sid].placement}
                    for sid in sorted(cmeans, key=lambda value: int(value.split()[-1]))
                },
                "aggregate": {
                    "sensor_count": len(deltas),
                    "mean_control_vol_pct": float(np.mean(list(cmeans.values()))),
                    "mean_obstacle_vol_pct": float(np.mean(list(omeans.values()))),
                    "mean_obstacle_minus_control_vol_pct": float(np.mean(list(deltas.values()))),
                    "mean_source_rate_adjusted_delta_vol_pct": float(np.mean(list(rate_adjusted.values()))),
                    "max_obstacle_minus_control_sensor_id": max(deltas, key=deltas.get),
                    "min_obstacle_minus_control_sensor_id": min(deltas, key=deltas.get),
                    "by_sensor_placement": by_placement,
                },
            })
        pair_reports.append({"control_case": control_id, "obstacle_case": obstacle_id,
                             "repeats": repeat_reports})
    return {
        "schema": SCHEMA,
        "window": {"start_s": start_s, "end_s": end_s,
                   "sensor_mean_operator": "piecewise-linear time-weighted trapezoidal mean; exact window edges interpolated",
                   "source_operator": "piecewise-linear flow interpolation and trapezoidal integral on same elapsed-time window"},
        "scope": "measured FFI open-channel H2 sensor baseline, paired control/obstacle repeats",
        "clock_caveat": "source-flow and sensor elapsed clocks are independently recorded; shared zero is not proof of drift-free synchronization",
        "source_coverage": "both source and sensor records must cover the entire requested window for every repeat",
        "concentration_unit": "vol%",
        "source_rate_adjustment": "linear concentration/source-mean-rate rescaling is a sensitivity diagnostic only; raw deltas are primary",
        "geometry_quantitative_allowed": benchmark.geometry_quantitative_allowed,
        "obstacle_h2_validation_allowed": benchmark.obstacle_h2_validation_allowed,
        "degali_score_allowed": benchmark.degali_score_allowed,
        "gate_codes": list(benchmark.gate_codes),
        "promotion_allowed": False,
        "comparisons": pair_reports,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=Path("reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"))
    parser.add_argument("--data-root", type=Path, default=Path("reference/hydrogen/open_channel_usn/extracted"))
    parser.add_argument("--start", type=float, default=60.0)
    parser.add_argument("--end", type=float, default=120.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError(args.output)
    payload = build_report(args.registry, args.data_root, args.start, args.end)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(f"FFI paired baseline: {args.output.resolve()}")
    print(f"Window: {args.start:g}--{args.end:g} s; comparisons: {len(payload['comparisons'])}; gate remains withheld")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
