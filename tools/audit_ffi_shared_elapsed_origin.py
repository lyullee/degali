"""Audit source/sensor elapsed-time origins in the public FFI channel runs.

This is a descriptive audit only. It does not fit a clock offset from H2
concentrations, interpolate either history, or assert that independent clocks
remain synchronized after their common recorded zero.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


SCHEMA = "degali.ffi-shared-elapsed-origin-audit.v2"
LES_SOURCE_CHECK_WINDOW_S = (60.0, 120.0)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _number(value: str | None) -> float | None:
    if value is None or not value.strip():
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def audit_run(
    path: Path,
    case_id: str,
    trial_id: str,
    article_les_source_rate_g_s: float | None = None,
) -> dict[str, object]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = set(reader.fieldnames or ())
        required = {"flow time [s]", "h2 sensor time"}
        if not required <= columns:
            raise ValueError(f"{path.name}: missing time columns {sorted(required-columns)}")
        flow_rate_column = next(
            (name for name in ("mass flow meter 1 [g/s]", "mass flow meter 2 [g/s]") if name in columns),
            None,
        )
        if flow_rate_column is None:
            raise ValueError(f"{path.name}: no mass-flow meter column")
        flow_times: list[float] = []
        sensor_times: list[float] = []
        flow_rates: list[tuple[float, float]] = []
        row_count = 0
        for row in reader:
            row_count += 1
            flow_t = _number(row.get("flow time [s]"))
            sensor_t = _number(row.get("h2 sensor time"))
            flow = _number(row.get(flow_rate_column))
            if flow_t is not None:
                flow_times.append(flow_t)
                if flow is not None:
                    flow_rates.append((flow_t, flow))
            if sensor_t is not None:
                sensor_times.append(sensor_t)

    if len(flow_times) < 2 or len(sensor_times) < 2:
        raise ValueError(f"{path.name}: fewer than two timestamps on one of the clocks")
    if any(b < a for a, b in zip(flow_times, flow_times[1:])):
        raise ValueError(f"{path.name}: flow clock is not nondecreasing")
    if any(b < a for a, b in zip(sensor_times, sensor_times[1:])):
        raise ValueError(f"{path.name}: sensor clock is not nondecreasing")
    flow_steps = [b - a for a, b in zip(flow_times, flow_times[1:]) if b > a]
    sensor_steps = [b - a for a, b in zip(sensor_times, sensor_times[1:]) if b > 0.0]
    flow_values = [q for _, q in flow_rates]
    raw_integral_g = math.fsum(
        0.5 * (left_q + right_q) * (right_t - left_t)
        for (left_t, left_q), (right_t, right_q) in zip(flow_rates, flow_rates[1:])
    )
    positive_integral_g = math.fsum(
        0.5 * (max(left_q, 0.0) + max(right_q, 0.0)) * (right_t - left_t)
        for (left_t, left_q), (right_t, right_q) in zip(flow_rates, flow_rates[1:])
    )
    final_rate = flow_rates[-1][1] if flow_rates else None
    max_rate = max(flow_values) if flow_values else None
    check_start_s, check_end_s = LES_SOURCE_CHECK_WINDOW_S
    check_rates = [
        q for time_s, q in flow_rates if check_start_s <= time_s < check_end_s
    ]
    check_mean = statistics.mean(check_rates) if check_rates else None
    return {
        "article_case": case_id,
        "article_trial": trial_id,
        "filename": path.name,
        "sha256": _sha256(path),
        "row_count": row_count,
        "flow": {
            "column": flow_rate_column,
            "finite_clock_samples": len(flow_times),
            "first_elapsed_s": flow_times[0],
            "last_elapsed_s": flow_times[-1],
            "median_positive_sample_interval_s": statistics.median(flow_steps),
            "finite_rate_samples": len(flow_rates),
            "minimum_recorded_rate_g_s": min(flow_values) if flow_values else None,
            "maximum_recorded_rate_g_s": max_rate,
            "last_recorded_rate_g_s": final_rate,
            "last_rate_as_fraction_of_recorded_max": (
                final_rate / max_rate if final_rate is not None and max_rate and max_rate > 0.0 else None
            ),
            "raw_trapezoidal_integral_g": raw_integral_g,
            "positive_part_trapezoidal_integral_g": positive_integral_g,
            "mean_rate_60_to_120_s_g_s": check_mean,
            "mean_rate_60_to_120_s_sample_count": len(check_rates),
            "mean_rate_60_to_120_s_minus_article_les_rate_g_s": (
                check_mean - article_les_source_rate_g_s
                if check_mean is not None and article_les_source_rate_g_s is not None
                else None
            ),
            "mean_rate_60_to_120_s_relative_difference_from_article_les": (
                check_mean / article_les_source_rate_g_s - 1.0
                if check_mean is not None and article_les_source_rate_g_s is not None
                and article_les_source_rate_g_s > 0.0
                else None
            ),
        },
        "article_les_source_rate_g_s": article_les_source_rate_g_s,
        "sensor": {
            "finite_clock_samples": len(sensor_times),
            "first_elapsed_s": sensor_times[0],
            "last_elapsed_s": sensor_times[-1],
            "median_positive_sample_interval_s": statistics.median(sensor_steps),
        },
        "recorded_elapsed_zero_matches": flow_times[0] == 0.0 and sensor_times[0] == 0.0,
        "source_window_covers_sensor_window": flow_times[-1] >= sensor_times[-1],
        "source_log_ends_with_nonzero_recorded_rate": bool(final_rate is not None and final_rate != 0.0),
        "sensor_window_extends_beyond_source_window_s": max(0.0, sensor_times[-1] - flow_times[-1]),
        "interpretation": (
            "Both logged elapsed-time columns begin at zero, but their recorded coverage differs. "
            "The sensor history extends beyond the source-flow record and the latter's endpoint "
            "must be checked before the pair can represent a complete transient. No concentration "
            "lag, baseline, drift, or sensor-response parameter is fitted."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    registry = json.loads(args.registry.read_text(encoding="utf-8"))
    runs = []
    for case in registry["case_crosswalk"]:
        if case["article_case"] not in {"O1", "O2"}:
            continue
        for trial, filename in zip(case["article_trials"], case["dataset_files"]):
            matches = list(args.data_root.rglob(filename))
            if len(matches) != 1:
                raise ValueError(f"expected one {filename}; found {len(matches)}")
            runs.append(audit_run(
                matches[0], case["article_case"], trial,
                float(case["article_les_source_rate_g_s"]),
            ))
    zero_matches = all(run["recorded_elapsed_zero_matches"] for run in runs)
    source_covers_sensor = all(run["source_window_covers_sensor_window"] for run in runs)
    payload = {
        "schema": SCHEMA,
        "scope": "FFI obstacle-present O1/O2 experiments only",
        "status": (
            "elapsed_origins_mismatch" if not zero_matches else
            "shared_elapsed_zero_and_full_coverage" if source_covers_sensor else
            "shared_elapsed_zero_but_source_coverage_incomplete"
        ),
        "concentration_values_used_for_alignment": False,
        "common_absolute_clock_proven": False,
        "source_sensor_zero_matches_in_every_run": zero_matches,
        "source_flow_history_covers_sensor_history_in_every_run": source_covers_sensor,
        "article_les_source_rate_check_window_s": list(LES_SOURCE_CHECK_WINDOW_S),
        "run_count": len(runs),
        "runs": runs,
        "limits": [
            "Independent timestamps are retained at their native sampling rates.",
            "No lag, drift, baseline, or sensor-response parameter is fitted.",
            "A matching first zero is not evidence that independent acquisition clocks remain drift-free.",
            "This audit does not resolve obstacle geometry, facility boundary conditions, or calibration.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "run_count": payload["run_count"],
        "source_sensor_zero_matches_in_every_run": payload["source_sensor_zero_matches_in_every_run"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
