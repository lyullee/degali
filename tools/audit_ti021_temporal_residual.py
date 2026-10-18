"""Audit TI021 timing and concentration-shape residuals without fitting them.

The workbook supplies dose plus arrival/departure summaries for the full sensor
array.  Eight independently raster-digitized report plots additionally supply
trace shape.  This audit uses those records only after a candidate has been
run; it never changes a model parameter or rescales an observed trace.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]


def _summary(values: list[float]) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if not data.size:
        return {"count": 0, "mean": None, "median": None,
                "mean_absolute": None, "p90_absolute": None}
    return {
        "count": int(data.size),
        "mean": float(np.mean(data)),
        "median": float(np.median(data)),
        "mean_absolute": float(np.mean(np.abs(data))),
        "p90_absolute": float(np.percentile(np.abs(data), 90.0)),
    }


def _relative_window(time_s: np.ndarray, trace: np.ndarray, fraction: float):
    peak = float(np.max(trace))
    if peak <= 0.0:
        return None
    hits = np.flatnonzero(trace >= fraction * peak)
    if not hits.size:
        return None
    return float(time_s[hits[0]]), float(time_s[hits[-1]])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True,
                        help="candidate JSON; the adjacent NPZ is required")
    parser.add_argument(
        "--digitized",
        type=Path,
        default=ROOT / "outputs/ti021-gas-traces-digitized-strict33p597-2026-10-17.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--relative-threshold", type=float, default=0.05)
    args = parser.parse_args()

    if not 0.0 < args.relative_threshold < 1.0:
        raise ValueError("relative threshold must lie strictly between zero and one")
    report = json.loads(args.candidate.read_text(encoding="utf-8"))
    if report.get("observation_species_mode") != "released_mixture":
        raise ValueError("temporal production audit requires released_mixture mode")
    candidate_npz = np.load(args.candidate.with_suffix(".npz"))
    time_s = np.asarray(candidate_npz["time_s"], dtype=float)
    concentration = np.asarray(
        candidate_npz["concentration_vol_percent"], dtype=float
    )
    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")

    arrival_errors: list[float] = []
    departure_errors: list[float] = []
    duration_errors: list[float] = []
    missed = 0
    timing_rows = []
    for index, sensor in enumerate(dataset.dose_sensors):
        if sensor.dose_percent_s <= 0.0 or sensor.arrival_time_s is None:
            continue
        window = _relative_window(
            time_s, concentration[index], args.relative_threshold
        )
        if window is None:
            missed += 1
            continue
        arrival, departure = window
        arrival_error = arrival - float(sensor.arrival_time_s)
        arrival_errors.append(arrival_error)
        departure_error = duration_error = None
        if sensor.departure_time_s is not None:
            departure_error = departure - float(sensor.departure_time_s)
            duration_error = (
                (departure - arrival)
                - (float(sensor.departure_time_s) - float(sensor.arrival_time_s))
            )
            departure_errors.append(departure_error)
            duration_errors.append(duration_error)
        timing_rows.append({
            "sensor_index_one_based": index + 1,
            "local_x_m": float(sensor.local_x_m),
            "local_y_m": float(sensor.local_y_m),
            "height_m": float(sensor.height_m),
            "observed_arrival_s": float(sensor.arrival_time_s),
            "predicted_relative_threshold_arrival_s": arrival,
            "arrival_error_s": arrival_error,
            "observed_departure_s": (
                None if sensor.departure_time_s is None
                else float(sensor.departure_time_s)
            ),
            "predicted_relative_threshold_departure_s": departure,
            "departure_error_s": departure_error,
            "duration_error_s": duration_error,
        })

    digitized_report = json.loads(args.digitized.read_text(encoding="utf-8"))
    digitized_npz = np.load(args.digitized.with_suffix(".npz"))
    trace_time = np.asarray(digitized_npz["time_s"], dtype=float)
    observed_traces = np.asarray(
        digitized_npz["observed_concentration_vol_percent"], dtype=float
    )
    selected_indices = np.asarray(
        digitized_npz["workbook_sensor_index_one_based"], dtype=int
    ) - 1
    trace_rows = []
    peak_time_errors = []
    log_peak_ratios = []
    equivalent_duration_ratios = []
    for row_index, sensor_index in enumerate(selected_indices):
        observed = observed_traces[row_index]
        predicted = np.interp(trace_time, time_s, concentration[sensor_index])
        observed_peak = float(np.max(observed))
        predicted_peak = float(np.max(predicted))
        observed_dose = float(np.trapezoid(observed, trace_time))
        predicted_dose = float(np.trapezoid(predicted, trace_time))
        observed_peak_time = float(trace_time[int(np.argmax(observed))])
        predicted_peak_time = float(trace_time[int(np.argmax(predicted))])
        peak_ratio = predicted_peak / observed_peak if observed_peak > 0.0 else None
        observed_equivalent_duration = (
            observed_dose / observed_peak if observed_peak > 0.0 else None
        )
        predicted_equivalent_duration = (
            predicted_dose / predicted_peak if predicted_peak > 0.0 else None
        )
        duration_ratio = (
            predicted_equivalent_duration / observed_equivalent_duration
            if observed_equivalent_duration and predicted_equivalent_duration is not None
            else None
        )
        peak_time_errors.append(predicted_peak_time - observed_peak_time)
        if peak_ratio and peak_ratio > 0.0:
            log_peak_ratios.append(float(np.log(peak_ratio)))
        if duration_ratio is not None:
            equivalent_duration_ratios.append(duration_ratio)
        source_row = digitized_report["sensor_rows"][row_index]
        trace_rows.append({
            "gas_id": source_row["gas_id"],
            "sensor_index_one_based": int(sensor_index + 1),
            "height_m": float(dataset.dose_sensors[sensor_index].height_m),
            "observed_peak_vol_percent": observed_peak,
            "predicted_peak_vol_percent": predicted_peak,
            "predicted_to_observed_peak_ratio": peak_ratio,
            "observed_peak_time_s": observed_peak_time,
            "predicted_peak_time_s": predicted_peak_time,
            "peak_time_error_s": predicted_peak_time - observed_peak_time,
            "observed_equivalent_duration_s": observed_equivalent_duration,
            "predicted_equivalent_duration_s": predicted_equivalent_duration,
            "predicted_to_observed_equivalent_duration_ratio": duration_ratio,
        })

    result = {
        "schema": "degali.ti021-temporal-residual-audit.v1",
        "status": "diagnostic_only_not_a_fitted_model",
        "candidate": str(args.candidate.resolve()),
        "candidate_dose_mape_percent": float(report["dose_mape_percent"]),
        "relative_threshold_fraction_of_local_peak": args.relative_threshold,
        "relative_threshold_is_diagnostic_not_confirmed_workbook_operator": True,
        "positive_timed_sensor_count": int(len(timing_rows) + missed),
        "missed_relative_threshold_sensor_count": int(missed),
        "arrival_error_s": _summary(arrival_errors),
        "departure_error_s": _summary(departure_errors),
        "duration_error_s": _summary(duration_errors),
        "digitized_trace_count": int(len(trace_rows)),
        "digitized_peak_time_error_s": _summary(peak_time_errors),
        "digitized_geometric_mean_peak_ratio": (
            None if not log_peak_ratios else float(np.exp(np.mean(log_peak_ratios)))
        ),
        "digitized_median_equivalent_duration_ratio": (
            None if not equivalent_duration_ratios
            else float(np.median(equivalent_duration_ratios))
        ),
        "interpretation": (
            "Negative timing errors mean the candidate is early; an equivalent-"
            "duration ratio below one means its concentration pulse is too narrow."
        ),
        "timing_rows": timing_rows,
        "digitized_trace_rows": trace_rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {"timing_rows", "digitized_trace_rows"}}, indent=2))


if __name__ == "__main__":
    main()
