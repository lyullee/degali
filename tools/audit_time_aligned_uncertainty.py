"""Build a time-aligned PRESLHY E3.5 screen and deterministic input envelope.

The primary operator uses the central 20 s of each frozen sustained-flow
window. Flowmeter samples are selected on their recorded clock. Xensor samples
are selected on their recorded clock after the documented 18 s sampling-line
delay. A zero-delay result is retained only as an observation-operator
sensitivity. Input envelopes use independently recorded mean/maximum wind and
the pre-registered sustained-window peak source-rate bound; they are not
probability or confidence intervals.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import math
import platform
from collections import defaultdict
from pathlib import Path

import numpy as np

from degali.validation import nearfield, preslhy, spadeadam
from degali.validation.statistics import statistics


ROOT = Path(__file__).resolve().parents[1]
REDUCED = ROOT / "reference/preslhy/e35_reduced.json"
CONDITIONS = ROOT / "reference/preslhy/conditions.csv"
RAW = ROOT / "reference/preslhy/raw"
FFI = ROOT / "reference/spadeadam"
TRIALS = (10, 11, 12, 22, 23, 24, 25)
COMMON_WINDOW_S = 20.0
PRIMARY_DELAY_S = 18.0
DETECTION_FLOOR_VOL_PCT = 0.05


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _end_time(clock: np.ndarray, last: int) -> float:
    if last < len(clock):
        return float(clock[last])
    increments = np.diff(clock)
    positive = increments[increments > 0]
    step = float(np.median(positive)) if positive.size else 1.0
    return float(clock[last - 1] + step)


def _read_raw(trial: dict, delay_s: float) -> tuple[dict, list[dict]]:
    """Read one common clock window without changing the source workbook."""
    import openpyxl

    number = int(trial["trial"])
    paths = list(RAW.glob(f"trial_{number}_*alldata.xlsx"))
    if len(paths) != 1:
        raise FileNotFoundError(f"trial {number}: expected one raw workbook, got {paths}")
    path = paths[0]
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    flow_rows = list(book["Flowmeter"].iter_rows(values_only=True))
    gas_rows = list(book["Xensor"].iter_rows(values_only=True))
    flow_header, gas_header = flow_rows[0], gas_rows[0]
    flow_col = next(i for i, heading in enumerate(flow_header)
                    if heading and "MassFlow" in str(heading))
    flow_clock = preslhy._clock_seconds([row[0] for row in flow_rows[1:]])
    gas_clock = preslhy._clock_seconds([row[0] for row in gas_rows[1:]])
    if not (np.all(np.isfinite(flow_clock)) and np.all(np.isfinite(gas_clock))):
        raise ValueError(f"trial {number}: recorded clocks are not fully parseable")
    first, last = map(int, trial["window"])
    sustained_start = float(flow_clock[first])
    sustained_end = _end_time(flow_clock, last)
    duration = sustained_end - sustained_start
    if duration + 1e-9 < COMMON_WINDOW_S:
        raise ValueError(f"trial {number}: sustained window {duration:g} s is too short")
    source_start = sustained_start + 0.5 * (duration - COMMON_WINDOW_S)
    source_end = source_start + COMMON_WINDOW_S
    flow_mask = (flow_clock >= source_start) & (flow_clock < source_end)
    flow = np.asarray([
        float(row[flow_col]) if row[flow_col] is not None else float("nan")
        for row in flow_rows[1:]
    ])
    rate = float(np.nanmean(flow[flow_mask]))
    if not math.isfinite(rate) or rate <= 0:
        raise ValueError(f"trial {number}: invalid common-window source rate {rate}")
    gas_start, gas_end = source_start + delay_s, source_end + delay_s
    gas_mask = (gas_clock >= gas_start) & (gas_clock < gas_end)
    columns = {
        match.group(1): index
        for index, heading in enumerate(gas_header)
        if heading and (match := preslhy._SERIAL.search(str(heading)))
    }
    sensors = []
    for sensor in trial["sensors"]:
        serial = sensor["serial"]
        if serial not in columns:
            continue
        values = np.asarray([
            float(row[columns[serial]]) if row[columns[serial]] is not None else float("nan")
            for row in gas_rows[1:]
        ])
        selected = values[gas_mask]
        finite = selected[np.isfinite(selected)]
        if not finite.size:
            continue
        sensors.append({
            "trial": number, "delay_s": delay_s, "serial": serial,
            "x_m": float(sensor["x"]), "y_m": float(sensor["y"]),
            "z_m": float(sensor["z"]), "samples": int(finite.size),
            "observed_mean_vol_pct": float(np.mean(finite)),
            "observed_sd_vol_pct": float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0,
            "observed_peak_vol_pct": float(np.max(finite)),
        })
    window = {
        "trial": number, "workbook": path.name,
        "sustained_start_clock_s": sustained_start,
        "sustained_end_clock_s": sustained_end,
        "sustained_duration_s": duration,
        "common_source_start_clock_s": source_start,
        "common_source_end_clock_s": source_end,
        "gas_delay_s": delay_s, "gas_start_clock_s": gas_start,
        "gas_end_clock_s": gas_end, "flow_samples": int(np.sum(flow_mask)),
        "gas_samples": int(np.sum(gas_mask)), "flow_mean_g_s": rate,
        "flow_sd_g_s": float(np.nanstd(flow[flow_mask], ddof=1)),
        "flow_cv": float(np.nanstd(flow[flow_mask], ddof=1) / rate),
        "sensor_count": len(sensors),
    }
    return window, sensors


def _trajectory(trial: dict, *, rate_g_s: float, wind_m_s: float):
    jet, initial = nearfield.hydrogen_jet(
        rate=rate_g_s / 1000.0,
        diameter=float(trial["orifice_mm"]) / 1000.0,
        wind=wind_m_s, height=float(trial["release_height_m"]),
        ambient_temperature=float(trial["T_C"]) + 273.15,
        relative_humidity=float(trial["RH_pct"]),
        storage_pressure_barg=float(trial["tanker_barg"]),
        storage_temperature=nearfield._liquid_source_temperature(trial, "tank_saturation"),
        wind_reference_height=float(trial["wind_ref_m"]),
        roughness=0.001, stability="D", averaging=COMMON_WINDOW_S,
        corrections=True, ground_effect=False,
    )
    result = nearfield.Trajectory(
        jet.th.table, jet.run(initial, distmx=nearfield.STEP, smax=nearfield.REACH).rows
    )
    if not result.ok:
        raise RuntimeError(f"trial {trial['trial']} integration failed")
    return result


def _arcs(sensors: list[dict], trajectory) -> list[dict]:
    grouped = defaultdict(list)
    for row in sensors:
        grouped[round(float(row["x_m"]), 3)].append(row)
    out = []
    for x, rows in sorted(grouped.items()):
        observed = max(float(row["observed_mean_vol_pct"]) for row in rows)
        predicted = max(trajectory.concentration_at(
            x, float(row["y_m"]), float(row["z_m"])) for row in rows)
        if observed <= DETECTION_FLOOR_VOL_PCT or predicted <= 0:
            continue
        out.append({"x_m": x, "n_sensors": len(rows),
                    "observed_arc_mean_max_vol_pct": observed,
                    "predicted_arc_max_vol_pct": predicted,
                    "observed_over_predicted": observed / predicted})
    return out


def _metric(rows: list[dict]) -> dict:
    score = statistics(
        [row["observed_arc_mean_max_vol_pct"] for row in rows],
        [row["predicted_arc_max_vol_pct"] for row in rows],
    )
    return dataclasses.asdict(score)


def _lfl_reach(trajectory, height_m: float = 1.5) -> tuple[float | None, str]:
    xs = np.geomspace(0.05, nearfield.REACH, 1000)
    points = [(float(x), trajectory.concentration_at(float(x), 0.0, height_m))
              for x in xs if trajectory.at(float(x)) is not None]
    above = [i for i, (_, concentration) in enumerate(points) if concentration >= 4.0]
    if not above:
        return None, "never_reaches_lfl"
    index = above[-1]
    if index == len(points) - 1:
        return None, "right_censored_at_40m"
    x0, c0 = points[index]
    x1, c1 = points[index + 1]
    return float(x0 + (4.0 - c0) * (x1 - x0) / (c1 - c0)), "resolved"


def _ffi_test6_envelope() -> list[dict]:
    trial = next(item for item in spadeadam.load(FFI) if item.test == 6)
    rows = []
    for label, wind in (("low_mast", trial.wind_low),
                        ("mast_midpoint", 0.5 * (trial.wind_low + trial.wind_high)),
                        ("high_mast", trial.wind_high)):
        jet, initial = nearfield.hydrogen_jet(
            rate=trial.rate, diameter=trial.orifice, wind=wind,
            height=spadeadam.RELEASE_HEIGHT,
            ambient_temperature=spadeadam.AMBIENT_TEMPERATURE,
            relative_humidity=spadeadam.RELATIVE_HUMIDITY,
            storage_pressure_barg=trial.line_pressure,
            roughness=0.001, stability="D", averaging=60.0,
            wind_reference_height=10.0, corrections=True, ground_effect=False,
        )
        trajectory = nearfield.Trajectory(
            jet.th.table, jet.run(initial, distmx=nearfield.STEP, smax=250.0).rows
        )
        # Same fixed-height definition used in the application example.
        xs = np.geomspace(0.05, 250.0, 1400)
        points = [(float(x), trajectory.concentration_at(float(x), 0.0, 1.5))
                  for x in xs if trajectory.at(float(x)) is not None]
        above = [i for i, (_, c) in enumerate(points) if c >= 4.0]
        reach = None
        status = "never_reaches_lfl"
        if above and above[-1] < len(points) - 1:
            i = above[-1]
            x0, c0 = points[i]
            x1, c1 = points[i + 1]
            reach = x0 + (4.0 - c0) * (x1 - x0) / (c1 - c0)
            status = "resolved"
        elif above:
            status = "right_censored_at_250m"
        rows.append({"case": "FFI_Test_6", "wind_label": label,
                     "wind_m_s": wind, "rate_kg_s": trial.rate,
                     "receptor_height_m": 1.5, "lfl_distance_m": reach,
                     "status": status,
                     "scope": "wind-only deterministic sensitivity; source-state uncertainty unavailable"})
    return rows


def _summary_svg(path: Path, trial_summary: list[dict], ffi_rows: list[dict]) -> None:
    width, height = 1000, 560
    left, top, plot_w, plot_h = 90, 65, 820, 380
    trials = [row["trial"] for row in trial_summary]
    finite = [float(row["lfl_upper_reporting_m"]) for row in trial_summary
              if row["lfl_upper_reporting_m"] not in (None, "")]
    ymax = max([40.0] + finite)
    def sx(index): return left + (index + 0.5) / len(trials) * plot_w
    def sy(value): return top + plot_h - float(value) / ymax * plot_h
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:20px;font-weight:bold}.label{font-size:14px}</style>',
             '<text x="90" y="34" class="title">E3.5 deterministic source/wind envelope: 4% reach at z = 1.5 m</text>',
             f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#222"/>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#222"/>']
    for index, row in enumerate(trial_summary):
        x = sx(index)
        lo, hi = row["lfl_min_m"], row["lfl_upper_reporting_m"]
        if lo not in (None, "") and hi not in (None, ""):
            parts.append(f'<line x1="{x:.1f}" y1="{sy(lo):.1f}" x2="{x:.1f}" y2="{sy(hi):.1f}" stroke="#0069aa" stroke-width="6"/>')
            parts.append(f'<circle cx="{x:.1f}" cy="{sy(lo):.1f}" r="5" fill="#174263"/>')
            parts.append(f'<circle cx="{x:.1f}" cy="{sy(hi):.1f}" r="5" fill="#174263"/>')
        parts.append(f'<text x="{x-10:.1f}" y="{top+plot_h+25}" class="label">{row["trial"]}</text>')
    parts += ['<text x="90" y="510" class="label">Bars use common-20-s mean / peak-rate bound and mean / maximum wind. Censored cases remain flagged in CSV.</text>',
              '</svg>']
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    reduced = {int(row["trial"]): row for row in
               json.loads(REDUCED.read_text(encoding="utf-8"))["trials"]}
    with CONDITIONS.open(encoding="utf-8", newline="") as handle:
        conditions = {int(row["trial"]): row for row in csv.DictReader(handle)}

    windows, sensor_rows, primary_arcs, scenario_rows = [], [], [], []
    scenario_arc_rows = []
    trial_metrics = []
    trial_envelopes = []
    for number in TRIALS:
        trial = reduced[number]
        delay_data = {}
        for delay in (PRIMARY_DELAY_S, 0.0):
            window, sensors = _read_raw(trial, delay)
            windows.append(window)
            sensor_rows.extend(sensors)
            delay_data[delay] = (window, sensors)
        primary_window, primary_sensors = delay_data[PRIMARY_DELAY_S]
        wind_mean = float(conditions[number]["wind_mean_ms"])
        wind_max = float(conditions[number]["wind_max_ms"])
        rate_mean = float(primary_window["flow_mean_g_s"])
        rate_peak = float(trial["flow_peak_gs"])

        primary_trajectory = _trajectory(trial, rate_g_s=rate_mean, wind_m_s=wind_mean)
        arcs = _arcs(primary_sensors, primary_trajectory)
        for row in arcs:
            primary_arcs.append({"trial": number, "window_s": COMMON_WINDOW_S,
                                 "gas_delay_s": PRIMARY_DELAY_S, **row})
        score = _metric(arcs)
        trial_metrics.append({"trial": number, "arc_count": len(arcs), **score,
                              "flow_mean_g_s": rate_mean, "wind_mean_m_s": wind_mean})

        reaches = []
        censored = 0
        for rate_label, rate in (("common_20s_mean", rate_mean),
                                 ("sustained_window_peak_bound", rate_peak)):
            for wind_label, wind in (("reported_mean", wind_mean),
                                     ("reported_max", wind_max)):
                trajectory = _trajectory(trial, rate_g_s=rate, wind_m_s=wind)
                reach, status = _lfl_reach(trajectory)
                if reach is not None:
                    reaches.append(reach)
                else:
                    censored += 1
                for delay in (PRIMARY_DELAY_S, 0.0):
                    scenario_arcs = _arcs(delay_data[delay][1], trajectory)
                    score_row = _metric(scenario_arcs)
                    for arc in scenario_arcs:
                        scenario_arc_rows.append({
                            "trial": number, "rate_label": rate_label,
                            "wind_label": wind_label, "gas_delay_s": delay,
                            **arc,
                        })
                    scenario_rows.append({
                        "trial": number, "rate_label": rate_label, "rate_g_s": rate,
                        "wind_label": wind_label, "wind_m_s": wind,
                        "gas_delay_s": delay, "arc_count": len(scenario_arcs),
                        **score_row, "lfl_distance_at_1p5m_m": reach,
                        "lfl_status": status,
                    })
        upper_censored = censored > 0
        trial_envelopes.append({
            "trial": number, "rate_min_g_s": min(rate_mean, rate_peak),
            "rate_max_g_s": max(rate_mean, rate_peak),
            "wind_min_m_s": min(wind_mean, wind_max),
            "wind_max_m_s": max(wind_mean, wind_max),
            "lfl_min_m": min(reaches) if reaches else None,
            "lfl_max_resolved_m": max(reaches) if reaches else None,
            "lfl_upper_reporting_m": nearfield.REACH if upper_censored else (
                max(reaches) if reaches else None
            ),
            "lfl_upper_status": "greater_than_40m" if upper_censored else "resolved",
            "right_or_left_censored_scenarios": censored,
            "scope": "deterministic evidence envelope; not a probability interval",
        })
        print(f"trial {number} complete", flush=True)

    aggregate = _metric(primary_arcs)
    aggregate_scenarios = []
    for key in sorted({(row["rate_label"], row["wind_label"], row["gas_delay_s"])
                       for row in scenario_arc_rows}):
        rate_label, wind_label, delay = key
        rows = [row for row in scenario_arc_rows if
                (row["rate_label"], row["wind_label"], row["gas_delay_s"]) == key]
        aggregate_scenarios.append({
            "rate_label": rate_label, "wind_label": wind_label,
            "gas_delay_s": delay, "trial_count": len({row["trial"] for row in rows}),
            "arc_count": len(rows), **_metric(rows),
        })
    ffi_rows = _ffi_test6_envelope()
    write_csv(out / "time_alignment_windows.csv", windows)
    write_csv(out / "time_aligned_sensors.csv", sensor_rows)
    write_csv(out / "time_aligned_arcs.csv", primary_arcs)
    write_csv(out / "time_aligned_trial_metrics.csv", trial_metrics)
    write_csv(out / "uncertainty_scenarios.csv", scenario_rows)
    write_csv(out / "uncertainty_scenario_arcs.csv", scenario_arc_rows)
    write_csv(out / "uncertainty_aggregate_scenarios.csv", aggregate_scenarios)
    write_csv(out / "uncertainty_trial_summary.csv", trial_envelopes)
    write_csv(out / "ffi_test6_wind_distance_envelope.csv", ffi_rows)
    _summary_svg(out / "e35_lfl_uncertainty_envelope.svg", trial_envelopes, ffi_rows)

    manifest = {
        "primary_time_operator": {
            "source_window": "central 20 s of frozen sustained-flow window",
            "flow": "recorded Flowmeter clock; arithmetic mean over same 20 s",
            "concentration": "recorded Xensor clock; arithmetic mean over source window shifted +18 s",
            "spatial": "predict every physical sensor, then take distance-wise arc maximum",
            "detection_floor_vol_pct": DETECTION_FLOOR_VOL_PCT,
        },
        "primary_aggregate_statistics": aggregate,
        "aggregate_sensitivity_file": "uncertainty_aggregate_scenarios.csv",
        "included_trials": list(TRIALS),
        "exclusions": "Trials 20 and 21 are absent because their raw workbooks are not present in this transfer; Trial 24 is retained because its sustained window exceeds 20 s.",
        "uncertainty_interpretation": {
            "source": "20 s mean to pre-registered sustained-window peak-rate bound",
            "wind": "reported campaign mean to reported maximum; gust maximum excluded from steady envelope",
            "delay": "18 s documented correction versus zero-delay operator sensitivity",
            "not_included": ["sensor calibration uncertainty", "source pressure/temperature uncertainty",
                             "wind-direction/meander uncertainty", "model-form probability",
                             "parameter covariance"],
            "label": "deterministic evidence/sensitivity envelope, not 95% confidence or credible interval",
        },
        "ffi_test6": "wind-only distance sensitivity; no defensible atmospheric source-state uncertainty is available",
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "platform": platform.platform()},
        "input_sha256": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
            for path in [REDUCED, CONDITIONS, ROOT / "src/degali/validation/preslhy.py",
                         ROOT / "src/degali/validation/nearfield.py", Path(__file__),
                         *[next(RAW.glob(f"trial_{n}_*alldata.xlsx")) for n in TRIALS]]
        },
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(out), "aggregate": aggregate,
                      "ffi_test6": ffi_rows}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
