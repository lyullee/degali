"""Paired, no-fit LH2 comparison and an explicitly conditional transfer screen.

Run from the repository root with ``PYTHONPATH=src``.  Derived files only are
written to --output-dir.  The legacy comparator is the local
``corrections=False`` reconstruction, *not* an independent DEGADIS executable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics as statlib
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from degali.validation import nearfield, spadeadam
from degali.validation.statistics import statistics

ROOT = Path(__file__).resolve().parents[1]
E35 = ROOT / "reference/preslhy/e35_reduced.json"
FFI = ROOT / "reference/spadeadam"
VARIANTS = (("historical_reconstruction", False), ("corrected", True))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"empty output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_e35(trial: dict, corrected: bool, repeats: int):
    timings = []
    trajectory = None
    for _ in range(repeats):
        start = time.perf_counter()
        jet, initial = nearfield.hydrogen_jet(
            rate=trial["flow_mean_gs"] / 1000.0,
            diameter=trial["orifice_mm"] / 1000.0,
            wind=trial["wind_ms"],
            height=trial["release_height_m"],
            ambient_temperature=trial["T_C"] + 273.15,
            relative_humidity=trial["RH_pct"],
            storage_pressure_barg=trial["tanker_barg"],
            storage_temperature=nearfield._liquid_source_temperature(
                trial, "tank_saturation"
            ),
            wind_reference_height=trial["wind_ref_m"],
            corrections=corrected,
            ground_effect=False,
        )
        candidate = nearfield.Trajectory(
            jet.th.table,
            jet.run(initial, distmx=nearfield.STEP, smax=nearfield.REACH).rows,
        )
        if not candidate.ok:
            raise RuntimeError(f"E3.5 trial {trial['trial']} integration failed")
        # Include the exact sensor/arc operator in the timing.
        for sensor in trial["sensors"]:
            candidate.concentration_at(sensor["x"], sensor["y"], sensor["z"])
        timings.append(time.perf_counter() - start)
        trajectory = candidate
    return trajectory, timings


def run_ffi(trial, corrected: bool, repeats: int):
    timings = []
    trajectory = None
    for _ in range(repeats):
        start = time.perf_counter()
        jet, initial = nearfield.hydrogen_jet(
            rate=trial.rate,
            diameter=trial.orifice,
            wind=trial.wind_low,
            height=spadeadam.RELEASE_HEIGHT,
            ambient_temperature=spadeadam.AMBIENT_TEMPERATURE,
            relative_humidity=spadeadam.RELATIVE_HUMIDITY,
            storage_pressure_barg=trial.line_pressure,
            roughness=0.001,
            wind_reference_height=10.0,
            corrections=corrected,
            ground_effect=False,
        )
        candidate = nearfield.Trajectory(
            jet.th.table, jet.run(initial, distmx=nearfield.STEP, smax=250.0).rows
        )
        if not candidate.ok:
            raise RuntimeError(f"FFI test {trial.test} integration failed")
        for radius in trial.radii():
            for height in trial.arc(radius):
                candidate.concentration_at(radius, 0.0, height)
        timings.append(time.perf_counter() - start)
        trajectory = candidate
    return trajectory, timings


def bootstrap_mg(ratios: list[float], rng, draws: int):
    logs = np.log(np.asarray(ratios, dtype=float))
    indexes = rng.integers(0, len(logs), size=(draws, len(logs)))
    samples = np.exp(np.mean(logs[indexes], axis=1))
    return tuple(float(v) for v in np.quantile(samples, [0.025, 0.975]))


def metric(rows: list[dict], variant: str) -> dict:
    pairs = [(r["observed_vol_pct"], r[f"{variant}_vol_pct"]) for r in rows]
    score = statistics([p[0] for p in pairs], [p[1] for p in pairs])
    return {"n_arcs": score.n, "mg": score.mg, "vg": score.vg,
            "fac2": score.fac2, "fb": score.fb, "nmse": score.nmse}


def cluster_interval(rows: list[dict], variant: str, rng, draws: int):
    groups = defaultdict(list)
    for row in rows:
        groups[row["test_id"]].append(row)
    keys = sorted(groups)
    mg, vg, fac2 = [], [], []
    for _ in range(draws):
        sample = [row for key in rng.choice(keys, len(keys), replace=True)
                  for row in groups[key]]
        result = metric(sample, variant)
        mg.append(result["mg"])
        vg.append(result["vg"])
        fac2.append(result["fac2"])
    return (float(np.quantile(mg, 0.025)), float(np.quantile(mg, 0.975)),
            float(np.quantile(vg, 0.025)), float(np.quantile(vg, 0.975)),
            float(np.quantile(fac2, 0.025)), float(np.quantile(fac2, 0.975)))


def lfl_reach(trajectory, *, height: float, limit: float = 250.0):
    """First sustained downwind exit below 4 vol% on a fixed-height axis."""
    xs = np.geomspace(0.05, limit, 1200)
    points = [(float(x), trajectory.concentration_at(float(x), 0.0, height))
              for x in xs if trajectory.at(float(x)) is not None]
    above = [i for i, (_, concentration) in enumerate(points)
             if concentration >= 4.0]
    if not above:
        return None, "never_reaches_lfl"
    last = above[-1]
    if last == len(points) - 1:
        return None, "right_censored"
    left_x, left_c = points[last]
    right_x, right_c = points[last + 1]
    crossing = left_x + (4.0 - left_c) * (right_x - left_x) / (right_c - left_c)
    return crossing, "resolved"


def _svg_plot(path: Path, rows: list[dict], case_rows: list[dict],
              trial_rows: list[dict]) -> None:
    """Compact, dependency-free manuscript figures."""
    width, height = 980, 560
    left, top, pw, ph = 90, 65, 820, 390
    def sx(x): return left + (x / 100.0) * pw
    def sy(y): return top + ph - min(max(y, -0.5), 3.5) / 4.0 * ph
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-weight:bold;font-size:20px}.label{font-size:14px}</style>',
             '<text x="90" y="34" class="title">FFI arc-maximum error, matched spatial operator</text>',
             f'<line x1="{left}" y1="{sy(1)}" x2="{left+pw}" y2="{sy(1)}" stroke="#777" stroke-dasharray="5 5"/>',
             f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#222"/>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#222"/>']
    for tick in (-0.5, 0, 1, 2, 3, 3.5):
        parts.append(f'<text x="48" y="{sy(tick)+5:.1f}" class="label">{tick:g}</text>')
    for tick in (0, 30, 50, 100):
        parts.append(f'<text x="{sx(tick)-10:.1f}" y="{top+ph+25}" class="label">{tick}</text>')
    for row in rows:
        if row["dataset"] != "FFI":
            continue
        for variant, color, offset in (("historical_reconstruction", "#777777", -6), ("corrected", "#0069aa", 6)):
            ratio = row["observed_vol_pct"] / row[f"{variant}_vol_pct"]
            x, y = sx(row["x_m"]) + offset, sy(math.log10(ratio) + 1)
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6" fill="{color}"/>')
            parts.append(f'<text x="{x+7:.1f}" y="{y-7:.1f}" class="label">T{row["test_id"]}</text>')
    parts += [f'<text x="{left+pw/2-70}" y="{height-48}" class="label">Arc radius (m)</text>',
              '<text x="90" y="520" class="label">Y axis: 1 + log10(observed / predicted); 1 = parity. Grey = historical, blue = corrected.</text>',
              '</svg>']
    path.write_text("\n".join(parts), encoding="utf-8")

    left, top, pw, ph = 90, 65, 820, 390
    max_x = max(row["x_m"] for row in case_rows)
    max_c = max(max(row["historical_reconstruction_vol_pct"], row["corrected_vol_pct"])
                for row in case_rows)
    max_c = max(8.0, min(100.0, math.ceil(max_c / 5.0) * 5.0))
    def px(x): return left + (x / max_x) * pw
    def py(y): return top + ph - min(max(y, 0.0), max_c) / max_c * ph
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-weight:bold;font-size:20px}.label{font-size:14px}</style>',
             '<text x="90" y="34" class="title">FFI Test 6 transfer analogue: 1.5 m axis, conditional screen</text>',
             f'<line x1="{left}" y1="{py(4):.1f}" x2="{left+pw}" y2="{py(4):.1f}" stroke="#aa3333" stroke-dasharray="6 4"/>',
             f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#222"/>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#222"/>']
    for variant, color in (("historical_reconstruction", "#777777"), ("corrected", "#0069aa")):
        points = " ".join(f'{px(r["x_m"]):.1f},{py(r[f"{variant}_vol_pct"]):.1f}' for r in case_rows)
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="3"/>')
    parts += [f'<text x="{left+pw/2-70}" y="{height-48}" class="label">Downwind distance (m)</text>',
              f'<text x="90" y="520" class="label">Hydrogen (vol %); dashed red = 4% LFL. Grey = historical, blue = corrected.</text>',
              '</svg>']
    path.with_name("transfer_lfl_screen.svg").write_text("\n".join(parts), encoding="utf-8")

    selected = [row for row in trial_rows if row["dataset"] == "E35"]
    trial_ids = sorted({row["test_id"] for row in selected})
    left, top, pw, ph = 90, 65, 820, 390
    ymin, ymax = -0.8, 1.65
    def tx(test): return left + (trial_ids.index(test) + 0.5) / len(trial_ids) * pw
    def ty(mg): return top + ph - (math.log10(mg) - ymin) / (ymax - ymin) * ph
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-weight:bold;font-size:20px}.label{font-size:14px}</style>',
             '<text x="90" y="34" class="title">PRESLHY E3.5: trial-level geometric mean bias</text>',
             f'<line x1="{left}" y1="{ty(1):.1f}" x2="{left+pw}" y2="{ty(1):.1f}" stroke="#777" stroke-dasharray="5 5"/>',
             f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#222"/>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#222"/>']
    for tick in (0.2, 0.5, 1, 2, 5, 10, 30):
        parts.append(f'<text x="45" y="{ty(tick)+5:.1f}" class="label">{tick:g}</text>')
    for number in trial_ids:
        parts.append(f'<text x="{tx(number)-10:.1f}" y="{top+ph+25}" class="label">{number}</text>')
    for row in selected:
        color, offset = (("#777777", -7) if row["variant"] == "historical_reconstruction"
                         else ("#0069aa", 7))
        x = tx(row["test_id"]) + offset
        lo, hi = ty(row["mg_resample_low_95"]), ty(row["mg_resample_high_95"])
        parts.append(f'<line x1="{x:.1f}" y1="{lo:.1f}" x2="{x:.1f}" y2="{hi:.1f}" stroke="{color}" stroke-width="2"/>')
        parts.append(f'<circle cx="{x:.1f}" cy="{ty(row["mg"]):.1f}" r="5" fill="{color}"/>')
    parts += [f'<text x="{left+pw/2-45}" y="{height-48}" class="label">Trial number</text>',
              '<text x="90" y="520" class="label">Log MG axis (observed/predicted); bars = descriptive within-trial arc resampling.</text>',
              '</svg>']
    path.with_name("e35_trial_bias.svg").write_text("\n".join(parts), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    if args.repeats < 1 or args.bootstrap < 100:
        raise ValueError("repeats >= 1 and bootstrap >= 100 required")
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"output directory not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261003)
    e35 = json.loads(E35.read_text(encoding="utf-8"))["trials"]
    ffi = {trial.test: trial for trial in spadeadam.load(FFI)}
    # Avoid charging one model with one-time process/table initialization.
    warmup = next(trial for trial in e35 if trial["trial"] == 10)
    for _, corrected in VARIANTS:
        run_e35(warmup, corrected, 1)
    arc_rows, times, trajectories = [], [], {}
    for trial in e35:
        if trial["orientation"] != "horizontal" or trial["flow_mean_gs"] <= 5:
            continue
        if nearfield._exit_ratio(trial, trial["flow_mean_gs"] / 1000) < nearfield.MOMENTUM_RATIO:
            continue
        number = int(trial["trial"])
        arcs = defaultdict(list)
        for sensor in trial["sensors"]:
            arcs[round(sensor["x"], 3)].append(sensor)
        by_variant = {}
        for name, corrected in VARIANTS:
            trajectory, samples = run_e35(trial, corrected, args.repeats)
            trajectories[("E35", number, name)] = trajectory
            by_variant[name] = trajectory
            times.append({"dataset": "E35", "test_id": number, "variant": name,
                          "repeats": args.repeats, "median_runtime_s": statlib.median(samples),
                          "min_runtime_s": min(samples), "max_runtime_s": max(samples),
                          "integration_limit_m": nearfield.REACH})
        for x, sensors in sorted(arcs.items()):
            predictions = {name: max(by_variant[name].concentration_at(x, s["y"], s["z"])
                                     for s in sensors) for name, _ in VARIANTS}
            if min(predictions.values()) <= 0:
                continue
            arc_rows.append({"dataset": "E35", "test_id": number, "x_m": x,
                             "n_sensors": len(sensors), "observed_vol_pct": max(s["peak"] for s in sensors),
                             **{f"{name}_vol_pct": predictions[name] for name, _ in VARIANTS}})
        print(f"E35 trial {number} complete", flush=True)
    for number in (4, 6):
        trial = ffi[number]
        by_variant = {}
        for name, corrected in VARIANTS:
            trajectory, samples = run_ffi(trial, corrected, args.repeats)
            trajectories[("FFI", number, name)] = trajectory
            by_variant[name] = trajectory
            times.append({"dataset": "FFI", "test_id": number, "variant": name,
                          "repeats": args.repeats, "median_runtime_s": statlib.median(samples),
                          "min_runtime_s": min(samples), "max_runtime_s": max(samples),
                          "integration_limit_m": 250.0})
        for x in trial.radii():
            observed = trial.arc(x)
            predictions = {name: max(by_variant[name].concentration_at(x, 0, z)
                                     for z in observed) for name, _ in VARIANTS}
            if min(predictions.values()) <= 0:
                continue
            arc_rows.append({"dataset": "FFI", "test_id": number, "x_m": x,
                             "n_sensors": len(observed), "observed_vol_pct": max(observed.values()),
                             **{f"{name}_vol_pct": predictions[name] for name, _ in VARIANTS}})
        print(f"FFI test {number} complete", flush=True)
    # Pair every model on exactly the same retained arcs.
    trial_rows = []
    for (dataset, test), group in sorted(defaultdict(list, {
        key: [r for r in arc_rows if (r["dataset"], r["test_id"]) == key]
        for key in {(r["dataset"], r["test_id"]) for r in arc_rows}
    }).items()):
        for name, _ in VARIANTS:
            result = metric(group, name)
            low, high = bootstrap_mg([r["observed_vol_pct"] / r[f"{name}_vol_pct"]
                                      for r in group], rng, args.bootstrap)
            time_row = next(t for t in times if t["dataset"] == dataset and
                            t["test_id"] == test and t["variant"] == name)
            trial_rows.append({"dataset": dataset, "test_id": test, "variant": name,
                               **result, "mg_resample_low_95": low,
                               "mg_resample_high_95": high,
                               "median_runtime_s": time_row["median_runtime_s"]})
    aggregate = []
    for dataset in ("E35", "FFI"):
        group = [r for r in arc_rows if r["dataset"] == dataset]
        for name, _ in VARIANTS:
            low_mg, high_mg, low_vg, high_vg, low_fac2, high_fac2 = cluster_interval(
                group, name, rng, args.bootstrap)
            aggregate.append({"dataset": dataset, "variant": name,
                              "n_tests": len({r["test_id"] for r in group}),
                              **metric(group, name), "mg_cluster_low_95": low_mg,
                              "mg_cluster_high_95": high_mg,
                              "vg_cluster_low_95": low_vg,
                              "vg_cluster_high_95": high_vg,
                              "fac2_cluster_low_95": low_fac2,
                              "fac2_cluster_high_95": high_fac2})
    # Test 6 is a real bunkering/transfer release, but its 30 m validation
    # mismatch makes the distance a *diagnostic*, not an operational setback.
    case = []
    for x in np.linspace(1.0, 250.0, 500):
        predictions = {name: trajectories[("FFI", 6, name)].concentration_at(
            float(x), 0.0, 1.5) for name, _ in VARIANTS}
        case.append({"x_m": float(x), **{f"{name}_vol_pct": predictions[name]
                                          for name, _ in VARIANTS}})
    reaches = {}
    for name, _ in VARIANTS:
        reaches[name] = dict(zip(("distance_m", "status"), lfl_reach(
            trajectories[("FFI", 6, name)], height=1.5)))
    write_csv(out / "paired_arcs.csv", arc_rows)
    write_csv(out / "trial_metrics.csv", trial_rows)
    write_csv(out / "aggregate_metrics.csv", aggregate)
    write_csv(out / "runtime_by_trial.csv", times)
    write_csv(out / "transfer_test6_lfl_profile.csv", case)
    _svg_plot(out / "ffi_paired_arc_errors.svg", arc_rows, case, trial_rows)
    manifest = {
        "comparator": "local corrections=False historical reconstruction; not original Fortran or external model",
        "matching": "identical source records, meteorology, model averaging=60 s, sensor coordinates and arc maxima",
        "temporal_limit": "observations are release-window peaks, not synchronized 60 s means; strict observation time matching unavailable",
        "uncertainty": "95% percentile bootstrap: arcs within test for trial rows; tests as clusters for aggregate rows; excludes source/meteorology/model-form uncertainty",
        "timing": "wall clock includes thermodynamics, source, integration, receptor extraction; trial 10 warmup of each variant excluded; warm single process; no concurrency",
        "lfl_case": {"source": "FFI Test 6 horizontal LH2 bunkering/transfer analogue",
                     "rate_kg_s": ffi[6].rate, "diameter_m": ffi[6].orifice,
                     "wind_low_m_s": ffi[6].wind_low,
                     "receptor": "centreline at 1.5 m height", "threshold_vol_pct": 4.0,
                     "distances": reaches,
                     "qualification": "conditional screening diagnostic, not an approved setback; Test 6 30 m model underprediction is unresolved"},
        "runtime_environment": {"python": platform.python_version(), "platform": platform.platform(),
                                "processor": platform.processor(), "numpy": np.__version__,
                                "timing_repeats": args.repeats, "bootstrap_draws": args.bootstrap},
        "input_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in
                         (E35, FFI / "conditions.csv", FFI / "sensors.csv",
                          ROOT / "src/degali/validation/nearfield.py", Path(__file__))},
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False,
                                                   allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out), "aggregate": aggregate, "lfl": reaches}, indent=2))


if __name__ == "__main__":
    main()
