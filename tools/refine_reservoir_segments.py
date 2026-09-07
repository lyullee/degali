"""Preserve failed 8-step records; separately check 16/32/64-step numerics."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from audit_reservoir_thermal_segments import ROOT, actual_source_model, serial, summarize
from degali.addons.reservoir_thermal import ReservoirShortSegment, decode_section


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("screen", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite refinement")
    screen = json.loads(args.screen.read_text(encoding="utf-8"))
    if not screen["completed"] or not screen["all_boundary_numerics_passed"]:
        raise RuntimeError("requires completed boundary checks")
    for relative, digest in screen["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"screen input/code changed: {relative}")
    paths = [ROOT/p for p in screen["sha256"]]
    paths += [args.screen.resolve(), Path(__file__).resolve(), ROOT/"docs/reservoir-segment-refinement-note.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    ref = ROOT/"reference/preslhy"
    trials = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    rows = []
    for original in screen["rows"]:
        if original.get("short_segment_passed", False):
            continue
        if "segment_length_m" not in original or "segment_failure" in original:
            raise RuntimeError("this is not a simple finite-step convergence failure")
        n, start = original["trial"], time.perf_counter()
        jp, th, _ = actual_source_model(next(t for t in trials if t["trial"] == n), measured_rows)
        driver = ReservoirShortSegment(jp, th, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
        parameters = np.array(original["initial_parameters"])
        values = np.array(original["initial_fluxes"])
        scales = np.maximum(abs(values), [1e-12, 1e-12, 1., 1., 1., 1.])
        results = []
        for divisions in (16, 32, 64):
            part = driver.integrate(parameters, original["segment_length_m"], divisions)
            end = driver.evaluate(part["parameters"], order=16, probes=1025)
            end_values = driver.flux_values(end["mixing"])
            balance = (end_values-values-part["cumulative_sources"])/scales
            part.update(divisions=divisions, endpoint=summarize(end), end_fluxes=end_values,
                        scaled_balance=balance, maximum_balance=float(max(abs(balance))))
            results.append(part)
            print(f"trial {n}: {divisions} steps, balance={part['maximum_balance']:.3e}", flush=True)
        change = float(np.max(abs(results[-1]["parameters"]-results[-2]["parameters"])/np.maximum(abs(results[-1]["parameters"]), 1.)))
        _, beta = decode_section(results[-1]["parameters"])
        row = dict(trial=n, original_short_segment_passed=False, segment_length_m=original["segment_length_m"],
            segments=results, endpoint_32_64_change=change, end_beta=beta,
            beta_fractional_change=beta/np.exp(parameters[7])-1.,
            passed=bool(results[-1]["maximum_balance"] <= 1e-5 and change <= 1e-5), elapsed_seconds=time.perf_counter()-start)
        rows.append(row)
    if hashes() != initial:
        raise RuntimeError("code/input changed during refinement")
    passed = {r["trial"] for r in screen["rows"] if r.get("short_segment_passed", False)}
    passed.update(r["trial"] for r in rows if r["passed"])
    payload = dict(phase="separate_reservoir_segment_refinement", rows=rows,
        original_all_short_segments_passed=screen["all_short_segments_passed"],
        combined_passed_trials=sorted(passed), combined_all_passed=set(screen["selected_trials"]) == passed,
        pointwise_boundary_closed=False, field_scored=False, promoted=False, sha256=initial)
    args.output.write_text(json.dumps(payload, default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
