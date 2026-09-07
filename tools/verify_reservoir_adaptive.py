"""Independent adaptive check after retained nonmonotonic midpoint failures."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp

from audit_reservoir_thermal_segments import ROOT, actual_source_model, serial, summarize
from degali.addons.reservoir_thermal import ReservoirShortSegment, decode_section


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("refinement", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite adaptive verification")
    refinement = json.loads(args.refinement.read_text(encoding="utf-8"))
    for relative, digest in refinement["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"refinement input/code changed: {relative}")
    screen_path = ROOT/"reference/preslhy/reservoir_thermal_segments_2026-09-05.json"
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    paths = [ROOT/p for p in refinement["sha256"]]
    paths += [args.refinement.resolve(), Path(__file__).resolve(), ROOT/"docs/reservoir-adaptive-verification-note.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    ref = ROOT/"reference/preslhy"
    trials = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    rows = []
    for previous in refinement["rows"]:
        if previous["passed"]:
            continue
        n, start = previous["trial"], time.perf_counter()
        original = next(r for r in screen["rows"] if r["trial"] == n)
        jp, th, _ = actual_source_model(next(t for t in trials if t["trial"] == n), measured_rows)
        driver = ReservoirShortSegment(jp, th, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
        parameters, values = np.array(original["initial_parameters"]), np.array(original["initial_fluxes"])
        scales = np.maximum(abs(values), [1e-12, 1e-12, 1., 1., 1., 1.])
        length, results = original["segment_length_m"], []
        for tolerance in (1e-8, 1e-10):
            calls, worst_error, largest_defect = 0, 0., 0.
            def rhs(s, state):
                nonlocal calls, worst_error, largest_defect
                calls += 1
                if calls > 2000:
                    raise RuntimeError("adaptive verification exceeded its 2000-call numerical budget")
                out = driver.evaluate(state[:8])
                worst_error = max(worst_error, out["weak_budget_scaled_error"])
                largest_defect = max(largest_defect, out["edge_gradient_defects"]["heat"])
                if calls % 128 == 0:
                    print(f"trial {n}, rtol={tolerance}: {calls} evaluations, s={1000*s:.3f} mm", flush=True)
                return np.r_[out["rates"], out["ledger"]["sources"], out["moment_rate"]]
            atol = tolerance*.01*np.r_[np.maximum(abs(parameters), 1.), scales]
            result = solve_ivp(rhs, (0., length), np.r_[parameters, np.zeros(6)], method="DOP853",
                               rtol=tolerance, atol=atol, max_step=length/16)
            if not result.success:
                raise RuntimeError(result.message)
            end_parameters, cumulative = result.y[:8, -1], result.y[8:, -1]
            end = driver.evaluate(end_parameters, order=16, probes=1025)
            end_values = driver.flux_values(end["mixing"])
            balance = (end_values-values-cumulative)/scales
            results.append(dict(rtol=tolerance, atol=atol, rhs_calls=calls, steps=len(result.t)-1,
                parameters=end_parameters, cumulative_sources=cumulative, end_fluxes=end_values, endpoint=summarize(end),
                maximum_weak_residual=worst_error, maximum_edge_heat_defect=largest_defect,
                scaled_balance=balance, maximum_balance=float(max(abs(balance)))))
            print(f"trial {n}, rtol={tolerance}: complete, {calls} evaluations, balance={max(abs(balance)):.3e}", flush=True)
        change = float(np.max(abs(results[1]["parameters"]-results[0]["parameters"])/np.maximum(abs(results[1]["parameters"]), 1.)))
        _, beta = decode_section(results[-1]["parameters"])
        rows.append(dict(trial=n, results=results, endpoint_tolerance_change=change,
            end_beta=beta, beta_fractional_change=beta/np.exp(parameters[7])-1.,
            passed=bool(results[-1]["maximum_balance"] <= 1e-5 and change <= 1e-5), elapsed_seconds=time.perf_counter()-start))
    if hashes() != initial:
        raise RuntimeError("input or code changed during verification")
    passed = set(refinement["combined_passed_trials"])
    passed.update(r["trial"] for r in rows if r["passed"])
    payload = dict(phase="independent_adaptive_reservoir_segment_verification", rows=rows,
        original_all_short_segments_passed=screen["all_short_segments_passed"],
        midpoint_refinement_all_passed=refinement["combined_all_passed"], combined_passed_trials=sorted(passed),
        combined_all_passed=set(screen["selected_trials"]) == passed,
        pointwise_boundary_closed=False, field_scored=False, promoted=False, sha256=initial)
    args.output.write_text(json.dumps(payload, default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
