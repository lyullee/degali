"""Frozen-boundary weak reservoir model and preregistered millimetre segments."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from audit_transverse_mixing import ROOT, actual_source_model
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.reservoir_thermal import (
    ReservoirThermalMoments, ReservoirShortSegment, encode_section, decode_section,
)


def summarize(evaluation):
    out = evaluation.copy()
    for key in ("mixing", "family"):
        out.pop(key, None)
    out["linear_flux_residual"] = evaluation["family"].maximum_scaled_residual
    return out


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"cannot serialize {type(value)}")


def actual_values(jp, th, parameters):
    state, beta = decode_section(parameters)
    section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=beta)
    mixing = ConservativeTransverseMixing(section, state)
    return ReservoirShortSegment.flux_values(mixing)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", nargs="+", type=int)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite a final or partial reservoir audit")
    ref = ROOT/"reference/preslhy"
    upstream_path = ref/"shear_thermal_boundary_diagnostic_2026-09-05.json"
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    if not upstream["all_numerics_passed"]:
        raise RuntimeError("upstream diagnosis failed")
    for relative, digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"upstream source/input changed: {relative}")
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    old_control = json.loads((ref/"phase_ambient_consistency_field_complete_2026-09-05.json").read_text(encoding="utf-8"))
    old_screen = json.loads((ref/"shear_thermal_compatibility_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    selected = args.trials or [10, 11, 12, 22, 23, 24, 25]
    if not selected or len(set(selected)) != len(selected) or not set(selected).issubset(boundary["selected_trials"]):
        parser.error("select unique frozen trials")
    paths = [ROOT/p for p in upstream["sha256"]]
    paths += [upstream_path, Path(__file__).resolve(), ROOT/"src/degali/addons/reservoir_thermal.py",
              ROOT/"tests/test_reservoir_thermal.py", ROOT/"docs/prereg-reservoir-thermal-moments.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes, rows, caches = hashes(), [], {}
    def payload(completed):
        return dict(phase="weak_ambient_reservoir_thermal_moments", completed=completed, selected_trials=selected,
            rows=rows, all_boundary_numerics_passed=len(rows)==len(selected) and all(r["boundary_numerics_passed"] for r in rows),
            all_short_segments_passed=len(rows)==len(selected) and all(r.get("short_segment_passed", False) for r in rows),
            thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work", pointwise_boundary_closed=False,
            field_scored=False, promoted=False, sha256=initial_hashes,
            limitations=["Natural-flux weak two-moment closure; edge gradient defects are retained, not solved.",
                "Instantaneous shear heating and reduced buoyancy work are research assumptions, not a full TKE/pressure/gravity model.",
                "Millimetre short-segment convergence is not full downstream validation or evidence of lower measured error."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    checkpoint()
    for n in selected:
        start = time.perf_counter()
        print(f"trial {n}: replay frozen source/boundary and weak reservoir closure", flush=True)
        trial = next(t for t in reduced if t["trial"] == n)
        jp, th, source = actual_source_model(trial, measured_rows)
        old = old_control["interfaces"][str(n)]["downstream"]
        replay = IndependentEnergyCrosswind(jp, th).source_terms(np.array(old["states"][0]))
        source_error = float(np.max(abs(replay-old["sources"][0])/np.maximum(abs(np.array(old["sources"][0])), 1.)))
        entry = boundary["interfaces"][str(n)]
        state, beta = np.array(entry["state"]), entry["thermal_width_ratio"]
        section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=beta, quadrature_points=entry["quadrature_points"])
        section._quadrature_cache = caches.setdefault(section.k.delta, {})
        moments = section.moments(state)
        moment_error = float(np.max(abs(moments-entry["moments"])/section.moment_scales(entry["moments"])))
        if source_error > 1e-8 or moment_error > 1e-9:
            raise RuntimeError("actual-source/frozen-boundary replay failed")
        mixing = ConservativeTransverseMixing(section, state)
        model = ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
        coarse, fine = model.evaluate(order=8, probes=1025), model.evaluate(order=16, probes=1025)
        old_sources = next(r for r in old_screen["rows"] if r["trial"] == n)["sources"]
        source_quad_error = float(np.max(abs(fine["ledger"]["original"]-old_sources)/np.maximum(abs(np.array(old_sources)), 1.)))
        refinement = max(float(np.max(abs(coarse["budgets"][k]-v)/np.maximum(abs(v), 1.))) for k, v in fine["budgets"].items())
        root_change = abs(coarse["gamma"]-fine["gamma"])/max(abs(fine["gamma"]), 1.)
        parameters = encode_section(state, beta)
        values = ReservoirShortSegment.flux_values(mixing)
        scales = np.maximum(abs(values), [1e-12, 1e-12, 1., 1., 1., 1.])
        steps = [ds/max(np.max(abs(fine["rates"])), 1.) for ds in (2e-4, 1e-4)]
        derivative = [(actual_values(jp, th, parameters+ds*fine["rates"])-actual_values(jp, th, parameters-ds*fine["rates"]))/ (2*ds) for ds in steps]
        target = np.r_[fine["ledger"]["sources"], fine["moment_rate"]]
        fd_error = float(np.max(abs(derivative[1]-target)/scales))
        fd_step = float(np.max(abs(derivative[1]-derivative[0])/scales))
        numerical = bool(fine["family"].maximum_scaled_residual <= 1e-8 and max(source_quad_error, refinement, root_change, fd_error, fd_step) <= 1e-5)
        row = dict(trial=n, source=source, source_replay_error=source_error, boundary_moment_replay_error=moment_error,
            source_quadrature_replay_error=source_quad_error, initial_parameters=parameters, initial_fluxes=values,
            boundary=summarize(fine), budget_8_16_scaled_difference=refinement, root_8_16_change=root_change,
            independent_flux_change=dict(steps_m=steps, derivatives=derivative, expected=target, scaled_error=fd_error, scaled_step_change=fd_step),
            boundary_numerics_passed=numerical, segments=[])
        rows.append(row)
        checkpoint()
        if numerical and fine["valid"]:
            driver = ReservoirShortSegment(jp, th, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
            length = min(.01, .05*min(mixing.sy, mixing.sn))
            row["segment_length_m"] = length
            for divisions in (2, 4, 8):
                try:
                    segment = driver.integrate(parameters, length, divisions)
                    endpoint = driver.evaluate(segment["parameters"], order=16, probes=1025)
                    end_coarse = driver.evaluate(segment["parameters"], order=8, probes=1025)
                    end_values = driver.flux_values(endpoint["mixing"])
                    balance = (end_values-values-segment["cumulative_sources"])/scales
                    end_refinement = max(float(np.max(abs(end_coarse["budgets"][k]-v)/np.maximum(abs(v), 1.))) for k, v in endpoint["budgets"].items())
                    segment.update(divisions=divisions, end_fluxes=end_values, scaled_balance=balance,
                        maximum_balance=float(max(abs(balance))), endpoint=summarize(endpoint), endpoint_8_16_difference=end_refinement,
                        passed=bool(max(abs(balance)) <= 1e-5 and end_refinement <= 1e-5))
                    row["segments"].append(segment)
                    checkpoint()
                    print(f"trial {n}: {divisions} steps over {1000*length:.3f} mm, balance={segment['maximum_balance']:.3e}", flush=True)
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                    row["segment_failure"] = str(exc)
                    checkpoint()
                    break
            if len(row["segments"]) == 3:
                ends = [s["parameters"] for s in row["segments"]]
                changes = [float(np.max(abs(ends[i+1]-ends[i])/np.maximum(abs(ends[i+1]), 1.))) for i in (0, 1)]
                row["endpoint_2_4_and_4_8_changes"] = changes
                row["short_segment_passed"] = bool(changes[1] <= 1e-5 and row["segments"][-1]["passed"])
                end_state, end_beta = decode_section(ends[-1])
                row["end_state"], row["end_beta"], row["beta_fractional_change"] = end_state, end_beta, end_beta/beta-1.
        row["elapsed_seconds"] = time.perf_counter()-start
        checkpoint()
        print(f"trial {n}: boundary numeric={numerical}, weak valid={fine['valid']}, short passed={row.get('short_segment_passed', False)}, "
              f"pointwise heat defect={fine['edge_gradient_defects']['heat']:.3f}, elapsed={row['elapsed_seconds']:.1f}s", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("code/input changed during the audit")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
