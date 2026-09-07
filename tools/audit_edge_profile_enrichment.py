"""Preregistered fixed-transport square profile enrichment; not a plume solve."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from audit_transverse_mixing import ROOT, actual_source_model
from audit_reservoir_thermal_segments import serial
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.reservoir_thermal import ReservoirThermalMoments
from degali.addons.edge_enrichment import FixedTransportEdgeProjection


def scalar_maxima(residual):
    return dict(zip(("hydrogen", "heat"), [float(max(abs(a))) for a in np.split(residual, 2)]))


def assess(projection, fitted):
    parameters = fitted["parameters"]
    orders = (96, 192, 384)
    values = [projection.moments(parameters, order) for order in orders]
    errors = [float(max(abs(v-projection.target)/projection.scales)) for v in values]
    changes = [float(max(abs(values[i+1]-values[i])/projection.scales)) for i in (0, 1)]
    baseline = scalar_maxima(projection.edge_residual(np.zeros(projection.count), 4097))
    enriched = scalar_maxima(projection.edge_residual(parameters, 4097))
    d = projection.fields(parameters, projection.grid(384))
    edge = projection.fields(parameters, projection.edge(4097))
    ranges = {k: [float(min(np.min(d[k]), np.min(edge[k]))),
                  float(max(np.max(d[k]), np.max(edge[k])))]
              for k in ("rho", "y", "temperature", "log_c", "log_h")}
    physical = bool(ranges["rho"][0] > 0. and 0. < ranges["y"][0] <= ranges["y"][1] < 1.
                    and max(abs(x) for k in ("log_c", "log_h") for x in ranges[k]) <= .1000000001)
    edge_passed = all(enriched[k] <= .05 and enriched[k] < baseline[k] for k in baseline)
    moment_passed = errors[0] <= 1e-8 and errors[-1] <= 1e-8
    quadrature_passed = changes[-1] <= 1e-5
    return dict(**fitted, quadrature_orders=orders, actual_moments=values,
        moment_scaled_errors=errors, adjacent_quadrature_changes=changes,
        independent_face_samples=4097, baseline_edge_defects=baseline,
        enriched_edge_defects=enriched, sampled_ranges=ranges,
        edge_passed=edge_passed, moment_passed=moment_passed,
        quadrature_passed=quadrature_passed, physical_passed=physical,
        passed=bool(edge_passed and moment_passed and quadrature_passed and physical))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", nargs="+", type=int)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite final or partial enrichment evidence")
    ref = ROOT/"reference/preslhy"
    upstream_path = ref/"reservoir_thermal_adaptive_verification_2026-09-05.json"
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    if not upstream["combined_all_passed"]:
        raise RuntimeError("upstream weak reservoir verification has not passed")
    for relative, digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"upstream source/input changed: {relative}")
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    old_control = json.loads((ref/"phase_ambient_consistency_field_complete_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    selected = args.trials or [10, 11, 12, 22, 23, 24, 25]
    if not selected or len(set(selected)) != len(selected) or not set(selected).issubset(upstream["combined_passed_trials"]):
        parser.error("select unique frozen, passed trials")
    paths = [ROOT/p for p in upstream["sha256"]]
    paths += [upstream_path, Path(__file__).resolve(), ROOT/"src/degali/addons/edge_enrichment.py",
              ROOT/"tests/test_edge_enrichment.py", ROOT/"docs/prereg-edge-profile-enrichment.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes, rows, caches = hashes(), [], {}
    def payload(completed):
        return dict(phase="fixed_transport_square_scalar_enrichment", completed=completed,
            selected_trials=selected, rows=rows, sha256=initial_hashes,
            all_passed=len(rows)==len(selected) and all(any(f.get("passed", False) for f in r["fits"]) for r in rows),
            transport_reconstructed=False, downstream_solved=False, field_scored=False, promoted=False,
            limitations=["Boundary flux improvement is conditional on the frozen Gaussian transport field.",
                "No transport equations for the additional coefficients or observed-error reduction are established.",
                "Trust-region and physical checks are sampled, not formal global polynomial bounds."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    checkpoint()
    for n in selected:
        start = time.perf_counter()
        print(f"trial {n}: replay actual source and frozen boundary", flush=True)
        trial = next(t for t in reduced if t["trial"] == n)
        jp, th, source = actual_source_model(trial, measured_rows)
        old = old_control["interfaces"][str(n)]["downstream"]
        replay = IndependentEnergyCrosswind(jp, th).source_terms(np.array(old["states"][0]))
        source_error = float(max(abs(replay-old["sources"][0])/np.maximum(abs(np.array(old["sources"][0])), 1.)))
        entry = boundary["interfaces"][str(n)]
        section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=entry["thermal_width_ratio"],
                                                    quadrature_points=entry["quadrature_points"])
        section._quadrature_cache = caches.setdefault(section.k.delta, {})
        state = np.array(entry["state"])
        moments = section.moments(state)
        moment_error = float(max(abs(moments-entry["moments"])/section.moment_scales(entry["moments"])))
        if source_error > 1e-8 or moment_error > 1e-9:
            raise RuntimeError("actual-source/frozen-boundary replay failed")
        mixing = ConservativeTransverseMixing(section, state)
        reservoir = ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
        baseline = reservoir.evaluate(order=16, probes=1025)
        row = dict(trial=n, source=source, source_replay_error=source_error,
            boundary_replay_error=moment_error, baseline_weak_valid=bool(baseline["valid"]), fits=[])
        rows.append(row)
        checkpoint()
        for degree in (4, 6):
            begin = time.perf_counter()
            fit = dict(degree=degree, passed=False)
            row["fits"].append(fit)
            try:
                projection = FixedTransportEdgeProjection(reservoir, baseline, degree=degree)
                t = projection.target
                replay_six = np.array([t[0], t[1], t[2]*np.cos(state[3]), t[2]*np.sin(state[3]), t[3], t[4]])
                physical_replay = abs(replay_six-moments)/section.moment_scales(moments)
                # Force is separately scaled because near-neutral net force can be small.
                force_replay = abs(t[4]-moments[5])/max(abs(moments[5]), 1e-3)
                fit.update(target=t, scales=projection.scales, original_flux_force_replay=physical_replay,
                           original_force_replay=float(force_replay))
                if max(physical_replay) > 1e-5 or force_replay > 1e-5:
                    raise RuntimeError("phase-partition target differs from frozen boundary")
                def progress(h):
                    print(f"trial {n}, degree {degree}, order {h['order']}, iteration {h['iteration']}: "
                          f"edge max={h['edge_max']:.4g}, moments={h['moment_error']:.3e}", flush=True)
                fitted = projection.fit(callback=progress)
                fit.update(fitted)
                checkpoint()
                fit.update(assess(projection, fitted))
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                fit["failure"] = str(exc)
            fit["elapsed_seconds"] = time.perf_counter()-begin
            checkpoint()
            print(f"trial {n}, degree {degree}: passed={fit['passed']}, "
                  f"edge={fit.get('enriched_edge_defects')}, moment={fit.get('moment_scaled_errors')}, "
                  f"failure={fit.get('failure')}, elapsed={fit['elapsed_seconds']:.1f}s", flush=True)
        row["elapsed_seconds"] = time.perf_counter()-start
    if hashes() != initial_hashes:
        raise RuntimeError("code/input changed during the audit")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
