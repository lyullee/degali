"""Frozen seven-case conservative shape refit with face/ray phase splitting."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from audit_transverse_mixing import ROOT, actual_source_model
from audit_reservoir_thermal_segments import serial
from audit_edge_profile_enrichment import scalar_maxima
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.reservoir_thermal import ReservoirThermalMoments
from degali.addons.edge_enrichment import FixedTransportEdgeProjection
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit


def construct_projection(old, boundary, reduced, measured, caches):
    n = old["trial"]
    trial = next(t for t in reduced if t["trial"] == n)
    jp, th, source = actual_source_model(trial, measured)
    entry = boundary["interfaces"][str(n)]
    section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=entry["thermal_width_ratio"],
                                               quadrature_points=entry["quadrature_points"])
    section._quadrature_cache = caches.setdefault(section.k.delta, {})
    mixing = ConservativeTransverseMixing(section, np.array(entry["state"]))
    reservoir = ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
    baseline = reservoir.evaluate(order=16, probes=1025)
    projection = FixedTransportEdgeProjection(reservoir, baseline, degree=old["degree"])
    if max(abs(projection.target-old["target"])/projection.scales) > 1e-12:
        raise RuntimeError("frozen moment target failed to replay")
    return projection, source


def verify(projection, parameters):
    quad = FaceSplitSquareMoments(projection)
    v8 = quad.moments(parameters)
    v16 = quad.moments(parameters, order=16, angular_order=16, probes=2049)
    errors = [float(max(abs(v-projection.target)/projection.scales)) for v in (v8, v16)]
    change = float(max(abs(v16-v8)/projection.scales))
    old_edge = scalar_maxima(projection.edge_residual(np.zeros(projection.count), 4097))
    new_edge = scalar_maxima(projection.edge_residual(parameters, 4097))
    t = np.linspace(0., 1., 129)
    a, b = np.meshgrid(t, t, indexing="ij")
    d = projection.fields(parameters, projection.prepare(a.ravel(), b.ravel()))
    edge = projection.fields(parameters, projection.edge(4097))
    ranges = {k: [float(min(np.min(d[k]), np.min(edge[k]))),
                  float(max(np.max(d[k]), np.max(edge[k])))]
              for k in ("rho", "y", "temperature", "log_c", "log_h")}
    physical = bool(ranges["rho"][0] > 0. and 0. < ranges["y"][0] <= ranges["y"][1] < 1.
                    and max(abs(x) for k in ("log_c", "log_h") for x in ranges[k]) <= .1000000001)
    edge_passed = bool(all(new_edge[k] <= .05 and new_edge[k] < old_edge[k] for k in new_edge))
    moments_passed = bool(max(errors) <= 1e-8 and change <= 1e-8)
    return dict(actual_moments=[v8, v16], moment_scaled_errors=errors, adjacent_difference=change,
        baseline_edge_defects=old_edge, enriched_edge_defects=new_edge, sampled_ranges=ranges,
        edge_passed=edge_passed, moment_passed=moments_passed, physical_passed=physical,
        passed=bool(edge_passed and moments_passed and physical))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", nargs="+", type=int)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite a final or partial refit audit")
    ref = ROOT/"reference/preslhy"
    upstream_path = ref/"edge_enrichment_phase_diagnosis_2026-09-05.json"
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    if not upstream["completed"]:
        raise RuntimeError("upstream diagnosis is incomplete")
    for relative, digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"upstream source/input changed: {relative}")
    selected = args.trials or [10, 11, 12, 22, 23, 24, 25]
    if len(set(selected)) != len(selected) or not set(selected).issubset({r["trial"] for r in upstream["rows"]}):
        parser.error("select unique frozen trials")
    paths = [ROOT/p for p in upstream["sha256"]]
    paths += [upstream_path, Path(__file__).resolve(), ROOT/"src/degali/addons/edge_conservative_refit.py",
              ROOT/"tests/test_edge_conservative_refit.py", ROOT/"docs/prereg-edge-conservative-refit.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes, rows, caches = hashes(), [], {}
    def payload(completed):
        return dict(phase="face_ray_phase_split_conservative_edge_refit", completed=completed,
            selected_trials=selected, rows=rows, sha256=initial_hashes,
            all_passed=len(rows)==len(selected) and all(r.get("passed", False) for r in rows),
            transport_reconstructed=False, shape_ode_solved=False, field_scored=False, promoted=False,
            limitations=["A fixed-transport cross-section projection, not a new downstream solution.",
                "Prior selected degrees are fixed; no observation fit or original failure reclassification.",
                "Face-event discovery, monotonicity and trust bounds are sampled, not globally proven."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))["trials"]
    measured = {r["trial"]: r for r in measured}
    checkpoint()
    for n in selected:
        begin = time.perf_counter()
        old = next(r for r in upstream["rows"] if r["trial"] == n)
        row = dict(trial=n, degree=old["degree"], initial_parameters=old["parameters"], passed=False)
        rows.append(row)
        print(f"trial {n}, degree {old['degree']}: independent face/ray splitting", flush=True)
        try:
            p, row["source"] = construct_projection(old, boundary, reduced, measured, caches)
            quad = FaceSplitSquareMoments(p)
            z = np.zeros(p.count)
            v8, v16 = quad.moments(z), quad.moments(z, order=16, angular_order=16, probes=2049)
            zero_errors = [float(max(abs(v-p.target)/p.scales)) for v in (v8, v16)]
            row.update(target=p.target, scales=p.scales, zero_scaled_errors=zero_errors,
                       zero_adjacent_difference=float(max(abs(v8-v16)/p.scales)))
            if max(zero_errors+[row["zero_adjacent_difference"]]) > 1e-8:
                raise RuntimeError("face/ray quadrature zero-reference check failed")
            row["quadrature_reference_passed"] = True
            checkpoint()
            print(f"trial {n}: zero reference errors={zero_errors}", flush=True)
            def progress(h):
                row["latest_progress"] = h
                checkpoint()
                print(f"trial {n}: {h}", flush=True)
            fitted = ConservativeEdgeRefit(p).fit(np.array(old["parameters"]), callback=progress)
            row.update(fitted)
            checkpoint()
            row.update(verify(p, fitted["parameters"]))
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
        row["elapsed_seconds"] = time.perf_counter()-begin
        checkpoint()
        print(f"trial {n}: passed={row['passed']}, moments={row.get('moment_scaled_errors')}, "
              f"edge={row.get('enriched_edge_defects')}, failure={row.get('failure')}, "
              f"elapsed={row['elapsed_seconds']:.1f}s", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("code/input changed during refit audit")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
