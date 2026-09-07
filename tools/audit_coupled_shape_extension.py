"""Unchanged pilot algorithm applied to the five remaining preregistered cases."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK, compact, independent_rays
from degali.addons.enriched_transport import PrescribedRadialMixing, EnrichedModalTransport
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit
from degali.addons.coupled_shape_initialization import FixedMeshCoupledTransport, SimultaneousShapeInitializer


def verify_candidate(p, par, supplied, row, checkpoint):
    row["stage"] = "independent_moving_phase_retraction"
    checkpoint()
    par, row["independent_retraction"] = ConservativeEdgeRefit(p).retract(par)
    row["parameters"] = par.copy()
    quad = FaceSplitSquareMoments(p)
    moments = [quad.moments(par), quad.moments(par, order=16, angular_order=16, probes=2049)]
    row["moment_errors"] = [float(max(abs(v-p.target)/p.scales)) for v in moments]
    row["moment_refinement"] = float(max(abs(moments[1]-moments[0])/p.scales))
    row["evaluations"] = []
    for order, angular in ((4, 8), (8, 16)):
        row["stage"] = f"independent_transport_{order}_{angular}"
        checkpoint()
        model = EnrichedModalTransport(p, par, scalar_mixing=supplied, thermal_species_ratio=1., mechanical_work=WORK)
        out = model.assemble(order=order, angular_order=angular)
        row["evaluations"].append(out)
        checkpoint()
        print(f"trial {row['trial']}: independent {order}/{angular}, edge={out['edge_defects']}, chi_P={out['minimum_chi_momentum']:.4g}", flush=True)
    coarse, fine = row["evaluations"]
    row["rate_refinement"] = float(max(abs(fine["rates"]-coarse["rates"])/np.maximum(abs(fine["rates"]), 1.)))
    row["matrix_refinement"] = float(np.max(abs(fine["matrix"]-coarse["matrix"])/np.maximum(abs(fine["matrix"]), 1.)))
    row["stage"] = "independent_65_rays_and_fields"
    checkpoint()
    row["independent_rays"] = independent_rays(model, fine["rates"])
    t = np.linspace(0., 1., 129)
    a, b = np.meshgrid(t, t, indexing="ij")
    fields = p.fields(par, p.prepare(a.ravel(), b.ravel()))
    ranges = {k: [float(min(fields[k])), float(max(fields[k]))] for k in ("rho", "y", "temperature", "log_c", "log_h")}
    row["ranges"] = ranges
    row["numerics_passed"] = bool(max(row["moment_errors"]+[row["moment_refinement"]]) <= 1e-8
        and max(row["rate_refinement"], row["matrix_refinement"]) <= 1e-5
        and max(fine["linear_scaled_error"], fine["weak_heat_scaled_error"]) <= 1e-8
        and abs(fine["mass_boundary_error"])/max(abs(fine["source"][0]), 1.) <= 1e-8)
    rays = row["independent_rays"]
    row["constitutive_passed"] = bool(max(list(fine["edge_defects"].values())+list(rays["edge_defects"].values())) <= .05
        and min(fine["minimum_chi_momentum"], rays["minimum_chi_momentum"]) >= 0.
        and max(fine["maximum_outward_mass"], rays["maximum_outward_mass"]) < 0.
        and fine["curvature_half_width"] < .1 and ranges["rho"][0] > 0.
        and 0. < ranges["y"][0] <= ranges["y"][1] < 1.
        and max(abs(v) for k in ("log_c", "log_h") for v in ranges[k]) <= .1000000001)
    row["passed"] = row["numerics_passed"] and row["constitutive_passed"]
    return par


def run_case(n, frozen, old_transport, boundary, reduced, measured, caches, row, checkpoint):
    old = next(r for r in frozen["rows"] if r["trial"] == n)
    row["degree"] = old["degree"]
    row["previous_defects"] = next(r for r in old_transport["rows"] if r["trial"] == n)["evaluations"][-1]["edge_defects"]
    row["stage"] = "source_reconstruction"
    checkpoint()
    print(f"trial {n}: reconstruct actual source, degree{old['degree']}", flush=True)
    p, row["source_reconstruction"] = construct_projection(old, boundary, reduced, measured, caches)
    supplied = PrescribedRadialMixing.from_weak_baseline(p, order=32)
    par = np.array(old["parameters"])
    row["initial_parameters"] = par.copy()
    for order, angular, iterations in ((4, 24, 6), (8, 48, 3)):
        row["stage"] = f"fitting_{order}_{angular}"
        stage = dict(radial_order=order, angular_order=angular, maximum_iterations=iterations, history=[])
        row["stages"].append(stage)
        checkpoint()
        mesh = FixedMeshCoupledTransport(p, par, scalar_mixing=supplied, thermal_species_ratio=1.,
            mechanical_work=WORK, order=order, angular_order=angular)
        mark = time.perf_counter()
        stage["before"] = compact(mesh.evaluate(par))
        stage["one_evaluation_seconds"] = time.perf_counter()-mark
        print(f"trial {n}: fit mesh {order}/{angular}, evaluation {stage['one_evaluation_seconds']:.2f}s", flush=True)
        def progress(item):
            stage["history"].append(item)
            checkpoint()
            print(f"trial {n}, mesh {order}/{angular}, iter {item['iteration']}: "
                f"edge={item['edge_defects']}, chi_P={item['minimum_chi_momentum']:.4g}, norm={item['objective_norm']:.4g}", flush=True)
        fitted = SimultaneousShapeInitializer(mesh).fit(par, maximum_iterations=iterations, callback=progress)
        par = fitted["parameters"]
        stage.update(history=fitted["history"], after=compact(fitted["output"]), parameters=par.copy())
        checkpoint()
    row["fit_parameters"] = par.copy()
    verify_candidate(p, par, supplied, row, checkpoint)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, nargs="+", default=[10, 12, 22, 23, 24])
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite extension evidence")
    if len(set(args.trials)) != len(args.trials) or not set(args.trials).issubset({10, 12, 22, 23, 24}):
        parser.error("select unique remaining trials10/12/22/23/24")
    ref = ROOT/"reference/preslhy"
    previous = ref/"coupled_shape_initialization_pilots_2026-09-05.json"
    pilot = json.loads(previous.read_text(encoding="utf-8"))
    if not pilot["completed"] or set(pilot["passed_trials"]) != {25, 11}:
        raise ValueError("requires completed two-pilot baseline")
    for name, digest in pilot["sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"frozen pilot changed: {name}")
    paths = [ROOT/n for n in pilot["sha256"]]+[previous, Path(__file__).resolve(),
        ROOT/"tests/test_coupled_shape_degree6.py", ROOT/"docs/prereg-coupled-shape-extension.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding="utf-8"))
    frozen = read("edge_conservative_refit_combined_2026-09-05.json")
    transport = read("enriched_transport_combined_2026-09-05.json")
    boundary = read("buoyancy_constrained_enthalpy_width_interface_2026-09-05.json")
    reduced = read("e35_reduced.json")["trials"]
    measured = {r["trial"]: r for r in read("measured_pipe_source_2026-09-05.json")["trials"]}
    rows, caches = [], {}
    began = datetime.now(timezone.utc).isoformat()
    def payload(completed):
        return dict(phase="simultaneous_shape_rate_initialization_extension", completed=completed,
            process_id=os.getpid(), started_utc=began, updated_utc=datetime.now(timezone.utc).isoformat(),
            selected_trials=args.trials, sha256=initial, rows=rows,
            passed_trials=[r["trial"] for r in rows if r.get("passed", False)],
            downstream_integrated=False, field_scored=False, promoted=False,
            limitations=["Same pilot method/gates; extension cases are not observations or downstream runs.",
                "Sampled positivity and fixed chi_C formula, not a full turbulence closure."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    checkpoint()
    for n in args.trials:
        started = time.perf_counter()
        row = dict(trial=n, passed=False, stages=[])
        rows.append(row)
        try:
            run_case(n, frozen, transport, boundary, reduced, measured, caches, row, checkpoint)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
        row["stage"] = "completed"
        row["elapsed_seconds"] = time.perf_counter()-started
        checkpoint()
        print(f"trial {n}: passed={row['passed']}, failure={row.get('failure')}, seconds={row['elapsed_seconds']:.1f}", flush=True)
    if hashes() != initial:
        raise RuntimeError("source/input changed during extension")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
