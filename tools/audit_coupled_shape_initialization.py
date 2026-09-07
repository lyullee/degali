"""Two preregistered pilots; fit meshes are never the acceptance authority."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from degali.addons.enriched_transport import PrescribedRadialMixing, EnrichedModalTransport
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit
from degali.addons.coupled_shape_initialization import FixedMeshCoupledTransport, SimultaneousShapeInitializer


WORK = "reduced_buoyancy_work_immediate_shear_heat"


def compact(out):
    return {k: out[k] for k in ("rates", "edge_defects", "minimum_chi_momentum", "moments",
        "maximum_outward_mass", "matrix_condition", "linear_scaled_error", "weak_heat_scaled_error",
        "mass_boundary_error", "curvature_half_width")}


def independent_rays(model, rates):
    p, m = model.projection, model.projection.mixing
    angles = np.linspace(0., math.pi/4, 65)
    at = np.r_[1., rates]
    minimum_chi, maximum_mass = math.inf, -math.inf
    defects = np.zeros(3)
    for start in range(0, len(angles), 8):
        subset = angles[start:start+8]
        for angle, k in zip(subset, model.quadrature.partitions(model.parameters, subset)):
            extra = model.mixing.knots[(model.mixing.knots > 0.) & (model.mixing.knots < k[-1])]
            k = np.sort(np.r_[k, extra])
            k = k[np.r_[True, np.diff(k) > k[-1]*2e-13]]
            d, fm, fp, fe, _ = model._ray(angle, k, 16)
            chi = -((fp-d["u"][:, None]*fm)@at)/(model.area*d["rho"]*d["uq"])
            minimum_chi = min(minimum_chi, float(min(chi)))
            t = math.tan(angle)
            ed = model.local(np.array([1.]), np.array([t]))
            em, ep = fe@at
            ue, rho, y, h = [ed[s][0] for s in ("u", "rho", "y", "specific_h")]
            qe = np.array([p.q0*(1+t*t)])
            ce, oldfm = model.mixing.evaluate(qe)[0]
            old = m.local(qe)
            fh = .5*m.wind**2*em-(ue*ep-.5*ue*ue*em)
            ec = y*em-model.area*rho*ce*ed["grad_y"][0, 0]/(2*p.q0)
            eh = h*em-model.area*rho*ce*ed["grad_h"][0, 0]/(2*p.q0)-fh
            epx = ep-m.wind*math.cos(m.state[3])*em
            scales = np.array([max(abs(old["y"][0]*oldfm), 1e-12),
                max(abs(old["h"][0]/old["rho"][0]*oldfm), 1.), max(abs(old["u"][0]*oldfm), 1.)])
            defects = np.maximum(defects, abs(np.array([ec, eh, epx]))/scales)
            minimum_chi = min(minimum_chi, float(-(ep-ue*em)/(model.area*rho*ed["uq"][0])))
            maximum_mass = max(maximum_mass, float(em))
    return dict(angles=65, radial_order=16, edge_defects=dict(zip(("hydrogen", "heat", "momentum"), defects)),
        minimum_chi_momentum=minimum_chi, maximum_outward_mass=maximum_mass)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, nargs="+", default=[25, 11])
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite initialization evidence")
    if len(set(args.trials)) != len(args.trials) or not set(args.trials).issubset({25, 11}):
        parser.error("only the two preregistered pilots 25 and 11")
    ref = ROOT/"reference/preslhy"
    previous = ref/"enriched_transport_combined_2026-09-05.json"
    upstream = json.loads(previous.read_text(encoding="utf-8"))
    for name, digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"frozen upstream changed: {name}")
    paths = [ROOT/n for n in upstream["sha256"]]+[previous, Path(__file__).resolve(),
        ROOT/"src/degali/addons/coupled_shape_initialization.py",
        ROOT/"tests/test_coupled_shape_initialization.py", ROOT/"docs/prereg-coupled-shape-initialization.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes = hashes()
    frozen = json.loads((ref/"edge_conservative_refit_combined_2026-09-05.json").read_text(encoding="utf-8"))
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))["trials"]
    measured = {r["trial"]: r for r in measured}
    rows, caches = [], {}
    def payload(completed):
        return dict(phase="simultaneous_shape_rate_initialization_pilots", completed=completed,
            selected_trials=args.trials, sha256=initial_hashes, rows=rows,
            passed_trials=[r["trial"] for r in rows if r.get("passed", False)],
            downstream_integrated=False, field_scored=False, promoted=False,
            limitations=["Fixed-mesh iterations require separate moving-phase acceptance checks.",
                "Only two preregistered pilots, not seven-case or observed validation.",
                "Prescribed chi_C, radial conservative flow and immediate shear heating remain conditional."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    checkpoint()
    for n in args.trials:
        started = time.perf_counter()
        row = dict(trial=n, passed=False, stages=[])
        rows.append(row)
        checkpoint()
        try:
            old = next(r for r in frozen["rows"] if r["trial"] == n)
            row["previous_defects"] = next(r for r in upstream["rows"] if r["trial"] == n)["evaluations"][-1]["edge_defects"]
            print(f"trial {n}: reconstruct actual source and initial shape", flush=True)
            p, row["source_reconstruction"] = construct_projection(old, boundary, reduced, measured, caches)
            supplied = PrescribedRadialMixing.from_weak_baseline(p, order=32)
            par = np.array(old["parameters"])
            row["initial_parameters"] = par.copy()
            for order, angular, iterations in ((4, 24, 6), (8, 48, 3)):
                stage = dict(radial_order=order, angular_order=angular, maximum_iterations=iterations, history=[])
                row["stages"].append(stage)
                mesh = FixedMeshCoupledTransport(p, par, scalar_mixing=supplied, thermal_species_ratio=1.,
                    mechanical_work=WORK, order=order, angular_order=angular)
                mark = time.perf_counter()
                stage["before"] = compact(mesh.evaluate(par))
                stage["one_evaluation_seconds"] = time.perf_counter()-mark
                print(f"trial {n}: fit mesh {order}/{angular}, one evaluation {stage['one_evaluation_seconds']:.2f}s", flush=True)
                def progress(item):
                    stage["history"].append(item)
                    checkpoint()
                    print(f"trial {n}, mesh {order}/{angular}, iter {item['iteration']}: "
                        f"edge={item['edge_defects']}, chi_P={item['minimum_chi_momentum']:.4g}, "
                        f"norm={item['objective_norm']:.4g}", flush=True)
                fitted = SimultaneousShapeInitializer(mesh).fit(par, maximum_iterations=iterations, callback=progress)
                par = fitted["parameters"]
                stage.update(history=fitted["history"], after=compact(fitted["output"]), parameters=par.copy())
                checkpoint()
            row["fit_parameters"] = par.copy()
            print(f"trial {n}: independent moving-phase conservation retraction", flush=True)
            par, row["independent_retraction"] = ConservativeEdgeRefit(p).retract(par)
            row["parameters"] = par.copy()
            quad = FaceSplitSquareMoments(p)
            moments = [quad.moments(par), quad.moments(par, order=16, angular_order=16, probes=2049)]
            row["moment_errors"] = [float(max(abs(v-p.target)/p.scales)) for v in moments]
            row["moment_refinement"] = float(max(abs(moments[1]-moments[0])/p.scales))
            row["evaluations"] = []
            for order, angular in ((4, 8), (8, 16)):
                model = EnrichedModalTransport(p, par, scalar_mixing=supplied, thermal_species_ratio=1., mechanical_work=WORK)
                out = model.assemble(order=order, angular_order=angular)
                row["evaluations"].append(out)
                checkpoint()
                print(f"trial {n}: independent {order}/{angular}, edge={out['edge_defects']}, chi_P={out['minimum_chi_momentum']:.4g}", flush=True)
            coarse, fine = row["evaluations"]
            row["rate_refinement"] = float(max(abs(fine["rates"]-coarse["rates"])/np.maximum(abs(fine["rates"]), 1.)))
            row["matrix_refinement"] = float(np.max(abs(fine["matrix"]-coarse["matrix"])/np.maximum(abs(fine["matrix"]), 1.)))
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
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
        row["elapsed_seconds"] = time.perf_counter()-started
        checkpoint()
        print(f"trial {n}: passed={row['passed']}, failure={row.get('failure')}, seconds={row['elapsed_seconds']:.1f}", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("source/input changed during pilot experiment")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
