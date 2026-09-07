"""Short, guarded continuation after independently passed initializations."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK, independent_rays
from degali.addons.enriched_transport import PrescribedRadialMixing, EnrichedModalTransport
from degali.addons.enriched_segments import EnrichedShortSegment, forward_shape_horizon, assert_transport_domain


def verify_endpoint(driver, segment, initial_fluxes):
    value = segment["parameters"]
    outputs = [driver.evaluate(value, order=o, angular_order=a, exact=True, enforce=False) for o, a in ((4, 8), (8, 16))]
    coarse, fine = outputs
    rate_change = float(max(abs(fine["rates"]-coarse["rates"])/np.maximum(abs(fine["rates"]), 1.)))
    matrix_change = float(np.max(abs(fine["matrix"]-coarse["matrix"])/np.maximum(abs(fine["matrix"]), 1.)))
    p, par, supplied, factor = driver.context(value)
    model = EnrichedModalTransport(p, par, scalar_mixing=supplied, thermal_species_ratio=1., mechanical_work=WORK)
    rays = independent_rays(model, fine["rates"])
    final_fluxes = driver.flux_values(value)
    balance = (final_fluxes-initial_fluxes-segment["cumulative_sources"])/np.maximum(abs(initial_fluxes), 1.)
    physical, failure = True, None
    try:
        assert_transport_domain(fine)
        if (max(rays["edge_defects"].values()) > .05 or rays["minimum_chi_momentum"] < 0.
                or rays["maximum_outward_mass"] >= 0.):
            raise ValueError("independent fixed rays failed their physical gate")
    except ValueError as exc:
        physical, failure = False, str(exc)
    return dict(evaluations=outputs, rate_refinement=rate_change, matrix_refinement=matrix_change,
        independent_rays=rays, final_fluxes=final_fluxes, scaled_balance=balance,
        maximum_balance=float(max(abs(balance))), mixing_amplitude_ratio=factor,
        numerics_passed=bool(max(rate_change, matrix_change, float(max(abs(balance)))) <= 1e-5
            and max(fine["linear_scaled_error"], fine["weak_heat_scaled_error"]) <= 1e-8
            and abs(fine["mass_boundary_error"])/max(abs(fine["source"][0]), 1.) <= 1e-8),
        physical_passed=physical, failure=failure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, nargs="+", default=[10, 11, 12, 22, 23, 24, 25])
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite segment evidence")
    ref = ROOT/"reference/preslhy"
    origin = ref/"coupled_initialization_combined_2026-09-06.json"
    data = json.loads(origin.read_text(encoding="utf-8"))
    if not data["completed"] or not data["all_passed"]:
        raise ValueError("first finish and inspect all seven initialization cases")
    if len(set(args.trials)) != len(args.trials) or not set(args.trials).issubset(data["selected_trials"]):
        parser.error("select unique initialized trials")
    for name, digest in data["sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"frozen initialization dependency changed: {name}")
    paths = [ROOT/n for n in data["sha256"]]+[origin, Path(__file__).resolve(),
        ROOT/"src/degali/addons/enriched_segments.py", ROOT/"tests/test_enriched_segments.py",
        ROOT/"docs/prereg-enriched-short-segments.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes = hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding="utf-8"))
    frozen = read("edge_conservative_refit_combined_2026-09-05.json")
    boundary = read("buoyancy_constrained_enthalpy_width_interface_2026-09-05.json")
    reduced = read("e35_reduced.json")["trials"]
    measured = {r["trial"]: r for r in read("measured_pipe_source_2026-09-05.json")["trials"]}
    rows, caches = [], {}
    began = datetime.now(timezone.utc).isoformat()
    def payload(completed):
        return dict(phase="guarded_enriched_short_segments", completed=completed, selected_trials=args.trials,
            rows=rows, passed_trials=[r["trial"] for r in rows if r.get("primary_passed", False)],
            sha256=initial_hashes, process_id=os.getpid(), started_utc=began,
            updated_utc=datetime.now(timezone.utc).isoformat(), full_downstream_integrated=False,
            field_scored=False, promoted=False,
            limitations=["Local segments with explicit conditional similarity mixing, not measured validation.",
                "No endpoint refit, clipped stress or silently relaxed shape/BC bounds.",
                "Alternative is a sensitivity run only, not replacement of a failed primary policy."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    checkpoint()
    for n in args.trials:
        began_case = time.perf_counter()
        row = dict(trial=n, primary_passed=False, segments=[], stage="source_reconstruction")
        rows.append(row)
        checkpoint()
        print(f"trial {n}: guarded new-shape downstream start", flush=True)
        try:
            seed = next(r for r in frozen["rows"] if r["trial"] == n)
            initial = next(r for r in data["rows"] if r["trial"] == n)
            p, row["source_reconstruction"] = construct_projection(seed, boundary, reduced, measured, caches)
            par = np.array(initial["parameters"])
            supplied = PrescribedRadialMixing.from_weak_baseline(p, order=32)
            horizon = forward_shape_horizon(p, par, np.array(initial["evaluations"][-1]["rates"]))
            length = min(.005, .05*min(p.mixing.sy, p.mixing.sn), .2*horizon)
            if not math.isfinite(length) or length <= 0.:
                raise ValueError("no positive admissible starting segment")
            row.update(length_m=length, linear_shape_horizon_m=horizon if math.isfinite(horizon) else None)
            kwargs = dict(scalar_mixing=supplied, thermal_species_ratio=1., mechanical_work=WORK)
            driver = EnrichedShortSegment(p, par, mixing_update="equilibrium_velocity_width", **kwargs)
            initial_fluxes = driver.flux_values(driver.initial)
            row["initial_state"], row["initial_fluxes"] = driver.initial.copy(), initial_fluxes
            for steps in (2, 4, 8):
                row["stage"] = f"primary_{steps}_steps"
                checkpoint()
                def progress(item):
                    row["latest_step"] = item
                    checkpoint()
                    print(f"trial {n}, {steps} divisions: {item['step']} at {item['reached_m']:.6g}m, "
                        f"edge={item['edge_defects']}, chi={item['minimum_chi_momentum']:.4g}", flush=True)
                out = driver.integrate(length, steps, callback=progress)
                row["segments"].append(out)
                checkpoint()
                print(f"trial {n}, {steps}: completed={out['completed']}, failure={out.get('failure')}", flush=True)
            if all(s["completed"] for s in row["segments"]):
                ends = [s["parameters"] for s in row["segments"]]
                changes = [float(max(abs(ends[i+1]-ends[i])/np.maximum(abs(ends[i+1]), 1.))) for i in (0, 1)]
                row["step_refinement"] = changes
                row["stage"] = "primary_independent_endpoint"
                checkpoint()
                verification = verify_endpoint(driver, row["segments"][-1], initial_fluxes)
                row["endpoint"] = verification
                row["primary_passed"] = bool(changes[-1] <= 1e-5 and verification["numerics_passed"] and verification["physical_passed"])
                checkpoint()
            if row["primary_passed"]:
                row["stage"] = "alternative_8_steps"
                checkpoint()
                alternative = EnrichedShortSegment(p, par, mixing_update="constant_geometric_diffusivity", **kwargs)
                alt = alternative.integrate(length, 8)
                row["alternative"] = alt
                checkpoint()
                if alt["completed"]:
                    row["stage"] = "alternative_independent_endpoint"
                    checkpoint()
                    alt["verification"] = verify_endpoint(alternative, alt, initial_fluxes)
                    primary = row["segments"][-1]["parameters"]
                    alt["relative_endpoint_difference"] = float(max(abs(alt["parameters"]-primary)/np.maximum(abs(primary), 1.)))
                    alt["step_refinement_checked"] = False
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
        row["stage"] = "completed"
        row["elapsed_seconds"] = time.perf_counter()-began_case
        checkpoint()
        print(f"trial {n}: primary={row['primary_passed']}, failure={row.get('failure')}, seconds={row['elapsed_seconds']:.1f}", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("code/input changed during guarded segments")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
