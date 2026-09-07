"""Independent finite changes of flux VALUES, not re-integration of slopes."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from audit_transverse_mixing import ROOT, actual_source_model
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.phase_radial_quadrature import PhaseRadialQuadrature


def displaced_flux(jp, th, state, beta, direction, ds):
    par = np.array([math.log(state[0]), math.log(state[0]*state[1]), math.log(state[2]), state[3],
                    math.log(state[4]), state[5], state[6], math.log(beta)])+ds*direction
    rho, c, area, uc, width = np.exp(par[[0, 1, 2, 4, 7]])
    physical = np.array([rho, c/rho, area, par[3], uc, par[5], par[6]])
    section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=width)
    partition = PhaseRadialQuadrature(section, physical)
    v = section._wind(physical)*math.cos(par[3])
    def values(q):
        rr, yy, _, hh = section.thermodynamic_profile(physical, np.exp(-q))
        u = v+uc*np.exp(-section.velocity_shape_exponent*q)
        return area*np.column_stack([rr*u, rr*yy*u, rr*u*u*math.cos(par[3]),
                                     rr*u*u*math.sin(par[3]), hh*u+.5*rr*u**3])
    return partition.integrate(values, 16, square=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("screen", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite verification")
    screen = json.loads(args.screen.read_text(encoding="utf-8"))
    for relative, digest in screen["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"screen code/input changed: {relative}")
    reference = ROOT/"reference/preslhy"
    boundary = json.loads((reference/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    trials = json.loads((reference/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((reference/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    paths = [args.screen.resolve(), Path(__file__).resolve(), ROOT/"docs/transverse-mixing-verification-note.md"]
    paths += [ROOT/p for p in screen["sha256"]]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    rows = []
    for original in screen["rows"]:
        n = original["trial"]
        start = time.perf_counter()
        jp, th, _ = actual_source_model(next(t for t in trials if t["trial"] == n), measured_rows)
        entry = boundary["interfaces"][str(n)]
        state, beta = np.array(entry["state"]), entry["thermal_width_ratio"]
        scale = np.maximum(abs(np.array(entry["moments"])[:5]), [1e-12, 1e-12, 1., 1., 1.])
        directions = []
        for name, target in (("origin", np.array(original["sources"])), ("response", np.zeros(5))):
            direction = np.array(original[f"tangent_{name}"])
            derivatives = []
            steps = [size/max(np.max(abs(direction)), 1.) for size in (2e-4, 1e-4)]
            for ds in steps:
                plus = displaced_flux(jp, th, state, beta, direction, ds)
                minus = displaced_flux(jp, th, state, beta, direction, -ds)
                derivatives.append((plus-minus)/(2*ds))
            error = abs(derivatives[1]-target)/scale
            change = abs(derivatives[1]-derivatives[0])/scale
            directions.append(dict(name=name, steps_m=steps, flux_derivative=derivatives[1].tolist(),
                target=target.tolist(), scaled_error=error.tolist(), scaled_step_change=change.tolist(),
                passed=bool(max(error) <= 1e-5 and max(change) <= 1e-5)))
        retained_checks = (original["original_source_replay_scaled_error"] <= 1e-8
            and original["boundary_moment_replay_error"] <= 1e-9 and original["linear_flux_residual"] <= 1e-8
            and max(original[k] for k in ("jacobian_8_16_scaled_difference", "transverse_8_16_scaled_difference",
                "mass_species_boundary_scaled_residual", "frozen_thermal_response_scaled_difference")) <= 1e-5)
        passed = bool(retained_checks and all(d["passed"] for d in directions))
        row = dict(trial=n, directions=directions, retained_checks_passed=bool(retained_checks), passed=passed,
                   elapsed_seconds=time.perf_counter()-start)
        rows.append(row)
        print(f"trial {n}: independent flux-change passed={passed}, maximum error "
              f"{max(max(d['scaled_error']) for d in directions):.3e}, elapsed {row['elapsed_seconds']:.1f}s", flush=True)
    if hashes() != initial:
        raise RuntimeError("verification input or code changed during run")
    payload = dict(phase="transverse_mixing_independent_flux_change_verification", rows=rows,
        all_passed=len(rows)==len(screen["selected_trials"]) and all(r["passed"] for r in rows),
        original_dense_gradient_audit_passed=screen["all_numerics_passed"],
        gamma_selected=False, thermal_species_ratio_selected=False, downstream_transport_closed=False,
        field_scored=False, promoted=False, sha256=initial)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
