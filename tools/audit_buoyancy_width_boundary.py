"""Independent adaptive Gauss--Kronrod integration of the frozen six moments."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cubature

from run_preslhy_ambient_profile_audit import ROOT, replay
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("boundary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite an independent boundary audit")
    source = json.loads(args.boundary.read_text(encoding="utf-8"))
    control_path = ROOT / "reference/preslhy/phase_ambient_consistency_field_complete_2026-09-05.json"
    reduced_path = ROOT / "reference/preslhy/e35_reduced.json"
    control = json.loads(control_path.read_text(encoding="utf-8"))
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    rows = []
    for trial in trials:
        n = trial["trial"]
        if n not in source["selected_trials"]:
            continue
        old, check = replay(control, trial)
        entry = source["interfaces"][str(n)]
        model = BuoyancyConstrainedEnthalpySection(old.model.jetplume, old.model.thermodynamics,
            thermal_width_ratio=entry["thermal_width_ratio"])
        state = np.array(entry["state"])
        theta, area, uc = state[3], state[2], state[4]
        v = model._wind(state)*math.cos(theta)
        q0 = math.pi*model.k.delta**2/8.
        expected = np.array(entry["moments"])
        scales = model.moment_scales(np.array(entry["target_moments"]))

        evaluations = [0]
        def integrand(coordinates):
            # Direct q-coordinate integration with the square's angular
            # measure, independent of the solver's two mapped GL node sets.
            q = coordinates[:, 0]
            evaluations[0] += len(q)
            angular = 2*math.pi-8*np.arccos(np.sqrt(np.minimum(q0/q, 1.)))
            rho, y, _, h = model.thermodynamic_profile(state, np.exp(-q))
            u = v+uc*np.exp(-model.velocity_shape_exponent*q)
            return area*angular[:, None]*np.column_stack([rho*u, rho*y*u, rho*u*u*math.cos(theta),
                rho*u*u*math.sin(theta), h*u+.5*rho*u**3, 9.81*(model.rhoa-rho)])/scales

        info = cubature(integrand, [0.], [2*q0], points=[[q0]], atol=1e-8,
            rtol=0., rule="gk21", max_subdivisions=5000)
        actual = info.estimate*scales
        difference = np.abs(actual-expected)/scales
        row = {"trial": n, "control_replay": check, "adaptive_moments": actual.tolist(),
               "difference_from_frozen_scaled": difference.tolist(),
               "estimated_scaled_error": np.asarray(info.error).tolist(), "evaluations": evaluations[0],
               "integrator_success": info.status == "converged", "message": info.status,
               "passed": bool(info.status == "converged" and max(difference) <= 1e-5)}
        rows.append(row)
        print(f"trial {n}: adaptive six-moment difference {max(difference):.3e}, passed={row['passed']}", flush=True)
    paths = [args.boundary, control_path, reduced_path, Path(__file__),
        ROOT / "tools/run_preslhy_ambient_profile_audit.py", ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/enthalpy_profile.py", ROOT / "src/degali/addons/axisymmetric_jet.py"]
    payload = {"method": "independent_vectorized_adaptive_GK21_direct_q_square", "threshold": 1e-5,
        "passed": len(rows) == 7 and all(r["passed"] for r in rows), "rows": rows,
        "sha256": {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):
                   hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "limitations": ["This checks quadrature of fixed states, not downstream transport or measured accuracy.",
                        "The near-field target moments are unchanged; their original discretization is not re-solved."]}
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
