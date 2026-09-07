"""Fixed-state thermal-moment audit; no ODE, diffusivity fit, or field score."""

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import cubature

from run_preslhy_ambient_profile_audit import ROOT, replay
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.thermal_moments import EnthalpyMomentOperators


def q_integral(p, delta):
    z = math.sqrt(p*math.pi*delta**2/8.)
    return (2*math.pi*math.erf(z)**2-4*math.sqrt(math.pi)*math.erf(z)*z*math.exp(-z*z))/p**2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("boundary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite a thermal-moment audit")
    control_path = ROOT / "reference/preslhy/phase_ambient_consistency_field_complete_2026-09-05.json"
    reduced_path = ROOT / "reference/preslhy/e35_reduced.json"
    source = json.loads(args.boundary.read_text(encoding="utf-8"))
    control = json.loads(control_path.read_text(encoding="utf-8"))
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    if not source["all_interfaces_accepted"] or source["selected_trials"] != [10, 11, 12, 22, 23, 24, 25]:
        parser.error("requires the accepted frozen seven-interface boundary")
    paths = [args.boundary, control_path, reduced_path, Path(__file__),
        ROOT / "docs/prereg-thermal-moment-operators.md",
        ROOT / "tools/run_preslhy_ambient_profile_audit.py",
        ROOT / "src/degali/addons/thermal_moments.py",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/enthalpy_profile.py",
        ROOT / "src/degali/addons/energy_crosswind.py",
        ROOT / "src/degali/addons/axisymmetric_jet.py",
        ROOT / "src/degali/core/jetplume.py",
        ROOT / "requirements-research.txt"]
    def hashes():
        return {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):
                hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes = hashes()
    rows = []
    for trial in trials:
        number = trial["trial"]
        if number not in source["selected_trials"]:
            continue
        old, replay_check = replay(control, trial)
        entry = source["interfaces"][str(number)]
        model = BuoyancyConstrainedEnthalpySection(old.model.jetplume, old.model.thermodynamics,
            thermal_width_ratio=entry["thermal_width_ratio"], quadrature_points=2048)
        state = np.array(entry["state"])
        op = EnthalpyMomentOperators(model)
        moment = op.enthalpy_second_moment(state)
        coarse = op.unit_diffusion_response(state, quadrature_points=1024)
        fine = op.unit_diffusion_response(state)
        sy, sn = model.section_widths(state)
        area, theta, uc = state[2:5]
        v = model._wind(state)*math.cos(theta)
        hc = float(model.phase_inverse.enthalpy_and_slope(state[0], state[0]*state[1])[0])
        p = 1./model.thermal_width_ratio**2
        analytic = area*(sy*sy+sn*sn)*hc*(v*q_integral(p, model.k.delta)
                   + uc*q_integral(p+model.velocity_shape_exponent, model.k.delta))
        gl = np.array([moment, fine.transverse_volume, fine.outward_boundary_flux])
        scales = np.maximum(abs(gl), 1.)
        length = math.sqrt(math.pi)/2.*model.k.delta
        q0 = .5*length**2
        evaluations = [0]

        def integrand(coordinate):
            t = coordinate[:, 0]
            q = 2.*q0*t
            angular = 2*math.pi-8*np.arccos(np.sqrt(np.minimum(q0/q, 1.)))
            u = v+uc*np.exp(-model.velocity_shape_exponent*q)
            h = hc*np.exp(-p*q)
            coefficient = op.radial_diffusion_coefficient(state, q)
            edge_coefficient = op.radial_diffusion_coefficient(state, q0*(1.+t*t))
            evaluations[0] += 2*len(t)
            return np.column_stack([
                area*2*q0*angular*(sy*sy+sn*sn)*q*h*u,
                area*2*q0*angular*4*q*coefficient,
                4*length**4*area*(2+((sn/sy)**2+(sy/sn)**2)*t*t)*edge_coefficient,
            ])/scales

        info = cubature(integrand, [0.], [1.], points=[[.5]], atol=1e-7, rtol=0.,
                        rule="gk21", max_subdivisions=10000)
        adaptive = info.estimate*scales
        errors = abs(gl-adaptive)/np.maximum(abs(adaptive), 1.)
        analytic_error = abs(moment-analytic)/max(abs(analytic), 1e-12)
        # Counterfactual -grad(H), not a second selectable transport model.
        q, weight = model._quadrature(2048)
        wrong_k = p*hc*np.exp(-p*q)
        wrong_volume = float(area*np.sum(weight*4*q*wrong_k))
        row = {
            "trial": number, "control_replay": replay_check,
            "mean_advective_enthalpy_second_moment_W_m2": moment,
            "analytic_second_moment_W_m2": analytic, "analytic_relative_error": analytic_error,
            "unit_diffusion_response_GL1024": {**asdict(coarse), "net": coarse.derivative},
            "unit_diffusion_response_GL2048": {**asdict(fine), "net": fine.derivative},
            "adaptive_moment_volume_boundary": adaptive.tolist(),
            "scaled_GL_adaptive_differences": errors.tolist(),
            "estimated_scaled_adaptive_error": np.asarray(info.error).tolist(),
            "adaptive_status": info.status, "phase_evaluations": evaluations[0],
            "incorrect_volumetric_gradient_volume_response": wrong_volume,
            "density_gradient_volume_correction_fraction": float((adaptive[1]-wrong_volume)/max(abs(adaptive[1]), 1.)),
            "boundary_to_volume_magnitude_ratio": float(abs(adaptive[2])/max(abs(adaptive[1]), 1.)),
            "passed": bool(info.status == "converged" and analytic_error <= 1e-10 and max(errors) <= 1e-5),
        }
        rows.append(row)
        print(f"trial {number}: M2 analytic {analytic_error:.3e}; GL/adaptive {max(errors):.3e}; "
              f"edge/volume {row['boundary_to_volume_magnitude_ratio']:.3f}; passed={row['passed']}", flush=True)
    final_hashes = hashes()
    if initial_hashes != final_hashes:
        raise RuntimeError("inputs or model changed while audit was running")
    payload = {
        "phase": "fixed_section_thermal_moment_operators_only",
        "moment_order": ["mean_advective_M2_W_m2", "diffusion_volume_per_Dh", "outward_diffusion_boundary_per_Dh"],
        "response_units": "(W m)/(m^2/s); no diffusivity selected",
        "thresholds": {"analytic_moment": 1e-10, "GL_adaptive_each_component": 1e-5},
        "selected_trials": source["selected_trials"],
        "all_passed": len(rows) == 7 and all(r["passed"] for r in rows), "rows": rows,
        "downstream_transport_closed": False, "field_scored": False, "promoted": False,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "sha256": final_hashes, "input_hashes_unchanged_during_run": True,
        "limitations": [
            "Fixed-state operators are not an LH2 downstream transport closure or validation against observations.",
            "Only mean advective axial enthalpy is in M2; axial turbulent flux and mechanical exchange are not supplied.",
            "The eddy specific-enthalpy-gradient response is not molecular conduction or a selected diffusivity.",
            "No ground-cut or torsional geometry is evaluated.",
            "Phase lookup derivative knots may require adaptive integration even when the advective M2 is smooth.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
