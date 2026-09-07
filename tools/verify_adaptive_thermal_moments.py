"""Verify the numerical repair against the preserved independent audit."""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from run_preslhy_ambient_profile_audit import ROOT, thermodynamics, hydrogen_gas_jet
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.thermal_moments import EnthalpyMomentOperators
from degali.addons.thermal_moment_quadrature import adaptive_unit_diffusion_response


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("audit", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite verification")
    independent = json.loads(args.audit.read_text(encoding="utf-8"))
    if independent["phase"] != "fixed_section_thermal_moment_operators_only":
        parser.error("requires the frozen thermal-operator audit")
    boundary_path = ROOT / "reference/preslhy/buoyancy_constrained_enthalpy_width_interface_2026-09-05.json"
    reduced_path = ROOT / "reference/preslhy/e35_reduced.json"
    boundary = json.loads(boundary_path.read_text(encoding="utf-8"))
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    for relative, digest in independent["sha256"].items():
        if hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"independent audit input/code changed: {relative}")
    paths = [args.audit, boundary_path, reduced_path, Path(__file__),
             ROOT / "src/degali/addons/thermal_moment_quadrature.py",
             ROOT / "docs/thermal-moment-quadrature-note.md"]
    def hashes():
        return {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):
                hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    rows, cache = [], {}
    for trial in trials:
        number = trial["trial"]
        if number not in independent["selected_trials"]:
            continue
        expected = next(r for r in independent["rows"] if r["trial"] == number)
        entry = boundary["interfaces"][str(number)]
        th = thermodynamics(trial, consistent=True)
        src = th.source  # Only property/geometry reconstruction, never an ODE source.
        jp, _ = hydrogen_gas_jet(
            rate=src.fuel_mass_flow, diameter=src.diameter, velocity=src.velocity,
            wind=trial["wind_ms"], height=trial["release_height_m"], source_temperature=src.temperature,
            source_density=src.density, theta=0., relative_humidity=trial["RH_pct"],
            ambient_temperature=th.ambient_temperature, wind_reference_height=trial["wind_ref_m"],
            sc=1.16**2, density_scaled_entrainment=False, momentum_entrainment_beta=.28, houf_entrainment=True)
        jp.th.ambient.humid, jp.rhoa = th.ambient_absolute_humidity, th.ambient_density
        section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=entry["thermal_width_ratio"])
        section._quadrature_cache = cache.setdefault(section.k.delta, {})
        state = np.array(entry["state"])
        op = EnthalpyMomentOperators(section)
        # Check reconstruction before computing any gradient response.
        moment_error = abs(op.enthalpy_second_moment(state)/expected["analytic_second_moment_W_m2"]-1.)
        if moment_error > 1e-10:
            raise RuntimeError("fixed-section geometry/enthalpy replay failed")
        result = adaptive_unit_diffusion_response(op, state)
        old = np.array(expected["adaptive_moment_volume_boundary"])[1:]
        actual = np.array([result.response.transverse_volume, result.response.outward_boundary_flux])
        error = abs(actual-old)/np.maximum(abs(old), 1.)
        net_error = abs(result.response.derivative-(old[0]-old[1]))/max(abs(old[0]-old[1]), 1.)
        passed = bool(result.converged and expected["adaptive_status"] == "converged"
                      and max(error) <= 1e-5 and net_error <= 1e-5)
        row = {"trial": number, "mean_moment_replay_relative_error": moment_error,
            "adaptive": asdict(result), "net_response": result.response.derivative,
            "independent_component_scaled_differences": error.tolist(),
            "independent_net_scaled_difference": net_error, "passed": passed}
        rows.append(row)
        print(f"trial {number}: adaptive response difference {max(error):.3e}, net {net_error:.3e}, passed={passed}", flush=True)
    if hashes() != initial:
        raise RuntimeError("verification files changed during run")
    payload = {"phase": "adaptive_thermal_moment_numerical_repair", "rows": rows,
        "all_passed": len(rows) == 7 and all(r["passed"] for r in rows),
        "downstream_transport_closed": False, "field_scored": False, "promoted": False,
        "sha256": initial, "upstream_sha256": independent["sha256"],
        "limitations": ["Numerical operator repair, not a fitted diffusivity or completed thermal-width ODE.",
                        "Original failed fixed-order records remain unchanged."]}
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
