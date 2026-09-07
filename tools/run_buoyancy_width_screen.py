"""Freeze a seven-case six-moment boundary screen; never score field data."""

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from importlib.metadata import version
import platform
from pathlib import Path

import numpy as np

from run_preslhy_ambient_profile_audit import ROOT, thermodynamics, hydrogen_gas_jet
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path,
        default=ROOT / "reference/preslhy/gaussian_enthalpy_profile_allflux_interface_2026-09-05.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite a boundary result")
    checkpoint = json.loads(args.checkpoint.read_text(encoding="utf-8"))
    selected = [10, 11, 12, 22, 23, 24, 25]
    required = {"phase": "interface", "mode": "full", "droplet_equilibrium_bound": True,
                "hydrogen_spin_isomer": "normal", "energy_quadrature": "polar_square_all_flux_v1",
                "phase_ambient_closure": "consistent_explicit_ideal_v1",
                "downstream_thermodynamic_profile": "enthalpy", "interfaces_accepted": True}
    if (any(checkpoint.get(k) != v for k, v in required.items())
            or checkpoint.get("selected_trials") != selected or checkpoint.get("failures")):
        parser.error("checkpoint does not match the frozen seven-case experiment")
    reduced_path = ROOT / "reference/preslhy/e35_reduced.json"
    expected_hash = checkpoint["provenance_sha256"][str(reduced_path.relative_to(ROOT))]
    if hashlib.sha256(reduced_path.read_bytes()).hexdigest() != expected_hash:
        parser.error("reduced input changed since the source checkpoint")
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    names = ["total_mass", "hydrogen", "momentum_x", "momentum_z", "energy"]
    rows, quadrature_caches = {}, {}
    for trial in trials:
        n = trial["trial"]
        if n not in selected:
            continue
        print(f"trial {n}: match heat, species and buoyancy with independent thermal width", flush=True)
        old = checkpoint["interfaces"][str(n)]
        th = thermodynamics(trial, consistent=True)
        src = th.source  # Property/geometry factory only; no near-field or ODE solve.
        jp, _ = hydrogen_gas_jet(
            rate=src.fuel_mass_flow, diameter=src.diameter, velocity=src.velocity,
            wind=trial["wind_ms"], height=trial["release_height_m"], source_temperature=src.temperature,
            source_density=src.density, theta=0., relative_humidity=trial["RH_pct"],
            ambient_temperature=th.ambient_temperature, wind_reference_height=trial["wind_ref_m"],
            sc=1.16**2, density_scaled_entrainment=False, momentum_entrainment_beta=.28, houf_entrainment=True)
        jp.th.ambient.humid = th.ambient_absolute_humidity
        jp.rhoa = th.ambient_density
        model = BuoyancyConstrainedEnthalpySection(jp, th, quadrature_points=256)
        model._quadrature_cache = quadrature_caches.setdefault(model.k.delta, {})
        initial = np.array(old["state"])
        represented = model.moments(initial, quadrature_points=old["energy_quadrature_points"])
        expected = np.array([old["projected_fluxes"][k] for k in names] + [old["projected_buoyancy_force"]])
        replay_error = float(np.max(np.abs(represented-expected)/model.moment_scales(expected)))
        if replay_error > 1e-9:
            raise RuntimeError(f"trial {n} property/geometry replay mismatch {replay_error}")
        target = np.array([old["target_fluxes"][k] for k in names] + [old["target_buoyancy_force"]])
        try:
            result = model.project_buoyancy(target, initial)
            row = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in asdict(result).items()}
            width = math.sqrt(2*math.log(2)*result.state[2])
            temperature = model.centre_temperature(result.state)
            width_error = abs(width/old["target_halfwidth_m"]-1.)
            temperature_error = abs(temperature-old["target_temperature_K"])
            reasons = []
            if not result.success:
                reasons.append(result.message)
            if width_error > .05:
                reasons.append("hydrogen half-width mismatch exceeds 5%")
            if temperature_error > 2.:
                reasons.append("centre-temperature mismatch exceeds 2 K")
            if result.at_width_bound:
                reasons.append("thermal-width solution touches a search boundary")
            row.update({"target_moments": target.tolist(), "target_halfwidth_m": old["target_halfwidth_m"],
                "projected_halfwidth_m": width, "halfwidth_residual": width_error,
                "target_temperature_K": old["target_temperature_K"], "projected_temperature_K": temperature,
                "temperature_residual_K": temperature_error, "accepted": not reasons,
                "failure_reasons": reasons, "same_width_replay_error": replay_error})
        except (ValueError, RuntimeError, FloatingPointError) as error:
            row = {"accepted": False, "failure_reasons": [str(error)], "same_width_replay_error": replay_error,
                   "target_moments": target.tolist()}
        rows[str(n)] = row
        print(f"trial {n}: accepted={row['accepted']}, beta={row.get('thermal_width_ratio')}, "
              f"temperature residual={row.get('temperature_residual_K')} K, reasons={row['failure_reasons']}", flush=True)
    hashed = [args.checkpoint, reduced_path, Path(__file__),
        ROOT / "requirements-research.txt",
        ROOT / "tools/run_preslhy_ambient_profile_audit.py",
        ROOT / "src/degali/addons/buoyancy_profile.py", ROOT / "src/degali/addons/enthalpy_profile.py",
        ROOT / "src/degali/addons/energy_crosswind.py", ROOT / "src/degali/addons/axisymmetric_jet.py",
        ROOT / "docs/prereg-buoyancy-constrained-enthalpy-width.md"]
    payload = {"phase": "boundary_only", "selected_trials": selected,
        "moment_order": names+["buoyancy_force_N_per_m"], "interfaces": rows,
        "all_interfaces_accepted": len(rows) == 7 and all(r["accepted"] for r in rows.values()),
        "downstream_transport_implemented": False, "field_scored": False, "promoted": False,
        "runtime": {"python": platform.python_version(), "platform": platform.platform(),
                    **{name: version(name) for name in ("numpy", "scipy", "CoolProp")}},
        "sha256": {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):
                   hashlib.sha256(p.read_bytes()).hexdigest() for p in hashed},
        "limitations": ["Buoyancy targets are frozen near-field moments, not measured forces.",
                        "Boundary agreement supplies no downstream thermal-width transport law.",
                        "No observed temperature/concentration is fitted or scored."]}
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"saved {args.output}; accepted {sum(r['accepted'] for r in rows.values())}/7", flush=True)


if __name__ == "__main__":
    main()
