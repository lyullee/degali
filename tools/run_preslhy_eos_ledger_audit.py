"""Freeze pure-H2 EOS screening and three source-energy ablations.

No observed concentration is loaded or scored. Profile changes are not an
experimental accuracy claim. Existing output artifacts cannot be overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from CoolProp.CoolProp import PropsSI
from degali.addons.cryogenic_air import multiphase_hydrogen_source_plane
from degali.addons.hydrogen_eos import hydrogen_gas_departure
from degali.addons.lh2_droplets import homogeneous_equilibrium_hydrogen_source
from degali.lh2 import run_lh2_near_field_research


def profile_record(result, distance):
    state = result.solution.state_at_s(float(distance))
    u, rho, fraction, _mw, temperature, rho_h = result.model._profiles(state)
    flux = result.model._fluxes(state)
    comparison = hydrogen_gas_departure(float(temperature[0]), 101325.)
    return {
        "S_m": float(distance), "state": state.tolist(),
        "radial_coordinate_m": (result.model._eta * state[1]).tolist(),
        "velocity_m_s": u.tolist(), "density_kg_m3": rho.tolist(),
        "hydrogen_mass_fraction": fraction.tolist(),
        "temperature_k": temperature.tolist(), "enthalpy_density_j_m3": rho_h.tolist(),
        "fluxes_mass_momentum_x_momentum_y_hydrogen_energy": flux.tolist(),
        "pure_H2_at_total_pressure_comparison_only": asdict(comparison),
    }


def run(output):
    if output.exists():
        raise FileExistsError(f"refusing to overwrite audit artifact: {output}")
    reduced_path = ROOT / "reference/preslhy/e35_reduced.json"
    measured_path = ROOT / "reference/preslhy/measured_pipe_source_2026-09-05.json"
    reduced = json.loads(reduced_path.read_text(encoding="utf-8"))
    measured = json.loads(measured_path.read_text(encoding="utf-8"))
    conditions = {r["trial"]: r for r in reduced["trials"]}
    records = []
    for row in measured["trials"]:
        if row["topology"] != "nozzle":
            continue
        number = row["trial"]
        condition = conditions[number]
        ta = condition["T_C"] + 273.15
        bound = homogeneous_equilibrium_hydrogen_source(
            mass_flow=row["pressure_loss_mass_flow_g_s"] / 1000.,
            orifice_diameter=condition["orifice_mm"] / 1000.,
            upstream_temperature=row["tc3_temperature_k"]["median"],
            upstream_pressure=101325. + row["pt2_barg"]["median"] * 1e5,
            ambient_temperature=ta, y=condition["release_height_m"],
        )
        pipe = bound.postflash
        old_plane = multiphase_hydrogen_source_plane(
            hydrogen_flow=pipe.mass_flow, orifice_diameter=pipe.orifice_diameter,
            orifice_density=pipe.upstream_density, storage_temperature=pipe.upstream_temperature,
            storage_pressure=pipe.upstream_pressure, ambient_temperature=ta,
            specific_momentum=pipe.postflash_velocity, include_kinetic_energy=True,
        )
        old_source = replace(
            bound.source, density=old_plane.density, velocity=old_plane.velocity,
            diameter=old_plane.diameter, mass_fraction=old_plane.endpoint.hydrogen_mass_fraction,
            x=old_plane.formation_distance, enthalpy_boundary=None,
        )
        variants = {
            "legacy_allgas_source_enthalpy": old_source,
            "incoming_kinetic_energy_only": replace(bound.source, enthalpy_boundary=None),
            "phase_energy_ledger": bound.source,
        }
        trial = {"trial": number, "corrected_source": asdict(bound), "variants": {}}
        for name, source in variants.items():
            print(f"trial {number}: {name}", flush=True)
            result = run_lh2_near_field_research(
                source, ambient_temperature=ta,
                relative_humidity=condition["RH_pct"], radial_points=41,
                establishment="entrained_mass", maximum_distance=10. * source.diameter,
                maximum_step=min(0.001, source.diameter / 20.), relative_tolerance=2e-6,
            )
            model = result.model
            grid, tables = model._phase_enthalpy_tables
            h_reference = float(np.interp(ta, grid, tables["Hydrogen"]))
            h_upstream = PropsSI("H", "T", pipe.upstream_temperature,
                                "P", pipe.upstream_pressure, "Hydrogen")
            expected_energy = pipe.mass_flow * (
                h_upstream - h_reference + 0.5 * pipe.upstream_velocity**2
            )
            source_energy = float(model._established_target_fluxes(0.)[3])
            solution = result.solution
            profile = [profile_record(result, distance) for distance in
                       np.linspace(solution.S[0], solution.S[-1], 7)]
            residual = abs(solution.energy_flux[0] - expected_energy) / abs(expected_energy)
            phase_temperature, _ = model._condensed_air_state_exact(
                np.array([source.density]), np.array([source.mass_fraction]),
            )
            trial["variants"][name] = {
                "source": asdict(source),
                "source_energy_flux_w": source_energy,
                "independent_pipe_energy_flux_w": float(expected_energy),
                "pipe_to_establishment_relative_energy_error": float(residual),
                "model_internal_boundary_residual": result.maximum_boundary_residual,
                "model_energy_drift": result.maximum_energy_drift,
                "source_temperature_reconstructed_by_existing_ideal_phase_eos_k": float(phase_temperature[0]),
                "profiles": profile,
            }
            print(f"  pipe energy error {residual:.6g}; Tc {solution.temperature[0]:.4f} -> {solution.temperature[-1]:.4f} K", flush=True)
        records.append(trial)
    eos_states = []
    for species in ("Hydrogen", "ParaHydrogen", "OrthoHydrogen"):
        saturation = float(PropsSI("T", "P", 101325., "Q", 1, species))
        for temperature in (saturation, 25., 33.2, 40., 60., 100., 200., 295.):
            state = hydrogen_gas_departure(temperature, 101325., species)
            eos_states.append({**asdict(state), "relative_volume_departure": state.relative_volume_departure})
    files = [reduced_path, measured_path, Path(__file__),
             ROOT / "docs/prereg-preslhy-source-eos-ledger.md"]
    files += [ROOT / f"src/degali/addons/{name}.py" for name in
              ("axisymmetric_jet", "cryogenic_air", "lh2_droplets", "hydrogen_eos")]
    corrected_errors = [r["variants"]["phase_energy_ledger"]["pipe_to_establishment_relative_energy_error"] for r in records]
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": "docs/prereg-preslhy-source-eos-ledger.md",
        "meaning": "source conservation and pure-component EOS screening, NOT experimental validation or mixture EOS",
        "nearfield_eos_changed": False, "default_dispersion_model_promoted": False,
        "source_ledger_gate_passed": all(e < 1e-8 for e in corrected_errors),
        "max_corrected_pipe_to_establishment_energy_error": max(corrected_errors),
        "file_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "pure_component_eos_states": eos_states, "trials": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"saved: {output}; source ledger gate: {payload['source_ledger_gate_passed']}", flush=True)
    return 0 if payload["source_ledger_gate_passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    raise SystemExit(run(parser.parse_args().output))
