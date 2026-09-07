"""Replay frozen fields at thermocouples without integrating the ODE again.

No raw workbook is modified. Reconstructed thermodynamics/geometry must
reproduce stored temperatures, five fluxes and sensor-arc concentrations
before any comparison is emitted. Missing old trajectories stay missing.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from CoolProp.CoolProp import PropsSI
from degali.addons.axisymmetric_jet import (
    AxisymmetricJetSource, ConservedGaussianJet, phase_ambient_from_rh,
)
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.validation.nearfield import IndependentEnergyTrajectory, hydrogen_gas_jet
from degali.validation.preslhy import TemperatureReading, read_nearfield_temperatures


def thermodynamics(trial, *, consistent, argon=False):
    ta, pa, rh = trial["T_C"] + 273.15, 101325., trial["RH_pct"]
    mw_air, mw_water = PropsSI("M", "Air"), PropsSI("M", "Water")
    pv = rh / 100. * PropsSI("P", "T", ta, "Q", 0, "Water")
    humidity = mw_water / mw_air * pv / (pa - pv)
    mw_humid = (1 + humidity) / (1 / mw_air + humidity / mw_water)
    rho = PropsSI("D", "T", ta, "P", pa, "Air") * mw_humid / mw_air
    if consistent:
        mw_air, humidity, rho = phase_ambient_from_rh(ta, pa, rh, include_argon=argon)
    # The source below only instantiates a thermodynamic property evaluator.
    # It is NOT a replacement release, and never enters a source/ODE solve.
    source = AxisymmetricJetSource(
        diameter=.01, velocity=100., temperature=100., density=PropsSI("D", "T", 100., "P", pa, "Hydrogen"),
        theta=0., y=trial["release_height_m"],
    )
    return ConservedGaussianJet(
        source, ambient_temperature=ta, ambient_pressure=pa, ambient_density=rho,
        fuel_molecular_weight=PropsSI("M", "Hydrogen"), ambient_molecular_weight=mw_air,
        fuel_heat_capacity=PropsSI("C", "T", ta, "P", pa, "Hydrogen"),
        ambient_heat_capacity=PropsSI("C", "T", ta, "P", pa, "Air"),
        ambient_absolute_humidity=humidity, equilibrium_air_condensation=True,
        temperature_dependent_phase_enthalpy=True, equilibrium_argon_condensation=argon,
        conservative_establishment="entrained_mass", radial_points=41,
        consistent_phase_ambient=consistent,
    )


def replay(field, trial):
    profile = field.get("downstream_thermodynamic_profile", "density")
    if profile not in {"density", "enthalpy"}:
        raise ValueError("replay cannot reinterpret an unknown thermodynamic profile")
    if (field.get("mode") != "full" or not field.get("droplet_equilibrium_bound")
            or field.get("hydrogen_spin_isomer", "normal") != "normal"
            or field.get("source_energy_ledger") != "moving_pipe_phase_enthalpy_v2"):
        raise ValueError("replay supports only the frozen normal-H2/full/HEM energy-ledger experiment")
    closure = field.get("phase_ambient_closure", "legacy_air_eos")
    if closure not in {"legacy_air_eos", "consistent_explicit_ideal_v1"}:
        raise ValueError(f"unsupported ambient closure: {closure}")
    entry = field["interfaces"].get(str(trial["trial"]), {})
    stored = entry.get("downstream")
    if stored is None:
        return None, {"status": "missing_stored_trajectory"}
    thermo = thermodynamics(trial, consistent=closure == "consistent_explicit_ideal_v1")
    src = thermo.source
    jp, _ = hydrogen_gas_jet(
        rate=src.fuel_mass_flow, diameter=src.diameter, velocity=src.velocity,
        wind=trial["wind_ms"], height=trial["release_height_m"], source_temperature=src.temperature,
        source_density=src.density, theta=0., relative_humidity=trial["RH_pct"],
        ambient_temperature=thermo.ambient_temperature, wind_reference_height=trial["wind_ref_m"],
        sc=1.16**2, density_scaled_entrainment=False, momentum_entrainment_beta=.28,
        houf_entrainment=True,
    )
    jp.th.ambient.humid = thermo.ambient_absolute_humidity
    jp.rhoa = thermo.ambient_density
    model_class = IndependentEnergyCrosswind
    if profile == "enthalpy":
        from degali.addons.enthalpy_profile import GaussianEnthalpyCrosswind
        model_class = GaussianEnthalpyCrosswind
    model = model_class(jp, thermo, quadrature_points=entry["energy_quadrature_points"])
    model.velocity_shape_exponent = entry["lambda"]**2
    result = SimpleNamespace(states=np.asarray(stored["states"]))
    trajectory = IndependentEnergyTrajectory(SimpleNamespace(model=model), result)
    states = result.states
    t, _ = thermo._condensed_air_state(states[:, 0], states[:, 1])
    t_error = float(np.max(np.abs(t - stored["temperatures"])))
    indices = np.unique(np.linspace(0, len(states)-1, 7, dtype=int))
    flux_errors = []
    for index in indices:
        actual = model._as_array(model.integral_fluxes(states[index]))
        expected = np.asarray(stored["fluxes"][index])
        flux_errors.append(float(np.max(np.abs(actual-expected) / np.maximum(np.abs(expected), 1.))))
    arc_errors = []
    for pair in field["candidate"]["pairs"]:
        if pair["trial"] != trial["trial"]:
            continue
        sensors = [s for s in trial["sensors"] if round(s["x"], 3) == pair["x"]]
        predicted = max(trajectory.concentration_at(pair["x"], s["y"], s["z"]) for s in sensors)
        arc_errors.append(abs(predicted / pair["predicted"] - 1.))
    checks = {"temperature_max_abs_K": t_error, "sampled_flux_max_scaled_error": max(flux_errors),
              "arc_max_relative_error": max(arc_errors), "flux_sections_checked": len(indices)}
    if t_error > 1e-6 or max(flux_errors) > 1e-9 or max(arc_errors) > 1e-9:
        raise RuntimeError(f"trial {trial['trial']} replay does not match stored fields: {checks}")
    return trajectory, checks


def statistics(rows, model_key, observed_key):
    errors = np.array([r[model_key] - r[observed_key] for r in rows if r.get(model_key) is not None])
    if not errors.size:
        return {"n": 0}
    return {"n": len(errors), "mean_bias_K": float(np.mean(errors)),
            "median_bias_K": float(np.median(errors)), "mae_K": float(np.mean(np.abs(errors))),
            "median_absolute_error_K": float(np.median(np.abs(errors))),
            "rmse_K": float(np.sqrt(np.mean(errors**2)))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--control", type=Path, default=ROOT / "reference/preslhy/source_phase_energy_ledger_field_complete_2026-09-05.json")
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--temperature-reference", type=Path,
                        help="reuse frozen sensor reductions after verifying original workbook checksums")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite an existing audit")
    reduced = ROOT / "reference/preslhy/e35_reduced.json"
    trials = json.loads(reduced.read_text(encoding="utf-8"))["trials"]
    fields = {"control": json.loads(args.control.read_text(encoding="utf-8")),
              "candidate": json.loads(args.candidate.read_text(encoding="utf-8"))}
    frozen_temperature = (json.loads(args.temperature_reference.read_text(encoding="utf-8"))
                          if args.temperature_reference is not None else None)
    checks, trajectories, rows, profiles, radial, raw = {}, {}, [], [], [], {}
    for trial in trials:
        number = trial["trial"]
        if number not in fields["candidate"]["selected_trials"]:
            continue
        trajectories[number] = {}
        for kind, field in fields.items():
            trajectory, check = replay(field, trial)
            checks[f"{kind}_{number}"] = check
            trajectories[number][kind] = trajectory
            print(f"replay {kind} trial {number}: {check}", flush=True)
            if trajectory is None:
                continue
            for pair in field["candidate"]["pairs"]:
                if pair["trial"] != number:
                    continue
                state = trajectory.state_at(pair["x"])
                sy, sz = trajectory.model.section_widths(state)
                profiles.append({"variant": kind, "trial": number, "x_m": pair["x"],
                    "predicted_over_observed": pair["predicted"] / pair["observed"],
                    "centre_T_K": trajectory.model.centre_temperature(state),
                    "centre_density_kg_m3": float(state[0]), "centre_Y_H2": float(state[1]),
                    "sigma_y_m": sy, "sigma_z_m": sz, "centre_z_m": float(state[6])})
                if number in (11, 12, 23, 24) and pair["x"] in (1.19, 4., 6.):
                    shape = np.array([1., .5, .1, .01, 1e-4, 1e-6])
                    th = trajectory.model.thermodynamics
                    density = th.ambient_density + (state[0] - th.ambient_density) * shape
                    fuel_density = state[0] * state[1] * shape
                    temp, rho_h = th._condensed_air_state(density, fuel_density / density)
                    if field.get("downstream_thermodynamic_profile", "density") == "enthalpy":
                        density, fraction, temp, rho_h = trajectory.model.thermodynamic_profile(state, shape)
                    radial.append({"variant": kind, "trial": number, "x_m": pair["x"],
                        "Gaussian_scalar_shape": shape.tolist(), "density_kg_m3": density.tolist(),
                        "T_K": temp.tolist(), "h_per_H2_J_kg": (rho_h / fuel_density).tolist()})
        if number not in (10, 23):
            continue
        workbook = next((ROOT / "reference/preslhy/raw").glob(f"trial_{number}_*alldata.xlsx"))
        digest_before = hashlib.sha256(workbook.read_bytes()).hexdigest()
        if frozen_temperature is None:
            window, readings = read_nearfield_temperatures(workbook, release_height=trial["release_height_m"])
        else:
            original = frozen_temperature["raw_workbooks"][str(number)]
            if original["file"] != workbook.name or original["sha256"].lower() != digest_before:
                raise RuntimeError("temperature reference does not match the original workbook checksum")
            window = original["data_window_zero_based_half_open"]
            readings = [TemperatureReading(**{key: row[key] for key in TemperatureReading.__dataclass_fields__})
                        for row in frozen_temperature["temperature_rows"] if row["trial"] == number]
            if len({r.channel for r in readings}) != len(readings):
                raise ValueError("temperature reference contains duplicated channels")
        if hashlib.sha256(workbook.read_bytes()).hexdigest() != digest_before:
            raise RuntimeError("raw workbook changed during read-only audit")
        raw[str(number)] = {"file": workbook.name, "sha256": digest_before,
                            "sheet": "Flexlogger", "data_window_zero_based_half_open": window,
                            "excel_rows_inclusive": [window[0]+2, window[1]+1]}
        for reading in readings:
            available = {kind: tr for kind, tr in trajectories[number].items()
                         if tr is not None and tr.state_at(reading.x) is not None}
            if "candidate" not in available:
                continue
            rows.append({**asdict(reading), "trial": number,
                "observed_minimum_K": reading.minimum_c + 273.15,
                "observed_p05_K": reading.percentile_05_c + 273.15,
                "observed_median_K": reading.median_c + 273.15,
                **{f"{kind}_K": (available[kind].temperature_at(reading.x, reading.y, reading.z)
                                  if kind in available else None) for kind in fields}})
    summaries = {}
    for label, subset in (("all_available", rows),
                          ("paired", [r for r in rows if r["control_K"] is not None]),
                          ("trial_10", [r for r in rows if r["trial"] == 10]),
                          ("trial_23", [r for r in rows if r["trial"] == 23])):
        summaries[label] = {basis: {kind: statistics(subset, f"{kind}_K", f"observed_{basis}_K")
                                    for kind in fields} for basis in ("minimum", "p05", "median")}
    fixed_points = []
    trial23 = next(t for t in trials if t["trial"] == 23)
    for argon in (False, True):
        for consistent in (False, True):
            th = thermodynamics(trial23, consistent=consistent, argon=argon)
            for name in ("exact", "interpolated"):
                inverse = th._condensed_air_state_exact if name == "exact" else th._condensed_air_state
                t, rho_h = inverse(np.array([th.ambient_density]), np.array([0.]))
                fixed_points.append({"argon": argon, "consistent": consistent, "inverse": name,
                    "ambient_density_kg_m3": th.ambient_density, "ambient_humidity_kg_kg": th.ambient_absolute_humidity,
                    "delta_T_K": float(t[0] - th.ambient_temperature), "rho_h_J_m3": float(rho_h[0])})
    hashed = [args.control, args.candidate, reduced, Path(__file__),
              ROOT / "src/degali/addons/axisymmetric_jet.py", ROOT / "src/degali/addons/energy_crosswind.py",
              ROOT / "src/degali/addons/enthalpy_profile.py",
              ROOT / "src/degali/validation/preslhy.py", ROOT / "src/degali/validation/nearfield.py"]
    if args.temperature_reference is not None:
        hashed.append(args.temperature_reference)
    payload = {"status": "research_only", "replay_checks": checks, "fixed_points": fixed_points,
        "temperature_summaries": summaries, "temperature_rows": rows, "axial_profiles": profiles,
        "radial_profiles": radial,
        "raw_workbooks": raw, "sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in hashed},
        "limitations": ["No trajectory ODE was rerun or calibrated by this audit.",
            "Missing control trajectories are not imputed; compare only paired temperature rows.",
            "Temperature minima, p05 and medians share the frozen release window but are not synchronized with gas-sampling extrema.",
            "No observed mixture density is inferred from non-simultaneous temperature/concentration extrema."]}
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"temperature": summaries, "saved": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
