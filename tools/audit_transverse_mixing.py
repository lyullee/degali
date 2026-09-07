"""Actual-source boundary tangent/mixing audit, with gamma left unknown."""

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from run_preslhy_ambient_profile_audit import ROOT, thermodynamics, hydrogen_gas_jet
from degali.addons.axisymmetric_jet import ConservedGaussianJet
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.thermal_moments import EnthalpyMomentOperators
from degali.addons.phase_radial_quadrature import gauss_rule
from degali.addons.transverse_mixing import ConservativeTransverseMixing, nonnegative_affine_interval
from degali.validation.nearfield import hydrogen_jet, _trial_source_inputs, _measured_lh2_equilibrium_source
from degali.core.jetplume import J_UC


def actual_source_model(trial, measured_rows):
    rate, temperature, pressure = _trial_source_inputs(trial, source="flow_mean_gs",
        liquid_temperature_mode="tank_saturation", measured_rows=measured_rows, measured_source_mode="full")
    setup, _ = hydrogen_jet(rate=rate, diameter=trial["orifice_mm"]/1000., wind=trial["wind_ms"],
        height=trial["release_height_m"], ambient_temperature=trial["T_C"]+273.15,
        relative_humidity=trial["RH_pct"], storage_pressure_barg=trial["tanker_barg"],
        storage_temperature=temperature, source_upstream_pressure_barg=pressure,
        wind_reference_height=trial["wind_ref_m"], corrections=True, hydrogen_spin_isomer="normal")
    src = _measured_lh2_equilibrium_source(enabled=True, setup=setup, trial=trial, rate=rate,
        source_temperature=temperature, source_pressure_barg=pressure, measured_rows=measured_rows,
        measured_source_mode="full", hydrogen_spin_isomer="normal")
    jp, initial = hydrogen_gas_jet(rate=src.fuel_mass_flow, diameter=src.diameter, velocity=src.velocity,
        wind=trial["wind_ms"], height=trial["release_height_m"], source_temperature=src.temperature,
        source_density=src.density, source_mass_fraction=src.mass_fraction, theta=0.,
        relative_humidity=trial["RH_pct"], ambient_temperature=trial["T_C"]+273.15,
        wind_reference_height=trial["wind_ref_m"], sc=1.16**2, density_scaled_entrainment=False,
        momentum_entrainment_beta=.28, houf_entrainment=True)
    local_wind = src.velocity-float(initial[J_UC])
    properties = thermodynamics(trial, consistent=True)
    # Construct source-dependent coefficients from the ACTUAL source/coflow.
    th = ConservedGaussianJet(src, ambient_temperature=properties.ambient_temperature,
        ambient_pressure=properties.ambient_pressure, ambient_density=properties.ambient_density,
        fuel_molecular_weight=properties.fuel_molecular_weight, ambient_molecular_weight=properties.ambient_molecular_weight,
        fuel_heat_capacity=properties.fuel_heat_capacity, ambient_heat_capacity=properties.ambient_heat_capacity,
        ambient_absolute_humidity=properties.ambient_absolute_humidity, ambient_coflow_velocity=local_wind,
        equilibrium_air_condensation=True, temperature_dependent_phase_enthalpy=True,
        consistent_phase_ambient=True, conservative_establishment="entrained_mass", radial_points=41)
    jp.th.ambient.humid, jp.rhoa = th.ambient_absolute_humidity, th.ambient_density
    return jp, th, dict(source=asdict(src), local_wind=local_wind,
                       source_froude=th._source_froude, buoyancy_coefficient=th._buoyancy_coefficient)


def dense_square_nodes(q0, panels=4096):
    nodes, weight = gauss_rule(4)
    edges = np.linspace(0., 1., panels+1)
    s = (edges[:-1, None]+np.diff(edges)[:, None]*nodes).ravel()
    w = (np.diff(edges)[:, None]*weight).ravel()
    return np.r_[q0*s, q0*(1+s*s)], np.r_[2*math.pi*q0*w, (2*math.pi-8*np.arctan(s))*2*q0*s*w]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, nargs="+")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite a transverse-mixing audit")
    names = ["e35_reduced.json", "measured_pipe_source_2026-09-05.json",
        "buoyancy_constrained_enthalpy_width_interface_2026-09-05.json",
        "phase_ambient_consistency_field_complete_2026-09-05.json",
        "thermal_moment_adaptive_verification_2026-09-05.json"]
    paths = [ROOT/"reference/preslhy"/n for n in names]
    reduced, measured, boundary, control, thermal_reference = [json.loads(p.read_text(encoding="utf-8")) for p in paths]
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    selected = args.trials or [10, 11, 12, 22, 23, 24, 25]
    if not selected or len(set(selected)) != len(selected) or not set(selected).issubset(boundary["selected_trials"]):
        parser.error("select unique frozen boundary trials")
    paths += [Path(__file__), ROOT/"tools/run_preslhy_ambient_profile_audit.py",
        ROOT/"docs/prereg-transverse-conservative-mixing.md", ROOT/"requirements-research.txt"]
    paths += [ROOT/"src/degali"/p for p in ["addons/phase_radial_quadrature.py", "addons/transverse_mixing.py",
        "addons/axisymmetric_jet.py", "addons/lh2_droplets.py", "addons/buoyancy_profile.py", "addons/enthalpy_profile.py",
        "addons/energy_crosswind.py", "addons/thermal_moments.py", "core/jetplume.py", "validation/nearfield.py"]]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes = hashes()
    rows, caches = [], {}
    for trial in reduced["trials"]:
        n = trial["trial"]
        if n not in selected:
            continue
        start = time.perf_counter()
        print(f"trial {n}: rebuild actual source and retain unknown thermal-width rate", flush=True)
        jp, th, source = actual_source_model(trial, measured_rows)
        old_entry = control["interfaces"][str(n)]
        old = IndependentEnergyCrosswind(jp, th)
        old_state = np.array(old_entry["downstream"]["states"][0])
        expected_source = np.array(old_entry["downstream"]["sources"][0])
        replay_source = old.source_terms(old_state)
        source_error = float(np.max(abs(replay_source-expected_source)/np.maximum(abs(expected_source), 1.)))
        if source_error > 1e-8:
            raise RuntimeError(f"trial {n}: original actual-source RHS replay failed: {source_error}")
        entry = boundary["interfaces"][str(n)]
        section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=entry["thermal_width_ratio"],
                                                   quadrature_points=entry["quadrature_points"])
        section._quadrature_cache = caches.setdefault(section.k.delta, {})
        state = np.array(entry["state"])
        moments = section.moments(state)
        moment_replay = float(np.max(abs(moments-np.array(entry["moments"]))/section.moment_scales(entry["moments"])))
        if moment_replay > 1e-9:
            raise RuntimeError("actual-source section does not replay the frozen boundary moments")
        # Deliberate base evaluation for a FIXED-STATE diagnostic only; this
        # does not expose or enable the boundary class's blocked ODE methods.
        sources = IndependentEnergyCrosswind.source_terms(section, state)
        op = ConservativeTransverseMixing(section, state)
        coarse, fine = op.tangent_family(sources, order=8), op.tangent_family(sources, order=16)
        jac_error = float(np.max(abs(coarse.jacobian-fine.jacobian)/np.maximum(abs(fine.jacobian), 1.)))
        dq, dw = dense_square_nodes(op.partition.q0)
        independent_jac = np.tensordot(dw, op.flux_rate_density(dq), axes=(0, 0))
        independent_jac_error = float(np.max(abs(independent_jac-fine.jacobian)/np.maximum(abs(fine.jacobian), 1.)))
        q = np.linspace(0., op.partition.qmax, 513)
        a, b = op.mixing_family(q, coarse, order=8), op.mixing_family(q, fine, order=16)
        transverse_error = max(float(np.max(abs(a[key]-b[key])/np.maximum(abs(b[key]), 1.))) for key in ("mass", "species"))
        interval = nonnegative_affine_interval(b["chi"][:, 0], b["chi"][:, 1])
        # Outward ALE flux on four square faces versus axial M/H2 change.
        outward = 16.*op.partition.q0*op.partition.edge_integrate(lambda qq: op.transverse_basis(qq), 16)
        boundary_error = float(np.max(abs(outward+fine.jacobian[:2])/np.maximum(abs(fine.jacobian[:2]), 1.)))
        thermal = EnthalpyMomentOperators(section)
        volume = op.area*op.partition.integrate(lambda qq: 4*qq*thermal.radial_diffusion_coefficient(state, qq), 16, square=True)
        length = math.sqrt(2*op.partition.q0)
        shape_ratio = (op.sn/op.sy)**2+(op.sy/op.sn)**2
        edge = 4*length**4*op.area*op.partition.edge_integrate(
            lambda qq: (2+shape_ratio*(qq/op.partition.q0-1))*thermal.radial_diffusion_coefficient(state, qq), 16)
        frozen = next(r for r in thermal_reference["rows"] if r["trial"] == n)["adaptive"]["response"]
        expected = np.array([frozen["transverse_volume"], frozen["outward_boundary_flux"]])
        thermal_error = float(np.max(abs(np.array([volume, edge])-expected)/np.maximum(abs(expected), 1.)))
        passed = (fine.maximum_scaled_residual <= 1e-8 and max(jac_error, independent_jac_error,
                  transverse_error, boundary_error, thermal_error) <= 1e-5)
        row = dict(trial=n, source=source, original_source_replay_scaled_error=source_error,
            boundary_moment_replay_error=moment_replay, sources=sources.tolist(),
            tangent_origin=fine.origin.tolist(), tangent_response=fine.response.tolist(),
            linear_flux_residual=fine.maximum_scaled_residual, phase_panels=len(op.partition.knots)-1,
            jacobian_8_16_scaled_difference=jac_error, dense_independent_jacobian_scaled_difference=independent_jac_error,
            transverse_8_16_scaled_difference=transverse_error, mass_species_boundary_scaled_residual=boundary_error,
            frozen_thermal_response_scaled_difference=thermal_error,
            sampled_gamma_interval=[v if math.isfinite(v) else None for v in interval],
            sampled_positive_diffusion_feasible=bool(interval[0] <= interval[1]),
            sample_q=q.tolist(), chi_origin=b["chi"][:, 0].tolist(), chi_response=b["chi"][:, 1].tolist(),
            numerics_passed=bool(passed), elapsed_seconds=time.perf_counter()-start)
        rows.append(row)
        print(f"trial {n}: numerical={passed}, gamma interval={interval}, "
              f"independent Jacobian={independent_jac_error:.3e}, boundary={boundary_error:.3e}, "
              f"elapsed={row['elapsed_seconds']:.1f}s", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("input or code changed during the audit")
    payload = dict(phase="conditional_conservative_transverse_mixing", selected_trials=selected, rows=rows,
        all_numerics_passed=len(rows)==len(selected) and all(r["numerics_passed"] for r in rows),
        gamma_selected=False, thermal_species_diffusivity_ratio_selected=False,
        downstream_transport_closed=False, field_scored=False, promoted=False,
        sha256=initial_hashes, limitations=[
            "Radial flux in deforming coordinates is an explicit shape ansatz, not a unique continuity result.",
            "Gamma intervals are sampled necessary conditions, not a selected ODE or all-radius positivity proof.",
            "No thermal/species diffusion ratio, mechanical heating distribution or ground closure was selected."])
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
