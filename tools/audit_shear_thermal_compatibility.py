"""Frozen-boundary falsification screen of unity mixing plus immediate heating."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from audit_transverse_mixing import ROOT, actual_source_model
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.phase_radial_quadrature import PhaseRadialQuadrature
from degali.addons.transverse_mixing import ConservativeTransverseMixing, nonnegative_affine_interval
from degali.addons.shear_thermal import ReducedShearThermal, affine_root


def displaced_thermal_values(jp, th, state, beta, direction, ds):
    """Actual H, mean KE and physical r²H fluxes; no Jacobian integration."""
    par = np.array([math.log(state[0]), math.log(state[0]*state[1]), math.log(state[2]), state[3],
                    math.log(state[4]), state[5], state[6], math.log(beta)])+ds*direction
    rho, c, area, uc, width = np.exp(par[[0, 1, 2, 4, 7]])
    physical = np.array([rho, c/rho, area, par[3], uc, par[5], par[6]])
    section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=width)
    partition = PhaseRadialQuadrature(section, physical)
    v = section._wind(physical)*math.cos(par[3])
    spread = sum(w*w for w in section.section_widths(physical))
    def values(q):
        rr, _, _, hh = section.thermodynamic_profile(physical, np.exp(-q))
        u = v+uc*np.exp(-section.velocity_shape_exponent*q)
        return area*np.column_stack([hh*u, .5*rr*u**3, spread*q*hh*u])
    return partition.integrate(values, 16, square=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, nargs="+")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite a physical compatibility screen")
    ref = ROOT/"reference/preslhy"
    upstream_path = ref/"transverse_mixing_flux_change_verification_2026-09-05.json"
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    if not upstream["all_passed"]:
        raise RuntimeError("upstream physical flux-change verification did not pass")
    for relative, digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"upstream source/input changed: {relative}")
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    trials = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    control = json.loads((ref/"phase_ambient_consistency_field_complete_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    selected = args.trials or [10, 11, 12, 22, 23, 24, 25]
    if not selected or len(set(selected)) != len(selected) or not set(selected).issubset(boundary["selected_trials"]):
        parser.error("select unique frozen boundary trials")
    paths = [ROOT/p for p in upstream["sha256"]]
    paths += [upstream_path, Path(__file__).resolve(), ROOT/"src/degali/addons/shear_thermal.py",
              ROOT/"tests/test_shear_thermal.py", ROOT/"docs/prereg-shear-thermal-compatibility.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes = hashes()
    rows, caches = [], {}
    for n in selected:
        start = time.perf_counter()
        print(f"trial {n}: actual-source replay and explicit equilibrium hypothesis", flush=True)
        trial = next(t for t in trials if t["trial"] == n)
        jp, th, source = actual_source_model(trial, measured_rows)
        old = control["interfaces"][str(n)]["downstream"]
        replay = IndependentEnergyCrosswind(jp, th).source_terms(np.array(old["states"][0]))
        source_error = float(np.max(abs(replay-old["sources"][0])/np.maximum(abs(np.array(old["sources"][0])), 1.)))
        entry = boundary["interfaces"][str(n)]
        state, beta = np.array(entry["state"]), entry["thermal_width_ratio"]
        section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=beta,
                                                    quadrature_points=entry["quadrature_points"])
        section._quadrature_cache = caches.setdefault(section.k.delta, {})
        replay_moments = section.moments(state)
        moment_error = float(np.max(abs(replay_moments-entry["moments"])/section.moment_scales(entry["moments"])))
        if source_error > 1e-8 or moment_error > 1e-9:
            raise RuntimeError("frozen actual source/boundary replay failed")
        sources = IndependentEnergyCrosswind.source_terms(section, state)
        mixing = ConservativeTransverseMixing(section, state)
        model = ReducedShearThermal(mixing, axial_force_density=lambda q, d: 9.81*(section.rhoa-d["rho"])*math.sin(state[3]))
        families = {order: mixing.tangent_family(sources, order=order) for order in (8, 16)}
        budgets = {order: model.equilibrium_budgets(families[order], thermal_species_ratio=1., order=order) for order in (8, 16)}
        b, family = budgets[16], families[16]
        refinement = max(float(np.max(abs(budgets[8][key]-value)/np.maximum(abs(value), 1.))) for key, value in b.items())
        kinetic_scale = np.maximum.reduce([np.ones(2)]+[abs(b[key]) for key in ("kinetic_axial", "kinetic_outward", "body_work", "production")])
        kinetic_error = float(np.max(abs(b["kinetic_identity"])/kinetic_scale))
        energy_identity_error = float(np.max(abs(b["total_energy_discrepancy"]-b["zeroth_residual"])/np.maximum(abs(b["zeroth_residual"]), 1.)))
        gamma2, gamma0 = affine_root(b["second_residual"]), affine_root(b["zeroth_residual"])
        root_refinement = abs(gamma2-affine_root(budgets[8]["second_residual"]))/max(abs(gamma2), 1.)
        at = np.array([1., gamma2])
        evaluated = {key: float(value@at) for key, value in b.items()}
        scale0 = max(1., *(abs(evaluated[key]) for key in ("enthalpy_axial", "enthalpy_outward", "production")))
        scale2 = max(1., *(abs(evaluated[key]) for key in ("weighted_enthalpy_axial", "weighted_enthalpy_outward", "weighted_transverse", "weighted_production")))
        error0, error2 = abs(evaluated["zeroth_residual"])/scale0, abs(evaluated["second_residual"])/scale2
        knots = mixing.partition.knots
        q = np.unique(np.r_[np.linspace(0., mixing.partition.qmax, 1025), .5*(knots[:-1]+knots[1:])])
        d = model.radial_fields(q, family, thermal_species_ratio=1.)
        chi_c, chi_p = d["chi_species"]@at, d["chi_momentum"]@at
        positive = bool(np.all(chi_c >= 0.) and np.all(chi_p >= 0.))
        interval = nonnegative_affine_interval(np.r_[d["chi_species"][:, 0], d["chi_momentum"][:, 0]],
                                               np.r_[d["chi_species"][:, 1], d["chi_momentum"][:, 1]])
        rates = family.at(gamma2)
        curvature = abs(rates[3])*math.sqrt(2*mixing.partition.q0)*mixing.sn
        iq = mixing.partition.integrate(lambda qq: qq*mixing.area*mixing.local(qq)["h"]*mixing.local(qq)["u"], 16, square=True)
        spread_rate = 2*mixing.sy**2*(mixing.log_sy_partials@rates)+2*mixing.sn**2*(mixing.log_sn_partials@rates)
        expected = np.array([evaluated["enthalpy_axial"], evaluated["kinetic_axial"],
                             evaluated["weighted_enthalpy_axial"]+spread_rate*iq])
        flux_values = displaced_thermal_values(jp, th, state, beta, rates, 0.)
        fd_scale = np.maximum(abs(flux_values), 1.)  # current flux / 1 m, not heat-source scale
        steps = [s/max(np.max(abs(rates)), 1.) for s in (2e-4, 1e-4)]
        derivatives = [(displaced_thermal_values(jp, th, state, beta, rates, ds)
                       -displaced_thermal_values(jp, th, state, beta, rates, -ds))/(2*ds) for ds in steps]
        fd_error = float(np.max(abs(derivatives[1]-expected)/fd_scale))
        fd_change = float(np.max(abs(derivatives[1]-derivatives[0])/fd_scale))
        numeric = bool(family.maximum_scaled_residual <= 1e-8 and max(refinement, kinetic_error,
            energy_identity_error, root_refinement, fd_error, fd_change) <= 1e-5)
        compatible = bool(error0 <= 1e-5 and error2 <= 1e-5 and positive and curvature < 1.)
        row = dict(trial=n, source=source, original_source_replay_error=source_error, boundary_moment_replay_error=moment_error,
            sources=sources.tolist(), linear_flux_residual=family.maximum_scaled_residual,
            budgets_affine={key: value.tolist() for key, value in b.items()},
            budget_8_16_scaled_difference=refinement, kinetic_identity_scaled_error=kinetic_error,
            energy_identity_scaled_error=energy_identity_error, root_8_16_scaled_difference=root_refinement,
            gamma_from_second_moment=gamma2, gamma_from_zeroth_moment=gamma0,
            at_second_moment_root=evaluated, zeroth_budget_scale_w_per_m=scale0, second_budget_scale_w_m=scale2,
            zeroth_scaled_residual=error0, second_scaled_residual=error2,
            positive_diffusion=positive, minimum_chi_species_per_s=float(min(chi_c)), minimum_chi_momentum_per_s=float(min(chi_p)),
            sample_q=q.tolist(), chi_species=chi_c.tolist(), chi_momentum=chi_p.tolist(),
            combined_sampled_gamma_interval=[v if math.isfinite(v) else None for v in interval],
            curvature_half_width=curvature, thermal_width_rate_candidate=rates.tolist(),
            flux_change=dict(steps_m=steps, expected_derivative=expected.tolist(), derivatives=[v.tolist() for v in derivatives],
                scaled_error=fd_error, scaled_step_change=fd_change, normalization="current individual flux / 1 m; floor 1"),
            numerics_passed=numeric, local_hypothesis_compatible=compatible,
            eligible_for_short_segment=bool(numeric and compatible and curvature < .1), elapsed_seconds=time.perf_counter()-start)
        rows.append(row)
        print(f"trial {n}: numeric={numeric}, compatible={compatible}, gamma2={gamma2:.6g}, gamma0={gamma0:.6g}, "
              f"heat residual={error0:.3e}, positive={positive}, elapsed={row['elapsed_seconds']:.1f}s", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("input or code changed during the audit")
    payload = dict(phase="reduced_shear_unity_diffusion_local_equilibrium_screen", selected_trials=selected, rows=rows,
        thermal_species_ratio=1., heating_hypothesis="instantaneous conversion of reconstructed axial shear production",
        all_numerics_passed=len(rows)==len(selected) and all(r["numerics_passed"] for r in rows),
        all_hypotheses_compatible=len(rows)==len(selected) and all(r["local_hypothesis_compatible"] for r in rows),
        downstream_transport_closed=False, short_segment_run=False, field_scored=False, promoted=False, sha256=initial_hashes,
        limitations=["Reduced even axial shear balance, not full curved Favre RANS or a TKE transport state.",
            "Unity turbulent scalar diffusion and instantaneous shear heating are explicit hypotheses, not fitted or validated LH2 laws.",
            "R2 root is only a candidate until independent R0, positivity and geometry gates pass.",
            "No existing total-energy source, input source, observation score or default model was changed."])
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
