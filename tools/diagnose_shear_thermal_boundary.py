"""Post-screen decomposition, not a redefined gate or selected heat coefficient."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_transverse_mixing import ROOT, actual_source_model
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.shear_thermal import ReducedShearThermal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("screen", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite a boundary diagnosis")
    screen = json.loads(args.screen.read_text(encoding="utf-8"))
    if not screen["all_numerics_passed"]:
        raise RuntimeError("screen numerics must pass before physical interpretation")
    for relative, digest in screen["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"screen code/input changed: {relative}")
    paths = [ROOT/p for p in screen["sha256"]]
    paths += [args.screen.resolve(), Path(__file__).resolve(), ROOT/"docs/shear-thermal-boundary-diagnostic-note.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    ref = ROOT/"reference/preslhy"
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    trials = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))
    measured_rows = {r["trial"]: r for r in measured["trials"]}
    rows = []
    for original in screen["rows"]:
        n = original["trial"]
        jp, th, _ = actual_source_model(next(t for t in trials if t["trial"] == n), measured_rows)
        entry = boundary["interfaces"][str(n)]
        state = np.array(entry["state"])
        section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=entry["thermal_width_ratio"])
        mixing = ConservativeTransverseMixing(section, state)
        family = mixing.tangent_family(np.array(original["sources"]))
        model = ReducedShearThermal(mixing, axial_force_density=lambda q, d: 9.81*(section.rhoa-d["rho"])*math.sin(state[3]))
        at = np.array([1., original["gamma_from_second_moment"]])
        def components(q, order):
            d = model.radial_fields(q, family, thermal_species_ratio=1., order=order)
            local = d["local"]
            advective = (local["h"]/local["rho"])*(d["mass_flux"]@at)
            diffusive = -mixing.area*local["rho"]*local["hq"]*(d["chi_species"]@at)
            return np.column_stack([advective, diffusive])
        edge = {order: 16*mixing.partition.q0*mixing.partition.edge_integrate(lambda q: components(q, order), order) for order in (8, 16)}
        expected = original["at_second_moment_root"]["enthalpy_outward"]
        replay_error = abs(float(sum(edge[16]))-expected)/max(abs(expected), 1.)
        refinement = float(np.max(abs(edge[16]-edge[8])/np.maximum(abs(edge[16]), 1.)))
        q = np.linspace(mixing.partition.q0, mixing.partition.qmax, 1025)
        _, _, temperature, h = section.thermodynamic_profile(state, np.exp(-q))
        hc = section.phase_inverse.enthalpy_and_slope(state[0], state[0]*state[1])[0]
        lower, upper = original["combined_sampled_gamma_interval"]
        if lower is None or upper is None or lower > upper:
            raise RuntimeError("this diagnostic requires a finite nonempty necessary interval")
        affine = original["budgets_affine"]
        r0, production = np.array(affine["zeroth_residual"]), np.array(affine["production"])
        corners = [dict(gamma=gamma, heating_fraction=eta, residual_w_per_m=float((r0+(1-eta)*production)@np.array([1., gamma])))
                   for gamma in (lower, upper) for eta in (0., 1.)]
        bounds = [min(c["residual_w_per_m"] for c in corners), max(c["residual_w_per_m"] for c in corners)]
        row = dict(trial=n, boundary_advective_w_per_m=float(edge[16][0]), boundary_diffusive_w_per_m=float(edge[16][1]),
            boundary_net_w_per_m=float(sum(edge[16])), replay_scaled_error=replay_error, order_refinement=refinement,
            modeled_edge_temperature_k_range=[float(min(temperature)), float(max(temperature))],
            ambient_temperature_k=th.ambient_temperature, modeled_edge_h_over_hc_range=[float(min(h/hc)), float(max(h/hc))],
            combined_necessary_gamma_interval=[lower, upper], reduced_heating_corners=corners,
            residual_bounds_over_gamma_and_heating_fraction=bounds,
            reduced_heating_cannot_close_zeroth_balance=bool(bounds[0] > 0. or bounds[1] < 0.),
            all_numerics_passed=bool(max(replay_error, refinement) <= 1e-5))
        rows.append(row)
        print(f"trial {n}: edge advection={edge[16][0]:.1f}, diffusion={edge[16][1]:.1f} W/m; "
              f"minimum compatible-interval heat defect={bounds[0]:.1f} W/m", flush=True)
    if hashes() != initial:
        raise RuntimeError("input or code changed during diagnosis")
    payload = dict(phase="post_screen_thermal_boundary_decomposition", rows=rows,
        all_numerics_passed=all(r["all_numerics_passed"] for r in rows),
        heating_fraction_selected=False, thermal_species_ratio_fitted=False,
        original_screen_reclassified=False, downstream_run=False, promoted=False, sha256=initial)
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
