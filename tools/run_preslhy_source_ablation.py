"""Run and freeze the pre-registered PRESLHY measured-source ablations."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.validation.nearfield import (  # noqa: E402
    CoupledFieldValidation,
    independent_energy_from_reduced,
    independent_energy_interfaces_from_reduced,
)


DEFAULT_REDUCED = ROOT / "reference" / "preslhy" / "e35_reduced.json"
DEFAULT_SOURCE = (
    ROOT / "reference" / "preslhy" / "measured_pipe_source_2026-09-05.json"
)


def _interfaces(result) -> dict[str, dict]:
    return {
        str(trial): {
            "station_m": interface.streamline_distance,
            "lambda": interface.velocity_spreading_ratio,
            "width_residual": interface.halfwidth_residual,
            "temperature_residual_k": interface.temperature_residual,
            "max_flux_residual": max(interface.relative_residuals.values()),
            "accepted": interface.accepted,
            "failure_reasons": interface.failure_reasons,
            "state": interface.state.tolist(),
            "target_fluxes": interface.target_fluxes,
            "projected_fluxes": interface.projected_fluxes,
            "target_temperature_K": interface.target_temperature,
            "projected_temperature_K": interface.projected_temperature,
            "target_halfwidth_m": interface.target_halfwidth,
            "projected_halfwidth_m": interface.projected_halfwidth,
            "energy_quadrature_points": interface.energy_quadrature_points,
            "energy_quadrature_residual": interface.energy_quadrature_residual,
            "target_buoyancy_force": interface.target_buoyancy_force,
            "projected_buoyancy_force": interface.projected_buoyancy_force,
        }
        for trial, interface in sorted(result.interfaces.items())
    }


def _near_field(result) -> dict:
    if not result.pairs:
        return {
            "arcs": 0,
            "MG": math.nan,
            "VG": math.nan,
            "FAC2": math.nan,
        }
    statistics = result.statistics()
    return {
        "arcs": len(result.pairs),
        "MG": statistics.mg,
        "VG": statistics.vg,
        "FAC2": statistics.fac2,
    }


def _field(result: CoupledFieldValidation) -> dict:
    def trajectory(interface):
        downstream = getattr(interface, "downstream_result", None)
        if downstream is None:
            return None
        return {
            name: value.tolist() if hasattr(value, "tolist") else value
            for name, value in vars(downstream).items()
        }

    baseline = _near_field(result.baseline)
    baseline["geometry"] = result.geometry(result.baseline_vertical)
    baseline["pairs"] = [asdict(pair) for pair in result.baseline.pairs]
    baseline["vertical_profiles"] = result.baseline_vertical
    candidate = _near_field(result.candidate)
    candidate["geometry"] = result.geometry(result.candidate_vertical)
    candidate["pairs"] = [asdict(pair) for pair in result.candidate.pairs]
    candidate["vertical_profiles"] = result.candidate_vertical
    return {
        "selected_trials": result.selected_trials,
        "interfaces_accepted": result.all_interfaces_accepted,
        "failures": result.failures,
        "baseline": baseline,
        "candidate": candidate,
        "vertical_rows": len(result.candidate_vertical),
        "promoted": result.promoted,
        "interfaces": {
            str(trial): {
                "station_m": interface.streamline_distance,
                "lambda": interface.velocity_spreading_ratio,
                "width_residual": interface.halfwidth_residual,
                "temperature_residual_k": interface.temperature_residual,
                "max_flux_residual": max(interface.relative_residuals.values()),
                "accepted": interface.accepted,
                "failure_reasons": interface.failure_reasons,
                "energy_quadrature_points": interface.energy_quadrature_points,
                "downstream": trajectory(interface),
            }
            for trial, interface in sorted(result.handoffs.items())
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("interface", "field"), required=True)
    parser.add_argument(
        "--mode", choices=("flow_only", "nozzle_only", "full"), required=True
    )
    parser.add_argument(
        "--droplet-equilibrium-bound", action="store_true",
        help=(
            "start measured nozzles at the collective all-H2-vapour "
            "equilibrium plane"
        ),
    )
    parser.add_argument("--reduced", type=Path, default=DEFAULT_REDUCED)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--maximum-step", type=float, default=0.02)
    parser.add_argument("--trials", nargs="+", type=int)
    parser.add_argument(
        "--hydrogen-spin-isomer",
        choices=("normal", "para", "ortho"),
        default="normal",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--consistent-phase-ambient", action="store_true")
    parser.add_argument("--downstream-profile", choices=("density", "enthalpy"), default="density")
    parser.add_argument("--interface-checkpoint", type=Path)
    args = parser.parse_args()
    if args.downstream_profile == "enthalpy":
        if not args.consistent_phase_ambient:
            parser.error("enthalpy profile requires --consistent-phase-ambient")
        if args.phase == "field":
            if args.interface_checkpoint is None:
                parser.error("enthalpy field requires a passed seven-case --interface-checkpoint")
            checkpoint = json.loads(args.interface_checkpoint.read_text(encoding="utf-8"))
            if (not checkpoint.get("interfaces_accepted") or checkpoint.get("failures")
                    or checkpoint.get("selected_trials") != [10, 11, 12, 22, 23, 24, 25]
                    or checkpoint.get("downstream_thermodynamic_profile") != "enthalpy"
                    or checkpoint.get("energy_quadrature") != "polar_square_all_flux_v1"
                    or checkpoint.get("phase_ambient_closure") != "consistent_explicit_ideal_v1"
                    or checkpoint.get("mode") != args.mode
                    or checkpoint.get("droplet_equilibrium_bound") != args.droplet_equilibrium_bound
                    or checkpoint.get("hydrogen_spin_isomer") != args.hydrogen_spin_isomer):
                parser.error("interface checkpoint does not match this seven-case physics experiment")
            for path in (args.reduced.resolve(), args.source.resolve()):
                key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
                if checkpoint.get("provenance_sha256", {}).get(key) != hashlib.sha256(path.read_bytes()).hexdigest():
                    parser.error("interface checkpoint input checksum does not match")
    if args.output is not None and args.output.exists():
        parser.error(f"refusing to overwrite existing result: {args.output}")

    common = {
        "measured_source_table": args.source,
        "measured_source_mode": args.mode,
        "measured_lh2_equilibrium_bound": args.droplet_equilibrium_bound,
        "hydrogen_spin_isomer": args.hydrogen_spin_isomer,
        "fit_velocity_spreading": args.downstream_profile == "density",
        "consistent_phase_ambient": args.consistent_phase_ambient,
        "downstream_thermodynamic_profile": args.downstream_profile,
        "progress": lambda message: print(message, flush=True),
    }
    if args.phase == "interface":
        result = independent_energy_interfaces_from_reduced(args.reduced, **common)
        payload = {
            "phase": "interface",
            "mode": args.mode,
            "droplet_equilibrium_bound": args.droplet_equilibrium_bound,
            "hydrogen_spin_isomer": args.hydrogen_spin_isomer,
            "selected_trials": result.selected_trials,
            "interfaces_accepted": result.all_interfaces_accepted,
            "failures": result.failures,
            "interfaces": _interfaces(result),
        }
        print(result.report(), flush=True)
    else:
        result = independent_energy_from_reduced(
            args.reduced,
            maximum_step=args.maximum_step,
            trial_filter=args.trials,
            **common,
        )
        payload = {
            "phase": "field",
            "mode": args.mode,
            "droplet_equilibrium_bound": args.droplet_equilibrium_bound,
            "hydrogen_spin_isomer": args.hydrogen_spin_isomer,
            **_field(result),
        }
        print(result.report(), flush=True)

    payload["source_energy_ledger"] = (
        "moving_pipe_phase_enthalpy_v2" if args.droplet_equilibrium_bound else "not_applicable"
    )
    payload["energy_quadrature"] = (
        "polar_square_all_flux_v1" if args.downstream_profile == "enthalpy" else "exact_symmetric_nodes_v1"
    )
    if args.interface_checkpoint is not None:
        payload["interface_checkpoint_sha256"] = hashlib.sha256(args.interface_checkpoint.read_bytes()).hexdigest()
    payload["downstream_thermodynamic_profile"] = args.downstream_profile
    payload["phase_ambient_closure"] = (
        "consistent_explicit_ideal_v1" if args.consistent_phase_ambient else "legacy_air_eos"
    )
    payload["provenance_sha256"] = {
        str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path):
        hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (
            args.reduced.resolve(), args.source.resolve(), Path(__file__).resolve(),
            ROOT / "src/degali/addons/axisymmetric_jet.py",
            ROOT / "src/degali/addons/energy_crosswind.py",
            ROOT / "src/degali/addons/enthalpy_profile.py",
            ROOT / "src/degali/lh2.py", ROOT / "src/degali/validation/nearfield.py",
            ROOT / "docs/prereg-phase-ambient-consistency.md",
            ROOT / "docs/prereg-gaussian-enthalpy-profile.md",
        )
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
        print(f"saved: {args.output}", flush=True)
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
