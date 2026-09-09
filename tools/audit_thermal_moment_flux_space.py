"""Audit Trial 10 with six physical fluxes as the primary RK4 state."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_reservoir_thermal_segments import actual_source_model, serial
from run_preslhy_ambient_profile_audit import replay
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.reservoir_thermal import ReservoirShortSegment, decode_section, encode_section
from degali.addons.thermal_moment_flux_march import (
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
)


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_parameters(parameters):
    state, beta = decode_section(parameters)
    return np.r_[state, beta]


def independent_values(driver, parameters):
    evaluation = driver.evaluate(parameters, order=16, probes=1025)
    return driver.flux_values(evaluation["mixing"], order=16), evaluation


def result_record(result, driver, initial_values_16):
    data = asdict(result)
    terminal_values_16, terminal = independent_values(driver, result.parameters[-1])
    scales_8 = np.maximum(np.abs(result.fluxes[0]), [1e-12, 1e-12, 1., 1., 1., 1.])
    scales_16 = np.maximum(np.abs(initial_values_16), [1e-12, 1e-12, 1., 1., 1., 1.])
    direct = (result.fluxes[-1]-result.fluxes[0]-result.cumulative_sources[-1])/scales_8
    independent = (terminal_values_16-initial_values_16-result.cumulative_sources[-1])/scales_16
    state, beta = decode_section(result.parameters[-1])
    data.update(
        final_arc_length_m=float(result.arc_length[-1]),
        final_downwind_m=float(state[5]),
        final_elevation_m=float(state[6]),
        final_thermal_width_ratio=beta,
        terminal_values_order16=terminal_values_16,
        direct_scaled_balance=direct,
        maximum_direct_scaled_balance=float(np.max(np.abs(direct))),
        independent_order16_scaled_balance=independent,
        maximum_independent_order16_scaled_balance=float(np.max(np.abs(independent))),
        endpoint_weak_budget_scaled_error=terminal["weak_budget_scaled_error"],
        endpoint_minimum_chi_species=terminal["minimum_chi_species"],
        endpoint_minimum_chi_momentum=terminal["minimum_chi_momentum"],
        balance_passed=bool(np.max(np.abs(independent)) <= 1e-5),
    )
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=REF / "thermal_moment_flux_space_2026-09-09.json",
    )
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite final or partial flux-space evidence")

    boundary_path = REF / "buoyancy_constrained_enthalpy_width_interface_2026-09-05.json"
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    control_path = REF / "phase_ambient_consistency_field_complete_2026-09-05.json"
    dependencies = [
        boundary_path,
        reduced_path,
        measured_path,
        control_path,
        ROOT / "docs/prereg-thermal-moment-flux-space-march.md",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        ROOT / "src/degali/addons/reservoir_thermal.py",
        ROOT / "src/degali/addons/transverse_mixing.py",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/thermal_moments.py",
        ROOT / "tests/test_thermal_moment_flux_march.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}

    boundary = json.loads(boundary_path.read_text(encoding="utf-8"))
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    measured = json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    control = json.loads(control_path.read_text(encoding="utf-8"))
    trial = next(row for row in trials if row["trial"] == 10)
    measured_rows = {row["trial"]: row for row in measured}
    jp, thermodynamics, source = actual_source_model(trial, measured_rows)
    _trajectory, control_check = replay(control, trial)

    entry = boundary["interfaces"]["10"]
    initial_state = np.asarray(entry["state"], float)
    initial_beta = float(entry["thermal_width_ratio"])
    section = BuoyancyConstrainedEnthalpySection(
        jp, thermodynamics, thermal_width_ratio=initial_beta,
        quadrature_points=entry["quadrature_points"],
    )
    replayed = section.moments(initial_state)
    boundary_error = float(np.max(
        np.abs(replayed-np.asarray(entry["moments"]))/section.moment_scales(entry["moments"])
    ))
    if boundary_error > 1e-9:
        raise RuntimeError("frozen six-moment boundary replay failed")

    initial_parameters = encode_section(initial_state, initial_beta)
    driver = ReservoirShortSegment(
        jp, thermodynamics, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work",
    )
    initial_values_16, initial_evaluation = independent_values(driver, initial_parameters)
    if not initial_evaluation["valid"]:
        raise RuntimeError("frozen initial local closure is invalid")

    def run(label, step):
        inverter = ThermalMomentFluxInverter(
            jp, thermodynamics, order=8, probes=257,
            quadrature_points=128, tolerance=1e-8,
        )
        marcher = FluxSpaceThermalMomentMarch(
            inverter, thermal_species_ratio=1.,
            mechanical_work="reduced_buoyancy_work",
        )
        def progress(count, arc, primary):
            if count == 1 or count % 25 == 0:
                print(f"{label}: step={count}, arc={arc:.6f} m, x={primary[6]:.6f} m", flush=True)
        return marcher.march(
            initial_parameters, 1.25, step=step, target_x=1.78,
            maximum_steps=4000, target_tolerance=1e-8, progress=progress,
        )

    payload = dict(
        completed=False,
        trial=10,
        source=source,
        target_downwind_m=1.78,
        thermal_species_ratio=1.,
        thermal_species_ratio_fitted=False,
        mechanical_work="reduced_buoyancy_work",
        primary_state=["mass", "hydrogen", "momentum_x", "momentum_z", "total_energy", "enthalpy_second_moment"],
        march_method="classical_rk4",
        inverse_tolerance=1e-8,
        independent_balance_limit=1e-5,
        convergence_limit=.005,
        yaw_included=False,
        ground_interaction_included=False,
        source_files_distributed=False,
        sensor_score_calculated=False,
        candidate_promoted=False,
        control_replay=control_check,
        boundary_moment_replay_error=boundary_error,
        initial_values_order16=initial_values_16,
        hashes=before,
    )

    coarse = run("coarse", .005)
    payload["coarse"] = result_record(coarse, driver, initial_values_16)
    partial.write_text(json.dumps(payload, default=serial, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(
        f"coarse complete: reached={coarse.reached_target}, "
        f"balance16={payload['coarse']['maximum_independent_order16_scaled_balance']:.3e}",
        flush=True,
    )

    refined = None
    if coarse.reached_target and payload["coarse"]["balance_passed"]:
        refined = run("refined", .0025)
        payload["refined"] = result_record(refined, driver, initial_values_16)
        coarse_physical = physical_parameters(coarse.parameters[-1])
        refined_physical = physical_parameters(refined.parameters[-1])
        parameter_scale = np.maximum(np.maximum(abs(coarse_physical), abs(refined_physical)), 1.)
        flux_scale = np.maximum(
            np.maximum(abs(payload["coarse"]["terminal_values_order16"]),
                       abs(payload["refined"]["terminal_values_order16"])),
            [1e-12, 1e-12, 1., 1., 1., 1.],
        )
        parameter_difference = float(np.max(abs(coarse_physical-refined_physical)/parameter_scale))
        flux_difference = float(np.max(
            abs(np.asarray(payload["coarse"]["terminal_values_order16"])
                -np.asarray(payload["refined"]["terminal_values_order16"]))/flux_scale
        ))
        payload["convergence"] = dict(
            maximum_physical_parameter_scaled_difference=parameter_difference,
            maximum_flux_scaled_difference=flux_difference,
            limit=.005,
            passed=bool(refined.reached_target and max(parameter_difference, flux_difference) <= .005),
        )
        partial.write_text(json.dumps(payload, default=serial, indent=2, allow_nan=False)+"\n", encoding="utf-8")

    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during the flux-space audit")
    payload["hashes"] = after
    payload["completed"] = True
    payload["finished_utc"] = datetime.now(timezone.utc).isoformat()
    payload["all_gates_passed"] = bool(
        coarse.reached_target
        and payload["coarse"]["balance_passed"]
        and refined is not None
        and payload["refined"]["balance_passed"]
        and payload["convergence"]["passed"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    partial.unlink()
    print(json.dumps({
        "completed": payload["completed"],
        "all_gates_passed": payload["all_gates_passed"],
        "coarse_stop": coarse.stop_reason,
        "refined_stop": None if refined is None else refined.stop_reason,
        "coarse_balance16": payload["coarse"]["maximum_independent_order16_scaled_balance"],
        "refined_balance16": None if refined is None else payload["refined"]["maximum_independent_order16_scaled_balance"],
        "convergence": payload.get("convergence"),
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
