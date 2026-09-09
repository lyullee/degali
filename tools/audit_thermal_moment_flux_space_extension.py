"""Extend the direct six-flux thermal march to frozen PRESLHY trials."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_reservoir_thermal_segments import actual_source_model, serial
from audit_thermal_moment_flux_space import (
    independent_values,
    physical_parameters,
    result_record,
)
from run_preslhy_ambient_profile_audit import replay
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.reservoir_thermal import ReservoirShortSegment, encode_section
from degali.addons.thermal_moment_flux_march import (
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
)


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"
TRIALS = (11, 12, 22, 23, 24, 25)
FROZEN_TARGETS = {11: .79, 12: 1.78, 22: 1.78, 23: .79, 24: 1.78, 25: 1.78}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def first_downstream_profile_target(trial, boundary_x):
    fits = trial.get("vertical_fits", [])
    candidates = sorted({
        float(row["x"])
        for row in fits
        if row.get("well_constrained") and float(row["x"]) > boundary_x+1e-9
    })
    if not candidates:
        raise ValueError("no well-constrained vertical profile lies downstream of the boundary")
    return candidates[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--trials", nargs="+", type=int, default=list(TRIALS))
    args = parser.parse_args()
    selected = tuple(args.trials)
    if len(set(selected)) != len(selected) or not set(selected).issubset(TRIALS):
        parser.error("select unique trials from 11,12,22,23,24,25")
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite final or partial extension evidence")

    boundary_path = REF / "buoyancy_constrained_enthalpy_width_interface_2026-09-05.json"
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    control_path = REF / "phase_ambient_consistency_field_complete_2026-09-05.json"
    dependencies = [
        boundary_path,
        reduced_path,
        measured_path,
        control_path,
        ROOT / "docs/prereg-thermal-moment-flux-space-extension.md",
        ROOT / "docs/prereg-flux-space-rk-stage-manifold.md",
        ROOT / "docs/prereg-flux-space-endpoint-step-rejection.md",
        ROOT / "docs/prereg-flux-space-cell-transition-adaptivity.md",
        ROOT / "docs/prereg-analytic-enthalpy-width-split.md",
        ROOT / "docs/prereg-analytic-width-six-trial-extension.md",
        ROOT / "src/degali/addons/energy_crosswind.py",
        ROOT / "src/degali/addons/enthalpy_profile.py",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        ROOT / "src/degali/addons/reservoir_thermal.py",
        ROOT / "src/degali/addons/transverse_mixing.py",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/thermal_moments.py",
        ROOT / "tools/audit_thermal_moment_flux_space.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    boundary = json.loads(boundary_path.read_text(encoding="utf-8"))
    trials = {row["trial"]: row for row in json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]}
    measured = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    control = json.loads(control_path.read_text(encoding="utf-8"))
    rows = []
    payload = dict(
        completed=False,
        selected_trials=list(selected),
        frozen_targets=FROZEN_TARGETS,
        thermal_species_ratio=1.,
        thermal_species_ratio_fitted=False,
        mechanical_work="reduced_buoyancy_work",
        coarse_step_m=.005,
        refined_step_m=.0025,
        stage_inverse_tolerance=2e-5,
        accepted_endpoint_inverse_tolerance=1e-8,
        minimum_step_m=.000078125,
        successful_reduced_steps_before_growth=8,
        independent_balance_limit=1e-5,
        convergence_limit=.005,
        sensor_score_calculated=False,
        candidate_promoted=False,
        rows=rows,
        hashes=before,
    )

    def checkpoint():
        partial.write_text(
            json.dumps(payload, default=serial, indent=2, allow_nan=False)+"\n",
            encoding="utf-8",
        )

    checkpoint()
    for number in selected:
        trial = trials[number]
        entry = boundary["interfaces"][str(number)]
        stored_state = np.asarray(entry["state"], float)
        stored_beta = float(entry["thermal_width_ratio"])
        target_x = first_downstream_profile_target(trial, stored_state[5])
        if target_x != FROZEN_TARGETS[number]:
            raise RuntimeError(f"trial {number} target no longer matches the frozen table")
        jp, thermodynamics, source = actual_source_model(trial, measured)
        _trajectory, control_check = replay(control, trial)
        section = BuoyancyConstrainedEnthalpySection(
            jp, thermodynamics, thermal_width_ratio=stored_beta,
            quadrature_points=entry["quadrature_points"],
        )
        replayed = section.moments(stored_state)
        stored_state_replay_error = float(np.max(
            abs(replayed-np.asarray(entry["moments"]))/section.moment_scales(entry["moments"])
        ))
        projection = section.project_buoyancy(
            entry["moments"], stored_state, maximum_quadrature_points=4096,
        )
        projection_passed = bool(
            projection.success and max(projection.relative_residuals) <= 1e-8
            and max(projection.quadrature_residuals) <= 1e-5
            and not projection.at_width_bound
        )
        if not projection_passed:
            raise RuntimeError(f"trial {number} analytic-width boundary projection failed")
        initial_state = projection.state
        initial_beta = projection.thermal_width_ratio
        initial_parameters = encode_section(initial_state, initial_beta)
        driver = ReservoirShortSegment(
            jp, thermodynamics, thermal_species_ratio=1.,
            mechanical_work="reduced_buoyancy_work",
        )
        initial_values_16, initial_evaluation = independent_values(driver, initial_parameters)
        if not initial_evaluation["valid"]:
            raise RuntimeError(f"trial {number} initial local closure is invalid")

        def run(label, step):
            inverter = ThermalMomentFluxInverter(
                jp, thermodynamics, order=8, probes=257,
                quadrature_points=128, tolerance=1e-8,
            )
            marcher = FluxSpaceThermalMomentMarch(
                inverter, thermal_species_ratio=1.,
                mechanical_work="reduced_buoyancy_work",
                stage_inverse_tolerance=2e-5,
            )
            def progress(count, arc, primary):
                if count == 1 or count % 25 == 0:
                    print(
                        f"trial {number} {label}: step={count}, arc={arc:.6f} m, x={primary[6]:.6f} m",
                        flush=True,
                    )
            return marcher.march(
                initial_parameters, 2., step=step, target_x=target_x,
                maximum_steps=4000, target_tolerance=1e-8,
                minimum_step=.000078125, progress=progress,
            )

        row = dict(
            trial=number,
            target_downwind_m=target_x,
            source=source,
            control_replay=control_check,
            stored_state_replay_error=stored_state_replay_error,
            boundary_projection=asdict(projection),
            boundary_projection_passed=projection_passed,
            initial_values_order16=initial_values_16,
        )
        rows.append(row)
        try:
            coarse = run("coarse", .005)
            row["coarse"] = result_record(coarse, driver, initial_values_16)
            checkpoint()
            if coarse.reached_target and row["coarse"]["balance_passed"]:
                refined = run("refined", .0025)
                row["refined"] = result_record(refined, driver, initial_values_16)
                coarse_physical = physical_parameters(coarse.parameters[-1])
                refined_physical = physical_parameters(refined.parameters[-1])
                parameter_scale = np.maximum(np.maximum(abs(coarse_physical), abs(refined_physical)), 1.)
                flux_scale = np.maximum(
                    np.maximum(abs(row["coarse"]["terminal_values_order16"]),
                               abs(row["refined"]["terminal_values_order16"])),
                    [1e-12, 1e-12, 1., 1., 1., 1.],
                )
                parameter_difference = float(max(abs(coarse_physical-refined_physical)/parameter_scale))
                flux_difference = float(max(
                    abs(np.asarray(row["coarse"]["terminal_values_order16"])
                        -np.asarray(row["refined"]["terminal_values_order16"]))/flux_scale
                ))
                row["convergence"] = dict(
                    maximum_physical_parameter_scaled_difference=parameter_difference,
                    maximum_flux_scaled_difference=flux_difference,
                    limit=.005,
                    passed=bool(refined.reached_target and max(parameter_difference, flux_difference) <= .005),
                )
            row["all_gates_passed"] = bool(
                row["coarse"]["reached_target"]
                and row["coarse"]["balance_passed"]
                and row.get("refined", {}).get("reached_target", False)
                and row.get("refined", {}).get("balance_passed", False)
                and row.get("convergence", {}).get("passed", False)
            )
        except (ArithmeticError, FloatingPointError, RuntimeError, ValueError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
            row["all_gates_passed"] = False
        checkpoint()
        print(f"trial {number}: all_gates_passed={row['all_gates_passed']}", flush=True)

    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during the extension audit")
    payload["completed"] = True
    payload["finished_utc"] = datetime.now(timezone.utc).isoformat()
    payload["all_selected_gates_passed"] = all(row["all_gates_passed"] for row in rows)
    payload["full_six_trial_claim"] = bool(set(selected) == set(TRIALS))
    payload["hashes"] = after
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    partial.unlink()
    print(json.dumps({
        "completed": True,
        "selected_trials": selected,
        "all_selected_gates_passed": payload["all_selected_gates_passed"],
        "full_six_trial_claim": payload["full_six_trial_claim"],
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
