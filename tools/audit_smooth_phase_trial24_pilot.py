"""Audit C1 phase lookup continuation beyond the retained trial 24 corner."""

from dataclasses import asdict
from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_reservoir_thermal_segments import actual_source_model, serial
from audit_thermal_moment_flux_space import independent_values, physical_parameters, result_record
from run_preslhy_ambient_profile_audit import replay
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.reservoir_thermal import ReservoirShortSegment, encode_section
from degali.addons.thermal_moment_flux_march import (
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
)


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite smooth trial 24 pilot evidence")
    boundary_path = REF / "buoyancy_constrained_enthalpy_width_interface_2026-09-05.json"
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    control_path = REF / "phase_ambient_consistency_field_complete_2026-09-05.json"
    screen_path = REF / "smooth_phase_hermite_lookup_2026-09-10.json"
    dependencies = [
        boundary_path, reduced_path, measured_path, control_path, screen_path,
        ROOT / "docs/prereg-smooth-phase-hermite-lookup.md",
        ROOT / "docs/prereg-smooth-phase-trial24-pilot.md",
        ROOT / "docs/prereg-analytic-enthalpy-width-split.md",
        ROOT / "src/degali/addons/smooth_phase_lookup.py",
        ROOT / "src/degali/addons/energy_crosswind.py",
        ROOT / "src/degali/addons/enthalpy_profile.py",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        ROOT / "src/degali/addons/reservoir_thermal.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    boundary = json.loads(boundary_path.read_text(encoding="utf-8"))
    trial = next(
        row for row in json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
        if row["trial"] == 24
    )
    measured = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    control = json.loads(control_path.read_text(encoding="utf-8"))
    jp, thermodynamics, source = actual_source_model(trial, measured)
    _, control_check = replay(control, trial)
    entry = boundary["interfaces"]["24"]
    section = BuoyancyConstrainedEnthalpySection(
        jp, thermodynamics, thermal_width_ratio=entry["thermal_width_ratio"],
        quadrature_points=entry["quadrature_points"],
        phase_interpolation="c1_hermite",
    )
    projection = section.project_buoyancy(
        entry["moments"], np.asarray(entry["state"], float),
        maximum_quadrature_points=4096,
    )
    projection_passed = bool(
        projection.success and max(projection.relative_residuals) <= 1e-8
        and max(projection.quadrature_residuals) <= 1e-5
        and not projection.at_width_bound
    )
    payload = {
        "completed": False, "trial": 24, "target_downwind_m": .43,
        "phase_interpolation": "c1_hermite", "default_changed": False,
        "source": source, "control_replay": control_check,
        "projection": asdict(projection), "projection_passed": projection_passed,
        "sensor_score_calculated": False, "candidate_promoted": False,
        "hashes": before,
    }
    partial.write_text(json.dumps(payload, default=serial, indent=2)+"\n", encoding="utf-8")
    if not projection_passed:
        raise RuntimeError("smooth trial 24 boundary projection failed")
    initial = encode_section(projection.state, projection.thermal_width_ratio)
    driver = ReservoirShortSegment(
        jp, thermodynamics, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", phase_interpolation="c1_hermite",
    )
    initial_values, initial_evaluation = independent_values(driver, initial)
    if not initial_evaluation["valid"]:
        raise RuntimeError("smooth projected boundary has invalid local closure")

    def run(label, step):
        inverter = ThermalMomentFluxInverter(
            jp, thermodynamics, order=8, probes=257, quadrature_points=128,
            tolerance=1e-8, phase_interpolation="c1_hermite",
        )
        marcher = FluxSpaceThermalMomentMarch(
            inverter, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work",
            stage_inverse_tolerance=2e-5,
        )
        def progress(count, arc, primary):
            print(f"trial 24 smooth {label}: step={count}, arc={arc:.8f}, x={primary[6]:.8f}", flush=True)
        return marcher.march(
            initial, .1, step=step, target_x=.43, maximum_steps=4000,
            target_tolerance=1e-8, minimum_step=.000078125, progress=progress,
        )

    coarse = run("coarse", .005)
    refined = run("refined", .0025)
    payload["coarse"] = result_record(coarse, driver, initial_values)
    payload["refined"] = result_record(refined, driver, initial_values)
    coarse_physical = physical_parameters(coarse.parameters[-1])
    refined_physical = physical_parameters(refined.parameters[-1])
    parameter_scale = np.maximum(np.maximum(abs(coarse_physical), abs(refined_physical)), 1.)
    flux_scale = np.maximum(
        np.maximum(abs(payload["coarse"]["terminal_values_order16"]),
                   abs(payload["refined"]["terminal_values_order16"])),
        [1e-12, 1e-12, 1., 1., 1., 1.],
    )
    parameter_difference = float(max(abs(coarse_physical-refined_physical)/parameter_scale))
    flux_difference = float(max(
        abs(np.asarray(payload["coarse"]["terminal_values_order16"])
            -np.asarray(payload["refined"]["terminal_values_order16"]))/flux_scale
    ))
    payload["convergence"] = {
        "maximum_physical_parameter_scaled_difference": parameter_difference,
        "maximum_flux_scaled_difference": flux_difference,
        "limit": .005,
        "passed": bool(max(parameter_difference, flux_difference) <= .005),
    }
    payload["all_gates_passed"] = bool(
        coarse.reached_target and refined.reached_target
        and payload["coarse"]["balance_passed"]
        and payload["refined"]["balance_passed"]
        and payload["convergence"]["passed"]
    )
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("a smooth trial 24 pilot input changed during execution")
    payload["hashes"] = after
    payload["completed"] = True
    payload["finished_utc"] = datetime.now(timezone.utc).isoformat()
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
        stream.write("\n")
    partial.unlink()
    print(json.dumps({
        "completed": True, "all_gates_passed": payload["all_gates_passed"],
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
