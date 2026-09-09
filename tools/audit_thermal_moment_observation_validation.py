"""Compare the direct-flux thermal-moment path at paired PRESLHY sensors."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_model_thermal_profile_observation import aggregate, model_summary
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
    ThermalMomentFluxTrajectory,
)
from degali.validation.thermal_profile_moments import project_model_profile


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"
TRIALS = (10, 23)
STATIONS = (1.78, 4.0)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def directionally_better(candidate: dict, control: dict) -> bool:
    """Apply the frozen componentwise amplitude/variance decision."""
    fields = ("centre_amplitude", "variance")
    candidate_thermal = candidate["thermal_deficit"]
    control_thermal = control["thermal_deficit"]
    complete = all(
        candidate_thermal[field]["complete"]
        and control_thermal[field]["complete"]
        for field in fields
    )
    if not complete:
        return False
    no_worse = all(
        candidate_thermal[field]["median_absolute_log_ratio"]
        <= control_thermal[field]["median_absolute_log_ratio"]
        for field in fields
    )
    better = any(
        candidate_thermal[field]["median_absolute_log_ratio"]
        < control_thermal[field]["median_absolute_log_ratio"]
        for field in fields
    )
    return bool(no_worse and better)


def sensor_refinement_difference(coarse_rows: list[dict], refined_rows: list[dict]) -> float:
    """Return the frozen maximum scaled difference over paired sensor values."""
    key = lambda row: (int(row["trial"]), float(row["station_m"]))
    coarse = {key(row): row for row in coarse_rows}
    refined = {key(row): row for row in refined_rows}
    if set(coarse) != set(refined) or len(coarse) != len(TRIALS)*len(STATIONS):
        raise ValueError("coarse and refined sensor populations do not match")
    worst = 0.0
    for item in sorted(coarse):
        left = coarse[item]["projection"]
        right = refined[item]["projection"]
        for field, floor in (("thermal_deficit_K", 1.0), ("hydrogen_vol_pct", .1)):
            a = np.asarray(left[field], float)
            b = np.asarray(right[field], float)
            if a.shape != (5,) or b.shape != (5,):
                raise ValueError("the frozen observation operator requires five sensors")
            scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), floor)
            worst = max(worst, float(np.max(np.abs(a-b)/scale)))
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=REF / "thermal_moment_observation_validation_2026-09-10.json",
    )
    parser.add_argument(
        "--positivity-domain", choices=("full_square", "radial_core"),
        default="full_square",
    )
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite final or partial observation evidence")

    boundary_path = REF / "buoyancy_constrained_enthalpy_width_interface_2026-09-05.json"
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    control_path = REF / "phase_ambient_consistency_field_complete_2026-09-05.json"
    observed_path = REF / "model_thermal_profile_observation_2026-09-08.json"
    dependencies = [
        boundary_path, reduced_path, measured_path, control_path, observed_path,
        ROOT / "docs/prereg-model-thermal-profile-observation.md",
        ROOT / "docs/prereg-thermal-moment-observation-validation.md",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/enthalpy_profile.py",
        ROOT / "src/degali/addons/energy_crosswind.py",
        ROOT / "src/degali/addons/reservoir_thermal.py",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        ROOT / "src/degali/addons/transverse_mixing.py",
        ROOT / "src/degali/validation/thermal_profile_moments.py",
        ROOT / "tools/audit_model_thermal_profile_observation.py",
        ROOT / "tools/audit_thermal_moment_flux_space.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    boundary = json.loads(boundary_path.read_text(encoding="utf-8"))
    trials = {
        row["trial"]: row
        for row in json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    }
    measured = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    control_evidence = json.loads(observed_path.read_text(encoding="utf-8"))
    control_field = json.loads(control_path.read_text(encoding="utf-8"))
    observed = control_evidence["observed"]
    control_rows = []
    candidate_rows = {"coarse": [], "refined": []}
    rows = []
    payload = {
        "completed": False,
        "selected_trials": list(TRIALS),
        "stations_m": list(STATIONS),
        "thermal_species_ratio": 1.0,
        "thermal_species_ratio_fitted": False,
        "mechanical_work": "reduced_buoyancy_work",
        "phase_interpolation": "linear",
        "positivity_domain": args.positivity_domain,
        "coarse_step_m": .005,
        "refined_step_m": .0025,
        "target_downwind_m": 4.0,
        "maximum_arc_length_m": 5.0,
        "stage_inverse_tolerance": 2e-5,
        "endpoint_inverse_tolerance": 1e-8,
        "minimum_step_m": .000078125,
        "balance_limit": 1e-5,
        "terminal_convergence_limit": .005,
        "sensor_convergence_limit": .005,
        "coefficient_fitted": False,
        "source_files_distributed": False,
        "candidate_promoted": False,
        "rows": rows,
        "control_rows": control_rows,
        "candidate_rows": candidate_rows,
        "hashes": before,
    }

    def checkpoint():
        partial.write_text(
            json.dumps(payload, default=serial, indent=2, allow_nan=False)+"\n",
            encoding="utf-8",
        )

    checkpoint()
    for number in TRIALS:
        trial = trials[number]
        jp, thermodynamics, source = actual_source_model(trial, measured)
        control_trajectory, control_check = replay(control_field, trial)
        if control_trajectory is None:
            raise RuntimeError(f"stored control trajectory is missing Trial {number}")
        for station in STATIONS:
            projection = project_model_profile(
                control_trajectory, station, trial["release_height_m"],
            )
            control_rows.append({
                "variant": "density_profile_control",
                "trial": number,
                "station_m": station,
                **model_summary(projection, observed[str(number)][str(station)]),
            })

        entry = boundary["interfaces"][str(number)]
        stored_state = np.asarray(entry["state"], float)
        section = BuoyancyConstrainedEnthalpySection(
            jp, thermodynamics,
            thermal_width_ratio=float(entry["thermal_width_ratio"]),
            quadrature_points=entry["quadrature_points"],
        )
        projection = section.project_buoyancy(
            entry["moments"], stored_state, maximum_quadrature_points=4096,
        )
        projection_passed = bool(
            projection.success and max(projection.relative_residuals) <= 1e-8
            and max(projection.quadrature_residuals) <= 1e-5
            and not projection.at_width_bound
        )
        row = {
            "trial": number,
            "source": source,
            "control_replay": control_check,
            "boundary_projection": asdict(projection),
            "boundary_projection_passed": projection_passed,
        }
        rows.append(row)
        if not projection_passed:
            row["failure"] = "analytic-width boundary projection failed"
            checkpoint()
            continue
        initial_parameters = encode_section(
            projection.state, projection.thermal_width_ratio,
        )
        driver = ReservoirShortSegment(
            jp, thermodynamics, thermal_species_ratio=1.,
            mechanical_work="reduced_buoyancy_work",
            positivity_domain=args.positivity_domain,
        )
        initial_values_16, initial_evaluation = independent_values(
            driver, initial_parameters,
        )
        if not initial_evaluation["valid"]:
            row["failure"] = "initial local closure is invalid"
            checkpoint()
            continue
        row["initial_values_order16"] = initial_values_16

        def run(label, step):
            inverter = ThermalMomentFluxInverter(
                jp, thermodynamics, order=8, probes=257,
                quadrature_points=128, tolerance=1e-8,
            )
            marcher = FluxSpaceThermalMomentMarch(
                inverter, thermal_species_ratio=1.,
                mechanical_work="reduced_buoyancy_work",
                stage_inverse_tolerance=2e-5,
                positivity_domain=args.positivity_domain,
            )

            def progress(count, arc, primary):
                if count == 1 or count % 50 == 0:
                    print(
                        f"trial {number} {label}: step={count}, "
                        f"arc={arc:.6f} m, x={primary[6]:.6f} m",
                        flush=True,
                    )

            return marcher.march(
                initial_parameters, 5., step=step, target_x=4.,
                maximum_steps=12000, target_tolerance=1e-8,
                minimum_step=.000078125, progress=progress,
            )

        try:
            for label, step in (("coarse", .005), ("refined", .0025)):
                result = run(label, step)
                row[label] = result_record(result, driver, initial_values_16)
                if result.reached_target and row[label]["balance_passed"]:
                    trajectory = ThermalMomentFluxTrajectory(
                        jp, thermodynamics, result.parameters,
                    )
                    for station in STATIONS:
                        receptor = project_model_profile(
                            trajectory, station, trial["release_height_m"],
                        )
                        candidate_rows[label].append({
                            "variant": f"thermal_moment_{label}",
                            "trial": number,
                            "station_m": station,
                            **model_summary(
                                receptor, observed[str(number)][str(station)],
                            ),
                        })
                checkpoint()
                if not result.reached_target or not row[label]["balance_passed"]:
                    break
            if "refined" in row:
                coarse_physical = physical_parameters(
                    np.asarray(row["coarse"]["parameters"])[-1]
                )
                refined_physical = physical_parameters(
                    np.asarray(row["refined"]["parameters"])[-1]
                )
                parameter_scale = np.maximum(
                    np.maximum(abs(coarse_physical), abs(refined_physical)), 1.,
                )
                parameter_difference = float(np.max(
                    abs(coarse_physical-refined_physical)/parameter_scale
                ))
                coarse_flux = np.asarray(row["coarse"]["terminal_values_order16"])
                refined_flux = np.asarray(row["refined"]["terminal_values_order16"])
                flux_scale = np.maximum(
                    np.maximum(abs(coarse_flux), abs(refined_flux)),
                    [1e-12, 1e-12, 1., 1., 1., 1.],
                )
                flux_difference = float(np.max(
                    abs(coarse_flux-refined_flux)/flux_scale
                ))
                row["terminal_convergence"] = {
                    "maximum_physical_parameter_scaled_difference": parameter_difference,
                    "maximum_flux_scaled_difference": flux_difference,
                    "limit": .005,
                    "passed": bool(max(parameter_difference, flux_difference) <= .005),
                }
        except (ArithmeticError, FloatingPointError, RuntimeError, ValueError,
                np.linalg.LinAlgError) as error:
            row["failure"] = str(error)
        checkpoint()
        print(
            f"trial {number}: coarse={row.get('coarse', {}).get('reached_target')}, "
            f"refined={row.get('refined', {}).get('reached_target')}",
            flush=True,
        )

    control_aggregate = aggregate(control_rows)
    candidate_aggregates = {
        label: aggregate(candidate_rows[label])
        for label in ("coarse", "refined")
        if len(candidate_rows[label]) == len(TRIALS)*len(STATIONS)
    }
    payload["aggregates"] = {
        "control": control_aggregate,
        **candidate_aggregates,
    }
    payload["directionally_better_than_control"] = {
        label: directionally_better(value, control_aggregate)
        for label, value in candidate_aggregates.items()
    }
    if all(
        len(candidate_rows[label]) == len(TRIALS)*len(STATIONS)
        for label in ("coarse", "refined")
    ):
        difference = sensor_refinement_difference(
            candidate_rows["coarse"], candidate_rows["refined"],
        )
    else:
        # Keep an interrupted/failed audit serializable with ``allow_nan=False``.
        # A missing population is a failed gate, not an infinite measurement.
        difference = None
    payload["sensor_refinement"] = {
        "maximum_scaled_difference": difference,
        "limit": .005,
        "passed": bool(difference is not None and difference <= .005),
    }
    numerical = all(
        row.get("coarse", {}).get("reached_target", False)
        and row.get("coarse", {}).get("balance_passed", False)
        and row.get("refined", {}).get("reached_target", False)
        and row.get("refined", {}).get("balance_passed", False)
        and row.get("terminal_convergence", {}).get("passed", False)
        for row in rows
    )
    observation = all(
        payload["directionally_better_than_control"].get(label, False)
        for label in ("coarse", "refined")
    )
    payload["all_numerical_gates_passed"] = bool(
        numerical and payload["sensor_refinement"]["passed"]
    )
    payload["observation_gate_passed"] = bool(observation)
    payload["research_candidate_may_advance"] = bool(
        payload["all_numerical_gates_passed"] and observation
    )
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during the observation audit")
    payload["hashes"] = after
    payload["completed"] = True
    payload["finished_utc"] = datetime.now(timezone.utc).isoformat()
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    partial.unlink()
    print(json.dumps({
        "completed": True,
        "all_numerical_gates_passed": payload["all_numerical_gates_passed"],
        "observation_gate_passed": payload["observation_gate_passed"],
        "research_candidate_may_advance": payload["research_candidate_may_advance"],
        "saved": str(args.output),
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
