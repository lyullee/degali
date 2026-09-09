"""Audit guarded downstream reach of the conservative thermal-moment state."""

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
from degali.addons.reservoir_thermal import (
    ReservoirShortSegment,
    decode_section,
    encode_section,
)
from degali.addons.thermal_moment_march import GuardedThermalMomentMarch
from degali.addons.transverse_mixing import ConservativeTransverseMixing


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def result_record(result) -> dict:
    record = asdict(result)
    final_state, final_beta = decode_section(result.states[-1, :8])
    record.update(
        accepted_steps=len(result.arc_length) - 1,
        final_arc_length_m=float(result.arc_length[-1]),
        final_downwind_m=float(final_state[5]),
        final_elevation_m=float(final_state[6]),
        final_thermal_width_ratio=final_beta,
    )
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REF / "thermal_moment_downstream_reach_2026-09-08.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite thermal-moment reach evidence")

    boundary_path = REF / "buoyancy_constrained_enthalpy_width_interface_2026-09-05.json"
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    control_path = REF / "phase_ambient_consistency_field_complete_2026-09-05.json"
    dependencies = [
        boundary_path,
        reduced_path,
        measured_path,
        control_path,
        ROOT / "docs/prereg-thermal-moment-downstream-reach.md",
        ROOT / "src/degali/addons/thermal_moment_march.py",
        ROOT / "src/degali/addons/reservoir_thermal.py",
        ROOT / "src/degali/addons/transverse_mixing.py",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "tests/test_thermal_moment_march.py",
        Path(__file__).resolve(),
    ]
    hashes_before = {
        str(path.relative_to(ROOT)): digest(path) for path in dependencies
    }
    boundary = json.loads(boundary_path.read_text(encoding="utf-8"))
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    measured = json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    control = json.loads(control_path.read_text(encoding="utf-8"))
    trial = next(row for row in trials if row["trial"] == 10)
    measured_rows = {row["trial"]: row for row in measured}
    jp, thermodynamics, source = actual_source_model(trial, measured_rows)

    _trajectory, control_check = replay(control, trial)
    entry = boundary["interfaces"]["10"]
    initial_state = np.asarray(entry["state"], dtype=float)
    initial_beta = float(entry["thermal_width_ratio"])
    replay_section = BuoyancyConstrainedEnthalpySection(
        jp,
        thermodynamics,
        thermal_width_ratio=initial_beta,
        quadrature_points=entry["quadrature_points"],
    )
    represented_boundary = replay_section.moments(initial_state)
    boundary_error = float(
        np.max(
            np.abs(represented_boundary - np.asarray(entry["moments"]))
            / replay_section.moment_scales(entry["moments"])
        )
    )
    if boundary_error > 1.0e-9:
        raise RuntimeError("frozen six-moment boundary replay failed")

    driver = ReservoirShortSegment(
        jp,
        thermodynamics,
        thermal_species_ratio=1.0,
        mechanical_work="reduced_buoyancy_work",
    )
    initial_parameters = encode_section(initial_state, initial_beta)
    initial_evaluation = driver.evaluate(initial_parameters, order=16, probes=1025)
    initial_values = driver.flux_values(initial_evaluation["mixing"])
    scales = np.maximum(
        np.abs(initial_values), np.asarray([1e-12, 1e-12, 1.0, 1.0, 1.0, 1.0])
    )
    target_x = 1.78
    target = lambda parameters: decode_section(parameters)[0][5] >= target_x
    marcher = GuardedThermalMomentMarch(driver.evaluate)
    coarse = marcher.march(
        initial_parameters,
        1.25,
        initial_step=0.002,
        maximum_step=0.005,
        minimum_step=1.0e-6,
        tolerance=2.0e-5,
        maximum_steps=4000,
        target=target,
    )
    coarse_record = result_record(coarse)
    final_evaluation = driver.evaluate(coarse.states[-1, :8], order=16, probes=1025)
    final_values = driver.flux_values(final_evaluation["mixing"])
    cumulative = coarse.states[-1, 8:]
    balance = (final_values - initial_values - cumulative) / scales
    coarse_record.update(
        final_values=final_values,
        cumulative_sources=cumulative,
        scaled_balance=balance,
        maximum_scaled_balance=float(np.max(np.abs(balance))),
        balance_passed=bool(np.max(np.abs(balance)) <= 5.0e-4),
    )

    refined_record = None
    convergence = None
    if coarse.reached_target:
        refined = marcher.march(
            initial_parameters,
            1.25,
            initial_step=0.001,
            maximum_step=0.0025,
            minimum_step=1.0e-6,
            tolerance=2.0e-5,
            maximum_steps=4000,
            target=target,
        )
        refined_record = result_record(refined)
        refined_evaluation = driver.evaluate(
            refined.states[-1, :8], order=16, probes=1025
        )
        refined_values = driver.flux_values(refined_evaluation["mixing"])
        parameter_scale = np.maximum(
            np.maximum(np.abs(coarse.states[-1, :8]), np.abs(refined.states[-1, :8])),
            1.0,
        )
        value_scale = np.maximum(
            np.maximum(np.abs(final_values), np.abs(refined_values)),
            np.asarray([1e-12, 1e-12, 1.0, 1.0, 1.0, 1.0]),
        )
        parameter_difference = float(
            np.max(np.abs(coarse.states[-1, :8] - refined.states[-1, :8]) / parameter_scale)
        )
        value_difference = float(
            np.max(np.abs(final_values - refined_values) / value_scale)
        )
        convergence = {
            "maximum_parameter_relative_difference": parameter_difference,
            "maximum_value_relative_difference": value_difference,
            "limit": 0.005,
            "passed": bool(
                refined.reached_target
                and max(parameter_difference, value_difference) <= 0.005
            ),
        }

    hashes_after = {
        str(path.relative_to(ROOT)): digest(path) for path in dependencies
    }
    if hashes_before != hashes_after:
        raise RuntimeError("an input changed during the reach audit")
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "trial": 10,
        "source": source,
        "target_downwind_m": target_x,
        "thermal_species_ratio": 1.0,
        "thermal_species_ratio_fitted": False,
        "mechanical_work": "reduced_buoyancy_work",
        "yaw_included": False,
        "ground_interaction_included": False,
        "source_files_distributed": False,
        "hashes": hashes_after,
        "control_replay": control_check,
        "boundary_moment_replay_error": boundary_error,
        "coarse": coarse_record,
        "refined": refined_record,
        "convergence": convergence,
        "sensor_score_calculated": False,
        "candidate_promoted": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    print(json.dumps({
        "reached_target": coarse.reached_target,
        "stop_reason": coarse.stop_reason,
        "final_arc_length_m": coarse_record["final_arc_length_m"],
        "final_downwind_m": coarse_record["final_downwind_m"],
        "accepted_steps": coarse_record["accepted_steps"],
        "rejected_steps": coarse.rejected_steps,
        "rhs_calls": coarse.rhs_calls,
        "maximum_scaled_balance": coarse_record["maximum_scaled_balance"],
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
