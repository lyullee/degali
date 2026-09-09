"""Apply the paired five-sensor operator to stored Trial-10 wind fields."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from degali.addons.observed_wind_profile import with_observed_wind
from degali.addons.yawed_crosswind import YawedCrosswind, YawedTrajectory
from degali.validation.thermal_profile_moments import (
    project_model_profile,
    read_paired_profiles,
)
from run_preslhy_ambient_profile_audit import replay
from audit_model_thermal_profile_observation import (
    model_summary,
    observed_profile_summary,
)


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"
FIELDS = {
    "yaw_direction_step002": REF / "yawed_trial10_step002_2026-09-06",
    "magnitude_step002": REF / "observed_wind_magnitude_step002_2026-09-06",
    "measured_vector_step002": REF / "observed_wind_vector_step002_2026-09-06",
    "measured_vector_step001": REF / "observed_wind_vector_step001_2026-09-06",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class YawedPercentTrajectory:
    """Expose a yawed 0--1 mole-fraction trajectory as the common vol-% API."""

    def __init__(self, trajectory: YawedTrajectory):
        self.trajectory = trajectory
        self.model = SimpleNamespace(
            thermodynamics=trajectory.model.base.thermodynamics
        )

    def state_at(self, station):
        try:
            return self.trajectory.section(station, 0.0)[0]
        except ValueError:
            return None

    def temperature_at(self, station, lateral, height):
        return self.trajectory.temperature_at(station, lateral, height)

    def concentration_at(self, station, lateral, height):
        return 100.0 * self.trajectory.concentration_at(station, lateral, height)


def group_errors(rows: list[dict]) -> dict:
    result = {}
    for scalar in ("thermal_deficit", "hydrogen"):
        result[scalar] = {}
        for field in ("centre_amplitude", "zeroth_moment", "variance"):
            values = [
                row["comparison"][scalar][field]["absolute_log_ratio"]
                for row in rows
            ]
            if any(value is None or not math.isfinite(value) for value in values):
                raise ValueError("all amplitude and moment ratios must be positive")
            result[scalar][field] = {
                "profiles": len(values),
                "median_absolute_log_ratio": float(np.median(values)),
                "maximum_absolute_log_ratio": float(np.max(values)),
                "within_observed_iqr": sum(
                    row["comparison"][scalar][field]["within_observed_iqr"]
                    for row in rows
                ),
            }
    return result


def displacement_gate(control: dict, candidate: dict) -> dict:
    amplitude_fields = (
        ("thermal_deficit", "centre_amplitude"),
        ("thermal_deficit", "zeroth_moment"),
        ("hydrogen", "centre_amplitude"),
        ("hydrogen", "zeroth_moment"),
    )
    variance_fields = (
        ("thermal_deficit", "variance"),
        ("hydrogen", "variance"),
    )
    comparisons = {}
    for scalar, field in amplitude_fields + variance_fields:
        old = control[scalar][field]["median_absolute_log_ratio"]
        new = candidate[scalar][field]["median_absolute_log_ratio"]
        comparisons[f"{scalar}_{field}"] = {
            "control": old,
            "candidate": new,
            "not_worse": bool(new <= old),
            "strictly_better": bool(new < old),
        }
    amplitude_pass = all(
        comparisons[f"{scalar}_{field}"]["not_worse"]
        for scalar, field in amplitude_fields
    )
    variance_pass = all(
        comparisons[f"{scalar}_{field}"]["not_worse"]
        for scalar, field in variance_fields
    )
    strict = any(entry["strictly_better"] for entry in comparisons.values())
    return {
        "passed": bool(amplitude_pass and variance_pass and strict),
        "amplitude_passed": bool(amplitude_pass),
        "variance_passed": bool(variance_pass),
        "at_least_one_strict_improvement": bool(strict),
        "comparisons": comparisons,
    }


def step_convergence(coarse_rows: list[dict], fine_rows: list[dict]) -> dict:
    coarse = {row["station_m"]: row for row in coarse_rows}
    fine = {row["station_m"]: row for row in fine_rows}
    if set(coarse) != set(fine) or not coarse:
        raise ValueError("step comparison requires matching stations")
    profile_errors = []
    moment_errors = []
    details = []
    for station in sorted(coarse):
        station_profile = []
        station_moment = []
        for key in ("thermal_deficit_K", "hydrogen_vol_pct"):
            a = np.asarray(coarse[station]["projection"][key], dtype=float)
            b = np.asarray(fine[station]["projection"][key], dtype=float)
            scale = np.maximum.reduce([np.abs(a), np.abs(b), np.full_like(a, 1e-6)])
            station_profile.extend((np.abs(a - b) / scale).tolist())
        for key in ("thermal_moment", "hydrogen_moment"):
            for field in ("zeroth", "centre", "variance"):
                a = float(coarse[station]["projection"][key][field])
                b = float(fine[station]["projection"][key][field])
                scale = max(abs(a), abs(b), 1e-6)
                station_moment.append(abs(a - b) / scale)
        profile_errors.extend(station_profile)
        moment_errors.extend(station_moment)
        details.append({
            "station_m": station,
            "maximum_profile_relative_difference": max(station_profile),
            "maximum_moment_relative_difference": max(station_moment),
        })
    limit = 0.005
    return {
        "maximum_profile_relative_difference": max(profile_errors),
        "maximum_moment_relative_difference": max(moment_errors),
        "limit": limit,
        "passed": bool(max(profile_errors) <= limit and max(moment_errors) <= limit),
        "stations": details,
    }


def reconstruct_yawed(original, directory: Path):
    inputs = json.loads((directory / "inputs.json").read_text(encoding="utf-8"))
    complete = json.loads((directory / "complete.json").read_text(encoding="utf-8"))
    field = json.loads((directory / "field.json").read_text(encoding="utf-8"))
    if not complete.get("completed") or not complete.get("accepted"):
        raise ValueError(f"stored wind field is incomplete: {directory.name}")
    if field["maximum_relative_balance_residual"] > 1.0e-5:
        raise ValueError(f"stored wind field fails balance: {directory.name}")
    if inputs.get("mode") is None:
        base = original.model
        angle = math.radians(inputs["wind_angle_deg"])
    else:
        base = with_observed_wind(
            original.model,
            speed=inputs["wind_speed"],
            reference_height=inputs["reference_height"],
        )
        angle = math.radians(inputs["relative_wind_angle_deg"])
    model = YawedCrosswind(base, angle)
    states = np.asarray(field["states"], dtype=float)
    stored_fluxes = np.asarray(field["fluxes"], dtype=float)
    stored_sources = np.asarray(field["sources"], dtype=float)
    indices = np.unique(np.linspace(0, len(states) - 1, 9, dtype=int))
    flux_error = source_error = 0.0
    for index in indices:
        flux = model.fluxes(states[index])
        source = model.sources(states[index])
        flux_error = max(
            flux_error,
            float(np.max(np.abs(flux - stored_fluxes[index])
                         / np.maximum(np.abs(stored_fluxes[index]), 1.0))),
        )
        source_error = max(
            source_error,
            float(np.max(np.abs(source - stored_sources[index])
                         / np.maximum(np.abs(stored_sources[index]), 1.0))),
        )
    if max(flux_error, source_error) > 1.0e-8:
        raise ValueError(f"stored wind field does not replay: {directory.name}")
    check = {
        "completed": True,
        "handoff_accepted": True,
        "states": len(states),
        "step_m": inputs["step_m"],
        "wind_speed_m_s": inputs.get("wind_speed"),
        "wind_reference_height_m": inputs.get("reference_height"),
        "relative_wind_angle_deg": math.degrees(angle),
        "maximum_relative_balance_residual": field["maximum_relative_balance_residual"],
        "sampled_flux_max_scaled_difference": flux_error,
        "sampled_source_max_scaled_difference": source_error,
        "sampled_sections": len(indices),
    }
    return YawedPercentTrajectory(YawedTrajectory(model, states)), check, inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REF / "wind_vector_profile_observation_2026-09-08.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite wind-vector observation evidence")

    control_path = REF / "phase_ambient_consistency_field_complete_2026-09-05.json"
    reduced_path = REF / "e35_reduced.json"
    report = ROOT / "tmp/pdfs/PRESLHY_D3.6_Summary_of_Rainout_Experiments_V1.20.pdf"
    workbook = REF / "raw/trial_10_13-09-2019alldata.xlsx"
    previous_audit = REF / "model_thermal_profile_observation_2026-09-08.json"
    dependencies = [
        control_path,
        reduced_path,
        report,
        workbook,
        previous_audit,
        ROOT / "docs/prereg-wind-vector-profile-observation.md",
        ROOT / "docs/prereg-model-thermal-profile-observation.md",
        ROOT / "src/degali/addons/yawed_crosswind.py",
        ROOT / "src/degali/addons/observed_wind_profile.py",
        ROOT / "src/degali/validation/thermal_profile_moments.py",
        ROOT / "tests/test_wind_vector_profile_observation.py",
        ROOT / "tools/audit_model_thermal_profile_observation.py",
        ROOT / "tools/run_preslhy_ambient_profile_audit.py",
        Path(__file__).resolve(),
    ]
    for directory in FIELDS.values():
        dependencies.extend(
            directory / name for name in ("inputs.json", "complete.json", "field.json")
        )
    hashes_before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    reduced = json.loads(reduced_path.read_text(encoding="utf-8"))
    trial = next(trial for trial in reduced["trials"] if trial["trial"] == 10)
    control_field = json.loads(control_path.read_text(encoding="utf-8"))
    original, control_check = replay(control_field, trial)
    paired = read_paired_profiles(workbook, report, flexlogger_rows=(30, 78))
    observed = {}
    for station in (1.78, 4.0):
        profiles = paired["profiles"][str(station)]
        observed[str(station)] = {
            "thermal_deficit": observed_profile_summary(
                paired["coordinate"], profiles["thermal_deficit"], centre_gate=5.0
            ),
            "hydrogen": observed_profile_summary(
                paired["coordinate"], profiles["hydrogen"], centre_gate=1.0
            ),
        }

    trajectories = {"planar_control": original}
    integrity = {"planar_control": control_check}
    inputs = {}
    for name, directory in FIELDS.items():
        trajectories[name], integrity[name], inputs[name] = reconstruct_yawed(
            original, directory
        )
    vector_inputs = [inputs[name] for name in (
        "measured_vector_step002", "measured_vector_step001"
    )]
    for key in ("mode", "wind_speed", "reference_height", "relative_wind_angle_deg"):
        if vector_inputs[0][key] != vector_inputs[1][key]:
            raise ValueError("measured-vector step inputs differ")

    rows = {}
    summaries = {}
    for name, trajectory in trajectories.items():
        rows[name] = []
        for station in (1.78, 4.0):
            projection = project_model_profile(
                trajectory, station, trial["release_height_m"]
            )
            rows[name].append({
                "variant": name,
                "trial": 10,
                "station_m": station,
                **model_summary(projection, observed[str(station)]),
            })
        summaries[name] = group_errors(rows[name])
    gates = {
        name: displacement_gate(summaries["planar_control"], summaries[name])
        for name in FIELDS
    }
    convergence = step_convergence(
        rows["measured_vector_step002"], rows["measured_vector_step001"]
    )
    hashes_after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if hashes_before != hashes_after:
        raise RuntimeError("an input changed while calculating")
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "trial": 10,
        "source_dataset": "PRESLHY E3.5",
        "source_doi": "10.35097/1481",
        "source_files_distributed": False,
        "trajectory_reintegrated": False,
        "wind_fitted_to_temperature": False,
        "trial23_excluded_faulty_wind": True,
        "hashes": hashes_after,
        "integrity": integrity,
        "observed": observed,
        "rows": rows,
        "summaries": summaries,
        "displacement_gates": gates,
        "measured_vector_step_convergence": convergence,
        "candidate_promoted": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
    print(json.dumps({
        "summaries": summaries,
        "displacement_gates": {name: gate["passed"] for name, gate in gates.items()},
        "step_convergence": convergence,
        "saved": str(args.output),
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
