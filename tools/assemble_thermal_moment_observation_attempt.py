"""Recover sensor diagnostics from an immutable stopped observation march."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_model_thermal_profile_observation import aggregate, model_summary
from audit_reservoir_thermal_segments import actual_source_model, serial
from audit_thermal_moment_observation_validation import directionally_better
from degali.addons.thermal_moment_flux_march import ThermalMomentFluxTrajectory
from degali.validation.thermal_profile_moments import project_model_profile


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"
TRIALS = (10, 23)
STATIONS = (1.78, 4.0)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--attempt", type=Path,
        default=REF / "thermal_moment_observation_validation_first_attempt_2026-09-10.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REF / "thermal_moment_observation_first_attempt_analysis_2026-09-10.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite recovered observation evidence")
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    observed_path = REF / "model_thermal_profile_observation_2026-09-08.json"
    dependencies = [
        args.attempt, reduced_path, measured_path, observed_path,
        ROOT / "docs/prereg-thermal-moment-observation-validation.md",
        ROOT / "src/degali/addons/buoyancy_profile.py",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        ROOT / "src/degali/validation/thermal_profile_moments.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    attempt = json.loads(args.attempt.read_text(encoding="utf-8"))
    trials = {
        row["trial"]: row
        for row in json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    }
    measured = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    observed = json.loads(observed_path.read_text(encoding="utf-8"))["observed"]
    candidate_rows = []
    unavailable = []
    for row in attempt["rows"]:
        number = int(row["trial"])
        trial = trials[number]
        jp, thermodynamics, _source = actual_source_model(trial, measured)
        coarse = row.get("coarse")
        if coarse is None:
            unavailable.extend({"trial": number, "station_m": x, "reason": "no coarse march"}
                               for x in STATIONS)
            continue
        trajectory = ThermalMomentFluxTrajectory(
            jp, thermodynamics, np.asarray(coarse["parameters"], float),
        )
        for station in STATIONS:
            if trajectory.parameters_at(station) is None:
                unavailable.append({
                    "trial": number,
                    "station_m": station,
                    "reason": coarse["stop_reason"],
                    "last_downwind_m": coarse["final_downwind_m"],
                })
                continue
            projection = project_model_profile(
                trajectory, station, trial["release_height_m"],
            )
            candidate_rows.append({
                "variant": "thermal_moment_coarse_stopped_attempt",
                "trial": number,
                "station_m": station,
                **model_summary(projection, observed[str(number)][str(station)]),
            })

    control_rows = attempt["control_rows"]
    candidate_keys = {
        (int(row["trial"]), float(row["station_m"]))
        for row in candidate_rows
    }
    matched_control_rows = [
        row for row in control_rows
        if (int(row["trial"]), float(row["station_m"])) in candidate_keys
    ]
    control_aggregate = aggregate(control_rows)
    matched_control_aggregate = aggregate(matched_control_rows)
    candidate_aggregate = aggregate(candidate_rows)
    full_population = len(candidate_rows) == len(TRIALS)*len(STATIONS)
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_attempt": str(args.attempt.relative_to(ROOT)),
        "source_attempt_completed": bool(attempt.get("completed", False)),
        "selected_trials": list(TRIALS),
        "stations_m": list(STATIONS),
        "control_rows": control_rows,
        "matched_control_rows": matched_control_rows,
        "candidate_rows": candidate_rows,
        "unavailable_profiles": unavailable,
        "complete_observation_population": full_population,
        "control_aggregate": control_aggregate,
        "matched_control_aggregate": matched_control_aggregate,
        "candidate_partial_aggregate": candidate_aggregate,
        "partial_directionally_better_on_matched_profiles": bool(
            directionally_better(candidate_aggregate, matched_control_aggregate)
        ),
        "directionally_better_than_control": bool(
            full_population
            and directionally_better(candidate_aggregate, control_aggregate)
        ),
        "observation_gate_passed": False,
        "research_candidate_may_advance": False,
        "interpretation": (
            "Partial sensor values are diagnostics only. Trial 10 stopped before "
            "4 m, so the frozen four-profile observation decision is unavailable."
        ),
        "hashes": before,
    }
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed while recovering the stopped attempt")
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    print(json.dumps({
        "candidate_profiles": len(candidate_rows),
        "unavailable_profiles": unavailable,
        "observation_gate_passed": False,
        "saved": str(args.output),
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
