"""Hash and summarize immutable subset results for the six-trial claim."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


TRIALS = (11, 12, 22, 23, 24, 25)
FROZEN_SETTINGS = {
    "thermal_species_ratio": 1.,
    "thermal_species_ratio_fitted": False,
    "mechanical_work": "reduced_buoyancy_work",
    "coarse_step_m": .005,
    "refined_step_m": .0025,
    "stage_inverse_tolerance": 2e-5,
    "accepted_endpoint_inverse_tolerance": 1e-8,
    "minimum_step_m": .000078125,
    "successful_reduced_steps_before_growth": 8,
    "independent_balance_limit": 1e-5,
    "convergence_limit": .005,
    "sensor_score_calculated": False,
    "candidate_promoted": False,
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite combined extension evidence")
    rows = {}
    sources = []
    for path in args.inputs:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not payload.get("completed") or not payload.get("all_selected_gates_passed"):
            raise RuntimeError(f"subset evidence is incomplete or failed: {path}")
        mismatched = {
            key: payload.get(key) for key, value in FROZEN_SETTINGS.items()
            if payload.get(key) != value
        }
        if mismatched:
            raise RuntimeError(f"subset evidence changed frozen settings: {mismatched}")
        if set(payload.get("selected_trials", ())) != {
                int(row["trial"]) for row in payload.get("rows", ())}:
            raise RuntimeError(f"subset trial declaration does not match its rows: {path}")
        sources.append({"path": str(path), "sha256": digest(path)})
        for row in payload["rows"]:
            trial = int(row["trial"])
            if trial in rows:
                raise RuntimeError(f"duplicate trial {trial} in subset evidence")
            rows[trial] = row
    if tuple(sorted(rows)) != TRIALS:
        raise RuntimeError(f"combined evidence does not contain exactly {TRIALS}")
    summary = []
    for trial in TRIALS:
        row = rows[trial]
        if not row.get("all_gates_passed"):
            raise RuntimeError(f"trial {trial} did not pass every registered gate")
        summary.append({
            "trial": trial,
            "stored_state_replay_error": row["stored_state_replay_error"],
            "boundary_projection_passed": row["boundary_projection_passed"],
            "coarse_reached_target": row["coarse"]["reached_target"],
            "refined_reached_target": row["refined"]["reached_target"],
            "coarse_balance_passed": row["coarse"]["balance_passed"],
            "refined_balance_passed": row["refined"]["balance_passed"],
            "coarse_accepted_steps": row["coarse"]["accepted_steps"],
            "refined_accepted_steps": row["refined"]["accepted_steps"],
            "coarse_rejected_steps": row["coarse"]["rejected_steps"],
            "refined_rejected_steps": row["refined"]["rejected_steps"],
            "coarse_maximum_inverse_residual": row["coarse"]["maximum_inverse_residual"],
            "refined_maximum_inverse_residual": row["refined"]["maximum_inverse_residual"],
            "minimum_sampled_diffusivity": min(
                row["coarse"]["minimum_sampled_diffusivity"],
                row["refined"]["minimum_sampled_diffusivity"],
            ),
            "maximum_physical_parameter_scaled_difference": row["convergence"][
                "maximum_physical_parameter_scaled_difference"
            ],
            "maximum_flux_scaled_difference": row["convergence"][
                "maximum_flux_scaled_difference"
            ],
            "all_gates_passed": True,
        })
    result = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "selected_trials": list(TRIALS),
        "full_six_trial_claim": True,
        "all_gates_passed": True,
        "source_evidence": sources,
        "summary": summary,
    }
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
