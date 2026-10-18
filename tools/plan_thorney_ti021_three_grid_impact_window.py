"""Freeze a concentration-blind TI021 three-grid obstacle-impact window."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

try:
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-three-grid-impact-window-plan.v3"


def build_plan(
    checkpoint_journal: Path,
    *,
    failed_fine_attempt: Path | None = None,
    interrupted_fine_attempt: Path | None = None,
) -> dict[str, object]:
    source = checkpoint_journal.resolve()
    journal = json.loads(source.read_text(encoding="utf-8"))
    if journal.get("schema") != "degali.thorney-ti021-inner-checkpoint.v1":
        raise ValueError("unsupported TI021 inner checkpoint schema")
    entries = [
        entry for entry in journal["segments"]
        if float(entry["stop_time_s"]) <= 20.0 + 1.0e-12
    ]
    if not entries or float(entries[-1]["stop_time_s"]) < 20.0 - 1.0e-12:
        raise ValueError("at least 20 committed model seconds are required")
    measured_dx2p5_wall_s = sum(float(entry["wall_clock_s"]) for entry in entries)
    source_configuration = journal.get("configuration", {})
    source_cell_size_m = float(source_configuration.get("cell_size_m", 2.5))
    source_time_step_s = float(source_configuration.get("time_step_s", 0.0625))
    source_advective_number = max(
        float(entry["diagnostics"]["maximum_advective_number"])
        for entry in entries
    )
    source_diffusion_number = max(
        float(entry["diagnostics"]["maximum_diffusion_number"])
        for entry in entries
    )
    source_thermal_diffusion_number = max(
        float(entry["diagnostics"]["maximum_thermal_diffusion_number"])
        for entry in entries
    )

    wall_radius_m = 50.0
    initial_cloud_flat_to_flat_m = 14.0
    event_wind_speed_m_s = 3.9
    nominal_contact_time_s = (
        wall_radius_m - initial_cloud_flat_to_flat_m / 2.0
    ) / event_wind_speed_m_s
    window_duration_s = 20.0
    grids = (
        (5.0, 0.0625, 12_800),
        (2.5, 0.03125, 102_400),
        (1.25, 0.015625, 819_200),
    )
    dx2p5_cell_steps = 102_400 * int(round(window_duration_s / 0.0625))
    measured_rate = measured_dx2p5_wall_s / dx2p5_cell_steps
    grid_records = []
    for cell_size, time_step, cells in grids:
        steps = int(round(window_duration_s / time_step))
        estimate = measured_rate * cells * steps
        projected_advective_number = (
            source_advective_number
            * (time_step / source_time_step_s)
            * (source_cell_size_m / cell_size)
        )
        projected_diffusion_number = (
            source_diffusion_number
            * (time_step / source_time_step_s)
            * (source_cell_size_m / cell_size) ** 2
        )
        projected_thermal_diffusion_number = (
            source_thermal_diffusion_number
            * (time_step / source_time_step_s)
            * (source_cell_size_m / cell_size) ** 2
        )
        grid_records.append({
            "cell_size_m": cell_size,
            "time_step_s": time_step,
            "cells": cells,
            "steps": steps,
            "wall_height_cells": 5.0 / cell_size,
            "estimated_wall_clock_s": estimate,
            "estimated_wall_clock_h": estimate / 3600.0,
            "projected_maximum_advective_number": projected_advective_number,
            "projected_maximum_diffusion_number": projected_diffusion_number,
            "projected_maximum_thermal_diffusion_number": (
                projected_thermal_diffusion_number
            ),
        })
    gates = {
        "window_reaches_nominal_cloud_wall_contact": (
            window_duration_s >= nominal_contact_time_s
        ),
        "window_contains_at_least_8s_post_contact_response": (
            window_duration_s - nominal_contact_time_s >= 8.0
        ),
        "all_grids_resolve_wall_height_by_integer_cells": all(
            math.isclose(record["wall_height_cells"], round(record["wall_height_cells"]))
            for record in grid_records
        ),
        "grids_have_exact_two_to_one_refinement": (
            [record["cell_size_m"] for record in grid_records]
            == [5.0, 2.5, 1.25]
        ),
        "revised_time_steps_have_exact_two_to_one_refinement": (
            [record["time_step_s"] for record in grid_records]
            == [0.0625, 0.03125, 0.015625]
        ),
        "selection_uses_no_concentration_or_dose_observation": True,
        "projected_diffusion_numbers_at_most_0p5": all(
            max(
                record["projected_maximum_diffusion_number"],
                record["projected_maximum_thermal_diffusion_number"],
            ) <= 0.5
            for record in grid_records
        ),
        "projected_explicit_advection_diffusion_sum_at_most_one": all(
            record["projected_maximum_advective_number"]
            + max(
                record["projected_maximum_diffusion_number"],
                record["projected_maximum_thermal_diffusion_number"],
            ) <= 1.0
            for record in grid_records
        ),
    }
    revision_basis = {
        "reason": (
            "the superseded 1.25 m dt=0.03125 s member exceeded the explicit "
            "species and thermal diffusion-number limit"
        ),
        "policy": (
            "halve every accepted member time step so the revised three-grid "
            "path retains exact 2:1 temporal refinement"
        ),
        "failed_fine_attempt": None,
        "interrupted_fine_attempt": None,
    }
    if failed_fine_attempt is not None:
        failed_path = failed_fine_attempt.resolve()
        failure = json.loads(failed_path.read_text(encoding="utf-8"))
        if failure.get("schema") != "degali.execution-failure.v1":
            raise ValueError("unsupported fine-member failure schema")
        if failure.get("status") != "rejected_explicit_stability_limit":
            raise ValueError("fine-member failure is not an explicit stability rejection")
        revision_basis["failed_fine_attempt"] = {
            "path": str(failed_path),
            "sha256": _sha256(failed_path),
            "maximum_species_diffusion_number": failure["failure"][
                "maximum_species_diffusion_number"
            ],
            "maximum_thermal_diffusion_number": failure["failure"][
                "maximum_thermal_diffusion_number"
            ],
        }
    if interrupted_fine_attempt is not None:
        interrupted_path = interrupted_fine_attempt.resolve()
        interruption = json.loads(interrupted_path.read_text(encoding="utf-8"))
        if interruption.get("schema") != "degali.execution-interruption.v1":
            raise ValueError("unsupported fine-member interruption schema")
        if interruption.get("status") != (
            "intentionally_stopped_before_first_checkpoint"
        ):
            raise ValueError("fine-member interruption is not a pre-checkpoint stop")
        configuration = interruption.get("configuration", {})
        expected_configuration = {
            "cell_size_m": 1.25,
            "time_step_s": 0.015625,
            "duration_s": 20.0,
            "checkpoint_interval_s": 10.0,
            "corrector_tolerance": 0.0001,
            "observation_value_fields_accessed": [],
        }
        if configuration != expected_configuration:
            raise ValueError("fine-member interruption configuration is not frozen v2")
        decision = interruption.get("decision", {})
        if decision.get("physics_changed") is not False:
            raise ValueError("fine-member interruption must not change physics")
        if decision.get("result") != "not_scored":
            raise ValueError("fine-member interruption must remain unscored")
        revision_basis["interrupted_fine_attempt"] = {
            "path": str(interrupted_path),
            "sha256": _sha256(interrupted_path),
            "last_frame_time_s": interruption["partial_execution"][
                "last_frame_time_s"
            ],
            "frame_count": interruption["partial_execution"]["frame_count"],
            "reason": decision["reason"],
            "result": decision["result"],
            "physics_changed": decision["physics_changed"],
            "next_checkpoint_interval_s": decision[
                "next_checkpoint_interval_s"
            ],
        }
    return {
        "schema": SCHEMA,
        "status": (
            "accepted_concentration_blind_diffusion_stable_three_grid_impact_window_plan"
            if all(gates.values()) else "withheld_three_grid_impact_window_plan"
        ),
        "scope": (
            "three-grid numerical sensitivity of initial-cloud/fence impact; "
            "not a substitute for full-600s finest-grid convergence"
        ),
        "selection_basis": {
            "wall_radius_m": wall_radius_m,
            "initial_cloud_flat_to_flat_m": initial_cloud_flat_to_flat_m,
            "event_wind_speed_m_s": event_wind_speed_m_s,
            "nominal_advective_cloud_edge_to_wall_time_s": nominal_contact_time_s,
            "window_start_s": 0.0,
            "window_stop_s": window_duration_s,
            "initialization": (
                "fresh identical report-priority initial cloud on every grid; "
                "no cross-grid state prolongation"
            ),
            "observation_value_fields_accessed": [],
        },
        "revision_basis": revision_basis,
        "runtime_basis": {
            "checkpoint_journal_path": str(source),
            "checkpoint_journal_sha256": _sha256(source),
            "measured_dx2p5_first_20s_wall_clock_s": measured_dx2p5_wall_s,
            "measured_dx2p5_s_per_cell_step": measured_rate,
            "source_cell_size_m": source_cell_size_m,
            "source_time_step_s": source_time_step_s,
            "measured_maximum_advective_number": source_advective_number,
            "measured_maximum_diffusion_number": source_diffusion_number,
            "measured_maximum_thermal_diffusion_number": (
                source_thermal_diffusion_number
            ),
            "scaling": (
                "measured seconds per cell-step times target cell-steps; "
                "estimate excludes potentially superlinear AMG setup"
            ),
        },
        "grids": grid_records,
        "required_outputs": [
            "common-grid final and maximum vol% concentration fields",
            "front/top/side/wake regional concentration metrics",
            "4 vol% maximum-ever volume and downwind/crosswind/vertical extents",
            "60 receptor peak, dose and threshold-arrival metrics",
            "species mass, sensible energy, continuity and nonlinear convergence",
            "observed-order/GCI only where the three values are monotonic",
        ],
        "promotion_policy": {
            "impact_window_grid_sensitivity": (
                "eligible after all three runs and monotonic/GCI audit"
            ),
            "full_event_grid_convergence": "remains_withheld",
            "reason": (
                "the 1.25 m run covers the obstacle-impact window, not the full "
                "600 s sensor and handoff history"
            ),
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-journal", type=Path, required=True)
    parser.add_argument("--failed-fine-attempt", type=Path)
    parser.add_argument("--interrupted-fine-attempt", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    record = build_plan(
        args.checkpoint_journal,
        failed_fine_attempt=args.failed_fine_attempt,
        interrupted_fine_attempt=args.interrupted_fine_attempt,
    )
    destination = args.output_json.resolve()
    destination.parent.mkdir(parents=True, exist_ok=False)
    destination.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
    )
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
