import json

import pytest

from tools.plan_thorney_ti021_three_grid_impact_window import build_plan


def test_three_grid_window_plan_is_observation_blind_and_two_to_one(tmp_path):
    journal = tmp_path / "checkpoint.json"
    journal.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-checkpoint.v1",
        "segments": [
            {"stop_time_s": 10.0, "wall_clock_s": 300.0, "diagnostics": {
                "maximum_advective_number": 0.49,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
            {"stop_time_s": 20.0, "wall_clock_s": 180.0, "diagnostics": {
                "maximum_advective_number": 0.31,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
        ],
    }), encoding="utf-8")
    record = build_plan(journal)
    assert record["all_gates_pass"]
    assert record["selection_basis"]["observation_value_fields_accessed"] == []
    assert [item["cell_size_m"] for item in record["grids"]] == [5.0, 2.5, 1.25]
    assert [item["time_step_s"] for item in record["grids"]] == [
        0.0625, 0.03125, 0.015625,
    ]
    assert record["grids"][2]["estimated_wall_clock_s"] == pytest.approx(
        32.0 * 480.0
    )
    assert record["grids"][2]["projected_maximum_diffusion_number"] == (
        pytest.approx(0.27)
    )
    assert record["gates"]["projected_diffusion_numbers_at_most_0p5"]
    assert record["promotion_policy"]["full_event_grid_convergence"] == (
        "remains_withheld"
    )


def test_three_grid_window_plan_requires_twenty_committed_seconds(tmp_path):
    journal = tmp_path / "checkpoint.json"
    journal.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-checkpoint.v1",
        "segments": [{"stop_time_s": 10.0, "wall_clock_s": 300.0,
                      "diagnostics": {
                          "maximum_advective_number": 0.49,
                          "maximum_diffusion_number": 0.27,
                          "maximum_thermal_diffusion_number": 0.26,
                      }}],
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="20 committed"):
        build_plan(journal)


def test_three_grid_window_plan_binds_explicit_stability_failure(tmp_path):
    journal = tmp_path / "checkpoint.json"
    journal.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-checkpoint.v1",
        "segments": [
            {"stop_time_s": 10.0, "wall_clock_s": 300.0, "diagnostics": {
                "maximum_advective_number": 0.49,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
            {"stop_time_s": 20.0, "wall_clock_s": 180.0, "diagnostics": {
                "maximum_advective_number": 0.31,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
        ],
    }), encoding="utf-8")
    failure = tmp_path / "failure.json"
    failure.write_text(json.dumps({
        "schema": "degali.execution-failure.v1",
        "status": "rejected_explicit_stability_limit",
        "failure": {
            "maximum_species_diffusion_number": 0.553501,
            "maximum_thermal_diffusion_number": 0.523581,
        },
    }), encoding="utf-8")
    record = build_plan(journal, failed_fine_attempt=failure)
    bound = record["revision_basis"]["failed_fine_attempt"]
    assert bound["sha256"]
    assert bound["maximum_species_diffusion_number"] == pytest.approx(0.553501)


def test_three_grid_window_plan_binds_unscored_operational_interruption(tmp_path):
    journal = tmp_path / "checkpoint.json"
    journal.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-checkpoint.v1",
        "segments": [
            {"stop_time_s": 10.0, "wall_clock_s": 300.0, "diagnostics": {
                "maximum_advective_number": 0.49,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
            {"stop_time_s": 20.0, "wall_clock_s": 180.0, "diagnostics": {
                "maximum_advective_number": 0.31,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
        ],
    }), encoding="utf-8")
    interruption = tmp_path / "interruption.json"
    interruption.write_text(json.dumps({
        "schema": "degali.execution-interruption.v1",
        "status": "intentionally_stopped_before_first_checkpoint",
        "configuration": {
            "cell_size_m": 1.25,
            "time_step_s": 0.015625,
            "duration_s": 20.0,
            "checkpoint_interval_s": 10.0,
            "corrector_tolerance": 0.0001,
            "observation_value_fields_accessed": [],
        },
        "partial_execution": {
            "frame_count": 22,
            "last_frame_time_s": 0.328125,
        },
        "decision": {
            "result": "not_scored",
            "reason": "first atomic checkpoint is too far from process start",
            "physics_changed": False,
            "next_checkpoint_interval_s": 1.0,
        },
    }), encoding="utf-8")
    record = build_plan(journal, interrupted_fine_attempt=interruption)
    bound = record["revision_basis"]["interrupted_fine_attempt"]
    assert bound["sha256"]
    assert bound["last_frame_time_s"] == pytest.approx(0.328125)
    assert bound["frame_count"] == 22
    assert bound["result"] == "not_scored"
    assert bound["physics_changed"] is False
    assert bound["next_checkpoint_interval_s"] == pytest.approx(1.0)


def test_three_grid_window_plan_rejects_interruption_that_changes_physics(
    tmp_path,
):
    journal = tmp_path / "checkpoint.json"
    journal.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-checkpoint.v1",
        "segments": [
            {"stop_time_s": 10.0, "wall_clock_s": 300.0, "diagnostics": {
                "maximum_advective_number": 0.49,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
            {"stop_time_s": 20.0, "wall_clock_s": 180.0, "diagnostics": {
                "maximum_advective_number": 0.31,
                "maximum_diffusion_number": 0.27,
                "maximum_thermal_diffusion_number": 0.26,
            }},
        ],
    }), encoding="utf-8")
    interruption = tmp_path / "interruption.json"
    interruption.write_text(json.dumps({
        "schema": "degali.execution-interruption.v1",
        "status": "intentionally_stopped_before_first_checkpoint",
        "configuration": {
            "cell_size_m": 1.25,
            "time_step_s": 0.015625,
            "duration_s": 20.0,
            "checkpoint_interval_s": 10.0,
            "corrector_tolerance": 0.0001,
            "observation_value_fields_accessed": [],
        },
        "partial_execution": {"frame_count": 1, "last_frame_time_s": 0.0},
        "decision": {
            "result": "not_scored",
            "reason": "invalid physics mutation",
            "physics_changed": True,
            "next_checkpoint_interval_s": 1.0,
        },
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="must not change physics"):
        build_plan(journal, interrupted_fine_attempt=interruption)
