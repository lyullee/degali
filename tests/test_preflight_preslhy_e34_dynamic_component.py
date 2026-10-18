from pathlib import Path
import hashlib
import json
import runpy

import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = runpy.run_path(
    ROOT / "tools/preflight_preslhy_e34_dynamic_component.py"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _replay_input(tmp_path: Path) -> Path:
    source = tmp_path / "source-freeze.json"
    audit = tmp_path / "audit.json"
    source.write_text("{}\n", encoding="utf-8")
    audit.write_text("{}\n", encoding="utf-8")
    replay = {
        "schema": "degali.preslhy-e34-dynamic-replay-input.v1",
        "trial_id": "fixture-Sand03Wind",
        "window_id": "fixture-window-01",
        "duration_s": 0.1,
        "source_freeze": {"path": str(source), "sha256": _sha256(source)},
        "dynamic_evidence_audit": {
            "path": str(audit), "sha256": _sha256(audit),
        },
        "source_schedule": {
            "time_s": [0.0, 0.1],
            "rate_kg_s": [0.001, 0.0],
            "source_id": "fixture-source",
            "rate_operator": "piecewise_constant",
        },
        "inlet_velocity_history": {
            "time_s": [0.0, 0.1],
            "u_m_s": [0.5, 0.5],
            "v_m_s": [0.0, 0.0],
            "w_m_s": [0.0, 0.0],
            "history_id": "fixture-wind",
            "interpolation": "piecewise_constant",
        },
        "wind_reduction": {
            "coordinate_frame": (
                "instrument_window_mean_horizontal_wind_aligned_not_site_registered"
            ),
            "whole_inlet_assumption": (
                "one point-history block mean times the declared spatial inlet shape"
            ),
        },
        "sealed_observation_operator": {
            "values_exported": False,
            "values_summarized": False,
            "sealed_until_model_configuration_is_frozen": True,
            "usable_channels": ["CH45", "CH55"],
            "positions_pool_centre_m": {
                "CH45": [0.0, 0.0, 0.45],
                "CH55": [0.0, 0.0, 0.55],
            },
        },
        "evidence_decision": {
            "dynamic_component_replay_input_allowed": True,
            "obstacle_transport_model_qualification_allowed": False,
        },
    }
    path = tmp_path / "dynamic-replay-input.json"
    path.write_text(json.dumps(replay), encoding="utf-8")
    return path


def test_freeze_is_observation_blind_and_explicitly_not_obstacle_validation(tmp_path):
    replay = _replay_input(tmp_path)
    output = tmp_path / "freeze"
    record, config, receptors = TOOL["freeze_model_configuration"](
        replay, output,
    )

    assert record["observation_control"][
        "model_configuration_frozen_before_observation_values_opened"
    ] is True
    assert record["observation_control"]["hydrogen_observation_values_opened"] is False
    assert record["evidence_decision"][
        "obstacle_transport_model_qualification_allowed"
    ] is False
    assert config.active_momentum_solver == "finite_volume_projection"
    assert config.active_inlet_velocity_history is not None
    assert config.obstacles == ()
    assert [item.label for item in receptors] == ["CH45", "CH55"]
    saved = json.loads(
        (output / "model-configuration-freeze.json").read_text(encoding="utf-8")
    )
    assert saved["model_configuration_sha256"] == record[
        "model_configuration_sha256"
    ]


def test_preflight_writes_hash_bound_numerical_gate_record(tmp_path):
    replay = _replay_input(tmp_path)
    output = tmp_path / "preflight"
    result = TOOL["run_preflight"](
        replay, output, preflight_duration_s=0.1,
    )

    assert result["observation_control"]["hydrogen_observation_values_opened"] is False
    assert result["all_predeclared_numerical_gates_passed"] is True
    assert result["evidence_decision"]["full_dynamic_component_run_allowed"] is True
    assert result["evidence_decision"][
        "obstacle_transport_model_qualification_allowed"
    ] is False
    assert (output / "preflight-result.json").is_file()
    assert result["configuration_freeze"]["sha256"] == _sha256(
        output / "model-configuration-freeze.json"
    )


def test_full_component_reconstructs_exact_freeze_before_model_output(tmp_path):
    replay = _replay_input(tmp_path)
    output = tmp_path / "full"
    TOOL["run_preflight"](replay, output, preflight_duration_s=0.1)
    result = TOOL["run_full_component"](output)

    assert result["observation_control"][
        "model_output_completed_before_observation_values_opened"
    ] is True
    assert result["all_predeclared_numerical_gates_passed"] is True
    assert result["evidence_decision"]["sealed_observation_scoring_allowed"] is True
    assert result["evidence_decision"][
        "obstacle_transport_model_qualification_allowed"
    ] is False
    assert (output / "full-run-result.json").is_file()
    model_output = output / "full-run-model-output.npz"
    assert model_output.is_file()
    assert result["model_output"]["sha256"] == _sha256(model_output)


def test_freeze_rejects_value_bearing_observation_field(tmp_path):
    replay = _replay_input(tmp_path)
    value = json.loads(replay.read_text(encoding="utf-8"))
    value["sealed_observation_operator"]["observed_values"] = [1.0, 2.0]
    replay.write_text(json.dumps(value), encoding="utf-8")

    with pytest.raises(ValueError, match="prohibited H2 value field"):
        TOOL["freeze_model_configuration"](replay, tmp_path / "freeze")
