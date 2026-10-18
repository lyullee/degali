from pathlib import Path
import hashlib
import json
import runpy

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TOOL = runpy.run_path(ROOT / "tools/score_preslhy_e34_dynamic_component.py")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_channel_score_uses_frozen_peak_dose_arrival_and_four_vol_operators():
    time = np.asarray((0.0, 1.0, 2.0, 3.0))
    observed = np.asarray((0.0, 1.0, 5.0, 0.0))
    model = np.asarray((0.0, 0.0, 4.0, 1.0))

    score = TOOL["_channel_score"](
        time,
        observed,
        model,
        log_floor=0.1,
        arrival_threshold=1.0,
        flammability_threshold=4.0,
    )

    assert score["mae_vol_percent_point"] == 0.75
    assert score["absolute_peak_error_vol_percent_point"] == 1.0
    assert score["observed_dose_vol_percent_s"] == 6.0
    assert score["model_dose_vol_percent_s"] == 4.5
    assert score["observed_arrival_1vol_s"] == 1.0
    assert score["model_arrival_1vol_s"] == 2.0
    assert score["normalized_arrival_time_error"] == 1.0
    classification = score["four_vol_sample_classification"]
    assert classification["true_positive"] == 1
    assert classification["false_negative"] == 0


def test_scoring_policy_freeze_is_observation_blind_and_non_obstacle(tmp_path):
    replay = tmp_path / "replay.json"
    replay_record = {
        "sealed_observation_operator": {
            "usable_channels": ["CH45", "CH55"],
            "excluded_channels": ["CH35"],
            "physical_observation_clock": "physical clock",
            "declared_sample_transport_delay_s": 35.0,
        },
        "duration_s": 10.0,
    }
    replay.write_text(json.dumps(replay_record), encoding="utf-8")
    freeze = tmp_path / "model-configuration-freeze.json"
    freeze_record = {
        "schema": "degali.preslhy-e34-dynamic-component-configuration-freeze.v1",
        "replay_input": {"path": str(replay), "sha256": _sha256(replay)},
        "model_configuration_sha256": "fixture-configuration",
        "observation_control": {"hydrogen_observation_values_opened": False},
    }
    freeze.write_text(json.dumps(freeze_record), encoding="utf-8")

    policy = TOOL["freeze_scoring_policy"](tmp_path)

    assert policy["observation_control"][
        "scoring_policy_frozen_before_observation_values_opened"
    ] is True
    assert policy["policy"]["arrival_threshold_vol_percent"] == 1.0
    assert policy["policy"]["flammability_threshold_vol_percent"] == 4.0
    assert policy["scope"]["obstacle_transport_model_qualification_allowed"] is False
    assert (tmp_path / "scoring-policy-freeze.json").is_file()
