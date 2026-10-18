from pathlib import Path
import hashlib
import json
import runpy

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TOOL = runpy.run_path(
    ROOT / "tools/diagnose_preslhy_e34_dynamic_component_score.py"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_diagnosis_retains_failed_primary_and_detects_initial_mismatch(tmp_path):
    replay_path = tmp_path / "replay.json"
    replay_path.write_text(json.dumps({
        "inlet_velocity_history": {
            "time_s": [0.0, 2.0, 4.0],
            "u_m_s": [1.0, 0.1, 0.1],
        },
        "sealed_observation_operator": {"usable_channels": ["CH45", "CH55"]},
        "dynamic_evidence_audit": {},
        "source_schedule_mass_relative_difference": 0.04,
    }), encoding="utf-8")
    freeze_path = tmp_path / "model-configuration-freeze.json"
    freeze_path.write_text(json.dumps({
        "replay_input": {"path": str(replay_path)},
    }), encoding="utf-8")
    arrays_path = tmp_path / "score-arrays.npz"
    np.savez_compressed(
        arrays_path,
        time_s=np.asarray((0.0, 1.0, 2.0, 3.0)),
        observed_CH45_h2_vol_percent=np.asarray((5.0, 2.0, 1.0, 0.0)),
        model_CH45_h2_vol_percent=np.asarray((0.0, 0.0, 10.0, 9.0)),
        observed_CH55_h2_vol_percent=np.asarray((2.0, 1.0, 0.5, 0.0)),
        model_CH55_h2_vol_percent=np.asarray((0.0, 0.0, 8.0, 7.0)),
    )
    score_path = tmp_path / "dynamic-component-score.json"
    score_path.write_text(json.dumps({
        "schema": "degali.preslhy-e34-dynamic-component-score.v1",
        "score_arrays": {"path": str(arrays_path), "sha256": _sha256(arrays_path)},
        "configuration_freeze": {"path": str(freeze_path)},
        "aggregate_targets": {"median_absolute_relative_dose_error": 1.2},
    }), encoding="utf-8")

    record = TOOL["diagnose_score"](tmp_path)

    assert record["status"] == "diagnosed_failed_primary_without_calibration"
    assert record["diagnostic_controls"]["parameter_adjustment_performed"] is False
    assert record["diagnostic_controls"]["primary_frozen_score_replaced"] is False
    assert record["cross_channel_findings"][
        "initial_condition_mismatch_observed"
    ] is True
    assert record["cross_channel_findings"][
        "model_peaks_in_lowest_inlet_wind_quartile"
    ] is True
    assert (tmp_path / "dynamic-component-diagnosis.png").is_file()
