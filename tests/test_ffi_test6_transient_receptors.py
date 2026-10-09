import json
from pathlib import Path
import runpy

import pytest


@pytest.fixture(scope="module")
def replay_tool():
    return runpy.run_path("tools/audit_ffi_test6_transient_receptors.py")


def _history(path):
    path.write_text(
        "time_s,wind_speed_ms,wind_direction_from_deg\n"
        "20,2.3,245.4\n"
        "80,2.3,245.4\n"
        "140,2.3,245.4\n",
        encoding="utf-8",
    )


def test_transient_replay_history_contract_is_strict(tmp_path, replay_tool):
    history = tmp_path / "wind.csv"
    _history(history)
    parsed, row_count = replay_tool["read_wind_history_csv"](history)

    assert row_count == 3
    assert parsed.arrays()[1].tolist() == [2.3, 2.3, 2.3]

    history.write_text(
        "time_s,wind_speed_ms,wind_direction_from_deg,extra\n"
        "20,2.3,245.4,0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="header must be exactly"):
        replay_tool["read_wind_history_csv"](history)


@pytest.mark.skipif(
    not Path("reference/spadeadam/conditions.csv").is_file(),
    reason="requires the locally controlled FFI/SPADEADAM reference tables",
)
def test_transient_replay_emits_hash_pinned_compact_observation_artifact(
    tmp_path, replay_tool,
):
    history = tmp_path / "wind.csv"
    _history(history)
    output = tmp_path / "replay.json"

    assert replay_tool["main"]([
        str(history),
        "--table-winds", "2.3",
        "--response-t90", "6",
        "--start", "20",
        "--end", "140",
        "--output", str(output),
    ]) == 0

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["schema"] == "degali.ffi-test6-transient-receptor-execution.v1"
    assert len(payload["replay"]["sensor_rows"]) == 30
    assert payload["replay"]["summary"]["sensor_count"] == 30
    assert payload["replay"]["promotion_allowed"] is False
    assert payload["replay"]["validation_qualified"] is False
    assert payload["input"]["history_row_count"] == 3
    assert len(payload["input"]["history_sha256"]) == 64
    assert payload["input"]["history_diagnostics"] == {
        "sample_count": 3,
        "time_start_s": 20.0,
        "time_end_s": 140.0,
        "minimum_dt_s": 60.0,
        "maximum_dt_s": 60.0,
        "uniform_sampling": True,
        "wind_speed_min_m_s": 2.3,
        "wind_speed_max_m_s": 2.3,
        "direction_from_deg_first": 245.4,
        "direction_from_deg_last": 245.4,
    }
    assert payload["input"]["response_t90_s"] == 6.0
    assert payload["observation_operator"]["observed_values_are_lower_bounds"] is True
    assert all(
        row["window"]["sample_count"] == 3
        for row in payload["replay"]["sensor_rows"]
    )
    verification = replay_tool["verify_replay_artifact"](output)
    assert verification["schema"] == "degali.ffi-test6-transient-receptor-verification.v1"
    assert verification["history_sha256_verified"] is True
    assert verification["reference_provenance_verified"] is True
    assert verification["promotion_allowed"] is False

    tampered = json.loads(json.dumps(payload))
    tampered["replay"]["sensor_rows"][0]["response_t90_s"] = 0.0
    tampered_path = tmp_path / "tampered-replay.json"
    tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="response time"):
        replay_tool["verify_replay_artifact"](tampered_path)

    tampered["replay"]["sensor_rows"][0]["response_t90_s"] = 6.0
    tampered["input"]["history_diagnostics"]["maximum_dt_s"] = 30.0
    tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="history_diagnostics.maximum_dt_s"):
        replay_tool["verify_replay_artifact"](tampered_path)

    assert replay_tool["main"]([
        str(history), "--table-winds", "2.3", "--start", "20",
        "--end", "140", "--output", str(output),
    ]) == 1
    assert replay_tool["main"](["--verify", str(output)]) == 0
