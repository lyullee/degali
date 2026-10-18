import json

import pytest

from tools.project_thorney_ti021_runtime import project_runtime


def test_runtime_projection_uses_latest_progress_interval(tmp_path):
    path = tmp_path / "progress.jsonl"
    path.write_text("\n".join([
        json.dumps({"model_time_s": 0.0,
                    "elapsed_wall_clock_s_this_process": 25.0}),
        json.dumps({"model_time_s": 0.25,
                    "elapsed_wall_clock_s_this_process": 581.0}),
        json.dumps({"model_time_s": 0.5,
                    "elapsed_wall_clock_s_this_process": 1126.0}),
    ]) + "\n", encoding="utf-8")
    record = project_runtime(
        path, duration_s=20.0, checkpoint_interval_s=1.0,
    )
    assert record["sample_count"] == 3
    assert record["latest_completed_fraction"] == pytest.approx(0.025)
    assert record["latest_interval"]["wall_seconds_per_model_second"] == (
        pytest.approx(2180.0)
    )
    assert record["projection"]["first_checkpoint_elapsed_wall_clock_s"] == (
        pytest.approx(2216.0)
    )
    assert record["observation_value_fields_accessed"] == []


def test_runtime_projection_rejects_duplicate_or_nonmonotonic_samples(tmp_path):
    path = tmp_path / "progress.jsonl"
    path.write_text("\n".join([
        json.dumps({"model_time_s": 0.0,
                    "elapsed_wall_clock_s_this_process": 1.0}),
        json.dumps({"model_time_s": 0.5,
                    "elapsed_wall_clock_s_this_process": 2.0}),
        json.dumps({"model_time_s": 0.5,
                    "elapsed_wall_clock_s_this_process": 3.0}),
    ]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="increase strictly"):
        project_runtime(path, duration_s=20.0, checkpoint_interval_s=1.0)
