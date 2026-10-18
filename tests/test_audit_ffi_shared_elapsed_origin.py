import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "audit_ffi_shared_elapsed_origin",
    ROOT / "tools" / "audit_ffi_shared_elapsed_origin.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_audit_records_matching_elapsed_zero_without_claiming_clock_sync(tmp_path):
    csv_path = tmp_path / "run.csv"
    csv_path.write_text(
        "flow time [s],mass flow meter 1 [g/s],h2 sensor time,sensor 1 h2 concentration [%]\n"
        "0,0,0,not-used\n"
        "0.01,0.1,0.5,also-not-used\n"
        "0.02,0.2,1.0,ignored\n",
        encoding="utf-8",
    )

    record = MODULE.audit_run(csv_path, "O1", "test-02-1A")

    assert record["recorded_elapsed_zero_matches"] is True
    assert record["source_window_covers_sensor_window"] is False
    assert record["source_log_ends_with_nonzero_recorded_rate"] is True
    assert record["flow"]["median_positive_sample_interval_s"] == 0.01
    assert record["sensor"]["median_positive_sample_interval_s"] == 0.5
    assert "drift" in record["interpretation"]


def test_audit_rejects_decreasing_clocks(tmp_path):
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text(
        "flow time [s],mass flow meter 1 [g/s],h2 sensor time\n"
        "0,0,0\n"
        "1,0.1,0.5\n"
        "0.5,0.2,1.0\n",
        encoding="utf-8",
    )

    try:
        MODULE.audit_run(csv_path, "O1", "test-02-1A")
    except ValueError as error:
        assert "flow clock is not nondecreasing" in str(error)
    else:
        raise AssertionError("decreasing flow clock should be rejected")
