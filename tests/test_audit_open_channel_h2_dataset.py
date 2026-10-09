import json

import pytest

from tools.audit_open_channel_h2_dataset import (
    audit_open_channel_h2_dataset,
    main,
    write_audit_json,
)


def _write_csv(path):
    path.write_text(
        "flow time [s],mass flow meter 1 [g/s],h2 sensor time,"
        "sensor 1 h2 concentration [%],sensor 2 h2 concentration [%]\n"
        "0,0.1,0,0.1,0.2\n"
        "0.01,0.2,0.5,1.2,0.3\n"
        "0.02,0.3,1.0,2.1,4.2\n"
        "0.03,0.2,1.5,0.8,0.3\n",
        encoding="utf-8",
    )


def test_open_channel_audit_keeps_boundary_and_hashes_files(tmp_path):
    root = tmp_path / "extracted"
    root.mkdir()
    _write_csv(root / "test.csv")

    record = audit_open_channel_h2_dataset(root)

    assert record["status"] == "boundary_only"
    assert record["file_count"] == 1
    assert record["parsed_file_count"] == 1
    assert record["promotion_allowed"] is False
    assert record["quantitative_lh2_pool_validation_allowed"] is False
    row = record["files"][0]
    assert len(row["sha256"]) == 64
    assert row["sensor_count"] == 2
    assert row["threshold_sensor_count"] == 1
    assert row["threshold_peak_max_percent"] == pytest.approx(4.2)
    assert "weather_geometry_common_clock_missing" in record["gate_codes"]


def test_open_channel_audit_missing_data_is_not_a_failure(tmp_path, capsys):
    record = audit_open_channel_h2_dataset(tmp_path / "missing")
    assert record["status"] == "not_available"
    assert record["promotion_allowed"] is False

    assert main([str(tmp_path / "missing")]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "not_available"


def test_open_channel_audit_writer_is_exclusive(tmp_path):
    root = tmp_path / "extracted"
    root.mkdir()
    _write_csv(root / "test.csv")
    record = audit_open_channel_h2_dataset(root)
    output = tmp_path / "audit.json"
    write_audit_json(record, output)
    assert json.loads(output.read_text(encoding="utf-8"))["schema"] == record["schema"]
    with pytest.raises(FileExistsError):
        write_audit_json(record, output)
