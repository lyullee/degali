import json

from tools.audit_ffi_field_evidence_manifest import build_intake


def test_ffi_intake_is_hash_pinned_and_fail_closed(tmp_path):
    data_root = tmp_path / "extracted"
    csv_dir = data_root / "23_FFI_P101_T00027"
    csv_dir.mkdir(parents=True)
    csv_path = csv_dir / "23_FFI_P101_T00027.csv"
    csv_path.write_text(
        ",mass flow meter 2 [g/s],flow time [s],h2 sensor time,"
        "sensor 1 h2 concentration [%]\n"
        "0,0.09,0.0,0.0,0.1\n",
        encoding="utf-8",
    )
    registry = tmp_path / "registry.json"
    registry.write_text(json.dumps({
        "obstacle": {"quantitative_geometry_allowed": False},
        "sensor_instrumentation": {"sensor_count": 29},
        "time_operator": {"dataset_source_and_sensor_clocks_are_distinct": True},
        "promotion": {"gate_codes": ["geometry_unresolved"]},
        "case_crosswalk": [{
            "article_case": "O1",
            "article_trials": ["test-02-1A"],
            "dataset_files": ["23_FFI_P101_T00027.csv"],
        }],
    }), encoding="utf-8")

    output = tmp_path / "intake.json"
    record = build_intake(
        data_root, registry, article_case="O1", output=output,
    )

    assert record["status"] == "withheld"
    assert record["promotion_allowed"] is False
    assert record["channels"]["source_boundary"]["status"] == "conditional"
    assert record["channels"]["weather"]["status"] == "withheld"
    assert record["channels"]["common_clock"]["status"] == "withheld"
    assert len(record["selected_trials"][0]["dataset_sha256"]) == 64
