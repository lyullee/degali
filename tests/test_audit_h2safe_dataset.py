import csv
import hashlib
from pathlib import Path
from zipfile import ZipFile

from tools.audit_h2safe_dataset import audit_h2safe_dataset


def _write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle).writerows(rows)


def _fixture(raw: Path):
    extracted = raw / "extracted"
    for lab, tests, count in (
        ("Lab1_Small", ("test1",), 24),
        ("Lab2_Large", ("test1", "test2", "test3", "test4"), 37),
    ):
        directory = extracted / lab
        _write_csv(
            directory / "sensors.csv",
            [["id", "X", "Y", "Z"]]
            + [[f"p{i}", str(i), "1", "2"] for i in range(1, count + 1)],
        )
        for test in tests:
            _write_csv(
                directory / f"{test}.csv",
                [["time", "p1", "p2"], ["0", "0", "#N/A"], ["1", "3", "#N/A"]],
            )
        stem = "Lab1_small" if lab == "Lab1_Small" else "Lab2_Large"
        (directory / f"{stem}.STEP").write_text(
            "#1 = PRODUCT ( 'room_box', 'room_box', '', (#2) );\n"
            "#3 = PRODUCT ( 'east_supply_vent', 'east_supply_vent', '', (#4) );\n"
            "#5 = PRODUCT ( 'west_exhaust_vent', 'west_exhaust_vent', '', (#6) );\n",
            encoding="ascii",
        )
        (directory / f"{stem}.STL").write_text(
            "solid lab\nfacet normal 0 0 1\nouter loop\n"
            "vertex -30 -5 -20\nvertex 50 -5 -20\nvertex 50 50 50\n"
            "endloop\nendfacet\nendsolid lab\n",
            encoding="ascii",
        )
    (raw / "HVAC_information.docx").write_bytes(b"DOCX")
    with ZipFile(raw / "Dataset.zip", "w") as archive:
        for path in sorted(extracted.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(extracted).as_posix())
    return (
        hashlib.sha256((raw / "Dataset.zip").read_bytes()).hexdigest(),
        hashlib.sha256((raw / "HVAC_information.docx").read_bytes()).hexdigest(),
    )


def test_h2safe_audit_emits_five_conditional_fail_closed_manifests(tmp_path):
    archive_sha, hvac_sha = _fixture(tmp_path)
    record = audit_h2safe_dataset(
        tmp_path,
        expected_archive_sha256=archive_sha,
        expected_hvac_sha256=hvac_sha,
    )

    assert record["status"] == "conditional_component_evidence_helium_surrogate"
    assert record["event_count"] == 5
    assert not record["promotion_allowed"]
    assert all(event["manifest"]["status"] == "conditional" for event in record["events"])
    assert all(
        not event["manifest"]["quantitative_concentration_validation_allowed"]
        for event in record["events"]
    )
    assert all(
        "concentration_unit_declaration" in event["manifest"]["missing_requirements"]
        for event in record["events"]
    )
    assert all(
        member["extracted_matches_archive"]
        for member in record["source_integrity"]["archive_members"]
    )
    assert all(
        event["coordinate_frame_check"]["status"] == "accepted"
        for event in record["events"]
    )
    assert all(
        event["coordinate_axis_mapping"]["solver_mapping"]
        == {"x": "X", "y": "Z", "z": "Y"}
        for event in record["events"]
    )
    assert all(
        event["step_semantic_inventory"]["boundary_face_mapping_status"].startswith("withheld")
        for event in record["events"]
    )
    lab1 = next(event for event in record["events"] if event["lab_id"] == "Lab1_Small")
    lab2 = next(event for event in record["events"] if event["lab_id"] == "Lab2_Large")
    assert lab1["hvac_boundary_summary"]["reported_relative_flow_imbalance"] < 0.01
    assert lab1["hvac_boundary_summary"]["mass_balance_status"] == (
        "accepted_reported_flow_balance"
    )
    assert lab2["hvac_boundary_summary"]["reported_outlet_flow_m3_s"] == 0.0
    assert lab2["hvac_boundary_summary"]["mass_balance_status"] == (
        "withheld_external_outlet_not_reported"
    )
    assert "external_hvac_outlet_flow_boundary" in (
        lab2["manifest"]["missing_requirements"]
    )
    assert "external_hvac_outlet_flow_missing" in lab2["manifest"]["gate_codes"]


def test_h2safe_audit_withholds_digest_drift(tmp_path):
    _, hvac_sha = _fixture(tmp_path)
    record = audit_h2safe_dataset(
        tmp_path,
        expected_archive_sha256="0" * 64,
        expected_hvac_sha256=hvac_sha,
    )
    assert record["status"] == "withheld"
    assert any("SHA-256" in error for error in record["errors"])
