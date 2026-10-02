"""Black-box guard for the no-fit E3.4 source audit command."""

import json
import subprocess
import sys
from pathlib import Path

import openpyxl


def test_e34_audit_requires_declared_window_and_writes_no_raw_data(tmp_path):
    workbook = tmp_path / "20200320-Concrete01-Final.xlsx"
    book = openpyxl.Workbook()
    sheet = book.active
    for _ in range(6):
        sheet.append([])
    sheet.append(["Sync. Time [s]", "m(LH2) [g]"])
    sheet.append([0.0, 1000.0])
    sheet.append([1.0, 990.0])
    sheet.append([2.0, 980.0])
    sheet.append([3.0, 970.0])
    book.save(workbook)
    manifest = tmp_path / "declared.json"
    manifest.write_text(json.dumps({
        "workbook": str(workbook), "trial_id": "20200320-Concrete01",
        "pool_contact_start_s": 0.0,
        "window_start_s": 1.0, "window_end_s": 3.0, "area_m2": 0.25,
        "substrate": {
            "conductivity_w_m_k": 1.5, "density_kg_m3": 2200.0,
            "heat_capacity_j_kg_k": 900.0, "initial_temperature_k": 293.15,
            "depth_m": 1.0,
        },
    }), encoding="utf-8")
    output = tmp_path / "audit.json"
    root = Path(__file__).resolve().parents[1]
    subprocess.run([
        sys.executable, str(root / "tools" / "audit_e34_pool_source.py"),
        str(manifest), "--output", str(output),
    ], check=True, capture_output=True, text=True)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["source_files_distributed"] is False
    assert result["coefficient_fitted"] is False
    assert result["report_to_workbook_identity_confirmed"] is False
    assert result["observed"]["mass_loss_kg"] > 0.0


def test_e34_audit_refuses_a_report_trial_mapped_to_a_different_workbook(tmp_path):
    workbook = tmp_path / "20200324-Concrete02-Final.xlsx"
    book = openpyxl.Workbook()
    sheet = book.active
    for _ in range(6):
        sheet.append([])
    sheet.append(["Sync. Time [s]", "m(LH2) [g]"])
    for time, mass in [(0, 1000), (1, 990), (2, 980), (3, 970)]:
        sheet.append([time, mass])
    book.save(workbook)
    manifest = tmp_path / "declared.json"
    manifest.write_text(json.dumps({
        "workbook": str(workbook), "trial_id": "20200324-Concrete02",
        "report_trial_id": "20200423-Concrete02",
        "pool_contact_start_s": 0.0, "window_start_s": 1.0,
        "window_end_s": 3.0, "area_m2": 0.25,
        "substrate": {
            "conductivity_w_m_k": 1.5, "density_kg_m3": 2200.0,
            "heat_capacity_j_kg_k": 900.0, "initial_temperature_k": 293.15,
            "depth_m": 1.0,
        },
    }), encoding="utf-8")
    root = Path(__file__).resolve().parents[1]
    failed = subprocess.run([
        sys.executable, str(root / "tools" / "audit_e34_pool_source.py"),
        str(manifest), "--output", str(tmp_path / "audit.json"),
    ], capture_output=True, text=True)
    assert failed.returncode != 0
    assert "report-to-workbook numerical comparison" in failed.stderr
