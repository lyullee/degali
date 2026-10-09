import hashlib
import json
from pathlib import Path

import pytest

from degali.cli import main
from degali.addons.field_batch_io import (
    FIELD_BATCH_INPUT_SCHEMA,
    field_batch_input_from_mapping,
    read_field_batch_input_json,
)


ROOT = Path(__file__).resolve().parents[1]


def _write_case(tmp_path: Path) -> Path:
    case = json.loads(
        (ROOT / "docs" / "field-screening-case.example.json").read_text(encoding="utf-8")
    )
    case["scenario"]["source"]["duration_s"] = 0.1
    case["scenario"]["source"]["mass_flow_kg_s"] = {
        "nominal": 0.265, "unit": "kg/s", "source": "FT-source-calibration-A",
    }
    case["scenario"]["weather"]["speed_m_s"] = {
        "nominal": 2.0, "unit": "m/s", "source": "met-mast-A",
    }
    case["transport"].update({
        "length_m": 4.0, "height_m": 2.0, "nx": 10, "nz": 10,
        "time_step_s": 0.005, "duration_s": 0.1, "source_sigma_m": 0.2,
    })
    case["obstacles"] = []
    case["sensor_deployments"] = []
    path = tmp_path / "field-case.json"
    path.write_text(json.dumps(case, indent=2), encoding="utf-8")
    return path


def _write_source(tmp_path: Path) -> Path:
    csv_path = tmp_path / "source.csv"
    csv_path.write_text(
        "time_s,rate_kg_s,rate_lower_kg_s,rate_upper_kg_s\n"
        "0.0,0.12,0.10,0.14\n"
        "0.05,0.08,0.06,0.10\n"
        "0.1,0.0,0.0,0.0\n",
        encoding="utf-8",
        newline="",
    )
    source = {
        "schema": "degali.field-atmospheric-schedule-input.v1",
        "csv_path": csv_path.name,
        "evidence": {
            "dataset_id": "batch-source-A",
            "path": csv_path.name,
            "sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
            "row_count": 3,
            "source_boundary_id": "post-flash-A",
            "common_clock_id": "clock-A",
            "source_kind": "post_flash_atmospheric_vapour",
        },
    }
    path = tmp_path / "source-case.json"
    path.write_text(json.dumps(source, indent=2), encoding="utf-8")
    return path


def _write_batch(tmp_path: Path, *, source: bool = True) -> Path:
    case_path = _write_case(tmp_path)
    payload = {
        "schema": FIELD_BATCH_INPUT_SCHEMA,
        "options": {
            "include_refinement": False,
            "max_uncertainty_cases": 8,
            "table_nodes": 17,
        },
        "cases": [{
            "label": "source" if source else "nominal",
            "case_path": case_path.name,
            **({"atmospheric_source_case_path": _write_source(tmp_path).name} if source else {}),
        }],
    }
    path = tmp_path / "batch.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def test_batch_input_resolves_case_and_source_provenance(tmp_path):
    parsed = read_field_batch_input_json(_write_batch(tmp_path))

    assert parsed.cases[0].label == "source"
    assert parsed.cases[0].case_path.name == "field-case.json"
    assert parsed.cases[0].source_case_path.name == "source-case.json"
    assert parsed.cases[0].source_case_sha256 == hashlib.sha256(
        (tmp_path / "source-case.json").read_bytes()
    ).hexdigest()
    assert parsed.atmospheric_source_schedules["source"].evidence.dataset_id == "batch-source-A"
    assert parsed.as_record()["schema"] == FIELD_BATCH_INPUT_SCHEMA



def test_batch_input_rejects_unknown_keys_and_duplicate_labels(tmp_path):
    batch = json.loads(_write_batch(tmp_path).read_text(encoding="utf-8"))
    batch["unexpected"] = True
    with pytest.raises(ValueError, match="unknown=unexpected"):
        field_batch_input_from_mapping(batch, base_directory=tmp_path)

    batch.pop("unexpected")
    batch["cases"].append(dict(batch["cases"][0]))
    with pytest.raises(ValueError, match="duplicate"):
        field_batch_input_from_mapping(batch, base_directory=tmp_path)


def test_batch_input_json_rejects_normalized_duplicate_keys(tmp_path):
    source = tmp_path / "duplicate-keys.json"
    source.write_text(
        '{"schema": "degali.field-batch-input.v1", '
        '" SCHEMA ": "degali.field-batch-input.v1"}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate or normalization-colliding"):
        read_field_batch_input_json(source)


def test_batch_input_rejects_unsafe_labels_before_loading_artifacts(tmp_path):
    batch = {
        "schema": FIELD_BATCH_INPUT_SCHEMA,
        "cases": [{"label": "../escape", "case_path": "missing.json"}],
    }
    with pytest.raises(ValueError, match="safe report-name pattern"):
        field_batch_input_from_mapping(batch, base_directory=tmp_path)


def test_direct_batch_input_requires_immutable_case_collection(tmp_path):
    parsed = read_field_batch_input_json(_write_batch(tmp_path))
    with pytest.raises(TypeError, match="cases must be a tuple"):
        type(parsed)(list(parsed.cases), parsed.options)


def test_batch_input_rejects_historian_or_phase_routing_cases(tmp_path):
    case_path = tmp_path / "phase-case.json"
    case_path.write_text(
        (ROOT / "docs" / "field-phase-routing-case.example.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    batch = {
        "schema": FIELD_BATCH_INPUT_SCHEMA,
        "cases": [{"label": "history", "case_path": case_path.name}],
    }
    with pytest.raises(ValueError, match="historian or phase-routing"):
        from degali.addons.field_batch_io import field_batch_input_from_mapping
        field_batch_input_from_mapping(batch, base_directory=tmp_path)


def test_field_batch_cli_writes_execution_provenance_and_manifest(tmp_path):
    batch_path = _write_batch(tmp_path)
    output = tmp_path / "batch-output"

    status = main([
        "field-batch", str(batch_path),
        "--output-directory", str(output),
    ])

    assert status == 0
    execution = json.loads((output / "batch-execution.json").read_text(encoding="utf-8"))
    assert execution["schema"] == "degali.field-batch-execution.v1"
    assert execution["input"]["sha256"] == hashlib.sha256(batch_path.read_bytes()).hexdigest()
    assert execution["manifest_sha256"] == hashlib.sha256(
        (output / "manifest.json").read_bytes()
    ).hexdigest()
    assert execution["input"]["batch"]["cases"][0]["atmospheric_source_case"]["sha256"]
    assert execution["cases"][0]["report_sha256"] == hashlib.sha256(
        (output / "source.json").read_bytes()
    ).hexdigest()
    assert execution["cases"][0]["operational_decision"] == "withheld"
    assert execution["cases"][0]["gate_codes"]
    assert execution["summary"]["schema"] == "degali.field-batch-summary.v1"
    assert execution["summary"]["case_count"] == 1
    assert execution["summary"]["withheld_case_count"] == 1
    assert not execution["summary"]["all_screening_allowed"]
    assert (output / "manifest.json").is_file()
    assert (output / "source.json").is_file()


def test_field_batch_cli_applies_case_review_to_source_envelope(tmp_path):
    batch_path = _write_batch(tmp_path)
    batch = json.loads(batch_path.read_text(encoding="utf-8"))
    batch["options"].update({
        "include_refinement": True,
        "relative_tolerance": 1.0,
        "allow_conditional_operational_screening": True,
    })
    batch_path.write_text(json.dumps(batch, indent=2), encoding="utf-8")
    output = tmp_path / "conditional-source-output"

    assert main([
        "field-batch", str(batch_path),
        "--output-directory", str(output),
    ]) == 0

    execution = json.loads((output / "batch-execution.json").read_text(encoding="utf-8"))
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert execution["cases"][0]["operational_decision"] == "conditional_allowed"
    assert "conditional_review_required" in execution["cases"][0]["gate_codes"]
    assert manifest["cases"][0]["conditional_review"]["review_id"] == "field-review-A"
    assert manifest["cases"][0]["conditional_review"]["applied_to_this_execution"]


def test_field_batch_cli_require_screening_returns_two_for_withheld_batch(tmp_path):
    batch_path = _write_batch(tmp_path, source=False)
    output = tmp_path / "required-output"

    assert main([
        "field-batch", str(batch_path),
        "--output-directory", str(output), "--require-screening",
    ]) == 2
    assert (output / "batch-execution.json").is_file()
