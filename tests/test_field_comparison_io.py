import json
import hashlib

import pytest

from degali.addons.field_comparison import (
    compare_field_model_sensor_sets,
    field_model_decision_impact,
    field_model_comparison_report,
)
from degali.addons.field_comparison_io import (
    FIELD_MODEL_COMPARISON_INPUT_SCHEMA,
    field_model_comparison_case_from_mapping,
    write_field_model_comparison_case_json,
    read_field_model_comparison_case_json,
)
from degali.cli import main


def _model(path, model_id, *, temporal_mode="transient", operator="mean-60s", obstacle_id=None):
    basis = {
        "source_boundary_id": "source-event-17",
        "weather_id": "met-window-17",
        "sensor_geometry_id": "sensor-layout-17",
        "temporal_operator_id": operator,
        "averaging_time_s": 60.0,
    }
    if obstacle_id is not None:
        basis["obstacle_representation_id"] = obstacle_id
    return {
        "model_id": model_id,
        "temporal_mode": temporal_mode,
        "csv_path": path,
        "prediction_column": "h2",
        "concentration_unit": "mole_fraction",
        "runtime_s": 0.1,
        "comparison_basis": basis,
    }


def _case(*, right_mode="transient", right_operator="mean-60s", left_obstacle=None, right_obstacle=None):
    return {
        "schema": FIELD_MODEL_COMPARISON_INPUT_SCHEMA,
        "threshold_mole_fraction": 0.04,
        "left": _model("degali.csv", "DEGALI", obstacle_id=left_obstacle),
        "right": _model(
            "slabx.csv", "SLABx", temporal_mode=right_mode, operator=right_operator,
            obstacle_id=right_obstacle,
        ),
    }


def _write_csvs(tmp_path):
    header = "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,h2\n"
    (tmp_path / "degali.csv").write_text(
        header + "S1,10,0,1.5,60,0.03\nS2,20,0,1.5,60,0.05\n",
        encoding="utf-8",
    )
    (tmp_path / "slabx.csv").write_text(
        header + "S1,10,0,1.5,60,0.05\nS2,20,0,1.5,60,0.05\n",
        encoding="utf-8",
    )


def test_json_comparison_case_resolves_relative_csvs_and_preserves_provenance(tmp_path):
    _write_csvs(tmp_path)
    source = tmp_path / "comparison.json"
    source.write_text(json.dumps(_case()), encoding="utf-8")

    parsed = read_field_model_comparison_case_json(source)
    comparison = compare_field_model_sensor_sets(
        parsed.left,
        parsed.right,
        threshold_mole_fraction=parsed.threshold_mole_fraction,
        position_tolerance_m=parsed.position_tolerance_m,
    )

    assert parsed.left.csv_provenance is not None
    assert parsed.left.csv_provenance.path == str((tmp_path / "degali.csv").resolve())
    assert comparison.applicability.status == "accepted"
    assert comparison.classification_disagreement_count == 1


def test_json_comparison_case_can_pin_csv_hash_and_row_count(tmp_path):
    _write_csvs(tmp_path)
    case = _case()
    for label in ("left", "right"):
        artifact = tmp_path / case[label]["csv_path"]
        case[label]["csv_sha256"] = hashlib.sha256(artifact.read_bytes()).hexdigest()
        case[label]["csv_row_count"] = 2
    source = tmp_path / "comparison-pinned.json"
    source.write_text(json.dumps(case), encoding="utf-8")

    parsed = read_field_model_comparison_case_json(source)

    assert parsed.left.csv_provenance is not None
    assert parsed.left.csv_provenance.row_count == 2
    assert parsed.left.csv_provenance.integrity_pinned is True
    report = field_model_comparison_report(
        compare_field_model_sensor_sets(parsed.left, parsed.right)
    )
    assert report["left"]["csv_provenance"]["integrity_pinned"] is True
    (tmp_path / "degali.csv").write_text(
        (tmp_path / "degali.csv").read_text(encoding="utf-8")
        + "S3,30,0,1.5,60,0.01\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="SHA-256"):
        read_field_model_comparison_case_json(source)


def test_json_comparison_case_rejects_unpaired_csv_fingerprint_fields(tmp_path):
    _write_csvs(tmp_path)
    case = _case()
    case["left"]["csv_row_count"] = 2
    source = tmp_path / "comparison-unpaired.json"
    source.write_text(json.dumps(case), encoding="utf-8")

    with pytest.raises(ValueError, match="must be supplied together"):
        read_field_model_comparison_case_json(source)


def test_strict_comparison_json_rejects_normalized_duplicate_keys(tmp_path):
    source = tmp_path / "duplicate-keys.json"
    source.write_text(
        '{"schema": "degali.field-model-comparison-input.v1", '
        '" SCHEMA ": "degali.field-model-comparison-input.v1"}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate or normalization-colliding"):
        read_field_model_comparison_case_json(source)


def test_typed_comparison_case_writer_pins_and_round_trips_csv_artifacts(tmp_path):
    _write_csvs(tmp_path)
    source = tmp_path / "comparison.json"
    source.write_text(json.dumps(_case()), encoding="utf-8")
    parsed = read_field_model_comparison_case_json(source)
    exported_path = tmp_path / "exported" / "comparison.json"

    exported = write_field_model_comparison_case_json(parsed, exported_path)

    payload = json.loads(exported_path.read_text(encoding="utf-8"))
    assert payload["left"]["csv_path"] == "../degali.csv"
    assert len(payload["left"]["csv_sha256"]) == 64
    assert exported.left.csv_provenance is not None
    assert exported.left.csv_provenance.integrity_pinned is True
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_field_model_comparison_case_json(parsed, exported_path)


def test_in_memory_comparison_mapping_refuses_unrooted_csv_paths():
    with pytest.raises(ValueError, match="require a case file path"):
        field_model_comparison_case_from_mapping(_case())


def test_comparison_case_parses_obstacle_representation_basis_and_withholds_mismatch(tmp_path):
    _write_csvs(tmp_path)
    source = tmp_path / "comparison.json"
    source.write_text(
        json.dumps(_case(left_obstacle="mask-v1", right_obstacle="mask-v2")),
        encoding="utf-8",
    )

    parsed = read_field_model_comparison_case_json(source)
    comparison = compare_field_model_sensor_sets(parsed.left, parsed.right)

    assert parsed.left.basis.obstacle_representation_id == "mask-v1"
    assert parsed.right.basis.obstacle_representation_id == "mask-v2"
    assert comparison.applicability.status == "conditional"
    assert field_model_decision_impact(comparison).status == "withheld"


def test_comparison_case_fingerprints_qualification_manifest_and_rejects_drift(tmp_path):
    _write_csvs(tmp_path)
    manifest = tmp_path / "comparison-manifest.json"
    manifest.write_text('{"window":"60 s","scope":"conditional"}', encoding="utf-8")
    case = _case()
    case["comparison_evidence"] = {
        "manifest_path": manifest.name,
        "sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "qualification": "The manifest records a conditional source mapping.",
    }
    source = tmp_path / "comparison.json"
    source.write_text(json.dumps(case), encoding="utf-8")

    parsed = read_field_model_comparison_case_json(source)

    assert parsed.comparison_evidence is not None
    assert parsed.comparison_evidence.path == str(manifest.resolve())
    assert parsed.comparison_evidence.sha256 == case["comparison_evidence"]["sha256"]

    exported = write_field_model_comparison_case_json(
        parsed, tmp_path / "exported-with-evidence.json",
    )
    assert exported.comparison_evidence is not None
    exported_payload = json.loads(
        (tmp_path / "exported-with-evidence.json").read_text(encoding="utf-8")
    )
    assert exported_payload["comparison_evidence"]["manifest_path"] == manifest.name

    manifest.write_text('{"window":"61 s"}', encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256"):
        read_field_model_comparison_case_json(source)


def test_field_compare_cli_emits_auditable_record_and_fails_safe_for_mismatch(tmp_path, capsys):
    _write_csvs(tmp_path)
    source = tmp_path / "comparison.json"
    output = tmp_path / "comparison-output.json"
    source.write_text(json.dumps(_case()), encoding="utf-8")

    code = main(["field-compare", str(source), "--output", str(output), "--require-comparable"])
    captured = capsys.readouterr()
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert code == 0
    assert "comparison applicability: accepted" in captured.out
    assert payload["schema"] == "degali.field-model-comparison-execution.v1"
    report = payload["field_model_comparison"]
    assert report["left"]["csv_provenance"]["row_count"] == 2
    assert report["comparison_evidence"] is None
    assert report["model_selection_impact"]["status"] == "different"
    assert main(["field-compare", str(source), "--output", str(output)]) == 1
    assert "refusing to overwrite" in capsys.readouterr().err

    source.write_text(json.dumps(_case(right_mode="steady")), encoding="utf-8")
    assert main(["field-compare", str(source), "--require-comparable"]) == 2
    conditional = json.loads(capsys.readouterr().out)
    assert conditional["field_model_comparison"]["comparison"]["applicability"]["status"] == "conditional"
