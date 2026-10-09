import json

import pytest

from degali.cli import main
from degali.addons.field_evidence_audit import (
    audit_field_evidence,
    write_field_evidence_audit_json,
)
from degali.addons.field_evidence_manifest import (
    FIELD_EVIDENCE_MANIFEST_SCHEMA,
    FieldEvidenceManifest,
    field_evidence_manifest_from_mapping,
    read_field_evidence_manifest_json,
    write_field_evidence_manifest_json,
)


def _write_channels(root):
    root.mkdir(parents=True, exist_ok=True)
    source = root / "source.csv"
    source.write_text(
        "source_boundary_id,mass_flow_kg_s\nsource-a,0.1\n", encoding="utf-8",
    )
    weather = root / "weather.csv"
    weather.write_text(
        "weather_id,wind_speed_m_s\nweather-a,2.0\n", encoding="utf-8",
    )
    obstacle = root / "obstacle.json"
    obstacle.write_text(
        json.dumps({"obstacle_geometry_id": "obstacle-a", "vertices": []}),
        encoding="utf-8",
    )
    receptors = root / "receptors.csv"
    receptors.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,common_clock_id\n"
        "s1,10,0,1.5,0.04,clock-a\n",
        encoding="utf-8",
    )
    clock = root / "clock.csv"
    clock.write_text("common_clock_id\nclock-a\n", encoding="utf-8")
    return {
        "source_boundary": source,
        "weather": weather,
        "obstacle_geometry": obstacle,
        "receptor_observations": receptors,
        "common_clock": clock,
    }


def _manifest(root):
    paths = _write_channels(root)
    audit = audit_field_evidence(root)
    assert audit.status == "candidate_complete"
    return FieldEvidenceManifest.from_audit(
        audit,
        manifest_id="manifest-a",
        event_id="event-a",
        selected_paths={channel: str(path) for channel, path in paths.items()},
        dataset_id="dataset-a",
        observed_row_count=1,
        source_boundary_id="source-a",
        weather_id="weather-a",
        obstacle_geometry_id="obstacle-a",
        receptor_geometry_id="receptors-a",
        temporal_operator_id="mean-1s",
        common_clock_id="clock-a",
    )


def test_manifest_requires_explicit_cross_file_join_and_converts_after_hash_check(tmp_path):
    manifest = _manifest(tmp_path)

    assert manifest.as_record()["schema"] == FIELD_EVIDENCE_MANIFEST_SCHEMA
    assert manifest.as_record()["promotion_allowed"] is False
    assert manifest.readiness_record()["status"] == "conditional"
    assert manifest.readiness_record()["missing_requirements"] == [
        "sensor_set_id", "operator_id",
    ]
    evidence = manifest.as_validation_evidence()
    assert evidence.dataset_id == "dataset-a"
    assert evidence.common_clock_id == "clock-a"
    assert evidence.supports_obstacle_transport


def test_manifest_writer_and_reader_pin_selected_files(tmp_path, capsys):
    manifest = _manifest(tmp_path / "evidence")
    output = tmp_path / "manifest-execution.json"

    saved = write_field_evidence_manifest_json(manifest, output)
    loaded = read_field_evidence_manifest_json(output)

    assert saved == manifest
    assert loaded == manifest
    with pytest.raises(FileExistsError):
        write_field_evidence_manifest_json(manifest, output)

    receptor = tmp_path / "evidence" / "receptors.csv"
    receptor.write_text(receptor.read_text(encoding="utf-8") + "s2,11,0,1.5,0.02,clock-a\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed since its provenance"):
        read_field_evidence_manifest_json(output)

    assert main(["field-evidence-manifest-verify", str(output)]) == 2
    withheld = json.loads(capsys.readouterr().out)
    assert withheld["evidence_readiness"]["status"] == "withheld"
    assert withheld["verification"]["selected_files_verified"] is False


def test_manifest_explicit_sensor_and_operator_metadata_is_accepted(tmp_path):
    paths = _write_channels(tmp_path / "evidence")
    registry = tmp_path / "evidence" / "sensor_registry.csv"
    registry.write_text(
        "sensor_set_id,sensor_id,calibration_id,valid_from,valid_to\n"
        "sensor-set-a,s1,cal-a,2026-09-21T00:00:00Z,2026-09-22T00:00:00Z\n",
        encoding="utf-8",
    )
    audit = audit_field_evidence(tmp_path / "evidence")
    manifest = FieldEvidenceManifest.from_audit(
        audit,
        manifest_id="manifest-complete",
        event_id="event-complete",
        selected_paths={channel: str(path) for channel, path in paths.items()},
        dataset_id="dataset-complete",
        observed_row_count=1,
        source_boundary_id="source-a",
        weather_id="weather-a",
        obstacle_geometry_id="obstacle-a",
        receptor_geometry_id="receptors-a",
        temporal_operator_id="mean-1s",
        common_clock_id="clock-a",
        sensor_set_id="sensor-set-a",
        operator_id="operator-a",
        sensor_registry_path=str(registry),
        sensor_calibration_status="certified",
    )
    assert manifest.readiness_record() == {
        "status": "accepted",
        "missing_requirements": [],
        "promotion_allowed": False,
        "gate_codes": ["manifest_metadata_complete", "manifest_not_promoted"],
    }
    record = manifest.as_record()
    assert record["sensor_set_id"] == "sensor-set-a"
    assert record["operator_id"] == "operator-a"
    assert record["sensor_calibration_status"] == "certified"
    assert record["evidence_readiness"]["status"] == "accepted"


def test_manifest_nominal_sensor_specification_remains_conditional(tmp_path):
    paths = _write_channels(tmp_path / "evidence")
    registry = tmp_path / "evidence" / "sensor_registry.csv"
    registry.write_text(
        "sensor_set_id,sensor_id,calibration_status\n"
        "sensor-set-a,s1,specification_only\n",
        encoding="utf-8",
    )
    audit = audit_field_evidence(tmp_path / "evidence")
    manifest = FieldEvidenceManifest.from_audit(
        audit,
        manifest_id="manifest-spec-only",
        event_id="event-spec-only",
        selected_paths={channel: str(path) for channel, path in paths.items()},
        dataset_id="dataset-spec-only",
        observed_row_count=1,
        source_boundary_id="source-a",
        weather_id="weather-a",
        obstacle_geometry_id="obstacle-a",
        receptor_geometry_id="receptors-a",
        temporal_operator_id="mean-1s",
        common_clock_id="clock-a",
        sensor_set_id="sensor-set-a",
        operator_id="operator-a",
        sensor_registry_path=str(registry),
        sensor_calibration_status="specification_only",
    )
    readiness = manifest.readiness_record()
    assert readiness["status"] == "conditional"
    assert readiness["missing_requirements"] == ["sensor_calibration_certificate"]


def test_manifest_rejects_tampered_readiness_record(tmp_path):
    manifest = _manifest(tmp_path / "evidence")
    payload = manifest.as_record()
    payload["evidence_readiness"] = {
        "status": "accepted",
        "missing_requirements": [],
        "promotion_allowed": False,
        "gate_codes": ["manifest_metadata_complete", "manifest_not_promoted"],
    }
    with pytest.raises(ValueError, match="evidence_readiness"):
        field_evidence_manifest_from_mapping(payload)


def test_manifest_rejects_promotion_or_channel_path_drift(tmp_path):
    manifest = _manifest(tmp_path / "evidence")
    payload = manifest.as_record()
    payload["promotion_allowed"] = True
    with pytest.raises(ValueError, match="promotion_allowed"):
        field_evidence_manifest_from_mapping(payload)

    payload = manifest.as_record()
    payload["channel_artifacts"][0]["path"] = str(tmp_path / "outside.csv")
    drifted = field_evidence_manifest_from_mapping(payload)
    with pytest.raises(ValueError, match="escapes audit_root|missing"):
        drifted.verify_files()


def test_manifest_rejects_incomplete_audit_and_wrong_observed_digest(tmp_path):
    paths = _write_channels(tmp_path)
    audit = audit_field_evidence(tmp_path)
    # Remove the explicit clock channel from the selection; no implicit join is allowed.
    with pytest.raises(ValueError, match="exactly the five audit channels"):
        FieldEvidenceManifest.from_audit(
            audit,
            manifest_id="manifest-a",
            event_id="event-a",
            selected_paths={channel: str(path) for channel, path in paths.items() if channel != "common_clock"},
            dataset_id="dataset-a",
            observed_row_count=1,
            source_boundary_id="source-a",
            weather_id="weather-a",
            obstacle_geometry_id="obstacle-a",
            receptor_geometry_id="receptors-a",
            temporal_operator_id="mean-1s",
            common_clock_id="clock-a",
        )


def test_manifest_rejects_partial_audit_before_channel_selection(tmp_path):
    paths = _write_channels(tmp_path)
    paths["source_boundary"].unlink()
    audit = audit_field_evidence(tmp_path)
    assert audit.status == "partial"
    with pytest.raises(ValueError, match="candidate_complete audit"):
        FieldEvidenceManifest.from_audit(
            audit,
            manifest_id="manifest-partial",
            event_id="event-partial",
            selected_paths={channel: str(path) for channel, path in paths.items()},
            dataset_id="dataset-partial",
            observed_row_count=1,
            source_boundary_id="source-a",
            weather_id="weather-a",
            obstacle_geometry_id="obstacle-a",
            receptor_geometry_id="receptors-a",
            temporal_operator_id="mean-1s",
            common_clock_id="clock-a",
        )


def test_manifest_verify_cli_rechecks_and_withholds_promotion(tmp_path, capsys):
    manifest = _manifest(tmp_path / "evidence")
    output = tmp_path / "manifest-execution.json"
    write_field_evidence_manifest_json(manifest, output)

    assert main(["field-evidence-manifest-verify", str(output)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema"] == "degali.field-evidence-manifest-verification.v1"
    assert payload["evidence_readiness"]["status"] == "conditional"
    assert payload["verification"] == {
        "selected_files_verified": True,
        "promotion_allowed": False,
    }


def test_manifest_create_cli_pins_sensor_registry_and_reports_accepted(tmp_path, capsys):
    root = tmp_path / "evidence"
    paths = _write_channels(root)
    registry = root / "sensor_registry.csv"
    registry.write_text(
        "sensor_set_id,sensor_id,calibration_id\nsensor-set-cli,s1,cal-cli\n",
        encoding="utf-8",
    )
    audit_path = tmp_path / "audit-execution.json"
    write_field_evidence_audit_json(audit_field_evidence(root), audit_path)
    output = tmp_path / "manifest-execution.json"
    arguments = [
        "field-evidence-manifest-create", str(audit_path),
        "--manifest-id", "manifest-cli-complete",
        "--event-id", "event-cli-complete",
        "--dataset-id", "dataset-cli-complete",
        "--observed-row-count", "1",
        "--source-boundary-id", "source-a",
        "--weather-id", "weather-a",
        "--obstacle-geometry-id", "obstacle-a",
        "--receptor-geometry-id", "receptors-a",
        "--temporal-operator-id", "mean-1s",
        "--common-clock-id", "clock-a",
        "--sensor-set-id", "sensor-set-cli",
        "--operator-id", "operator-cli",
        "--sensor-registry-path", str(registry),
        "--sensor-calibration-status", "certified",
        "--output", str(output),
    ]
    for channel, path in paths.items():
        arguments.extend(["--selected-path", f"{channel}={path}"])

    assert main(arguments) == 0
    assert "evidence readiness: accepted" in capsys.readouterr().out
    loaded = read_field_evidence_manifest_json(output)
    assert loaded.sensor_registry_artifact is not None
    assert loaded.readiness_record()["status"] == "accepted"


def test_manifest_create_cli_requires_explicit_five_channel_selection(tmp_path, capsys):
    root = tmp_path / "evidence"
    paths = _write_channels(root)
    audit_path = tmp_path / "audit-execution.json"
    write_field_evidence_audit_json(audit_field_evidence(root), audit_path)
    output = tmp_path / "manifest-execution.json"

    arguments = [
        "field-evidence-manifest-create", str(audit_path),
        "--manifest-id", "manifest-cli",
        "--event-id", "event-cli",
        "--dataset-id", "dataset-cli",
        "--observed-row-count", "1",
        "--source-boundary-id", "source-a",
        "--weather-id", "weather-a",
        "--obstacle-geometry-id", "obstacle-a",
        "--receptor-geometry-id", "receptors-a",
        "--temporal-operator-id", "mean-1s",
        "--common-clock-id", "clock-a",
        "--output", str(output),
    ]
    for channel, path in paths.items():
        arguments.extend(["--selected-path", f"{channel}={path}"])

    assert main(arguments) == 0
    assert "promotion allowed: false" in capsys.readouterr().out
    loaded = read_field_evidence_manifest_json(output)
    assert loaded.manifest_id == "manifest-cli"
    assert loaded.event_id == "event-cli"

    assert main(arguments) == 1


def test_manifest_create_cli_refuses_partial_audit_without_output(tmp_path):
    root = tmp_path / "partial-evidence"
    paths = _write_channels(root)
    paths["common_clock"].unlink()
    audit_path = tmp_path / "partial-audit-execution.json"
    write_field_evidence_audit_json(audit_field_evidence(root), audit_path)
    output = tmp_path / "must-not-exist.json"
    arguments = [
        "field-evidence-manifest-create", str(audit_path),
        "--manifest-id", "manifest-partial",
        "--event-id", "event-partial",
        "--dataset-id", "dataset-partial",
        "--observed-row-count", "1",
        "--source-boundary-id", "source-a",
        "--weather-id", "weather-a",
        "--obstacle-geometry-id", "obstacle-a",
        "--receptor-geometry-id", "receptors-a",
        "--temporal-operator-id", "mean-1s",
        "--common-clock-id", "clock-a",
        "--output", str(output),
    ]
    for channel, path in paths.items():
        arguments.extend(["--selected-path", f"{channel}={path}"])

    assert main(arguments) == 1
    assert not output.exists()
