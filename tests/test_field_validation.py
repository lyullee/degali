import hashlib
import json
from dataclasses import replace

import pytest

from degali.cli import main
from degali.addons.field_comparison import (
    FieldComparisonBasis,
    FieldModelSensorSet,
    FieldSensorPrediction,
    field_model_sensor_set_from_csv,
)
from degali.addons.field_contracts import FieldValidationEvidence
from degali.addons.field_evidence_audit import audit_field_evidence
from degali.addons.field_evidence_manifest import (
    FieldEvidenceManifest,
    write_field_evidence_manifest_json,
)
from degali.addons.field_validation import (
    FieldValidationDataset,
    FieldValidationObservation,
    field_validation_dataset_from_csv,
    field_validation_score_record,
    score_field_model_against_validation,
)
from degali.addons.field_validation_io import (
    FIELD_VALIDATION_INPUT_SCHEMA,
    FieldValidationManifestReference,
    FieldValidationCase,
    field_validation_case_from_mapping,
    read_field_validation_case_json,
    write_field_validation_case_json,
)


def _basis(obstacle="mask-v1"):
    return FieldComparisonBasis(
        "source-a", "weather-a", "sensors-a", "mean-60s", 60.0, obstacle,
    )


def _model(basis=None):
    return FieldModelSensorSet(
        model_id="DEGALI",
        temporal_mode="transient",
        basis=basis or _basis(),
        predictions=(
            FieldSensorPrediction("s1", (10.0, 0.0, 1.5), 0.05),
            FieldSensorPrediction("s2", (20.0, 0.0, 1.5), 0.04),
        ),
    )


def _csv(tmp_path, *, wrong_clock=False, obstacle="obstacle-a"):
    path = tmp_path / "observed.csv"
    clock = "clock-wrong" if wrong_clock else "clock-a"
    path.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        f"s1,10,0,1.5,0.04,2,60,{clock},{obstacle}\n"
        f"s2,20,0,1.5,0.02,2,60,clock-a,{obstacle}\n",
        encoding="utf-8",
    )
    return path


def _evidence(path, *, scope="lh2_obstacle_transport", obstacle="obstacle-a"):
    return FieldValidationEvidence(
        dataset_id="trial-a",
        path=str(path),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        row_count=2,
        source_boundary_id="source-a",
        weather_id="weather-a",
        obstacle_geometry_id=obstacle,
        receptor_geometry_id="sensors-a",
        temporal_operator_id="mean-60s",
        common_clock_id="clock-a",
        scope=scope,
    )


def test_validation_csv_is_fingerprinted_and_scored_against_matched_model(tmp_path):
    path = _csv(tmp_path)
    evidence = _evidence(path)
    dataset = field_validation_dataset_from_csv(
        path, evidence=evidence, basis=_basis(), temporal_mode="transient",
    )
    score = score_field_model_against_validation(_model(), dataset)

    assert score.status == "qualified"
    assert score.mean_absolute_error_mole_fraction == pytest.approx(0.015)
    assert score.root_mean_square_error_mole_fraction == pytest.approx((0.00025) ** 0.5)
    assert score.mean_bias_model_minus_observed == pytest.approx(0.015)
    assert score.comparison.classification_disagreement_count == 1
    record = field_validation_score_record(score)
    assert record["schema"] == "degali.field-validation-score.v1"
    assert record["dataset_id"] == "trial-a"
    assert "validation_qualified" in record["gate_codes"]
    assert record["comparison_gate_codes"] == ["matched_basis"]
    assert "no spatial interpolation" in record["qualification"]


def test_typed_validation_case_writer_pins_and_round_trips_both_csvs(tmp_path):
    observed_path = _csv(tmp_path)
    dataset = field_validation_dataset_from_csv(
        observed_path, evidence=_evidence(observed_path), basis=_basis(),
        temporal_mode="transient",
    )
    model_path = tmp_path / "model.csv"
    model_path.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "s1,10,0,1.5,60,0.05\n"
        "s2,20,0,1.5,60,0.04\n",
        encoding="utf-8",
    )
    model = field_model_sensor_set_from_csv(
        model_path,
        model_id="DEGALI",
        temporal_mode="transient",
        basis=_basis(),
        prediction_column="prediction",
        concentration_unit="mole_fraction",
    )
    case = FieldValidationCase(model=model, dataset=dataset)
    output = tmp_path / "exported" / "validation.json"

    exported = write_field_validation_case_json(case, output)

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["model"]["csv_path"] == "../model.csv"
    assert payload["dataset"]["csv_path"] == "../observed.csv"
    assert payload["validation_evidence"]["path"] == "../observed.csv"
    assert exported.model.csv_provenance is not None
    assert exported.model.csv_provenance.integrity_pinned is True
    assert exported.dataset.observations == dataset.observations
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_field_validation_case_json(case, output)


def _manifest_for_validation(evidence_root, observed_path):
    evidence_root.mkdir(parents=True, exist_ok=True)
    (evidence_root / "source.csv").write_text(
        "source_boundary_id,mass_flow_kg_s\nsource-a,0.1\n", encoding="utf-8",
    )
    (evidence_root / "weather.csv").write_text(
        "weather_id,wind_speed_m_s\nweather-a,2.0\n", encoding="utf-8",
    )
    (evidence_root / "obstacle.json").write_text(
        json.dumps({"obstacle_geometry_id": "obstacle-a", "vertices": []}),
        encoding="utf-8",
    )
    (evidence_root / "clock.csv").write_text(
        "common_clock_id\nclock-a\n", encoding="utf-8",
    )
    audit = audit_field_evidence(evidence_root)
    assert audit.status == "candidate_complete"
    paths = {
        "source_boundary": evidence_root / "source.csv",
        "weather": evidence_root / "weather.csv",
        "obstacle_geometry": evidence_root / "obstacle.json",
        "receptor_observations": observed_path,
        "common_clock": evidence_root / "clock.csv",
    }
    return FieldEvidenceManifest.from_audit(
        audit,
        manifest_id="manifest-validation-a",
        event_id="event-validation-a",
        selected_paths={channel: str(path) for channel, path in paths.items()},
        dataset_id="trial-a",
        observed_row_count=2,
        source_boundary_id="source-a",
        weather_id="weather-a",
        obstacle_geometry_id="obstacle-a",
        receptor_geometry_id="sensors-a",
        temporal_operator_id="mean-60s",
        common_clock_id="clock-a",
    )


def test_validation_case_round_trip_rechecks_optional_evidence_manifest(tmp_path):
    evidence_root = tmp_path / "evidence"
    evidence_root.mkdir(parents=True, exist_ok=True)
    observed_path = _csv(evidence_root)
    dataset = field_validation_dataset_from_csv(
        observed_path, evidence=_evidence(observed_path), basis=_basis(),
        temporal_mode="transient",
    )
    model_path = tmp_path / "model.csv"
    model_path.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "s1,10,0,1.5,60,0.05\n"
        "s2,20,0,1.5,60,0.04\n",
        encoding="utf-8",
    )
    model = field_model_sensor_set_from_csv(
        model_path, model_id="DEGALI", temporal_mode="transient", basis=_basis(),
        prediction_column="prediction", concentration_unit="mole_fraction",
    )
    manifest = _manifest_for_validation(evidence_root, observed_path)
    manifest_path = tmp_path / "manifest.json"
    write_field_evidence_manifest_json(manifest, manifest_path)
    reference = FieldValidationManifestReference(
        path=str(manifest_path.resolve()),
        sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        manifest=manifest,
    )
    case_path = tmp_path / "case" / "validation.json"
    exported = write_field_validation_case_json(
        FieldValidationCase(model=model, dataset=dataset, manifest_reference=reference),
        case_path,
    )
    payload = json.loads(case_path.read_text(encoding="utf-8"))
    assert payload["validation_manifest"]["path"] == "../manifest.json"
    assert exported.manifest_reference is not None
    assert exported.manifest_reference.manifest.event_id == "event-validation-a"

    payload["validation_evidence"]["weather_id"] = "weather-tampered"
    tampered = tmp_path / "case" / "tampered.json"
    tampered.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        read_field_validation_case_json(tampered)

    manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_payload["field_evidence_manifest"]["weather_id"] = "weather-tampered"
    manifest_tampered_path = tmp_path / "manifest-tampered.json"
    manifest_tampered_path.write_text(json.dumps(manifest_payload), encoding="utf-8")
    manifest_case = json.loads(case_path.read_text(encoding="utf-8"))
    manifest_case["validation_manifest"] = {
        "path": "../manifest-tampered.json",
        "sha256": hashlib.sha256(manifest_tampered_path.read_bytes()).hexdigest(),
    }
    manifest_case_path = tmp_path / "case" / "manifest-mismatch.json"
    manifest_case_path.write_text(json.dumps(manifest_case), encoding="utf-8")
    with pytest.raises(ValueError, match="manifest evidence does not match"):
        read_field_validation_case_json(manifest_case_path)


def test_validation_case_json_rejects_normalized_duplicate_keys(tmp_path):
    source = tmp_path / "duplicate-keys.json"
    source.write_text(
        '{"schema": "degali.field-validation-input.v1", '
        '" SCHEMA ": "degali.field-validation-input.v1"}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate or normalization-colliding"):
        read_field_validation_case_json(source)


def test_validation_score_contract_rejects_inconsistent_metrics_and_status(tmp_path):
    path = _csv(tmp_path)
    dataset = field_validation_dataset_from_csv(
        path, evidence=_evidence(path), basis=_basis(), temporal_mode="transient",
    )
    score = score_field_model_against_validation(_model(), dataset)

    with pytest.raises(ValueError, match="mean_absolute_error"):
        replace(score, mean_absolute_error_mole_fraction=0.0)
    with pytest.raises(ValueError, match="accepted comparison"):
        replace(score, reasons=("fabricated qualification",))
    with pytest.raises(TypeError, match="boolean"):
        FieldValidationObservation("s1", (10.0, 0.0, 1.5), True, 2.0)
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        score_field_model_against_validation(
            _model(), dataset, threshold_mole_fraction=True,
        )
    with pytest.raises(ValueError, match="position_tolerance_m"):
        score_field_model_against_validation(
            _model(), dataset, position_tolerance_m=True,
        )
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        FieldValidationCase(_model(), dataset, threshold_mole_fraction=True)
    with pytest.raises(ValueError, match="position_tolerance_m"):
        FieldValidationCase(_model(), dataset, position_tolerance_m=True)
    with pytest.raises(TypeError, match="FieldModelSensorSet"):
        FieldValidationCase("model", dataset)
    with pytest.raises(TypeError, match="FieldValidationDataset"):
        FieldValidationCase(_model(), "dataset")
    with pytest.raises(TypeError, match="reasons must be a tuple"):
        replace(score, reasons=["mutable reason"])
    with pytest.raises(ValueError, match="conditional validation score"):
        replace(score, status="conditional", reasons=())


def test_validation_dataset_direct_contract_matches_fingerprinted_csv_boundary(tmp_path):
    path = _csv(tmp_path)
    dataset = field_validation_dataset_from_csv(
        path, evidence=_evidence(path), basis=_basis(), temporal_mode="transient",
    )
    with pytest.raises(ValueError, match="observation count"):
        replace(dataset, observations=dataset.observations[:1])
    mismatched_provenance = replace(
        dataset.csv_provenance, path=str(tmp_path / "other.csv"),
    )
    with pytest.raises(ValueError, match="evidence path"):
        replace(dataset, csv_provenance=mismatched_provenance)
    with pytest.raises(TypeError, match="observations must be a tuple"):
        FieldValidationDataset(
            evidence=dataset.evidence, basis=dataset.basis,
            temporal_mode=dataset.temporal_mode,
            observations=list(dataset.observations),
            csv_provenance=dataset.csv_provenance,
        )


def test_validation_score_direct_contract_matches_both_comparison_sides(tmp_path):
    path = _csv(tmp_path)
    dataset = field_validation_dataset_from_csv(
        path, evidence=_evidence(path), basis=_basis(), temporal_mode="transient",
    )
    score = score_field_model_against_validation(_model(), dataset)
    fabricated_right = replace(score.comparison.right, model_id="observed:other-dataset")
    fabricated_comparison = replace(score.comparison, right=fabricated_right)
    with pytest.raises(ValueError, match="right model"):
        replace(score, comparison=fabricated_comparison)
    fabricated_left = replace(score.comparison.left, runtime_s=123.0)
    fabricated_comparison = replace(score.comparison, left=fabricated_left)
    with pytest.raises(ValueError, match="complete comparison left model"):
        replace(score, comparison=fabricated_comparison)


@pytest.mark.parametrize(
    ("basis_kwargs", "expected_label"),
    [
        ({"source_boundary_id": "source-other"}, "source boundary"),
        ({"weather_id": "weather-other"}, "weather"),
        ({"sensor_geometry_id": "sensors-other"}, "sensor geometry"),
        ({"temporal_operator_id": "mean-10s"}, "temporal operator"),
    ],
)
def test_validation_dataset_rejects_basis_evidence_identity_mismatch(
    tmp_path, basis_kwargs, expected_label,
):
    path = _csv(tmp_path)
    basis_values = {
        "source_boundary_id": "source-a",
        "weather_id": "weather-a",
        "sensor_geometry_id": "sensors-a",
        "temporal_operator_id": "mean-60s",
        "averaging_time_s": 60.0,
        "obstacle_representation_id": "mask-v1",
    }
    basis_values.update(basis_kwargs)
    with pytest.raises(ValueError, match=expected_label):
        field_validation_dataset_from_csv(
            path,
            evidence=_evidence(path),
            basis=FieldComparisonBasis(**basis_values),
            temporal_mode="transient",
        )


def test_validation_csv_rejects_a_mismatched_common_clock(tmp_path):
    path = _csv(tmp_path, wrong_clock=True)
    with pytest.raises(ValueError, match="common_clock_id"):
        field_validation_dataset_from_csv(
            path, evidence=_evidence(path), basis=_basis(), temporal_mode="transient",
        )


def test_validation_csv_rejects_ambiguous_headers_extra_values_and_role_collisions(tmp_path):
    duplicate = tmp_path / "duplicate-header.csv"
    duplicate.write_text(
        "sensor,sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        "s1,s1,10,0,1.5,0.04,2,60,clock-a,obstacle-a\n"
        "s2,s2,20,0,1.5,0.02,2,60,clock-a,obstacle-a\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate columns"):
        field_validation_dataset_from_csv(
            duplicate, evidence=_evidence(duplicate), basis=_basis(),
            temporal_mode="transient",
        )

    extra = tmp_path / "extra-value.csv"
    extra.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        "s1,10,0,1.5,0.04,2,60,clock-a,obstacle-a,unexpected\n"
        "s2,20,0,1.5,0.02,2,60,clock-a,obstacle-a,unexpected\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="more values than header"):
        field_validation_dataset_from_csv(
            extra, evidence=_evidence(extra), basis=_basis(),
            temporal_mode="transient",
        )

    missing = tmp_path / "missing-value.csv"
    missing.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        "s1,10,0,1.5,0.04,2,60,clock-a\n"
        "s2,20,0,1.5,0.02,2,60,clock-a,obstacle-a\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="fewer values than header"):
        field_validation_dataset_from_csv(
            missing, evidence=_evidence(missing), basis=_basis(),
            temporal_mode="transient",
        )

    normalized = tmp_path / "normalized-header.csv"
    normalized.write_text(
        "sensor,x_downwind_m, x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        "s1,10,10,0,1.5,0.04,2,60,clock-a,obstacle-a\n"
        "s2,20,20,0,1.5,0.02,2,60,clock-a,obstacle-a\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate columns"):
        field_validation_dataset_from_csv(
            normalized, evidence=_evidence(normalized), basis=_basis(),
            temporal_mode="transient",
        )

    path = _csv(tmp_path)
    with pytest.raises(ValueError, match="column roles must be distinct"):
        field_validation_dataset_from_csv(
            path, evidence=_evidence(path), basis=_basis(), temporal_mode="transient",
            sensor_id_column="x_downwind_m",
        )


def test_validation_csv_rejects_evidence_path_mismatch_even_when_digest_matches(tmp_path):
    path = _csv(tmp_path)
    evidence = _evidence(path)
    mismatched = FieldValidationEvidence(
        dataset_id=evidence.dataset_id,
        path=str(tmp_path / "copied-observed.csv"),
        sha256=evidence.sha256,
        row_count=evidence.row_count,
        source_boundary_id=evidence.source_boundary_id,
        weather_id=evidence.weather_id,
        obstacle_geometry_id=evidence.obstacle_geometry_id,
        receptor_geometry_id=evidence.receptor_geometry_id,
        temporal_operator_id=evidence.temporal_operator_id,
        common_clock_id=evidence.common_clock_id,
        scope=evidence.scope,
    )
    with pytest.raises(ValueError, match="evidence path"):
        field_validation_dataset_from_csv(
            path, evidence=mismatched, basis=_basis(), temporal_mode="transient",
        )


def test_free_field_validation_cannot_qualify_an_obstacle_model(tmp_path):
    path = _csv(tmp_path, obstacle="none")
    evidence = _evidence(path, scope="lh2_free_field", obstacle="none")
    dataset = field_validation_dataset_from_csv(
        path, evidence=evidence, basis=_basis(), temporal_mode="transient",
    )
    score = score_field_model_against_validation(_model(), dataset)

    assert score.status == "withheld"
    assert any("does not cover obstacle transport" in reason for reason in score.reasons)
    assert "obstacle_validation_unsupported" in score.gate_codes


def test_obstacle_validation_requires_a_model_obstacle_representation(tmp_path):
    path = _csv(tmp_path)
    evidence = _evidence(path)
    dataset = field_validation_dataset_from_csv(
        path, evidence=evidence, basis=_basis(), temporal_mode="transient",
    )
    score = score_field_model_against_validation(_model(_basis(obstacle=None)), dataset)

    assert score.status == "withheld"
    assert any("model obstacle_representation_id" in reason for reason in score.reasons)
    assert "obstacle_validation_unsupported" in score.gate_codes


def test_strict_validation_case_resolves_model_and_observation_csvs(tmp_path):
    observed = _csv(tmp_path)
    model_csv = tmp_path / "model.csv"
    model_csv.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "s1,10,0,1.5,60,0.05\n"
        "s2,20,0,1.5,60,0.04\n",
        encoding="utf-8",
    )
    evidence = _evidence(observed)
    case = field_validation_case_from_mapping(
        {
            "schema": FIELD_VALIDATION_INPUT_SCHEMA,
            "model": {
                "model_id": "DEGALI",
                "temporal_mode": "transient",
                "csv_path": "model.csv",
                "prediction_column": "prediction",
                "concentration_unit": "mole_fraction",
                "comparison_basis": {
                    "source_boundary_id": "source-a",
                    "weather_id": "weather-a",
                    "sensor_geometry_id": "sensors-a",
                    "temporal_operator_id": "mean-60s",
                    "averaging_time_s": 60.0,
                    "obstacle_representation_id": "mask-v1",
                },
            },
            "dataset": {
                "csv_path": "observed.csv",
                "temporal_mode": "transient",
                "comparison_basis": {
                    "source_boundary_id": "source-a",
                    "weather_id": "weather-a",
                    "sensor_geometry_id": "sensors-a",
                    "temporal_operator_id": "mean-60s",
                    "averaging_time_s": 60.0,
                    "obstacle_representation_id": "mask-v1",
                },
            },
            "validation_evidence": {
                "dataset_id": "trial-a",
                "path": "observed.csv",
                "sha256": evidence.sha256,
                "row_count": 2,
                "source_boundary_id": "source-a",
                "weather_id": "weather-a",
                "obstacle_geometry_id": "obstacle-a",
                "receptor_geometry_id": "sensors-a",
                "temporal_operator_id": "mean-60s",
                "common_clock_id": "clock-a",
                "scope": "lh2_obstacle_transport",
            },
        },
        base_directory=tmp_path,
    )

    score = score_field_model_against_validation(
        case.model, case.dataset,
        threshold_mole_fraction=case.threshold_mole_fraction,
        position_tolerance_m=case.position_tolerance_m,
    )
    assert score.status == "qualified"


def test_field_validate_cli_emits_a_qualified_record(tmp_path, capsys):
    observed = _csv(tmp_path)
    model_csv = tmp_path / "model.csv"
    model_csv.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "s1,10,0,1.5,60,0.05\n"
        "s2,20,0,1.5,60,0.04\n",
        encoding="utf-8",
    )
    evidence = _evidence(observed)
    case_path = tmp_path / "validation-case.json"
    case_path.write_text(json.dumps({
        "schema": FIELD_VALIDATION_INPUT_SCHEMA,
        "model": {
            "model_id": "DEGALI", "temporal_mode": "transient",
            "csv_path": "model.csv", "prediction_column": "prediction",
            "concentration_unit": "mole_fraction",
            "comparison_basis": {
                "source_boundary_id": "source-a", "weather_id": "weather-a",
                "sensor_geometry_id": "sensors-a", "temporal_operator_id": "mean-60s",
                "averaging_time_s": 60.0, "obstacle_representation_id": "mask-v1",
            },
        },
        "dataset": {
            "csv_path": "observed.csv", "temporal_mode": "transient",
            "comparison_basis": {
                "source_boundary_id": "source-a", "weather_id": "weather-a",
                "sensor_geometry_id": "sensors-a", "temporal_operator_id": "mean-60s",
                "averaging_time_s": 60.0, "obstacle_representation_id": "mask-v1",
            },
        },
        "validation_evidence": {
            "dataset_id": "trial-a", "path": "observed.csv",
            "sha256": evidence.sha256, "row_count": 2,
            "source_boundary_id": "source-a", "weather_id": "weather-a",
            "obstacle_geometry_id": "obstacle-a", "receptor_geometry_id": "sensors-a",
            "temporal_operator_id": "mean-60s", "common_clock_id": "clock-a",
            "scope": "lh2_obstacle_transport",
        },
    }, ensure_ascii=False), encoding="utf-8")

    assert main(["field-validate", str(case_path), "--require-qualified"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["field_validation_score"]["status"] == "qualified"


def _censored_csv(tmp_path, *, first_kind="lower_bound", first_value=0.04, second_kind="lower_bound", second_value=0.02):
    path = tmp_path / "censored-observed.csv"
    path.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_kind,observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        f"s1,10,0,1.5,{first_value},{first_kind},2,60,clock-a,obstacle-a\n"
        f"s2,20,0,1.5,{second_value},{second_kind},2,60,clock-a,obstacle-a\n",
        encoding="utf-8",
    )
    return path


def test_lower_bound_observations_use_one_sided_metrics_and_withhold_qualification(tmp_path):
    path = _censored_csv(tmp_path)
    dataset = field_validation_dataset_from_csv(
        path,
        evidence=_evidence(path),
        basis=_basis(),
        temporal_mode="transient",
        observation_kind_column="observation_kind",
    )
    score = score_field_model_against_validation(_model(), dataset)

    assert score.status == "conditional"
    assert score.mean_absolute_error_mole_fraction is None
    assert score.lower_bound_observation_count == 2
    assert score.lower_bound_satisfied_count == 2
    assert score.lower_bound_satisfaction_fraction == pytest.approx(1.0)
    assert score.maximum_lower_bound_deficit_mole_fraction == pytest.approx(0.0)
    record = field_validation_score_record(score)
    assert record["observation_counts"] == {"exact": 0, "lower_bound": 2}
    assert record["lower_bound_constraint"]["status"] == "satisfied"
    assert "lower_bound_present" in score.gate_codes


def test_lower_bound_violation_withholds_score_and_reports_deficit(tmp_path):
    path = _censored_csv(tmp_path, first_value=0.06)
    dataset = field_validation_dataset_from_csv(
        path,
        evidence=_evidence(path),
        basis=_basis(),
        temporal_mode="transient",
        observation_kind_column="observation_kind",
    )
    score = score_field_model_against_validation(_model(), dataset)

    assert score.status == "withheld"
    assert score.lower_bound_constraint_status == "violated"
    assert score.lower_bound_satisfied_count == 1
    assert score.lower_bound_satisfaction_fraction == pytest.approx(0.5)
    assert score.maximum_lower_bound_deficit_mole_fraction == pytest.approx(0.01)
    assert any("violates lower-bound" in reason for reason in score.reasons)
    assert "lower_bound_violated" in score.gate_codes


def test_validation_score_gate_code_contract_rejects_unknown_or_duplicate_codes(tmp_path):
    path = _csv(tmp_path)
    dataset = field_validation_dataset_from_csv(
        path, evidence=_evidence(path), basis=_basis(), temporal_mode="transient",
    )
    score = score_field_model_against_validation(_model(), dataset)

    with pytest.raises(ValueError, match="unsupported code"):
        replace(score, gate_codes=("not-a-code",))
    with pytest.raises(ValueError, match="unique"):
        replace(score, gate_codes=("matched_basis", "matched_basis"))


def test_mixed_exact_and_lower_bound_rows_keep_symmetric_metrics_exact_only(tmp_path):
    path = _censored_csv(tmp_path, first_kind="exact", first_value=0.04, second_value=0.03)
    dataset = field_validation_dataset_from_csv(
        path,
        evidence=_evidence(path),
        basis=_basis(),
        temporal_mode="transient",
        observation_kind_column="observation_kind",
    )
    score = score_field_model_against_validation(_model(), dataset)

    assert score.status == "conditional"
    assert score.mean_absolute_error_mole_fraction == pytest.approx(0.01)
    assert score.root_mean_square_error_mole_fraction == pytest.approx(0.01)
    assert score.mean_bias_model_minus_observed == pytest.approx(0.01)
    assert score.lower_bound_satisfaction_fraction == pytest.approx(1.0)


def test_validation_observation_kind_is_strict(tmp_path):
    path = _censored_csv(tmp_path, first_kind="threshold", second_kind="lower_bound")
    with pytest.raises(ValueError, match="observation_kind must be exact or lower_bound"):
        field_validation_dataset_from_csv(
            path,
            evidence=_evidence(path),
            basis=_basis(),
            temporal_mode="transient",
            observation_kind_column="observation_kind",
        )


def test_strict_validation_case_accepts_censored_observation_column(tmp_path):
    observed = _censored_csv(tmp_path)
    model_csv = tmp_path / "model.csv"
    model_csv.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "s1,10,0,1.5,60,0.05\n"
        "s2,20,0,1.5,60,0.04\n",
        encoding="utf-8",
    )
    evidence = _evidence(observed)
    basis = {
        "source_boundary_id": "source-a",
        "weather_id": "weather-a",
        "sensor_geometry_id": "sensors-a",
        "temporal_operator_id": "mean-60s",
        "averaging_time_s": 60.0,
        "obstacle_representation_id": "mask-v1",
    }
    case = field_validation_case_from_mapping(
        {
            "schema": FIELD_VALIDATION_INPUT_SCHEMA,
            "model": {
                "model_id": "DEGALI",
                "temporal_mode": "transient",
                "csv_path": "model.csv",
                "prediction_column": "prediction",
                "concentration_unit": "mole_fraction",
                "comparison_basis": basis,
            },
            "dataset": {
                "csv_path": "censored-observed.csv",
                "temporal_mode": "transient",
                "observation_kind_column": "observation_kind",
                "comparison_basis": basis,
            },
            "validation_evidence": {
                "dataset_id": "trial-censored",
                "path": "censored-observed.csv",
                "sha256": evidence.sha256,
                "row_count": 2,
                "source_boundary_id": "source-a",
                "weather_id": "weather-a",
                "obstacle_geometry_id": "obstacle-a",
                "receptor_geometry_id": "sensors-a",
                "temporal_operator_id": "mean-60s",
                "common_clock_id": "clock-a",
                "scope": "lh2_obstacle_transport",
            },
        },
        base_directory=tmp_path,
    )
    assert all(item.is_lower_bound for item in case.dataset.observations)
    assert case.dataset.csv_provenance.observation_kind_column == "observation_kind"
