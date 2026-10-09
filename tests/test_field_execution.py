import hashlib
import json
from pathlib import Path

import pytest

from degali.cli import main


ROOT = Path(__file__).resolve().parents[1]


def _source_case(tmp_path, *, duration=2.0):
    csv_path = tmp_path / "source.csv"
    csv_path.write_text(
        f"time_s,rate_kg_s\n0,0.2\n{duration / 2:g},0.1\n{duration:g},0\n",
        encoding="utf-8",
    )
    case_path = tmp_path / "source-case.json"
    case_path.write_text(json.dumps({
        "schema": "degali.field-atmospheric-schedule-input.v1",
        "csv_path": csv_path.name,
        "evidence": {
            "dataset_id": "source-a",
            "path": csv_path.name,
            "sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
            "row_count": 3,
            "source_boundary_id": "boundary-a",
            "common_clock_id": "clock-a",
            "source_kind": "post_flash_atmospheric_vapour",
        },
    }), encoding="utf-8")
    return case_path, csv_path


def _comparison_case(tmp_path):
    header = "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
    left = tmp_path / "left.csv"
    right = tmp_path / "right.csv"
    left.write_text(header + "s1,10,0,1.5,60,0.03\ns2,20,0,1.5,60,0.05\n", encoding="utf-8")
    right.write_text(header + "s1,10,0,1.5,60,0.04\ns2,20,0,1.5,60,0.05\n", encoding="utf-8")
    basis = {
        "source_boundary_id": "source-a",
        "weather_id": "weather-a",
        "sensor_geometry_id": "sensors-a",
        "temporal_operator_id": "mean-60s",
        "averaging_time_s": 60.0,
    }
    case_path = tmp_path / "comparison-case.json"
    case_path.write_text(json.dumps({
        "schema": "degali.field-model-comparison-input.v1",
        "left": {
            "model_id": "left", "temporal_mode": "steady", "csv_path": left.name,
            "prediction_column": "prediction", "concentration_unit": "mole_fraction",
            "comparison_basis": basis,
        },
        "right": {
            "model_id": "right", "temporal_mode": "steady", "csv_path": right.name,
            "prediction_column": "prediction", "concentration_unit": "mole_fraction",
            "comparison_basis": basis,
        },
    }), encoding="utf-8")
    return case_path, left


def _validation_case(tmp_path):
    observed = tmp_path / "observed.csv"
    observed.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,"
        "observation_time_s,averaging_time_s,common_clock_id,obstacle_geometry_id\n"
        "s1,10,0,1.5,0.03,2,60,clock-a,obstacle-a\n"
        "s2,20,0,1.5,0.05,2,60,clock-a,obstacle-a\n",
        encoding="utf-8",
    )
    model = tmp_path / "model.csv"
    model.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "s1,10,0,1.5,60,0.03\n"
        "s2,20,0,1.5,60,0.05\n",
        encoding="utf-8",
    )
    basis = {
        "source_boundary_id": "source-a",
        "weather_id": "weather-a",
        "sensor_geometry_id": "sensors-a",
        "temporal_operator_id": "mean-60s",
        "averaging_time_s": 60.0,
        "obstacle_representation_id": "mask-a",
    }
    case_path = tmp_path / "validation-case.json"
    case_path.write_text(json.dumps({
        "schema": "degali.field-validation-input.v1",
        "model": {
            "model_id": "model", "temporal_mode": "transient", "csv_path": model.name,
            "prediction_column": "prediction", "concentration_unit": "mole_fraction",
            "comparison_basis": basis,
        },
        "dataset": {
            "csv_path": observed.name, "temporal_mode": "transient",
            "comparison_basis": basis,
        },
        "validation_evidence": {
            "dataset_id": "observed-a", "path": observed.name,
            "sha256": hashlib.sha256(observed.read_bytes()).hexdigest(), "row_count": 2,
            "source_boundary_id": "source-a", "weather_id": "weather-a",
            "obstacle_geometry_id": "obstacle-a", "receptor_geometry_id": "sensors-a",
            "temporal_operator_id": "mean-60s", "common_clock_id": "clock-a",
            "scope": "lh2_obstacle_transport",
        },
    }), encoding="utf-8")
    return case_path, observed


def _screening_case(tmp_path):
    case = json.loads(
        (ROOT / "docs" / "field-screening-case.example.json").read_text(encoding="utf-8")
    )
    case["scenario"]["source"]["duration_s"] = 0.1
    case["scenario"]["source"]["mass_flow_kg_s"] = {
        "nominal": 0.265, "unit": "kg/s", "source": "FT-source-a",
    }
    case["scenario"]["weather"]["speed_m_s"] = {
        "nominal": 2.0, "unit": "m/s", "source": "met-a",
    }
    case["transport"].update({
        "length_m": 4.0, "height_m": 2.0, "nx": 10, "nz": 10,
        "time_step_s": 0.005, "duration_s": 0.1, "source_sigma_m": 0.2,
    })
    case["obstacles"] = []
    case["sensor_deployments"] = []
    path = tmp_path / "screening-case.json"
    path.write_text(json.dumps(case), encoding="utf-8")
    return path


def _pressure_screening_case(tmp_path):
    path = _screening_case(tmp_path)
    case = json.loads(path.read_text(encoding="utf-8"))
    source = case["scenario"]["source"]
    source.pop("mass_flow_kg_s")
    source["pressure_driven_mass_flow"] = {
        "ambient_pressure_pa": {
            "nominal": 101325.0,
            "lower": 101325.0,
            "upper": 101325.0,
            "unit": "Pa",
            "source": "ambient-a",
        },
        "source_id": "orifice-a",
    }
    path.write_text(json.dumps(case), encoding="utf-8")
    return path


def _pressure_sensor_envelope_case(tmp_path):
    path = _pressure_screening_case(tmp_path)
    case = json.loads(path.read_text(encoding="utf-8"))
    sensor = case["scenario"]["sensor"]
    sensor["response_time_s"] = {
        "nominal": 0.25,
        "lower": 0.2,
        "upper": 0.3,
        "unit": "s",
        "source": "calibration-a",
    }
    sensor["gain"] = {
        "nominal": 1.0,
        "lower": 0.9,
        "upper": 1.1,
        "unit": "1",
        "source": "calibration-a",
    }
    path.write_text(json.dumps(case), encoding="utf-8")
    return path


def _batch_case(tmp_path):
    case_path = _screening_case(tmp_path)
    path = tmp_path / "batch.json"
    path.write_text(json.dumps({
        "schema": "degali.field-batch-input.v1",
        "options": {"include_refinement": False, "max_uncertainty_cases": 8},
        "cases": [{"label": "case-a", "case_path": case_path.name}],
    }), encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("builder", "command", "output_key"),
    [
        (_source_case, "field-source", "atmospheric_source_schedule"),
        (_comparison_case, "field-compare", "field_model_comparison"),
        (_validation_case, "field-validate", "field_validation_score"),
    ],
)
def test_field_verify_recomputes_supported_execution_artifacts(
    tmp_path, capsys, builder, command, output_key,
):
    case_path, _referenced = builder(tmp_path)
    execution = tmp_path / f"{command}-execution.json"
    assert main([command, str(case_path), "--output", str(execution)]) == 0
    capsys.readouterr()

    assert main(["field-verify", str(execution)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema"] == "degali.field-execution-verification.v1"
    assert payload["execution_schema"].startswith("degali.field-")
    assert payload["input_sha256_verified"] is True
    assert payload["referenced_artifacts_verified"] is True
    assert payload["report_recomputed"] is True
    assert payload["promotion_allowed"] is False


def test_field_verify_fails_when_a_referenced_source_changes(tmp_path, capsys):
    case_path, csv_path = _source_case(tmp_path)
    execution = tmp_path / "source-execution.json"
    assert main(["field-source", str(case_path), "--output", str(execution)]) == 0
    capsys.readouterr()
    csv_path.write_text(csv_path.read_text(encoding="utf-8").replace("0.2", "0.3"), encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 1
    assert "SHA-256 does not match evidence" in capsys.readouterr().err


def test_field_verify_upgrades_a_legacy_source_schedule_operator(tmp_path, capsys):
    case_path, _csv_path = _source_case(tmp_path)
    execution = tmp_path / "legacy-source-execution.json"
    assert main(["field-source", str(case_path), "--output", str(execution)]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    payload["atmospheric_source_schedule"]["schedule"].pop("rate_operator")
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True


def test_field_verify_checks_screening_input_and_batch_manifest_artifacts(tmp_path, capsys):
    screening = _screening_case(tmp_path)
    screening_execution = tmp_path / "screening-execution.json"
    assert main([
        "field-screen", str(screening), "--no-refinement",
        "--output", str(screening_execution),
    ]) == 0
    capsys.readouterr()
    assert main(["field-verify", str(screening_execution)]) == 0
    screening_verification = json.loads(capsys.readouterr().out)
    assert screening_verification["execution_schema"] == "degali.field-screening-execution.v1"
    assert screening_verification["report_recomputed"] is True

    original_execution_text = screening_execution.read_text(encoding="utf-8")
    legacy_payload = json.loads(original_execution_text)
    legacy_payload["input"].pop("max_cell_steps")
    legacy_payload["field_report"]["schema"] = "degali.untrusted-field-report.v1"
    screening_execution.write_text(
        json.dumps(legacy_payload, indent=2) + "\n", encoding="utf-8",
    )
    assert main(["field-verify", str(screening_execution)]) == 1
    assert "unsupported schema" in capsys.readouterr().err
    screening_execution.write_text(original_execution_text, encoding="utf-8")

    refined_execution = tmp_path / "screening-refined-execution.json"
    assert main([
        "field-screen", str(screening), "--output", str(refined_execution),
    ]) == 0
    capsys.readouterr()
    assert main(["field-verify", str(refined_execution)]) == 0
    refined_verification = json.loads(capsys.readouterr().out)
    assert refined_verification["report_recomputed"] is True

    source_case, _source_csv = _source_case(tmp_path, duration=0.1)
    source_execution = tmp_path / "screening-source-execution.json"
    assert main([
        "field-screen", str(screening), "--no-refinement",
        "--uncertainty-envelope", "--atmospheric-source-case", str(source_case),
        "--output", str(source_execution),
    ]) == 0
    capsys.readouterr()
    assert main(["field-verify", str(source_execution)]) == 0
    source_verification = json.loads(capsys.readouterr().out)
    assert source_verification["report_recomputed"] is True

    batch = _batch_case(tmp_path)
    output_directory = tmp_path / "batch-output"
    assert main([
        "field-batch", str(batch), "--output-directory", str(output_directory),
    ]) == 0
    capsys.readouterr()
    execution = output_directory / "batch-execution.json"
    assert main(["field-verify", str(execution)]) == 0
    batch_verification = json.loads(capsys.readouterr().out)
    assert batch_verification["execution_schema"] == "degali.field-batch-execution.v1"
    assert batch_verification["referenced_artifacts_verified"] is True
    assert batch_verification["report_recomputed"] is False
    assert main([
        "field-verify", str(execution), "--require-recomputed",
    ]) == 2
    json.loads(capsys.readouterr().out)

    execution_payload = json.loads(execution.read_text(encoding="utf-8"))
    original_summary = dict(execution_payload["summary"])
    original_allowed = execution_payload["cases"][0]["screening_allowed"]
    original_gate_codes = list(execution_payload["cases"][0]["gate_codes"])
    execution_payload["summary"]["all_screening_allowed"] = not original_summary[
        "all_screening_allowed"
    ]
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "summary does not match" in capsys.readouterr().err
    execution_payload["summary"] = original_summary
    execution_payload["cases"][0]["screening_allowed"] = not execution_payload[
        "cases"
    ][0]["screening_allowed"]
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "screening_allowed mismatch" in capsys.readouterr().err
    execution_payload["cases"][0]["screening_allowed"] = original_allowed
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")

    execution_payload["cases"][0]["gate_codes"] = original_gate_codes + ["corner_withheld"]
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "gate_codes mismatch" in capsys.readouterr().err
    execution_payload["cases"][0]["gate_codes"] = original_gate_codes
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")

    manifest_path = output_directory / "manifest.json"
    manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    original_manifest_summary = json.loads(json.dumps(manifest_payload["summary"]))
    manifest_payload["summary"]["gate_code_counts"]["refinement_missing"] += 1
    manifest_path.write_text(json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    execution_payload["manifest_sha256"] = hashlib.sha256(
        manifest_path.read_bytes()
    ).hexdigest()
    execution_payload["summary"] = manifest_payload["summary"]
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "summary does not match its case decisions" in capsys.readouterr().err
    manifest_payload["summary"] = original_manifest_summary
    manifest_path.write_text(json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    execution_payload["manifest_sha256"] = hashlib.sha256(
        manifest_path.read_bytes()
    ).hexdigest()
    execution_payload["summary"] = original_manifest_summary
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")

    report = output_directory / "case-a.json"
    report.write_text(json.dumps({"schema": "degali.untrusted-report.v1"}), encoding="utf-8")
    report_sha = hashlib.sha256(report.read_bytes()).hexdigest()
    manifest = json.loads((output_directory / "manifest.json").read_text(encoding="utf-8"))
    manifest["cases"][0]["report_sha256"] = report_sha
    (output_directory / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
    )
    execution_payload = json.loads(execution.read_text(encoding="utf-8"))
    execution_payload["manifest_sha256"] = hashlib.sha256(
        (output_directory / "manifest.json").read_bytes()
    ).hexdigest()
    execution_payload["cases"][0]["report_sha256"] = report_sha
    execution.write_text(json.dumps(execution_payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "unsupported schema" in capsys.readouterr().err

    report.write_text(report.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "batch report changed since execution" in capsys.readouterr().err


def test_field_verify_rejects_unblocked_source_schedule_residual(
    tmp_path, capsys,
):
    """Integrity-only verification must still enforce transport fail-safe gates."""
    screening = _screening_case(tmp_path)
    execution = tmp_path / "screening-execution.json"
    assert main([
        "field-screen", str(screening), "--no-refinement",
        "--output", str(execution),
    ]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    # Removing one replay option intentionally selects the verifier's
    # integrity-only path; the input digest is not an execution self-hash.
    payload["input"].pop("max_cell_steps")
    diagnostics = payload["field_report"]["transport_result"]["diagnostics"]
    label = diagnostics["source_mass_schedule_residual_kg"][0][0]
    diagnostics["source_mass_schedule_residual_kg"] = [[label, 1.0]]
    diagnostics["maximum_source_mass_schedule_residual_kg"] = 1.0
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 1
    error = capsys.readouterr().err
    assert "source schedule residual" in error
    assert "blocked applicability reason" in error

    # Even a blocked report cannot be paired with a decision that omits the
    # machine-readable transport cause.
    payload["field_report"]["applicability"] = {
        "status": "blocked",
        "reasons": ["semi-FV source-mass schedule residual 1 kg exceeds the numerical gate"],
        "warnings": ["transport diagnostic was tampered for this verifier test"],
        "uncertainty_complete": False,
    }
    payload["operational_screening"]["gate_codes"] = []
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "missing transport gate_codes" in capsys.readouterr().err

    payload["operational_screening"]["gate_codes"] = [
        "refinement_missing",
    ]
    payload["field_report"]["sensor_result"]["maximum_true_mole_fraction"] = 2.0
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "outside [0, 1]" in capsys.readouterr().err

    payload["field_report"]["sensor_result"]["maximum_true_mole_fraction"] = 0.5
    payload["field_report"]["transport_result"]["diagnostics"][
        "source_mass_injected_kg"
    ][0][1] = 99.0
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "does not match injected source values" in capsys.readouterr().err


def test_field_verify_rejects_operational_decision_physical_status_mismatch(
    tmp_path, capsys,
):
    """Integrity-only verification must bind allowance to physical status."""
    screening_case = _screening_case(tmp_path)
    execution = tmp_path / "screening-execution.json"
    assert main([
        "field-screen", str(screening_case), "--no-refinement",
        "--output", str(execution),
    ]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    # Omit one replay option so the verifier takes its integrity-only path.
    payload["input"].pop("max_cell_steps")
    payload["field_report"]["applicability"] = {
        "status": "blocked",
        "reasons": ["tampered physical applicability"],
        "warnings": [],
        "uncertainty_complete": False,
    }
    decision = payload["operational_screening"]
    decision["status"] = "screening_allowed"
    decision["screening_allowed"] = True
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 1
    assert "accepted physical applicability" in capsys.readouterr().err


def test_field_verify_rejects_declared_schedule_mass_mismatch(
    tmp_path, capsys,
):
    """Integrity-only verification must bind source ledgers to input provenance."""
    screening = _screening_case(tmp_path)
    execution = tmp_path / "screening-execution.json"
    assert main([
        "field-screen", str(screening), "--no-refinement",
        "--output", str(execution),
    ]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    payload["input"].pop("max_cell_steps")
    diagnostics = payload["field_report"]["transport_result"]["diagnostics"]
    primary_mass = dict(diagnostics["source_mass_injected_kg"])["primary"]
    payload["field_report"]["transport_input"]["declared_direct_vapour_schedule"] = {
        "source_id": "tamper-test-schedule",
        "released_mass_kg": primary_mass,
    }
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    # The added provenance is internally consistent and therefore still
    # passes the integrity-only semantic verifier.
    assert main(["field-verify", str(execution)]) == 0
    capsys.readouterr()

    payload["field_report"]["transport_input"][
        "declared_direct_vapour_schedule"
    ]["released_mass_kg"] = primary_mass + 0.01
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "does not match its declared schedule mass" in capsys.readouterr().err

    payload["field_report"]["transport_input"][
        "declared_direct_vapour_schedule"
    ]["released_mass_kg"] = primary_mass
    payload["field_report"]["transport_input"]["distributed_vapour_sources"] = [{
        "label": "unlisted-source",
        "released_mass_kg": 0.0,
    }]
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "labels do not match declared transport inputs" in capsys.readouterr().err

    payload["field_report"]["transport_input"].pop("distributed_vapour_sources")
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "missing distributed_vapour_sources" in capsys.readouterr().err


def test_field_verify_rejects_closed_form_mass_total_tampering(
    tmp_path, capsys,
):
    """Integrity-only verification must recompute the final inventory identity."""
    screening = _screening_case(tmp_path)
    execution = tmp_path / "screening-execution.json"
    assert main([
        "field-screen", str(screening), "--no-refinement",
        "--output", str(execution),
    ]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    payload["input"].pop("max_cell_steps")
    diagnostics = payload["field_report"]["transport_result"]["diagnostics"]
    original_domain_mass = diagnostics["mass_domain_kg"]
    diagnostics["mass_domain_kg"] += 0.01
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 1
    assert "mass totals do not close" in capsys.readouterr().err

    diagnostics["mass_domain_kg"] = -1.0
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "mass totals must be non-negative" in capsys.readouterr().err

    diagnostics["mass_domain_kg"] = original_domain_mass
    diagnostics["final_mass_residual_kg"] = 0.1
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "final mass residual does not match mass totals" in capsys.readouterr().err


def test_field_verify_binds_pressure_derived_source_provenance(
    tmp_path, capsys,
):
    """Integrity-only verification must retain the typed pressure boundary."""
    screening = _pressure_screening_case(tmp_path)
    execution = tmp_path / "pressure-screening-execution.json"
    assert main([
        "field-screen", str(screening), "--no-refinement",
        "--output", str(execution),
    ]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    payload["input"].pop("max_cell_steps")
    source_record = payload["field_report"]["scenario"]["source"]
    source_record["pressure_driven_mass_flow"]["source_id"] = "tampered-source"
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 1
    error = capsys.readouterr().err
    assert "pressure-driven source identity does not match strict input" in error


def test_field_verify_binds_pressure_provenance_in_nested_sensor_envelope(
    tmp_path, capsys,
):
    """Nested envelope reports cannot hide a changed pressure source ID."""
    screening = _pressure_sensor_envelope_case(tmp_path)
    execution = tmp_path / "pressure-sensor-envelope-execution.json"
    assert main([
        "field-screen", str(screening), "--sensor-array-envelope",
        "--allow-conditional", "--relative-tolerance", "1.0",
        "--max-uncertainty-cases", "8", "--output", str(execution),
    ]) == 0
    capsys.readouterr()

    payload = json.loads(execution.read_text(encoding="utf-8"))
    payload["input"].pop("max_cell_steps")
    nested_source = payload["field_report"]["sensor_calibration_envelope"][
        "base_field_screening"
    ]["scenario"]["source"]
    nested_source["pressure_driven_mass_flow"]["source_id"] = "tampered-source"
    execution.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    assert main(["field-verify", str(execution)]) == 1
    error = capsys.readouterr().err
    assert "pressure-driven source identity does not match strict input" in error
