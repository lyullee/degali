import hashlib
import json
import math
from dataclasses import replace

import pytest

from degali.addons.field_batch import (
    FIELD_BATCH_MANIFEST_SCHEMA,
    FieldBatchCase,
    export_field_screening_batch,
)
import degali.addons.field_batch as field_batch_module
from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_decision import (
    FieldConditionalReviewAuthorization,
    FieldOperationalScreeningDecision,
)
from degali.addons.field_source_io import (
    FieldAtmosphericSourceSchedule,
    FieldSourceScheduleEvidence,
)
from degali.addons.field_workflow import FieldSemiFVRequest
from degali.addons.field_workflow import FieldSensorDeployment
from degali.addons.semi_fv_obstacle import SemiFVConfig, SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


def _request():
    return FieldSemiFVRequest(
        FieldScenario(
            source=ReleaseSource(
                fluid="lh2", location_m=(0.0, 0.0, 0.5),
                upstream_pressure=BoundedValue(0.4e6),
                upstream_temperature=BoundedValue(26.084),
                mass_flow_kg_s=BoundedValue(0.265),
                opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8),
                discharge_coefficient=BoundedValue(0.8),
                liquid_fraction=BoundedValue(0.922),
                flash_model="homogeneous_equilibrium", duration_s=0.1,
            ),
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((0.6, 0.0, 0.5)), temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=10, nz=10,
            duration_s=0.1, time_step_s=0.005, source_sigma_m=0.2,
        ),
    )


def _review():
    return FieldConditionalReviewAuthorization(
        review_id="batch-review-A",
        reviewer_id="process-safety-reviewer-01",
        reviewer_role="process_safety_engineer",
        reviewed_at_utc="2026-10-05T00:00:00Z",
        evidence_id="batch-review-record-A",
    )


def _atmospheric_source_schedule():
    time_s = (0.0, 0.05, 0.1)
    return FieldAtmosphericSourceSchedule(
        schedule=SourceRateSchedule(time_s, (0.12, 0.08, 0.0), source_id="post-flash-A"),
        lower_schedule=SourceRateSchedule(time_s, (0.10, 0.06, 0.0), source_id="post-flash-A"),
        upper_schedule=SourceRateSchedule(time_s, (0.14, 0.10, 0.0), source_id="post-flash-A"),
        evidence=FieldSourceScheduleEvidence(
            dataset_id="batch-atmospheric-A",
            path="C:/evidence/atmospheric-A.csv",
            sha256="0" * 64,
            row_count=3,
            source_boundary_id="post-flash-A",
            common_clock_id="clock-A",
            source_kind="post_flash_atmospheric_vapour",
        ),
        time_column="time_s",
        rate_column="rate_kg_s",
        lower_rate_column="rate_lower_kg_s",
        upper_rate_column="rate_upper_kg_s",
    )


def test_batch_export_writes_complete_named_reports_without_refinement(tmp_path):
    exported = export_field_screening_batch(
        {"nominal": _request()}, tmp_path / "batch", include_refinement=False,
    )

    manifest = json.loads(exported.manifest_path.read_text(encoding="utf-8"))
    assert manifest["schema"] == FIELD_BATCH_MANIFEST_SCHEMA
    assert manifest["cases"][0]["report_file"] == "nominal.json"
    assert not manifest["cases"][0]["refinement_requested"]
    assert manifest["cases"][0]["operational_decision"]["status"] == "withheld"
    assert manifest["cases"][0]["operational_decision"]["refinement_required"]
    assert manifest["operational_screening_policy"]["require_resolved_uncertainty"]
    assert manifest["summary"] == {
        "schema": "degali.field-batch-summary.v1",
        "case_count": 1,
        "decision_status_counts": {
            "screening_allowed": 0,
            "conditional_allowed": 0,
            "withheld": 1,
        },
        "screening_allowed_case_count": 0,
        "conditional_allowed_case_count": 0,
        "withheld_case_count": 1,
        "all_screening_allowed": False,
            "gate_code_counts": {
                "obstacle_transport_conditional": 1,
                "physical_applicability_conditional": 1,
                "refinement_missing": 1,
            },
    }
    assert not exported.cases[0].operational_decision.screening_allowed
    assert json.loads((tmp_path / "batch" / "nominal.json").read_text(encoding="utf-8"))["schema"] == (
        "degali.field-screening.v1"
    )
    with pytest.raises(FileExistsError, match="non-empty"):
        export_field_screening_batch({"new": _request()}, tmp_path / "batch")


def test_batch_report_render_failure_creates_no_partial_output_directory(tmp_path, monkeypatch):
    def fail_render(_result):
        raise ValueError("synthetic report-render failure")

    monkeypatch.setattr(field_batch_module, "field_screening_report", fail_render)
    output = tmp_path / "render-failure"
    with pytest.raises(ValueError, match="synthetic report-render failure"):
        export_field_screening_batch(
            {"nominal": _request()}, output, include_refinement=False,
        )
    assert not output.exists()


def test_default_batch_preserves_a_blocked_refinement_instead_of_omitting_it(tmp_path):
    exported = export_field_screening_batch(
        {"screen": _request()}, tmp_path / "refinement-batch", max_cell_steps=1,
    )

    report = json.loads((tmp_path / "refinement-batch" / "screen.json").read_text(encoding="utf-8"))
    assert exported.cases[0].refinement_requested
    assert not exported.cases[0].refinement_completed
    assert report["refinement_applicability"]["status"] == "blocked"
    assert report["numerical_refinement"] is None
    manifest = json.loads((tmp_path / "refinement-batch" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["cases"][0]["operational_decision"]["status"] == "withheld"
    assert "refinement" in " ".join(manifest["cases"][0]["operational_decision"]["reasons"])


def test_batch_conditional_opt_in_requires_and_preserves_per_case_review(tmp_path):
    with pytest.raises(ValueError, match="missing: nominal"):
        export_field_screening_batch(
            {"nominal": _request()},
            tmp_path / "missing-review",
            include_refinement=False,
            allow_conditional_operational_screening=True,
            require_refinement_for_operational_screening=False,
        )
    assert not (tmp_path / "missing-review").exists()

    exported = export_field_screening_batch(
        {"nominal": _request()},
        tmp_path / "reviewed-batch",
        include_refinement=False,
        allow_conditional_operational_screening=True,
        conditional_review_authorizations={"nominal": _review()},
        require_refinement_for_operational_screening=False,
    )
    manifest = json.loads(exported.manifest_path.read_text(encoding="utf-8"))

    assert manifest["operational_screening_policy"]["conditional_review_required"]
    assert manifest["cases"][0]["conditional_review"]["review_id"] == "batch-review-A"
    assert manifest["cases"][0]["conditional_review"]["applied_to_this_execution"]
    assert manifest["cases"][0]["operational_decision"]["status"] == "conditional_allowed"
    assert manifest["summary"]["decision_status_counts"] == {
        "screening_allowed": 0,
        "conditional_allowed": 1,
        "withheld": 0,
    }
    assert "conditional_review_required" in manifest["summary"]["gate_code_counts"]
    assert manifest["summary"]["all_screening_allowed"]
    assert exported.cases[0].conditional_review_applied

    with pytest.raises(ValueError, match="requires an applied conditional review"):
        replace(
            exported.cases[0],
            conditional_review=None,
            conditional_review_applied=False,
        )
    with pytest.raises(ValueError, match="conditional physical applicability"):
        replace(exported.cases[0], applicability_status="accepted")


def test_batch_case_rejects_status_and_refinement_metadata_mismatches(tmp_path):
    exported = export_field_screening_batch(
        {"nominal": _request()}, tmp_path / "metadata-batch", include_refinement=False,
    )
    case = exported.cases[0]
    allowed = FieldOperationalScreeningDecision(
        "screening_allowed", True, False, (), (),
        refinement_required=False, refinement_available=False,
        uncertainty_resolution_required=False, uncertainty_resolved=True,
    )
    with pytest.raises(ValueError, match="accepted physical applicability"):
        replace(case, operational_decision=allowed)
    with pytest.raises(ValueError, match="completed refinement"):
        replace(case, refinement_completed=True)
    with pytest.raises(TypeError, match="FieldOperationalScreeningDecision"):
        FieldBatchCase(
            case.label, case.report_path, case.applicability_status,
            case.refinement_requested, case.refinement_completed,
            "not-a-decision",
        )


def test_batch_export_rejects_missing_or_out_of_directory_artifacts(tmp_path):
    exported = export_field_screening_batch(
        {"nominal": _request()}, tmp_path / "artifact-batch", include_refinement=False,
    )
    empty_output = tmp_path / "empty-output"
    empty_output.mkdir()
    with pytest.raises(ValueError, match="manifest_path must point to an existing"):
        replace(
            exported,
            output_directory=empty_output,
            manifest_path=empty_output / "manifest.json",
        )
    with pytest.raises(ValueError, match="report paths must be inside"):
        replace(
            exported,
            cases=(replace(exported.cases[0], report_path=tmp_path / "nominal.json"),),
        )
    with pytest.raises(ValueError, match="report path must match"):
        replace(
            exported,
            cases=(replace(exported.cases[0], report_path=exported.output_directory / "other.json"),),
        )


def test_batch_export_rejects_tampered_manifest_and_report_content(tmp_path):
    exported = export_field_screening_batch(
        {"nominal": _request()}, tmp_path / "tampered-batch", include_refinement=False,
    )
    manifest_path = exported.manifest_path
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["cases"][0]["report_file"] = "other.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="report_file does not match"):
        replace(exported, cases=exported.cases)

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["cases"][0]["report_file"] = "nominal.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    report_path = exported.cases[0].report_path
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report["schema"] = "tampered.schema.v1"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    with pytest.raises(ValueError, match="report_sha256"):
        replace(exported, cases=exported.cases)

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["cases"][0]["report_sha256"] = hashlib.sha256(
        report_path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported schema"):
        replace(exported, cases=exported.cases)


def test_batch_atmospheric_source_schedule_uses_complete_operational_envelope(tmp_path):
    exported = export_field_screening_batch(
        {"source": _request()},
        tmp_path / "source-batch",
        include_refinement=False,
        atmospheric_source_schedules={"source": _atmospheric_source_schedule()},
        max_uncertainty_cases=8,
        table_nodes=17,
    )

    report = json.loads((tmp_path / "source-batch" / "source.json").read_text(encoding="utf-8"))
    manifest = json.loads(exported.manifest_path.read_text(encoding="utf-8"))
    assert report["schema"] == "degali.field-operational-uncertainty-envelope.v1"
    assert report["field_uncertainty_envelope"]["atmospheric_source_schedule"]["evidence"]["dataset_id"] == (
        "batch-atmospheric-A"
    )
    assert len(report["cases"]) == 6
    assert {case["selection"]["source_schedule"] for case in report["cases"]} == {
        "lower", "nominal", "upper",
    }
    assert manifest["cases"][0]["execution_kind"] == "operational_uncertainty_envelope"
    assert manifest["cases"][0]["operational_decision"]["status"] == "withheld"
    assert exported.cases[0].applicability_status in {"accepted", "conditional", "blocked"}


def test_batch_atmospheric_source_schedule_cannot_waive_operational_gates(tmp_path):
    with pytest.raises(ValueError, match="require refinement"):
        export_field_screening_batch(
            {"source": _request()},
            tmp_path / "invalid-source-policy",
            include_refinement=False,
            atmospheric_source_schedules={"source": _atmospheric_source_schedule()},
            require_refinement_for_operational_screening=False,
        )
    assert not (tmp_path / "invalid-source-policy").exists()

    exported = export_field_screening_batch(
        {"source": _request()},
        tmp_path / "conditional-source-policy",
        atmospheric_source_schedules={"source": _atmospheric_source_schedule()},
        allow_conditional_operational_screening=True,
        conditional_review_authorizations={"source": _review()},
        relative_tolerance=1.0,
        max_uncertainty_cases=8,
        table_nodes=17,
    )
    manifest = json.loads(
        (tmp_path / "conditional-source-policy" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["cases"][0]["conditional_review"]["review_id"] == "batch-review-A"
    assert manifest["cases"][0]["conditional_review"]["applied_to_this_execution"]
    assert manifest["cases"][0]["operational_decision"]["status"] == "conditional_allowed"
    assert exported.cases[0].conditional_review_applied


def test_batch_manifest_keeps_unrepresented_obstacle_and_withheld_sensor_in_sync(tmp_path):
    request = replace(
        _request(),
        obstacles=(AxisAlignedCuboid(
            1.5, 2.0, 2.0, 3.0, 0.0, 1.0, "off-plane-obstacle",
        ),),
        sensor_deployments=(FieldSensorDeployment(
            "off-plane-detector", SensorModel((1.0, 1.0, 0.5)),
        ),),
    )

    exported = export_field_screening_batch(
        {"geometry": request}, tmp_path / "geometry-batch", include_refinement=False,
    )
    report = json.loads(
        (tmp_path / "geometry-batch" / "geometry.json").read_text(encoding="utf-8")
    )
    manifest = json.loads(
        (tmp_path / "geometry-batch" / "manifest.json").read_text(encoding="utf-8")
    )
    decision = manifest["cases"][0]["operational_decision"]

    assert exported.cases[0].operational_decision.status == "withheld"
    assert report["obstacle_projections"][0]["local_obstacle"] is None
    assert report["obstacle_projections"][0]["status"] == "accepted"
    assert report["sensor_results"][1]["label"] == "off-plane-detector"
    assert report["sensor_results"][1]["withheld"]
    assert any("off-plane-detector" in reason for reason in decision["reasons"])
    assert any("off-plane-obstacle" in reason for reason in decision["reasons"])
    assert decision["status"] == "withheld"
