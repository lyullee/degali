from dataclasses import replace

import pytest

from degali.addons.field_comparison import (
    FieldComparisonBasis,
    FieldComparisonEvidence,
    FieldModelComparisonRow,
    FieldModelDecisionImpact,
    field_model_sensor_set_from_screening,
    FieldModelSensorSet,
    FieldSensorCsvProvenance,
    FieldSensorPrediction,
    compare_field_model_sensor_sets,
    field_model_decision_impact,
    field_model_decision_impact_record,
    field_model_comparison_report,
    field_model_sensor_set_from_csv,
    write_field_model_sensor_set_csv,
)
from degali.addons.field_comparison_io import FieldModelComparisonCase
from degali.addons.field_contracts import (
    BoundedValue,
    FieldApplicability,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _set(model_id, *, temporal_mode="steady", basis=None, value=0.03):
    return FieldModelSensorSet(
        model_id=model_id,
        temporal_mode=temporal_mode,
        basis=basis or FieldComparisonBasis("source-a", "weather-a", "sensors-a", "mean-60s", 60.0),
        predictions=(
            FieldSensorPrediction("s1", (10.0, 0.0, 1.5), value),
            FieldSensorPrediction("s2", (20.0, 0.0, 1.5), 0.05),
        ),
        runtime_s=0.1,
    )


def _screening_for_comparison():
    source = ReleaseSource(
        fluid="lh2", location_m=(0.0, 0.0, 0.5),
        upstream_pressure=BoundedValue(0.4e6),
        upstream_temperature=BoundedValue(26.084),
        mass_flow_kg_s=BoundedValue(0.265),
        opening_area_m2=BoundedValue(3.141592653589793 * 0.012**2 / 4.0 / 0.8),
        discharge_coefficient=BoundedValue(0.8), liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium", duration_s=1.0,
    )
    scenario = FieldScenario(
        source=source,
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        sensor=SensorModel((1.0, 0.0, 0.5), response_time_s=BoundedValue(0.2)),
        temporal_mode="transient",
    )
    transport = SemiFVConfig(
        length_m=10.0, height_m=5.0, nx=40, nz=20,
        duration_s=3.0, time_step_s=0.01, source_sigma_m=0.25,
    )
    return run_field_semi_fv_screening(FieldSemiFVRequest(scenario, transport=transport))


def test_matched_sensor_sets_report_model_form_difference_and_lfl_disagreement():
    comparison = compare_field_model_sensor_sets(
        _set("DEGALI", value=0.03), _set("SLABx", temporal_mode="transient", value=0.05)
    )
    assert comparison.applicability.status == "conditional"
    assert len(comparison.rows) == 2
    assert comparison.classification_disagreement_count == 1
    assert "temporal_mode_mismatch" in comparison.gate_codes
    assert comparison.mean_absolute_difference_mole_fraction == pytest.approx(0.01)
    impact = field_model_decision_impact(comparison)
    assert impact.status == "withheld"
    assert "model_selection_withheld" in impact.gate_codes
    assert "steady/transient" in impact.reasons[0]


def test_matched_basis_decision_impact_reports_only_monitored_threshold_change():
    comparison = compare_field_model_sensor_sets(
        _set("DEGALI", value=0.03), _set("SLABx", value=0.05),
    )
    impact = field_model_decision_impact(comparison)
    record = field_model_decision_impact_record(impact)

    assert comparison.applicability.status == "accepted"
    assert impact.status == "different"
    assert impact.left_alert_sensor_ids == ("s2",)
    assert impact.right_alert_sensor_ids == ("s1", "s2")
    assert impact.left_farthest_monitored_exceedance_x_m == pytest.approx(20.0)
    assert impact.right_farthest_monitored_exceedance_x_m == pytest.approx(20.0)
    assert impact.gate_codes == ("classification_disagreement",)
    assert record["qualification"].startswith("farthest monitored")
    assert record["gate_codes"] == ["classification_disagreement"]


def test_matched_basis_aligned_decision_has_explicit_alignment_code():
    comparison = compare_field_model_sensor_sets(_set("left"), _set("right"))
    impact = field_model_decision_impact(comparison)

    assert comparison.gate_codes == ("matched_basis",)
    assert impact.status == "aligned"
    assert impact.gate_codes == ("classification_aligned",)


def test_comparison_report_preserves_operator_and_input_scope():
    comparison = compare_field_model_sensor_sets(
        _set("DEGALI", value=0.03), _set("SLABx", value=0.05),
    )
    report = field_model_comparison_report(comparison)

    assert report["schema"] == "degali.field-model-comparison-report.v1"
    assert report["left"]["comparison_basis"]["temporal_operator_id"] == "mean-60s"
    assert report["comparison"]["position_tolerance_m"] == pytest.approx(1.0e-9)
    assert report["comparison"]["gate_codes"] == ["matched_basis"]
    assert report["comparison"]["rows"][0]["sensor_id"] == "s1"
    assert report["model_selection_impact"]["status"] == "different"
    assert report["comparison"]["applicability"]["gate_codes"] == ["matched_basis"]
    assert "not an independent accuracy validation" in report["scope"]


def test_screening_adapter_declares_trace_operator_and_preserves_conditional_execution():
    screening = _screening_for_comparison()
    assert screening.completed
    basis = FieldComparisonBasis(
        "field-source", "weather-a", "sensors-a", "peak-indicated",
    )
    model = field_model_sensor_set_from_screening(
        screening, model_id="DEGALI-field", basis=basis,
        prediction_operator="peak_indicated",
    )
    assert model.temporal_mode == "transient"
    assert model.predictions[0].mole_fraction == pytest.approx(
        max(screening.sensor_trace.indicated_mole_fraction)
    )
    assert model.execution_provenance is not None
    assert model.execution_provenance.prediction_operator == "peak_indicated"
    assert model.execution_provenance.source_applicability.status == "conditional"

    comparison = compare_field_model_sensor_sets(
        model, replace(model, model_id="external", execution_provenance=None),
    )
    assert comparison.applicability.status == "conditional"
    assert any("field execution is conditional" in item for item in comparison.applicability.warnings)
    assert field_model_decision_impact(comparison).status == "withheld"
    report = field_model_comparison_report(comparison)
    assert report["left"]["execution_provenance"]["prediction_operator"] == "peak_indicated"
    assert report["left"]["execution_provenance"]["source_applicability"]["status"] == "conditional"


def test_screening_adapter_uses_explicit_trace_operator_and_rejects_average_mismatch():
    screening = _screening_for_comparison()
    basis = FieldComparisonBasis(
        "field-source", "weather-a", "sensors-a", "full-trace-mean",
    )
    peak = field_model_sensor_set_from_screening(
        screening, model_id="peak", basis=basis, prediction_operator="peak_true",
    )
    final = field_model_sensor_set_from_screening(
        screening, model_id="final", basis=basis, prediction_operator="final_true",
    )
    mean = field_model_sensor_set_from_screening(
        screening, model_id="mean", basis=basis, prediction_operator="time_average_true",
    )
    assert peak.predictions[0].mole_fraction >= final.predictions[0].mole_fraction
    assert min(screening.sensor_trace.true_mole_fraction) <= mean.predictions[0].mole_fraction <= peak.predictions[0].mole_fraction
    with pytest.raises(ValueError, match="averaging time"):
        field_model_sensor_set_from_screening(
            screening,
            model_id="bad-average",
            basis=FieldComparisonBasis(
                "field-source", "weather-a", "sensors-a", "mean-60s", 60.0,
            ),
            prediction_operator="time_average_indicated",
        )


def test_screening_adapter_rejects_blocked_or_invalid_operator():
    screening = _screening_for_comparison()
    blocked = replace(
        screening,
        applicability=FieldApplicability(
            "blocked", reasons=("test block",),
        ),
        transport=None,
        sensor_trace=None,
        sensor_results=(),
    )
    with pytest.raises(ValueError, match="completed"):
        field_model_sensor_set_from_screening(
            blocked, model_id="bad", basis=_set("basis").basis,
        )
    with pytest.raises(ValueError, match="prediction_operator"):
        field_model_sensor_set_from_screening(
            screening, model_id="bad", basis=_set("basis").basis,
            prediction_operator="not-an-operator",
        )


def test_comparison_report_preserves_checked_external_qualification():
    comparison = compare_field_model_sensor_sets(
        _set("DEGALI", value=0.03), _set("SLABx", value=0.05),
    )
    evidence = FieldComparisonEvidence(
        path="C:/evidence/manifest.json",
        sha256="a" * 64,
        qualification="Native source models differ; do not rank them.",
    )

    report = field_model_comparison_report(comparison, comparison_evidence=evidence)

    assert report["comparison_evidence"] == {
        "path": "C:/evidence/manifest.json",
        "sha256": "a" * 64,
        "qualification": "Native source models differ; do not rank them.",
    }


def test_comparison_blocks_incompatible_source_or_sensor_geometry():
    incompatible = _set(
        "SLABx",
        basis=FieldComparisonBasis("different-source", "weather-a", "sensors-a", "mean-60s", 60.0),
    )
    comparison = compare_field_model_sensor_sets(_set("DEGALI"), incompatible)
    assert comparison.applicability.status == "conditional"
    assert "source boundary differs" in comparison.applicability.warnings[0]

    different_sensor = FieldModelSensorSet(
        model_id="SLABx",
        temporal_mode="steady",
        basis=_set("DEGALI").basis,
        predictions=(
            FieldSensorPrediction("s1", (10.1, 0.0, 1.5), 0.03),
            FieldSensorPrediction("s2", (20.0, 0.0, 1.5), 0.05),
        ),
    )
    blocked = compare_field_model_sensor_sets(_set("DEGALI"), different_sensor)
    assert blocked.applicability.status == "blocked"
    assert blocked.gate_codes == ("sensor_geometry_mismatch",)
    assert "s1" in blocked.applicability.reasons[0]


def test_obstacle_representation_mismatch_withholds_model_selection_impact():
    left_basis = FieldComparisonBasis(
        "source-a", "weather-a", "sensors-a", "mean-60s", 60.0, "mask-v1",
    )
    right_basis = FieldComparisonBasis(
        "source-a", "weather-a", "sensors-a", "mean-60s", 60.0, "mask-v2",
    )
    comparison = compare_field_model_sensor_sets(
        _set("DEGALI", basis=left_basis),
        _set("SLABx", basis=right_basis),
    )

    assert comparison.applicability.status == "conditional"
    assert any("obstacle representation" in warning for warning in comparison.applicability.warnings)
    assert field_model_decision_impact(comparison).status == "withheld"
    assert "basis_mismatch" in comparison.gate_codes


def test_undeclared_obstacle_representation_does_not_match_declared_one():
    declared = FieldComparisonBasis(
        "source-a", "weather-a", "sensors-a", "mean-60s", 60.0, "mask-v1",
    )
    comparison = compare_field_model_sensor_sets(
        _set("DEGALI"), _set("SLABx", basis=declared),
    )

    assert comparison.applicability.status == "conditional"
    assert "obstacle representation" in comparison.applicability.warnings[0]


def test_csv_sensor_import_preserves_geometry_units_and_file_fingerprint(tmp_path):
    source = tmp_path / "external.csv"
    source.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,pred_vol_pct\n"
        "A,10,0,1.5,60,4\n"
        "B,20,1,2.0,60,0.25\n",
        encoding="utf-8",
    )
    basis = FieldComparisonBasis("source-a", "weather-a", "sensors-a", "mean-60s", 60.0)

    imported = field_model_sensor_set_from_csv(
        source, model_id="external", temporal_mode="transient", basis=basis,
        prediction_column="pred_vol_pct", concentration_unit="volume_percent",
    )

    assert [item.mole_fraction for item in imported.predictions] == pytest.approx([0.04, 0.0025])
    assert imported.predictions[1].position_m == (20.0, 1.0, 2.0)
    assert imported.csv_provenance is not None
    assert imported.csv_provenance.row_count == 2
    assert len(imported.csv_provenance.sha256) == 64


def test_typed_model_sensor_set_writer_round_trips_with_fingerprint(tmp_path):
    model = _set("DEGALI")
    path = tmp_path / "degali-predictions.csv"

    exported = write_field_model_sensor_set_csv(model, path)

    assert exported.csv_provenance is not None
    assert exported.csv_provenance.path == str(path.resolve())
    assert exported.csv_provenance.row_count == 2
    assert exported.csv_provenance.averaging_time_column == "averaging_time_s"
    parsed = field_model_sensor_set_from_csv(
        path,
        model_id="DEGALI",
        temporal_mode="steady",
        basis=model.basis,
        prediction_column="mole_fraction",
        concentration_unit="mole_fraction",
    )
    assert parsed.predictions == model.predictions
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_field_model_sensor_set_csv(model, path)


def test_csv_sensor_import_rejects_a_different_declared_averaging_operator(tmp_path):
    source = tmp_path / "wrong-window.csv"
    source.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,averaging_time_s,prediction\n"
        "A,10,0,1.5,59,0.04\n",
        encoding="utf-8",
    )
    basis = FieldComparisonBasis("source-a", "weather-a", "sensors-a", "mean-60s", 60.0)

    with pytest.raises(ValueError, match="averaging time"):
        field_model_sensor_set_from_csv(
            source, model_id="external", temporal_mode="transient", basis=basis,
            prediction_column="prediction", concentration_unit="mole_fraction",
        )


def test_csv_sensor_import_rejects_ambiguous_headers_and_extra_values(tmp_path):
    basis = FieldComparisonBasis("source-a", "weather-a", "sensors-a", "snapshot")
    duplicate = tmp_path / "duplicate-header.csv"
    duplicate.write_text(
        "sensor,sensor,x_downwind_m,y_crosswind_m,height_m,prediction\n"
        "A,A,10,0,1.5,0.04\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate columns"):
        field_model_sensor_set_from_csv(
            duplicate, model_id="external", temporal_mode="steady", basis=basis,
            prediction_column="prediction", concentration_unit="mole_fraction",
            averaging_time_column=None,
        )

    extra = tmp_path / "extra-value.csv"
    extra.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,prediction\n"
        "A,10,0,1.5,0.04,unexpected\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="more values than header"):
        field_model_sensor_set_from_csv(
            extra, model_id="external", temporal_mode="steady", basis=basis,
            prediction_column="prediction", concentration_unit="mole_fraction",
            averaging_time_column=None,
        )

    missing = tmp_path / "missing-value.csv"
    missing.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,prediction\n"
        "A,10,0,1.5\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="fewer values than header"):
        field_model_sensor_set_from_csv(
            missing, model_id="external", temporal_mode="steady", basis=basis,
            prediction_column="prediction", concentration_unit="mole_fraction",
            averaging_time_column=None,
        )

    normalized = tmp_path / "normalized-header.csv"
    normalized.write_text(
        "sensor,x_downwind_m, x_downwind_m,y_crosswind_m,height_m,prediction\n"
        "A,10,10,0,1.5,0.04\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate columns"):
        field_model_sensor_set_from_csv(
            normalized, model_id="external", temporal_mode="steady", basis=basis,
            prediction_column="prediction", concentration_unit="mole_fraction",
            averaging_time_column=None,
        )


def test_csv_sensor_import_rejects_colliding_column_roles(tmp_path):
    source = tmp_path / "colliding-columns.csv"
    source.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,prediction\n"
        "A,10,0,1.5,0.04\n",
        encoding="utf-8",
    )
    basis = FieldComparisonBasis("source-a", "weather-a", "sensors-a", "snapshot")

    with pytest.raises(ValueError, match="column roles must be distinct"):
        field_model_sensor_set_from_csv(
            source,
            model_id="external",
            temporal_mode="steady",
            basis=basis,
            prediction_column="prediction",
            concentration_unit="mole_fraction",
            sensor_id_column="x_downwind_m",
            averaging_time_column=None,
        )


def test_comparison_withholds_unverified_csv_averaging_operator(tmp_path):
    source = tmp_path / "no-average-column.csv"
    source.write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,prediction\n"
        "s1,10,0,1.5,0.03\n"
        "s2,20,0,1.5,0.05\n",
        encoding="utf-8",
    )
    basis = FieldComparisonBasis(
        "source-a", "weather-a", "sensors-a", "mean-60s", 60.0,
    )
    imported = field_model_sensor_set_from_csv(
        source, model_id="external", temporal_mode="steady", basis=basis,
        prediction_column="prediction", concentration_unit="mole_fraction",
        averaging_time_column=None,
    )

    assert imported.csv_provenance is not None
    assert imported.csv_provenance.averaging_time_column is None
    comparison = compare_field_model_sensor_sets(imported, _set("DEGALI", basis=basis))

    assert comparison.applicability.status == "conditional"
    assert any("no averaging-time column" in warning for warning in comparison.applicability.warnings)
    assert field_model_decision_impact(comparison).status == "withheld"
    assert "averaging_operator_unverified" in comparison.gate_codes


def test_sensor_prediction_rejects_impossible_mole_fraction():
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        FieldSensorPrediction("sensor", (0.0, 0.0, 0.0), 1.001)


def test_model_comparison_input_boundary_rejects_mutable_or_incomplete_records():
    basis = _set("basis").basis
    with pytest.raises(TypeError, match="basis must be a FieldComparisonBasis"):
        FieldModelSensorSet(
            model_id="bad", temporal_mode="steady", basis="basis",
            predictions=_set("basis").predictions,
        )
    with pytest.raises(TypeError, match="predictions must be a tuple"):
        FieldModelSensorSet(
            model_id="bad", temporal_mode="steady", basis=basis,
            predictions=list(_set("basis").predictions),
        )
    with pytest.raises(ValueError, match="positive integer"):
        FieldSensorCsvProvenance(
            path="C:/evidence/model.csv", sha256="a" * 64, row_count=True,
            prediction_column="prediction", concentration_unit="mole_fraction",
        )
    comparison = compare_field_model_sensor_sets(_set("left"), _set("right"))
    blocked = replace(
        comparison,
        applicability=FieldApplicability("blocked", reasons=("manual block",)),
        rows=(), mean_absolute_difference_mole_fraction=None,
        classification_disagreement_count=0,
    )
    with pytest.raises(ValueError, match="blocked comparison"):
        replace(blocked, rows=comparison.rows)
    with pytest.raises(ValueError, match="classification disagreements"):
        replace(blocked, classification_disagreement_count=1)
    with pytest.raises(ValueError, match="mean difference must be None"):
        replace(blocked, mean_absolute_difference_mole_fraction=0.1)


def test_comparison_contract_rejects_boolean_numeric_inputs():
    with pytest.raises(TypeError, match="boolean"):
        FieldSensorPrediction("sensor", (0.0, 0.0, 0.0), True)
    with pytest.raises(TypeError, match="boolean"):
        FieldComparisonBasis("source", "weather", "sensor", "snapshot", True)
    with pytest.raises(TypeError, match="boolean"):
        FieldModelSensorSet(
            model_id="model", temporal_mode="steady", basis=_set("model").basis,
            predictions=_set("model").predictions, runtime_s=True,
        )
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        compare_field_model_sensor_sets(_set("left"), _set("right"), threshold_mole_fraction=True)
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        compare_field_model_sensor_sets(_set("left"), _set("right"), threshold_mole_fraction=1.01)
    with pytest.raises(ValueError, match="position_tolerance_m"):
        compare_field_model_sensor_sets(_set("left"), _set("right"), position_tolerance_m=True)
    with pytest.raises(TypeError, match="FieldModelSensorSet"):
        compare_field_model_sensor_sets("left", _set("right"))
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        FieldModelComparisonCase(
            _set("left"), _set("right"), threshold_mole_fraction=True,
        )
    with pytest.raises(ValueError, match="threshold_mole_fraction"):
        FieldModelComparisonCase(
            _set("left"), _set("right"), threshold_mole_fraction=1.01,
        )
    with pytest.raises(ValueError, match="position_tolerance_m"):
        FieldModelComparisonCase(
            _set("left"), _set("right"), position_tolerance_m=True,
        )
    with pytest.raises(TypeError, match="FieldModelSensorSet"):
        FieldModelComparisonCase("left", _set("right"))


def test_comparison_result_contract_rejects_inconsistent_rows():
    # The row constructor itself must reject a fabricated difference before a
    # malformed comparison can be serialised.
    with pytest.raises(ValueError, match="difference"):
        FieldModelComparisonRow(
            sensor_id="s1", position_m=(10.0, 0.0, 1.5),
            left_mole_fraction=0.03, right_mole_fraction=0.05,
            difference_mole_fraction=0.0, ratio_right_to_left=0.05 / 0.03,
            left_above_threshold=False, right_above_threshold=True,
        )

    comparison = compare_field_model_sensor_sets(_set("DEGALI"), _set("SLABx"))
    with pytest.raises(ValueError, match="disagreement count"):
        replace(comparison, classification_disagreement_count=1)


def test_model_decision_impact_contract_rejects_malformed_direct_records():
    comparison = compare_field_model_sensor_sets(_set("DEGALI"), _set("SLABx"))
    valid = field_model_decision_impact(comparison)
    with pytest.raises(TypeError, match="reasons must be a tuple"):
        replace(valid, reasons=["mutable reason"])
    with pytest.raises(TypeError, match="left_alert_sensor_ids must be a tuple"):
        replace(valid, left_alert_sensor_ids=list(valid.left_alert_sensor_ids))
    with pytest.raises(ValueError, match="withheld impact"):
        FieldModelDecisionImpact(
            comparison, "withheld", ("bad status",), (), (), None, None, None,
        )
    with pytest.raises(ValueError, match="alert identifiers"):
        FieldModelDecisionImpact(
            comparison, "different", ("bad ids",), ("missing",), (), None, None, None,
        )
    with pytest.raises(ValueError, match="runtime_ratio"):
        FieldModelDecisionImpact(
            comparison, "different", ("bad ratio",), (), (), None, None, -1.0,
        )
    with pytest.raises(ValueError, match="unsupported code"):
        replace(comparison, gate_codes=("not-a-code",))
    with pytest.raises(ValueError, match="unsupported code"):
        replace(valid, gate_codes=("not-a-code",))
    accepted = compare_field_model_sensor_sets(_set("left"), _set("right"))
    accepted_impact = field_model_decision_impact(accepted)
    with pytest.raises(ValueError, match="monitored distance"):
        replace(
            accepted_impact,
            left_farthest_monitored_exceedance_x_m=(
                accepted_impact.left_farthest_monitored_exceedance_x_m + 1.0
            ),
        )


def test_model_sensor_set_rejects_mismatched_csv_provenance_row_count():
    provenance = FieldSensorCsvProvenance(
        path="C:/evidence/model.csv", sha256="a" * 64, row_count=2,
        prediction_column="prediction", concentration_unit="mole_fraction",
    )
    with pytest.raises(ValueError, match="row_count"):
        FieldModelSensorSet(
            model_id="external", temporal_mode="steady",
            basis=FieldComparisonBasis("source-a", "weather-a", "sensors-a", "snapshot"),
            predictions=(FieldSensorPrediction("s1", (1.0, 0.0, 1.0), 0.01),),
            csv_provenance=provenance,
        )
