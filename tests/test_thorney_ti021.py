from pathlib import Path
import json

import pytest

from degali.validation.thorney_ti021 import (
    THORNEY_TI021_EVIDENCE_SCHEMA,
    _parse_rows,
    thorney_ti021_evidence_record,
)


def _rows():
    rows = [
        ["# Dataset_Reference", "Thorney Island"],
        ["# Trial_Identification", "Trial 021"],
        ["# Release_Conditions", ""],
        ["substance", "R12"], ["Initial concentration", 0.32],
        ["release type (puff/jet/cyclone)", "puff"],
        ["release point x", 400.0], ["               y", 200.0],
        ["               z", 0.0], ["exit temperature(C)", 18.0],
        ["exit pressure (mb)", 1005.0], ["effective mol wt (g/mol)", 58.5],
        ["release rate  (kg/s)", 4962.0], ["release duration (s)", 1.0],
        ["#Meteorological Conditions", ""],
        ["site av wind direction (deg from N)", 194.2],
        ["ideal wind direction (deg from N)", 210.0],
        ["y-axis orientation (deg from N)", 30.0],
        ["site average windspeed (m/s)", 3.4],
        ["reference height for wind (m)", 10.0],
        ["friction velocity (m/s)", 0.205], ["surface roughness (m)", 0.005],
        ["stability class (from data prior to & during rel)", "D/E"],
        ["#Obstacle details", ""], ["radius (m)", 50.0],
        ["centre coordinates x, y (m)", 400200.0], ["height (m)", 6.0],
        ["porosity", 1.0], ["#Arc Positions", ""],
        ["X(m)", "height", "dose", "max", "peak", "width", "arr", "dep"],
    ]
    rows.extend([[75.0, float(i), 10.0, 1.0, 20.0, 4.0, 5.0, 30.0]
                 for i in range(6)])
    rows.extend([
        ["# Concentration_Sensor", ""],
        ["X(m)", "Y(m)", "height", "Dose", "arrival", "departure"],
    ])
    rows.extend([[400.0, 250.0, 0.4, float(i), 10.0, -999.0]
                 for i in range(145)])
    return rows


def test_ti021_parser_preserves_dose_semantics_and_site_coordinates():
    dataset = _parse_rows(
        _rows(), workbook_path=Path("ti21.xls"), workbook_sha256="a" * 64
    )
    assert dataset.workbook_wall_center_site_m == (400.0, 200.0)
    assert dataset.dose_sensors[0].local_x_m == pytest.approx(0.0)
    assert dataset.dose_sensors[0].local_y_m == pytest.approx(50.0)
    assert dataset.dose_sensors[-1].dose_percent_s == pytest.approx(144.0)
    assert dataset.dose_sensors[-1].departure_time_s is None
    assert dataset.arc_observations[0].maximum_concentration_percent == 1.0
    assert dataset.workbook_release_rate_kg_s == pytest.approx(4962.0)
    assert dataset.workbook_release_duration_s == pytest.approx(1.0)
    assert dataset.workbook_release_mass_proxy_kg == pytest.approx(4962.0)


def test_ti021_audit_fails_closed_on_report_workbook_conflicts(tmp_path):
    dataset = _parse_rows(
        _rows(), workbook_path=Path("ti21.xls"), workbook_sha256="a" * 64
    )
    report = tmp_path / "report.pdf"
    report.write_bytes(b"report")
    record = thorney_ti021_evidence_record(dataset, report_path=report)
    assert record["schema"] == THORNEY_TI021_EVIDENCE_SCHEMA
    assert record["decision"]["obstacle_dense_gas_component"] == "withheld"
    assert not record["promotion_gates"]["eligible_for_no_fit_model_run"]
    assert "wall_height_report_5m_vs_workbook_6m" in record[
        "conflicts_and_missing_evidence"
    ]
    assert "full_sensor_concentration_time_series_not_in_workbook" in record[
        "conflicts_and_missing_evidence"
    ]
    assert "sensor_arrival_peak_relative_operator_not_recoverable_from_workbook" in record[
        "conflicts_and_missing_evidence"
    ]
    assert record["observations"]["arrival_operator"]["definition"] == (
        "first concentration reaching 5 percent of the local peak"
    )


def test_ti021_audit_promotes_only_to_full_scale_preflight_after_component_audits(
    tmp_path,
):
    dataset = _parse_rows(
        _rows(), workbook_path=Path("ti21.xls"), workbook_sha256="a" * 64
    )
    report = tmp_path / "report.pdf"
    report.write_bytes(b"report")
    geometry = tmp_path / "geometry.json"
    geometry.write_text(json.dumps({
        "schema": "degali.thorney-ti021-curved-wall-geometry-audit.v1",
        "status": "accepted_geometry_component",
        "all_gates_pass": True,
    }), encoding="utf-8")
    initial = tmp_path / "initial.json"
    initial.write_text(json.dumps({
        "schema": "degali.thorney-ti021-initial-cloud-audit.v1",
        "status": "accepted_initial_state_component",
        "gates": {"geometry": True, "mixture": True, "canonical_step": True},
    }), encoding="utf-8")

    record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
    )

    assert record["status"] == (
        "conditional_pending_property_evidence_and_full_scale_no_fit_run"
    )
    assert record["promotion_gates"][
        "full_binary_r12_air_mixture_solver_available"
    ]
    assert record["promotion_gates"][
        "finite_volume_initial_cloud_state_supported"
    ]
    assert record["promotion_gates"]["eligible_for_full_scale_no_fit_preflight"]
    assert not record["promotion_gates"]["eligible_for_no_fit_model_run"]
    assert record["decision"]["obstacle_dense_gas_component"] == "conditional"

    properties = tmp_path / "properties.json"
    properties.write_text(json.dumps({
        "schema": "degali.thorney-ti021-r12-property-audit.v1",
        "status": "accepted_molecular_property_component",
        "gates": {"diffusivity": True, "cp": True},
    }), encoding="utf-8")
    property_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
    )
    assert property_record["status"] == (
        "conditional_pending_full_scale_no_fit_preflight"
    )
    assert property_record["promotion_gates"][
        "r12_molecular_property_inputs_frozen"
    ]
    assert not property_record["promotion_gates"][
        "turbulent_scalar_mixing_closure_frozen"
    ]

    preflight = tmp_path / "full-scale-preflight.json"
    preflight.write_text(json.dumps({
        "schema": "degali.thorney-ti021-full-scale-preflight.v1",
        "status": "accepted_full_scale_single_step_preflight",
        "all_gates_pass": True,
    }), encoding="utf-8")
    preflight_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
    )
    assert preflight_record["status"] == (
        "conditional_pending_turbulent_closure_and_no_fit_run"
    )
    assert preflight_record["promotion_gates"][
        "full_scale_single_step_preflight_passed"
    ]
    assert not preflight_record["promotion_gates"][
        "eligible_for_no_fit_model_run"
    ]

    mixing = tmp_path / "mixing.json"
    mixing.write_text(json.dumps({
        "schema": "degali.thorney-ti021-turbulent-mixing-policy.v1",
        "status": (
            "accepted_concentration_blind_turbulent_mixing_envelope_component"
        ),
        "all_gates_pass": True,
        "frozen_policy": {
            "ensemble": [
                {"member_id": f"member-{index}"} for index in range(6)
            ],
        },
    }), encoding="utf-8")
    ready_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
    )
    assert ready_record["status"] == (
        "conditional_pending_full_transient_convergence_and_no_fit_run"
    )
    assert ready_record["promotion_gates"][
        "turbulent_scalar_mixing_closure_frozen"
    ]
    assert ready_record["promotion_gates"]["eligible_for_inner_arc_no_fit_run"]
    assert not ready_record["promotion_gates"]["eligible_for_no_fit_model_run"]

    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({
        "schema": "degali.thorney-ti021-full-transient-plan.v2",
        "status": "accepted_coordinate_consistent_nested_execution_plan",
        "all_gates_pass": True,
        "coverage": {"nested_union": {"combined_unique_sensor_count": 145}},
    }), encoding="utf-8")
    planned_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
    )
    assert planned_record["status"] == (
        "conditional_pending_conservative_nested_handoff_and_full_transient"
    )
    assert planned_record["promotion_gates"][
        "sensor_complete_execution_plan_frozen"
    ]
    assert not planned_record["promotion_gates"][
        "conservative_nested_handoff_verified"
    ]
    assert not planned_record["promotion_gates"][
        "eligible_for_full_145_sensor_no_fit_run"
    ]

    handoff = tmp_path / "handoff.json"
    handoff.write_text(json.dumps({
        "schema": "degali.thorney-ti021-nested-handoff-audit.v1",
        "status": "accepted_conservative_handoff_remap_operator_only",
        "all_gates_pass": True,
    }), encoding="utf-8")
    remap_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
        nested_handoff_audit_path=handoff,
    )
    assert remap_record["status"] == (
        "conditional_pending_outer_boundary_coupling_and_full_transient"
    )
    assert remap_record["promotion_gates"][
        "conservative_handoff_remap_operator_verified"
    ]
    assert not remap_record["promotion_gates"][
        "conservative_nested_handoff_verified"
    ]

    outer_inlet = tmp_path / "outer-inlet.json"
    outer_inlet.write_text(json.dumps({
        "schema": "degali.thorney-ti021-outer-inlet-coupling-audit.v1",
        "status": (
            "accepted_outer_binary_boundary_component_pending_momentum_overlap"
        ),
        "all_gates_pass": True,
    }), encoding="utf-8")
    boundary_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
        nested_handoff_audit_path=handoff,
        outer_inlet_coupling_audit_path=outer_inlet,
    )
    assert boundary_record["status"] == (
        "conditional_pending_outer_momentum_overlap_and_full_transient"
    )
    assert boundary_record["promotion_gates"][
        "outer_binary_species_energy_continuity_boundary_verified"
    ]
    assert not boundary_record["promotion_gates"][
        "outer_momentum_and_overlap_verified"
    ]

    outer_momentum = tmp_path / "outer-momentum.json"
    outer_momentum.write_text(json.dumps({
        "schema": "degali.thorney-ti021-outer-momentum-overlap-audit.v1",
        "status": "accepted_outer_momentum_and_overlap_component",
        "all_gates_pass": True,
        "configuration": {"observation_value_fields_accessed": []},
    }), encoding="utf-8")
    nested_ready_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
        nested_handoff_audit_path=handoff,
        outer_inlet_coupling_audit_path=outer_inlet,
        outer_momentum_overlap_audit_path=outer_momentum,
    )
    assert nested_ready_record["status"] == (
        "conditional_pending_streaming_handoff_and_full_transient"
    )
    assert not nested_ready_record["promotion_gates"][
        "coordinate_consistent_streaming_handoff_verified"
    ]
    streaming = tmp_path / "streaming.json"
    streaming.write_text(json.dumps({
        "schema": "degali.thorney-ti021-handoff-streaming-preflight.v1",
        "status": (
            "accepted_coordinate_consistent_streaming_handoff_preflight"
        ),
        "all_gates_pass": True,
        "configuration": {"observation_value_fields_accessed": []},
        "gates": {
            "durable_hdf5_reopened_with_expected_frames": True,
            "all_five_fluxes_conserved_by_area_remap": True,
        },
    }), encoding="utf-8")
    nested_ready_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
        nested_handoff_audit_path=handoff,
        outer_inlet_coupling_audit_path=outer_inlet,
        outer_momentum_overlap_audit_path=outer_momentum,
        handoff_streaming_preflight_path=streaming,
    )
    assert nested_ready_record["status"] == (
        "conditional_pending_full_transient_convergence_and_no_fit_scoring"
    )
    assert nested_ready_record["promotion_gates"][
        "outer_momentum_and_overlap_verified"
    ]
    assert nested_ready_record["promotion_gates"][
        "conservative_nested_handoff_verified"
    ]
    assert nested_ready_record["promotion_gates"][
        "coordinate_consistent_streaming_handoff_verified"
    ]
    assert nested_ready_record["promotion_gates"][
        "eligible_for_full_145_sensor_no_fit_run"
    ]
    assert nested_ready_record["promotion_gates"][
        "eligible_for_no_fit_model_run"
    ]

    completed_inner = tmp_path / "completed-inner.json"
    completed_inner.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-transient.v1",
        "status": "accepted_completed_inner_transient",
        "all_gates_pass": True,
        "configuration": {
            "duration_s": 600.0,
            "inner_receptor_count": 60,
            "initial_binary_cloud_momentum": "quiescent",
            "observation_values_used_to_configure_run": [],
        },
    }), encoding="utf-8")
    inner_score = tmp_path / "inner-score.json"
    inner_score.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-transient-score.v1",
        "status": "scored_no_fit_inner_dx5_not_grid_converged",
        "protocol": {"configuration_accessed_observation_values": []},
        "resolution_limit": {"promotion_allowed": False},
        "decision": {"result": "withheld_as_quantitative_validation"},
    }), encoding="utf-8")
    scored_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
        nested_handoff_audit_path=handoff,
        outer_inlet_coupling_audit_path=outer_inlet,
        outer_momentum_overlap_audit_path=outer_momentum,
        handoff_streaming_preflight_path=streaming,
        completed_inner_transient_path=completed_inner,
        inner_transient_score_path=inner_score,
    )
    assert scored_record["status"] == (
        "conditional_pending_outer_transient_and_grid_convergence"
    )
    assert scored_record["promotion_gates"][
        "primary_dx5_inner_600s_transient_completed"
    ]
    assert scored_record["promotion_gates"][
        "coarse_dx5_inner_no_fit_score_available"
    ]
    assert not scored_record["promotion_gates"][
        "coarse_dx5_inner_quantitative_validation_promoted"
    ]
    assert not scored_record["promotion_gates"][
        "dx2p5_corrector_tolerance_verified_with_peak_postcheck_pending"
    ]

    tolerance = tmp_path / "corrector-tolerance.json"
    tolerance.write_text(json.dumps({
        "schema": "degali.thorney-ti021-corrector-tolerance-audit.v1",
        "status": (
            "accepted_1e4_corrector_tolerance_for_dx2p5_production_with_postcheck"
        ),
        "all_gates_pass": True,
        "configuration": {
            "cell_size_m": 2.5,
            "reference_tolerance": 1.0e-6,
            "candidate_tolerance": 1.0e-4,
            "observation_value_fields_accessed": [],
        },
        "decision": {"mandatory_postcheck": "strict peak-window comparison"},
    }), encoding="utf-8")
    tolerance_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        corrector_tolerance_audit_path=tolerance,
    )
    assert tolerance_record["promotion_gates"][
        "dx2p5_corrector_tolerance_verified_with_peak_postcheck_pending"
    ]
    assert tolerance_record["provenance"][
        "corrector_tolerance_audit"
    ]["status"].startswith("accepted_1e4")

    refined_completed = tmp_path / "refined-completed-inner.json"
    refined_completed.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-transient.v1",
        "status": "accepted_completed_inner_transient",
        "all_gates_pass": True,
        "configuration": {
            "duration_s": 600.0,
            "cell_size_m": 2.5,
            "inner_receptor_count": 60,
            "initial_binary_cloud_momentum": "quiescent",
            "observation_values_used_to_configure_run": [],
        },
    }), encoding="utf-8")
    refined_score = tmp_path / "refined-inner-score.json"
    refined_score.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-transient-score.v1",
        "status": "scored_no_fit_inner_dx2p5_pending_finest_grid",
        "protocol": {"configuration_accessed_observation_values": []},
        "resolution_limit": {
            "cell_size_m": 2.5, "promotion_allowed": False,
        },
        "decision": {"result": "withheld_as_quantitative_validation"},
    }), encoding="utf-8")
    strict = tmp_path / "strict-peak-window.json"
    strict.write_text(json.dumps({
        "schema": "degali.thorney-ti021-strict-peak-window-audit.v1",
        "status": "accepted_1e4_corrector_tolerance_at_model_peak",
        "configuration": {
            "cell_size_m": 2.5,
            "candidate_tolerance": 1.0e-4,
            "strict_tolerance": 1.0e-6,
        },
        "selection": {"observation_value_fields_accessed": []},
        "gates": {"representative_gate": True},
    }), encoding="utf-8")
    comparison = tmp_path / "full-event-grid-comparison.json"
    comparison.write_text(json.dumps({
        "schema": "degali.thorney-ti021-inner-grid-comparison.v1",
        "status": "diagnostic_two_grid_comparison_pending_1p25m_window",
        "grid_convergence_promoted": False,
        "configuration": {
            "coarse_cell_size_m": 5.0,
            "fine_cell_size_m": 2.5,
            "duration_s": 600.0,
        },
        "provenance": {"observation_value_fields_accessed": []},
    }), encoding="utf-8")
    refined_record = thorney_ti021_evidence_record(
        dataset,
        report_path=report,
        curved_wall_geometry_audit_path=geometry,
        initial_cloud_audit_path=initial,
        r12_property_audit_path=properties,
        full_scale_preflight_path=preflight,
        turbulent_mixing_policy_path=mixing,
        full_transient_plan_path=plan,
        nested_handoff_audit_path=handoff,
        outer_inlet_coupling_audit_path=outer_inlet,
        outer_momentum_overlap_audit_path=outer_momentum,
        handoff_streaming_preflight_path=streaming,
        coarse_inner_transient_path=completed_inner,
        coarse_inner_score_path=inner_score,
        completed_inner_transient_path=refined_completed,
        inner_transient_score_path=refined_score,
        corrector_tolerance_audit_path=tolerance,
        strict_peak_window_audit_path=strict,
        full_event_grid_comparison_path=comparison,
    )
    assert refined_record["status"] == (
        "conditional_pending_outer_transient_and_grid_convergence"
    )
    gates = refined_record["promotion_gates"]
    assert gates["refined_dx2p5_inner_600s_transient_completed"]
    assert gates["refined_dx2p5_inner_no_fit_score_available"]
    assert gates["dx2p5_strict_peak_window_postcheck_passed"]
    assert gates["full_event_5m_to_2p5m_grid_diagnostic_available"]
    assert gates["primary_dx5_inner_600s_transient_completed"]
    assert gates["coarse_dx5_inner_no_fit_score_available"]
    assert not gates["full_event_grid_convergence_promoted"]
