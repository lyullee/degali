import pytest

from tools.build_goal_completion_audit import build


def test_goal_audit_is_explicitly_conditional_and_lists_withheld_gates():
    record = build()
    assert record["schema"] == "degali.goal-completion-audit.v1"
    assert record["promotion_allowed"] is False
    assert record["completion_status"] == "conditional_research_and_screening_system"
    assert "grid_and_time_convergence" in record["blocking_requirements"]
    assert "public_obstacle_h2_validation" in record["blocking_requirements"]
    assert "embedded_surface_wall_shear" in record["blocking_requirements"]
    case_j_component = next(
        item for item in record["requirements"]
        if item["id"] == "nonisothermal_obstacle_inlet_component"
    )
    assert case_j_component["status"] == "conditional"
    assert case_j_component["details"]["height_varying_temperature_profile_supported"] is True
    assert case_j_component["details"]["heated_ground_boundary_supported"] is True
    assert case_j_component["details"]["heated_obstacle_surface_boundary_supported"] is True
    assert case_j_component["details"]["spatially_resolved_measured_wall_temperature_applied"] is False
    assert case_j_component["details"]["binary_mixture_source_species_and_enthalpy_coupled_to_projection"] is True
    assert case_j_component["details"]["finite_window_transient_mean_output_supported"] is True
    assert case_j_component["details"]["same_condition_model_comparison"] == "development_run_scored_validation_withheld"
    aij = next(
        item for item in record["requirements"]
        if item["id"] == "aij_neutral_wake_component"
    )
    assert aij["status"] == "conditional"
    assert aij["details"]["time_averaged_aligned_grid_status"] == "accepted_component_grid"
    assert aij["details"]["time_averaged_aligned_grid_passed"] is True
    assert aij["details"]["common_points"] == 176
    assert aij["details"]["finest_flow_converged"] is True
    assert aij["details"]["finest_obstacle_grid_aligned"] is True
    assert aij["details"]["quantitative_hydrogen_validation"] is False
    assert aij["details"]["concentration_validation_withheld"] is True
    tracer = aij["details"]["passive_tracer_resolution_sensitivity"]
    assert tracer["status"] == "joint_grid_timestep_source_regularization_sensitivity_not_convergence"
    assert tracer["common_measurement_point_count"] == 582
    assert tracer["coarse_flow_converged"] is True
    assert tracer["fine_flow_converged"] is True
    assert tracer["calibration_used"] is False
    assert tracer["promotion_allowed"] is False
    assert tracer["fine_rmse_ppmv"] > tracer["coarse_rmse_ppmv"]
    surface = next(
        item for item in record["requirements"]
        if item["id"] == "embedded_surface_wall_shear"
    )
    assert surface["status"] == "withheld"
    assert surface["details"]["existing_600s_screening_uses_updated_closure"] is False
    assert "20-s coarse execution witness" in surface["details"]["interpretation"]
    assert surface["details"]["embedded_arc_tangential_no_slip_stress_applied"] is True
    assert surface["details"]["wall_stress_benchmark_validation"] is False
    assert surface["details"]["one_step_smoke"]["duration_s"] == 0.02
    assert surface["details"]["one_step_smoke"]["all_numerical_gates_pass"] is True
    assert surface["details"]["one_step_smoke"]["maximum_wall_shear_stress_pa"] > 0.0
    transient = surface["details"]["two_second_transient_smoke_not_validation"]
    assert transient["duration_s"] == 2.0
    assert transient["cell_size_m"] == 5.0
    assert transient["all_numerical_gates_pass"] is True
    assert transient["wall_model"]["tangential_no_slip_stress_resolved"] is True
    assert transient["maximum_wall_shear_stress_pa"] > 0.0
    assert transient["maximum_relative_species_mass_residual"] < 1.0e-12
    assert transient["maximum_relative_sensible_enthalpy_residual"] < 1.0e-12
    assert "not wall-stress calibration" in transient["interpretation"]
    wall_20s = surface["details"]["twenty_second_coarse_transient_not_validation"]
    assert wall_20s["duration_s"] == 20.0
    assert wall_20s["cell_size_m"] == 5.0
    assert wall_20s["all_numerical_gates_pass"] is True
    assert 0.007 < wall_20s["maximum_wall_shear_stress_pa"] < 0.008
    assert wall_20s["maximum_relative_species_mass_residual"] < 1.0e-12
    comparison = surface["details"]["pre_post_wall_shear_20s_comparison"]
    assert comparison["common_stop_time_s"] == 20.0
    assert comparison["historical_4vol_percent_volume_difference_m3"] == 0.0
    assert comparison["sensor_trace_rmse_vol_percentage_points"] < 1.0e-4
    regions = comparison["regional_differences"]
    assert set(regions) == {
        "windward_inner_face_shell", "roof_shear_band",
        "side_tip_bypass", "exterior_lee_wake_shell",
    }
    assert regions["windward_inner_face_shell"]["common_fluid_cell_count"] > 0
    assert regions["windward_inner_face_shell"]["velocity_u_rmse_m_s"] > regions[
        "side_tip_bypass"
    ]["velocity_u_rmse_m_s"]
    sensitivity = surface["details"]["one_step_grid_sensitivity_not_validation"]
    assert sensitivity["cell_sizes_m"] == [5.0, 2.5, 1.25]
    assert sensitivity["cell_counts"] == [12800, 102400, 819200]
    assert sensitivity["wall_clock_s"][0] < sensitivity["wall_clock_s"][1] < sensitivity["wall_clock_s"][2]
    assert sensitivity["fine_to_coarse_stress_ratio"] > 6.0
    assert "not benchmark validation" in sensitivity["interpretation"]
    h2 = next(
        item for item in record["requirements"]
        if item["id"] == "public_obstacle_h2_validation"
    )
    assert h2["details"]["elvhys_tcs"] == (
        "conditional_cryogenic_h2_enclosure_vent_restriction_component"
    )
    assert h2["details"]["elvhys_tcs_open_field_wake_qualified"] is False
    grid = next(
        item for item in record["requirements"]
        if item["id"] == "grid_and_time_convergence"
    )
    assert grid["details"]["matched_window_temporal_sensitivity"]["stop_time_s"] == 20.0
    assert grid["details"]["receptor_gate"] is False
    fixed_surface_20s = grid["details"]["fixed_surface_20s_grid_diagnostic"]
    assert fixed_surface_20s["duration_s"] == 20.0
    assert fixed_surface_20s["field_norms_contract"] is True
    assert fixed_surface_20s["hazard_volume_gci_gate"] is False
    assert fixed_surface_20s["peak_monotonic_fraction"] == 0.24
    assert fixed_surface_20s["dose_monotonic_fraction"] == 0.36
    assert fixed_surface_20s["physical_geometry_held_constant"] is False
    assert fixed_surface_20s["full_event_promoted"] is False
    assert grid["details"]["matched_dt_obstacle_arc_length_relative_error"] == {
        "coarse": 0.10873231868538613,
        "medium": 0.04507034144862799,
        "fine": 0.06633811871569874,
    }
    sensor = next(
        item for item in record["requirements"]
        if item["id"] == "sensor_time_series_and_error_metrics"
    )
    nested = sensor["details"]["current_nested_145_sensor_600s_score"]
    assert nested["unique_sensor_count"] == 145
    assert nested["observed_positive_dose_sensor_count"] == 84
    assert nested["predicted_detected_sensor_count"] == 54
    assert nested["overlap_gate_passed"] is False
    assert nested["maximum_overlap_dose_change"] == pytest.approx(0.19616119135125276)
    assert nested["interior_overlap_location_count"] == 6
    assert nested["boundary_coincident_sensor_count"] == 9
    assert nested["dose_MAPE"] == pytest.approx(0.8498489814900475)
    assert nested["normalized_dose_bias"] == pytest.approx(-0.8985249120550491)
    assert nested["arrival_MAE_s"] == pytest.approx(31.764150943396228)
    assert nested["receptors_at_or_below_bottom_cell_centre_count"] > 0
    assert nested["decision"] == "withheld_as_unconditional_quantitative_validation"
    interior = sensor["details"]["short_true_interior_overlap_diagnostic"]
    assert interior["model_horizon_s"] == 60.0
    assert interior["interior_overlap_horizontal_location_count"] == 6
    assert interior["shared_sensor_count"] == 33
    assert interior["union_sensor_count"] == 145
    assert interior["expected_sensor_count"] == 145
    assert interior["maximum_peak_relative_change"] < 0.10
    assert interior["maximum_dose_relative_change"] < 0.10
    assert interior["inner_runtime_s"] > 0
    assert interior["outer_runtime_s"] > interior["inner_runtime_s"]
    assert interior["experimental_model_promotion_allowed"] is False
    latest_score = sensor["details"]["latest_nested_145_sensor_60s_score"]
    assert latest_score["model_horizon_s"] == 60.0
    assert latest_score["unique_sensor_count"] == 145
    assert latest_score["interior_overlap_sensor_count"] == 24
    assert latest_score["interior_overlap_location_count"] == 6
    assert latest_score["interior_overlap_gate_passed"] is True
    assert latest_score["decision"] == "withheld_as_unconditional_quantitative_validation"
    corrected_inner = sensor["details"]["corrected_nested_600s_inner_run"]
    assert corrected_inner["status"] == "accepted_completed_inner_transient"
    assert corrected_inner["all_internal_gates_pass"] is True
    assert corrected_inner["model_horizon_s"] == 600.0
    assert corrected_inner["domain_m"] == [300.0, 200.0, 40.0]
    assert corrected_inner["handoff_site_downwind_y_m"] == 100.0
    assert corrected_inner["outer_600s_run_completed"] is True
    assert corrected_inner["outer_runtime_s"] > 0
    assert corrected_inner["maximum_relative_species_mass_residual"] < 1.0e-14
