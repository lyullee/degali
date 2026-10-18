"""Build a compact, hash-pinned index of the current DEGALI validation evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


EVIDENCE = (
    ("resolved_600s_screening", "outputs/thorney-ti021-coarse-600s-dx2p5m-dt0p0625-checkpoint1s-tol1e3-2026-10-11/ti021-inner-transient.json", "accepted_full_horizon_screening_grid"),
    ("resolved_600s_ti021_score", "outputs/thorney-ti021-inner-transient-score-dx2p5m-2026-10-11-v1/ti021-inner-transient-score.json", "full_horizon_score_quantitative_promotion_withheld_grid_and_error_gates"),
    ("prior_nested_145_sensor_600s_score", "outputs/thorney-ti021-nested-145-sensor-score-v2-600s-2026-10-13-v3/ti021-nested-145-sensor-score.json", "prior_full_horizon_nested_score_old_handoff_has_boundary_only_overlap_not_current_geometry"),
    ("latest_nested_145_sensor_60s_score", "outputs/thorney-ti021-nested-145-sensor-score-60s-dx5m-handoff100-2026-10-13-v1/ti021-nested-145-sensor-score.json", "all_145_sensor_union_and_24_interior_overlap_sensors_pass_10pct_component_gate_60s_only_full_event_withheld"),
    ("current_nested_600s_inner_run", "outputs/thorney-ti021-inner-handoff-v3-600s-dx5m-2026-10-13-v1/ti021-inner-transient.json", "completed_600s_inner_run_all_internal_mass_energy_continuity_and_corrector_gates_pass"),
    ("current_nested_600s_inner_checkpoint", "outputs/thorney-ti021-inner-handoff-v3-600s-dx5m-2026-10-13-v1/ti021-inner-checkpoint.json", "finalized_600s_inner_checkpoint_journal"),
    ("current_nested_600s_outer_run", "outputs/thorney-ti021-outer-handoff-v3-600s-dx5m-dt0p125-2026-10-13-v1/ti021-outer-transient.json", "completed_600s_outer_run_all_internal_mass_energy_continuity_remap_and_handoff_gates_pass_coupled_dose_gate_failed"),
    ("current_nested_600s_outer_checkpoint", "outputs/thorney-ti021-outer-handoff-v3-600s-dx5m-dt0p125-2026-10-13-v1/ti021-outer-checkpoint.json", "finalized_600s_outer_checkpoint_journal"),
    ("current_nested_145_sensor_600s_score", "outputs/thorney-ti021-nested-145-sensor-score-v3-600s-2026-10-13-v2/ti021-nested-145-sensor-score.json", "full_event_all145_covered_peak_coupling_passes_dose_coupling_19p6pct_fails_experiment_dose_mape_85pct_validation_withheld"),
    ("ti021_corrected_wide_domain_120s_inner", "outputs/thorney-ti021-wide-inner-domain-120s-x450y640-centered-2026-10-13-v1/inner/ti021-inner-transient.json", "120s_diagnostic_centered_crosswind_origin_320m_x450_y640_domain_completed_internal_gates_pass_not_full_event_or_validation"),
    ("ti021_corrected_wide_domain_120s_outer", "outputs/thorney-ti021-wide-outer-domain-120s-x450y640-centered-2026-10-13-v1/outer/ti021-outer-transient.json", "120s_diagnostic_correct_q_bounds_plusminus320m_completed_internal_gates_pass_not_full_event_or_validation"),
    ("ti021_corrected_wide_domain_120s_nested_score", "outputs/thorney-ti021-wide-domain-overlap-score-120s-centered-2026-10-13-v1/ti021-nested-145-sensor-score.json", "54_interior_overlap_sensors_15_locations_coordinate_match_zero_but_all_sensor_peak_17p4pct_dose_116pct_low_signal_driven_withheld"),
    ("ti021_crosswind_span_prefix_comparison", "outputs/thorney-ti021-crosswind-span-prefix-comparison-2026-10-13-v1/crosswind-span-prefix-comparison.json", "120s_common_prefix_36_shared_interior_sensors_active12_peak_max_4p87_to_2p30pct_dose_12p97_to_7p20pct_diagnostic_only"),
    ("ti021_crosswind_span_prefix_comparison_tool", "tools/compare_thorney_ti021_crosswind_span_prefixes.py", "reproducible_hashed_same_prefix_same_coordinates_crosswind_domain_coupling_comparator"),
    ("ti021_invalid_wide_domain_200s_pilot", "outputs/thorney-ti021-wide-domain-overlap-score-200s-coordinate-check-v2/ti021-nested-145-sensor-score.json", "invalidated_due_to_crosswind_origin_and_handoff_edge_mapping_bug_do_not_use_as_physical_evidence"),
    ("prior_nested_600s_outer_checkpoint", "outputs/thorney-ti021-outer-handoff-v2-600s-dx5m-dt0p125-2026-10-12-v1/ti021-outer-checkpoint.json", "prior_finalized_600s_outer_run_superseded_by_expanded_inner_domain"),
    ("prior_nested_600s_inner_checkpoint", "outputs/thorney-ti021-inner-handoff-v2-600s-dx5m-2026-10-12-v1/ti021-inner-checkpoint.json", "prior_finalized_600s_inner_run_superseded_for_shared_interior_coverage"),
    ("thorney_interior_overlap_audit", "outputs/thorney-ti021-interior-overlap-60s-dx5m-handoff100-2026-10-13-v1/audit/ti021-interior-overlap-audit.json", "all_145_sensor_union_24_interior_overlap_sensors_6_locations_60s_diagnostic_only"),
    ("thorney_interior_overlap_audit_tool", "tools/audit_thorney_ti021_interior_overlap.py", "observation_blind_shared_interior_sensor_coupling_diagnostic"),
    ("resolved_600s_screening_package", "outputs/thorney-ti021-600s-screening-package-2026-10-12-v1/ti021-600s-screening-package.json", "front_top_side_wake_screening_package"),
    ("aij_case_h_resolved_3d_grid", "outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_resolved_3d_grid_audit.json", "conditional_neutral_wake_component"),
    ("aij_case_h_time_averaged_aligned_3d_grid", "outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_time_averaged_aligned_grid_audit.json", "accepted_neutral_wake_mean_flow_component_grid_concentration_withheld"),
    ("aij_case_h_lateral_domain_sensitivity", "outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_lateral_domain_sensitivity.json", "matched_lateral_boundary_sensitivity_flow_component_validation_withheld"),
    ("aij_case_h_lateral_domain_sensitivity_tool", "tools/audit_aij_case_h_lateral_domain.py", "reproducible_lateral_domain_boundary_sensitivity_audit"),
    ("aij_case_h_lateral_domain_tracer_comparison", "outputs/aij-case-h-lateral-domain-tracer-comparison-2026-10-14-v1/aij-case-h-lateral-domain-tracer-comparison.json", "paired_582_point_fv_tracer_lateral_boundary_sensitivity_validation_withheld"),
    ("aij_case_h_lateral_domain_tracer_comparator", "tools/compare_aij_case_h_lateral_domain_tracer.py", "reproducible_paired_lateral_domain_tracer_comparator"),
    ("aij_case_h_tracer_scoring_tool", "tools/score_aij_case_h_resolved_tracer.py", "no_fit_passive_tracer_component_scorer"),
    ("aij_case_h_tracer_coarse_run", "outputs/aij-case-h-resolved-tracer-coarse-2026-10-12-v1/aij-case-h-tracer-score.json", "20s_coarse_passive_tracer_gross_overprediction_not_validation"),
    ("aij_case_h_tracer_matched_coarse_run", "outputs/aij-case-h-resolved-tracer-matched-dx2p5cm-2026-10-12-v1/aij-case-h-tracer-score.json", "8to12s_common_window_coarse_passive_tracer_score"),
    ("aij_case_h_tracer_fine_run", "outputs/aij-case-h-resolved-tracer-fine-2026-10-12-v1/aij-case-h-tracer-score.json", "8to12s_fine_passive_tracer_gross_overprediction_not_validation"),
    ("aij_case_h_tracer_resolution_sensitivity", "outputs/aij-case-h-tracer-resolution-sensitivity-2026-10-12-v1/aij-case-h-tracer-resolution-sensitivity.json", "joint_grid_timestep_source_regularization_not_convergence"),
    ("aij_case_h_tracer_sgs_scorer", "tools/compare_aij_case_h_tracer_sgs.py", "paired_no_fit_sgs_scalar_sensitivity_scorer"),
    ("aij_case_h_tracer_sgs_run", "outputs/aij-case-h-tracer-sgs-sensitivity-dx2p5cm-2026-10-12-v1/aij-case-h-tracer-score.json", "uncalibrated_sgs_mixing_component_run_still_gross_overprediction"),
    ("aij_case_h_tracer_sgs_comparison", "outputs/aij-case-h-tracer-sgs-sensitivity-2026-10-12-v1/aij-case-h-tracer-sgs-sensitivity.json", "modest_rmse_improvement_not_validation_promotion_withheld"),
    ("aij_case_h_source_momentum_comparator", "tools/compare_aij_case_h_source_momentum.py", "paired_ethylene_air_fv_source_momentum_smoke_comparator"),
    ("aij_case_h_source_momentum_run", "outputs/aij-case-h-tracer-fv-momentum-smoke-2026-10-12-v1/aij-case-h-tracer-score.json", "1s_ethylene_air_fv_nozzle_momentum_execution_and_conservation_smoke"),
    ("aij_case_h_source_momentum_sensitivity", "outputs/aij-case-h-tracer-source-momentum-sensitivity-2026-10-12-v1/aij-case-h-source-momentum-sensitivity.json", "startup_window_momentum_effect_small_120s_validation_withheld"),
    ("aij_case_h_momentum_tier_comparator", "tools/compare_aij_case_h_tracer_momentum_tiers.py", "matched_mesh_time_source_sgs_d3q19_vs_ethylene_air_fv_no_fit_comparator"),
    ("aij_case_h_matched_d3q19_sgs_12s", "outputs/aij-case-h-tracer-d3q19-sgs-matcheddt-12s-dx2p5cm-dt1ms-2026-10-12-v1/aij-case-h-tracer-score.json", "matched_8to12s_d3q19_sgs_scalar_baseline"),
    ("aij_case_h_ethylene_air_fv_momentum_sgs_12s", "outputs/aij-case-h-tracer-fv-momentum-sgs-12s-dx2p5cm-dt1ms-2026-10-12-v1/aij-case-h-tracer-score.json", "matched_8to12s_ethylene_air_fv_and_inferred_source_momentum_candidate"),
    ("aij_case_h_matched_momentum_tier_comparison", "outputs/aij-case-h-momentum-tier-comparison-2026-10-12-v1/aij-case-h-momentum-tier-comparison.json", "55p9pct_rmse_reduction_candidate_120s_validation_and_promotion_withheld"),
    ("aij_case_h_momentum_tier_figure", "outputs/aij-case-h-momentum-tier-comparison-2026-10-12-v1/aij-case-h-momentum-tier-comparison.png", "matched_8to12s_sensor_parity_by_solver_tier"),
    ("aij_case_h_fv_zero_momentum_sgs_12s", "outputs/aij-case-h-tracer-fv-zero-momentum-sgs-12s-dx2p5cm-dt1ms-2026-10-12-v1/aij-case-h-tracer-score.json", "matched_fv_eos_pressure_sgs_control_without_source_momentum"),
    ("aij_case_h_fv_source_momentum_8to12s_comparator", "tools/compare_aij_case_h_source_momentum.py", "strict_matched_fv_case_source_momentum_toggle_comparator"),
    ("aij_case_h_fv_source_momentum_8to12s_sensitivity", "outputs/aij-case-h-fv-source-momentum-effect-12s-2026-10-12-v1/aij-case-h-source-momentum-sensitivity.json", "source_momentum_rmse_effect_below_0p04pct_validation_withheld"),
    ("aij_case_h_fv_source_momentum_8to12s_figure", "outputs/aij-case-h-fv-source-momentum-effect-12s-2026-10-12-v1/aij-case-h-source-momentum-sensitivity.png", "paired_field_collapse_shows_small_configured_momentum_effect"),
    ("active_thermal_buoyancy_witness", "outputs/obstacle-physics-rebuild-2026-10-10/obstacle_thermal_buoyancy_comparison.json", "conditional_development_witness"),
    ("dataverse_no_obstructed_release_audit", "outputs/dataverse-no-wskbij-audit-2026-10-12.json", "component_source_pressure_only"),
    ("resolved_20s_transient", "outputs/thorney-ti021-checkpoint-prefix-20s-dt0p02-2026-10-12-v1/ti021-inner-transient.json", "conditional_computed_window"),
    ("three_grid_convergence", "outputs/thorney-ti021-three-grid-impact-convergence-2026-10-12-v1/ti021-three-grid-impact-convergence.json", "hazard_and_receptor_promotion_withheld"),
    ("zero_surface_3s_grid_diagnostic", "outputs/thorney-ti021-zero-surface-matched-3s-convergence-2026-10-12-v2/ti021-three-grid-impact-convergence.json", "short_prefix_only_hazard_and_receptor_gates_withheld"),
    ("zero_surface_18s_prefix", "outputs/thorney-ti021-zero-surface-prefix-18s-dt0p02-2026-10-12-v2/ti021-inner-transient.json", "18s_conservative_computed_window_not_validation"),
    ("zero_surface_18s_checkpoint_snapshot", "outputs/thorney-ti021-zero-surface-prefix-18s-dt0p02-2026-10-12-v2/source-ti021-inner-checkpoint.json", "immutable_source_checkpoint_journal_snapshot"),
    ("zero_surface_18s_prefix_manifest", "outputs/thorney-ti021-zero-surface-prefix-18s-dt0p02-2026-10-12-v2/ti021-checkpoint-prefix.json", "segment_and_output_hash_manifest_with_all_gates"),
    ("zero_surface_18s_field_figure", "outputs/thorney-ti021-zero-surface-prefix-18s-dt0p02-2026-10-12-v2/ti021-zero-surface-prefix-18s-contours.png", "time_maximum_field_front_top_side_views_4volpct_contour"),
    ("zero_surface_2s_to18s_evolution", "outputs/thorney-ti021-zero-surface-prefix-evolution-2s-to18s-2026-10-12-v2/ti021-2s-to18s.json", "same_geometry_hazard_and_receptor_evolution_not_convergence"),
    ("zero_surface_18s_ti021_score", "outputs/thorney-ti021-zero-surface-prefix-18s-score-2026-10-12-v1/ti021-inner-transient-score.json", "horizon_limited_independent_sensor_score_withheld_not_validation"),
    ("zero_surface_18s_ti021_score_figure", "outputs/thorney-ti021-zero-surface-prefix-18s-score-2026-10-12-v1/ti021-inner-transient-score.png", "horizon_limited_sensor_score_visual"),
    ("thorney_18s_wall_geometry_sensitivity", "outputs/thorney-ti021-geometry-representation-prefix-18s-2026-10-12-v2.json", "matched_18s_voxel_vs_zero_surface_sensitivity_not_validation"),
    ("zero_surface_19s_prefix", "outputs/thorney-ti021-zero-surface-prefix-19s-dt0p02-2026-10-12-v1/ti021-inner-transient.json", "19s_conservative_computed_window_not_validation"),
    ("zero_surface_19s_checkpoint_snapshot", "outputs/thorney-ti021-zero-surface-prefix-19s-dt0p02-2026-10-12-v1/source-ti021-inner-checkpoint.json", "immutable_source_checkpoint_journal_snapshot"),
    ("zero_surface_19s_prefix_manifest", "outputs/thorney-ti021-zero-surface-prefix-19s-dt0p02-2026-10-12-v1/ti021-checkpoint-prefix.json", "segment_and_output_hash_manifest_with_all_gates"),
    ("zero_surface_19s_field_figure", "outputs/thorney-ti021-zero-surface-prefix-19s-dt0p02-2026-10-12-v1/ti021-zero-surface-prefix-19s-contours.png", "time_maximum_field_front_top_side_views_4volpct_contour"),
    ("zero_surface_19s_ti021_score", "outputs/thorney-ti021-zero-surface-prefix-19s-score-2026-10-12-v1/ti021-inner-transient-score.json", "horizon_limited_independent_sensor_score_withheld_not_validation"),
    ("thorney_19s_wall_geometry_sensitivity", "outputs/thorney-ti021-geometry-representation-prefix-19s-2026-10-12-v1.json", "matched_19s_voxel_vs_zero_surface_sensitivity_not_validation"),
    ("zero_surface_20s_prefix", "outputs/thorney-ti021-zero-surface-prefix-20s-dt0p02-2026-10-12-v1/ti021-inner-transient.json", "20s_conservative_computed_window_not_validation"),
    ("zero_surface_20s_checkpoint_snapshot", "outputs/thorney-ti021-zero-surface-prefix-20s-dt0p02-2026-10-12-v1/source-ti021-inner-checkpoint.json", "immutable_source_checkpoint_journal_snapshot"),
    ("zero_surface_20s_prefix_manifest", "outputs/thorney-ti021-zero-surface-prefix-20s-dt0p02-2026-10-12-v1/ti021-checkpoint-prefix.json", "segment_and_output_hash_manifest_with_all_gates"),
    ("zero_surface_20s_field_figure", "outputs/thorney-ti021-zero-surface-prefix-20s-dt0p02-2026-10-12-v1/ti021-zero-surface-prefix-20s-contours.png", "time_maximum_field_front_top_side_views_4volpct_contour"),
    ("zero_surface_20s_ti021_score", "outputs/thorney-ti021-zero-surface-prefix-20s-score-2026-10-12-v1/ti021-inner-transient-score.json", "horizon_limited_independent_sensor_score_withheld_not_validation"),
    ("zero_surface_20s_ti021_score_figure", "outputs/thorney-ti021-zero-surface-prefix-20s-score-2026-10-12-v1/ti021-inner-transient-score.png", "horizon_limited_sensor_score_visual"),
    ("thorney_20s_wall_geometry_sensitivity", "outputs/thorney-ti021-geometry-representation-prefix-20s-2026-10-12-v1.json", "matched_20s_voxel_vs_zero_surface_sensitivity_not_validation"),
    ("zero_surface_20s_three_grid_sensitivity", "outputs/thorney-ti021-zero-surface-matched-dt-grid-convergence-20s-2026-10-12-v2/ti021-three-grid-impact-convergence.json", "20s_fixed_ideal_surface_field_contracts_but_hazard_and_receptor_gates_fail_physical_scaffold_omitted"),
    ("zero_surface_20s_5m_vs_2p5m", "outputs/thorney-ti021-zero-surface-grid-5m-vs-2p5m-20s-2026-10-12-v1/ti021-inner-grid-comparison.json", "20s_two_grid_concentration_diagnostic_not_convergence"),
    ("zero_surface_20s_2p5m_vs_1p25m", "outputs/thorney-ti021-zero-surface-grid-2p5m-vs-1p25m-20s-2026-10-12-v1/ti021-inner-grid-comparison.json", "20s_two_grid_concentration_diagnostic_not_convergence"),
    ("embedded_wall_stress_smoke", "outputs/thorney-ti021-zero-surface-wall-stress-smoke-dx5m-1step-2026-10-12-v1/ti021-inner-transient.json", "one_step_smooth_wall_closure_smoke_not_benchmark_validated"),
    ("embedded_wall_stress_grid_smoke_2p5m", "outputs/thorney-ti021-zero-surface-wall-stress-smoke-2p5m-1step-2026-10-12-v1/ti021-inner-transient.json", "one_step_wall_stress_grid_sensitivity_not_validation"),
    ("embedded_wall_stress_grid_smoke_1p25m", "outputs/thorney-ti021-zero-surface-wall-stress-smoke-1p25m-1step-2026-10-12-v1/ti021-inner-transient.json", "one_step_wall_stress_grid_sensitivity_not_validation"),
    ("embedded_wall_stress_updated_2s_transient", "outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-2s-dt0p02-2026-10-12-v1/ti021-inner-transient.json", "updated_smooth_wall_closure_2s_coarse_transient_conservative_smoke_not_validation"),
    ("embedded_wall_stress_updated_20s_transient", "outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-20s-dt0p02-2026-10-12-v1/ti021-inner-transient.json", "updated_smooth_wall_closure_20s_coarse_transient_conservative_smoke_not_validation"),
    ("embedded_wall_stress_updated_20s_prefix", "outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-prefix-20s-2026-10-12-v1/ti021-checkpoint-prefix.json", "immutable_20s_materialized_prefix_with_all_gates"),
    ("embedded_wall_stress_20s_pre_post_comparison", "outputs/thorney-ti021-wall-shear-closure-20s-comparison-2026-10-12-v1.json", "same_grid_20s_multi_field_and_sensor_comparison_not_validation"),
    ("field_manifest_trial10", "outputs/preslhy-trial10-manifest-verification-2026-10-12.json", "conditional_calibration_missing"),
    ("h2safe_facility_audit", "outputs/h2safe-evidence-audit-2026-10-12-v1/h2safe-audit.json", "conditional_helium_surrogate"),
    ("exjet_benchmark_audit", "outputs/exjet-benchmark-audit-2026-10-12.json", "withheld_missing_machine_data"),
    ("bos_obstacle_benchmark_audit", "outputs/bos-obstacle-benchmark-audit-2026-10-12.json", "withheld_literature_aggregate_only"),
    ("public_obstacle_evidence_triage", "docs/public-obstacle-evidence-triage-2026-10-12.md", "direct_h2_experimental_and_embargoed_cfd_evidence_scope"),
    ("ffi_obstacle_registry", "reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json", "conditional_geometry_review_required_two_top_cross_beams_recorded"),
    ("ffi_paired_obstacle_effect_delta", "outputs/ffi-o1-control-obstacle-delta-2026-10-14.json", "conditional_paired_obstacle_effect_diagnostic_delta_pattern_mismatch_promotion_withheld"),
    ("ffi_paired_obstacle_effect_delta_tool", "tools/compare_ffi_control_obstacle_delta.py", "reproducible_control_vs_obstacle_delta_diagnostic_no_fit"),
    ("ffi_dataset_readme", "reference/hydrogen/open_channel_usn/ReadMe.txt", "public_experimental_dataset_metadata"),
    ("ffi_dataset_provenance", "reference/hydrogen/open_channel_usn/PROVENANCE.md", "download_and_hash_provenance"),
    ("ffi_control_case_a_run_1", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00012.zip", "raw_sensor_data_control_case_A"),
    ("ffi_control_case_a_run_2", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00013.zip", "raw_sensor_data_control_case_A"),
    ("ffi_obstacle_case_o1_run_1", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00027.zip", "raw_sensor_data_obstacle_case_O1"),
    ("ffi_obstacle_case_o1_run_2", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00028.zip", "raw_sensor_data_obstacle_case_O1"),
    ("ffi_control_case_e_run_1", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00003.zip", "raw_sensor_data_control_case_E"),
    ("ffi_control_case_e_run_2", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00004.zip", "raw_sensor_data_control_case_E"),
    ("ffi_obstacle_case_o2_run_1", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00029.zip", "raw_sensor_data_obstacle_case_O2"),
    ("ffi_obstacle_case_o2_run_2", "reference/hydrogen/open_channel_usn/23_FFI_P101_T00030.zip", "raw_sensor_data_obstacle_case_O2"),
    ("ffi_field_evidence_audit", "outputs/ffi-o1-field-evidence-audit-2026-10-14-v3.json", "wide_csv_parser_fixed_sha_pinned_channels_still_withheld"),
    ("ffi_field_evidence_intake_tool", "tools/audit_ffi_field_evidence_manifest.py", "reproducible_fail_closed_ffi_manifest_intake"),
    ("ffi_o1_field_evidence_intake", "outputs/ffi-o1-field-evidence-intake-2026-10-14-v2.json", "o1_source_receptor_conditional_weather_clock_obstacle_withheld"),
    ("ffi_o2_field_evidence_intake", "outputs/ffi-o2-field-evidence-intake-2026-10-14.json", "o2_source_receptor_conditional_weather_clock_obstacle_withheld"),
    ("ffi_hybrid_patch_60s_thermal", "outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-2026-10-14.json", "active_thermal_fast_resolved_fast_interface_witness_validation_withheld"),
    ("ffi_hybrid_patch_60s_flow20", "outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-flow20-2026-10-14.json", "active_thermal_flow_iteration_speed_sensitivity_validation_withheld"),
    ("ffi_hybrid_temporal_convergence", "outputs/ffi-o1-hybrid-temporal-convergence-3s-2026-10-14.json", "explicit_isothermal_time_interface_convergence_audit_validation_withheld"),
    ("ffi_hybrid_runner", "tools/run_ffi_o1_hybrid_patch_prefix.py", "reproducible_fast_resolved_fast_conservative_bridge_runner"),
    ("aij_case_j_source_provenance", "reference/obstacle_wake/aij_case_j_zenodo_15447888/PROVENANCE.md", "source_metadata_and_scope_boundary"),
    ("aij_case_j_approach_flow_data", "reference/obstacle_wake/aij_case_j_zenodo_15447888/AF_caseJ.csv", "raw_benchmark_used_in_development_comparison"),
    ("aij_case_j_surface_temperature_data", "reference/obstacle_wake/aij_case_j_zenodo_15447888/BD-ST_caseJ.csv", "raw_benchmark_used_as_uniform_wall_temperature_mean"),
    ("aij_case_j_wake_temperature_concentration_data", "reference/obstacle_wake/aij_case_j_zenodo_15447888/RS_caseJ.csv", "raw_benchmark_scored_440_of_444_points_validation_withheld"),
    ("aij_case_j_pointwise_scorer", "src/degali/validation/aij_case_j_comparison.py", "applied_development_run_scored_validation_withheld"),
    ("aij_case_j_development_comparison", "outputs/aij-case-j-resolved-development-0p60s-dx0p02-dt0p001-2026-10-12-v2/aij-case-j-resolved-development-comparison.json", "finite_transient_ethylene_score_validation_withheld"),
    ("aij_case_j_resolved_setup", "src/degali/validation/aij_case_j_setup.py", "physical_property_and_source_mapping_development_setup"),
    ("aij_case_j_resolved_runner", "tools/run_aij_case_j_resolved_comparison.py", "finite_transient_comparison_runner_validation_withheld"),
    ("aij_case_j_structural_dataset_audit", "outputs/aij-case-j-dataset-audit-2026-10-12.json", "checksum_verified_data_and_scored_development_comparison_validation_withheld"),
    ("elvhys_tcs_component_audit", "outputs/elvhys-tcs-component-audit-2026-10-12.json", "conditional_cryogenic_h2_vent_restriction_component"),
    ("completion_audit", "docs/degali-goal-completion-audit-2026-10-12.md", "requirement_matrix"),
)


def build() -> dict[str, Any]:
    records = []
    for evidence_id, relative, status in EVIDENCE:
        path = ROOT / relative
        records.append({
            "id": evidence_id,
            "path": relative,
            "absolute_path": str(path.resolve()),
            "exists": path.is_file(),
            "sha256": sha256(path),
            "status": status,
        })
    return {
        "schema": "degali.validation-evidence-index.v1",
        "root": str(ROOT.resolve()),
        "records": records,
        "all_files_present": all(item["exists"] for item in records),
        "promotion_allowed": False,
        "promotion_boundary": "independent calibrated obstacle-H2 observations and accepted hazard/receptor convergence remain required",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.output} all_files_present={record['all_files_present']}")
    return 0 if record["all_files_present"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
