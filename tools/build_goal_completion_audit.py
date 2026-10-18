"""Build a machine-readable requirement audit for the DEGALI reconstruction goal."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def _read(relative: str) -> dict[str, Any]:
    path = ROOT / relative
    return json.loads(path.read_text(encoding="utf-8"))


def build() -> dict[str, Any]:
    index = _read("outputs/degali-validation-evidence-index-2026-10-13.json")
    screening = _read("outputs/thorney-ti021-600s-screening-package-2026-10-12-v1/ti021-600s-screening-package.json")
    screening_score = _read("outputs/thorney-ti021-inner-transient-score-dx2p5m-2026-10-11-v1/ti021-inner-transient-score.json")
    nested_145_score = _read("outputs/thorney-ti021-nested-145-sensor-score-v3-600s-2026-10-13-v2/ti021-nested-145-sensor-score.json")
    interior_overlap = _read("outputs/thorney-ti021-interior-overlap-60s-dx5m-handoff100-2026-10-13-v1/audit/ti021-interior-overlap-audit.json")
    latest_nested_60s_score = _read("outputs/thorney-ti021-nested-145-sensor-score-60s-dx5m-handoff100-2026-10-13-v1/ti021-nested-145-sensor-score.json")
    corrected_nested_600s_inner = _read("outputs/thorney-ti021-inner-handoff-v3-600s-dx5m-2026-10-13-v1/ti021-inner-transient.json")
    corrected_nested_600s_outer = _read("outputs/thorney-ti021-outer-handoff-v3-600s-dx5m-dt0p125-2026-10-13-v1/ti021-outer-transient.json")
    grid = _read("outputs/thorney-ti021-matched-dt-grid-convergence-2026-10-12-v4/ti021-three-grid-impact-convergence.json")
    grid_zero_surface_20s = _read("outputs/thorney-ti021-zero-surface-matched-dt-grid-convergence-20s-2026-10-12-v2/ti021-three-grid-impact-convergence.json")
    grid_short = _read("outputs/thorney-ti021-zero-surface-matched-3s-convergence-2026-10-12-v2/ti021-three-grid-impact-convergence.json")
    temporal = _read("outputs/thorney-ti021-temporal-prefix-comparison-20s-2026-10-12-v2/ti021-temporal-prefix-comparison-20s.json")
    aij = _read("outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_resolved_3d_grid_audit.json")
    aij_mean = _read("outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_time_averaged_aligned_grid_audit.json")
    aij_tracer_resolution = _read("outputs/aij-case-h-tracer-resolution-sensitivity-2026-10-12-v1/aij-case-h-tracer-resolution-sensitivity.json")
    aij_tracer_sgs = _read("outputs/aij-case-h-tracer-sgs-sensitivity-2026-10-12-v1/aij-case-h-tracer-sgs-sensitivity.json")
    aij_source_momentum = _read("outputs/aij-case-h-fv-source-momentum-effect-12s-2026-10-12-v1/aij-case-h-source-momentum-sensitivity.json")
    aij_tier_comparison = _read("outputs/aij-case-h-momentum-tier-comparison-2026-10-12-v1/aij-case-h-momentum-tier-comparison.json")
    aij_lateral_domain = _read("outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_lateral_domain_sensitivity.json")
    aij_lateral_tracer = _read("outputs/aij-case-h-lateral-domain-tracer-comparison-2026-10-14-v1/aij-case-h-lateral-domain-tracer-comparison.json")
    ffi_hybrid_thermal = _read("outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-2026-10-14.json")
    ffi_hybrid_flow20 = _read("outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-flow20-2026-10-14.json")
    ffi_hybrid_temporal = _read("outputs/ffi-o1-hybrid-temporal-convergence-3s-2026-10-14.json")
    thermal = _read("outputs/obstacle-physics-rebuild-2026-10-10/obstacle_thermal_buoyancy_comparison.json")
    h2safe = _read("outputs/h2safe-evidence-audit-2026-10-12-v1/h2safe-audit.json")
    exjet = _read("outputs/exjet-benchmark-audit-2026-10-12.json")
    bos = _read("outputs/bos-obstacle-benchmark-audit-2026-10-12.json")
    dataverse = _read("outputs/dataverse-no-wskbij-audit-2026-10-12.json")
    elvhys = _read("outputs/elvhys-tcs-component-audit-2026-10-12.json")
    ffi_o1_intake = _read("outputs/ffi-o1-field-evidence-intake-2026-10-14-v2.json")
    ffi_o2_intake = _read("outputs/ffi-o2-field-evidence-intake-2026-10-14.json")
    wall_smoke = _read(
        "outputs/thorney-ti021-zero-surface-wall-stress-smoke-dx5m-1step-2026-10-12-v1/ti021-inner-transient.json"
    )
    wall_smoke_2p5m = _read(
        "outputs/thorney-ti021-zero-surface-wall-stress-smoke-2p5m-1step-2026-10-12-v1/ti021-inner-transient.json"
    )
    wall_smoke_1p25m = _read(
        "outputs/thorney-ti021-zero-surface-wall-stress-smoke-1p25m-1step-2026-10-12-v1/ti021-inner-transient.json"
    )
    wall_stress_2s = _read(
        "outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-2s-dt0p02-2026-10-12-v1/ti021-inner-transient.json"
    )
    wall_stress_20s = _read(
        "outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-20s-dt0p02-2026-10-12-v1/ti021-inner-transient.json"
    )
    wall_stress_20s_comparison = _read(
        "outputs/thorney-ti021-wall-shear-closure-20s-comparison-2026-10-12-v1.json"
    )
    case_j = _read("outputs/aij-case-j-dataset-audit-2026-10-12.json")
    zero_surface_18s = _read("outputs/thorney-ti021-zero-surface-prefix-18s-dt0p02-2026-10-12-v2/ti021-inner-transient.json")
    zero_surface_18s_evolution = _read("outputs/thorney-ti021-zero-surface-prefix-evolution-2s-to18s-2026-10-12-v2/ti021-2s-to18s.json")
    zero_surface_18s_score = _read("outputs/thorney-ti021-zero-surface-prefix-18s-score-2026-10-12-v1/ti021-inner-transient-score.json")
    geometry_18s = _read("outputs/thorney-ti021-geometry-representation-prefix-18s-2026-10-12-v2.json")
    zero_surface_19s = _read("outputs/thorney-ti021-zero-surface-prefix-19s-dt0p02-2026-10-12-v1/ti021-inner-transient.json")
    zero_surface_19s_score = _read("outputs/thorney-ti021-zero-surface-prefix-19s-score-2026-10-12-v1/ti021-inner-transient-score.json")
    geometry_19s = _read("outputs/thorney-ti021-geometry-representation-prefix-19s-2026-10-12-v1.json")
    zero_surface_20s = _read("outputs/thorney-ti021-zero-surface-prefix-20s-dt0p02-2026-10-12-v1/ti021-inner-transient.json")
    zero_surface_20s_score = _read("outputs/thorney-ti021-zero-surface-prefix-20s-score-2026-10-12-v1/ti021-inner-transient-score.json")
    geometry_20s = _read("outputs/thorney-ti021-geometry-representation-prefix-20s-2026-10-12-v1.json")

    requirements = [
        {
            "id": "embedded_surface_wall_shear",
            "status": "withheld",
            "evidence": (
                "resolved_obstacle_flow_3d._solid_mask, "
                "advance_variable_density_momentum_3d call, and "
                "_increment_projected_faces_from_cell_delta"
            ),
            "details": {
                "zero_thickness_arc_excluded_from_solid_mask": True,
                "finite_volume_predictor_receives_solid_mask_only": True,
                "blocked_arc_faces_enforce_zero_normal_flux": True,
                "embedded_arc_tangential_no_slip_stress_applied": True,
                "implemented_closure": "smooth-wall Reichardt immersed stress",
                "existing_600s_screening_uses_updated_closure": False,
                "wall_stress_benchmark_validation": False,
                "one_step_smoke": {
                    "status": wall_smoke.get("status"),
                    "duration_s": wall_smoke.get("configuration", {}).get("duration_s"),
                    "all_numerical_gates_pass": wall_smoke.get("all_gates_pass"),
                    "maximum_wall_shear_stress_pa": wall_smoke.get(
                        "embedded_arc_wall_stress_summary", {}
                    ).get("maximum_shear_stress_pa"),
                    "maximum_momentum_balance_residual_kg_m_s": wall_smoke.get(
                        "numerics", {}
                    ).get("maximum_momentum_balance_residual_kg_m_s"),
                },
                "one_step_grid_sensitivity_not_validation": {
                    "cell_sizes_m": [5.0, 2.5, 1.25],
                    "cell_counts": [
                        wall_smoke.get("configuration", {}).get("cell_count"),
                        wall_smoke_2p5m.get("configuration", {}).get("cell_count"),
                        wall_smoke_1p25m.get("configuration", {}).get("cell_count"),
                    ],
                    "wall_clock_s": [
                        wall_smoke.get("runtime", {}).get("wall_clock_s"),
                        wall_smoke_2p5m.get("runtime", {}).get("wall_clock_s"),
                        wall_smoke_1p25m.get("runtime", {}).get("wall_clock_s"),
                    ],
                    "maximum_wall_shear_stress_pa": [
                        wall_smoke.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa"),
                        wall_smoke_2p5m.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa"),
                        wall_smoke_1p25m.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa"),
                    ],
                    "fine_to_coarse_stress_ratio": (
                        wall_smoke_1p25m.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa")
                        / wall_smoke.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa")
                    ),
                    "interpretation": "single 0.02 s startup step; strong grid sensitivity; not benchmark validation or grid convergence",
                },
                "two_second_transient_smoke_not_validation": {
                    "status": wall_stress_2s.get("status"),
                    "duration_s": wall_stress_2s.get("configuration", {}).get("duration_s"),
                    "cell_size_m": wall_stress_2s.get("configuration", {}).get("cell_size_m"),
                    "wall_model": wall_stress_2s.get("obstacle_wall_model"),
                    "all_numerical_gates_pass": wall_stress_2s.get("all_gates_pass"),
                    "maximum_wall_shear_stress_pa": wall_stress_2s.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa"),
                    "maximum_friction_velocity_m_s": wall_stress_2s.get("embedded_arc_wall_stress_summary", {}).get("maximum_friction_velocity_m_s"),
                    "cumulative_wall_impulse_n_s": wall_stress_2s.get("embedded_arc_wall_stress_summary", {}).get("cumulative_wall_impulse_n_s"),
                    "maximum_relative_species_mass_residual": wall_stress_2s.get("numerics", {}).get("maximum_relative_species_mass_residual"),
                    "maximum_relative_sensible_enthalpy_residual": wall_stress_2s.get("numerics", {}).get("maximum_relative_sensible_enthalpy_residual"),
                    "maximum_momentum_balance_residual_kg_m_s": wall_stress_2s.get("numerics", {}).get("maximum_momentum_balance_residual_kg_m_s"),
                    "wall_clock_s": wall_stress_2s.get("runtime", {}).get("wall_clock_s"),
                    "interpretation": "updated closure active in a 2-s coarse transient; verifies execution and conservation only, not wall-stress calibration, separated-wake validation, concentration accuracy, or full-event convergence",
                },
                "twenty_second_coarse_transient_not_validation": {
                    "status": wall_stress_20s.get("status"),
                    "duration_s": wall_stress_20s.get("configuration", {}).get("duration_s"),
                    "cell_size_m": wall_stress_20s.get("configuration", {}).get("cell_size_m"),
                    "all_numerical_gates_pass": wall_stress_20s.get("all_gates_pass"),
                    "maximum_wall_shear_stress_pa": wall_stress_20s.get("embedded_arc_wall_stress_summary", {}).get("maximum_shear_stress_pa"),
                    "maximum_friction_velocity_m_s": wall_stress_20s.get("embedded_arc_wall_stress_summary", {}).get("maximum_friction_velocity_m_s"),
                    "cumulative_wall_impulse_n_s": wall_stress_20s.get("embedded_arc_wall_stress_summary", {}).get("cumulative_wall_impulse_n_s"),
                    "maximum_relative_species_mass_residual": wall_stress_20s.get("numerics", {}).get("maximum_relative_species_mass_residual"),
                    "maximum_relative_sensible_enthalpy_residual": wall_stress_20s.get("numerics", {}).get("maximum_relative_sensible_enthalpy_residual"),
                    "maximum_momentum_balance_residual_kg_m_s": wall_stress_20s.get("numerics", {}).get("maximum_momentum_balance_residual_kg_m_s"),
                    "wall_clock_s": wall_stress_20s.get("runtime", {}).get("wall_clock_s"),
                    "interpretation": "updated closure active over a 20-s coarse impact window; numerical execution/conservation only, not wall-stress or wake validation",
                },
                "pre_post_wall_shear_20s_comparison": {
                    "common_stop_time_s": wall_stress_20s_comparison.get("receptors", {}).get("common_stop_time_s"),
                    "max_concentration_field_rmse_vol_percentage_points": wall_stress_20s_comparison.get("physical_field_components", {}).get("compared_components", {}).get("concentration_vol_percent", {}).get("rmse"),
                    "max_u_velocity_difference_m_s": wall_stress_20s_comparison.get("physical_field_components", {}).get("compared_components", {}).get("velocity_u_m_s", {}).get("maximum_absolute_difference"),
                    "sensor_trace_rmse_vol_percentage_points": wall_stress_20s_comparison.get("receptors", {}).get("rmse"),
                    "historical_4vol_percent_volume_difference_m3": wall_stress_20s_comparison.get("field", {}).get("hazard_volume_difference_m3"),
                    "regional_differences": {
                        label: {
                            "common_fluid_cell_count": region.get("common_fluid_cell_count"),
                            "concentration_rmse_vol_percentage_points": region.get("components", {}).get("concentration_vol_percent", {}).get("rmse"),
                            "velocity_u_rmse_m_s": region.get("components", {}).get("velocity_u_m_s", {}).get("rmse"),
                            "pressure_rmse_pa": region.get("components", {}).get("pressure_pa_gauge", {}).get("rmse"),
                        }
                        for label, region in wall_stress_20s_comparison.get("physical_field_components", {}).get("obstacle_region_components", {}).items()
                    },
                    "interpretation": "pre/post closure scalar changes are negligible on this 5-m 20-s case; this is a numerical effect study, not experimental validation or proof the closure is unnecessary",
                },
                "research_basis": (
                    "Liu and Stevens, Wall modeled immersed boundary method "
                    "for high Reynolds number flow over complex terrain, "
                    "Computers and Fluids 208 (2020) 104604, "
                    "doi:10.1016/j.compfluid.2020.104604"
                ),
                "next_required": (
                    "Run the updated closure through grid- and time-refined "
                    "flat-wall friction and wall-mounted-cube wake benchmarks; "
                    "quantify surface-roughness sensitivity before promotion."
                ),
                "interpretation": (
                    "A geometry-aware smooth-wall shear model is implemented "
                    "and has a conservative 20-s coarse execution witness, "
                    "but the existing 600-s screening run predates this "
                    "implementation. Wall shear and separated-wake behavior "
                    "remain unvalidated; do not use the screening wake as "
                    "physical validation."
                ),
            },
        },
        {
            "id": "conservative_obstacle_solver",
            "status": "pass",
            "evidence": "screening.physics.conservation_residuals",
            "details": screening["physics"]["conservation_residuals"],
        },
        {
            "id": "nonisothermal_obstacle_inlet_component",
            "status": "conditional",
            "evidence": "AIJ Case J checksum-pinned package, measured thermal inlet, source-aware binary-mixture FV smoke, and pointwise scorer",
            "details": {
                "dataset_status": case_j["status"],
                "dataset_scope": case_j["dataset_scope"],
                "row_counts": case_j["row_counts"],
                "height_varying_temperature_profile_supported": True,
                "initial_thermal_stratification_supported": True,
                "conservative_inlet_enthalpy_flux_supported": True,
                "binary_mixture_source_species_and_enthalpy_coupled_to_projection": True,
                "source_total_to_species_mass_ratio_supported": True,
                "finite_window_transient_mean_output_supported": True,
                "heated_ground_boundary_supported": True,
                "heated_obstacle_surface_boundary_supported": True,
                "obstacle_wall_temperature_representation": "uniform_scalar_temperature_first_order_half_cell_conduction",
                "spatially_resolved_measured_wall_temperature_applied": False,
                "same_condition_model_comparison": "development_run_scored_validation_withheld",
                "hydrogen_validation": False,
            },
        },
        {
            "id": "front_top_side_wake_vol_percent_fields",
            "status": "conditional",
            "evidence": "screening.obstacle_regions, front/top/side/wake figure, and matched 18-s wall-representation sensitivity",
            "details": {
                "regions": list(screening["obstacle_regions"]),
                "unit": screening["concentration_unit"],
                "18s_wall_representation_sensitivity": {
                    "status": geometry_18s.get("status"),
                    "common_fluid_endpoint_rmse_vol_percentage_points": geometry_18s.get("common_fluid_field_difference_vol_percentage_points", {}).get("endpoint_rmse"),
                    "endpoint_max_abs_difference_vol_percentage_points": geometry_18s.get("common_fluid_field_difference_vol_percentage_points", {}).get("endpoint_max_abs"),
                    "endpoint_4vol_percent_volume_m3": geometry_18s.get("four_vol_percent_volume_m3"),
                    "interpretation": "wall representation measurably changes the computed-window field/hazard; no accuracy inference without observations",
                },
                "19s_wall_representation_sensitivity": {
                    "status": geometry_19s.get("status"),
                    "common_fluid_endpoint_rmse_vol_percentage_points": geometry_19s.get("common_fluid_field_difference_vol_percentage_points", {}).get("endpoint_rmse"),
                    "endpoint_4vol_percent_volume_m3": geometry_19s.get("four_vol_percent_volume_m3"),
                    "interpretation": "same-grid wall-representation sensitivity only; not independent validation",
                },
                "20s_wall_representation_sensitivity": {
                    "status": geometry_20s.get("status"),
                    "common_fluid_endpoint_rmse_vol_percentage_points": geometry_20s.get("common_fluid_field_difference_vol_percentage_points", {}).get("endpoint_rmse"),
                    "endpoint_4vol_percent_volume_m3": geometry_20s.get("four_vol_percent_volume_m3"),
                    "interpretation": "same-grid 20-s wall-representation sensitivity only; not independent validation",
                },
            },
        },
        {
            "id": "four_vol_percent_hazard_region",
            "status": "conditional",
            "evidence": "screening.hazard and committed zero-surface 18-s prefix",
            "details": {
                "screening": screening["hazard"],
                "zero_surface_18s": {
                    "prefix_status": zero_surface_18s.get("status"),
                    "stop_time_s": zero_surface_18s.get("configuration", {}).get("duration_s"),
                    "final_4vol_percent_volume_m3": zero_surface_18s_evolution.get("later", {}).get("hazard_volume_final_m3"),
                    "historical_4vol_percent_envelope_m3": zero_surface_18s_evolution.get("later", {}).get("hazard_volume_ever_m3"),
                    "maximum_vol_percent": zero_surface_18s.get("numerics", {}).get("maximum_vol_percent"),
                    "evolution_2s_to18s": zero_surface_18s_evolution.get("delta_later_minus_earlier"),
                    "interpretation": "computed-window evolution only; not full-event hazard validation or convergence",
                },
                "zero_surface_20s": {
                    "prefix_status": zero_surface_20s.get("status"),
                    "stop_time_s": zero_surface_20s.get("configuration", {}).get("duration_s"),
                    "final_4vol_percent_volume_m3": geometry_20s.get("four_vol_percent_volume_m3", {}).get("zero_thickness_surface_endpoint"),
                    "historical_4vol_percent_envelope_m3": geometry_20s.get("four_vol_percent_volume_m3", {}).get("zero_thickness_surface_ever"),
                    "maximum_vol_percent": zero_surface_20s.get("numerics", {}).get("maximum_vol_percent"),
                    "interpretation": "completed 20-s computed window; full TI021 event horizon, hazard convergence, and independent validation remain unmet",
                },
            },
        },
        {
            "id": "sensor_time_series_and_error_metrics",
            "status": "conditional",
            "evidence": "screening.receptors and screening.independent_error_metrics",
            "details": {
                "full_horizon_screening": screening["independent_error_metrics"],
                "full_horizon_600s_TI021_score": {
                    "status": screening_score.get("status"),
                    "decision": screening_score.get("decision"),
                    "normalized_dose_bias": screening_score.get("dose_sensor_score", {}).get("normalized_mean_bias"),
                    "dose_MAPE": screening_score.get("dose_sensor_score", {}).get("mean_absolute_percentage_error"),
                    "detected_positive_sensors": screening_score.get("dose_sensor_score", {}).get("predicted_detected_count_at_0p1vol_percent"),
                    "observed_positive_sensors": screening_score.get("dose_sensor_score", {}).get("observed_positive_count"),
                    "missed_arrivals": screening_score.get("arrival_score", {}).get("missed_count"),
                    "arrival_MAE_s": screening_score.get("arrival_score", {}).get("mean_absolute_error_s"),
                    "six_arc_peak_MAPE": screening_score.get("six_arc_score", {}).get("peak_mean_absolute_percentage_error"),
                    "six_arc_observed_peak_vol_percent": screening_score.get("six_arc_score", {}).get("observed_peak_vol_percent"),
                    "six_arc_predicted_peak_vol_percent": screening_score.get("six_arc_score", {}).get("predicted_peak_0p6s_average_vol_percent"),
                    "four_vol_percent_sensitivity_tested": False,
                    "interpretation": "full observation horizon is covered, but coarse grid/receptor aliasing and poor dose/peak errors withhold quantitative validation",
                },
                "current_nested_145_sensor_600s_score": {
                    "status": nested_145_score.get("status"),
                    "decision": nested_145_score.get("decision", {}).get("result"),
                    "inner_sensor_count": nested_145_score.get("coverage", {}).get("inner_sensor_count"),
                    "outer_sensor_count": nested_145_score.get("coverage", {}).get("outer_sensor_count"),
                    "unique_sensor_count": nested_145_score.get("coverage", {}).get("combined_unique_sensor_count"),
                    "observed_positive_dose_sensor_count": nested_145_score.get("dose_sensor_score", {}).get("observed_positive_count"),
                    "predicted_detected_sensor_count": nested_145_score.get("dose_sensor_score", {}).get("predicted_detected_count_at_0p1vol_percent"),
                    "normalized_dose_bias": nested_145_score.get("dose_sensor_score", {}).get("normalized_mean_bias"),
                    "dose_MAPE": nested_145_score.get("dose_sensor_score", {}).get("mean_absolute_percentage_error"),
                    "arrival_MAE_s": nested_145_score.get("arrival_score", {}).get("mean_absolute_error_s"),
                    "departure_MAE_s": nested_145_score.get("departure_score", {}).get("mean_absolute_error_s"),
                    "overlap_sensor_count": nested_145_score.get("coverage", {}).get("overlap_sensor_count"),
                    "overlap_gate_passed": nested_145_score.get("overlap_audit", {}).get("all_gates_pass"),
                    "maximum_overlap_dose_change": nested_145_score.get("overlap_audit", {}).get("maximum_dose_absolute_relative_change"),
                    "interior_overlap_location_count": nested_145_score.get("overlap_audit", {}).get("geometry", {}).get("interior_overlap_horizontal_location_count"),
                    "boundary_coincident_sensor_count": len(nested_145_score.get("overlap_audit", {}).get("geometry", {}).get("boundary_coincident_sensor_labels", [])),
                    "bottom_cell_centre_m": nested_145_score.get("vertical_sensor_resolution_diagnostic", {}).get("bottom_cell_centre_m"),
                    "receptors_at_or_below_bottom_cell_centre_count": nested_145_score.get("vertical_sensor_resolution_diagnostic", {}).get("receptors_at_or_below_bottom_cell_centre_count"),
                    "interpretation": "corrected full-horizon nested diagnostic covers all sensors and has real interior overlap; peak coupling passes but integrated-dose coupling fails and experiment dose error remains gross, so quantitative promotion is withheld",
                },
                "short_true_interior_overlap_diagnostic": {
                    "status": interior_overlap.get("status"),
                    "model_horizon_s": interior_overlap.get("configuration", {}).get("model_horizon_s"),
                    "handoff_site_downwind_y_m": interior_overlap.get("configuration", {}).get("handoff_site_downwind_y_m"),
                    "inner_sensor_count": interior_overlap.get("coverage", {}).get("inner_sensor_count"),
                    "outer_sensor_count": interior_overlap.get("coverage", {}).get("outer_sensor_count"),
                    "shared_sensor_count": interior_overlap.get("coverage", {}).get("shared_sensor_count"),
                    "union_sensor_count": interior_overlap.get("coverage", {}).get("union_sensor_count"),
                    "expected_sensor_count": interior_overlap.get("coverage", {}).get("expected_sensor_count"),
                    "interior_overlap_horizontal_location_count": interior_overlap.get("geometry", {}).get("interior_overlap_horizontal_location_count"),
                    "boundary_coincident_sensor_count": len(interior_overlap.get("geometry", {}).get("boundary_coincident_sensor_labels", [])),
                    "maximum_peak_relative_change": interior_overlap.get("comparison", {}).get("maximum_peak_relative_change_using_0p1vol_floor"),
                    "maximum_dose_relative_change": interior_overlap.get("comparison", {}).get("maximum_dose_relative_change_using_0p1vol_s_floor"),
                    "inner_runtime_s": interior_overlap.get("runtime", {}).get("inner_wall_clock_s"),
                    "outer_runtime_s": interior_overlap.get("runtime", {}).get("outer_wall_clock_s"),
                    "experimental_model_promotion_allowed": interior_overlap.get("gates", {}).get("experimental_model_promotion_allowed"),
                    "interpretation": "true shared interior sensor comparison is now executable and within the diagnostic 10% limit over 60 s; the complete 145-sensor union is covered, but the full event and physical-model validation remain unproven",
                },
                "latest_nested_145_sensor_60s_score": {
                    "status": latest_nested_60s_score.get("status"),
                    "model_horizon_s": latest_nested_60s_score.get("model_horizon_s"),
                    "unique_sensor_count": latest_nested_60s_score.get("coverage", {}).get("combined_unique_sensor_count"),
                    "interior_overlap_sensor_count": latest_nested_60s_score.get("overlap_audit", {}).get("geometry", {}).get("interior_overlap_sensor_count"),
                    "interior_overlap_location_count": latest_nested_60s_score.get("overlap_audit", {}).get("geometry", {}).get("interior_overlap_horizontal_location_count"),
                    "interior_overlap_gate_passed": latest_nested_60s_score.get("overlap_audit", {}).get("all_gates_pass"),
                    "maximum_interior_peak_change": latest_nested_60s_score.get("overlap_audit", {}).get("maximum_peak_absolute_relative_change"),
                    "maximum_interior_dose_change": latest_nested_60s_score.get("overlap_audit", {}).get("maximum_dose_absolute_relative_change"),
                    "observed_positive_dose_count": latest_nested_60s_score.get("dose_sensor_score", {}).get("observed_positive_count"),
                    "dose_MAPE": latest_nested_60s_score.get("dose_sensor_score", {}).get("mean_absolute_percentage_error"),
                    "normalized_dose_bias": latest_nested_60s_score.get("dose_sensor_score", {}).get("normalized_mean_bias"),
                    "decision": latest_nested_60s_score.get("decision", {}).get("result"),
                    "interpretation": "60 s diagnostic only; this does not score the full experimental event or establish grid convergence",
                },
                "corrected_nested_600s_inner_run": {
                    "status": corrected_nested_600s_inner.get("status"),
                    "all_internal_gates_pass": corrected_nested_600s_inner.get("all_gates_pass"),
                    "model_horizon_s": corrected_nested_600s_inner.get("configuration", {}).get("duration_s"),
                    "inner_receptor_count": corrected_nested_600s_inner.get("configuration", {}).get("inner_receptor_count"),
                    "domain_m": corrected_nested_600s_inner.get("configuration", {}).get("domain_m"),
                    "handoff_site_downwind_y_m": corrected_nested_600s_inner.get("configuration", {}).get("handoff_site_downwind_y_m"),
                    "maximum_relative_species_mass_residual": corrected_nested_600s_inner.get("numerics", {}).get("maximum_relative_species_mass_residual"),
                    "maximum_relative_sensible_enthalpy_residual": corrected_nested_600s_inner.get("numerics", {}).get("maximum_relative_sensible_enthalpy_residual"),
                    "maximum_continuity_residual_kg_m3_s": corrected_nested_600s_inner.get("numerics", {}).get("maximum_continuity_residual_kg_m3_s"),
                    "maximum_global_mass_residual_kg_s": corrected_nested_600s_inner.get("numerics", {}).get("maximum_global_mass_residual_kg_s"),
                    "corrector_nonconverged_steps": corrected_nested_600s_inner.get("numerics", {}).get("corrector_nonconverged_steps"),
                    "runtime_s": corrected_nested_600s_inner.get("runtime", {}).get("wall_clock_s"),
                    "outer_600s_run_completed": corrected_nested_600s_outer.get("status") == "accepted_completed_outer_transient" and corrected_nested_600s_outer.get("all_gates_pass") is True,
                    "outer_runtime_s": corrected_nested_600s_outer.get("runtime", {}).get("wall_clock_s"),
                    "outer_relative_species_mass_residual": corrected_nested_600s_outer.get("numerics", {}).get("cumulative_relative_species_mass_residual"),
                    "outer_maximum_relative_remap_residual": corrected_nested_600s_outer.get("numerics", {}).get("maximum_relative_integral_remap_residual"),
                    "interpretation": "full-horizon conservative inner and outer calculations completed; the nested dose-coupling gate still fails and this does not establish experimental validation",
                },
                "zero_surface_18s_horizon_limited_score": {
                    "status": zero_surface_18s_score.get("status"),
                    "decision": zero_surface_18s_score.get("decision"),
                    "sensor_count": zero_surface_18s_score.get("dose_sensor_score", {}).get("sensor_count"),
                    "observed_positive_count": zero_surface_18s_score.get("dose_sensor_score", {}).get("observed_positive_count"),
                    "predicted_detected_count": zero_surface_18s_score.get("dose_sensor_score", {}).get("predicted_detected_count_at_0p1vol_percent"),
                    "missed_arrival_count": zero_surface_18s_score.get("arrival_score", {}).get("missed_count"),
                    "normalized_mean_bias": zero_surface_18s_score.get("dose_sensor_score", {}).get("normalized_mean_bias"),
                    "dose_mape": zero_surface_18s_score.get("dose_sensor_score", {}).get("mean_absolute_percentage_error"),
                    "observed_six_arc_peak_time_range_s": [min(zero_surface_18s_score["six_arc_score"]["observed_peak_time_s"]), max(zero_surface_18s_score["six_arc_score"]["observed_peak_time_s"])],
                    "model_horizon_s": zero_surface_18s.get("configuration", {}).get("duration_s"),
                    "interpretation": "Run stops before most observed TI021 peaks; dose/arrival score is horizon-limited and withheld, not a valid final accuracy estimate.",
                },
                "zero_surface_19s_horizon_limited_score": {
                    "status": zero_surface_19s_score.get("status"),
                    "model_horizon_s": zero_surface_19s.get("configuration", {}).get("duration_s"),
                    "observed_positive_count": zero_surface_19s_score.get("dose_sensor_score", {}).get("observed_positive_count"),
                    "predicted_detected_count": zero_surface_19s_score.get("dose_sensor_score", {}).get("predicted_detected_count_at_0p1vol_percent"),
                    "missed_arrival_count": zero_surface_19s_score.get("arrival_score", {}).get("missed_count"),
                    "arrival_MAE_s": zero_surface_19s_score.get("arrival_score", {}).get("mean_absolute_error_s"),
                    "decision": zero_surface_19s_score.get("decision"),
                },
                "zero_surface_20s_horizon_limited_score": {
                    "status": zero_surface_20s_score.get("status"),
                    "decision": zero_surface_20s_score.get("decision"),
                    "model_horizon_s": zero_surface_20s.get("configuration", {}).get("duration_s"),
                    "observed_positive_count": zero_surface_20s_score.get("dose_sensor_score", {}).get("observed_positive_count"),
                    "predicted_detected_count": zero_surface_20s_score.get("dose_sensor_score", {}).get("predicted_detected_count_at_0p1vol_percent"),
                    "normalized_mean_bias": zero_surface_20s_score.get("dose_sensor_score", {}).get("normalized_mean_bias"),
                    "dose_mape": zero_surface_20s_score.get("dose_sensor_score", {}).get("mean_absolute_percentage_error"),
                    "arrival_mae_s": zero_surface_20s_score.get("arrival_score", {}).get("mean_absolute_error_s"),
                    "decision_boundary": zero_surface_20s_score.get("decision"),
                    "interpretation": "score is horizon-limited and withheld; observed TI021 peak window extends well beyond 20 s",
                },
            },
        },
        {
            "id": "grid_and_time_convergence",
            "status": "withheld",
            "evidence": "three_grid_impact_convergence",
            "details": {
                "matched_window_temporal_sensitivity": {
                    "stop_time_s": temporal["receptors"]["common_stop_time_s"],
                    "maximum_field_rmse_vol_percent": temporal["field"]["rmse"],
                    "receptor_rmse_vol_percent": temporal["receptors"]["rmse"],
                    "four_vol_percent_envelope_volume_difference_m3": abs(temporal["field"]["hazard_volume_difference_m3"]),
                    "input_sha256_recorded": bool(temporal.get("provenance")),
                },
                "matched_dt_obstacle_arc_length_relative_error": {
                    label: value["relative_arc_length_error"]
                    for label, value in grid["voxelized_obstacle_geometry"]["by_grid"].items()
                },
                "field_contract": grid["gates"].get("field_l1_and_rmse_contract"),
                "hazard_gate": grid["gates"].get("hazard_volume_monotonic_gci_le_10_percent"),
                "receptor_gate": grid["gates"].get("receptor_peak_and_dose_monotonic_fraction_ge_70_percent"),
                "full_event_promoted": grid.get("full_event_grid_convergence_promoted"),
                "short_fixed_surface_prefix_diagnostic": {
                    "window_classification": grid_short["configuration"]["window_classification"],
                    "duration_s": grid_short["configuration"]["duration_s"],
                    "fixed_numerical_surface": grid_short["gates"]["numerical_obstacle_representation_held_constant"],
                    "field_norms_contract": grid_short["gates"]["field_l1_and_rmse_contract"],
                    "hazard_volume_gci_gate": grid_short["gates"]["hazard_volume_monotonic_gci_le_10_percent"],
                    "receptor_gate": grid_short["gates"]["receptor_peak_and_dose_monotonic_fraction_ge_70_percent"],
                    "full_20s_gate": grid_short["gates"]["frozen_observation_blind_20s_window"],
                },
                "fixed_surface_20s_grid_diagnostic": {
                    "window_classification": grid_zero_surface_20s["configuration"]["window_classification"],
                    "duration_s": grid_zero_surface_20s["configuration"]["duration_s"],
                    "field_norms_contract": grid_zero_surface_20s["gates"]["field_l1_and_rmse_contract"],
                    "hazard_volume_gci_gate": grid_zero_surface_20s["gates"]["hazard_volume_monotonic_gci_le_10_percent"],
                    "peak_monotonic_fraction": grid_zero_surface_20s["receptor_convergence"]["peak_vol_percent"]["monotonic_convergent_fraction"],
                    "dose_monotonic_fraction": grid_zero_surface_20s["receptor_convergence"]["dose_vol_percent_s"]["monotonic_convergent_fraction"],
                    "peak_median_fine_gci_percent": grid_zero_surface_20s["receptor_convergence"]["peak_vol_percent"]["fine_gci_relative_percent_median"],
                    "dose_median_fine_gci_percent": grid_zero_surface_20s["receptor_convergence"]["dose_vol_percent_s"]["fine_gci_relative_percent_median"],
                    "physical_geometry_held_constant": grid_zero_surface_20s["gates"]["physical_obstacle_geometry_identified_and_held_constant"],
                    "full_event_promoted": grid_zero_surface_20s["full_event_grid_convergence_promoted"],
                    "interpretation": "full 0-20 s numerical grid sensitivity; field norms contract but hazard and receptor convergence gates fail, and represented physical geometry is not held constant",
                },
                "ffi_hybrid_patch_interface_audit": {
                    "status": "completed_interface_and_conservation_witness_not_facility_validation",
                    "duration_s": ffi_hybrid_thermal.get("duration_s"),
                    "wall_clock_s": ffi_hybrid_thermal.get("runtime_s"),
                    "time_step_s": ffi_hybrid_thermal.get("time_step_s"),
                    "thermal_enthalpy_exchange": ffi_hybrid_thermal.get("interface_relaxation", {}).get("thermal_enthalpy_exchange"),
                    "incoming_bridge_exchanges": ffi_hybrid_thermal.get("bridge", {}).get("incoming_exchanged_steps"),
                    "outgoing_bridge_exchanges": ffi_hybrid_thermal.get("bridge", {}).get("outgoing_exchanged_steps"),
                    "shared_sensor_count": ffi_hybrid_thermal.get("sensor_count"),
                    "sensor_mae_vol_percent_points": ffi_hybrid_thermal.get("diagnostic_sensor_error_not_validation_score", {}).get("mae_vol_percent"),
                    "sensor_rmse_vol_percent_points": ffi_hybrid_thermal.get("diagnostic_sensor_error_not_validation_score", {}).get("rmse_vol_percent"),
                    "maximum_concentration_vol_percent": ffi_hybrid_thermal.get("maximum_predicted_vol_percent"),
                    "relative_species_mass_residual_kg": ffi_hybrid_thermal.get("maximum_partitioned_mass_balance_residual_kg"),
                    "relative_enthalpy_residual_j": ffi_hybrid_thermal.get("resolved_patch_energy_residual_j"),
                    "maximum_continuity_residual_kg_m3_s": ffi_hybrid_thermal.get("resolved_patch_continuity_residual_kg_m3_s"),
                    "maximum_pressure_correction_pa": ffi_hybrid_thermal.get("resolved_patch_max_abs_pressure_correction_pa"),
                    "maximum_mach": ffi_hybrid_thermal.get("resolved_patch_diagnostics", {}).get("maximum_lattice_mach_number"),
                    "interpretation": "fast-resolved-fast interface and conservation witness; FFI geometry remains conditional and independent obstacle-H2 validation remains withheld",
                },
                "ffi_isothermal_time_step_audit": {
                    "status": ffi_hybrid_temporal.get("status"),
                    "duration_s": ffi_hybrid_temporal.get("duration_s"),
                    "coarse_dt_s": ffi_hybrid_temporal.get("coarse_time_step_s"),
                    "fine_dt_s": ffi_hybrid_temporal.get("fine_time_step_s"),
                    "sensor_mean_absolute_change_vol_percent_points": ffi_hybrid_temporal.get("sensor_mean_delta_vol_percent", {}).get("mean_abs"),
                    "sensor_rmse_change_vol_percent_points": ffi_hybrid_temporal.get("sensor_mean_delta_vol_percent", {}).get("rmse"),
                    "field_rmse_change_vol_percent": ffi_hybrid_temporal.get("field_delta_vol_percent", {}).get("rmse"),
                    "coarse_runtime_s": ffi_hybrid_temporal.get("runtime_s", {}).get("coarse"),
                    "fine_runtime_s": ffi_hybrid_temporal.get("runtime_s", {}).get("fine"),
                    "interpretation": "explicit isothermal audit because active thermal D3Q19 path requires the physical step; does not promote facility prediction",
                },
                "ffi_flow_iteration_sensitivity": {
                    "status": "completed_execution_sensitivity_not_facility_validation",
                    "baseline_runtime_s": ffi_hybrid_thermal.get("runtime_s"),
                    "reduced_runtime_s": ffi_hybrid_flow20.get("runtime_s"),
                    "baseline_flow_steps_per_scalar_step": ffi_hybrid_thermal.get("resolved_flow_steps_per_scalar_step"),
                    "reduced_flow_steps_per_scalar_step": ffi_hybrid_flow20.get("resolved_flow_steps_per_scalar_step"),
                    "runtime_reduction_fraction": 1.0 - ffi_hybrid_flow20.get("runtime_s", 0.0) / ffi_hybrid_thermal.get("runtime_s", 1.0),
                    "interpretation": "optional speed mode; same conditional geometry and not a validation score",
                },
            },
        },
        {
            "id": "aij_neutral_wake_component",
            "status": "conditional",
            "evidence": "aij_case_h_time_averaged_aligned_grid_audit plus instantaneous resolved 3-D diagnostic",
            "details": {
                "time_averaged_aligned_grid_status": aij_mean.get("status"),
                "time_averaged_aligned_grid_passed": aij_mean.get("passed"),
                "common_points": aij_mean.get("common_point_count"),
                "finest_u_correlation": aij_mean["cases"][-1]["matched_metrics"]["u_correlation"],
                "finest_w_correlation": aij_mean["cases"][-1]["matched_metrics"]["w_correlation"],
                "finest_flow_converged": aij_mean["cases"][-1]["flow_converged"],
                "finest_obstacle_grid_aligned": aij_mean["cases"][-1]["obstacle_grid_alignment_exact"],
                "finest_pair_u_mae_relative_change": aij_mean["finest_pair_metric_changes"]["u_mae_normalized"],
                "finest_pair_w_rmse_relative_change": aij_mean["finest_pair_metric_changes"]["w_rmse_normalized"],
                "instantaneous_grid_audit_status": aij.get("status"),
                "instantaneous_grid_audit_passed": aij.get("passed"),
                "quantitative_hydrogen_validation": aij_mean.get("quantitative_hydrogen_validation"),
                "concentration_validation_withheld": aij_mean.get("concentration_validation_withheld"),
                "lateral_domain_sensitivity": {
                    "status": aij_lateral_domain.get("status"),
                    "base_lateral_domain_m": aij_lateral_domain.get("base_domain", {}).get("lateral_domain_m"),
                    "wide_lateral_domain_m": aij_lateral_domain.get("wide_domain", {}).get("lateral_domain_m"),
                    "base_flow_converged": aij_lateral_domain.get("base_domain", {}).get("flow_converged"),
                    "wide_flow_converged": aij_lateral_domain.get("wide_domain", {}).get("flow_converged"),
                    "wide_minus_base_metric_change": aij_lateral_domain.get("wide_minus_base_metric_change"),
                    "calibration_used": aij_lateral_domain.get("calibration_used"),
                    "promotion_allowed": aij_lateral_domain.get("promotion_allowed"),
                },
                "lateral_domain_tracer_sensitivity": {
                    "status": aij_lateral_tracer.get("status"),
                    "common_point_count": aij_lateral_tracer.get("common_point_count"),
                    "base_rmse_ppmv": aij_lateral_tracer.get("base_common_point_metrics", {}).get("rmse_ppmv"),
                    "wide_rmse_ppmv": aij_lateral_tracer.get("wide_common_point_metrics", {}).get("rmse_ppmv"),
                    "rmse_reduction_fraction": aij_lateral_tracer.get("rmse_reduction_fraction"),
                    "base_pearson_r": aij_lateral_tracer.get("base_common_point_metrics", {}).get("pearson_r"),
                    "wide_pearson_r": aij_lateral_tracer.get("wide_common_point_metrics", {}).get("pearson_r"),
                    "calibration_used": aij_lateral_tracer.get("calibration_used"),
                    "promotion_allowed": aij_lateral_tracer.get("promotion_allowed"),
                },
                "passive_tracer_resolution_sensitivity": {
                    "status": aij_tracer_resolution.get("status"),
                    "common_measurement_point_count": aij_tracer_resolution.get("common_measurement_point_count"),
                    "time_window_s": aij_tracer_resolution.get("comparison_window_s"),
                    "coarse_rmse_ppmv": aij_tracer_resolution["coarse"]["metrics_on_common_points"]["rmse_ppmv"],
                    "fine_rmse_ppmv": aij_tracer_resolution["fine"]["metrics_on_common_points"]["rmse_ppmv"],
                    "coarse_bias_ppmv": aij_tracer_resolution["coarse"]["metrics_on_common_points"]["bias_ppmv"],
                    "fine_bias_ppmv": aij_tracer_resolution["fine"]["metrics_on_common_points"]["bias_ppmv"],
                    "coarse_flow_converged": aij_tracer_resolution["coarse"]["flow_converged"],
                    "fine_flow_converged": aij_tracer_resolution["fine"]["flow_converged"],
                    "calibration_used": aij_tracer_resolution.get("calibration_used"),
                    "promotion_allowed": aij_tracer_resolution.get("promotion_allowed"),
                },
                "passive_tracer_sgs_sensitivity": {
                    "status": aij_tracer_sgs.get("status"),
                    "common_measurement_point_count": aij_tracer_sgs.get("common_measurement_point_count"),
                    "comparison_window_s": aij_tracer_sgs.get("comparison_window_s"),
                    "molecular_only_rmse_ppmv": aij_tracer_sgs["baseline_molecular_only"]["scores"]["rmse_ppmv"],
                    "sgs_mixed_rmse_ppmv": aij_tracer_sgs["sgs_gradient_diffusion_sensitivity"]["scores"]["rmse_ppmv"],
                    "rmse_reduction_fraction": aij_tracer_sgs.get("rmse_reduction_fraction"),
                    "molecular_only_bias_ppmv": aij_tracer_sgs["baseline_molecular_only"]["scores"]["bias_ppmv"],
                    "sgs_mixed_bias_ppmv": aij_tracer_sgs["sgs_gradient_diffusion_sensitivity"]["scores"]["bias_ppmv"],
                    "sgs_effective_schmidt_number": aij_tracer_sgs["sgs_gradient_diffusion_sensitivity"]["configuration"].get("effective_schmidt_number"),
                    "sgs_mass_residual": aij_tracer_sgs["sgs_gradient_diffusion_sensitivity"].get("mass_residual"),
                    "calibration_used": aij_tracer_sgs.get("calibration_used"),
                    "promotion_allowed": aij_tracer_sgs.get("promotion_allowed"),
                },
                "ethylene_air_fv_source_momentum_smoke": {
                    "status": aij_source_momentum.get("status"),
                    "averaging_window_s": aij_source_momentum.get("averaging_window_s"),
                    "source_exit_speed_m_s_inferred": aij_source_momentum.get("source_exit_speed_m_s"),
                    "continuity_residual_max_kg_m3_s": aij_source_momentum.get("finite_volume_continuity_residual_max_kg_m3_s"),
                    "global_mass_residual_max_kg_s": aij_source_momentum.get("finite_volume_global_mass_residual_max_kg_s"),
                    "momentum_balance_residual_max_kg_m_s": aij_source_momentum.get("momentum_balance_residual_max_kg_m_s"),
                    "paired_mean_field_change_ppmv": aij_source_momentum.get("paired_field_change", {}).get("mean_delta_ppmv"),
                    "paired_field_rmse_change_ppmv": aij_source_momentum.get("paired_field_change", {}).get("rmse_delta_ppmv"),
                    "without_momentum_rmse_ppmv": aij_source_momentum.get("without_momentum", {}).get("rmse_ppmv"),
                    "with_momentum_rmse_ppmv": aij_source_momentum.get("with_inferred_momentum", {}).get("rmse_ppmv"),
                    "relative_rmse_change_fraction": (
                        (aij_source_momentum["without_momentum"]["rmse_ppmv"]
                         - aij_source_momentum["with_inferred_momentum"]["rmse_ppmv"])
                        / aij_source_momentum["without_momentum"]["rmse_ppmv"]
                    ),
                    "time_window_is_experimental_120s": aij_source_momentum.get("averaging_window_s") == [0.0, 120.0],
                    "score_is_120s_validation": False,
                    "promotion_allowed": aij_source_momentum.get("promotion_allowed"),
                },
                "matched_d3q19_vs_binary_fv_tracer_tiers": {
                    "status": aij_tier_comparison.get("status"),
                    "common_sensor_count": aij_tier_comparison.get("matched_conditions", {}).get("comparison_ids"),
                    "mean_window_s": aij_tier_comparison.get("matched_conditions", {}).get("mean_window_s"),
                    "D3Q19_RMSE_ppmv": aij_tier_comparison.get("D3Q19_molecular_air_carrier_with_SGS_scalar", {}).get("scores", {}).get("rmse_ppmv"),
                    "FV_RMSE_ppmv": aij_tier_comparison.get("finite_volume_binary_ethylene_air_with_inferred_nozzle_momentum", {}).get("scores", {}).get("rmse_ppmv"),
                    "RMSE_reduction_fraction": aij_tier_comparison.get("rmse_reduction_fraction"),
                    "FV_continuity_residual_kg_m3_s": aij_tier_comparison.get("finite_volume_binary_ethylene_air_with_inferred_nozzle_momentum", {}).get("continuity_residual_kg_m3_s"),
                    "FV_global_mass_residual_kg_s": aij_tier_comparison.get("finite_volume_binary_ethylene_air_with_inferred_nozzle_momentum", {}).get("global_mass_residual_kg_s"),
                    "D3Q19_runtime_s": aij_tier_comparison.get("D3Q19_molecular_air_carrier_with_SGS_scalar", {}).get("runtime_wall_s"),
                    "FV_runtime_s": aij_tier_comparison.get("finite_volume_binary_ethylene_air_with_inferred_nozzle_momentum", {}).get("runtime_wall_s"),
                    "measured_window_matched": False,
                    "calibration_used": aij_tier_comparison.get("calibration_used"),
                    "promotion_allowed": aij_tier_comparison.get("promotion_allowed"),
                },
            },
        },
        {
            "id": "low_temperature_density_buoyancy_path",
            "status": "conditional",
            "evidence": "obstacle_thermal_buoyancy_comparison and hydrogen_air_buoyancy_transition_audit",
            "details": {
                "maximum_h2_vol_percent": thermal["active_thermal"]["maximum_h2_vol_percent"],
                "minimum_temperature_k": thermal["active_thermal"]["minimum_temperature_k"],
                "maximum_density_departure": thermal["active_thermal"]["diagnostics"]["maximum_relative_mixture_density_departure"],
                "quantitative_prediction_allowed": thermal["active_thermal"]["diagnostics"]["quantitative_prediction_allowed"],
            },
        },
        {
            "id": "public_obstacle_h2_validation",
            "status": "withheld",
            "evidence": "h2safe, exjet, bos and ELVHYS fail-closed audits",
            "details": {
                "h2safe": h2safe.get("status"),
                "exjet": exjet.get("status"),
                "bos": bos.get("status"),
                "dataverse_no": dataverse.get("status"),
                "elvhys_tcs": elvhys.get("status"),
                "elvhys_tcs_open_field_wake_qualified": False,
                "ffi_field_evidence_intake": {
                    "O1": {
                        "status": ffi_o1_intake.get("status"),
                        "gate_codes": ffi_o1_intake.get("gate_codes", []),
                        "source_channel": ffi_o1_intake.get("channels", {}).get("source_boundary", {}).get("status"),
                        "receptor_channel": ffi_o1_intake.get("channels", {}).get("receptor_observations", {}).get("status"),
                        "weather_channel": ffi_o1_intake.get("channels", {}).get("weather", {}).get("status"),
                        "common_clock_channel": ffi_o1_intake.get("channels", {}).get("common_clock", {}).get("status"),
                    },
                    "O2": {
                        "status": ffi_o2_intake.get("status"),
                        "gate_codes": ffi_o2_intake.get("gate_codes", []),
                        "source_channel": ffi_o2_intake.get("channels", {}).get("source_boundary", {}).get("status"),
                        "receptor_channel": ffi_o2_intake.get("channels", {}).get("receptor_observations", {}).get("status"),
                        "weather_channel": ffi_o2_intake.get("channels", {}).get("weather", {}).get("status"),
                        "common_clock_channel": ffi_o2_intake.get("channels", {}).get("common_clock", {}).get("status"),
                    },
                },
            },
        },
        {
            "id": "reproducible_evidence_index",
            "status": "pass" if index.get("all_files_present") else "withheld",
            "evidence": "degali-validation-evidence-index",
            "details": {"all_files_present": index.get("all_files_present"), "record_count": len(index.get("records", []))},
        },
    ]
    status_counts = {}
    for item in requirements:
        status_counts[item["status"]] = status_counts.get(item["status"], 0) + 1
    return {
        "schema": "degali.goal-completion-audit.v1",
        "objective": "validated multifidelity transient hydrogen concentration fields around real-facility obstacles",
        "requirements": requirements,
        "status_counts": status_counts,
        "promotion_allowed": False,
        "completion_status": "conditional_research_and_screening_system",
        "blocking_requirements": [item["id"] for item in requirements if item["status"] == "withheld"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status_counts": result["status_counts"], "promotion_allowed": result["promotion_allowed"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
