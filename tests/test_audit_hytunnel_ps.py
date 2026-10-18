import json
from pathlib import Path
import runpy

import pytest


ROOT = Path("reference/datasets/hytunnel-ps-4668554/raw")


@pytest.mark.skipif(
    not (ROOT / "HTE2440PS001MIXED_300321.xlsx").is_file(),
    reason="requires the locally controlled public HyTunnel data",
)
def test_public_hytunnel_audit_preserves_metadata_gaps(tmp_path):
    tool = runpy.run_path("tools/audit_hytunnel_ps.py")
    record = tool["audit_hytunnel_ps"](ROOT)

    assert record["status"] == "component_evidence_conditional"
    assert record["all_published_md5_verified"]
    assert record["concentration"]["workbook_sheet_count"] == 42
    assert record["concentration"]["hydrogen_labeled_sheet_count"] == 41
    assert record["concentration"]["helium_labeled_sheet_count"] == 1
    assert record["concentration"]["point_count"] == 15_768
    assert record["concentration"]["metadata_configuration_count_discrepancy"]
    reconciliation = record["concentration"]["test_matrix_reconciliation"]
    assert reconciliation["reconciled"]
    assert reconciliation["matched_declared_hydrogen_case_count"] == 39
    assert reconciliation["missing_declared_hydrogen_cases"] == [{
        "nozzle_diameter_mm": 1.0,
        "nominal_mass_flow_g_s": 1.5,
        "nominal_ventilation_speed_m_s": 3.5,
        "flow_direction": "cross-flow",
        "gas_label": "H2",
    }]
    assert len(reconciliation["additional_hydrogen_baselines"]) == 2
    assert reconciliation["non_hydrogen_cases"][0]["gas_label"] == "He"
    assert record["ventilation_flow"]["point_count"] == 1_260
    assert record["promotion_allowed"] is False
    matched = record["matched_reference_case"]
    assert matched["official_cfd_domain_m"] == [2.5, 0.9, 0.9]
    assert matched["published_experiment_4volpct_length_m"] == 1.65
    assert matched["degali_score_status"].startswith("withheld_")
    adapter = matched["degali_source_adapter"]
    assert adapter["status"] == (
        "implemented_profile_preflight_nominal_failed_0p06_candidate_only"
    )
    assert adapter["source_model_match_to_official_cfd"] is False
    assert adapter["degali_notional_nozzle"]["diameter_m"] == pytest.approx(
        0.004363832229758349
    )
    assert adapter["official_cfd_notional_nozzle"]["diameter_m"] == 0.00725
    assert adapter["low_mach_handoff"]["mach_number"] < 0.30
    assert adapter["low_mach_handoff"][
        "relative_density_departure"
    ] == pytest.approx(-0.10)
    assert adapter["low_mach_handoff"][
        "maximum_relative_conservation_residual"
    ] < 1.0e-12
    assert adapter["resolved_source_mapping"]["spatially_resolved"] is True
    assert adapter["resolved_source_mapping"][
        "cells_across_equivalent_diameter"
    ] > 29.0
    assert adapter["resolved_source_mapping"][
        "mapped_gaussian_profile_density_admissibility"
    ] == "nominal_failed_0p06_candidate_only_not_promoted"
    assert "nominal_resolved_source_profile_preflight_failed" in record[
        "gate_codes"
    ]
    assert "matched_high_pressure_source_adapter_not_implemented" not in record[
        "gate_codes"
    ]
    boundary = matched["degali_boundary_configuration"]
    assert boundary["flow_boundary_scheme"] == (
        "characteristic_multidimensional_extrapolated"
    )
    assert boundary["nonstreamwise_wall_boundary"] == "all_pressure_open"
    assert matched["official_cfd_boundary_configuration"][
        "remaining_five_faces"
    ] == "open non-reflecting"
    assert "matched_pressure_open_unified_boundary_candidate_failed" in record[
        "gate_codes"
    ]
    assert "matched_full_duration_domain_validation_required" not in record[
        "gate_codes"
    ]
    assert "matched_official_duration_runtime_target_failed" in record[
        "gate_codes"
    ]
    assert "matched_resolved_concentration_score_failed" in record[
        "gate_codes"
    ]
    assert boundary["streamwise_linear_acoustic_reflection_gate_passed"] is True
    assert boundary[
        "lateral_normal_linear_acoustic_reflection_gate_passed"
    ] is True
    assert boundary[
        "oblique_30deg_linear_acoustic_reflection_gate_passed"
    ] is True
    assert boundary["oblique_angle_wavelength_sweep_gate_passed"] is True
    assert boundary["angle_wavelength_worst_reflection_ratio"] < 0.05
    assert boundary["canonical_xy_corner_reflection_gate_passed"] is True
    assert boundary["canonical_xy_corner_reflection_ratio"] < 0.05
    assert boundary[
        "three_dimensional_corner_orientation_validation_required"
    ] is False
    assert boundary["extrapolated_wavelength8_3d_corner_gate_passed"] is True
    assert boundary["extrapolated_wavelength10_3d_corner_gate_passed"] is True
    assert boundary["multiwavelength_3d_corner_validation_required"] is False
    assert boundary[
        "extrapolated_wavelength8_xyz_corner_reflection_ratio"
    ] < 0.05
    assert boundary[
        "extrapolated_wavelength10_xyz_corner_reflection_ratio"
    ] < 0.05
    assert boundary["predecessor_wavelength10_3d_corner_gate_passed"] is False
    assert boundary["vortical_outflow_boundary_validation_required"] is True
    assert boundary["vortical_limited_candidate_stable"] is True
    assert boundary["vortical_limited_candidate_gate_passed"] is False
    assert boundary["vortical_limited_maximum_normalised_velocity_l2"] > 0.05
    assert boundary["vortical_limited_maximum_normalised_vorticity_l2"] < 0.05
    assert boundary["vortical_limited_wavelength8_3d_corner_gate_passed"] is True
    assert boundary["vortical_limited_wavelength10_3d_corner_gate_passed"] is False
    assert boundary["unified_pressure_open_boundary_gate_passed"] is False
    assert boundary["short_horizon_plume_domain_sensitivity_passed"] is True
    assert boundary["short_horizon_maximum_primary_relative_change"] < 0.10
    assert boundary["short_horizon_runtime_target_passed"] is True
    assert boundary["short_horizon_expanded_runtime_s"] < 300.0
    assert boundary["short_horizon_projection_linear_solver"] == "cached_direct"
    assert boundary["official_0p707s_domain_sensitivity_required"] is False
    assert boundary["official_0p707s_domain_sensitivity_passed"] is True
    assert boundary["official_0p707s_maximum_primary_relative_change"] < 0.10
    assert boundary["official_0p707s_common_field_normalised_l1"] < 0.10
    assert boundary["official_0p707s_expanded_boundary_clearance_passed"] is True
    assert boundary["official_0p707s_runtime_target_passed"] is False
    assert boundary["official_0p707s_expanded_runtime_s"] > 300.0
    score = matched["degali_resolved_overlap_score"]
    assert score["supported_downstream_point_count"] == 54
    assert score["operator_match"] is False
    assert score["qualified_steady_endpoint_proxy"] is True
    assert score["maximum_slice_replay_normalised_l1"] == 0.0
    assert score["mae_vol_pct_point"] == pytest.approx(2.1792203023)
    assert score["four_vol_sensitivity"] == 1.0
    assert score["four_vol_specificity"] < 0.90
    sensitivity = matched["degali_source_mixing_sensitivity"]
    assert sensitivity["common_supported_point_count"] == 45
    assert sensitivity["calibration_or_best_variant_selection_allowed"] is False
    assert sensitivity["handoff0p06_sc0p5_mae_vol_pct_point"] < sensitivity[
        "baseline_handoff0p06_sc0p7_mae_vol_pct_point"
    ]
    assert sensitivity["handoff0p06_sc0p5_four_vol_specificity"] < 0.90
    local_sgs = matched["degali_local_sgs_closure_candidate"]
    assert local_sgs["common_supported_point_count"] == 45
    assert local_sgs["maximum_relative_density_departure"] < 0.10
    assert local_sgs["maximum_relative_mass_residual"] < 1.0e-6
    assert local_sgs["mae_vol_pct_point"] > sensitivity[
        "handoff0p05_sc0p7_mae_vol_pct_point"
    ]
    assert local_sgs["four_vol_specificity"] == pytest.approx(0.25)
    assert local_sgs["calibration_or_best_variant_selection_allowed"] is False
    assert "matched_scalar_mixing_external_basis_required" in record[
        "gate_codes"
    ]
    assert "local_sgs_scalar_closure_accuracy_failed" in record["gate_codes"]

    output = tmp_path / "hytunnel-audit.json"
    tool["write_audit"](record, output)
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded["cad"]["unsafe_path_count"] == 0
    with pytest.raises(FileExistsError):
        tool["write_audit"](record, output)
