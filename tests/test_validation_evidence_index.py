from __future__ import annotations

from tools.build_validation_evidence_index import build


def test_validation_evidence_index_is_complete_and_hash_pinned():
    record = build()
    assert record["schema"] == "degali.validation-evidence-index.v1"
    assert record["all_files_present"] is True
    assert record["promotion_allowed"] is False
    assert len(record["records"]) >= 6
    assert all(item["exists"] for item in record["records"])
    assert all(item["sha256"] and len(item["sha256"]) == 64 for item in record["records"])
    records = {item["id"]: item for item in record["records"]}
    assert records["ffi_obstacle_registry"]["exists"] is True
    assert records["ffi_obstacle_case_o1_run_1"]["exists"] is True
    assert records["ffi_obstacle_case_o2_run_2"]["exists"] is True
    assert records["aij_case_h_time_averaged_aligned_3d_grid"]["exists"] is True
    assert records["current_nested_145_sensor_600s_score"]["exists"] is True
    assert records["ti021_corrected_wide_domain_120s_nested_score"]["exists"] is True
    assert records["ti021_invalid_wide_domain_200s_pilot"]["exists"] is True
    assert "invalidated_due_to_crosswind_origin" in records["ti021_invalid_wide_domain_200s_pilot"]["status"]
    assert records["ti021_crosswind_span_prefix_comparison"]["exists"] is True
    assert records["ti021_crosswind_span_prefix_comparison_tool"]["exists"] is True
    assert records["current_nested_600s_inner_checkpoint"]["exists"] is True
    assert records["current_nested_600s_outer_checkpoint"]["exists"] is True
    assert records["thorney_interior_overlap_audit"]["exists"] is True
    assert records["thorney_interior_overlap_audit_tool"]["exists"] is True
    assert "goal_completion_audit_machine" not in records
