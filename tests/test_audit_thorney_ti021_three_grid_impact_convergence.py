import json

import numpy as np
import pytest

from tools.audit_thorney_ti021_three_grid_impact_convergence import (
    audit_three_grid_runs,
    hazard_geometry,
    solid_geometry_metrics,
    scalar_triplet_convergence,
)


def test_scalar_triplet_reports_monotonic_order_and_gci():
    result = scalar_triplet_convergence(5.0, 3.0, 2.5)
    assert result["classification"] == "monotonic_convergent"
    assert result["observed_order"] == pytest.approx(2.0)
    assert result["fine_gci_absolute"] == pytest.approx(1.25 * 0.5 / 3.0)
    assert result["richardson_extrapolated"] == pytest.approx(2.5 - 0.5 / 3.0)


def test_scalar_triplet_withholds_gci_for_oscillation():
    result = scalar_triplet_convergence(5.0, 3.0, 4.0)
    assert result["classification"] == "oscillatory_or_stationary"
    assert result["observed_order"] is None
    assert result["fine_gci_relative"] is None


def test_scalar_triplet_exact_agreement_has_zero_gci():
    result = scalar_triplet_convergence(4.0, 4.0, 4.0)
    assert result["classification"] == "exact_agreement"
    assert result["fine_gci_absolute"] == 0.0


def test_hazard_geometry_uses_cell_faces_and_excludes_solids():
    field = {
        "x_m": np.asarray((0.5, 1.5, 2.5)),
        "y_m": np.asarray((0.5, 1.5)),
        "z_m": np.asarray((0.5, 1.5)),
        "solid_mask": np.zeros((2, 2, 3), dtype=bool),
        "maximum_concentration_vol_percent": np.zeros((2, 2, 3)),
    }
    field["maximum_concentration_vol_percent"][0, 0, 1:] = 4.0
    field["solid_mask"][0, 0, 2] = True
    result = hazard_geometry(field, cell_size_m=1.0)
    assert result["cell_count"] == 1
    assert result["volume_m3"] == pytest.approx(1.0)
    assert result["x_minimum_face_m"] == pytest.approx(1.0)
    assert result["x_maximum_face_m"] == pytest.approx(2.0)
    assert result["x_extent_m"] == pytest.approx(1.0)


def test_solid_geometry_metrics_reports_voxel_volume_and_bounds():
    field = {
        "solid_mask": np.zeros((2, 2, 2), dtype=bool),
        "x_m": np.asarray((0.5, 1.5)),
        "y_m": np.asarray((0.5, 1.5)),
        "z_m": np.asarray((0.5, 1.5)),
    }
    field["solid_mask"][0, 1, 1] = True
    result = solid_geometry_metrics(field, cell_size_m=1.0)
    assert result["solid_cell_count"] == 1
    assert result["represented_solid_voxel_volume_m3"] == pytest.approx(1.0)
    assert result["x_bounds_m"] == pytest.approx([1.0, 2.0])
    assert result["y_bounds_m"] == pytest.approx([1.0, 2.0])
    assert result["z_bounds_m"] == pytest.approx([0.0, 1.0])
    assert result["horizontal_solid_cell_count"] == 1
    assert result["estimated_arc_length_m"] == pytest.approx(1.0)
    assert result["relative_arc_length_error"] > 0.0


def _write_synthetic_run(
    root, *, cell_size_m, time_step_s,
    circular_arc_representation="finite_thickness_voxel",
    duration_s=20.0,
):
    root.mkdir()
    x = np.arange(0.5 * cell_size_m, 200.0, cell_size_m)
    y = np.arange(0.5 * cell_size_m, 200.0, cell_size_m)
    z = np.arange(0.5 * cell_size_m, 40.0, cell_size_m)
    shape = (z.size, y.size, x.size)
    value = 5.0 + 0.01 * cell_size_m**2
    record = {
        "all_gates_pass": True,
        "configuration": {
            "cell_size_m": cell_size_m,
            "time_step_s": time_step_s,
            "duration_s": duration_s,
            "circular_arc_representation": circular_arc_representation,
        },
    }
    (root / "ti021-inner-transient.json").write_text(
        json.dumps(record), encoding="utf-8",
    )
    np.savez_compressed(
        root / "ti021-inner-final-field.npz",
        x_m=x, y_m=y, z_m=z,
        solid_mask=np.zeros(shape, dtype=bool),
        concentration_vol_percent=np.full(shape, value),
        maximum_concentration_vol_percent=np.full(shape, value),
    )
    np.savez_compressed(
        root / "ti021-inner-receptor-traces.npz",
        labels=np.asarray(("R1", "R2")),
        solver_positions_m=np.asarray(((60.0, 100.0, 1.0), (80.0, 100.0, 1.0))),
        time_s=np.asarray((0.0, 0.5 * duration_s, duration_s)),
        concentration_vol_percent=np.full((2, 3), value),
    )
    return root


def test_three_grid_audit_runs_end_to_end_and_preserves_full_event_holdback(tmp_path):
    coarse = _write_synthetic_run(tmp_path / "coarse", cell_size_m=5.0, time_step_s=0.125)
    medium = _write_synthetic_run(tmp_path / "medium", cell_size_m=2.5, time_step_s=0.0625)
    fine = _write_synthetic_run(tmp_path / "fine", cell_size_m=1.25, time_step_s=0.03125)
    record, paths = audit_three_grid_runs(
        coarse_run=coarse, medium_run=medium, fine_run=fine,
        output_directory=tmp_path / "audit",
    )
    assert not record["impact_window_grid_sensitivity_supported"]
    assert not record["gates"][
        "physical_obstacle_geometry_identified_and_held_constant"
    ]
    assert "one-cell wall-regularisation path" in record[
        "grid_comparison_interpretation"
    ]
    assert not record["full_event_grid_convergence_promoted"]
    assert record["field_convergence"]["final_concentration"][
        "effective_convergence"
    ]["rmse_vol_percent"]["contracts_with_refinement"]
    assert paths["json"].is_file()
    assert paths["plot"].is_file()


def test_matched_time_step_is_accepted_and_receptor_gci_uses_spatial_ratio(tmp_path):
    coarse = _write_synthetic_run(
        tmp_path / "coarse", cell_size_m=5.0, time_step_s=0.02,
    )
    medium = _write_synthetic_run(
        tmp_path / "medium", cell_size_m=2.5, time_step_s=0.02,
    )
    fine = _write_synthetic_run(
        tmp_path / "fine", cell_size_m=1.25, time_step_s=0.02,
    )
    record, _ = audit_three_grid_runs(
        coarse_run=coarse, medium_run=medium, fine_run=fine,
        output_directory=tmp_path / "matched-dt-audit",
    )
    assert record["configuration"]["temporal_comparison_mode"] == (
        "matched_time_step_isolates_spatial_grid_effects"
    )
    assert record["configuration"]["temporal_refinement_ratios"] == [1.0, 1.0]
    assert record["configuration"]["receptor_gci_refinement_ratios"] == [2.0, 2.0]
    assert record["voxelized_obstacle_geometry"][
        "geometry_comparison_classification"
    ] == "one_cell_surface_regularisation_refinement_not_fixed_physical_geometry"


def test_zero_thickness_surface_is_separated_from_physical_thickness_claim(tmp_path):
    runs = [
        _write_synthetic_run(
            tmp_path / label, cell_size_m=dx, time_step_s=0.02,
            circular_arc_representation="zero_thickness_surface",
        )
        for label, dx in (("coarse", 5.0), ("medium", 2.5), ("fine", 1.25))
    ]
    record, _ = audit_three_grid_runs(
        coarse_run=runs[0], medium_run=runs[1], fine_run=runs[2],
        output_directory=tmp_path / "surface-audit",
    )

    assert record["gates"]["numerical_obstacle_representation_held_constant"]
    assert not record["gates"]["physical_fence_thickness_known"]
    assert not record["gates"][
        "physical_obstacle_geometry_identified_and_held_constant"
    ]
    geometry = record["voxelized_obstacle_geometry"]
    assert geometry["geometry_comparison_classification"] == (
        "grid_independent_zero_thickness_embedded_surface"
    )
    reference = geometry["reference_geometry"]
    assert reference["published_construction"] == (
        "tarpaulin sheets hung over a scaffold frame"
    )
    assert reference["numerical_idealization"] == (
        "impermeable zero-thickness surface; scaffold frame omitted"
    )
    assert record["gates"]["field_l1_and_rmse_contract"] is True
    assert isinstance(record["gates"]["hazard_volume_monotonic_gci_le_10_percent"], bool)
    assert isinstance(record["gates"]["receptor_peak_and_dose_monotonic_fraction_ge_70_percent"], bool)
    assert geometry["declared_wall_thickness_by_grid_m"] == [0.0, 0.0, 0.0]
    assert "unknown physical fence thickness" in record[
        "grid_comparison_interpretation"
    ]
    assert "matched-dt 0--20 s spatial audit is complete" in record["next_required"]


def test_short_prefix_is_explicitly_diagnostic_and_cannot_pass_20s_gate(tmp_path):
    runs = [
        _write_synthetic_run(
            tmp_path / label, cell_size_m=dx, time_step_s=0.02,
            circular_arc_representation="zero_thickness_surface", duration_s=3.0,
        )
        for label, dx in (("coarse", 5.0), ("medium", 2.5), ("fine", 1.25))
    ]
    record, _ = audit_three_grid_runs(
        coarse_run=runs[0], medium_run=runs[1], fine_run=runs[2],
        output_directory=tmp_path / "short-prefix-audit",
        allow_short_prefix_diagnostic=True,
    )
    assert record["status"] == "short_prefix_diagnostic_only"
    assert record["configuration"]["window_classification"] == (
        "matched_short_prefix_diagnostic_only"
    )
    assert not record["gates"]["frozen_observation_blind_20s_window"]
    assert not record["impact_window_grid_sensitivity_supported"]
