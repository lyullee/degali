import numpy as np
import pytest
import json
from types import SimpleNamespace

import tools.score_thorney_ti021_nested_145_sensors as nested_score
from tools.score_thorney_ti021_nested_145_sensors import (
    _interpolate_rows,
    merge_inner_outer_traces,
)


def test_nested_merge_covers_145_and_prefers_inner_on_overlap():
    time = np.asarray((0.0, 1.0, 2.0))
    inner_labels = np.asarray([f"TI021-S{i:03d}" for i in range(1, 61)])
    outer_labels = np.asarray([f"TI021-S{i:03d}" for i in range(52, 146)])
    inner = np.tile(np.asarray((0.0, 1.0, 2.0)), (60, 1))
    outer = np.tile(np.asarray((0.0, 10.0, 20.0)), (94, 1))
    merged, coverage = merge_inner_outer_traces(
        inner_time_s=time, inner_labels=inner_labels, inner_traces=inner,
        outer_time_s=time, outer_labels=outer_labels, outer_traces=outer,
    )
    assert coverage["combined_unique_sensor_count"] == 145
    assert coverage["overlap_sensor_count"] == 9
    assert merged["concentration_vol_percent"][51] == pytest.approx(inner[51])
    assert merged["concentration_vol_percent"][60] == pytest.approx(outer[9])
    assert merged["selected_source"].tolist().count("inner") == 60
    assert merged["selected_source"].tolist().count("outer") == 85


def test_outer_rows_are_interpolated_onto_inner_clock():
    result = _interpolate_rows(
        np.asarray((0.0, 1.0, 2.0)), np.asarray(((0.0, 2.0, 4.0),)),
        np.asarray((0.0, 0.5, 1.0, 1.5, 2.0)),
    )
    assert np.allclose(result, ((0.0, 1.0, 2.0, 3.0, 4.0),))


def test_nested_merge_rejects_incomplete_union():
    time = np.asarray((0.0, 1.0))
    with pytest.raises(ValueError, match="does not cover 145"):
        merge_inner_outer_traces(
            inner_time_s=time, inner_labels=np.asarray(("TI021-S001",)),
            inner_traces=np.zeros((1, 2)), outer_time_s=time,
            outer_labels=np.asarray(("TI021-S002",)), outer_traces=np.zeros((1, 2)),
        )


def test_nested_score_runs_after_both_observation_blind_executions(tmp_path, monkeypatch):
    inner = tmp_path / "inner"; outer = tmp_path / "outer"
    inner.mkdir(); outer.mkdir()
    workbook = tmp_path / "ti21.xls"; workbook.write_bytes(b"synthetic")
    common_configuration = {
        "observation_values_used_to_configure_run": [],
        "cells": [40, 40, 8], "domain_m": [200.0, 200.0, 40.0],
        "cell_size_m": 1.0, "duration_s": 600.0,
    }
    def receptor_rows(first: int, last: int, *, inner: bool) -> list[dict[str, object]]:
        rows = []
        for index in range(first, last + 1):
            site_x = float(index)
            site_y = 0.0 if index < 52 or index == 52 else 10.0
            row = {
                "label": f"TI021-S{index:03d}",
                "site_local_xyz_m": [site_x, site_y, 1.0],
            }
            if inner:
                row["solver_xyz_m"] = [60.0 + site_y, 100.0 - site_x, 1.0]
            else:
                row["outer_solver_xyz_m"] = [site_y, 100.0 - site_x, 1.0]
            rows.append(row)
        return rows
    inner_configuration = {
        **common_configuration,
        "receptors": receptor_rows(1, 60, inner=True),
    }
    outer_configuration = {
        **common_configuration,
        "receptors": receptor_rows(52, 145, inner=False),
    }
    (inner / "ti021-inner-transient.json").write_text(json.dumps({
        "status": "accepted_completed_inner_transient",
        "configuration": inner_configuration,
        "coordinates": {
            "handoff_site_downwind_y_m": 0.0,
            "source_solver_xyz_m": [60.0, 100.0, 0.0],
        },
    }), encoding="utf-8")
    (outer / "ti021-outer-transient.json").write_text(json.dumps({
        "status": "accepted_completed_outer_transient",
        "configuration": outer_configuration,
    }), encoding="utf-8")
    time = np.asarray((0.0, 1.0, 2.0))
    inner_labels = np.asarray([f"TI021-S{i:03d}" for i in range(1, 61)])
    outer_labels = np.asarray([f"TI021-S{i:03d}" for i in range(52, 146)])
    scale_inner = np.arange(1.0, 61.0)[:, None] / 100.0
    scale_inner[1] = scale_inner[0]
    scale_outer = np.arange(52.0, 146.0)[:, None] / 100.0
    scale_outer[0] *= 2.0  # boundary monitor must not contaminate interior acceptance
    shape = np.asarray((0.0, 1.0, 0.0))[None, :]
    np.savez_compressed(
        inner / "ti021-inner-receptor-traces.npz", time_s=time,
        labels=inner_labels, concentration_vol_percent=scale_inner * shape,
    )
    np.savez_compressed(
        outer / "ti021-outer-receptor-traces.npz", time_s=time,
        labels=outer_labels, concentration_vol_percent=scale_outer * shape,
    )
    sensors = tuple(SimpleNamespace(
        local_x_m=0.0 if index < 2 else float(index),
        local_y_m=0.0,
        height_m=0.4 if index == 0 else 2.4 if index == 1 else 4.4,
        dose_percent_s=float(index + 1) / 100.0,
        arrival_time_s=1.0, departure_time_s=1.0,
    ) for index in range(145))
    monkeypatch.setattr(
        nested_score, "read_thorney_ti021_workbook",
        lambda path: SimpleNamespace(dose_sensors=sensors),
    )
    score, arrays = nested_score.score_nested_run(
        workbook=workbook, inner_run_directory=inner,
        outer_run_directory=outer,
    )
    assert score["coverage"]["combined_unique_sensor_count"] == 145
    assert score["dose_sensor_score"]["sensor_count"] == 145
    assert arrays["concentration_vol_percent"].shape == (145, 3)
    resolution = score["vertical_sensor_resolution_diagnostic"]
    assert resolution["bottom_cell_centre_m"] == pytest.approx(2.5)
    assert resolution["same_horizontal_position_pairs_on_bottom_cell"] == 1
    assert resolution["pairs_with_identical_predicted_traces"] == 1
    assert resolution["positive_observed_dose_ratio_median"] == pytest.approx(2.0)
    geometry = score["overlap_audit"]["geometry"]
    assert geometry["boundary_coincident_sensor_labels"] == ["TI021-S052"]
    assert geometry["interior_overlap_sensor_count"] == 8
    assert geometry["interior_overlap_horizontal_location_count"] == 8
    assert score["overlap_audit"]["gates"]["at_least_three_interior_overlap_locations"] is True
    assert score["overlap_audit"]["all_gates_pass"] is True
    assert score["overlap_audit"]["maximum_boundary_monitor_peak_absolute_relative_change"] > 0.10


def test_overlap_geometry_rejects_boundary_only_sensors():
    from tools.score_thorney_ti021_nested_145_sensors import _overlap_geometry_diagnostic

    labels = np.asarray(["TI021-S040", "TI021-S047", "TI021-S048"])
    rows = [
        {
            "label": str(label), "site_local_xyz_m": [0.0, 100.0, float(index)],
            "solver_xyz_m": [160.0, 100.0, float(index)],
        }
        for index, label in enumerate(labels)
    ]
    inner = {
        "coordinates": {"handoff_site_downwind_y_m": 100.0},
        "configuration": {
            "cell_size_m": 5.0, "domain_m": [300.0, 200.0, 40.0],
            "receptors": rows,
        },
    }
    outer = {"configuration": {
        "cell_size_m": 5.0, "domain_m": [440.0, 640.0, 40.0],
        "receptors": [
            {
                **row,
                "outer_solver_xyz_m": [0.0, 320.0, row["site_local_xyz_m"][2]],
            }
            for row in rows
        ],
    }}
    diagnostic = _overlap_geometry_diagnostic(inner, outer, labels)
    assert diagnostic["boundary_coincident_sensor_labels"] == labels.tolist()
    assert diagnostic["interior_overlap_sensor_count"] == 0
    assert diagnostic["interior_overlap_horizontal_location_count"] == 0
    assert diagnostic["at_least_three_interior_horizontal_locations"] is False

    shifted_rows = [
        {
            **row,
            "solver_xyz_m": [160.0, 100.0, row["site_local_xyz_m"][2]],
        }
        for row in rows
    ]
    wide_inner = {
        "coordinates": {
            "handoff_site_downwind_y_m": 100.0,
            "source_solver_xyz_m": [60.0, 100.0, 0.0],
        },
        "configuration": {
            "domain_m": [300.0, 640.0, 40.0], "cell_size_m": 5.0,
            "receptors": shifted_rows,
        },
    }
    wide_outer = {
        "configuration": {
            "domain_m": [440.0, 640.0, 40.0], "cell_size_m": 5.0,
            "receptors": [
                {**row, "outer_solver_xyz_m": [0.0, 320.0, row["site_local_xyz_m"][2]]}
                for row in rows
            ],
        },
    }
    bad_mapping = _overlap_geometry_diagnostic(wide_inner, wide_outer, labels)
    assert bad_mapping["inner_solver_coordinate_mismatch_sensor_labels"] == labels.tolist()
    assert bad_mapping["inner_crosswind_origin_mismatch"] is True
    assert bad_mapping["all_shared_sensor_coordinates_match"] is False
