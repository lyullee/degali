from __future__ import annotations

import json

import numpy as np

from tools.compare_thorney_ti021_temporal_prefix import (
    _common_fluid_component_metrics,
    compare_temporal_prefixes,
)


def _write_prefix(path, *, time_s, field, receptor, flow_offset=0.0):
    path.mkdir()
    (path / "ti021-inner-transient.json").write_text(
        json.dumps({"all_gates_pass": True}), encoding="utf-8",
    )
    np.savez_compressed(
        path / "ti021-inner-final-field.npz",
        x_m=np.array([0.0, 1.0]),
        y_m=np.array([0.0, 1.0]),
        z_m=np.array([0.0, 1.0]),
        maximum_concentration_vol_percent=np.asarray(field, dtype=float),
        solid_mask=np.zeros((2, 2, 2), dtype=bool),
        temperature_k=np.full((2, 2, 2), 290.0 + flow_offset),
        velocity_u_m_s=np.full((2, 2, 2), flow_offset),
        pressure_pa_gauge=np.full((2, 2, 2), 2.0 * flow_offset),
    )
    np.savez_compressed(
        path / "ti021-inner-receptor-traces.npz",
        time_s=np.asarray(time_s, dtype=float),
        labels=np.array(["S001", "S002"]),
        concentration_vol_percent=np.asarray(receptor, dtype=float),
    )


def test_temporal_prefix_comparison_interpolates_receptors(tmp_path):
    reference = tmp_path / "reference"
    comparison = tmp_path / "comparison"
    _write_prefix(
        reference,
        time_s=[0.0, 1.0, 2.0],
        field=np.array([[[0.0, 4.0], [2.0, 8.0]], [[1.0, 3.0], [5.0, 7.0]]]),
        receptor=np.array([[0.0, 2.0, 4.0], [1.0, 3.0, 5.0]]),
        flow_offset=0.0,
    )
    _write_prefix(
        comparison,
        time_s=[0.0, 2.0],
        field=np.array([[[0.0, 4.0], [2.0, 8.0]], [[1.0, 3.0], [5.0, 7.0]]]),
        receptor=np.array([[0.0, 4.0], [1.0, 5.0]]),
        flow_offset=1.0,
    )
    output = tmp_path / "comparison.json"
    record = compare_temporal_prefixes(reference, comparison, output)
    assert record["field"]["rmse"] == 0.0
    assert record["receptors"]["rmse"] == 0.0
    assert record["receptors"]["common_frame_count"] == 3
    components = record["physical_field_components"]["compared_components"]
    assert components["velocity_u_m_s"]["bias"] == -1.0
    assert components["temperature_k"]["rmse"] == 1.0
    assert components["pressure_pa_gauge"]["maximum_absolute_difference"] == 2.0
    assert len(record["provenance"]["reference_field_sha256"]) == 64
    assert len(record["provenance"]["comparison_trace_sha256"]) == 64
    assert json.loads(output.read_text(encoding="utf-8"))["status"] == (
        "computed_window_temporal_comparison"
    )


def test_component_comparison_reports_front_roof_side_and_wake_regions():
    shape = (2, 2, 3)
    reference = {
        "x_m": np.array([105.0, 115.0, 60.0]),
        "y_m": np.array([100.0, 55.0]),
        "z_m": np.array([2.0, 7.0]),
        "solid_mask": np.zeros(shape, dtype=bool),
        "concentration_vol_percent": np.zeros(shape),
        "velocity_u_m_s": np.zeros(shape),
    }
    comparison = {
        **reference,
        "concentration_vol_percent": np.ones(shape),
        "velocity_u_m_s": np.full(shape, 0.25),
    }
    result = _common_fluid_component_metrics(reference, comparison)
    regions = result["obstacle_region_components"]
    assert set(regions) == {
        "windward_inner_face_shell", "roof_shear_band",
        "side_tip_bypass", "exterior_lee_wake_shell",
    }
    for region in regions.values():
        assert region["common_fluid_cell_count"] > 0
        assert region["components"]["concentration_vol_percent"]["rmse"] == 1.0
        assert region["components"]["velocity_u_m_s"]["rmse"] == 0.25
