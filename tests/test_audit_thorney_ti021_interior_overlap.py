from __future__ import annotations

import json

import numpy as np

from tools.audit_thorney_ti021_interior_overlap import build_audit


def test_short_interior_overlap_audit_requires_real_shared_geometry(tmp_path):
    inner = tmp_path / "inner"
    outer = tmp_path / "outer"
    inner.mkdir(); outer.mkdir()
    time = np.asarray((0.0, 1.0, 2.0))
    inner_labels = ["TI021-S001", "TI021-S002", "TI021-S003", "TI021-S004", "TI021-S047"]
    outer_labels = ["TI021-S002", "TI021-S003", "TI021-S004", "TI021-S005", "TI021-S047"]

    def receptors(labels, *, inner):
        rows = []
        for label in labels:
            site_x = float(int(label[-3:]))
            site_y = 100.0 if label.endswith("047") else 110.0
            row = {
                "label": label,
                "site_local_xyz_m": [site_x, site_y, 1.0],
            }
            if inner:
                row["solver_xyz_m"] = [60.0 + site_y, 100.0 - site_x, 1.0]
            else:
                row["outer_solver_xyz_m"] = [site_y - 100.0, 320.0 - site_x, 1.0]
            rows.append(row)
        return rows

    inner_record = {
        "status": "accepted_inner_transient_qualification",
        "all_gates_pass": True,
        "runtime": {"wall_clock_s": 1.0},
        "coordinates": {
            "handoff_site_downwind_y_m": 100.0,
            "source_solver_xyz_m": [60.0, 100.0, 0.0],
        },
        "configuration": {
            "duration_s": 2.0, "domain_m": [200.0, 200.0, 40.0],
            "cell_size_m": 5.0, "receptors": receptors(inner_labels, inner=True),
        },
    }
    outer_record = {
        "status": "accepted_completed_outer_transient",
        "all_gates_pass": True,
        "runtime": {"wall_clock_s": 2.0},
        "configuration": {
            "duration_s": 2.0, "domain_m": [440.0, 640.0, 40.0],
            "cell_size_m": 5.0, "receptors": receptors(outer_labels, inner=False),
        },
    }
    (inner / "ti021-inner-transient.json").write_text(json.dumps(inner_record), encoding="utf-8")
    (outer / "ti021-outer-transient.json").write_text(json.dumps(outer_record), encoding="utf-8")
    common_trace = np.tile(np.asarray((0.0, 0.2, 0.0)), (5, 1))
    np.savez_compressed(
        inner / "ti021-inner-receptor-traces.npz", time_s=time,
        labels=np.asarray(inner_labels), concentration_vol_percent=common_trace,
    )
    np.savez_compressed(
        outer / "ti021-outer-receptor-traces.npz", time_s=time,
        labels=np.asarray(outer_labels),
        concentration_vol_percent=common_trace * np.asarray((1.0, 1.0, 1.0, 1.0, 2.0))[:, None],
    )

    record, _ = build_audit(inner_run_directory=inner, outer_run_directory=outer)
    assert record["status"] == "accepted_short_interior_overlap_diagnostic_only"
    assert record["coverage"]["shared_sensor_count"] == 4
    assert record["geometry"]["interior_overlap_horizontal_location_count"] == 3
    assert record["comparison"]["maximum_peak_relative_change_using_0p1vol_floor"] == 0.0
    boundary = record["comparison"]["boundary_interface_monitoring"]
    assert boundary["sensor_count"] == 1
    assert boundary["maximum_dose_relative_change"] > 0.10
    assert record["gates"]["full_observation_horizon_scored"] is False
    assert record["gates"]["experimental_model_promotion_allowed"] is False
