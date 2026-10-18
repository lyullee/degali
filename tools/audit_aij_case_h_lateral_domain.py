"""Audit AIJ Case H lateral-domain truncation without concentration fitting.

The earlier 3-D AIJ component used a 0.30 m lateral domain while the public
measurement grid spans approximately 0.40 m.  This audit keeps the obstacle,
inlet profile, lattice spacing and physical flow horizon fixed, extends only
the lateral domain to 0.50 m, and compares the center-plane mean U/W metrics.
It is a boundary-sensitivity diagnostic; it is not H2 validation.
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.validation.aij_case_h import (
    compare_aij_case_h_centerplane_flow_3d,
    read_aij_case_h,
)
from tools.audit_aij_resolved_3d_grid import _config, _metrics


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reference/aij_case_h"
OUTPUT = ROOT / "outputs/obstacle-physics-rebuild-2026-10-10"
CELL_SIZE_M = 0.025
FLOW_HORIZON_S = 0.75
AVERAGE_LAST_HORIZON_S = 0.25
BASE_LATERAL_DOMAIN_M = 0.30
WIDE_LATERAL_DOMAIN_M = 0.50
X_ORIGIN_M = -0.25
BASE_Y_ORIGIN_M = -0.15
WIDE_Y_ORIGIN_M = -0.25


def _case_config(dataset, lateral_domain_m: float, y_origin_m: float):
    base = _config(dataset, CELL_SIZE_M)
    cells = (
        int(round(base.domain_m[0] / CELL_SIZE_M)),
        int(round(lateral_domain_m / CELL_SIZE_M)),
        int(round(base.domain_m[2] / CELL_SIZE_M)),
    )
    shift = BASE_Y_ORIGIN_M - y_origin_m
    flow_dt = base.lattice_inlet_speed * CELL_SIZE_M / base.wind_speed_m_s
    flow_steps = int(round(FLOW_HORIZON_S / flow_dt))
    average_steps = max(2, int(round(AVERAGE_LAST_HORIZON_S / flow_dt)))
    return replace(
        base,
        domain_m=(base.domain_m[0], lateral_domain_m, base.domain_m[2]),
        cells=cells,
        flow_steps=flow_steps,
        minimum_flow_steps=flow_steps,
        flow_average_last_steps=average_steps,
        flow_statistical_convergence_tolerance=0.01,
        subgrid_model="none",
        smagorinsky_constant=0.0,
        flow_boundary_scheme="equilibrium_copy",
        source_position_m=(base.source_position_m[0], base.source_position_m[1] + shift, base.source_position_m[2]),
        obstacles=tuple(
            replace(obstacle, y_min_m=obstacle.y_min_m + shift, y_max_m=obstacle.y_max_m + shift)
            for obstacle in base.obstacles
        ),
    )


def _run(dataset, lateral_domain_m: float, y_origin_m: float) -> dict[str, object]:
    config = _case_config(dataset, lateral_domain_m, y_origin_m)
    result = solve_resolved_obstacle_flow_3d(config)
    comparison = compare_aij_case_h_centerplane_flow_3d(
        dataset,
        result,
        model_x_origin_in_experiment_m=X_ORIGIN_M,
        model_y_origin_in_experiment_m=y_origin_m,
    )
    mapping = {item.point_id: item for item in comparison.points}
    metrics = _metrics([mapping[key] for key in sorted(mapping)])
    return {
        "lateral_domain_m": lateral_domain_m,
        "lateral_cells": config.cells[1],
        "model_y_origin_in_experiment_m": y_origin_m,
        "common_point_count": len(mapping),
        "skipped_point_count": len(comparison.skipped_point_ids),
        "metrics": metrics,
        "runtime_s": result.diagnostics["total_runtime_s"],
        "flow_physical_horizon_s": (
            result.diagnostics["flow_steps_completed"]
            * result.diagnostics["flow_physical_time_step_s"]
        ),
        "flow_converged": result.diagnostics["flow_converged"],
        "flow_statistical_residual": result.diagnostics["flow_statistical_residual"],
        "final_flow_residual": result.diagnostics["final_flow_residual"],
        "recirculating_lee_cell_count": result.diagnostics["recirculating_lee_cell_count"],
        "minimum_lee_streamwise_velocity_m_s": result.diagnostics["minimum_lee_streamwise_velocity_m_s"],
    }


def main() -> int:
    dataset = read_aij_case_h(DATA)
    base = _run(dataset, BASE_LATERAL_DOMAIN_M, BASE_Y_ORIGIN_M)
    wide = _run(dataset, WIDE_LATERAL_DOMAIN_M, WIDE_Y_ORIGIN_M)
    metric_changes = {
        key: float(wide["metrics"][key] - base["metrics"][key])
        for key in ("u_mae_normalized", "u_rmse_normalized", "w_mae_normalized", "w_rmse_normalized")
    }
    output = {
        "schema": "degali.aij-case-h-lateral-domain-sensitivity.v1",
        "status": "matched_flow_boundary_sensitivity_not_hydrogen_validation",
        "calibration_used": False,
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "fixed_configuration": {
            "cell_size_m": CELL_SIZE_M,
            "physical_flow_horizon_s": base["flow_physical_horizon_s"],
            "obstacle_and_inlet_profile_held_fixed": True,
            "center_plane_comparison": True,
        },
        "base_domain": base,
        "wide_domain": wide,
        "wide_minus_base_metric_change": metric_changes,
        "interpretation": (
            "The lateral domain was expanded without fitting concentration or changing the obstacle/source physics. "
            "Any change diagnoses side-boundary sensitivity; it does not establish facility or hydrogen validation."
        ),
        "promotion_allowed": False,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_lateral_domain_sensitivity.json"
    figure_path = OUTPUT / "aij_case_h_lateral_domain_sensitivity.png"
    json_path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    labels = ["0.30 m\nbase", "0.50 m\nwide"]
    cases = [base, wide]
    x = np.arange(2)
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.8), constrained_layout=True)
    for axis, metric, title in zip(
        axes,
        ("u_rmse_normalized", "w_rmse_normalized"),
        ("Center-plane U RMSE", "Center-plane W RMSE"),
    ):
        axis.bar(x, [case["metrics"][metric] for case in cases], color=("#64748b", "#0f766e"))
        axis.set_xticks(x, labels)
        axis.set_title(title)
        axis.set_ylabel("normalized error")
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle("AIJ Case H — lateral-domain sensitivity (no concentration fit)")
    figure.savefig(figure_path, dpi=200)
    plt.close(figure)
    print(json.dumps({"json": str(json_path.resolve()), "figure": str(figure_path.resolve()), **output}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
