"""Three-grid AIJ mean-flow audit for the resolved D3Q19 obstacle prototype."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.validation.aij_case_h import (
    compare_aij_case_h_centerplane_flow_3d,
    read_aij_case_h,
)


DATA = Path("reference/aij_case_h")
OUTPUT = Path("outputs/obstacle-physics-rebuild-2026-10-10")
GRID_SPACINGS_M = (0.025, 0.020, 0.0125)
FLOW_PHYSICAL_HORIZON_S = 0.2786377708978328
X_ORIGIN_M = -0.25
Y_ORIGIN_M = -0.15


def _config(dataset, spacing: float) -> ResolvedObstacleFlow3DConfig:
    cells = tuple(int(round(length / spacing)) for length in (0.70, 0.30, 0.40))
    profile_z = (0.0,) + tuple(item.z_m for item in dataset.approach_points)
    profile_u = (0.0,) + tuple(item.u_m_s for item in dataset.approach_points)
    flow_dt = 0.030 * spacing / dataset.reference_speed_m_s
    flow_steps = int(round(FLOW_PHYSICAL_HORIZON_S / flow_dt))
    return ResolvedObstacleFlow3DConfig(
        domain_m=(0.70, 0.30, 0.40),
        cells=cells,
        wind_speed_m_s=dataset.reference_speed_m_s,
        kinematic_viscosity_m2_s=0.020,
        flow_steps=flow_steps,
        minimum_flow_steps=flow_steps,
        lattice_inlet_speed=0.030,
        duration_s=1.0e-5,
        scalar_time_step_s=1.0e-5,
        scalar_diffusivity_m2_s=1.0e-5,
        source_schedule=SourceRateSchedule((0.0, 1.0e-5), (0.0, 0.0)),
        source_position_m=(0.30, 0.15, spacing / 2.0),
        source_sigma_m=spacing,
        obstacles=(AxisAlignedCuboid(
            0.15, 0.25, 0.10, 0.20, 0.0, 0.20, "AIJ Case H finite block"
        ),),
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
        smagorinsky_constant=0.12,
    )


def _metrics(points) -> dict[str, float | int]:
    measured_u = np.asarray([point.measured_u_over_uh for point in points])
    predicted_u = np.asarray([point.predicted_u_over_uh for point in points])
    measured_w = np.asarray([point.measured_w_over_uh for point in points])
    predicted_w = np.asarray([point.predicted_w_over_uh for point in points])
    sign_mask = np.abs(measured_u) > 0.02
    return {
        "point_count": len(points),
        "u_mae_normalized": float(np.mean(np.abs(predicted_u - measured_u))),
        "u_rmse_normalized": float(np.sqrt(np.mean((predicted_u - measured_u) ** 2))),
        "u_correlation": float(np.corrcoef(measured_u, predicted_u)[0, 1]),
        "w_mae_normalized": float(np.mean(np.abs(predicted_w - measured_w))),
        "w_rmse_normalized": float(np.sqrt(np.mean((predicted_w - measured_w) ** 2))),
        "w_correlation": float(np.corrcoef(measured_w, predicted_w)[0, 1]),
        "u_recirculation_sign_accuracy": float(np.mean(
            np.sign(measured_u[sign_mask]) == np.sign(predicted_u[sign_mask])
        )),
    }


def _relative_change(left: float, right: float, floor: float = 1.0e-12) -> float:
    return abs(right - left) / max(abs(left), abs(right), floor)


def main() -> int:
    dataset = read_aij_case_h(DATA)
    runs = []
    for spacing in GRID_SPACINGS_M:
        base_config = _config(dataset, spacing)
        config = replace(
            base_config,
            flow_average_last_steps=max(
                2, int(round(0.75 * base_config.flow_steps))
            ),
            flow_statistical_convergence_tolerance=0.01,
        )
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset, result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((spacing, config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _spacing, _config_value, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    records = []
    for (spacing, config, result, comparison), mapping in zip(runs, point_maps):
        metrics = _metrics([mapping[identifier] for identifier in common_ids])
        records.append({
            "cell_size_m": spacing,
            "cells": list(config.cells),
            "flow_steps": config.flow_steps,
            "flow_physical_horizon_s": (
                config.flow_steps
                * config.lattice_inlet_speed * spacing
                / config.wind_speed_m_s
            ),
            "evaluated_point_count_before_matching": len(comparison.points),
            "skipped_point_count": len(comparison.skipped_point_ids),
            "matched_metrics": metrics,
            "runtime_s": result.diagnostics["total_runtime_s"],
            "final_flow_residual": result.diagnostics["final_flow_residual"],
            "flow_converged": result.diagnostics["flow_converged"],
            "flow_statistical_residual": result.diagnostics[
                "flow_statistical_residual"
            ],
            "maximum_relative_density_variation": result.diagnostics[
                "maximum_relative_lattice_density_variation"
            ],
            "maximum_absolute_raw_face_divergence_s_1": result.diagnostics[
                "maximum_absolute_raw_face_divergence_s_1"
            ],
            "maximum_absolute_divergence_s_1": result.diagnostics[
                "maximum_absolute_divergence_s_1"
            ],
            "face_projection_iterations": result.diagnostics[
                "face_projection_iterations"
            ],
            "recirculating_lee_cell_count": result.diagnostics[
                "recirculating_lee_cell_count"
            ],
            "minimum_lee_streamwise_velocity_m_s": result.diagnostics[
                "minimum_lee_streamwise_velocity_m_s"
            ],
            "maximum_eddy_viscosity_m2_s": result.diagnostics[
                "maximum_eddy_viscosity_m2_s"
            ],
        })

    middle, finest = records[-2], records[-1]
    metric_changes = {}
    for name in (
        "u_mae_normalized", "u_rmse_normalized", "u_correlation",
        "w_mae_normalized", "w_rmse_normalized", "w_correlation",
        "u_recirculation_sign_accuracy",
    ):
        metric_changes[name] = _relative_change(
            float(middle["matched_metrics"][name]),
            float(finest["matched_metrics"][name]),
        )
    lee_velocity_change = _relative_change(
        float(middle["minimum_lee_streamwise_velocity_m_s"]),
        float(finest["minimum_lee_streamwise_velocity_m_s"]),
    )
    tolerance = 0.10
    gate_reasons = []
    for name, change in metric_changes.items():
        if change > tolerance:
            gate_reasons.append(
                f"finest-pair {name} change {change:.3%} exceeds {tolerance:.1%}"
            )
    if lee_velocity_change > tolerance:
        gate_reasons.append(
            "finest-pair minimum lee velocity change "
            f"{lee_velocity_change:.3%} exceeds {tolerance:.1%}"
        )
    if not bool(finest["flow_converged"]):
        statistical_residual = finest["flow_statistical_residual"]
        if statistical_residual is None or float(statistical_residual) > 0.01:
            gate_reasons.append(
                "finest-grid flow did not meet instantaneous or statistical "
                "residual tolerance"
            )
    if float(finest["maximum_absolute_divergence_s_1"]) > 1.0e-7:
        gate_reasons.append("finest-grid projected face divergence exceeds 1e-7 s-1")
    passed = not gate_reasons
    record = {
        "schema": "degali.aij-case-h-resolved-3d-grid-audit.v1",
        "status": "accepted" if passed else "withheld",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "comparison_scope": "AIJ center-plane mean U/W on common observation IDs",
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "relative_tolerance": tolerance,
        "cases": records,
        "finest_pair_relative_changes": metric_changes,
        "finest_pair_minimum_lee_velocity_change": lee_velocity_change,
        "passed": passed,
        "gate_reasons": gate_reasons,
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_resolved_3d_grid_audit.json"
    figure_path = OUTPUT / "aij_case_h_resolved_3d_grid_audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    spacings = np.asarray([item["cell_size_m"] for item in records])
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    axes[0, 0].plot(
        spacings,
        [item["matched_metrics"]["u_rmse_normalized"] for item in records],
        "o-", label="$U/U_H$ RMSE",
    )
    axes[0, 0].plot(
        spacings,
        [item["matched_metrics"]["w_rmse_normalized"] for item in records],
        "s-", label="$W/U_H$ RMSE",
    )
    axes[0, 0].set(ylabel="Normalised RMSE", title="Matched AIJ error")
    axes[0, 0].legend()
    axes[0, 1].plot(
        spacings,
        [item["matched_metrics"]["u_correlation"] for item in records],
        "o-", label="$U$ correlation",
    )
    axes[0, 1].plot(
        spacings,
        [item["matched_metrics"]["w_correlation"] for item in records],
        "s-", label="$W$ correlation",
    )
    axes[0, 1].set(ylabel="Pearson r", title="Matched AIJ correlation")
    axes[0, 1].legend()
    axes[1, 0].plot(
        spacings,
        [item["minimum_lee_streamwise_velocity_m_s"] for item in records],
        "o-", color="#7c3aed",
    )
    axes[1, 0].set(ylabel="Minimum lee U (m/s)", title="Wake-grid sensitivity")
    axes[1, 1].plot(
        spacings, [item["runtime_s"] for item in records], "o-", label="runtime"
    )
    residual_axis = axes[1, 1].twinx()
    residual_axis.plot(
        spacings,
        [item["final_flow_residual"] for item in records],
        "s--", color="#dc2626", label="flow residual",
    )
    axes[1, 1].set(ylabel="Runtime (s)", title="Cost and residual")
    residual_axis.set_ylabel("Final flow residual")
    for axis in axes.flat:
        axis.invert_xaxis()
        axis.set_xlabel("Cell size (m), finer →")
        axis.grid(alpha=0.25)
    figure.suptitle(
        "AIJ Case H resolved 3-D grid audit — "
        + ("PASS" if passed else "WITHHELD"),
        fontsize=14,
    )
    figure.savefig(figure_path, dpi=220)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "common_point_count": len(common_ids),
        "gate_reasons": gate_reasons,
        "finest_runtime_s": finest["runtime_s"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
