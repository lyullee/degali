"""Three-grid AIJ audit using a fixed spin-up plus time-mean observation."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.validation.aij_case_h import (
    compare_aij_case_h_centerplane_flow_3d,
    read_aij_case_h,
)
from tools.audit_aij_resolved_3d_grid import (
    DATA,
    OUTPUT,
    X_ORIGIN_M,
    Y_ORIGIN_M,
    _config,
    _metrics,
    _relative_change,
)


LEGACY_GRID_SPACINGS_M = (0.025, 0.020, 0.0125)
ALIGNED_GRID_SPACINGS_M = (0.025, 1.0 / 60.0, 0.0125)
TOTAL_FLOW_HORIZON_S = 0.5572755417956656
AVERAGING_DURATION_S = 0.18575851393188854
RELATIVE_TOLERANCE = 0.10
STATISTICAL_TOLERANCE = 0.01


def _relative_l2(left: np.ndarray, right: np.ndarray) -> float:
    return float(
        np.linalg.norm(right - left)
        / max(np.linalg.norm(left), np.linalg.norm(right), 1.0e-15)
    )


def build(
    grid_spacings_m: tuple[float, ...],
    *,
    output_stem: str,
    grid_family: str,
) -> int:
    dataset = read_aij_case_h(DATA)
    runs = []
    for spacing in grid_spacings_m:
        base = _config(dataset, spacing)
        flow_dt = base.lattice_inlet_speed * spacing / base.wind_speed_m_s
        flow_steps = int(round(TOTAL_FLOW_HORIZON_S / flow_dt))
        average_steps = int(round(AVERAGING_DURATION_S / flow_dt))
        config = replace(
            base,
            flow_steps=flow_steps,
            minimum_flow_steps=flow_steps,
            flow_average_last_steps=average_steps,
            flow_statistical_convergence_tolerance=STATISTICAL_TOLERANCE,
            subgrid_model="none",
            smagorinsky_constant=0.0,
            flow_boundary_scheme="equilibrium_copy",
        )
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset,
            result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((spacing, config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _spacing, _config_value, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    cases = []
    predictions_u = []
    predictions_w = []
    for (spacing, config, result, comparison), mapping in zip(runs, point_maps):
        points = [mapping[identifier] for identifier in common_ids]
        predictions_u.append(np.asarray([point.predicted_u_over_uh for point in points]))
        predictions_w.append(np.asarray([point.predicted_w_over_uh for point in points]))
        diagnostics = result.diagnostics
        cases.append({
            "cell_size_m": spacing,
            "cells": list(config.cells),
            "flow_steps": config.flow_steps,
            "flow_physical_horizon_s": (
                config.flow_steps * diagnostics["flow_physical_time_step_s"]
            ),
            "flow_average_last_steps": config.flow_average_last_steps,
            "flow_average_duration_s": (
                config.flow_average_last_steps
                * diagnostics["flow_physical_time_step_s"]
            ),
            "matched_metrics": _metrics(points),
            "runtime_s": diagnostics["total_runtime_s"],
            "flow_converged": diagnostics["flow_converged"],
            "flow_statistical_residual": diagnostics["flow_statistical_residual"],
            "final_instantaneous_flow_residual": diagnostics[
                "final_instantaneous_flow_residual"
            ],
            "minimum_lee_streamwise_velocity_m_s": diagnostics[
                "minimum_lee_streamwise_velocity_m_s"
            ],
            "recirculating_lee_cell_count": diagnostics[
                "recirculating_lee_cell_count"
            ],
            "maximum_absolute_divergence_s_1": diagnostics[
                "maximum_absolute_divergence_s_1"
            ],
            "obstacle_grid_alignment_exact": diagnostics[
                "obstacle_grid_alignment_exact"
            ],
            "maximum_obstacle_face_error_m": diagnostics[
                "maximum_obstacle_face_error_m"
            ],
            "maximum_obstacle_face_error_cells": diagnostics[
                "maximum_obstacle_face_error_cells"
            ],
            "maximum_obstacle_size_error_m": diagnostics[
                "maximum_obstacle_size_error_m"
            ],
        })

    pair_labels = ("coarse_to_medium", "medium_to_fine")
    pair_field_changes = {}
    for index, label in enumerate(pair_labels):
        pair_field_changes[label] = {
            "u_relative_l2": _relative_l2(
                predictions_u[index], predictions_u[index + 1]
            ),
            "w_relative_l2": _relative_l2(
                predictions_w[index], predictions_w[index + 1]
            ),
        }

    middle, finest = cases[-2:]
    finest_pair_metric_changes = {
        name: _relative_change(
            float(middle["matched_metrics"][name]),
            float(finest["matched_metrics"][name]),
        )
        for name in (
            "u_mae_normalized", "u_rmse_normalized", "u_correlation",
            "w_mae_normalized", "w_rmse_normalized", "w_correlation",
            "u_recirculation_sign_accuracy",
        )
    }
    finest_pair_lee_change = _relative_change(
        float(middle["minimum_lee_streamwise_velocity_m_s"]),
        float(finest["minimum_lee_streamwise_velocity_m_s"]),
    )
    gate_reasons = []
    for name, change in finest_pair_metric_changes.items():
        if change > RELATIVE_TOLERANCE:
            gate_reasons.append(
                f"finest-pair {name} change {change:.3%} exceeds "
                f"{RELATIVE_TOLERANCE:.1%}"
            )
    if finest_pair_lee_change > RELATIVE_TOLERANCE:
        gate_reasons.append(
            "finest-pair minimum lee velocity change "
            f"{finest_pair_lee_change:.3%} exceeds {RELATIVE_TOLERANCE:.1%}"
        )
    for component in ("u_relative_l2", "w_relative_l2"):
        coarse_change = pair_field_changes["coarse_to_medium"][component]
        fine_change = pair_field_changes["medium_to_fine"][component]
        if fine_change > RELATIVE_TOLERANCE:
            gate_reasons.append(
                f"finest-pair {component} {fine_change:.3%} exceeds "
                f"{RELATIVE_TOLERANCE:.1%}"
            )
        if fine_change >= coarse_change:
            gate_reasons.append(
                f"{component} does not decrease with refinement "
                f"({coarse_change:.3%} -> {fine_change:.3%})"
            )
    for case in cases:
        if not bool(case["flow_converged"]):
            gate_reasons.append(
                f"{case['cell_size_m']:.4f} m time-mean failed the 1% "
                "statistical-stationarity gate"
            )
        if (
            grid_family.startswith("all obstacle faces aligned")
            and not bool(case["obstacle_grid_alignment_exact"])
        ):
            gate_reasons.append(
                f"{case['cell_size_m']:.4f} m declared aligned grid has nonzero "
                "obstacle voxelization error"
            )
    if float(finest["maximum_absolute_divergence_s_1"]) > 1.0e-7:
        gate_reasons.append("finest-grid projected divergence exceeds 1e-7 s-1")
    passed = not gate_reasons
    record = {
        "schema": "degali.aij-case-h-time-averaged-grid-audit.v1",
        "status": "accepted_component_grid" if passed else "withheld_grid_gate",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "comparison_scope": "AIJ center-plane mean U/W on common observation IDs",
        "grid_family": grid_family,
        "total_flow_horizon_s": TOTAL_FLOW_HORIZON_S,
        "averaging_duration_s": AVERAGING_DURATION_S,
        "subgrid_model": "none",
        "flow_boundary_scheme": "equilibrium_copy",
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "relative_tolerance": RELATIVE_TOLERANCE,
        "statistical_tolerance": STATISTICAL_TOLERANCE,
        "cases": cases,
        "model_pair_field_changes": pair_field_changes,
        "finest_pair_metric_changes": finest_pair_metric_changes,
        "finest_pair_minimum_lee_velocity_change": finest_pair_lee_change,
        "passed": passed,
        "gate_reasons": gate_reasons,
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / f"{output_stem}.json"
    figure_path = OUTPUT / f"{output_stem}.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    spacings = np.asarray([item["cell_size_m"] for item in cases])
    figure, axes = plt.subplots(2, 2, figsize=(11.8, 8.0), constrained_layout=True)
    axes[0, 0].plot(
        spacings,
        [item["matched_metrics"]["u_rmse_normalized"] for item in cases],
        "o-", label="$U/U_H$ RMSE",
    )
    axes[0, 0].plot(
        spacings,
        [item["matched_metrics"]["w_rmse_normalized"] for item in cases],
        "s-", label="$W/U_H$ RMSE",
    )
    axes[0, 0].set(title="Matched AIJ error", ylabel="Normalised RMSE")
    axes[0, 0].legend()

    pair_x = np.arange(2)
    axes[0, 1].plot(
        pair_x,
        [pair_field_changes[label]["u_relative_l2"] for label in pair_labels],
        "o-", label="$U$ field change",
    )
    axes[0, 1].plot(
        pair_x,
        [pair_field_changes[label]["w_relative_l2"] for label in pair_labels],
        "s-", label="$W$ field change",
    )
    axes[0, 1].axhline(
        RELATIVE_TOLERANCE, linestyle="--", color="black", linewidth=1,
        label="10% tolerance",
    )
    axes[0, 1].set(
        title="Model-to-model observation-field change",
        ylabel="Relative L2",
        xticks=pair_x,
        xticklabels=("coarse→medium", "medium→fine"),
    )
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(
        spacings,
        [item["flow_statistical_residual"] for item in cases],
        "o-", color="#059669",
    )
    axes[1, 0].axhline(
        STATISTICAL_TOLERANCE, linestyle="--", color="black", linewidth=1,
    )
    axes[1, 0].set(title="Time-mean stationarity", ylabel="Half-window mean change")

    axes[1, 1].plot(
        spacings, [item["runtime_s"] for item in cases], "o-", color="#475569",
    )
    axes[1, 1].set(title="Computational cost", ylabel="Runtime (s)")
    for axis in (axes[0, 0], axes[1, 0], axes[1, 1]):
        axis.invert_xaxis()
        axis.set_xlabel("Cell size (m), finer →")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    figure.suptitle(
        "AIJ Case H — time-averaged three-grid audit — "
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
        "model_pair_field_changes": pair_field_changes,
        "finest_pair_metric_changes": finest_pair_metric_changes,
        "gate_reasons": gate_reasons,
        "runtimes_s": [item["runtime_s"] for item in cases],
    }, indent=2, allow_nan=False))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--grid-family",
        choices=("legacy", "aligned"),
        default="legacy",
        help="use the historical 0.020 m middle grid or an obstacle-face-aligned family",
    )
    args = parser.parse_args()
    if args.grid_family == "aligned":
        return build(
            ALIGNED_GRID_SPACINGS_M,
            output_stem="aij_case_h_time_averaged_aligned_grid_audit",
            grid_family="all obstacle faces aligned with cell faces",
        )
    return build(
        LEGACY_GRID_SPACINGS_M,
        output_stem="aij_case_h_time_averaged_grid_audit",
        grid_family="historical middle grid; obstacle x faces cut cell centres",
    )


if __name__ == "__main__":
    raise SystemExit(main())
