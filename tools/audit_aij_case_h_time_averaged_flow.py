"""Evaluate a spin-up plus finite-window mean against AIJ Case H mean flow."""

from __future__ import annotations

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
)


CELL_SIZE_M = 0.020
CASES = (
    ("instantaneous", 1_500, 0),
    ("time_mean", 3_000, 1_000),
)


def _percent_change(candidate: float, reference: float) -> float:
    return 100.0 * (candidate - reference) / max(abs(reference), 1.0e-12)


def main() -> int:
    dataset = read_aij_case_h(DATA)
    base = replace(
        _config(dataset, CELL_SIZE_M),
        subgrid_model="none",
        smagorinsky_constant=0.0,
        flow_boundary_scheme="equilibrium_copy",
        flow_statistical_convergence_tolerance=0.01,
    )
    runs = []
    for label, steps, average_steps in CASES:
        config = replace(
            base,
            flow_steps=steps,
            minimum_flow_steps=steps,
            flow_average_last_steps=average_steps,
        )
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset,
            result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((label, config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _label, _config_value, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    cases = []
    for (label, config, result, comparison), mapping in zip(runs, point_maps):
        diagnostics = result.diagnostics
        cases.append({
            "label": label,
            "flow_steps": config.flow_steps,
            "flow_physical_horizon_s": (
                config.flow_steps * diagnostics["flow_physical_time_step_s"]
            ),
            "flow_average_last_steps": config.flow_average_last_steps,
            "flow_average_duration_s": (
                config.flow_average_last_steps
                * diagnostics["flow_physical_time_step_s"]
            ),
            "flow_field_statistic": diagnostics["flow_field_statistic"],
            "metrics": _metrics([mapping[identifier] for identifier in common_ids]),
            "runtime_s": diagnostics["total_runtime_s"],
            "flow_converged": diagnostics["flow_converged"],
            "final_flow_residual": diagnostics["final_flow_residual"],
            "final_instantaneous_flow_residual": diagnostics[
                "final_instantaneous_flow_residual"
            ],
            "flow_statistical_residual": diagnostics[
                "flow_statistical_residual"
            ],
            "minimum_lee_streamwise_velocity_m_s": diagnostics[
                "minimum_lee_streamwise_velocity_m_s"
            ],
            "recirculating_lee_cell_count": diagnostics[
                "recirculating_lee_cell_count"
            ],
            "residual_history": result.flow_residual_history.tolist(),
        })

    instantaneous, averaged = cases
    metric_changes = {
        name: _percent_change(
            float(averaged["metrics"][name]),
            float(instantaneous["metrics"][name]),
        )
        for name in (
            "u_mae_normalized", "u_rmse_normalized", "u_correlation",
            "w_mae_normalized", "w_rmse_normalized", "w_correlation",
            "u_recirculation_sign_accuracy",
        )
    }
    gate_reasons = []
    if not bool(averaged["flow_converged"]):
        gate_reasons.append(
            "time-mean first-half/second-half residual exceeds 1%"
        )
    # Accuracy is reported, not used to select the averaging window.  The
    # window is accepted only as statistically stationary component evidence.
    stationarity_passed = not gate_reasons
    record = {
        "schema": "degali.aij-case-h-time-averaged-flow-audit.v1",
        "status": (
            "accepted_component_stationarity"
            if stationarity_passed else "withheld_statistical_stationarity"
        ),
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "comparison_scope": "AIJ center-plane mean U/W on common observation IDs",
        "cell_size_m": CELL_SIZE_M,
        "cells": list(base.cells),
        "subgrid_model": "none",
        "flow_boundary_scheme": "equilibrium_copy",
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "cases": cases,
        "time_mean_vs_instantaneous_change_percent": metric_changes,
        "statistical_tolerance": base.flow_statistical_convergence_tolerance,
        "stationarity_passed": stationarity_passed,
        "gate_reasons": gate_reasons,
        "grid_convergence_passed": False,
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_time_averaged_flow_audit.json"
    figure_path = OUTPUT / "aij_case_h_time_averaged_flow_audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    labels = ("Final instantaneous", "Spin-up + time mean")
    x = np.arange(2)
    width = 0.36
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), constrained_layout=True)
    axes[0, 0].bar(
        x - width / 2,
        [item["metrics"]["u_rmse_normalized"] for item in cases],
        width, label="$U/U_H$", color="#2563eb",
    )
    axes[0, 0].bar(
        x + width / 2,
        [item["metrics"]["w_rmse_normalized"] for item in cases],
        width, label="$W/U_H$", color="#f59e0b",
    )
    axes[0, 0].set(title="Matched AIJ RMSE", ylabel="Normalised RMSE")
    axes[0, 0].legend()
    axes[0, 1].bar(
        x - width / 2,
        [item["metrics"]["u_correlation"] for item in cases],
        width, label="$U$ correlation", color="#0f766e",
    )
    axes[0, 1].bar(
        x + width / 2,
        [item["metrics"]["w_correlation"] for item in cases],
        width, label="$W$ correlation", color="#7c3aed",
    )
    axes[0, 1].set(title="Matched AIJ correlation", ylabel="Pearson r")
    axes[0, 1].legend()
    for item, color, label in zip(cases, ("#64748b", "#dc2626"), labels):
        history = np.asarray(item["residual_history"])
        time_s = history[:, 0] * base.lattice_inlet_speed * CELL_SIZE_M / base.wind_speed_m_s
        axes[1, 0].semilogy(time_s, history[:, 1], "o-", color=color, label=label)
    axes[1, 0].set(
        title="Instantaneous residual history",
        xlabel="Physical horizon (s)",
        ylabel="100-step velocity-change residual",
    )
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].bar(
        (0,), (averaged["flow_statistical_residual"],), color="#059669",
        label="first-half vs second-half mean",
    )
    axes[1, 1].axhline(
        base.flow_statistical_convergence_tolerance,
        color="black", linestyle="--", label="1% tolerance",
    )
    axes[1, 1].set(
        title="Time-mean stationarity gate", ylabel="Normalised mean-field change",
        xticks=(0,), xticklabels=("1000-step window",),
    )
    axes[1, 1].legend(fontsize=8)
    for axis in (axes[0, 0], axes[0, 1]):
        axis.set_xticks(x, labels, rotation=8)
    for axis in axes.flat:
        axis.grid(axis="y", alpha=0.22)
    figure.suptitle(
        "AIJ Case H — instantaneous versus time-averaged carrier flow\n"
        + ("STATIONARY" if stationarity_passed else "WITHHELD"),
        fontsize=14,
    )
    figure.savefig(figure_path, dpi=220)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "flow_statistical_residual": averaged["flow_statistical_residual"],
        "time_mean_metrics": averaged["metrics"],
        "time_mean_vs_instantaneous_change_percent": metric_changes,
        "runtimes_s": [item["runtime_s"] for item in cases],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
