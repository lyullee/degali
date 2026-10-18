"""Compare legacy and non-equilibrium-extrapolation D3Q19 open boundaries."""

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
BOUNDARY_CASES = (
    ("equilibrium_copy", "Equilibrium inlet / copy outlet"),
    ("nonequilibrium_extrapolation", "NEE velocity / pressure boundary"),
)


def _percent_change(candidate: float, reference: float) -> float:
    return 100.0 * (candidate - reference) / max(abs(reference), 1.0e-12)


def main() -> int:
    dataset = read_aij_case_h(DATA)
    base = replace(
        _config(dataset, CELL_SIZE_M),
        subgrid_model="none",
        smagorinsky_constant=0.0,
    )
    runs = []
    for scheme, label in BOUNDARY_CASES:
        config = replace(base, flow_boundary_scheme=scheme)
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset,
            result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((scheme, label, config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _scheme, _label, _config_value, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    cases = []
    for (scheme, label, config, result, comparison), mapping in zip(runs, point_maps):
        inlet_flux = float(np.sum(result.velocity_u_m_s[:, :, 0]))
        outlet_flux = float(np.sum(result.velocity_u_m_s[:, :, -1]))
        cases.append({
            "scheme": scheme,
            "label": label,
            "metrics": _metrics([mapping[identifier] for identifier in common_ids]),
            "runtime_s": result.diagnostics["total_runtime_s"],
            "flow_steps": result.diagnostics["flow_steps_completed"],
            "flow_physical_horizon_s": (
                result.diagnostics["flow_steps_completed"]
                * result.diagnostics["flow_physical_time_step_s"]
            ),
            "final_flow_residual": result.diagnostics["final_flow_residual"],
            "flow_converged": result.diagnostics["flow_converged"],
            "maximum_relative_lattice_density_variation": result.diagnostics[
                "maximum_relative_lattice_density_variation"
            ],
            "minimum_lee_streamwise_velocity_m_s": result.diagnostics[
                "minimum_lee_streamwise_velocity_m_s"
            ],
            "recirculating_lee_cell_count": result.diagnostics[
                "recirculating_lee_cell_count"
            ],
            "inlet_raw_cell_velocity_sum_m_s": inlet_flux,
            "outlet_raw_cell_velocity_sum_m_s": outlet_flux,
            "raw_cell_velocity_flux_imbalance_fraction": (
                abs(outlet_flux - inlet_flux) / max(abs(inlet_flux), 1.0e-12)
            ),
            "maximum_absolute_divergence_s_1": result.diagnostics[
                "maximum_absolute_divergence_s_1"
            ],
            "residual_history": result.flow_residual_history.tolist(),
        })

    legacy, candidate = cases
    error_changes = {
        name: _percent_change(
            float(candidate["metrics"][name]), float(legacy["metrics"][name])
        )
        for name in (
            "u_mae_normalized", "u_rmse_normalized",
            "w_mae_normalized", "w_rmse_normalized",
        )
    }
    correlation_changes = {
        name: float(candidate["metrics"][name]) - float(legacy["metrics"][name])
        for name in (
            "u_correlation", "w_correlation", "u_recirculation_sign_accuracy",
        )
    }
    residual_change = _percent_change(
        float(candidate["final_flow_residual"]),
        float(legacy["final_flow_residual"]),
    )
    density_change = _percent_change(
        float(candidate["maximum_relative_lattice_density_variation"]),
        float(legacy["maximum_relative_lattice_density_variation"]),
    )
    # Promotion requires a material accuracy gain, no correlation regression,
    # and a converged candidate.  This is intentionally stricter than merely
    # obtaining a smaller residual.
    promotion_reasons = []
    if not all(value <= -5.0 for value in error_changes.values()):
        promotion_reasons.append(
            "candidate does not reduce all four primary AIJ errors by at least 5%"
        )
    if any(value < -0.01 for value in correlation_changes.values()):
        promotion_reasons.append(
            "candidate reduces a correlation/sign score by more than 0.01"
        )
    if not bool(candidate["flow_converged"]):
        promotion_reasons.append("candidate flow residual did not converge")
    promoted = not promotion_reasons

    record = {
        "schema": "degali.aij-case-h-flow-boundary-comparison.v1",
        "status": "accepted" if promoted else "withheld_candidate_comparison",
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
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "cases": cases,
        "candidate_vs_legacy_error_change_percent": error_changes,
        "candidate_vs_legacy_correlation_change_absolute": correlation_changes,
        "candidate_vs_legacy_residual_change_percent": residual_change,
        "candidate_vs_legacy_density_variation_change_percent": density_change,
        "candidate_promoted": promoted,
        "promotion_reasons": promotion_reasons,
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_flow_boundary_comparison.json"
    figure_path = OUTPUT / "aij_case_h_flow_boundary_comparison.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    labels = [item["label"] for item in cases]
    x = np.arange(len(cases))
    width = 0.36
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), constrained_layout=True)
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
    axes[0, 0].set_title("Matched AIJ RMSE (lower is better)")
    axes[0, 0].set_ylabel("Normalised RMSE")
    axes[0, 0].legend()

    for item, color in zip(cases, ("#64748b", "#dc2626")):
        history = np.asarray(item["residual_history"])
        physical_time = (
            history[:, 0] * base.lattice_inlet_speed * CELL_SIZE_M
            / base.wind_speed_m_s
        )
        axes[0, 1].semilogy(
            physical_time, history[:, 1], "o-", label=item["label"], color=color,
        )
    axes[0, 1].axhline(
        base.flow_convergence_tolerance,
        color="black", linestyle="--", linewidth=1, label="declared tolerance",
    )
    axes[0, 1].set_title("Flow residual history")
    axes[0, 1].set_xlabel("Physical horizon (s)")
    axes[0, 1].set_ylabel("100-step velocity-change residual")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].bar(
        x,
        [item["maximum_relative_lattice_density_variation"] for item in cases],
        color=("#64748b", "#0f766e"),
    )
    axes[1, 0].set_title("Weak-compressibility density variation")
    axes[1, 0].set_ylabel("(max ρ − min ρ) / mean ρ")

    axes[1, 1].bar(
        x - width / 2,
        [item["runtime_s"] for item in cases],
        width, color="#475569", label="runtime",
    )
    flux_axis = axes[1, 1].twinx()
    flux_axis.bar(
        x + width / 2,
        [item["raw_cell_velocity_flux_imbalance_fraction"] for item in cases],
        width, color="#e11d48", alpha=0.75, label="raw flux imbalance",
    )
    axes[1, 1].set_title("Cost and raw boundary flux balance")
    axes[1, 1].set_ylabel("Runtime (s)")
    flux_axis.set_ylabel("Absolute imbalance fraction")

    for axis in (axes[0, 0], axes[1, 0], axes[1, 1]):
        axis.set_xticks(x, labels, rotation=8)
    for axis in axes.flat:
        axis.grid(axis="y", alpha=0.22)
    figure.suptitle(
        "AIJ Case H — D3Q19 open-boundary candidate comparison\n"
        + ("PROMOTED" if promoted else "WITHHELD"),
        fontsize=14,
    )
    figure.savefig(figure_path, dpi=220)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "candidate_promoted": promoted,
        "candidate_vs_legacy_error_change_percent": error_changes,
        "candidate_vs_legacy_residual_change_percent": residual_change,
        "promotion_reasons": promotion_reasons,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
