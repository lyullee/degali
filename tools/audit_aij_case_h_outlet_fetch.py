"""Audit downstream-domain sensitivity of the AIJ resolved 3-D flow case."""

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
    _relative_change,
)


CELL_SIZE_M = 0.020
DOMAIN_LENGTHS_M = (0.70, 0.90, 1.10)
BASE_FLOW_STEPS = 1_500


def main() -> int:
    dataset = read_aij_case_h(DATA)
    base = replace(
        _config(dataset, CELL_SIZE_M),
        subgrid_model="none",
        smagorinsky_constant=0.0,
        flow_boundary_scheme="equilibrium_copy",
    )
    runs = []
    for length in DOMAIN_LENGTHS_M:
        nx = int(round(length / CELL_SIZE_M))
        # Preserve the number of nominal domain flow-through times rather than
        # giving a longer domain less development time.
        flow_steps = int(round(BASE_FLOW_STEPS * length / DOMAIN_LENGTHS_M[0]))
        config = replace(
            base,
            domain_m=(length, base.domain_m[1], base.domain_m[2]),
            cells=(nx, base.cells[1], base.cells[2]),
            flow_steps=flow_steps,
            minimum_flow_steps=flow_steps,
        )
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset,
            result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _config_value, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    cases = []
    for (config, result, comparison), mapping in zip(runs, point_maps):
        cases.append({
            "domain_length_m": config.domain_m[0],
            "downstream_fetch_from_block_m": config.domain_m[0] - 0.25,
            "cells": list(config.cells),
            "flow_steps": config.flow_steps,
            "nominal_domain_flow_throughs": (
                config.flow_steps * config.lattice_inlet_speed / config.cells[0]
            ),
            "metrics": _metrics([mapping[identifier] for identifier in common_ids]),
            "runtime_s": result.diagnostics["total_runtime_s"],
            "final_flow_residual": result.diagnostics["final_flow_residual"],
            "flow_converged": result.diagnostics["flow_converged"],
            "minimum_lee_streamwise_velocity_m_s": result.diagnostics[
                "minimum_lee_streamwise_velocity_m_s"
            ],
            "recirculating_lee_cell_count": result.diagnostics[
                "recirculating_lee_cell_count"
            ],
            "maximum_relative_lattice_density_variation": result.diagnostics[
                "maximum_relative_lattice_density_variation"
            ],
        })

    medium, longest = cases[-2:]
    finest_pair_changes = {
        name: _relative_change(
            float(medium["metrics"][name]), float(longest["metrics"][name])
        )
        for name in (
            "u_mae_normalized", "u_rmse_normalized", "u_correlation",
            "w_mae_normalized", "w_rmse_normalized", "w_correlation",
            "u_recirculation_sign_accuracy",
        )
    }
    finest_pair_lee_change = _relative_change(
        float(medium["minimum_lee_streamwise_velocity_m_s"]),
        float(longest["minimum_lee_streamwise_velocity_m_s"]),
    )
    tolerance = 0.05
    gate_reasons = [
        f"0.90-to-1.10 m {name} change {change:.3%} exceeds {tolerance:.1%}"
        for name, change in finest_pair_changes.items() if change > tolerance
    ]
    if finest_pair_lee_change > tolerance:
        gate_reasons.append(
            "0.90-to-1.10 m minimum lee velocity change "
            f"{finest_pair_lee_change:.3%} exceeds {tolerance:.1%}"
        )
    if not bool(longest["flow_converged"]):
        gate_reasons.append("longest-domain flow residual did not converge")
    passed = not gate_reasons
    record = {
        "schema": "degali.aij-case-h-outlet-fetch-audit.v1",
        "status": "accepted" if passed else "withheld",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "comparison_scope": "AIJ center-plane mean U/W on common observation IDs",
        "cell_size_m": CELL_SIZE_M,
        "subgrid_model": "none",
        "flow_boundary_scheme": "equilibrium_copy",
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "relative_tolerance": tolerance,
        "cases": cases,
        "longest_pair_relative_changes": finest_pair_changes,
        "longest_pair_minimum_lee_velocity_change": finest_pair_lee_change,
        "passed": passed,
        "gate_reasons": gate_reasons,
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_outlet_fetch_audit.json"
    figure_path = OUTPUT / "aij_case_h_outlet_fetch_audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    lengths = np.asarray([item["domain_length_m"] for item in cases])
    figure, axes = plt.subplots(2, 2, figsize=(11.8, 8.0), constrained_layout=True)
    axes[0, 0].plot(
        lengths, [item["metrics"]["u_rmse_normalized"] for item in cases],
        "o-", label="$U/U_H$ RMSE",
    )
    axes[0, 0].plot(
        lengths, [item["metrics"]["w_rmse_normalized"] for item in cases],
        "s-", label="$W/U_H$ RMSE",
    )
    axes[0, 0].set(title="Matched AIJ error", ylabel="Normalised RMSE")
    axes[0, 0].legend()
    axes[0, 1].plot(
        lengths, [item["metrics"]["u_correlation"] for item in cases],
        "o-", label="$U$ correlation",
    )
    axes[0, 1].plot(
        lengths, [item["metrics"]["w_correlation"] for item in cases],
        "s-", label="$W$ correlation",
    )
    axes[0, 1].set(title="Matched AIJ correlation", ylabel="Pearson r")
    axes[0, 1].legend()
    axes[1, 0].semilogy(
        lengths, [item["final_flow_residual"] for item in cases],
        "o-", color="#dc2626", label="final residual",
    )
    axes[1, 0].axhline(
        base.flow_convergence_tolerance,
        linestyle="--", color="black", linewidth=1, label="declared tolerance",
    )
    axes[1, 0].set(title="Flow residual", ylabel="100-step velocity-change residual")
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].plot(
        lengths,
        [item["minimum_lee_streamwise_velocity_m_s"] for item in cases],
        "o-", color="#7c3aed",
    )
    axes[1, 1].set(title="Wake response", ylabel="Minimum lee U (m/s)")
    for axis in axes.flat:
        axis.set_xlabel("Streamwise domain length (m)")
        axis.grid(alpha=0.25)
    figure.suptitle(
        "AIJ Case H — downstream-fetch sensitivity — "
        + ("PASS" if passed else "WITHHELD"),
        fontsize=14,
    )
    figure.savefig(figure_path, dpi=220)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "gate_reasons": gate_reasons,
        "runtimes_s": [item["runtime_s"] for item in cases],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
