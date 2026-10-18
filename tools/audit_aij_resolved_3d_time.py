"""Flow-horizon convergence audit for the resolved 3-D AIJ prototype."""

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


FLOW_STEPS = (1_500, 3_000, 6_000)
CELL_SIZE_M = 0.020


def main() -> int:
    dataset = read_aij_case_h(DATA)
    # The obstacle wake is inherently unsteady; use a trailing-window
    # statistical average in addition to the instantaneous residual so that a
    # persistent wake is not incorrectly treated as a failed transient.
    base = replace(
        _config(dataset, CELL_SIZE_M),
        flow_average_last_steps=1_000,
        flow_statistical_convergence_tolerance=0.01,
    )
    runs = []
    for steps in FLOW_STEPS:
        config = replace(base, flow_steps=steps, minimum_flow_steps=steps)
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset, result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _config_value, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    records = []
    for (config, result, comparison), mapping in zip(runs, point_maps):
        records.append({
            "flow_steps": config.flow_steps,
            "flow_physical_horizon_s": (
                config.flow_steps * config.lattice_inlet_speed * CELL_SIZE_M
                / config.wind_speed_m_s
            ),
            "matched_metrics": _metrics([
                mapping[identifier] for identifier in common_ids
            ]),
            "runtime_s": result.diagnostics["total_runtime_s"],
            "final_flow_residual": result.diagnostics["final_flow_residual"],
            "flow_converged": result.diagnostics["flow_converged"],
            "flow_statistical_residual": result.diagnostics[
                "flow_statistical_residual"
            ],
            "minimum_lee_streamwise_velocity_m_s": result.diagnostics[
                "minimum_lee_streamwise_velocity_m_s"
            ],
            "recirculating_lee_cell_count": result.diagnostics[
                "recirculating_lee_cell_count"
            ],
            "maximum_absolute_divergence_s_1": result.diagnostics[
                "maximum_absolute_divergence_s_1"
            ],
        })
    middle, longest = records[-2], records[-1]
    metric_changes = {}
    for name in (
        "u_mae_normalized", "u_rmse_normalized", "u_correlation",
        "w_mae_normalized", "w_rmse_normalized", "w_correlation",
        "u_recirculation_sign_accuracy",
    ):
        metric_changes[name] = _relative_change(
            float(middle["matched_metrics"][name]),
            float(longest["matched_metrics"][name]),
        )
    lee_change = _relative_change(
        float(middle["minimum_lee_streamwise_velocity_m_s"]),
        float(longest["minimum_lee_streamwise_velocity_m_s"]),
    )
    tolerance = 0.10
    gate_reasons = [
        f"longest-pair {name} change {change:.3%} exceeds {tolerance:.1%}"
        for name, change in metric_changes.items() if change > tolerance
    ]
    if lee_change > tolerance:
        gate_reasons.append(
            f"longest-pair minimum lee velocity change {lee_change:.3%} "
            f"exceeds {tolerance:.1%}"
        )
    if not bool(longest["flow_converged"]):
        statistical_residual = longest["flow_statistical_residual"]
        if statistical_residual is None or float(statistical_residual) > 0.01:
            gate_reasons.append(
                "longest flow horizon did not meet instantaneous or statistical "
                "residual tolerance"
            )
    passed = not gate_reasons
    record = {
        "schema": "degali.aij-case-h-resolved-3d-time-audit.v1",
        "status": "accepted" if passed else "withheld",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "cell_size_m": CELL_SIZE_M,
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "relative_tolerance": tolerance,
        "cases": records,
        "longest_pair_relative_changes": metric_changes,
        "longest_pair_minimum_lee_velocity_change": lee_change,
        "passed": passed,
        "gate_reasons": gate_reasons,
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_resolved_3d_time_audit.json"
    figure_path = OUTPUT / "aij_case_h_resolved_3d_time_audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    horizons = np.asarray([item["flow_physical_horizon_s"] for item in records])
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    axes[0, 0].plot(
        horizons,
        [item["matched_metrics"]["u_rmse_normalized"] for item in records],
        "o-", label="$U/U_H$ RMSE",
    )
    axes[0, 0].plot(
        horizons,
        [item["matched_metrics"]["w_rmse_normalized"] for item in records],
        "s-", label="$W/U_H$ RMSE",
    )
    axes[0, 0].set(ylabel="Normalised RMSE", title="Matched AIJ error")
    axes[0, 0].legend()
    axes[0, 1].plot(
        horizons,
        [item["matched_metrics"]["u_correlation"] for item in records],
        "o-", label="$U$ correlation",
    )
    axes[0, 1].plot(
        horizons,
        [item["matched_metrics"]["w_correlation"] for item in records],
        "s-", label="$W$ correlation",
    )
    axes[0, 1].set(ylabel="Pearson r", title="Matched AIJ correlation")
    axes[0, 1].legend()
    axes[1, 0].plot(
        horizons,
        [item["minimum_lee_streamwise_velocity_m_s"] for item in records],
        "o-", color="#7c3aed",
    )
    axes[1, 0].set(ylabel="Minimum lee U (m/s)", title="Wake time sensitivity")
    axes[1, 1].semilogy(
        horizons,
        [item["final_flow_residual"] for item in records],
        "s--", color="#dc2626", label="flow residual",
    )
    axes[1, 1].set(ylabel="Final flow residual", title="Residual decay")
    for axis in axes.flat:
        axis.set_xlabel("Flow physical horizon (s)")
        axis.grid(alpha=0.25)
    figure.suptitle(
        "AIJ Case H resolved 3-D flow-horizon audit — "
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
        "longest_runtime_s": longest["runtime_s"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
