"""Compare no-SGS, Smagorinsky, and WALE closures against AIJ Case H.

This is an evaluation-only component check.  The official hash-pinned AIJ
mean-flow observations are never used to tune the closure constants, and the
result cannot qualify hydrogen concentration predictions.
"""

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
CELL_SIZE_M = 0.020
FLOW_PHYSICAL_HORIZON_S = 0.2786377708978328
X_ORIGIN_M = -0.25
Y_ORIGIN_M = -0.15

MODEL_CASES = (
    ("none", "No SGS", 0.0, 0.0),
    ("smagorinsky", "Smagorinsky $C_s=0.12$", 0.12, 0.0),
    ("wale", "WALE $C_w=0.50$", 0.0, 0.50),
)


def _base_config(dataset) -> ResolvedObstacleFlow3DConfig:
    spacing = CELL_SIZE_M
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


def _percent_change(candidate: float, reference: float) -> float:
    """Signed change; negative means the candidate is smaller."""
    return 100.0 * (candidate - reference) / max(abs(reference), 1.0e-12)


def main() -> int:
    dataset = read_aij_case_h(DATA)
    base_config = _base_config(dataset)
    runs = []
    for model, label, smagorinsky_constant, wale_constant in MODEL_CASES:
        config = replace(
            base_config,
            subgrid_model=model,
            smagorinsky_constant=smagorinsky_constant,
            wale_constant=wale_constant,
        )
        result = solve_resolved_obstacle_flow_3d(config)
        comparison = compare_aij_case_h_centerplane_flow_3d(
            dataset,
            result,
            model_x_origin_in_experiment_m=X_ORIGIN_M,
            model_y_origin_in_experiment_m=Y_ORIGIN_M,
        )
        runs.append((model, label, config, result, comparison))

    point_maps = [
        {point.point_id: point for point in comparison.points}
        for _model, _label, _config, _result, comparison in runs
    ]
    common_ids = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    records = []
    for (model, label, config, result, comparison), mapping in zip(runs, point_maps):
        diagnostics = result.diagnostics
        records.append({
            "model": model,
            "label": label,
            "smagorinsky_constant": config.smagorinsky_constant,
            "wale_constant": config.wale_constant,
            "metrics": _metrics([mapping[identifier] for identifier in common_ids]),
            "runtime_s": diagnostics["total_runtime_s"],
            "final_flow_residual": diagnostics["final_flow_residual"],
            "flow_converged": diagnostics["flow_converged"],
            "maximum_lattice_relaxation_time": diagnostics[
                "maximum_lattice_relaxation_time"
            ],
            "maximum_eddy_viscosity_m2_s": diagnostics[
                "maximum_eddy_viscosity_m2_s"
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
        })

    by_model = {item["model"]: item for item in records}
    smag = by_model["smagorinsky"]
    wale = by_model["wale"]
    metric_changes = {
        name: _percent_change(float(wale["metrics"][name]), float(smag["metrics"][name]))
        for name in ("u_mae_normalized", "u_rmse_normalized", "w_mae_normalized", "w_rmse_normalized")
    }
    correlation_changes = {
        name: float(wale["metrics"][name]) - float(smag["metrics"][name])
        for name in ("u_correlation", "w_correlation", "u_recirculation_sign_accuracy")
    }
    lower_error_count = sum(value < 0.0 for value in metric_changes.values())
    higher_correlation_count = sum(value > 0.0 for value in correlation_changes.values())
    wale_improves_all_primary_errors = lower_error_count == len(metric_changes)

    record = {
        "schema": "degali.aij-case-h-subgrid-model-comparison.v1",
        "status": "withheld_candidate_comparison",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "comparison_scope": "AIJ center-plane mean U/W on common observation IDs",
        "cell_size_m": CELL_SIZE_M,
        "cells": list(base_config.cells),
        "flow_steps": base_config.flow_steps,
        "flow_physical_horizon_s": (
            base_config.flow_steps
            * base_config.lattice_inlet_speed * CELL_SIZE_M
            / base_config.wind_speed_m_s
        ),
        "calibration_used": False,
        "common_point_count": len(common_ids),
        "common_point_ids": common_ids,
        "cases": records,
        "wale_vs_smagorinsky_primary_error_change_percent": metric_changes,
        "wale_vs_smagorinsky_correlation_change_absolute": correlation_changes,
        "wale_lower_primary_error_count": lower_error_count,
        "wale_higher_correlation_or_sign_count": higher_correlation_count,
        "wale_improves_all_primary_errors": wale_improves_all_primary_errors,
        "default_model_changed": False,
        "decision": (
            "retain WALE as opt-in pending grid/time and wall-flow verification"
        ),
        "quantitative_hydrogen_validation": False,
        "concentration_validation_withheld": True,
        "limitations": [
            "single grid and a deliberately fixed closure constant per candidate",
            "flow did not receive case-specific closure calibration",
            "mean velocity component evidence does not validate scalar dispersion",
            "AIJ Case H uses neutral ethylene, not cryogenic hydrogen",
        ],
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_subgrid_model_comparison.json"
    figure_path = OUTPUT / "aij_case_h_subgrid_model_comparison.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    labels = [item["label"] for item in records]
    colors = ("#64748b", "#2563eb", "#dc2626")
    x = np.arange(len(records))
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), constrained_layout=True)
    width = 0.36
    axes[0, 0].bar(
        x - width / 2,
        [item["metrics"]["u_rmse_normalized"] for item in records],
        width,
        label="$U/U_H$",
        color="#2563eb",
    )
    axes[0, 0].bar(
        x + width / 2,
        [item["metrics"]["w_rmse_normalized"] for item in records],
        width,
        label="$W/U_H$",
        color="#f59e0b",
    )
    axes[0, 0].set_title("Matched mean-velocity RMSE (lower is better)")
    axes[0, 0].set_ylabel("Normalised RMSE")
    axes[0, 0].legend()

    axes[0, 1].bar(
        x - width / 2,
        [item["metrics"]["u_correlation"] for item in records],
        width,
        label="$U$ correlation",
        color="#0f766e",
    )
    axes[0, 1].bar(
        x + width / 2,
        [item["metrics"]["w_correlation"] for item in records],
        width,
        label="$W$ correlation",
        color="#7c3aed",
    )
    axes[0, 1].set_title("Matched correlation (higher is better)")
    axes[0, 1].set_ylabel("Pearson r")
    axes[0, 1].legend()

    axes[1, 0].bar(
        x,
        [item["minimum_lee_streamwise_velocity_m_s"] for item in records],
        color=colors,
    )
    axes[1, 0].set_title("Resolved reverse flow behind block")
    axes[1, 0].set_ylabel("Minimum lee $U$ (m/s)")
    axes[1, 0].axhline(0.0, color="black", linewidth=0.8)

    axes[1, 1].bar(
        x - width / 2,
        [item["runtime_s"] for item in records],
        width,
        color="#475569",
        label="runtime",
    )
    eddy_axis = axes[1, 1].twinx()
    eddy_axis.bar(
        x + width / 2,
        [item["maximum_eddy_viscosity_m2_s"] for item in records],
        width,
        color="#e11d48",
        alpha=0.75,
        label="max eddy viscosity",
    )
    axes[1, 1].set_title("Cost and closure strength")
    axes[1, 1].set_ylabel("Runtime (s)")
    eddy_axis.set_ylabel("Maximum eddy viscosity (m²/s)")

    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=10)
        axis.grid(axis="y", alpha=0.22)
    figure.suptitle(
        "AIJ Case H — SGS candidate comparison (evaluation only)", fontsize=14
    )
    figure.savefig(figure_path, dpi=220)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "common_point_count": len(common_ids),
        "wale_vs_smagorinsky_primary_error_change_percent": metric_changes,
        "wale_vs_smagorinsky_correlation_change_absolute": correlation_changes,
        "runtimes_s": {item["model"]: item["runtime_s"] for item in records},
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
