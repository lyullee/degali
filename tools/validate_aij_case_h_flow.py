"""Run the conditional AIJ Case H mean-flow component comparison."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from degali.addons.resolved_obstacle_flow_2d import (
    ResolvedObstacleFlow2DConfig,
    solve_resolved_obstacle_flow_2d,
)
from degali.addons.semi_fv_obstacle import (
    RectangularObstacle2D,
    SourceRateSchedule,
)
from degali.validation.aij_case_h import (
    aij_case_h_flow_comparison_record,
    compare_aij_case_h_centerplane_flow,
    read_aij_case_h,
)


DEFAULT_DATA = Path("reference/aij_case_h")
DEFAULT_OUTPUT = Path("outputs/obstacle-physics-rebuild-2026-10-10")
MODEL_X_ORIGIN_IN_EXPERIMENT_M = -0.25


def _config(dataset) -> tuple[ResolvedObstacleFlow2DConfig, dict[str, object]]:
    # The wind-tunnel profile starts at z=0.01 m.  A no-slip point at the
    # ground is prepended as an explicit boundary-completion assumption.
    profile_z = (0.0,) + tuple(item.z_m for item in dataset.approach_points)
    profile_u = (0.0,) + tuple(item.u_m_s for item in dataset.approach_points)
    config = ResolvedObstacleFlow2DConfig(
        length_m=0.70,
        height_m=0.40,
        nx=140,
        nz=80,
        wind_speed_m_s=dataset.reference_speed_m_s,
        kinematic_viscosity_m2_s=0.020,
        flow_steps=12_000,
        minimum_flow_steps=12_000,
        lattice_inlet_speed=0.040,
        duration_s=1.0e-5,
        scalar_time_step_s=1.0e-5,
        scalar_diffusivity_m2_s=1.0e-5,
        source_schedule=SourceRateSchedule((0.0, 1.0e-5), (0.0, 0.0)),
        source_x_m=0.30,
        source_height_m=0.005,
        source_sigma_m=0.01,
        obstacle=RectangularObstacle2D(
            0.20, 0.10, dataset.building_height_m, label="AIJ Case H building"
        ),
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
    )
    assumptions = {
        "model_x_origin_in_experiment_m": MODEL_X_ORIGIN_IN_EXPERIMENT_M,
        "ground_profile_completion": {"z_m": 0.0, "u_m_s": 0.0},
        "effective_kinematic_viscosity_m2_s": 0.020,
        "calibration_used": False,
        "selection_basis": (
            "declared stability-qualified effective viscosity; not fitted to Case H"
        ),
        "concentration_source_enabled": False,
    }
    return config, assumptions


def _plot(dataset, result, comparison, path: Path) -> None:
    height = dataset.reference_height_m
    measured_x = np.asarray([item.x_over_h for item in comparison.points])
    measured_z = np.asarray([item.z_over_h for item in comparison.points])
    measured_u = np.asarray([item.measured_u_over_uh for item in comparison.points])
    predicted_u = np.asarray([item.predicted_u_over_uh for item in comparison.points])
    measured_w = np.asarray([item.measured_w_over_uh for item in comparison.points])
    predicted_w = np.asarray([item.predicted_w_over_uh for item in comparison.points])

    model_x = (result.x_m + MODEL_X_ORIGIN_IN_EXPERIMENT_M) / height
    model_z = result.z_m / height
    model_u = np.ma.masked_where(
        result.solid_mask, result.velocity_u_m_s / dataset.reference_speed_m_s
    )
    limits = (-0.8, 1.8)
    levels = np.linspace(*limits, 27)
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.7), constrained_layout=True)

    measurement = axes[0, 0].scatter(
        measured_x, measured_z, c=measured_u, cmap="turbo", vmin=limits[0],
        vmax=limits[1], s=27, marker="s", linewidths=0,
    )
    axes[0, 0].set_title("AIJ wind-tunnel measurements: $U/U_H$")
    figure.colorbar(measurement, ax=axes[0, 0], label="$U/U_H$")

    contour = axes[0, 1].contourf(
        model_x, model_z, model_u, levels=levels, cmap="turbo", extend="both"
    )
    axes[0, 1].scatter(
        measured_x, measured_z, s=8, facecolors="none", edgecolors="white",
        linewidths=0.45, alpha=0.65,
    )
    axes[0, 1].set_title("DEGALI Tier-2 field: $U/U_H$")
    figure.colorbar(contour, ax=axes[0, 1], label="$U/U_H$")

    for axis in axes[0]:
        axis.add_patch(Rectangle((-0.5, 0.0), 0.5, 1.0, color="0.25", zorder=5))
        axis.set(xlim=(-1.1, 2.1), ylim=(0.0, 1.55), xlabel="$x/H$", ylabel="$z/H$")
        axis.grid(alpha=0.15)

    def parity(axis, measured, predicted, label, rmse, correlation):
        low = float(min(np.min(measured), np.min(predicted)))
        high = float(max(np.max(measured), np.max(predicted)))
        pad = 0.08 * max(high - low, 0.1)
        axis.scatter(measured, predicted, c=measured_x, cmap="coolwarm", s=23,
                     alpha=0.78, edgecolors="none")
        axis.plot((low - pad, high + pad), (low - pad, high + pad), "k--", lw=1)
        axis.set(
            xlim=(low - pad, high + pad), ylim=(low - pad, high + pad),
            xlabel=f"measured {label}", ylabel=f"predicted {label}",
        )
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.25)
        axis.text(
            0.03, 0.97, f"RMSE={rmse:.3f}\nr={correlation:.3f}",
            transform=axis.transAxes, va="top",
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
        )

    parity(
        axes[1, 0], measured_u, predicted_u, "$U/U_H$",
        comparison.u_rmse_normalized, comparison.u_correlation,
    )
    axes[1, 0].set_title("Streamwise velocity parity (220 points)")
    parity(
        axes[1, 1], measured_w, predicted_w, "$W/U_H$",
        comparison.w_rmse_normalized, comparison.w_correlation,
    )
    axes[1, 1].set_title("Vertical velocity parity (220 points)")

    figure.suptitle(
        "AIJ Case H conditional component check — concentration withheld\n"
        f"experiment Re$_H$={dataset.reynolds_number:,.0f}; "
        f"model effective Re$_H$={result.diagnostics['reynolds_number_obstacle_height']:.1f}; "
        f"runtime={result.diagnostics['total_runtime_s']:.1f} s",
        fontsize=13,
    )
    figure.savefig(path, dpi=220)
    plt.close(figure)


def run(data_directory: Path, output_directory: Path) -> dict[str, object]:
    dataset = read_aij_case_h(data_directory)
    config, assumptions = _config(dataset)
    result = solve_resolved_obstacle_flow_2d(config)
    comparison = compare_aij_case_h_centerplane_flow(
        dataset, result,
        model_x_origin_in_experiment_m=MODEL_X_ORIGIN_IN_EXPERIMENT_M,
    )
    record = aij_case_h_flow_comparison_record(comparison)
    record["dataset"] = {
        "record_doi": dataset.record_doi,
        "declared_readme_doi": dataset.declared_readme_doi,
        "approach_flow_sha256": dataset.approach_flow_sha256,
        "result_sha256": dataset.result_sha256,
        "readme_sha256": dataset.readme_sha256,
        "approach_point_count": len(dataset.approach_points),
        "observation_count": len(dataset.observations),
        "centerplane_flow_point_count": len(dataset.centerplane_flow_observations),
        "tracer": dataset.tracer,
        "experimental_reynolds_number": dataset.reynolds_number,
    }
    record["model_setup"] = assumptions
    record["numerics"] = result.diagnostics
    record["decision"] = {
        "status": "withheld_for_quantitative_hydrogen_use",
        "accepted_scope": "qualitative/conditional mean-flow component evidence",
        "reasons": [
            "model effective Reynolds number does not match the experiment",
            "flow residual did not meet the declared convergence tolerance",
            "two-dimensional solver omits finite-width lateral bypass",
            "ethylene concentration is not comparable to the per-span 2-D source",
            "cryogenic hydrogen heat transfer and buoyancy are absent",
        ],
    }

    output_directory.mkdir(parents=True, exist_ok=True)
    json_path = output_directory / "aij_case_h_flow_component_comparison.json"
    figure_path = output_directory / "aij_case_h_flow_component_comparison.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _plot(dataset, result, comparison, figure_path)
    return {
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": comparison.status,
        "metrics": record["metrics"],
        "runtime_s": result.diagnostics["total_runtime_s"],
        "quantitative_hydrogen_use": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.data, args.output), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
