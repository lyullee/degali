"""Compare 2-D viscous and finite-width 3-D potential baselines on AIJ Case H."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.resolved_obstacle_flow_2d import (
    ResolvedObstacleFlow2DConfig,
    solve_resolved_obstacle_flow_2d,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import RectangularObstacle2D, SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.validation.aij_case_h import (
    aij_case_h_flow_comparison_record,
    compare_aij_case_h_centerplane_flow,
    compare_aij_case_h_centerplane_flow_3d,
    read_aij_case_h,
)
from degali.validation.obstacle_flow_canonical import (
    D3Q19PoiseuilleVerificationConfig,
    verify_d3q19_poiseuille_channel,
)


DATA = Path("reference/aij_case_h")
OUTPUT = Path("outputs/obstacle-physics-rebuild-2026-10-10")
X_ORIGIN_M = -0.25
Y_ORIGIN_M = -0.15


def _profile(dataset) -> tuple[tuple[float, ...], tuple[float, ...]]:
    return (
        (0.0,) + tuple(item.z_m for item in dataset.approach_points),
        (0.0,) + tuple(item.u_m_s for item in dataset.approach_points),
    )


def _zero_schedule() -> SourceRateSchedule:
    return SourceRateSchedule((0.0, 1.0e-5), (0.0, 0.0))


def _solve_2d(dataset):
    profile_z, profile_u = _profile(dataset)
    return solve_resolved_obstacle_flow_2d(ResolvedObstacleFlow2DConfig(
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
        source_schedule=_zero_schedule(),
        source_x_m=0.30,
        source_height_m=0.005,
        source_sigma_m=0.01,
        obstacle=RectangularObstacle2D(0.20, 0.10, 0.20, label="AIJ Case H 2-D block"),
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
    ))


def _solve_3d(dataset):
    profile_z, profile_u = _profile(dataset)
    return solve_fast_obstacle_transport_3d(FastObstacleTransport3DConfig(
        domain_m=(0.70, 0.30, 0.40),
        cells=(70, 30, 40),
        duration_s=1.0e-5,
        time_step_s=1.0e-5,
        wind_speed_m_s=dataset.reference_speed_m_s,
        diffusivity_m2_s=1.0e-5,
        source_schedule=_zero_schedule(),
        source_position_m=(0.30, 0.15, 0.005),
        source_sigma_m=0.01,
        obstacles=(AxisAlignedCuboid(
            0.15, 0.25, 0.10, 0.20, 0.0, 0.20, "AIJ Case H finite block"
        ),),
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
    ))


def _solve_resolved_3d(dataset):
    profile_z, profile_u = _profile(dataset)
    return solve_resolved_obstacle_flow_3d(ResolvedObstacleFlow3DConfig(
        domain_m=(0.70, 0.30, 0.40),
        cells=(35, 15, 20),
        wind_speed_m_s=dataset.reference_speed_m_s,
        kinematic_viscosity_m2_s=0.020,
        flow_steps=1_500,
        minimum_flow_steps=1_500,
        lattice_inlet_speed=0.030,
        duration_s=1.0e-5,
        scalar_time_step_s=1.0e-5,
        scalar_diffusivity_m2_s=1.0e-5,
        source_schedule=_zero_schedule(),
        source_position_m=(0.30, 0.15, 0.01),
        source_sigma_m=0.02,
        obstacles=(AxisAlignedCuboid(
            0.15, 0.25, 0.10, 0.20, 0.0, 0.20, "AIJ Case H finite block"
        ),),
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
        smagorinsky_constant=0.12,
    ))


def _matched_metrics(*comparisons):
    point_maps = [
        {point.point_id: point for point in comparison.points}
        for comparison in comparisons
    ]
    common = sorted(set.intersection(*(set(mapping) for mapping in point_maps)))
    output = []
    for mapping in point_maps:
        points = [mapping[identifier] for identifier in common]
        measured_u = np.asarray([point.measured_u_over_uh for point in points])
        predicted_u = np.asarray([point.predicted_u_over_uh for point in points])
        measured_w = np.asarray([point.measured_w_over_uh for point in points])
        predicted_w = np.asarray([point.predicted_w_over_uh for point in points])
        sign_mask = np.abs(measured_u) > 0.02
        output.append({
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
        })
    return common, output


def _add_building(axis) -> None:
    axis.add_patch(Rectangle((-0.5, 0.0), 0.5, 1.0, color="0.23", zorder=5))
    axis.set(xlim=(-1.1, 2.1), ylim=(0.0, 1.55), xlabel="$x/H$", ylabel="$z/H$")
    axis.grid(alpha=0.14)


def _plot(
    dataset, result_2d, comparison_2d, result_3d, comparison_3d,
    result_resolved_3d, comparison_resolved_3d, path,
):
    height = dataset.reference_height_m
    x_measure = np.asarray([item.x_over_h for item in comparison_2d.points])
    z_measure = np.asarray([item.z_over_h for item in comparison_2d.points])
    measured_u = np.asarray([item.measured_u_over_uh for item in comparison_2d.points])
    measured_w = np.asarray([item.measured_w_over_uh for item in comparison_2d.points])
    x_2d = (result_2d.x_m + X_ORIGIN_M) / height
    z_2d = result_2d.z_m / height
    y_mid = int(np.argmin(np.abs(result_3d.y_m - 0.15)))
    y_mid_resolved = int(np.argmin(np.abs(result_resolved_3d.y_m - 0.15)))
    x_3d = (result_3d.x_m + X_ORIGIN_M) / height
    z_3d = result_3d.z_m / height
    x_resolved_3d = (result_resolved_3d.x_m + X_ORIGIN_M) / height
    z_resolved_3d = result_resolved_3d.z_m / height

    fields = (
        np.ma.masked_where(
            result_2d.solid_mask,
            result_2d.velocity_u_m_s / dataset.reference_speed_m_s,
        ),
        np.ma.masked_where(
            result_3d.solid_mask[:, y_mid, :],
            result_3d.velocity_u_m_s[:, y_mid, :] / dataset.reference_speed_m_s,
        ),
        np.ma.masked_where(
            result_resolved_3d.solid_mask[:, y_mid_resolved, :],
            result_resolved_3d.velocity_u_m_s[:, y_mid_resolved, :]
            / dataset.reference_speed_m_s,
        ),
        np.ma.masked_where(
            result_2d.solid_mask,
            result_2d.velocity_w_m_s / dataset.reference_speed_m_s,
        ),
        np.ma.masked_where(
            result_3d.solid_mask[:, y_mid, :],
            result_3d.velocity_w_m_s[:, y_mid, :] / dataset.reference_speed_m_s,
        ),
        np.ma.masked_where(
            result_resolved_3d.solid_mask[:, y_mid_resolved, :],
            result_resolved_3d.velocity_w_m_s[:, y_mid_resolved, :]
            / dataset.reference_speed_m_s,
        ),
    )
    figure, axes = plt.subplots(2, 4, figsize=(20.0, 8.4), constrained_layout=True)
    u_limits = (-0.8, 1.8)
    w_limits = (-0.8, 0.8)

    measured = axes[0, 0].scatter(
        x_measure, z_measure, c=measured_u, cmap="turbo",
        vmin=u_limits[0], vmax=u_limits[1], marker="s", s=24, linewidths=0,
    )
    axes[0, 0].set_title("AIJ measured $U/U_H$")
    figure.colorbar(measured, ax=axes[0, 0], label="$U/U_H$")
    for axis, field, coordinates, comparison, name in (
        (axes[0, 1], fields[0], (x_2d, z_2d), comparison_2d, "2-D viscous"),
        (axes[0, 2], fields[1], (x_3d, z_3d), comparison_3d, "3-D finite-width potential"),
        (
            axes[0, 3], fields[2], (x_resolved_3d, z_resolved_3d),
            comparison_resolved_3d, "3-D finite-width viscous",
        ),
    ):
        contour = axis.contourf(
            *coordinates, field, levels=np.linspace(*u_limits, 27),
            cmap="turbo", extend="both",
        )
        axis.set_title(
            f"{name} $U/U_H$\nRMSE={comparison.u_rmse_normalized:.3f}, "
            f"r={comparison.u_correlation:.3f}"
        )
        figure.colorbar(contour, ax=axis, label="$U/U_H$")

    measured = axes[1, 0].scatter(
        x_measure, z_measure, c=measured_w, cmap="coolwarm",
        vmin=w_limits[0], vmax=w_limits[1], marker="s", s=24, linewidths=0,
    )
    axes[1, 0].set_title("AIJ measured $W/U_H$")
    figure.colorbar(measured, ax=axes[1, 0], label="$W/U_H$")
    for axis, field, coordinates, comparison, name in (
        (axes[1, 1], fields[3], (x_2d, z_2d), comparison_2d, "2-D viscous"),
        (axes[1, 2], fields[4], (x_3d, z_3d), comparison_3d, "3-D finite-width potential"),
        (
            axes[1, 3], fields[5], (x_resolved_3d, z_resolved_3d),
            comparison_resolved_3d, "3-D finite-width viscous",
        ),
    ):
        contour = axis.contourf(
            *coordinates, field, levels=np.linspace(*w_limits, 25),
            cmap="coolwarm", extend="both",
        )
        axis.set_title(
            f"{name} $W/U_H$\nRMSE={comparison.w_rmse_normalized:.3f}, "
            f"r={comparison.w_correlation:.3f}"
        )
        figure.colorbar(contour, ax=axis, label="$W/U_H$")

    for axis in axes.flat:
        _add_building(axis)
    figure.suptitle(
        "AIJ Case H: finite width plus viscous wake improves the velocity field\n"
        "All tiers remain conditional; tracer concentration and H2 validation are withheld",
        fontsize=14,
    )
    figure.savefig(path, dpi=220)
    plt.close(figure)


def main() -> int:
    dataset = read_aij_case_h(DATA)
    d3q19_channel = verify_d3q19_poiseuille_channel(
        D3Q19PoiseuilleVerificationConfig(
            nx=12, ny=4, nz=12, steps=6_000,
            convergence_check_interval=500,
        )
    )
    result_2d = _solve_2d(dataset)
    comparison_2d = compare_aij_case_h_centerplane_flow(
        dataset, result_2d, model_x_origin_in_experiment_m=X_ORIGIN_M,
    )
    result_3d = _solve_3d(dataset)
    comparison_3d = compare_aij_case_h_centerplane_flow_3d(
        dataset, result_3d,
        model_x_origin_in_experiment_m=X_ORIGIN_M,
        model_y_origin_in_experiment_m=Y_ORIGIN_M,
    )
    result_resolved_3d = _solve_resolved_3d(dataset)
    comparison_resolved_3d = compare_aij_case_h_centerplane_flow_3d(
        dataset, result_resolved_3d,
        model_x_origin_in_experiment_m=X_ORIGIN_M,
        model_y_origin_in_experiment_m=Y_ORIGIN_M,
    )
    record_2d = aij_case_h_flow_comparison_record(comparison_2d)
    record_3d = aij_case_h_flow_comparison_record(comparison_3d)
    record_resolved_3d = aij_case_h_flow_comparison_record(comparison_resolved_3d)
    common_ids, matched = _matched_metrics(
        comparison_2d, comparison_3d, comparison_resolved_3d,
    )
    matched_2d, matched_3d, matched_resolved_3d = matched
    reduction_u = 1.0 - (
        matched_resolved_3d["u_rmse_normalized"]
        / matched_2d["u_rmse_normalized"]
    )
    reduction_w = 1.0 - (
        matched_resolved_3d["w_rmse_normalized"]
        / matched_2d["w_rmse_normalized"]
    )
    record = {
        "schema": "degali.aij-case-h-fidelity-comparison.v2",
        "status": "withheld_for_quantitative_hydrogen_use",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
            "readme_sha256": dataset.readme_sha256,
        },
        "d3q19_canonical_poiseuille": {
            "schema": d3q19_channel.schema,
            "passed": d3q19_channel.passed,
            "relative_l2_error": d3q19_channel.relative_l2_error,
            "relative_linf_error": d3q19_channel.relative_linf_error,
            "relative_mass_flux_error": d3q19_channel.relative_mass_flux_error,
            "maximum_horizontal_nonuniformity": (
                d3q19_channel.maximum_horizontal_nonuniformity
            ),
            "maximum_relative_density_variation": (
                d3q19_channel.maximum_relative_density_variation
            ),
            "runtime_s": d3q19_channel.runtime_s,
            "diagnostics": d3q19_channel.diagnostics,
        },
        "two_dimensional_viscous": {
            **record_2d,
            "runtime_s": result_2d.diagnostics["total_runtime_s"],
            "numerics": result_2d.diagnostics,
        },
        "three_dimensional_finite_width_potential": {
            **record_3d,
            "runtime_s": result_3d.diagnostics["total_runtime_s"],
            "numerics": result_3d.diagnostics,
        },
        "three_dimensional_finite_width_viscous": {
            **record_resolved_3d,
            "runtime_s": result_resolved_3d.diagnostics["total_runtime_s"],
            "numerics": result_resolved_3d.diagnostics,
        },
        "matched_point_comparison": {
            "point_count": len(common_ids),
            "point_ids": common_ids,
            "two_dimensional_viscous": matched_2d,
            "three_dimensional_finite_width_potential": matched_3d,
            "three_dimensional_finite_width_viscous": matched_resolved_3d,
        },
        "matched_rmse_reduction_from_2d_to_resolved_3d": {
            "u_fraction": reduction_u,
            "u_percent": 100.0 * reduction_u,
            "w_fraction": reduction_w,
            "w_percent": 100.0 * reduction_w,
        },
        "interpretation": [
            "finite lateral bypass materially reduces the 2-D blockage magnitude error",
            "D3Q19 no-slip momentum further improves the matched mean velocity comparison",
            "the D3Q19 boundary and Smagorinsky wall treatment remain unverified",
            "neither comparison validates ethylene concentration or cryogenic hydrogen",
            "the next gate is grid/time convergence plus energy-density-buoyancy coupling",
        ],
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_2d_3d_fidelity_comparison.json"
    figure_path = OUTPUT / "aij_case_h_2d_3d_fidelity_comparison.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _plot(
        dataset, result_2d, comparison_2d, result_3d, comparison_3d,
        result_resolved_3d, comparison_resolved_3d, figure_path,
    )
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "u_rmse_2d": comparison_2d.u_rmse_normalized,
        "u_rmse_3d": comparison_3d.u_rmse_normalized,
        "u_rmse_resolved_3d": comparison_resolved_3d.u_rmse_normalized,
        "matched_point_count": len(common_ids),
        "matched_u_rmse_reduction_percent": 100.0 * reduction_u,
        "w_rmse_2d": comparison_2d.w_rmse_normalized,
        "w_rmse_3d": comparison_3d.w_rmse_normalized,
        "w_rmse_resolved_3d": comparison_resolved_3d.w_rmse_normalized,
        "matched_w_rmse_reduction_percent": 100.0 * reduction_w,
        "runtime_2d_s": result_2d.diagnostics["total_runtime_s"],
        "runtime_3d_s": result_3d.diagnostics["total_runtime_s"],
        "runtime_resolved_3d_s": result_resolved_3d.diagnostics["total_runtime_s"],
        "d3q19_poiseuille_passed": d3q19_channel.passed,
        "status": record["status"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
