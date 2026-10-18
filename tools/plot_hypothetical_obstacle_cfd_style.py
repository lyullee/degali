"""Render a CFD-style obstacle screening comparison from DEGALI semi-FV.

The figure uses DEGALI's conservative x-z semi-finite-volume obstacle
transport operator.  It is a reduced-order screening result, not an
obstacle-resolved RANS/LES calculation and not a validated wake prediction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.semi_fv_obstacle import (
    RectangularObstacle2D,
    SemiFVConfig,
    SemiFVResult,
    solve_semi_fv_obstacle,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "hypothetical-obstacle-cfd-style-2026-10-10"

R_GAS_J_MOL_K = 8.314462618
TEMPERATURE_K = 293.15
PRESSURE_PA = 101_325.0
MOLAR_MASS_H2_KG_MOL = 2.01588e-3


def _base_config(duration_s: float, *, obstacle: bool) -> SemiFVConfig:
    declared_obstacle = None
    if obstacle:
        declared_obstacle = RectangularObstacle2D(
            distance_m=9.0,
            width_m=2.5,
            height_m=3.0,
            base_height_m=0.0,
            label="hypothetical-cuboid-section",
        )
    return SemiFVConfig(
        length_m=24.0,
        height_m=8.0,
        nx=144,
        nz=64,
        time_step_s=0.025,
        duration_s=duration_s,
        wind_speed_m_s=1.8,
        diffusivity_m2_s=0.08,
        gravitational_settling_m_s=0.0,
        source_rate_kg_s=0.01,
        source_height_m=0.60,
        source_sigma_m=0.35,
        obstacle=declared_obstacle,
    )


def _kg_m3_to_vol_percent(concentration_kg_m3: np.ndarray) -> np.ndarray:
    """Ideal-gas screening conversion at the declared ambient state."""
    factor = (
        R_GAS_J_MOL_K
        * TEMPERATURE_K
        / (PRESSURE_PA * MOLAR_MASS_H2_KG_MOL)
        * 100.0
    )
    return np.clip(concentration_kg_m3 * factor, 0.0, 100.0)


def _timed_solve(config: SemiFVConfig) -> tuple[SemiFVResult, float]:
    started = time.perf_counter()
    result = solve_semi_fv_obstacle(config)
    return result, time.perf_counter() - started


def _threshold_extent(
    field_vol_percent: np.ndarray,
    x_m: np.ndarray,
    threshold: float = 4.0,
) -> float | None:
    vertical_maximum = np.nanmax(field_vol_percent, axis=0)
    mask = vertical_maximum >= threshold
    if not np.any(mask):
        return None
    return float(np.max(x_m[mask]))


def _draw_obstacle(axis, obstacle: RectangularObstacle2D) -> None:
    from matplotlib.patches import Rectangle

    patch = Rectangle(
        (obstacle.x_min_m, obstacle.base_height_m),
        obstacle.width_m,
        obstacle.height_m,
        facecolor="#263238",
        edgecolor="white",
        linewidth=1.4,
        hatch="////",
        zorder=8,
    )
    axis.add_patch(patch)


def _draw_field(
    axis,
    result: SemiFVResult,
    config: SemiFVConfig,
    *,
    title: str,
    levels: np.ndarray,
    cmap,
    norm,
    show_y_label: bool,
):
    field = _kg_m3_to_vol_percent(result.concentration_kg_m3)
    plotted = np.ma.masked_less_equal(np.ma.masked_invalid(field), 0.0)
    image = axis.contourf(
        result.x_m,
        result.z_m,
        plotted,
        levels=levels,
        cmap=cmap,
        norm=norm,
        extend="max",
    )
    if float(np.nanmax(field)) >= 4.0:
        axis.contour(
            result.x_m,
            result.z_m,
            plotted,
            levels=[4.0],
            colors="white",
            linewidths=3.6,
            zorder=5,
        )
        axis.contour(
            result.x_m,
            result.z_m,
            plotted,
            levels=[4.0],
            colors="#111827",
            linewidths=1.4,
            zorder=6,
        )
    if config.obstacle is not None:
        _draw_obstacle(axis, config.obstacle)
    axis.scatter(
        [result.x_m[0]],
        [config.source_height_m],
        marker="*",
        s=115,
        color="white",
        edgecolor="#111827",
        linewidth=0.9,
        zorder=10,
    )
    axis.annotate(
        "",
        xy=(4.0, 7.15),
        xytext=(1.0, 7.15),
        arrowprops={"arrowstyle": "-|>", "lw": 2.0, "color": "#111827"},
    )
    axis.text(1.0, 7.42, "wind 1.8 m s$^{-1}$", fontsize=8.5, weight="bold")
    axis.set_title(title, weight="bold")
    axis.set_xlim(0.0, config.length_m)
    axis.set_ylim(0.0, config.height_m)
    axis.set_xlabel("Downwind distance, x (m)")
    if show_y_label:
        axis.set_ylabel("Height, z (m)")
    axis.set_facecolor("#f8fafc")
    axis.grid(False)
    return image, field


def build(output_dir: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    wall_started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)

    snapshot_times = (3.0, 6.0, 10.0)
    obstacle_cases: list[tuple[float, SemiFVConfig, SemiFVResult, float]] = []
    for snapshot_time in snapshot_times:
        config = _base_config(snapshot_time, obstacle=True)
        result, runtime = _timed_solve(config)
        obstacle_cases.append((snapshot_time, config, result, runtime))

    open_config = _base_config(snapshot_times[-1], obstacle=False)
    open_result, open_runtime = _timed_solve(open_config)
    solve_runtime = sum(item[3] for item in obstacle_cases) + open_runtime

    levels = np.asarray((0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_under("#f8fafc")
    cmap.set_bad("#263238")
    norm = colors.LogNorm(vmin=float(levels[0]), vmax=float(levels[-1]))

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.titlesize": 16,
        "figure.titleweight": "bold",
        "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.0), constrained_layout=False)

    last_image = None
    obstacle_fields = []
    for column, (snapshot_time, config, result, runtime) in enumerate(obstacle_cases):
        last_image, field = _draw_field(
            axes[0, column],
            result,
            config,
            title=f"Obstacle case: t = {snapshot_time:.0f} s  |  solve {runtime:.2f} s",
            levels=levels,
            cmap=cmap,
            norm=norm,
            show_y_label=column == 0,
        )
        obstacle_fields.append(field)

    last_time, last_config, last_result, last_runtime = obstacle_cases[-1]
    last_obstacle_field = obstacle_fields[-1]
    _, open_field = _draw_field(
        axes[1, 0],
        open_result,
        open_config,
        title=f"Open field reference: t = {last_time:.0f} s  |  solve {open_runtime:.2f} s",
        levels=levels,
        cmap=cmap,
        norm=norm,
        show_y_label=True,
    )

    difference_axis = axes[1, 1]
    difference = last_obstacle_field - open_field
    difference_masked = np.ma.masked_invalid(difference)
    finite_difference = np.abs(difference[np.isfinite(difference)])
    difference_limit = max(1.0, float(np.percentile(finite_difference, 99.5)))
    difference_levels = np.linspace(-difference_limit, difference_limit, 17)
    difference_image = difference_axis.contourf(
        last_result.x_m,
        last_result.z_m,
        difference_masked,
        levels=difference_levels,
        cmap="coolwarm",
        extend="both",
    )
    assert isinstance(last_config.obstacle, RectangularObstacle2D)
    _draw_obstacle(difference_axis, last_config.obstacle)
    difference_axis.axhline(0.0, color="#111827", linewidth=1.3)
    difference_axis.set_xlim(0.0, last_config.length_m)
    difference_axis.set_ylim(0.0, last_config.height_m)
    difference_axis.set_xlabel("Downwind distance, x (m)")
    difference_axis.set_title("Obstacle effect at 10 s: obstacle minus open field", weight="bold")
    difference_axis.text(
        0.03,
        0.94,
        "red: local increase   blue: local decrease",
        transform=difference_axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 3},
    )
    difference_bar = fig.colorbar(
        difference_image,
        ax=difference_axis,
        orientation="horizontal",
        pad=0.18,
        fraction=0.065,
    )
    difference_bar.set_label("Concentration change (vol percentage points)")

    profile_axis = axes[1, 2]
    open_envelope = np.nanmax(open_field, axis=0)
    obstacle_envelope = np.nanmax(last_obstacle_field, axis=0)
    profile_axis.plot(
        open_result.x_m,
        open_envelope,
        color="#64748b",
        linewidth=2.0,
        linestyle="--",
        label="Open field",
    )
    profile_axis.plot(
        last_result.x_m,
        obstacle_envelope,
        color="#c026d3",
        linewidth=2.4,
        label="With obstacle",
    )
    profile_axis.axhline(4.0, color="#111827", linewidth=1.5, linestyle=":", label="4 vol%")
    profile_axis.axvspan(
        last_config.obstacle.x_min_m,
        last_config.obstacle.x_max_m,
        color="#263238",
        alpha=0.18,
        label="Obstacle",
    )
    profile_axis.set_xlim(0.0, last_config.length_m)
    profile_axis.set_ylim(0.0, max(36.0, float(np.nanmax([open_envelope, obstacle_envelope])) * 1.08))
    profile_axis.set_xlabel("Downwind distance, x (m)")
    profile_axis.set_ylabel("Vertical-maximum H$_2$ (vol%)")
    profile_axis.set_title("Hazard-envelope comparison at 10 s", weight="bold")
    profile_axis.grid(True, color="#cbd5e1", linewidth=0.6, alpha=0.7)
    profile_axis.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=8.5)

    assert last_image is not None
    concentration_bar_axis = fig.add_axes((0.16, 0.075, 0.68, 0.022))
    concentration_bar = fig.colorbar(
        last_image,
        cax=concentration_bar_axis,
        orientation="horizontal",
        ticks=levels,
    )
    concentration_bar.ax.set_xticklabels([f"{value:g}" for value in levels])
    concentration_bar.set_label("Hydrogen concentration (vol%, ideal-gas screening conversion)")

    obstacle_extent = _threshold_extent(last_obstacle_field, last_result.x_m)
    open_extent = _threshold_extent(open_field, open_result.x_m)
    total_cell_steps = sum(
        config.nx * config.nz * int(np.ceil(config.duration_s / config.time_step_s))
        for _, config, _, _ in obstacle_cases
    ) + open_config.nx * open_config.nz * int(
        np.ceil(open_config.duration_s / open_config.time_step_s)
    )
    per_million_cell_steps = solve_runtime / (total_cell_steps / 1_000_000.0)

    legend_items = [
        Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
               markeredgecolor="#111827", markersize=10, label="H$_2$ source"),
        Patch(facecolor="#263238", edgecolor="white", hatch="////", label="Solid obstacle"),
        Line2D([0], [0], color="#111827", linewidth=1.6, label="4 vol% isoline"),
    ]
    fig.legend(
        handles=legend_items,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.928),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        "DEGALI hypothetical obstacle case — transient concentration screening",
        y=0.992,
    )
    fig.text(
        0.5,
        0.958,
        "x–z section | cuboid surrogate: 2.5 m long × 3.0 m high at x = 9.0 m | "
        "source 0.01 kg s$^{-1}$ at z = 0.60 m",
        ha="center",
        va="center",
        fontsize=10,
    )
    fig.text(
        0.5,
        0.008,
        f"Reduced-order conservative semi-FV screen — not obstacle-resolved RANS/LES CFD; "
        f"wake recirculation is not resolved. Total solver time {solve_runtime:.2f} s "
        f"({total_cell_steps / 1_000_000.0:.2f} million cell-steps; "
        f"{per_million_cell_steps:.3f} s per million cell-steps).",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#334155",
    )
    concentration_bar_axis.set_position((0.16, 0.090, 0.68, 0.020))
    fig.subplots_adjust(left=0.065, right=0.975, top=0.875, bottom=0.165, hspace=0.43, wspace=0.24)

    png_path = output_dir / "degali_hypothetical_obstacle_cfd_style.png"
    svg_path = output_dir / "degali_hypothetical_obstacle_cfd_style.svg"
    pdf_path = output_dir / "degali_hypothetical_obstacle_cfd_style.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    wall_runtime = time.perf_counter() - wall_started
    metadata = {
        "classification": "hypothetical_reduced_order_obstacle_screen",
        "not_claimed": [
            "obstacle-resolved CFD",
            "RANS or LES",
            "validated wake recirculation prediction",
            "experimental agreement",
        ],
        "coordinate_plane": "x-z wind-aligned vertical section",
        "ambient_conversion": {
            "basis": "ideal gas",
            "temperature_K": TEMPERATURE_K,
            "pressure_Pa": PRESSURE_PA,
            "hydrogen_molar_mass_kg_mol": MOLAR_MASS_H2_KG_MOL,
        },
        "scenario": {
            "domain_m": [last_config.length_m, last_config.height_m],
            "cells": [last_config.nx, last_config.nz],
            "requested_time_step_s": last_config.time_step_s,
            "wind_speed_m_s": last_config.wind_speed_m_s,
            "diffusivity_m2_s": last_config.diffusivity_m2_s,
            "source_rate_kg_s": last_config.source_rate_kg_s,
            "source_height_m": last_config.source_height_m,
            "source_sigma_m": last_config.source_sigma_m,
            "obstacle": {
                "centre_x_m": last_config.obstacle.distance_m,
                "width_m": last_config.obstacle.width_m,
                "height_m": last_config.obstacle.height_m,
                "base_height_m": last_config.obstacle.base_height_m,
            },
            "snapshot_times_s": list(snapshot_times),
        },
        "runtime": {
            "obstacle_solve_s_by_duration": {
                f"{duration:g}": runtime
                for duration, _config, _result, runtime in obstacle_cases
            },
            "open_field_10_s_solve_s": open_runtime,
            "total_solver_s": solve_runtime,
            "total_wall_including_render_and_export_s": wall_runtime,
            "total_cell_steps": total_cell_steps,
            "seconds_per_million_cell_steps": per_million_cell_steps,
        },
        "diagnostics": {
            "applicability": last_result.diagnostics.applicability,
            "obstacle_contact": last_result.diagnostics.obstacle_contact,
            "diverted_mass_fraction": last_result.diagnostics.diverted_mass_fraction,
            "mass_injected_kg": last_result.diagnostics.mass_injected_kg,
            "mass_domain_kg": last_result.diagnostics.mass_domain_kg,
            "mass_outflow_kg": last_result.diagnostics.mass_outflow_kg,
            "final_mass_residual_kg": last_result.diagnostics.final_mass_residual_kg,
            "maximum_mass_residual_kg": last_result.diagnostics.maximum_mass_residual_kg,
            "maximum_concentration_vol_percent": float(np.nanmax(last_obstacle_field)),
            "four_vol_percent_extent_m": obstacle_extent,
            "open_field_four_vol_percent_extent_m": open_extent,
            "warnings": list(last_result.diagnostics.warnings),
        },
        "outputs": {
            "png": str(png_path),
            "svg": str(svg_path),
            "pdf": str(pdf_path),
        },
    }
    metadata_path = output_dir / "figure_metadata.json"
    metadata["outputs"]["metadata"] = str(metadata_path)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    metadata = build(args.output_dir.resolve())
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
