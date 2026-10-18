"""Compare legacy non-local routing with the Tier-1 local obstacle model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport import (
    FastObstacleReceptor,
    FastObstacleTransportConfig,
    solve_fast_obstacle_transport,
)
from degali.addons.semi_fv_obstacle import (
    RectangularObstacle2D,
    SemiFVConfig,
    SemiFVReceptor,
    SourceRateSchedule,
    solve_semi_fv_obstacle,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"
R_GAS_J_MOL_K = 8.314462618
TEMPERATURE_K = 293.15
PRESSURE_PA = 101_325.0
MOLAR_MASS_H2_KG_MOL = 2.01588e-3


def _vol_percent(concentration_kg_m3: np.ndarray) -> np.ndarray:
    molar_density = PRESSURE_PA / (R_GAS_J_MOL_K * TEMPERATURE_K)
    return np.clip(
        concentration_kg_m3 / MOLAR_MASS_H2_KG_MOL / molar_density * 100.0,
        0.0,
        100.0,
    )


def _configs():
    obstacle = RectangularObstacle2D(
        distance_m=9.0,
        width_m=2.5,
        height_m=3.0,
        label="ground-cuboid-section",
    )
    legacy = SemiFVConfig(
        length_m=24.0,
        height_m=8.0,
        nx=144,
        nz=64,
        duration_s=10.0,
        time_step_s=0.025,
        wind_speed_m_s=1.8,
        diffusivity_m2_s=0.08,
        source_rate_kg_s=0.01,
        source_height_m=0.60,
        source_sigma_m=0.35,
        obstacle=obstacle,
    )
    replacement = FastObstacleTransportConfig(
        length_m=24.0,
        height_m=8.0,
        nx=144,
        nz=64,
        duration_s=10.0,
        time_step_s=0.0025,
        wind_speed_m_s=1.8,
        diffusivity_m2_s=0.08,
        source_schedule=SourceRateSchedule(
            (0.0, 10.0),
            (0.01, 0.0),
            source_id="matched-finite-atmospheric-h2-per-metre-span",
        ),
        source_x_m=0.0,
        source_height_m=0.60,
        source_sigma_m=0.35,
        obstacle=obstacle,
    )
    return obstacle, legacy, replacement


def _draw_concentration(
    axis,
    x,
    z,
    field,
    *,
    title,
    obstacle,
    levels,
    cmap,
    norm,
):
    from matplotlib.patches import Rectangle

    plotted = np.ma.masked_less_equal(np.ma.masked_invalid(field), 0.0)
    image = axis.contourf(
        x, z, plotted, levels=levels, cmap=cmap, norm=norm, extend="max",
    )
    if float(np.nanmax(field)) >= 4.0:
        axis.contour(x, z, plotted, levels=[4.0], colors="white", linewidths=3.5)
        axis.contour(x, z, plotted, levels=[4.0], colors="#111827", linewidths=1.3)
    axis.add_patch(Rectangle(
        (obstacle.x_min_m, 0.0), obstacle.width_m, obstacle.height_m,
        facecolor="#263238", edgecolor="white", linewidth=1.3, hatch="////",
        zorder=8,
    ))
    axis.set_xlim(0.0, 24.0)
    axis.set_ylim(0.0, 8.0)
    axis.set_xlabel("Downwind x (m)")
    axis.set_ylabel("Height z (m)")
    axis.set_title(title, weight="bold")
    return image


def build(output_dir: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.patches import Rectangle

    output_dir.mkdir(parents=True, exist_ok=True)
    obstacle, legacy_config, replacement_config = _configs()
    legacy_receptors = (
        SemiFVReceptor("windward", 7.5, 0.60),
        SemiFVReceptor("roof", 9.0, 3.20),
        SemiFVReceptor("lee", 12.0, 0.60),
    )
    replacement_receptors = tuple(
        FastObstacleReceptor(item.label, item.x_m, item.z_m)
        for item in legacy_receptors
    )

    started = time.perf_counter()
    legacy_started = time.perf_counter()
    legacy = solve_semi_fv_obstacle(
        legacy_config, receptors=legacy_receptors,
    )
    legacy_runtime = time.perf_counter() - legacy_started
    replacement = solve_fast_obstacle_transport(
        replacement_config, receptors=replacement_receptors,
    )
    calculation_runtime = time.perf_counter() - started

    legacy_vol = _vol_percent(legacy.concentration_kg_m3)
    replacement_vol = replacement.concentration_vol_percent
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
    fig, axes = plt.subplots(2, 3, figsize=(16.4, 9.1))
    image = _draw_concentration(
        axes[0, 0], legacy.x_m, legacy.z_m, legacy_vol,
        title=f"Legacy routing | solve {legacy_runtime:.2f} s",
        obstacle=obstacle, levels=levels, cmap=cmap, norm=norm,
    )
    _draw_concentration(
        axes[0, 1], replacement.x_m, replacement.z_m, replacement_vol,
        title=(
            "Tier 1: local flux + divergence-free flow | solve "
            f"{replacement.diagnostics['total_runtime_s']:.2f} s"
        ),
        obstacle=obstacle, levels=levels, cmap=cmap, norm=norm,
    )
    speed = np.sqrt(
        replacement.velocity_u_m_s**2 + replacement.velocity_w_m_s**2
    )
    linewidth = 0.45 + 1.2 * speed / max(float(np.nanmax(speed)), 1.0e-12)
    axes[0, 1].streamplot(
        replacement.x_m,
        replacement.z_m,
        np.ma.masked_where(replacement.solid_mask, replacement.velocity_u_m_s),
        np.ma.masked_where(replacement.solid_mask, replacement.velocity_w_m_s),
        color="white",
        density=1.1,
        linewidth=linewidth,
        arrowsize=0.7,
        minlength=0.2,
        zorder=7,
    )

    difference = replacement_vol - legacy_vol
    finite = np.abs(difference[np.isfinite(difference)])
    limit = max(1.0, float(np.percentile(finite, 99.0)))
    difference_levels = np.linspace(-limit, limit, 17)
    difference_image = axes[0, 2].contourf(
        replacement.x_m,
        replacement.z_m,
        np.ma.masked_invalid(difference),
        levels=difference_levels,
        cmap="coolwarm",
        extend="both",
    )
    axes[0, 2].add_patch(Rectangle(
        (obstacle.x_min_m, 0.0), obstacle.width_m, obstacle.height_m,
        facecolor="#263238", edgecolor="white", linewidth=1.3, hatch="////",
        zorder=8,
    ))
    axes[0, 2].set_xlim(0.0, 24.0)
    axes[0, 2].set_ylim(0.0, 8.0)
    axes[0, 2].set_xlabel("Downwind x (m)")
    axes[0, 2].set_ylabel("Height z (m)")
    axes[0, 2].set_title("Tier 1 minus legacy (vol percentage points)", weight="bold")
    fig.colorbar(difference_image, ax=axes[0, 2], orientation="horizontal", pad=0.17).set_label(
        "Concentration change (vol percentage points)"
    )

    cp_axis = axes[1, 0]
    cp_image = cp_axis.contourf(
        replacement.x_m,
        replacement.z_m,
        np.ma.masked_invalid(replacement.pressure_coefficient),
        levels=np.linspace(-3.0, 1.0, 17),
        cmap="viridis",
        extend="min",
    )
    cp_axis.streamplot(
        replacement.x_m,
        replacement.z_m,
        np.ma.masked_where(replacement.solid_mask, replacement.velocity_u_m_s),
        np.ma.masked_where(replacement.solid_mask, replacement.velocity_w_m_s),
        color="white",
        density=1.1,
        linewidth=0.7,
        arrowsize=0.75,
        minlength=0.2,
    )
    cp_axis.add_patch(Rectangle(
        (obstacle.x_min_m, 0.0), obstacle.width_m, obstacle.height_m,
        facecolor="#263238", edgecolor="white", linewidth=1.3, hatch="////",
        zorder=8,
    ))
    cp_axis.set_xlim(0.0, 24.0)
    cp_axis.set_ylim(0.0, 8.0)
    cp_axis.set_xlabel("Downwind x (m)")
    cp_axis.set_ylabel("Height z (m)")
    cp_axis.set_title("Tier-1 pressure coefficient and streamlines", weight="bold")
    fig.colorbar(cp_image, ax=cp_axis, orientation="horizontal", pad=0.17).set_label(
        "Potential-flow pressure coefficient, $C_p$"
    )

    envelope_axis = axes[1, 1]
    envelope_axis.plot(
        legacy.x_m,
        np.nanmax(legacy_vol, axis=0),
        color="#ef4444",
        linewidth=2.2,
        label="Legacy routing",
    )
    envelope_axis.plot(
        replacement.x_m,
        np.nanmax(replacement_vol, axis=0),
        color="#2563eb",
        linewidth=2.2,
        label="Tier 1 local flux",
    )
    envelope_axis.axhline(4.0, color="#111827", linestyle=":", linewidth=1.5, label="4 vol%")
    envelope_axis.axvspan(obstacle.x_min_m, obstacle.x_max_m, color="#263238", alpha=0.16)
    envelope_axis.set_xlim(0.0, 24.0)
    envelope_axis.set_xlabel("Downwind x (m)")
    envelope_axis.set_ylabel("Vertical-maximum H$_2$ (vol%)")
    envelope_axis.set_title("Hazard envelope: routing spike removed", weight="bold")
    envelope_axis.grid(True, color="#cbd5e1", linewidth=0.6)
    envelope_axis.legend(fontsize=8.5)

    bar_axis = axes[1, 2]
    names = ("windward", "roof", "lee")
    legacy_traces = {item.receptor.label: item for item in legacy.receptor_traces}
    legacy_peaks = np.asarray([
        float(np.max(legacy_traces[name].concentration_kg_m3)) for name in names
    ])
    replacement_peaks = np.asarray([
        float(np.max(replacement.receptor_traces_kg_m3[name])) for name in names
    ])
    conversion = (
        R_GAS_J_MOL_K * TEMPERATURE_K
        / (PRESSURE_PA * MOLAR_MASS_H2_KG_MOL)
        * 100.0
    )
    legacy_peaks *= conversion
    replacement_peaks *= conversion
    positions = np.arange(len(names))
    width = 0.36
    bar_axis.bar(
        positions - width / 2,
        legacy_peaks,
        width,
        color="#ef4444",
        label="Legacy routing",
    )
    bar_axis.bar(
        positions + width / 2,
        replacement_peaks,
        width,
        color="#2563eb",
        label="Tier 1 local flux",
    )
    bar_axis.set_xticks(positions, ("Windward", "Roof", "Lee ground"))
    bar_axis.set_ylabel("Peak H$_2$ (vol%)")
    bar_axis.set_title("Matched receptor peaks", weight="bold")
    bar_axis.grid(True, axis="y", color="#cbd5e1", linewidth=0.6)
    bar_axis.legend(fontsize=8.5)

    color_axis = fig.add_axes((0.17, 0.072, 0.66, 0.020))
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal", ticks=levels)
    colorbar.ax.set_xticklabels([f"{value:g}" for value in levels])
    colorbar.set_label("Hydrogen concentration (vol%, ideal-gas screening conversion)")
    fig.suptitle("DEGALI obstacle-physics rebuild: matched legacy vs local-flux screen", y=0.985)
    fig.text(
        0.5,
        0.947,
        "Same 24 m × 8 m wind plane, 2.5 m × 3.0 m obstacle, 1.8 m/s wind and 0.10 kg/m injected H$_2$",
        ha="center",
    )
    fig.text(
        0.5,
        0.010,
        "Tier 1 corrects local conservation and windward impact, but remains conditional: "
        "potential flow does not resolve no-slip separation, recirculating wake, thermal buoyancy or 3-D lateral bypass.",
        ha="center",
        fontsize=8.7,
        color="#334155",
    )
    fig.subplots_adjust(left=0.065, right=0.975, top=0.895, bottom=0.165, hspace=0.45, wspace=0.24)

    png_path = output_dir / "obstacle_physics_legacy_vs_tier1.png"
    svg_path = output_dir / "obstacle_physics_legacy_vs_tier1.svg"
    pdf_path = output_dir / "obstacle_physics_legacy_vs_tier1.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    legacy_max_index = np.unravel_index(np.nanargmax(legacy_vol), legacy_vol.shape)
    replacement_max_index = np.unravel_index(
        np.nanargmax(replacement_vol), replacement_vol.shape,
    )
    metadata = {
        "classification": "model-rebuild-diagnostic-not-validation",
        "scenario": {
            "domain_m": [24.0, 8.0],
            "cells": [144, 64],
            "duration_s": 10.0,
            "wind_speed_m_s": 1.8,
            "injected_h2_mass_per_span_kg_m": 0.1,
            "obstacle": {
                "centre_x_m": obstacle.distance_m,
                "width_m": obstacle.width_m,
                "height_m": obstacle.height_m,
            },
        },
        "legacy": {
            "runtime_s": legacy_runtime,
            "maximum_vol_percent": float(np.nanmax(legacy_vol)),
            "maximum_location_m": [
                float(legacy.x_m[legacy_max_index[1]]),
                float(legacy.z_m[legacy_max_index[0]]),
            ],
            "maximum_mass_residual_kg_m": legacy.diagnostics.maximum_mass_residual_kg,
            "nonlocal_scalar_routing": True,
        },
        "tier_1": {
            "runtime_s": replacement.diagnostics["total_runtime_s"],
            "maximum_vol_percent": float(np.nanmax(replacement_vol)),
            "maximum_location_m": [
                float(replacement.x_m[replacement_max_index[1]]),
                float(replacement.z_m[replacement_max_index[0]]),
            ],
            "maximum_mass_residual_kg_m": (
                replacement.maximum_mass_residual_per_span_kg_m
            ),
            "maximum_absolute_divergence_s_1": replacement.diagnostics[
                "maximum_absolute_divergence_s_1"
            ],
            "nonlocal_scalar_routing": False,
            "warnings": replacement.diagnostics["warnings"],
        },
        "receptor_peak_vol_percent": {
            name: {
                "legacy": float(left),
                "tier_1": float(right),
            }
            for name, left, right in zip(names, legacy_peaks, replacement_peaks)
        },
        "combined_calculation_runtime_s": calculation_runtime,
        "outputs": {
            "png": str(png_path),
            "svg": str(svg_path),
            "pdf": str(pdf_path),
        },
    }
    metadata_path = output_dir / "comparison_metadata.json"
    metadata["outputs"]["metadata"] = str(metadata_path)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output_dir.resolve()), indent=2))


if __name__ == "__main__":
    main()
