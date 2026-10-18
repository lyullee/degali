"""Build a matched legacy/Tier-1/Tier-2 obstacle-physics comparison figure."""

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
from degali.addons.resolved_obstacle_flow_2d import (
    ResolvedObstacleFlow2DConfig,
    solve_resolved_obstacle_flow_2d,
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
R_GAS = 8.314462618
TEMPERATURE_K = 293.15
PRESSURE_PA = 101_325.0
MOLAR_MASS_H2 = 2.01588e-3


def _vol_percent(concentration: np.ndarray) -> np.ndarray:
    molar_density = PRESSURE_PA / (R_GAS * TEMPERATURE_K)
    return np.clip(concentration / MOLAR_MASS_H2 / molar_density * 100.0, 0.0, 100.0)


def _draw_field(axis, x, z, field, obstacle, title, levels, cmap, norm):
    from matplotlib.patches import Rectangle

    plotted = np.ma.masked_less_equal(np.ma.masked_invalid(field), 0.0)
    image = axis.contourf(x, z, plotted, levels=levels, cmap=cmap, norm=norm, extend="max")
    if float(np.nanmax(field)) >= 4.0:
        axis.contour(x, z, plotted, levels=[4.0], colors="white", linewidths=3.2)
        axis.contour(x, z, plotted, levels=[4.0], colors="#111827", linewidths=1.1)
    axis.add_patch(Rectangle(
        (obstacle.x_min_m, 0.0), obstacle.width_m, obstacle.height_m,
        facecolor="#263238", edgecolor="white", hatch="////", linewidth=1.2,
        zorder=8,
    ))
    axis.set(xlim=(0.0, 16.0), ylim=(0.0, 6.0), xlabel="Downwind x (m)", ylabel="Height z (m)")
    axis.set_title(title, weight="bold")
    return image


def build(output_dir: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.patches import Rectangle

    output_dir.mkdir(parents=True, exist_ok=True)
    obstacle = RectangularObstacle2D(6.0, 2.0, 2.5, label="ground-cuboid-section")
    schedule = SourceRateSchedule(
        (0.0, 8.0), (0.01, 0.0),
        source_id="matched-finite-atmospheric-h2-per-metre-span",
    )
    legacy_config = SemiFVConfig(
        length_m=16.0, height_m=6.0, nx=128, nz=48,
        duration_s=8.0, time_step_s=0.02, wind_speed_m_s=1.8,
        diffusivity_m2_s=0.06, source_rate_kg_s=0.01,
        source_height_m=0.6, source_sigma_m=0.30, obstacle=obstacle,
    )
    tier1_config = FastObstacleTransportConfig(
        length_m=16.0, height_m=6.0, nx=128, nz=48,
        duration_s=8.0, time_step_s=0.003, wind_speed_m_s=1.8,
        diffusivity_m2_s=0.06, source_schedule=schedule,
        source_x_m=0.6, source_height_m=0.6, source_sigma_m=0.30,
        obstacle=obstacle,
    )
    tier2_config = ResolvedObstacleFlow2DConfig(
        length_m=16.0, height_m=6.0, nx=128, nz=48,
        wind_speed_m_s=1.8, kinematic_viscosity_m2_s=0.25,
        flow_steps=3000, lattice_inlet_speed=0.03,
        duration_s=8.0, scalar_time_step_s=0.003,
        scalar_diffusivity_m2_s=0.06, source_schedule=schedule,
        source_x_m=0.6, source_height_m=0.6, source_sigma_m=0.30,
        obstacle=obstacle, minimum_flow_steps=1000,
    )
    legacy_receptors = (
        SemiFVReceptor("windward", 4.8, 0.75),
        SemiFVReceptor("roof", 6.0, 2.75),
        SemiFVReceptor("lee", 7.5, 0.75),
    )
    receptors = tuple(
        FastObstacleReceptor(item.label, item.x_m, item.z_m)
        for item in legacy_receptors
    )

    started = time.perf_counter()
    legacy_started = time.perf_counter()
    legacy = solve_semi_fv_obstacle(legacy_config, receptors=legacy_receptors)
    legacy_runtime = time.perf_counter() - legacy_started
    tier1 = solve_fast_obstacle_transport(tier1_config, receptors=receptors)
    tier2 = solve_resolved_obstacle_flow_2d(tier2_config, receptors=receptors)
    combined_runtime = time.perf_counter() - started

    fields = (
        _vol_percent(legacy.concentration_kg_m3),
        tier1.concentration_vol_percent,
        tier2.concentration_vol_percent,
    )
    levels = np.asarray((0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_under("#f8fafc")
    cmap.set_bad("#263238")
    norm = colors.LogNorm(vmin=levels[0], vmax=levels[-1])
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "figure.titlesize": 16, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.2))
    titles = (
        f"Legacy non-local routing | {legacy_runtime:.2f} s",
        f"Tier 1 local potential flow | {tier1.diagnostics['total_runtime_s']:.2f} s",
        f"Tier 2 viscous resolved flow | {tier2.diagnostics['total_runtime_s']:.2f} s",
    )
    models = (legacy, tier1, tier2)
    image = None
    for axis, model, field, title in zip(axes[0], models, fields, titles):
        image = _draw_field(
            axis, model.x_m, model.z_m, field, obstacle, title, levels, cmap, norm,
        )
    for axis, model in zip(axes[0, 1:], (tier1, tier2)):
        speed = np.sqrt(model.velocity_u_m_s**2 + model.velocity_w_m_s**2)
        axis.streamplot(
            model.x_m, model.z_m,
            np.ma.masked_where(model.solid_mask, model.velocity_u_m_s),
            np.ma.masked_where(model.solid_mask, model.velocity_w_m_s),
            color="white", density=1.0,
            linewidth=0.45 + speed / max(float(np.nanmax(speed)), 1.0e-12),
            arrowsize=0.65, minlength=0.2, zorder=7,
        )

    velocity_axis = axes[1, 0]
    velocity_limit = max(1.0, float(np.nanpercentile(np.abs(tier2.velocity_u_m_s), 99)))
    velocity_image = velocity_axis.contourf(
        tier2.x_m, tier2.z_m,
        np.ma.masked_where(tier2.solid_mask, tier2.velocity_u_m_s),
        levels=np.linspace(-velocity_limit, velocity_limit, 21), cmap="coolwarm", extend="both",
    )
    velocity_axis.streamplot(
        tier2.x_m, tier2.z_m,
        np.ma.masked_where(tier2.solid_mask, tier2.velocity_u_m_s),
        np.ma.masked_where(tier2.solid_mask, tier2.velocity_w_m_s),
        color="#111827", density=1.2, linewidth=0.65, arrowsize=0.7,
    )
    velocity_axis.add_patch(Rectangle(
        (obstacle.x_min_m, 0.0), obstacle.width_m, obstacle.height_m,
        facecolor="#263238", edgecolor="white", hatch="////", zorder=8,
    ))
    velocity_axis.set(xlim=(0, 16), ylim=(0, 6), xlabel="Downwind x (m)", ylabel="Height z (m)")
    velocity_axis.set_title("Tier 2 streamwise velocity: negative wake is resolved", weight="bold")
    fig.colorbar(velocity_image, ax=velocity_axis, orientation="horizontal", pad=0.17).set_label("u velocity (m/s)")

    envelope_axis = axes[1, 1]
    colors_by_model = ("#ef4444", "#2563eb", "#059669")
    labels = ("Legacy", "Tier 1", "Tier 2")
    for model, field, color, label in zip(models, fields, colors_by_model, labels):
        envelope_axis.plot(model.x_m, np.nanmax(field, axis=0), color=color, linewidth=2.1, label=label)
    envelope_axis.axhline(4.0, color="#111827", linestyle=":", linewidth=1.5, label="4 vol%")
    envelope_axis.axvspan(obstacle.x_min_m, obstacle.x_max_m, color="#263238", alpha=0.16)
    envelope_axis.set(xlim=(0, 16), xlabel="Downwind x (m)", ylabel="Vertical-maximum H$_2$ (vol%)")
    envelope_axis.set_title("Hazard envelope across fidelity tiers", weight="bold")
    envelope_axis.grid(True, color="#cbd5e1", linewidth=0.6)
    envelope_axis.legend(fontsize=8.5)

    peak_axis = axes[1, 2]
    names = ("windward", "roof", "lee")
    conversion = R_GAS * TEMPERATURE_K / (PRESSURE_PA * MOLAR_MASS_H2) * 100.0
    legacy_by_name = {trace.receptor.label: trace for trace in legacy.receptor_traces}
    peaks = np.asarray((
        [np.max(legacy_by_name[name].concentration_kg_m3) * conversion for name in names],
        [np.max(tier1.receptor_traces_kg_m3[name]) * conversion for name in names],
        [np.max(tier2.receptor_traces_kg_m3[name]) * conversion for name in names],
    ))
    positions = np.arange(len(names))
    width = 0.24
    for index, (values, color, label) in enumerate(zip(peaks, colors_by_model, labels)):
        peak_axis.bar(positions + (index - 1) * width, values, width, color=color, label=label)
    peak_axis.set_xticks(positions, ("Windward", "Roof", "Lee ground"))
    peak_axis.set_ylabel("Peak H$_2$ (vol%)")
    peak_axis.set_title("Matched receptor peaks", weight="bold")
    peak_axis.grid(True, axis="y", color="#cbd5e1", linewidth=0.6)
    peak_axis.legend(fontsize=8.5)

    assert image is not None
    color_axis = fig.add_axes((0.17, 0.071, 0.66, 0.020))
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal", ticks=levels)
    colorbar.ax.set_xticklabels([f"{value:g}" for value in levels])
    colorbar.set_label("Hydrogen concentration (vol%, ideal-gas conversion)")
    fig.suptitle("DEGALI obstacle physics: routing rule → local flow → resolved viscous wake", y=0.985, weight="bold")
    fig.text(
        0.5, 0.946,
        "Matched 16 m × 6 m plane, 2.0 m × 2.5 m obstacle, 1.8 m/s wind and 0.08 kg/m injected H$_2$",
        ha="center",
    )
    fig.text(
        0.5, 0.010,
        "Diagnostic model comparison, not experimental validation. Tier 2 is 2-D/isothermal with prescribed effective viscosity; quantitative facility use remains withheld.",
        ha="center", fontsize=8.7, color="#334155",
    )
    fig.subplots_adjust(left=0.065, right=0.975, top=0.895, bottom=0.165, hspace=0.45, wspace=0.24)

    stem = output_dir / "obstacle_fidelity_legacy_tier1_tier2"
    paths = {suffix: stem.with_suffix(f".{suffix}") for suffix in ("png", "svg", "pdf")}
    fig.savefig(paths["png"], dpi=220, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)

    model_metadata = {}
    for label, model, field in zip(("legacy", "tier_1", "tier_2"), models, fields):
        index = np.unravel_index(np.nanargmax(field), field.shape)
        diagnostics = model.diagnostics
        model_metadata[label] = {
            "runtime_s": (
                legacy_runtime if label == "legacy" else diagnostics["total_runtime_s"]
            ),
            "maximum_vol_percent": float(np.nanmax(field)),
            "maximum_location_m": [float(model.x_m[index[1]]), float(model.z_m[index[0]])],
            "nonlocal_scalar_routing": label == "legacy",
        }
    model_metadata["tier_2"].update({
        "minimum_lee_streamwise_velocity_m_s": tier2.diagnostics["minimum_lee_streamwise_velocity_m_s"],
        "recirculating_lee_cell_count": tier2.diagnostics["recirculating_lee_cell_count"],
        "maximum_absolute_divergence_s_1": tier2.diagnostics["maximum_absolute_divergence_s_1"],
        "maximum_mass_residual_per_span_kg_m": tier2.maximum_mass_residual_per_span_kg_m,
        "quantitative_prediction_allowed": False,
    })
    metadata = {
        "classification": "model-rebuild-diagnostic-not-validation",
        "scenario": {
            "domain_m": [16.0, 6.0], "cells": [128, 48], "duration_s": 8.0,
            "wind_speed_m_s": 1.8, "injected_h2_mass_per_span_kg_m": 0.08,
            "obstacle": {"centre_x_m": 6.0, "width_m": 2.0, "height_m": 2.5},
        },
        "models": model_metadata,
        "receptor_peak_vol_percent": {
            name: {label: float(peaks[j, i]) for j, label in enumerate(("legacy", "tier_1", "tier_2"))}
            for i, name in enumerate(names)
        },
        "combined_calculation_runtime_s": combined_runtime,
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    metadata_path = stem.with_name(stem.name + "_metadata.json")
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
