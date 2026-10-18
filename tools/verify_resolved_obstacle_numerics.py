"""Run canonical and three-grid numerical audits for Tier-2 obstacle flow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from degali.addons.fast_obstacle_transport import FastObstacleReceptor
from degali.addons.resolved_obstacle_flow_2d import ResolvedObstacleFlow2DConfig
from degali.addons.semi_fv_obstacle import RectangularObstacle2D, SourceRateSchedule
from degali.validation.obstacle_flow_canonical import (
    PoiseuilleChannelVerificationConfig,
    verify_poiseuille_channel,
)
from degali.validation.obstacle_flow_convergence import (
    resolved_obstacle_grid_convergence_record,
    run_resolved_obstacle_grid_convergence,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"


def _grid_configs() -> tuple[ResolvedObstacleFlow2DConfig, ...]:
    schedule = SourceRateSchedule(
        (0.0, 1.0, 6.0), (0.02, 0.0, 0.0),
        source_id="tier2-grid-convergence-finite-source",
    )
    obstacle = RectangularObstacle2D(4.0, 1.2, 2.0, label="grid-audit-block")
    configs = []
    for nx, nz in ((60, 20), (90, 30), (120, 40), (150, 50)):
        cell_size = 12.0 / nx
        flow_steps = int(round(4.5 / (0.03 * cell_size / 2.0)))
        configs.append(ResolvedObstacleFlow2DConfig(
            length_m=12.0,
            height_m=4.0,
            nx=nx,
            nz=nz,
            wind_speed_m_s=2.0,
            kinematic_viscosity_m2_s=0.3,
            flow_steps=flow_steps,
            lattice_inlet_speed=0.03,
            duration_s=6.0,
            scalar_time_step_s=0.02 * cell_size,
            scalar_diffusivity_m2_s=0.03,
            source_schedule=schedule,
            source_x_m=0.8,
            source_height_m=0.5,
            source_sigma_m=0.25,
            obstacle=obstacle,
            minimum_flow_steps=flow_steps,
        ))
    return tuple(configs)


def _receptors() -> tuple[FastObstacleReceptor, ...]:
    return (
        FastObstacleReceptor("windward", 3.2, 1.0),
        FastObstacleReceptor("roof", 4.0, 2.25),
        FastObstacleReceptor("lee-low", 5.2, 1.0),
        FastObstacleReceptor("lee-shear", 5.2, 2.25),
    )


def build(output_dir: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)
    channel = verify_poiseuille_channel(
        PoiseuilleChannelVerificationConfig(
            nx=32, nz=20, steps=5_000, convergence_check_interval=500,
        )
    )
    convergence = run_resolved_obstacle_grid_convergence(
        _grid_configs(), receptors=_receptors(), relative_tolerance=0.10,
    )
    record = resolved_obstacle_grid_convergence_record(convergence)
    record["canonical_poiseuille"] = {
        "schema": channel.schema,
        "passed": channel.passed,
        "relative_l2_error": channel.relative_l2_error,
        "relative_linf_error": channel.relative_linf_error,
        "relative_mass_flux_error": channel.relative_mass_flux_error,
        "maximum_streamwise_nonuniformity": channel.maximum_streamwise_nonuniformity,
        "final_convergence_residual": channel.final_convergence_residual,
        "maximum_relative_density_variation": channel.maximum_relative_density_variation,
        "runtime_s": channel.runtime_s,
        "diagnostics": channel.diagnostics,
    }
    record["overall_numerical_gate_passed"] = bool(
        channel.passed and convergence.passed
    )

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "figure.titlesize": 15,
        "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2))
    axes[0, 0].plot(
        channel.numerical_velocity_lattice,
        channel.wall_normal_coordinate_lattice,
        "o", color="#2563eb", markersize=4, label="D2Q9",
    )
    axes[0, 0].plot(
        channel.analytic_velocity_lattice,
        channel.wall_normal_coordinate_lattice,
        color="#111827", linewidth=2, label="Analytic",
    )
    axes[0, 0].set(
        xlabel="Streamwise velocity (lattice units)",
        ylabel="Wall-normal coordinate (lattice units)",
        title=(
            "Canonical Poiseuille profile\n"
            f"L2={channel.relative_l2_error:.3%}, "
            f"L∞={channel.relative_linf_error:.3%}"
        ),
    )
    axes[0, 0].grid(True, color="#cbd5e1", linewidth=0.6)
    axes[0, 0].legend()

    cell_sizes = np.asarray([case.cell_size_m for case in convergence.cases])
    colors = ("#2563eb", "#f97316", "#dc2626", "#059669")
    labels = tuple(metric.label for metric in convergence.cases[0].receptor_metrics)
    for color, label in zip(colors, labels):
        peaks = [
            next(item.peak_vol_percent for item in case.receptor_metrics if item.label == label)
            for case in convergence.cases
        ]
        axes[0, 1].plot(cell_sizes, peaks, "o-", color=color, label=label)
    axes[0, 1].invert_xaxis()
    axes[0, 1].set(
        xlabel="Cell size Δ (m), finer →",
        ylabel="Peak H$_2$ (vol%)",
        title="Receptor peak grid convergence",
    )
    axes[0, 1].grid(True, color="#cbd5e1", linewidth=0.6)
    axes[0, 1].legend(fontsize=8)

    for color, label in zip(colors, labels):
        doses = [
            next(item.dose_vol_percent_s for item in case.receptor_metrics if item.label == label)
            for case in convergence.cases
        ]
        axes[1, 0].plot(cell_sizes, doses, "o-", color=color, label=label)
    axes[1, 0].invert_xaxis()
    axes[1, 0].set(
        xlabel="Cell size Δ (m), finer →",
        ylabel="H$_2$ dose (vol%·s)",
        title="Receptor dose grid convergence",
    )
    axes[1, 0].grid(True, color="#cbd5e1", linewidth=0.6)
    axes[1, 0].legend(fontsize=8)

    area_axis = axes[1, 1]
    areas = [case.four_vol_percent_envelope_area_per_span_m2 for case in convergence.cases]
    runtimes = [case.runtime_s for case in convergence.cases]
    area_axis.plot(cell_sizes, areas, "o-", color="#7c3aed", linewidth=2, label="4 vol% area")
    area_axis.invert_xaxis()
    area_axis.set(
        xlabel="Cell size Δ (m), finer →",
        ylabel="4 vol% envelope area per span (m²)",
        title=(
            "Finest-pair gate: "
            + ("PASS" if convergence.passed else "WITHHELD")
        ),
    )
    runtime_axis = area_axis.twinx()
    runtime_axis.plot(cell_sizes, runtimes, "s--", color="#475569", label="Runtime")
    runtime_axis.set_ylabel("Runtime (s)")
    area_axis.grid(True, color="#cbd5e1", linewidth=0.6)
    handles_left, labels_left = area_axis.get_legend_handles_labels()
    handles_right, labels_right = runtime_axis.get_legend_handles_labels()
    area_axis.legend(handles_left + handles_right, labels_left + labels_right, loc="best")
    if convergence.gate_reasons:
        area_axis.text(
            0.03, 0.04,
            "\n".join(convergence.gate_reasons[:4]),
            transform=area_axis.transAxes, fontsize=8, va="bottom",
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "#fff7ed", "alpha": 0.9},
        )

    fig.suptitle("DEGALI Tier-2 numerical verification and grid-convergence gate", y=0.985, weight="bold")
    fig.text(
        0.5, 0.018,
        "Numerical verification only. Passing does not establish obstacle-dispersion accuracy; thermal, buoyancy, turbulence, 3-D and external-data gates remain separate.",
        ha="center", fontsize=8.7, color="#334155",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.955))
    stem = output_dir / "tier2_numerical_verification"
    paths = {suffix: stem.with_suffix(f".{suffix}") for suffix in ("png", "svg", "pdf")}
    fig.savefig(paths["png"], dpi=220, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)

    record["outputs"] = {key: str(value) for key, value in paths.items()}
    json_path = stem.with_suffix(".json")
    record["outputs"]["json"] = str(json_path)
    json_path.write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output_dir.resolve()), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
