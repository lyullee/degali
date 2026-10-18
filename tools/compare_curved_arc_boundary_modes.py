"""Paired small-domain sensitivity run for curved-arc boundary treatments.

This is a numerical implementation comparison, not an independent validation
or a production-scale prediction. The two cases share mesh, flow/scalar time
steps, source, obstacle centreline and receptors; the baseline uses an
explicitly declared finite-thickness voxel wall while the alternate case uses
the zero-thickness surface mode.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import CircularArcWall


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "curved-arc-boundary-mode-comparison-2026-10-12-v1.json"


def _case(representation: str):
    wall = CircularArcWall(
        center_m=(2.0, 1.2), radius_m=0.5,
        thickness_m=0.2 if representation == "finite_thickness_voxel" else 0.0,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=-90.0,
        arc_end_deg=90.0, label="semicircular-fence",
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(4.0, 2.4, 2.0), cells=(40, 24, 20),
        wind_speed_m_s=2.0, kinematic_viscosity_m2_s=0.5,
        flow_steps=500, minimum_flow_steps=500, lattice_inlet_speed=0.03,
        duration_s=1.0, scalar_time_step_s=0.0015,
        scalar_diffusivity_m2_s=0.03,
        source_schedule=SourceRateSchedule(
            (0.0, 0.6, 1.0), (0.02, 0.0, 0.0),
            source_id="paired-curved-arc-boundary-mode-source",
        ),
        source_position_m=(0.5, 1.2, 0.4), source_sigma_m=0.15,
        obstacles=(wall,), circular_arc_representation=representation,
        smagorinsky_constant=0.10,
    )
    receptors = (
        FastObstacleReceptor3D("windward", (1.2, 1.2, 0.5)),
        FastObstacleReceptor3D("roof", (2.0, 1.2, 1.2)),
        FastObstacleReceptor3D("lee", (2.7, 1.2, 0.5)),
        FastObstacleReceptor3D("side", (2.6, 0.6, 0.5)),
    )
    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(config, receptors=receptors)
    elapsed = time.perf_counter() - started
    dx = config.domain_m[0] / config.cells[0]
    maximum = result.maximum_concentration_vol_percent
    return result, {
        "representation": representation,
        "grid_cells_xyz": list(config.cells),
        "cell_size_m": dx,
        "sample_count": int(result.time_s.size),
        "simulation_time_s": float(result.time_s[-1]),
        "runtime_wall_s": elapsed,
        "maximum_ever_h2_vol_percent": float(np.nanmax(maximum)),
        "maximum_ever_4vol_percent_volume_m3": float(
            np.count_nonzero(maximum >= 4.0) * dx**3
        ),
        "endpoint_4vol_percent_volume_m3": float(
            np.count_nonzero(result.concentration_vol_percent >= 4.0) * dx**3
        ),
        "maximum_species_mass_residual_kg": result.maximum_mass_residual_kg,
        "maximum_energy_residual_j": result.maximum_energy_residual_j,
        "flow_steps": result.diagnostics["flow_steps_completed"],
        "d3q19_zero_thickness_arc_link_count": result.diagnostics[
            "d3q19_zero_thickness_arc_link_count"
        ],
        "d3q19_zero_thickness_arc_link_rule": result.diagnostics[
            "d3q19_zero_thickness_arc_link_rule"
        ],
        "receptor_peak_vol_percent": {
            name: float(np.max(values))
            for name, values in result.receptor_traces_vol_percent.items()
        },
        "receptor_dose_vol_percent_s": {
            name: float(np.trapezoid(values, result.time_s))
            for name, values in result.receptor_traces_vol_percent.items()
        },
    }


def build(output_path: Path) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    baseline, baseline_record = _case("finite_thickness_voxel")
    surface, surface_record = _case("zero_thickness_surface")
    if not np.array_equal(baseline.time_s, surface.time_s):
        raise RuntimeError("paired runs produced different time axes")
    baseline_endpoint = baseline.concentration_vol_percent
    surface_endpoint = surface.concentration_vol_percent
    baseline_maximum = baseline.maximum_concentration_vol_percent
    surface_maximum = surface.maximum_concentration_vol_percent
    common_fluid = ~(baseline.solid_mask | surface.solid_mask)
    delta_endpoint = (
        surface_endpoint[common_fluid] - baseline_endpoint[common_fluid]
    )
    delta_maximum = (
        surface_maximum[common_fluid] - baseline_maximum[common_fluid]
    )
    shared_setup = {
        "domain_m": [4.0, 2.4, 2.0],
        "grid_cells_xyz": [40, 24, 20],
        "cell_size_m": 0.1,
        "wind_speed_m_s": 2.0,
        "kinematic_viscosity_m2_s": 0.5,
        "lattice_inlet_speed": 0.03,
        "d3q19_flow_steps": 500,
        "duration_s": 1.0,
        "scalar_time_step_s": 0.0015,
        "scalar_diffusivity_m2_s": 0.03,
        "subgrid_model": "smagorinsky",
        "smagorinsky_constant": 0.1,
        "source_position_m": [0.5, 1.2, 0.4],
        "source_sigma_m": 0.15,
        "source_schedule_time_s": [0.0, 0.6, 1.0],
        "source_schedule_rate_kg_s": [0.02, 0.0, 0.0],
        "source_schedule_rate_operator": "piecewise_constant",
        "source_total_mass_kg": 0.012,
        "obstacle_center_m": [2.0, 1.2],
        "obstacle_radius_m": 0.5,
        "obstacle_arc_deg": [-90.0, 90.0],
        "wall_height_m": 1.0,
        "voxel_mode_radial_thickness_m": 0.2,
        "surface_mode_radial_thickness_m": 0.0,
        "receptors_m": {
            "windward": [1.2, 1.2, 0.5],
            "roof": [2.0, 1.2, 1.2],
            "lee": [2.7, 1.2, 0.5],
            "side": [2.6, 0.6, 0.5],
        },
        "common_fluid_comparison_cell_count": int(np.count_nonzero(common_fluid)),
        "excluded_solid_or_voxelized_cells": int(
            common_fluid.size - np.count_nonzero(common_fluid)
        ),
    }
    canonical_setup = json.dumps(
        shared_setup, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    record: dict[str, object] = {
        "schema": "degali.curved-arc-boundary-mode-comparison.v1",
        "status": "numerical_sensitivity_only",
        "scope_warning": (
            "The voxel and zero-thickness modes represent different wall "
            "geometries. This paired synthetic case is not independent "
            "curved-wall or hydrogen-dispersion validation."
        ),
        "shared_setup": shared_setup,
        "shared_setup_sha256": hashlib.sha256(canonical_setup).hexdigest(),
        "runner_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": [baseline_record, surface_record],
        "paired_field_difference": {
            "endpoint_h2_vol_percent_mae": float(np.mean(np.abs(delta_endpoint))),
            "endpoint_h2_vol_percent_rmse": float(np.sqrt(np.mean(delta_endpoint**2))),
            "endpoint_h2_vol_percent_max_abs": float(np.max(np.abs(delta_endpoint))),
            "maximum_ever_h2_vol_percent_mae": float(np.mean(np.abs(delta_maximum))),
            "maximum_ever_h2_vol_percent_rmse": float(np.sqrt(np.mean(delta_maximum**2))),
            "maximum_ever_h2_vol_percent_max_abs": float(np.max(np.abs(delta_maximum))),
            "endpoint_four_percent_volume_delta_m3": (
                surface_record["endpoint_4vol_percent_volume_m3"]
                - baseline_record["endpoint_4vol_percent_volume_m3"]
            ),
            "maximum_ever_four_percent_volume_delta_m3": (
                surface_record["maximum_ever_4vol_percent_volume_m3"]
                - baseline_record["maximum_ever_4vol_percent_volume_m3"]
            ),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    x = baseline.x_m
    y = baseline.y_m
    z = baseline.z_m
    z_index = int(np.argmin(np.abs(z - 0.5)))
    y_index = int(np.argmin(np.abs(y - 1.2)))
    wake_x_index = int(np.argmin(np.abs(x - 2.7)))
    theta = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 300)
    figure, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    image = None
    for row, (case_result, label) in enumerate((
        (baseline, "Finite-thickness voxel"),
        (surface, "Zero-thickness link surface"),
    )):
        concentration = case_result.concentration_vol_percent
        planes = (
            (x, y, concentration[z_index], "Plan at z≈0.5 m", "x (m)", "y (m)"),
            (x, z, concentration[:, y_index, :], "Centre vertical plane", "x (m)", "z (m)"),
            (y, z, concentration[:, :, wake_x_index], "Wake cross-plane at x≈2.7 m", "y (m)", "z (m)"),
        )
        for column, (horizontal, vertical, values, title, xlabel, ylabel) in enumerate(planes):
            axis = axes[row, column]
            visible = np.ma.masked_where(~np.isfinite(values) | (values < 0.05), values)
            image = axis.pcolormesh(
                horizontal, vertical, visible,
                norm=LogNorm(vmin=0.05, vmax=50.0), cmap="turbo", shading="auto",
            )
            axis.set_title(title)
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)
            axis.set_aspect("equal", adjustable="box")
            if column == 0:
                arc_x = 2.0 + 0.5 * np.cos(theta)
                arc_y = 1.2 + 0.5 * np.sin(theta)
                axis.plot(arc_x, arc_y, color="white", linewidth=1.2, linestyle="--")
                axis.plot(0.5, 1.2, marker="*", color="black", markersize=8)
            elif column == 1:
                axis.axvline(2.0, color="white", linewidth=0.8, linestyle=":")
                axis.axhline(1.0, color="white", linewidth=0.8, linestyle=":")
            else:
                axis.axhline(1.0, color="white", linewidth=0.8, linestyle=":")
        axes[row, 0].text(
            0.02, 0.98, label, transform=axes[row, 0].transAxes,
            ha="left", va="top", color="white", weight="bold",
            bbox={"facecolor": "black", "alpha": 0.55, "pad": 3},
        )
    figure.colorbar(image, ax=axes, label="H₂ concentration (vol%)", shrink=0.88)
    figure.suptitle(
        "Curved-arc boundary sensitivity at t=1.0 s\n"
        "Common scale; values below 0.05 vol% masked",
        weight="bold",
    )
    figure_path = output_path.with_suffix(".png")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    record["figure_path"] = str(figure_path.resolve())
    record["figure_sha256"] = hashlib.sha256(figure_path.read_bytes()).hexdigest()
    record["output_path"] = str(output_path.resolve())
    output_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output.resolve()), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
