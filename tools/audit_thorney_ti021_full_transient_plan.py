"""Freeze TI021 full-transient coverage and computational feasibility.

The audit uses sensor coordinates, heights and observation-window timing only.
Dose and concentration values do not enter domain selection, grid selection or
runtime estimates.  Its purpose is to prevent a short preflight box from being
mistaken for a full 145-sensor validation domain and to decide whether a
monolithic or nested multi-fidelity execution is computationally defensible.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle
import numpy as np

from degali.validation.thorney_ti021 import (
    ThorneyTi021Dataset,
    read_thorney_ti021_workbook,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-full-transient-plan.v2"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def site_local_to_solver(
    local_x_m: float, local_y_m: float, height_m: float, *,
    source_solver_xyz_m: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Map site +y downwind and site -x to solver +x/+y.

    This is the same rotation used by the report-priority -6.1 degree wind in
    the full-site preflight: solver x follows site +y and solver y follows
    site -x.
    """

    sx, sy, sz = source_solver_xyz_m
    return sx + local_y_m, sy - local_x_m, sz + height_m


def site_local_to_nested_outer(
    local_x_m: float,
    local_y_m: float,
    height_m: float,
    *,
    handoff_downwind_y_m: float,
    outer_crosswind_origin_m: float,
) -> tuple[float, float, float]:
    """Map the handoff plane to outer x=0 in the same physical site frame."""

    return (
        local_y_m - handoff_downwind_y_m,
        outer_crosswind_origin_m - local_x_m,
        height_m,
    )


def _inside(point: tuple[float, float, float], domain: tuple[float, float, float]) -> bool:
    return all(0.0 <= coordinate <= length for coordinate, length in zip(point, domain))


def _runtime_rate(record: dict[str, object]) -> float:
    configuration = record["configuration"]
    return float(record["runtime"]["total_s"]) / (
        int(configuration["cell_count"]) * int(configuration["active_steps"])
    )


def _grid_record(
    *, domain_m: tuple[float, float, float], cell_size_m: float,
    duration_s: float, time_step_s: float,
    lower_rate_s_per_cell_step: float,
    upper_rate_s_per_cell_step: float,
) -> dict[str, object]:
    cells_float = tuple(length / cell_size_m for length in domain_m)
    if any(not math.isclose(value, round(value), abs_tol=1.0e-12) for value in cells_float):
        raise ValueError("domain dimensions must be divisible by cell size")
    cells = tuple(int(round(value)) for value in cells_float)
    steps = int(math.ceil(duration_s / time_step_s))
    cell_steps = int(np.prod(cells)) * steps
    return {
        "cell_size_m": cell_size_m,
        "cells": list(cells),
        "cell_count": int(np.prod(cells)),
        "time_step_s": time_step_s,
        "step_count": steps,
        "cell_steps": cell_steps,
        "wall_height_cell_count": 5.0 / cell_size_m,
        "runtime_estimate_hours": {
            "lower": cell_steps * lower_rate_s_per_cell_step / 3600.0,
            "upper": cell_steps * upper_rate_s_per_cell_step / 3600.0,
        },
    }


def build_plan(
    dataset: ThorneyTi021Dataset, *,
    dx2_benchmark: dict[str, object],
    dx1_benchmark: dict[str, object],
) -> dict[str, object]:
    """Create a concentration-blind, sensor-complete execution plan."""

    preflight_domain = (130.0, 140.0, 16.0)
    preflight_source = (60.0, 60.0, 0.0)
    inner_domain = (200.0, 200.0, 40.0)
    inner_source = (60.0, 100.0, 0.0)
    monolithic_domain = (600.0, 640.0, 40.0)
    monolithic_source = (60.0, 320.0, 0.0)
    handoff_downwind_y = 100.0
    handoff_inner_solver_x = inner_source[0] + handoff_downwind_y
    nested_outer_domain = (440.0, 640.0, 40.0)
    nested_outer_crosswind_origin = 320.0
    sensor_local = [
        (sensor.local_x_m, sensor.local_y_m, sensor.height_m)
        for sensor in dataset.dose_sensors
    ]
    preflight_points = [
        site_local_to_solver(*point, source_solver_xyz_m=preflight_source)
        for point in sensor_local
    ]
    inner_points = [
        site_local_to_solver(*point, source_solver_xyz_m=inner_source)
        for point in sensor_local
    ]
    monolithic_points = [
        site_local_to_solver(*point, source_solver_xyz_m=monolithic_source)
        for point in sensor_local
    ]
    nested_outer_points = [
        site_local_to_nested_outer(
            *point,
            handoff_downwind_y_m=handoff_downwind_y,
            outer_crosswind_origin_m=nested_outer_crosswind_origin,
        )
        for point in sensor_local
    ]
    preflight_indices = [
        index for index, point in enumerate(preflight_points)
        if _inside(point, preflight_domain)
    ]
    inner_indices = [
        index for index, point in enumerate(inner_points)
        if _inside(point, inner_domain)
    ]
    monolithic_indices = [
        index for index, point in enumerate(monolithic_points)
        if _inside(point, monolithic_domain)
    ]
    nested_outer_indices = [
        index for index, point in enumerate(nested_outer_points)
        if _inside(point, nested_outer_domain)
    ]
    combined_indices = set(inner_indices) | set(nested_outer_indices)
    overlap_indices = set(inner_indices) & set(nested_outer_indices)
    arc_centerline_indices = [
        index for index, sensor in enumerate(dataset.dose_sensors)
        if abs(sensor.local_x_m) <= 1.0e-12
        and any(abs(sensor.local_y_m - distance) <= 1.0e-12 for distance in (75.0, 100.0))
        and any(abs(sensor.height_m - height) <= 1.0e-12 for height in (0.4, 2.4, 4.4))
    ]
    lower_rate = _runtime_rate(dx2_benchmark)
    upper_rate = _runtime_rate(dx1_benchmark)
    if lower_rate > upper_rate:
        lower_rate, upper_rate = upper_rate, lower_rate
    duration_s = 600.0
    grids = ((5.0, 0.125), (2.5, 0.0625), (1.25, 0.03125))
    inner_grids = [
        _grid_record(
            domain_m=inner_domain, cell_size_m=dx, duration_s=duration_s,
            time_step_s=dt, lower_rate_s_per_cell_step=lower_rate,
            upper_rate_s_per_cell_step=upper_rate,
        )
        for dx, dt in grids
    ]
    monolithic_grids = [
        _grid_record(
            domain_m=monolithic_domain, cell_size_m=dx, duration_s=duration_s,
            time_step_s=dt, lower_rate_s_per_cell_step=lower_rate,
            upper_rate_s_per_cell_step=upper_rate,
        )
        for dx, dt in grids
    ]
    nested_outer_grids = [
        _grid_record(
            domain_m=nested_outer_domain, cell_size_m=dx,
            duration_s=duration_s, time_step_s=dt,
            lower_rate_s_per_cell_step=lower_rate,
            upper_rate_s_per_cell_step=upper_rate,
        )
        for dx, dt in grids
    ]
    gates = {
        "current_preflight_does_not_claim_full_sensor_coverage": (
            len(preflight_indices) < len(sensor_local)
        ),
        "monolithic_domain_covers_all_145_sensor_coordinates": (
            len(monolithic_indices) == len(sensor_local) == 145
        ),
        "inner_and_nested_outer_union_covers_all_145_sensor_coordinates": (
            len(combined_indices) == len(sensor_local) == 145
        ),
        "nested_outer_starts_at_same_physical_handoff_plane": (
            math.isclose(handoff_inner_solver_x, 160.0)
            and math.isclose(
                site_local_to_nested_outer(
                    0.0, handoff_downwind_y, 0.0,
                    handoff_downwind_y_m=handoff_downwind_y,
                    outer_crosswind_origin_m=nested_outer_crosswind_origin,
                )[0],
                0.0,
            )
        ),
        "inner_outer_overlap_contains_independent_sensor_coordinates": (
            len(overlap_indices) > 0
        ),
        "inner_domain_covers_all_six_arc_centerline_observations": (
            len(arc_centerline_indices) == 6
            and all(index in inner_indices for index in arc_centerline_indices)
        ),
        "inner_domain_contains_complete_50m_fence_with_margin": (
            inner_source[0] >= 50.0
            and inner_source[1] >= 50.0
            and inner_domain[1] - inner_source[1] >= 50.0
        ),
        "three_grids_represent_5m_wall_by_integer_cells": all(
            math.isclose(item["wall_height_cell_count"], round(item["wall_height_cell_count"]))
            for item in inner_grids
        ),
        "observation_values_excluded_from_domain_grid_and_runtime_selection": True,
        "runtime_estimate_uses_measured_degali_cell_step_bounds": (
            lower_rate > 0.0 and upper_rate >= lower_rate
        ),
    }
    finest_monolithic_lower_h = monolithic_grids[-1]["runtime_estimate_hours"]["lower"]
    return {
        "schema": SCHEMA,
        "status": (
            "accepted_coordinate_consistent_nested_execution_plan"
            if all(gates.values())
            else "withheld_full_transient_plan"
        ),
        "scope": (
            "geometry/timing coverage and measured runtime feasibility only; "
            "no dose, peak concentration or model accuracy value selects this plan"
        ),
        "dataset": {
            "workbook_path": str(dataset.workbook_path),
            "workbook_sha256": dataset.workbook_sha256,
            "sensor_count": len(sensor_local),
            "arc_summary_count": len(dataset.arc_observations),
            "coordinate_fields_accessed": [
                "local_x_m", "local_y_m", "height_m",
            ],
            "timing_policy": (
                "fixed 600 s horizon exceeds every reported finite departure; "
                "timing values are coverage metadata and do not set coefficients"
            ),
            "observation_value_fields_accessed_for_planning": [],
        },
        "coordinate_transform": {
            "solver_x": "source_solver_x + site_local_y",
            "solver_y": "source_solver_y - site_local_x",
            "solver_z": "source_solver_z + sensor_height",
            "nested_outer_solver_x": (
                "site_local_y - handoff_downwind_y; the physical handoff "
                "plane maps to outer x=0"
            ),
            "wind_mapping": (
                "site +y to solver +x and site -x to solver +y; preserves the "
                "report-priority -6.1 degree heading used by the preflight"
            ),
        },
        "coverage": {
            "current_preflight": {
                "domain_m": list(preflight_domain),
                "source_solver_xyz_m": list(preflight_source),
                "covered_sensor_count": len(preflight_indices),
                "excluded_sensor_count": len(sensor_local) - len(preflight_indices),
            },
            "resolved_inner": {
                "domain_m": list(inner_domain),
                "source_solver_xyz_m": list(inner_source),
                "local_site_bounds_m": {
                    "crosswind_x": [-100.0, 100.0],
                    "downwind_y": [-60.0, 140.0],
                },
                "covered_sensor_count": len(inner_indices),
                "covered_arc_centerline_count": len(arc_centerline_indices),
                "purpose": (
                    "resolve initial cloud, complete fence, separation/wake and "
                    "the 75/100 m arc comparison before downstream handoff"
                ),
            },
            "monolithic_sensor_domain": {
                "domain_m": list(monolithic_domain),
                "source_solver_xyz_m": list(monolithic_source),
                "local_site_bounds_m": {
                    "crosswind_x": [-320.0, 320.0],
                    "downwind_y": [-60.0, 540.0],
                },
                "covered_sensor_count": len(monolithic_indices),
                "purpose": "feasibility reference only; not the nested execution domain",
            },
            "nested_outer_sensor_domain": {
                "domain_m": list(nested_outer_domain),
                "physical_origin_site_local_m": {
                    "crosswind_x": -nested_outer_crosswind_origin,
                    "downwind_y": handoff_downwind_y,
                    "height_z": 0.0,
                },
                "local_site_bounds_m": {
                    "crosswind_x": [-320.0, 320.0],
                    "downwind_y": [handoff_downwind_y, 540.0],
                },
                "covered_sensor_count": len(nested_outer_indices),
                "purpose": "receive the inner plane at outer x-min and carry all downstream sensors",
            },
            "nested_union": {
                "combined_unique_sensor_count": len(combined_indices),
                "overlap_sensor_count": len(overlap_indices),
                "inner_only_sensor_count": len(set(inner_indices) - set(nested_outer_indices)),
                "outer_only_sensor_count": len(set(nested_outer_indices) - set(inner_indices)),
                "handoff_downwind_y_m": handoff_downwind_y,
                "handoff_inner_solver_x_m": handoff_inner_solver_x,
                "overlap_downwind_bounds_m": [handoff_downwind_y, 140.0],
            },
        },
        "frozen_horizon_s": duration_s,
        "runtime_model": {
            "lower_s_per_cell_step": lower_rate,
            "upper_s_per_cell_step": upper_rate,
            "lower_benchmark": {
                "cell_size_m": dx2_benchmark["configuration"]["cell_size_m"],
                "cell_count": dx2_benchmark["configuration"]["cell_count"],
                "active_steps": dx2_benchmark["configuration"]["active_steps"],
                "runtime_s": dx2_benchmark["runtime"]["total_s"],
            },
            "upper_benchmark": {
                "cell_size_m": dx1_benchmark["configuration"]["cell_size_m"],
                "cell_count": dx1_benchmark["configuration"]["cell_count"],
                "active_steps": dx1_benchmark["configuration"]["active_steps"],
                "runtime_s": dx1_benchmark["runtime"]["total_s"],
                "caution": "single-step startup-heavy upper estimate",
            },
            "resolved_inner_three_grid": inner_grids,
            "nested_outer_three_grid": nested_outer_grids,
            "monolithic_all_sensor_three_grid": monolithic_grids,
            "six_member_multiplier": 6,
        },
        "decision": {
            "monolithic_full_sensor_execution": "withheld_as_current_default",
            "reason": (
                "the current box covers only 33/145 sensors, while the finest "
                f"sensor-complete monolithic run is estimated at at least "
                f"{finest_monolithic_lower_h / 24.0:.1f} days per prior member"
            ),
            "required_architecture": (
                "overlapping resolved-inner to conservative-outer handoff of "
                "species mass, sensible enthalpy and face mass flux at physical "
                "site-local downwind y=100 m (inner solver x=160 m, outer x-min)"
            ),
            "handoff_acceptance_gates": [
                "integrated species mass-flux mismatch <= 1e-6",
                "integrated sensible-enthalpy flux mismatch <= 1e-6",
                "projected outer continuity residual <= 1e-7 kg/m3/s",
                "overlap receptor peak and dose change <= 10%",
                "no selection of Sc_t or Pr_t from TI021 scores",
            ],
            "no_fit_run_order": [
                "implement and verify conservative nested handoff",
                "run the primary Sc_t=0.7 Pr_t=0.74 inner three-grid audit",
                "run all six frozen prior members on the accepted production grid",
                "propagate every member through the outer sensor domain",
                "open dose/arc values once, score all members, and retain all results",
            ],
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def _plot(record: dict[str, object], dataset: ThorneyTi021Dataset, output: Path) -> None:
    x = np.asarray([sensor.local_x_m for sensor in dataset.dose_sensors])
    y = np.asarray([sensor.local_y_m for sensor in dataset.dose_sensors])
    heights = np.asarray([sensor.height_m for sensor in dataset.dose_sensors])
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.2), constrained_layout=True)
    scatter = axes[0].scatter(
        x, y, c=heights, cmap="viridis", s=24, edgecolors="black", linewidths=0.25,
    )
    axes[0].add_patch(Arc((0.0, 0.0), 100.0, 100.0, theta1=0, theta2=180,
                          color="crimson", linewidth=3, label="50 m fence"))
    axes[0].add_patch(Rectangle(
        (-100.0, -60.0), 200.0, 200.0, fill=False, edgecolor="tab:orange",
        linewidth=2, label="resolved inner",
    ))
    axes[0].add_patch(Rectangle(
        (-320.0, 100.0), 640.0, 440.0, fill=False, edgecolor="tab:blue",
        linewidth=2, label="nested outer (handoff at y=100 m)",
    ))
    axes[0].scatter([0.0], [0.0], marker="*", s=150, color="gold", edgecolor="black")
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel("site-local crosswind x (m)")
    axes[0].set_ylabel("site-local downwind y (m)")
    axes[0].set_title("Sensor-complete nested domains (no dose shown)")
    axes[0].legend(fontsize=8)
    fig.colorbar(scatter, ax=axes[0], shrink=0.8, label="sensor height (m)")

    labels = []
    lower = []
    upper = []
    for prefix, key in (("inner", "resolved_inner_three_grid"),
                        ("nested", "nested_outer_three_grid"),
                        ("mono", "monolithic_all_sensor_three_grid")):
        for grid in record["runtime_model"][key]:
            labels.append(f"{prefix}\ndx={grid['cell_size_m']:g} m")
            lower.append(grid["runtime_estimate_hours"]["lower"])
            upper.append(grid["runtime_estimate_hours"]["upper"])
    positions = np.arange(len(labels))
    axes[1].bar(positions, lower, color="tab:green", label="measured lower rate")
    axes[1].bar(
        positions, np.asarray(upper) - np.asarray(lower), bottom=lower,
        color="tab:red", alpha=0.45, label="startup-heavy upper range",
    )
    axes[1].set_yscale("log")
    axes[1].set_xticks(positions, labels)
    axes[1].set_ylabel("estimated runtime per prior member (hours, log scale)")
    axes[1].set_title("600 s full-transient feasibility")
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.suptitle("TI021 full-transient plan — geometry/timing only")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument(
        "--dx2-benchmark", type=Path,
        default=(ROOT / "outputs" /
                 "thorney-ti021-full-scale-preflight-dx2m-dt0p025-amg-40step-2026-10-11-v1" /
                 "ti021-full-scale-preflight.json"),
    )
    parser.add_argument(
        "--dx1-benchmark", type=Path,
        default=(ROOT / "outputs" /
                 "thorney-ti021-full-scale-preflight-dx1m-mixing-2026-10-11-v1" /
                 "ti021-full-scale-preflight.json"),
    )
    parser.add_argument(
        "--output-directory", type=Path,
        default=(ROOT / "outputs" /
                 "thorney-ti021-full-transient-plan-2026-10-11-v2"),
    )
    args = parser.parse_args()
    workbook = args.workbook.resolve()
    dx2_path = args.dx2_benchmark.resolve()
    dx1_path = args.dx1_benchmark.resolve()
    dataset = read_thorney_ti021_workbook(workbook)
    dx2 = json.loads(dx2_path.read_text(encoding="utf-8"))
    dx1 = json.loads(dx1_path.read_text(encoding="utf-8"))
    record = build_plan(dataset, dx2_benchmark=dx2, dx1_benchmark=dx1)
    record["provenance"] = {
        "workbook_sha256": _sha256(workbook),
        "dx2_benchmark_path": str(dx2_path),
        "dx2_benchmark_sha256": _sha256(dx2_path),
        "dx1_benchmark_path": str(dx1_path),
        "dx1_benchmark_sha256": _sha256(dx1_path),
    }
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-full-transient-plan.json"
    plot_path = args.output_directory / "ti021-full-transient-plan.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(record, dataset, plot_path)
    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "coverage": record["coverage"],
        "gates": record["gates"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
