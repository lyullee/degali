"""Grid audit for the TI021 zero-thickness curved-fence representation.

The final report gives an impermeable membrane geometry but no physical wall
thickness.  The numerical policy is therefore a one-cell radial solid band at
each grid level.  Its thickness shrinks with the mesh and is not presented as
a measured dimension.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.site_geometry import CircularArcWall


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "thorney-ti021-curved-wall-geometry-audit-2026-10-11-v1"
SCHEMA = "degali.thorney-ti021-curved-wall-geometry-audit.v1"


def _component_count(mask: np.ndarray) -> int:
    remaining = set(map(tuple, np.argwhere(mask)))
    components = 0
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            row, column = stack.pop()
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == dc == 0:
                        continue
                    neighbour = (row + dr, column + dc)
                    if neighbour in remaining:
                        remaining.remove(neighbour)
                        stack.append(neighbour)
    return components


def _case(cell_size_m: float) -> tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray]:
    domain_m = (130.0, 140.0, 10.0)
    cells = tuple(int(round(length / cell_size_m)) for length in domain_m)
    x = (np.arange(cells[0]) + 0.5) * cell_size_m
    y = (np.arange(cells[1]) + 0.5) * cell_size_m
    z = (np.arange(cells[2]) + 0.5) * cell_size_m
    wall = CircularArcWall(
        center_m=(60.0, 60.0), radius_m=50.0,
        thickness_m=cell_size_m, z_min_m=0.0, z_max_m=5.0,
        arc_start_deg=-90.0, arc_end_deg=90.0,
        label="TI021 report-priority 5 m curved fence",
    )
    mask = wall.voxel_mask(x, y, z)
    horizontal = np.any(mask, axis=0)
    indices = np.argwhere(horizontal)
    radii = np.hypot(
        x[indices[:, 1]] - wall.center_m[0],
        y[indices[:, 0]] - wall.center_m[1],
    )
    estimated_length = float(np.count_nonzero(horizontal) * cell_size_m)
    target_length = math.pi * wall.radius_m
    active_z = np.flatnonzero(np.any(mask, axis=(1, 2)))
    effective_height = float(active_z.size * cell_size_m)
    source_index = (
        int(np.argmin(abs(z - 0.5 * cell_size_m))),
        int(np.argmin(abs(y - wall.center_m[1]))),
        int(np.argmin(abs(x - wall.center_m[0]))),
    )
    record = {
        "cell_size_m": cell_size_m,
        "cells": list(cells),
        "numerical_wall_thickness_m": cell_size_m,
        "thickness_interpretation": "one-cell numerical surface regularisation",
        "horizontal_solid_cell_count": int(np.count_nonzero(horizontal)),
        "solid_cell_count": int(np.count_nonzero(mask)),
        "eight_connected_component_count": _component_count(horizontal),
        "minimum_solid_center_radius_m": float(np.min(radii)),
        "maximum_solid_center_radius_m": float(np.max(radii)),
        "maximum_solid_center_radial_error_m": float(
            np.max(abs(radii - wall.radius_m))
        ),
        "maximum_solid_center_radial_error_cells": float(
            np.max(abs(radii - wall.radius_m)) / cell_size_m
        ),
        "estimated_arc_length_m": estimated_length,
        "target_arc_length_m": target_length,
        "relative_arc_length_error": abs(estimated_length - target_length) / target_length,
        "effective_height_m": effective_height,
        "height_error_m": effective_height - wall.z_max_m,
        "source_cell_is_solid": bool(mask[source_index]),
    }
    return record, x, y, horizontal


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    cases_with_fields = [_case(dx) for dx in (1.0, 0.5, 0.25)]
    cases = [item[0] for item in cases_with_fields]
    length_errors = [float(case["relative_arc_length_error"]) for case in cases]
    gates = {
        "single_connected_wall_each_grid": all(
            case["eight_connected_component_count"] == 1 for case in cases
        ),
        "source_remains_fluid_each_grid": all(
            not case["source_cell_is_solid"] for case in cases
        ),
        "radial_center_error_at_most_half_cell": all(
            case["maximum_solid_center_radial_error_cells"] <= 0.5 + 1.0e-12
            for case in cases
        ),
        "height_exact_each_grid": all(
            abs(case["height_error_m"]) <= 1.0e-12 for case in cases
        ),
        "arc_length_error_decreases": all(
            right < left for left, right in zip(length_errors, length_errors[1:])
        ),
        "finest_arc_length_error_at_most_2_percent": length_errors[-1] <= 0.02,
    }
    passed = all(gates.values())
    record = {
        "schema": SCHEMA,
        "status": "accepted_geometry_component" if passed else "withheld_geometry_component",
        "geometry_source": {
            "dataset_id": "thorney-island-phase2-trial-021",
            "report_priority_height_m": 5.0,
            "radius_m": 50.0,
            "shape": "semi-circular impermeable fence",
        },
        "policy": {
            "physical_thickness_known": False,
            "representation": "one-cell annular-sector solid band",
            "thickness_shrinks_with_grid": True,
            "prohibited_inference": "cell-band thickness is not a measured fence thickness",
            "solver_coordinate_transform": (
                "solver +x is site +y (ideal downwind); solver +y is site -x; "
                "site arc 0..180 deg becomes solver -90..90 deg"
            ),
        },
        "cases": cases,
        "gates": gates,
        "all_gates_pass": passed,
        "scope": (
            "voxel geometry only; does not validate R12 thermodynamics, wind, "
            "turbulence, concentration, dose or hydrogen accuracy"
        ),
    }
    json_path = OUTPUT / "curved-wall-geometry-audit.json"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.4), constrained_layout=True)
    for axis, (case, x, y, horizontal) in zip(axes, cases_with_fields):
        axis.pcolormesh(x, y, horizontal, shading="nearest", cmap="Greys")
        theta = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 361)
        axis.plot(60.0 + 50.0 * np.cos(theta), 60.0 + 50.0 * np.sin(theta),
                  color="crimson", linewidth=1.2, label="declared centreline")
        axis.scatter([60.0], [60.0], marker="*", color="orange", edgecolor="black")
        axis.set_aspect("equal")
        axis.set_xlim(55.0, 115.0)
        axis.set_ylim(5.0, 115.0)
        axis.set_title(
            f"dx={case['cell_size_m']:.2f} m\n"
            f"arc error={100*case['relative_arc_length_error']:.2f}%"
        )
        axis.set_xlabel("solver x / site downwind (m)")
    axes[0].set_ylabel("solver y / site crosswind (m)")
    axes[-1].legend(loc="upper right", fontsize=8)
    plot_path = OUTPUT / "curved-wall-grid-comparison.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    print(json.dumps({
        "status": record["status"], "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()), "gates": gates,
    }, indent=2, ensure_ascii=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
