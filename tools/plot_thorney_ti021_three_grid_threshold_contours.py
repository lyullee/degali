"""Plot matched TI021 maximum-concentration threshold contours by grid.

This is a diagnostic visualization only. It does not score experiment data or
promote the R12 surrogate calculation as hydrogen validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "degali.thorney-ti021-three-grid-threshold-contours.v1"
EXPECTED_GRIDS_M = (5.0, 2.5, 1.25)
LEVELS_VOL_PERCENT = (0.0, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
THRESHOLD_VOL_PERCENT = 4.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_run(root: Path) -> dict[str, object]:
    record_path = root / "ti021-inner-transient.json"
    field_path = root / "ti021-inner-final-field.npz"
    if not record_path.is_file() or not field_path.is_file():
        raise FileNotFoundError(f"missing transient record or field in {root}")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("status") not in {
        "accepted_inner_transient_qualification",
        "accepted_completed_inner_transient",
    }:
        raise ValueError(f"run is not a completed accepted run: {root}")
    with np.load(field_path, allow_pickle=False) as source:
        field = {key: source[key].copy() for key in source.files}
    required = {
        "x_m", "y_m", "z_m", "solid_mask",
        "maximum_concentration_vol_percent",
    }
    if not required.issubset(field):
        raise ValueError(f"field in {root} lacks required arrays")
    concentration = np.asarray(field["maximum_concentration_vol_percent"])
    solid = np.asarray(field["solid_mask"], dtype=bool)
    if concentration.shape != solid.shape or concentration.ndim != 3:
        raise ValueError(f"concentration/solid array shapes are invalid in {root}")
    return {
        "root": root,
        "record": record,
        "field": field,
        "record_path": record_path,
        "field_path": field_path,
    }


def _masked(array: np.ndarray, solid: np.ndarray) -> np.ma.MaskedArray:
    values = np.asarray(array, dtype=float)
    mask = np.asarray(solid, dtype=bool) | ~np.isfinite(values)
    return np.ma.array(values, mask=mask)


def _interpolate_cell_center_slice(
    values: np.ndarray,
    solid: np.ndarray,
    coordinates: np.ndarray,
    *,
    axis: int,
    target: float,
) -> np.ma.MaskedArray:
    """Interpolate a 3-D cell-centred field to one shared physical plane."""

    array = np.asarray(values, dtype=float)
    obstacle = np.asarray(solid, dtype=bool)
    coord = np.asarray(coordinates, dtype=float)
    if array.shape != obstacle.shape or array.shape[axis] != coord.size:
        raise ValueError("field, mask, and slice coordinate dimensions differ")
    if coord.size == 0 or np.any(np.diff(coord) <= 0.0):
        raise ValueError("slice coordinates must be nonempty and strictly increasing")
    if target < coord[0] or target > coord[-1]:
        raise ValueError(f"target {target:g} lies outside cell-centre range")
    upper = int(np.searchsorted(coord, target, side="left"))
    if upper < coord.size and math.isclose(float(coord[upper]), target,
                                           rel_tol=0.0, abs_tol=1.0e-10):
        result = np.take(array, upper, axis=axis)
        result_solid = np.take(obstacle, upper, axis=axis)
        return _masked(result, result_solid)
    if upper == 0 or upper == coord.size:
        raise ValueError("target must be bracketed by cell centres")
    lower = upper - 1
    fraction = float((target - coord[lower]) / (coord[upper] - coord[lower]))
    low_values = np.take(array, lower, axis=axis)
    high_values = np.take(array, upper, axis=axis)
    low_valid = ~np.take(obstacle, lower, axis=axis) & np.isfinite(low_values)
    high_valid = ~np.take(obstacle, upper, axis=axis) & np.isfinite(high_values)
    both = low_valid & high_valid
    output = np.zeros_like(low_values, dtype=float)
    output[both] = (
        (1.0 - fraction) * low_values[both] + fraction * high_values[both]
    )
    output[low_valid & ~high_valid] = low_values[low_valid & ~high_valid]
    output[high_valid & ~low_valid] = high_values[high_valid & ~low_valid]
    return np.ma.array(output, mask=~(low_valid | high_valid))


def plot_three_grid_threshold_contours(
    coarse_run: Path,
    medium_run: Path,
    fine_run: Path,
    output_directory: Path,
) -> tuple[dict[str, object], dict[str, Path]]:
    runs = [_load_run(path.resolve()) for path in (coarse_run, medium_run, fine_run)]
    configurations = [item["record"]["configuration"] for item in runs]
    cells = [float(item["cell_size_m"]) for item in configurations]
    if not np.allclose(cells, EXPECTED_GRIDS_M, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"expected grid order {EXPECTED_GRIDS_M}, got {cells}")
    for key in (
        "duration_s", "time_step_s", "corrector_tolerance",
        "circular_arc_representation", "workbook_sha256", "property_audit_sha256",
    ):
        values = [item.get(key) for item in configurations]
        if any(value != values[0] for value in values[1:]):
            raise ValueError(f"runs differ in required matched field {key}: {values}")
    if configurations[0].get("circular_arc_representation") != "zero_thickness_surface":
        raise ValueError("contour comparison requires the matched zero-thickness surface")

    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    output_png = output / "ti021-three-grid-maximum-concentration-contours.png"
    output_json = output / "ti021-three-grid-maximum-concentration-contours.json"
    if output_png.exists() or output_json.exists():
        raise FileExistsError("threshold contour outputs already exist; choose a new output directory")

    fig, axes = plt.subplots(3, 3, figsize=(15.5, 12.0), constrained_layout=True)
    norm = matplotlib.colors.BoundaryNorm(LEVELS_VOL_PERCENT, ncolors=256, clip=False)
    cmap = matplotlib.colormaps["turbo"]
    mappable = None
    panels: list[dict[str, object]] = []

    for col, run in enumerate(runs):
        field = run["field"]
        record = run["record"]
        conf = record["configuration"]
        x = np.asarray(field["x_m"], dtype=float)
        y = np.asarray(field["y_m"], dtype=float)
        z = np.asarray(field["z_m"], dtype=float)
        solid = np.asarray(field["solid_mask"], dtype=bool)
        max_c = np.asarray(field["maximum_concentration_vol_percent"], dtype=float)
        if max_c.shape != (z.size, y.size, x.size):
            raise ValueError(f"field axis order is invalid for {run['root']}")
        dx = float(conf["cell_size_m"])
        axes[0, col].set_title(f"{dx:g} m grid", fontsize=13)

        # Interpolate to the same physical plane on every grid.
        plan = _interpolate_cell_center_slice(
            max_c, solid, z, axis=0, target=2.5,
        )
        ax = axes[0, col]
        mappable = ax.contourf(x, y, plan, levels=LEVELS_VOL_PERCENT, norm=norm,
                               cmap=cmap, extend="max")
        if plan.count() and float(plan.min()) <= THRESHOLD_VOL_PERCENT <= float(plan.max()):
            ax.contour(x, y, plan, levels=[THRESHOLD_VOL_PERCENT], colors="black",
                       linewidths=1.1)
        theta = np.linspace(-0.5 * math.pi, 0.5 * math.pi, 240)
        ax.plot(60.0 + 50.0 * np.cos(theta), 100.0 + 50.0 * np.sin(theta),
                color="white", linewidth=1.0, linestyle="--")
        ax.plot(60.0, 100.0, marker="+", color="white", markersize=7)
        ax.set(xlim=(0, 200), ylim=(0, 200), ylabel=("y (m)" if col == 0 else ""))
        ax.set_aspect("equal")
        above = int(np.count_nonzero(plan.filled(-np.inf) >= THRESHOLD_VOL_PERCENT))
        panels.append({"grid_m": dx, "view": "plan_xy", "slice_z_m": 2.5,
                       "cells_at_or_above_4vol_percent": above,
                       "maximum_vol_percent": float(plan.max()) if plan.count() else None})

        # Streamwise vertical section through the arc apex and release centreline.
        xz = _interpolate_cell_center_slice(
            max_c, solid, y, axis=1, target=100.0,
        )
        ax = axes[1, col]
        mappable = ax.contourf(x, z, xz, levels=LEVELS_VOL_PERCENT, norm=norm,
                               cmap=cmap, extend="max")
        if xz.count() and float(xz.min()) <= THRESHOLD_VOL_PERCENT <= float(xz.max()):
            ax.contour(x, z, xz, levels=[THRESHOLD_VOL_PERCENT], colors="black",
                       linewidths=1.1)
        ax.plot((110.0, 110.0), (0.0, 5.0), color="white", linewidth=2.0)
        ax.set(xlim=(0, 200), ylim=(0, 20), ylabel=("z (m)" if col == 0 else ""))
        above = int(np.count_nonzero(xz.filled(-np.inf) >= THRESHOLD_VOL_PERCENT))
        panels.append({"grid_m": dx, "view": "streamwise_xz", "slice_y_m": 100.0,
                       "cells_at_or_above_4vol_percent": above,
                       "maximum_vol_percent": float(xz.max()) if xz.count() else None})

        # Crosswind vertical section through the release centre; exposes the two arc tips.
        yz = _interpolate_cell_center_slice(
            max_c, solid, x, axis=2, target=60.0,
        )
        ax = axes[2, col]
        mappable = ax.contourf(y, z, yz, levels=LEVELS_VOL_PERCENT, norm=norm,
                               cmap=cmap, extend="max")
        if yz.count() and float(yz.min()) <= THRESHOLD_VOL_PERCENT <= float(yz.max()):
            ax.contour(y, z, yz, levels=[THRESHOLD_VOL_PERCENT], colors="black",
                       linewidths=1.1)
        for wall_y in (50.0, 150.0):
            ax.plot((wall_y, wall_y), (0.0, 5.0), color="white", linewidth=2.0)
        ax.set(xlim=(0, 200), ylim=(0, 20), xlabel=("y (m)" if col == 1 else ""),
               ylabel=("z (m)" if col == 0 else ""))
        above = int(np.count_nonzero(yz.filled(-np.inf) >= THRESHOLD_VOL_PERCENT))
        panels.append({"grid_m": dx, "view": "crosswind_yz", "slice_x_m": 60.0,
                       "cells_at_or_above_4vol_percent": above,
                       "maximum_vol_percent": float(yz.max()) if yz.count() else None})

    axes[0, 0].set_ylabel("Plan, z=2.5 m\ny (m)")
    axes[1, 0].set_ylabel("Streamwise, y=100 m\nz (m)")
    axes[2, 0].set_ylabel("Crosswind, x=60 m\nz (m)")
    for ax in axes[0, :]:
        ax.set_xlabel("x (m)")
    for ax in axes[1, :]:
        ax.set_xlabel("x (m)")
    for ax in axes[2, :]:
        ax.set_xlabel("y (m)")
    colorbar = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.86, pad=0.015)
    colorbar.set_label("Maximum R12 tracer concentration (vol%)")
    fig.suptitle(
        "TI021 R12 surrogate — maximum concentration over 0–20 s\n"
        "Matched ideal fence; black contour = 4 vol%; visual diagnostic, not H₂ validation",
        fontsize=15,
    )
    fig.savefig(output_png, dpi=180)
    plt.close(fig)

    record = {
        "schema": SCHEMA,
        "status": "visual_diagnostic_not_validation",
        "scope": "matched 5/2.5/1.25 m zero-thickness TI021 R12 0-20 s concentration-envelope slices",
        "threshold_vol_percent": THRESHOLD_VOL_PERCENT,
        "field": "maximum_concentration_vol_percent",
        "slice_interpolation": (
            "linear interpolation from cell centres to shared physical planes "
            "z=2.5m, y=100m, and x=60m; at solid cells use the adjacent fluid value"
        ),
        "shared_configuration": {
            key: configurations[0].get(key) for key in (
                "duration_s", "time_step_s", "corrector_tolerance",
                "circular_arc_representation", "workbook_sha256", "property_audit_sha256",
            )
        },
        "contour_levels_vol_percent": list(LEVELS_VOL_PERCENT),
        "panels": panels,
        "provenance": [{
            "run_record": str(run["record_path"]),
            "run_record_sha256": _sha256(run["record_path"]),
            "field": str(run["field_path"]),
            "field_sha256": _sha256(run["field_path"]),
        } for run in runs],
        "interpretation_warning": (
            "R12 surrogate, idealized fence and 20 s numerical window only; "
            "do not call this a hydrogen flammability prediction or validated field."
        ),
    }
    output_json.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")
    return record, {"plot": output_png, "record": output_json}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse-run", type=Path, required=True)
    parser.add_argument("--medium-run", type=Path, required=True)
    parser.add_argument("--fine-run", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = plot_three_grid_threshold_contours(
        args.coarse_run, args.medium_run, args.fine_run, args.output_directory,
    )
    print(json.dumps({"status": record["status"], **{k: str(v) for k, v in paths.items()}},
                     indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
