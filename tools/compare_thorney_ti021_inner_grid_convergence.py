"""Compare matched TI021 5 m and 2.5 m inner-domain model outputs.

This audit reads only model outputs.  It is a two-grid diagnostic and cannot
promote grid convergence without the declared 1.25 m representative-window
third level.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-inner-grid-comparison.v1"
GRID_AND_TIME_CONFIGURATION_KEYS = {
    "duration_s", "time_step_s", "checkpoint_interval_s", "cell_size_m",
    "cells", "active_steps",
}


def physics_configuration_differences(
    coarse: dict[str, object], fine: dict[str, object],
) -> list[str]:
    """Return non-grid/time configuration differences between the two runs."""
    coarse_physics = {
        key: value for key, value in coarse.items()
        if key not in GRID_AND_TIME_CONFIGURATION_KEYS
    }
    fine_physics = {
        key: value for key, value in fine.items()
        if key not in GRID_AND_TIME_CONFIGURATION_KEYS
    }
    differences = []
    for key in sorted(set(coarse_physics) | set(fine_physics)):
        if coarse_physics.get(key) != fine_physics.get(key):
            differences.append(key)
    return differences


def obstacle_region_masks(
    x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray,
    solid_mask: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return fixed physical masks around the TI021 semicircular fence.

    The definitions use only report geometry: centre (60, 100) m in solver
    coordinates, radius 50 m and height 5 m.  They are diagnostic sampling
    volumes, not fitted flow zones or substitutes for streamline topology.
    """

    x = np.asarray(x_m, dtype=float)[None, None, :]
    y = np.asarray(y_m, dtype=float)[None, :, None]
    z = np.asarray(z_m, dtype=float)[:, None, None]
    solid = np.asarray(solid_mask, dtype=bool)
    if solid.shape != (z.size, y.size, x.size):
        raise ValueError("solid mask shape must match z-y-x coordinates")
    fluid = ~solid
    radius = np.sqrt((x - 60.0) ** 2 + (y - 100.0) ** 2)
    downstream_half = x >= 60.0
    front = (
        downstream_half & (radius >= 40.0) & (radius < 50.0)
        & (z < 5.0) & fluid
    )
    roof = (
        downstream_half & (radius >= 45.0) & (radius <= 55.0)
        & (z >= 5.0) & (z < 10.0) & fluid
    )
    tip_distance = np.minimum(
        np.sqrt((x - 60.0) ** 2 + (y - 50.0) ** 2),
        np.sqrt((x - 60.0) ** 2 + (y - 150.0) ** 2),
    )
    side = (tip_distance <= 10.0) & (z < 5.0) & fluid
    wake = (
        downstream_half & (radius > 50.0) & (radius <= 70.0)
        & (z < 7.5) & fluid
    )
    return {
        "windward_inner_face_shell": front,
        "roof_shear_band": roof,
        "side_tip_bypass": side,
        "exterior_lee_wake_shell": wake,
    }


def obstacle_region_metrics(
    field: dict[str, np.ndarray], *, cell_size_m: float,
) -> dict[str, dict[str, float | int]]:
    masks = obstacle_region_masks(
        field["x_m"], field["y_m"], field["z_m"], field["solid_mask"],
    )
    final = np.asarray(field["concentration_vol_percent"], dtype=float)
    maximum = np.asarray(
        field["maximum_concentration_vol_percent"], dtype=float,
    )
    result = {}
    for label, mask in masks.items():
        if not np.any(mask):
            raise ValueError(f"obstacle region {label} contains no fluid cells")
        result[label] = {
            "fluid_cell_count": int(np.count_nonzero(mask)),
            "sampling_volume_m3": float(np.count_nonzero(mask) * cell_size_m**3),
            "final_mean_vol_percent": float(np.mean(final[mask])),
            "final_maximum_vol_percent": float(np.max(final[mask])),
            "maximum_ever_mean_vol_percent": float(np.mean(maximum[mask])),
            "maximum_ever_peak_vol_percent": float(np.max(maximum[mask])),
            "maximum_ever_four_vol_percent_volume_m3": float(
                np.count_nonzero(mask & (maximum >= 4.0)) * cell_size_m**3
            ),
        }
    return result


def coarsen_two_to_one(field: np.ndarray) -> np.ndarray:
    """Volume-average a cell-centred 3-D field by exactly two per axis."""

    values = np.asarray(field, dtype=float)
    if values.ndim != 3 or any(size % 2 for size in values.shape):
        raise ValueError("fine field must be 3-D with even dimensions")
    nz, ny, nx = values.shape
    return np.nanmean(
        values.reshape(nz // 2, 2, ny // 2, 2, nx // 2, 2),
        axis=(1, 3, 5),
    )


def coarsen_fluid_mask_two_to_one(solid_mask: np.ndarray) -> np.ndarray:
    """Require all eight fine cells to be fluid for an unbiased comparison."""

    solid = np.asarray(solid_mask, dtype=bool)
    if solid.ndim != 3 or any(size % 2 for size in solid.shape):
        raise ValueError("fine solid mask must be 3-D with even dimensions")
    nz, ny, nx = solid.shape
    return ~np.any(
        solid.reshape(nz // 2, 2, ny // 2, 2, nx // 2, 2),
        axis=(1, 3, 5),
    )


def trace_metrics(
    time_s: np.ndarray, traces_vol_percent: np.ndarray, *, threshold: float,
) -> dict[str, np.ndarray]:
    """Return per-receptor peak, dose and first threshold arrival."""

    time_values = np.asarray(time_s, dtype=float)
    traces = np.asarray(traces_vol_percent, dtype=float)
    if (
        time_values.ndim != 1 or traces.ndim != 2
        or traces.shape[1] != time_values.size
        or np.any(np.diff(time_values) <= 0.0)
    ):
        raise ValueError("trace array must be receptor-by-increasing-time")
    peak = np.max(traces, axis=1)
    dose = np.trapezoid(traces, time_values, axis=1)
    arrival = np.full(traces.shape[0], np.nan)
    for index, trace in enumerate(traces):
        hit = np.flatnonzero(trace >= threshold)
        if hit.size:
            arrival[index] = time_values[int(hit[0])]
    return {"peak_vol_percent": peak, "dose_vol_percent_s": dose, "arrival_s": arrival}


def _field_difference(
    coarse: np.ndarray, fine_coarsened: np.ndarray, common_fluid: np.ndarray,
) -> dict[str, float]:
    left = np.asarray(coarse, dtype=float)[common_fluid]
    right = np.asarray(fine_coarsened, dtype=float)[common_fluid]
    difference = left - right
    return {
        "common_cell_count": int(left.size),
        "maximum_absolute_difference_vol_percent": float(np.max(np.abs(difference))),
        "root_mean_square_difference_vol_percent": float(np.sqrt(np.mean(difference**2))),
        "relative_l1_difference": float(
            np.sum(np.abs(difference)) / max(np.sum(np.abs(right)), 1.0e-30)
        ),
    }


def _relative_summary(coarse: np.ndarray, fine: np.ndarray, eligible: np.ndarray):
    absolute = np.abs(np.asarray(coarse) - np.asarray(fine))
    relative = absolute[eligible] / np.maximum(np.abs(np.asarray(fine)[eligible]), 1.0e-30)
    return {
        "eligible_receptor_count": int(np.count_nonzero(eligible)),
        "mean_absolute_difference": float(np.mean(absolute[eligible])) if np.any(eligible) else None,
        "maximum_absolute_difference": float(np.max(absolute[eligible])) if np.any(eligible) else None,
        "mean_absolute_percentage_difference": float(100.0 * np.mean(relative)) if np.any(eligible) else None,
        "maximum_absolute_percentage_difference": float(100.0 * np.max(relative)) if np.any(eligible) else None,
    }


def compare_grid_runs(
    *, coarse_run: Path, fine_run: Path, output_directory: Path,
) -> tuple[dict[str, object], dict[str, Path]]:
    coarse_directory = coarse_run.resolve()
    fine_directory = fine_run.resolve()
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    coarse_record_path = coarse_directory / "ti021-inner-transient.json"
    fine_record_path = fine_directory / "ti021-inner-transient.json"
    coarse_record = json.loads(coarse_record_path.read_text(encoding="utf-8"))
    fine_record = json.loads(fine_record_path.read_text(encoding="utf-8"))
    if not (coarse_record.get("all_gates_pass") and fine_record.get("all_gates_pass")):
        raise ValueError("both grid runs must pass their numerical gates")
    coarse_dx = float(coarse_record["configuration"]["cell_size_m"])
    fine_dx = float(fine_record["configuration"]["cell_size_m"])
    if not np.isclose(coarse_dx / fine_dx, 2.0):
        raise ValueError("grid comparison requires an exact 2:1 cell-size ratio")

    with np.load(
        coarse_directory / "ti021-inner-final-field.npz", allow_pickle=False,
    ) as source:
        coarse_field = {key: np.asarray(source[key]).copy() for key in source.files}
    with np.load(
        fine_directory / "ti021-inner-final-field.npz", allow_pickle=False,
    ) as source:
        fine_field = {key: np.asarray(source[key]).copy() for key in source.files}
    if tuple(2 * size for size in coarse_field["solid_mask"].shape) != fine_field["solid_mask"].shape:
        raise ValueError("fine grid dimensions must be exactly twice coarse dimensions")
    fine_fluid_coarse = coarsen_fluid_mask_two_to_one(fine_field["solid_mask"])
    common_fluid = (~coarse_field["solid_mask"]) & fine_fluid_coarse
    fine_final_coarse = coarsen_two_to_one(fine_field["concentration_vol_percent"])
    fine_maximum_coarse = coarsen_two_to_one(
        fine_field["maximum_concentration_vol_percent"]
    )
    field_differences = {
        "final_concentration": _field_difference(
            coarse_field["concentration_vol_percent"], fine_final_coarse, common_fluid,
        ),
        "maximum_concentration": _field_difference(
            coarse_field["maximum_concentration_vol_percent"],
            fine_maximum_coarse, common_fluid,
        ),
    }
    coarse_hazard_volume = float(np.count_nonzero(
        (~coarse_field["solid_mask"])
        & (coarse_field["maximum_concentration_vol_percent"] >= 4.0)
    ) * coarse_dx**3)
    fine_hazard_volume = float(np.count_nonzero(
        (~fine_field["solid_mask"])
        & (fine_field["maximum_concentration_vol_percent"] >= 4.0)
    ) * fine_dx**3)

    with np.load(
        coarse_directory / "ti021-inner-receptor-traces.npz", allow_pickle=False,
    ) as source:
        coarse_trace = {key: np.asarray(source[key]).copy() for key in source.files}
    with np.load(
        fine_directory / "ti021-inner-receptor-traces.npz", allow_pickle=False,
    ) as source:
        fine_trace = {key: np.asarray(source[key]).copy() for key in source.files}
    if not np.array_equal(coarse_trace["labels"], fine_trace["labels"]):
        raise ValueError("coarse and fine receptor labels must match exactly")
    if not np.allclose(
        coarse_trace["solver_positions_m"], fine_trace["solver_positions_m"],
        rtol=0.0, atol=1.0e-12,
    ):
        raise ValueError("coarse and fine receptor coordinates must match")
    coarse_metrics = trace_metrics(
        coarse_trace["time_s"], coarse_trace["concentration_vol_percent"],
        threshold=0.1,
    )
    fine_metrics = trace_metrics(
        fine_trace["time_s"], fine_trace["concentration_vol_percent"],
        threshold=0.1,
    )
    peak_eligible = fine_metrics["peak_vol_percent"] >= 0.1
    dose_floor_vol_percent_s = 0.1
    dose_eligible = (
        fine_metrics["dose_vol_percent_s"] >= dose_floor_vol_percent_s
    )
    both_arrived = np.isfinite(coarse_metrics["arrival_s"]) & np.isfinite(
        fine_metrics["arrival_s"]
    )
    arrival_difference = np.abs(
        coarse_metrics["arrival_s"][both_arrived]
        - fine_metrics["arrival_s"][both_arrived]
    )
    receptor_differences = {
        "threshold_vol_percent": 0.1,
        "peak": _relative_summary(
            coarse_metrics["peak_vol_percent"],
            fine_metrics["peak_vol_percent"], peak_eligible,
        ),
        "dose": _relative_summary(
            coarse_metrics["dose_vol_percent_s"],
            fine_metrics["dose_vol_percent_s"], dose_eligible,
        ),
        "arrival": {
            "both_detected_count": int(np.count_nonzero(both_arrived)),
            "coarse_only_count": int(np.count_nonzero(
                np.isfinite(coarse_metrics["arrival_s"])
                & ~np.isfinite(fine_metrics["arrival_s"])
            )),
            "fine_only_count": int(np.count_nonzero(
                ~np.isfinite(coarse_metrics["arrival_s"])
                & np.isfinite(fine_metrics["arrival_s"])
            )),
            "mean_absolute_difference_s": (
                None if arrival_difference.size == 0
                else float(np.mean(arrival_difference))
            ),
            "maximum_absolute_difference_s": (
                None if arrival_difference.size == 0
                else float(np.max(arrival_difference))
            ),
        },
    }
    hazard_relative = abs(coarse_hazard_volume - fine_hazard_volume) / max(
        fine_hazard_volume, 1.0e-30,
    )
    regional_metrics = {
        "definitions": {
            "windward_inner_face_shell": "40 <= radial distance < 50 m, x >= 60 m, z < 5 m",
            "roof_shear_band": "45 <= radial distance <= 55 m, x >= 60 m, 5 <= z < 10 m",
            "side_tip_bypass": "within 10 m of either fence endpoint, z < 5 m",
            "exterior_lee_wake_shell": "50 < radial distance <= 70 m, x >= 60 m, z < 7.5 m",
            "interpretation": (
                "fixed report-geometry sampling volumes; labels do not infer "
                "separation or recirculation from concentration alone"
            ),
        },
        "coarse": obstacle_region_metrics(coarse_field, cell_size_m=coarse_dx),
        "fine": obstacle_region_metrics(fine_field, cell_size_m=fine_dx),
    }
    physics_differences = physics_configuration_differences(
        coarse_record["configuration"], fine_record["configuration"],
    )
    gates = {
        "matched_time_horizon": (
            coarse_record["configuration"]["duration_s"]
            == fine_record["configuration"]["duration_s"]
        ),
        "matched_physics_except_grid_and_time_step": not physics_differences,
        "both_numerically_accepted": True,
        "third_grid_or_representative_window_available": False,
    }
    record = {
        "schema": SCHEMA,
        "status": (
            "diagnostic_two_grid_physics_mismatch"
            if physics_differences
            else "diagnostic_two_grid_comparison_pending_1p25m_window"
        ),
        "scope": (
            "observation-blind 5 m versus 2.5 m model-output comparison; a "
            "third 1.25 m representative window is mandatory for convergence"
        ),
        "provenance": {
            "coarse_record": {"path": str(coarse_record_path), "sha256": _sha256(coarse_record_path)},
            "fine_record": {"path": str(fine_record_path), "sha256": _sha256(fine_record_path)},
            "observation_value_fields_accessed": [],
        },
        "configuration": {
            "coarse_cell_size_m": coarse_dx,
            "fine_cell_size_m": fine_dx,
            "grid_ratio": coarse_dx / fine_dx,
            "coarse_time_step_s": coarse_record["configuration"]["time_step_s"],
            "fine_time_step_s": fine_record["configuration"]["time_step_s"],
            "duration_s": coarse_record["configuration"]["duration_s"],
            "receptor_count": int(coarse_trace["labels"].size),
            "physics_configuration_differences": physics_differences,
            "peak_eligibility_floor_vol_percent": 0.1,
            "dose_eligibility_floor_vol_percent_s": dose_floor_vol_percent_s,
            "common_fully_fluid_coarse_cell_count": int(np.count_nonzero(common_fluid)),
        },
        "field_differences": field_differences,
        "receptor_differences": receptor_differences,
        "obstacle_region_metrics": regional_metrics,
        "four_vol_percent_maximum_hazard_volume": {
            "coarse_m3": coarse_hazard_volume,
            "fine_m3": fine_hazard_volume,
            "absolute_difference_m3": abs(coarse_hazard_volume - fine_hazard_volume),
            "relative_difference": hazard_relative,
        },
        "gates": {key: bool(value) for key, value in gates.items()},
        "grid_convergence_promoted": False,
        "next_required": (
            "run a matched 1.25 m representative plume/obstacle window and "
            "evaluate monotonicity plus an observed order or GCI"
        ),
    }
    json_path = output / "ti021-inner-grid-comparison.json"
    plot_path = output / "ti021-inner-grid-comparison.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
    )

    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), constrained_layout=True)
    axes[0].scatter(
        fine_metrics["peak_vol_percent"], coarse_metrics["peak_vol_percent"],
        s=22, alpha=0.8,
    )
    limit = max(float(np.max(fine_metrics["peak_vol_percent"])), float(np.max(coarse_metrics["peak_vol_percent"])), 0.1)
    axes[0].plot((0.0, limit), (0.0, limit), "k--", linewidth=1.0)
    axes[0].set(xlabel="2.5 m peak (vol%)", ylabel="5 m peak (vol%)", title="60 receptor peaks")
    axes[0].grid(alpha=0.25)
    axes[1].scatter(
        fine_metrics["dose_vol_percent_s"], coarse_metrics["dose_vol_percent_s"],
        s=22, alpha=0.8,
    )
    dose_limit = max(float(np.max(fine_metrics["dose_vol_percent_s"])), float(np.max(coarse_metrics["dose_vol_percent_s"])), 1.0)
    axes[1].plot((0.0, dose_limit), (0.0, dose_limit), "k--", linewidth=1.0)
    axes[1].set(xlabel="2.5 m dose (vol% s)", ylabel="5 m dose (vol% s)", title="60 receptor doses")
    axes[1].grid(alpha=0.25)
    labels = ("5 m", "2.5 m")
    axes[2].bar(labels, (coarse_hazard_volume, fine_hazard_volume), color=("#4C78A8", "#F58518"))
    axes[2].set(ylabel="Maximum-ever C >= 4 vol% volume (m3)", title="Flammable-volume sensitivity")
    axes[2].grid(axis="y", alpha=0.25)
    figure.suptitle("TI021 inner-domain two-grid diagnostic (not a convergence claim)")
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    return record, {"json": json_path, "plot": plot_path}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse-run", type=Path, required=True)
    parser.add_argument("--fine-run", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = compare_grid_runs(
        coarse_run=args.coarse_run,
        fine_run=args.fine_run,
        output_directory=args.output_directory,
    )
    print(json.dumps({
        "status": record["status"],
        "grid_convergence_promoted": record["grid_convergence_promoted"],
        "field_differences": record["field_differences"],
        "receptor_differences": record["receptor_differences"],
        "four_vol_percent_maximum_hazard_volume": record[
            "four_vol_percent_maximum_hazard_volume"
        ],
        "outputs": {key: str(path) for key, path in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
