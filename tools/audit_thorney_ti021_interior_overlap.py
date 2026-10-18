"""Compare shared TI021 receptors inside an actual nested-domain overlap.

This audit is observation-blind and intentionally supports short transient
windows. It diagnoses spatial handoff consistency only; it does not score
experimental doses or concentrations and cannot promote an incomplete event.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from tools.score_thorney_ti021_nested_145_sensors import (
    _interpolate_rows,
    _overlap_geometry_diagnostic,
    _trace_metrics,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-interior-overlap-audit.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_audit(
    *, inner_run_directory: Path, outer_run_directory: Path,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    inner_dir = inner_run_directory.resolve()
    outer_dir = outer_run_directory.resolve()
    inner_record_path = inner_dir / "ti021-inner-transient.json"
    outer_record_path = outer_dir / "ti021-outer-transient.json"
    inner_trace_path = inner_dir / "ti021-inner-receptor-traces.npz"
    outer_trace_path = outer_dir / "ti021-outer-receptor-traces.npz"
    inner_record = json.loads(inner_record_path.read_text(encoding="utf-8"))
    outer_record = json.loads(outer_record_path.read_text(encoding="utf-8"))
    if not inner_record.get("all_gates_pass") or not outer_record.get("all_gates_pass"):
        raise ValueError("both transient records must pass their numerical gates")
    inner_duration = float(inner_record["configuration"]["duration_s"])
    outer_duration = float(outer_record["configuration"]["duration_s"])
    if not np.isclose(inner_duration, outer_duration):
        raise ValueError("inner and outer durations must match for this audit")

    with np.load(inner_trace_path, allow_pickle=False) as loaded:
        inner = {key: np.asarray(loaded[key]).copy() for key in loaded.files}
    with np.load(outer_trace_path, allow_pickle=False) as loaded:
        outer = {key: np.asarray(loaded[key]).copy() for key in loaded.files}
    inner_labels = [str(value) for value in inner["labels"]]
    outer_labels = [str(value) for value in outer["labels"]]
    common_labels = sorted(set(inner_labels) & set(outer_labels))
    if not common_labels:
        raise ValueError("nested runs have no shared receptor labels")
    geometry = _overlap_geometry_diagnostic(
        inner_record, outer_record, np.asarray(common_labels),
    )

    inner_index = {label: index for index, label in enumerate(inner_labels)}
    outer_index = {label: index for index, label in enumerate(outer_labels)}
    inner_time = inner["time_s"]
    outer_on_inner_clock = _interpolate_rows(
        outer["time_s"], outer["concentration_vol_percent"], inner_time,
    )
    inner_common = np.asarray([inner["concentration_vol_percent"][inner_index[label]] for label in common_labels])
    outer_common = np.asarray([outer_on_inner_clock[outer_index[label]] for label in common_labels])
    inner_metrics = _trace_metrics(inner_time, inner_common)
    outer_metrics = _trace_metrics(inner_time, outer_common)
    peak_delta = np.abs(outer_metrics["peak"] - inner_metrics["peak"])
    dose_delta = np.abs(outer_metrics["dose"] - inner_metrics["dose"])
    peak_relative = peak_delta / np.maximum(np.abs(inner_metrics["peak"]), 0.1)
    dose_relative = dose_delta / np.maximum(np.abs(inner_metrics["dose"]), 0.1)
    interior_labels = set(geometry["interior_overlap_sensor_labels"])
    interior_mask = np.asarray([label in interior_labels for label in common_labels])
    interior_indices = np.flatnonzero(interior_mask)
    if interior_indices.size == 0:
        raise ValueError("no shared sensors lie inside both domains beyond the interface buffer")
    boundary_indices = np.flatnonzero(~interior_mask)
    interior_peak_relative = peak_relative[interior_indices]
    interior_dose_relative = dose_relative[interior_indices]
    overlap_records = []
    for index, label in enumerate(common_labels):
        receptor = next(
            row for row in inner_record["configuration"]["receptors"]
            if row["label"] == label
        )
        overlap_records.append({
            "label": label,
            "site_local_xyz_m": receptor["site_local_xyz_m"],
            "location_class": "interior_overlap" if interior_mask[index] else "boundary_or_noninterior_monitor",
            "inner_peak_vol_percent": float(inner_metrics["peak"][index]),
            "outer_peak_vol_percent": float(outer_metrics["peak"][index]),
            "absolute_peak_difference_vol_percentage_points": float(peak_delta[index]),
            "inner_dose_vol_percent_s": float(inner_metrics["dose"][index]),
            "outer_dose_vol_percent_s": float(outer_metrics["dose"][index]),
            "absolute_dose_difference_vol_percent_s": float(dose_delta[index]),
            "peak_relative_change_using_0p1vol_floor": float(peak_relative[index]),
            "dose_relative_change_using_0p1vol_s_floor": float(dose_relative[index]),
        })

    inner_set, outer_set = set(inner_labels), set(outer_labels)
    expected_labels = {f"TI021-S{index:03d}" for index in range(1, 146)}
    union = inner_set | outer_set
    locations_gate = bool(geometry["at_least_three_interior_horizontal_locations"])
    coordinates_gate = bool(geometry["all_shared_sensor_coordinates_match"])
    peak_gate = bool(np.max(interior_peak_relative) <= 0.10)
    dose_gate = bool(np.max(interior_dose_relative) <= 0.10)
    record: dict[str, object] = {
        "schema": SCHEMA,
        "status": (
            "accepted_short_interior_overlap_diagnostic_only"
            if locations_gate and coordinates_gate and peak_gate and dose_gate else
            "withheld_interior_overlap_gate"
        ),
        "scope": (
            "no-observation spatial coupling diagnostic; not an experimental score, "
            "full-event validation, grid-convergence result, or hydrogen prediction"
        ),
        "provenance": {
            "inner_record": {"path": str(inner_record_path), "sha256": _sha256(inner_record_path)},
            "inner_trace": {"path": str(inner_trace_path), "sha256": _sha256(inner_trace_path)},
            "outer_record": {"path": str(outer_record_path), "sha256": _sha256(outer_record_path)},
            "outer_trace": {"path": str(outer_trace_path), "sha256": _sha256(outer_trace_path)},
        },
        "configuration": {
            "model_horizon_s": inner_duration,
            "inner_domain_m": inner_record["configuration"]["domain_m"],
            "outer_domain_m": outer_record["configuration"]["domain_m"],
            "inner_cell_size_m": inner_record["configuration"]["cell_size_m"],
            "outer_cell_size_m": outer_record["configuration"]["cell_size_m"],
            "handoff_site_downwind_y_m": geometry["handoff_site_downwind_y_m"],
            "observation_values_used_to_configure_runs": [],
        },
        "runtime": {
            "inner_wall_clock_s": inner_record.get("runtime", {}).get("wall_clock_s"),
            "outer_wall_clock_s": outer_record.get("runtime", {}).get("wall_clock_s"),
            "combined_wall_clock_s": (
                float(inner_record["runtime"]["wall_clock_s"])
                + float(outer_record["runtime"]["wall_clock_s"])
            ),
        },
        "numerics": {
            "inner_maximum_continuity_residual_kg_m3_s": inner_record.get("numerics", {}).get("maximum_continuity_residual_kg_m3_s"),
            "outer_maximum_continuity_residual_kg_m3_s": outer_record.get("numerics", {}).get("maximum_continuity_residual_kg_m3_s"),
            "inner_cumulative_relative_species_mass_residual": inner_record.get("numerics", {}).get("maximum_relative_species_mass_residual"),
            "outer_cumulative_relative_species_mass_residual": outer_record.get("numerics", {}).get("cumulative_relative_species_mass_residual"),
        },
        "coverage": {
            "inner_sensor_count": len(inner_set),
            "outer_sensor_count": len(outer_set),
            "shared_sensor_count": len(common_labels),
            "union_sensor_count": len(union),
            "expected_sensor_count": len(expected_labels),
            "uncovered_sensor_labels": sorted(expected_labels - union),
        },
        "geometry": geometry,
        "comparison": {
            "comparison_horizon_s": [float(inner_time[0]), float(inner_time[-1])],
            "outer_trace_interpolation": "linear onto inner event clock",
            "acceptance_scope": "interior_overlap_sensor_locations_only",
            "interior_overlap_sensor_count": int(interior_indices.size),
            "maximum_peak_absolute_difference_vol_percentage_points": float(np.max(peak_delta[interior_indices])),
            "maximum_dose_absolute_difference_vol_percent_s": float(np.max(dose_delta[interior_indices])),
            "maximum_peak_relative_change_using_0p1vol_floor": float(np.max(interior_peak_relative)),
            "maximum_dose_relative_change_using_0p1vol_s_floor": float(np.max(interior_dose_relative)),
            "relative_change_limit": 0.10,
            "peak_10_percent_gate": peak_gate,
            "dose_10_percent_gate": dose_gate,
            "boundary_interface_monitoring": {
                "sensor_count": int(boundary_indices.size),
                "labels": [common_labels[index] for index in boundary_indices],
                "maximum_peak_relative_change": (
                    None if boundary_indices.size == 0 else float(np.max(peak_relative[boundary_indices]))
                ),
                "maximum_dose_relative_change": (
                    None if boundary_indices.size == 0 else float(np.max(dose_relative[boundary_indices]))
                ),
                "used_for_interior_acceptance": False,
            },
            "receptors": overlap_records,
        },
        "gates": {
            "both_transients_passed_numerical_gates": True,
            "shared_sensor_coordinates_match": coordinates_gate,
            "at_least_three_interior_horizontal_locations": locations_gate,
            "interior_peak_difference_at_most_10_percent": peak_gate,
            "interior_dose_difference_at_most_10_percent": dose_gate,
            "full_145_sensor_union_covered": union == expected_labels,
            "full_observation_horizon_scored": False,
            "experimental_model_promotion_allowed": False,
        },
        "interpretation": (
            "This short matched-grid run establishes shared sensor locations inside "
            "the nested domains and measures transient coupling differences only; "
            "the sensor union and full event are not covered by this diagnostic."
        ),
    }
    arrays = {
        "time_s": inner_time,
        "labels": np.asarray(common_labels),
        "inner_traces": inner_common,
        "outer_traces": outer_common,
        "peak_relative": peak_relative,
        "dose_relative": dose_relative,
        "interior_mask": interior_mask,
    }
    return record, arrays


def _plot(record: dict[str, object], arrays: dict[str, np.ndarray], path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 8.0), constrained_layout=True)
    for label in ("TI021-S058", "TI021-S059", "TI021-S060", "TI021-S061",
                  "TI021-S066", "TI021-S067", "TI021-S068", "TI021-S069"):
        members = np.flatnonzero(arrays["labels"] == label)
        if members.size == 0:
            continue
        index = int(members[0])
        axes[0].plot(arrays["time_s"], arrays["inner_traces"][index], label=f"{label} inner")
        axes[0].plot(arrays["time_s"], arrays["outer_traces"][index], "--", label=f"{label} outer")
    axes[0].set(xlabel="time (s)", ylabel="R12 concentration (vol%)",
                title="Shared interior sensor histories (60 s diagnostic; no experimental data)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(ncol=4, fontsize=7)
    positions = np.flatnonzero(arrays["interior_mask"])
    width = 0.38
    axes[1].bar(positions - width / 2, 100.0 * arrays["peak_relative"][positions], width, label="0.6 s peak change")
    axes[1].bar(positions + width / 2, 100.0 * arrays["dose_relative"][positions], width, label="60 s dose change")
    axes[1].axhline(10.0, color="black", linestyle="--", linewidth=1, label="10% diagnostic limit")
    axes[1].set(xticks=positions, xticklabels=arrays["labels"][positions], ylabel="relative change (%)",
                title="Interior overlap sensitivity by shared sensor")
    axes[1].tick_params(axis="x", labelrotation=75, labelsize=6)
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(ncol=3, fontsize=8)
    fig.suptitle("TI021 nested-domain interior overlap — diagnostic only, not experimental validation")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inner-run-directory", type=Path, required=True)
    parser.add_argument("--outer-run-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, arrays = build_audit(
        inner_run_directory=args.inner_run_directory,
        outer_run_directory=args.outer_run_directory,
    )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-interior-overlap-audit.json"
    png_path = args.output_directory / "ti021-interior-overlap-audit.png"
    json_path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _plot(record, arrays, png_path)
    print(json.dumps({
        "status": record["status"],
        "coverage": record["coverage"],
        "interior_overlap_horizontal_location_count": record["geometry"]["interior_overlap_horizontal_location_count"],
        "maximum_peak_relative_change": record["comparison"]["maximum_peak_relative_change_using_0p1vol_floor"],
        "maximum_dose_relative_change": record["comparison"]["maximum_dose_relative_change_using_0p1vol_s_floor"],
        "json": str(json_path.resolve()),
        "png": str(png_path.resolve()),
    }, indent=2))
    return 0 if record["status"] == "accepted_short_interior_overlap_diagnostic_only" else 2


if __name__ == "__main__":
    raise SystemExit(main())
