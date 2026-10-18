"""Merge no-fit inner/outer histories and score all 145 TI021 sensors."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

try:
    from tools.run_thorney_ti021_inner_transient import _sha256
    from tools.score_thorney_ti021_inner_transient import (
        _ratio_metrics,
        trailing_boxcar_average,
    )
except ModuleNotFoundError:
    from run_thorney_ti021_inner_transient import _sha256
    from score_thorney_ti021_inner_transient import _ratio_metrics, trailing_boxcar_average


SCHEMA = "degali.thorney-ti021-nested-145-sensor-score.v1"


def _interpolate_rows(
    source_time_s: np.ndarray, source_values: np.ndarray, target_time_s: np.ndarray,
) -> np.ndarray:
    source_time = np.asarray(source_time_s, dtype=float)
    target_time = np.asarray(target_time_s, dtype=float)
    values = np.asarray(source_values, dtype=float)
    if values.ndim != 2 or values.shape[1] != source_time.size:
        raise ValueError("source trace array must be receptor-by-time")
    if np.any(np.diff(source_time) <= 0.0) or np.any(np.diff(target_time) <= 0.0):
        raise ValueError("source and target clocks must be strictly increasing")
    if target_time[0] < source_time[0] or target_time[-1] > source_time[-1]:
        raise ValueError("target clock must lie within source clock")
    return np.asarray([
        np.interp(target_time, source_time, row) for row in values
    ])


def merge_inner_outer_traces(
    *, inner_time_s: np.ndarray, inner_labels: np.ndarray,
    inner_traces: np.ndarray, outer_time_s: np.ndarray,
    outer_labels: np.ndarray, outer_traces: np.ndarray,
    sensor_count: int = 145,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Prefer the resolved inner value and retain outer overlap diagnostics."""

    inner_labels_text = [str(value) for value in np.asarray(inner_labels)]
    outer_labels_text = [str(value) for value in np.asarray(outer_labels)]
    if len(set(inner_labels_text)) != len(inner_labels_text) or len(set(outer_labels_text)) != len(outer_labels_text):
        raise ValueError("inner and outer labels must be unique")
    inner = np.asarray(inner_traces, dtype=float)
    outer = _interpolate_rows(outer_time_s, outer_traces, inner_time_s)
    inner_map = dict(zip(inner_labels_text, inner))
    outer_map = dict(zip(outer_labels_text, outer))
    expected = [f"TI021-S{index + 1:03d}" for index in range(sensor_count)]
    union = set(inner_map) | set(outer_map)
    if union != set(expected):
        missing = sorted(set(expected) - union)
        extra = sorted(union - set(expected))
        raise ValueError(f"inner/outer union does not cover 145 sensors; missing={missing}, extra={extra}")
    overlap = sorted(set(inner_map) & set(outer_map))
    merged = np.asarray([
        inner_map[label] if label in inner_map else outer_map[label]
        for label in expected
    ])
    source = np.asarray([
        "inner" if label in inner_map else "outer" for label in expected
    ])
    overlap_inner = np.asarray([inner_map[label] for label in overlap])
    overlap_outer = np.asarray([outer_map[label] for label in overlap])
    return {
        "time_s": np.asarray(inner_time_s, dtype=float),
        "labels": np.asarray(expected),
        "concentration_vol_percent": merged,
        "selected_source": source,
        "overlap_labels": np.asarray(overlap),
        "overlap_inner_vol_percent": overlap_inner,
        "overlap_outer_vol_percent": overlap_outer,
    }, {
        "inner_sensor_count": len(inner_map),
        "outer_sensor_count": len(outer_map),
        "overlap_sensor_count": len(overlap),
        "combined_unique_sensor_count": len(union),
        "selection_policy": "prefer resolved inner where both domains cover a sensor",
    }


def _timing_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, object]:
    available = np.isfinite(observed)
    paired = available & np.isfinite(predicted)
    difference = predicted[paired] - observed[paired]
    return {
        "observed_available_count": int(np.count_nonzero(available)),
        "model_detected_among_available_count": int(np.count_nonzero(paired)),
        "missed_count": int(np.count_nonzero(available & ~np.isfinite(predicted))),
        "mean_error_s": None if difference.size == 0 else float(np.mean(difference)),
        "mean_absolute_error_s": None if difference.size == 0 else float(np.mean(np.abs(difference))),
        "maximum_absolute_error_s": None if difference.size == 0 else float(np.max(np.abs(difference))),
    }


def _trace_metrics(time_s: np.ndarray, traces: np.ndarray) -> dict[str, np.ndarray]:
    filtered = trailing_boxcar_average(time_s, traces, window_s=0.6)
    peak = np.max(filtered, axis=1)
    peak_time = time_s[np.argmax(filtered, axis=1)]
    dose = np.trapezoid(traces, time_s, axis=1)
    arrival = np.full(traces.shape[0], np.nan)
    departure = np.full(traces.shape[0], np.nan)
    for index, trace in enumerate(filtered):
        hit = np.flatnonzero(trace >= 0.1)
        if hit.size:
            arrival[index] = time_s[int(hit[0])]
            departure[index] = time_s[int(hit[-1])]
    return {
        "filtered": filtered, "peak": peak, "peak_time": peak_time,
        "dose": dose, "arrival": arrival, "departure": departure,
    }


def _overlap_geometry_diagnostic(
    inner_record: dict[str, object], outer_record: dict[str, object],
    overlap_labels: np.ndarray,
) -> dict[str, object]:
    """Check that overlap sensors sample shared fluid, not the inlet plane.

    Coincident observations on the coupling plane are useful boundary monitors,
    but they do not independently validate the solution in a common interior
    region. Require three distinct horizontal sensor locations at least one
    coarse cell downstream of the plane before the overlap score is considered
    a spatial coupling check.
    """
    inner_coordinates = inner_record.get("coordinates", {})
    handoff_y = inner_coordinates.get("handoff_site_downwind_y_m")
    if handoff_y is None:
        raise ValueError("inner transient is missing the physical handoff location")
    inner_config = inner_record["configuration"]
    outer_config = outer_record["configuration"]
    buffer_m = max(
        float(inner_config["cell_size_m"]),
        float(outer_config["cell_size_m"]),
    )
    inner_receptors = {
        item["label"]: item
        for item in inner_config["receptors"]
    }
    outer_receptors = {
        item["label"]: item
        for item in outer_config["receptors"]
    }
    mismatched = []
    inner_solver_mismatched = []
    outer_solver_mismatched = []
    coincident = []
    interior = []
    horizontal_positions: set[tuple[float, float]] = set()
    tolerance = 1.0e-8
    source_solver = np.asarray(
        inner_coordinates.get("source_solver_xyz_m", (60.0, 0.5 * float(inner_config["domain_m"][1]), 0.0)),
        dtype=float,
    )
    inner_crosswind_origin = 0.5 * float(inner_config["domain_m"][1])
    declared_inner_crosswind_origin = float(
        inner_config.get("inner_crosswind_origin_m", source_solver[1])
    )
    inner_frame_origin_mismatch = not (
        np.isclose(source_solver[1], inner_crosswind_origin, rtol=0.0, atol=tolerance)
        and np.isclose(
            declared_inner_crosswind_origin, inner_crosswind_origin,
            rtol=0.0, atol=tolerance,
        )
    )
    outer_crosswind_origin = 0.5 * float(outer_config["domain_m"][1])
    for raw_label in overlap_labels:
        label = str(raw_label)
        inner_row = inner_receptors[label]
        outer_row = outer_receptors[label]
        inner_xyz = np.asarray(inner_row["site_local_xyz_m"], dtype=float)
        outer_xyz = np.asarray(outer_row["site_local_xyz_m"], dtype=float)
        if inner_xyz.shape != (3,) or outer_xyz.shape != (3,):
            raise ValueError(f"overlap sensor {label} must have 3-D site coordinates")
        if not np.allclose(inner_xyz, outer_xyz, rtol=0.0, atol=tolerance):
            mismatched.append(label)
            continue
        expected_inner_solver = np.asarray((
            float(source_solver[0]) + float(inner_xyz[1]),
            inner_crosswind_origin - float(inner_xyz[0]),
            float(inner_xyz[2]),
        ))
        expected_outer_solver = np.asarray((
            float(inner_xyz[1]) - float(handoff_y),
            outer_crosswind_origin - float(inner_xyz[0]),
            float(inner_xyz[2]),
        ))
        actual_inner_solver = inner_row.get("solver_xyz_m")
        actual_outer_solver = outer_row.get("outer_solver_xyz_m")
        if inner_frame_origin_mismatch or actual_inner_solver is None or not np.allclose(
            np.asarray(actual_inner_solver, dtype=float), expected_inner_solver,
            rtol=0.0, atol=tolerance,
        ):
            inner_solver_mismatched.append(label)
        if actual_outer_solver is None or not np.allclose(
            np.asarray(actual_outer_solver, dtype=float), expected_outer_solver,
            rtol=0.0, atol=tolerance,
        ):
            outer_solver_mismatched.append(label)
        if abs(float(inner_xyz[1]) - float(handoff_y)) <= tolerance:
            coincident.append(label)
        if float(inner_xyz[1]) > float(handoff_y) + buffer_m + tolerance:
            interior.append(label)
            horizontal_positions.add((round(float(inner_xyz[0]), 8), round(float(inner_xyz[1]), 8)))
    return {
        "handoff_site_downwind_y_m": float(handoff_y),
        "minimum_interior_buffer_m": buffer_m,
        "boundary_coincident_sensor_labels": coincident,
        "interior_overlap_sensor_labels": interior,
        "interior_overlap_sensor_count": len(interior),
        "interior_overlap_horizontal_location_count": len(horizontal_positions),
        "minimum_required_interior_horizontal_locations": 3,
        "coordinate_mismatch_sensor_labels": mismatched,
        "inner_solver_coordinate_mismatch_sensor_labels": inner_solver_mismatched,
        "outer_solver_coordinate_mismatch_sensor_labels": outer_solver_mismatched,
        "inner_crosswind_origin_mismatch": inner_frame_origin_mismatch,
        "all_shared_sensor_coordinates_match": not (
            mismatched or inner_solver_mismatched or outer_solver_mismatched
        ),
        "at_least_three_interior_horizontal_locations": len(horizontal_positions) >= 3,
    }


def _vertical_sensor_resolution_diagnostic(
    sensors: tuple[object, ...], labels: np.ndarray, traces: np.ndarray,
    *, bottom_cell_centre_m: float,
) -> dict[str, object]:
    """Expose receptor pairs collapsed onto the first cell-centred z value."""
    label_index = {str(label): index for index, label in enumerate(labels)}
    by_horizontal_position: dict[tuple[float, float], list[tuple[int, object]]] = {}
    for sensor_index, sensor in enumerate(sensors):
        if sensor.height_m <= bottom_cell_centre_m + 1.0e-12:
            key = (round(float(sensor.local_x_m), 9), round(float(sensor.local_y_m), 9))
            by_horizontal_position.setdefault(key, []).append((sensor_index, sensor))

    pairs: list[dict[str, object]] = []
    positive_observed_dose_ratios: list[float] = []
    for members in by_horizontal_position.values():
        for left_position in range(len(members)):
            for right_position in range(left_position + 1, len(members)):
                left_index, left = members[left_position]
                right_index, right = members[right_position]
                left_label = f"TI021-S{left_index + 1:03d}"
                right_label = f"TI021-S{right_index + 1:03d}"
                left_trace = np.asarray(traces[label_index[left_label]])
                right_trace = np.asarray(traces[label_index[right_label]])
                identical = bool(np.array_equal(left_trace, right_trace))
                left_dose = float(left.dose_percent_s)
                right_dose = float(right.dose_percent_s)
                ratio = None
                if left_dose > 0.0 and right_dose > 0.0:
                    ratio = max(left_dose, right_dose) / min(left_dose, right_dose)
                    positive_observed_dose_ratios.append(ratio)
                pairs.append({
                    "labels": [left_label, right_label],
                    "heights_m": [float(left.height_m), float(right.height_m)],
                    "observed_doses_vol_percent_s": [left_dose, right_dose],
                    "observed_dose_ratio": ratio,
                    "predicted_traces_identical": identical,
                })

    ratios = np.asarray(positive_observed_dose_ratios, dtype=float)
    return {
        "scope": "post-run diagnostic; does not alter configuration or score",
        "inner_grid_vertical_spacing_m": 2.0 * bottom_cell_centre_m,
        "bottom_cell_centre_m": float(bottom_cell_centre_m),
        "receptors_at_or_below_bottom_cell_centre_count": int(sum(
            sensor.height_m <= bottom_cell_centre_m + 1.0e-12
            for sensor in sensors
        )),
        "same_horizontal_position_pairs_on_bottom_cell": len(pairs),
        "pairs_with_identical_predicted_traces": int(sum(
            pair["predicted_traces_identical"] for pair in pairs
        )),
        "pairs_with_positive_observed_dose_at_both_heights": int(ratios.size),
        "positive_observed_dose_ratio_median": (
            None if ratios.size == 0 else float(np.median(ratios))
        ),
        "positive_observed_dose_ratio_maximum": (
            None if ratios.size == 0 else float(np.max(ratios))
        ),
        "pairs": pairs,
        "interpretation": (
            "sensor height variation within the bottom cell is not resolved; "
            "identical predictions must not be interpreted as validated vertical profiles"
        ),
    }


def score_nested_run(
    *, workbook: Path, inner_run_directory: Path, outer_run_directory: Path,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    workbook_path = workbook.resolve()
    inner = inner_run_directory.resolve()
    outer = outer_run_directory.resolve()
    inner_record_path = inner / "ti021-inner-transient.json"
    outer_record_path = outer / "ti021-outer-transient.json"
    inner_trace_path = inner / "ti021-inner-receptor-traces.npz"
    outer_trace_path = outer / "ti021-outer-receptor-traces.npz"
    inner_record = json.loads(inner_record_path.read_text(encoding="utf-8"))
    outer_record = json.loads(outer_record_path.read_text(encoding="utf-8"))
    if inner_record.get("status") not in {
        "accepted_completed_inner_transient",
        "accepted_inner_transient_qualification",
    }:
        raise ValueError("inner transient is not completed or explicitly qualified")
    if outer_record.get("status") != "accepted_completed_outer_transient":
        raise ValueError("outer transient is not completed and accepted")
    if inner_record["configuration"].get("observation_values_used_to_configure_run") != [] or outer_record["configuration"].get("observation_values_used_to_configure_run") != []:
        raise ValueError("nested executions must be observation-blind")
    with np.load(inner_trace_path, allow_pickle=False) as source:
        inner_trace = {key: np.asarray(source[key]).copy() for key in source.files}
    with np.load(outer_trace_path, allow_pickle=False) as source:
        outer_trace = {key: np.asarray(source[key]).copy() for key in source.files}
    merged, coverage = merge_inner_outer_traces(
        inner_time_s=inner_trace["time_s"], inner_labels=inner_trace["labels"],
        inner_traces=inner_trace["concentration_vol_percent"],
        outer_time_s=outer_trace["time_s"], outer_labels=outer_trace["labels"],
        outer_traces=outer_trace["concentration_vol_percent"],
    )
    if coverage["combined_unique_sensor_count"] != 145:
        raise ValueError(f"nested sensor union must cover all 145 receptors: {coverage}")
    time_s = merged["time_s"]
    traces = merged["concentration_vol_percent"]
    metrics = _trace_metrics(time_s, traces)
    overlap_inner_metrics = _trace_metrics(time_s, merged["overlap_inner_vol_percent"])
    overlap_outer_metrics = _trace_metrics(time_s, merged["overlap_outer_vol_percent"])
    overlap_peak_relative = np.abs(
        overlap_outer_metrics["peak"] - overlap_inner_metrics["peak"]
    ) / np.maximum(np.abs(overlap_inner_metrics["peak"]), 0.1)
    overlap_dose_relative = np.abs(
        overlap_outer_metrics["dose"] - overlap_inner_metrics["dose"]
    ) / np.maximum(np.abs(overlap_inner_metrics["dose"]), 0.1)
    overlap_geometry = _overlap_geometry_diagnostic(
        inner_record, outer_record, merged["overlap_labels"],
    )
    label_to_index = {str(label): index for index, label in enumerate(merged["overlap_labels"])}
    interior_indices = np.asarray([
        label_to_index[label]
        for label in overlap_geometry["interior_overlap_sensor_labels"]
    ], dtype=int)
    boundary_indices = np.asarray([
        label_to_index[label]
        for label in overlap_geometry["boundary_coincident_sensor_labels"]
    ], dtype=int)
    interior_peak_relative = overlap_peak_relative[interior_indices]
    interior_dose_relative = overlap_dose_relative[interior_indices]
    boundary_peak_relative = overlap_peak_relative[boundary_indices]
    boundary_dose_relative = overlap_dose_relative[boundary_indices]
    overlap_gates = {
        "at_least_three_interior_overlap_sensors": interior_indices.size >= 3,
        "at_least_three_interior_overlap_locations": overlap_geometry[
            "at_least_three_interior_horizontal_locations"
        ],
        "shared_sensor_coordinates_match": overlap_geometry[
            "all_shared_sensor_coordinates_match"
        ],
        "maximum_interior_overlap_peak_change_at_most_10_percent": bool(
            interior_peak_relative.size and float(np.max(interior_peak_relative)) <= 0.10
        ),
        "maximum_interior_overlap_dose_change_at_most_10_percent": bool(
            interior_dose_relative.size and float(np.max(interior_dose_relative)) <= 0.10
        ),
    }
    dataset = read_thorney_ti021_workbook(workbook_path)
    sensors = dataset.dose_sensors
    if len(sensors) != 145:
        raise ValueError("TI021 workbook must contain exactly 145 dose sensors")
    observed_dose = np.asarray([sensor.dose_percent_s for sensor in sensors])
    observed_arrival = np.asarray([
        np.nan if sensor.arrival_time_s is None else sensor.arrival_time_s
        for sensor in sensors
    ])
    observed_departure = np.asarray([
        np.nan if sensor.departure_time_s is None else sensor.departure_time_s
        for sensor in sensors
    ])
    inner_configuration = inner_record["configuration"]
    inner_horizon_s = float(inner_configuration["duration_s"])
    outer_horizon_s = float(outer_record["configuration"]["duration_s"])
    if not np.isclose(inner_horizon_s, outer_horizon_s):
        raise ValueError("inner and outer transient horizons must match before scoring")
    full_event_horizon = np.isclose(inner_horizon_s, 600.0)
    inner_vertical_spacing_m = (
        float(inner_configuration["domain_m"][2])
        / int(inner_configuration["cells"][2])
    )
    vertical_resolution = _vertical_sensor_resolution_diagnostic(
        sensors, merged["labels"], traces,
        bottom_cell_centre_m=0.5 * inner_vertical_spacing_m,
    )
    dose_score = _ratio_metrics(observed_dose, metrics["dose"])
    score = {
        "schema": SCHEMA,
        "status": (
            f"scored_nested_145_sensors_{'full_event' if full_event_horizon else 'short_horizon_diagnostic'}_"
            f"overlap_{'pass' if all(overlap_gates.values()) else 'failed'}_grid_withheld"
        ),
        "scope": (
            "no-fit inner/outer union of all 145 TI021 dose sensors; model-form "
            "score is independent of configuration; horizon and finest-grid status are reported separately"
        ),
        "provenance": {
            "workbook_path": str(workbook_path), "workbook_sha256": _sha256(workbook_path),
            "inner_record": {"path": str(inner_record_path), "sha256": _sha256(inner_record_path)},
            "inner_trace": {"path": str(inner_trace_path), "sha256": _sha256(inner_trace_path)},
            "outer_record": {"path": str(outer_record_path), "sha256": _sha256(outer_record_path)},
            "outer_trace": {"path": str(outer_trace_path), "sha256": _sha256(outer_trace_path)},
        },
        "protocol": {
            "configuration_accessed_observation_values": [],
            "scoring_opened_after_both_completed_runs": True,
            "concentration_unit": "vol%", "dose_unit": "vol% s",
            "sensor_response_average_s": 0.6,
            "arrival_and_departure_threshold_vol_percent": 0.1,
            "outer_trace_interpolation": "linear onto the finer inner event clock before union",
            "overlap_acceptance_gate": (
                "at least three distinct shared interior horizontal sensor locations "
                "more than one coarse cell from the interface, matching physical "
                "coordinates, and maximum peak and dose change each <=10%"
            ),
        },
        "coverage": coverage,
        "model_horizon_s": inner_horizon_s,
        "overlap_audit": {
            "labels": merged["overlap_labels"].tolist(),
            "interior_labels": overlap_geometry["interior_overlap_sensor_labels"],
            "inner_peak_vol_percent": overlap_inner_metrics["peak"].tolist(),
            "outer_peak_vol_percent": overlap_outer_metrics["peak"].tolist(),
            "peak_absolute_relative_change": overlap_peak_relative.tolist(),
            "inner_dose_vol_percent_s": overlap_inner_metrics["dose"].tolist(),
            "outer_dose_vol_percent_s": overlap_outer_metrics["dose"].tolist(),
            "dose_absolute_relative_change": overlap_dose_relative.tolist(),
            "maximum_peak_absolute_relative_change": (
                float(np.max(interior_peak_relative)) if interior_peak_relative.size else None
            ),
            "maximum_dose_absolute_relative_change": (
                float(np.max(interior_dose_relative)) if interior_dose_relative.size else None
            ),
            "maximum_boundary_monitor_peak_absolute_relative_change": (
                float(np.max(boundary_peak_relative)) if boundary_peak_relative.size else None
            ),
            "maximum_boundary_monitor_dose_absolute_relative_change": (
                float(np.max(boundary_dose_relative)) if boundary_dose_relative.size else None
            ),
            "geometry": overlap_geometry,
            "gates": overlap_gates,
            "all_gates_pass": bool(all(overlap_gates.values())),
        },
        "dose_sensor_score": {
            "sensor_count": 145,
            "observed_positive_count": int(np.count_nonzero(observed_dose > 0.0)),
            "predicted_detected_count_at_0p1vol_percent": int(np.count_nonzero(metrics["peak"] >= 0.1)),
            **dose_score,
        },
        "arrival_score": _timing_metrics(observed_arrival, metrics["arrival"]),
        "departure_score": _timing_metrics(observed_departure, metrics["departure"]),
        "four_vol_percent_sensor_summary": {
            "predicted_sensor_count_at_or_above_4vol_percent": int(np.count_nonzero(metrics["peak"] >= 4.0)),
            "maximum_predicted_0p6s_peak_vol_percent": float(np.max(metrics["peak"])),
            "note": "dose-sheet sensors do not provide 145 independent measured peak concentrations",
        },
        "vertical_sensor_resolution_diagnostic": vertical_resolution,
        "decision": {
            "result": "withheld_as_unconditional_quantitative_validation",
            "overlap_coupling": "accepted" if all(overlap_gates.values()) else "withheld",
            "reason": (
                (
                    "nested-domain coupling failed its independent overlap gate: "
                    f"maximum interior peak change={float(np.max(interior_peak_relative)) if interior_peak_relative.size else float('nan'):.6g}, "
                    f"maximum interior dose change={float(np.max(interior_dose_relative)) if interior_dose_relative.size else float('nan'):.6g}; "
                    f"interior overlap locations={overlap_geometry['interior_overlap_horizontal_location_count']}, "
                    "so the shared sensor set does not establish in-domain overlap coupling; "
                    "full-event finest-grid convergence also remains withheld"
                )
                if not all(overlap_gates.values()) else
                (
                    (
                        "overlap coupling passed, but this is a 60 s diagnostic "
                        "rather than a full-event score; finest-grid convergence "
                        "remains withheld"
                        if not full_event_horizon else
                        "overlap coupling passed, but full-event finest-grid "
                        "convergence remains withheld"
                    )
                )
            ),
        },
    }
    arrays = {
        **merged,
        "filtered_vol_percent": metrics["filtered"],
        "predicted_peak": metrics["peak"], "predicted_peak_time": metrics["peak_time"],
        "predicted_dose": metrics["dose"], "predicted_arrival": metrics["arrival"],
        "predicted_departure": metrics["departure"],
        "observed_dose": observed_dose, "observed_arrival": observed_arrival,
        "observed_departure": observed_departure,
        "overlap_peak_relative": overlap_peak_relative,
        "overlap_dose_relative": overlap_dose_relative,
    }
    return score, arrays


def _plot(score: dict[str, object], arrays: dict[str, np.ndarray], output: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(13.0, 10.0), constrained_layout=True)
    positive = arrays["observed_dose"] > 0.0
    floor = 0.1
    observed = np.maximum(arrays["observed_dose"][positive], floor)
    predicted = np.maximum(arrays["predicted_dose"][positive], floor)
    axes[0, 0].scatter(observed, predicted, s=20, alpha=0.75)
    limit = max(float(np.max(observed)), float(np.max(predicted)), 1.0)
    axes[0, 0].plot((floor, limit), (floor, limit), "k-")
    axes[0, 0].plot((floor, limit), (2 * floor, 2 * limit), "k--", alpha=0.5)
    axes[0, 0].plot((floor, limit), (0.5 * floor, 0.5 * limit), "k--", alpha=0.5)
    axes[0, 0].set(xscale="log", yscale="log", xlabel="Observed dose (vol% s)", ylabel="Predicted dose (vol% s)", title="All 145 sensors")
    interior_labels = score["overlap_audit"]["interior_labels"]
    overlap_indices = [score["overlap_audit"]["labels"].index(label) for label in interior_labels]
    overlap_x = np.arange(len(overlap_indices))
    inner_peaks = _trace_metrics(arrays["time_s"], arrays["overlap_inner_vol_percent"])["peak"][overlap_indices]
    outer_peaks = _trace_metrics(arrays["time_s"], arrays["overlap_outer_vol_percent"])["peak"][overlap_indices]
    axes[0, 1].bar(overlap_x - 0.2, inner_peaks, 0.4, label="inner")
    axes[0, 1].bar(overlap_x + 0.2, outer_peaks, 0.4, label="outer")
    axes[0, 1].set(xticks=overlap_x, xticklabels=interior_labels, ylabel="0.6 s peak (vol%)", title="Interior overlap sensors (boundary monitors excluded)")
    axes[0, 1].tick_params(axis="x", rotation=45); axes[0, 1].legend()
    available = np.isfinite(arrays["observed_arrival"]) & np.isfinite(arrays["predicted_arrival"])
    axes[1, 0].scatter(arrays["observed_arrival"][available], arrays["predicted_arrival"][available], s=20, alpha=0.75)
    if np.any(available):
        timing_limit = max(float(np.max(arrays["observed_arrival"][available])), float(np.max(arrays["predicted_arrival"][available])))
        axes[1, 0].plot((0, timing_limit), (0, timing_limit), "k-")
    axes[1, 0].set(xlabel="Observed arrival (s)", ylabel="Predicted arrival (s)", title="0.1 vol% arrival")
    order = np.argsort(arrays["observed_dose"])[::-1]
    axes[1, 1].plot(arrays["observed_dose"][order], label="observed", color="black")
    axes[1, 1].plot(arrays["predicted_dose"][order], label="predicted", color="tab:orange")
    axes[1, 1].set(yscale="symlog", xlabel="Sensor rank by observed dose", ylabel="Dose (vol% s)", title="Dose distribution")
    axes[1, 1].legend()
    for axis in axes.flat:
        axis.grid(alpha=0.2)
    figure.suptitle(
        f"TI021 nested inner/outer no-fit score — {score['model_horizon_s']:g} s; "
        "full-event grid convergence withheld"
    )
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--inner-run-directory", type=Path, required=True)
    parser.add_argument("--outer-run-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    score, arrays = score_nested_run(
        workbook=args.workbook, inner_run_directory=args.inner_run_directory,
        outer_run_directory=args.outer_run_directory,
    )
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    json_path = output / "ti021-nested-145-sensor-score.json"
    npz_path = output / "ti021-nested-145-sensor-traces.npz"
    plot_path = output / "ti021-nested-145-sensor-score.png"
    json_path.write_text(json.dumps(score, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez_compressed(npz_path, **arrays)
    _plot(score, arrays, plot_path)
    print(json.dumps({
        "status": score["status"], "decision": score["decision"],
        "coverage": score["coverage"], "overlap_audit": score["overlap_audit"],
        "dose_sensor_score": score["dose_sensor_score"],
        "arrival_score": score["arrival_score"],
        "outputs": {"json": str(json_path), "npz": str(npz_path), "plot": str(plot_path)},
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
