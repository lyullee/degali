"""Compare predicted and measured obstacle-induced sensor deltas for FFI O1.

The absolute FFI concentration score is intentionally withheld because the
published O1 beam assembly and source/sensor clock operator are unresolved.
This paired diagnostic runs the same fast transport path for the unobstructed
control A and the conditional O1 two-top-crossbeam surrogate, then compares
sensor-wise ``obstacle - control`` changes.  It is useful for separating an
obstacle-response error from a common source/absolute-level error; it is not a
facility validation score.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleReceptor3D,
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.validation.open_channel_h2_obstacle import (
    conditional_ffi_o1_crossbeam_assembly,
    load_open_channel_h2_obstacle_benchmark,
)


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"
DATA_ROOT = ROOT / "reference/hydrogen/open_channel_usn/extracted"
H2_MOLAR_MASS = 2.01588e-3
R_GAS = 8.314462618
T_AMBIENT_K = 293.0
P_AMBIENT_PA = 101325.0


def _window_mean(time_s: np.ndarray, values: np.ndarray, start: float, end: float) -> float:
    keep = (time_s > start) & (time_s < end)
    local_t = np.concatenate(([start], time_s[keep], [end]))
    local_y = np.interp(local_t, time_s, values)
    return float(np.trapezoid(local_y, local_t) / (end - start))


def _source_schedule(run, start: float, end: float) -> SourceRateSchedule:
    keep = (run.flow_time_s >= start) & (run.flow_time_s <= end)
    source_t = np.asarray(run.flow_time_s[keep] - start, dtype=float)
    source_q = np.asarray(run.mass_flow_g_s[keep] * 1.0e-3, dtype=float)
    if source_t.size < 2:
        raise ValueError("FFI source history does not cover the requested window")
    source_t[0] = 0.0
    source_t[-1] = end - start
    if source_q[-1] != 0.0:
        source_t = np.append(source_t, end - start + 1.0e-6)
        source_q = np.append(source_q, 0.0)
    return SourceRateSchedule(tuple(source_t), tuple(source_q), source_id="FFI-paired-delta")


def _conditional_beams() -> tuple[AxisAlignedCuboid, ...]:
    return tuple(
        AxisAlignedCuboid(
            beam.y_min_m, beam.y_max_m,
            beam.x_min_m, beam.x_max_m,
            beam.z_min_m, beam.z_max_m,
            beam.label,
        )
        for beam in conditional_ffi_o1_crossbeam_assembly()
    )


def _run_case(run, *, obstacles, receptors, cell_size: float, start: float, end: float):
    domain = (5.8, 0.9, 0.8)
    cells = tuple(int(round(length / cell_size)) for length in domain)
    config = FastObstacleTransport3DConfig(
        domain_m=domain,
        cells=cells,
        duration_s=end - start + 1.0e-6,
        time_step_s=0.02,
        wind_speed_m_s=0.1,
        diffusivity_m2_s=0.01,
        source_schedule=_source_schedule(run, start, end),
        source_position_m=(0.5, 0.45, 0.76),
        source_sigma_m=0.08,
        obstacles=obstacles,
        store_fields=False,
        open_boundary_faces=("x_max",),
    )
    started = time.perf_counter()
    result = solve_fast_obstacle_transport_3d(config, receptors=receptors)
    runtime = time.perf_counter() - started
    molar_density = P_AMBIENT_PA / (R_GAS * T_AMBIENT_K)
    means = {
        label: float(np.mean(trace[1:]) / H2_MOLAR_MASS / molar_density * 100.0)
        for label, trace in result.receptor_traces_kg_m3.items()
    }
    return result, means, runtime


def build_report(*, registry: Path, data_root: Path, start: float, end: float,
                 cell_size: float) -> dict[str, object]:
    benchmark = load_open_channel_h2_obstacle_benchmark(registry, data_root)
    cases = {case.article_case: case for case in benchmark.cases}
    control, obstacle = cases["A"].runs[0], cases["O1"].runs[0]
    receptors = tuple(
        FastObstacleReceptor3D(
            sensor.sensor_id,
            (sensor.position_m[1], sensor.position_m[0], sensor.position_m[2]),
        )
        for sensor in benchmark.sensor_receptors
    )
    # Sensors intersecting the conditional beam stencil are omitted from both
    # model runs, so the paired delta has identical receptor support.
    beams = _conditional_beams()
    def near_beam(sensor):
        x, y, z = sensor.position_m
        return any(
            beam.x_min_m - cell_size <= x <= beam.x_max_m + cell_size
            and beam.y_min_m - cell_size <= y <= beam.y_max_m + cell_size
            and beam.z_min_m - cell_size <= z <= beam.z_max_m + cell_size
            for beam in beams
        )
    receptors = tuple(item for item in receptors if not near_beam(item))
    control_result, control_pred, control_runtime = _run_case(
        control, obstacles=(), receptors=receptors, cell_size=cell_size,
        start=start, end=end,
    )
    obstacle_result, obstacle_pred, obstacle_runtime = _run_case(
        obstacle, obstacles=beams, receptors=receptors, cell_size=cell_size,
        start=start, end=end,
    )
    predicted_delta = {
        label: obstacle_pred[label] - control_pred[label]
        for label in sorted(set(control_pred) & set(obstacle_pred))
    }
    measured_delta = {}
    for sensor in benchmark.sensor_receptors:
        if sensor.sensor_id not in predicted_delta:
            continue
        c = _window_mean(control.sensor_time_s, control.sensor_percent[sensor.sensor_id], start, end)
        o = _window_mean(obstacle.sensor_time_s, obstacle.sensor_percent[sensor.sensor_id], start, end)
        measured_delta[sensor.sensor_id] = o - c
    labels = sorted(predicted_delta, key=lambda value: int(value.split()[-1]))
    errors = np.asarray([predicted_delta[label] - measured_delta[label] for label in labels])
    measured = np.asarray([measured_delta[label] for label in labels])
    return {
        "schema": "degali.ffi-o1-control-obstacle-delta.v1",
        "status": "conditional_paired_obstacle_effect_diagnostic",
        "window_s": [start, end],
        "cases": {"control": "A", "obstacle": "O1"},
        "geometry_mode": "conditional_two_top_crossbeam_surrogate",
        "quantitative_geometry_allowed": False,
        "cell_size_m": cell_size,
        "cells": list(tuple(int(round(length / cell_size)) for length in (5.8, 0.9, 0.8))),
        "runtime_s": {"control": control_runtime, "obstacle": obstacle_runtime,
                      "total": control_runtime + obstacle_runtime},
        "sensor_count": len(labels),
        "predicted_sensor_mean_vol_percent": {"control": control_pred, "obstacle": obstacle_pred},
        "predicted_obstacle_minus_control_vol_percent": predicted_delta,
        "measured_obstacle_minus_control_vol_percent_reference": measured_delta,
        "diagnostic_error": {
            "mae_vol_percent_points": float(np.mean(np.abs(errors))),
            "rmse_vol_percent_points": float(np.sqrt(np.mean(errors ** 2))),
            "bias_vol_percent_points": float(np.mean(errors)),
            "signed_correlation": float(np.corrcoef(errors + measured, measured)[0, 1])
            if len(labels) > 1 and np.std(measured) > 0.0 else None,
        },
        "conservation": {
            "control_mass_residual_kg": control_result.maximum_mass_residual_kg,
            "obstacle_mass_residual_kg": obstacle_result.maximum_mass_residual_kg,
        },
        "gate_codes": list(benchmark.gate_codes) + [
            "conditional_beam_geometry",
            "source_sensor_clock_alignment_unverified",
        ],
        "promotion_allowed": False,
        "interpretation": (
            "Paired obstacle-effect diagnostic only. It separates a common absolute "
            "source/level error from the conditional obstacle response, but does not "
            "validate the published beam assembly or promote H2 prediction."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cell-size-m", type=float, default=0.05)
    parser.add_argument("--start", type=float, default=60.0)
    parser.add_argument("--end", type=float, default=120.0)
    args = parser.parse_args()
    report = build_report(
        registry=REGISTRY, data_root=DATA_ROOT, start=args.start,
        end=args.end, cell_size=args.cell_size_m,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
