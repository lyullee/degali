"""Run a conditional FFI O1 obstacle replay through the fast 3-D solver.

The default geometry is a two-top-crossbeam surrogate that preserves the
published 36.8% frontal blockage ratio.  ``solid_proxy`` remains available as
the historical rectangular sensitivity case.  Neither is machine-registered
facility CAD, so sensor histories and mass-balance evidence are emitted while
the facility accuracy score remains withheld.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleReceptor3D, FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.validation.open_channel_h2_obstacle import (
    conditional_ffi_o1_crossbeam_assembly,
    conditional_ffi_o1_thin_top_crossbeam_assembly,
    load_open_channel_h2_obstacle_benchmark,
)


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"
DATA_ROOT = ROOT / "reference/hydrogen/open_channel_usn"
H2_MOLAR_MASS = 2.01588e-3
R_GAS = 8.314462618
AMBIENT_TEMPERATURE_K = 293.0
AMBIENT_PRESSURE_PA = 101325.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "outputs/ffi-o1-conditional-replay.json")
    parser.add_argument("--cell-size-m", type=float, default=0.1)
    parser.add_argument("--proxy-width-m", type=float, default=0.6)
    parser.add_argument("--wind-speed-m-s", type=float, default=0.1)
    parser.add_argument("--diffusivity-m2-s", type=float, default=0.01)
    parser.add_argument("--source-sigma-m", type=float, default=0.08)
    parser.add_argument("--thin-beam-height-m", type=float, default=0.15)
    parser.add_argument(
        "--geometry-mode", choices=(
            "solid_proxy", "crossbeam_assembly", "thin_top_crossbeam",
        ),
        default="crossbeam_assembly",
    )
    args = parser.parse_args()
    benchmark = load_open_channel_h2_obstacle_benchmark(REGISTRY, DATA_ROOT)
    case = next(item for item in benchmark.cases if item.article_case == "O1")
    run = case.runs[0]
    t0, t1 = 60.0, 120.0
    inside = (run.flow_time_s >= t0) & (run.flow_time_s <= t1)
    source_t = np.asarray(run.flow_time_s[inside] - t0, dtype=float)
    source_q = np.asarray(run.mass_flow_g_s[inside] * 1.0e-3, dtype=float)
    if source_t.size < 2:
        raise ValueError("O1 source history does not cover the replay window")
    source_t[0] = 0.0
    source_t[-1] = t1 - t0
    if source_q[-1] != 0.0:
        source_t = np.append(source_t, t1 - t0 + 1.0e-6)
        source_q = np.append(source_q, 0.0)
    size = float(args.cell_size_m)
    domain = (5.8, 0.9, 0.8)  # solver x/y/z = FFI length/width/height
    cells = tuple(int(round(length / size)) for length in domain)
    if any(abs(length / count - size) > 1e-12 for length, count in zip(domain, cells)):
        raise ValueError("cell size must divide all FFI dimensions")
    all_sensors = tuple(FastObstacleReceptor3D(
        sensor.sensor_id,
        (sensor.position_m[1], sensor.position_m[0], sensor.position_m[2]),
    ) for sensor in benchmark.sensor_receptors)
    # Conditional proxy: reported location/length and an assumed 0.6 m beam
    # projection. The registry explicitly rejects this as quantitative geometry.
    proxy_width = float(args.proxy_width_m)
    if not 0.0 < proxy_width < domain[1]:
        raise ValueError("proxy-width-m must lie inside the 0.9 m channel width")
    if args.geometry_mode in {"crossbeam_assembly", "thin_top_crossbeam"}:
        # Registry coordinates are (channel width x, length y, height z),
        # while the solver uses (streamwise x, lateral y, vertical z).
        channel_beams = (
            conditional_ffi_o1_thin_top_crossbeam_assembly(
                beam_height_m=float(args.thin_beam_height_m),
            )
            if args.geometry_mode == "thin_top_crossbeam"
            else conditional_ffi_o1_crossbeam_assembly()
        )
        obstacles = tuple(
            AxisAlignedCuboid(
                beam.y_min_m, beam.y_max_m,
                beam.x_min_m, beam.x_max_m,
                beam.z_min_m, beam.z_max_m,
                beam.label,
            )
            for beam in channel_beams
        )
        geometry_assumption = (
            "two top crossbeam surrogate preserving reported O1 frontal blockage "
            "ratio 0.368; section and datum remain conditional"
            if args.geometry_mode == "crossbeam_assembly" else
            "thin full-span upper-beam hypothesis inferred from reported 0.368 "
            f"blockage with beam height {args.thin_beam_height_m:.3f} m; "
            "section and datum remain conditional"
        )
    else:
        obstacles = (AxisAlignedCuboid(
            2.1, 3.1,
            0.45 - proxy_width / 2.0,
            0.45 + proxy_width / 2.0,
            0.0, 0.72,
            "FFI O1 conditional solid rectangular proxy",
        ),)
        geometry_assumption = "solid rectangular proxy for unresolved beam assembly"
    sensors = tuple(item for item in all_sensors if not (
        any(
            obstacle.x_min_m - size <= item.position_m[0] <= obstacle.x_max_m + size
            and obstacle.y_min_m - size <= item.position_m[1] <= obstacle.y_max_m + size
            and obstacle.z_min_m - size <= item.position_m[2] <= obstacle.z_max_m + size
            for obstacle in obstacles
        )
    ))
    config = FastObstacleTransport3DConfig(
        domain_m=domain, cells=cells, duration_s=60.000001, time_step_s=0.02,
        # The fast solver requires a positive carrier speed; this conditional
        # proxy uses the channel's weak mean transport and does not claim the
        # article's unresolved expanded-end boundary.
        wind_speed_m_s=float(args.wind_speed_m_s),
        diffusivity_m2_s=float(args.diffusivity_m2_s),
        source_schedule=SourceRateSchedule(tuple(source_t), tuple(source_q),
                                            source_id="FFI-O1-measured-60-120s"),
        source_position_m=(0.5, 0.45, 0.76), source_sigma_m=float(args.source_sigma_m),
        obstacles=obstacles, store_fields=False,
        open_boundary_faces=("x_max",),
    )
    started = time.perf_counter()
    result = solve_fast_obstacle_transport_3d(config, receptors=sensors)
    runtime = time.perf_counter() - started
    def window_mean(sensor_id: str) -> float:
        sensor_time = np.asarray(run.sensor_time_s, dtype=float)
        sensor_values = np.asarray(run.sensor_percent[sensor_id], dtype=float)
        keep = (sensor_time > t0) & (sensor_time < t1)
        local_t = np.concatenate(([t0], sensor_time[keep], [t1]))
        local_y = np.interp(local_t, sensor_time, sensor_values)
        return float(np.trapezoid(local_y, local_t) / (t1 - t0))
    measured = {sensor.sensor_id: window_mean(sensor.sensor_id)
                for sensor in benchmark.sensor_receptors}
    predicted = {label: float(np.mean(values[1:]))
                 for label, values in result.receptor_traces_kg_m3.items()}
    molar_density = AMBIENT_PRESSURE_PA / (R_GAS * AMBIENT_TEMPERATURE_K)
    predicted_vol = {
        label: value / H2_MOLAR_MASS / molar_density * 100.0
        for label, value in predicted.items()
    }
    def exceedance_duration(times: np.ndarray, values: np.ndarray, threshold: float = 4.0) -> float:
        if len(times) < 2:
            return 0.0
        dt = np.diff(times)
        return float(np.sum(dt[np.asarray(values[:-1]) >= threshold]))
    predicted_exceedance_s = {
        label: exceedance_duration(np.asarray(result.time_s), np.asarray(trace))
        for label, trace in (
            (label, np.asarray(values) / H2_MOLAR_MASS / molar_density * 100.0)
            for label, values in result.receptor_traces_kg_m3.items()
        )
    }
    measured_exceedance_s = {}
    for sensor in benchmark.sensor_receptors:
        sensor_time = np.asarray(run.sensor_time_s, dtype=float)
        sensor_values = np.asarray(run.sensor_percent[sensor.sensor_id], dtype=float)
        keep = (sensor_time >= t0) & (sensor_time <= t1)
        measured_exceedance_s[sensor.sensor_id] = exceedance_duration(
            sensor_time[keep] - t0, sensor_values[keep]
        )
    paired = [(predicted_vol[label], measured[label])
              for label in predicted_vol if label in measured]
    errors = np.asarray([p - o for p, o in paired], dtype=float)
    observed = np.asarray([o for _, o in paired], dtype=float)
    diagnostic_error = {
        "sensor_count": int(len(paired)),
        "mae_vol_percent": float(np.mean(np.abs(errors))),
        "rmse_vol_percent": float(np.sqrt(np.mean(errors ** 2))),
        "bias_vol_percent": float(np.mean(errors)),
        "mape_percent_floor_0p1_vol_percent": float(
            np.mean(np.abs(errors) / np.maximum(np.abs(observed), 0.1)) * 100.0
        ),
    }
    for cutoff in (1.0, 4.0):
        keep = np.abs(observed) >= cutoff
        if np.any(keep):
            selected = errors[keep]
            diagnostic_error[f"observed_ge_{cutoff:g}_vol_percent"] = {
                "sensor_count": int(np.count_nonzero(keep)),
                "mae_vol_percent": float(np.mean(np.abs(selected))),
                "rmse_vol_percent": float(np.sqrt(np.mean(selected ** 2))),
                "bias_vol_percent": float(np.mean(selected)),
            }
    report = {
        "schema": "degali.ffi-o1-conditional-replay.v1",
        "status": "conditional_geometry_replay_not_validation_score",
        "case": "O1",
        "window_s": [t0, t1],
        "geometry_mode": args.geometry_mode,
        "geometry_assumption": geometry_assumption,
        "proxy_width_m": proxy_width,
        "wind_speed_m_s": float(args.wind_speed_m_s),
        "diffusivity_m2_s": float(args.diffusivity_m2_s),
        "source_sigma_m": float(args.source_sigma_m),
        "reported_blockage_ratio": 0.368,
        "quantitative_geometry_allowed": False,
        "cells": list(cells), "cell_size_m": size,
        "open_boundary_faces": ["x_max"],
        "runtime_s": runtime,
        "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        "maximum_predicted_vol_percent": float(np.nanmax(result.maximum_concentration_vol_percent)),
        "sensor_count": len(sensors),
        "sensor_ids_omitted_inside_or_near_proxy_stencil": sorted(
            set(item.label for item in all_sensors) - set(predicted)
        ),
        "predicted_sensor_mean_kg_m3": predicted,
        "predicted_sensor_mean_vol_percent": predicted_vol,
        "predicted_sensor_time_s": [float(value) for value in result.time_s],
        "predicted_sensor_trace_vol_percent": {
            label: [float(value / H2_MOLAR_MASS / molar_density * 100.0)
                    for value in trace]
            for label, trace in result.receptor_traces_kg_m3.items()
        },
        "predicted_4vol_exceedance_duration_s": predicted_exceedance_s,
        "measured_4vol_exceedance_duration_s_reference": measured_exceedance_s,
        "measured_sensor_mean_vol_percent_reference": measured,
        "diagnostic_sensor_error_not_validation_score": diagnostic_error,
        "gate_codes": list(benchmark.gate_codes),
        "promotion_allowed": False,
        "warnings": [
            "conditional obstacle surrogate is not machine-registered published CAD",
            "source history is replayed on elapsed time but event-clock alignment remains unverified",
            "predicted sensor values are kg/m3; no facility score is computed",
            "error strata are descriptive only because geometry and clock gates remain withheld",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
