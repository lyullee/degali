"""TI021 fast-3-D obstacle witness with the workbook sensor operator.

This is a conservative scalar R12 witness, not a promoted facility model.  It
exists to test whether a resolved obstacle route changes the 600 s dose shape
relative to the single-layer hybrid without fitting sensor gains.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleReceptor3D, FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import CircularArcWall
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

ROOT = Path(__file__).resolve().parents[1]
R_GAS = 8.31446261815324
R12_MOLAR_MASS = 0.120913


def run(output: Path, *, cell_size_m: float = 5.0, vertical_cell_size_m: float = 2.0,
        duration_s: float = 600.0, diffusivity_m2_s: float = 1.0,
        source_sigma_m: float = 7.0) -> dict:
    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    source_xy = np.array((60.0, 320.0))
    points = np.array([
        (60.0 + item.local_y_m, 320.0 - item.local_x_m, item.height_m)
        for item in dataset.dose_sensors
    ], dtype=float)
    anchors = np.vstack((
        points[:, :2],
        source_xy - np.array((60.0, 60.0)),
        source_xy + np.array((60.0, 60.0)),
    ))
    lower = np.floor((np.min(anchors, axis=0) - 5.0) / cell_size_m) * cell_size_m
    upper = np.ceil((np.max(anchors, axis=0) + 5.0) / cell_size_m) * cell_size_m
    lengths = upper - lower
    nx, ny = (int(round(value / cell_size_m)) for value in lengths)
    nz = int(round(20.0 / vertical_cell_size_m))
    domain = (float(lengths[0]), float(lengths[1]), float(nz * vertical_cell_size_m))
    source_local = tuple(float(value) for value in (*source_xy - lower, 6.5))
    arc = CircularArcWall(
        center_m=source_local[:2], radius_m=dataset.workbook_wall_radius_m,
        thickness_m=0.0, z_min_m=0.0, z_max_m=dataset.workbook_wall_height_m,
        arc_start_deg=0.0, arc_end_deg=180.0, label="TI021-180deg-arc",
    )
    closure = BinaryIdealGasMixtureClosure(
        "R12", R12_MOLAR_MASS, 606.0749,
        ambient_temperature_k=291.15, ambient_pressure_pa=100500.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        (60.0, 320.0), 14.0, 0.0, 13.0,
        float(dataset.initial_r12_mole_fraction), 291.15,
    )
    source_state = closure.state_from_mole_fraction_and_temperature(
        cloud.species_mole_fraction, cloud.temperature_k,
    )
    source_total_mass = cloud.exact_volume_m3 * float(source_state.mixture_density_kg_m3)
    source_species_mass = source_total_mass * float(source_state.species_mass_fraction)
    release_end = min(1.0, 0.5 * duration_s)
    schedule = SourceRateSchedule(
        time_s=(0.0, release_end, duration_s),
        rate_kg_s=(source_species_mass / max(release_end, 1.0e-30), 0.0, 0.0),
        source_id="TI021-R12-1s-release",
    )
    config = FastObstacleTransport3DConfig(
        domain_m=domain, cells=(nx, ny, nz), duration_s=duration_s,
        time_step_s=0.1, wind_speed_m_s=3.4, diffusivity_m2_s=diffusivity_m2_s,
        source_schedule=schedule, source_position_m=source_local,
        source_sigma_m=source_sigma_m, source_sigma_xyz_m=(source_sigma_m, source_sigma_m, 3.75),
        source_kernel="gaussian", obstacles=(arc,), store_fields=False,
    )
    receptors = tuple(
        FastObstacleReceptor3D(f"TI021-S{i + 1:03d}", tuple(point - np.array((*lower, 0.0))))
        for i, point in enumerate(points)
    )
    started = time.perf_counter()
    result = solve_fast_obstacle_transport_3d(config, receptors=receptors)
    runtime = time.perf_counter() - started
    molar_density = 100500.0 / (R_GAS * 291.15)
    conversion = 100.0 / (R12_MOLAR_MASS * molar_density)
    dose = np.array([
        np.trapezoid(result.receptor_traces_kg_m3[item.label], result.time_s) * conversion
        for item in receptors
    ])
    observed = np.array([item.dose_percent_s for item in dataset.dose_sensors])
    positive = observed > 0.0
    report = {
        "schema": "degali.ti021-fast-3d-witness.v1",
        "status": "conditional_fast_3d_scalar_witness",
        "quantitative_prediction_allowed": False,
        "duration_s": duration_s,
        "dose_mape_percent": float(np.mean(np.abs(dose[positive] / observed[positive] - 1.0)) * 100.0),
        "dose_mae_vol_percent_s": float(np.mean(np.abs(dose[positive] - observed[positive]))),
        "normalized_mean_bias": float(np.sum(dose[positive] - observed[positive]) / np.sum(observed[positive])),
        "positive_observed_sensor_count": int(np.sum(positive)),
        "positive_observed_dose_sum_vol_percent_s": float(np.sum(observed[positive])),
        "positive_predicted_dose_sum_vol_percent_s": float(np.sum(dose[positive])),
        "runtime_s": runtime,
        "cells": list(config.cells),
        "cell_size_m": [cell_size_m, cell_size_m, vertical_cell_size_m],
        "source_total_mass_kg": source_total_mass,
        "source_species_mass_kg": source_species_mass,
        "source_sigma_xyz_m": list(config.source_sigma_xyz_m),
        "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        "embedded_surface_blocked_face_count": result.diagnostics["embedded_surface_blocked_face_count"],
        "sensor_rows": [
            {"label": item.label, "observed_dose_vol_percent_s": float(observed[i]),
             "predicted_dose_vol_percent_s": float(dose[i])}
            for i, item in enumerate(receptors)
        ],
        "diagnostics": result.diagnostics,
        "warnings": [
            "passive scalar R12, not a binary variable-density mixture",
            "potential-flow carrier omits separated recirculating wake and buoyancy feedback",
            "witness is not promoted as TI021 validation",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "sensor_rows"}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cell-size-m", type=float, default=5.0)
    parser.add_argument("--vertical-cell-size-m", type=float, default=2.0)
    parser.add_argument("--duration-s", type=float, default=600.0)
    parser.add_argument("--diffusivity-m2-s", type=float, default=1.0)
    parser.add_argument("--source-sigma-m", type=float, default=7.0)
    args = parser.parse_args()
    run(args.output, cell_size_m=args.cell_size_m,
        vertical_cell_size_m=args.vertical_cell_size_m, duration_s=args.duration_s,
        diffusivity_m2_s=args.diffusivity_m2_s, source_sigma_m=args.source_sigma_m)
