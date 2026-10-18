"""Run a reproducible H2 obstacle screening case with vol% and sensor outputs.

This case is a solver/IO qualification for the H2 path.  It is deliberately
not labelled as facility validation: the release and obstacle geometry are
declared inputs, while the external evidence gate is kept separate.
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


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "h2-obstacle-fast-screening-600s-2026-10-14.json"


def run_case(*, output_path: Path, duration_s: float = 600.0) -> dict[str, object]:
    if duration_s <= 0.0:
        raise ValueError("duration_s must be positive")
    release_end_s = min(10.0, 0.5 * duration_s)
    config = FastObstacleTransport3DConfig(
        domain_m=(100.0, 40.0, 20.0),
        cells=(50, 20, 10),
        duration_s=duration_s,
        time_step_s=0.1,
        wind_speed_m_s=3.0,
        diffusivity_m2_s=0.08,
        source_schedule=SourceRateSchedule(
            time_s=(0.0, release_end_s, duration_s),
            rate_kg_s=(0.20, 0.0, 0.0),
            source_id="declared H2 transient screening release",
        ),
        source_position_m=(5.0, 20.0, 1.0),
        source_sigma_m=1.5,
        obstacles=(AxisAlignedCuboid(
            30.0, 40.0, 14.0, 26.0, 0.0, 10.0,
            label="declared facility obstacle",
        ),),
        store_fields=False,
    )
    receptors = (
        FastObstacleReceptor3D("front", (27.0, 20.0, 3.0)),
        FastObstacleReceptor3D("top", (40.0, 20.0, 12.0)),
        FastObstacleReceptor3D("side", (40.0, 10.0, 3.0)),
        FastObstacleReceptor3D("wake", (55.0, 20.0, 3.0)),
    )
    started = time.perf_counter()
    result = solve_fast_obstacle_transport_3d(config, receptors=receptors)
    runtime_s = time.perf_counter() - started
    vol = result.maximum_concentration_vol_percent
    finite_fluid = np.isfinite(vol) & (~result.solid_mask)
    threshold_mask = finite_fluid & (vol >= 4.0)
    cell_volume = float(np.prod(result.diagnostics["cell_size_m"]))
    threshold_count = int(np.count_nonzero(threshold_mask))
    if threshold_count:
        zz, yy, xx = np.where(threshold_mask)
        centroid = [
            float(np.mean(result.x_m[xx])),
            float(np.mean(result.y_m[yy])),
            float(np.mean(result.z_m[zz])),
        ]
    else:
        centroid = None
    summary = {
        "schema": "degali.h2-obstacle-fast-screening.v1",
        "status": "conditional_h2_obstacle_screening",
        "quantitative_prediction_allowed": False,
        "external_validation_gate": "withheld_no_facility_manifest_attached",
        "duration_s": duration_s,
        "runtime_s": runtime_s,
        "runtime_per_simulated_second_s": runtime_s / duration_s,
        "cells": list(config.cells),
        "cell_size_m": list(result.diagnostics["cell_size_m"]),
        "time_step_s": result.diagnostics["time_step_s"],
        "source": {
            "species": "H2",
            "total_released_mass_kg": config.source_schedule.released_mass_kg,
            "schedule_source_id": config.source_schedule.source_id,
        },
        "obstacles": [item.label for item in config.obstacles],
        "maximum_h2_vol_percent": float(np.nanmax(vol[finite_fluid])),
        "four_vol_percent_volume_m3": threshold_count * cell_volume,
        "four_vol_percent_cell_count": threshold_count,
        "four_vol_percent_centroid_m": centroid,
        "receptor_peak_h2_vol_percent": {
            label: float(np.nanmax(trace))
            for label, trace in result.receptor_traces_kg_m3.items()
        },
        "receptor_trace_kg_m3": {
            label: trace.tolist()
            for label, trace in result.receptor_traces_kg_m3.items()
        },
        "mass_conservation": {
            "maximum_absolute_residual_kg": result.maximum_mass_residual_kg,
            "passed_screening_tolerance": result.maximum_mass_residual_kg <= 1.0e-9,
        },
        "physics": result.diagnostics,
        "warnings": [
            "potential-flow carrier omits separated recirculating wake and buoyancy feedback",
            "declared synthetic H2 release and obstacle geometry are not a facility event",
            "external sensor calibration, common clock and FieldEvidenceManifest are required for promotion",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez_compressed(
        output_path.with_name(output_path.stem + "-maximum-field.npz"),
        x_m=result.x_m, y_m=result.y_m, z_m=result.z_m,
        maximum_h2_vol_percent=vol,
        solid_mask=result.solid_mask,
        threshold_mask=threshold_mask,
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--duration-s", type=float, default=600.0)
    args = parser.parse_args()
    print(json.dumps(run_case(output_path=args.output, duration_s=args.duration_s), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
