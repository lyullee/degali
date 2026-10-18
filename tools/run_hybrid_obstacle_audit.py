"""Run a small reproducible coarse/local obstacle audit.

This is an execution and bookkeeping harness, not a claim that two overlapping
domains are already a validated facility prediction.  It emits explicit
conditional status when the pair is not suitable for quantitative use.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleTransport3DConfig, solve_fast_obstacle_transport_3d)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.addons.hybrid_field_audit import (
    audit_partitioned_solver_results_with_plan, infer_net_internal_flux_rate,
    threshold_region_metrics)
from degali.addons.hybrid_amr import make_obstacle_patch_bounds
from degali.addons.hybrid_amr_plan import HybridAMRRunPlan
from degali.addons.hybrid_amr import overlap_mass_remap_to_patch
from degali.addons.hydrogen_air_thermodynamics import HydrogenAirActiveScalarClosure

def _run(cells, dt, duration):
    config = FastObstacleTransport3DConfig(
        domain_m=(4.0, 2.4, 2.0), cells=cells, duration_s=duration,
        time_step_s=dt, wind_speed_m_s=2.0, diffusivity_m2_s=0.03,
        source_schedule=SourceRateSchedule((0.0, 0.05, duration), (0.02, 0.0, 0.0)),
        source_position_m=(0.5, 1.2, 0.4), source_sigma_m=0.15,
        obstacles=(AxisAlignedCuboid(1.5, 2.0, 0.8, 1.6, 0.0, 1.0, "audit block"),),
        store_fields=False,
    )
    started = time.perf_counter()
    result = solve_fast_obstacle_transport_3d(config)
    return result, time.perf_counter() - started

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("outputs/hybrid_obstacle_audit.json"))
    args = parser.parse_args()
    duration = 0.2
    coarse, coarse_wall = _run((16, 10, 8), 0.01, duration)
    patch, patch_wall = _run((32, 20, 16), 0.0025, duration)
    patch_bounds = make_obstacle_patch_bounds(
        (1.5, 2.0, 0.8, 1.6, 0.0, 1.0), upstream_m=0.5,
        downstream_m=0.5, lateral_m=0.25, vertical_m=0.25,
        domain_m=(4.0, 2.4, 2.0))
    plan = HybridAMRRunPlan(
        domain_m=(4.0, 2.4, 2.0), coarse_cell_size_m=0.25,
        patch_cell_size_m=0.125, patch=patch_bounds, duration_s=duration,
        coarse_time_step_s=0.01, patch_time_step_s=0.0025)
    audit = audit_partitioned_solver_results_with_plan(
        coarse, patch, plan, mass_tolerance_kg=1e-8)
    coarse_spacing = tuple(float(v) for v in (
        coarse.x_m[1] - coarse.x_m[0], coarse.y_m[1] - coarse.y_m[0],
        coarse.z_m[1] - coarse.z_m[0]))
    coarse_mass_field = np.nan_to_num(coarse.concentration_kg_m3, nan=0.0) * np.prod(coarse_spacing)
    _, overlap_mass = overlap_mass_remap_to_patch(
        coarse_mass_field, x_m=coarse.x_m, y_m=coarse.y_m, z_m=coarse.z_m,
        patch=patch_bounds, patch_shape=(32, 20, 16))
    coarse_internal_flux = infer_net_internal_flux_rate(coarse)
    patch_internal_flux = infer_net_internal_flux_rate(patch)
    concentration = patch.maximum_concentration_kg_m3
    closure = HydrogenAirActiveScalarClosure()
    thermo = closure.state_from_conserved(
        np.nan_to_num(concentration, nan=0.0),
        np.zeros_like(concentration))
    patch_metrics = threshold_region_metrics(
        100.0 * thermo.hydrogen_mole_fraction,
        threshold_vol_percent=4.0,
        cell_volume_m3=float(np.prod(tuple(float(v) for v in (
            patch.x_m[1] - patch.x_m[0], patch.y_m[1] - patch.y_m[0],
            patch.z_m[1] - patch.z_m[0])))),
        coordinates_m=(patch.z_m, patch.y_m, patch.x_m),
    )
    report = {
        "schema": "degali.hybrid-obstacle-audit.v1",
        "status": "conditional",
        "quantitative_ready": bool(audit.quantitative_ready),
        "coarse_runtime_s": coarse_wall,
        "patch_runtime_s": patch_wall,
        "coarse_cells": [16, 10, 8], "patch_cells": [32, 20, 16],
        "duration_s": duration, "patch_substeps_per_coarse_step": 4,
        "partition_weights": {"coarse": plan.partition_weights[0],
                               "patch": plan.partition_weights[1]},
        "maximum_total_mass_drift_kg": audit.maximum_total_mass_drift_kg,
        "coarse_to_patch_overlap_mass_kg": overlap_mass,
        "inferred_internal_flux_rate_kg_s": {
            "coarse_max_abs": float(np.max(np.abs(coarse_internal_flux))),
            "patch_max_abs": float(np.max(np.abs(patch_internal_flux))),
        },
        "threshold": {"basis": "hydrogen mole fraction (vol%)",
                       **patch_metrics},
        "warnings": [
            "coarse and patch fields are independent diagnostic runs",
            "overlap remap is diagnostic only; dynamic interface flux coupling is not claimed",
            "4 vol% conversion uses ambient-pressure/zero-enthalpy closure in this audit",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
