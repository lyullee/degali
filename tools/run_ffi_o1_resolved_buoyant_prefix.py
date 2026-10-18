"""Run a resolved, composition-buoyant FFI O1 sensitivity replay.

The public O1 package does not contain machine-registered beam CAD.  This
runner therefore uses a separately labelled lower-blockage beam placement that
leaves the resolved solver's side/top fetch intact.  It is a physics-path
experiment, not a facility score.  The purpose is to test the missing
quiescent H2/air density-buoyancy and wake path against the same sensor bundle
without promoting a surrogate geometry.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    ResolvedSourceMomentumHistory3D,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.validation.open_channel_h2_obstacle import (
    conditional_ffi_o1_crossbeam_assembly,
    load_registered_ffi_obstacle_geometry,
    load_open_channel_h2_obstacle_benchmark,
)
from degali.validation.field_manifest_runtime import verify_runtime_manifest


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"
DATA_ROOT = ROOT / "reference/hydrogen/open_channel_usn"
H2_MOLAR_MASS = 2.01588e-3
R_GAS = 8.314462618
FFI_INLET_DIAMETER_M = 4.6e-3


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "outputs/ffi-o1-resolved-buoyant-prefix-2026-10-14.json",
    )
    parser.add_argument("--duration-s", type=float, default=60.0)
    parser.add_argument("--cell-size-m", type=float, default=0.1)
    parser.add_argument("--edge-margin-m", type=float, default=0.21)
    parser.add_argument("--wind-speed-m-s", type=float, default=0.1)
    parser.add_argument("--diffusivity-m2-s", type=float, default=0.01)
    parser.add_argument(
        "--source-sigma-m", type=float, default=0.08,
        help="Gaussian source spread used by the resolved sensitivity (m)",
    )
    parser.add_argument(
        "--source-momentum-z-m-s", type=float, default=0.0,
        help=(
            "downward source-jet momentum velocity (negative z); zero keeps "
            "the legacy scalar-only sensitivity path"
        ),
    )
    parser.add_argument(
        "--maximum-lattice-mach-number", type=float, default=0.15,
        help="strict active-flow Mach gate; relaxing it is sensitivity-only",
    )
    parser.add_argument(
        "--subgrid-model", choices=("none", "smagorinsky", "wale"),
        default="smagorinsky",
        help="resolved-flow SGS closure; C=0 preserves the prior molecular-only path",
    )
    parser.add_argument(
        "--smagorinsky-constant", type=float, default=0.0,
        help="Smagorinsky C_s (0 disables eddy viscosity)",
    )
    parser.add_argument(
        "--wale-constant", type=float, default=0.5,
        help="WALE constant when --subgrid-model=wale",
    )
    parser.add_argument(
        "--effective-schmidt-number", type=float, default=None,
        help="optional scalar SGS Schmidt number; requires an explicit scalar closure",
    )
    parser.add_argument(
        "--effective-prandtl-number", type=float, default=None,
        help="optional thermal SGS Prandtl number; requires thermal buoyancy",
    )
    parser.add_argument(
        "--scalar-turbulent-diffusivity-model",
        choices=("legacy_total_kinematic_viscosity_ratio", "molecular_plus_subgrid_eddy_viscosity"),
        default="legacy_total_kinematic_viscosity_ratio",
    )
    parser.add_argument(
        "--thermal-turbulent-diffusivity-model",
        choices=("legacy_total_kinematic_viscosity_ratio", "molecular_plus_subgrid_eddy_viscosity"),
        default="legacy_total_kinematic_viscosity_ratio",
    )
    parser.add_argument(
        "--flow-boundary-scheme",
        choices=(
            "equilibrium_copy", "nonequilibrium_extrapolation",
            "characteristic_nonreflecting", "characteristic_multidimensional",
            "characteristic_multidimensional_extrapolated",
            "characteristic_multidimensional_vortical_limited",
            "characteristic_nscbc_convective",
        ),
        default="equilibrium_copy",
    )
    parser.add_argument(
        "--nonstreamwise-wall-boundary",
        choices=(
            "ground_no_slip_lateral_top_free_slip", "all_no_slip_tunnel",
            "all_symmetry_free_slip", "all_pressure_open",
        ),
        default="ground_no_slip_lateral_top_free_slip",
    )
    parser.add_argument("--flow-steps", type=int, default=600)
    parser.add_argument("--field-output", type=Path, default=None)
    parser.add_argument(
        "--obstacle-geometry-json", type=Path, default=None,
        help=(
            "machine-registered exact O1 cuboid assembly; when omitted, use "
            "the explicitly conditional two-beam surrogate"
        ),
    )
    parser.add_argument("--field-evidence-manifest", type=Path, default=None)
    parser.add_argument(
        "--require-quantitative-manifest", action="store_true",
        help="fail before solving unless the attached manifest clock/calibration gate is accepted",
    )
    args = parser.parse_args()

    field_manifest = None
    if args.field_evidence_manifest is not None:
        field_manifest = verify_runtime_manifest(
            args.field_evidence_manifest,
            require_quantitative=bool(args.require_quantitative_manifest),
        )
    elif args.require_quantitative_manifest:
        raise ValueError(
            "--require-quantitative-manifest requires --field-evidence-manifest"
        )

    if args.duration_s <= 0.0:
        raise ValueError("duration-s must be positive")
    benchmark = load_open_channel_h2_obstacle_benchmark(REGISTRY, DATA_ROOT)
    case = next(item for item in benchmark.cases if item.article_case == "O1")
    run = case.runs[0]
    t0 = 60.0
    t1 = t0 + args.duration_s
    inside = (run.flow_time_s >= t0) & (run.flow_time_s <= t1)
    source_t = np.asarray(run.flow_time_s[inside] - t0, dtype=float)
    source_q = np.asarray(run.mass_flow_g_s[inside] * 1.0e-3, dtype=float)
    if source_t.size < 2:
        raise ValueError("O1 source history does not cover the requested prefix")
    source_t[0] = 0.0
    source_t[-1] = args.duration_s
    if source_q[-1] != 0.0:
        # Keep the measured terminal rate just inside the replay window, then
        # close the schedule with an explicit zero endpoint.
        source_t[-1] = max(0.0, args.duration_s - 1.0e-6)
        source_t = np.append(source_t, args.duration_s)
        source_q = np.append(source_q, 0.0)
    measured_source_q = source_q[source_q > 0.0]
    inlet_area_m2 = np.pi * (FFI_INLET_DIAMETER_M / 2.0) ** 2
    hydrogen_density_kg_m3 = (
        101_325.0 * H2_MOLAR_MASS / (R_GAS * 293.15)
    )
    inferred_orifice_velocity_m_s = (
        measured_source_q / (hydrogen_density_kg_m3 * inlet_area_m2)
        if measured_source_q.size else np.asarray([], dtype=float)
    )
    source_momentum_history = None
    if abs(float(args.source_momentum_z_m_s)) > 0.0:
        peak_source_q = max(float(np.max(measured_source_q)), 1.0e-12)
        source_momentum_history = ResolvedSourceMomentumHistory3D(
            time_s=tuple(float(value) for value in source_t),
            u_m_s=tuple(0.0 for _ in source_t),
            v_m_s=tuple(0.0 for _ in source_t),
            w_m_s=tuple(
                float(args.source_momentum_z_m_s)
                * max(float(rate), 0.0) / peak_source_q
                for rate in source_q
            ),
            history_id="FFI-O1-source-flow-scaled-jet-momentum",
        )

    size = float(args.cell_size_m)
    domain = (5.8, 0.9, 0.8)
    cells = tuple(int(round(length / size)) for length in domain)
    if any(abs(length / count - size) > 1.0e-12 for length, count in zip(domain, cells)):
        raise ValueError("cell size must divide all FFI dimensions")
    physical_step_s = 0.03 * size / float(args.wind_speed_m_s)
    if abs(args.duration_s / physical_step_s - round(args.duration_s / physical_step_s)) > 1.0e-10:
        raise ValueError(
            f"duration-s must be a multiple of the D3Q19 physical step {physical_step_s:g} s"
        )
    if args.obstacle_geometry_json is None:
        channel_beams = conditional_ffi_o1_crossbeam_assembly(
            edge_margin_m=float(args.edge_margin_m),
        )
        geometry_mode = "lower_blockage_resolved_sensitivity"
        geometry_assumption = (
            "two top crossbeam surrogate with edge margin %.3f m; "
            "not the machine-registered 36.8%% facility assembly"
        ) % args.edge_margin_m
        geometry_quantitative_allowed = False
    else:
        channel_beams = load_registered_ffi_obstacle_geometry(
            args.obstacle_geometry_json, article_case="O1",
        )
        geometry_mode = "machine_registered_exact"
        geometry_assumption = (
            "machine-registered axis-aligned cuboid assembly loaded from "
            f"{Path(args.obstacle_geometry_json).resolve()}"
        )
        geometry_quantitative_allowed = True
    registered_projected_blockage_ratio = float(
        sum(
            (beam.x_max_m - beam.x_min_m) * (beam.z_max_m - beam.z_min_m)
            for beam in channel_beams
        ) / (0.9 * 0.8)
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
    all_sensors = tuple(
        FastObstacleReceptor3D(
            sensor.sensor_id,
            (sensor.position_m[1], sensor.position_m[0], sensor.position_m[2]),
        )
        for sensor in benchmark.sensor_receptors
    )
    sensors = tuple(
        sensor for sensor in all_sensors
        if all(0.0 <= value <= bound for value, bound in zip(sensor.position_m, domain))
        and not any(
            obstacle.x_min_m - size <= sensor.position_m[0] <= obstacle.x_max_m + size
            and obstacle.y_min_m - size <= sensor.position_m[1] <= obstacle.y_max_m + size
            and obstacle.z_min_m - size <= sensor.position_m[2] <= obstacle.z_max_m + size
            for obstacle in obstacles
        )
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=domain,
        cells=cells,
        wind_speed_m_s=float(args.wind_speed_m_s),
        kinematic_viscosity_m2_s=0.01,
        flow_steps=int(args.flow_steps),
        minimum_flow_steps=int(args.flow_steps),
        lattice_inlet_speed=0.03,
        duration_s=float(args.duration_s),
        scalar_time_step_s=physical_step_s,
        scalar_diffusivity_m2_s=float(args.diffusivity_m2_s),
        thermal_diffusivity_m2_s=float(args.diffusivity_m2_s),
        source_schedule=SourceRateSchedule(
            tuple(float(value) for value in source_t),
            tuple(float(value) for value in source_q),
            source_id="FFI-O1-measured-60s-resolved-prefix",
        ),
        source_position_m=(0.5, 0.45, 0.76),
        source_sigma_m=float(args.source_sigma_m),
        source_momentum_velocity_m_s=(0.0, 0.0, float(args.source_momentum_z_m_s)),
        source_momentum_velocity_history=source_momentum_history,
        obstacles=obstacles,
        source_temperature_k=293.0,
        thermal_buoyancy_coupling=True,
        subgrid_model=str(args.subgrid_model),
        smagorinsky_constant=float(args.smagorinsky_constant),
        wale_constant=float(args.wale_constant),
        effective_schmidt_number=(
            None if args.effective_schmidt_number is None
            else float(args.effective_schmidt_number)
        ),
        effective_prandtl_number=(
            None if args.effective_prandtl_number is None
            else float(args.effective_prandtl_number)
        ),
        scalar_turbulent_diffusivity_model=(
            str(args.scalar_turbulent_diffusivity_model)
        ),
        thermal_turbulent_diffusivity_model=(
            str(args.thermal_turbulent_diffusivity_model)
        ),
        flow_boundary_scheme=str(args.flow_boundary_scheme),
        nonstreamwise_wall_boundary=str(args.nonstreamwise_wall_boundary),
        active_momentum_solver=(
            "finite_volume_projection"
            if abs(float(args.source_momentum_z_m_s)) > 0.0
            else "d3q19_bgk"
        ),
        carrier_flow_initialization=(
            "projected_inlet"
            if abs(float(args.source_momentum_z_m_s)) > 0.0
            else "d3q19_spinup"
        ),
        maximum_lattice_mach_number=float(args.maximum_lattice_mach_number),
        maximum_active_scalar_density_departure=0.20,
        store_scalar_fields=False,
    )
    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(config, receptors=sensors)
    runtime = time.perf_counter() - started
    molar_density = result.ambient_pressure_pa / (R_GAS * result.ambient_temperature_k)
    predicted = {
        label: [
            float(value / H2_MOLAR_MASS / molar_density * 100.0)
            for value in trace
        ]
        for label, trace in result.receptor_traces_kg_m3.items()
    }

    def window_mean(sensor_id: str) -> float:
        sensor_time = np.asarray(run.sensor_time_s, dtype=float)
        sensor_values = np.asarray(run.sensor_percent[sensor_id], dtype=float)
        keep = (sensor_time >= t0) & (sensor_time <= t1)
        local_t = np.concatenate(([t0], sensor_time[keep], [t1]))
        local_y = np.interp(local_t, sensor_time, sensor_values)
        return float(np.trapezoid(local_y, local_t) / args.duration_s)

    measured = {sensor.label: window_mean(sensor.label) for sensor in sensors}
    predicted_mean = {
        label: float(np.mean(values[1:])) for label, values in predicted.items()
    }
    paired = [
        (predicted_mean[label], measured[label])
        for label in predicted_mean if label in measured
    ]
    errors = np.asarray([predicted - observed for predicted, observed in paired])
    observed = np.asarray([observed for _, observed in paired])
    diagnostics = {
        "sensor_count": len(paired),
        "mae_vol_percent": float(np.mean(np.abs(errors))),
        "rmse_vol_percent": float(np.sqrt(np.mean(errors ** 2))),
        "bias_vol_percent": float(np.mean(errors)),
        "mape_percent_floor_0p1_vol_percent": float(
            np.mean(np.abs(errors) / np.maximum(np.abs(observed), 0.1)) * 100.0
        ),
    }
    finite_max = float(np.nanmax(result.maximum_concentration_vol_percent))
    field_output = args.field_output or args.output.with_name(
        args.output.stem + "-fields.npz"
    )
    maximum_field = np.asarray(result.maximum_concentration_vol_percent, dtype=float)
    final_field = np.asarray(result.concentration_vol_percent, dtype=float)
    np.savez_compressed(
        field_output,
        x_m=np.asarray(result.x_m),
        y_m=np.asarray(result.y_m),
        z_m=np.asarray(result.z_m),
        maximum_concentration_vol_percent=maximum_field,
        final_concentration_vol_percent=final_field,
        solid_mask=np.asarray(result.solid_mask, dtype=bool),
    )
    four_vol_volume_m3 = float(
        np.count_nonzero(np.isfinite(maximum_field) & (maximum_field >= 4.0))
        * size**3
    )
    report = {
        "schema": "degali.ffi-o1-resolved-buoyant-prefix.v1",
        "status": "conditional_resolved_physics_witness_not_validation_score",
        "case": "O1",
        "window_s": [t0, t1],
        "grid": {"domain_m": list(domain), "cells": list(cells), "cell_size_m": size,
                 "physical_step_s": physical_step_s},
        "geometry_mode": geometry_mode,
        "geometry_assumption": geometry_assumption,
        "reported_blockage_ratio": 0.368,
        "registered_geometry_projected_blockage_ratio": (
            registered_projected_blockage_ratio
        ),
        "quantitative_geometry_allowed": geometry_quantitative_allowed,
        "physics": {
            "resolved_flow": True,
            "thermal_buoyancy_coupling": True,
            "composition_density_buoyancy": True,
            "scalar_diffusivity_m2_s": float(args.diffusivity_m2_s),
            "source_temperature_k": 293.0,
            "source_sigma_m": float(args.source_sigma_m),
            "source_momentum_velocity_m_s": [
                0.0, 0.0, float(args.source_momentum_z_m_s)
            ],
            "source_momentum_velocity_history_id": (
                None if source_momentum_history is None
                else source_momentum_history.history_id
            ),
            "maximum_lattice_mach_number": float(args.maximum_lattice_mach_number),
            "subgrid_model": str(args.subgrid_model),
            "smagorinsky_constant": float(args.smagorinsky_constant),
            "wale_constant": float(args.wale_constant),
            "effective_schmidt_number": (
                None if args.effective_schmidt_number is None
                else float(args.effective_schmidt_number)
            ),
            "effective_prandtl_number": (
                None if args.effective_prandtl_number is None
                else float(args.effective_prandtl_number)
            ),
            "scalar_turbulent_diffusivity_model": (
                str(args.scalar_turbulent_diffusivity_model)
            ),
            "thermal_turbulent_diffusivity_model": (
                str(args.thermal_turbulent_diffusivity_model)
            ),
            "flow_boundary_scheme": str(args.flow_boundary_scheme),
            "nonstreamwise_wall_boundary": str(args.nonstreamwise_wall_boundary),
            "source_orifice_diagnostic": {
                "diameter_m": FFI_INLET_DIAMETER_M,
                "positive_mass_flow_range_g_s": [
                    float(np.min(measured_source_q) * 1.0e3),
                    float(np.max(measured_source_q) * 1.0e3),
                ] if measured_source_q.size else [0.0, 0.0],
                "ideal_gas_h2_density_kg_m3": hydrogen_density_kg_m3,
                "inferred_jet_velocity_range_m_s": [
                    float(np.min(inferred_orifice_velocity_m_s)),
                    float(np.max(inferred_orifice_velocity_m_s)),
                ] if inferred_orifice_velocity_m_s.size else [0.0, 0.0],
                "interpretation": (
                    "diagnostic only: the resolved Gaussian source is not a "
                    "4.6 mm nozzle mesh"
                ),
            },
        },
        "runtime_s": runtime,
        "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        "maximum_energy_residual_j": result.maximum_energy_residual_j,
        "maximum_predicted_vol_percent": finite_max,
        "four_vol_percent_max_volume_m3": four_vol_volume_m3,
        "field_artifact": str(field_output),
        "sensor_count": len(sensors),
        "sensor_ids_omitted_at_boundary_or_solid_stencil": sorted(
            set(sensor.label for sensor in all_sensors) - set(predicted)
        ),
        "predicted_sensor_mean_vol_percent": predicted_mean,
        "predicted_sensor_time_s": [float(value) for value in result.time_s],
        "predicted_sensor_trace_vol_percent": predicted,
        "measured_sensor_mean_vol_percent_reference": measured,
        "diagnostic_sensor_error_not_validation_score": diagnostics,
        "promotion_allowed": False,
        "field_evidence_manifest": field_manifest,
        "warnings": [
            (
                "machine-registered geometry was loaded, but source/sensor clock "
                "and calibration gates remain unresolved"
                if geometry_quantitative_allowed else
                "conditional beam placement is a resolved sensitivity geometry"
            ),
            "FFI source/sensor clock and calibration gates remain unresolved",
            "positive carrier speed is a numerical transport condition; the facility is not claimed to have a measured wind",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
