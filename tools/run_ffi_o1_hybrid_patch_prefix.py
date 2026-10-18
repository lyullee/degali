"""Run a three-block FFI O1 fast/resolved/fast obstacle-patch prefix.

The outer blocks use the fast conservative transport model.  Only the middle
block contains the conditional O1 beam surrogate and uses the resolved
obstacle-flow solver.  Two shared-clock bridges exchange species mass,
enthalpy (when present), and momentum packets at the patch faces.  This is a
hybrid coupling witness: the beam geometry and FFI clocks remain conditional,
but unlike the older overlap smoke it uses real O1 sensors and non-overlapping
streamwise blocks.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleReceptor3D,
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.hybrid_flux_packets import CartesianInterfacePlane, boundary_cell_mask
from degali.addons.hybrid_solver_bridge import SynchronizedSolverBridge
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
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
AMBIENT_TEMPERATURE_K = 293.15
AMBIENT_PRESSURE_PA = 101_325.0
PATCH_X_MIN_M = 1.8
PATCH_X_MAX_M = 3.4
CELL_SIZE_M = 0.1


def _chain_one(*callbacks):
    def chained(*args):
        value = args[0]
        rest = args[1:]
        for callback in callbacks:
            value = callback(value, *rest)
        return value
    return chained


def _chain_conserved_state(*callbacks):
    """Compose resolved callbacks carrying concentration and enthalpy."""
    def chained(concentration, sensible_enthalpy, start, dt):
        values = (concentration, sensible_enthalpy)
        for callback in callbacks:
            values = callback(*values, start, dt)
        return values
    return chained


def _chain_six(*callbacks):
    def chained(u, v, w, u_face, v_face, w_face, start, dt):
        values = (u, v, w, u_face, v_face, w_face)
        for callback in callbacks:
            values = callback(*values, start, dt)
        return values
    return chained


def _chain_flux(*callbacks):
    def chained(flux_x, flux_y, flux_z):
        for callback in callbacks:
            callback(flux_x, flux_y, flux_z)
    return chained


def _window_source(run, duration_s: float) -> SourceRateSchedule:
    t0 = 60.0
    t1 = t0 + duration_s
    keep = (run.flow_time_s >= t0) & (run.flow_time_s <= t1)
    time_s = np.asarray(run.flow_time_s[keep] - t0, dtype=float)
    flow_kg_s = np.asarray(run.mass_flow_g_s[keep] * 1.0e-3, dtype=float)
    if time_s.size < 2:
        raise ValueError("O1 source history does not cover requested prefix")
    time_s[0] = 0.0
    time_s[-1] = duration_s
    if flow_kg_s[-1] != 0.0:
        time_s[-1] = max(0.0, duration_s - 1.0e-9)
        time_s = np.append(time_s, duration_s)
        flow_kg_s = np.append(flow_kg_s, 0.0)
    return SourceRateSchedule(
        tuple(float(value) for value in time_s),
        tuple(float(value) for value in flow_kg_s),
        source_id="FFI-O1-measured-60s-hybrid-prefix",
    )


def _receptors(
    benchmark, *, block: str, patch_x_min_m: float, patch_x_max_m: float,
    obstacles: tuple[AxisAlignedCuboid, ...] = (),
) -> tuple[FastObstacleReceptor3D, ...]:
    selected = []
    for sensor in benchmark.sensor_receptors:
        streamwise = float(sensor.position_m[1])
        if block == "upstream" and not streamwise < patch_x_min_m:
            continue
        if block == "patch" and not patch_x_min_m <= streamwise <= patch_x_max_m:
            continue
        if block == "downstream" and not streamwise > patch_x_max_m:
            continue
        local_x = streamwise if block == "upstream" else (
            streamwise - patch_x_min_m if block == "patch" else streamwise - patch_x_max_m
        )
        if any(
            obstacle.x_min_m - CELL_SIZE_M <= local_x <= obstacle.x_max_m + CELL_SIZE_M
            and obstacle.y_min_m - CELL_SIZE_M <= float(sensor.position_m[0]) <= obstacle.y_max_m + CELL_SIZE_M
            and obstacle.z_min_m - CELL_SIZE_M <= float(sensor.position_m[2]) <= obstacle.z_max_m + CELL_SIZE_M
            for obstacle in obstacles
        ):
            continue
        selected.append(FastObstacleReceptor3D(
            sensor.sensor_id, (local_x, float(sensor.position_m[0]), float(sensor.position_m[2]))
        ))
    return tuple(selected)


def _bridge_pair(
    *, duration_s: float,
    time_step_s: float,
    expects_enthalpy: bool,
    cells_upstream: tuple[int, int, int],
    cells_patch: tuple[int, int, int],
    cells_downstream: tuple[int, int, int],
):
    area_a = CELL_SIZE_M * CELL_SIZE_M
    patch_dy = 0.9 / cells_patch[1]
    patch_dz = 0.8 / cells_patch[2]
    area_b = patch_dy * patch_dz
    volume_a = CELL_SIZE_M ** 3
    volume_b = CELL_SIZE_M * patch_dy * patch_dz
    def array_shape(cells: tuple[int, int, int]) -> tuple[int, int, int]:
        # Solver fields are stored as (z, y, x); configs are (x, y, z).
        return (cells[2], cells[1], cells[0])
    incoming = SynchronizedSolverBridge(time_step_s=time_step_s, timeout_s=120.0)
    incoming.register_side(
        name="upstream", plane=CartesianInterfacePlane(0, cells_upstream[0], 1.0, area_a),
        target_mask=boundary_cell_mask(array_shape(cells_upstream), axis=2, side="max"),
        cell_volume_m3=volume_a,
    )
    incoming.register_side(
        name="patch_in", plane=CartesianInterfacePlane(0, 0, -1.0, area_b),
        target_mask=boundary_cell_mask(array_shape(cells_patch), axis=2, side="min"),
        cell_volume_m3=volume_b,
        expects_enthalpy=expects_enthalpy,
        momentum_application="characteristic",
    )
    outgoing = SynchronizedSolverBridge(time_step_s=time_step_s, timeout_s=120.0)
    outgoing.register_side(
        name="patch_out", plane=CartesianInterfacePlane(0, cells_patch[0], 1.0, area_b),
        target_mask=boundary_cell_mask(array_shape(cells_patch), axis=2, side="max"),
        cell_volume_m3=volume_b,
        expects_enthalpy=expects_enthalpy,
        momentum_application="characteristic",
    )
    outgoing.register_side(
        name="downstream", plane=CartesianInterfacePlane(0, 0, -1.0, area_a),
        target_mask=boundary_cell_mask(array_shape(cells_downstream), axis=2, side="min"),
        cell_volume_m3=volume_a,
    )
    return incoming, outgoing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration-s", type=float, default=3.0)
    parser.add_argument("--time-step-s", type=float, default=0.02)
    parser.add_argument(
        "--flow-steps", type=int, default=40,
        help="resolved patch D3Q19 flow iterations per scalar step",
    )
    parser.add_argument(
        "--isothermal", action="store_true",
        help="disable active thermal buoyancy for numerical time-convergence audits",
    )
    parser.add_argument(
        "--patch-x-min-m", type=float, default=PATCH_X_MIN_M,
        help="global streamwise start of the resolved obstacle patch",
    )
    parser.add_argument(
        "--patch-x-max-m", type=float, default=PATCH_X_MAX_M,
        help="global streamwise end of the resolved obstacle patch",
    )
    parser.add_argument(
        "--subgrid-model", choices=("none", "smagorinsky", "wale"),
        default="smagorinsky",
    )
    parser.add_argument("--smagorinsky-constant", type=float, default=0.0)
    parser.add_argument("--wale-constant", type=float, default=0.5)
    parser.add_argument("--effective-schmidt-number", type=float, default=None)
    parser.add_argument("--effective-prandtl-number", type=float, default=None)
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
    parser.add_argument(
        "--obstacle-geometry-json", type=Path, default=None,
        help=(
            "machine-registered exact O1 cuboid assembly; omitted means the "
            "conditional two-beam surrogate"
        ),
    )
    parser.add_argument("--field-evidence-manifest", type=Path, default=None)
    parser.add_argument(
        "--require-quantitative-manifest", action="store_true",
        help="fail before solving unless the attached manifest clock/calibration gate is accepted",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "outputs/ffi-o1-hybrid-patch-prefix-2026-10-14.json",
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
    if (
        args.duration_s <= 0.0 or args.time_step_s <= 0.0
        or abs(args.duration_s / args.time_step_s - round(args.duration_s / args.time_step_s)) > 1e-10
    ):
        raise ValueError("duration-s must be a positive multiple of time-step-s")
    if args.flow_steps < 1:
        raise ValueError("flow-steps must be positive")
    patch_x_min_m = float(args.patch_x_min_m)
    patch_x_max_m = float(args.patch_x_max_m)
    if not np.isfinite(patch_x_min_m) or not np.isfinite(patch_x_max_m):
        raise ValueError("patch bounds must be finite")
    if not 0.0 < patch_x_min_m < patch_x_max_m < 5.8:
        raise ValueError("patch bounds must lie strictly inside the 5.8 m channel")

    benchmark = load_open_channel_h2_obstacle_benchmark(REGISTRY, DATA_ROOT)
    thermal_enabled = not args.isothermal
    case = next(item for item in benchmark.cases if item.article_case == "O1")
    run = case.runs[0]
    source = _window_source(run, args.duration_s)
    upstream_length = patch_x_min_m
    patch_length = patch_x_max_m - patch_x_min_m
    downstream_length = 5.8 - patch_x_max_m
    upstream_cells = (int(round(upstream_length / CELL_SIZE_M)), 9, 8)
    patch_cells = (int(round(patch_length / CELL_SIZE_M)), 9, 8)
    downstream_cells = (int(round(downstream_length / CELL_SIZE_M)), 9, 8)
    # Keep the same conservative edge margin used by the standalone resolved
    # witness so the patch leaves at least one fetch cell above and beside the
    # surrogate beams.  The nominal zero-margin geometry touches the top wall
    # at this coarse 0.1 m patch resolution and is rejected by the resolved
    # solver's boundary-fetch guard.
    if args.obstacle_geometry_json is None:
        beams = conditional_ffi_o1_crossbeam_assembly(edge_margin_m=0.21)
        geometry_mode = "conditional_two_crossbeam_surrogate"
        geometry_quantitative_allowed = False
    else:
        beams = load_registered_ffi_obstacle_geometry(
            args.obstacle_geometry_json, article_case="O1",
        )
        geometry_mode = "machine_registered_exact"
        geometry_quantitative_allowed = True
    patch_obstacles_list = []
    for beam in beams:
        if geometry_quantitative_allowed:
            fetch_tolerance = 1.0e-10
            if not (
                beam.y_min_m >= patch_x_min_m + 3.0 * CELL_SIZE_M - fetch_tolerance
                and beam.y_max_m <= patch_x_max_m - 3.0 * CELL_SIZE_M + fetch_tolerance
            ):
                raise ValueError(
                    "machine-registered obstacle must lie wholly inside the "
                    "resolved patch with a three-cell streamwise fetch; "
                    "refusing geometry clipping"
                )
            local_y_min = beam.y_min_m - patch_x_min_m
            local_y_max = beam.y_max_m - patch_x_min_m
        else:
            local_y_min = beam.y_min_m - patch_x_min_m
            local_y_max = min(
                beam.y_max_m - patch_x_min_m,
                patch_length - 3.0 * CELL_SIZE_M - 1.0e-6,
            )
        patch_obstacles_list.append(AxisAlignedCuboid(
            local_y_min,
            local_y_max,
            beam.x_min_m,
            beam.x_max_m,
            beam.z_min_m,
            beam.z_max_m,
            beam.label,
        ))
    patch_obstacles = tuple(patch_obstacles_list)
    zero_source = SourceRateSchedule((0.0, args.duration_s), (0.0, 0.0), source_id="hybrid-patch-zero")
    upstream_config = FastObstacleTransport3DConfig(
        domain_m=(upstream_length, 0.9, 0.8), cells=upstream_cells,
        duration_s=args.duration_s, time_step_s=args.time_step_s, wind_speed_m_s=0.15,
        diffusivity_m2_s=0.01, source_schedule=source,
        source_position_m=(0.5, 0.45, 0.76), source_sigma_m=0.08,
        obstacles=(), open_boundary_faces=("x_max",), store_fields=True,
    )
    patch_config = ResolvedObstacleFlow3DConfig(
        domain_m=(patch_length, 0.9, 0.8), cells=patch_cells,
        wind_speed_m_s=0.15, kinematic_viscosity_m2_s=0.01,
        flow_steps=args.flow_steps, minimum_flow_steps=args.flow_steps,
        lattice_inlet_speed=0.03,
        duration_s=args.duration_s, scalar_time_step_s=args.time_step_s,
        scalar_diffusivity_m2_s=0.01, source_schedule=zero_source,
        source_position_m=(0.1, 0.45, 0.2), source_sigma_m=0.08,
        obstacles=patch_obstacles, store_scalar_fields=True,
        thermal_buoyancy_coupling=thermal_enabled,
        source_temperature_k=(293.15 if thermal_enabled else None),
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
        # The conservative interface momentum packet adds a small numerical
        # overshoot above the standalone 0.15 Mach witness; retain the gate
        # explicitly at 0.30 and report it as a conditional coupling limit.
        maximum_lattice_mach_number=0.30,
    )
    downstream_config = FastObstacleTransport3DConfig(
        domain_m=(downstream_length, 0.9, 0.8), cells=downstream_cells,
        duration_s=args.duration_s, time_step_s=args.time_step_s, wind_speed_m_s=0.15,
        diffusivity_m2_s=0.01, source_schedule=zero_source,
        source_position_m=(0.1, 0.45, 0.2), source_sigma_m=0.08,
        obstacles=(), open_boundary_faces=("x_min", "x_max"), store_fields=True,
    )
    incoming, outgoing = _bridge_pair(
        duration_s=args.duration_s,
        time_step_s=args.time_step_s,
        expects_enthalpy=thermal_enabled,
        cells_upstream=upstream_cells,
        cells_patch=patch_cells,
        cells_downstream=downstream_cells,
    )
    up_cb = incoming.callbacks_for("upstream")
    patch_in_cb = incoming.callbacks_for("patch_in")
    patch_out_cb = outgoing.callbacks_for("patch_out")
    down_cb = outgoing.callbacks_for("downstream")
    patch_state = _chain_conserved_state(
        patch_in_cb["interface_conserved_state_callback"],
        patch_out_cb["interface_conserved_state_callback"],
    )
    def patch_mass(start, dt):
        incoming_mass, incoming_enthalpy = patch_in_cb[
            "interface_conserved_mass_callback"
        ](start, dt)
        outgoing_mass, outgoing_enthalpy = patch_out_cb[
            "interface_conserved_mass_callback"
        ](start, dt)
        return (
            incoming_mass + outgoing_mass,
            incoming_enthalpy + outgoing_enthalpy,
        )
    patch_momentum = _chain_six(
        patch_in_cb["interface_momentum_callback"],
        patch_out_cb["interface_momentum_callback"],
    )
    patch_enthalpy_flux = _chain_flux(
        patch_in_cb["enthalpy_flux_callback"],
        patch_out_cb["enthalpy_flux_callback"],
    )
    patch_flux = _chain_flux(
        patch_in_cb["face_flux_callback"], patch_out_cb["face_flux_callback"],
    )
    started = time.perf_counter()

    def run_upstream():
        return solve_fast_obstacle_transport_3d(
            upstream_config, receptors=_receptors(
                benchmark, block="upstream",
                patch_x_min_m=patch_x_min_m, patch_x_max_m=patch_x_max_m,
            ),
            face_flux_callback=up_cb["face_flux_callback"],
            interface_state_callback=up_cb["interface_state_callback"],
            interface_state_mass_callback=up_cb["interface_state_mass_callback"],
            interface_momentum_callback=up_cb["interface_momentum_callback"],
        )

    def run_patch():
        return solve_resolved_obstacle_flow_3d(
            patch_config, receptors=_receptors(
                benchmark, block="patch", patch_x_min_m=patch_x_min_m,
                patch_x_max_m=patch_x_max_m, obstacles=patch_obstacles,
            ),
            face_flux_callback=patch_flux,
            enthalpy_flux_callback=(patch_enthalpy_flux if thermal_enabled else None),
            interface_conserved_state_callback=patch_state,
            interface_conserved_mass_callback=patch_mass,
            interface_momentum_callback=patch_momentum,
        )

    def run_downstream():
        return solve_fast_obstacle_transport_3d(
            downstream_config, receptors=_receptors(
                benchmark, block="downstream",
                patch_x_min_m=patch_x_min_m, patch_x_max_m=patch_x_max_m,
            ),
            face_flux_callback=down_cb["face_flux_callback"],
            interface_state_callback=down_cb["interface_state_callback"],
            interface_state_mass_callback=down_cb["interface_state_mass_callback"],
            interface_momentum_callback=down_cb["interface_momentum_callback"],
        )

    with ThreadPoolExecutor(max_workers=3) as executor:
        # Collect the resolved patch first so a patch-side interface failure is
        # reported directly instead of being masked by an outer block waiting
        # on the shared-clock exchange.
        futures = [executor.submit(run_patch), executor.submit(run_upstream), executor.submit(run_downstream)]
        results = [future.result() for future in futures]
    patch, upstream, downstream = results
    final_step = round(args.duration_s / args.time_step_s) - 1
    incoming.finalize_step(final_step)
    outgoing.finalize_step(final_step)
    runtime = time.perf_counter() - started

    def total_history(field: str) -> np.ndarray:
        return sum(np.asarray(getattr(result, field), dtype=float) for result in results)

    total_residual = total_history("domain_mass_kg") + total_history("outflow_mass_kg") - total_history("injected_mass_kg")
    fields = [
        np.nan_to_num(upstream.maximum_concentration_kg_m3, nan=0.0),
        np.nan_to_num(patch.maximum_concentration_kg_m3, nan=0.0),
        np.nan_to_num(downstream.maximum_concentration_kg_m3, nan=0.0),
    ]
    assembled_max = np.concatenate(fields, axis=2)
    molar_density = AMBIENT_PRESSURE_PA / (R_GAS * AMBIENT_TEMPERATURE_K)
    assembled_vol = assembled_max / H2_MOLAR_MASS / molar_density * 100.0
    sensor_means = {}
    sensor_means.update({key: float(np.mean(value[1:]) / H2_MOLAR_MASS / molar_density * 100.0) for key, value in upstream.receptor_traces_kg_m3.items()})
    # Resolved results already expose this trace in vol%; do not convert it a
    # second time from kg/m3 (that would inflate the patch sensors by ~400x).
    sensor_means.update({key: float(np.mean(value[1:])) for key, value in patch.receptor_traces_vol_percent.items()})
    sensor_means.update({key: float(np.mean(value[1:]) / H2_MOLAR_MASS / molar_density * 100.0) for key, value in downstream.receptor_traces_kg_m3.items()})
    measured = {}
    sensor_time = np.asarray(run.sensor_time_s, dtype=float)
    for sensor_id, values in run.sensor_percent.items():
        keep = (sensor_time >= 60.0) & (sensor_time <= 60.0 + args.duration_s)
        local_t = np.concatenate(([60.0], sensor_time[keep], [60.0 + args.duration_s]))
        local_y = np.interp(local_t, sensor_time, np.asarray(values, dtype=float))
        measured[sensor_id] = float(np.trapezoid(local_y, local_t) / args.duration_s)
    paired = [(sensor_means[key], measured[key]) for key in sensor_means if key in measured]
    errors = np.asarray([a - b for a, b in paired], dtype=float)
    report = {
        "schema": "degali.ffi-o1-hybrid-patch-prefix.v1",
        "status": "conditional_nonoverlap_fast_resolved_fast_coupling_witness",
        "quantitative_ready": False,
        "field_evidence_manifest": field_manifest,
        "geometry_mode": geometry_mode,
        "geometry_quantitative_allowed": geometry_quantitative_allowed,
        "duration_s": args.duration_s,
        "time_step_s": args.time_step_s,
        "resolved_flow_steps_per_scalar_step": args.flow_steps,
        "blocks": {
            "upstream": {"domain_m": [upstream_length, 0.9, 0.8], "cells": list(upstream_cells), "solver": "fast"},
            "patch": {
                "domain_m": [patch_length, 0.9, 0.8],
                "cells": list(patch_cells),
                "solver": "resolved",
                "global_x_m": [patch_x_min_m, patch_x_max_m],
                "thermal_buoyancy_coupling": thermal_enabled,
                "subgrid_model": str(args.subgrid_model),
                "smagorinsky_constant": float(args.smagorinsky_constant),
                "wale_constant": float(args.wale_constant),
                "effective_schmidt_number": args.effective_schmidt_number,
                "effective_prandtl_number": args.effective_prandtl_number,
                "scalar_turbulent_diffusivity_model": str(
                    args.scalar_turbulent_diffusivity_model
                ),
                "thermal_turbulent_diffusivity_model": str(
                    args.thermal_turbulent_diffusivity_model
                ),
                "flow_boundary_scheme": str(args.flow_boundary_scheme),
                "nonstreamwise_wall_boundary": str(
                    args.nonstreamwise_wall_boundary
                ),
            },
            "downstream": {"domain_m": [downstream_length, 0.9, 0.8], "cells": list(downstream_cells), "solver": "fast"},
        },
        "bridge": {"incoming_exchanged_steps": incoming.exchanged_steps, "outgoing_exchanged_steps": outgoing.exchanged_steps},
        "interface_relaxation": {
            "patch_momentum_application": "characteristic_packet_velocity",
            "species_mass_exchange": "unrelaxed_conservative_packet",
            "thermal_enthalpy_exchange": thermal_enabled,
        },
        "runtime_s": runtime,
        "sensor_count": len(paired),
        "predicted_sensor_mean_vol_percent": sensor_means,
        "measured_sensor_mean_vol_percent_reference": measured,
        "diagnostic_sensor_error_not_validation_score": {
            "mae_vol_percent": float(np.mean(np.abs(errors))) if errors.size else None,
            "rmse_vol_percent": float(np.sqrt(np.mean(errors ** 2))) if errors.size else None,
            "bias_vol_percent": float(np.mean(errors)) if errors.size else None,
        },
        "maximum_predicted_vol_percent": float(np.nanmax(assembled_vol)),
        "four_vol_percent_max_volume_m3": float(np.count_nonzero(assembled_vol >= 4.0) * CELL_SIZE_M ** 3),
        "maximum_partitioned_mass_balance_residual_kg": float(np.max(np.abs(total_residual))),
        "maximum_block_mass_residual_kg": float(max(result.maximum_mass_residual_kg for result in results)),
        "resolved_patch_energy_residual_j": float(
            patch.maximum_energy_residual_j
            if patch.maximum_energy_residual_j is not None else 0.0
        ),
        "resolved_patch_continuity_residual_kg_m3_s": float(
            patch.diagnostics.get(
                "maximum_variable_density_continuity_residual_kg_m3_s", 0.0,
            )
        ),
        "resolved_patch_max_abs_pressure_correction_pa": float(
            np.nanmax(np.abs(patch.pressure_correction_pa))
            if patch.pressure_correction_pa is not None else 0.0
        ),
        "resolved_patch_diagnostics": {
            "maximum_lattice_mach_number": float(
                patch.diagnostics.get("maximum_lattice_mach_number", float("nan"))
            ),
            "maximum_allowed_lattice_mach_number": float(
                patch.diagnostics.get("maximum_allowed_lattice_mach_number", float("nan"))
            ),
            "maximum_relative_lattice_density_variation": float(
                patch.diagnostics.get("maximum_relative_lattice_density_variation", float("nan"))
            ),
        },
        "field_shape_assembled": list(assembled_vol.shape),
        "warnings": [
            (
                "machine-registered geometry was loaded; source/sensor clock and "
                "calibration gates remain unresolved"
                if geometry_quantitative_allowed else
                "beam geometry is the conditional O1 surrogate, not registered facility CAD"
            ),
            "patch receives conservative packets through boundary strips; full six-face AMR remains future work",
            (
                "resolved patch has active thermal buoyancy but zero local source; upstream enthalpy remains ambient in this prefix"
                if thermal_enabled else
                "isothermal mode is used only for the declared time-step/interface convergence audit"
            ),
            "resolved patch interface momentum uses characteristic packet velocity with a conditional Mach gate of 0.30; this is not a qualified CFD closure",
            "sensor errors are diagnostic and not a validation score",
        ],
    }
    field_path = args.output.with_name(args.output.stem + "-fields.npz")
    x_m = np.concatenate((
        np.asarray(upstream.x_m, dtype=float),
        np.asarray(patch.x_m, dtype=float) + patch_x_min_m,
        np.asarray(downstream.x_m, dtype=float) + patch_x_max_m,
    ))
    y_m = np.asarray(upstream.y_m, dtype=float)
    z_m = np.asarray(upstream.z_m, dtype=float)
    solid_mask = np.concatenate((upstream.solid_mask, patch.solid_mask, downstream.solid_mask), axis=2)
    np.savez_compressed(
        field_path,
        maximum_h2_vol_percent=assembled_vol,
        maximum_concentration_vol_percent=assembled_vol,
        x_m=x_m, y_m=y_m, z_m=z_m, solid_mask=solid_mask,
    )
    report["field_artifact"] = str(field_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
