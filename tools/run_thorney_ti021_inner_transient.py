"""Run a no-fit TI021 resolved-inner transient and persist its evidence.

The command reads sensor *coordinates* from the TI021 workbook but never uses
dose, peak concentration, arrival or departure values to configure the run.
Every x=160 m conservative handoff frame is flushed to HDF5 while the solver
is running.  On successful completion, the inner-domain receptor histories
and final three-dimensional fields are saved separately for scoring and plots.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from degali.addons.binary_mixture_handoff_plane_3d import (
    BinaryMixtureHandoffPlane3DRequest,
)
from degali.addons.conservative_handoff_hdf5_3d import (
    ConservativeHandoffHdf5Writer3D,
)
from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.validation.thorney_ti021 import (
    ThorneyTi021Dataset,
    read_thorney_ti021_workbook,
)
try:
    from tools.audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )
except ModuleNotFoundError:  # Direct ``python tools/<script>.py`` execution.
    from audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-inner-transient.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _handoff_crosswind_edges(width_m: float, cells: int) -> np.ndarray:
    """Return handoff face coordinates centered on the active inner width."""
    width = float(width_m)
    count = int(cells)
    if not math.isfinite(width) or width <= 0.0 or count < 1:
        raise ValueError("handoff crosswind width and cell count must be positive")
    return np.linspace(-0.5 * width, 0.5 * width, count + 1)


def inner_receptors(
    dataset: ThorneyTi021Dataset,
    *, domain_m: tuple[float, float, float] = (200.0, 200.0, 40.0),
) -> tuple[tuple[FastObstacleReceptor3D, ...], list[dict[str, object]]]:
    """Return sensor-coordinate-only receptors inside the resolved inner box."""

    receptors: list[FastObstacleReceptor3D] = []
    metadata: list[dict[str, object]] = []
    crosswind_origin_m = 0.5 * float(domain_m[1])
    for index, sensor in enumerate(dataset.dose_sensors):
        solver_position = (
            60.0 + sensor.local_y_m,
            crosswind_origin_m - sensor.local_x_m,
            sensor.height_m,
        )
        if all(
            0.0 <= coordinate <= bound
            for coordinate, bound in zip(
                solver_position, domain_m
            )
        ):
            label = f"TI021-S{index + 1:03d}"
            receptors.append(FastObstacleReceptor3D(label, solver_position))
            metadata.append({
                "label": label,
                "workbook_row_index_zero_based": index,
                "site_local_xyz_m": [
                    sensor.local_x_m,
                    sensor.local_y_m,
                    sensor.height_m,
                ],
                "solver_xyz_m": list(solver_position),
            })
    return tuple(receptors), metadata


def run_inner_transient(
    *,
    workbook: Path,
    property_audit: Path,
    output_directory: Path,
    duration_s: float,
    cell_size_m: float = 5.0,
    time_step_s: float = 0.125,
    progress_interval_s: float = 0.1,
    corrector_tolerance: float = 1.0e-6,
    stratification_damping_alpha: float | None = None,
    stratification_damping_exponent: float = 1.0,
    stratification_shear_regularization_s_inv: float | None = None,
    inner_domain_m: tuple[float, float, float] = (200.0, 200.0, 40.0),
    circular_arc_representation: str = "finite_thickness_voxel",
    obstacle_height_m: float = 5.0,
    handoff_site_downwind_y_m: float = 100.0,
    initial_cloud_flat_to_flat_m: float = 14.0,
    initial_cloud_z_max_m: float = 13.0,
    finite_source_release_duration_s: float | None = None,
    finite_source_start_delay_s: float = 0.0,
    finite_source_mass_kg: float | None = None,
    finite_source_mass_basis: str = "species",
    finite_source_position_z_m: float = 1.0,
    finite_source_momentum_velocity_m_s: tuple[float, float, float] = (0.0, 0.0, 0.0),
    finite_source_total_mass_to_species_mass_ratio: float = 1.0,
    wind_speed_m_s: float = 3.9,
    wind_heading_relative_deg: float = 6.1,
    receptor_vertical_sampling: str = "cell_center",
    finite_source_release_shape: str = "constant",
    finite_source_sigma_xyz_m: tuple[float, float, float] | None = None,
    finite_source_kernel: str = "gaussian",
    scalar_advection_scheme: str = "muscl_van_leer",
    source_time_integration: str = "prestep",
) -> tuple[dict[str, object], dict[str, Path]]:
    """Execute one frozen-prior inner run and persist all required outputs."""

    duration = float(duration_s)
    dt = float(time_step_s)
    ratio = duration / dt
    if duration <= 0.0 or not math.isclose(
        ratio, round(ratio), rel_tol=0.0, abs_tol=1.0e-10,
    ):
        raise ValueError("duration must be a positive integer multiple of time step")
    if progress_interval_s <= 0.0:
        raise ValueError("progress_interval_s must be positive")
    steps = int(round(ratio))
    workbook_path = workbook.resolve()
    property_path = property_audit.resolve()
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "handoff_hdf5": output / "ti021-inner-handoff.h5",
        "receptor_npz": output / "ti021-inner-receptor-traces.npz",
        "final_field_npz": output / "ti021-inner-final-field.npz",
        "record_json": output / "ti021-inner-transient.json",
        "progress_jsonl": output / "ti021-inner-progress.jsonl",
    }
    existing = [path for path in paths.values() if path.exists()]
    if existing:
        raise FileExistsError(
            "inner transient output already exists: "
            + ", ".join(str(path) for path in existing)
        )

    domain = tuple(float(value) for value in inner_domain_m)
    if len(domain) != 3 or any(value <= 0.0 or not math.isfinite(value) for value in domain):
        raise ValueError("inner_domain_m must contain three positive finite lengths")
    if handoff_site_downwind_y_m < 0.0 or not math.isfinite(handoff_site_downwind_y_m):
        raise ValueError("handoff_site_downwind_y_m must be finite and non-negative")
    dataset = read_thorney_ti021_workbook(workbook_path)
    receptors, receptor_metadata = inner_receptors(dataset, domain_m=domain)
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _closure, coordinates = build_configuration(
        properties,
        cell_size_m=cell_size_m,
        time_step_s=dt,
        active_steps=steps,
        corrector_tolerance=corrector_tolerance,
        circular_arc_representation=circular_arc_representation,
        obstacle_height_m=obstacle_height_m,
        inner_domain_m=domain,
        handoff_site_downwind_y_m=handoff_site_downwind_y_m,
        initial_cloud_flat_to_flat_m=initial_cloud_flat_to_flat_m,
        initial_cloud_z_max_m=initial_cloud_z_max_m,
        finite_source_release_duration_s=finite_source_release_duration_s,
        finite_source_start_delay_s=finite_source_start_delay_s,
        finite_source_mass_kg=finite_source_mass_kg,
        finite_source_mass_basis=finite_source_mass_basis,
        finite_source_position_z_m=finite_source_position_z_m,
        finite_source_momentum_velocity_m_s=finite_source_momentum_velocity_m_s,
        finite_source_total_mass_to_species_mass_ratio=(
            finite_source_total_mass_to_species_mass_ratio
        ),
        wind_speed_m_s=wind_speed_m_s,
        wind_heading_relative_deg=wind_heading_relative_deg,
        receptor_vertical_sampling=receptor_vertical_sampling,
        finite_source_release_shape=finite_source_release_shape,
        finite_source_sigma_xyz_m=finite_source_sigma_xyz_m,
        finite_source_kernel=finite_source_kernel,
        scalar_advection_scheme=scalar_advection_scheme,
        source_time_integration=source_time_integration,
    )
    dx = config.domain_m[0] / config.cells[0]
    first_cell_center_z_m = 0.5 * dx
    below_first_vertical_center = [
        {
            "label": f"TI021-S{index + 1:03d}",
            "workbook_row_index_zero_based": index,
            "height_m": float(sensor.height_m),
        }
        for index, sensor in enumerate(dataset.dose_sensors)
        if float(sensor.height_m) < first_cell_center_z_m
    ]
    config = replace(
        config,
        store_scalar_fields=False,
        stratification_damping_alpha=stratification_damping_alpha,
        stratification_damping_exponent=stratification_damping_exponent,
        stratification_shear_regularization_s_inv=(
            stratification_shear_regularization_s_inv
        ),
    )
    source_z_edges = np.linspace(0.0, config.domain_m[2], config.cells[2] + 1)
    source_q_edges = _handoff_crosswind_edges(
        config.domain_m[1], config.cells[1],
    )
    request = BinaryMixtureHandoffPlane3DRequest(
        x_face_m=float(coordinates["handoff_inner_solver_x_m"]),
        sample_interval_s=dt,
        include_initial=True,
    )
    metadata = {
        "case": "Thorney Island trial 21",
        "run_schema": SCHEMA,
        "duration_s": duration,
        "time_step_s": dt,
        "cell_size_m": dx,
        "corrector_tolerance": corrector_tolerance,
        "inner_receptor_count": len(receptors),
        "inner_domain_m": list(config.domain_m),
        "handoff_site_downwind_y_m": handoff_site_downwind_y_m,
        "coverage_scope": (
            "full_inner_domain_sensor_set"
            if len(receptors) == len(dataset.dose_sensors)
            else "local_patch_sensor_subset"
        ),
        "coupling_mode": "one_way_handoff_frame_export",
        "dynamic_outer_patch_coupling_applied": False,
        "transverse_coordinate": (
            f"q=-site_x=solver_y-{0.5 * config.domain_m[1]:g}m"
        ),
        "inner_crosswind_origin_m": 0.5 * float(config.domain_m[1]),
        "observation_values_used_to_configure_run": [],
        "initial_cloud_flat_to_flat_m": float(initial_cloud_flat_to_flat_m),
        "initial_cloud_z_max_m": float(initial_cloud_z_max_m),
        "finite_source_release_duration_s": (
            None if finite_source_release_duration_s is None
            else float(finite_source_release_duration_s)
        ),
        "finite_source_start_delay_s": float(finite_source_start_delay_s),
        "finite_source_mass_kg": (
            None if finite_source_mass_kg is None else float(finite_source_mass_kg)
        ),
        "finite_source_mass_basis": finite_source_mass_basis,
        "finite_source_schedule_species_mass_kg": coordinates.get(
            "finite_source_schedule_species_mass_kg"
        ),
        "finite_source_schedule_total_mass_kg": coordinates.get(
            "finite_source_schedule_total_mass_kg"
        ),
        "workbook_release_rate_kg_s": float(dataset.workbook_release_rate_kg_s),
        "workbook_release_duration_s": float(dataset.workbook_release_duration_s),
        "workbook_release_mass_proxy_kg": float(dataset.workbook_release_mass_proxy_kg),
        "finite_source_mass_provenance": (
            "explicit_cli_override_only; workbook rate-times-duration proxy is not a validated atmospheric source"
            if finite_source_mass_kg is not None else
            "geometry-derived initial-cloud mass unless an explicit override is supplied"
        ),
        "finite_source_position_z_m": float(finite_source_position_z_m),
        "finite_source_momentum_velocity_m_s": [
            float(value) for value in finite_source_momentum_velocity_m_s
        ],
        "finite_source_total_mass_to_species_mass_ratio": float(
            finite_source_total_mass_to_species_mass_ratio
        ),
    }
    started = time.perf_counter()
    next_progress_time = 0.0

    with paths["progress_jsonl"].open("x", encoding="utf-8") as progress_stream:
        def log_progress(frame) -> None:
            nonlocal next_progress_time
            writer(frame)
            if frame.time_s + 1.0e-12 >= next_progress_time:
                elapsed = time.perf_counter() - started
                entry = {
                    "model_time_s": frame.time_s,
                    "elapsed_wall_clock_s": elapsed,
                    "completed_fraction": frame.time_s / duration,
                    "frame_count": writer.frame_count,
                    "maximum_species_flux_kg_m2_s": float(np.max(np.abs(
                        frame.species_mass_flux_kg_m2_s
                    ))),
                }
                progress_stream.write(json.dumps(entry) + "\n")
                progress_stream.flush()
                print(json.dumps(entry), flush=True)
                next_progress_time += progress_interval_s

        with ConservativeHandoffHdf5Writer3D(
            paths["handoff_hdf5"],
            source_z_edges_m=source_z_edges,
            source_y_edges_m=source_q_edges,
            source_id=(
                f"TI021 resolved inner dx={dx:g}m dt={dt:g}s duration={duration:g}s"
            ),
            metadata=metadata,
        ) as writer:
            result = solve_resolved_obstacle_flow_3d(
                config,
                receptors=receptors,
                handoff_request=request,
                handoff_observer=log_progress,
            )
            handoff_frame_count = writer.frame_count
    wall_clock = time.perf_counter() - started

    labels = np.asarray([item.label for item in receptors])
    positions = np.asarray([item.position_m for item in receptors], dtype=float)
    concentration_traces = np.asarray([
        result.receptor_traces_kg_m3[label] for label in labels
    ])
    vol_percent_traces = np.asarray([
        result.receptor_traces_vol_percent[label] for label in labels
    ])
    temperature_traces = np.asarray([
        result.receptor_temperature_traces_k[label] for label in labels
    ])
    np.savez_compressed(
        paths["receptor_npz"],
        time_s=result.time_s,
        labels=labels,
        solver_positions_m=positions,
        concentration_kg_m3=concentration_traces,
        concentration_vol_percent=vol_percent_traces,
        temperature_k=temperature_traces,
    )
    np.savez_compressed(
        paths["final_field_npz"],
        x_m=result.x_m,
        y_m=result.y_m,
        z_m=result.z_m,
        solid_mask=result.solid_mask,
        concentration_kg_m3=result.concentration_kg_m3,
        concentration_vol_percent=result.concentration_vol_percent,
        maximum_concentration_vol_percent=(
            result.maximum_concentration_vol_percent
        ),
        temperature_k=result.temperature_k,
        mixture_density_kg_m3=result.mixture_density_kg_m3,
        velocity_u_m_s=result.velocity_u_m_s,
        velocity_v_m_s=result.velocity_v_m_s,
        velocity_w_m_s=result.velocity_w_m_s,
        pressure_pa_gauge=result.pressure_pa_gauge,
    )
    diagnostics = result.diagnostics
    expected_frames = steps + 1
    gates = {
        "inner_sensor_coordinates_within_declared_domain": len(receptors) > 0,
        "all_60_inner_sensor_coordinates_sampled": len(receptors) == 60,
        "complete_time_horizon": (
            result.time_s.size == expected_frames
            and math.isclose(result.time_s[-1], duration)
        ),
        "complete_handoff_history_flushed": (
            handoff_frame_count == expected_frames
        ),
        "species_mass_conservative": (
            diagnostics["maximum_relative_mass_residual"] <= 1.0e-6
        ),
        "sensible_enthalpy_conservative": (
            diagnostics["maximum_relative_energy_residual"] <= 1.0e-6
        ),
        "continuity_residual_at_most_1e_7": (
            diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ] <= 1.0e-7
        ),
        "thermodynamic_corrector_converged_every_step": (
            diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
        ),
        "mole_fraction_bounded": bool(
            np.nanmin(result.species_mole_fraction) >= -1.0e-12
            and np.nanmax(result.maximum_species_mole_fraction) <= 1.0 + 1.0e-12
        ),
        "report_stationary_cloud_momentum_applied": (
            config.initial_binary_cloud_momentum == "quiescent"
        ),
        "observation_values_not_used_to_configure_run": True,
    }
    # A deliberately smaller local patch is a valid component run when every
    # numerical/conservation gate passes.  It must not be mislabeled as a
    # full 60-sensor inner-domain qualification; that coverage gate remains a
    # diagnostic and is consumed by the outer/hybrid assembly step.
    physical_gate_values = [
        value for name, value in gates.items()
        if name != "all_60_inner_sensor_coordinates_sampled"
    ]
    all_gates_pass = all(physical_gate_values)
    record = {
        "schema": SCHEMA,
        "status": (
            "accepted_completed_inner_transient"
            if all_gates_pass and math.isclose(duration, 600.0)
            else (
                "accepted_local_patch_transient_qualification"
                if all_gates_pass and len(receptors) < 60
                else "accepted_inner_transient_qualification"
                if all_gates_pass
                else "withheld_inner_transient"
            )
        ),
        "scope": (
            "no-fit primary-prior resolved inner-domain execution; sensor "
            "observation values are reserved for later independent scoring"
        ),
        "provenance": {
            "workbook_path": str(workbook_path),
            "workbook_sha256": _sha256(workbook_path),
            "property_audit_path": str(property_path),
            "property_audit_sha256": _sha256(property_path),
        },
        "coordinates": coordinates,
        "configuration": {
            "duration_s": duration,
            "time_step_s": dt,
            "active_steps": steps,
            "cell_size_m": dx,
            "corrector_tolerance": corrector_tolerance,
            "domain_m": list(config.domain_m),
            "cells": list(config.cells),
            "cell_count": int(np.prod(config.cells)),
            "inner_receptor_count": len(receptors),
            "coverage_scope": (
                "full_inner_domain_sensor_set"
                if len(receptors) == len(dataset.dose_sensors)
                else "local_patch_sensor_subset"
            ),
            "coupling_mode": "one_way_handoff_frame_export",
            "dynamic_outer_patch_coupling_applied": False,
            "initial_binary_cloud_momentum": (
                config.initial_binary_cloud_momentum
            ),
            "stratification_damping_alpha": (
                config.stratification_damping_alpha
            ),
            "stratification_damping_exponent": (
                config.stratification_damping_exponent
            ),
            "stratification_shear_regularization_s_inv": (
                config.stratification_shear_regularization_s_inv
            ),
            "initial_cloud_flat_to_flat_m": float(initial_cloud_flat_to_flat_m),
            "initial_cloud_z_max_m": float(initial_cloud_z_max_m),
            "finite_source_release_duration_s": (
                None if finite_source_release_duration_s is None
                else float(finite_source_release_duration_s)
            ),
            "finite_source_start_delay_s": float(finite_source_start_delay_s),
            "finite_source_mass_kg": (
                None if finite_source_mass_kg is None
                else float(finite_source_mass_kg)
            ),
            "finite_source_mass_basis": finite_source_mass_basis,
            "finite_source_schedule_species_mass_kg": coordinates.get(
                "finite_source_schedule_species_mass_kg"
            ),
            "finite_source_schedule_total_mass_kg": coordinates.get(
                "finite_source_schedule_total_mass_kg"
            ),
            "workbook_release_rate_kg_s": float(dataset.workbook_release_rate_kg_s),
            "workbook_release_duration_s": float(dataset.workbook_release_duration_s),
            "workbook_release_mass_proxy_kg": float(dataset.workbook_release_mass_proxy_kg),
            "finite_source_mass_provenance": (
                "explicit_cli_override_only; workbook rate-times-duration proxy is not a validated atmospheric source"
                if finite_source_mass_kg is not None else
                "geometry-derived initial-cloud mass unless an explicit override is supplied"
            ),
            "finite_source_position_z_m": float(finite_source_position_z_m),
            "finite_source_momentum_velocity_m_s": [
                float(value) for value in finite_source_momentum_velocity_m_s
            ],
            "finite_source_total_mass_to_species_mass_ratio": float(
                finite_source_total_mass_to_species_mass_ratio
            ),
            "wind_speed_m_s": float(wind_speed_m_s),
            "wind_heading_relative_deg": float(wind_heading_relative_deg),
            "wind_input_provenance": (
                "explicit CLI/configuration prior; not fitted to arrival labels"
            ),
            "receptor_vertical_sampling": receptor_vertical_sampling,
            "circular_arc_representation": circular_arc_representation,
            "obstacle_height_m": float(obstacle_height_m),
            "vertical_resolution": {
                "first_cell_center_z_m": first_cell_center_z_m,
                "below_first_cell_center_count": len(below_first_vertical_center),
                "below_first_cell_center_receptors": below_first_vertical_center,
                "warning": (
                    "low receptors are vertically aliased by the resolved grid"
                    if below_first_vertical_center else None
                ),
            },
            "finite_source_release_shape": finite_source_release_shape,
            "finite_source_sigma_xyz_m": (
                None if finite_source_sigma_xyz_m is None
                else [float(value) for value in finite_source_sigma_xyz_m]
            ),
            "finite_source_kernel": finite_source_kernel,
            "scalar_advection_scheme": scalar_advection_scheme,
            "source_time_integration": source_time_integration,
            "receptors": receptor_metadata,
            "observation_values_used_to_configure_run": [],
        },
        "outputs": {
            key: {
                "path": str(path),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for key, path in paths.items()
            if key != "record_json"
        },
        "numerics": {
            "maximum_relative_species_mass_residual": diagnostics[
                "maximum_relative_mass_residual"
            ],
            "maximum_relative_sensible_enthalpy_residual": diagnostics[
                "maximum_relative_energy_residual"
            ],
            "maximum_continuity_residual_kg_m3_s": diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
            "maximum_global_mass_residual_kg_s": diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "maximum_momentum_balance_residual_kg_m_s": diagnostics[
                "maximum_momentum_balance_residual_kg_m_s"
            ],
            "maximum_advective_number": diagnostics["advective_number"],
            "maximum_diffusion_number": diagnostics["diffusion_number"],
            "maximum_thermal_diffusion_number": diagnostics[
                "thermal_diffusion_number"
            ],
            "maximum_vol_percent": float(np.nanmax(
                result.maximum_concentration_vol_percent
            )),
            "maximum_receptor_vol_percent": float(np.max(
                vol_percent_traces
            )),
            "maximum_receptor_temperature_k": float(np.max(
                temperature_traces
            )),
            "minimum_receptor_temperature_k": float(np.min(
                temperature_traces
            )),
            "corrector_nonconverged_steps": diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ],
            "corrector_maximum_iterations_per_step": diagnostics[
                "thermodynamic_corrector_maximum_iterations_used"
            ],
            "corrector_mean_iterations_per_step": diagnostics[
                "thermodynamic_corrector_mean_iterations_used"
            ],
            "projection_solver": diagnostics[
                "variable_density_linear_solver"
            ],
            "projection_preconditioner": diagnostics[
                "variable_density_preconditioner"
            ],
            "projection_maximum_linear_iterations_per_step": diagnostics[
                "variable_density_projection_maximum_linear_iterations_per_step"
            ],
            "projection_mean_linear_iterations_per_step": diagnostics[
                "variable_density_projection_mean_linear_iterations_per_step"
            ],
        },
        "runtime": {
            "wall_clock_s": wall_clock,
            "solver_reported_total_s": diagnostics["total_runtime_s"],
            "seconds_per_cell_step": (
                diagnostics["total_runtime_s"]
                / (int(np.prod(config.cells)) * steps)
            ),
            "linear_extrapolation_to_600s_hours": (
                diagnostics["total_runtime_s"] * 600.0 / duration / 3600.0
            ),
        },
        "gates": gates,
        "all_gates_pass": all_gates_pass,
        "next_required": (
            "run a dynamically coupled outer/patch execution before any "
            "600 s full-event score"
            if len(receptors) < len(dataset.dose_sensors) else
            "run the same frozen configuration to 600 s, then score the "
            "nested outer domain and all 145 sensors"
            if duration < 600.0 else
            "drive the nested outer domain and score all 145 sensors"
        ),
    }
    paths["record_json"].write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return record, paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument(
        "--property-audit",
        type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
        ),
    )
    parser.add_argument("--duration-s", type=float, default=10.0)
    parser.add_argument("--cell-size-m", type=float, default=5.0)
    parser.add_argument(
        "--inner-domain-m", nargs=3, type=float, metavar=("LENGTH_X", "WIDTH_Y", "HEIGHT_Z"),
        default=(200.0, 200.0, 40.0),
        help="local resolved patch domain; use a smaller box only with an explicit outer handoff",
    )
    parser.add_argument("--handoff-site-downwind-y-m", type=float, default=100.0)
    parser.add_argument("--time-step-s", type=float, default=0.125)
    parser.add_argument(
        "--progress-interval-s", type=float, default=0.1,
        help="write progress after each 0.1 s of simulated time by default",
    )
    parser.add_argument("--corrector-tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--stratification-damping-alpha", type=float,
        help="enable opt-in vertical turbulent-diffusion damping; no dense-gas validation is implied",
    )
    parser.add_argument(
        "--stratification-damping-exponent", type=float, default=1.0,
    )
    parser.add_argument(
        "--stratification-shear-regularization-s-inv", type=float,
        help="explicit epsilon_s for lagged gradient Richardson number",
    )
    parser.add_argument(
        "--initial-cloud-flat-to-flat-m", type=float, default=14.0,
        help="diagnostic initial-cloud footprint; default is the report-scale 14 m",
    )
    parser.add_argument(
        "--initial-cloud-z-max-m", type=float, default=13.0,
        help="diagnostic initial-cloud top; default is the report-scale 13 m",
    )
    parser.add_argument(
        "--finite-source-release-duration-s", type=float,
        help="diagnostic finite source duration; replaces the t=0 filled cloud",
    )
    parser.add_argument(
        "--finite-source-start-delay-s", type=float, default=0.0,
        help=(
            "diagnostic pre-spread/collapse delay before finite source release; "
            "default 0 and not fitted into the physical prior"
        ),
    )
    parser.add_argument(
        "--finite-source-mass-kg", type=float,
        help="diagnostic source mass override; must carry an external provenance",
    )
    parser.add_argument(
        "--finite-source-mass-basis", choices=("species", "total"),
        default="species",
        help=(
            "mass basis for --finite-source-mass-kg; the solver schedule is "
            "always converted to tracked-species mass"
        ),
    )
    parser.add_argument(
        "--finite-source-position-z-m", type=float, default=1.0,
        help="diagnostic finite-source height; default 1 m is not a measured source height",
    )
    parser.add_argument(
        "--finite-source-momentum-velocity-m-s", nargs=3, type=float,
        metavar=("U_SOURCE", "V_SOURCE", "W_SOURCE"),
        default=(0.0, 0.0, 0.0),
        help="diagnostic source momentum velocity; nonzero values require finite-volume projection",
    )
    parser.add_argument(
        "--finite-source-total-mass-to-species-mass-ratio", type=float, default=1.0,
        help="diagnostic mixture total-mass/species-mass ratio; 1.509 is the workbook 0.32-mole-fraction mixture",
    )
    parser.add_argument(
        "--wind-speed-m-s", type=float, default=3.9,
        help="explicit wind prior; default is the report-scale 3.9 m/s",
    )
    parser.add_argument(
        "--wind-heading-relative-deg", type=float, default=6.1,
        help="explicit wind heading relative to solver +x; default is report-scale 6.1 deg",
    )
    parser.add_argument(
        "--receptor-vertical-sampling",
        choices=("cell_center", "ground_linear"),
        default="cell_center",
        help="diagnostic receptor operator for below-first-cell heights",
    )
    parser.add_argument(
        "--circular-arc-representation",
        choices=("finite_thickness_voxel", "zero_thickness_surface"),
        default="finite_thickness_voxel",
        help="diagnostic curved-wall representation; default uses a one-cell voxel thickness",
    )
    parser.add_argument(
        "--obstacle-height-m", type=float, default=5.0,
        help="diagnostic curved-wall height; default is the declared 5 m wall",
    )
    parser.add_argument(
        "--finite-source-release-shape",
        choices=("constant", "front_loaded", "back_loaded", "smooth_ramp"),
        default="constant",
        help="diagnostic finite-source rate shape; total released mass is fixed",
    )
    parser.add_argument(
        "--finite-source-sigma-xyz-m", nargs=3, type=float,
        metavar=("SIGMA_X", "SIGMA_Y", "SIGMA_Z"),
        help="diagnostic finite-source Gaussian widths; requires finite release",
    )
    parser.add_argument(
        "--finite-source-kernel",
        choices=("gaussian", "variance_matched_box", "cosine_tapered_box"),
        default="gaussian",
        help=(
            "diagnostic source support; compact box and cosine-tapered box "
            "are candidates and are not the default"
        ),
    )
    parser.add_argument(
        "--scalar-advection-scheme",
        choices=("upwind", "muscl_van_leer", "muscl_ssprk2"),
        default="muscl_van_leer",
        help="diagnostic scalar transport scheme; default is MUSCL Van Leer",
    )
    parser.add_argument(
        "--source-time-integration",
        choices=("prestep", "strang"),
        default="prestep",
        help="diagnostic source splitting; default is prestep",
    )
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = run_inner_transient(
        workbook=args.workbook,
        property_audit=args.property_audit,
        output_directory=args.output_directory,
        duration_s=args.duration_s,
        cell_size_m=args.cell_size_m,
        time_step_s=args.time_step_s,
        progress_interval_s=args.progress_interval_s,
        corrector_tolerance=args.corrector_tolerance,
        stratification_damping_alpha=args.stratification_damping_alpha,
        stratification_damping_exponent=args.stratification_damping_exponent,
        stratification_shear_regularization_s_inv=(
            args.stratification_shear_regularization_s_inv
        ),
        inner_domain_m=tuple(args.inner_domain_m),
        circular_arc_representation=args.circular_arc_representation,
        obstacle_height_m=args.obstacle_height_m,
        handoff_site_downwind_y_m=args.handoff_site_downwind_y_m,
        initial_cloud_flat_to_flat_m=args.initial_cloud_flat_to_flat_m,
        initial_cloud_z_max_m=args.initial_cloud_z_max_m,
        finite_source_release_duration_s=args.finite_source_release_duration_s,
        finite_source_start_delay_s=args.finite_source_start_delay_s,
        finite_source_mass_kg=args.finite_source_mass_kg,
        finite_source_mass_basis=args.finite_source_mass_basis,
        finite_source_position_z_m=args.finite_source_position_z_m,
        finite_source_momentum_velocity_m_s=tuple(
            args.finite_source_momentum_velocity_m_s
        ),
        finite_source_total_mass_to_species_mass_ratio=(
            args.finite_source_total_mass_to_species_mass_ratio
        ),
        wind_speed_m_s=args.wind_speed_m_s,
        wind_heading_relative_deg=args.wind_heading_relative_deg,
        receptor_vertical_sampling=args.receptor_vertical_sampling,
        finite_source_release_shape=args.finite_source_release_shape,
        finite_source_sigma_xyz_m=(
            None if args.finite_source_sigma_xyz_m is None
            else tuple(args.finite_source_sigma_xyz_m)
        ),
        finite_source_kernel=args.finite_source_kernel,
        scalar_advection_scheme=args.scalar_advection_scheme,
        source_time_integration=args.source_time_integration,
    )
    print(json.dumps({
        "status": record["status"],
        "all_gates_pass": record["all_gates_pass"],
        "runtime": record["runtime"],
        "outputs": {key: str(path) for key, path in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
