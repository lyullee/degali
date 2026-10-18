"""Run the TI021 resolved inner transient as durable restartable segments.

Each completed segment atomically commits three mutually checked artifacts:
the conservative solver restart, the segment histories/final fields, and a
small JSON commit journal.  Handoff frames are flushed continuously; after an
interruption, frames newer than the journal are truncated before continuation.
Sensor observation values are never read for configuration or execution.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import time

import numpy as np

from degali.addons.binary_mixture_handoff_plane_3d import (
    BinaryMixtureHandoffPlane3DRequest,
)
from degali.addons.conservative_handoff_hdf5_3d import (
    ConservativeHandoffHdf5Writer3D,
)
from degali.addons.resolved_obstacle_flow_3d import (
    read_resolved_obstacle_flow_3d_restart,
    solve_resolved_obstacle_flow_3d,
    write_resolved_obstacle_flow_3d_restart,
)
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

try:
    from tools.audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )
    from tools.run_thorney_ti021_inner_transient import (
        ROOT, SCHEMA, _handoff_crosswind_edges, _sha256, inner_receptors,
    )
except ModuleNotFoundError:  # Direct ``python tools/<script>.py`` execution.
    from audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )
    from run_thorney_ti021_inner_transient import (
        ROOT, SCHEMA, _handoff_crosswind_edges, _sha256, inner_receptors,
    )


CHECKPOINT_SCHEMA = "degali.thorney-ti021-inner-checkpoint.v1"


def _atomic_json(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _atomic_npz(path: Path, **arrays: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("xb") as stream:
            np.savez_compressed(stream, **arrays)
        os.replace(temporary, path)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def _segment_diagnostics(result) -> dict[str, object]:
    diagnostics = result.diagnostics
    return {
        "step_count": int(result.time_s.size - 1),
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
        "corrector_nonconverged_steps": diagnostics[
            "thermodynamic_corrector_nonconverged_steps"
        ],
        "corrector_maximum_iterations_per_step": diagnostics[
            "thermodynamic_corrector_maximum_iterations_used"
        ],
        "corrector_mean_iterations_per_step": diagnostics[
            "thermodynamic_corrector_mean_iterations_used"
        ],
        "projection_solver": diagnostics["variable_density_linear_solver"],
        "projection_preconditioner": diagnostics[
            "variable_density_preconditioner"
        ],
        "projection_maximum_linear_iterations_per_step": diagnostics[
            "variable_density_projection_maximum_linear_iterations_per_step"
        ],
        "projection_mean_linear_iterations_per_step": diagnostics[
            "variable_density_projection_mean_linear_iterations_per_step"
        ],
        "embedded_arc_wall_stress_closure": diagnostics[
            "embedded_arc_wall_stress_closure"
        ],
        "maximum_embedded_arc_wall_shear_stress_pa": diagnostics[
            "maximum_embedded_arc_wall_shear_stress_pa"
        ],
        "maximum_embedded_arc_wall_friction_velocity_m_s": diagnostics[
            "maximum_embedded_arc_wall_friction_velocity_m_s"
        ],
        "maximum_embedded_arc_wall_tangential_speed_m_s": diagnostics[
            "maximum_embedded_arc_wall_tangential_speed_m_s"
        ],
        "cumulative_embedded_arc_wall_impulse_x_n_s": diagnostics[
            "cumulative_embedded_arc_wall_impulse_x_n_s"
        ],
        "cumulative_embedded_arc_wall_impulse_y_n_s": diagnostics[
            "cumulative_embedded_arc_wall_impulse_y_n_s"
        ],
        "cumulative_embedded_arc_wall_impulse_z_n_s": diagnostics[
            "cumulative_embedded_arc_wall_impulse_z_n_s"
        ],
        "solver_reported_total_s": diagnostics["total_runtime_s"],
    }


def _write_segment(path: Path, result, labels: np.ndarray) -> None:
    assert result.restart_state is not None
    _atomic_npz(
        path,
        time_s=result.time_s,
        labels=labels,
        concentration_kg_m3=np.asarray([
            result.receptor_traces_kg_m3[label] for label in labels
        ]),
        concentration_vol_percent=np.asarray([
            result.receptor_traces_vol_percent[label] for label in labels
        ]),
        receptor_temperature_k=np.asarray([
            result.receptor_temperature_traces_k[label] for label in labels
        ]),
        injected_mass_kg=result.injected_mass_kg,
        domain_mass_kg=result.domain_mass_kg,
        outflow_mass_kg=result.outflow_mass_kg,
        injected_sensible_enthalpy_j=result.injected_sensible_enthalpy_j,
        domain_sensible_enthalpy_j=result.domain_sensible_enthalpy_j,
        outflow_sensible_enthalpy_j=result.outflow_sensible_enthalpy_j,
        x_m=result.x_m,
        y_m=result.y_m,
        z_m=result.z_m,
        solid_mask=result.solid_mask,
        final_concentration_kg_m3=result.concentration_kg_m3,
        final_concentration_vol_percent=result.concentration_vol_percent,
        maximum_concentration_vol_percent=(
            result.maximum_concentration_vol_percent
        ),
        final_temperature_k=result.temperature_k,
        final_mixture_density_kg_m3=result.mixture_density_kg_m3,
        final_velocity_u_m_s=result.velocity_u_m_s,
        final_velocity_v_m_s=result.velocity_v_m_s,
        final_velocity_w_m_s=result.velocity_w_m_s,
        final_pressure_pa_gauge=result.pressure_pa_gauge,
        diagnostics_json=np.asarray(json.dumps(_segment_diagnostics(result))),
    )


def _read_committed_segments(
    output: Path, entries: list[dict[str, object]],
) -> list[dict[str, object]]:
    segments: list[dict[str, object]] = []
    for entry in entries:
        segment_path = output / str(entry["segment_file"])
        restart_path = output / str(entry["restart_file"])
        if _sha256(segment_path) != entry["segment_sha256"]:
            raise ValueError(f"segment SHA-256 mismatch: {segment_path}")
        if _sha256(restart_path) != entry["restart_sha256"]:
            raise ValueError(f"restart SHA-256 mismatch: {restart_path}")
        with np.load(segment_path, allow_pickle=False) as payload:
            segments.append({
                key: np.asarray(payload[key]).copy() for key in payload.files
            })
    return segments


def _weighted_mean(
    entries: list[dict[str, object]], key: str,
) -> float | None:
    weighted = 0.0
    count = 0
    for entry in entries:
        diagnostics = entry["diagnostics"]
        value = diagnostics[key]
        steps = int(diagnostics["step_count"])
        if value is not None:
            weighted += float(value) * steps
            count += steps
    return None if count == 0 else weighted / count


def run_checkpointed_inner_transient(
    *,
    workbook: Path,
    property_audit: Path,
    output_directory: Path,
    duration_s: float,
    checkpoint_interval_s: float,
    cell_size_m: float = 2.5,
    time_step_s: float = 0.0625,
    progress_interval_s: float = 0.1,
    corrector_tolerance: float = 1.0e-4,
    circular_arc_representation: str = "finite_thickness_voxel",
    domain_m: tuple[float, float, float] = (200.0, 200.0, 40.0),
    handoff_site_downwind_y_m: float = 100.0,
    resume: bool = False,
) -> tuple[dict[str, object], dict[str, Path]]:
    """Execute or resume a journaled no-fit TI021 inner-domain run."""

    duration = float(duration_s)
    dt = float(time_step_s)
    checkpoint_interval = float(checkpoint_interval_s)
    for value, label in (
        (duration / dt, "duration"),
        (checkpoint_interval / dt, "checkpoint interval"),
    ):
        if value <= 0.0 or not math.isclose(
            value, round(value), rel_tol=0.0, abs_tol=1.0e-10,
        ):
            raise ValueError(
                f"{label} must be a positive integer multiple of time step"
            )
    if progress_interval_s <= 0.0:
        raise ValueError("progress_interval_s must be positive")
    steps = int(round(duration / dt))
    workbook_path = workbook.resolve()
    property_path = property_audit.resolve()
    output = output_directory.resolve()
    paths = {
        "handoff_hdf5": output / "ti021-inner-handoff.h5",
        "receptor_npz": output / "ti021-inner-receptor-traces.npz",
        "final_field_npz": output / "ti021-inner-final-field.npz",
        "record_json": output / "ti021-inner-transient.json",
        "progress_jsonl": output / "ti021-inner-progress.jsonl",
        "checkpoint_json": output / "ti021-inner-checkpoint.json",
    }
    segments_directory = output / "checkpoint-segments"
    restarts_directory = output / "checkpoint-restarts"
    if resume:
        required = (
            paths["handoff_hdf5"], paths["progress_jsonl"],
            paths["checkpoint_json"], segments_directory, restarts_directory,
        )
        if not all(path.exists() for path in required):
            raise FileNotFoundError("resume requires a complete checkpoint journal")
        if paths["record_json"].exists():
            raise FileExistsError("the checkpointed run is already finalized")
    else:
        existing = [path for path in (*paths.values(), segments_directory, restarts_directory) if path.exists()]
        if existing:
            raise FileExistsError(
                "checkpointed inner output already exists: "
                + ", ".join(str(path) for path in existing)
            )
        output.mkdir(parents=True, exist_ok=True)
        segments_directory.mkdir()
        restarts_directory.mkdir()

    dataset = read_thorney_ti021_workbook(workbook_path)
    receptors, receptor_metadata = inner_receptors(dataset, domain_m=domain_m)
    labels = np.asarray([item.label for item in receptors])
    positions = np.asarray([item.position_m for item in receptors], dtype=float)
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _closure, coordinates = build_configuration(
        properties,
        cell_size_m=cell_size_m,
        time_step_s=dt,
        active_steps=steps,
        corrector_tolerance=corrector_tolerance,
        circular_arc_representation=circular_arc_representation,
        inner_domain_m=domain_m,
        handoff_site_downwind_y_m=handoff_site_downwind_y_m,
    )
    dx = config.domain_m[0] / config.cells[0]
    source_z_edges = np.linspace(0.0, config.domain_m[2], config.cells[2] + 1)
    source_q_edges = _handoff_crosswind_edges(
        config.domain_m[1], config.cells[1],
    )
    source_id = f"TI021 resolved inner dx={dx:g}m dt={dt:g}s duration={duration:g}s"
    configuration_record = {
        "duration_s": duration,
        "time_step_s": dt,
        "checkpoint_interval_s": checkpoint_interval,
        "cell_size_m": dx,
        "corrector_tolerance": corrector_tolerance,
        "circular_arc_representation": config.circular_arc_representation,
        "domain_m": list(config.domain_m),
        "handoff_site_downwind_y_m": float(handoff_site_downwind_y_m),
        "handoff_inner_solver_x_m": float(coordinates["handoff_inner_solver_x_m"]),
        "inner_crosswind_origin_m": 0.5 * float(config.domain_m[1]),
        "cells": list(config.cells),
        "workbook_sha256": _sha256(workbook_path),
        "property_audit_sha256": _sha256(property_path),
    }
    writer_metadata = {
        "case": "Thorney Island trial 21",
        "run_schema": SCHEMA,
        "execution": "atomic segmented restart",
        **configuration_record,
        "inner_receptor_count": len(receptors),
        "observation_values_used_to_configure_run": [],
    }

    if resume:
        journal = json.loads(paths["checkpoint_json"].read_text(encoding="utf-8"))
        if journal.get("schema") != CHECKPOINT_SCHEMA:
            raise ValueError("unsupported TI021 checkpoint journal schema")
        if journal.get("configuration") != configuration_record:
            raise ValueError("resume configuration does not match checkpoint")
        entries = list(journal["segments"])
        if not entries:
            raise ValueError("checkpoint journal has no committed segment")
        latest_restart = output / str(entries[-1]["restart_file"])
        if _sha256(latest_restart) != entries[-1]["restart_sha256"]:
            raise ValueError("latest restart SHA-256 does not match journal")
        restart_state = read_resolved_obstacle_flow_3d_restart(latest_restart)
        completed_time = float(entries[-1]["stop_time_s"])
        if not math.isclose(restart_state.time_s, completed_time):
            raise ValueError("restart time does not match checkpoint journal")
    else:
        entries = []
        restart_state = None
        completed_time = 0.0

    progress_mode = "a" if resume else "x"
    writer_resume = resume
    started = time.perf_counter()
    next_progress_time = (
        completed_time if resume else 0.0
    )
    with paths["progress_jsonl"].open(
        progress_mode, encoding="utf-8",
    ) as progress_stream:
        with ConservativeHandoffHdf5Writer3D(
            paths["handoff_hdf5"],
            source_z_edges_m=source_z_edges,
            source_y_edges_m=source_q_edges,
            source_id=source_id,
            metadata=writer_metadata,
            resume=writer_resume,
        ) as writer:
            if resume:
                removed = writer.truncate_after(completed_time)
                progress_stream.write(json.dumps({
                    "event": "resume",
                    "model_time_s": completed_time,
                    "discarded_uncommitted_handoff_frames": removed,
                }) + "\n")
                progress_stream.flush()

            while completed_time < duration - 1.0e-12:
                stop_time = min(duration, completed_time + checkpoint_interval)
                segment_started = time.perf_counter()

                def log_progress(frame) -> None:
                    nonlocal next_progress_time
                    writer(frame)
                    if frame.time_s + 1.0e-12 >= next_progress_time:
                        entry = {
                            "model_time_s": frame.time_s,
                            "elapsed_wall_clock_s_this_process": (
                                time.perf_counter() - started
                            ),
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

                request = BinaryMixtureHandoffPlane3DRequest(
                    x_face_m=float(coordinates["handoff_inner_solver_x_m"]),
                    sample_interval_s=dt,
                    include_initial=(completed_time == 0.0),
                )
                result = solve_resolved_obstacle_flow_3d(
                    config,
                    receptors=receptors,
                    handoff_request=request,
                    handoff_observer=log_progress,
                    restart_state=restart_state,
                    stop_time_s=stop_time,
                )
                assert result.restart_state is not None
                segment_index = len(entries) + 1
                segment_relative = Path("checkpoint-segments") / (
                    f"segment-{segment_index:04d}.npz"
                )
                restart_relative = Path("checkpoint-restarts") / (
                    f"restart-{segment_index:04d}.npz"
                )
                segment_path = output / segment_relative
                restart_path = output / restart_relative
                _write_segment(segment_path, result, labels)
                write_resolved_obstacle_flow_3d_restart(
                    restart_path, result.restart_state,
                )
                diagnostics = _segment_diagnostics(result)
                entry = {
                    "segment_index": segment_index,
                    "start_time_s": completed_time,
                    "stop_time_s": stop_time,
                    "segment_file": segment_relative.as_posix(),
                    "segment_sha256": _sha256(segment_path),
                    "restart_file": restart_relative.as_posix(),
                    "restart_sha256": _sha256(restart_path),
                    "wall_clock_s": time.perf_counter() - segment_started,
                    "diagnostics": diagnostics,
                }
                entries.append(entry)
                journal = {
                    "schema": CHECKPOINT_SCHEMA,
                    "status": (
                        "completed_not_finalized"
                        if math.isclose(stop_time, duration)
                        else "running_checkpoint_committed"
                    ),
                    "configuration": configuration_record,
                    "last_completed_time_s": stop_time,
                    "segments": entries,
                }
                _atomic_json(paths["checkpoint_json"], journal)
                restart_state = result.restart_state
                completed_time = stop_time
            handoff_frame_count = writer.frame_count

    committed = _read_committed_segments(output, entries)
    times: list[np.ndarray] = []
    concentration: list[np.ndarray] = []
    vol_percent: list[np.ndarray] = []
    receptor_temperature: list[np.ndarray] = []
    domain_mass: list[np.ndarray] = []
    cumulative_injected: list[np.ndarray] = []
    cumulative_outflow: list[np.ndarray] = []
    domain_energy: list[np.ndarray] = []
    cumulative_injected_energy: list[np.ndarray] = []
    cumulative_outflow_energy: list[np.ndarray] = []
    injected_offset = outflow_offset = 0.0
    injected_energy_offset = outflow_energy_offset = 0.0
    for index, segment in enumerate(committed):
        keep = slice(None) if index == 0 else slice(1, None)
        times.append(segment["time_s"][keep])
        concentration.append(segment["concentration_kg_m3"][:, keep])
        vol_percent.append(segment["concentration_vol_percent"][:, keep])
        receptor_temperature.append(segment["receptor_temperature_k"][:, keep])
        domain_mass.append(segment["domain_mass_kg"][keep])
        cumulative_injected.append(
            injected_offset + segment["injected_mass_kg"][keep]
        )
        cumulative_outflow.append(
            outflow_offset + segment["outflow_mass_kg"][keep]
        )
        domain_energy.append(segment["domain_sensible_enthalpy_j"][keep])
        cumulative_injected_energy.append(
            injected_energy_offset
            + segment["injected_sensible_enthalpy_j"][keep]
        )
        cumulative_outflow_energy.append(
            outflow_energy_offset
            + segment["outflow_sensible_enthalpy_j"][keep]
        )
        injected_offset += float(segment["injected_mass_kg"][-1])
        outflow_offset += float(segment["outflow_mass_kg"][-1])
        injected_energy_offset += float(
            segment["injected_sensible_enthalpy_j"][-1]
        )
        outflow_energy_offset += float(
            segment["outflow_sensible_enthalpy_j"][-1]
        )
    time_values = np.concatenate(times)
    concentration_values = np.concatenate(concentration, axis=1)
    vol_percent_values = np.concatenate(vol_percent, axis=1)
    temperature_values = np.concatenate(receptor_temperature, axis=1)
    domain_mass_values = np.concatenate(domain_mass)
    injected_values = np.concatenate(cumulative_injected)
    outflow_values = np.concatenate(cumulative_outflow)
    domain_energy_values = np.concatenate(domain_energy)
    injected_energy_values = np.concatenate(cumulative_injected_energy)
    outflow_energy_values = np.concatenate(cumulative_outflow_energy)
    initial_mass = float(domain_mass_values[0])
    mass_residual = (
        domain_mass_values + outflow_values - initial_mass - injected_values
    )
    initial_energy = float(domain_energy_values[0])
    energy_residual = (
        domain_energy_values + outflow_energy_values
        - initial_energy - injected_energy_values
    )
    last = committed[-1]
    np.savez_compressed(
        paths["receptor_npz"],
        time_s=time_values,
        labels=labels,
        solver_positions_m=positions,
        concentration_kg_m3=concentration_values,
        concentration_vol_percent=vol_percent_values,
        temperature_k=temperature_values,
    )
    np.savez_compressed(
        paths["final_field_npz"],
        x_m=last["x_m"], y_m=last["y_m"], z_m=last["z_m"],
        solid_mask=last["solid_mask"],
        concentration_kg_m3=last["final_concentration_kg_m3"],
        concentration_vol_percent=last["final_concentration_vol_percent"],
        maximum_concentration_vol_percent=(
            last["maximum_concentration_vol_percent"]
        ),
        temperature_k=last["final_temperature_k"],
        mixture_density_kg_m3=last["final_mixture_density_kg_m3"],
        velocity_u_m_s=last["final_velocity_u_m_s"],
        velocity_v_m_s=last["final_velocity_v_m_s"],
        velocity_w_m_s=last["final_velocity_w_m_s"],
        pressure_pa_gauge=last["final_pressure_pa_gauge"],
    )

    diagnostic_entries = [entry["diagnostics"] for entry in entries]
    def maximum(key: str) -> float:
        return max(float(item[key]) for item in diagnostic_entries)

    maximum_relative_mass = float(np.max(np.abs(mass_residual))) / max(
        abs(initial_mass + injected_values[-1]), 1.0e-30,
    )
    maximum_relative_energy = float(np.max(np.abs(energy_residual))) / max(
        abs(initial_energy + injected_energy_values[-1]), 1.0e-30,
    )
    nonconverged = sum(
        int(item["corrector_nonconverged_steps"])
        for item in diagnostic_entries
    )
    expected_frames = steps + 1
    gates = {
        "all_declared_inner_sensor_coordinates_sampled": (
            len(receptors) == len(receptor_metadata) and len(receptors) > 0
        ),
        "complete_time_horizon": (
            time_values.size == expected_frames
            and math.isclose(float(time_values[-1]), duration)
        ),
        "complete_handoff_history_flushed": handoff_frame_count == expected_frames,
        "species_mass_conservative": maximum_relative_mass <= 1.0e-6,
        "sensible_enthalpy_conservative": maximum_relative_energy <= 1.0e-6,
        "continuity_residual_at_most_1e_7": (
            maximum("maximum_continuity_residual_kg_m3_s") <= 1.0e-7
        ),
        "thermodynamic_corrector_converged_every_step": nonconverged == 0,
        "mole_fraction_bounded": bool(
            np.nanmin(last["final_concentration_vol_percent"]) >= -1.0e-10
            and np.nanmax(last["maximum_concentration_vol_percent"])
            <= 100.0 + 1.0e-10
        ),
        "report_stationary_cloud_momentum_applied": (
            config.initial_binary_cloud_momentum == "quiescent"
        ),
        "embedded_obstacle_tangential_no_slip_stress_resolved": (
            config.circular_arc_representation != "zero_thickness_surface"
            or config.active_momentum_solver in {
                "d3q19_bgk", "finite_volume_projection",
            }
        ),
        "embedded_arc_wall_stress_diagnostics_complete": (
            config.circular_arc_representation != "zero_thickness_surface"
            or config.active_momentum_solver != "finite_volume_projection"
            or all(
                "maximum_embedded_arc_wall_shear_stress_pa" in item
                for item in diagnostic_entries
            )
        ),
        "observation_values_not_used_to_configure_run": True,
        "checkpoint_journal_complete": math.isclose(completed_time, duration),
    }
    gates = {key: bool(value) for key, value in gates.items()}
    all_gates_pass = all(gates.values())
    total_wall_clock = sum(float(entry["wall_clock_s"]) for entry in entries)
    total_solver_runtime = sum(
        float(item["solver_reported_total_s"])
        for item in diagnostic_entries
    )
    latest_restart = output / str(entries[-1]["restart_file"])
    journal["status"] = "finalized"
    _atomic_json(paths["checkpoint_json"], journal)
    record = {
        "schema": SCHEMA,
        "status": (
            "accepted_completed_inner_transient"
            if all_gates_pass and math.isclose(duration, 600.0)
            else (
                "accepted_inner_transient_qualification"
                if all_gates_pass else "withheld_inner_transient"
            )
        ),
        "scope": (
            "no-fit primary-prior resolved inner-domain execution with atomic "
            "segment restarts; sensor observation values are reserved for scoring"
        ),
        "provenance": {
            "workbook_path": str(workbook_path),
            "workbook_sha256": _sha256(workbook_path),
            "property_audit_path": str(property_path),
            "property_audit_sha256": _sha256(property_path),
        },
        "coordinates": coordinates,
        "configuration": {
            **configuration_record,
            "active_steps": steps,
            "domain_m": list(config.domain_m),
            "cell_count": int(np.prod(config.cells)),
            "inner_receptor_count": len(receptors),
            "initial_binary_cloud_momentum": config.initial_binary_cloud_momentum,
            "receptors": receptor_metadata,
            "observation_values_used_to_configure_run": [],
        },
        "obstacle_wall_model": {
            "representation": config.circular_arc_representation,
            "active_momentum_solver": config.active_momentum_solver,
            "tangential_no_slip_stress_resolved": gates[
                "embedded_obstacle_tangential_no_slip_stress_resolved"
            ],
            "interpretation": (
                "smooth-wall Reichardt immersed stress with blocked normal "
                "faces; roughness sensitivity and separated-wake validation "
                "are still required"
                if config.circular_arc_representation == "zero_thickness_surface"
                and config.active_momentum_solver == "finite_volume_projection"
                else "wall shear boundary is represented by the selected solver"
            ),
        },
        "embedded_arc_wall_stress_summary": {
            "closure": (
                "smooth_wall_reichardt"
                if (
                    config.circular_arc_representation == "zero_thickness_surface"
                    and config.active_momentum_solver == "finite_volume_projection"
                ) else None
            ),
            "maximum_shear_stress_pa": max(
                (
                    float(item.get("maximum_embedded_arc_wall_shear_stress_pa", 0.0))
                    for item in diagnostic_entries
                ),
                default=0.0,
            ),
            "maximum_friction_velocity_m_s": max(
                (
                    float(item.get("maximum_embedded_arc_wall_friction_velocity_m_s", 0.0))
                    for item in diagnostic_entries
                ),
                default=0.0,
            ),
            "cumulative_wall_impulse_n_s": [
                sum(float(item.get(f"cumulative_embedded_arc_wall_impulse_{axis}_n_s", 0.0))
                    for item in diagnostic_entries)
                for axis in ("x", "y", "z")
            ],
            "diagnostics_complete": gates[
                "embedded_arc_wall_stress_diagnostics_complete"
            ],
        },
        "checkpointing": {
            "schema": CHECKPOINT_SCHEMA,
            "segment_count": len(entries),
            "last_completed_time_s": completed_time,
            "latest_restart_path": str(latest_restart),
            "latest_restart_sha256": _sha256(latest_restart),
            "recovery_policy": (
                "truncate HDF5 frames after the last atomically journaled restart"
            ),
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
            "maximum_relative_species_mass_residual": maximum_relative_mass,
            "maximum_relative_sensible_enthalpy_residual": maximum_relative_energy,
            "maximum_continuity_residual_kg_m3_s": maximum(
                "maximum_continuity_residual_kg_m3_s"
            ),
            "maximum_global_mass_residual_kg_s": maximum(
                "maximum_global_mass_residual_kg_s"
            ),
            "maximum_momentum_balance_residual_kg_m_s": maximum(
                "maximum_momentum_balance_residual_kg_m_s"
            ),
            "maximum_advective_number": maximum("maximum_advective_number"),
            "maximum_diffusion_number": maximum("maximum_diffusion_number"),
            "maximum_thermal_diffusion_number": maximum(
                "maximum_thermal_diffusion_number"
            ),
            "maximum_vol_percent": float(np.nanmax(
                last["maximum_concentration_vol_percent"]
            )),
            "maximum_receptor_vol_percent": float(np.max(vol_percent_values)),
            "maximum_receptor_temperature_k": float(np.max(temperature_values)),
            "minimum_receptor_temperature_k": float(np.min(temperature_values)),
            "corrector_nonconverged_steps": nonconverged,
            "corrector_maximum_iterations_per_step": max(
                int(item["corrector_maximum_iterations_per_step"])
                for item in diagnostic_entries
            ),
            "corrector_mean_iterations_per_step": _weighted_mean(
                entries, "corrector_mean_iterations_per_step"
            ),
            "projection_solver": diagnostic_entries[-1]["projection_solver"],
            "projection_preconditioner": diagnostic_entries[-1][
                "projection_preconditioner"
            ],
            "projection_maximum_linear_iterations_per_step": max(
                int(item["projection_maximum_linear_iterations_per_step"])
                for item in diagnostic_entries
            ),
            "projection_mean_linear_iterations_per_step": _weighted_mean(
                entries, "projection_mean_linear_iterations_per_step"
            ),
        },
        "runtime": {
            "wall_clock_s": total_wall_clock,
            "solver_reported_total_s": total_solver_runtime,
            "seconds_per_cell_step": (
                total_solver_runtime / (int(np.prod(config.cells)) * steps)
            ),
            "linear_extrapolation_to_600s_hours": (
                total_solver_runtime * 600.0 / duration / 3600.0
            ),
        },
        "gates": gates,
        "all_gates_pass": all_gates_pass,
        "next_required": (
            "run the same frozen configuration to 600 s"
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
        "--property-audit", type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
        ),
    )
    parser.add_argument("--duration-s", type=float, default=600.0)
    parser.add_argument("--checkpoint-interval-s", type=float, default=10.0)
    parser.add_argument("--cell-size-m", type=float, default=2.5)
    parser.add_argument("--time-step-s", type=float, default=0.0625)
    parser.add_argument(
        "--progress-interval-s", type=float, default=0.1,
        help="write progress after each 0.1 s of simulated time by default",
    )
    parser.add_argument("--corrector-tolerance", type=float, default=1.0e-4)
    parser.add_argument(
        "--circular-arc-representation",
        choices=("finite_thickness_voxel", "zero_thickness_surface"),
        default="finite_thickness_voxel",
    )
    parser.add_argument(
        "--domain-m", nargs=3, type=float, default=(200.0, 200.0, 40.0),
        metavar=("X", "Y", "Z"),
    )
    parser.add_argument("--handoff-site-downwind-y-m", type=float, default=100.0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = run_checkpointed_inner_transient(
        workbook=args.workbook,
        property_audit=args.property_audit,
        output_directory=args.output_directory,
        duration_s=args.duration_s,
        checkpoint_interval_s=args.checkpoint_interval_s,
        cell_size_m=args.cell_size_m,
        time_step_s=args.time_step_s,
        progress_interval_s=args.progress_interval_s,
        corrector_tolerance=args.corrector_tolerance,
        circular_arc_representation=args.circular_arc_representation,
        domain_m=tuple(args.domain_m),
        handoff_site_downwind_y_m=args.handoff_site_downwind_y_m,
        resume=args.resume,
    )
    print(json.dumps({
        "status": record["status"],
        "all_gates_pass": record["all_gates_pass"],
        "checkpointing": record["checkpointing"],
        "runtime": record["runtime"],
        "outputs": {key: str(path) for key, path in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
