"""Run the TI021 145-sensor nested outer transient with atomic checkpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np

from degali.addons.atmospheric_surface_layer_mixing import NeutralSurfaceLayerMixing
from degali.addons.conservative_handoff_hdf5_3d import ConservativeHandoffHdf5Reader3D
from degali.addons.nested_outer_binary_flow_3d import (
    NESTED_OUTER_BINARY_FLOW_3D_SCHEMA,
    NestedOuterBinaryFlow3DState,
    advance_nested_outer_binary_flow_step_3d,
    initialise_nested_outer_binary_flow_3d,
)
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

try:
    from tools.audit_thorney_ti021_full_transient_plan import site_local_to_nested_outer
    from tools.audit_thorney_ti021_handoff_streaming_preflight import build_configuration
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from audit_thorney_ti021_full_transient_plan import site_local_to_nested_outer
    from audit_thorney_ti021_handoff_streaming_preflight import build_configuration
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-outer-transient.v1"
CHECKPOINT_SCHEMA = "degali.thorney-ti021-outer-checkpoint.v1"
RESTART_SCHEMA = "degali.thorney-ti021-outer-restart.v1"
DEFAULT_DOMAIN_M = (440.0, 640.0, 40.0)
DEFAULT_HANDOFF_DOWNWIND_Y_M = 100.0
OUTER_DOWNWIND_END_SITE_Y_M = 540.0
OUTER_CROSSWIND_ORIGIN_M = 320.0


def _atomic_json(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
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


def _configuration_sha256(value: dict[str, object]) -> str:
    payload = json.dumps(
        value, sort_keys=True, ensure_ascii=False, separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def outer_receptors(
    dataset, *, handoff_downwind_y_m: float = DEFAULT_HANDOFF_DOWNWIND_Y_M,
    domain_m: tuple[float, float, float] = DEFAULT_DOMAIN_M,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, object]]]:
    """Return the coordinate-only sensor subset covered by the nested outer box."""

    labels = []
    positions = []
    metadata = []
    for index, sensor in enumerate(dataset.dose_sensors):
        position = site_local_to_nested_outer(
            sensor.local_x_m, sensor.local_y_m, sensor.height_m,
            handoff_downwind_y_m=handoff_downwind_y_m,
            outer_crosswind_origin_m=OUTER_CROSSWIND_ORIGIN_M,
        )
        if all(0.0 <= coordinate <= length for coordinate, length in zip(position, domain_m)):
            label = f"TI021-S{index + 1:03d}"
            labels.append(label)
            positions.append(position)
            metadata.append({
                "label": label,
                "workbook_row_index_zero_based": index,
                "site_local_xyz_m": [sensor.local_x_m, sensor.local_y_m, sensor.height_m],
                "outer_solver_xyz_m": list(position),
            })
    return np.asarray(labels), np.asarray(positions, dtype=float), metadata


def sample_trilinear_many(
    field: np.ndarray, x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray,
    positions_m: np.ndarray,
) -> np.ndarray:
    """Sample fixed physical receptors with boundary-to-centre clamping."""

    values = np.asarray(field, dtype=float)
    coordinates = tuple(np.asarray(item, dtype=float) for item in (x_m, y_m, z_m))
    positions = np.asarray(positions_m, dtype=float)
    if values.shape != (coordinates[2].size, coordinates[1].size, coordinates[0].size):
        raise ValueError("field shape must match z-y-x coordinates")
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions_m must be an n-by-3 array")
    result = np.empty(positions.shape[0], dtype=float)
    for receptor_index, position in enumerate(positions):
        bounds = []
        fractions = []
        for axis, point in zip(coordinates, position):
            clipped = float(np.clip(point, axis[0], axis[-1]))
            right = int(np.clip(np.searchsorted(axis, clipped), 1, axis.size - 1))
            left = right - 1
            bounds.append((left, right))
            fractions.append((clipped - axis[left]) / (axis[right] - axis[left]))
        sampled = 0.0
        for iz in (0, 1):
            for iy in (0, 1):
                for ix in (0, 1):
                    weight = (
                        (fractions[2] if iz else 1.0 - fractions[2])
                        * (fractions[1] if iy else 1.0 - fractions[1])
                        * (fractions[0] if ix else 1.0 - fractions[0])
                    )
                    sampled += weight * values[
                        bounds[2][iz], bounds[1][iy], bounds[0][ix]
                    ]
        result[receptor_index] = sampled
    return result


def _write_restart(
    path: Path, state: NestedOuterBinaryFlow3DState,
    maximum_concentration_vol_percent: np.ndarray, *, configuration_sha256: str,
) -> None:
    potential = state.correction_potential_kg_m_s
    _atomic_npz(
        path, schema=np.asarray(RESTART_SCHEMA),
        configuration_sha256=np.asarray(configuration_sha256),
        time_s=np.asarray(state.time_s),
        species_mass_concentration_kg_m3=state.species_mass_concentration_kg_m3,
        sensible_enthalpy_density_j_m3=state.sensible_enthalpy_density_j_m3,
        u_m_s=state.u_m_s, v_m_s=state.v_m_s, w_m_s=state.w_m_s,
        u_face_m_s=state.u_face_m_s, v_face_m_s=state.v_face_m_s,
        w_face_m_s=state.w_face_m_s,
        correction_potential_present=np.asarray(potential is not None),
        correction_potential_kg_m_s=(
            np.empty((0,), dtype=float) if potential is None else potential
        ),
        maximum_concentration_vol_percent=maximum_concentration_vol_percent,
        cumulative_species_boundary_outflow_mass_kg=np.asarray(state.cumulative_species_boundary_outflow_mass_kg),
        cumulative_sensible_enthalpy_boundary_outflow_j=np.asarray(state.cumulative_sensible_enthalpy_boundary_outflow_j),
        cumulative_species_mass_residual_kg=np.asarray(state.cumulative_species_mass_residual_kg),
        cumulative_sensible_enthalpy_residual_j=np.asarray(state.cumulative_sensible_enthalpy_residual_j),
        cumulative_prescribed_inlet_species_mass_kg=np.asarray(state.cumulative_prescribed_inlet_species_mass_kg),
        cumulative_prescribed_inlet_sensible_enthalpy_j=np.asarray(state.cumulative_prescribed_inlet_sensible_enthalpy_j),
    )


def _read_restart(path: Path, *, closure, configuration_sha256: str):
    with np.load(path, allow_pickle=False) as source:
        if str(source["schema"].item()) != RESTART_SCHEMA:
            raise ValueError("unsupported outer restart schema")
        if str(source["configuration_sha256"].item()) != configuration_sha256:
            raise ValueError("outer restart configuration SHA-256 mismatch")
        species = np.asarray(source["species_mass_concentration_kg_m3"]).copy()
        enthalpy = np.asarray(source["sensible_enthalpy_density_j_m3"]).copy()
        thermodynamic = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
            species, enthalpy,
        )
        potential = (
            np.asarray(source["correction_potential_kg_m_s"]).copy()
            if bool(source["correction_potential_present"].item()) else None
        )
        state = NestedOuterBinaryFlow3DState(
            schema=NESTED_OUTER_BINARY_FLOW_3D_SCHEMA,
            time_s=float(source["time_s"].item()),
            species_mass_concentration_kg_m3=species,
            sensible_enthalpy_density_j_m3=enthalpy,
            thermodynamic_state=thermodynamic,
            u_m_s=np.asarray(source["u_m_s"]).copy(),
            v_m_s=np.asarray(source["v_m_s"]).copy(),
            w_m_s=np.asarray(source["w_m_s"]).copy(),
            u_face_m_s=np.asarray(source["u_face_m_s"]).copy(),
            v_face_m_s=np.asarray(source["v_face_m_s"]).copy(),
            w_face_m_s=np.asarray(source["w_face_m_s"]).copy(),
            correction_potential_kg_m_s=potential,
            cumulative_species_boundary_outflow_mass_kg=float(source["cumulative_species_boundary_outflow_mass_kg"].item()),
            cumulative_sensible_enthalpy_boundary_outflow_j=float(source["cumulative_sensible_enthalpy_boundary_outflow_j"].item()),
            cumulative_species_mass_residual_kg=float(source["cumulative_species_mass_residual_kg"].item()),
            cumulative_sensible_enthalpy_residual_j=float(source["cumulative_sensible_enthalpy_residual_j"].item()),
            cumulative_prescribed_inlet_species_mass_kg=float(source["cumulative_prescribed_inlet_species_mass_kg"].item()),
            cumulative_prescribed_inlet_sensible_enthalpy_j=float(source["cumulative_prescribed_inlet_sensible_enthalpy_j"].item()),
        )
        maximum = np.asarray(source["maximum_concentration_vol_percent"]).copy()
    if maximum.shape != species.shape:
        raise ValueError("outer restart maximum field shape mismatch")
    return state, maximum


def _segment_diagnostics(values: dict[str, object]) -> dict[str, object]:
    return {
        "step_count": int(values["step_count"]),
        "maximum_continuity_residual_kg_m3_s": float(values["maximum_continuity"]),
        "maximum_momentum_balance_residual_kg_m_s": float(values["maximum_momentum"]),
        "maximum_absolute_species_mass_residual_kg": float(values["maximum_species_residual"]),
        "maximum_absolute_sensible_enthalpy_residual_j": float(values["maximum_energy_residual"]),
        "sum_absolute_species_mass_residual_kg": float(values["sum_absolute_species_residual"]),
        "sum_absolute_sensible_enthalpy_residual_j": float(values["sum_absolute_energy_residual"]),
        "corrector_nonconverged_steps": int(values["nonconverged"]),
        "maximum_corrector_iterations": int(values["maximum_iterations"]),
        "prescribed_inlet_mass_flux_exact_every_step": bool(values["inlet_exact"]),
        "maximum_relative_integral_remap_residual": float(values["maximum_remap_residual"]),
        "lazy_handoff_query_count": int(values["handoff_query_count"]),
        "lazy_source_dataset_frame_read_count": int(values["source_frame_read_count"]),
        "wall_clock_s": float(values["wall_clock_s"]),
    }


def run_outer_transient(
    *, inner_run_directory: Path, output_directory: Path,
    duration_s: float = 600.0, checkpoint_interval_s: float = 10.0,
    cell_size_m: float = 5.0, time_step_s: float = 0.125,
    progress_interval_s: float = 5.0, resume: bool = False,
) -> tuple[dict[str, object], dict[str, Path]]:
    inner = inner_run_directory.resolve()
    output = output_directory.resolve()
    inner_record_path = inner / "ti021-inner-transient.json"
    inner_journal_path = inner / "ti021-inner-checkpoint.json"
    handoff_path = inner / "ti021-inner-handoff.h5"
    inner_record = json.loads(inner_record_path.read_text(encoding="utf-8"))
    inner_journal = json.loads(inner_journal_path.read_text(encoding="utf-8"))
    if not inner_record.get("all_gates_pass") or inner_journal.get("status") != "finalized":
        raise ValueError("outer transient requires a finalized accepted inner run")
    if float(inner_record["configuration"]["duration_s"]) < duration_s:
        raise ValueError("inner handoff history is shorter than outer duration")
    workbook_path = Path(inner_record["provenance"]["workbook_path"]).resolve()
    property_path = Path(inner_record["provenance"]["property_audit_path"]).resolve()
    handoff_downwind_y_m = float(
        inner_record.get("coordinates", {}).get(
            "handoff_site_downwind_y_m",
            inner_record["configuration"].get(
                "handoff_site_downwind_y_m", DEFAULT_HANDOFF_DOWNWIND_Y_M,
            ),
        )
    )
    domain_m = (
        OUTER_DOWNWIND_END_SITE_Y_M - handoff_downwind_y_m,
        640.0,
        40.0,
    )
    if domain_m[0] <= 0.0:
        raise ValueError("handoff must lie upstream of the declared outer domain end")
    dataset = read_thorney_ti021_workbook(workbook_path)
    labels, positions, receptor_metadata = outer_receptors(
        dataset, handoff_downwind_y_m=handoff_downwind_y_m,
        domain_m=domain_m,
    )
    if labels.size == 0:
        raise ValueError("outer domain contains no sensor-coordinate receptors")
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    steps = int(round(duration_s / time_step_s))
    if not math.isclose(steps * time_step_s, duration_s):
        raise ValueError("duration must be an integer number of time steps")
    segment_steps = int(round(checkpoint_interval_s / time_step_s))
    if segment_steps < 1 or not math.isclose(segment_steps * time_step_s, checkpoint_interval_s):
        raise ValueError("checkpoint interval must be an integer number of time steps")
    if any(not math.isclose(length / cell_size_m, round(length / cell_size_m)) for length in domain_m):
        raise ValueError("cell size must divide the outer domain")
    cells = tuple(int(round(length / cell_size_m)) for length in domain_m)
    _inner_config, closure, _coordinates = build_configuration(
        properties, cell_size_m=cell_size_m, time_step_s=time_step_s,
        active_steps=steps, corrector_tolerance=1.0e-6,
        inner_domain_m=tuple(inner_record["configuration"].get(
            "domain_m", (200.0, 200.0, 40.0),
        )),
        handoff_site_downwind_y_m=handoff_downwind_y_m,
    )
    config_record = {
        "duration_s": float(duration_s), "checkpoint_interval_s": float(checkpoint_interval_s),
        "cell_size_m": float(cell_size_m), "time_step_s": float(time_step_s),
        "domain_m": list(domain_m), "cells": list(cells),
        "handoff_site_downwind_y_m": handoff_downwind_y_m,
        "inner_handoff_sha256": _sha256(handoff_path),
        "inner_record_sha256": _sha256(inner_record_path),
        "workbook_sha256": _sha256(workbook_path),
        "property_audit_sha256": _sha256(property_path),
        "projection_linear_solver": "amg_cg",
        "corrector_tolerance": 1.0e-6,
        "observation_values_used_to_configure_run": [],
    }
    config_hash = _configuration_sha256(config_record)
    paths = {
        "record_json": output / "ti021-outer-transient.json",
        "checkpoint_json": output / "ti021-outer-checkpoint.json",
        "progress_jsonl": output / "ti021-outer-progress.jsonl",
        "receptor_npz": output / "ti021-outer-receptor-traces.npz",
        "final_field_npz": output / "ti021-outer-final-field.npz",
    }
    segment_directory = output / "checkpoint-segments"
    restart_directory = output / "checkpoint-restarts"
    if resume:
        if not all(path.exists() for path in (paths["checkpoint_json"], paths["progress_jsonl"], segment_directory, restart_directory)):
            raise FileNotFoundError("resume requires a complete outer checkpoint journal")
        if paths["record_json"].exists():
            raise FileExistsError("outer transient is already finalized")
    else:
        existing = [path for path in (*paths.values(), segment_directory, restart_directory) if path.exists()]
        if existing:
            raise FileExistsError("outer output already exists: " + ", ".join(str(path) for path in existing))
        output.mkdir(parents=True, exist_ok=True)
        segment_directory.mkdir(); restart_directory.mkdir()

    nx, ny, nz = cells
    x_edges = np.linspace(0.0, domain_m[0], nx + 1)
    q_edges = np.linspace(-OUTER_CROSSWIND_ORIGIN_M, OUTER_CROSSWIND_ORIGIN_M, ny + 1)
    y_edges = np.linspace(0.0, domain_m[1], ny + 1)
    z_edges = np.linspace(0.0, domain_m[2], nz + 1)
    x_centres = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centres = 0.5 * (y_edges[:-1] + y_edges[1:])
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    solid = np.zeros((nz, ny, nx), dtype=bool)
    profile_speed = 0.205 / 0.4 * np.log((z_centres + 0.005) / 0.005)
    heading = math.radians(6.1)
    ambient_u = np.broadcast_to((profile_speed * math.cos(heading))[:, None], (nz, ny)).copy()
    ambient_v = np.broadcast_to((profile_speed * math.sin(heading))[:, None], (nz, ny)).copy()
    ambient_w = np.zeros((nz, ny), dtype=float)
    ambient_mass_flux = closure.ambient_density_kg_m3 * ambient_u
    mixing = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=0.205, roughness_length_m=0.005,
        source_id="TI021 u*=0.205 m/s neutral outer transient",
    )
    eddy = mixing.eddy_viscosity_m2_s(z_centres, solid)
    viscosity = 1.5e-5 + eddy
    species_diffusivity = float(properties["solver_policy"]["molecular_species_diffusivity_m2_s"]) + eddy / 0.7
    thermal_diffusivity = float(properties["pure_r12_ambient_properties"]["thermal_diffusivity_m2_s"]) + eddy / 0.74

    if resume:
        journal = json.loads(paths["checkpoint_json"].read_text(encoding="utf-8"))
        if journal.get("schema") != CHECKPOINT_SCHEMA or journal.get("configuration") != config_record:
            raise ValueError("outer resume configuration does not match checkpoint")
        entries = list(journal["segments"])
        last = entries[-1]
        restart_path = output / str(last["restart_file"])
        if _sha256(restart_path) != last["restart_sha256"]:
            raise ValueError("outer restart SHA-256 mismatch")
        state, maximum_field = _read_restart(
            restart_path, closure=closure, configuration_sha256=config_hash,
        )
        completed_time = float(last["stop_time_s"])
    else:
        journal = {"schema": CHECKPOINT_SCHEMA, "status": "running", "configuration": config_record, "segments": []}
        entries = []
        state = initialise_nested_outer_binary_flow_3d(
            cells=cells, closure=closure, ambient_u_m_s=ambient_u,
            ambient_v_m_s=ambient_v, ambient_w_m_s=ambient_w,
            solid_mask=solid,
        )
        maximum_field = np.zeros((nz, ny, nx), dtype=float)
        completed_time = 0.0

    started = time.perf_counter()
    next_progress = completed_time
    progress_mode = "a" if resume else "x"
    with ConservativeHandoffHdf5Reader3D(
        handoff_path, target_z_edges_m=z_edges, target_y_edges_m=q_edges,
    ) as handoff_reader:
        if handoff_reader.time_s[0] > 0.0 or handoff_reader.time_s[-1] < duration_s:
            raise ValueError("durable inner handoff does not cover the outer clock")
        with paths["progress_jsonl"].open(progress_mode, encoding="utf-8") as progress:
            while completed_time < duration_s - 1.0e-12:
                stop_time = min(duration_s, completed_time + checkpoint_interval_s)
                segment_started = time.perf_counter()
                query_count_start = handoff_reader.query_count
                frame_read_count_start = handoff_reader.source_frame_read_count
                remap_residual_start = handoff_reader.maximum_relative_integral_residual
                times = [completed_time]
                initial_vol = 100.0 * state.thermodynamic_state.species_mole_fraction
                traces = [sample_trilinear_many(initial_vol, x_centres, y_centres, z_centres, positions)]
                diagnostics = {
                    "step_count": 0, "maximum_continuity": 0.0, "maximum_momentum": 0.0,
                    "maximum_species_residual": 0.0, "maximum_energy_residual": 0.0,
                    "sum_absolute_species_residual": 0.0,
                    "sum_absolute_energy_residual": 0.0,
                    "nonconverged": 0, "maximum_iterations": 0, "inlet_exact": True,
                    "maximum_remap_residual": remap_residual_start,
                    "handoff_query_count": 0, "source_frame_read_count": 0,
                }
                while state.time_s < stop_time - 1.0e-12:
                    result = advance_nested_outer_binary_flow_step_3d(
                        state, handoff=handoff_reader,
                        ambient_total_mass_flux_kg_m2_s=ambient_mass_flux,
                        closure=closure, solid_mask=solid,
                        spacing_m=(cell_size_m, cell_size_m, cell_size_m),
                        time_step_s=time_step_s,
                        kinematic_viscosity_m2_s=viscosity,
                        species_diffusivity_m2_s=species_diffusivity,
                        thermal_diffusivity_m2_s=thermal_diffusivity,
                        projection_linear_solver="amg_cg",
                        maximum_corrector_iterations=32,
                        minimum_corrector_iterations=2,
                        relative_density_tolerance=1.0e-6,
                        relative_mass_flux_tolerance=1.0e-6,
                    )
                    state = result.state
                    vol = 100.0 * state.thermodynamic_state.species_mole_fraction
                    maximum_field = np.maximum(maximum_field, vol)
                    times.append(state.time_s)
                    traces.append(sample_trilinear_many(vol, x_centres, y_centres, z_centres, positions))
                    corrector = result.binary_corrector
                    diagnostics["step_count"] += 1
                    diagnostics["maximum_continuity"] = max(diagnostics["maximum_continuity"], corrector.mass_projection.maximum_absolute_residual_after_kg_m3_s)
                    diagnostics["maximum_momentum"] = max(diagnostics["maximum_momentum"], float(np.max(np.abs(result.momentum.momentum_balance_residual_kg_m_s))))
                    diagnostics["maximum_species_residual"] = max(diagnostics["maximum_species_residual"], abs(corrector.species_mass_residual_kg))
                    diagnostics["maximum_energy_residual"] = max(diagnostics["maximum_energy_residual"], abs(corrector.sensible_enthalpy_residual_j))
                    diagnostics["sum_absolute_species_residual"] += abs(corrector.species_mass_residual_kg)
                    diagnostics["sum_absolute_energy_residual"] += abs(corrector.sensible_enthalpy_residual_j)
                    diagnostics["nonconverged"] += int(not corrector.converged)
                    diagnostics["maximum_iterations"] = max(diagnostics["maximum_iterations"], corrector.iterations)
                    diagnostics["inlet_exact"] = diagnostics["inlet_exact"] and bool(np.allclose(
                        corrector.mass_projection.mass_flux_x_kg_m2_s[:, :, 0],
                        result.inlet_frame.total_mass_flux_kg_m2_s,
                        rtol=0.0, atol=1.0e-12,
                    ))
                    if state.time_s + 1.0e-12 >= next_progress:
                        item = {
                            "model_time_s": state.time_s,
                            "completed_fraction": state.time_s / duration_s,
                            "elapsed_wall_clock_s_this_process": time.perf_counter() - started,
                            "maximum_vol_percent": float(np.max(maximum_field)),
                            "handoff_query_count": handoff_reader.query_count,
                        }
                        progress.write(json.dumps(item) + "\n"); progress.flush()
                        print(json.dumps(item), flush=True)
                        next_progress += progress_interval_s
                diagnostics["wall_clock_s"] = time.perf_counter() - segment_started
                diagnostics["maximum_remap_residual"] = handoff_reader.maximum_relative_integral_residual
                diagnostics["handoff_query_count"] = handoff_reader.query_count - query_count_start
                diagnostics["source_frame_read_count"] = handoff_reader.source_frame_read_count - frame_read_count_start
                segment_index = len(entries) + 1
                segment_relative = Path("checkpoint-segments") / f"segment-{segment_index:04d}.npz"
                restart_relative = Path("checkpoint-restarts") / f"restart-{segment_index:04d}.npz"
                segment_path = output / segment_relative
                restart_path = output / restart_relative
                _atomic_npz(
                    segment_path, time_s=np.asarray(times), labels=labels,
                    concentration_vol_percent=np.asarray(traces).T,
                    diagnostics_json=np.asarray(json.dumps(_segment_diagnostics(diagnostics))),
                )
                _write_restart(
                    restart_path, state, maximum_field,
                    configuration_sha256=config_hash,
                )
                entry = {
                    "segment_index": segment_index,
                    "start_time_s": completed_time, "stop_time_s": stop_time,
                    "segment_file": segment_relative.as_posix(), "segment_sha256": _sha256(segment_path),
                    "restart_file": restart_relative.as_posix(), "restart_sha256": _sha256(restart_path),
                    "diagnostics": _segment_diagnostics(diagnostics),
                }
                entries.append(entry)
                completed_time = stop_time
                journal = {
                    "schema": CHECKPOINT_SCHEMA,
                    "status": "completed_not_finalized" if math.isclose(completed_time, duration_s) else "running_checkpoint_committed",
                    "configuration": config_record,
                    "last_completed_time_s": completed_time,
                    "segments": entries,
                }
                _atomic_json(paths["checkpoint_json"], journal)
        current_process_remap_residual = handoff_reader.maximum_relative_integral_residual

    all_times = []
    all_traces = []
    diagnostic_records = []
    for index, entry in enumerate(entries):
        segment_path = output / str(entry["segment_file"])
        if _sha256(segment_path) != entry["segment_sha256"]:
            raise ValueError("outer segment SHA-256 mismatch during finalization")
        with np.load(segment_path, allow_pickle=False) as segment:
            keep = slice(None) if index == 0 else slice(1, None)
            all_times.append(np.asarray(segment["time_s"])[keep])
            all_traces.append(np.asarray(segment["concentration_vol_percent"])[:, keep])
        diagnostic_records.append(entry["diagnostics"])
    time_values = np.concatenate(all_times)
    trace_values = np.concatenate(all_traces, axis=1)
    _atomic_npz(
        paths["receptor_npz"], time_s=time_values, labels=labels,
        outer_solver_positions_m=positions,
        concentration_vol_percent=trace_values,
    )
    final_vol = 100.0 * state.thermodynamic_state.species_mole_fraction
    _atomic_npz(
        paths["final_field_npz"], x_m=x_centres, y_m=y_centres, z_m=z_centres,
        crosswind_q_m=0.5 * (q_edges[:-1] + q_edges[1:]), solid_mask=solid,
        concentration_vol_percent=final_vol,
        maximum_concentration_vol_percent=maximum_field,
        temperature_k=state.thermodynamic_state.temperature_k,
        mixture_density_kg_m3=state.thermodynamic_state.mixture_density_kg_m3,
        velocity_u_m_s=state.u_m_s, velocity_v_m_s=state.v_m_s,
        velocity_w_m_s=state.w_m_s,
    )
    maximum_continuity = max(float(item["maximum_continuity_residual_kg_m3_s"]) for item in diagnostic_records)
    maximum_momentum = max(float(item["maximum_momentum_balance_residual_kg_m_s"]) for item in diagnostic_records)
    nonconverged = sum(int(item["corrector_nonconverged_steps"]) for item in diagnostic_records)
    exact_inlet = all(bool(item["prescribed_inlet_mass_flux_exact_every_step"]) for item in diagnostic_records)
    sum_absolute_species_residual = sum(float(item["sum_absolute_species_mass_residual_kg"]) for item in diagnostic_records)
    sum_absolute_energy_residual = sum(float(item["sum_absolute_sensible_enthalpy_residual_j"]) for item in diagnostic_records)
    species_relative = sum_absolute_species_residual / max(abs(state.cumulative_prescribed_inlet_species_mass_kg), 1.0e-30)
    energy_relative = sum_absolute_energy_residual / max(abs(state.cumulative_prescribed_inlet_sensible_enthalpy_j), 1.0e-30)
    remap_residual = max(
        current_process_remap_residual,
        max(float(item["maximum_relative_integral_remap_residual"]) for item in diagnostic_records),
    )
    handoff_queries = sum(int(item["lazy_handoff_query_count"]) for item in diagnostic_records)
    source_frame_reads = sum(int(item["lazy_source_dataset_frame_read_count"]) for item in diagnostic_records)
    gates = {
        "all_declared_outer_sensor_coordinates_sampled": (
            labels.size == len(receptor_metadata) and labels.size > 0
        ),
        "complete_requested_transient_clock": math.isclose(time_values[-1], duration_s) and time_values.size == steps + 1,
        "lazy_handoff_queries_match_steps": handoff_queries == steps,
        "conservative_remap_residual_at_most_1e_10": remap_residual <= 1.0e-10,
        "prescribed_inlet_mass_flux_exact_every_step": exact_inlet,
        "maximum_continuity_residual_at_most_1e_7": maximum_continuity <= 1.0e-7,
        "relative_species_mass_residual_at_most_1e_8": species_relative <= 1.0e-8,
        "relative_sensible_enthalpy_residual_at_most_1e_8": energy_relative <= 1.0e-8,
        "corrector_converged_every_step": nonconverged == 0,
        "observation_values_not_used_to_configure_run": True,
    }
    gates = {key: bool(value) for key, value in gates.items()}
    accepted = all(gates.values())
    total_wall = sum(float(item["wall_clock_s"]) for item in diagnostic_records)
    journal["status"] = "finalized"
    _atomic_json(paths["checkpoint_json"], journal)
    record = {
        "schema": SCHEMA,
        "status": "accepted_completed_outer_transient" if accepted else "withheld_outer_transient",
        "scope": "no-fit conservative nested outer execution; observation values reserved for scoring",
        "provenance": {
            "inner_record_path": str(inner_record_path), "inner_record_sha256": _sha256(inner_record_path),
            "inner_handoff_path": str(handoff_path), "inner_handoff_sha256": _sha256(handoff_path),
            "workbook_path": str(workbook_path), "workbook_sha256": _sha256(workbook_path),
            "property_audit_path": str(property_path), "property_audit_sha256": _sha256(property_path),
        },
        "configuration": {
            **config_record, "cell_count": int(np.prod(cells)), "active_steps": steps,
            "outer_receptor_count": int(labels.size), "receptors": receptor_metadata,
            "coordinate_contract": {
                "solver_x": "site_local_y - 100 m",
                "solver_y": "320 m - site_local_x",
                "handoff_target_crosswind_q_bounds_m": [-320.0, 320.0],
            },
        },
        "checkpointing": {
            "schema": CHECKPOINT_SCHEMA, "segment_count": len(entries),
            "last_completed_time_s": completed_time,
            "latest_restart_path": str(output / str(entries[-1]["restart_file"])),
        },
        "numerics": {
            "maximum_continuity_residual_kg_m3_s": maximum_continuity,
            "maximum_momentum_balance_residual_kg_m_s": maximum_momentum,
            "corrector_nonconverged_steps": nonconverged,
            "cumulative_relative_species_mass_residual": species_relative,
            "cumulative_relative_sensible_enthalpy_residual": energy_relative,
            "sum_absolute_species_mass_residual_kg": sum_absolute_species_residual,
            "sum_absolute_sensible_enthalpy_residual_j": sum_absolute_energy_residual,
            "signed_cumulative_species_mass_residual_kg": state.cumulative_species_mass_residual_kg,
            "signed_cumulative_sensible_enthalpy_residual_j": state.cumulative_sensible_enthalpy_residual_j,
            "maximum_vol_percent": float(np.max(maximum_field)),
            "maximum_receptor_vol_percent": float(np.max(trace_values)),
            "maximum_relative_integral_remap_residual": remap_residual,
            "lazy_handoff_query_count": handoff_queries,
            "lazy_source_dataset_frame_read_count": source_frame_reads,
        },
        "runtime": {
            "wall_clock_s": total_wall,
            "seconds_per_cell_step": total_wall / (int(np.prod(cells)) * steps),
        },
        "outputs": {
            key: {"path": str(path), "sha256": _sha256(path), "size_bytes": path.stat().st_size}
            for key, path in paths.items() if key != "record_json"
        },
        "gates": gates, "all_gates_pass": accepted,
        "next_required": "merge sensor histories, verify interior overlap geometry, and score only when the observation horizon is covered",
    }
    paths["record_json"].write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return record, paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inner-run-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--duration-s", type=float, default=600.0)
    parser.add_argument("--checkpoint-interval-s", type=float, default=10.0)
    parser.add_argument("--cell-size-m", type=float, default=5.0)
    parser.add_argument("--time-step-s", type=float, default=0.125)
    parser.add_argument("--progress-interval-s", type=float, default=5.0)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    record, paths = run_outer_transient(
        inner_run_directory=args.inner_run_directory,
        output_directory=args.output_directory,
        duration_s=args.duration_s, checkpoint_interval_s=args.checkpoint_interval_s,
        cell_size_m=args.cell_size_m, time_step_s=args.time_step_s,
        progress_interval_s=args.progress_interval_s, resume=args.resume,
    )
    print(json.dumps({
        "status": record["status"], "all_gates_pass": record["all_gates_pass"],
        "numerics": record["numerics"], "runtime": record["runtime"],
        "outputs": {key: str(value) for key, value in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
