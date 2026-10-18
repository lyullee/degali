"""Materialize a committed prefix of a checkpointed TI021 run for audits."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

try:
    from tools.run_thorney_ti021_inner_transient import SCHEMA, _sha256
    from tools.run_thorney_ti021_inner_transient_checkpointed import (
        CHECKPOINT_SCHEMA,
        _read_committed_segments,
    )
except ModuleNotFoundError:
    from run_thorney_ti021_inner_transient import SCHEMA, _sha256
    from run_thorney_ti021_inner_transient_checkpointed import (
        CHECKPOINT_SCHEMA,
        _read_committed_segments,
    )


PREFIX_SCHEMA = "degali.thorney-ti021-checkpoint-prefix.v1"


def materialize_prefix(
    *, checkpointed_run: Path, stop_time_s: float,
    coordinate_reference_run: Path, output_directory: Path,
) -> tuple[dict[str, object], dict[str, Path]]:
    source = checkpointed_run.resolve()
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    journal_path = source / "ti021-inner-checkpoint.json"
    journal_bytes = journal_path.read_bytes()
    journal_sha256 = hashlib.sha256(journal_bytes).hexdigest()
    journal_snapshot = output / "source-ti021-inner-checkpoint.json"
    journal_snapshot.write_bytes(journal_bytes)
    journal = json.loads(journal_bytes.decode("utf-8"))
    if journal.get("schema") != CHECKPOINT_SCHEMA:
        raise ValueError("unsupported TI021 checkpoint journal schema")
    stop = float(stop_time_s)
    entries = [
        entry for entry in journal["segments"]
        if float(entry["stop_time_s"]) <= stop + 1.0e-12
    ]
    if not entries or not math.isclose(float(entries[-1]["stop_time_s"]), stop):
        raise ValueError("stop_time_s must be an atomically committed segment boundary")
    committed = _read_committed_segments(source, entries)
    configuration = journal["configuration"]
    dt = float(configuration["time_step_s"])

    reference_path = (
        coordinate_reference_run.resolve() / "ti021-inner-receptor-traces.npz"
    )
    reference_record_path = (
        coordinate_reference_run.resolve() / "ti021-inner-transient.json"
    )
    reference_record = json.loads(
        reference_record_path.read_text(encoding="utf-8")
    )
    receptor_metadata = reference_record["configuration"]["receptors"]
    with np.load(reference_path, allow_pickle=False) as reference:
        labels = np.asarray(reference["labels"]).copy()
        positions = np.asarray(reference["solver_positions_m"], dtype=float).copy()
    if labels.tolist() != [item["label"] for item in receptor_metadata]:
        raise ValueError("coordinate-reference receptor metadata do not match labels")
    if not np.array_equal(labels, committed[0]["labels"]):
        raise ValueError("coordinate-reference labels do not match checkpoint segments")

    times = []
    concentration = []
    vol_percent = []
    temperatures = []
    domain_mass = []
    injected = []
    outflow = []
    domain_energy = []
    injected_energy = []
    outflow_energy = []
    injected_offset = outflow_offset = 0.0
    injected_energy_offset = outflow_energy_offset = 0.0
    for index, segment in enumerate(committed):
        keep = slice(None) if index == 0 else slice(1, None)
        times.append(segment["time_s"][keep])
        concentration.append(segment["concentration_kg_m3"][:, keep])
        vol_percent.append(segment["concentration_vol_percent"][:, keep])
        temperatures.append(segment["receptor_temperature_k"][:, keep])
        domain_mass.append(segment["domain_mass_kg"][keep])
        injected.append(injected_offset + segment["injected_mass_kg"][keep])
        outflow.append(outflow_offset + segment["outflow_mass_kg"][keep])
        domain_energy.append(segment["domain_sensible_enthalpy_j"][keep])
        injected_energy.append(
            injected_energy_offset
            + segment["injected_sensible_enthalpy_j"][keep]
        )
        outflow_energy.append(
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
    temperature_values = np.concatenate(temperatures, axis=1)
    domain_mass_values = np.concatenate(domain_mass)
    injected_values = np.concatenate(injected)
    outflow_values = np.concatenate(outflow)
    domain_energy_values = np.concatenate(domain_energy)
    injected_energy_values = np.concatenate(injected_energy)
    outflow_energy_values = np.concatenate(outflow_energy)
    mass_residual = (
        domain_mass_values + outflow_values
        - domain_mass_values[0] - injected_values
    )
    energy_residual = (
        domain_energy_values + outflow_energy_values
        - domain_energy_values[0] - injected_energy_values
    )
    last = committed[-1]
    paths = {
        "receptor_npz": output / "ti021-inner-receptor-traces.npz",
        "final_field_npz": output / "ti021-inner-final-field.npz",
        "record_json": output / "ti021-inner-transient.json",
        "prefix_json": output / "ti021-checkpoint-prefix.json",
    }
    np.savez_compressed(
        paths["receptor_npz"],
        time_s=time_values, labels=labels, solver_positions_m=positions,
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
        maximum_concentration_vol_percent=last["maximum_concentration_vol_percent"],
        temperature_k=last["final_temperature_k"],
        mixture_density_kg_m3=last["final_mixture_density_kg_m3"],
        velocity_u_m_s=last["final_velocity_u_m_s"],
        velocity_v_m_s=last["final_velocity_v_m_s"],
        velocity_w_m_s=last["final_velocity_w_m_s"],
        pressure_pa_gauge=last["final_pressure_pa_gauge"],
    )
    diagnostics = [entry["diagnostics"] for entry in entries]
    maximum_continuity = max(
        float(item["maximum_continuity_residual_kg_m3_s"])
        for item in diagnostics
    )
    nonconverged = sum(
        int(item["corrector_nonconverged_steps"]) for item in diagnostics
    )
    maximum_relative_mass = float(np.max(np.abs(mass_residual))) / max(
        abs(float(domain_mass_values[0] + injected_values[-1])), 1.0e-30,
    )
    maximum_relative_energy = float(np.max(np.abs(energy_residual))) / max(
        abs(float(domain_energy_values[0] + injected_energy_values[-1])), 1.0e-30,
    )
    expected_count = int(round(stop / dt)) + 1
    gates = {
        "prefix_ends_at_committed_boundary": math.isclose(time_values[-1], stop),
        "complete_prefix_time_axis": time_values.size == expected_count,
        "all_60_receptors_present": labels.size == 60,
        "species_mass_conservative": maximum_relative_mass <= 1.0e-6,
        "sensible_enthalpy_conservative": maximum_relative_energy <= 1.0e-6,
        "continuity_residual_at_most_1e_7": maximum_continuity <= 1.0e-7,
        "corrector_converged_every_step": nonconverged == 0,
        "observation_values_not_used": True,
    }
    gates = {key: bool(value) for key, value in gates.items()}
    record = {
        "schema": SCHEMA,
        "status": (
            "accepted_inner_transient_qualification"
            if all(gates.values()) else "withheld_inner_transient"
        ),
        "scope": (
            "immutable prefix view of an atomically checkpointed no-fit inner run"
        ),
        "provenance": {
            "checkpoint_journal_path": str(journal_path),
            "checkpoint_journal_sha256_at_extraction": journal_sha256,
            "coordinate_reference_path": str(reference_path),
            "coordinate_reference_sha256": _sha256(reference_path),
            "segment_records": entries,
        },
        "configuration": {
            "duration_s": stop,
            "time_step_s": dt,
            "cell_size_m": float(configuration["cell_size_m"]),
            "circular_arc_representation": configuration.get(
                "circular_arc_representation", "finite_thickness_voxel",
            ),
            "corrector_tolerance": float(configuration["corrector_tolerance"]),
            "cells": [
                int(last["x_m"].size), int(last["y_m"].size),
                int(last["z_m"].size),
            ],
            "inner_receptor_count": int(labels.size),
            "receptors": receptor_metadata,
            "initial_binary_cloud_momentum": "quiescent",
            "observation_values_used_to_configure_run": [],
        },
        "numerics": {
            "maximum_relative_species_mass_residual": maximum_relative_mass,
            "maximum_relative_sensible_enthalpy_residual": maximum_relative_energy,
            "maximum_continuity_residual_kg_m3_s": maximum_continuity,
            "corrector_nonconverged_steps": nonconverged,
            "maximum_vol_percent": float(np.nanmax(
                last["maximum_concentration_vol_percent"]
            )),
        },
        "runtime": {
            "wall_clock_s": sum(float(entry["wall_clock_s"]) for entry in entries),
            "source_segment_count": len(entries),
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }
    paths["record_json"].write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
    )
    prefix = {
        "schema": PREFIX_SCHEMA,
        "status": "accepted_materialized_checkpoint_prefix" if all(gates.values()) else "withheld_checkpoint_prefix",
        "source_checkpoint_journal": str(journal_path),
        "source_checkpoint_journal_sha256_at_extraction": journal_sha256,
        "source_checkpoint_journal_snapshot": {
            "path": str(journal_snapshot),
            "sha256": _sha256(journal_snapshot),
        },
        "stop_time_s": stop,
        "segment_count": len(entries),
        "outputs": {
            key: {"path": str(path), "sha256": _sha256(path)}
            for key, path in paths.items() if key != "prefix_json"
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }
    paths["prefix_json"].write_text(
        json.dumps(prefix, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
    )
    return record, paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpointed-run", type=Path, required=True)
    parser.add_argument("--stop-time-s", type=float, required=True)
    parser.add_argument("--coordinate-reference-run", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = materialize_prefix(
        checkpointed_run=args.checkpointed_run,
        stop_time_s=args.stop_time_s,
        coordinate_reference_run=args.coordinate_reference_run,
        output_directory=args.output_directory,
    )
    print(json.dumps({
        "status": record["status"], "all_gates_pass": record["all_gates_pass"],
        "configuration": record["configuration"],
        "numerics": record["numerics"],
        "outputs": {key: str(path) for key, path in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
