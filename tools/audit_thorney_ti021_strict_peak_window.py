"""Run the mandatory strict corrector postcheck near the model plume peak."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.binary_mixture_handoff_plane_3d import (
    BinaryMixtureHandoffPlane3DRequest,
)
from degali.addons.resolved_obstacle_flow_3d import (
    read_resolved_obstacle_flow_3d_restart,
    rebind_restart_for_corrector_tolerance_audit,
    solve_resolved_obstacle_flow_3d,
)
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

try:
    from tools.audit_thorney_ti021_corrector_tolerance import _difference
    from tools.audit_thorney_ti021_handoff_streaming_preflight import build_configuration
    from tools.run_thorney_ti021_inner_transient import _sha256, inner_receptors
except ModuleNotFoundError:
    from audit_thorney_ti021_corrector_tolerance import _difference
    from audit_thorney_ti021_handoff_streaming_preflight import build_configuration
    from run_thorney_ti021_inner_transient import _sha256, inner_receptors


SCHEMA = "degali.thorney-ti021-strict-peak-window-audit.v1"
FIELD_NAMES = (
    "concentration_vol_percent", "velocity_u_m_s", "velocity_v_m_s",
    "velocity_w_m_s", "pressure_pa_gauge", "temperature_k",
)
HANDOFF_NAMES = (
    "total_mass_flux_anomaly_kg_m2_s", "species_mass_flux_kg_m2_s",
    "sensible_enthalpy_flux_w_m2", "lateral_momentum_flux_anomaly_pa",
    "vertical_momentum_flux_anomaly_pa",
)


def select_model_peak_segment(
    time_s: np.ndarray, concentration_vol_percent: np.ndarray,
    segments: list[dict[str, object]],
) -> dict[str, float | int]:
    """Select the committed segment containing the model-only receptor peak."""

    time_values = np.asarray(time_s, dtype=float)
    traces = np.asarray(concentration_vol_percent, dtype=float)
    if time_values.ndim != 1 or traces.ndim != 2 or traces.shape[1] != time_values.size:
        raise ValueError("receptor traces must be receptor-by-time")
    if np.any(np.diff(time_values) <= 0.0) or np.any(~np.isfinite(traces)):
        raise ValueError("peak selection requires finite traces on an increasing clock")
    flattened = int(np.argmax(traces))
    receptor_index, time_index = np.unravel_index(flattened, traces.shape)
    peak_time = float(time_values[time_index])
    for index, entry in enumerate(segments):
        start = float(entry["start_time_s"])
        stop = float(entry["stop_time_s"])
        is_last = index == len(segments) - 1
        if start <= peak_time < stop or (is_last and math.isclose(peak_time, stop)):
            return {
                "segment_index": int(entry["segment_index"]),
                "start_time_s": start,
                "stop_time_s": stop,
                "model_peak_time_s": peak_time,
                "model_peak_receptor_index": int(receptor_index),
                "model_peak_vol_percent": float(traces[receptor_index, time_index]),
            }
    raise ValueError("model peak does not lie in any committed segment")


def _result_fields(result) -> dict[str, np.ndarray]:
    return {
        "concentration_vol_percent": np.asarray(result.concentration_vol_percent),
        "velocity_u_m_s": np.asarray(result.velocity_u_m_s),
        "velocity_v_m_s": np.asarray(result.velocity_v_m_s),
        "velocity_w_m_s": np.asarray(result.velocity_w_m_s),
        "pressure_pa_gauge": np.asarray(result.pressure_pa_gauge),
        "temperature_k": np.asarray(result.temperature_k),
    }


def _trace_array(result, labels: list[str]) -> np.ndarray:
    return np.asarray([
        result.receptor_traces_vol_percent[label] for label in labels
    ])


def _frame_array(frames: list[object], name: str) -> np.ndarray:
    return np.asarray([np.asarray(getattr(frame, name)) for frame in frames])


def _atomic_npz(path: Path, **arrays: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("xb") as stream:
            np.savez_compressed(stream, **arrays)
        temporary.replace(path)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def run_strict_peak_window_audit(
    *, candidate_directory: Path, output_directory: Path,
) -> tuple[dict[str, object], dict[str, Path]]:
    candidate = candidate_directory.resolve()
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    record_path = candidate / "ti021-inner-transient.json"
    journal_path = candidate / "ti021-inner-checkpoint.json"
    trace_path = candidate / "ti021-inner-receptor-traces.npz"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    journal = json.loads(journal_path.read_text(encoding="utf-8"))
    if not record.get("all_gates_pass") or record["configuration"]["duration_s"] != 600.0:
        raise ValueError("candidate must be a completed accepted 600 s run")
    if journal.get("status") != "finalized":
        raise ValueError("candidate checkpoint journal must be finalized")
    if record["configuration"]["corrector_tolerance"] != 1.0e-4:
        raise ValueError("candidate production tolerance must be 1e-4")
    with np.load(trace_path, allow_pickle=False) as source:
        full_trace = {key: np.asarray(source[key]).copy() for key in source.files}
    selected = select_model_peak_segment(
        full_trace["time_s"], full_trace["concentration_vol_percent"],
        list(journal["segments"]),
    )
    segment_index = int(selected["segment_index"])
    selected_entry = journal["segments"][segment_index - 1]
    selected_segment_path = candidate / str(selected_entry["segment_file"])
    if _sha256(selected_segment_path) != selected_entry["segment_sha256"]:
        raise ValueError("selected candidate segment SHA-256 mismatch")

    workbook_path = Path(record["provenance"]["workbook_path"]).resolve()
    property_path = Path(record["provenance"]["property_audit_path"]).resolve()
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    configuration = record["configuration"]
    duration = float(configuration["duration_s"])
    dt = float(configuration["time_step_s"])
    dx = float(configuration["cell_size_m"])
    active_steps = int(round(duration / dt))
    candidate_config, _closure, _coordinates = build_configuration(
        properties, cell_size_m=dx, time_step_s=dt,
        active_steps=active_steps, corrector_tolerance=1.0e-4,
    )
    strict_config, _closure, _coordinates = build_configuration(
        properties, cell_size_m=dx, time_step_s=dt,
        active_steps=active_steps, corrector_tolerance=1.0e-6,
    )
    start_time = float(selected["start_time_s"])
    stop_time = float(selected["stop_time_s"])
    source_restart_path = None
    candidate_restart = None
    strict_restart = None
    if start_time > 0.0:
        previous_entry = journal["segments"][segment_index - 2]
        source_restart_path = candidate / str(previous_entry["restart_file"])
        if _sha256(source_restart_path) != previous_entry["restart_sha256"]:
            raise ValueError("source restart SHA-256 mismatch")
        candidate_restart = read_resolved_obstacle_flow_3d_restart(source_restart_path)
        if not math.isclose(candidate_restart.time_s, start_time):
            raise ValueError("source restart does not match selected window start")
        strict_restart = rebind_restart_for_corrector_tolerance_audit(
            candidate_restart,
            source_config=candidate_config,
            target_config=strict_config,
        )

    dataset = read_thorney_ti021_workbook(workbook_path)
    receptors, _metadata = inner_receptors(dataset)
    labels = [item.label for item in receptors]
    request = BinaryMixtureHandoffPlane3DRequest(
        x_face_m=160.0, sample_interval_s=dt, include_initial=True,
    )
    candidate_frames: list[object] = []
    started = time.perf_counter()
    candidate_result = solve_resolved_obstacle_flow_3d(
        candidate_config, receptors=receptors, handoff_request=request,
        handoff_observer=candidate_frames.append,
        restart_state=candidate_restart, stop_time_s=stop_time,
    )
    candidate_wall = time.perf_counter() - started
    strict_frames: list[object] = []
    started = time.perf_counter()
    strict_result = solve_resolved_obstacle_flow_3d(
        strict_config, receptors=receptors, handoff_request=request,
        handoff_observer=strict_frames.append,
        restart_state=strict_restart, stop_time_s=stop_time,
    )
    strict_wall = time.perf_counter() - started
    candidate_fields = _result_fields(candidate_result)
    strict_fields = _result_fields(strict_result)
    field_differences = {
        name: _difference(strict_fields[name], candidate_fields[name])
        for name in FIELD_NAMES
    }
    candidate_traces = _trace_array(candidate_result, labels)
    strict_traces = _trace_array(strict_result, labels)
    receptor_difference = _difference(strict_traces, candidate_traces)
    if not np.array_equal(candidate_result.time_s, strict_result.time_s):
        raise ValueError("strict and candidate clocks differ")
    if len(candidate_frames) != len(strict_frames):
        raise ValueError("strict and candidate handoff frame counts differ")
    handoff_differences = {
        name: _difference(
            _frame_array(strict_frames, name), _frame_array(candidate_frames, name),
        )
        for name in HANDOFF_NAMES
    }
    with np.load(selected_segment_path, allow_pickle=False) as committed:
        committed_final_difference = _difference(
            committed["final_concentration_vol_percent"],
            candidate_fields["concentration_vol_percent"],
        )
        committed_trace_difference = _difference(
            committed["concentration_vol_percent"], candidate_traces,
        )

    candidate_diag = candidate_result.diagnostics
    strict_diag = strict_result.diagnostics
    velocity_linf = max(
        field_differences[name]["maximum_absolute_difference"]
        for name in ("velocity_u_m_s", "velocity_v_m_s", "velocity_w_m_s")
    )
    gates = {
        "window_contains_model_receptor_peak": (
            start_time <= float(selected["model_peak_time_s"]) <= stop_time
        ),
        "window_selected_without_observation_values": True,
        "candidate_replay_final_linf_at_most_1e_10vol_percent": (
            committed_final_difference["maximum_absolute_difference"] <= 1.0e-10
        ),
        "candidate_replay_trace_linf_at_most_1e_10vol_percent": (
            committed_trace_difference["maximum_absolute_difference"] <= 1.0e-10
        ),
        "candidate_concentration_relative_l1_at_most_5e_4": (
            field_differences["concentration_vol_percent"]["relative_l1_difference"] <= 5.0e-4
        ),
        "candidate_concentration_linf_at_most_0p05vol_percent": (
            field_differences["concentration_vol_percent"]["maximum_absolute_difference"] <= 0.05
        ),
        "candidate_receptor_linf_at_most_0p01vol_percent": (
            receptor_difference["maximum_absolute_difference"] <= 0.01
        ),
        "candidate_velocity_component_linf_at_most_0p002m_s": velocity_linf <= 0.002,
        "candidate_handoff_mass_anomaly_relative_l1_at_most_5e_5": (
            handoff_differences["total_mass_flux_anomaly_kg_m2_s"]["relative_l1_difference"] <= 5.0e-5
        ),
        "both_continuity_residuals_at_most_1e_7": max(
            candidate_diag["maximum_variable_density_continuity_residual_kg_m3_s"],
            strict_diag["maximum_variable_density_continuity_residual_kg_m3_s"],
        ) <= 1.0e-7,
        "both_species_mass_residuals_at_most_1e_8": max(
            candidate_diag["maximum_relative_mass_residual"],
            strict_diag["maximum_relative_mass_residual"],
        ) <= 1.0e-8,
        "both_energy_residuals_at_most_1e_8": max(
            candidate_diag["maximum_relative_energy_residual"],
            strict_diag["maximum_relative_energy_residual"],
        ) <= 1.0e-8,
        "both_correctors_converged_every_step": (
            candidate_diag["thermodynamic_corrector_nonconverged_steps"] == 0
            and strict_diag["thermodynamic_corrector_nonconverged_steps"] == 0
        ),
    }
    gates = {key: bool(value) for key, value in gates.items()}
    accepted = all(gates.values())
    candidate_npz = output / "candidate-1e4-peak-window.npz"
    strict_npz = output / "reference-1e6-peak-window.npz"
    _atomic_npz(
        candidate_npz, time_s=candidate_result.time_s, labels=np.asarray(labels),
        receptor_concentration_vol_percent=candidate_traces,
        **candidate_fields,
    )
    _atomic_npz(
        strict_npz, time_s=strict_result.time_s, labels=np.asarray(labels),
        receptor_concentration_vol_percent=strict_traces,
        **strict_fields,
    )
    audit = {
        "schema": SCHEMA,
        "status": (
            "accepted_1e4_corrector_tolerance_at_model_peak"
            if accepted else "withheld_1e4_corrector_tolerance_at_model_peak"
        ),
        "scope": (
            "matched candidate/strict replay of the checkpoint segment containing "
            "the model-only maximum receptor concentration"
        ),
        "selection": {
            **selected,
            "selection_signal": "model receptor concentration only",
            "observation_value_fields_accessed": [],
        },
        "provenance": {
            "candidate_record": {"path": str(record_path), "sha256": _sha256(record_path)},
            "checkpoint_journal": {"path": str(journal_path), "sha256": _sha256(journal_path)},
            "selected_candidate_segment": {"path": str(selected_segment_path), "sha256": _sha256(selected_segment_path)},
            "source_restart": None if source_restart_path is None else {"path": str(source_restart_path), "sha256": _sha256(source_restart_path)},
            "candidate_replay_npz": {"path": str(candidate_npz), "sha256": _sha256(candidate_npz)},
            "strict_reference_npz": {"path": str(strict_npz), "sha256": _sha256(strict_npz)},
        },
        "configuration": {
            "cell_size_m": dx, "time_step_s": dt,
            "candidate_tolerance": 1.0e-4, "strict_tolerance": 1.0e-6,
            "window_start_s": start_time, "window_stop_s": stop_time,
            "identical_conserved_initial_state": True,
            "restart_rebind_scope": (
                "only thermodynamic corrector density and velocity tolerances"
            ),
        },
        "runtime": {
            "candidate_wall_clock_s": candidate_wall,
            "strict_wall_clock_s": strict_wall,
            "strict_to_candidate_ratio": strict_wall / candidate_wall,
        },
        "candidate_replay_differences": {
            "final_concentration": committed_final_difference,
            "receptor_trace": committed_trace_difference,
        },
        "strict_candidate_field_differences": field_differences,
        "strict_candidate_receptor_difference": receptor_difference,
        "strict_candidate_handoff_differences": handoff_differences,
        "conservation": {
            "candidate": {
                "maximum_relative_species_mass_residual": candidate_diag["maximum_relative_mass_residual"],
                "maximum_relative_sensible_enthalpy_residual": candidate_diag["maximum_relative_energy_residual"],
                "maximum_continuity_residual_kg_m3_s": candidate_diag["maximum_variable_density_continuity_residual_kg_m3_s"],
            },
            "strict": {
                "maximum_relative_species_mass_residual": strict_diag["maximum_relative_mass_residual"],
                "maximum_relative_sensible_enthalpy_residual": strict_diag["maximum_relative_energy_residual"],
                "maximum_continuity_residual_kg_m3_s": strict_diag["maximum_variable_density_continuity_residual_kg_m3_s"],
            },
        },
        "gates": gates,
        "all_gates_pass": accepted,
        "production_tolerance_promoted": accepted,
    }
    json_path = output / "ti021-strict-peak-window-audit.json"
    plot_path = output / "ti021-strict-peak-window-audit.png"
    json_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    difference = candidate_fields["concentration_vol_percent"] - strict_fields["concentration_vol_percent"]
    z_index = int(np.argmax(np.max(np.abs(difference), axis=(1, 2))))
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    image = axes[0].imshow(difference[z_index], origin="lower", cmap="coolwarm")
    axes[0].set(title=f"Candidate − strict concentration, z-index {z_index}", xlabel="x cell", ylabel="y cell")
    figure.colorbar(image, ax=axes[0], label="vol%")
    axes[1].scatter(strict_traces.ravel(), candidate_traces.ravel(), s=8, alpha=0.5)
    limit = max(float(np.max(strict_traces)), float(np.max(candidate_traces)), 0.01)
    axes[1].plot((0, limit), (0, limit), "k--", linewidth=1)
    axes[1].set(xlabel="Strict 1e-6 receptor concentration (vol%)", ylabel="Candidate 1e-4 (vol%)", title="Peak-window receptor replay")
    axes[1].grid(alpha=0.25)
    figure.suptitle("TI021 strict nonlinear-corrector peak-window postcheck")
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    return audit, {"json": json_path, "plot": plot_path, "candidate_npz": candidate_npz, "strict_npz": strict_npz}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = run_strict_peak_window_audit(
        candidate_directory=args.candidate_directory,
        output_directory=args.output_directory,
    )
    print(json.dumps({
        "status": record["status"], "selection": record["selection"],
        "runtime": record["runtime"], "gates": record["gates"],
        "outputs": {key: str(value) for key, value in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
