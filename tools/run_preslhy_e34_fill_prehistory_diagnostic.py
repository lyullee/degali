"""Run and score the post-primary PRESLHY fill-prehistory diagnostic.

This is explicitly an observation-informed diagnostic sequence, not a new
independent validation. The fill-onset rule itself is mass-only and the source
rate is not fitted, but the need for prehistory was diagnosed from the failed
primary. The failed primary remains the authoritative frozen score.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import runpy
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import (  # noqa: E402
    ResolvedInletVelocityHistory3D,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule  # noqa: E402


PRIMARY_TOOL = runpy.run_path(
    ROOT / "tools/preflight_preslhy_e34_dynamic_component.py"
)
SCORE_TOOL = runpy.run_path(
    ROOT / "tools/score_preslhy_e34_dynamic_component.py"
)
SCHEMA = "degali.preslhy-e34-fill-prehistory-diagnostic-run.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def _schedule(mapping: dict[str, Any]) -> SourceRateSchedule:
    return SourceRateSchedule(
        time_s=tuple(float(value) for value in mapping["time_s"]),
        rate_kg_s=tuple(float(value) for value in mapping["rate_kg_s"]),
        source_id=str(mapping["source_id"]),
        rate_operator=str(mapping["rate_operator"]),
    )


def _history(mapping: dict[str, Any]) -> ResolvedInletVelocityHistory3D:
    return ResolvedInletVelocityHistory3D(
        time_s=tuple(float(value) for value in mapping["time_s"]),
        u_m_s=tuple(float(value) for value in mapping["u_m_s"]),
        v_m_s=tuple(float(value) for value in mapping["v_m_s"]),
        w_m_s=tuple(float(value) for value in mapping["w_m_s"]),
        history_id=str(mapping["history_id"]),
        interpolation=str(mapping["interpolation"]),
    )


def _numerical_gates(diagnostics: dict[str, Any]) -> dict[str, bool]:
    return {
        "relative_hydrogen_mass_residual": bool(
            diagnostics["maximum_relative_mass_residual"] <= 1.0e-6
        ),
        "relative_sensible_energy_residual": bool(
            diagnostics["maximum_relative_energy_residual"] <= 1.0e-6
        ),
        "variable_density_continuity": bool(
            diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ] <= 1.0e-7
        ),
        "lattice_mach": bool(
            diagnostics["maximum_lattice_mach_number"]
            <= diagnostics["maximum_allowed_lattice_mach_number"]
        ),
        "active_scalar_applicability": bool(
            diagnostics["active_scalar_applicability_passed"]
        ),
        "thermodynamic_corrector_convergence": bool(
            diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
        ),
    }


def run_fill_prehistory_diagnostic(
    diagnostic_input_path: str | Path,
) -> dict[str, Any]:
    input_path = Path(diagnostic_input_path).resolve()
    directory = input_path.parent
    config_path = directory / "diagnostic-model-configuration.json"
    model_path = directory / "diagnostic-model-output.npz"
    result_path = directory / "diagnostic-run-score.json"
    plot_path = directory / "diagnostic-run-score.png"
    if any(path.exists() for path in (config_path, model_path, result_path, plot_path)):
        raise FileExistsError("one or more diagnostic run outputs already exist")
    package = _load_json(input_path)
    if package.get("schema") != (
        "degali.preslhy-e34-fill-prehistory-diagnostic-input.v1"
    ):
        raise ValueError("unsupported fill-prehistory diagnostic input schema")
    if not package["evidence_decision"]["diagnostic_run_allowed"]:
        raise ValueError("fill-prehistory diagnostic input withholds execution")
    primary_score_path = Path(package["failed_primary_score"]["path"]).resolve()
    if _sha256(primary_score_path) != package["failed_primary_score"]["sha256"]:
        raise ValueError("failed primary score fingerprint changed")
    primary_score = _load_json(primary_score_path)
    primary_freeze_path = Path(
        primary_score["configuration_freeze"]["path"]
    ).resolve()
    primary_freeze, freeze_record, primary_config, receptors = PRIMARY_TOOL[
        "_load_frozen_configuration"
    ](primary_freeze_path.parent)
    if primary_freeze != primary_freeze_path:
        raise ValueError("loaded primary freeze path differs from score")
    duration = float(package["mass_only_prehistory"]["full_duration_s"])
    diagnostic_config = replace(
        primary_config,
        duration_s=duration,
        source_schedule=_schedule(package["source_schedule"]),
        active_inlet_velocity_history=_history(package["inlet_velocity_history"]),
    )
    serialized = PRIMARY_TOOL["_serializable_config"](
        diagnostic_config, receptors,
    )
    config_record = {
        "schema": "degali.preslhy-e34-fill-prehistory-diagnostic-configuration.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "diagnostic_input": {"path": str(input_path), "sha256": _sha256(input_path)},
        "failed_primary_configuration": {
            "path": str(primary_freeze_path),
            "sha256": _sha256(primary_freeze_path),
            "model_configuration_sha256": freeze_record[
                "model_configuration_sha256"
            ],
        },
        "configuration_changes_from_failed_primary": [
            "prepend_mass_only_fill_prehistory",
            "extend_frozen_postfill_source_rate_backward",
            "prepend_measured_wind_only_history_in_the_same_primary_axes",
        ],
        "unchanged_model_choices": [
            "domain_grid_source_geometry_and_receptors",
            "transport_thermodynamics_and_momentum_closures",
            "time_step_boundary_scheme_and_numerical_gates",
        ],
        "observation_control": {
            "diagnostic_developed_after_failed_primary_observations_were_scored": True,
            "observation_used_to_select_fill_onset_fraction": False,
            "observation_used_to_select_source_rate": False,
            "eligible_as_independent_validation": False,
        },
        "model_configuration": serialized,
        "model_configuration_sha256": PRIMARY_TOOL["_canonical_sha256"](serialized),
    }
    config_path.write_text(
        json.dumps(config_record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    result = solve_resolved_obstacle_flow_3d(
        diagnostic_config, receptors=receptors,
    )
    diagnostics = result.diagnostics
    gates = _numerical_gates(diagnostics)
    arrays: dict[str, np.ndarray] = {
        "x_m": result.x_m,
        "y_m": result.y_m,
        "z_m": result.z_m,
        "time_s": result.time_s,
        "final_h2_vol_percent": result.concentration_vol_percent,
        "maximum_h2_vol_percent": result.maximum_concentration_vol_percent,
        "injected_mass_kg": result.injected_mass_kg,
        "domain_mass_kg": result.domain_mass_kg,
        "outflow_mass_kg": result.outflow_mass_kg,
        "mass_residual_kg": result.mass_residual_kg,
    }
    for label, values in result.receptor_traces_vol_percent.items():
        arrays[f"model_receptor_{label}_h2_vol_percent"] = values
    np.savez_compressed(model_path, **arrays)

    primary_arrays_path = Path(primary_score["score_arrays"]["path"]).resolve()
    if _sha256(primary_arrays_path) != primary_score["score_arrays"]["sha256"]:
        raise ValueError("primary score arrays fingerprint changed")
    offset = float(package["model_clock"]["score_start_model_time_s"])
    scoring_policy_path = Path(
        primary_score["scoring_policy_freeze"]["path"]
    ).resolve()
    scoring_policy = _load_json(scoring_policy_path)["policy"]
    channel_scores: dict[str, dict[str, Any]] = {}
    aligned_model: dict[str, np.ndarray] = {}
    observed_arrays: dict[str, np.ndarray] = {}
    primary_model_arrays: dict[str, np.ndarray] = {}
    with np.load(primary_arrays_path) as primary_arrays:
        observation_time = np.asarray(primary_arrays["time_s"], dtype=float)
        for label in package["sealed_observation_operator_metadata"]["usable_channels"]:
            observed = np.asarray(
                primary_arrays[f"observed_{label}_h2_vol_percent"], dtype=float,
            )
            model = np.interp(
                observation_time + offset,
                result.time_s,
                result.receptor_traces_vol_percent[label],
            )
            observed_arrays[label] = observed
            primary_model_arrays[label] = np.asarray(
                primary_arrays[f"model_{label}_h2_vol_percent"], dtype=float,
            )
            aligned_model[label] = model
            channel_scores[label] = SCORE_TOOL["_channel_score"](
                observation_time,
                observed,
                model,
                log_floor=float(scoring_policy["log_metric_floor_vol_percent"]),
                arrival_threshold=float(scoring_policy["arrival_threshold_vol_percent"]),
                flammability_threshold=float(
                    scoring_policy["flammability_threshold_vol_percent"]
                ),
            )
    eligible = [
        label for label, metric in channel_scores.items()
        if metric["observed_peak_vol_percent"] >= 1.0
    ]
    median_peak_error = float(np.median([
        channel_scores[label]["absolute_peak_error_vol_percent_point"]
        for label in eligible
    ]))
    median_dose_error = float(np.median([
        channel_scores[label]["absolute_relative_dose_error"] for label in eligible
    ]))
    arrival_values = [
        channel_scores[label]["normalized_arrival_time_error"] for label in eligible
    ]
    median_arrival_error = (
        None if any(value is None for value in arrival_values)
        else float(np.median(arrival_values))
    )
    pooled_observed = np.concatenate([observed_arrays[label] for label in eligible])
    pooled_model = np.concatenate([aligned_model[label] for label in eligible])
    pooled_four = SCORE_TOOL["_binary_counts"](
        pooled_observed, pooled_model, 4.0,
    )
    primary_aggregate = primary_score["aggregate_targets"]
    result_record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "diagnostic_input": {"path": str(input_path), "sha256": _sha256(input_path)},
        "diagnostic_configuration": {
            "path": str(config_path), "sha256": _sha256(config_path),
        },
        "diagnostic_model_output": {
            "path": str(model_path), "sha256": _sha256(model_path),
        },
        "numerical_diagnostics": {
            "duration_s": duration,
            "steps": int(round(duration / diagnostic_config.scalar_time_step_s)),
            "runtime_s": diagnostics["total_runtime_s"],
            "maximum_h2_vol_percent": float(
                np.nanmax(result.maximum_concentration_vol_percent)
            ),
            "maximum_relative_hydrogen_mass_residual": diagnostics[
                "maximum_relative_mass_residual"
            ],
            "maximum_relative_sensible_energy_residual": diagnostics[
                "maximum_relative_energy_residual"
            ],
            "maximum_relative_mixture_density_departure": diagnostics[
                "maximum_relative_mixture_density_departure"
            ],
            "maximum_lattice_mach_number": diagnostics[
                "maximum_lattice_mach_number"
            ],
            "maximum_variable_density_continuity_residual_kg_m3_s": diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
        },
        "numerical_gate_results": gates,
        "all_numerical_gates_passed": bool(all(gates.values())),
        "diagnostic_channel_scores": channel_scores,
        "diagnostic_aggregate": {
            "median_absolute_peak_error_vol_percent_point": median_peak_error,
            "median_absolute_relative_dose_error": median_dose_error,
            "median_normalized_arrival_time_error": median_arrival_error,
            "pooled_four_vol_sample_classification": pooled_four,
        },
        "change_from_failed_primary": {
            "median_peak_error_vol_percent_point": (
                median_peak_error
                - primary_aggregate[
                    "median_absolute_peak_error_vol_percent_point"
                ]
            ),
            "median_relative_dose_error": (
                median_dose_error
                - primary_aggregate["median_absolute_relative_dose_error"]
            ),
            "median_normalized_arrival_time_error": (
                None
                if median_arrival_error is None
                else median_arrival_error
                - primary_aggregate["median_normalized_arrival_time_error"]
            ),
            "pooled_four_vol_sensitivity": (
                None if pooled_four["sensitivity"] is None else
                pooled_four["sensitivity"]
                - primary_aggregate[
                    "pooled_four_vol_sample_classification"
                ]["sensitivity"]
            ),
            "pooled_four_vol_specificity": (
                None if pooled_four["specificity"] is None else
                pooled_four["specificity"]
                - primary_aggregate[
                    "pooled_four_vol_sample_classification"
                ]["specificity"]
            ),
        },
        "evidence_decision": {
            "status": "diagnostic_only_after_failed_primary",
            "failed_primary_replaced": False,
            "eligible_as_independent_validation": False,
            "obstacle_transport_model_qualification_allowed": False,
            "quantitative_facility_prediction_allowed": False,
        },
        "review_plot": str(plot_path),
    }
    figure, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for axis, label in zip(axes, eligible):
        axis.plot(observation_time, observed_arrays[label], label=f"{label} observed")
        axis.plot(
            observation_time, aligned_model[label],
            label=f"{label} fill-prehistory diagnostic",
        )
        axis.plot(
            observation_time,
            primary_model_arrays[label],
            label=f"{label} failed primary", alpha=0.65,
        )
        axis.axhline(4.0, color="tab:red", ls="--", lw=0.9)
        axis.set_ylabel("H2 (vol%)")
        axis.grid(alpha=0.25)
        axis.legend(loc="upper right")
    axes[-1].set_xlabel("Primary score-window time (s)")
    figure.suptitle(
        "PRESLHY E3.4 fill-prehistory diagnostic — failed primary retained"
    )
    figure.tight_layout()
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    result_path.write_text(
        json.dumps(result_record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return result_record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diagnostic-input", type=Path, required=True)
    arguments = parser.parse_args(argv)
    record = run_fill_prehistory_diagnostic(arguments.diagnostic_input)
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if record["all_numerical_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
