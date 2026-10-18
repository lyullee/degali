"""Freeze and preflight a concentration-blind PRESLHY E3.4 replay.

This command consumes the hash-bound replay input produced by
``build_preslhy_e34_dynamic_replay_package.py``.  It validates that no H2
values were exported, freezes every numerical and physical modelling choice,
writes that freeze to disk, and only then runs a short numerical preflight.

The E3.4 record has no machine-registered obstacle geometry.  Consequently
this artifact is dynamic inlet/source/sensor-height component evidence only;
it must never be promoted to obstacle-transport validation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.fast_obstacle_transport import (  # noqa: E402
    HYDROGEN_MOLAR_MASS_KG_MOL,
    R_GAS_J_MOL_K,
)
from degali.addons.fast_obstacle_transport_3d import (  # noqa: E402
    FastObstacleReceptor3D,
)
from degali.addons.resolved_obstacle_flow_3d import (  # noqa: E402
    ResolvedInletVelocityHistory3D,
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule  # noqa: E402


SCHEMA = "degali.preslhy-e34-dynamic-component-configuration-freeze.v1"
RESULT_SCHEMA = "degali.preslhy-e34-dynamic-component-preflight.v1"
SUPPORTED_REPLAY_SCHEMA = "degali.preslhy-e34-dynamic-replay-input.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _load_replay_input(path: Path) -> dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8")
    record = json.loads(raw_text)
    if not isinstance(record, dict):
        raise ValueError("dynamic replay input must be a JSON object")
    if record.get("schema") != SUPPORTED_REPLAY_SCHEMA:
        raise ValueError("unsupported dynamic replay input schema")
    decision = record.get("evidence_decision", {})
    if not decision.get("dynamic_component_replay_input_allowed", False):
        raise ValueError("dynamic component replay is not allowed by evidence gate")
    if decision.get("obstacle_transport_model_qualification_allowed") is not False:
        raise ValueError("E3.4 replay must not claim obstacle qualification")
    sealed = record.get("sealed_observation_operator", {})
    if sealed.get("values_exported") is not False:
        raise ValueError("H2 observation values must remain unexported before freeze")
    if sealed.get("values_summarized") is not False:
        raise ValueError("H2 observation values must remain unsummarized before freeze")
    if sealed.get("sealed_until_model_configuration_is_frozen") is not True:
        raise ValueError("observation operator is not declared sealed until freeze")
    # This rejects accidental future additions of a value-bearing channel even
    # if the booleans above are left unchanged. Channel names and positions are
    # intentionally permitted metadata.
    prohibited = (
        "hydrogen_vol_percent", "h2_vol_percent", "observed_values",
        "concentration_values",
    )
    folded = raw_text.casefold()
    if any(token in folded for token in prohibited):
        raise ValueError("dynamic replay input contains a prohibited H2 value field")

    for dependency_name in ("source_freeze", "dynamic_evidence_audit"):
        dependency = record.get(dependency_name, {})
        dependency_path = Path(dependency.get("path", "")).resolve()
        if not dependency_path.is_file():
            raise FileNotFoundError(dependency_path)
        if _sha256(dependency_path) != str(dependency.get("sha256", "")).lower():
            raise ValueError(f"{dependency_name} fingerprint changed")
    return record


def _history_from_mapping(mapping: dict[str, Any]) -> ResolvedInletVelocityHistory3D:
    return ResolvedInletVelocityHistory3D(
        time_s=tuple(float(value) for value in mapping["time_s"]),
        u_m_s=tuple(float(value) for value in mapping["u_m_s"]),
        v_m_s=tuple(float(value) for value in mapping["v_m_s"]),
        w_m_s=tuple(float(value) for value in mapping["w_m_s"]),
        history_id=str(mapping["history_id"]),
        interpolation=str(mapping["interpolation"]),
    )


def _schedule_from_mapping(mapping: dict[str, Any]) -> SourceRateSchedule:
    return SourceRateSchedule(
        time_s=tuple(float(value) for value in mapping["time_s"]),
        rate_kg_s=tuple(float(value) for value in mapping["rate_kg_s"]),
        source_id=str(mapping["source_id"]),
        rate_operator=str(mapping["rate_operator"]),
    )


def _model_objects(record: dict[str, Any]) -> tuple[
    ResolvedObstacleFlow3DConfig,
    tuple[FastObstacleReceptor3D, ...],
    dict[str, Any],
]:
    duration_s = float(record["duration_s"])
    source_schedule = _schedule_from_mapping(record["source_schedule"])
    inlet_history = _history_from_mapping(record["inlet_velocity_history"])
    if not math.isclose(source_schedule.duration_s, duration_s, abs_tol=1.0e-9):
        raise ValueError("source schedule and replay duration differ")
    if not math.isclose(inlet_history.time_s[-1], duration_s, abs_tol=1.0e-9):
        raise ValueError("inlet history and replay duration differ")

    source_temperature_k = 20.3
    source_area_m2 = 0.50 * 0.50
    source_rate_kg_s = float(source_schedule.rate_kg_s[0])
    ideal_hydrogen_density = (
        101_325.0 * HYDROGEN_MOLAR_MASS_KG_MOL
        / (R_GAS_J_MOL_K * source_temperature_k)
    )
    source_vertical_velocity = (
        source_rate_kg_s / (ideal_hydrogen_density * source_area_m2)
    )
    reference_wind_speed = float(np.mean(inlet_history.u_m_s[:-1]))
    dx_m = 0.10
    source_centre = (0.75, 0.80, 0.05)
    square_second_moment_sigma = 0.50 / math.sqrt(12.0)
    vertical_second_moment_sigma = dx_m / math.sqrt(12.0)
    receptors = tuple(
        FastObstacleReceptor3D(
            label,
            (source_centre[0], source_centre[1], float(position[2])),
        )
        for label, position in record["sealed_observation_operator"][
            "positions_pool_centre_m"
        ].items()
        if label in record["sealed_observation_operator"]["usable_channels"]
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(2.0, 1.6, 1.0),
        cells=(20, 16, 10),
        wind_speed_m_s=reference_wind_speed,
        kinematic_viscosity_m2_s=0.005,
        flow_steps=100,
        minimum_flow_steps=100,
        lattice_inlet_speed=0.03,
        duration_s=duration_s,
        scalar_time_step_s=0.01,
        scalar_diffusivity_m2_s=6.1e-5,
        thermal_diffusivity_m2_s=2.2e-5,
        source_schedule=source_schedule,
        source_position_m=source_centre,
        source_sigma_m=square_second_moment_sigma,
        source_sigma_xyz_m=(
            square_second_moment_sigma,
            square_second_moment_sigma,
            vertical_second_moment_sigma,
        ),
        source_momentum_velocity_m_s=(0.0, 0.0, source_vertical_velocity),
        obstacles=(),
        scalar_advection_scheme="muscl_van_leer",
        source_time_integration="prestep",
        subgrid_model="wale",
        wale_constant=0.50,
        flow_boundary_scheme="characteristic_multidimensional",
        nonstreamwise_wall_boundary="all_pressure_open",
        active_inlet_velocity_history=inlet_history,
        ambient_temperature_k=293.15,
        ambient_pressure_pa=101_325.0,
        ambient_air_density_kg_m3=1.204,
        thermal_buoyancy_coupling=True,
        source_temperature_k=source_temperature_k,
        effective_schmidt_number=0.70,
        effective_prandtl_number=0.85,
        scalar_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        thermal_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        maximum_active_scalar_density_departure=0.10,
        active_momentum_solver="finite_volume_projection",
        active_provisional_face_projection="variable_density_only",
        require_thermodynamic_corrector_convergence=True,
        store_scalar_fields=False,
    )
    assumptions = {
        "coordinate_frame": record["wind_reduction"]["coordinate_frame"],
        "site_registered_geometry_available": False,
        "obstacles": [],
        "source_geometry": (
            "0.50 m square pool mapped to a cell-integrated Gaussian with "
            "matching horizontal second moment; one-cell vertical slab"
        ),
        "source_phase_model": (
            "corrected liquid mass loss treated conditionally as atmospheric "
            "hydrogen mass flux; no flash or real-gas source calculation"
        ),
        "source_temperature_basis": (
            "predeclared 20.3 K approximate normal-boiling-point LH2 proxy; "
            "not a measured event-mean pool temperature"
        ),
        "source_momentum_basis": (
            "ideal-gas H2 density at 101325 Pa and 20.3 K through 0.25 m2"
        ),
        "whole_inlet_basis": record["wind_reduction"]["whole_inlet_assumption"],
        "seed_viscosity_limitation": (
            "0.005 m2/s is a numerical D3Q19 seed viscosity required by the "
            "relaxation-time gate on the 0.10 m grid; it is not a physical "
            "atmospheric viscosity or Reynolds-matched claim"
        ),
        "pressure_open_limitation": (
            "characteristic_multidimensional passes selected acoustic tests "
            "but not the unified nonlinear vortical boundary gate"
        ),
        "qualification_scope": "dynamic_component_only_not_obstacle_transport",
    }
    return config, receptors, assumptions


def _serializable_config(
    config: ResolvedObstacleFlow3DConfig,
    receptors: tuple[FastObstacleReceptor3D, ...],
) -> dict[str, Any]:
    value = asdict(config)
    value["receptors"] = [asdict(receptor) for receptor in receptors]
    return value


def freeze_model_configuration(
    replay_input_path: str | Path,
    output_directory: str | Path,
) -> tuple[dict[str, Any], ResolvedObstacleFlow3DConfig, tuple[FastObstacleReceptor3D, ...]]:
    """Freeze the full replay configuration without opening H2 observations."""

    replay_path = Path(replay_input_path).resolve()
    output = Path(output_directory).resolve()
    if output.exists():
        raise FileExistsError(output)
    replay = _load_replay_input(replay_path)
    config, receptors, assumptions = _model_objects(replay)
    serialized_config = _serializable_config(config, receptors)
    record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "replay_input": {
            "path": str(replay_path),
            "sha256": _sha256(replay_path),
            "trial_id": replay["trial_id"],
            "window_id": replay["window_id"],
        },
        "observation_control": {
            "hydrogen_observation_values_opened": False,
            "hydrogen_observation_values_used": False,
            "model_configuration_frozen_before_observation_values_opened": True,
            "usable_channel_names_only": replay["sealed_observation_operator"][
                "usable_channels"
            ],
        },
        "model_configuration": serialized_config,
        "model_configuration_sha256": _canonical_sha256(serialized_config),
        "assumptions_and_limitations": assumptions,
        "predeclared_numerical_gates": {
            "maximum_relative_hydrogen_mass_residual": 1.0e-6,
            "maximum_relative_sensible_energy_residual": 1.0e-6,
            "maximum_variable_density_continuity_residual_kg_m3_s": 1.0e-7,
            "maximum_lattice_mach_number": config.maximum_lattice_mach_number,
            "active_scalar_out_of_range_cell_steps": 0,
            "thermodynamic_corrector_nonconverged_steps": 0,
        },
        "evidence_decision": {
            "status": "conditional_configuration_frozen",
            "dynamic_component_preflight_allowed": True,
            "hydrogen_observation_scoring_allowed_after_preflight_passes": True,
            "obstacle_transport_model_qualification_allowed": False,
            "quantitative_facility_prediction_allowed": False,
        },
    }
    output.mkdir(parents=True)
    freeze_path = output / "model-configuration-freeze.json"
    freeze_path.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return record, config, receptors


def _truncate_for_preflight(
    config: ResolvedObstacleFlow3DConfig,
    duration_s: float,
) -> ResolvedObstacleFlow3DConfig:
    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise ValueError("preflight duration must be positive and finite")
    if duration_s > config.duration_s:
        raise ValueError("preflight duration cannot exceed frozen replay duration")
    step_ratio = duration_s / config.scalar_time_step_s
    steps = int(round(step_ratio))
    if steps < 1 or not math.isclose(step_ratio, steps, abs_tol=1.0e-10):
        raise ValueError("preflight duration must align with the frozen time step")
    rate = config.source_schedule.rate_kg_s[0]
    schedule = SourceRateSchedule(
        time_s=(0.0, duration_s),
        rate_kg_s=(rate, 0.0),
        source_id=config.source_schedule.source_id + ":preflight-truncation",
        rate_operator=config.source_schedule.rate_operator,
    )
    history = config.active_inlet_velocity_history
    assert history is not None
    initial_velocity = history.velocity_at(0.0)
    truncated_history = ResolvedInletVelocityHistory3D(
        time_s=(0.0, duration_s),
        u_m_s=(initial_velocity[0], initial_velocity[0]),
        v_m_s=(initial_velocity[1], initial_velocity[1]),
        w_m_s=(initial_velocity[2], initial_velocity[2]),
        history_id=history.history_id + ":preflight-truncation",
        interpolation=history.interpolation,
    )
    return replace(
        config,
        duration_s=duration_s,
        source_schedule=schedule,
        active_inlet_velocity_history=truncated_history,
    )


def _preflight_record(
    result: object,
    *,
    freeze_path: Path,
    config: ResolvedObstacleFlow3DConfig,
) -> dict[str, Any]:
    diagnostics = result.diagnostics
    gates = {
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
    maximum_field = result.maximum_concentration_vol_percent
    return {
        "schema": RESULT_SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "configuration_freeze": {
            "path": str(freeze_path),
            "sha256": _sha256(freeze_path),
        },
        "preflight_duration_s": config.duration_s,
        "preflight_steps": int(round(config.duration_s / config.scalar_time_step_s)),
        "observation_control": {
            "hydrogen_observation_values_opened": False,
            "hydrogen_observation_values_used": False,
        },
        "numerical_diagnostics": {
            "runtime_s": diagnostics["total_runtime_s"],
            "maximum_h2_vol_percent": float(np.nanmax(maximum_field)),
            "four_vol_percent_volume_m3": float(
                np.count_nonzero(maximum_field >= 4.0) * 0.1**3
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
            "maximum_allowed_active_scalar_density_departure": diagnostics[
                "maximum_active_scalar_density_departure"
            ],
            "active_scalar_out_of_range_cell_steps": diagnostics[
                "active_scalar_out_of_range_cell_steps"
            ],
            "maximum_lattice_mach_number": diagnostics[
                "maximum_lattice_mach_number"
            ],
            "maximum_allowed_lattice_mach_number": diagnostics[
                "maximum_allowed_lattice_mach_number"
            ],
            "maximum_variable_density_continuity_residual_kg_m3_s": diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
            "maximum_variable_density_global_residual_kg_s": diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "maximum_momentum_balance_residual_kg_m_s": diagnostics[
                "maximum_momentum_balance_residual_kg_m_s"
            ],
            "maximum_momentum_advective_courant_number": diagnostics[
                "maximum_momentum_advective_courant_number"
            ],
            "maximum_momentum_viscous_stability_number": diagnostics[
                "maximum_momentum_viscous_stability_number"
            ],
            "scalar_advective_number": diagnostics["advective_number"],
            "scalar_diffusion_number": diagnostics["diffusion_number"],
            "thermal_diffusion_number": diagnostics["thermal_diffusion_number"],
            "thermodynamic_corrector_maximum_iterations_used": diagnostics[
                "thermodynamic_corrector_maximum_iterations_used"
            ],
            "thermodynamic_corrector_nonconverged_steps": diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ],
            "model_receptor_peak_h2_vol_percent": {
                label: float(np.max(values))
                for label, values in result.receptor_traces_vol_percent.items()
            },
        },
        "gate_results": gates,
        "all_predeclared_numerical_gates_passed": bool(all(gates.values())),
        "evidence_decision": {
            "status": (
                "conditional_preflight_passed"
                if all(gates.values())
                else "withheld_preflight_failed"
            ),
            "full_dynamic_component_run_allowed": bool(all(gates.values())),
            "obstacle_transport_model_qualification_allowed": False,
            "quantitative_facility_prediction_allowed": False,
        },
    }


def run_preflight(
    replay_input_path: str | Path,
    output_directory: str | Path,
    *,
    preflight_duration_s: float = 0.10,
) -> dict[str, Any]:
    """Freeze configuration, run the short preflight, and save diagnostics."""

    output = Path(output_directory).resolve()
    _, full_config, receptors = freeze_model_configuration(
        replay_input_path, output,
    )
    freeze_path = output / "model-configuration-freeze.json"
    preflight_config = _truncate_for_preflight(full_config, preflight_duration_s)
    result = solve_resolved_obstacle_flow_3d(
        preflight_config, receptors=receptors,
    )
    record = _preflight_record(
        result, freeze_path=freeze_path, config=preflight_config,
    )
    (output / "preflight-result.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return record


def _load_frozen_configuration(
    freeze_directory: str | Path,
) -> tuple[Path, dict[str, Any], ResolvedObstacleFlow3DConfig, tuple[FastObstacleReceptor3D, ...]]:
    directory = Path(freeze_directory).resolve()
    freeze_path = directory / "model-configuration-freeze.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("schema") != SCHEMA:
        raise ValueError("unsupported model-configuration freeze schema")
    control = freeze.get("observation_control", {})
    if control.get("hydrogen_observation_values_opened") is not False:
        raise ValueError("configuration freeze indicates H2 observations were opened")
    if control.get("model_configuration_frozen_before_observation_values_opened") is not True:
        raise ValueError("configuration was not frozen before H2 observation access")
    replay_path = Path(freeze["replay_input"]["path"]).resolve()
    if _sha256(replay_path) != freeze["replay_input"]["sha256"]:
        raise ValueError("dynamic replay input changed after configuration freeze")
    replay = _load_replay_input(replay_path)
    config, receptors, _ = _model_objects(replay)
    reconstructed = _serializable_config(config, receptors)
    if _canonical_sha256(reconstructed) != freeze["model_configuration_sha256"]:
        raise ValueError("reconstructed model configuration differs from freeze")
    normalized_reconstructed = json.loads(json.dumps(reconstructed, allow_nan=False))
    if normalized_reconstructed != freeze["model_configuration"]:
        raise ValueError("serialized model configuration differs from freeze")
    return freeze_path, freeze, config, receptors


def run_full_component(
    freeze_directory: str | Path,
) -> dict[str, Any]:
    """Run the complete frozen component replay without opening H2 values."""

    directory = Path(freeze_directory).resolve()
    result_path = directory / "full-run-result.json"
    array_path = directory / "full-run-model-output.npz"
    if result_path.exists() or array_path.exists():
        raise FileExistsError("full-run output already exists")
    freeze_path, freeze, config, receptors = _load_frozen_configuration(directory)
    result = solve_resolved_obstacle_flow_3d(config, receptors=receptors)
    diagnostics = result.diagnostics
    gates = {
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
    maximum_field = result.maximum_concentration_vol_percent
    final_field = result.concentration_vol_percent
    arrays: dict[str, np.ndarray] = {
        "x_m": result.x_m,
        "y_m": result.y_m,
        "z_m": result.z_m,
        "time_s": result.time_s,
        "final_h2_vol_percent": final_field,
        "maximum_h2_vol_percent": maximum_field,
        "final_temperature_k": result.temperature_k,
        "velocity_u_m_s": result.velocity_u_m_s,
        "velocity_v_m_s": result.velocity_v_m_s,
        "velocity_w_m_s": result.velocity_w_m_s,
        "injected_mass_kg": result.injected_mass_kg,
        "domain_mass_kg": result.domain_mass_kg,
        "outflow_mass_kg": result.outflow_mass_kg,
        "mass_residual_kg": result.mass_residual_kg,
    }
    for label, values in result.receptor_traces_vol_percent.items():
        arrays[f"model_receptor_{label}_h2_vol_percent"] = values
    for label, values in result.receptor_temperature_traces_k.items():
        arrays[f"model_receptor_{label}_temperature_k"] = values
    np.savez_compressed(array_path, **arrays)
    record = {
        "schema": "degali.preslhy-e34-dynamic-component-full-run.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "configuration_freeze": {
            "path": str(freeze_path),
            "sha256": _sha256(freeze_path),
            "model_configuration_sha256": freeze["model_configuration_sha256"],
        },
        "model_output": {
            "path": str(array_path),
            "sha256": _sha256(array_path),
        },
        "duration_s": config.duration_s,
        "steps": int(round(config.duration_s / config.scalar_time_step_s)),
        "observation_control": {
            "hydrogen_observation_values_opened": False,
            "hydrogen_observation_values_used": False,
            "model_output_completed_before_observation_values_opened": True,
        },
        "numerical_diagnostics": {
            "runtime_s": diagnostics["total_runtime_s"],
            "maximum_h2_vol_percent": float(np.nanmax(maximum_field)),
            "final_maximum_h2_vol_percent": float(np.nanmax(final_field)),
            "four_vol_percent_maximum_field_volume_m3": float(
                np.count_nonzero(maximum_field >= 4.0) * 0.1**3
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
            "active_scalar_out_of_range_cell_steps": diagnostics[
                "active_scalar_out_of_range_cell_steps"
            ],
            "maximum_lattice_mach_number": diagnostics[
                "maximum_lattice_mach_number"
            ],
            "maximum_variable_density_continuity_residual_kg_m3_s": diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
            "maximum_variable_density_global_residual_kg_s": diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "maximum_momentum_balance_residual_kg_m_s": diagnostics[
                "maximum_momentum_balance_residual_kg_m_s"
            ],
            "maximum_momentum_advective_courant_number": diagnostics[
                "maximum_momentum_advective_courant_number"
            ],
            "maximum_momentum_viscous_stability_number": diagnostics[
                "maximum_momentum_viscous_stability_number"
            ],
            "scalar_advective_number": diagnostics["advective_number"],
            "scalar_diffusion_number": diagnostics["diffusion_number"],
            "thermal_diffusion_number": diagnostics["thermal_diffusion_number"],
            "thermodynamic_corrector_maximum_iterations_used": diagnostics[
                "thermodynamic_corrector_maximum_iterations_used"
            ],
            "thermodynamic_corrector_nonconverged_steps": diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ],
            "final_injected_mass_kg": float(result.injected_mass_kg[-1]),
            "final_domain_mass_kg": float(result.domain_mass_kg[-1]),
            "final_outflow_mass_kg": float(result.outflow_mass_kg[-1]),
            "model_receptor_peak_h2_vol_percent": {
                label: float(np.max(values))
                for label, values in result.receptor_traces_vol_percent.items()
            },
        },
        "gate_results": gates,
        "all_predeclared_numerical_gates_passed": bool(all(gates.values())),
        "evidence_decision": {
            "status": (
                "conditional_model_output_ready_for_sealed_scoring"
                if all(gates.values())
                else "withheld_full_run_numerical_gate_failed"
            ),
            "sealed_observation_scoring_allowed": bool(all(gates.values())),
            "obstacle_transport_model_qualification_allowed": False,
            "quantitative_facility_prediction_allowed": False,
        },
    }
    result_path.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--replay-input",
        type=Path,
        default=(
            ROOT / "outputs" / "preslhy-e34-sand03-window01-replay-input-2026-10-11-v2"
            / "dynamic-replay-input.json"
        ),
    )
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--preflight-duration-s", type=float, default=0.10)
    parser.add_argument(
        "--run-full-component", action="store_true",
        help=(
            "use an existing output directory with a verified configuration "
            "freeze and run the complete observation-sealed component case"
        ),
    )
    arguments = parser.parse_args(argv)
    if arguments.run_full_component:
        record = run_full_component(arguments.output_directory)
    else:
        record = run_preflight(
            arguments.replay_input,
            arguments.output_directory,
            preflight_duration_s=arguments.preflight_duration_s,
        )
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if record["all_predeclared_numerical_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
