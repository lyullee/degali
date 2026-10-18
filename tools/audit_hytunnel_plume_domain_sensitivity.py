"""Same-grid HyTunnel plume-domain sensitivity audit.

The numerically admissible 6% handoff candidate is run in the published-size
box and in a larger box at identical 45 mm spacing, time step, source physics
and common source-relative receptors.  The comparison is numerical boundary
evidence only: it neither fits nor scores measured concentration.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import time
from typing import Any

import numpy as np

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.high_pressure_source_handoff_3d import (
    low_mach_handoff_from_energy_nozzle,
)
from degali.addons.notional import energy_conserving_notional_nozzle
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.validation.hytunnel_ps import (
    DATASET_DOI,
    read_hytunnel_flow_workbook,
)


SCHEMA = "degali.hytunnel-plume-domain-sensitivity.v7"
CASE_SCHEMA = "degali.hytunnel-resolved-case.v2"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FLOW_WORKBOOK = (
    REPOSITORY_ROOT
    / "reference/datasets/hytunnel-ps-4668554/raw/HTE2440PS005FLMT00300321.xlsx"
)
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-plume-domain-sensitivity-2026-10-10-v3.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-plume-domain-sensitivity-2026-10-10-v3.png"
)
GRID_SPACING_M = 0.045
BASE_DOMAIN_M = (2.52, 0.90, 0.90)
BASE_CELLS = (56, 20, 20)
EXPANDED_DOMAIN_M = (3.15, 1.35, 1.35)
EXPANDED_CELLS = (70, 30, 30)
WIND_SPEED_M_S = 3.5
LATTICE_INLET_SPEED = 0.018
ACTIVE_TIME_STEP_S = (
    LATTICE_INLET_SPEED * GRID_SPACING_M / WIND_SPEED_M_S
)
TIME_STEPS = 520
DURATION_S = TIME_STEPS * ACTIVE_TIME_STEP_S
FLOW_STEPS = 100
PRIMARY_RELATIVE_CHANGE_LIMIT = 0.10
COMMON_FIELD_NORMALISED_L1_LIMIT = 0.10
DEVELOPMENT_RUNTIME_TARGET_S = 300.0
OPEN_FACE_MAXIMUM_CONCENTRATION_LIMIT_VOL_PERCENT = 0.10
FOUR_VOL_MINIMUM_OPEN_BOUNDARY_CLEARANCE_CELLS = 3
FAN_INLET_PROFILES = ("uniform", "measured_first_station")


def _handoff(maximum_density_departure: float = 0.06):
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen",
        storage_temperature=286.15,
        storage_pressure=11.515e6,
        mass_flow=0.005,
        orifice_diameter=0.001,
        ambient_pressure=98_500.0,
        allow_supercritical_gas=True,
    )
    return low_mach_handoff_from_energy_nozzle(
        nozzle,
        source_id="hytunnel-1mm-5g-3p5ms-domain-rho0.060",
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_axial_velocity_m_s=WIND_SPEED_M_S,
        maximum_mach_number=0.30,
        maximum_relative_density_departure=maximum_density_departure,
    )


def _case_specification(
    name: str, maximum_density_departure: float = 0.06,
) -> dict[str, Any]:
    if name == "published_box":
        domain = BASE_DOMAIN_M
        cells = BASE_CELLS
    elif name == "expanded_box":
        domain = EXPANDED_DOMAIN_M
        cells = EXPANDED_CELLS
    else:
        raise ValueError(f"unsupported domain case {name!r}")
    centre_y = 0.5 * domain[1]
    centre_z = 0.5 * domain[2]
    handoff = _handoff(maximum_density_departure)
    source_x = handoff.handoff_distance_m
    receptors = (
        ("centre_0p30", (source_x + 0.30, centre_y, centre_z)),
        ("centre_0p75", (source_x + 0.75, centre_y, centre_z)),
        ("centre_1p20", (source_x + 1.20, centre_y, centre_z)),
        ("side_0p75", (source_x + 0.75, centre_y + 0.225, centre_z)),
        ("upper_0p75", (source_x + 0.75, centre_y, centre_z + 0.225)),
    )
    return {
        "name": name,
        "domain_m": domain,
        "cells": cells,
        "centre_y_m": centre_y,
        "centre_z_m": centre_z,
        "source_x_m": source_x,
        "maximum_density_departure": maximum_density_departure,
        "receptors": receptors,
    }


def _fan_inlet_profile(
    specification: dict[str, Any],
    selection: str,
    workbook_path: str | Path = DEFAULT_FLOW_WORKBOOK,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return an active FV inlet declaration and its immutable provenance.

    The public workbook contains one-minute mean/minimum/maximum velocity
    matrices, not a velocity time series.  Only the first measured mean plane
    is transferred here; it does not define turbulence intensity or spectrum.
    """

    if selection not in FAN_INLET_PROFILES:
        raise ValueError(f"unsupported fan inlet profile {selection!r}")
    if selection == "uniform":
        return {}, {
            "selection": selection,
            "model": "uniform_declared_nominal_velocity",
            "nominal_velocity_m_s": WIND_SPEED_M_S,
            "turbulence_statistics_available": False,
        }

    source = Path(workbook_path)
    dataset = read_hytunnel_flow_workbook(source)
    configuration_id = dataset.configuration_ids[0]
    average = tuple(
        point for point in dataset.points
        if point.configuration_id == configuration_id
        and point.statistic == "average"
    )
    if not average:
        raise ValueError("HyTunnel first fan-flow sheet has no average points")
    first_x = min(point.x_m for point in average)
    first_plane = sorted(
        (
            point for point in average
            if math.isclose(point.x_m, first_x, abs_tol=1.0e-12)
        ),
        key=lambda point: point.y_m,
    )
    relative_y = tuple(float(point.y_m) for point in first_plane)
    velocity = tuple(float(point.velocity_x_m_s) for point in first_plane)
    if len(relative_y) < 2 or any(
        right <= left for left, right in zip(relative_y, relative_y[1:])
    ):
        raise ValueError("HyTunnel first fan-flow mean plane is not ordered")
    absolute_y = tuple(
        float(specification["centre_y_m"]) + value for value in relative_y
    )
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    return {
        "active_inlet_profile_y_m": absolute_y,
        "active_inlet_profile_u_m_s": velocity,
    }, {
        "selection": selection,
        "model": "first_measured_one_minute_mean_plane_transferred_to_x0",
        "dataset_doi": DATASET_DOI,
        "workbook": source.name,
        "workbook_sha256": digest,
        "configuration_id": configuration_id,
        "sheet_index": 1,
        "statistic": "average",
        "measured_x_m": first_x,
        "relative_y_m": list(relative_y),
        "velocity_x_m_s": list(velocity),
        "transfer_to_computational_x_m": 0.0,
        "outside_measured_lateral_range_policy": "constant_endpoint_extension",
        "turbulence_statistics_available": False,
        "missing_turbulence_inputs": [
            "velocity_time_series", "rms_velocity", "turbulence_intensity",
            "integral_length_scale", "spectrum", "dissipation_rate",
        ],
    }


def _case_config(
    specification: dict[str, Any], *, time_steps: int = TIME_STEPS,
    effective_schmidt_number: float = 0.7,
    turbulent_diffusivity_model: str = (
        "legacy_total_kinematic_viscosity_ratio"
    ),
    fan_inlet_profile: str = "uniform",
    flow_workbook_path: str | Path = DEFAULT_FLOW_WORKBOOK,
) -> ResolvedObstacleFlow3DConfig:
    if not isinstance(time_steps, int) or time_steps < 1:
        raise ValueError("time_steps must be a positive integer")
    handoff = _handoff(float(specification["maximum_density_departure"]))
    dt = ACTIVE_TIME_STEP_S
    duration_s = time_steps * dt
    source = handoff.resolved_source_fields(
        nozzle_position_m=(
            0.0,
            float(specification["centre_y_m"]),
            float(specification["centre_z_m"]),
        ),
        grid_spacing_m=GRID_SPACING_M,
        duration_s=duration_s,
    )
    fan_kwargs, _fan_evidence = _fan_inlet_profile(
        specification, fan_inlet_profile, flow_workbook_path,
    )
    return ResolvedObstacleFlow3DConfig(
        domain_m=tuple(specification["domain_m"]),
        cells=tuple(specification["cells"]),
        wind_speed_m_s=WIND_SPEED_M_S,
        kinematic_viscosity_m2_s=0.04,
        flow_steps=FLOW_STEPS,
        minimum_flow_steps=FLOW_STEPS,
        lattice_inlet_speed=LATTICE_INLET_SPEED,
        duration_s=duration_s,
        scalar_time_step_s=dt,
        scalar_diffusivity_m2_s=2.2e-5,
        thermal_diffusivity_m2_s=2.2e-5,
        obstacles=(),
        scalar_advection_scheme="muscl_van_leer",
        source_time_integration="prestep",
        subgrid_model="wale",
        wale_constant=0.5,
        flow_boundary_scheme="characteristic_multidimensional_extrapolated",
        nonstreamwise_wall_boundary="all_pressure_open",
        effective_schmidt_number=effective_schmidt_number,
        effective_prandtl_number=0.85,
        scalar_turbulent_diffusivity_model=turbulent_diffusivity_model,
        thermal_turbulent_diffusivity_model=turbulent_diffusivity_model,
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_air_density_kg_m3=handoff.ambient_density_kg_m3,
        maximum_active_scalar_density_departure=0.10,
        maximum_lattice_mach_number=0.25,
        variable_density_projection_linear_solver="cached_direct",
        store_scalar_fields=False,
        **fan_kwargs,
        **source.as_config_kwargs(),
    )


def _trace_metrics(time_s: np.ndarray, values: np.ndarray) -> dict[str, Any]:
    peak = float(np.max(values))
    dose = float(np.trapezoid(values, time_s))
    above = values >= 4.0
    arrival = None
    if np.any(above):
        arrival = float(time_s[int(np.argmax(above))])
    duration = float(np.sum(above[1:] * np.diff(time_s)))
    return {
        "peak_vol_percent": peak,
        "dose_vol_percent_s": dose,
        "arrival_4vol_s": arrival,
        "duration_above_4vol_s": duration,
    }


def _common_crop(array: np.ndarray, specification: dict[str, Any]) -> np.ndarray:
    if specification["name"] == "published_box":
        return np.asarray(array, dtype=float)
    # The expanded box adds five cells on either lateral/vertical side and
    # fourteen cells downstream.  This slice is exactly source-relative and
    # cell-centre aligned with the published box.
    return np.asarray(array, dtype=float)[5:25, 5:25, :56]


def _boundary_clearance(maximum: np.ndarray) -> dict[str, Any]:
    field = np.asarray(maximum, dtype=float)
    if field.ndim != 3 or np.any(~np.isfinite(field)):
        raise ValueError("maximum concentration field must be finite and 3-D")
    face_maximum = {
        "x_max": float(np.max(field[:, :, -1])),
        "y_min": float(np.max(field[:, 0, :])),
        "y_max": float(np.max(field[:, -1, :])),
        "z_min": float(np.max(field[0, :, :])),
        "z_max": float(np.max(field[-1, :, :])),
    }
    above = np.argwhere(field >= 4.0)
    clearances = {
        "x_max": None,
        "y_min": None,
        "y_max": None,
        "z_min": None,
        "z_max": None,
    }
    if above.size:
        nz, ny, nx = field.shape
        distances = {
            "x_max": nx - 1 - above[:, 2],
            "y_min": above[:, 1],
            "y_max": ny - 1 - above[:, 1],
            "z_min": above[:, 0],
            "z_max": nz - 1 - above[:, 0],
        }
        clearances = {
            face: int(np.min(values))
            for face, values in distances.items()
        }
    lateral_faces = ("y_min", "y_max", "z_min", "z_max")
    lateral_clearance = None
    if above.size:
        lateral_clearance = min(
            int(clearances[face]) for face in lateral_faces
        )
    face_passed = max(face_maximum[face] for face in lateral_faces) <= (
        OPEN_FACE_MAXIMUM_CONCENTRATION_LIMIT_VOL_PERCENT
    )
    four_vol_passed = bool(
        lateral_clearance is None
        or lateral_clearance
        >= FOUR_VOL_MINIMUM_OPEN_BOUNDARY_CLEARANCE_CELLS
    )
    return {
        "open_face_maximum_concentration_vol_percent": face_maximum,
        "maximum_open_face_concentration_vol_percent": max(
            face_maximum.values()
        ),
        "maximum_lateral_open_face_concentration_vol_percent": max(
            face_maximum[face] for face in lateral_faces
        ),
        "open_face_maximum_concentration_limit_vol_percent": (
            OPEN_FACE_MAXIMUM_CONCENTRATION_LIMIT_VOL_PERCENT
        ),
        "minimum_4vol_open_boundary_clearance_cells_by_face": clearances,
        "minimum_4vol_lateral_open_boundary_clearance_cells": (
            lateral_clearance
        ),
        "minimum_4vol_open_boundary_clearance_limit_cells": (
            FOUR_VOL_MINIMUM_OPEN_BOUNDARY_CLEARANCE_CELLS
        ),
        "open_face_concentration_gate_passed": face_passed,
        "four_vol_boundary_clearance_gate_passed": four_vol_passed,
        "boundary_clearance_gate_passed": bool(
            face_passed and four_vol_passed
        ),
        "x_max_policy": (
            "advective_outlet_evaluated_by_domain_sensitivity_and_mass_outflow_"
            "not_zero_concentration"
        ),
    }


def _run_case(
    name: str,
    time_steps: int = TIME_STEPS,
    maximum_density_departure: float = 0.06,
    effective_schmidt_number: float = 0.7,
    turbulent_diffusivity_model: str = (
        "legacy_total_kinematic_viscosity_ratio"
    ),
    fan_inlet_profile: str = "uniform",
    flow_workbook_path: str | Path = DEFAULT_FLOW_WORKBOOK,
) -> dict[str, Any]:
    specification = _case_specification(name, maximum_density_departure)
    _fan_kwargs, fan_evidence = _fan_inlet_profile(
        specification, fan_inlet_profile, flow_workbook_path,
    )
    config = _case_config(
        specification,
        time_steps=time_steps,
        effective_schmidt_number=effective_schmidt_number,
        turbulent_diffusivity_model=turbulent_diffusivity_model,
        fan_inlet_profile=fan_inlet_profile,
        flow_workbook_path=flow_workbook_path,
    )
    receptors = tuple(
        FastObstacleReceptor3D(label, tuple(position))
        for label, position in specification["receptors"]
    )
    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(config, receptors=receptors)
    full_maximum = np.asarray(
        result.maximum_concentration_vol_percent, dtype=float,
    )
    maximum = _common_crop(
        full_maximum, specification,
    )
    final = _common_crop(result.concentration_vol_percent, specification)
    centre_z = BASE_CELLS[2] // 2 - 1
    cell_volume = GRID_SPACING_M**3
    diagnostics = result.diagnostics
    return {
        "name": name,
        "domain_m": list(specification["domain_m"]),
        "cells": list(specification["cells"]),
        "source_position_m": list(config.source_position_m),
        "fan_inlet_evidence": fan_evidence,
        "time_s": result.time_s.tolist(),
        "receptor_traces_vol_percent": {
            label: values.tolist()
            for label, values in result.receptor_traces_vol_percent.items()
        },
        "receptor_metrics": {
            label: _trace_metrics(
                result.time_s, result.receptor_traces_vol_percent[label],
            )
            for label in result.receptor_traces_vol_percent
        },
        "common_maximum_field": maximum.tolist(),
        "common_final_field": final.tolist(),
        "common_maximum_centre_slice": maximum[centre_z].tolist(),
        "common_final_centre_slice": final[centre_z].tolist(),
        "common_maximum_hydrogen_vol_percent": float(np.max(maximum)),
        "common_4vol_maximum_volume_m3": float(
            np.count_nonzero(maximum >= 4.0) * cell_volume
        ),
        "final_domain_hydrogen_mass_kg": float(result.domain_mass_kg[-1]),
        "cumulative_outflow_hydrogen_mass_kg": float(
            result.outflow_mass_kg[-1]
        ),
        "injected_hydrogen_mass_kg": float(result.injected_mass_kg[-1]),
        "maximum_relative_mass_residual": float(
            diagnostics["maximum_relative_mass_residual"]
        ),
        "maximum_relative_energy_residual": float(
            diagnostics["maximum_relative_energy_residual"]
        ),
        "maximum_relative_mixture_density_departure": float(
            diagnostics["maximum_relative_mixture_density_departure"]
        ),
        "maximum_lattice_mach_number": float(
            diagnostics["maximum_lattice_mach_number"]
        ),
        "active_scalar_applicability_passed": bool(
            diagnostics["active_scalar_applicability_passed"]
        ),
        "scalar_diffusivity_model": diagnostics["scalar_diffusivity_model"],
        "applied_scalar_diffusivity_minimum_m2_s": diagnostics[
            "applied_scalar_diffusivity_minimum_m2_s"
        ],
        "applied_scalar_diffusivity_maximum_m2_s": diagnostics[
            "applied_scalar_diffusivity_maximum_m2_s"
        ],
        "applied_scalar_diffusivity_mean_m2_s": diagnostics[
            "applied_scalar_diffusivity_mean_m2_s"
        ],
        "active_inlet_profile_model": diagnostics[
            "active_inlet_profile_model"
        ],
        "active_inlet_plane_minimum_m_s": diagnostics[
            "active_inlet_plane_minimum_m_s"
        ],
        "active_inlet_plane_maximum_m_s": diagnostics[
            "active_inlet_plane_maximum_m_s"
        ],
        "active_inlet_plane_mean_m_s": diagnostics[
            "active_inlet_plane_mean_m_s"
        ],
        "boundary_clearance": _boundary_clearance(full_maximum),
        "variable_density_linear_solver": diagnostics[
            "variable_density_linear_solver"
        ],
        "variable_density_projection_warm_start": diagnostics[
            "variable_density_projection_warm_start"
        ],
        "variable_density_projection_mean_linear_iterations_per_step": (
            diagnostics[
                "variable_density_projection_mean_linear_iterations_per_step"
            ]
        ),
        "variable_density_projection_maximum_linear_iterations_per_step": (
            diagnostics[
                "variable_density_projection_maximum_linear_iterations_per_step"
            ]
        ),
        "runtime_s": time.perf_counter() - started,
    }


def _relative_change(a: float, b: float, floor: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), floor)


def classify_domain_sensitivity(
    published: dict[str, Any], expanded: dict[str, Any],
) -> dict[str, Any]:
    published_maximum = np.asarray(published["common_maximum_field"], dtype=float)
    expanded_maximum = np.asarray(expanded["common_maximum_field"], dtype=float)
    if published_maximum.shape != expanded_maximum.shape:
        raise ValueError("common maximum fields must have identical shapes")
    normaliser = max(float(np.sum(np.abs(expanded_maximum))), 1.0e-12)
    field_l1 = float(
        np.sum(np.abs(published_maximum - expanded_maximum)) / normaliser
    )
    primary: dict[str, float] = {
        "common_maximum_hydrogen_vol_percent": _relative_change(
            float(published["common_maximum_hydrogen_vol_percent"]),
            float(expanded["common_maximum_hydrogen_vol_percent"]),
            0.1,
        ),
        "common_4vol_maximum_volume_m3": _relative_change(
            float(published["common_4vol_maximum_volume_m3"]),
            float(expanded["common_4vol_maximum_volume_m3"]),
            GRID_SPACING_M**3,
        ),
    }
    for label, expanded_metrics in expanded["receptor_metrics"].items():
        published_metrics = published["receptor_metrics"][label]
        if float(expanded_metrics["peak_vol_percent"]) < 1.0:
            continue
        primary[f"{label}.peak_vol_percent"] = _relative_change(
            float(published_metrics["peak_vol_percent"]),
            float(expanded_metrics["peak_vol_percent"]),
            0.1,
        )
        primary[f"{label}.dose_vol_percent_s"] = _relative_change(
            float(published_metrics["dose_vol_percent_s"]),
            float(expanded_metrics["dose_vol_percent_s"]),
            1.0e-4,
        )
        if expanded_metrics["arrival_4vol_s"] is not None:
            if published_metrics["arrival_4vol_s"] is None:
                primary[f"{label}.arrival_4vol_s"] = math.inf
            else:
                primary[f"{label}.arrival_4vol_s"] = _relative_change(
                    float(published_metrics["arrival_4vol_s"]),
                    float(expanded_metrics["arrival_4vol_s"]),
                    float(expanded["time_s"][1]),
                )
            primary[f"{label}.duration_above_4vol_s"] = _relative_change(
                float(published_metrics["duration_above_4vol_s"]),
                float(expanded_metrics["duration_above_4vol_s"]),
                float(expanded["time_s"][1]),
            )
    worst_name, worst_change = max(primary.items(), key=lambda item: item[1])
    conservation_passed = all(
        float(case["maximum_relative_mass_residual"]) <= 1.0e-6
        and float(case["maximum_relative_energy_residual"]) <= 1.0e-6
        for case in (published, expanded)
    )
    applicability_passed = all(
        bool(case["active_scalar_applicability_passed"])
        for case in (published, expanded)
    )
    comparison_boundary_clearance_passed = bool(
        published["boundary_clearance"]["boundary_clearance_gate_passed"]
    )
    boundary_clearance_passed = bool(
        expanded["boundary_clearance"]["boundary_clearance_gate_passed"]
    )
    passed = bool(
        worst_change <= PRIMARY_RELATIVE_CHANGE_LIMIT
        and field_l1 <= COMMON_FIELD_NORMALISED_L1_LIMIT
        and conservation_passed
        and applicability_passed
        and boundary_clearance_passed
    )
    return {
        "status": (
            "same_grid_plume_domain_sensitivity_gate_passed"
            if passed else "same_grid_plume_domain_sensitivity_gate_failed"
        ),
        "passed": passed,
        "promotion_allowed": False,
        "primary_relative_changes": primary,
        "maximum_primary_relative_change": float(worst_change),
        "maximum_primary_relative_change_metric": worst_name,
        "primary_relative_change_limit": PRIMARY_RELATIVE_CHANGE_LIMIT,
        "common_maximum_field_normalised_l1": field_l1,
        "common_field_normalised_l1_limit": COMMON_FIELD_NORMALISED_L1_LIMIT,
        "conservation_gate_passed": conservation_passed,
        "active_scalar_applicability_gate_passed": applicability_passed,
        "boundary_clearance_gate_passed": boundary_clearance_passed,
        "smaller_comparison_domain_boundary_clearance_passed": (
            comparison_boundary_clearance_passed
        ),
        "boundary_clearance_policy": (
            "larger_reference_domain_must_clear_lateral_y_z_faces; "
            "x_max_is_an_intentional_advective_outlet_and_is_tested_by_"
            "same_grid_domain_sensitivity_and_mass_conservation"
        ),
    }


def run_audit(
    workers: int = 2, time_steps: int = TIME_STEPS,
    maximum_density_departure: float = 0.06,
    effective_schmidt_number: float = 0.7,
    turbulent_diffusivity_model: str = (
        "legacy_total_kinematic_viscosity_ratio"
    ),
    fan_inlet_profile: str = "uniform",
    flow_workbook_path: str | Path = DEFAULT_FLOW_WORKBOOK,
) -> dict[str, Any]:
    if not isinstance(time_steps, int) or time_steps < 1:
        raise ValueError("time_steps must be a positive integer")
    started = time.perf_counter()
    names = ("published_box", "expanded_box")
    if workers == 1:
        cases = [
            _run_case(
                name, time_steps, maximum_density_departure,
                effective_schmidt_number, turbulent_diffusivity_model,
                fan_inlet_profile, flow_workbook_path,
            )
            for name in names
        ]
    else:
        with ProcessPoolExecutor(max_workers=min(workers, 2)) as executor:
            futures = [
                executor.submit(
                    _run_case, name, time_steps, maximum_density_departure,
                    effective_schmidt_number, turbulent_diffusivity_model,
                    fan_inlet_profile, flow_workbook_path,
                )
                for name in names
            ]
            cases = [future.result() for future in futures]
    by_name = {case["name"]: case for case in cases}
    classification = classify_domain_sensitivity(
        by_name["published_box"], by_name["expanded_box"],
    )
    # Full 3-D fields served their comparison purpose and are omitted from the
    # persisted JSON; the aligned centre slices and scalar metrics remain.
    for case in cases:
        case.pop("common_maximum_field")
        case.pop("common_final_field")
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **classification,
        "configuration": {
            "grid_spacing_m": GRID_SPACING_M,
            "duration_s": time_steps * ACTIVE_TIME_STEP_S,
            "scalar_time_step_s": (
                ACTIVE_TIME_STEP_S
            ),
            "time_steps": time_steps,
            "flow_steps": FLOW_STEPS,
            "source_candidate_mean_density_departure_limit": (
                maximum_density_departure
            ),
            "effective_schmidt_number": effective_schmidt_number,
            "turbulent_diffusivity_model": turbulent_diffusivity_model,
            "fan_inlet_profile": fan_inlet_profile,
            "flow_boundary_scheme": (
                "characteristic_multidimensional_extrapolated"
            ),
            "primary_receptor_inclusion_peak_floor_vol_percent": 1.0,
            "common_region_domain_m": list(BASE_DOMAIN_M),
        },
        "runtime_evidence": {
            "development_target_s": DEVELOPMENT_RUNTIME_TARGET_S,
            "published_box_target_passed": bool(
                by_name["published_box"]["runtime_s"]
                <= DEVELOPMENT_RUNTIME_TARGET_S
            ),
            "expanded_box_target_passed": bool(
                by_name["expanded_box"]["runtime_s"]
                <= DEVELOPMENT_RUNTIME_TARGET_S
            ),
            "assessment": (
                "expanded_box_meets_5min_development_target_but_is_not_"
                "the_converged_reference_case"
                if by_name["expanded_box"]["runtime_s"]
                <= DEVELOPMENT_RUNTIME_TARGET_S
                else
                "expanded_box_misses_5min_development_target_but_is_not_"
                "the_converged_reference_case"
            ),
            "machine": {
                "platform": platform.platform(),
                "machine": platform.machine(),
                "processor": platform.processor(),
                "logical_cpu_count": os.cpu_count(),
            },
        },
        "cases": cases,
        "scope": (
            "same_grid_two_domain_numerical_sensitivity_not_external_accuracy_or_grid_convergence"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def _receptor_tail_stationarity(case: dict[str, Any]) -> dict[str, Any]:
    """Measure endpoint stability over the last 20% of each receptor trace."""

    by_receptor: dict[str, float] = {}
    for label, values in case["receptor_traces_vol_percent"].items():
        trace = np.asarray(values, dtype=float)
        start = int(math.floor(0.8 * trace.size))
        tail = trace[start:]
        scale = max(abs(float(np.mean(tail))), 0.10)
        by_receptor[label] = float(np.ptp(tail) / scale)
    maximum = max(by_receptor.values(), default=math.inf)
    return {
        "window": "last_20_percent_of_transient",
        "relative_peak_to_peak_by_receptor": by_receptor,
        "maximum_relative_peak_to_peak": maximum,
        "relative_peak_to_peak_limit": 0.01,
        "passed": bool(maximum <= 0.01),
    }


def run_single_case(
    name: str,
    time_steps: int = TIME_STEPS,
    maximum_density_departure: float = 0.06,
    effective_schmidt_number: float = 0.7,
    turbulent_diffusivity_model: str = (
        "legacy_total_kinematic_viscosity_ratio"
    ),
    fan_inlet_profile: str = "uniform",
    flow_workbook_path: str | Path = DEFAULT_FLOW_WORKBOOK,
) -> dict[str, Any]:
    """Run and persist one resolved case with its final observation slice."""

    started = time.perf_counter()
    case = _run_case(
        name,
        time_steps,
        maximum_density_departure,
        effective_schmidt_number,
        turbulent_diffusivity_model,
        fan_inlet_profile,
        flow_workbook_path,
    )
    stationarity = _receptor_tail_stationarity(case)
    case.pop("common_maximum_field")
    case.pop("common_final_field")
    return {
        "schema": CASE_SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "resolved_case_endpoint_stationarity_passed"
            if stationarity["passed"]
            else "resolved_case_endpoint_stationarity_failed"
        ),
        "configuration": {
            "grid_spacing_m": GRID_SPACING_M,
            "duration_s": time_steps * ACTIVE_TIME_STEP_S,
            "scalar_time_step_s": ACTIVE_TIME_STEP_S,
            "time_steps": time_steps,
            "flow_steps": FLOW_STEPS,
            "source_candidate_mean_density_departure_limit": (
                maximum_density_departure
            ),
            "effective_schmidt_number": effective_schmidt_number,
            "turbulent_diffusivity_model": turbulent_diffusivity_model,
            "fan_inlet_profile": fan_inlet_profile,
            "applied_scalar_diffusivity_m2_s": (
                0.04 / effective_schmidt_number
                if turbulent_diffusivity_model
                == "legacy_total_kinematic_viscosity_ratio"
                else None
            ),
            "flow_boundary_scheme": "characteristic_multidimensional_extrapolated",
            "common_region_domain_m": list(BASE_DOMAIN_M),
        },
        "endpoint_stationarity": stationarity,
        "case": case,
        "scope": (
            "single_case_endpoint_and_maximum_observation_slices_not_domain_"
            "sensitivity_or_external_validation"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    cases = {case["name"]: case for case in record["cases"]}
    published = np.asarray(
        cases["published_box"]["common_maximum_centre_slice"], dtype=float,
    )
    expanded = np.asarray(
        cases["expanded_box"]["common_maximum_centre_slice"], dtype=float,
    )
    difference = published - expanded
    extent = [0.0, BASE_DOMAIN_M[0], -0.45, 0.45]
    maximum = max(float(np.max(published)), float(np.max(expanded)), 4.0)
    figure, axes = plt.subplots(
        2, 2, figsize=(11.0, 7.6), constrained_layout=True,
    )
    for axis, field, title in (
        (axes[0, 0], published, "Published-size box"),
        (axes[0, 1], expanded, "Expanded box, common crop"),
    ):
        image = axis.imshow(
            field, origin="lower", extent=extent, aspect="auto",
            cmap="turbo", vmin=0.0, vmax=maximum,
        )
        axis.contour(
            field, levels=[4.0], colors="white", linewidths=1.0,
            origin="lower", extent=extent,
        )
        axis.set_title(title)
        axis.set_xlabel("x (m)")
        axis.set_ylabel("source-relative y (m)")
    difference_limit = max(float(np.max(np.abs(difference))), 1.0e-9)
    diff_image = axes[1, 0].imshow(
        difference, origin="lower", extent=extent, aspect="auto",
        cmap="coolwarm", vmin=-difference_limit, vmax=difference_limit,
    )
    axes[1, 0].set_title("Published - expanded (vol%p)")
    axes[1, 0].set_xlabel("x (m)")
    axes[1, 0].set_ylabel("source-relative y (m)")
    labels = list(cases["expanded_box"]["receptor_metrics"])
    locations = np.arange(len(labels))
    width = 0.38
    axes[1, 1].bar(
        locations - width / 2,
        [
            cases["published_box"]["receptor_metrics"][label][
                "peak_vol_percent"
            ]
            for label in labels
        ],
        width, label="published-size",
    )
    axes[1, 1].bar(
        locations + width / 2,
        [
            cases["expanded_box"]["receptor_metrics"][label][
                "peak_vol_percent"
            ]
            for label in labels
        ],
        width, label="expanded",
    )
    axes[1, 1].set_xticks(locations, labels, rotation=30, ha="right")
    axes[1, 1].set_ylabel("Peak H2 (vol%)")
    axes[1, 1].set_title("Common receptor peaks")
    axes[1, 1].legend(frameon=False)
    figure.colorbar(image, ax=axes[0, :], label="Maximum H2 (vol%)")
    figure.colorbar(diff_image, ax=axes[1, 0], label="Difference (vol%p)")
    figure.suptitle(
        "HyTunnel 6% candidate: same-grid plume-domain sensitivity\n"
        f"{record['status']}",
        fontsize=11,
    )
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--time-steps", type=int, default=TIME_STEPS)
    parser.add_argument(
        "--case-only", choices=("published_box", "expanded_box"),
        help="run one case and persist final/maximum centre slices",
    )
    parser.add_argument(
        "--handoff-density-departure-limit", type=float, default=0.06,
        help="mean-plane low-Mach handoff density-departure limit",
    )
    parser.add_argument(
        "--effective-schmidt-number", type=float, default=0.7,
        help="effective momentum-to-species diffusivity ratio",
    )
    parser.add_argument(
        "--turbulent-diffusivity-model",
        choices=(
            "legacy_total_kinematic_viscosity_ratio",
            "molecular_plus_subgrid_eddy_viscosity",
        ),
        default="legacy_total_kinematic_viscosity_ratio",
        help="legacy uniform surrogate or local SGS gradient diffusion",
    )
    parser.add_argument(
        "--fan-inlet-profile",
        choices=FAN_INLET_PROFILES,
        default="uniform",
        help=(
            "uniform nominal inlet or the public workbook's first measured "
            "one-minute mean fan plane"
        ),
    )
    parser.add_argument(
        "--flow-workbook", type=Path, default=DEFAULT_FLOW_WORKBOOK,
        help="HyTunnel fan-flow workbook used by measured_first_station",
    )
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.time_steps < 1:
        parser.error("--time-steps must be positive")
    if not 0.0 < args.handoff_density_departure_limit <= 0.10:
        parser.error("--handoff-density-departure-limit must lie in (0, 0.10]")
    if args.effective_schmidt_number <= 0.0:
        parser.error("--effective-schmidt-number must be positive")
    if args.case_only:
        record = run_single_case(
            args.case_only,
            args.time_steps,
            args.handoff_density_departure_limit,
            args.effective_schmidt_number,
            args.turbulent_diffusivity_model,
            args.fan_inlet_profile,
            args.flow_workbook,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(record, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({
            "status": record["status"],
            "maximum_tail_relative_peak_to_peak": record[
                "endpoint_stationarity"
            ]["maximum_relative_peak_to_peak"],
            "runtime_s": record["runtime_s"],
        }, indent=2, allow_nan=False))
        return 0 if record["endpoint_stationarity"]["passed"] else 1
    record = run_audit(
        args.workers,
        args.time_steps,
        args.handoff_density_departure_limit,
        args.effective_schmidt_number,
        args.turbulent_diffusivity_model,
        args.fan_inlet_profile,
        args.flow_workbook,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot_audit(record, args.plot)
    print(json.dumps({
        key: record[key]
        for key in (
            "status", "maximum_primary_relative_change",
            "maximum_primary_relative_change_metric",
            "common_maximum_field_normalised_l1", "runtime_s",
        )
    }, indent=2, allow_nan=False))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
