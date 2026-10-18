"""Grid-convergence audit for the resolved two-dimensional obstacle tier."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from degali.addons.fast_obstacle_transport import (
    HYDROGEN_MOLAR_MASS_KG_MOL,
    R_GAS_J_MOL_K,
    FastObstacleReceptor,
)
from degali.addons.resolved_obstacle_flow_2d import (
    ResolvedObstacleFlow2DConfig,
    ResolvedObstacleFlow2DResult,
    solve_resolved_obstacle_flow_2d,
)


RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA = (
    "degali.resolved-obstacle-grid-convergence.v1"
)


@dataclass(frozen=True)
class ObstacleReceptorConvergenceMetric:
    label: str
    peak_vol_percent: float
    dose_vol_percent_s: float
    arrival_time_4vol_s: float | None
    duration_above_4vol_s: float


@dataclass(frozen=True)
class ResolvedObstacleGridCase:
    nx: int
    nz: int
    cell_size_m: float
    flow_physical_duration_s: float
    runtime_s: float
    maximum_vol_percent: float
    four_vol_percent_envelope_area_per_span_m2: float
    four_vol_percent_maximum_downwind_x_m: float | None
    minimum_lee_streamwise_velocity_m_s: float
    maximum_absolute_divergence_s_1: float
    maximum_mass_residual_per_span_kg_m: float
    receptor_metrics: tuple[ObstacleReceptorConvergenceMetric, ...]
    result: ResolvedObstacleFlow2DResult


@dataclass(frozen=True)
class ResolvedObstacleGridConvergenceResult:
    schema: str
    cases: tuple[ResolvedObstacleGridCase, ...]
    relative_tolerance: float
    concentration_floor_vol_percent: float
    dose_floor_vol_percent_s: float
    arrival_time_tolerance_s: float
    threshold_duration_floor_s: float
    finest_adjacent_relative_changes: dict[str, object]
    passed: bool
    gate_reasons: tuple[str, ...]
    diagnostics: dict[str, object]


def _conversion(config: ResolvedObstacleFlow2DConfig) -> float:
    molar_density = config.ambient_pressure_pa / (
        R_GAS_J_MOL_K * config.ambient_temperature_k
    )
    return 100.0 / (HYDROGEN_MOLAR_MASS_KG_MOL * molar_density)


def _threshold_metrics(
    time_s: np.ndarray,
    values: np.ndarray,
    threshold: float,
) -> tuple[float | None, float]:
    above = values >= threshold
    arrival = None if not np.any(above) else float(time_s[np.flatnonzero(above)[0]])
    duration = float(np.trapezoid(above.astype(float), time_s))
    return arrival, duration


def _scenario_signature(config: ResolvedObstacleFlow2DConfig) -> tuple[object, ...]:
    return (
        config.length_m, config.height_m, config.wind_speed_m_s,
        config.kinematic_viscosity_m2_s, config.lattice_inlet_speed,
        config.duration_s, config.scalar_diffusivity_m2_s,
        config.source_schedule, config.source_x_m, config.source_height_m,
        config.source_sigma_m, config.obstacle,
        config.inlet_profile_z_m, config.inlet_profile_u_m_s,
        config.ambient_temperature_k,
        config.ambient_pressure_pa, config.ambient_air_density_kg_m3,
    )


def _relative_change(left: float, right: float, floor: float) -> float:
    return abs(right - left) / max(abs(left), abs(right), floor)


def run_resolved_obstacle_grid_convergence(
    configs: Sequence[ResolvedObstacleFlow2DConfig],
    *,
    receptors: Sequence[FastObstacleReceptor],
    relative_tolerance: float = 0.10,
    concentration_floor_vol_percent: float = 0.01,
    dose_floor_vol_percent_s: float = 0.01,
    arrival_time_tolerance_s: float = 0.10,
    threshold_duration_floor_s: float = 0.10,
) -> ResolvedObstacleGridConvergenceResult:
    """Run at least three matched grids and audit the finest adjacent pair."""
    checked_configs = tuple(configs)
    checked_receptors = tuple(receptors)
    if len(checked_configs) < 3:
        raise ValueError("grid convergence requires at least three configurations")
    if any(not isinstance(item, ResolvedObstacleFlow2DConfig) for item in checked_configs):
        raise TypeError("configs must contain ResolvedObstacleFlow2DConfig values")
    if not checked_receptors:
        raise ValueError("grid convergence requires at least one receptor")
    if any(not isinstance(item, FastObstacleReceptor) for item in checked_receptors):
        raise TypeError("receptors must contain FastObstacleReceptor values")
    labels = tuple(item.label for item in checked_receptors)
    if len(set(labels)) != len(labels):
        raise ValueError("grid convergence receptor labels must be unique")
    for name, value in (
        ("relative_tolerance", relative_tolerance),
        ("concentration_floor_vol_percent", concentration_floor_vol_percent),
        ("dose_floor_vol_percent_s", dose_floor_vol_percent_s),
        ("arrival_time_tolerance_s", arrival_time_tolerance_s),
        ("threshold_duration_floor_s", threshold_duration_floor_s),
    ):
        if isinstance(value, bool) or not math.isfinite(float(value)) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    signature = _scenario_signature(checked_configs[0])
    if any(_scenario_signature(item) != signature for item in checked_configs[1:]):
        raise ValueError("grid convergence configurations must share one physical scenario")
    cell_sizes = tuple(item.length_m / item.nx for item in checked_configs)
    if any(left <= right for left, right in zip(cell_sizes, cell_sizes[1:])):
        raise ValueError("grid convergence configurations must be ordered coarse to fine")
    flow_durations = tuple(
        item.flow_steps * item.lattice_inlet_speed * (item.length_m / item.nx)
        / item.wind_speed_m_s
        for item in checked_configs
    )
    if not all(math.isclose(value, flow_durations[0], rel_tol=1.0e-9, abs_tol=1.0e-12)
               for value in flow_durations[1:]):
        raise ValueError("grid convergence configurations must use one physical flow duration")

    cases: list[ResolvedObstacleGridCase] = []
    for config, cell_size, flow_duration in zip(
        checked_configs, cell_sizes, flow_durations,
    ):
        result = solve_resolved_obstacle_flow_2d(
            config, receptors=checked_receptors,
        )
        conversion = _conversion(config)
        envelope = result.maximum_concentration_kg_m3 * conversion
        above_four = envelope >= 4.0
        if np.any(above_four):
            downwind = float(np.max(result.x_m[np.where(above_four)[1]]))
        else:
            downwind = None
        receptor_metrics = []
        for receptor in checked_receptors:
            values = result.receptor_traces_kg_m3[receptor.label] * conversion
            arrival, duration = _threshold_metrics(result.time_s, values, 4.0)
            receptor_metrics.append(ObstacleReceptorConvergenceMetric(
                label=receptor.label,
                peak_vol_percent=float(np.max(values)),
                dose_vol_percent_s=float(np.trapezoid(values, result.time_s)),
                arrival_time_4vol_s=arrival,
                duration_above_4vol_s=duration,
            ))
        cases.append(ResolvedObstacleGridCase(
            nx=config.nx,
            nz=config.nz,
            cell_size_m=cell_size,
            flow_physical_duration_s=flow_duration,
            runtime_s=float(result.diagnostics["total_runtime_s"]),
            maximum_vol_percent=float(np.nanmax(envelope)),
            four_vol_percent_envelope_area_per_span_m2=float(
                np.count_nonzero(above_four) * cell_size**2
            ),
            four_vol_percent_maximum_downwind_x_m=downwind,
            minimum_lee_streamwise_velocity_m_s=float(
                result.diagnostics["minimum_lee_streamwise_velocity_m_s"]
            ),
            maximum_absolute_divergence_s_1=float(
                result.diagnostics["maximum_absolute_divergence_s_1"]
            ),
            maximum_mass_residual_per_span_kg_m=(
                result.maximum_mass_residual_per_span_kg_m
            ),
            receptor_metrics=tuple(receptor_metrics),
            result=result,
        ))

    coarser, finer = cases[-2], cases[-1]
    coarser_by_label = {item.label: item for item in coarser.receptor_metrics}
    finer_by_label = {item.label: item for item in finer.receptor_metrics}
    receptor_changes = {}
    reasons: list[str] = []
    for label in labels:
        left = coarser_by_label[label]
        right = finer_by_label[label]
        peak_change = _relative_change(
            left.peak_vol_percent, right.peak_vol_percent,
            concentration_floor_vol_percent,
        )
        dose_change = _relative_change(
            left.dose_vol_percent_s, right.dose_vol_percent_s,
            dose_floor_vol_percent_s,
        )
        if left.arrival_time_4vol_s is None and right.arrival_time_4vol_s is None:
            arrival_change: float | str = 0.0
        elif left.arrival_time_4vol_s is None or right.arrival_time_4vol_s is None:
            arrival_change = "threshold_status_changed"
            reasons.append(
                f"{label} 4 vol% arrival status changes on the finest adjacent grids"
            )
        else:
            arrival_change = abs(
                right.arrival_time_4vol_s - left.arrival_time_4vol_s
            )
            if arrival_change > arrival_time_tolerance_s:
                reasons.append(
                    f"{label} 4 vol% arrival changes {arrival_change:.6g} s on the finest adjacent grids"
                )
        duration_change = _relative_change(
            left.duration_above_4vol_s, right.duration_above_4vol_s,
            threshold_duration_floor_s,
        )
        receptor_changes[label] = {
            "peak_vol_percent": peak_change,
            "dose_vol_percent_s": dose_change,
            "arrival_time_4vol_absolute_s": arrival_change,
            "duration_above_4vol_s": duration_change,
        }
        if peak_change > relative_tolerance:
            reasons.append(
                f"{label} peak changes {peak_change:.3%} on the finest adjacent grids"
            )
        if dose_change > relative_tolerance:
            reasons.append(
                f"{label} dose changes {dose_change:.3%} on the finest adjacent grids"
            )
        if duration_change > relative_tolerance:
            reasons.append(
                f"{label} 4 vol% duration changes {duration_change:.3%} on the finest adjacent grids"
            )
    area_change = _relative_change(
        coarser.four_vol_percent_envelope_area_per_span_m2,
        finer.four_vol_percent_envelope_area_per_span_m2,
        min(cell_sizes[-2]**2, cell_sizes[-1]**2),
    )
    if area_change > relative_tolerance:
        reasons.append(
            f"4 vol% envelope area changes {area_change:.3%} on the finest adjacent grids"
        )
    left_distance = coarser.four_vol_percent_maximum_downwind_x_m
    right_distance = finer.four_vol_percent_maximum_downwind_x_m
    if left_distance is None and right_distance is None:
        distance_change: float | str = 0.0
    elif left_distance is None or right_distance is None:
        distance_change = "threshold_status_changed"
        reasons.append(
            "4 vol% maximum downwind extent status changes on the finest adjacent grids"
        )
    else:
        distance_change = _relative_change(left_distance, right_distance, cell_sizes[-1])
        if distance_change > relative_tolerance:
            reasons.append(
                f"4 vol% maximum downwind extent changes {distance_change:.3%} on the finest adjacent grids"
            )
    changes: dict[str, object] = {
        "coarser_cells": [coarser.nx, coarser.nz],
        "finer_cells": [finer.nx, finer.nz],
        "receptors": receptor_changes,
        "four_vol_percent_envelope_area_per_span_m2": area_change,
        "four_vol_percent_maximum_downwind_x_m": distance_change,
        "minimum_lee_streamwise_velocity_m_s": _relative_change(
            coarser.minimum_lee_streamwise_velocity_m_s,
            finer.minimum_lee_streamwise_velocity_m_s,
            1.0e-9,
        ),
    }
    return ResolvedObstacleGridConvergenceResult(
        schema=RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA,
        cases=tuple(cases),
        relative_tolerance=float(relative_tolerance),
        concentration_floor_vol_percent=float(concentration_floor_vol_percent),
        dose_floor_vol_percent_s=float(dose_floor_vol_percent_s),
        arrival_time_tolerance_s=float(arrival_time_tolerance_s),
        threshold_duration_floor_s=float(threshold_duration_floor_s),
        finest_adjacent_relative_changes=changes,
        passed=not reasons,
        gate_reasons=tuple(reasons),
        diagnostics={
            "classification": "numerical-convergence-not-experimental-validation",
            "quantitative_prediction_allowed": False,
            "grid_order": "coarse_to_fine",
            "physical_flow_duration_s": flow_durations[0],
            "warnings": [
                "passing grid convergence does not establish model-form accuracy",
                "Tier-2 thermal, buoyancy, turbulence and three-dimensional gates remain open",
            ],
        },
    )


def resolved_obstacle_grid_convergence_record(
    result: ResolvedObstacleGridConvergenceResult,
) -> dict[str, object]:
    """Return a JSON-native audit record without embedding numerical arrays."""
    if not isinstance(result, ResolvedObstacleGridConvergenceResult):
        raise TypeError("result must be ResolvedObstacleGridConvergenceResult")
    return {
        "schema": result.schema,
        "relative_tolerance": result.relative_tolerance,
        "concentration_floor_vol_percent": result.concentration_floor_vol_percent,
        "dose_floor_vol_percent_s": result.dose_floor_vol_percent_s,
        "arrival_time_tolerance_s": result.arrival_time_tolerance_s,
        "threshold_duration_floor_s": result.threshold_duration_floor_s,
        "passed": result.passed,
        "gate_reasons": list(result.gate_reasons),
        "finest_adjacent_relative_changes": result.finest_adjacent_relative_changes,
        "diagnostics": result.diagnostics,
        "cases": [
            {
                "cells": [case.nx, case.nz],
                "cell_size_m": case.cell_size_m,
                "flow_physical_duration_s": case.flow_physical_duration_s,
                "runtime_s": case.runtime_s,
                "maximum_vol_percent": case.maximum_vol_percent,
                "four_vol_percent_envelope_area_per_span_m2": (
                    case.four_vol_percent_envelope_area_per_span_m2
                ),
                "four_vol_percent_maximum_downwind_x_m": (
                    case.four_vol_percent_maximum_downwind_x_m
                ),
                "minimum_lee_streamwise_velocity_m_s": (
                    case.minimum_lee_streamwise_velocity_m_s
                ),
                "maximum_absolute_divergence_s_1": (
                    case.maximum_absolute_divergence_s_1
                ),
                "maximum_mass_residual_per_span_kg_m": (
                    case.maximum_mass_residual_per_span_kg_m
                ),
                "receptors": [
                    {
                        "label": metric.label,
                        "peak_vol_percent": metric.peak_vol_percent,
                        "dose_vol_percent_s": metric.dose_vol_percent_s,
                        "arrival_time_4vol_s": metric.arrival_time_4vol_s,
                        "duration_above_4vol_s": metric.duration_above_4vol_s,
                    }
                    for metric in case.receptor_metrics
                ],
            }
            for case in result.cases
        ],
    }


__all__ = [
    "RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA",
    "ObstacleReceptorConvergenceMetric",
    "ResolvedObstacleGridCase",
    "ResolvedObstacleGridConvergenceResult",
    "run_resolved_obstacle_grid_convergence",
    "resolved_obstacle_grid_convergence_record",
]
