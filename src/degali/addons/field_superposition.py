"""Failure-safe sensor superposition for separately screened field sources.

The local semi-FV scalar equation is linear only for a common declared wind
and transport boundary.  Direct flash vapour, jet-handoff vapour and a
separately routed pool-vapour schedule may therefore be added at a *shared
sensor* after their true H2 mass-concentration traces have been established.
Instrument gain, bias, response and averaging are then applied once to the
sum.  This module deliberately refuses convenient but non-equivalent cases.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from .field_contracts import FieldApplicability, SensorModel
from .field_observation import (
    FieldSensorTrace,
    apply_sensor_model,
    h2_mole_fraction_from_mass_concentration,
)
from .field_workflow import FieldSemiFVScreeningResult


@dataclass(frozen=True)
class FieldSensorBranch:
    """One named completed local-source calculation for sensor superposition."""

    label: str
    result: FieldSemiFVScreeningResult

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("field sensor branch label must be non-empty")
        if not isinstance(self.result, FieldSemiFVScreeningResult):
            raise TypeError("field sensor branch result must be a FieldSemiFVScreeningResult")


@dataclass(frozen=True)
class FieldSensorSuperposition:
    """One common-sensor sum, with raw scalar and indicated trace retained."""

    branches: tuple[FieldSensorBranch, ...]
    applicability: FieldApplicability
    time_s: np.ndarray
    concentration_kg_m3: np.ndarray
    sensor_trace: FieldSensorTrace
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if len(self.branches) < 2 or any(
            not isinstance(branch, FieldSensorBranch) for branch in self.branches
        ):
            raise ValueError("sensor superposition requires at least two valid branches")
        labels = tuple(branch.label for branch in self.branches)
        if len(set(labels)) != len(labels):
            raise ValueError("field sensor superposition branch labels must be unique")
        if not isinstance(self.applicability, FieldApplicability):
            raise TypeError("sensor superposition applicability must be FieldApplicability")
        if not isinstance(self.sensor_trace, FieldSensorTrace):
            raise TypeError("sensor superposition sensor_trace must be FieldSensorTrace")
        arrays = {
            "time_s": self.time_s,
            "concentration_kg_m3": self.concentration_kg_m3,
        }
        for name, value in arrays.items():
            array = np.asarray(value, dtype=float)
            if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
                raise ValueError(
                    f"sensor superposition {name} must be a finite non-empty one-dimensional array"
                )
        time = np.asarray(self.time_s, dtype=float)
        concentration = np.asarray(self.concentration_kg_m3, dtype=float)
        if time.size != concentration.size or time.size != self.sensor_trace.time_s.size:
            raise ValueError("sensor superposition arrays must have equal lengths")
        if time.size > 1 and np.any(np.diff(time) <= 0.0):
            raise ValueError("sensor superposition time_s must be strictly increasing")
        if np.any(concentration < 0.0):
            raise ValueError("sensor superposition concentration cannot be negative")
        if not np.array_equal(time, self.sensor_trace.time_s):
            raise ValueError("sensor superposition time_s must match the sensor trace")
        expected_true = h2_mole_fraction_from_mass_concentration(
            concentration,
            ambient_air_density_kg_m3=self.sensor_trace.ambient_air_density_kg_m3,
        )
        if not np.allclose(
            expected_true, self.sensor_trace.true_mole_fraction,
            rtol=1.0e-9, atol=1.0e-12,
        ):
            raise ValueError("sensor superposition concentration does not match the true sensor trace")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("sensor superposition warnings must be non-empty strings")

    @property
    def completed(self) -> bool:
        return self.applicability.status != "blocked"


def _common_sensor(branch: FieldSensorBranch) -> SensorModel:
    result = branch.result
    if not result.completed or result.transport is None or result.sensor_trace is None:
        raise ValueError(
            f"field branch {branch.label!r} must be completed with an in-plane sensor trace"
        )
    if len(result.transport.receptor_traces) != 1:
        raise ValueError(
            f"field branch {branch.label!r} must expose exactly one declared sensor receptor"
        )
    return result.sensor_trace.sensor


def _same_number(left: float, right: float, name: str) -> None:
    if not math.isclose(float(left), float(right), rel_tol=1.0e-12, abs_tol=1.0e-12):
        raise ValueError(f"field branches do not share {name}")


def _unresolved_request_uncertainty(request) -> tuple[str, ...]:
    """Return deterministic field bounds not resolved by this nominal branch."""
    labels: list[str] = [
        name for name, value in request.scenario.uncertainty_fields().items()
        if not value.is_exact
    ]
    labels.extend(
        name for name, value in request.ambient_uncertainty_fields().items()
        if not value.is_exact
    )
    for uncertainty in request.obstacle_geometry_uncertainty:
        labels.extend(
            f"obstacle.{uncertainty.obstacle.label}.{name}"
            for name, value in uncertainty.uncertainty_fields().items()
            if not value.is_exact
        )
    for source in request.distributed_vapour_sources:
        labels.extend(
            name for name, value in source.uncertainty_fields().items()
            if not value.is_exact
        )
        if source.has_schedule_uncertainty:
            labels.append(f"distributed_source.{source.label}.schedule")
    if request.stability_alternatives is not None and not request.stability_uncertainty_resolved:
        labels.append("weather.stability")
    if (
        any(source.source_kind == "pool_vapour" for source in request.distributed_vapour_sources)
        and not request.phase_routing_uncertainty_resolved
    ):
        labels.append("phase_routing")
    if (
        request.measured_history_quality is not None
        and not request.measured_history_source_uncertainty_resolved
    ):
        labels.append("measured_history.source")
    return tuple(dict.fromkeys(labels))


def _validate_common_boundary(
    branches: tuple[FieldSensorBranch, ...],
) -> tuple[SensorModel, float, float]:
    first = branches[0].result
    sensor = _common_sensor(branches[0])
    assert first.transport is not None
    base_trace = first.transport.receptor_traces[0]
    base_request = first.request
    start, end = float(base_trace.time_s[0]), float(base_trace.time_s[-1])
    for branch in branches[1:]:
        result = branch.result
        candidate_sensor = _common_sensor(branch)
        assert result.transport is not None
        trace = result.transport.receptor_traces[0]
        if candidate_sensor != sensor:
            raise ValueError("field branches must use one identical sensor boundary")
        if result.request.scenario.weather != base_request.scenario.weather:
            raise ValueError("field branches must use one identical weather boundary")
        if result.request.scenario.temporal_mode != base_request.scenario.temporal_mode:
            raise ValueError("field branches must use one identical temporal mode")
        if result.request.obstacle != base_request.obstacle or result.request.obstacles != base_request.obstacles:
            raise ValueError("field branches must use one identical global obstacle boundary")
        if result.request.stability_mixing_closure != base_request.stability_mixing_closure:
            raise ValueError("field branches must use one identical stability mixing closure")
        _same_number(
            result.request.ambient_air_density_kg_m3,
            base_request.ambient_air_density_kg_m3,
            "ambient_air_density_kg_m3",
        )
        _same_number(result.transport.diagnostics.scalar_diffusivity_m2_s,
                     first.transport.diagnostics.scalar_diffusivity_m2_s,
                     "resolved scalar diffusivity")
        _same_number(float(trace.time_s[0]), start, "trace start time")
        _same_number(float(trace.time_s[-1]), end, "trace end time")
        # Summation is only physically reproducible when every branch was
        # advanced on the same observation grid.  Matching endpoints alone is
        # insufficient: np.interp would otherwise silently insert a different
        # temporal operator (and potentially hide a source shutoff or peak)
        # into one branch before the sensor operator is applied once to the
        # total.
        if not np.array_equal(np.asarray(trace.time_s), np.asarray(base_trace.time_s)):
            raise ValueError(
                "field branches must use one identical sensor time grid"
            )
    return sensor, start, end


def superpose_field_sensor_branches(
    branches: Sequence[FieldSensorBranch],
) -> FieldSensorSuperposition:
    """Sum true scalar branch traces and apply one common sensor operator.

    Every branch must have the same fixed weather, obstacle, scalar-diffusion,
    sensor and complete time coverage.  A branch's indicated reading is never
    summed because gain/bias and first-order response must apply to the total
    H2 signal once.
    """
    declared = tuple(branches)
    if len(declared) < 2:
        raise ValueError("sensor superposition requires at least two field branches")
    if any(not isinstance(branch, FieldSensorBranch) for branch in declared):
        raise TypeError("every item must be a FieldSensorBranch")
    labels = tuple(branch.label for branch in declared)
    if len(set(labels)) != len(labels):
        raise ValueError("field sensor branch labels must be unique")
    sensor, _start, _end = _validate_common_boundary(declared)
    first = declared[0].result
    assert first.transport is not None
    time = np.unique(np.concatenate([
        branch.result.transport.receptor_traces[0].time_s
        for branch in declared
    ]))
    concentration = np.zeros_like(time, dtype=float)
    applicability_warnings: list[str] = []
    uncertainty_complete = True
    for branch in declared:
        result = branch.result
        assert result.transport is not None
        trace = result.transport.receptor_traces[0]
        concentration += np.interp(time, trace.time_s, trace.concentration_kg_m3)
        applicability_warnings.extend(
            f"{branch.label}: {warning}" for warning in result.applicability.warnings
        )
        unresolved = _unresolved_request_uncertainty(result.request)
        if unresolved:
            applicability_warnings.append(
                f"{branch.label}: nominal superposition leaves unresolved field uncertainty "
                + ", ".join(unresolved)
                + "; run the corresponding deterministic uncertainty envelope"
            )
        uncertainty_complete = (
            uncertainty_complete
            and result.applicability.uncertainty_complete
            and not unresolved
        )
    sensor_trace = apply_sensor_model(
        time,
        concentration,
        sensor,
        ambient_air_density_kg_m3=first.request.ambient_air_density_kg_m3,
    )
    warnings = (
        "separate field-source scalar traces are superposed only at the shared sensor; "
        "source-source interaction, three-dimensional lateral mixing and nonlinear cold-cloud thermodynamics are unresolved",
    ) + tuple(applicability_warnings)
    applicability = FieldApplicability(
        "conditional",
        warnings=warnings,
        uncertainty_complete=uncertainty_complete,
    )
    return FieldSensorSuperposition(
        declared, applicability, time, concentration, sensor_trace, warnings,
    )


__all__ = [
    "FieldSensorBranch", "FieldSensorSuperposition",
    "superpose_field_sensor_branches",
]
