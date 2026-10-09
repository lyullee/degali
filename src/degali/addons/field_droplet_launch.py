"""Explicit handoff of resolved LH2 droplet-evaporation paths to field transport.

The droplet model retains liquid-to-vapour mass on trajectories.  This module
does not collapse that mass to the nozzle or a rainout centroid: it converts
each recorded, mass-weighted trajectory section into its own atmospheric
scalar source.  It is still a two-dimensional field-screen handoff, so an
off-plane trajectory is refused by :mod:`field_workflow`.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import TYPE_CHECKING

import numpy as np

from .droplet_rainout import (
    DropletPopulationResult,
    DropletTransportInput,
    transport_droplet_population,
)
from .field_distributed_source import FieldDistributedVapourSource
from .semi_fv_obstacle import SourceRateSchedule

if TYPE_CHECKING:
    from .field_workflow import FieldSemiFVRequest, FieldSemiFVScreeningResult


@dataclass(frozen=True)
class DropletVapourLaunchBoundary:
    """Declared scalar closure for trajectory-resolved LH2 droplet vapour."""

    vertical_sigma_m: float
    evidence_id: str
    maximum_sources: int = 512

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.vertical_sigma_m)) or self.vertical_sigma_m <= 0.0:
            raise ValueError("droplet-vapour vertical_sigma_m must be positive and finite")
        if not isinstance(self.evidence_id, str) or not self.evidence_id.strip() or self.evidence_id == "unspecified":
            raise ValueError("droplet-vapour launch requires a declared evidence_id")
        if not isinstance(self.maximum_sources, int) or self.maximum_sources < 1:
            raise ValueError("droplet-vapour maximum_sources must be a positive integer")


def _delayed_constant_schedule(
    *,
    rate_kg_s: float,
    release_duration_s: float,
    age_midpoint_s: float,
    source_id: str,
) -> SourceRateSchedule:
    """Use the segment's midpoint age without moving evaporation to the nozzle."""
    if age_midpoint_s <= 1.0e-12:
        return SourceRateSchedule(
            (0.0, release_duration_s), (rate_kg_s, 0.0), source_id=source_id,
        )
    return SourceRateSchedule(
        (0.0, age_midpoint_s, age_midpoint_s + release_duration_s),
        (0.0, rate_kg_s, 0.0), source_id=source_id,
    )


def distributed_sources_from_droplet_evaporation(
    droplets: DropletPopulationResult,
    *,
    release_duration_s: float,
    launch: DropletVapourLaunchBoundary,
    label_prefix: str = "droplet-vapour",
) -> tuple[FieldDistributedVapourSource, ...]:
    """Map recorded droplet evaporation sections to delayed field sources.

    The mass for each source is exact for the recorded section.  Its location
    and onset age are the section's vapour-mass-weighted midpoint; reducing
    ``trajectory_segment_duration_s`` in the droplet calculation refines
    this representation.  It does not map unvaporised airborne liquid.
    """
    if not isinstance(droplets, DropletPopulationResult):
        raise TypeError("droplets must be a DropletPopulationResult")
    if not math.isfinite(float(release_duration_s)) or release_duration_s <= 0.0:
        raise ValueError("release_duration_s must be positive and finite")
    if not isinstance(launch, DropletVapourLaunchBoundary):
        raise TypeError("launch must be a DropletVapourLaunchBoundary")
    if not isinstance(label_prefix, str) or not label_prefix.strip():
        raise ValueError("label_prefix must be non-empty")

    pieces = [
        (class_index, outcome, segment)
        for class_index, outcome in enumerate(droplets.outcomes)
        for segment in outcome.evaporation_segments
    ]
    expected_rate = droplets.airborne_vapour_mass_flow_kg_s
    if expected_rate > 0.0 and not pieces:
        raise ValueError(
            "in-flight droplet vapour has no recorded trajectory segments; "
            "rerun droplet transport with trajectory_segment_duration_s"
        )
    if len(pieces) > launch.maximum_sources:
        raise ValueError(
            f"{len(pieces)} droplet-vapour sources exceed maximum_sources={launch.maximum_sources}; "
            "use a coarser recorded trajectory segment duration or raise the declared limit"
        )

    sources: list[FieldDistributedVapourSource] = []
    represented_rate = 0.0
    for class_index, outcome, segment in pieces:
        rate = outcome.inlet_mass_flow_kg_s * segment.vapour_mass_fraction
        if rate <= 0.0:
            continue
        age_midpoint_s = 0.5 * (segment.start_time_s + segment.end_time_s)
        label = f"{label_prefix}-{class_index:02d}-{len(sources):04d}"
        schedule = _delayed_constant_schedule(
            rate_kg_s=rate,
            release_duration_s=release_duration_s,
            age_midpoint_s=age_midpoint_s,
            source_id=f"{label}:trajectory-segment",
        )
        sources.append(FieldDistributedVapourSource(
            label=label,
            position_m=segment.vapour_centroid_m,
            vertical_sigma_m=launch.vertical_sigma_m,
            schedule=schedule,
            evidence_id=launch.evidence_id,
            source_kind="in_flight_droplet_evaporation",
            warnings=(
                "trajectory-section mass is conserved at its vapour-mass-weighted spatial and temporal midpoint",
                "declared d-squared evaporation and mean-wind droplet trajectory are used; turbulent droplet dispersion, cold-cloud thermodynamics and lateral dilution remain unresolved",
            ) + tuple(droplets.qualifications),
        ))
        represented_rate += rate
    if not math.isclose(represented_rate, expected_rate, rel_tol=1.0e-9, abs_tol=1.0e-12):
        raise ValueError("recorded droplet-evaporation source rates do not close the vapour ledger")
    return tuple(sources)


@dataclass(frozen=True)
class FieldDropletHandoffRefinementCase:
    """One trajectory-section duration evaluated at the declared field sensor."""

    trajectory_segment_duration_s: float
    source_count: int
    in_flight_vapour_mass_kg: float
    sensor_true_peak_mole_fraction: float
    sensor_true_dose_mole_fraction_s: float
    screening: "FieldSemiFVScreeningResult"

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.trajectory_segment_duration_s)) or self.trajectory_segment_duration_s <= 0.0:
            raise ValueError("droplet refinement trajectory duration must be positive and finite")
        if not isinstance(self.source_count, int) or isinstance(self.source_count, bool) or self.source_count < 1:
            raise ValueError("droplet refinement source_count must be a positive integer")
        for name, value in {
            "in_flight_vapour_mass_kg": self.in_flight_vapour_mass_kg,
            "sensor_true_peak_mole_fraction": self.sensor_true_peak_mole_fraction,
            "sensor_true_dose_mole_fraction_s": self.sensor_true_dose_mole_fraction_s,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)) or value < 0.0:
                raise ValueError(f"droplet refinement {name} must be finite and non-negative")
        from .field_workflow import FieldSemiFVScreeningResult

        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError("droplet refinement screening must be a FieldSemiFVScreeningResult")


@dataclass(frozen=True)
class FieldDropletHandoffRefinementStudy:
    """Field response to making the droplet-evaporation handoff finer."""

    droplet_boundary: DropletTransportInput
    launch: DropletVapourLaunchBoundary
    release_duration_s: float
    maximum_droplet_time_s: float
    cases: tuple[FieldDropletHandoffRefinementCase, ...]
    reference_trajectory_segment_duration_s: float
    peak_relative_changes: tuple[tuple[float, float], ...]
    dose_relative_changes: tuple[tuple[float, float], ...]
    relative_tolerance: float
    converged: bool
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.droplet_boundary, DropletTransportInput):
            raise TypeError("droplet refinement droplet_boundary must be DropletTransportInput")
        if not isinstance(self.launch, DropletVapourLaunchBoundary):
            raise TypeError("droplet refinement launch must be DropletVapourLaunchBoundary")
        for name, value in {
            "release_duration_s": self.release_duration_s,
            "maximum_droplet_time_s": self.maximum_droplet_time_s,
            "relative_tolerance": self.relative_tolerance,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)) or value <= 0.0:
                raise ValueError(f"droplet refinement {name} must be positive and finite")
        if not isinstance(self.cases, tuple) or len(self.cases) < 2:
            raise ValueError("droplet refinement study requires at least two cases")
        if any(not isinstance(case, FieldDropletHandoffRefinementCase) for case in self.cases):
            raise TypeError(
                "droplet refinement cases must contain only "
                "FieldDropletHandoffRefinementCase values"
            )
        durations = tuple(case.trajectory_segment_duration_s for case in self.cases)
        if len(set(durations)) != len(durations):
            raise ValueError("droplet refinement case durations must be unique")
        if not math.isfinite(float(self.reference_trajectory_segment_duration_s)) or self.reference_trajectory_segment_duration_s <= 0.0:
            raise ValueError("droplet refinement reference duration must be positive and finite")
        if self.reference_trajectory_segment_duration_s not in durations:
            raise ValueError("droplet refinement reference duration must identify a case")
        for name, values in {
            "peak_relative_changes": self.peak_relative_changes,
            "dose_relative_changes": self.dose_relative_changes,
        }.items():
            if not isinstance(values, tuple) or len(values) != len(self.cases):
                raise ValueError(f"droplet refinement {name} must match the case count")
            for duration, change in values:
                if not math.isfinite(float(duration)) or duration not in durations:
                    raise ValueError(f"droplet refinement {name} has an unknown duration")
                if not math.isfinite(float(change)) or change < 0.0:
                    raise ValueError(f"droplet refinement {name} changes must be finite and non-negative")
        if not isinstance(self.converged, bool):
            raise TypeError("droplet refinement converged must be boolean")
        if not isinstance(self.warnings, tuple) or any(
            not isinstance(item, str) or not item.strip() for item in self.warnings
        ):
            raise TypeError("droplet refinement warnings must be non-empty strings")


def _relative_change(value: float, reference: float, floor: float) -> float:
    return abs(value - reference) / max(abs(reference), floor)


def run_field_droplet_handoff_refinement_study(
    request: "FieldSemiFVRequest",
    droplet_boundary: DropletTransportInput,
    *,
    release_duration_s: float,
    launch: DropletVapourLaunchBoundary,
    trajectory_segment_durations_s: tuple[float, ...] = (0.2, 0.1, 0.05),
    maximum_droplet_time_s: float = 120.0,
    relative_tolerance: float = 0.05,
    absolute_floor_mole_fraction: float = 1.0e-12,
    maximum_cell_source_steps: int = 50_000_000,
) -> FieldDropletHandoffRefinementStudy:
    """Audit field-sensor sensitivity to the droplet handoff discretisation.

    The study is intentionally restricted to one declared transient field
    scenario and its in-plane sensor. It holds the droplet boundary, source
    duration, transport mesh and all direct-flash inputs fixed while only the
    trajectory-recorded evaporation section duration changes. Existing
    distributed sources are refused so a reported change cannot be attributed
    to an unrelated pool or secondary release.
    """
    from .field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(droplet_boundary, DropletTransportInput):
        raise TypeError("droplet_boundary must be a DropletTransportInput")
    if not isinstance(launch, DropletVapourLaunchBoundary):
        raise TypeError("launch must be a DropletVapourLaunchBoundary")
    durations = tuple(float(value) for value in trajectory_segment_durations_s)
    if len(durations) < 2 or len(set(durations)) != len(durations):
        raise ValueError("trajectory_segment_durations_s needs at least two unique values")
    if any(not math.isfinite(value) or value <= 0.0 for value in durations):
        raise ValueError("trajectory segment durations must be positive and finite")
    if not math.isfinite(float(release_duration_s)) or release_duration_s <= 0.0:
        raise ValueError("release_duration_s must be positive and finite")
    if not math.isfinite(float(maximum_droplet_time_s)) or maximum_droplet_time_s <= 0.0:
        raise ValueError("maximum_droplet_time_s must be positive and finite")
    if not math.isfinite(float(relative_tolerance)) or relative_tolerance < 0.0:
        raise ValueError("relative_tolerance must be non-negative and finite")
    if not math.isfinite(float(absolute_floor_mole_fraction)) or absolute_floor_mole_fraction <= 0.0:
        raise ValueError("absolute_floor_mole_fraction must be positive and finite")
    if not isinstance(maximum_cell_source_steps, int) or maximum_cell_source_steps < 1:
        raise ValueError("maximum_cell_source_steps must be a positive integer")
    scenario = request.scenario
    if scenario.sensor is None:
        raise ValueError("droplet handoff refinement requires an in-plane declared field sensor")
    if scenario.temporal_mode != "transient" or scenario.source.duration_s is None:
        raise ValueError("droplet handoff refinement requires a transient finite field release")
    if not math.isclose(release_duration_s, scenario.source.duration_s, rel_tol=1.0e-9, abs_tol=1.0e-12):
        raise ValueError("release_duration_s must match the field ReleaseSource duration")
    if request.direct_vapour_schedule is not None:
        raise ValueError("droplet handoff refinement requires the nominal direct-flash branch, not a replacement schedule")
    if request.distributed_vapour_sources:
        raise ValueError("droplet handoff refinement requires no pre-existing distributed vapour sources")
    source_location = np.asarray(scenario.source.location_m, dtype=float)
    droplet_location = np.asarray(droplet_boundary.release_position_m, dtype=float)
    if not np.allclose(source_location, droplet_location, rtol=0.0, atol=1.0e-9):
        raise ValueError("droplet boundary release_position_m must match the field ReleaseSource location")

    cases = []
    for duration in sorted(durations, reverse=True):
        droplets = transport_droplet_population(
            droplet_boundary,
            max_time_s=maximum_droplet_time_s,
            trajectory_segment_duration_s=duration,
        )
        sources = distributed_sources_from_droplet_evaporation(
            droplets, release_duration_s=release_duration_s, launch=launch,
            label_prefix=f"droplet-vapour-dt-{duration:.6g}",
        )
        estimated_work = (
            len(sources) * request.transport.nx * request.transport.nz
            * math.ceil((release_duration_s + request.post_release_duration_s) / request.transport.time_step_s)
        )
        if estimated_work > maximum_cell_source_steps:
            raise ValueError(
                "droplet handoff refinement exceeds maximum_cell_source_steps="
                f"{maximum_cell_source_steps} at trajectory interval {duration:g} s"
            )
        try:
            screening = run_field_semi_fv_screening(replace(
                request, distributed_vapour_sources=sources,
            ))
        except ValueError as error:
            raise ValueError(
                f"droplet handoff refinement cannot construct interval {duration:g} s: {error}"
            ) from error
        if not screening.completed or screening.sensor_trace is None:
            detail = "; ".join(screening.applicability.reasons)
            raise ValueError(
                f"droplet handoff refinement interval {duration:g} s did not produce a sensor field result: {detail}"
            )
        trace = screening.sensor_trace
        cases.append(FieldDropletHandoffRefinementCase(
            trajectory_segment_duration_s=duration,
            source_count=len(sources),
            in_flight_vapour_mass_kg=sum(source.schedule.released_mass_kg for source in sources),
            sensor_true_peak_mole_fraction=float(np.max(trace.true_mole_fraction)),
            sensor_true_dose_mole_fraction_s=float(np.trapezoid(trace.true_mole_fraction, trace.time_s)),
            screening=screening,
        ))
    ordered = tuple(sorted(cases, key=lambda case: case.trajectory_segment_duration_s, reverse=True))
    reference = ordered[-1]
    peak_changes = tuple((
        case.trajectory_segment_duration_s,
        _relative_change(
            case.sensor_true_peak_mole_fraction,
            reference.sensor_true_peak_mole_fraction,
            absolute_floor_mole_fraction,
        ),
    ) for case in ordered)
    dose_changes = tuple((
        case.trajectory_segment_duration_s,
        _relative_change(
            case.sensor_true_dose_mole_fraction_s,
            reference.sensor_true_dose_mole_fraction_s,
            absolute_floor_mole_fraction,
        ),
    ) for case in ordered)
    converged = all(
        change <= relative_tolerance
        for _duration, change in (*peak_changes, *dose_changes)
    )
    return FieldDropletHandoffRefinementStudy(
        droplet_boundary=droplet_boundary,
        launch=launch,
        release_duration_s=release_duration_s,
        maximum_droplet_time_s=maximum_droplet_time_s,
        cases=ordered,
        reference_trajectory_segment_duration_s=reference.trajectory_segment_duration_s,
        peak_relative_changes=peak_changes,
        dose_relative_changes=dose_changes,
        relative_tolerance=relative_tolerance,
        converged=converged,
        warnings=(
            "this is a trajectory-to-scalar-source discretisation study; it does not validate the d-squared evaporation coefficient, turbulence closure or off-plane droplet transport",
            "all cases retain the same field mesh, direct-flash source, weather and sensor; only the recorded droplet-evaporation section duration changes",
        ),
    )


__all__ = [
    "DropletVapourLaunchBoundary",
    "distributed_sources_from_droplet_evaporation",
    "FieldDropletHandoffRefinementCase", "FieldDropletHandoffRefinementStudy",
    "run_field_droplet_handoff_refinement_study",
]
