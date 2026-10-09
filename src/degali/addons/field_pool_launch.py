"""Explicit handoff of a conservative LH2 pool ledger to local scalar transport.

Pool evaporation mass and timing can be resolved by the dynamic-pool ledger
without defining an atmospheric launch condition. This adapter requires the
caller to choose a ground-level scalar closure and records that choice instead
of merging pool vapour into a nozzle jet.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import TYPE_CHECKING

from .semi_fv_obstacle import SourceRateSchedule
from .field_distributed_source import FieldDistributedVapourSource

if TYPE_CHECKING:
    from .field_phase_routing import FieldPhaseRoutingResult
    from .field_workflow import FieldSemiFVRequest


@dataclass(frozen=True)
class PoolVapourLaunchBoundary:
    """Declared atmospheric handoff for mass leaving a liquid-hydrogen pool.

    The only currently supported closure starts the scalar source at a
    specified ground-level height. It does not resolve pool-vapour temperature,
    vertical momentum, lateral footprint dilution or a 3-D obstacle wake.
    """

    source_height_m: float = 0.0
    closure_id: str = "ground_level_scalar"
    evidence_id: str = "unspecified"

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.source_height_m)) or self.source_height_m < 0.0:
            raise ValueError("pool-vapour source_height_m must be finite and non-negative")
        if self.closure_id != "ground_level_scalar":
            raise ValueError("only the explicit ground_level_scalar pool-vapour closure is supported")
        if not isinstance(self.evidence_id, str) or not self.evidence_id.strip() or self.evidence_id == "unspecified":
            raise ValueError("pool-vapour launch requires a declared evidence_id")


@dataclass(frozen=True)
class FieldPoolVapourSchedule:
    """Conservative pool evaporation history at an explicit impact location."""

    schedule: SourceRateSchedule
    position_m: tuple[float, float, float]
    maximum_wet_area_m2: float | None
    evaporated_mass_kg: float
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.schedule, SourceRateSchedule):
            raise TypeError("pool-vapour schedule must be a SourceRateSchedule")
        self.schedule.require_zero_endpoint("pool-vapour schedule")
        if (
            not isinstance(self.position_m, tuple)
            or len(self.position_m) != 3
            or any(isinstance(value, bool) for value in self.position_m)
            or not all(math.isfinite(float(value)) for value in self.position_m)
        ):
            raise ValueError("pool-vapour position_m must contain three finite values")
        if self.position_m[2] != 0.0:
            raise ValueError("pool-vapour impact position must be at ground height")
        if not math.isfinite(float(self.evaporated_mass_kg)) or self.evaporated_mass_kg < 0.0:
            raise ValueError("evaporated_mass_kg must be finite and non-negative")
        if self.maximum_wet_area_m2 is not None and (
            not math.isfinite(float(self.maximum_wet_area_m2)) or self.maximum_wet_area_m2 <= 0.0
        ):
            raise ValueError("maximum_wet_area_m2 must be positive and finite when supplied")
        if not math.isclose(self.schedule.released_mass_kg, self.evaporated_mass_kg, rel_tol=1.0e-8, abs_tol=1.0e-10):
            raise ValueError("pool-vapour schedule mass must match its evaporation ledger")
        if not isinstance(self.warnings, tuple):
            raise TypeError("pool-vapour warnings must be a tuple")
        if any(not isinstance(item, str) or not item.strip() for item in self.warnings):
            raise ValueError("pool-vapour warnings must contain non-empty strings")


def pool_vapour_schedule_from_phase_routing(
    result: "FieldPhaseRoutingResult",
    *,
    source_id: str = "field-pool:evaporation-vapour",
) -> FieldPoolVapourSchedule:
    """Convert an existing conservative dynamic-pool history into source bins.

    The extraction requires a resolved impact centroid and dynamic-pool step
    ledger. Fixed-footprint outputs without timed steps are withheld rather
    than treated as a transient source history.
    """
    from .field_phase_routing import FieldPhaseRoutingResult

    if not isinstance(result, FieldPhaseRoutingResult):
        raise TypeError("result must be a FieldPhaseRoutingResult")
    if result.coupled_result is None:
        raise ValueError("field phase routing did not produce a coupled pool result")
    coupled = result.coupled_result
    if not getattr(coupled, "conservative", False):
        raise ValueError("pool-vapour handoff requires a conservative coupled mass ledger")
    phase = getattr(coupled, "phase_routing", None)
    if phase is None:
        raise ValueError("coupled result contains no rainout/pool routing branch")
    coupling = phase.pool_coupling
    pool = coupling.pool
    if pool is None or not getattr(pool, "steps", ()):
        raise ValueError("pool-vapour handoff requires a time-resolved dynamic-pool ledger")
    centroid = phase.droplets.impact_centroid_m
    if centroid is None:
        raise ValueError("pool-vapour handoff requires a resolved rainout impact centroid")
    if not isinstance(source_id, str) or not source_id.strip():
        raise ValueError("source_id must be non-empty")

    times = [0.0]
    rates = []
    previous_time = 0.0
    previous_mass = 0.0
    wet_areas = []
    for step in pool.steps:
        elapsed = float(step.elapsed_s)
        cumulative = float(step.evaporated_mass_kg)
        if elapsed < previous_time:
            raise ValueError("dynamic-pool step times are not ordered")
        if math.isclose(elapsed, previous_time, abs_tol=1.0e-12):
            previous_mass = cumulative
            continue
        evaporated = cumulative - previous_mass
        if evaporated < -1.0e-10:
            raise ValueError("dynamic-pool cumulative evaporation decreases")
        rates.append(max(0.0, evaporated) / (elapsed - previous_time))
        times.append(elapsed)
        previous_time, previous_mass = elapsed, cumulative
        area = getattr(step, "wet_area_m2", None)
        if area is not None and math.isfinite(float(area)) and area > 0.0:
            wet_areas.append(float(area))
    if len(times) < 2:
        raise ValueError("dynamic-pool ledger contains no positive-duration source interval")
    schedule = SourceRateSchedule(tuple(times), tuple(rates) + (0.0,), source_id=source_id)
    reported_mass = float(coupling.pool_vapour_mass_kg)
    position = (float(centroid[0]), float(centroid[1]), 0.0)
    warnings = (
        "pool-vapour schedule is extracted from the conservative pool evaporation ledger only",
        "in-flight droplet vapour and direct flash vapour remain separate atmospheric source branches",
    )
    return FieldPoolVapourSchedule(
        schedule=schedule,
        position_m=position,
        maximum_wet_area_m2=max(wet_areas) if wet_areas else None,
        evaporated_mass_kg=reported_mass,
        warnings=warnings,
    )


def request_with_pool_vapour_schedule(
    request: "FieldSemiFVRequest",
    pool_schedule: FieldPoolVapourSchedule,
    launch: PoolVapourLaunchBoundary,
) -> "FieldSemiFVRequest":
    """Create a separate local scalar screen for the declared pool-vapour branch."""
    from .field_workflow import FieldSemiFVRequest

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(pool_schedule, FieldPoolVapourSchedule):
        raise TypeError("pool_schedule must be a FieldPoolVapourSchedule")
    if not isinstance(launch, PoolVapourLaunchBoundary):
        raise TypeError("launch must be a PoolVapourLaunchBoundary")
    if request.direct_vapour_schedule is not None:
        raise ValueError("request already declares another atmospheric direct-vapour schedule")
    source = replace(
        request.scenario.source,
        location_m=(
            pool_schedule.position_m[0], pool_schedule.position_m[1], launch.source_height_m,
        ),
        duration_s=pool_schedule.schedule.duration_s,
        duration_uncertainty=None,
    )
    scenario = replace(request.scenario, source=source, temporal_mode="transient")
    warnings = (
        "pool evaporation is injected through the declared ground_level_scalar closure; "
        "pool-vapour temperature, vertical momentum and three-dimensional footprint dilution are unresolved",
        f"pool-vapour launch evidence_id={launch.evidence_id!r}",
    ) + pool_schedule.warnings
    return replace(
        request,
        scenario=scenario,
        direct_vapour_schedule=pool_schedule.schedule,
        direct_vapour_warnings=tuple(request.direct_vapour_warnings) + warnings,
        direct_vapour_effective_area_m2=pool_schedule.maximum_wet_area_m2,
    )


def distributed_source_from_pool_vapour_schedule(
    pool_schedule: FieldPoolVapourSchedule,
    launch: PoolVapourLaunchBoundary,
    *,
    vertical_sigma_m: float,
    label: str = "pool-vapour",
) -> FieldDistributedVapourSource:
    """Create an internal field source without replacing direct flash vapour.

    This adapter preserves the dynamic-pool mass/timing ledger and its impact
    location. It does not claim that pool vapour is a nozzle jet: the caller
    must still supply the vertical scalar width and the field workflow will
    reject the source if its global position is outside the declared wind
    plane.
    """
    if not isinstance(pool_schedule, FieldPoolVapourSchedule):
        raise TypeError("pool_schedule must be a FieldPoolVapourSchedule")
    if not isinstance(launch, PoolVapourLaunchBoundary):
        raise TypeError("launch must be a PoolVapourLaunchBoundary")
    return FieldDistributedVapourSource(
        label=label,
        position_m=(
            pool_schedule.position_m[0], pool_schedule.position_m[1],
            launch.source_height_m,
        ),
        vertical_sigma_m=vertical_sigma_m,
        schedule=pool_schedule.schedule,
        evidence_id=launch.evidence_id,
        source_kind="pool_vapour",
        warnings=(
            "pool evaporation is injected as a declared internal scalar source; "
            "pool-vapour temperature, vertical momentum and three-dimensional footprint dilution are unresolved",
        ) + pool_schedule.warnings,
    )


__all__ = [
    "PoolVapourLaunchBoundary", "FieldPoolVapourSchedule",
    "pool_vapour_schedule_from_phase_routing", "request_with_pool_vapour_schedule",
    "distributed_source_from_pool_vapour_schedule",
]
