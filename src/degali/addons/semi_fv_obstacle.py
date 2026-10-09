"""Local conservative semi-finite-volume obstacle transport.

This is a reduced-order local transport operator, not a CFD solver.  It
advects a scalar H2 inventory on an x-z finite-volume mesh and routes blocked
downwind flux through available cells above and below a rectangular obstacle.
The routing is conservative by construction and exposes its residual and
applicability status.  It is intended for screening and sensitivity studies
until an obstacle-specific LH2 validation set is available.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class RectangularObstacle2D:
    """A wind-aligned obstacle in local downwind/height coordinates."""

    distance_m: float
    width_m: float
    height_m: float
    base_height_m: float = 0.0
    label: str = "obstacle"

    def __post_init__(self) -> None:
        for name, value in (("distance_m", self.distance_m), ("width_m", self.width_m),
                            ("height_m", self.height_m), ("base_height_m", self.base_height_m)):
            if isinstance(value, bool):
                raise ValueError(f"{name} must be finite and not boolean")
            if not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
        if self.distance_m < 0.0 or self.width_m <= 0.0 or self.height_m <= 0.0:
            raise ValueError("obstacle distance must be non-negative and size positive")
        if self.base_height_m < 0.0:
            raise ValueError("base_height_m cannot be negative")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("label must be non-empty")

    @property
    def x_min_m(self) -> float:
        return self.distance_m - self.width_m / 2.0

    @property
    def x_max_m(self) -> float:
        return self.distance_m + self.width_m / 2.0

    @property
    def z_max_m(self) -> float:
        return self.base_height_m + self.height_m


@dataclass(frozen=True)
class SourceRateSchedule:
    """Atmospheric source rate over an explicit duration.

    ``rate_operator='piecewise_constant'`` retains the legacy interval-start
    interpretation: ``rate_kg_s[i]`` applies from ``time_s[i]`` up to
    ``time_s[i + 1]``.  ``rate_operator='linear'`` treats the rates as values
    at the declared time nodes and integrates the linear segment between them.
    In both modes the final rate is retained as an endpoint record and is not
    integrated beyond ``time_s[-1]``. This class carries an already-
    atmospheric scalar source only; it does not turn an upstream liquid or
    two-phase historian record into a flash calculation.
    """

    time_s: tuple[float, ...]
    rate_kg_s: tuple[float, ...]
    source_id: str = "declared"
    rate_operator: str = "piecewise_constant"

    def __post_init__(self) -> None:
        if len(self.time_s) < 2 or len(self.time_s) != len(self.rate_kg_s):
            raise ValueError("source schedule requires equal time/rate sequences of at least two values")
        if not isinstance(self.source_id, str) or not self.source_id.strip():
            raise ValueError("source schedule source_id must be non-empty")
        if not isinstance(self.rate_operator, str) or self.rate_operator not in {
            "piecewise_constant", "linear"
        }:
            raise ValueError(
                "source schedule rate_operator must be piecewise_constant or linear"
            )
        if any(isinstance(value, bool) for value in (*self.time_s, *self.rate_kg_s)):
            raise ValueError("source schedule time and rate values must be numeric, not boolean")
        times = tuple(float(value) for value in self.time_s)
        rates = tuple(float(value) for value in self.rate_kg_s)
        if not all(math.isfinite(value) for value in (*times, *rates)):
            raise ValueError("source schedule time and rate values must be finite")
        if times[0] != 0.0:
            raise ValueError("source schedule must start at time zero")
        if any(second <= first for first, second in zip(times, times[1:])):
            raise ValueError("source schedule times must be strictly increasing")
        if any(value < 0.0 for value in rates):
            raise ValueError("source schedule rates must be non-negative")
        object.__setattr__(self, "time_s", times)
        object.__setattr__(self, "rate_kg_s", rates)

    @property
    def duration_s(self) -> float:
        return self.time_s[-1]

    @property
    def has_zero_endpoint(self) -> bool:
        """Whether the retained endpoint explicitly shuts the source off."""
        return self.rate_kg_s[-1] == 0.0

    def require_zero_endpoint(self, label: str = "source schedule") -> None:
        """Require an explicit zero endpoint before finite transport use.

        ``SourceRateSchedule`` remains a low-level schedule record and may be
        inspected with a non-zero retained endpoint.  A finite transport
        boundary must, however, state its shutoff explicitly so that a caller
        cannot accidentally interpret the endpoint as continued release.
        """
        if not isinstance(label, str) or not label.strip():
            raise ValueError("source schedule endpoint label must be non-empty")
        if not self.has_zero_endpoint:
            raise ValueError(f"{label} final endpoint rate must be zero")

    @property
    def released_mass_kg(self) -> float:
        if self.rate_operator == "linear":
            return math.fsum(
                0.5 * (left_rate + right_rate) * (end - start)
                for start, end, left_rate, right_rate in zip(
                    self.time_s,
                    self.time_s[1:],
                    self.rate_kg_s,
                    self.rate_kg_s[1:],
                )
            )
        return math.fsum(
            rate * (end - start)
            for start, end, rate in zip(self.time_s, self.time_s[1:], self.rate_kg_s)
        )

    def mass_between(self, start_s: float, end_s: float) -> float:
        """Integrate the declared source rate over one solver time interval."""
        start, end = float(start_s), float(end_s)
        if not (math.isfinite(start) and math.isfinite(end)) or end < start:
            raise ValueError("source-schedule interval must be finite and ordered")
        tolerance = 1.0e-12 * max(1.0, self.duration_s)
        if start < -tolerance or end > self.duration_s + tolerance:
            raise ValueError("source-schedule interval lies outside declared duration")
        start = max(start, 0.0)
        end = min(end, self.duration_s)
        if self.rate_operator == "piecewise_constant":
            return math.fsum(
                rate * max(0.0, min(end, right) - max(start, left))
                for left, right, rate in zip(self.time_s, self.time_s[1:], self.rate_kg_s)
            )
        masses = []
        for left, right, left_rate, right_rate in zip(
            self.time_s,
            self.time_s[1:],
            self.rate_kg_s,
            self.rate_kg_s[1:],
        ):
            clipped_start = max(start, left)
            clipped_end = min(end, right)
            if clipped_end <= clipped_start:
                continue
            span = right - left
            u0 = clipped_start - left
            u1 = clipped_end - left
            slope = (right_rate - left_rate) / span
            masses.append(
                left_rate * (u1 - u0)
                + 0.5 * slope * (u1 * u1 - u0 * u0)
            )
        return math.fsum(masses)


@dataclass(frozen=True)
class DistributedScalarSource:
    """One declared time-varying scalar source inside the local wind plane.

    The mass rate is already atmospheric H2 vapour. It is intentionally not a
    flash, droplet or pool model; callers must provide the physical handoff
    that determines its x/z location and schedule.
    """

    label: str
    x_m: float
    height_m: float
    sigma_m: float
    schedule: SourceRateSchedule

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("distributed scalar source label must be non-empty")
        for name, value, allow_zero in (
            ("x_m", self.x_m, True),
            ("height_m", self.height_m, True),
            ("sigma_m", self.sigma_m, False),
        ):
            if isinstance(value, bool):
                raise ValueError(f"distributed scalar source {name} must be numeric, not boolean")
            if not math.isfinite(float(value)) or (float(value) < 0.0 if allow_zero else float(value) <= 0.0):
                qualifier = "non-negative" if allow_zero else "positive"
                raise ValueError(f"distributed scalar source {name} must be finite and {qualifier}")
        if not isinstance(self.schedule, SourceRateSchedule):
            raise TypeError("distributed scalar source schedule must be a SourceRateSchedule")


@dataclass(frozen=True)
class SemiFVConfig:
    """Numerical boundary for the local wind-plane scalar transport.

    ``gravitational_settling_m_s`` is a signed vertical velocity in the local
    ``z`` coordinate: positive values settle downward (toward decreasing
    height), negative values represent an explicitly supplied upward drift.
    It is a transport sensitivity input, not an inferred LH2 droplet closure.
    """

    length_m: float = 100.0
    height_m: float = 30.0
    nx: int = 200
    nz: int = 60
    time_step_s: float = 0.05
    duration_s: float = 20.0
    wind_speed_m_s: float = 2.0
    diffusivity_m2_s: float = 0.5
    gravitational_settling_m_s: float = 0.0
    source_rate_kg_s: float = 1.0
    source_schedule: SourceRateSchedule | None = None
    source_height_m: float = 0.5
    source_sigma_m: float = 0.5
    obstacle: RectangularObstacle2D | tuple[RectangularObstacle2D, ...] | None = None
    distributed_sources: tuple[DistributedScalarSource, ...] = ()

    def validate(self) -> None:
        for name, value in (("length_m", self.length_m), ("height_m", self.height_m),
                            ("time_step_s", self.time_step_s), ("duration_s", self.duration_s),
                            ("diffusivity_m2_s", self.diffusivity_m2_s),
                            ("source_sigma_m", self.source_sigma_m)):
            if isinstance(value, bool):
                raise ValueError(f"{name} must be positive and not boolean")
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if isinstance(self.nx, bool) or isinstance(self.nz, bool):
            raise ValueError("nx and nz must be integers of at least four")
        if not isinstance(self.nx, int) or not isinstance(self.nz, int) or self.nx < 4 or self.nz < 4:
            raise ValueError("nx and nz must be at least four")
        if isinstance(self.wind_speed_m_s, bool):
            raise ValueError("wind_speed_m_s must be positive and not boolean")
        if not math.isfinite(float(self.wind_speed_m_s)) or self.wind_speed_m_s <= 0.0:
            raise ValueError("wind_speed_m_s must be positive and finite")
        if isinstance(self.gravitational_settling_m_s, bool):
            raise ValueError("gravitational_settling_m_s must be numeric, not boolean")
        if not math.isfinite(float(self.gravitational_settling_m_s)):
            raise ValueError("gravitational_settling_m_s must be finite")
        if isinstance(self.source_rate_kg_s, bool):
            raise ValueError("source_rate_kg_s must be numeric, not boolean")
        if not math.isfinite(float(self.source_rate_kg_s)) or self.source_rate_kg_s < 0.0:
            raise ValueError("source_rate_kg_s must be finite and non-negative")
        if self.source_schedule is not None:
            if not isinstance(self.source_schedule, SourceRateSchedule):
                raise TypeError("source_schedule must be a SourceRateSchedule")
            if self.source_rate_kg_s != 0.0:
                raise ValueError("source_rate_kg_s must be zero when source_schedule is supplied")
            self.source_schedule.require_zero_endpoint("source schedule")
            if self.source_schedule.duration_s > self.duration_s + 1.0e-12:
                raise ValueError("source_schedule duration cannot exceed duration_s")
        if isinstance(self.source_height_m, bool):
            raise ValueError("source_height_m must be numeric, not boolean")
        if not math.isfinite(float(self.source_height_m)):
            raise ValueError("source_height_m must be finite")
        if not 0.0 <= self.source_height_m <= self.height_m:
            raise ValueError("source_height_m must lie inside the domain")
        for obstacle in _obstacles(self):
            if (
                obstacle.x_min_m <= 0.0
                or obstacle.x_max_m >= self.length_m
                or obstacle.z_max_m > self.height_m
            ):
                raise ValueError(
                    "obstacle must lie strictly inside the local domain without "
                    "touching the inlet or outlet"
                )
        if any(not isinstance(source, DistributedScalarSource) for source in self.distributed_sources):
            raise TypeError("distributed_sources must contain only DistributedScalarSource values")
        labels = tuple(source.label for source in self.distributed_sources)
        if len(set(labels)) != len(labels):
            raise ValueError("distributed scalar source labels must be unique")
        for source in self.distributed_sources:
            source.schedule.require_zero_endpoint(
                f"distributed source {source.label!r} schedule"
            )
            if source.x_m > self.length_m or source.height_m > self.height_m:
                raise ValueError("distributed scalar source must lie inside the local domain")
            if source.schedule.duration_s > self.duration_s + 1.0e-12:
                raise ValueError("distributed scalar source schedule duration cannot exceed duration_s")


@dataclass(frozen=True)
class SemiFVDiagnostics:
    mass_injected_kg: float
    mass_domain_kg: float
    mass_outflow_kg: float
    maximum_mass_residual_kg: float
    resolved_time_step_s: float
    x_advection_courant: float
    scalar_diffusivity_m2_s: float
    vertical_diffusion_number: float
    vertical_settling_courant: float
    source_mode: str
    distributed_source_count: int
    obstacle_contact: bool
    obstacle_count: int
    diverted_mass_fraction: float
    applicability: str
    warnings: tuple[str, ...]
    # Added after the required fields so older direct constructors remain
    # source-compatible; solver-produced diagnostics always populate it.
    final_mass_residual_kg: float = 0.0
    source_mass_injected_kg: tuple[tuple[str, float], ...] = ()
    source_mass_ledger_residual_kg: float = 0.0
    source_mass_schedule_residual_kg: tuple[tuple[str, float], ...] = ()
    maximum_source_mass_schedule_residual_kg: float = 0.0

    def __post_init__(self) -> None:
        non_negative = {
            "mass_injected_kg": self.mass_injected_kg,
            "mass_domain_kg": self.mass_domain_kg,
            "mass_outflow_kg": self.mass_outflow_kg,
            "maximum_mass_residual_kg": self.maximum_mass_residual_kg,
            "resolved_time_step_s": self.resolved_time_step_s,
            "x_advection_courant": self.x_advection_courant,
            "scalar_diffusivity_m2_s": self.scalar_diffusivity_m2_s,
            "vertical_diffusion_number": self.vertical_diffusion_number,
            "vertical_settling_courant": self.vertical_settling_courant,
            "maximum_source_mass_schedule_residual_kg": (
                self.maximum_source_mass_schedule_residual_kg
            ),
        }
        for name, value in non_negative.items():
            if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"semi-FV diagnostic {name} must be finite and non-negative")
        if isinstance(self.final_mass_residual_kg, bool) or not math.isfinite(
            float(self.final_mass_residual_kg)
        ):
            raise ValueError("semi-FV diagnostic final_mass_residual_kg must be finite")
        if isinstance(self.source_mass_ledger_residual_kg, bool) or not math.isfinite(
            float(self.source_mass_ledger_residual_kg)
        ):
            raise ValueError("semi-FV diagnostic source_mass_ledger_residual_kg must be finite")
        if self.source_mode not in {
            "constant", "scheduled", "distributed", "primary_plus_distributed",
        }:
            raise ValueError("semi-FV diagnostic source_mode is unsupported")
        if isinstance(self.distributed_source_count, bool) or not isinstance(
            self.distributed_source_count, int
        ) or self.distributed_source_count < 0:
            raise ValueError("semi-FV diagnostic distributed_source_count must be non-negative")
        if not isinstance(self.obstacle_contact, bool):
            raise TypeError("semi-FV diagnostic obstacle_contact must be boolean")
        if isinstance(self.obstacle_count, bool) or not isinstance(
            self.obstacle_count, int
        ) or self.obstacle_count < 0:
            raise ValueError("semi-FV diagnostic obstacle_count must be non-negative")
        if (
            isinstance(self.diverted_mass_fraction, bool)
            or not math.isfinite(float(self.diverted_mass_fraction))
            or not 0.0 <= float(self.diverted_mass_fraction) <= 1.0
        ):
            raise ValueError("semi-FV diagnostic diverted_mass_fraction must lie in [0, 1]")
        if self.applicability not in {"accepted", "conditional", "blocked"}:
            raise ValueError("semi-FV diagnostic applicability is unsupported")
        if not isinstance(self.warnings, tuple) or any(
            not isinstance(item, str) or not item.strip() for item in self.warnings
        ):
            raise ValueError("semi-FV diagnostic warnings must be non-empty strings")
        for name, entries in (
            ("source_mass_injected_kg", self.source_mass_injected_kg),
            ("source_mass_schedule_residual_kg", self.source_mass_schedule_residual_kg),
        ):
            if not isinstance(entries, tuple):
                raise TypeError(f"semi-FV diagnostic {name} must be a tuple")
            labels: set[str] = set()
            for item in entries:
                if not isinstance(item, tuple) or len(item) != 2:
                    raise ValueError(f"semi-FV diagnostic {name} entries must be pairs")
                label, value = item
                if not isinstance(label, str) or not label.strip() or label in labels:
                    raise ValueError(f"semi-FV diagnostic {name} labels must be unique strings")
                if isinstance(value, bool) or not math.isfinite(float(value)):
                    raise ValueError(f"semi-FV diagnostic {name} values must be finite")
                if name == "source_mass_injected_kg" and float(value) < 0.0:
                    raise ValueError(
                        "semi-FV diagnostic source_mass_injected_kg values must be non-negative"
                    )
                labels.add(label)


@dataclass(frozen=True)
class SemiFVReceptor:
    """One local wind-plane receptor at downwind coordinate and height."""

    label: str
    x_m: float
    z_m: float

    def validate(self, config: SemiFVConfig) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("receptor label must be non-empty")
        if any(isinstance(value, bool) for value in (self.x_m, self.z_m)):
            raise ValueError("receptor coordinates must be numeric, not boolean")
        if not all(math.isfinite(float(value)) for value in (self.x_m, self.z_m)):
            raise ValueError("receptor coordinates must be finite")
        if not 0.0 <= self.x_m <= config.length_m or not 0.0 <= self.z_m <= config.height_m:
            raise ValueError("receptor must lie inside the local domain")


@dataclass(frozen=True)
class SemiFVReceptorTrace:
    """True scalar-concentration history at one local receptor."""

    receptor: SemiFVReceptor
    time_s: np.ndarray
    concentration_kg_m3: np.ndarray


@dataclass(frozen=True)
class SemiFVResult:
    x_m: np.ndarray
    z_m: np.ndarray
    concentration_kg_m3: np.ndarray
    diagnostics: SemiFVDiagnostics
    receptor_traces: tuple[SemiFVReceptorTrace, ...] = ()

    @property
    def maximum_concentration_kg_m3(self) -> float:
        return float(np.nanmax(self.concentration_kg_m3))

    def concentration_at(self, x_m: float, z_m: float) -> float:
        if not (math.isfinite(float(x_m)) and math.isfinite(float(z_m))):
            raise ValueError("receptor coordinates must be finite")
        i = int(np.clip(np.searchsorted(self.x_m, x_m), 0, len(self.x_m) - 1))
        j = int(np.clip(np.searchsorted(self.z_m, z_m), 0, len(self.z_m) - 1))
        return float(self.concentration_kg_m3[j, i])


@dataclass(frozen=True)
class SemiFVReceptorRefinement:
    """Change in peak and dose between adjacent numerical refinements."""

    receptor_label: str
    coarse_factor: int
    fine_factor: int
    coarse_peak_kg_m3: float
    fine_peak_kg_m3: float
    peak_relative_change: float
    coarse_dose_kg_s_m3: float
    fine_dose_kg_s_m3: float
    dose_relative_change: float


@dataclass(frozen=True)
class SemiFVRefinementStudy:
    """Explicit grid/time sensitivity result for declared local receptors.

    ``converged`` is only a numerical-resolution screen for the reduced-order
    operator. It does not validate the obstacle closure or make it CFD.
    """

    cases: tuple[tuple[int, SemiFVResult], ...]
    receptor_changes: tuple[SemiFVReceptorRefinement, ...]
    relative_tolerance: float
    converged: bool
    estimated_cell_steps: int
    warnings: tuple[str, ...] = ()


def _obstacle_mask(config: SemiFVConfig, x: np.ndarray, z: np.ndarray) -> np.ndarray:
    mask = np.zeros((len(z), len(x)), dtype=bool)
    for obstacle in _obstacles(config):
        mask |= (
            (x[None, :] >= obstacle.x_min_m) & (x[None, :] <= obstacle.x_max_m)
            & (z[:, None] >= obstacle.base_height_m)
            & (z[:, None] <= obstacle.z_max_m)
        )
    return mask


def _obstacles(config: SemiFVConfig) -> tuple[RectangularObstacle2D, ...]:
    """Normalise the backward-compatible one-or-many obstacle contract."""
    declared = config.obstacle
    if declared is None:
        return ()
    if isinstance(declared, RectangularObstacle2D):
        return (declared,)
    if not isinstance(declared, tuple) or not declared:
        raise TypeError("obstacle must be one RectangularObstacle2D or a non-empty tuple")
    if any(not isinstance(item, RectangularObstacle2D) for item in declared):
        raise TypeError("every declared obstacle must be a RectangularObstacle2D")
    labels = tuple(item.label for item in declared)
    if len(set(labels)) != len(labels):
        raise ValueError("declared obstacle labels must be unique")
    return declared


def _route_target(j: int, i: int, solid: np.ndarray) -> tuple[int, int] | None:
    """Find the first fluid cell downstream of a blocked horizontal face.

    The search is intentionally local and deterministic: it first looks in the
    next column and then walks downstream until it finds a clear cell, testing
    vertical offsets from the donor cell outwards.  This is a conservative
    obstacle-routing closure, not a claim to resolve a recirculating wake.
    """
    nz, nx = solid.shape
    for target_i in range(i + 1, nx):
        for offset in range(nz):
            candidates = (j,) if offset == 0 else (j - offset, j + offset)
            for target_j in candidates:
                if 0 <= target_j < nz and not solid[target_j, target_i]:
                    return target_j, target_i
    return None


def _schedule_mass(
    schedule: SourceRateSchedule,
    start_time_s: float,
    end_time_s: float,
) -> float:
    """Integrate a finite schedule over one solver step, with zero after end."""
    if start_time_s >= schedule.duration_s:
        return 0.0
    return schedule.mass_between(start_time_s, min(end_time_s, schedule.duration_s))


def _inject_gaussian_source(
    concentration: np.ndarray,
    *,
    solid: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    cell_volume_m3: float,
    mass_kg: float,
    x_m: float,
    height_m: float,
    sigma_m: float,
    label: str,
) -> None:
    """Inject one cell-centred vertical Gaussian without source-in-solid loss."""
    if mass_kg <= 0.0:
        return
    column = int(np.clip(np.searchsorted(x, x_m), 0, len(x) - 1))
    centre = int(np.clip(np.searchsorted(z, height_m), 0, len(z) - 1))
    if solid[centre, column]:
        raise ValueError(f"distributed scalar source {label!r} lies inside a declared solid obstacle")
    weights = np.exp(-0.5 * ((z - height_m) / sigma_m) ** 2)
    weights[solid[:, column]] = 0.0
    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError(f"distributed scalar source {label!r} has no fluid injection cell")
    concentration[:, column] += mass_kg * weights / total / cell_volume_m3


def _source_mode(config: SemiFVConfig) -> str:
    """Describe whether the primary and/or internal source contracts are used."""
    if config.distributed_sources:
        if config.source_rate_kg_s == 0.0 and config.source_schedule is None:
            return "distributed"
        return "primary_plus_distributed"
    return "constant" if config.source_schedule is None else "scheduled"


def solve_semi_fv_obstacle(
    config: SemiFVConfig,
    *,
    receptors: Sequence[SemiFVReceptor] = (),
) -> SemiFVResult:
    """Solve the local reduced-order transport problem."""
    config.validate()
    dx, dz = config.length_m / config.nx, config.height_m / config.nz
    x = (np.arange(config.nx) + 0.5) * dx
    z = (np.arange(config.nz) + 0.5) * dz
    solid = _obstacle_mask(config, x, z)
    if np.any(solid[:, 0]) or np.any(solid[:, -1]):
        raise ValueError("obstacle may not touch the inlet or outlet boundary")
    checked_receptors = tuple(receptors)
    labels = set()
    for receptor in checked_receptors:
        if not isinstance(receptor, SemiFVReceptor):
            raise TypeError("every receptor must be a SemiFVReceptor")
        receptor.validate(config)
        if receptor.label in labels:
            raise ValueError("receptor labels must be unique")
        labels.add(receptor.label)
        receptor_i = int(np.clip(np.searchsorted(x, receptor.x_m), 0, config.nx - 1))
        receptor_j = int(np.clip(np.searchsorted(z, receptor.z_m), 0, config.nz - 1))
        if solid[receptor_j, receptor_i]:
            raise ValueError("receptor lies inside declared solid obstacle")
    concentration = np.zeros((config.nz, config.nx), dtype=float)
    cell_volume = dx * dz * 1.0
    total_in = 0.0
    total_out = 0.0
    maximum_residual = 0.0
    diverted_total = 0.0
    source_mass_injected: dict[str, float] = {
        "primary": 0.0,
        **{f"distributed:{source.label}": 0.0 for source in config.distributed_sources},
    }
    expected_source_mass: dict[str, float] = {
        "primary": (
            config.source_rate_kg_s * config.duration_s
            if config.source_schedule is None
            else config.source_schedule.released_mass_kg
        ),
        **{
            f"distributed:{source.label}": source.schedule.released_mass_kg
            for source in config.distributed_sources
        },
    }
    steps = int(math.ceil(config.duration_s / config.time_step_s))
    dt = config.duration_s / steps
    x_courant = config.wind_speed_m_s * dt / dx
    vertical_diffusion = config.diffusivity_m2_s * dt / dz**2
    vertical_settling = abs(config.gravitational_settling_m_s) * dt / dz
    if x_courant > 1.0 + 1.0e-12:
        raise ValueError(
            "semi-FV x-advection Courant number exceeds one; reduce time_step_s or refine nx"
        )
    if 2.0 * vertical_diffusion + vertical_settling > 1.0 + 1.0e-12:
        raise ValueError(
            "semi-FV vertical diffusion/settling stability limit is exceeded; reduce time_step_s or refine nz"
        )
    trace_time = [0.0]
    trace_values = [[0.0] for _ in checked_receptors]

    def sample(receptor: SemiFVReceptor) -> float:
        i = int(np.clip(np.searchsorted(x, receptor.x_m), 0, config.nx - 1))
        j = int(np.clip(np.searchsorted(z, receptor.z_m), 0, config.nz - 1))
        return float(concentration[j, i])

    for step in range(steps):
        before = float(np.sum(concentration) * cell_volume)
        start_time = step * dt
        end_time = (step + 1) * dt
        source = (
            config.source_rate_kg_s * dt
            if config.source_schedule is None
            else _schedule_mass(config.source_schedule, start_time, end_time)
        )
        _inject_gaussian_source(
            concentration, solid=solid, x=x, z=z, cell_volume_m3=cell_volume,
            mass_kg=source, x_m=0.0, height_m=config.source_height_m,
            sigma_m=config.source_sigma_m, label="primary",
        )
        total_in += source
        source_mass_injected["primary"] += source
        distributed_mass_total = 0.0
        for distributed in config.distributed_sources:
            distributed_mass = _schedule_mass(distributed.schedule, start_time, end_time)
            _inject_gaussian_source(
                concentration, solid=solid, x=x, z=z, cell_volume_m3=cell_volume,
                mass_kg=distributed_mass, x_m=distributed.x_m,
                height_m=distributed.height_m, sigma_m=distributed.sigma_m,
                label=distributed.label,
            )
            total_in += distributed_mass
            source_mass_injected[f"distributed:{distributed.label}"] += distributed_mass
            distributed_mass_total += distributed_mass
        # x advection is evaluated from the old field and applied through a
        # conservative increment.  The explicit old-field evaluation prevents
        # a parcel from crossing multiple cells in one time step.
        old = concentration.copy()
        increment = np.zeros_like(old)
        diverted = 0.0
        out_step = 0.0
        for j in range(config.nz):
            for i in range(config.nx):
                if solid[j, i] or old[j, i] <= 0.0:
                    continue
                amount = config.wind_speed_m_s * old[j, i] * dz * dt
                if i == config.nx - 1:
                    increment[j, i] -= amount / cell_volume
                    total_out += amount
                    out_step += amount
                    continue
                if not solid[j, i + 1]:
                    increment[j, i] -= amount / cell_volume
                    increment[j, i + 1] += amount / cell_volume
                    continue
                target = _route_target(j, i, solid)
                if target is None:
                    # Keep the mass in the donor cell when the obstacle spans
                    # the complete available height; this is fail-closed.
                    continue
                target_j, target_i = target
                increment[j, i] -= amount / cell_volume
                increment[target_j, target_i] += amount / cell_volume
                diverted += amount
        concentration += increment
        diverted_total += diverted

        # Vertical settling and turbulent diffusion use antisymmetric internal
        # face fluxes.  Solid cells block the face; no wake concentration is
        # invented behind the obstacle.
        old = concentration.copy()
        vertical_increment = np.zeros_like(old)
        diff = config.diffusivity_m2_s * dt * dx / max(dz, 1.0e-30)
        # The public boundary uses the conventional positive-downward settling
        # sign, while the cell index increases with z.  The signed face flux
        # below therefore receives the negated velocity coefficient.
        settling = -config.gravitational_settling_m_s * dt * dx
        for j in range(config.nz - 1):
            for i in range(config.nx):
                if solid[j, i] or solid[j + 1, i]:
                    continue
                amount = diff * (old[j, i] - old[j + 1, i])
                if settling > 0.0:
                    amount += settling * old[j, i]
                elif settling < 0.0:
                    amount += settling * old[j + 1, i]
                vertical_increment[j, i] -= amount / cell_volume
                vertical_increment[j + 1, i] += amount / cell_volume
        concentration += vertical_increment
        concentration[solid] = 0.0
        concentration = np.maximum(concentration, 0.0)
        after = float(np.sum(concentration) * cell_volume)
        expected = before + source + distributed_mass_total - out_step
        maximum_residual = max(maximum_residual, abs(after - expected))
        trace_time.append((step + 1) * dt)
        for values, receptor in zip(trace_values, checked_receptors):
            values.append(sample(receptor))
    obstacles = _obstacles(config)
    obstacle_contact = bool(np.any(solid))
    fraction = diverted_total / max(total_in, 1.0e-30)
    source_mass_schedule_residual = tuple(
        (
            label,
            float(source_mass_injected[label] - expected_source_mass[label]),
        )
        for label in source_mass_injected
    )
    maximum_source_mass_schedule_residual = max(
        (abs(residual) for _label, residual in source_mass_schedule_residual),
        default=0.0,
    )
    warnings = (
        "reduced-order obstacle routing; not obstacle-resolved CFD",
        "no LH2 obstacle concentration validation is included",
    )
    applicability = "conditional" if obstacle_contact else "accepted"
    mass_domain = float(np.sum(concentration) * cell_volume)
    final_mass_residual = float(total_in - mass_domain - total_out)
    diagnostics = SemiFVDiagnostics(
        mass_injected_kg=total_in,
        mass_domain_kg=mass_domain,
        mass_outflow_kg=total_out,
        maximum_mass_residual_kg=maximum_residual,
        final_mass_residual_kg=final_mass_residual,
        resolved_time_step_s=dt,
        x_advection_courant=x_courant,
        scalar_diffusivity_m2_s=config.diffusivity_m2_s,
        vertical_diffusion_number=vertical_diffusion,
        vertical_settling_courant=vertical_settling,
        source_mode=_source_mode(config),
        distributed_source_count=len(config.distributed_sources),
        obstacle_contact=obstacle_contact,
        obstacle_count=len(obstacles),
        diverted_mass_fraction=float(np.clip(fraction, 0.0, 1.0)),
        applicability=applicability,
        warnings=warnings,
        source_mass_injected_kg=tuple(
            (label, float(mass)) for label, mass in source_mass_injected.items()
        ),
        source_mass_ledger_residual_kg=float(
            total_in - math.fsum(source_mass_injected.values())
        ),
        source_mass_schedule_residual_kg=source_mass_schedule_residual,
        maximum_source_mass_schedule_residual_kg=float(
            maximum_source_mass_schedule_residual
        ),
    )
    concentration[solid] = np.nan
    traces = tuple(
        SemiFVReceptorTrace(
            receptor, np.asarray(trace_time), np.asarray(values),
        )
        for receptor, values in zip(checked_receptors, trace_values)
    )
    return SemiFVResult(x, z, concentration, diagnostics, traces)


def run_semi_fv_refinement_study(
    config: SemiFVConfig,
    *,
    receptors: Sequence[SemiFVReceptor],
    refinement_factors: Sequence[int] = (1, 2),
    relative_tolerance: float = 0.05,
    absolute_floor_kg_m3: float = 1.0e-12,
    max_cell_steps: int = 20_000_000,
) -> SemiFVRefinementStudy:
    """Compare local receptor peak and dose across explicit grid/time refinements.

    Each factor multiplies both grid dimensions while dividing the requested
    time step by factor-squared. This preserves the explicit vertical-diffusion
    stability scale. A cell-step budget prevents a screening request from
    silently becoming an impractically large calculation.
    """
    config.validate()
    factors = tuple(refinement_factors)
    if len(factors) < 2:
        raise ValueError("at least two refinement_factors are required")
    if any(isinstance(factor, bool) or not isinstance(factor, int) or factor < 1 for factor in factors):
        raise ValueError("refinement_factors must be positive integers")
    if tuple(sorted(factors)) != factors or len(set(factors)) != len(factors):
        raise ValueError("refinement_factors must be unique and strictly increasing")
    if not receptors:
        raise ValueError("at least one receptor is required for a refinement study")
    if isinstance(relative_tolerance, bool) or not math.isfinite(float(relative_tolerance)) or relative_tolerance < 0.0:
        raise ValueError("relative_tolerance must be finite and non-negative")
    if isinstance(absolute_floor_kg_m3, bool) or not math.isfinite(float(absolute_floor_kg_m3)) or absolute_floor_kg_m3 <= 0.0:
        raise ValueError("absolute_floor_kg_m3 must be positive and finite")
    if not isinstance(max_cell_steps, int) or max_cell_steps <= 0:
        raise ValueError("max_cell_steps must be a positive integer")

    refined_configs = tuple(
        (
            factor,
            replace(
                config,
                nx=config.nx * factor,
                nz=config.nz * factor,
                time_step_s=config.time_step_s / factor**2,
            ),
        )
        for factor in factors
    )
    estimated_cell_steps = sum(
        item.nx * item.nz * int(math.ceil(item.duration_s / item.time_step_s))
        for _, item in refined_configs
    )
    if estimated_cell_steps > max_cell_steps:
        raise ValueError(
            "semi-FV refinement study exceeds max_cell_steps; reduce factors/domain/duration "
            "or explicitly raise the budget"
        )

    cases = tuple(
        (factor, solve_semi_fv_obstacle(item, receptors=receptors))
        for factor, item in refined_configs
    )
    changes = []
    for (coarse_factor, coarse), (fine_factor, fine) in zip(cases, cases[1:]):
        coarse_traces = {trace.receptor.label: trace for trace in coarse.receptor_traces}
        fine_traces = {trace.receptor.label: trace for trace in fine.receptor_traces}
        for label in sorted(coarse_traces):
            first, second = coarse_traces[label], fine_traces[label]
            first_peak = float(np.max(first.concentration_kg_m3))
            second_peak = float(np.max(second.concentration_kg_m3))
            first_dose = float(np.trapezoid(first.concentration_kg_m3, first.time_s))
            second_dose = float(np.trapezoid(second.concentration_kg_m3, second.time_s))
            peak_change = abs(second_peak - first_peak) / max(
                abs(first_peak), abs(second_peak), absolute_floor_kg_m3,
            )
            dose_change = abs(second_dose - first_dose) / max(
                abs(first_dose), abs(second_dose), absolute_floor_kg_m3,
            )
            changes.append(SemiFVReceptorRefinement(
                label, coarse_factor, fine_factor,
                first_peak, second_peak, peak_change,
                first_dose, second_dose, dose_change,
            ))
    converged = all(
        max(item.peak_relative_change, item.dose_relative_change) <= relative_tolerance
        for item in changes
    )
    warnings = () if converged else (
        "declared receptor peak or dose did not meet the requested numerical refinement tolerance",
    )
    return SemiFVRefinementStudy(
        cases, tuple(changes), float(relative_tolerance), converged,
        estimated_cell_steps, warnings,
    )


__all__ = [
    "RectangularObstacle2D", "SourceRateSchedule", "DistributedScalarSource",
    "SemiFVConfig", "SemiFVDiagnostics",
    "SemiFVReceptor", "SemiFVReceptorTrace", "SemiFVResult",
    "SemiFVReceptorRefinement", "SemiFVRefinementStudy",
    "solve_semi_fv_obstacle", "run_semi_fv_refinement_study",
]
