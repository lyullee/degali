"""Geometry screens for using an integral plume near site structures.

This module deliberately does *not* invent a building-wake dilution factor.
The Gaussian/integral closures in :mod:`degali.addons` assume an unobstructed
control volume.  A solid obstacle can split that volume, generate separation
and recirculation, and exchange heat with a cryogenic cloud.  Those effects
cannot be determined from a cuboid's dimensions alone.  The functions here
therefore make the geometry explicit and fail closed when a calculated plume
centre trajectory enters a solid object.

It is useful before an obstacle-resolved calculation in two ways: it accepts
arbitrary horizontal wind bearings, and it reports the exact first contact of
the modelled centre trajectory with a cuboid or a cross-wind wall.  It changes
no source flux, entrainment, temperature, or concentration.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Protocol, Sequence

import numpy as np


_GEOMETRY_TOLERANCE = 1.0e-12


@dataclass(frozen=True)
class WindFrame:
    """Horizontal coordinate system directed *towards* the mean wind.

    ``direction_rad`` is measured counter-clockwise from the global ``+x``
    axis.  It is a *to* direction, not a meteorological ``from`` direction.
    This explicit convention avoids silently reversing an obstacle layout.
    ``s`` increases downwind and ``n`` increases to the left of the wind.
    """

    origin_x_m: float = 0.0
    origin_y_m: float = 0.0
    direction_rad: float = 0.0

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (
            self.origin_x_m, self.origin_y_m, self.direction_rad
        )):
            raise ValueError("wind-frame origin and direction must be finite")

    @property
    def downwind_unit(self) -> tuple[float, float]:
        return math.cos(self.direction_rad), math.sin(self.direction_rad)

    @property
    def lateral_unit(self) -> tuple[float, float]:
        return -math.sin(self.direction_rad), math.cos(self.direction_rad)

    def local(self, x_m: float, y_m: float) -> tuple[float, float]:
        """Return downwind/lateral coordinates for a global horizontal point."""
        if not math.isfinite(x_m) or not math.isfinite(y_m):
            raise ValueError("global coordinates must be finite")
        dx, dy = x_m-self.origin_x_m, y_m-self.origin_y_m
        ex, ey = self.downwind_unit
        nx, ny = self.lateral_unit
        return dx*ex + dy*ey, dx*nx + dy*ny

    def global_(self, downwind_m: float, lateral_m: float) -> tuple[float, float]:
        """Return a global horizontal point from downwind/lateral coordinates."""
        if not math.isfinite(downwind_m) or not math.isfinite(lateral_m):
            raise ValueError("wind-frame coordinates must be finite")
        ex, ey = self.downwind_unit
        nx, ny = self.lateral_unit
        return (
            self.origin_x_m + downwind_m*ex + lateral_m*nx,
            self.origin_y_m + downwind_m*ey + lateral_m*ny,
        )


@dataclass(frozen=True)
class SegmentIntersection:
    """Closed line-segment contact with a solid, parameterised on ``[0, 1]``."""

    entry_fraction: float
    exit_fraction: float
    entry_point_m: tuple[float, float, float]
    exit_point_m: tuple[float, float, float]


class _SegmentObstacle(Protocol):
    label: str

    def segment_intersection(
        self, start_m: Sequence[float], end_m: Sequence[float]
    ) -> SegmentIntersection | None:
        ...


def _point3(value: Sequence[float], *, name: str) -> np.ndarray:
    point = np.asarray(value, dtype=float)
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError(f"{name} must contain exactly three finite coordinates")
    return point


def _segment_aabb_intersection(
    start_m: Sequence[float], end_m: Sequence[float], lower_m: Sequence[float],
    upper_m: Sequence[float],
) -> SegmentIntersection | None:
    """Exact slab intersection of a closed segment and axis-aligned box.

    Infinite lateral bounds are accepted for an idealised wall.  The method is
    geometric only: touching a face counts as contact because even a grazing
    trajectory invalidates an unobstructed control volume.
    """
    start, end = _point3(start_m, name="segment start"), _point3(end_m, name="segment end")
    lower, upper = np.asarray(lower_m, float), np.asarray(upper_m, float)
    if (lower.shape != (3,) or upper.shape != (3,)
            or np.any(np.isnan(lower)) or np.any(np.isnan(upper))
            or np.any(lower > upper)):
        raise ValueError("box bounds must be ordered real or infinite values")
    delta = end-start
    entry, exit_ = 0.0, 1.0
    for index in range(3):
        if abs(delta[index]) <= _GEOMETRY_TOLERANCE:
            if start[index] < lower[index]-_GEOMETRY_TOLERANCE or start[index] > upper[index]+_GEOMETRY_TOLERANCE:
                return None
            continue
        first = (lower[index]-start[index])/delta[index]
        second = (upper[index]-start[index])/delta[index]
        entry, exit_ = max(entry, min(first, second)), min(exit_, max(first, second))
        if entry > exit_ + _GEOMETRY_TOLERANCE:
            return None
    entry, exit_ = max(0.0, entry), min(1.0, exit_)
    if entry > exit_ + _GEOMETRY_TOLERANCE:
        return None
    entry_point, exit_point = start+entry*delta, start+exit_*delta
    return SegmentIntersection(
        float(entry), float(exit_), tuple(float(value) for value in entry_point),
        tuple(float(value) for value in exit_point),
    )


@dataclass(frozen=True)
class AxisAlignedCuboid:
    """Closed, solid global-coordinate cuboid representing a site obstacle."""

    x_min_m: float
    x_max_m: float
    y_min_m: float
    y_max_m: float
    z_min_m: float
    z_max_m: float
    label: str = "cuboid"

    def __post_init__(self) -> None:
        bounds = (self.x_min_m, self.x_max_m, self.y_min_m, self.y_max_m,
                  self.z_min_m, self.z_max_m)
        if not all(math.isfinite(value) for value in bounds):
            raise ValueError("cuboid bounds must be finite")
        if self.x_min_m >= self.x_max_m or self.y_min_m >= self.y_max_m or self.z_min_m >= self.z_max_m:
            raise ValueError("a cuboid requires strictly positive extent on every axis")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("obstacle label must be a non-empty string")

    def segment_intersection(
        self, start_m: Sequence[float], end_m: Sequence[float]
    ) -> SegmentIntersection | None:
        return _segment_aabb_intersection(
            start_m, end_m,
            (self.x_min_m, self.y_min_m, self.z_min_m),
            (self.x_max_m, self.y_max_m, self.z_max_m),
        )


@dataclass(frozen=True)
class TransverseWall:
    """Solid wall normal to a :class:`WindFrame` downwind axis.

    ``half_width_m=None`` is an ideal laterally unbounded wall.  If the plume
    centre hits such a wall, any continuous centre path would have to clear
    ``top_height_m``; this is a geometric necessity, *not* a prediction that
    the cloud climbs to that height.  A finite wall may be passed around as
    well as over, so geometry alone cannot select a path.
    """

    frame: WindFrame
    downwind_m: float
    base_height_m: float
    height_m: float
    half_width_m: float | None = None
    thickness_m: float = 0.0
    label: str = "transverse wall"

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (
            self.downwind_m, self.base_height_m, self.height_m, self.thickness_m
        )):
            raise ValueError("wall location, height and thickness must be finite")
        if self.height_m <= 0.0 or self.thickness_m < 0.0:
            raise ValueError("wall height must be positive and thickness non-negative")
        if self.half_width_m is not None and (
            not math.isfinite(self.half_width_m) or self.half_width_m <= 0.0
        ):
            raise ValueError("wall half-width must be positive and finite, or None")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("wall label must be a non-empty string")

    @property
    def top_height_m(self) -> float:
        return self.base_height_m+self.height_m

    @property
    def is_laterally_unbounded(self) -> bool:
        return self.half_width_m is None

    def segment_intersection(
        self, start_m: Sequence[float], end_m: Sequence[float]
    ) -> SegmentIntersection | None:
        start, end = _point3(start_m, name="segment start"), _point3(end_m, name="segment end")
        start_s, start_n = self.frame.local(start[0], start[1])
        end_s, end_n = self.frame.local(end[0], end[1])
        width = math.inf if self.half_width_m is None else self.half_width_m
        local_hit = _segment_aabb_intersection(
            (start_s, start_n, start[2]), (end_s, end_n, end[2]),
            (self.downwind_m-self.thickness_m/2.0, -width, self.base_height_m),
            (self.downwind_m+self.thickness_m/2.0, width, self.top_height_m),
        )
        if local_hit is None:
            return None
        # Return global points so all geometry screens share one convention.
        delta = end-start
        entry = start+local_hit.entry_fraction*delta
        exit_ = start+local_hit.exit_fraction*delta
        return SegmentIntersection(
            local_hit.entry_fraction, local_hit.exit_fraction,
            tuple(float(value) for value in entry), tuple(float(value) for value in exit_),
        )


@dataclass(frozen=True)
class ObstacleEncounter:
    """One contiguous trajectory encounter with a named solid obstacle."""

    obstacle_label: str
    obstacle_type: str
    first_segment_index: int
    last_segment_index: int
    entry_arc_length_m: float
    exit_arc_length_m: float
    entry_point_m: tuple[float, float, float]
    exit_point_m: tuple[float, float, float]
    laterally_unbounded_wall: bool
    minimum_overflight_height_m: float | None


@dataclass(frozen=True)
class ObstacleScreenResult:
    """Applicability result; it contains no obstacle-corrected concentration."""

    encounters: tuple[ObstacleEncounter, ...]
    free_plume_prediction_applicable: bool
    requires_obstacle_resolved_model: bool
    required_overflight_height_m: float | None
    notes: tuple[str, ...]

    @property
    def is_clear(self) -> bool:
        return not self.encounters

    @property
    def first_contact_arc_length_m(self) -> float | None:
        """Return first solid contact along the supplied trajectory, if any."""
        return self.encounters[0].entry_arc_length_m if self.encounters else None

    def free_plume_applicable_through(self, arc_length_m: float) -> bool:
        """Whether the trajectory is clear from its start through this arclength.

        This permits upstream receptors to retain the unobstructed prediction.
        It does not restore applicability after an obstacle is encountered.
        """
        if not math.isfinite(arc_length_m) or arc_length_m < 0.0:
            raise ValueError("arclength must be finite and non-negative")
        return all(
            encounter.entry_arc_length_m > arc_length_m+_GEOMETRY_TOLERANCE
            for encounter in self.encounters
        )


def _obstacle_metadata(obstacle: _SegmentObstacle) -> tuple[str, bool, float | None]:
    if isinstance(obstacle, TransverseWall):
        return "transverse_wall", obstacle.is_laterally_unbounded, (
            obstacle.top_height_m if obstacle.is_laterally_unbounded else None
        )
    if isinstance(obstacle, AxisAlignedCuboid):
        return "axis_aligned_cuboid", False, None
    return type(obstacle).__name__, False, None


def screen_trajectory(
    points_m: Sequence[Sequence[float]], obstacles: Sequence[_SegmentObstacle]
) -> ObstacleScreenResult:
    """Screen a centre-trajectory polyline against declared solid geometry.

    The coordinates are global ``(x, y, z)`` metres.  The resulting flag is
    intentionally strict: any contact makes a free-plume integral prediction
    inapplicable after the contact.  It does not overwrite a baseline field,
    because neither bypass partition nor wake mixing is identifiable here.
    """
    points = np.asarray(points_m, dtype=float)
    if points.ndim != 2 or points.shape[0] < 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("trajectory requires at least two finite (x, y, z) points")
    checked = tuple(obstacles)
    if any(not isinstance(getattr(obstacle, "label", None), str) or not obstacle.label.strip()
           or not callable(getattr(obstacle, "segment_intersection", None)) for obstacle in checked):
        raise ValueError("each obstacle requires a non-empty label and segment_intersection")

    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if np.any(segment_lengths <= _GEOMETRY_TOLERANCE):
        raise ValueError("trajectory points must be separated by a positive distance")
    arc_start = np.r_[0.0, np.cumsum(segment_lengths)]
    encounters: list[ObstacleEncounter] = []

    for obstacle in checked:
        obstacle_type, full_wall, overflight = _obstacle_metadata(obstacle)
        current: ObstacleEncounter | None = None
        for index, (start, end, length) in enumerate(zip(points[:-1], points[1:], segment_lengths)):
            hit = obstacle.segment_intersection(start, end)
            if hit is None:
                if current is not None:
                    encounters.append(current)
                    current = None
                continue
            entry_arc = float(arc_start[index]+hit.entry_fraction*length)
            exit_arc = float(arc_start[index]+hit.exit_fraction*length)
            contiguous = current is not None and entry_arc <= current.exit_arc_length_m+_GEOMETRY_TOLERANCE
            if contiguous:
                current = ObstacleEncounter(
                    current.obstacle_label, current.obstacle_type,
                    current.first_segment_index, index,
                    current.entry_arc_length_m, exit_arc,
                    current.entry_point_m, hit.exit_point_m,
                    current.laterally_unbounded_wall, current.minimum_overflight_height_m,
                )
            else:
                if current is not None:
                    encounters.append(current)
                current = ObstacleEncounter(
                    obstacle.label, obstacle_type, index, index, entry_arc, exit_arc,
                    hit.entry_point_m, hit.exit_point_m, full_wall, overflight,
                )
        if current is not None:
            encounters.append(current)

    encounters.sort(key=lambda encounter: (encounter.entry_arc_length_m, encounter.obstacle_label))
    heights = [encounter.minimum_overflight_height_m for encounter in encounters
               if encounter.minimum_overflight_height_m is not None]
    required_height = max(heights) if heights else None
    if encounters:
        notes = (
            "A calculated plume-centre trajectory enters declared solid geometry; "
            "the unobstructed integral-plume field is not applicable downstream of first contact.",
            "No obstacle wake, bypass partition, wall heat transfer, reflection, or concentration correction was applied.",
        )
    else:
        notes = (
            "No declared solid intersects the calculated centre trajectory.",
            "This screen does not prove that a finite Gaussian plume envelope clears an obstacle.",
        )
    return ObstacleScreenResult(
        tuple(encounters), not bool(encounters), bool(encounters), required_height, notes
    )
