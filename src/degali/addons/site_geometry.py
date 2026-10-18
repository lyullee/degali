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


def segment_circular_arc_surface_intersections(
    start_m: Sequence[float],
    end_m: Sequence[float],
    *,
    center_xy_m: Sequence[float],
    radius_m: float,
    z_min_m: float,
    z_max_m: float,
    arc_start_deg: float,
    arc_end_deg: float,
) -> tuple[SegmentIntersection, ...]:
    """Return exact segment contacts with a zero-thickness vertical arc.

    This is a geometric primitive for embedded/immersed wall treatments. It
    deliberately does not voxelise the arc or assign it a numerical thickness.
    Each result is a point contact (equal entry/exit fraction); a chord may
    therefore produce two ordered contacts and a tangent produces one.
    """
    start = _point3(start_m, name="segment start")
    end = _point3(end_m, name="segment end")
    center = np.asarray(center_xy_m, dtype=float)
    values = (
        radius_m, z_min_m, z_max_m, arc_start_deg, arc_end_deg,
    )
    if center.shape != (2,) or not np.all(np.isfinite(center)):
        raise ValueError("center_xy_m must contain exactly two finite coordinates")
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("arc dimensions and angles must be finite")
    if radius_m <= 0.0 or z_min_m > z_max_m:
        raise ValueError("arc radius must be positive and vertical bounds ordered")
    span = (arc_end_deg - arc_start_deg) % 360.0
    if span <= _GEOMETRY_TOLERANCE:
        raise ValueError("arc angular span must lie in (0, 360) degrees")

    delta = end - start
    quadratic = float(delta[0] ** 2 + delta[1] ** 2)
    if quadratic <= _GEOMETRY_TOLERANCE**2:
        return ()
    relative = start[:2] - center
    linear = 2.0 * float(np.dot(relative, delta[:2]))
    constant = float(np.dot(relative, relative) - radius_m**2)
    discriminant = linear**2 - 4.0 * quadratic * constant
    discriminant_tolerance = _GEOMETRY_TOLERANCE * max(
        1.0, linear**2, abs(4.0 * quadratic * constant), quadratic * radius_m**2,
    )
    if discriminant < -discriminant_tolerance:
        return ()
    root = math.sqrt(max(0.0, discriminant))
    fractions = sorted((
        (-linear - root) / (2.0 * quadratic),
        (-linear + root) / (2.0 * quadratic),
    ))
    contacts: list[SegmentIntersection] = []
    for fraction in fractions:
        if fraction < -_GEOMETRY_TOLERANCE or fraction > 1.0 + _GEOMETRY_TOLERANCE:
            continue
        fraction = min(1.0, max(0.0, fraction))
        point = start + fraction * delta
        if not z_min_m - _GEOMETRY_TOLERANCE <= point[2] <= z_max_m + _GEOMETRY_TOLERANCE:
            continue
        angle = math.degrees(math.atan2(point[1] - center[1], point[0] - center[0])) % 360.0
        angle_offset = (angle - arc_start_deg) % 360.0
        if angle_offset > span + _GEOMETRY_TOLERANCE:
            continue
        point_tuple = tuple(float(value) for value in point)
        if contacts and abs(fraction - contacts[-1].entry_fraction) <= _GEOMETRY_TOLERANCE:
            continue
        contacts.append(SegmentIntersection(
            float(fraction), float(fraction), point_tuple, point_tuple,
        ))
    return tuple(contacts)


def circular_arc_surface_face_masks(
    x_m: Sequence[float],
    y_m: Sequence[float],
    z_m: Sequence[float],
    *,
    center_xy_m: Sequence[float],
    radius_m: float,
    z_min_m: float,
    z_max_m: float,
    arc_start_deg: float,
    arc_end_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mark Cartesian faces crossed by an exact zero-thickness vertical arc.

    The returned boolean masks have the standard finite-volume shapes
    ``(nz, ny, nx+1)``, ``(nz, ny+1, nx)``, and ``(nz+1, ny, nx)``. A face is
    blocked when the mathematical arc intersects its open horizontal span and
    the cell-height interval overlaps the wall. The z-face mask is empty: flow
    normal to a horizontal face is tangent to this vertical wall. Face masks
    are binary subgrid geometry, not fractional cut-cell areas.
    """
    coordinates = tuple(np.asarray(axis, dtype=float) for axis in (x_m, y_m, z_m))
    names = ("x_m", "y_m", "z_m")
    spacings: list[float] = []
    for name, axis in zip(names, coordinates):
        if axis.ndim != 1 or axis.size < 2 or not np.all(np.isfinite(axis)):
            raise ValueError(f"{name} must be a finite 1-D cell-centre array")
        differences = np.diff(axis)
        if np.any(differences <= 0.0) or not np.allclose(
            differences, differences[0], rtol=1.0e-10, atol=1.0e-12,
        ):
            raise ValueError(f"{name} must be strictly increasing and uniform")
        spacings.append(float(differences[0]))
    x, y, z = coordinates
    dx, dy, dz = spacings
    center = np.asarray(center_xy_m, dtype=float)
    if center.shape != (2,) or not np.all(np.isfinite(center)):
        raise ValueError("center_xy_m must contain exactly two finite coordinates")
    if not all(math.isfinite(float(value)) for value in (
        radius_m, z_min_m, z_max_m, arc_start_deg, arc_end_deg,
    )):
        raise ValueError("arc dimensions and angles must be finite")
    if radius_m <= 0.0 or z_min_m > z_max_m:
        raise ValueError("arc radius must be positive and vertical bounds ordered")
    span = (arc_end_deg - arc_start_deg) % 360.0
    if span <= _GEOMETRY_TOLERANCE:
        raise ValueError("arc angular span must lie in (0, 360) degrees")

    def in_arc(angle_degrees: np.ndarray) -> np.ndarray:
        offset = (angle_degrees - arc_start_deg) % 360.0
        return offset <= span + _GEOMETRY_TOLERANCE

    def overlaps_vertical_cells() -> np.ndarray:
        cell_lower = z - 0.5 * dz
        cell_upper = z + 0.5 * dz
        return (
            (cell_upper > z_min_m + _GEOMETRY_TOLERANCE)
            & (cell_lower < z_max_m - _GEOMETRY_TOLERANCE)
        )

    x_faces = x[0] - 0.5 * dx + np.arange(x.size + 1) * dx
    x_offset = x_faces - center[0]
    x_radial_square = radius_m**2 - x_offset**2
    x_intersects_circle = x_radial_square >= -_GEOMETRY_TOLERANCE
    x_root = np.sqrt(np.maximum(x_radial_square, 0.0))
    y_candidates = np.stack((center[1] - x_root, center[1] + x_root), axis=1)
    x_angles = np.degrees(np.arctan2(
        y_candidates - center[1], x_offset[:, None],
    )) % 360.0
    x_candidate_valid = x_intersects_circle[:, None] & in_arc(x_angles)
    x_horizontal = np.any(
        (np.abs(y[None, :, None] - y_candidates[:, None, :])
         <= 0.5 * dy + _GEOMETRY_TOLERANCE)
        & x_candidate_valid[:, None, :],
        axis=2,
    ).T
    x_blocked = overlaps_vertical_cells()[:, None, None] & x_horizontal[None, :, :]

    y_faces = y[0] - 0.5 * dy + np.arange(y.size + 1) * dy
    y_offset = y_faces - center[1]
    y_radial_square = radius_m**2 - y_offset**2
    y_intersects_circle = y_radial_square >= -_GEOMETRY_TOLERANCE
    y_root = np.sqrt(np.maximum(y_radial_square, 0.0))
    x_candidates = np.stack((center[0] - y_root, center[0] + y_root), axis=1)
    y_angles = np.degrees(np.arctan2(
        y_offset[:, None], x_candidates - center[0],
    )) % 360.0
    y_candidate_valid = y_intersects_circle[:, None] & in_arc(y_angles)
    y_horizontal = np.any(
        (np.abs(x[None, :, None] - x_candidates[:, None, :])
         <= 0.5 * dx + _GEOMETRY_TOLERANCE)
        & y_candidate_valid[:, None, :],
        axis=2,
    )
    y_blocked = overlaps_vertical_cells()[:, None, None] & y_horizontal[None, :, :]
    z_blocked = np.zeros((z.size + 1, y.size, x.size), dtype=bool)
    return x_blocked, y_blocked, z_blocked


def circular_arc_surface_cell_metrics(
    x_m: Sequence[float],
    y_m: Sequence[float],
    z_m: Sequence[float],
    *,
    center_xy_m: Sequence[float],
    radius_m: float,
    z_min_m: float,
    z_max_m: float,
    arc_start_deg: float,
    arc_end_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Integrate a vertical circular-arc surface over Cartesian cells.

    Returns ``(area, area_normal_x, area_normal_y, wall_distance)`` arrays with shape
    ``(nz, ny, nx)``.  The two normal arrays are *area-weighted* outward
    normal components, not unit normals.  Arc intervals are split at every
    Cartesian cell edge before integration, so the summed area is invariant
    to grid alignment up to floating-point roundoff.  This supplies geometric
    measures for embedded wall-stress models; it does not itself impose a
    momentum boundary condition.
    """
    coordinates = tuple(np.asarray(axis, dtype=float) for axis in (x_m, y_m, z_m))
    if any(
        axis.ndim != 1 or axis.size < 2 or not np.all(np.isfinite(axis))
        for axis in coordinates
    ):
        raise ValueError("cell-centre coordinates must be finite 1-D arrays")
    spacings = tuple(np.diff(axis) for axis in coordinates)
    if any(
        np.any(delta <= 0.0)
        or not np.allclose(delta, delta[0], rtol=1.0e-10, atol=1.0e-12)
        for delta in spacings
    ):
        raise ValueError("cell-centre coordinates must be uniformly increasing")
    x, y, z = coordinates
    dx, dy, dz = (float(delta[0]) for delta in spacings)
    center = np.asarray(center_xy_m, dtype=float)
    values = (radius_m, z_min_m, z_max_m, arc_start_deg, arc_end_deg)
    if center.shape != (2,) or not np.all(np.isfinite(center)):
        raise ValueError("center_xy_m must contain exactly two finite coordinates")
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("arc dimensions and angles must be finite")
    if radius_m <= 0.0 or z_min_m > z_max_m:
        raise ValueError("radius must be positive and z_min_m <= z_max_m")
    span_deg = (float(arc_end_deg) - float(arc_start_deg)) % 360.0
    if span_deg <= _GEOMETRY_TOLERANCE:
        raise ValueError("arc angular span must lie in (0, 360) degrees")

    theta_start = math.radians(float(arc_start_deg) % 360.0)
    span = math.radians(span_deg)
    breaks = [0.0, span]
    x_edges = x[0] - 0.5 * dx + np.arange(x.size + 1) * dx
    y_edges = y[0] - 0.5 * dy + np.arange(y.size + 1) * dy

    def add_angle(angle_rad: float) -> None:
        offset = (angle_rad - theta_start) % (2.0 * math.pi)
        if _GEOMETRY_TOLERANCE < offset < span - _GEOMETRY_TOLERANCE:
            breaks.append(offset)

    for edge in x_edges:
        cosine = (float(edge) - center[0]) / radius_m
        if -1.0 <= cosine <= 1.0:
            angle = math.acos(cosine)
            add_angle(angle)
            add_angle(-angle)
    for edge in y_edges:
        sine = (float(edge) - center[1]) / radius_m
        if -1.0 <= sine <= 1.0:
            angle = math.asin(sine)
            add_angle(angle)
            add_angle(math.pi - angle)
    breaks.sort()
    unique_breaks = [breaks[0]]
    for value in breaks[1:]:
        if value - unique_breaks[-1] > _GEOMETRY_TOLERANCE:
            unique_breaks.append(value)

    area = np.zeros((z.size, y.size, x.size), dtype=float)
    area_normal_x = np.zeros_like(area)
    area_normal_y = np.zeros_like(area)
    wall_distance = np.full_like(area, np.inf)
    z_edges = z[0] - 0.5 * dz + np.arange(z.size + 1) * dz
    z_overlap = np.maximum(
        0.0,
        np.minimum(z_edges[1:], z_max_m) - np.maximum(z_edges[:-1], z_min_m),
    )
    for lower, upper in zip(unique_breaks[:-1], unique_breaks[1:]):
        midpoint = theta_start + 0.5 * (lower + upper)
        px = center[0] + radius_m * math.cos(midpoint)
        py = center[1] + radius_m * math.sin(midpoint)
        i = int(np.searchsorted(x_edges, px, side="right") - 1)
        j = int(np.searchsorted(y_edges, py, side="right") - 1)
        if not (0 <= i < x.size and 0 <= j < y.size):
            continue
        delta = upper - lower
        nx_integral = math.sin(theta_start + upper) - math.sin(theta_start + lower)
        ny_integral = math.cos(theta_start + lower) - math.cos(theta_start + upper)
        theta_cell = math.atan2(y[j] - center[1], x[i] - center[0])
        cell_offset = (theta_cell - theta_start) % (2.0 * math.pi)
        cell_radius = math.hypot(x[i] - center[0], y[j] - center[1])
        if lower <= cell_offset <= upper:
            horizontal_distance = abs(cell_radius - radius_m)
        else:
            endpoint_a = theta_start + lower
            endpoint_b = theta_start + upper
            horizontal_distance = min(
                math.hypot(
                    x[i] - (center[0] + radius_m * math.cos(endpoint_a)),
                    y[j] - (center[1] + radius_m * math.sin(endpoint_a)),
                ),
                math.hypot(
                    x[i] - (center[0] + radius_m * math.cos(endpoint_b)),
                    y[j] - (center[1] + radius_m * math.sin(endpoint_b)),
                ),
            )
        for k, height in enumerate(z_overlap):
            if height <= 0.0:
                continue
            area[k, j, i] += radius_m * delta * height
            area_normal_x[k, j, i] += radius_m * nx_integral * height
            area_normal_y[k, j, i] += radius_m * ny_integral * height
            vertical_distance = max(z_min_m - z[k], 0.0, z[k] - z_max_m)
            wall_distance[k, j, i] = min(
                wall_distance[k, j, i],
                math.hypot(horizontal_distance, vertical_distance),
            )
    return area, area_normal_x, area_normal_y, wall_distance


def circular_arc_surface_link_fractions(
    x_m: Sequence[float],
    y_m: Sequence[float],
    z_m: Sequence[float],
    directions: Sequence[Sequence[int]],
    *,
    center_xy_m: Sequence[float],
    radius_m: float,
    z_min_m: float,
    z_max_m: float,
    arc_start_deg: float,
    arc_end_deg: float,
) -> np.ndarray:
    """Return first wall-hit fraction for each D3Q19 link, NaN when clear.

    Fractions lie in [0, 1] measured from the starting cell centre along the
    link. Contacts are exact for the declared vertical arc, including diagonal
    D3Q19 links; if a segment has two contacts, the first is returned.
    Coordinates must be uniform and cubic. Geometry is kept separate from the
    selected kinetic boundary rule.
    """
    coordinates = tuple(np.asarray(axis, dtype=float) for axis in (x_m, y_m, z_m))
    names = ("x_m", "y_m", "z_m")
    spacings: list[float] = []
    for name, axis in zip(names, coordinates):
        if axis.ndim != 1 or axis.size < 2 or not np.all(np.isfinite(axis)):
            raise ValueError(f"{name} must be a finite 1-D cell-centre array")
        differences = np.diff(axis)
        if np.any(differences <= 0.0) or not np.allclose(
            differences, differences[0], rtol=1.0e-10, atol=1.0e-12,
        ):
            raise ValueError(f"{name} must be strictly increasing and uniform")
        spacings.append(float(differences[0]))
    if not np.allclose(spacings, spacings[0], rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("D3Q19 link masks require uniform cubic cells")
    direction_array = np.asarray(directions)
    if (
        direction_array.ndim != 2 or direction_array.shape[1] != 3
        or not np.issubdtype(direction_array.dtype, np.integer)
        or np.any(np.abs(direction_array) > 1)
        or np.any(np.all(direction_array == 0, axis=1))
    ):
        raise ValueError("directions must be non-zero D3Q19 integer links")
    center = np.asarray(center_xy_m, dtype=float)
    if center.shape != (2,) or not np.all(np.isfinite(center)):
        raise ValueError("center_xy_m must contain exactly two finite coordinates")
    if not all(math.isfinite(float(value)) for value in (
        radius_m, z_min_m, z_max_m, arc_start_deg, arc_end_deg,
    )):
        raise ValueError("arc dimensions and angles must be finite")
    if radius_m <= 0.0 or z_min_m > z_max_m:
        raise ValueError("arc radius must be positive and vertical bounds ordered")
    span = (arc_end_deg - arc_start_deg) % 360.0
    if span <= _GEOMETRY_TOLERANCE:
        raise ValueError("arc angular span must lie in (0, 360) degrees")

    x, y, z = coordinates
    xx, yy, zz = np.meshgrid(x, y, z, indexing="xy")
    # meshgrid(indexing='xy') is (ny,nx,nz); move z to the leading axis to
    # match the flow solver's (nz,ny,nx) field layout.
    start_x = np.moveaxis(xx, -1, 0) - center[0]
    start_y = np.moveaxis(yy, -1, 0) - center[1]
    start_z = np.moveaxis(zz, -1, 0)
    result = np.full((len(direction_array), z.size, y.size, x.size), np.nan)
    for index, (cx, cy, cz) in enumerate(direction_array):
        if cx == 0 and cy == 0:
            continue  # a vertical link is tangent to a vertical cylinder
        delta_x, delta_y, delta_z = cx * spacings[0], cy * spacings[0], cz * spacings[0]
        quadratic = delta_x**2 + delta_y**2
        linear = 2.0 * (start_x * delta_x + start_y * delta_y)
        constant = start_x**2 + start_y**2 - radius_m**2
        discriminant = linear**2 - 4.0 * quadratic * constant
        tolerance = _GEOMETRY_TOLERANCE * np.maximum.reduce((
            np.ones_like(discriminant), linear**2,
            np.abs(4.0 * quadratic * constant),
            np.full_like(discriminant, quadratic * radius_m**2),
        ))
        root = np.sqrt(np.maximum(discriminant, 0.0))
        fractions = (
            (-linear - root) / (2.0 * quadratic),
            (-linear + root) / (2.0 * quadratic),
        )
        for fraction in fractions:
            valid = (
                (discriminant >= -tolerance)
                & (fraction >= -_GEOMETRY_TOLERANCE)
                & (fraction <= 1.0 + _GEOMETRY_TOLERANCE)
            )
            hit_x = start_x + fraction * delta_x
            hit_y = start_y + fraction * delta_y
            hit_z = start_z + fraction * delta_z
            angle = np.degrees(np.arctan2(hit_y, hit_x)) % 360.0
            offset = (angle - arc_start_deg) % 360.0
            valid &= (
                (hit_z >= z_min_m - _GEOMETRY_TOLERANCE)
                & (hit_z <= z_max_m + _GEOMETRY_TOLERANCE)
                & (offset <= span + _GEOMETRY_TOLERANCE)
            )
            clipped_fraction = np.clip(fraction, 0.0, 1.0)
            result[index] = np.where(
                valid, np.fmin(result[index], clipped_fraction), result[index],
            )
    return result


def circular_arc_surface_link_masks(*args, **kwargs) -> np.ndarray:
    """Return a boolean view of :func:`circular_arc_surface_link_fractions`."""
    return np.isfinite(circular_arc_surface_link_fractions(*args, **kwargs))


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

    def horizontal_corners(self) -> tuple[tuple[float, float], ...]:
        """Return the four global horizontal footprint corners."""
        return (
            (self.x_min_m, self.y_min_m), (self.x_min_m, self.y_max_m),
            (self.x_max_m, self.y_min_m), (self.x_max_m, self.y_max_m),
        )

    @property
    def bounding_box_m(self) -> tuple[float, float, float, float, float, float]:
        return (
            self.x_min_m, self.x_max_m, self.y_min_m, self.y_max_m,
            self.z_min_m, self.z_max_m,
        )

    def contains_point(self, point_m: Sequence[float]) -> bool:
        point = _point3(point_m, name="point")
        return bool(
            self.x_min_m <= point[0] <= self.x_max_m
            and self.y_min_m <= point[1] <= self.y_max_m
            and self.z_min_m <= point[2] <= self.z_max_m
        )

    def voxel_mask(
        self, x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray,
    ) -> np.ndarray:
        """Cell-centre, half-open Cartesian solid mask."""
        return (
            (x_m[None, None, :] >= self.x_min_m)
            & (x_m[None, None, :] < self.x_max_m)
            & (y_m[None, :, None] >= self.y_min_m)
            & (y_m[None, :, None] < self.y_max_m)
            & (z_m[:, None, None] >= self.z_min_m)
            & (z_m[:, None, None] < self.z_max_m)
        )


@dataclass(frozen=True)
class CircularArcWall:
    """Finite-height annular-sector wall represented on a Cartesian grid.

    Angles use the mathematical convention: degrees counter-clockwise from
    global ``+x``.  ``arc_end_deg`` is reached counter-clockwise from
    ``arc_start_deg`` and the span must lie in ``(0, 360)``.  The declared
    thickness is radial. A positive thickness defines an annular solid band;
    zero thickness denotes an ideal impermeable surface and has no voxel
    volume. The resolved solver must select its matching surface-boundary mode
    for a zero-thickness wall. A one-cell value may be used as an explicitly
    named numerical regularisation when physical membrane thickness is
    unavailable, but that choice requires its own grid-sensitivity evidence.
    """

    center_m: tuple[float, float]
    radius_m: float
    thickness_m: float
    z_min_m: float
    z_max_m: float
    arc_start_deg: float
    arc_end_deg: float
    label: str = "circular_arc_wall"

    def __post_init__(self) -> None:
        if len(self.center_m) != 2 or not all(
            math.isfinite(float(value)) for value in self.center_m
        ):
            raise ValueError("arc-wall center_m must contain two finite values")
        for name in (
            "radius_m", "thickness_m", "z_min_m", "z_max_m",
            "arc_start_deg", "arc_end_deg",
        ):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"arc-wall {name} must be finite")
        if self.radius_m <= 0.0 or self.thickness_m < 0.0:
            raise ValueError("arc-wall radius must be positive and thickness non-negative")
        if self.thickness_m >= 2.0 * self.radius_m:
            raise ValueError("arc-wall thickness must be less than its diameter")
        if self.z_min_m >= self.z_max_m:
            raise ValueError("arc-wall requires positive vertical extent")
        span = (self.arc_end_deg - self.arc_start_deg) % 360.0
        if span <= 1.0e-12:
            raise ValueError("arc-wall angular span must lie in (0, 360) degrees")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("arc-wall label must be a non-empty string")
        object.__setattr__(
            self, "center_m", (float(self.center_m[0]), float(self.center_m[1]))
        )
        object.__setattr__(self, "radius_m", float(self.radius_m))
        object.__setattr__(self, "thickness_m", float(self.thickness_m))
        object.__setattr__(self, "z_min_m", float(self.z_min_m))
        object.__setattr__(self, "z_max_m", float(self.z_max_m))
        object.__setattr__(self, "arc_start_deg", float(self.arc_start_deg) % 360.0)
        object.__setattr__(self, "arc_end_deg", float(self.arc_end_deg) % 360.0)

    @property
    def angular_span_deg(self) -> float:
        return (self.arc_end_deg - self.arc_start_deg) % 360.0

    @property
    def inner_radius_m(self) -> float:
        return self.radius_m - 0.5 * self.thickness_m

    @property
    def outer_radius_m(self) -> float:
        return self.radius_m + 0.5 * self.thickness_m

    def _angle_in_arc(self, angle_deg: np.ndarray | float) -> np.ndarray | bool:
        delta = (np.asarray(angle_deg) - self.arc_start_deg) % 360.0
        result = delta <= self.angular_span_deg + 1.0e-12
        return bool(result) if np.ndim(result) == 0 else result

    @property
    def bounding_box_m(self) -> tuple[float, float, float, float, float, float]:
        candidate_angles = [self.arc_start_deg, self.arc_end_deg]
        candidate_angles.extend(
            angle for angle in (0.0, 90.0, 180.0, 270.0)
            if self._angle_in_arc(angle)
        )
        points = []
        for angle in candidate_angles:
            radians = math.radians(angle)
            for radius in (self.inner_radius_m, self.outer_radius_m):
                points.append((
                    self.center_m[0] + radius * math.cos(radians),
                    self.center_m[1] + radius * math.sin(radians),
                ))
        return (
            min(point[0] for point in points), max(point[0] for point in points),
            min(point[1] for point in points), max(point[1] for point in points),
            self.z_min_m, self.z_max_m,
        )

    def contains_point(self, point_m: Sequence[float]) -> bool:
        point = _point3(point_m, name="point")
        if point[2] < self.z_min_m or point[2] > self.z_max_m:
            return False
        dx = point[0] - self.center_m[0]
        dy = point[1] - self.center_m[1]
        radius = math.hypot(dx, dy)
        angle = math.degrees(math.atan2(dy, dx)) % 360.0
        return bool(
            self.inner_radius_m <= radius <= self.outer_radius_m
            and self._angle_in_arc(angle)
        )

    def voxel_mask(
        self, x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray,
    ) -> np.ndarray:
        dx = x_m[None, :] - self.center_m[0]
        dy = y_m[:, None] - self.center_m[1]
        radius = np.hypot(dx, dy)
        angle = np.degrees(np.arctan2(dy, dx)) % 360.0
        horizontal = (
            (radius >= self.inner_radius_m)
            & (radius < self.outer_radius_m)
            & self._angle_in_arc(angle)
        )
        vertical = (
            (z_m >= self.z_min_m) & (z_m < self.z_max_m)
        )
        return vertical[:, None, None] & horizontal[None, :, :]


@dataclass(frozen=True)
class CrosswindCylinder:
    """Finite horizontal cylinder with its axis along the global y direction.

    This geometry directly represents transverse rods or pipe arrays in a
    streamwise-x / crosswind-y / vertical-z site coordinate system. The caps
    are flat at ``y_min_m`` and ``y_max_m``; the circular section is in x-z.
    """

    center_x_m: float
    center_z_m: float
    radius_m: float
    y_min_m: float
    y_max_m: float
    label: str = "crosswind_cylinder"

    def __post_init__(self) -> None:
        values = (
            self.center_x_m, self.center_z_m, self.radius_m,
            self.y_min_m, self.y_max_m,
        )
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("crosswind-cylinder dimensions must be finite")
        if self.radius_m <= 0.0 or self.y_min_m >= self.y_max_m:
            raise ValueError("crosswind cylinder requires positive radius and span")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("crosswind-cylinder label must be non-empty")
        for name in (
            "center_x_m", "center_z_m", "radius_m", "y_min_m", "y_max_m",
        ):
            object.__setattr__(self, name, float(getattr(self, name)))

    @property
    def bounding_box_m(self) -> tuple[float, float, float, float, float, float]:
        return (
            self.center_x_m - self.radius_m,
            self.center_x_m + self.radius_m,
            self.y_min_m,
            self.y_max_m,
            self.center_z_m - self.radius_m,
            self.center_z_m + self.radius_m,
        )

    def contains_point(self, point_m: Sequence[float]) -> bool:
        x, y, z = _point3(point_m, name="point")
        return bool(
            self.y_min_m <= y <= self.y_max_m
            and (x - self.center_x_m)**2 + (z - self.center_z_m)**2
            <= self.radius_m**2
        )

    def voxel_mask(
        self, x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray,
    ) -> np.ndarray:
        radial_square = (
            (x_m[None, :] - self.center_x_m)**2
            + (z_m[:, None] - self.center_z_m)**2
        )
        cross_section = radial_square <= self.radius_m**2
        span = (y_m >= self.y_min_m) & (y_m < self.y_max_m)
        return span[None, :, None] & cross_section[:, None, :]


SolidObstacle3D = AxisAlignedCuboid | CircularArcWall | CrosswindCylinder


def solid_obstacle_bounds_m(
    obstacle: SolidObstacle3D,
) -> tuple[float, float, float, float, float, float]:
    if not isinstance(obstacle, (AxisAlignedCuboid, CircularArcWall, CrosswindCylinder)):
        raise TypeError("unsupported 3-D solid obstacle type")
    return obstacle.bounding_box_m


def solid_obstacle_contains_point(
    obstacle: SolidObstacle3D, point_m: Sequence[float],
) -> bool:
    if not isinstance(obstacle, (AxisAlignedCuboid, CircularArcWall, CrosswindCylinder)):
        raise TypeError("unsupported 3-D solid obstacle type")
    return obstacle.contains_point(point_m)


def solid_obstacle_voxel_mask(
    obstacle: SolidObstacle3D,
    x_m: np.ndarray,
    y_m: np.ndarray,
    z_m: np.ndarray,
) -> np.ndarray:
    if not isinstance(obstacle, (AxisAlignedCuboid, CircularArcWall, CrosswindCylinder)):
        raise TypeError("unsupported 3-D solid obstacle type")
    return obstacle.voxel_mask(x_m, y_m, z_m)


@dataclass(frozen=True)
class OrientedCuboid:
    """A rectangular global-footprint cuboid with a declared long-axis bearing.

    ``long_axis_bearing_deg`` is a mathematical *to* direction measured
    counter-clockwise from global ``+x``. This convention is deliberately
    separate from meteorological wind-*from* bearing. It preserves an actual
    facility orientation before projection into the local wind plane; it does
    not introduce a 3-D wake or bypass closure.
    """

    center_m: tuple[float, float]
    length_m: float
    width_m: float
    z_min_m: float
    z_max_m: float
    long_axis_bearing_deg: float
    label: str = "oriented_cuboid"

    def __post_init__(self) -> None:
        if len(self.center_m) != 2 or not all(math.isfinite(float(value)) for value in self.center_m):
            raise ValueError("oriented cuboid center_m must contain two finite coordinates")
        values = {
            "length_m": self.length_m,
            "width_m": self.width_m,
            "z_min_m": self.z_min_m,
            "z_max_m": self.z_max_m,
            "long_axis_bearing_deg": self.long_axis_bearing_deg,
        }
        if not all(math.isfinite(float(value)) for value in values.values()):
            raise ValueError("oriented cuboid dimensions and bearing must be finite")
        if self.length_m <= 0.0 or self.width_m <= 0.0 or self.z_min_m >= self.z_max_m:
            raise ValueError("an oriented cuboid requires positive extent on every axis")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("obstacle label must be a non-empty string")
        object.__setattr__(self, "center_m", (float(self.center_m[0]), float(self.center_m[1])))
        object.__setattr__(self, "long_axis_bearing_deg", float(self.long_axis_bearing_deg) % 360.0)

    @property
    def long_axis_direction_rad(self) -> float:
        return math.radians(self.long_axis_bearing_deg)

    def horizontal_corners(self) -> tuple[tuple[float, float], ...]:
        """Return the four actual global footprint corners in stable order."""
        ex, ey = math.cos(self.long_axis_direction_rad), math.sin(self.long_axis_direction_rad)
        nx, ny = -ey, ex
        half_length, half_width = 0.5 * self.length_m, 0.5 * self.width_m
        return tuple(
            (
                self.center_m[0] + along * ex + across * nx,
                self.center_m[1] + along * ey + across * ny,
            )
            for along, across in (
                (-half_length, -half_width), (-half_length, half_width),
                (half_length, -half_width), (half_length, half_width),
            )
        )

    def segment_intersection(
        self, start_m: Sequence[float], end_m: Sequence[float]
    ) -> SegmentIntersection | None:
        """Intersect a global segment after an exact horizontal rotation."""
        start, end = _point3(start_m, name="segment start"), _point3(end_m, name="segment end")
        ex, ey = math.cos(self.long_axis_direction_rad), math.sin(self.long_axis_direction_rad)
        nx, ny = -ey, ex

        def local(point: np.ndarray) -> tuple[float, float, float]:
            dx, dy = point[0] - self.center_m[0], point[1] - self.center_m[1]
            return dx * ex + dy * ey, dx * nx + dy * ny, float(point[2])

        hit = _segment_aabb_intersection(
            local(start), local(end),
            (-0.5 * self.length_m, -0.5 * self.width_m, self.z_min_m),
            (0.5 * self.length_m, 0.5 * self.width_m, self.z_max_m),
        )
        if hit is None:
            return None
        delta = end - start
        entry = start + hit.entry_fraction * delta
        exit_ = start + hit.exit_fraction * delta
        return SegmentIntersection(
            hit.entry_fraction, hit.exit_fraction,
            tuple(float(value) for value in entry), tuple(float(value) for value in exit_),
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
    if isinstance(obstacle, OrientedCuboid):
        return "oriented_cuboid", False, None
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
