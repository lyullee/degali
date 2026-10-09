"""Directional site-geometry bridge for the local semi-FV obstacle plane."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import Mapping

from .field_contracts import BoundedValue, CircularBoundedValue
from .semi_fv_obstacle import RectangularObstacle2D
from .site_geometry import AxisAlignedCuboid, OrientedCuboid, WindFrame


_GlobalObstacle = AxisAlignedCuboid | OrientedCuboid
_GeometryBound = BoundedValue | CircularBoundedValue


@dataclass(frozen=True)
class FieldObstacleGeometryUncertainty:
    """Explicit deterministic bounds for one declared global obstacle.

    The nominal geometry remains the obstacle used by a single screening run.
    Every non-exact bound is expanded into Cartesian lower/upper geometry
    corners by the field envelope; no probability, wake correction, or
    lateral bypass is inferred from a dimensional tolerance.
    """

    obstacle: _GlobalObstacle
    fields: tuple[tuple[str, _GeometryBound], ...]
    evidence_id: str

    def __post_init__(self) -> None:
        if not isinstance(self.obstacle, (AxisAlignedCuboid, OrientedCuboid)):
            raise TypeError("obstacle must be an AxisAlignedCuboid or OrientedCuboid")
        if (
            not isinstance(self.evidence_id, str)
            or not self.evidence_id.strip()
            or self.evidence_id.strip().lower() == "unspecified"
        ):
            raise ValueError("obstacle geometry uncertainty requires an explicit evidence_id")
        expected = self._expected_field_names()
        if len(self.fields) != len(expected):
            raise ValueError(
                f"obstacle geometry uncertainty for {self.obstacle.label!r} requires "
                f"exactly {len(expected)} fields"
            )
        names = tuple(name for name, _value in self.fields)
        if names != expected or len(set(names)) != len(names):
            raise ValueError(
                "obstacle geometry uncertainty fields must use the ordered names "
                + ", ".join(expected)
            )
        for name, bound in self.fields:
            if not isinstance(bound, (BoundedValue, CircularBoundedValue)):
                raise TypeError(f"{name} must be a BoundedValue or CircularBoundedValue")
            if (
                not isinstance(bound.source, str)
                or not bound.source.strip()
                or bound.source.strip().lower() == "unspecified"
            ):
                raise ValueError(f"{name} requires an explicit source")
            expected_unit = "deg" if name == "long_axis_bearing_deg" else "m"
            if bound.unit != expected_unit:
                raise ValueError(f"{name} must use unit {expected_unit!r}")
            if name == "long_axis_bearing_deg":
                if not isinstance(bound, (BoundedValue, CircularBoundedValue)):
                    raise TypeError(f"{name} must be an angular bounded value")
            elif not isinstance(bound, BoundedValue):
                raise TypeError(f"{name} must use a linear BoundedValue")
        nominal = self._nominal_values()
        for name, bound in self.fields:
            expected_nominal = nominal[name]
            if name == "long_axis_bearing_deg":
                difference = (float(bound.nominal) - expected_nominal) % 360.0
                if not (math.isclose(difference, 0.0, abs_tol=1.0e-9)
                        or math.isclose(difference, 360.0, abs_tol=1.0e-9)):
                    raise ValueError(f"{name}.nominal must match the declared obstacle")
            elif not math.isclose(
                float(bound.nominal), expected_nominal,
                rel_tol=1.0e-9, abs_tol=1.0e-12,
            ):
                raise ValueError(f"{name}.nominal must match the declared obstacle")
        # Validate every declared lower/upper combination before a workflow
        # can advertise an envelope. Invalid dimensions are rejected here,
        # rather than becoming a blocked transport corner later.
        choices = tuple(
            ((bound.nominal,) if bound.is_exact else bound.corners())
            for _name, bound in self.fields
        )
        for selected in product(*choices):
            self._build(dict(zip(expected, selected)))

    def _expected_field_names(self) -> tuple[str, ...]:
        if isinstance(self.obstacle, AxisAlignedCuboid):
            return (
                "x_min_m", "x_max_m", "y_min_m", "y_max_m", "z_min_m", "z_max_m",
            )
        return (
            "center_x_m", "center_y_m", "length_m", "width_m",
            "z_min_m", "z_max_m", "long_axis_bearing_deg",
        )

    def _nominal_values(self) -> dict[str, float]:
        if isinstance(self.obstacle, AxisAlignedCuboid):
            return {
                "x_min_m": self.obstacle.x_min_m,
                "x_max_m": self.obstacle.x_max_m,
                "y_min_m": self.obstacle.y_min_m,
                "y_max_m": self.obstacle.y_max_m,
                "z_min_m": self.obstacle.z_min_m,
                "z_max_m": self.obstacle.z_max_m,
            }
        return {
            "center_x_m": self.obstacle.center_m[0],
            "center_y_m": self.obstacle.center_m[1],
            "length_m": self.obstacle.length_m,
            "width_m": self.obstacle.width_m,
            "z_min_m": self.obstacle.z_min_m,
            "z_max_m": self.obstacle.z_max_m,
            "long_axis_bearing_deg": self.obstacle.long_axis_bearing_deg,
        }

    def _build(self, values: Mapping[str, float]) -> _GlobalObstacle:
        if isinstance(self.obstacle, AxisAlignedCuboid):
            return AxisAlignedCuboid(
                float(values["x_min_m"]), float(values["x_max_m"]),
                float(values["y_min_m"]), float(values["y_max_m"]),
                float(values["z_min_m"]), float(values["z_max_m"]),
                self.obstacle.label,
            )
        return OrientedCuboid(
            (float(values["center_x_m"]), float(values["center_y_m"])),
            float(values["length_m"]), float(values["width_m"]),
            float(values["z_min_m"]), float(values["z_max_m"]),
            float(values["long_axis_bearing_deg"]), self.obstacle.label,
        )

    def uncertainty_fields(self) -> dict[str, _GeometryBound]:
        return {
            f"obstacle.{self.obstacle.label}.{name}": bound
            for name, bound in self.fields
        }

    def corner_cases(
        self,
    ) -> tuple[tuple[_GlobalObstacle, tuple[tuple[str, float], ...]], ...]:
        names = self._expected_field_names()
        choices = tuple(
            ((bound.nominal,) if bound.is_exact else bound.corners())
            for _name, bound in self.fields
        )
        cases = []
        for selected in product(*choices):
            raw = dict(zip(names, (float(value) for value in selected)))
            cases.append((
                self._build(raw),
                tuple(
                    (f"obstacle.{self.obstacle.label}.{name}", float(raw[name]))
                    for name in names
                ),
            ))
        return tuple(cases)

    def at_corner(self, values: Mapping[str, float]) -> _GlobalObstacle:
        """Fix every declared geometry bound at one deterministic corner."""
        fields = self.uncertainty_fields()
        if set(values) != set(fields):
            missing = sorted(set(fields) - set(values))
            extra = sorted(set(values) - set(fields))
            raise ValueError(
                "obstacle geometry corner keys must match uncertainty fields; "
                f"missing={missing}, extra={extra}"
            )
        prefix = f"obstacle.{self.obstacle.label}."
        raw = {
            name[len(prefix):]: float(value)
            for name, value in values.items()
        }
        return self._build(raw)

    def as_record(self) -> dict[str, object]:
        return {
            "geometry_type": type(self.obstacle).__name__,
            "label": self.obstacle.label,
            "evidence_id": self.evidence_id,
            "uncertainty": {
                name: bound.as_dict() for name, bound in self.fields
            },
            "scope": (
                "deterministic obstacle-geometry corners; no probability, "
                "lateral-bypass or 3-D wake closure"
            ),
        }


@dataclass(frozen=True)
class WindPlaneObstacleProjection:
    """A global cuboid projected into the wind-aligned x-z transport plane.

    The semi-FV operator is two dimensional. A finite-width cuboid that
    crosses the wind-plane centreline is represented as a laterally unbounded
    obstacle in that plane, so the result is necessarily ``conditional``. A
    cuboid wholly away from the centreline is not inserted into the local
    plane; that is not evidence that its three-dimensional wake is harmless.
    """

    obstacle: RectangularObstacle2D | None
    downwind_bounds_m: tuple[float, float]
    crosswind_bounds_m: tuple[float, float]
    status: str
    reasons: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.status not in {"accepted", "conditional", "blocked"}:
            raise ValueError("projection status must be accepted, conditional, or blocked")
        if self.status == "blocked" and not self.reasons:
            raise ValueError("blocked projection requires a reason")

    @property
    def centreline_intersects(self) -> bool:
        lower, upper = self.crosswind_bounds_m
        return lower <= 0.0 <= upper


def project_cuboid_to_wind_plane(
    cuboid: AxisAlignedCuboid | OrientedCuboid,
    frame: WindFrame,
) -> WindPlaneObstacleProjection:
    """Project one global cuboid using an explicit wind-*to* frame.

    A building overlapping the source plane cannot be represented by a local
    downwind inlet and is blocked. An upstream building is outside this local
    downwind operator and is reported but not inserted. The function never
    derives a wake correction from geometry alone.
    """
    local = tuple(frame.local(x_m, y_m) for x_m, y_m in cuboid.horizontal_corners())
    downwind = tuple(value[0] for value in local)
    crosswind = tuple(value[1] for value in local)
    s_min, s_max = min(downwind), max(downwind)
    n_min, n_max = min(crosswind), max(crosswind)
    bounds_s, bounds_n = (s_min, s_max), (n_min, n_max)

    if s_min <= 0.0 <= s_max:
        return WindPlaneObstacleProjection(
            None, bounds_s, bounds_n, "blocked",
            reasons=(
                "obstacle overlaps the source wind plane; local downwind inlet is undefined",
            ),
        )
    if s_max <= 0.0:
        return WindPlaneObstacleProjection(
            None, bounds_s, bounds_n, "accepted",
            warnings=(
                "obstacle lies fully upwind of the local semi-FV domain and was not inserted",
            ),
        )
    if not n_min <= 0.0 <= n_max:
        return WindPlaneObstacleProjection(
            None, bounds_s, bounds_n, "accepted",
            warnings=(
                "obstacle does not intersect the wind-plane centreline; "
                "two-dimensional calculation does not resolve its lateral wake",
            ),
        )
    # The local semi-FV plane uses z=0 as its lower boundary.  Do not pass a
    # cuboid that crosses that datum to RectangularObstacle2D: clipping its
    # negative base would silently change the declared geometry and could
    # otherwise escape the workflow as an uncaught ValueError.
    if cuboid.z_min_m < 0.0 < cuboid.z_max_m:
        return WindPlaneObstacleProjection(
            None, bounds_s, bounds_n, "blocked",
            reasons=(
                "obstacle vertical bounds cross the local z=0 datum; "
                "declare a compatible vertical reference before semi-FV transport",
            ),
        )
    if cuboid.z_max_m <= 0.0:
        return WindPlaneObstacleProjection(
            None, bounds_s, bounds_n, "accepted",
            warnings=(
                "obstacle lies at or below the local z=0 datum and was not inserted "
                "into the atmospheric semi-FV plane",
            ),
        )
    obstacle = RectangularObstacle2D(
        distance_m=0.5 * (s_min + s_max),
        width_m=s_max - s_min,
        height_m=cuboid.z_max_m - cuboid.z_min_m,
        base_height_m=cuboid.z_min_m,
        label=cuboid.label,
    )
    return WindPlaneObstacleProjection(
        obstacle, bounds_s, bounds_n, "conditional",
        warnings=(
            "finite crosswind building width is represented as a full wind-plane obstacle",
            "no three-dimensional lateral-bypass or wake-turbulence correction is applied",
        ),
    )


__all__ = [
    "FieldObstacleGeometryUncertainty",
    "WindPlaneObstacleProjection", "project_cuboid_to_wind_plane",
]
