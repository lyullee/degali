import math

import pytest

from degali.addons.field_geometry import project_cuboid_to_wind_plane
from degali.addons.site_geometry import AxisAlignedCuboid, OrientedCuboid, WindFrame


def test_global_cuboid_is_projected_with_an_explicit_rotated_wind_direction():
    cuboid = AxisAlignedCuboid(-2.0, 2.0, 8.0, 10.0, 0.0, 4.0, "tank")
    projection = project_cuboid_to_wind_plane(
        cuboid, WindFrame(direction_rad=math.pi / 2.0)
    )
    assert projection.status == "conditional"
    assert projection.centreline_intersects
    assert projection.obstacle is not None
    assert projection.obstacle.distance_m == pytest.approx(9.0)
    assert projection.obstacle.width_m == pytest.approx(2.0)
    assert projection.obstacle.height_m == pytest.approx(4.0)


def test_off_centre_and_source_overlapping_cuboids_fail_safely_by_scope():
    frame = WindFrame()
    off_centre = project_cuboid_to_wind_plane(
        AxisAlignedCuboid(5.0, 7.0, 3.0, 4.0, 0.0, 3.0), frame
    )
    assert off_centre.status == "accepted"
    assert off_centre.obstacle is None
    assert "lateral wake" in off_centre.warnings[0]

    source_overlap = project_cuboid_to_wind_plane(
        AxisAlignedCuboid(-1.0, 1.0, -1.0, 1.0, 0.0, 3.0), frame
    )
    assert source_overlap.status == "blocked"
    assert source_overlap.obstacle is None
    assert "source wind plane" in source_overlap.reasons[0]


def test_oriented_cuboid_projects_its_rotated_footprint_into_the_wind_plane():
    cuboid = OrientedCuboid(
        center_m=(8.0, 0.0), length_m=4.0, width_m=2.0,
        z_min_m=0.0, z_max_m=3.0, long_axis_bearing_deg=90.0,
        label="rotated-skid",
    )
    projection = project_cuboid_to_wind_plane(cuboid, WindFrame())

    assert projection.status == "conditional"
    assert projection.centreline_intersects
    assert projection.downwind_bounds_m == pytest.approx((7.0, 9.0))
    assert projection.obstacle is not None
    assert projection.obstacle.distance_m == pytest.approx(8.0)
    assert projection.obstacle.width_m == pytest.approx(2.0)


def test_vertical_datum_crossing_is_blocked_instead_of_escaping_as_an_exception():
    projection = project_cuboid_to_wind_plane(
        AxisAlignedCuboid(5.0, 7.0, -1.0, 1.0, -1.0, 2.0, "buried-base"),
        WindFrame(),
    )
    assert projection.status == "blocked"
    assert projection.obstacle is None
    assert "z=0 datum" in projection.reasons[0]


def test_geometry_below_vertical_datum_is_reported_without_inserting_a_wall():
    projection = project_cuboid_to_wind_plane(
        AxisAlignedCuboid(5.0, 7.0, -1.0, 1.0, -3.0, -1.0, "buried-base"),
        WindFrame(),
    )
    assert projection.status == "accepted"
    assert projection.obstacle is None
    assert "below the local z=0 datum" in projection.warnings[0]
