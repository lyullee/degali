import math
from types import SimpleNamespace

import pytest

from degali.addons.site_geometry import (
    AxisAlignedCuboid,
    OrientedCuboid,
    TransverseWall,
    WindFrame,
    screen_trajectory,
)
from degali.addons.yawed_crosswind import YawedTrajectory


def test_wind_frame_is_an_exact_horizontal_rotation_and_translation():
    frame = WindFrame(10.0, -2.0, math.pi/2.0)
    assert frame.local(7.0, 2.0) == pytest.approx((4.0, 3.0))
    assert frame.global_(4.0, 3.0) == pytest.approx((7.0, 2.0))


def test_cuboid_screen_reports_one_contiguous_centreline_encounter():
    building = AxisAlignedCuboid(1.0, 2.0, -0.5, 0.5, 0.0, 3.0, "building A")
    result = screen_trajectory(
        [(0.0, 0.0, 1.0), (1.5, 0.0, 1.0), (3.0, 0.0, 1.0)], [building]
    )
    assert not result.free_plume_prediction_applicable
    assert result.requires_obstacle_resolved_model
    assert len(result.encounters) == 1
    hit = result.encounters[0]
    assert hit.obstacle_label == "building A"
    assert hit.obstacle_type == "axis_aligned_cuboid"
    assert hit.first_segment_index == 0
    assert hit.last_segment_index == 1
    assert hit.entry_arc_length_m == pytest.approx(1.0)
    assert hit.exit_arc_length_m == pytest.approx(2.0)
    assert hit.entry_point_m == pytest.approx((1.0, 0.0, 1.0))
    assert hit.exit_point_m == pytest.approx((2.0, 0.0, 1.0))
    assert result.first_contact_arc_length_m == pytest.approx(1.0)
    assert result.free_plume_applicable_through(.99)
    assert not result.free_plume_applicable_through(1.0)
    with pytest.raises(ValueError, match="non-negative"):
        result.free_plume_applicable_through(-.1)
    assert "No obstacle wake" in result.notes[1]


def test_clear_centreline_is_not_mistaken_for_an_obstacle_wake_prediction():
    result = screen_trajectory(
        [(0.0, 0.0, 1.0), (4.0, 0.0, 1.0)],
        [AxisAlignedCuboid(1.0, 2.0, 2.0, 3.0, 0.0, 3.0)],
    )
    assert result.is_clear
    assert result.first_contact_arc_length_m is None
    assert result.free_plume_prediction_applicable
    assert not result.requires_obstacle_resolved_model
    assert "does not prove" in result.notes[1]


def test_oriented_cuboid_uses_its_actual_rotated_footprint_for_contact():
    building = OrientedCuboid(
        center_m=(2.0, 0.0), length_m=4.0, width_m=2.0,
        z_min_m=0.0, z_max_m=3.0, long_axis_bearing_deg=90.0,
        label="north-south-skid",
    )
    result = screen_trajectory(
        [(0.0, 0.0, 1.0), (4.0, 0.0, 1.0)], [building]
    )

    assert len(result.encounters) == 1
    hit = result.encounters[0]
    assert hit.obstacle_type == "oriented_cuboid"
    assert hit.entry_point_m == pytest.approx((1.0, 0.0, 1.0))
    assert hit.exit_point_m == pytest.approx((3.0, 0.0, 1.0))


def test_unbounded_wall_is_screened_in_an_arbitrarily_rotated_wind_frame():
    frame = WindFrame(direction_rad=math.pi/2.0)
    wall = TransverseWall(frame, downwind_m=5.0, base_height_m=0.0, height_m=2.0)
    # For a +y wind, the global y-coordinate is the downwind coordinate.
    result = screen_trajectory([(0.0, 0.0, 1.0), (0.0, 10.0, 1.0)], [wall])
    assert len(result.encounters) == 1
    assert result.encounters[0].entry_point_m == pytest.approx((0.0, 5.0, 1.0))
    assert result.encounters[0].laterally_unbounded_wall
    assert result.encounters[0].minimum_overflight_height_m == pytest.approx(2.0)
    assert result.required_overflight_height_m == pytest.approx(2.0)


def test_finite_wall_can_be_screened_but_geometry_does_not_choose_over_or_around():
    wall = TransverseWall(
        WindFrame(), downwind_m=2.0, base_height_m=0.0, height_m=2.0,
        half_width_m=1.0, thickness_m=0.2,
    )
    result = screen_trajectory([(0.0, 0.0, 1.0), (4.0, 0.0, 1.0)], [wall])
    assert len(result.encounters) == 1
    assert not result.encounters[0].laterally_unbounded_wall
    assert result.required_overflight_height_m is None


def test_zero_thickness_wall_and_grazing_cuboid_are_solid_contacts():
    wall = TransverseWall(WindFrame(), 1.0, 0.0, 2.0)
    cuboid = AxisAlignedCuboid(2.0, 3.0, 0.0, 1.0, 0.0, 2.0)
    result = screen_trajectory([(0.0, 0.0, 1.0), (4.0, 0.0, 1.0)], [wall, cuboid])
    assert [hit.obstacle_label for hit in result.encounters] == ["transverse wall", "cuboid"]


def test_geometry_rejects_degenerate_trajectory_and_invalid_structure():
    with pytest.raises(ValueError, match="positive distance"):
        screen_trajectory([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [])
    with pytest.raises(ValueError, match="strictly positive extent"):
        AxisAlignedCuboid(0.0, 0.0, 0.0, 1.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="half-width"):
        TransverseWall(WindFrame(), 1.0, 0.0, 1.0, half_width_m=0.0)


def test_yawed_trajectory_wrapper_uses_global_path_and_never_mutates_it():
    states = []
    for x in (0.0, 2.0):
        state = [0.0]*10
        state[7], state[8], state[9] = x, 0.0, 1.0
        states.append(state)
    trajectory = YawedTrajectory(SimpleNamespace(), states)
    before = trajectory.states.copy()
    result = trajectory.obstacle_screen([AxisAlignedCuboid(.5, 1.5, -.5, .5, 0., 2.)])
    assert not result.is_clear
    assert trajectory.states == pytest.approx(before)
