import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.site_geometry import (
    AxisAlignedCuboid,
    CircularArcWall,
    CrosswindCylinder,
    OrientedCuboid,
    TransverseWall,
    WindFrame,
    circular_arc_surface_face_masks,
    circular_arc_surface_cell_metrics,
    circular_arc_surface_link_masks,
    circular_arc_surface_link_fractions,
    segment_circular_arc_surface_intersections,
    screen_trajectory,
)


def test_circular_arc_wall_has_exact_sector_semantics_and_voxel_mask():
    wall = CircularArcWall(
        center_m=(0.0, 0.0), radius_m=5.0, thickness_m=1.0,
        z_min_m=0.0, z_max_m=2.0, arc_start_deg=0.0, arc_end_deg=180.0,
        label="semicircular fence",
    )
    assert wall.contains_point((5.0, 0.0, 1.0))
    assert wall.contains_point((0.0, 5.0, 1.0))
    assert wall.contains_point((-5.0, 0.0, 1.0))
    assert not wall.contains_point((0.0, -5.0, 1.0))
    assert not wall.contains_point((0.0, 0.0, 1.0))
    assert wall.bounding_box_m == pytest.approx((-5.5, 5.5, 0.0, 5.5, 0.0, 2.0))

    coordinates = np.arange(-5.0, 5.1, 1.0)
    z = np.asarray((0.5, 1.5, 2.5))
    mask = wall.voxel_mask(coordinates, coordinates, z)
    assert mask.shape == (3, 11, 11)
    assert mask[0, 10, 5]  # (x=0, y=5)
    assert not mask[0, 0, 5]  # (x=0, y=-5)
    assert not np.any(mask[2])


def test_zero_thickness_circular_arc_is_an_ideal_surface_without_voxel_volume():
    wall = CircularArcWall(
        center_m=(0.0, 0.0), radius_m=2.0, thickness_m=0.0,
        z_min_m=0.0, z_max_m=1.0, arc_start_deg=0.0, arc_end_deg=180.0,
    )
    coordinates = np.arange(-2.5, 3.0, 1.0)
    z = np.asarray((0.5, 1.5))

    assert wall.contains_point((2.0, 0.0, 0.5))
    assert not np.any(wall.voxel_mask(coordinates, coordinates, z))


def test_crosswind_cylinder_has_exact_round_section_and_finite_span():
    rod = CrosswindCylinder(
        center_x_m=2.0, center_z_m=0.5, radius_m=0.1,
        y_min_m=0.2, y_max_m=1.7, label="ExJet horizontal rod",
    )
    x = np.asarray((1.85, 1.95, 2.05, 2.15))
    y = np.asarray((0.15, 0.25, 1.65, 1.75))
    z = np.asarray((0.35, 0.45, 0.55, 0.65))
    mask = rod.voxel_mask(x, y, z)

    assert rod.bounding_box_m == pytest.approx((1.9, 2.1, 0.2, 1.7, 0.4, 0.6))
    assert rod.contains_point((2.0, 0.3, 0.5))
    assert not rod.contains_point((2.11, 0.3, 0.5))
    assert not rod.contains_point((2.0, 1.71, 0.5))
    assert mask.shape == (4, 4, 4)
    assert np.any(mask[:, 1:3, :])
    assert not np.any(mask[:, 0, :]) and not np.any(mask[:, 3, :])


def test_zero_thickness_arc_intersection_returns_chord_contacts_and_tangent():
    semicircle = dict(
        center_xy_m=(0.0, 0.0), radius_m=5.0, z_min_m=0.0, z_max_m=2.0,
        arc_start_deg=0.0, arc_end_deg=180.0,
    )
    chord = segment_circular_arc_surface_intersections(
        (-6.0, 2.0, 1.0), (6.0, 2.0, 1.0), **semicircle,
    )
    assert len(chord) == 2
    assert chord[0].entry_point_m == pytest.approx((-math.sqrt(21.0), 2.0, 1.0))
    assert chord[1].entry_point_m == pytest.approx((math.sqrt(21.0), 2.0, 1.0))
    assert all(hit.entry_fraction == pytest.approx(hit.exit_fraction) for hit in chord)

    tangent = segment_circular_arc_surface_intersections(
        (-6.0, 5.0, 1.0), (6.0, 5.0, 1.0), **semicircle,
    )
    assert len(tangent) == 1
    assert tangent[0].entry_point_m == pytest.approx((0.0, 5.0, 1.0))


def test_zero_thickness_arc_intersection_respects_arc_and_height_bounds():
    arc = dict(
        center_xy_m=(0.0, 0.0), radius_m=5.0, z_min_m=0.0, z_max_m=2.0,
        arc_start_deg=0.0, arc_end_deg=180.0,
    )
    assert segment_circular_arc_surface_intersections(
        (-6.0, -2.0, 1.0), (6.0, -2.0, 1.0), **arc,
    ) == ()
    assert segment_circular_arc_surface_intersections(
        (-6.0, 2.0, 3.0), (6.0, 2.0, 3.0), **arc,
    ) == ()
    lower_half = segment_circular_arc_surface_intersections(
        (-6.0, -2.0, 1.0), (6.0, -2.0, 1.0),
        **{**arc, "arc_start_deg": 180.0, "arc_end_deg": 0.0},
    )
    assert len(lower_half) == 2


def test_zero_thickness_arc_face_masks_match_finite_volume_shapes_and_height():
    coordinates = np.arange(-4.5, 5.0, 1.0)
    z = np.asarray((0.5, 1.5, 2.5))
    x_faces, y_faces, z_faces = circular_arc_surface_face_masks(
        coordinates, coordinates, z,
        center_xy_m=(0.0, 0.0), radius_m=5.0,
        z_min_m=0.0, z_max_m=2.0,
        arc_start_deg=0.0, arc_end_deg=180.0,
    )
    assert x_faces.shape == (3, 10, 11)
    assert y_faces.shape == (3, 11, 10)
    assert z_faces.shape == (4, 10, 10)
    assert np.any(x_faces[0]) and np.any(y_faces[0])
    assert not np.any(x_faces[2]) and not np.any(y_faces[2])
    assert not np.any(z_faces)
    # An upper semicircle cannot block faces at the corresponding lower arc.
    assert not x_faces[0, 0, 5]
    assert not y_faces[0, 5, 5]


def test_zero_thickness_arc_cell_surface_area_is_grid_alignment_invariant():
    z = np.asarray((0.25, 0.75, 1.25, 1.75))
    totals = []
    for spacing, first_center in ((1.0, -9.5), (2.0, -9.0)):
        coordinates = np.arange(first_center, 10.0, spacing)
        area, normal_x, normal_y, wall_distance = circular_arc_surface_cell_metrics(
            coordinates, coordinates, z,
            center_xy_m=(0.0, 0.0), radius_m=5.0,
            z_min_m=0.0, z_max_m=2.0,
            arc_start_deg=0.0, arc_end_deg=180.0,
        )
        assert area.shape == (z.size, coordinates.size, coordinates.size)
        assert np.all(np.isfinite(wall_distance[area > 0.0]))
        assert np.all(np.isinf(wall_distance[area == 0.0]))
        totals.append((area.sum(), normal_x.sum(), normal_y.sum()))

    expected_area = math.pi * 5.0 * 2.0
    assert totals[0][0] == pytest.approx(expected_area, rel=1.0e-12)
    assert totals[1][0] == pytest.approx(expected_area, rel=1.0e-12)
    assert totals[0][1] == pytest.approx(0.0, abs=1.0e-12)
    assert totals[1][1] == pytest.approx(0.0, abs=1.0e-12)
    assert totals[0][2] == pytest.approx(20.0, rel=1.0e-12)
    assert totals[1][2] == pytest.approx(20.0, rel=1.0e-12)
def test_zero_thickness_arc_link_masks_find_d3q19_chord_contacts():
    coordinates = np.asarray((-1.5, -0.5, 0.5, 1.5))
    z = np.asarray((0.5, 1.5))
    directions = np.asarray(((1, 0, 0), (-1, 0, 0), (1, -1, 0), (0, 0, 1)))
    masks = circular_arc_surface_link_masks(
        coordinates, coordinates, z, directions,
        center_xy_m=(0.0, 0.0), radius_m=1.0,
        z_min_m=0.0, z_max_m=1.0,
        arc_start_deg=0.0, arc_end_deg=180.0,
    )
    assert masks.shape == (4, 2, 4, 4)
    # From (-1.5, 0.5) to (-0.5, 0.5), the link crosses the upper arc.
    assert masks[0, 0, 2, 0]
    # The reverse link sees the same contact from its opposite endpoint.
    assert masks[1, 0, 2, 1]
    # A diagonal chord crosses the arc at this height as well.
    assert masks[2, 0, 2, 0]
    # Vertical links are tangent to the vertical wall and do not cross it.
    assert not np.any(masks[3])


def test_zero_thickness_arc_link_fractions_return_first_hit_along_link():
    coordinates = np.asarray((-1.5, -0.5, 0.5, 1.5))
    z = np.asarray((0.5, 1.5))
    fractions = circular_arc_surface_link_fractions(
        coordinates, coordinates, z, ((1, 0, 0),),
        center_xy_m=(0.0, 0.0), radius_m=1.0,
        z_min_m=0.0, z_max_m=1.0,
        arc_start_deg=0.0, arc_end_deg=180.0,
    )
    assert fractions.shape == (1, 2, 4, 4)
    assert fractions[0, 0, 2, 0] == pytest.approx(1.5 - math.sqrt(0.75))
    assert np.isnan(fractions[0, 1, 2, 0])


def test_zero_thickness_arc_link_masks_respect_height_and_circumferential_sector():
    coordinates = np.asarray((-1.5, -0.5, 0.5, 1.5))
    z = np.asarray((0.5, 1.5))
    mask = circular_arc_surface_link_masks(
        coordinates, coordinates, z, ((1, 0, 0),),
        center_xy_m=(0.0, 0.0), radius_m=1.0,
        z_min_m=1.0, z_max_m=2.0,
        arc_start_deg=0.0, arc_end_deg=90.0,
    )
    assert not np.any(mask[0, 0])  # below the wall's lower edge
    assert not mask[0, 1, 2, 0]  # intersection angle is outside this quadrant


def test_zero_thickness_arc_link_masks_require_cubic_grid_and_valid_links():
    coordinates = np.asarray((-1.5, -0.5, 0.5, 1.5))
    kwargs = dict(
        center_xy_m=(0.0, 0.0), radius_m=1.0,
        z_min_m=0.0, z_max_m=2.0,
        arc_start_deg=0.0, arc_end_deg=180.0,
    )
    with pytest.raises(ValueError, match="cubic"):
        circular_arc_surface_link_masks(
            coordinates, coordinates, np.asarray((0.25, 0.75, 1.25, 1.75)),
            ((1, 0, 0),), **kwargs,
        )
    with pytest.raises(ValueError, match="D3Q19 integer links"):
        circular_arc_surface_link_masks(
            coordinates, coordinates, coordinates,
            ((2, 0, 0),), **kwargs,
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
