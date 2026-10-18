import numpy as np
import pytest

from degali.addons.hybrid_amr import (
    conservative_patch_remap, make_obstacle_patch_bounds, overlap_mass_remap_to_patch,
    extract_patch_boundary_flux, PatchFluxAccumulator, cartesian_face_coordinates,
)


def test_patch_bounds_are_buffered_and_clipped():
    patch = make_obstacle_patch_bounds(
        (10, 20, 5, 15, 0, 8), upstream_m=4, downstream_m=12,
        lateral_m=3, vertical_m=2, domain_m=(30, 25, 20),
    )
    assert (patch.x_min, patch.x_max) == (6, 30)
    assert (patch.y_min, patch.y_max) == (2, 18)
    assert (patch.z_min, patch.z_max) == (0, 10)


def test_patch_remap_preserves_mass():
    coarse = np.arange(8, dtype=float).reshape(2, 2, 2)
    patch = conservative_patch_remap(coarse, coarse_cell_volume_m3=2.0, patch_shape=(4, 4, 4))
    assert np.sum(patch) == pytest.approx(np.sum(coarse))

def test_overlap_remap_selects_only_patch_mass():
    bounds = make_obstacle_patch_bounds(
        (1, 2, 1, 2, 0, 1), upstream_m=0, downstream_m=0,
        lateral_m=0, vertical_m=0, domain_m=(3, 3, 2))
    field, selected = overlap_mass_remap_to_patch(
        np.ones((2, 3, 3)), x_m=np.array([.5, 1.5, 2.5]),
        y_m=np.array([.5, 1.5, 2.5]), z_m=np.array([.5, 1.5]),
        patch=bounds, patch_shape=(4, 4, 4))
    assert selected == pytest.approx(1.0)
    assert np.sum(field) == pytest.approx(selected)

def test_patch_boundary_flux_integrates_faces_and_net_outward_rate():
    patch = make_obstacle_patch_bounds((1, 2, 1, 2, 0, 1), upstream_m=0,
                                       downstream_m=0, lateral_m=0, vertical_m=0,
                                       domain_m=(3, 3, 2))
    flux = extract_patch_boundary_flux(
        flux_x=np.ones((2, 3, 4)), flux_y=np.ones((2, 4, 3)),
        flux_z=np.ones((3, 3, 3)), x_face_m=np.array([0, 1, 2, 3.]),
        y_face_m=np.array([0, 1, 2, 3.]), z_face_m=np.array([0, 1, 2.]), patch=patch)
    assert flux["x_min_kg_s"] == pytest.approx(1.0)
    assert flux["net_outward_kg_s"] == pytest.approx(0.0)

def test_patch_flux_accumulator_is_transport_callback_compatible():
    patch = make_obstacle_patch_bounds((1, 2, 1, 2, 0, 1), upstream_m=0,
                                       downstream_m=0, lateral_m=0, vertical_m=0,
                                       domain_m=(3, 3, 2))
    collector = PatchFluxAccumulator(np.array([0, 1, 2, 3.]),
                                     np.array([0, 1, 2, 3.]),
                                     np.array([0, 1, 2.]), patch)
    collector(np.ones((2, 3, 4)), np.ones((2, 4, 3)), np.ones((3, 3, 3)))
    assert len(collector.net_outward_kg_s) == 1
    assert collector.rates_kg_s.shape == (1,)

def test_face_coordinates_match_solver_grid_convention():
    x, y, z = cartesian_face_coordinates((4., 2.4, 2.), (16, 10, 8))
    assert len(x) == 17 and x[0] == 0 and x[-1] == 4
    assert len(y) == 11 and len(z) == 9
