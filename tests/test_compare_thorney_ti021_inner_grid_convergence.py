import numpy as np
import pytest

from tools.compare_thorney_ti021_inner_grid_convergence import (
    coarsen_fluid_mask_two_to_one,
    coarsen_two_to_one,
    obstacle_region_masks,
    physics_configuration_differences,
    trace_metrics,
)


def test_two_to_one_coarsening_is_a_volume_average():
    field = np.arange(64.0).reshape(4, 4, 4)
    result = coarsen_two_to_one(field)
    assert result.shape == (2, 2, 2)
    assert result[0, 0, 0] == pytest.approx(np.mean(field[:2, :2, :2]))
    assert result[1, 1, 1] == pytest.approx(np.mean(field[2:, 2:, 2:]))


def test_common_fluid_coarsening_rejects_any_cut_block():
    solid = np.zeros((4, 4, 4), dtype=bool)
    solid[1, 1, 1] = True
    fluid = coarsen_fluid_mask_two_to_one(solid)
    assert fluid.shape == (2, 2, 2)
    assert not fluid[0, 0, 0]
    assert np.count_nonzero(fluid) == 7


def test_trace_metrics_preserve_peak_dose_and_threshold_arrival():
    time_s = np.asarray((0.0, 1.0, 2.0, 3.0))
    traces = np.asarray(((0.0, 0.05, 0.2, 0.1), (0.0, 0.0, 0.0, 0.0)))
    metrics = trace_metrics(time_s, traces, threshold=0.1)
    assert metrics["peak_vol_percent"] == pytest.approx((0.2, 0.0))
    assert metrics["dose_vol_percent_s"] == pytest.approx((0.3, 0.0))
    assert metrics["arrival_s"][0] == pytest.approx(2.0)
    assert np.isnan(metrics["arrival_s"][1])


def test_physics_configuration_gate_ignores_only_grid_and_time_discretization():
    coarse = {
        "duration_s": 600.0, "time_step_s": 0.125, "cell_size_m": 5.0,
        "cells": [40, 40, 8], "active_steps": 4800,
        "circular_arc_representation": "finite_thickness_voxel",
        "corrector_tolerance": 1.0e-6, "initial_binary_cloud_momentum": "quiescent",
    }
    fine = {
        "duration_s": 60.0, "time_step_s": 0.0625, "cell_size_m": 2.5,
        "cells": [80, 80, 16], "active_steps": 960,
        "circular_arc_representation": "zero_thickness_surface",
        "corrector_tolerance": 1.0e-4, "initial_binary_cloud_momentum": "quiescent",
    }
    assert physics_configuration_differences(coarse, fine) == [
        "circular_arc_representation", "corrector_tolerance",
    ]


def test_obstacle_regions_sample_front_roof_side_and_wake():
    coordinates = np.arange(2.5, 200.0, 5.0)
    z = np.arange(2.5, 40.0, 5.0)
    solid = np.zeros((z.size, coordinates.size, coordinates.size), dtype=bool)
    masks = obstacle_region_masks(coordinates, coordinates, z, solid)
    assert set(masks) == {
        "windward_inner_face_shell", "roof_shear_band",
        "side_tip_bypass", "exterior_lee_wake_shell",
    }
    assert all(np.any(mask) for mask in masks.values())
