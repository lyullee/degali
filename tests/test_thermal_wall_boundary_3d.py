from __future__ import annotations

import numpy as np
import pytest

from degali.addons.thermal_wall_boundary_3d import isothermal_wall_heat_source_3d


def test_isothermal_obstacle_wall_heat_is_added_once_per_exposed_face():
    shape = (3, 3, 3)
    solid = np.zeros(shape, dtype=bool)
    solid[1, 1, 1] = True
    temperature = np.full(shape, 290.0)

    result = isothermal_wall_heat_source_3d(
        temperature, np.ones(shape), np.full(shape, 1000.0), 1.0e-5,
        solid, spacing_m=(0.1, 0.1, 0.1),
        obstacle_wall_temperature_k=300.0,
    )

    assert result.obstacle_wall_face_count == 6
    assert result.ground_wall_cell_count == 0
    assert result.obstacle_wall_heat_rate_w == pytest.approx(0.12)
    assert result.wall_heat_rate_w == pytest.approx(
        np.sum(result.sensible_enthalpy_source_rate_w_m3[~solid]) * 0.001
    )
    assert np.all(result.sensible_enthalpy_source_rate_w_m3[solid] == 0.0)


def test_ground_and_obstacle_wall_heat_use_distinct_boundary_areas():
    shape = (3, 2, 2)
    solid = np.zeros(shape, dtype=bool)
    solid[1, 0, 0] = True
    temperature = np.full(shape, 290.0)

    result = isothermal_wall_heat_source_3d(
        temperature, np.ones(shape), np.full(shape, 1000.0), 1.0e-5,
        solid, spacing_m=(0.1, 0.1, 0.1),
        obstacle_wall_temperature_k=290.0,
        ground_temperature_k=300.0,
    )

    assert result.obstacle_wall_heat_rate_w == pytest.approx(0.0)
    assert result.ground_wall_cell_count == 4
    assert result.ground_wall_heat_rate_w == pytest.approx(0.08)
    assert result.wall_heat_rate_w == pytest.approx(0.08)


def test_wall_conduction_uses_diffusivity_normal_to_each_face():
    shape = (3, 3, 3)
    solid = np.zeros(shape, dtype=bool)
    solid[1, 1, 1] = True
    temperature = np.full(shape, 290.0)

    result = isothermal_wall_heat_source_3d(
        temperature, np.ones(shape), np.ones(shape),
        (0.3 * np.ones(shape), 0.2 * np.ones(shape), 0.1 * np.ones(shape)),
        solid, spacing_m=(2.0, 1.0, 0.5),
        obstacle_wall_temperature_k=300.0,
    )

    # Two x, y, and z faces use their respective normal diffusivity.
    assert result.obstacle_wall_heat_rate_w == pytest.approx(3.0 + 8.0 + 16.0)
    assert result.obstacle_wall_face_count == 6
