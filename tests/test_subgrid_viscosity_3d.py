import numpy as np
import pytest

from degali.addons.subgrid_viscosity_3d import (
    smagorinsky_eddy_viscosity_lattice,
    subgrid_eddy_viscosity_physical,
    subgrid_relaxation_time_lattice,
    velocity_gradient_tensor_lattice,
    wale_eddy_viscosity_lattice,
)


def _coordinates(shape):
    z, y, x = np.meshgrid(
        np.arange(shape[0], dtype=float),
        np.arange(shape[1], dtype=float),
        np.arange(shape[2], dtype=float),
        indexing="ij",
    )
    return x, y, z


def test_velocity_gradient_tensor_recovers_affine_field_exactly():
    shape = (7, 8, 9)
    x, y, z = _coordinates(shape)
    solid = np.zeros(shape, dtype=bool)
    ux = 2.0 * x - 3.0 * y + 0.5 * z
    uy = -x + 4.0 * y + 2.0 * z
    uz = 0.25 * x - 0.75 * y + 1.5 * z

    gradient = velocity_gradient_tensor_lattice(ux, uy, uz, solid)
    expected = np.asarray((
        (2.0, -3.0, 0.5),
        (-1.0, 4.0, 2.0),
        (0.25, -0.75, 1.5),
    ))

    assert np.allclose(gradient, expected)


def test_wale_switches_off_for_uniform_flow_and_pure_shear():
    shape = (8, 9, 10)
    x, y, _z = _coordinates(shape)
    solid = np.zeros(shape, dtype=bool)
    zero = np.zeros(shape)
    uniform = np.full(shape, 0.04)
    uniform_wale = wale_eddy_viscosity_lattice(
        uniform, zero, zero, solid, constant=0.5,
    )
    shear_rate = 0.01
    shear_u = shear_rate * y
    shear_wale = wale_eddy_viscosity_lattice(
        shear_u, zero, zero, solid, constant=0.5,
    )
    shear_smagorinsky = smagorinsky_eddy_viscosity_lattice(
        shear_u, zero, zero, solid, constant=0.1,
    )

    assert np.max(np.abs(uniform_wale)) == 0.0
    assert np.max(np.abs(shear_wale)) < 1.0e-15
    assert np.allclose(shear_smagorinsky, 0.1**2 * shear_rate)


def test_wale_responds_to_rotation_while_smagorinsky_does_not():
    shape = (8, 9, 10)
    x, y, _z = _coordinates(shape)
    solid = np.zeros(shape, dtype=bool)
    rotation_rate = 0.02
    ux = -rotation_rate * y
    uy = rotation_rate * x
    uz = np.zeros(shape)
    constant = 0.5

    wale = wale_eddy_viscosity_lattice(
        ux, uy, uz, solid, constant=constant,
    )
    smagorinsky = smagorinsky_eddy_viscosity_lattice(
        ux, uy, uz, solid, constant=0.1,
    )
    expected = constant**2 * (2.0 / 3.0)**0.25 * rotation_rate

    assert np.allclose(wale, expected)
    assert np.max(np.abs(smagorinsky)) < 1.0e-15


def test_subgrid_relaxation_caps_tau_and_zeroes_solid_cells():
    shape = (8, 9, 10)
    x, y, _z = _coordinates(shape)
    solid = np.zeros(shape, dtype=bool)
    solid[2:5, 3:6, 4:7] = True
    ux = -0.3 * y
    uy = 0.3 * x
    uz = np.zeros(shape)

    tau, viscosity = subgrid_relaxation_time_lattice(
        ux, uy, uz, solid,
        base_lattice_viscosity=0.01,
        model="wale",
        smagorinsky_constant=0.1,
        wale_constant=0.8,
        maximum_relaxation_time=0.56,
    )

    assert np.all(viscosity[solid] == 0.0)
    assert np.max(tau) == pytest.approx(0.56)
    assert np.all(tau >= 0.53)


def test_subgrid_model_rejects_unknown_name():
    shape = (3, 3, 3)
    zero = np.zeros(shape)
    with pytest.raises(ValueError, match="model"):
        subgrid_relaxation_time_lattice(
            zero, zero, zero, np.zeros(shape, dtype=bool),
            base_lattice_viscosity=0.01,
            model="invented",
            smagorinsky_constant=0.1,
            wale_constant=0.5,
            maximum_relaxation_time=1.9,
        )


def test_physical_subgrid_viscosity_applies_filter_width_scaling():
    shape = (5, 6, 7)
    z, y, x = np.meshgrid(
        np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]),
        indexing="ij",
    )
    spacing = 0.04
    u = 0.8 * x - 0.3 * y
    v = 0.2 * x + 0.5 * z
    w = -0.4 * y + 0.1 * z
    solid = np.zeros(shape, dtype=bool)
    lattice = wale_eddy_viscosity_lattice(
        u, v, w, solid, constant=0.5,
    )
    physical = subgrid_eddy_viscosity_physical(
        u, v, w, solid,
        grid_spacing_m=spacing,
        model="wale",
        smagorinsky_constant=0.1,
        wale_constant=0.5,
    )

    assert physical == pytest.approx(lattice * spacing)
