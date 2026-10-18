import numpy as np
import pytest

from degali.addons.variable_density_momentum_3d import (
    VARIABLE_DENSITY_MOMENTUM_3D_SCHEMA,
    _upwind_momentum_fluxes,
    advance_variable_density_momentum_3d,
    smooth_wall_log_law_force_density_3d,
)


def _faces(shape, density=1.2, u=0.0):
    nz, ny, nx = shape
    return (
        np.full((nz, ny, nx + 1), density * u),
        np.zeros((nz, ny + 1, nx)),
        np.zeros((nz + 1, ny, nx)),
    )


def test_pressure_outlet_backflow_uses_ambient_rest_momentum_state():
    scalar = np.array([[[0.7]]])
    mass_flux_x = np.array([[[0.0, -1.0]]])
    mass_flux_y = np.zeros((1, 2, 1))
    mass_flux_z = np.zeros((2, 1, 1))
    solid = np.zeros((1, 1, 1), dtype=bool)

    flux_x, _, _ = _upwind_momentum_fluxes(
        scalar, mass_flux_x, mass_flux_y, mass_flux_z, solid,
        component_axis=0, inlet_value=np.zeros((1, 1)),
        nonstreamwise_wall_boundary="all_no_slip_tunnel",
    )

    assert flux_x[0, 0, -1] == 0.0


def test_pressure_outlet_outflow_uses_interior_momentum_state():
    scalar = np.array([[[0.7]]])
    mass_flux_x = np.array([[[0.0, 1.0]]])
    mass_flux_y = np.zeros((1, 2, 1))
    mass_flux_z = np.zeros((2, 1, 1))
    solid = np.zeros((1, 1, 1), dtype=bool)

    flux_x, _, _ = _upwind_momentum_fluxes(
        scalar, mass_flux_x, mass_flux_y, mass_flux_z, solid,
        component_axis=0, inlet_value=np.zeros((1, 1)),
        nonstreamwise_wall_boundary="all_no_slip_tunnel",
    )

    assert flux_x[0, 0, -1] == pytest.approx(0.7)


def test_uniform_free_stream_and_global_momentum_are_preserved():
    shape = (5, 6, 8)
    density = np.full(shape, 1.2)
    u = np.full(shape, 0.7)
    zero = np.zeros(shape)
    fluxes = _faces(shape, density=1.2, u=0.7)

    result = advance_variable_density_momentum_3d(
        density, u, zero, zero,
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.01,
        kinematic_viscosity_m2_s=0.0,
        inlet_u_m_s=0.7,
    )

    assert result.schema == VARIABLE_DENSITY_MOMENTUM_3D_SCHEMA
    assert result.provisional_u_m_s == pytest.approx(u)
    assert result.provisional_v_m_s == pytest.approx(zero)
    assert result.provisional_w_m_s == pytest.approx(zero)
    assert result.momentum_provisional_kg_m_s == pytest.approx(
        result.momentum_initial_kg_m_s
    )
    assert np.max(np.abs(result.momentum_balance_residual_kg_m_s)) < 1.0e-14


def test_uniform_body_acceleration_has_exact_discrete_impulse():
    shape = (4, 5, 6)
    density = np.full(shape, 1.1)
    zero = np.zeros(shape)
    fluxes = _faces(shape)
    dt = 0.02

    result = advance_variable_density_momentum_3d(
        density, zero, zero, zero,
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.25, 0.25, 0.25),
        time_step_s=dt,
        kinematic_viscosity_m2_s=0.0,
        acceleration_z_m_s2=0.3,
    )

    assert result.provisional_w_m_s == pytest.approx(0.3 * dt)
    assert result.body_force_n[2] > 0.0
    assert np.max(np.abs(result.momentum_balance_residual_kg_m_s)) < 1.0e-14


def test_no_slip_ground_and_obstacle_viscosity_dissipate_momentum():
    shape = (6, 6, 8)
    density = np.ones(shape)
    u = np.full(shape, 0.2)
    zero = np.zeros(shape)
    solid = np.zeros(shape, dtype=bool)
    solid[1:4, 2:4, 3:5] = True
    u[solid] = 0.0
    fluxes = _faces(shape, density=1.0, u=0.0)

    result = advance_variable_density_momentum_3d(
        density, u, zero, zero,
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.002,
        kinematic_viscosity_m2_s=0.05,
        inlet_u_m_s=0.2,
    )

    fluid = ~solid
    assert np.all(result.provisional_u_m_s[solid] == 0.0)
    assert np.sum(result.provisional_u_m_s[fluid]) < np.sum(u[fluid])
    assert result.viscous_boundary_force_n[0] < 0.0
    assert np.max(np.abs(result.momentum_balance_residual_kg_m_s)) < 1.0e-13


def test_tunnel_sidewalls_and_roof_add_tangential_viscous_drag():
    shape = (6, 7, 8)
    density = np.ones(shape)
    u = np.full(shape, 0.2)
    zero = np.zeros(shape)
    fluxes = _faces(shape, density=1.0, u=0.0)
    common = dict(
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.002,
        kinematic_viscosity_m2_s=0.05,
        inlet_u_m_s=0.2,
    )
    free = advance_variable_density_momentum_3d(
        density, u, zero, zero, **common,
    )
    tunnel = advance_variable_density_momentum_3d(
        density, u, zero, zero,
        nonstreamwise_wall_boundary="all_no_slip_tunnel",
        **common,
    )

    assert tunnel.viscous_boundary_force_n[0] < free.viscous_boundary_force_n[0]
    assert np.sum(tunnel.provisional_u_m_s) < np.sum(free.provisional_u_m_s)
    assert np.max(np.abs(tunnel.momentum_balance_residual_kg_m_s)) < 1.0e-13


def test_all_symmetry_mode_removes_lower_tangential_wall_drag():
    shape = (6, 7, 8)
    density = np.ones(shape)
    u = np.full(shape, 0.2)
    zero = np.zeros(shape)
    fluxes = _faces(shape, density=1.0, u=0.0)
    common = dict(
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.002,
        kinematic_viscosity_m2_s=0.05,
        inlet_u_m_s=0.2,
    )
    ground = advance_variable_density_momentum_3d(
        density, u, zero, zero, **common,
    )
    symmetry = advance_variable_density_momentum_3d(
        density, u, zero, zero,
        nonstreamwise_wall_boundary="all_symmetry_free_slip",
        **common,
    )

    assert symmetry.viscous_boundary_force_n[0] > ground.viscous_boundary_force_n[0]
    assert np.sum(symmetry.provisional_u_m_s) > np.sum(ground.provisional_u_m_s)
    assert np.max(np.abs(symmetry.momentum_balance_residual_kg_m_s)) < 1.0e-13


def test_pressure_open_faces_remove_normal_and_tangential_wall_drag():
    shape = (6, 7, 8)
    density = np.ones(shape)
    zero = np.zeros(shape)
    w = np.full(shape, 0.2)
    fluxes = _faces(shape, density=1.0, u=0.0)
    common = dict(
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.002,
        kinematic_viscosity_m2_s=0.05,
        inlet_w_m_s=0.2,
    )
    symmetry = advance_variable_density_momentum_3d(
        density, zero, zero, w,
        nonstreamwise_wall_boundary="all_symmetry_free_slip",
        **common,
    )
    pressure_open = advance_variable_density_momentum_3d(
        density, zero, zero, w,
        nonstreamwise_wall_boundary="all_pressure_open",
        **common,
    )

    assert symmetry.viscous_boundary_force_n[2] < 0.0
    assert pressure_open.viscous_boundary_force_n[2] == pytest.approx(
        0.0, abs=1.0e-15,
    )
    assert pressure_open.provisional_w_m_s == pytest.approx(w)
    assert np.max(
        np.abs(pressure_open.momentum_balance_residual_kg_m_s)
    ) < 1.0e-13


def test_declared_source_velocity_supplies_matching_momentum():
    shape = (4, 4, 5)
    density = np.ones(shape)
    zero = np.zeros(shape)
    source = np.zeros(shape)
    source[1, 2, 2] = 0.5
    fluxes = _faces(shape)

    result = advance_variable_density_momentum_3d(
        density, zero, zero, zero,
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.1, 0.1, 0.1),
        time_step_s=0.01,
        kinematic_viscosity_m2_s=0.0,
        mass_source_rate_kg_m3_s=source,
        source_velocity_m_s=(2.0, 0.0, 0.0),
    )

    assert result.source_momentum_force_n[0] == pytest.approx(0.001)
    assert result.provisional_u_m_s[1, 2, 2] == pytest.approx(0.01)
    assert np.max(np.abs(result.momentum_balance_residual_kg_m_s)) < 1.0e-14


def test_complete_zy_inlet_plane_drives_local_boundary_momentum_flux():
    shape = (3, 4, 5)
    density = np.ones(shape)
    zero = np.zeros(shape)
    fluxes = _faces(shape, density=1.0, u=0.5)
    inlet = np.asarray((
        (0.1, 0.2, 0.3, 0.4),
        (0.2, 0.3, 0.4, 0.5),
        (0.3, 0.4, 0.5, 0.6),
    ))

    result = advance_variable_density_momentum_3d(
        density, zero, zero, zero,
        mass_flux_x_kg_m2_s=fluxes[0],
        mass_flux_y_kg_m2_s=fluxes[1],
        mass_flux_z_kg_m2_s=fluxes[2],
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.1, 0.1, 0.1),
        time_step_s=0.01,
        kinematic_viscosity_m2_s=0.0,
        inlet_u_m_s=inlet,
    )

    expected_first_cell = 0.01 * 0.5 * inlet / 0.1
    assert result.provisional_u_m_s[:, :, 0] == pytest.approx(
        expected_first_cell
    )
    assert np.max(np.abs(result.momentum_balance_residual_kg_m_s)) < 1.0e-14


def test_smooth_wall_log_law_damps_only_tangential_velocity_and_is_conservative():
    shape = (1, 1, 1)
    density = np.full(shape, 1.2)
    zeros = np.zeros(shape)
    area = np.full(shape, 2.0)
    normal_x_area = area.copy()
    normal_y_area = zeros.copy()
    distance = np.full(shape, 0.5)
    force, diagnostics = smooth_wall_log_law_force_density_3d(
        zeros, np.full(shape, 10.0), zeros, density,
        surface_area_m2=area,
        area_normal_x_m2=normal_x_area,
        area_normal_y_m2=normal_y_area,
        wall_distance_m=distance,
        spacing_m=(1.0, 1.0, 1.0),
        molecular_kinematic_viscosity_m2_s=1.5e-5,
    )

    assert force[0] == pytest.approx(zeros)
    assert force[1][0, 0, 0] < 0.0
    assert force[2] == pytest.approx(zeros)
    assert diagnostics["surface_cell_count"] == 1
    assert diagnostics["maximum_friction_velocity_m_s"] > 0.0
    assert diagnostics["maximum_wall_shear_stress_pa"] > 0.0
    assert diagnostics["net_wall_force_y_n"] == pytest.approx(
        force[1].sum(), rel=1e-12,
    )

    normal_force, normal_diagnostics = smooth_wall_log_law_force_density_3d(
        np.full(shape, 10.0), zeros, zeros, density,
        surface_area_m2=area,
        area_normal_x_m2=normal_x_area,
        area_normal_y_m2=normal_y_area,
        wall_distance_m=distance,
        spacing_m=(1.0, 1.0, 1.0),
        molecular_kinematic_viscosity_m2_s=1.5e-5,
    )
    assert all(np.allclose(item, 0.0) for item in normal_force)
    assert normal_diagnostics["maximum_wall_shear_stress_pa"] == 0.0


def test_reichardt_wall_law_recovers_prescribed_friction_velocity():
    shape = (1, 1, 1)
    density = np.full(shape, 1.2)
    zeros = np.zeros(shape)
    friction_velocity = 0.32
    distance = 0.4
    viscosity = 1.5e-5
    y_plus = friction_velocity * distance / viscosity
    kappa = 0.41
    scaled_y = y_plus / 11.0
    u_plus = (
        np.log1p(kappa * y_plus) / kappa
        + 7.8 * (
            1.0 - np.exp(-scaled_y)
            - scaled_y * np.exp(-y_plus / 3.0)
        )
    )
    tangential_speed = friction_velocity * u_plus
    _force, diagnostics = smooth_wall_log_law_force_density_3d(
        zeros, np.full(shape, tangential_speed), zeros, density,
        surface_area_m2=np.ones(shape),
        area_normal_x_m2=np.ones(shape),
        area_normal_y_m2=zeros,
        wall_distance_m=np.full(shape, distance),
        spacing_m=(1.0, 1.0, 1.0),
        molecular_kinematic_viscosity_m2_s=viscosity,
    )

    assert diagnostics["maximum_friction_velocity_m_s"] == pytest.approx(
        friction_velocity, rel=1.0e-10,
    )
    assert diagnostics["maximum_wall_shear_stress_pa"] == pytest.approx(
        1.2 * friction_velocity**2, rel=1.0e-10,
    )
