import numpy as np
import pytest

from degali.addons.variable_density_projection_3d import (
    VARIABLE_DENSITY_MASS_FLUX_PROJECTION_SCHEMA,
    _projection_component_labels,
    project_variable_density_mass_flux_3d,
)


def _zero_flux(shape):
    nz, ny, nx = shape
    return (
        np.zeros((nz, ny, nx + 1)),
        np.zeros((nz, ny + 1, nx)),
        np.zeros((nz + 1, ny, nx)),
    )


def test_projection_component_labels_are_cached_for_fixed_geometry():
    shape = (3, 4, 5)
    solid = np.zeros(shape, dtype=bool)
    blocked = tuple(
        np.zeros(face_shape, dtype=bool)
        for face_shape in ((3, 4, 6), (3, 5, 5), (4, 4, 5))
    )
    blocked_bytes = tuple(mask.tobytes() for mask in blocked)
    args = (
        shape,
        solid.tobytes(),
        (1.0, 1.0, 1.0),
        ("x_max",),
        blocked_bytes,
    )
    _projection_component_labels.cache_clear()

    first = _projection_component_labels(*args)
    second = _projection_component_labels(*args)

    assert first[0] == 1
    assert np.all(first[1] == 0)
    assert not first[1].flags.writeable
    assert second[1] is first[1]
    assert _projection_component_labels.cache_info().hits == 1


def test_projection_component_cache_separates_distinct_obstacle_topologies():
    shape = (5, 5, 5)
    open_solid = np.zeros(shape, dtype=bool)
    pocket_solid = open_solid.copy()
    pocket_solid[:, :, 2] = True
    blocked_bytes = tuple(
        mask.tobytes()
        for mask in (
            np.zeros((5, 5, 6), dtype=bool),
            np.zeros((5, 6, 5), dtype=bool),
            np.zeros((6, 5, 5), dtype=bool),
        )
    )
    common = (
        shape, (1.0, 1.0, 1.0), ("x_min", "x_max"), blocked_bytes,
    )
    _projection_component_labels.cache_clear()

    connected = _projection_component_labels(
        common[0], open_solid.tobytes(), *common[1:],
    )
    pocket = _projection_component_labels(
        common[0], pocket_solid.tobytes(), *common[1:],
    )

    assert connected[0] == 1
    assert pocket[0] == 2
    assert connected[1] is not pocket[1]
    assert _projection_component_labels.cache_info().misses == 2


def test_projection_component_cache_separates_blocked_surface_topologies():
    shape = (4, 5, 6)
    solid = np.zeros(shape, dtype=bool)
    face_shapes = ((4, 5, 7), (4, 6, 6), (5, 5, 6))
    open_faces = ("x_min", "x_max")
    no_surface = tuple(np.zeros(face_shape, dtype=bool) for face_shape in face_shapes)
    surface = tuple(mask.copy() for mask in no_surface)
    surface[0][:, :, 3] = True
    _projection_component_labels.cache_clear()

    connected = _projection_component_labels(
        shape, solid.tobytes(), (1.0, 1.0, 1.0), open_faces,
        tuple(mask.tobytes() for mask in no_surface),
    )
    split_by_surface = _projection_component_labels(
        shape, solid.tobytes(), (1.0, 1.0, 1.0), open_faces,
        tuple(mask.tobytes() for mask in surface),
    )

    assert connected[0] == 1
    assert split_by_surface[0] == 2
    assert connected[1] is not split_by_surface[1]
    assert _projection_component_labels.cache_info().misses == 2


def test_projection_reduces_random_uniform_density_continuity_residual():
    shape = (7, 8, 10)
    rng = np.random.default_rng(42)
    solid = np.zeros(shape, dtype=bool)
    fx, fy, fz = _zero_flux(shape)
    fx[:, :, 1:-1] = rng.normal(0.0, 0.1, fx[:, :, 1:-1].shape)
    fy[:, 1:-1, :] = rng.normal(0.0, 0.1, fy[:, 1:-1, :].shape)
    fz[1:-1, :, :] = rng.normal(0.0, 0.1, fz[1:-1, :, :].shape)
    density = np.full(shape, 1.2)

    result = project_variable_density_mass_flux_3d(
        fx, fy, fz,
        density_old_kg_m3=density,
        density_new_kg_m3=density,
        mass_source_rate_kg_m3_s=0.0,
        solid_mask=solid,
        spacing_m=(0.2, 0.25, 0.3),
        time_step_s=0.05,
    )

    assert result.schema == VARIABLE_DENSITY_MASS_FLUX_PROJECTION_SCHEMA
    assert result.maximum_absolute_residual_before_kg_m3_s > 0.1
    assert result.maximum_absolute_residual_after_kg_m3_s < 1.0e-9
    assert abs(result.global_mass_balance_residual_kg_s) < 1.0e-10
    assert np.all(np.isfinite(result.pressure_correction_pa))


def test_uniform_density_decrease_generates_matching_outlet_mass_flux():
    shape = (6, 7, 9)
    solid = np.zeros(shape, dtype=bool)
    fluxes = _zero_flux(shape)
    old = np.full(shape, 1.20)
    new = np.full(shape, 1.19)
    spacing = (0.4, 0.3, 0.2)
    dt = 0.1

    result = project_variable_density_mass_flux_3d(
        *fluxes,
        density_old_kg_m3=old,
        density_new_kg_m3=new,
        mass_source_rate_kg_m3_s=0.0,
        solid_mask=solid,
        spacing_m=spacing,
        time_step_s=dt,
    )

    assert result.domain_mass_change_rate_kg_s < 0.0
    assert result.net_boundary_outflow_rate_kg_s == pytest.approx(
        -result.domain_mass_change_rate_kg_s, abs=1.0e-10,
    )
    assert result.maximum_absolute_residual_after_kg_m3_s < 1.0e-9


def test_localized_mass_source_closes_cellwise_and_global_budget_around_obstacle():
    shape = (8, 9, 12)
    solid = np.zeros(shape, dtype=bool)
    solid[:5, 3:6, 5:8] = True
    fx, fy, fz = _zero_flux(shape)
    density = np.full(shape, 1.2)
    source = np.zeros(shape)
    source[1:3, 1:3, 2:4] = 0.08

    result = project_variable_density_mass_flux_3d(
        fx, fy, fz,
        density_old_kg_m3=density,
        density_new_kg_m3=density,
        mass_source_rate_kg_m3_s=source,
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.02,
    )

    assert result.integrated_source_rate_kg_s > 0.0
    assert result.net_boundary_outflow_rate_kg_s == pytest.approx(
        result.integrated_source_rate_kg_s, abs=1.0e-10,
    )
    assert result.maximum_absolute_residual_after_kg_m3_s < 1.0e-9
    assert np.all(result.mass_flux_x_kg_m2_s[:, :, 1:-1][
        solid[:, :, :-1] | solid[:, :, 1:]
    ] == 0.0)
    assert np.all(result.mass_flux_y_kg_m2_s[:, 1:-1, :][
        solid[:, :-1, :] | solid[:, 1:, :]
    ] == 0.0)
    assert np.all(result.mass_flux_z_kg_m2_s[1:-1, :, :][
        solid[:-1, :, :] | solid[1:, :, :]
    ] == 0.0)


def test_nonuniform_density_front_is_balanced_without_fabricated_source():
    shape = (6, 7, 11)
    solid = np.zeros(shape, dtype=bool)
    fx, fy, fz = _zero_flux(shape)
    old = np.full(shape, 1.2)
    new = old.copy()
    new[:, :, :5] -= 0.015
    new[:, :, 5:8] += 0.005

    result = project_variable_density_mass_flux_3d(
        fx, fy, fz,
        density_old_kg_m3=old,
        density_new_kg_m3=new,
        mass_source_rate_kg_m3_s=0.0,
        solid_mask=solid,
        spacing_m=(0.15, 0.20, 0.25),
        time_step_s=0.04,
    )

    assert result.maximum_absolute_residual_before_kg_m3_s > 0.1
    assert result.maximum_absolute_residual_after_kg_m3_s < 1.0e-8
    assert abs(result.global_mass_balance_residual_kg_s) < 1.0e-9


def test_cached_direct_and_cg_projection_are_equivalent():
    shape = (6, 7, 10)
    rng = np.random.default_rng(7)
    solid = np.zeros(shape, dtype=bool)
    solid[:4, 3:5, 5:7] = True
    fx, fy, fz = _zero_flux(shape)
    fx[:, :, 1:-1] = rng.normal(0.0, 0.02, fx[:, :, 1:-1].shape)
    fy[:, 1:-1, :] = rng.normal(0.0, 0.02, fy[:, 1:-1, :].shape)
    fz[1:-1, :, :] = rng.normal(0.0, 0.02, fz[1:-1, :, :].shape)
    fx[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    fy[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    fz[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    old = np.full(shape, 1.2)
    new = old.copy()
    new[1:3, 1:3, 2:4] -= 0.01
    kwargs = dict(
        density_old_kg_m3=old,
        density_new_kg_m3=new,
        mass_source_rate_kg_m3_s=0.0,
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.02,
    )

    direct = project_variable_density_mass_flux_3d(
        fx, fy, fz, linear_solver="cached_direct", **kwargs,
    )
    iterative = project_variable_density_mass_flux_3d(
        fx, fy, fz, linear_solver="cg", **kwargs,
    )
    amg = project_variable_density_mass_flux_3d(
        fx, fy, fz, linear_solver="amg_cg", **kwargs,
    )
    automatic = project_variable_density_mass_flux_3d(fx, fy, fz, **kwargs)

    assert direct.solver_method == "cached_direct"
    assert direct.solver_iterations == 1
    assert iterative.solver_method == "cg"
    assert iterative.solver_iterations > 1
    assert amg.solver_method == "amg_cg"
    assert amg.preconditioner == "symmetric_smoothed_aggregation_v_cycle"
    assert amg.preconditioner_level_count is not None
    assert amg.preconditioner_level_count >= 1
    assert amg.preconditioner_operator_complexity is not None
    assert amg.preconditioner_operator_complexity >= 1.0
    assert amg.pyamg_version is not None
    assert automatic.solver_method == "cached_direct"
    assert np.allclose(
        direct.pressure_correction_pa,
        iterative.pressure_correction_pa,
        rtol=1.0e-9, atol=1.0e-9,
    )
    assert np.allclose(
        direct.mass_flux_x_kg_m2_s,
        iterative.mass_flux_x_kg_m2_s,
        rtol=1.0e-9, atol=1.0e-10,
    )
    assert np.allclose(
        direct.mass_flux_x_kg_m2_s,
        amg.mass_flux_x_kg_m2_s,
        rtol=1.0e-9, atol=1.0e-10,
    )
    assert amg.maximum_absolute_residual_after_kg_m3_s < 1.0e-8

    warm = project_variable_density_mass_flux_3d(
        fx, fy, fz,
        linear_solver="cg",
        initial_correction_potential_kg_m_s=(
            iterative.correction_potential_kg_m_s
        ),
        **kwargs,
    )
    assert warm.solver_iterations <= iterative.solver_iterations
    assert np.allclose(
        warm.correction_potential_kg_m_s,
        iterative.correction_potential_kg_m_s,
        rtol=1.0e-9, atol=1.0e-10,
    )


def test_projection_rejects_solid_penetration():
    shape = (5, 5, 7)
    solid = np.zeros(shape, dtype=bool)
    solid[1:4, 1:4, 3:5] = True
    fx, fy, fz = _zero_flux(shape)
    fx[2, 2, 3] = 0.1
    density = np.ones(shape)
    with pytest.raises(ValueError, match="penetrates"):
        project_variable_density_mass_flux_3d(
            fx, fy, fz,
            density_old_kg_m3=density,
            density_new_kg_m3=density,
            mass_source_rate_kg_m3_s=0.0,
            solid_mask=solid,
            spacing_m=(1.0, 1.0, 1.0),
            time_step_s=0.1,
        )


def test_projection_rejects_fluid_pocket_disconnected_from_pressure_outlet():
    shape = (7, 7, 9)
    solid = np.zeros(shape, dtype=bool)
    centre = (3, 3, 3)
    for dk, dj, di in (
        (-1, 0, 0), (1, 0, 0), (0, -1, 0),
        (0, 1, 0), (0, 0, -1), (0, 0, 1),
    ):
        solid[centre[0] + dk, centre[1] + dj, centre[2] + di] = True
    density = np.ones(shape)
    with pytest.raises(ValueError, match="must reach"):
        project_variable_density_mass_flux_3d(
            *_zero_flux(shape),
            density_old_kg_m3=density,
            density_new_kg_m3=density,
            mass_source_rate_kg_m3_s=0.0,
            solid_mask=solid,
            spacing_m=(1.0, 1.0, 1.0),
            time_step_s=0.1,
        )


def test_six_pressure_open_faces_share_a_centered_source_and_close_mass():
    shape = (7, 7, 7)
    solid = np.zeros(shape, dtype=bool)
    density = np.ones(shape)
    source = np.zeros(shape)
    source[3, 3, 3] = 0.7

    result = project_variable_density_mass_flux_3d(
        *_zero_flux(shape),
        density_old_kg_m3=density,
        density_new_kg_m3=density,
        mass_source_rate_kg_m3_s=source,
        solid_mask=solid,
        spacing_m=(0.1, 0.1, 0.1),
        time_step_s=0.01,
        open_boundary_faces=(
            "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
        ),
    )

    face_outflows = np.asarray((
        -np.sum(result.mass_flux_x_kg_m2_s[:, :, 0]),
        np.sum(result.mass_flux_x_kg_m2_s[:, :, -1]),
        -np.sum(result.mass_flux_y_kg_m2_s[:, 0, :]),
        np.sum(result.mass_flux_y_kg_m2_s[:, -1, :]),
        -np.sum(result.mass_flux_z_kg_m2_s[0, :, :]),
        np.sum(result.mass_flux_z_kg_m2_s[-1, :, :]),
    )) * 0.1**2
    assert result.open_boundary_faces == (
        "x_max", "x_min", "y_max", "y_min", "z_max", "z_min",
    )
    assert np.max(face_outflows) - np.min(face_outflows) < 1.0e-13
    assert np.sum(face_outflows) == pytest.approx(
        result.integrated_source_rate_kg_s, abs=1.0e-12,
    )
    assert result.maximum_absolute_residual_after_kg_m3_s < 1.0e-10
    assert abs(result.global_mass_balance_residual_kg_s) < 1.0e-12


def test_component_may_reach_a_lateral_pressure_open_face_only():
    shape = (5, 7, 9)
    solid = np.ones(shape, dtype=bool)
    solid[2, :4, 2] = False
    density = np.ones(shape)

    result = project_variable_density_mass_flux_3d(
        *_zero_flux(shape),
        density_old_kg_m3=density,
        density_new_kg_m3=density,
        mass_source_rate_kg_m3_s=0.0,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.1,
        open_boundary_faces=("y_min",),
    )

    assert result.open_boundary_faces == ("y_min",)
    assert result.maximum_absolute_residual_after_kg_m3_s == 0.0


@pytest.mark.parametrize(
    "faces", [(), ("x_max", "x_max"), ("diagonal",)],
)
def test_projection_rejects_invalid_pressure_open_face_sets(faces):
    shape = (3, 3, 3)
    density = np.ones(shape)
    with pytest.raises(ValueError, match="open_boundary_faces"):
        project_variable_density_mass_flux_3d(
            *_zero_flux(shape),
            density_old_kg_m3=density,
            density_new_kg_m3=density,
            mass_source_rate_kg_m3_s=0.0,
            solid_mask=np.zeros(shape, dtype=bool),
            spacing_m=(1.0, 1.0, 1.0),
            time_step_s=0.1,
            open_boundary_faces=faces,
        )
