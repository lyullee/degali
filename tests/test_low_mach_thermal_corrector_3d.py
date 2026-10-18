import numpy as np

from degali.addons.hydrogen_air_thermodynamics import (
    HydrogenAirActiveScalarClosure,
)
from degali.addons.low_mach_thermal_corrector_3d import (
    LOW_MACH_THERMAL_CORRECTOR_3D_SCHEMA,
    solve_low_mach_thermal_corrector_3d,
)


def _faces(shape):
    nz, ny, nx = shape
    u = np.full((nz, ny, nx + 1), 0.25)
    v = np.zeros((nz, ny + 1, nx))
    w = np.zeros((nz + 1, ny, nx))
    return u, v, w


def test_corrector_preserves_ambient_state_and_uniform_flow():
    shape = (6, 7, 10)
    solid = np.zeros(shape, dtype=bool)
    zeros = np.zeros(shape)
    u, v, w = _faces(shape)
    closure = HydrogenAirActiveScalarClosure()

    result = solve_low_mach_thermal_corrector_3d(
        zeros, zeros,
        hydrogen_source_increment_kg_m3=zeros,
        sensible_enthalpy_source_increment_j_m3=zeros,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.01,
        species_diffusivity_m2_s=(0.01, 0.01, 0.005),
        thermal_diffusivity_m2_s=(0.01, 0.01, 0.005),
        closure=closure,
        projection_linear_solver="cached_direct",
    )

    assert result.schema == LOW_MACH_THERMAL_CORRECTOR_3D_SCHEMA
    assert result.converged
    assert result.iterations == 2
    assert result.projection_iterations_history.shape == (2,)
    assert np.all(result.projection_iterations_history == 1)
    assert result.mass_projection.solver_method == "cached_direct"
    assert np.max(np.abs(result.hydrogen_mass_concentration_kg_m3)) == 0.0
    assert np.max(np.abs(result.sensible_enthalpy_density_j_m3)) == 0.0
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1.0e-10
    assert abs(result.hydrogen_mass_residual_kg) < 1.0e-14
    assert abs(result.sensible_enthalpy_residual_j) < 1.0e-12


def test_corrector_accepts_spatial_species_and_thermal_diffusivity():
    shape = (4, 5, 6)
    solid = np.zeros(shape, dtype=bool)
    zeros = np.zeros(shape)
    u, v, w = _faces(shape)
    species = np.linspace(0.001, 0.01, np.prod(shape)).reshape(shape)
    thermal = np.linspace(0.002, 0.02, np.prod(shape)).reshape(shape)

    result = solve_low_mach_thermal_corrector_3d(
        zeros, zeros,
        hydrogen_source_increment_kg_m3=zeros,
        sensible_enthalpy_source_increment_j_m3=zeros,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.01,
        species_diffusivity_m2_s=species,
        thermal_diffusivity_m2_s=thermal,
        closure=HydrogenAirActiveScalarClosure(),
        projection_linear_solver="cached_direct",
    )

    assert result.converged
    assert abs(result.hydrogen_mass_residual_kg) < 1.0e-14
    assert abs(result.sensible_enthalpy_residual_j) < 1.0e-12


def test_corrector_couples_local_cold_source_to_density_pressure_and_mass_flux():
    shape = (7, 8, 12)
    solid = np.zeros(shape, dtype=bool)
    solid[:5, 3:6, 6:8] = True
    zeros = np.zeros(shape)
    source = np.zeros(shape)
    source[1:3, 1:3, 2:4] = 2.0e-4
    closure = HydrogenAirActiveScalarClosure(
        maximum_relative_density_departure=0.10,
    )
    source_energy = source * closure.source_specific_sensible_enthalpy_j_kg(80.0)
    u, v, w = _faces(shape)
    u[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    v[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    w[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0

    result = solve_low_mach_thermal_corrector_3d(
        zeros, zeros,
        hydrogen_source_increment_kg_m3=source,
        sensible_enthalpy_source_increment_j_m3=source_energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.01,
        species_diffusivity_m2_s=0.01,
        thermal_diffusivity_m2_s=0.01,
        closure=closure,
        advection_scheme="muscl_van_leer",
        source_time_integration="strang",
        maximum_iterations=12,
        relative_density_tolerance=1.0e-7,
        relative_velocity_tolerance=1.0e-7,
    )

    assert result.converged
    assert np.min(result.thermodynamic_state.temperature_k[~solid]) < 293.15
    assert np.max(np.abs(result.pressure_correction_pa[~solid])) > 0.0
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1.0e-9
    assert abs(result.mass_projection.global_mass_balance_residual_kg_s) < 1.0e-9
    assert abs(result.hydrogen_mass_residual_kg) < 1.0e-12
    assert abs(result.sensible_enthalpy_residual_j) < 1.0e-8
    assert np.all(result.u_face_m_s[:, :, 1:-1][
        solid[:, :, :-1] | solid[:, :, 1:]
    ] == 0.0)


def test_corrector_uses_declared_lateral_pressure_open_faces():
    shape = (7, 7, 7)
    zeros = np.zeros(shape)
    source = np.zeros(shape)
    source[3, 3, 3] = 1.0e-4
    closure = HydrogenAirActiveScalarClosure()
    u = np.zeros((7, 7, 8))
    v = np.zeros((7, 8, 7))
    w = np.zeros((8, 7, 7))
    open_faces = (
        "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
    )

    result = solve_low_mach_thermal_corrector_3d(
        zeros, zeros,
        hydrogen_source_increment_kg_m3=source,
        sensible_enthalpy_source_increment_j_m3=(
            source * closure.source_specific_sensible_enthalpy_j_kg(100.0)
        ),
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=np.zeros(shape, dtype=bool),
        spacing_m=(0.1, 0.1, 0.1),
        time_step_s=0.01,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        maximum_iterations=12,
        open_boundary_faces=open_faces,
    )

    assert result.mass_projection.open_boundary_faces == (
        "x_max", "x_min", "y_max", "y_min", "z_max", "z_min",
    )
    assert np.any(np.abs(result.v_face_m_s[:, (0, -1), :]) > 0.0)
    assert np.any(np.abs(result.w_face_m_s[(0, -1), :, :]) > 0.0)
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1.0e-9
    assert abs(result.mass_projection.global_mass_balance_residual_kg_s) < 1.0e-10
    assert abs(result.hydrogen_mass_residual_kg) < 1.0e-12


def test_corrector_applies_prescribed_thermal_inlet_and_closes_energy_ledger():
    shape = (3, 3, 4)
    zeros = np.zeros(shape)
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _faces(shape)
    closure = HydrogenAirActiveScalarClosure()
    inlet_enthalpy = np.full((shape[0], shape[1]), 12_000.0)

    result = solve_low_mach_thermal_corrector_3d(
        zeros, zeros,
        hydrogen_source_increment_kg_m3=zeros,
        sensible_enthalpy_source_increment_j_m3=zeros,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(0.2, 0.2, 0.2),
        time_step_s=0.01,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        inlet_sensible_enthalpy_density_j_m3=inlet_enthalpy,
        maximum_iterations=12,
        relative_density_tolerance=1.0e-8,
        relative_velocity_tolerance=1.0e-8,
        open_boundary_faces=("x_min", "x_max"),
    )

    assert result.converged
    assert np.max(result.sensible_enthalpy_density_j_m3) > 0.0
    assert result.sensible_enthalpy_outflow_j < 0.0
    assert abs(result.sensible_enthalpy_residual_j) < 1.0e-11
