import numpy as np
import pytest

from degali.addons.binary_mixture_corrector_3d import (
    BINARY_MIXTURE_CORRECTOR_3D_SCHEMA,
    solve_binary_mixture_initial_cloud_step_3d,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)


def _closure() -> BinaryIdealGasMixtureClosure:
    return BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )


def _zero_faces(shape):
    nz, ny, nx = shape
    return (
        np.zeros((nz, ny, nx + 1)),
        np.zeros((nz, ny + 1, nx)),
        np.zeros((nz + 1, ny, nx)),
    )


def test_uniform_dense_binary_mixture_is_an_exact_zero_flow_identity():
    closure = _closure()
    shape = (4, 4, 6)
    prescribed = closure.state_from_mole_fraction_and_temperature(0.32, 291.15)
    species = np.full(shape, float(prescribed.species_mass_concentration_kg_m3))
    energy = np.zeros(shape)
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)

    result = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.1,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
    )

    assert result.schema == BINARY_MIXTURE_CORRECTOR_3D_SCHEMA
    assert result.converged
    assert result.species_mass_concentration_kg_m3 == pytest.approx(species)
    assert np.all(result.sensible_enthalpy_density_j_m3 == 0.0)
    assert np.max(np.abs(result.u_face_m_s)) == 0.0
    assert result.species_mass_residual_kg == pytest.approx(0.0, abs=1e-14)
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1e-12


def test_dense_cloud_advection_preserves_species_and_total_continuity():
    closure = _closure()
    shape = (4, 4, 8)
    mole_fraction = np.zeros(shape)
    mole_fraction[:, :, 1:4] = 0.32
    prescribed = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, np.full(shape, 291.15),
    )
    species = prescribed.species_mass_concentration_kg_m3
    energy = prescribed.sensible_enthalpy_density_j_m3
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)
    u[:] = 0.08

    result = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.05,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        maximum_iterations=40,
        relative_density_tolerance=1e-8,
        relative_mass_flux_tolerance=1e-8,
    )

    assert result.converged
    assert np.min(result.species_mass_concentration_kg_m3) >= 0.0
    assert np.max(result.thermodynamic_state.species_mole_fraction) <= 1.0
    assert abs(result.species_mass_residual_kg) < 2e-12
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 2e-10
    assert abs(result.mass_projection.global_mass_balance_residual_kg_s) < 2e-10

    warm = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.05,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        maximum_iterations=40,
        relative_density_tolerance=1e-8,
        relative_mass_flux_tolerance=1e-8,
        initial_correction_potential_kg_m_s=(
            result.mass_projection.correction_potential_kg_m_s
        ),
    )
    assert warm.species_mass_concentration_kg_m3 == pytest.approx(
        result.species_mass_concentration_kg_m3, rel=1e-10, abs=1e-12,
    )
    assert np.sum(warm.projection_iterations_history) <= np.sum(
        result.projection_iterations_history
    )


def test_species_diffusion_uses_mass_fraction_gradient_and_is_conservative():
    closure = _closure()
    shape = (4, 4, 8)
    mole_fraction = np.zeros(shape)
    mole_fraction[1:3, 1:3, 2:5] = 0.32
    prescribed = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, np.full(shape, 291.15),
    )
    species = prescribed.species_mass_concentration_kg_m3
    energy = prescribed.sensible_enthalpy_density_j_m3
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)

    result = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.05,
        species_diffusivity_m2_s=0.02,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        maximum_iterations=40,
        relative_density_tolerance=1e-8,
        relative_mass_flux_tolerance=1e-8,
    )

    assert result.converged
    assert np.max(result.species_mass_concentration_kg_m3) < np.max(species)
    assert abs(result.species_mass_residual_kg) < 2e-12
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 2e-10


def test_binary_corrector_accepts_directional_diffusion_and_conserves_species_energy():
    closure = _closure()
    shape = (4, 4, 8)
    mole_fraction = np.zeros(shape)
    mole_fraction[1:3, :, :] = 0.32
    prescribed = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, np.full(shape, 291.15),
    )
    species = prescribed.species_mass_concentration_kg_m3.copy()
    energy = prescribed.sensible_enthalpy_density_j_m3.copy()
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)
    zero = np.zeros(shape)
    vertical = np.full(shape, 0.02)

    result = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.02,
        species_diffusivity_m2_s=(zero, zero, vertical),
        thermal_diffusivity_m2_s=(zero, zero, vertical),
        closure=closure,
        maximum_iterations=40,
        relative_density_tolerance=1e-8,
        relative_mass_flux_tolerance=1e-8,
    )

    assert result.converged
    assert result.species_mass_concentration_kg_m3[0, 1, 1] > species[0, 1, 1]
    assert abs(result.species_mass_residual_kg) < 2e-10
    assert abs(result.sensible_enthalpy_residual_j) < 2e-7
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 2e-8


def test_prescribed_outer_inlet_flux_closes_species_energy_and_total_continuity():
    closure = _closure()
    shape = (4, 4, 8)
    species = np.zeros(shape)
    energy = np.zeros(shape)
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)
    u[:] = 0.2
    total_flux = np.full(shape[:2], closure.ambient_density_kg_m3 * 0.2)
    species_flux = np.full(shape[:2], 0.01)
    enthalpy_flux = np.full(shape[:2], -1_000.0)
    dt = 0.02

    result = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=dt,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        maximum_iterations=40,
        relative_density_tolerance=1e-8,
        relative_mass_flux_tolerance=1e-8,
        open_boundary_faces=("x_max",),
        inlet_total_mass_flux_kg_m2_s=total_flux,
        inlet_species_mass_flux_kg_m2_s=species_flux,
        inlet_sensible_enthalpy_flux_w_m2=enthalpy_flux,
    )

    expected_species = float(np.sum(species_flux) * dt)
    expected_energy = float(np.sum(enthalpy_flux) * dt)
    assert result.converged
    assert result.prescribed_inlet_species_mass_kg == pytest.approx(expected_species)
    assert result.prescribed_inlet_sensible_enthalpy_j == pytest.approx(expected_energy)
    assert abs(result.species_mass_residual_kg) < 1e-10
    assert abs(result.sensible_enthalpy_residual_j) < 1e-7
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1e-8
    assert np.sum(result.species_mass_concentration_kg_m3) == pytest.approx(
        expected_species, abs=1e-10,
    )


def test_binary_corrector_preserves_species_and_heat_across_blocked_surface():
    closure = _closure()
    shape = (4, 4, 8)
    mole_fraction = np.zeros(shape)
    mole_fraction[:, :, :4] = 0.32
    temperature = np.full(shape, 280.0)
    temperature[:, :, :4] = 300.0
    initial = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, temperature,
    )
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)
    blocked_x = np.zeros_like(u, dtype=bool)
    blocked_y = np.zeros_like(v, dtype=bool)
    blocked_z = np.zeros_like(w, dtype=bool)
    blocked_x[:, :, 4] = True

    result = solve_binary_mixture_initial_cloud_step_3d(
        initial.species_mass_concentration_kg_m3,
        initial.sensible_enthalpy_density_j_m3,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.01,
        species_diffusivity_m2_s=0.02,
        thermal_diffusivity_m2_s=0.02,
        closure=closure,
        open_boundary_faces=("x_min", "x_max"),
        blocked_faces=(blocked_x, blocked_y, blocked_z),
    )

    assert result.converged
    assert result.species_mass_concentration_kg_m3 == pytest.approx(
        initial.species_mass_concentration_kg_m3,
    )
    assert result.sensible_enthalpy_density_j_m3 == pytest.approx(
        initial.sensible_enthalpy_density_j_m3,
    )
    assert np.all(result.mass_projection.mass_flux_x_kg_m2_s[blocked_x] == 0.0)
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1.0e-10
    assert abs(result.species_mass_residual_kg) < 1.0e-12
    assert abs(result.sensible_enthalpy_residual_j) < 1.0e-10


@pytest.mark.parametrize("source_time_integration", ("prestep", "strang"))
def test_binary_corrector_couples_species_and_energy_sources_to_total_projection(
    source_time_integration,
):
    closure = _closure()
    shape = (4, 4, 6)
    species = np.zeros(shape)
    energy = np.zeros(shape)
    source_species = np.zeros(shape)
    source_energy = np.zeros(shape)
    source_species[1, 1, 1] = 1.0e-3
    source_energy[1, 1, 1] = 25.0
    solid = np.zeros(shape, dtype=bool)
    u, v, w = _zero_faces(shape)

    result = solve_binary_mixture_initial_cloud_step_3d(
        species, energy,
        provisional_u_face_m_s=u,
        provisional_v_face_m_s=v,
        provisional_w_face_m_s=w,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.1,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        source_time_integration=source_time_integration,
        species_source_increment_kg_m3=source_species,
        total_mass_source_increment_kg_m3=20.0 * source_species,
        sensible_enthalpy_source_increment_j_m3=source_energy,
        maximum_iterations=40,
        relative_density_tolerance=1.0e-8,
        relative_mass_flux_tolerance=1.0e-8,
    )

    assert result.converged
    assert np.sum(result.species_mass_concentration_kg_m3) == pytest.approx(1.0e-3)
    assert np.sum(result.sensible_enthalpy_density_j_m3) == pytest.approx(25.0)
    assert result.species_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)
    assert result.sensible_enthalpy_residual_j == pytest.approx(0.0, abs=1.0e-9)
    assert result.mass_projection.maximum_absolute_residual_after_kg_m3_s < 1.0e-9
    assert abs(result.mass_projection.global_mass_balance_residual_kg_s) < 1.0e-9
