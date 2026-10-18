import numpy as np
import pytest

from degali.addons.binary_mixture_handoff_plane_3d import (
    BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
    extract_binary_mixture_handoff_plane_3d,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)


def _closure():
    return BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0749,
        ambient_temperature_k=291.15,
        ambient_pressure_pa=100_500.0,
        minimum_temperature_k=150.0,
    )


def test_uniform_interior_plane_recovers_all_declared_fluxes():
    closure = _closure()
    shape = (3, 4, 5)
    state = closure.state_from_mole_fraction_and_temperature(
        np.full(shape, 0.08), np.full(shape, 280.0),
    )
    solid = np.zeros(shape, dtype=bool)
    u_face = np.full((3, 4, 6), 2.0)
    v = np.full(shape, 0.3)
    w = np.full(shape, -0.1)
    ambient_total = closure.ambient_density_kg_m3 * 2.0
    frame = extract_binary_mixture_handoff_plane_3d(
        state,
        state.sensible_enthalpy_density_j_m3,
        u_face_m_s=u_face,
        cell_v_m_s=v,
        cell_w_m_s=w,
        solid_mask=solid,
        x_face_index=3,
        spacing_m=(1.0, 2.0, 3.0),
        time_s=4.0,
        time_step_s=0.1,
        species_diffusivity_m2_s=1.0e-5,
        thermal_diffusivity_m2_s=2.0e-5,
        ambient_total_mass_flux_kg_m2_s=ambient_total,
        advection_scheme="muscl_van_leer",
    )

    rho = state.mixture_density_kg_m3[0, 0, 0]
    mdot = rho * 2.0
    specific_h = (
        state.sensible_enthalpy_density_j_m3[0, 0, 0] / rho
    )
    assert frame.schema == BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA
    assert frame.total_mass_flux_anomaly_kg_m2_s == pytest.approx(
        mdot - ambient_total
    )
    assert frame.species_mass_flux_kg_m2_s == pytest.approx(
        mdot * state.species_mass_fraction[0, 0, 0]
    )
    assert frame.advective_species_mass_flux_kg_m2_s == pytest.approx(
        frame.species_mass_flux_kg_m2_s
    )
    assert frame.advective_sensible_enthalpy_flux_w_m2 == pytest.approx(
        frame.sensible_enthalpy_flux_w_m2
    )
    assert frame.sensible_enthalpy_flux_w_m2 == pytest.approx(mdot * specific_h)
    assert frame.lateral_momentum_flux_anomaly_pa == pytest.approx(mdot * 0.3)
    assert frame.vertical_momentum_flux_anomaly_pa == pytest.approx(mdot * -0.1)
    assert frame.diagnostics["backflow_face_count"] == 0
    assert not frame.diagnostics["clipping_or_posthoc_renormalisation"]
    packet=frame.dense_layer_packet
    assert packet['momentum_x_flux_pa'] == pytest.approx(2.*packet['mixture_mass_flux_kg_m2_s'])
    assert packet['momentum_y_flux_pa'] == pytest.approx(.3*packet['mixture_mass_flux_kg_m2_s'])
    assert packet['momentum_z_flux_pa'] == pytest.approx(-.1*packet['mixture_mass_flux_kg_m2_s'])


def test_interior_plane_retains_diffusion_and_zeroes_solid_faces():
    closure = _closure()
    shape = (2, 3, 4)
    mole_fraction = np.zeros(shape)
    mole_fraction[:, :, 2:] = 0.1
    temperature = np.full(shape, 291.15)
    temperature[:, :, 2:] = 280.0
    state = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, temperature,
    )
    solid = np.zeros(shape, dtype=bool)
    solid[0, 1, 2] = True
    u_face = np.zeros((2, 3, 5))
    zero = np.zeros(shape)
    frame = extract_binary_mixture_handoff_plane_3d(
        state,
        state.sensible_enthalpy_density_j_m3,
        u_face_m_s=u_face,
        cell_v_m_s=zero,
        cell_w_m_s=zero,
        solid_mask=solid,
        x_face_index=2,
        spacing_m=(1.0, 1.0, 1.0),
        time_s=0.0,
        time_step_s=0.1,
        species_diffusivity_m2_s=0.02,
        thermal_diffusivity_m2_s=0.03,
        ambient_total_mass_flux_kg_m2_s=1.0,
        advection_scheme="upwind",
    )

    fluid = np.ones((2, 3), dtype=bool)
    fluid[0, 1] = False
    assert np.all(frame.species_mass_flux_kg_m2_s[fluid] < 0.0)
    assert np.all(frame.advective_species_mass_flux_kg_m2_s[fluid] == 0.0)
    assert np.all(frame.sensible_enthalpy_flux_w_m2[fluid] > 0.0)
    assert np.all(frame.advective_sensible_enthalpy_flux_w_m2[fluid] == 0.0)
    assert frame.species_mass_flux_kg_m2_s[0, 1] == 0.0
    assert frame.sensible_enthalpy_flux_w_m2[0, 1] == 0.0

