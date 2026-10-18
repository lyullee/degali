import numpy as np

from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.conservative_nested_handoff_3d import (
    nested_outer_inlet_frame_at,
    remap_conservative_nested_handoff_3d,
)
from degali.addons.nested_outer_binary_flow_3d import (
    NESTED_OUTER_BINARY_FLOW_3D_SCHEMA,
    advance_nested_outer_binary_flow_step_3d,
    initialise_nested_outer_binary_flow_3d,
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


def _constant_handoff():
    z_edges = np.arange(0.0, 4.0, 1.0)
    y_edges = np.arange(0.0, 5.0, 1.0)
    shape = (2, 3, 4)
    species = np.full(shape, 0.015)
    return remap_conservative_nested_handoff_3d(
        time_s=np.asarray((0.0, 1.0)),
        source_z_edges_m=z_edges,
        source_y_edges_m=y_edges,
        target_z_edges_m=z_edges,
        target_y_edges_m=y_edges,
        total_mass_flux_anomaly_kg_m2_s=0.2 * species,
        species_mass_flux_kg_m2_s=species,
        sensible_enthalpy_flux_w_m2=-2.0e5 * species,
        lateral_momentum_flux_anomaly_pa=np.full(shape, 0.04),
        vertical_momentum_flux_anomaly_pa=np.full(shape, 0.02),
        source_id="constant three-component test inlet",
    )


def test_nested_outer_step_applies_tangential_momentum_and_closes_ledgers():
    closure = _closure()
    handoff = _constant_handoff()
    solid = np.zeros((3, 4, 4), dtype=bool)
    state = initialise_nested_outer_binary_flow_3d(
        cells=(4, 4, 3), closure=closure, ambient_u_m_s=1.0,
        solid_mask=solid,
    )
    ambient_mass_flux = closure.ambient_density_kg_m3
    result = advance_nested_outer_binary_flow_step_3d(
        state,
        handoff=handoff,
        ambient_total_mass_flux_kg_m2_s=ambient_mass_flux,
        closure=closure,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.1,
        kinematic_viscosity_m2_s=0.0,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        projection_linear_solver="cached_direct",
        relative_density_tolerance=1.0e-8,
        relative_mass_flux_tolerance=1.0e-8,
    )

    assert result.schema == NESTED_OUTER_BINARY_FLOW_3D_SCHEMA
    assert result.diagnostics["three_component_inlet_momentum_applied"]
    assert np.all(result.inlet_frame.velocity_v_m_s > 0.0)
    assert np.all(result.inlet_frame.velocity_w_m_s > 0.0)
    assert np.all(result.momentum.provisional_v_m_s[:, :, 0] > 0.0)
    assert np.all(result.momentum.provisional_w_m_s[:, :, 0] > 0.0)
    assert np.allclose(
        result.binary_corrector.mass_projection.mass_flux_x_kg_m2_s[:, :, 0],
        result.inlet_frame.total_mass_flux_kg_m2_s,
        rtol=0.0,
        atol=1.0e-12,
    )
    assert result.binary_corrector.converged
    assert abs(result.binary_corrector.species_mass_residual_kg) < 1.0e-10
    assert abs(result.binary_corrector.sensible_enthalpy_residual_j) < 1.0e-6
    assert (
        result.binary_corrector.mass_projection
        .maximum_absolute_residual_after_kg_m3_s < 1.0e-7
    )
    assert np.max(np.abs(
        result.momentum.momentum_balance_residual_kg_m_s
    )) < 1.0e-12


def test_nested_outer_step_carries_projected_state_and_cumulative_ledgers():
    closure = _closure()
    handoff = _constant_handoff()
    solid = np.zeros((3, 4, 4), dtype=bool)
    state = initialise_nested_outer_binary_flow_3d(
        cells=(4, 4, 3), closure=closure, ambient_u_m_s=1.0,
        solid_mask=solid,
    )
    common = dict(
        handoff=handoff,
        ambient_total_mass_flux_kg_m2_s=closure.ambient_density_kg_m3,
        closure=closure,
        solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0),
        time_step_s=0.1,
        kinematic_viscosity_m2_s=0.0,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        projection_linear_solver="cached_direct",
    )
    first = advance_nested_outer_binary_flow_step_3d(state, **common)
    second = advance_nested_outer_binary_flow_step_3d(first.state, **common)

    assert second.state.time_s == 0.2
    assert second.state.correction_potential_kg_m_s is not None
    assert second.state.cumulative_prescribed_inlet_species_mass_kg == (
        first.state.cumulative_prescribed_inlet_species_mass_kg
        + second.binary_corrector.prescribed_inlet_species_mass_kg
    )
    assert abs(second.state.cumulative_species_mass_residual_kg) < 1.0e-9
    assert abs(second.state.cumulative_sensible_enthalpy_residual_j) < 1.0e-5


def test_nested_outer_step_accepts_a_lazy_inlet_provider():
    closure = _closure()
    handoff = _constant_handoff()
    ambient_mass_flux = closure.ambient_density_kg_m3

    class Provider:
        calls = 0

        def frame_at(self, *, time_s, ambient_total_mass_flux_kg_m2_s, closure):
            self.calls += 1
            return nested_outer_inlet_frame_at(
                handoff, time_s=time_s,
                ambient_total_mass_flux_kg_m2_s=ambient_total_mass_flux_kg_m2_s,
                closure=closure,
            )

    provider = Provider()
    solid = np.zeros((3, 4, 4), dtype=bool)
    state = initialise_nested_outer_binary_flow_3d(
        cells=(4, 4, 3), closure=closure, ambient_u_m_s=1.0,
        solid_mask=solid,
    )
    result = advance_nested_outer_binary_flow_step_3d(
        state, handoff=provider,
        ambient_total_mass_flux_kg_m2_s=ambient_mass_flux,
        closure=closure, solid_mask=solid,
        spacing_m=(1.0, 1.0, 1.0), time_step_s=0.1,
        kinematic_viscosity_m2_s=0.0,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        projection_linear_solver="cached_direct",
    )
    assert provider.calls == 1
    assert result.diagnostics["three_component_inlet_momentum_applied"]
