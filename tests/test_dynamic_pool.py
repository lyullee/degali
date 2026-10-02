import pytest

from degali.addons.dynamic_pool import (
    ConstantHeatFluxSurface,
    DynamicPoolNumerics,
    simulate_axisymmetric_spreading_pool,
)
from degali.addons.pool_evaporation import SolidSubstrate


def _substrate(temperature_k: float) -> SolidSubstrate:
    return SolidSubstrate(
        conductivity_w_m_k=1.4,
        density_kg_m3=2200.0,
        heat_capacity_j_kg_k=850.0,
        initial_temperature_k=temperature_k,
        depth_m=0.2,
        cells=24,
    )


def _numerics(**changes) -> DynamicPoolNumerics:
    values = dict(
        domain_radius_m=3.0,
        radial_step_m=0.02,
        maximum_time_step_s=0.01,
        source_radius_m=0.2,
        reported_front_depth_m=1.0e-5,
    )
    values.update(changes)
    return DynamicPoolNumerics(**values)


def test_isothermal_pool_spreads_and_conserves_all_inflow():
    result = simulate_axisymmetric_spreading_pool(
        _substrate(20.27),
        liquid_inflow_rate_kg_s=0.2,
        inflow_duration_s=2.0,
        duration_s=5.0,
        output_time_step_s=0.5,
        liquid_density_kg_m3=70.0,
        liquid_kinematic_viscosity_m2_s=2.0e-7,
        numerics=_numerics(),
    )
    radii = [step.reported_radius_m for step in result.steps]
    assert result.total_liquid_inflow_kg == pytest.approx(0.4)
    assert result.total_evaporated_mass_kg == pytest.approx(0.0)
    assert result.remaining_liquid_mass_kg == pytest.approx(0.4)
    assert max(radii) > result.numerics.source_radius_m
    assert all(right >= left for left, right in zip(radii[1:], radii[2:]))
    assert result.maximum_absolute_mass_residual_kg < 1.0e-12


def test_hot_substrate_evaporation_is_positive_and_depth_stays_nonnegative():
    result = simulate_axisymmetric_spreading_pool(
        _substrate(293.15),
        liquid_inflow_rate_kg_s=0.2,
        inflow_duration_s=2.0,
        duration_s=3.0,
        output_time_step_s=0.25,
        liquid_density_kg_m3=70.0,
        liquid_kinematic_viscosity_m2_s=2.0e-7,
        numerics=_numerics(),
    )
    assert result.total_evaporated_mass_kg > 0.0
    assert min(result.final_depth_m) >= 0.0
    assert result.maximum_absolute_mass_residual_kg < 1.0e-10


def test_domain_escape_is_reported_and_closes_the_mass_ledger():
    numerics = _numerics(
        domain_radius_m=0.5,
        radial_step_m=0.01,
        source_radius_m=0.1,
        reported_front_depth_m=1.0e-6,
    )
    result = simulate_axisymmetric_spreading_pool(
        _substrate(20.27),
        liquid_inflow_rate_kg_s=1.0,
        inflow_duration_s=1.0,
        duration_s=3.0,
        output_time_step_s=0.25,
        liquid_density_kg_m3=70.0,
        liquid_kinematic_viscosity_m2_s=0.0,
        numerics=numerics,
    )
    assert result.escaped_domain_liquid_mass_kg > 0.0
    final = result.steps[-1]
    assert (
        final.liquid_mass_kg + final.escaped_domain_liquid_mass_kg
    ) == pytest.approx(final.cumulative_liquid_inflow_kg, abs=1.0e-9)
    assert result.maximum_absolute_mass_residual_kg < 1.0e-9


def test_unresolved_annular_footprint_is_rejected():
    numerics = DynamicPoolNumerics(
        domain_radius_m=2.0,
        radial_step_m=0.2,
        source_inner_radius_m=0.01,
        source_radius_m=0.02,
    )
    with pytest.raises(ValueError, match="unresolved"):
        simulate_axisymmetric_spreading_pool(
            _substrate(20.27),
            liquid_inflow_rate_kg_s=0.1,
            inflow_duration_s=1.0,
            duration_s=2.0,
            output_time_step_s=0.5,
            liquid_density_kg_m3=70.0,
            liquid_kinematic_viscosity_m2_s=2.0e-7,
            numerics=numerics,
        )


def test_reported_front_is_stable_under_radial_grid_refinement():
    radii = []
    for radial_step in (0.04, 0.02, 0.01):
        result = simulate_axisymmetric_spreading_pool(
            _substrate(20.27),
            liquid_inflow_rate_kg_s=0.2,
            inflow_duration_s=2.0,
            duration_s=3.0,
            output_time_step_s=0.5,
            liquid_density_kg_m3=70.0,
            liquid_kinematic_viscosity_m2_s=2.0e-7,
            numerics=_numerics(radial_step_m=radial_step),
        )
        radii.append(result.steps[-1].reported_radius_m)
    assert abs(radii[1] - radii[2]) <= 0.02
    assert abs(radii[0] - radii[2]) / radii[2] < 0.05


def test_declared_constant_heat_flux_surface_closes_mass_and_evaporates():
    result = simulate_axisymmetric_spreading_pool(
        ConstantHeatFluxSurface(200_000.0, "synthetic water boundary"),
        liquid_inflow_rate_kg_s=0.35,
        inflow_duration_s=10.0,
        duration_s=12.0,
        output_time_step_s=1.0,
        liquid_density_kg_m3=70.0,
        liquid_kinematic_viscosity_m2_s=2.0e-7,
        numerics=_numerics(
            domain_radius_m=1.2,
            source_inner_radius_m=0.18,
            source_radius_m=0.22,
            reported_front_depth_m=7.5e-4,
        ),
    )
    assert result.total_evaporated_mass_kg > 0.0
    assert result.evaporation_momentum_closure \
        == "zero_radial_momentum_vapor"
    assert result.maximum_absolute_mass_residual_kg < 1.0e-10


def test_evaporation_momentum_closure_is_explicit_and_changes_only_dynamics():
    common = dict(
        liquid_inflow_rate_kg_s=0.42,
        inflow_duration_s=10.0,
        duration_s=12.0,
        output_time_step_s=1.0,
        liquid_density_kg_m3=70.0,
        liquid_kinematic_viscosity_m2_s=2.0e-7,
        numerics=_numerics(
            domain_radius_m=1.2,
            source_inner_radius_m=0.18,
            source_radius_m=0.22,
            reported_front_depth_m=7.5e-4,
        ),
        solid_heat_flux_multiplier=0.4,
    )
    substrate = SolidSubstrate(204.0, 2700.0, 879.0, 277.9, 0.2, 24)
    zero_vapour_momentum = simulate_axisymmetric_spreading_pool(
        substrate, **common,
    )
    liquid_carryoff = simulate_axisymmetric_spreading_pool(
        substrate, **common,
        evaporation_momentum_closure="liquid_velocity_carryoff",
    )
    assert zero_vapour_momentum.total_liquid_inflow_kg \
        == pytest.approx(liquid_carryoff.total_liquid_inflow_kg)
    assert max(step.reported_radius_m for step in zero_vapour_momentum.steps) \
        != max(step.reported_radius_m for step in liquid_carryoff.steps)
    assert zero_vapour_momentum.maximum_absolute_mass_residual_kg < 1.0e-10
    assert liquid_carryoff.maximum_absolute_mass_residual_kg < 1.0e-10
