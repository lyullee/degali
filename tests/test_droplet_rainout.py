import math
from types import SimpleNamespace

import pytest

from degali.addons.droplet_rainout import (
    DropletClass,
    DropletTransportInput,
    concurrent_rainout_pool,
    droplet_transport_input_from_flash,
    post_release_rainout_pool,
    transport_droplet_population,
)
from degali.addons.pool_evaporation import SolidSubstrate


def _boundary():
    return DropletTransportInput(
        liquid_mass_flow_kg_s=0.2,
        classes=(
            DropletClass(diameter_m=1.0e-4, mass_fraction=0.25),
            DropletClass(diameter_m=1.0e-3, mass_fraction=0.75),
        ),
        liquid_density_kg_m3=70.0,
        initial_velocity_m_s=(0.0, 0.0, 0.0),
        release_position_m=(0.0, 0.0, 0.1),
        air_density_kg_m3=1.2,
        air_kinematic_viscosity_m2_s=1.5e-5,
        evaporation_coefficient_m2_s=1.0e-7,
        schmidt_number=0.7,
        wind_velocity_m_s=(5.0, 0.0, 0.0),
        latent_heat_j_kg=4.46e5,
        environmental_heat_input_w=5.0e4,
    )


def test_population_closes_mass_and_resolves_evaporation_and_rainout():
    result = transport_droplet_population(_boundary(), max_time_s=2.0)
    assert result.mass_residual_kg_s == pytest.approx(0.0, abs=1.0e-12)
    assert result.airborne_vapour_mass_flow_kg_s > 0.0
    assert result.ground_liquid_mass_flow_kg_s > 0.0
    assert {item.status for item in result.outcomes} == {
        "complete_evaporation_before_ground", "ground_impact"
    }
    assert result.impact_centroid_m is not None
    assert result.impact_centroid_m[0] > 0.0
    assert result.impact_centroid_m[2] == pytest.approx(0.0)
    assert result.evaporation_heat_requirement_w == pytest.approx(
        result.airborne_vapour_mass_flow_kg_s * 4.46e5
    )


def test_post_release_pool_preserves_end_to_end_hydrogen_inventory():
    droplets = transport_droplet_population(_boundary(), max_time_s=2.0)
    substrate = SolidSubstrate(
        conductivity_w_m_k=1.4,
        density_kg_m3=2200.0,
        heat_capacity_j_kg_k=850.0,
        initial_temperature_k=293.15,
        depth_m=0.2,
        cells=24,
    )
    coupled = post_release_rainout_pool(
        droplets,
        direct_vapour_mass_flow_kg_s=0.1,
        release_duration_s=2.0,
        pool_area_m2=0.5,
        pool_duration_s=10.0,
        pool_time_step_s=0.5,
        substrate=substrate,
    )
    assert coupled.pool is not None
    assert coupled.hydrogen_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)
    assert coupled.pool_vapour_mass_kg + coupled.remaining_pool_liquid_mass_kg \
        == pytest.approx(coupled.deposited_liquid_mass_kg)


def test_concurrent_rainout_forms_pool_during_release_and_closes_mass():
    droplets = transport_droplet_population(_boundary(), max_time_s=2.0)
    substrate = SolidSubstrate(
        conductivity_w_m_k=1.4, density_kg_m3=2200.0,
        heat_capacity_j_kg_k=850.0, initial_temperature_k=293.15,
        depth_m=0.2, cells=24,
    )
    coupled = concurrent_rainout_pool(
        droplets, direct_vapour_mass_flow_kg_s=0.1,
        release_duration_s=2.0, post_release_duration_s=8.0,
        pool_area_m2=0.5, pool_time_step_s=0.5, substrate=substrate,
    )
    assert coupled.pool is not None
    assert coupled.pool.total_liquid_inflow_kg == pytest.approx(
        coupled.deposited_liquid_mass_kg
    )
    assert coupled.hydrogen_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)


def test_population_requires_explicit_normalized_classes():
    boundary = _boundary()
    invalid = DropletTransportInput(
        **{**boundary.__dict__, "classes": (DropletClass(1.0e-3, 0.9),)}
    )
    with pytest.raises(ValueError, match="sum to one"):
        transport_droplet_population(invalid)


def test_time_limit_keeps_unfinished_droplets_as_airborne_liquid():
    result = transport_droplet_population(_boundary(), max_time_s=1.0e-5)
    assert result.airborne_liquid_mass_flow_kg_s > 0.0
    assert result.mass_residual_kg_s == pytest.approx(0.0, abs=1.0e-12)


def test_postflash_adapter_preserves_liquid_rate_and_vector_jet_velocity():
    source = SimpleNamespace(
        liquid_mass_flow=0.3,
        droplet_diameter=2.0e-4,
        liquid_density=70.0,
        postflash_velocity=12.0,
    )
    boundary = droplet_transport_input_from_flash(
        source,
        release_position_m=(0.0, 0.0, 1.0),
        jet_direction=(3.0, 4.0, 0.0),
        wind_velocity_m_s=(1.0, 0.0, 0.0),
        air_density_kg_m3=1.2,
        air_kinematic_viscosity_m2_s=1.5e-5,
        evaporation_coefficient_m2_s=1.0e-8,
        schmidt_number=0.7,
    )
    assert boundary.liquid_mass_flow_kg_s == pytest.approx(0.3)
    assert boundary.initial_velocity_m_s == pytest.approx((7.2, 9.6, 0.0))
    assert boundary.classes == (DropletClass(2.0e-4, 1.0),)
