from types import SimpleNamespace

import pytest

from degali.addons.droplet_rainout import DropletClass
from degali.addons.dynamic_pool import DynamicPoolResult
from degali.addons.pool_evaporation import SolidSubstrate
from degali.lh2 import run_lh2_rainout_pool_research


def test_one_call_postflash_droplet_rainout_and_pool_path_closes_hydrogen():
    source = SimpleNamespace(
        mass_flow=0.3,
        vapour_mass_flow=0.1,
        liquid_mass_flow=0.2,
        droplet_diameter=1.0e-3,
        liquid_density=70.0,
        postflash_velocity=10.0,
    )
    substrate = SolidSubstrate(
        conductivity_w_m_k=1.4, density_kg_m3=2200.0,
        heat_capacity_j_kg_k=850.0, initial_temperature_k=293.15,
        depth_m=0.2, cells=24,
    )
    result = run_lh2_rainout_pool_research(
        source,
        release_duration_s=2.0,
        post_release_duration_s=8.0,
        release_position_m=(0.0, 0.0, 0.1),
        release_azimuth_rad=0.0,
        release_elevation_rad=0.0,
        wind_speed_m_s=5.0,
        wind_to_angle_rad=0.0,
        evaporation_coefficient_m2_s=1.0e-7,
        pool_area_m2=0.5,
        pool_time_step_s=0.5,
        substrate=substrate,
        droplet_classes=(
            DropletClass(1.0e-4, 0.25),
            DropletClass(1.0e-3, 0.75),
        ),
        maximum_droplet_time_s=2.0,
    )
    assert result.accepted
    assert result.pool_coupling.direct_vapour_mass_kg == pytest.approx(0.2)
    assert result.pool_coupling.hydrogen_mass_residual_kg \
        == pytest.approx(0.0, abs=1.0e-12)
    assert result.droplets.airborne_vapour_mass_flow_kg_s > 0.0
    assert result.droplets.ground_liquid_mass_flow_kg_s > 0.0
    assert isinstance(result.pool_coupling.pool, DynamicPoolResult)
    assert "execution screen            : pass" in result.report()
