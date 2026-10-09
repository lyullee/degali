from types import SimpleNamespace

import pytest
from CoolProp.CoolProp import PropsSI

import degali.lh2 as lh2
from degali.addons import SolidSubstrate, flashing_hydrogen_droplet_source
from degali.addons.lh2_droplets import direct_vapour_jet_source_from_flash


def _flash_source(*, vapour=0.1, liquid=0.2):
    pressure = 101325.0
    temperature = float(PropsSI("T", "P", pressure, "Q", 0, "Hydrogen"))
    quality = vapour / (vapour + liquid)
    return SimpleNamespace(
        hydrogen_species="Hydrogen",
        mass_flow=vapour + liquid,
        vapour_mass_flow=vapour,
        liquid_mass_flow=liquid,
        ambient_pressure=pressure,
        postflash_quality=quality,
        postflash_density=1.0,
        postflash_specific_enthalpy=1.0e5,
        postflash_temperature=temperature,
        postflash_velocity=20.0,
    )


def test_direct_flash_vapour_adapter_preserves_vapour_flux_and_momentum():
    source = _flash_source()
    gas = direct_vapour_jet_source_from_flash(
        source, ambient_temperature=295.0, theta=0.0, x=3.0, y=0.5
    )

    assert gas.fuel_mass_flow == pytest.approx(source.vapour_mass_flow)
    assert gas.mass_flow * gas.velocity == pytest.approx(
        source.vapour_mass_flow * source.postflash_velocity
    )
    assert gas.x == pytest.approx(3.0)
    assert gas.y == pytest.approx(0.5)
    assert gas.enthalpy_boundary is not None


def test_coupled_transient_closes_mass_and_exposes_unresolved_vapour(
    monkeypatch,
):
    source = _flash_source()
    pool_steps = (
        SimpleNamespace(
            elapsed_s=1.0, evaporated_mass_kg=0.04, wet_area_m2=0.2
        ),
        SimpleNamespace(
            elapsed_s=2.0, evaporated_mass_kg=0.10, wet_area_m2=0.3
        ),
    )
    coupling = SimpleNamespace(
        direct_vapour_mass_kg=0.20,
        airborne_droplet_vapour_mass_kg=0.08,
        airborne_liquid_mass_kg=0.04,
        deposited_liquid_mass_kg=0.28,
        pool_vapour_mass_kg=0.10,
        remaining_pool_liquid_mass_kg=0.15,
        escaped_pool_liquid_mass_kg=0.03,
        hydrogen_mass_residual_kg=0.0,
        pool=SimpleNamespace(steps=pool_steps),
    )
    droplets = SimpleNamespace(
        outcomes=(SimpleNamespace(terminal_time_s=0.5),),
        impact_centroid_m=(1.0, 2.0, 0.0),
    )
    phase = SimpleNamespace(
        accepted=True,
        pool_coupling=coupling,
        droplets=droplets,
        warnings=["phase qualification"],
    )
    gas = SimpleNamespace(
        accepted=True,
        handoff=SimpleNamespace(
            status="transition_ready", hydrogen_mass_kg=0.20
        ),
        puff=SimpleNamespace(
            maximum_relative_mass_residual=0.0,
            maximum_relative_hydrogen_residual=0.0,
        ),
    )
    calls = {}

    def fake_phase(_source, **kwargs):
        calls["phase"] = kwargs
        return phase

    def fake_gas(gas_source, **kwargs):
        calls["gas_source"] = gas_source
        calls["gas"] = kwargs
        return gas

    monkeypatch.setattr(lh2, "run_lh2_rainout_pool_research", fake_phase)
    monkeypatch.setattr(lh2, "run_lh2_finite_release_research", fake_gas)

    result = lh2.run_lh2_coupled_transient_research(
        source,
        release_duration_s=2.0,
        post_release_duration_s=3.0,
        puff_duration_s=4.0,
        release_position_m=(5.0, 6.0, 0.5),
        release_azimuth_rad=0.2,
        release_elevation_rad=0.0,
        wind_speed_m_s=2.0,
        wind_to_angle_rad=0.1,
        evaporation_coefficient_m2_s=1.0e-7,
        pool_area_m2=0.5,
        pool_time_step_s=0.1,
        substrate=object(),
    )

    assert calls["gas"]["source_lateral_offset_m"] == pytest.approx(6.0)
    assert calls["gas_source"].fuel_mass_flow == pytest.approx(0.1)
    assert result.released_hydrogen_mass_kg == pytest.approx(0.6)
    assert result.hydrogen_mass_residual_kg == pytest.approx(0.0)
    assert result.conservative
    assert result.accepted
    assert not result.atmospherically_complete
    assert result.unresolved_atmospheric_mass_kg == pytest.approx(0.18)
    pool_terms = [
        item for item in result.atmospheric_sources
        if item.mechanism == "pool_evaporation"
    ]
    assert sum(item.hydrogen_mass_kg for item in pool_terms) == pytest.approx(0.1)
    assert "atmospheric completion      : partial" in result.report()


def test_coupled_transient_refuses_to_hide_nonhorizontal_gas_launch():
    source = _flash_source()
    with pytest.raises(ValueError, match="horizontal"):
        lh2.run_lh2_coupled_transient_research(
            source,
            release_duration_s=1.0,
            post_release_duration_s=1.0,
            puff_duration_s=1.0,
            release_position_m=(0.0, 0.0, 0.5),
            release_azimuth_rad=0.0,
            release_elevation_rad=0.1,
            wind_speed_m_s=2.0,
            wind_to_angle_rad=0.0,
            evaporation_coefficient_m2_s=1.0e-7,
            pool_area_m2=0.5,
            pool_time_step_s=0.1,
            substrate=object(),
        )


@pytest.mark.slow
def test_real_flash_runs_end_to_end_through_both_transient_branches():
    source = flashing_hydrogen_droplet_source(
        mass_flow=0.265,
        orifice_diameter=0.012,
        upstream_temperature=26.084,
        upstream_pressure=0.4e6,
        upstream_quality=0.078,
        ambient_temperature=293.15,
    )
    substrate = SolidSubstrate(
        conductivity_w_m_k=1.4,
        density_kg_m3=2200.0,
        heat_capacity_j_kg_k=850.0,
        initial_temperature_k=293.15,
        depth_m=0.2,
        cells=12,
    )
    result = lh2.run_lh2_coupled_transient_research(
        source,
        release_duration_s=0.002,
        post_release_duration_s=0.1,
        puff_duration_s=0.05,
        release_position_m=(1.0, 2.0, 0.5),
        release_azimuth_rad=0.0,
        release_elevation_rad=0.0,
        wind_speed_m_s=2.5,
        wind_to_angle_rad=0.0,
        evaporation_coefficient_m2_s=1.0e-7,
        pool_area_m2=0.5,
        pool_time_step_s=0.05,
        substrate=substrate,
        gas_model_options={
            "maximum_nearfield_distance": 0.2,
            "radial_points": 41,
            "nearfield_maximum_step": 0.001,
            "nearfield_relative_tolerance": 2.0e-6,
            "crosswind_maximum_distance": 2.0,
            "crosswind_maximum_step": 0.02,
            "puff_time_step": 0.01,
        },
        phase_model_options={"maximum_droplet_time_s": 0.1},
    )

    assert result.conservative
    assert result.hydrogen_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)
    assert result.gas_dispersion.transition_position_m[1] == pytest.approx(2.0)
    assert result.gas_dispersion.puff.maximum_relative_hydrogen_residual < 1.0e-10
