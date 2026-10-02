"""Phase-domain and balance tests for the cryogenic-air add-on."""

import math

import pytest

from degali.addons.cryogenic_air import (
    air_saturation_pressure,
    advance_two_velocity_phase_step,
    cunningham_corrected_particle_relaxation_time,
    davies_cunningham_slip_correction,
    equilibrium_air_phase_split,
    n2o2_condensed_phase_scope,
    li2026_zone3,
    multiphase_hydrogen_evaporation_endpoint,
    multiphase_hydrogen_source_plane,
    minimum_heat_limited_sublimation_time,
    particle_relaxation_time,
    particle_terminal_velocity,
    particle_knudsen_number,
    ranz_marshall_transfer_number,
    relax_two_velocity_drag,
    relax_two_velocity_schiller_naumann_drag,
    schiller_naumann_relaxation_time,
    transported_condensed_air_source,
)
from degali.addons.notional import (
    energy_conserving_notional_nozzle,
    expand_measured_throat_to_ambient,
    isentropic_throat,
)
from degali.addons.lh2_droplets import (
    critical_droplet_diameter_for_phase_delay,
    critical_diffusivity_for_phase_delay,
    flashing_hydrogen_droplet_source,
    gasflow_phase_relaxation_coefficient,
    homogeneous_equilibrium_hydrogen_source,
    homogeneous_equilibrium_hydrogen_source_from_postflash,
    minimum_heat_limited_hydrogen_evaporation_time,
)
from degali.validation.hecht_panda import load_data as load_hecht_panda_data
from degali.validation.nearfield import hydrogen_jet


def test_flashing_droplet_source_reproduces_published_test_3_scale():
    """ESREL 2025 table 3.1: 0.606e-3 mm, quality 0.153."""
    pytest.importorskip("CoolProp")
    source = flashing_hydrogen_droplet_source(
        mass_flow=0.265,
        orifice_diameter=0.012,
        upstream_temperature=26.084,
        upstream_pressure=0.4e6,
        upstream_quality=0.078,
        ambient_temperature=293.15,
    )
    # The paper's thermodynamic property implementation gives 0.153.  The
    # current CoolProp value is 0.185; both preserve a predominantly liquid
    # post-flash source and the independently tabulated atomisation scale.
    assert source.postflash_quality == pytest.approx(0.153, abs=0.04)
    assert source.droplet_diameter == pytest.approx(0.606e-6, rel=0.03)
    assert source.shattered
    assert source.liquid_mass_flow > source.vapour_mass_flow
    assert source.mass_residual < 1e-12
    assert source.momentum_residual < 1e-12
    assert source.energy_residual < 1e-12


def test_measured_lh2_pipe_source_retains_residual_liquid():
    pytest.importorskip("CoolProp")
    common = dict(
        mass_flow=0.237153394213,
        orifice_diameter=0.012,
        upstream_temperature=20.22341379,
        upstream_pressure=101325.0 + 2.249465067105e5,
        ambient_temperature=288.15,
    )
    small = flashing_hydrogen_droplet_source(
        **common, droplet_size_coefficient=10.0
    )
    central = flashing_hydrogen_droplet_source(
        **common, droplet_size_coefficient=15.0
    )
    large = flashing_hydrogen_droplet_source(
        **common, droplet_size_coefficient=20.0
    )
    assert central.upstream_quality is None
    assert central.postflash_quality < 0.1
    assert central.liquid_mass_fraction > 0.9
    assert central.droplet_number_flux > 0.0
    assert central.droplet_diameter / small.droplet_diameter == pytest.approx(1.5)
    assert large.droplet_diameter / small.droplet_diameter == pytest.approx(2.0)
    assert central.mass_residual < 1e-12
    assert central.momentum_residual < 1e-12
    assert central.energy_residual < 1e-12


def test_flash_rejects_a_mass_flow_and_pressure_thrust_that_overextract_energy():
    pytest.importorskip("CoolProp")
    with pytest.raises(ValueError, match="mass-flow/throat-momentum"):
        flashing_hydrogen_droplet_source(
            mass_flow=0.2,
            orifice_diameter=0.004,
            upstream_temperature=80.0,
            upstream_pressure=20.0e6,
            ambient_temperature=293.0,
        )


def test_lh2_fastest_evaporation_lifetime_is_a_diameter_squared_bound():
    pytest.importorskip("CoolProp")
    common = dict(
        gas_temperature=288.15,
        droplet_temperature=20.2,
        gas_thermal_conductivity=0.025,
        gas_prandtl=0.71,
    )
    one = minimum_heat_limited_hydrogen_evaporation_time(
        diameter=1.0e-6, **common
    )
    two = minimum_heat_limited_hydrogen_evaporation_time(
        diameter=2.0e-6, **common
    )
    forced = minimum_heat_limited_hydrogen_evaporation_time(
        diameter=1.0e-6, particle_reynolds=100.0, **common
    )
    assert one > 0.0
    assert two / one == pytest.approx(4.0)
    assert forced < one
    assert math.isinf(minimum_heat_limited_hydrogen_evaporation_time(
        diameter=1.0e-6,
        gas_temperature=20.2,
        droplet_temperature=20.2,
        gas_thermal_conductivity=0.025,
        gas_prandtl=0.71,
    ))


def test_collective_lh2_evaporation_bound_starts_after_liquid_is_spent():
    pytest.importorskip("CoolProp")
    bound = homogeneous_equilibrium_hydrogen_source(
        mass_flow=0.26527468341375925,
        orifice_diameter=0.012,
        upstream_temperature=20.186769311178438,
        upstream_pressure=101325.0 + 1.760862291691e5,
        ambient_temperature=288.15,
        theta=0.0,
        y=0.5,
    )
    assert bound.postflash.liquid_mass_fraction > 0.98
    assert bound.formation_distance > 0.0
    assert bound.source.x == pytest.approx(bound.formation_distance)
    assert bound.source.y == pytest.approx(0.5)
    assert 0.0 < bound.source.mass_fraction < 1.0
    assert bound.source.temperature == pytest.approx(20.368903539, rel=1e-9)
    assert bound.hydrogen_mass_residual < 1e-12
    assert bound.total_mass_residual < 1e-12
    assert bound.momentum_residual < 1e-12
    assert bound.energy_residual < 1e-12


def test_supercritical_hydrogen_postflash_reaches_the_same_phase_safe_handoff():
    """A measured 80 K compressed-gas state must not masquerade as LH2 storage."""
    pytest.importorskip("CoolProp")
    postflash = flashing_hydrogen_droplet_source(
        mass_flow=0.35033891417451984,
        orifice_diameter=0.004,
        upstream_temperature=80.0,
        upstream_pressure=201.72256977519947e5,
        ambient_temperature=293.0,
    )
    source = homogeneous_equilibrium_hydrogen_source_from_postflash(
        postflash, ambient_temperature=293.0,
    )
    assert 0.0 < postflash.postflash_quality < 1.0
    assert postflash.liquid_mass_flow > 0.0
    assert source.formation_distance > 0.0
    assert source.source.temperature == pytest.approx(20.368903539, rel=1e-9)
    assert 0.0 < source.source.mass_fraction < 1.0
    assert max(
        source.hydrogen_mass_residual,
        source.total_mass_residual,
        source.momentum_residual,
        source.energy_residual,
    ) < 1e-12


def test_collective_equilibrium_endpoint_does_not_fit_droplet_size():
    pytest.importorskip("CoolProp")
    common = dict(
        mass_flow=0.12654769978438706,
        orifice_diameter=0.006,
        upstream_temperature=21.89929295627428,
        upstream_pressure=101325.0 + 4.276562515975e5,
        ambient_temperature=288.15,
    )
    small = homogeneous_equilibrium_hydrogen_source(
        **common, droplet_size_coefficient=10.0
    )
    large = homogeneous_equilibrium_hydrogen_source(
        **common, droplet_size_coefficient=20.0
    )
    assert large.postflash.droplet_diameter == pytest.approx(
        2.0 * small.postflash.droplet_diameter
    )
    assert large.source == small.source
    assert large.phase_plane == small.phase_plane


def test_para_hydrogen_source_is_explicit_consistent_and_conservative():
    pytest.importorskip("CoolProp")
    common = dict(
        mass_flow=0.26527468341375925,
        orifice_diameter=0.012,
        upstream_temperature=20.186769311178438,
        upstream_pressure=101325.0 + 1.760862291691e5,
        ambient_temperature=288.15,
    )
    normal = homogeneous_equilibrium_hydrogen_source(
        **common, hydrogen_species="Hydrogen"
    )
    para = homogeneous_equilibrium_hydrogen_source(
        **common, hydrogen_species="ParaHydrogen"
    )
    assert normal.postflash.hydrogen_species == "Hydrogen"
    assert para.postflash.hydrogen_species == "ParaHydrogen"
    assert para.phase_plane.endpoint.hydrogen_species == "ParaHydrogen"
    assert para.source.temperature == pytest.approx(20.271250661, rel=1e-9)
    assert para.postflash.upstream_density != pytest.approx(
        normal.postflash.upstream_density, rel=1e-4
    )
    assert para.phase_plane.endpoint.air_ratio != pytest.approx(
        normal.phase_plane.endpoint.air_ratio, rel=1e-4
    )
    assert max(
        para.hydrogen_mass_residual,
        para.total_mass_residual,
        para.momentum_residual,
        para.energy_residual,
    ) < 1e-12


def test_gasflow_phase_relaxation_uses_physical_sherwood_scaling():
    one = gasflow_phase_relaxation_coefficient(
        droplet_diameter=1.0e-6,
        diffusion_coefficient=1.0e-7,
        schmidt_number=1.0,
    )
    two = gasflow_phase_relaxation_coefficient(
        droplet_diameter=2.0e-6,
        diffusion_coefficient=1.0e-7,
        schmidt_number=1.0,
    )
    forced = gasflow_phase_relaxation_coefficient(
        droplet_diameter=1.0e-6,
        diffusion_coefficient=1.0e-7,
        schmidt_number=1.0,
        particle_reynolds=100.0,
    )
    assert one == pytest.approx(1.2e6)
    assert one / two == pytest.approx(4.0)
    assert forced > one
    recovered = critical_diffusivity_for_phase_delay(
        droplet_diameter=1.0e-6,
        delay_time=1.0 / one,
        schmidt_number=1.0,
    )
    assert recovered == pytest.approx(1.0e-7)
    recovered_diameter = critical_droplet_diameter_for_phase_delay(
        diffusion_coefficient=1.0e-7,
        delay_time=1.0 / one,
        schmidt_number=1.0,
    )
    assert recovered_diameter == pytest.approx(1.0e-6)


@pytest.mark.parametrize("coefficient", [9.99, 20.01])
def test_flashing_droplet_source_rejects_unpublished_size_coefficient(
    coefficient,
):
    pytest.importorskip("CoolProp")
    with pytest.raises(ValueError, match="published 10--20 range"):
        flashing_hydrogen_droplet_source(
            mass_flow=0.1,
            orifice_diameter=0.012,
            upstream_temperature=20.5,
            upstream_pressure=2.0e5,
            droplet_size_coefficient=coefficient,
        )


def _hecht_panda_data_or_skip():
    try:
        return load_hecht_panda_data()
    except FileNotFoundError:
        pytest.skip(
            "external Hecht-Panda benchmark data are not redistributed; "
            "configure the controlled local validation material to run this test"
        )


def test_hecht_panda_uses_final_journal_aggregate_fits():
    data = _hecht_panda_data_or_skip()
    assert data["active_benchmark"] == "journal_2019"
    assert data["published_fits"][
        "centerline_inverse_mass_fraction_slope"
    ] == pytest.approx(0.2771)
    assert data["published_fits"][
        "mass_fraction_half_width_slope_mm"
    ] == pytest.approx(0.07069)
    assert data["conference_2017_published_fits"][
        "centerline_inverse_mass_fraction_slope"
    ] == pytest.approx(0.2626)


def test_hecht_panda_condition_count_ambiguity_is_not_silently_filled():
    data = _hecht_panda_data_or_skip()
    assert len(data["conditions"]) == 9
    assert not any(
        condition["P_nozzle_bar_abs"] == 4.0
        and condition["T_nozzle_K"] == 45.0
        for condition in data["conditions"]
    )
    assert any("inclusion" in note and "unresolved" in note for note in data["notes"])


def test_hecht_panda_throat_density_supports_real_gas_normalization():
    from CoolProp.CoolProp import PropsSI

    data = _hecht_panda_data_or_skip()
    relative_errors = [
        PropsSI(
            "D",
            "T",
            condition["T_throat_K"],
            "P",
            condition["P_throat_bar_abs"] * 1.0e5,
            "Hydrogen",
        )
        / condition["rho_throat_kg_m3"]
        - 1.0
        for condition in data["conditions"]
    ]
    assert max(map(abs, relative_errors)) < 0.006


def test_nbs_solid_nitrogen_pressure_is_monotone_and_matches_table():
    pressures = [
        air_saturation_pressure("Nitrogen", temperature)
        for temperature in (20.0, 30.0, 40.0, 52.0, 60.0)
    ]
    assert pressures == sorted(pressures)
    assert pressures[3] == pytest.approx(5.7 * 133.322368, rel=0.005)
    assert pressures[4] == pytest.approx(47.2 * 133.322368, rel=0.005)


def test_nbs_solid_oxygen_pressure_matches_table():
    assert air_saturation_pressure("Oxygen", 50.74) == pytest.approx(
        0.291 * 133.322368, rel=1e-12
    )
    assert air_saturation_pressure("Oxygen", 39.41) == pytest.approx(
        0.00092 * 133.322368, rel=1e-12
    )


def test_li_zone3_rejects_liquid_saturation_below_nitrogen_triple_point():
    with pytest.raises(ValueError, match="below the .* nitrogen triple point"):
        li2026_zone3(
            hydrogen_flow=0.00059,
            station2_temperature=51.0,
            station2_velocity=599.48,
            station2_density=0.485,
            station2_diameter=0.00161,
            ambient_temperature=298.0,
        )


def test_li_zone3_reproduces_paper_test2_and_conservation_balances():
    result = li2026_zone3(
        hydrogen_flow=0.00059,
        station2_temperature=51.0,
        station2_velocity=599.48,
        station2_density=0.485,
        station2_diameter=0.00161,
        ambient_temperature=298.0,
        allow_subtriple_liquid_extrapolation=True,
    )

    # Fig. 4 is about 0.79 at 51 K.  The tighter number below follows from
    # the printed equations and Table-2 Station-2 values.
    assert result.hydrogen_mass_fraction == pytest.approx(0.7904812, rel=1e-6)
    assert result.velocity == pytest.approx(473.8777, rel=1e-6)
    assert result.gas_mass_fraction_sum == pytest.approx(1.0, abs=1e-12)
    assert result.used_subtriple_liquid_extrapolation

    # Equation (13): inlet H2 + lateral air = outlet gas + stationary LN2.
    assert 0.00059 + result.entrained_air_flow == pytest.approx(
        result.gas_flow + result.condensed_nitrogen_flow, rel=1e-12
    )
    # Equation (14), with v_LN2 = 0.
    assert 0.00059 * 599.48 == pytest.approx(
        result.gas_flow * result.velocity, rel=1e-12
    )
    # Equation (23), continuity at Station 3.
    area = math.pi * result.diameter**2 / 4.0
    assert result.density * result.velocity * area == pytest.approx(
        result.gas_flow, rel=1e-12
    )


def test_li_zone3_is_identity_above_its_condensation_threshold():
    result = li2026_zone3(
        hydrogen_flow=0.00331,
        station2_temperature=80.0,
        station2_velocity=924.48,
        station2_density=0.308,
        station2_diameter=0.00385,
    )
    assert result.hydrogen_mass_fraction == 1.0
    assert result.velocity == 924.48
    assert result.diameter == 0.00385
    assert result.entrained_air_flow == 0.0
    assert not result.used_subtriple_liquid_extrapolation


def test_equilibrium_split_freezes_both_bulk_air_species_at_lh2_temperature():
    split = equilibrium_air_phase_split(
        temperature=23.0,
        pressure=101325.0,
        hydrogen_flow=1.0,
        nitrogen_flow=3.025,
        oxygen_flow=0.975,
    )
    assert split.nitrogen_condensed_flow > 3.024999
    assert split.oxygen_condensed_flow > 0.974999
    assert split.nitrogen_partial_pressure == pytest.approx(
        air_saturation_pressure("Nitrogen", 23.0), rel=1e-12
    )
    assert split.oxygen_partial_pressure == pytest.approx(
        air_saturation_pressure("Oxygen", 23.0), rel=1e-12
    )
    assert split.nitrogen_gas_flow + split.nitrogen_condensed_flow == (
        pytest.approx(3.025, rel=1e-12)
    )
    assert split.oxygen_gas_flow + split.oxygen_condensed_flow == (
        pytest.approx(0.975, rel=1e-12)
    )
    scope = n2o2_condensed_phase_scope(split)
    assert scope.classification == "co_condensed_subtriple_n2o2_mixture_unsupported"


def test_equilibrium_split_makes_oxygen_the_last_bulk_air_phase_boundary():
    cold = equilibrium_air_phase_split(
        temperature=68.0, pressure=101325.0, hydrogen_flow=1.0,
        nitrogen_flow=3.025, oxygen_flow=0.975,
    )
    warm = equilibrium_air_phase_split(
        temperature=77.35, pressure=101325.0, hydrogen_flow=1.0,
        nitrogen_flow=3.025, oxygen_flow=0.975,
    )
    assert cold.nitrogen_condensed_flow == pytest.approx(0.0, abs=1e-12)
    assert cold.oxygen_condensed_flow > 0.1
    assert warm.condensed_flow == pytest.approx(0.0, abs=1e-12)


def test_n2o2_scope_does_not_promote_separate_condensates_to_a_mixture_eos():
    all_gas = equilibrium_air_phase_split(
        temperature=90.0, pressure=101325.0, hydrogen_flow=1.0,
        nitrogen_flow=3.025, oxygen_flow=0.975,
    )
    all_gas_scope = n2o2_condensed_phase_scope(all_gas)
    assert all_gas_scope.classification == "uncondensed_bulk_air"
    assert all_gas_scope.separate_pure_component_closure_usable

    oxygen_only = equilibrium_air_phase_split(
        temperature=68.0, pressure=101325.0, hydrogen_flow=1.0,
        nitrogen_flow=3.025, oxygen_flow=0.975,
    )
    oxygen_scope = n2o2_condensed_phase_scope(oxygen_only)
    assert oxygen_scope.classification == "single_component_condensate_bound"
    assert oxygen_scope.oxygen_condensed and not oxygen_scope.nitrogen_condensed
    assert not oxygen_scope.separate_pure_component_closure_usable

    mixed = equilibrium_air_phase_split(
        temperature=23.0, pressure=101325.0, hydrogen_flow=1.0,
        nitrogen_flow=3.025, oxygen_flow=0.975,
    )
    mixed_scope = n2o2_condensed_phase_scope(mixed)
    assert mixed_scope.classification == (
        "co_condensed_subtriple_n2o2_mixture_unsupported"
    )
    assert mixed_scope.mixture_thermodynamics_required
    assert 0.0 < mixed_scope.condensed_nitrogen_mole_fraction < 1.0
    assert not mixed_scope.separate_pure_component_closure_usable


def test_hydrogen_evaporation_endpoint_carries_its_mixed_phase_scope():
    endpoint = multiphase_hydrogen_evaporation_endpoint(
        storage_temperature=27.0, ambient_temperature=288.0,
    )
    assert endpoint.n2o2_phase_scope.classification == (
        "co_condensed_subtriple_n2o2_mixture_unsupported"
    )
    assert endpoint.n2o2_phase_scope.mixture_thermodynamics_required


def test_particle_response_bounds_span_carried_to_settling_regimes():
    sizes = (1e-6, 1e-5, 1e-4)
    relaxation = [
        particle_relaxation_time(
            diameter=d, particle_density=800.0, gas_viscosity=4e-6
        )
        for d in sizes
    ]
    settling = [
        particle_terminal_velocity(
            diameter=d, particle_density=800.0, gas_density=0.4,
            gas_viscosity=4e-6,
        )
        for d in sizes
    ]
    assert relaxation[1] / relaxation[0] == pytest.approx(100.0)
    assert relaxation[2] / relaxation[1] == pytest.approx(100.0)
    assert settling == sorted(settling)
    assert settling[0] < 0.001
    assert settling[2] > 0.5
    # Nonlinear drag must reduce the 100-um result below its Stokes estimate.
    assert settling[2] < relaxation[2] * 9.80665


def test_declared_rarefaction_slip_corrects_only_the_low_re_response_time():
    diameter = 1e-6
    mean_free_path = .05e-6
    knudsen = particle_knudsen_number(
        mean_free_path=mean_free_path, diameter=diameter,
    )
    correction = davies_cunningham_slip_correction(
        mean_free_path=mean_free_path, diameter=diameter,
    )
    stokes = particle_relaxation_time(
        diameter=diameter, particle_density=800.0, gas_viscosity=4e-6,
    )
    corrected = cunningham_corrected_particle_relaxation_time(
        diameter=diameter, particle_density=800.0, gas_viscosity=4e-6,
        mean_free_path=mean_free_path,
    )
    assert knudsen == pytest.approx(.1)
    assert correction == pytest.approx(
        1.0 + .1 * (1.257 + .4 * math.exp(-11.0)), rel=1e-15,
    )
    assert corrected == pytest.approx(stokes * correction, rel=1e-15)
    assert davies_cunningham_slip_correction(
        mean_free_path=.05e-6, diameter=1e-3,
    ) < 1.001
    with pytest.raises(ValueError, match="positive"):
        particle_knudsen_number(mean_free_path=0.0, diameter=diameter)


def test_schiller_naumann_relaxation_is_a_finite_re_local_stokes_limit():
    common = dict(
        diameter=1e-4, particle_density=800.0, gas_density=.4,
        gas_viscosity=4e-6,
    )
    stokes = particle_relaxation_time(
        diameter=common["diameter"],
        particle_density=common["particle_density"],
        gas_viscosity=common["gas_viscosity"],
    )
    at_rest = schiller_naumann_relaxation_time(
        **common, relative_velocity=0.0,
    )
    finite_re = schiller_naumann_relaxation_time(
        **common, relative_velocity=20.0,
    )
    assert at_rest == pytest.approx(stokes, rel=1e-15)
    assert finite_re < stokes
    # Reynolds number here is 200: use the declared finite-Re multiplier,
    # rather than silently applying a terminal-velocity relation.
    assert finite_re == pytest.approx(
        stokes / (1.0 + .15 * 200.0**.687), rel=1e-14,
    )
    assert schiller_naumann_relaxation_time(
        **common, relative_velocity=-20.0,
    ) == pytest.approx(finite_re, rel=1e-14)

    with pytest.raises(ValueError, match="Reynolds"):
        schiller_naumann_relaxation_time(
            **common, relative_velocity=100.1,
        )
    with pytest.raises(ValueError, match="finite"):
        schiller_naumann_relaxation_time(
            **common, relative_velocity=math.nan,
        )


def test_exact_finite_re_two_velocity_drag_conserves_and_exceeds_stokes_drag():
    properties = dict(
        diameter=1e-4, particle_density=800.0, gas_density=.4,
        gas_viscosity=4e-6,
    )
    stokes_time = particle_relaxation_time(
        diameter=properties["diameter"],
        particle_density=properties["particle_density"],
        gas_viscosity=properties["gas_viscosity"],
    )
    finite_re = relax_two_velocity_schiller_naumann_drag(
        gas_mass_flow=2.0, particle_mass_flow=.5,
        gas_velocity=80.0, particle_velocity=60.0,
        duration=.02, **properties,
    )
    stokes = relax_two_velocity_drag(
        gas_mass_flow=2.0, particle_mass_flow=.5,
        gas_velocity=80.0, particle_velocity=60.0,
        duration=.02, particle_relaxation=stokes_time,
    )
    assert finite_re.momentum_after == pytest.approx(finite_re.momentum_before)
    assert finite_re.kinetic_energy_after + finite_re.thermalised_kinetic_energy == (
        pytest.approx(finite_re.kinetic_energy_before)
    )
    assert 60.0 < finite_re.particle_velocity < finite_re.gas_velocity < 80.0
    # At the declared initial Re=200, finite-Re drag reduces slip more than
    # the Stokes-only update over the same constant-property interval.
    assert (
        finite_re.gas_velocity - finite_re.particle_velocity
        < stokes.gas_velocity - stokes.particle_velocity
    )

    with pytest.raises(ValueError, match="Reynolds"):
        relax_two_velocity_schiller_naumann_drag(
            gas_mass_flow=2.0, particle_mass_flow=.5,
            gas_velocity=160.1, particle_velocity=60.0,
            duration=.02, **properties,
        )


def test_two_velocity_drag_relaxation_conserves_momentum_and_heats_slip():
    result = relax_two_velocity_drag(
        gas_mass_flow=2.0, particle_mass_flow=.5,
        gas_velocity=80.0, particle_velocity=20.0,
        duration=.1, particle_relaxation=.05,
    )
    assert result.momentum_after == pytest.approx(result.momentum_before)
    assert result.kinetic_energy_after < result.kinetic_energy_before
    assert result.thermalised_kinetic_energy > 0.0
    assert result.kinetic_energy_after + result.thermalised_kinetic_energy == (
        pytest.approx(result.kinetic_energy_before)
    )
    assert 20.0 < result.particle_velocity < result.gas_velocity < 80.0


def test_two_velocity_drag_has_exact_single_phase_and_long_time_limits():
    no_particle = relax_two_velocity_drag(
        gas_mass_flow=2.0, particle_mass_flow=0.0,
        gas_velocity=80.0, particle_velocity=-500.0,
        duration=10.0, particle_relaxation=.01,
    )
    assert no_particle.gas_velocity == 80.0
    assert no_particle.particle_velocity == -500.0
    assert no_particle.thermalised_kinetic_energy == 0.0

    long_time = relax_two_velocity_drag(
        gas_mass_flow=2.0, particle_mass_flow=.5,
        gas_velocity=80.0, particle_velocity=20.0,
        duration=10.0, particle_relaxation=.01,
    )
    assert long_time.gas_velocity == pytest.approx(68.0)
    assert long_time.particle_velocity == pytest.approx(68.0)
    with pytest.raises(ValueError, match="relaxation"):
        relax_two_velocity_drag(
            gas_mass_flow=2.0, particle_mass_flow=.5,
            gas_velocity=80.0, particle_velocity=20.0,
            duration=.1, particle_relaxation=0.0,
        )


def test_two_velocity_phase_step_separates_entrainment_and_drag_heating():
    result = advance_two_velocity_phase_step(
        gas_mass_flow=2.0, particle_mass_flow=.5,
        gas_velocity=80.0, particle_velocity=20.0,
        entrained_gas_mass_flow=1.0, entrained_gas_velocity=0.0,
        final_particle_mass_flow=1.0,
        duration=.1, particle_relaxation=.05,
    )
    assert result.gas_mass_flow == pytest.approx(2.5)
    assert result.particle_mass_flow == pytest.approx(1.0)
    assert result.momentum_after == pytest.approx(result.momentum_before)
    assert result.entrainment_thermalised_kinetic_energy > 0.0
    assert result.phase_transfer_thermalised_kinetic_energy > 0.0
    assert result.drag_thermalised_kinetic_energy > 0.0
    assert result.total_thermalised_kinetic_energy == pytest.approx(
        result.entrainment_thermalised_kinetic_energy
        + result.phase_transfer_thermalised_kinetic_energy
        + result.drag_thermalised_kinetic_energy
    )
    assert result.kinetic_energy_after + result.total_thermalised_kinetic_energy == (
        pytest.approx(result.kinetic_energy_before)
    )


def test_two_velocity_phase_step_handles_re_evaporation_without_impulse():
    result = advance_two_velocity_phase_step(
        gas_mass_flow=2.0, particle_mass_flow=1.0,
        gas_velocity=60.0, particle_velocity=10.0,
        entrained_gas_mass_flow=.5, entrained_gas_velocity=0.0,
        final_particle_mass_flow=.25,
        duration=.0, particle_relaxation=.05,
    )
    assert result.gas_mass_flow == pytest.approx(3.25)
    assert result.particle_mass_flow == pytest.approx(.25)
    assert result.momentum_after == pytest.approx(result.momentum_before)
    assert result.drag_thermalised_kinetic_energy == 0.0
    assert result.total_thermalised_kinetic_energy == pytest.approx(
        result.entrainment_thermalised_kinetic_energy
        + result.phase_transfer_thermalised_kinetic_energy
    )
    with pytest.raises(ValueError, match="positive gas"):
        advance_two_velocity_phase_step(
            gas_mass_flow=2.0, particle_mass_flow=.5,
            gas_velocity=80.0, particle_velocity=20.0,
            entrained_gas_mass_flow=0.0, entrained_gas_velocity=0.0,
            final_particle_mass_flow=2.5,
            duration=.1, particle_relaxation=.05,
        )


def test_two_velocity_phase_step_handles_complete_re_evaporation_limit():
    result = advance_two_velocity_phase_step(
        gas_mass_flow=2.0, particle_mass_flow=1.0,
        gas_velocity=60.0, particle_velocity=10.0,
        entrained_gas_mass_flow=.5, entrained_gas_velocity=0.0,
        final_particle_mass_flow=0.0,
        duration=.1, particle_relaxation=.05,
    )
    # Re-evaporated mass joins the gas at the particle velocity.  With no
    # particle inventory left, the drag term must vanish rather than using the
    # stale particle velocity in a fictitious force.
    assert result.gas_mass_flow == pytest.approx(3.5)
    assert result.particle_mass_flow == 0.0
    assert result.gas_velocity == pytest.approx(
        result.momentum_before / result.gas_mass_flow
    )
    assert result.momentum_after == pytest.approx(result.momentum_before)
    assert result.phase_transfer_thermalised_kinetic_energy > 0.0
    assert result.drag_thermalised_kinetic_energy == 0.0
    assert result.total_thermalised_kinetic_energy == pytest.approx(
        result.entrainment_thermalised_kinetic_energy
        + result.phase_transfer_thermalised_kinetic_energy
    )
    assert result.kinetic_energy_after + result.total_thermalised_kinetic_energy == (
        pytest.approx(result.kinetic_energy_before)
    )


def test_two_velocity_phase_step_conserves_declared_coflow_entrainment():
    result = advance_two_velocity_phase_step(
        gas_mass_flow=2.0, particle_mass_flow=0.0,
        gas_velocity=80.0, particle_velocity=-1.0,
        entrained_gas_mass_flow=1.0, entrained_gas_velocity=30.0,
        final_particle_mass_flow=0.0,
        duration=.1, particle_relaxation=.05,
    )
    # The ambient stream supplies real incoming axial momentum/KE.  The
    # resolved loss is the inelastic two-gas mixing loss, not a spurious loss
    # from treating a 30 m/s coflow as a zero-velocity reservoir.
    assert result.momentum_before == pytest.approx(190.0)
    assert result.momentum_after == pytest.approx(190.0)
    assert result.entrained_gas_mass_flow == 1.0
    assert result.entrained_gas_velocity == 30.0
    assert result.gas_velocity == pytest.approx(190.0 / 3.0)
    assert result.entrainment_thermalised_kinetic_energy == pytest.approx(
        .5 * (2.0 / 3.0) * (80.0 - 30.0)**2
    )
    assert result.phase_transfer_thermalised_kinetic_energy == 0.0
    assert result.drag_thermalised_kinetic_energy == 0.0
    assert result.kinetic_energy_after + result.total_thermalised_kinetic_energy == (
        pytest.approx(result.kinetic_energy_before)
    )


@pytest.mark.parametrize(
    "gas,particle,ug,up,entrained,ue,final_particle",
    [
        (1.0, .3, 60.0, 10.0, .5, -20.0, .7),
        (2.0, 1.0, -20.0, 30.0, .7, 40.0, .2),
    ],
)
def test_two_velocity_phase_step_preserves_co_and_counterflow_ledgers(
    gas, particle, ug, up, entrained, ue, final_particle,
):
    result = advance_two_velocity_phase_step(
        gas_mass_flow=gas, particle_mass_flow=particle,
        gas_velocity=ug, particle_velocity=up,
        entrained_gas_mass_flow=entrained, entrained_gas_velocity=ue,
        final_particle_mass_flow=final_particle,
        duration=.04, particle_relaxation=.07,
    )
    incoming_momentum = gas * ug + particle * up + entrained * ue
    assert result.momentum_before == pytest.approx(incoming_momentum)
    assert result.momentum_after == pytest.approx(incoming_momentum)
    assert result.entrainment_thermalised_kinetic_energy >= 0.0
    assert result.phase_transfer_thermalised_kinetic_energy >= 0.0
    assert result.drag_thermalised_kinetic_energy >= 0.0
    assert result.kinetic_energy_after + result.total_thermalised_kinetic_energy == (
        pytest.approx(result.kinetic_energy_before)
    )


def test_two_velocity_phase_step_recovers_the_existing_single_gas_entrainment_limit():
    result = advance_two_velocity_phase_step(
        gas_mass_flow=2.0, particle_mass_flow=0.0,
        gas_velocity=80.0, particle_velocity=-1.0,
        entrained_gas_mass_flow=1.0, entrained_gas_velocity=0.0,
        final_particle_mass_flow=0.0,
        duration=.1, particle_relaxation=.05,
    )
    assert result.gas_mass_flow == pytest.approx(3.0)
    assert result.gas_velocity == pytest.approx(2.0 * 80.0 / 3.0)
    assert result.particle_velocity == -1.0
    assert result.phase_transfer_thermalised_kinetic_energy == 0.0
    assert result.drag_thermalised_kinetic_energy == 0.0
    assert result.total_thermalised_kinetic_energy == pytest.approx(
        result.entrainment_thermalised_kinetic_energy
    )


@pytest.mark.parametrize(
    "gas,particle,ug,up,entrained,final_particle,duration,tau",
    [
        (1.0, 0.0, 40.0, -20.0, .5, .2, .03, .2),
        (2.0, .7, -15.0, 30.0, .1, 1.1, .2, .04),
        (3.0, 1.5, 4.0, -8.0, 2.0, .3, .0, .3),
        (.2, .4, 150.0, 20.0, .8, 1.0, .4, .01),
    ],
)
def test_two_velocity_phase_step_has_nonnegative_partitioned_losses(
    gas, particle, ug, up, entrained, final_particle, duration, tau,
):
    result = advance_two_velocity_phase_step(
        gas_mass_flow=gas, particle_mass_flow=particle,
        gas_velocity=ug, particle_velocity=up,
        entrained_gas_mass_flow=entrained, entrained_gas_velocity=0.0,
        final_particle_mass_flow=final_particle,
        duration=duration, particle_relaxation=tau,
    )
    assert result.momentum_after == pytest.approx(result.momentum_before)
    assert result.entrainment_thermalised_kinetic_energy >= 0.0
    assert result.phase_transfer_thermalised_kinetic_energy >= 0.0
    assert result.drag_thermalised_kinetic_energy >= 0.0
    assert result.total_thermalised_kinetic_energy == pytest.approx(
        result.entrainment_thermalised_kinetic_energy
        + result.phase_transfer_thermalised_kinetic_energy
        + result.drag_thermalised_kinetic_energy,
        abs=1e-10,
    )


def test_ranz_marshall_transfer_number_has_conduction_and_forced_limits():
    assert ranz_marshall_transfer_number(
        particle_reynolds=0.0, transport_number=0.7
    ) == pytest.approx(2.0)
    assert ranz_marshall_transfer_number(
        particle_reynolds=100.0, transport_number=0.7
    ) > 2.0
    with pytest.raises(ValueError, match="Reynolds"):
        ranz_marshall_transfer_number(
            particle_reynolds=-1.0, transport_number=0.7
        )


def test_heat_limited_sublimation_is_a_fast_d_squared_lower_bound():
    common = dict(
        species="Nitrogen", gas_temperature=68.0,
        particle_temperature=20.4, gas_thermal_conductivity=0.04,
        gas_prandtl=0.7,
    )
    one = minimum_heat_limited_sublimation_time(
        **common, diameter=1.0e-6
    )
    ten = minimum_heat_limited_sublimation_time(
        **common, diameter=1.0e-5
    )
    forced = minimum_heat_limited_sublimation_time(
        **common, diameter=1.0e-5, particle_reynolds=100.0
    )
    assert ten / one == pytest.approx(100.0)
    assert forced < ten
    assert minimum_heat_limited_sublimation_time(
        **{**common, "gas_temperature": 20.4}, diameter=1.0e-6
    ) == math.inf


def test_multiphase_lh2_endpoint_closes_pressure_energy_and_component_mass():
    endpoint = multiphase_hydrogen_evaporation_endpoint(
        storage_temperature=27.0,
        ambient_temperature=288.0,
        ambient_pressure=101325.0,
    )
    assert endpoint.temperature == pytest.approx(20.36890354, rel=1e-7)
    assert endpoint.air_ratio == pytest.approx(0.70008134, rel=1e-7)
    assert endpoint.hydrogen_mass_fraction == pytest.approx(
        1.0 / (1.0 + endpoint.air_ratio), rel=1e-12
    )
    assert (
        endpoint.hydrogen_partial_pressure
        + endpoint.nitrogen_partial_pressure
        + endpoint.oxygen_partial_pressure
    ) == pytest.approx(101325.0, rel=1e-11)
    assert endpoint.nitrogen_gas_ratio < 3e-12
    assert endpoint.oxygen_gas_ratio < 1e-15
    assert endpoint.relative_energy_residual < 1e-12


def test_multiphase_endpoint_closes_total_energy_with_source_momentum():
    static = multiphase_hydrogen_evaporation_endpoint(
        storage_temperature=27.0,
        ambient_temperature=288.0,
        ambient_pressure=101325.0,
    )
    moving = multiphase_hydrogen_evaporation_endpoint(
        storage_temperature=27.0,
        ambient_temperature=288.0,
        ambient_pressure=101325.0,
        specific_momentum=100.0,
    )
    assert moving.air_ratio > static.air_ratio
    assert moving.specific_kinetic_energy > 0.0
    assert moving.relative_energy_residual < 1e-12

    stationary = multiphase_hydrogen_evaporation_endpoint(
        storage_temperature=27.0,
        ambient_temperature=288.0,
        ambient_pressure=101325.0,
        specific_momentum=100.0,
        particle_velocity_fraction=0.0,
    )
    assert stationary.air_ratio > moving.air_ratio
    assert stationary.specific_kinetic_energy > moving.specific_kinetic_energy
    assert stationary.relative_energy_residual < 1e-12


def test_hydrogen_hem_throat_is_a_local_mass_flux_maximum():
    throat = isentropic_throat(
        fluid="Hydrogen", storage_temperature=27.0,
        ambient_pressure=101325.0,
    )
    from CoolProp.CoolProp import PropsSI

    for factor in (0.999, 1.001):
        pressure = throat.pressure * factor
        enthalpy = PropsSI(
            "H", "P", pressure, "S", throat.storage_entropy, "Hydrogen"
        )
        density = PropsSI(
            "D", "P", pressure, "S", throat.storage_entropy, "Hydrogen"
        )
        velocity = math.sqrt(
            2.0 * (throat.storage_enthalpy - enthalpy)
        )
        assert density * velocity <= throat.mass_flux * (1.0 + 1.0e-8)
    assert throat.choked
    assert throat.relative_energy_residual < 1e-12


def test_subcooled_hydrogen_throat_uses_independent_temperature_and_pressure():
    """A measured pipe state must not be projected onto saturation.

    Trial 23's corrected TC3 and PT2 medians are representative of the
    PRESLHY source plane: about 20.2 K at 2.25 barg.  That is compressed
    liquid, not saturated hydrogen at either of those two measurements.
    """
    throat = isentropic_throat(
        fluid="Hydrogen",
        storage_temperature=20.2234,
        storage_pressure=326246.5,
        ambient_pressure=101325.0,
    )
    assert throat.storage_temperature == pytest.approx(20.2234)
    assert throat.storage_pressure == pytest.approx(326246.5)
    assert throat.pressure == pytest.approx(101325.0, rel=1e-8)
    assert not throat.choked
    assert throat.relative_energy_residual < 1e-12


def test_supercritical_hydrogen_throat_requires_explicit_pure_gas_scope():
    """Cold compressed H2 has a measured P/T state, not a saturated state."""
    with pytest.raises(ValueError, match="allow_supercritical_gas"):
        isentropic_throat(
            fluid="Hydrogen", storage_temperature=80.0,
            storage_pressure=20.0e6,
        )

    throat = isentropic_throat(
        fluid="Hydrogen", storage_temperature=80.0,
        storage_pressure=20.0e6, allow_supercritical_gas=True,
    )
    from CoolProp.CoolProp import PropsSI

    for factor in (0.999, 1.001):
        pressure = throat.pressure * factor
        enthalpy = PropsSI("H", "P", pressure, "S", throat.storage_entropy, "Hydrogen")
        density = PropsSI("D", "P", pressure, "S", throat.storage_entropy, "Hydrogen")
        velocity = math.sqrt(2.0 * (throat.storage_enthalpy - enthalpy))
        assert density * velocity <= throat.mass_flux * (1.0 + 1.0e-8)
    assert throat.choked
    assert throat.relative_energy_residual < 1e-12


def test_supercritical_ideal_throat_can_reach_the_notional_source_plane():
    throat = isentropic_throat(
        fluid="Hydrogen", storage_temperature=80.0,
        storage_pressure=20.0e6, allow_supercritical_gas=True,
    )
    area = math.pi * 0.004**2 / 4.0
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen", storage_temperature=80.0,
        storage_pressure=20.0e6,
        mass_flow=throat.mass_flux * area,
        orifice_diameter=0.004, allow_supercritical_gas=True,
    )
    assert nozzle.admissible
    assert nozzle.discharge_coefficient == pytest.approx(1.0, abs=1e-9)
    assert nozzle.relative_mass_residual < 1e-12
    assert nozzle.relative_momentum_residual < 1e-12
    assert nozzle.relative_energy_residual < 1e-12


def test_subcooled_hydrogen_notional_nozzle_is_admissible():
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen",
        storage_temperature=20.2234,
        storage_pressure=326246.5,
        mass_flow=0.236,
        orifice_diameter=0.012,
        ambient_pressure=101325.0,
    )
    assert nozzle.admissible
    assert 0.0 < nozzle.discharge_coefficient < 1.0
    assert nozzle.relative_mass_residual < 1e-12
    assert nozzle.relative_momentum_residual < 1e-12
    assert nozzle.relative_energy_residual < 1e-12


def test_measured_pressure_momentum_reaches_the_axisymmetric_source():
    rate = 0.2371533942133789
    jet, _initial = hydrogen_jet(
        rate=rate,
        diameter=0.012,
        wind=1.7,
        height=1.5,
        ambient_temperature=288.65,
        relative_humidity=53.7,
        corrections=True,
        storage_temperature=20.22341379150174,
        source_upstream_pressure_barg=2.249465067105,
    )
    source = jet.axisymmetric_source
    nozzle = jet.source_notional_nozzle
    total_flow = source.density * source.velocity * source.area
    momentum = source.density * source.velocity**2 * source.area
    assert source.mass_fraction * total_flow == pytest.approx(rate, rel=1e-12)
    assert momentum == pytest.approx(
        # The legacy DEGADIS PI constant changes the reconstructed area by
        # about three parts per trillion.
        rate * nozzle.specific_momentum, rel=1e-10
    )
    assert nozzle.admissible


def test_hydrogen_jet_spin_selector_preserves_normal_default():
    common = dict(
        rate=0.2371533942133789,
        diameter=0.012,
        wind=1.7,
        height=1.5,
        ambient_temperature=288.65,
        relative_humidity=53.7,
        corrections=True,
        storage_temperature=20.22341379150174,
        source_upstream_pressure_barg=2.249465067105,
    )
    default, default_state = hydrogen_jet(**common)
    normal, normal_state = hydrogen_jet(
        **common, hydrogen_spin_isomer="normal"
    )
    para, _para_state = hydrogen_jet(
        **common, hydrogen_spin_isomer="para"
    )
    assert normal.axisymmetric_source == default.axisymmetric_source
    assert normal_state == pytest.approx(default_state)
    assert para.axisymmetric_source != default.axisymmetric_source
    assert para.axisymmetric_source.temperature != pytest.approx(
        default.axisymmetric_source.temperature, rel=1e-4
    )


def test_energy_conserving_notional_nozzle_closes_three_invariants():
    throat = isentropic_throat(
        fluid="Hydrogen", storage_temperature=23.0,
        ambient_pressure=101325.0,
    )
    diameter = 0.01
    area = math.pi * diameter**2 / 4.0
    # A conventional 0.8 discharge coefficient keeps the test independent
    # of any field campaign's measured source rate.
    mass_flow = 0.8 * throat.mass_flux * area
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen", storage_temperature=23.0,
        mass_flow=mass_flow, orifice_diameter=diameter,
        ambient_pressure=101325.0,
    )
    assert nozzle.admissible
    assert nozzle.discharge_coefficient == pytest.approx(0.8, rel=1e-9)
    assert nozzle.pressure_thrust > 0.0
    assert nozzle.relative_mass_residual < 1e-12
    assert nozzle.relative_momentum_residual < 1e-12
    assert nozzle.relative_energy_residual < 1e-12


def test_measured_throat_expansion_closes_three_invariants():
    diameter = 0.001
    density = 0.55
    velocity = 544.5
    mass_flow = density * velocity * math.pi * diameter**2 / 4.0
    outlet = expand_measured_throat_to_ambient(
        fluid="Hydrogen",
        mass_flow=mass_flow,
        throat_diameter=diameter,
        throat_pressure=0.972e5,
        throat_temperature=43.5,
        throat_density=density,
        throat_velocity=velocity,
        ambient_pressure=101325.0,
    )
    # The reported throat pressure is slightly sub-ambient, so no artificial
    # negative pressure thrust is allowed.
    assert outlet.pressure_thrust == 0.0
    assert outlet.velocity == pytest.approx(velocity)
    # DEGADIS's historical PI constant differs from ``math.pi`` in the last
    # few parts per trillion; this remains far below the frozen 1e-8 limit.
    assert outlet.relative_input_mass_residual < 1e-10
    assert outlet.relative_mass_residual < 1e-12
    assert outlet.relative_momentum_residual < 1e-12
    assert outlet.relative_energy_residual < 1e-12


def test_multiphase_source_plane_closes_mass_momentum_and_continuity():
    plane = multiphase_hydrogen_source_plane(
        hydrogen_flow=0.1,
        orifice_diameter=0.01,
        orifice_density=61.0,
        storage_temperature=27.0,
        ambient_temperature=288.0,
    )
    area = math.pi * plane.diameter**2 / 4.0
    assert plane.total_flow == pytest.approx(
        plane.hydrogen_flow + plane.nitrogen_flow + plane.oxygen_flow,
        rel=1e-12,
    )
    assert plane.density * plane.velocity * area == pytest.approx(
        plane.total_flow, rel=1e-12
    )
    orifice_velocity = 0.1 / (61.0 * math.pi * 0.01**2 / 4.0)
    assert plane.total_flow * plane.velocity == pytest.approx(
        0.1 * orifice_velocity, rel=1e-12
    )


def test_multiphase_source_plane_can_close_storage_total_energy():
    legacy = multiphase_hydrogen_source_plane(
        hydrogen_flow=0.1,
        orifice_diameter=0.01,
        orifice_density=61.0,
        storage_temperature=27.0,
        ambient_temperature=288.0,
    )
    conserved = multiphase_hydrogen_source_plane(
        hydrogen_flow=0.1,
        orifice_diameter=0.01,
        orifice_density=61.0,
        storage_temperature=27.0,
        ambient_temperature=288.0,
        include_kinetic_energy=True,
    )
    assert conserved.endpoint.air_ratio > legacy.endpoint.air_ratio
    assert conserved.endpoint.relative_energy_residual < 1e-12
    assert conserved.total_flow * conserved.velocity == pytest.approx(
        conserved.momentum_flux, rel=1e-12
    )
    assert conserved.formation_distance > 0.0

    from CoolProp.CoolProp import PropsSI

    rho_ambient = PropsSI("D", "T", 288.0, "P", 101325.0, "Air")
    entrainment = 0.281 * math.sqrt(
        conserved.momentum_flux / rho_ambient
    )
    expected = (
        (1.0 - conserved.endpoint.hydrogen_mass_fraction)
        * conserved.endpoint.air_ratio * conserved.hydrogen_flow
        / (entrainment * rho_ambient)
    )
    assert conserved.formation_distance == pytest.approx(expected, rel=1e-12)

    stationary = multiphase_hydrogen_source_plane(
        hydrogen_flow=0.1,
        orifice_diameter=0.01,
        orifice_density=61.0,
        storage_temperature=27.0,
        ambient_temperature=288.0,
        include_kinetic_energy=True,
        particle_velocity_fraction=0.0,
    )
    assert stationary.particle_velocity == 0.0
    assert stationary.gas_velocity > stationary.velocity
    assert stationary.momentum_flux == pytest.approx(
        stationary.gas_velocity
        * stationary.hydrogen_flow
        * (
            1.0 + stationary.endpoint.nitrogen_gas_ratio
            + stationary.endpoint.oxygen_gas_ratio
        ),
        rel=1e-12,
    )
    assert stationary.kinetic_energy_flow == pytest.approx(
        stationary.hydrogen_flow
        * stationary.endpoint.specific_kinetic_energy,
        rel=1e-12,
    )


def test_stationary_condensate_bound_deposits_with_zero_momentum():
    plane = multiphase_hydrogen_source_plane(
        hydrogen_flow=0.1,
        orifice_diameter=0.01,
        orifice_density=61.0,
        storage_temperature=27.0,
        ambient_temperature=288.0,
        include_kinetic_energy=True,
        particle_velocity_fraction=0.0,
    )
    results = []
    for diameter in (1e-6, 1e-4):
        results.append(transported_condensed_air_source(
            hydrogen_flow=plane.hydrogen_flow,
            station2_temperature=plane.endpoint.temperature,
            station2_velocity=plane.velocity,
            station2_density=plane.density,
            station2_diameter=plane.diameter,
            particle_diameter=diameter,
            initial_nitrogen_flow=plane.nitrogen_flow,
            initial_oxygen_flow=plane.oxygen_flow,
            ambient_temperature=288.0,
            station2_kinetic_energy=plane.kinetic_energy_flow,
            stationary_condensate=True,
        ))
    first, second = results
    assert first.handoff_reached
    assert first.hydrogen_mass_fraction > 0.999999
    assert first.nitrogen_dropped_flow + first.oxygen_dropped_flow > (
        plane.nitrogen_flow + plane.oxygen_flow
    )
    assert first.maximum_mass_residual < 1e-12
    assert first.maximum_momentum_residual < 1e-12
    assert first.maximum_energy_residual < 1e-12
    # Particle size cannot affect the declared zero-velocity limit.
    assert first == second

    with pytest.raises(ValueError, match="two-velocity kinetic energy"):
        transported_condensed_air_source(
            hydrogen_flow=plane.hydrogen_flow,
            station2_temperature=plane.endpoint.temperature,
            station2_velocity=plane.velocity,
            station2_density=plane.density,
            station2_diameter=plane.diameter,
            particle_diameter=1e-6,
            initial_nitrogen_flow=plane.nitrogen_flow,
            initial_oxygen_flow=plane.oxygen_flow,
            ambient_temperature=288.0,
            stationary_condensate=True,
        )


def test_transported_condensed_air_particle_bound_conserves_and_settles_more():
    plane = multiphase_hydrogen_source_plane(
        hydrogen_flow=0.1,
        orifice_diameter=0.01,
        orifice_density=61.0,
        storage_temperature=27.0,
        ambient_temperature=288.0,
    )
    dropped = []
    for diameter in (1e-6, 1e-5, 1e-4):
        result = transported_condensed_air_source(
            hydrogen_flow=plane.hydrogen_flow,
            station2_temperature=plane.endpoint.temperature,
            station2_velocity=plane.velocity,
            station2_density=plane.density,
            station2_diameter=plane.diameter,
            particle_diameter=diameter,
            initial_nitrogen_flow=plane.nitrogen_flow,
            initial_oxygen_flow=plane.oxygen_flow,
            ambient_temperature=288.0,
        )
        assert result.handoff_reached
        assert result.distance < 2.0
        assert result.maximum_mass_residual < 1e-6
        assert result.maximum_momentum_residual < 1e-6
        assert result.maximum_energy_residual < 1e-6
        assert result.axial_kinematic_limit == "no_slip_retained_condensate_limit"
        assert not result.finite_particle_slip_modelled
        assert not result.quantitative_lh2_prediction_allowed
        dropped.append(
            result.nitrogen_dropped_flow + result.oxygen_dropped_flow
        )
    assert dropped == sorted(dropped)
    assert dropped[0] < 1e-3 * plane.total_flow
    assert dropped[2] > 0.25 * (plane.nitrogen_flow + plane.oxygen_flow)


def test_transported_source_is_wired_to_jet_only_as_an_explicit_option():
    pytest.importorskip("CoolProp")
    from degali.core.jetplume import J_X
    from degali.validation.nearfield import hydrogen_jet

    plume, initial = hydrogen_jet(
        rate=0.1,
        diameter=0.01,
        wind=5.0,
        height=0.5,
        ambient_temperature=288.0,
        corrections=True,
        condensed_air_particle_diameter=1e-6,
    )
    source = plume.condensed_air_source
    assert source.handoff_reached
    assert initial[J_X] > source.distance
    assert source.maximum_mass_residual < 1e-6
    assert source.maximum_momentum_residual < 1e-6
    assert source.maximum_energy_residual < 1e-6


def test_source_energy_and_thrust_options_require_transport_source():
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import hydrogen_jet

    common = dict(
        rate=0.1, diameter=0.01, wind=5.0, height=0.5,
        ambient_temperature=288.0, corrections=True,
    )
    with pytest.raises(ValueError, match="require the transported"):
        hydrogen_jet(**common, source_total_energy_consistency=True)
    with pytest.raises(ValueError, match="require the transported"):
        hydrogen_jet(**common, source_pressure_thrust=True)
    with pytest.raises(ValueError, match="require the transported"):
        hydrogen_jet(**common, evaporation_zone_distance=True)
    with pytest.raises(ValueError, match="require the transported"):
        hydrogen_jet(**common, condensed_air_stationary_bound=True)

    with pytest.raises(ValueError, match="requires total-energy"):
        hydrogen_jet(
            **common,
            condensed_air_particle_diameter=1e-6,
            condensed_air_stationary_bound=True,
        )


def test_stationary_condensate_bound_is_wired_as_an_explicit_option():
    pytest.importorskip("CoolProp")
    from degali.validation.nearfield import hydrogen_jet

    plume, _initial = hydrogen_jet(
        rate=0.1, diameter=0.01, wind=5.0, height=0.5,
        ambient_temperature=288.0, corrections=True,
        condensed_air_particle_diameter=1e-6,
        source_total_energy_consistency=True,
        condensed_air_stationary_bound=True,
    )
    source = plume.condensed_air_source
    assert source.handoff_reached
    assert source.hydrogen_mass_fraction > 0.999999
    assert source.maximum_mass_residual < 1e-6
    assert source.maximum_momentum_residual < 1e-6
    assert source.maximum_energy_residual < 1e-6


@pytest.mark.slow
def test_evaporation_zone_distance_moves_only_the_source_coordinate():
    pytest.importorskip("CoolProp")
    from degali.core.jetplume import J_X
    from degali.validation.nearfield import hydrogen_jet

    common = dict(
        rate=0.1, diameter=0.01, wind=5.0, height=0.5,
        ambient_temperature=288.0, corrections=True,
        condensed_air_particle_diameter=1e-6,
        source_total_energy_consistency=True,
    )
    ordinary, y_ordinary = hydrogen_jet(**common)
    shifted, y_shifted = hydrogen_jet(
        **common, evaporation_zone_distance=True
    )
    assert shifted.evaporation_zone_distance > 0.0
    assert y_shifted[J_X] - y_ordinary[J_X] == pytest.approx(
        shifted.evaporation_zone_distance, rel=1e-12
    )
    assert shifted.condensed_air_source == ordinary.condensed_air_source
