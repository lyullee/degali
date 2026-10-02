import pytest

from degali.addons.cryogenic_blowdown import (
    CryogenicBlowdownConfig,
    TransientPipeWall,
    TransientTankWall,
    blowdown_state_to_flashing_droplet_source,
    blowdown_state_to_homogeneous_evaporation_source,
    blowdown_state_to_lh2_source,
    hem_blowdown_state_to_flashing_droplet_source,
    hem_blowdown_state_to_homogeneous_evaporation_source,
    run_cryogenic_blowdown,
)


def _config(**overrides):
    values = dict(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=5.0e-4,
        initial_temperature_k=80.0,
        initial_pressure_pa=2.0e6,
        discharge_coefficient=0.7,
        duration_s=0.04,
        time_step_s=0.01,
    )
    values.update(overrides)
    return CryogenicBlowdownConfig(**values)


def test_adiabatic_blowdown_releases_mass_and_closes_its_discrete_energy_balance():
    result = run_cryogenic_blowdown(_config())

    assert result.states[0].mass_flow_kg_s > 0.0
    assert result.released_mass_kg > 0.0
    assert result.states[-1].pressure_pa < result.states[0].pressure_pa
    assert result.states[-1].temperature_k < result.states[0].temperature_k
    assert result.maximum_discrete_energy_residual_j < 1e-8


def test_prescribed_warm_wall_is_an_explicit_opt_in_heat_path():
    adiabatic = run_cryogenic_blowdown(_config())
    warmed = run_cryogenic_blowdown(
        _config(tank_wall_temperature_k=293.0, tank_wall_ua_w_k=30.0)
    )

    assert warmed.states[0].heat_transfer_w > 0.0
    assert warmed.states[-1].temperature_k > adiabatic.states[-1].temperature_k


def test_transient_wall_keeps_its_heat_storage_separate_from_a_lumped_ua():
    wall = TransientTankWall(
        thickness_m=0.02,
        area_m2=0.1,
        density_kg_m3=8000.0,
        specific_heat_j_kg_k=200.0,
        conductivity_w_m_k=9.0,
        initial_temperature_k=293.0,
        external_temperature_k=293.0,
        inner_heat_transfer_w_m2_k=20.0,
        outer_heat_transfer_w_m2_k=6.0,
    )
    result = run_cryogenic_blowdown(_config(transient_tank_wall=wall))

    assert result.states[0].tank_wall_inner_temperature_k == pytest.approx(293.0)
    assert result.states[0].heat_transfer_w > 0.0
    assert result.final_tank_wall_temperature_profile_k is not None
    assert result.final_tank_wall_temperature_profile_k[-1] < 293.0
    with pytest.raises(ValueError, match="either"):
        run_cryogenic_blowdown(_config(transient_tank_wall=wall, tank_wall_ua_w_k=1.0, tank_wall_temperature_k=293.0))


def test_prescribed_pipe_wall_changes_the_nozzle_inlet_but_not_tank_energy_balance():
    baseline = run_cryogenic_blowdown(_config())
    heated_pipe = run_cryogenic_blowdown(
        _config(
            pipe_wall_temperature_k=293.0,
            pipe_wall_ua_w_k=2.0,
            pipe_inner_diameter_m=0.01,
        )
    )

    assert heated_pipe.states[0].pipe_heat_transfer_w > 0.0
    assert (
        heated_pipe.states[0].nozzle_inlet_temperature_k
        > baseline.states[0].nozzle_inlet_temperature_k
    )
    assert heated_pipe.maximum_discrete_energy_residual_j < 1e-8


def test_transient_pipe_wall_is_an_explicit_thermal_store():
    pipe = TransientPipeWall(
        length_m=0.055,
        inner_diameter_m=0.01,
        thickness_m=0.001,
        density_kg_m3=8000.0,
        specific_heat_j_kg_k=200.0,
        conductivity_w_m_k=9.0,
        initial_temperature_k=293.0,
        external_temperature_k=293.0,
        inner_heat_transfer_w_m2_k=2.0,
        outer_heat_transfer_w_m2_k=6.0,
    )
    result = run_cryogenic_blowdown(_config(transient_pipe_wall=pipe))

    assert result.states[0].pipe_heat_transfer_w > 0.0
    assert result.final_pipe_wall_temperature_profile_k is not None
    assert result.final_pipe_wall_temperature_profile_k[-1] < 293.0
    with pytest.raises(ValueError, match="either"):
        run_cryogenic_blowdown(
            _config(
                transient_pipe_wall=pipe,
                pipe_wall_temperature_k=293.0,
                pipe_wall_ua_w_k=1.0,
                pipe_inner_diameter_m=0.01,
            )
        )


def test_single_phase_blowdown_snapshot_maps_to_a_conserved_ambient_source_plane():
    config = _config(initial_temperature_k=300.0)
    state = run_cryogenic_blowdown(config).states[0]

    expanded = blowdown_state_to_lh2_source(config, state, theta=0.0)

    assert expanded.source.mass_flow == pytest.approx(state.mass_flow_kg_s, rel=1e-8)
    assert expanded.source.theta == 0.0
    assert expanded.expansion.relative_input_mass_residual < 1e-8
    assert expanded.expansion.relative_momentum_residual < 1e-8
    assert expanded.expansion.relative_energy_residual < 1e-8


def test_e31_single_phase_tank_state_refuses_a_two_phase_ambient_gas_adapter():
    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
        duration_s=2.0,
        time_step_s=0.02,
    )
    state = run_cryogenic_blowdown(config).states[0]

    with pytest.raises(ValueError, match="ambient-pressure expansion is two-phase"):
        blowdown_state_to_lh2_source(config, state)


def test_e31_single_phase_tank_state_can_use_the_explicit_fast_evaporation_bound():
    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
        duration_s=2.0,
        time_step_s=0.02,
    )
    state = run_cryogenic_blowdown(config).states[0]
    homogeneous = blowdown_state_to_homogeneous_evaporation_source(config, state)

    assert 0.0 < homogeneous.postflash.postflash_quality < 1.0
    assert homogeneous.source.fuel_mass_flow == pytest.approx(state.mass_flow_kg_s)
    assert homogeneous.source.enthalpy_boundary is not None
    assert homogeneous.source.enthalpy_boundary.ambient_temperature == 295.0
    assert homogeneous.hydrogen_mass_residual < 1e-8
    assert homogeneous.momentum_residual < 1e-8
    assert homogeneous.energy_residual < 1e-8


def test_e31_flash_adapter_retains_the_liquid_inventory_before_evaporation():
    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
    )
    state = run_cryogenic_blowdown(config).states[0]
    flashing = blowdown_state_to_flashing_droplet_source(config, state)

    assert flashing.mass_flow == pytest.approx(state.mass_flow_kg_s)
    assert flashing.liquid_mass_flow > 0.0
    assert flashing.vapour_mass_flow > 0.0
    assert flashing.liquid_mass_flow + flashing.vapour_mass_flow == pytest.approx(
        state.mass_flow_kg_s
    )
    assert flashing.mass_residual < 1e-8
    assert flashing.momentum_residual < 1e-8
    assert flashing.energy_residual < 1e-8


def test_blowdown_snapshot_refuses_two_phase_and_terminal_states():
    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
        duration_s=1.2,
        time_step_s=0.02,
        two_phase_withdrawal="vapour",
        tank_internal_diameter_m=0.16,
        outlet_height_from_bottom_m=0.03,
    )
    result = run_cryogenic_blowdown(config)
    two_phase = next(state for state in result.states if state.vapour_quality is not None)
    with pytest.raises(ValueError, match="two-phase"):
        blowdown_state_to_lh2_source(config, two_phase)
    with pytest.raises(ValueError, match="positive-flow"):
        blowdown_state_to_lh2_source(config, result.states[-1])


def test_transient_pipe_wall_is_available_from_the_public_addons_namespace():
    from degali.addons import TransientPipeWall as public_pipe_wall

    assert public_pipe_wall is TransientPipeWall


def test_blowdown_rejects_an_undefined_wall_boundary_or_overunity_cd():
    with pytest.raises(ValueError, match="requires"):
        run_cryogenic_blowdown(_config(tank_wall_ua_w_k=1.0))
    with pytest.raises(ValueError, match="discharge_coefficient"):
        run_cryogenic_blowdown(_config(discharge_coefficient=1.01))
    with pytest.raises(ValueError, match="pipe_wall_temperature"):
        run_cryogenic_blowdown(_config(pipe_wall_ua_w_k=1.0))


def test_blowdown_stops_at_a_two_phase_tank_boundary_instead_of_choosing_an_outflow_phase():
    result = run_cryogenic_blowdown(
        _config(
            nozzle_diameter_m=4.0e-3,
            initial_pressure_pa=20.0e6,
            duration_s=2.0,
            time_step_s=0.02,
        )
    )

    assert result.termination == "two_phase_tank_boundary"
    assert result.states[-1].pressure_pa > result.config.ambient_pressure_pa


def test_equilibrium_vapour_withdrawal_requires_geometry_and_reports_phase_state():
    with pytest.raises(ValueError, match="requires tank_internal_diameter"):
        run_cryogenic_blowdown(
            _config(two_phase_withdrawal="vapour", duration_s=0.1)
        )

    result = run_cryogenic_blowdown(
        _config(
            nozzle_diameter_m=4.0e-3,
            initial_pressure_pa=20.0e6,
            duration_s=1.2,
            time_step_s=0.02,
            two_phase_withdrawal="vapour",
            tank_internal_diameter_m=0.16,
            outlet_height_from_bottom_m=0.03,
        )
    )

    assert result.termination in {"duration", "liquid_reaches_outlet"}
    two_phase = [state for state in result.states if state.vapour_quality is not None]
    assert two_phase
    assert two_phase[0].liquid_height_m is not None


def test_homogeneous_equilibrium_withdrawal_is_an_explicit_two_phase_option():
    from CoolProp.CoolProp import PropsSI

    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
        duration_s=1.2,
        time_step_s=0.02,
        two_phase_withdrawal="homogeneous",
    )
    result = run_cryogenic_blowdown(config)
    two_phase = [state for state in result.states if state.vapour_quality is not None]

    assert result.termination == "duration"
    assert two_phase
    assert two_phase[0].mass_flow_kg_s > 0.0
    assert two_phase[0].specific_enthalpy_j_kg == pytest.approx(
        PropsSI(
            "H", "D", two_phase[0].density_kg_m3,
            "U", two_phase[0].specific_internal_energy_j_kg, "Hydrogen",
        )
    )
    assert result.maximum_discrete_energy_residual_j < 1e-8


def test_hem_two_phase_state_maps_to_a_pressure_thrust_flash_plane():
    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
        duration_s=1.2,
        time_step_s=0.02,
        two_phase_withdrawal="homogeneous",
    )
    result = run_cryogenic_blowdown(config)
    state = next(item for item in result.states if item.vapour_quality is not None)
    flashing = hem_blowdown_state_to_flashing_droplet_source(config, state)

    assert flashing.mass_flow == pytest.approx(state.mass_flow_kg_s)
    assert 0.0 <= flashing.postflash_quality <= 1.0
    assert flashing.mass_residual < 1e-8
    assert flashing.momentum_residual < 1e-8
    assert flashing.energy_residual < 1e-8


def test_hem_two_phase_state_can_reach_the_explicit_fast_evaporation_bound():
    config = _config(
        nozzle_diameter_m=4.0e-3,
        initial_pressure_pa=20.0e6,
        duration_s=1.2,
        time_step_s=0.02,
        two_phase_withdrawal="homogeneous",
    )
    state = next(
        item for item in run_cryogenic_blowdown(config).states
        if item.vapour_quality is not None
    )
    homogeneous = hem_blowdown_state_to_homogeneous_evaporation_source(
        config, state
    )

    assert homogeneous.source.fuel_mass_flow == pytest.approx(state.mass_flow_kg_s)
    assert homogeneous.hydrogen_mass_residual < 1e-8
    assert homogeneous.momentum_residual < 1e-8
    assert homogeneous.energy_residual < 1e-8
