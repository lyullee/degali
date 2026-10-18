import numpy as np
import pytest

from degali.addons.hydrogen_air_thermodynamics import (
    HydrogenAirActiveScalarClosure,
)


def test_ambient_state_is_neutral_and_energy_reference_is_exact():
    closure = HydrogenAirActiveScalarClosure()
    state = closure.state_from_conserved(0.0, 0.0)

    assert state.temperature_k == pytest.approx(closure.ambient_temperature_k)
    assert state.hydrogen_mole_fraction == pytest.approx(0.0)
    assert state.mixture_density_kg_m3 == pytest.approx(
        closure.ambient_pressure_pa * 28.9652e-3
        / (8.31446261815324 * closure.ambient_temperature_k),
    )
    assert float(state.buoyancy_acceleration_m_s2) == pytest.approx(0.0, abs=1.0e-14)


def test_cold_hydrogen_air_mixture_crosses_neutral_buoyancy_temperature():
    closure = HydrogenAirActiveScalarClosure()
    mole_fraction = 0.04
    neutral = float(closure.neutral_buoyancy_temperature_k(mole_fraction))
    air_moles = closure.ambient_air_density_kg_m3 / 28.9652e-3
    hydrogen_concentration = (
        mole_fraction / (1.0 - mole_fraction) * air_moles * 2.01588e-3
    )
    heat_capacity = (
        closure.ambient_air_density_kg_m3 * closure.dry_air_heat_capacity_j_kg_k
        + hydrogen_concentration * closure.hydrogen_heat_capacity_j_kg_k
    )

    dense = closure.state_from_conserved(
        hydrogen_concentration,
        heat_capacity * ((neutral - 5.0) - closure.ambient_temperature_k),
    )
    light = closure.state_from_conserved(
        hydrogen_concentration,
        heat_capacity * ((neutral + 5.0) - closure.ambient_temperature_k),
    )

    assert float(dense.hydrogen_mole_fraction) == pytest.approx(mole_fraction)
    assert float(dense.buoyancy_acceleration_m_s2) < 0.0
    assert float(light.buoyancy_acceleration_m_s2) > 0.0


def test_source_sensible_enthalpy_reconstructs_declared_temperature_for_h2_dominance():
    closure = HydrogenAirActiveScalarClosure()
    concentration = 100.0
    source_temperature = 80.0
    enthalpy = concentration * closure.source_specific_sensible_enthalpy_j_kg(
        source_temperature
    )
    state = closure.state_from_conserved(concentration, enthalpy)

    assert source_temperature < float(state.temperature_k) < source_temperature + 2.0
    assert float(state.hydrogen_mole_fraction) > 0.99
    assert float(state.buoyancy_acceleration_m_s2) > 0.0


def test_closure_flags_strong_density_departure_and_rejects_nonphysical_temperature():
    closure = HydrogenAirActiveScalarClosure(maximum_relative_density_departure=0.05)
    state = closure.state_from_conserved(np.array([0.0, 0.1]), np.array([0.0, 0.0]))

    assert bool(state.active_scalar_applicability_mask[0])
    assert not bool(state.active_scalar_applicability_mask[1])
    with pytest.raises(ValueError, match="outside closure bounds"):
        closure.state_from_conserved(0.0, -1.0e9)


def test_zero_gravity_is_valid_for_buoyancy_isolation():
    closure = HydrogenAirActiveScalarClosure(gravity_m_s2=0.0)
    state = closure.state_from_conserved(
        np.asarray([[[0.01]]]), np.asarray([[[0.0]]]),
    )

    assert np.all(state.buoyancy_acceleration_m_s2 == 0.0)
