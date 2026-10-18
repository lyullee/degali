import numpy as np
import pytest

from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)


def _r12_closure() -> BinaryIdealGasMixtureClosure:
    return BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )


def test_ti021_composition_reproduces_effective_molecular_weight_and_density_ratio():
    closure = _r12_closure()
    state = closure.state_from_mole_fraction_and_temperature(
        0.32, closure.ambient_temperature_k,
    )

    assert float(state.mixture_molar_mass_kg_mol) == pytest.approx(
        58.3887e-3, rel=2.0e-4,
    )
    assert float(state.mixture_density_kg_m3 / closure.ambient_density_kg_m3) == (
        pytest.approx(2.01580, rel=2.0e-4)
    )
    assert float(state.relative_density_departure) > 1.0
    assert float(state.buoyancy_acceleration_m_s2) < 0.0


def test_conserved_species_and_enthalpy_inversion_recovers_binary_state():
    closure = _r12_closure()
    mole_fraction = np.asarray([0.0, 0.02, 0.32, 0.8])
    temperature = np.asarray([291.15, 275.0, 291.15, 320.0])
    source = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, temperature,
    )
    recovered = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        source.species_mass_concentration_kg_m3,
        source.sensible_enthalpy_density_j_m3,
    )

    assert recovered.species_mole_fraction == pytest.approx(mole_fraction, abs=2e-13)
    assert recovered.temperature_k == pytest.approx(temperature, rel=2e-13)
    assert recovered.mixture_density_kg_m3 == pytest.approx(
        source.mixture_density_kg_m3, rel=2e-13,
    )


def test_analytic_conserved_inversion_recovers_random_physical_states():
    closure = _r12_closure()
    rng = np.random.default_rng(20261013)
    mole_fraction = rng.uniform(0.0, 1.0, size=512)
    temperature = rng.uniform(
        closure.minimum_temperature_k, closure.maximum_temperature_k, size=512,
    )
    source = closure.state_from_mole_fraction_and_temperature(
        mole_fraction, temperature,
    )

    recovered = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        source.species_mass_concentration_kg_m3,
        source.sensible_enthalpy_density_j_m3,
    )

    assert recovered.species_mole_fraction == pytest.approx(mole_fraction, abs=2e-12)
    assert recovered.temperature_k == pytest.approx(temperature, rel=2e-12)
    assert recovered.mixture_density_kg_m3 == pytest.approx(
        source.mixture_density_kg_m3, rel=2e-12,
    )


def test_binary_closure_rejects_species_concentration_above_pure_gas_bound():
    closure = _r12_closure()
    pure = closure.state_from_mole_fraction_and_temperature(
        1.0, closure.ambient_temperature_k,
    )
    with pytest.raises(ValueError, match="pure-species bound"):
        closure.state_from_species_mass_concentration_and_sensible_enthalpy(
            float(pure.mixture_density_kg_m3) * 1.01, 0.0,
        )


def test_zero_species_is_ambient_carrier_without_fixed_inventory_approximation():
    closure = _r12_closure()
    state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        np.zeros((2, 3)), np.zeros((2, 3)),
    )

    assert np.all(state.species_mass_fraction == 0.0)
    assert np.all(state.species_mole_fraction == 0.0)
    assert state.mixture_density_kg_m3 == pytest.approx(
        np.full((2, 3), closure.ambient_density_kg_m3),
    )


def test_conserved_inversion_handles_pure_species_and_heated_carrier_limits():
    closure = _r12_closure()
    source = closure.state_from_mole_fraction_and_temperature(
        np.asarray([0.0, 1.0]), np.asarray([320.0, 320.0]),
    )

    recovered = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        source.species_mass_concentration_kg_m3,
        source.sensible_enthalpy_density_j_m3,
    )

    assert recovered.species_mole_fraction == pytest.approx([0.0, 1.0], abs=2e-13)
    assert recovered.temperature_k == pytest.approx([320.0, 320.0], rel=2e-13)


def test_specific_enthalpy_boundary_state_recovers_flux_ratio_temperature():
    closure = _r12_closure()
    mass_fraction = np.asarray((0.0, 0.2, 0.6))
    specific_enthalpy = np.asarray((0.0, -8_000.0, 12_000.0))

    state = closure.state_from_mass_fraction_and_specific_sensible_enthalpy(
        mass_fraction, specific_enthalpy,
    )

    represented = (
        state.sensible_enthalpy_density_j_m3
        / state.mixture_density_kg_m3
    )
    assert state.species_mass_fraction == pytest.approx(mass_fraction)
    assert represented == pytest.approx(specific_enthalpy, rel=1e-12, abs=1e-12)
