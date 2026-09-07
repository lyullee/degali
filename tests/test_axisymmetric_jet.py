"""Tests for the conserved-energy axisymmetric near-field model."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.axisymmetric_jet import (
    AxisymmetricJetSource,
    ConservedGaussianJet,
    entrain_and_heat_initial_plug,
)
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons import axisymmetric_jet
from degali.core.constants import R_UNIVERSAL
from degali.core.jetplume import JetCoefficients
from degali.lh2 import (
    audit_lh2_independent_energy_interface,
    handoff_lh2_near_field_to_crosswind,
    lh2_source_from_measured_throat,
    run_lh2_crosswind_research,
    run_lh2_near_field_research,
)


def model(
    *, hydrogen_nonideal_volume_correction: bool = False,
) -> ConservedGaussianJet:
    source = AxisymmetricJetSource(
        diameter=9.81128e-4,
        velocity=544.5,
        density=0.571362,
        temperature=43.5243,
    )
    return ConservedGaussianJet(
        source,
        ambient_temperature=295.0,
        ambient_pressure=101325.0,
        ambient_density=1.197,
        fuel_molecular_weight=0.00201588,
        ambient_molecular_weight=0.02896546,
        fuel_heat_capacity=14294.8,
        ambient_heat_capacity=1006.2,
        radial_points=81,
        hydrogen_nonideal_volume_correction=hydrogen_nonideal_volume_correction,
    )


def test_plug_to_gaussian_establishment_is_physical():
    jet = model()
    assert not jet.equilibrium_air_condensation
    distance, state = jet.established_initial_state()
    assert distance == pytest.approx(6.2 * jet.source.diameter)
    assert state[0] == jet.source.velocity
    assert state[1] > 0.0
    assert 0.0 < state[3] < jet.source.mass_fraction
    assert state[6] == pytest.approx(distance)
    assert jet.centreline_temperature(state) > jet.source.temperature
    assert jet.centreline_temperature(state) < jet.ambient_temperature
    # The published algebra is retained for oracle reproduction, but its
    # boundary loss is explicit rather than mistaken for ODE conservation.
    assert jet.establishment_residuals["species"] < -0.05


def test_nonideal_volume_correction_is_optional_and_defaults_to_ideal():
    jet = model()
    density = 1.1
    fraction = 0.45
    expected = (
        jet.ambient_pressure
        / (R_UNIVERSAL * density)
        / (
            fraction / jet.fuel_molecular_weight
            + (1.0 - fraction) / jet._humid_ambient_molecular_weight
        )
    )

    temperature = jet._temperature_from_density(density, fraction)
    recomputed_density = jet._density_from_temperature(temperature, fraction)

    assert temperature == pytest.approx(expected, rel=1.0e-12)
    assert recomputed_density == pytest.approx(density, rel=1.0e-12)


def test_nonideal_volume_correction_activates_hydrogen_eos(monkeypatch):
    def fake_departure(*_args, **_kwargs):
        return SimpleNamespace(compressibility=0.88)

    monkeypatch.setattr(axisymmetric_jet, "hydrogen_gas_departure", fake_departure)
    jet = model(hydrogen_nonideal_volume_correction=True)
    density = 1.1
    fraction = 0.72
    molecular_weight = 1.0 / (
        fraction / jet.fuel_molecular_weight
        + (1.0 - fraction) / jet._humid_ambient_molecular_weight
    )
    ideal_temperature = jet.ambient_pressure * molecular_weight / (
        R_UNIVERSAL * density
    )
    corrected_temperature = jet._temperature_from_density(density, fraction)

    assert corrected_temperature > ideal_temperature
    assert jet._density_from_temperature(corrected_temperature, fraction) == pytest.approx(
        density, rel=1.0e-8
    )


def test_nonideal_volume_correction_ignored_when_disabled(monkeypatch):
    def forbidden_departure(*_args, **_kwargs):
        raise AssertionError("EOS should not be called")

    monkeypatch.setattr(axisymmetric_jet, "hydrogen_gas_departure", forbidden_departure)
    jet = model(hydrogen_nonideal_volume_correction=False)
    density = 1.1
    fraction = 0.4
    _ = jet._temperature_from_density(density, fraction)


def test_initial_entrainment_heating_closes_published_plug_balances():
    base = model()
    heated = entrain_and_heat_initial_plug(
        base.source,
        minimum_temperature=82.15083337712093,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
    )
    assert heated.source.temperature == pytest.approx(82.15083337712093)
    assert 0.0 < heated.source.mass_fraction < 1.0
    assert heated.source.mass_flow > heated.initial_mass_flow
    assert heated.source.velocity < base.source.velocity
    assert heated.source.diameter > 0.0
    assert heated.length > 0.0
    assert max(
        heated.fuel_mass_residual,
        heated.momentum_residual,
        heated.energy_residual,
        heated.geometry_mass_residual,
    ) < 1.0e-10


def test_scalar_peak_establishment_closes_all_three_invariants():
    base = model()
    jet = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment="scalar_peak",
    )
    _distance, state = jet.established_initial_state()
    assert min(state[:4]) > 0.0
    assert max(abs(value) for value in jet.establishment_residuals.values()) < 1e-10


def test_source_flux_establishment_does_not_add_development_zone_air():
    source = AxisymmetricJetSource(
        diameter=0.07966117529039163,
        velocity=30.75253187296913,
        density=2.8097694013401364,
        temperature=20.086235774350637,
        mass_fraction=0.43984781562378894,
        y=0.5,
    )
    common = dict(
        ambient_temperature=291.25,
        ambient_pressure=101325.0,
        ambient_density=1.219048590943567,
        fuel_molecular_weight=0.00201588,
        ambient_molecular_weight=0.02896546,
        fuel_heat_capacity=14294.8,
        ambient_heat_capacity=1006.2,
        radial_points=41,
    )
    source_flux = ConservedGaussianJet(
        source, conservative_establishment="source_flux", **common,
    )
    entrained = ConservedGaussianJet(
        source, conservative_establishment="entrained_mass", **common,
    )
    distance, source_state = source_flux.established_initial_state()
    _distance, entrained_state = entrained.established_initial_state()
    source_mass = source_flux._fluxes(source_state)[0]
    entrained_mass = entrained._fluxes(entrained_state)[0]
    assert distance == pytest.approx(6.2 * source.diameter)
    assert source_mass == pytest.approx(source.mass_flow, rel=1e-10)
    assert entrained_mass > source_mass
    assert max(
        abs(value) for value in source_flux.establishment_residuals.values()
    ) < 1e-8


def test_li_enthalpy_transport_is_explicit_and_fails_closed_at_bad_boundary():
    base = model()
    common = dict(
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment="source_flux",
    )
    total = ConservedGaussianJet(base.source, **common)
    enthalpy = ConservedGaussianJet(
        base.source, energy_transport="enthalpy", **common
    )
    with pytest.raises(ValueError, match="energy transport"):
        ConservedGaussianJet(base.source, energy_transport="turbulence", **common)

    assert total.energy_transport == "total"
    assert enthalpy.energy_transport == "enthalpy"
    _distance, state = enthalpy._published_initial_state()
    assert total._fluxes(state)[4] > enthalpy._fluxes(state)[4]
    with pytest.raises(
        RuntimeError, match="conservative Gaussian establishment did not close"
    ):
        enthalpy.established_initial_state()


def test_independent_energy_ground_bounds_are_geometric_and_fail_closed():
    coefficients = JetCoefficients()

    class CrosswindStub:
        th = None
        k = coefficients
        deltay = 0.0
        betay = 1.0
        deltaz = 0.0
        betaz = 1.0
        gammaz = 0.0
        ustar = 0.2

        @staticmethod
        def _split(sysz, _sya, _sza):
            width = math.sqrt(sysz)
            return width, width

        @staticmethod
        def _wind_profile(_z, _sz, _ct):
            return 2.5, 0.1

        @classmethod
        def _wind(cls, z, sz, ct):
            return cls._wind_profile(z, sz, ct)[0]

    thermodynamics = SimpleNamespace(
        ambient_density=1.2,
        spreading_ratio=1.16,
        _buoyancy_coefficient=0.05,
        plume_entrainment_limit=0.082,
    )
    common = dict(
        jetplume=CrosswindStub(),
        thermodynamics=thermodynamics,
        ground_interaction="free",
    )
    free = IndependentEnergyCrosswind(**common)
    geometry = IndependentEnergyCrosswind(
        common["jetplume"], thermodynamics, ground_interaction="geometry"
    )
    surface = IndependentEnergyCrosswind(
        common["jetplume"], thermodynamics, ground_interaction="surface_layer"
    )
    with pytest.raises(ValueError, match="ground interaction"):
        IndependentEnergyCrosswind(
            common["jetplume"], thermodynamics, ground_interaction="fitted"
        )

    contact = np.array([0.5, 0.1, 1.0, 0.1, 1.0, 2.0, 0.2])
    free_geometry = free._geometry(contact)
    contact_geometry = geometry._geometry(contact)
    assert 0.0 < contact_geometry[5] < 1.0
    assert contact_geometry[3] < free_geometry[3]
    assert contact_geometry[4] > 0.0
    assert geometry.source_terms(contact)[3] < free.source_terms(contact)[3]
    assert surface.source_terms(contact)[0] > geometry.source_terms(contact)[0]

    airborne = contact.copy()
    airborne[6] = 5.0
    assert geometry._geometry(airborne)[3:] == pytest.approx(
        free._geometry(airborne)[3:]
    )
    assert surface.source_terms(airborne) == pytest.approx(
        free.source_terms(airborne)
    )


def test_independent_energy_receptor_temperature_uses_local_profile():
    """A temperature comparison must use the sensor, not hidden centre state."""
    base = model()

    class CrosswindStub:
        th = None
        k = JetCoefficients()
        deltay = 0.0
        betay = 1.0
        deltaz = 0.0
        betaz = 1.0
        gammaz = 0.0
        ustar = 0.2

        @staticmethod
        def _split(sysz, _sya, _sza):
            width = math.sqrt(sysz)
            return width, width

        @staticmethod
        def _wind_profile(_z, _sz, _ct):
            return 2.5, 0.1

        @classmethod
        def _wind(cls, z, sz, ct):
            return cls._wind_profile(z, sz, ct)[0]

    crosswind = IndependentEnergyCrosswind(CrosswindStub(), base)
    state = np.array([0.8, 0.3, 0.25, 0.0, 5.0, 2.0, 1.0])
    centre = crosswind.point_temperature(state, 0.0, 1.0)
    edge = crosswind.point_temperature(state, 3.0, 1.0)
    assert centre < edge < base.ambient_temperature


def test_differential_scalar_spreading_requires_conservative_boundary():
    base = model()
    with pytest.raises(ValueError, match="requires a conservative establishment"):
        ConservedGaussianJet(
            base.source,
            ambient_temperature=base.ambient_temperature,
            ambient_pressure=base.ambient_pressure,
            ambient_density=base.ambient_density,
            fuel_molecular_weight=base.fuel_molecular_weight,
            ambient_molecular_weight=base.ambient_molecular_weight,
            fuel_heat_capacity=base.fuel_heat_capacity,
            ambient_heat_capacity=base.ambient_heat_capacity,
            thermodynamic_spreading_ratio=1.05,
        )


@pytest.mark.slow
def test_differential_scalar_spreading_conserves_all_fluxes():
    base = model()
    jet = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment=True,
        thermodynamic_spreading_ratio=1.16 * math.sqrt(0.70 / 0.85),
    )
    result = jet.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    assert max(abs(value) for value in jet.establishment_residuals.values()) < 1e-8
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )
    state = result.state_at_y(0.04)
    radii = np.linspace(0.0, 10.0 * state[1], 101)
    radial = np.array([jet.radial_state(state, radius) for radius in radii])
    assert np.all((radial[:, 0] >= 0.0) & (radial[:, 0] <= 1.0))
    assert np.all(radial[:, 1:] > 0.0)


@pytest.mark.slow
def test_independent_mass_fraction_temperature_conserves_all_fluxes():
    base = model()
    jet = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment=True,
        thermodynamic_spreading_ratio=1.16 * math.sqrt(0.70 / 0.85),
        thermodynamic_profile="temperature_mass_fraction",
    )
    result = jet.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    assert max(abs(value) for value in jet.establishment_residuals.values()) < 1e-8
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )
    state = result.state_at_y(0.04)
    assert jet.half_width(state, "temperature") == pytest.approx(
        jet.thermodynamic_spreading_ratio
        * state[1]
        * math.sqrt(math.log(2.0)),
        rel=1.0e-10,
    )
    radii = np.linspace(0.0, 10.0 * state[1], 101)
    radial = np.array([jet.radial_state(state, radius) for radius in radii])
    assert np.all((radial[:, 0] >= 0.0) & (radial[:, 0] <= 1.0))
    assert np.all(radial[:, 1:] > 0.0)


@pytest.mark.slow
def test_temperature_dependent_enthalpy_is_conservative():
    base = model()
    jet = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment=True,
        thermodynamic_spreading_ratio=1.16 * math.sqrt(0.70 / 0.85),
        thermodynamic_profile="temperature_mass_fraction",
        ideal_gas_enthalpy_fluids=("Hydrogen", "Air"),
    )
    variable_delta_h = (
        jet._mixture_enthalpy(295.0, 1.0)
        - jet._mixture_enthalpy(50.0, 1.0)
    )
    constant_delta_h = base.fuel_heat_capacity * (295.0 - 50.0)
    assert 0.8 < variable_delta_h / constant_delta_h < 1.0
    result = jet.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    assert max(abs(value) for value in jet.establishment_residuals.values()) < 1e-8
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )


@pytest.mark.slow
def test_equilibrium_condensation_is_an_explicit_conservative_bound():
    base = model()
    gas = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
    )
    phase = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
    )
    delayed = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
        equilibrium_dry_air_condensation=False,
    )
    _distance, gas_state = gas.established_initial_state()
    _distance, phase_state = phase.established_initial_state()
    _distance, delayed_state = delayed.established_initial_state()
    assert max(
        abs(value) for value in phase.establishment_residuals.values()
    ) < 1e-10
    assert phase.centreline_temperature(phase_state) > gas.centreline_temperature(
        gas_state
    )
    assert phase.centreline_temperature(
        phase_state
    ) > delayed.centreline_temperature(delayed_state)


@pytest.mark.slow
def test_humid_frost_adds_warming_without_breaking_flux_conservation():
    base = model()
    dry = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
    )
    humidity = 0.0165
    humid_mw = (1.0 + humidity) / (
        1.0 / base.ambient_molecular_weight + humidity / 0.01801528
    )
    humid_density = (
        base.ambient_pressure * humid_mw
        / (8.31446261815324 * base.ambient_temperature)
    )
    humid = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=humid_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
        ambient_absolute_humidity=humidity,
    )
    _distance, dry_state = dry.established_initial_state()
    _distance, humid_state = humid.established_initial_state()
    assert max(abs(v) for v in humid.establishment_residuals.values()) < 1e-8
    assert humid.centreline_temperature(humid_state) > dry.centreline_temperature(
        dry_state
    )
    result = humid.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=5.0e-6,
        method="LSODA",
    )
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )


@pytest.mark.slow
def test_measured_coflow_adds_ambient_momentum_and_energy_conservatively():
    base = model()
    coflow = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment=True,
        ambient_coflow_velocity=0.3,
    )
    result = coflow.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    assert max(abs(v) for v in coflow.establishment_residuals.values()) < 1e-8
    assert result.momentum_flux[-1] > result.momentum_flux[0]
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    relative_energy = (
        result.energy_flux - 0.5 * coflow.ambient_coflow_velocity**2
        * result.mass_flux
    )
    assert relative_energy[-1] == pytest.approx(
        relative_energy[0], rel=2.0e-4
    )


@pytest.mark.slow
def test_temperature_dependent_phase_enthalpy_is_physical_and_conservative():
    base = model()
    phase = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
        temperature_dependent_phase_enthalpy=True,
    )
    assert phase._mixture_enthalpy(50.0, 1.0) > (
        base._mixture_enthalpy(50.0, 1.0)
        - base._mixture_enthalpy(295.0, 1.0)
    )
    result = phase.solve(
        maximum_distance=0.04,
        maximum_step=0.0005,
        relative_tolerance=2.0e-7,
        method="LSODA",
    )
    assert max(abs(v) for v in phase.establishment_residuals.values()) < 1e-8
    assert np.all((result.temperature >= 14.1) & (result.temperature <= 400.0))
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )


@pytest.mark.slow
def test_para_hydrogen_enthalpy_is_explicit_and_conservative():
    base = model()
    normal = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
        temperature_dependent_phase_enthalpy=True,
    )
    para = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
        temperature_dependent_phase_enthalpy=True,
        hydrogen_enthalpy_species="ParaHydrogen",
    )
    assert para.hydrogen_enthalpy_species == "ParaHydrogen"
    assert para._mixture_enthalpy(80.0, 1.0) != pytest.approx(
        normal._mixture_enthalpy(80.0, 1.0), rel=1.0e-3
    )
    result = para.solve(
        maximum_distance=0.04,
        maximum_step=0.0005,
        relative_tolerance=2.0e-7,
        method="LSODA",
    )
    assert max(abs(v) for v in para.establishment_residuals.values()) < 1e-8
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )


def test_hydrogen_spin_isomer_requires_supported_phase_enthalpy():
    base = model()
    common = dict(
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
    )
    with pytest.raises(ValueError, match="spin-isomer"):
        ConservedGaussianJet(
            base.source, hydrogen_enthalpy_species="ParaHydrogen", **common
        )
    with pytest.raises(ValueError, match="hydrogen enthalpy species"):
        ConservedGaussianJet(
            base.source,
            temperature_dependent_phase_enthalpy=True,
            hydrogen_enthalpy_species="EquilibriumHydrogen",
            **common,
        )


def test_two_scalar_phase_profile_is_rejected_before_integration():
    base = model()
    with pytest.raises(ValueError, match="requires four-flux"):
        ConservedGaussianJet(
            base.source,
            ambient_temperature=base.ambient_temperature,
            ambient_pressure=base.ambient_pressure,
            ambient_density=base.ambient_density,
            fuel_molecular_weight=base.fuel_molecular_weight,
            ambient_molecular_weight=base.ambient_molecular_weight,
            fuel_heat_capacity=base.fuel_heat_capacity,
            ambient_heat_capacity=base.ambient_heat_capacity,
            radial_points=41,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            thermodynamic_profile="temperature_mass_fraction",
            thermodynamic_spreading_ratio=1.16 * math.sqrt(0.70 / 0.85),
        )


@pytest.mark.slow
def test_four_flux_two_scalar_phase_profile_is_conservative():
    base = model()
    phase = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment=True,
        equilibrium_air_condensation=True,
        temperature_dependent_phase_enthalpy=True,
        thermodynamic_profile="temperature_mass_fraction",
        thermodynamic_spreading_ratio=1.16 * math.sqrt(0.70 / 0.85),
    )
    result = phase.solve(
        maximum_distance=0.04,
        maximum_step=0.0005,
        relative_tolerance=2.0e-7,
        method="LSODA",
    )
    assert max(abs(v) for v in phase.establishment_residuals.values()) < 1e-8
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )


@pytest.mark.slow
def test_argon_phase_is_explicit_physical_and_conservative():
    base = model()
    argon = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=41,
        conservative_establishment="scalar_peak",
        equilibrium_air_condensation=True,
        temperature_dependent_phase_enthalpy=True,
        equilibrium_argon_condensation=True,
    )
    result = argon.solve(
        maximum_distance=0.04,
        maximum_step=0.0005,
        relative_tolerance=2.0e-7,
        method="LSODA",
    )
    assert argon._dry_argon_mass_fraction == pytest.approx(0.0129, rel=0.02)
    assert max(abs(v) for v in argon.establishment_residuals.values()) < 1e-8
    assert np.all((result.temperature >= 14.1) & (result.temperature <= 400.0))
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )


@pytest.mark.slow
def test_public_lh2_research_path_uses_the_accepted_configuration():
    expanded = lh2_source_from_measured_throat(
        throat_diameter=0.001,
        throat_pressure=2.42e5,
        throat_temperature=37.4,
        throat_density=1.65,
        throat_velocity=498.2,
    )
    assert expanded.expansion.relative_mass_residual < 1.0e-12
    assert expanded.expansion.relative_momentum_residual < 1.0e-12
    assert expanded.expansion.relative_energy_residual < 1.0e-12
    result = run_lh2_near_field_research(
        expanded,
        maximum_distance=0.025,
        radial_points=41,
        maximum_step=0.001,
        relative_tolerance=2.0e-6,
    )
    assert result.model.establishment == "scalar_peak"
    assert result.model.equilibrium_air_condensation
    assert result.model.temperature_dependent_phase_enthalpy
    assert result.model.hydrogen_enthalpy_species == "Hydrogen"
    assert result.model.ambient_absolute_humidity == 0.0
    assert result.conservative
    assert not result.warnings
    assert "conservation screen         : pass" in result.report()


def test_public_lh2_research_path_rejects_unknown_spin_isomer():
    with pytest.raises(ValueError, match="spin isomer"):
        run_lh2_near_field_research(object(), hydrogen_spin_isomer="mixed")


def test_handoff_search_rejects_nonpositive_step_before_running_source():
    source = AxisymmetricJetSource(
        diameter=0.001, velocity=500.0, density=0.5,
        temperature=45.0, theta=0.0, y=0.5,
    )
    with pytest.raises(ValueError, match="search step"):
        run_lh2_crosswind_research(
            source, wind=2.5, maximum_nearfield_distance=0.02,
            handoff_search_start=0.01, handoff_search_step=0.0,
            radial_points=41, maximum_step=0.001,
            relative_tolerance=2.0e-6,
        )


def test_public_lh2_research_path_rejects_two_phase_hydrogen_source():
    source = AxisymmetricJetSource(
        diameter=0.0276,
        velocity=122.0,
        density=11.3,
        temperature=20.37,
        theta=0.0,
        y=0.5,
    )
    with pytest.raises(ValueError, match="single-phase gas source"):
        run_lh2_near_field_research(source, maximum_distance=0.01)


def test_axisymmetric_result_rejects_streamline_extrapolation():
    jet = model()
    result = jet.solve(
        maximum_distance=0.01,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    state = result.state_at_s(float(result.S[-1]))
    assert state[0] == pytest.approx(result.velocity[-1])
    with pytest.raises(ValueError, match="outside the computed range"):
        result.state_at_s(float(result.S[-1]) + 0.01)


@pytest.mark.slow
def test_phase_profile_handoff_passes_without_relaxing_failed_candidate():
    expanded = lh2_source_from_measured_throat(
        throat_diameter=0.001,
        throat_pressure=2.422e5,
        throat_temperature=37.4,
        throat_density=1.65,
        throat_velocity=498.2,
        theta=0.0,
        y=0.5,
    )
    coupled = run_lh2_crosswind_research(
        expanded,
        wind=2.5,
        height=0.5,
        relative_humidity=0.0,
        ambient_temperature=295.0,
        ambient_pressure=101325.0,
        maximum_nearfield_distance=0.08,
        handoff_distance=0.08,
        minimum_mass_fraction=1.0e-5,
        radial_points=41,
        maximum_step=0.001,
        relative_tolerance=2.0e-6,
    )
    near = coupled.near_field
    crosswind = coupled.handoff.model

    rejected = handoff_lh2_near_field_to_crosswind(
        near,
        crosswind,
        streamline_distance=0.025,
        thermodynamic_closure="centre_state",
    )
    assert not rejected.accepted
    with pytest.raises(RuntimeError, match="handoff screen failed"):
        rejected.run(distmx=1.0)

    handoff = coupled.handoff
    assert coupled.accepted
    assert near.model.establishment == "entrained_mass"
    assert handoff.accepted
    assert max(
        handoff.relative_residuals[name]
        for name in ("total_mass", "hydrogen", "momentum_x", "momentum_z")
    ) < 1.0e-8
    assert handoff.relative_residuals["energy"] < 0.02
    assert handoff.halfwidth_residual < 0.05
    assert handoff.temperature_residual < 2.0
    assert handoff.energy_quadrature_residual < 1.0e-5
    assert handoff.energy_quadrature_points in (64, 128, 256, 512)
    assert handoff.crosswind_entrainment == "source_momentum"
    assert crosswind.k.sc == pytest.approx(near.model.spreading_ratio**2)
    energy_interface = audit_lh2_independent_energy_interface(
        near, crosswind, streamline_distance=0.08,
    )
    assert energy_interface.accepted
    assert max(energy_interface.relative_residuals.values()) < 1.0e-8
    assert energy_interface.energy_quadrature_residual < 1.0e-5
    assert energy_interface.halfwidth_residual < 0.05
    assert energy_interface.temperature_residual < 2.0
    assert energy_interface.energy_quadrature_points in (64, 128, 256, 512)
    enthalpy_interface = audit_lh2_independent_energy_interface(
        near,
        crosswind,
        streamline_distance=0.08,
        energy_transport="enthalpy",
    )
    assert enthalpy_interface.accepted
    assert enthalpy_interface.energy_transport == "enthalpy"
    assert max(enthalpy_interface.relative_residuals.values()) < 1.0e-8
    assert enthalpy_interface.energy_quadrature_residual < 1.0e-5
    assert enthalpy_interface.halfwidth_residual < 0.05
    assert enthalpy_interface.temperature_residual < 2.0
    assert math.isfinite(enthalpy_interface.target_buoyancy_force)
    assert math.isfinite(enthalpy_interface.projected_buoyancy_force)
    assert math.isfinite(enthalpy_interface.buoyancy_force_ratio)
    scalar_width = math.sqrt(2.0 * enthalpy_interface.state[2])
    assert enthalpy_interface.houf_width_mapping == "velocity"
    assert enthalpy_interface.model.houf_velocity_width(
        enthalpy_interface.state
    ) == pytest.approx(
        scalar_width / near.model.spreading_ratio
    )
    scalar_interface = audit_lh2_independent_energy_interface(
        near,
        crosswind,
        streamline_distance=0.08,
        energy_transport="enthalpy",
        houf_width_mapping="scalar",
    )
    assert scalar_interface.model.houf_velocity_width(
        scalar_interface.state
    ) == pytest.approx(math.sqrt(2.0 * scalar_interface.state[2]))
    derivatives = energy_interface.model.derivatives(
        0.0, energy_interface.state
    )
    assert np.all(np.isfinite(derivatives))
    assert derivatives[1] < 0.0
    assert derivatives[2] > 0.0
    coarse = energy_interface.model.solve(
        energy_interface.state,
        maximum_distance=0.01,
        maximum_step=0.002,
        relative_tolerance=1.0e-5,
    )
    fine = energy_interface.model.solve(
        energy_interface.state,
        maximum_distance=0.01,
        maximum_step=0.001,
        relative_tolerance=1.0e-5,
    )
    assert coarse.maximum_relative_balance_residual < 1.0e-4
    assert fine.maximum_relative_balance_residual < 1.0e-4
    assert coarse.states[-1] == pytest.approx(fine.states[-1], rel=2.0e-4)
    assert fine.states[-1, 1] < fine.states[0, 1]
    assert fine.states[-1, 2] > fine.states[0, 2]
    assert fine.temperatures[-1] > fine.temperatures[0]
    downstream = handoff.run(distmx=0.5, smax=1.0, tol=2.0e-4)
    assert len(downstream.rows) > 2


def test_ground_image_geometry_is_refitted_at_sensor_heights():
    """The geometry audit compares fitted observables, not hidden state."""
    from degali.validation.nearfield import _fit_model_vertical_profile

    class ImageProfile:
        @staticmethod
        def concentration_at(_x, _y, height):
            centre = 0.5
            sigma = 0.6
            return 50.0 * (
                math.exp(-0.5 * ((height - centre) / sigma) ** 2)
                + math.exp(-0.5 * ((height + centre) / sigma) ** 2)
            )

    fitted = _fit_model_vertical_profile(
        ImageProfile(),
        6.0,
        np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
    )
    assert fitted["r_squared"] > 0.99
    assert fitted["centre"] < 0.2
    assert fitted["sigma_z"] > 0.7


@pytest.mark.slow
def test_five_balance_jet_conserves_species_and_total_energy():
    jet = model()
    result = jet.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    assert len(result.S) > 10
    assert np.all(np.diff(result.width) > 0.0)
    assert np.all(np.diff(result.mass_fraction) < 0.0)
    assert np.all(np.diff(result.temperature) > 0.0)
    assert result.mass_flux[-1] > result.mass_flux[0]
    assert result.species_flux[-1] == pytest.approx(
        result.species_flux[0], rel=2.0e-4
    )
    assert result.energy_flux[-1] == pytest.approx(
        result.energy_flux[0], rel=2.0e-4
    )

    state = result.state_at_y(0.04)
    mass_width = jet.half_width(state, "mass")
    temperature_width = jet.half_width(state, "temperature")
    assert 0.0 < mass_width < temperature_width < 20.0 * state[1]
    assert math.isfinite(jet.radial_state(state, mass_width)[2])


def test_bad_axisymmetric_source_is_rejected():
    source = AxisymmetricJetSource(
        diameter=0.001, velocity=500.0, density=0.5,
        temperature=45.0, mass_fraction=0.0,
    )
    with pytest.raises(ValueError, match="mass fraction"):
        ConservedGaussianJet(
            source,
            ambient_temperature=295.0,
            ambient_pressure=101325.0,
            ambient_density=1.197,
            fuel_molecular_weight=0.00201588,
            ambient_molecular_weight=0.02896546,
            fuel_heat_capacity=14294.8,
            ambient_heat_capacity=1006.2,
        )


def test_radiative_heat_has_a_separate_external_energy_ledger():
    base = model()
    jet = ConservedGaussianJet(
        base.source,
        ambient_temperature=base.ambient_temperature,
        ambient_pressure=base.ambient_pressure,
        ambient_density=base.ambient_density,
        fuel_molecular_weight=base.fuel_molecular_weight,
        ambient_molecular_weight=base.ambient_molecular_weight,
        fuel_heat_capacity=base.fuel_heat_capacity,
        ambient_heat_capacity=base.ambient_heat_capacity,
        radial_points=81,
        conservative_establishment="scalar_peak",
        radiative_absorptivity=1.0,
    )
    result = jet.solve(
        maximum_distance=0.04,
        maximum_step=0.002,
        relative_tolerance=2.0e-5,
    )
    assert result.radiative_heat_added[0] == 0.0
    assert result.radiative_heat_added[-1] > 0.0
    corrected = result.energy_flux - result.radiative_heat_added
    assert corrected[-1] == pytest.approx(corrected[0], rel=2.0e-4)


def test_radiative_absorptivity_is_bounded():
    base = model()
    with pytest.raises(ValueError, match="absorptivity"):
        ConservedGaussianJet(
            base.source,
            ambient_temperature=base.ambient_temperature,
            ambient_pressure=base.ambient_pressure,
            ambient_density=base.ambient_density,
            fuel_molecular_weight=base.fuel_molecular_weight,
            ambient_molecular_weight=base.ambient_molecular_weight,
            fuel_heat_capacity=base.fuel_heat_capacity,
            ambient_heat_capacity=base.ambient_heat_capacity,
            radiative_absorptivity=1.01,
        )
