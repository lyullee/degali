"""Independent EOS identities and pipe-to-Gaussian energy-ledger tests."""

import math
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("CoolProp")
from CoolProp.CoolProp import PropsSI

from degali.addons.axisymmetric_jet import ConservedGaussianJet
from degali.addons.cryogenic_air import multiphase_hydrogen_evaporation_endpoint
from degali.addons.hydrogen_eos import hydrogen_gas_departure
from degali.addons.lh2_droplets import homogeneous_equilibrium_hydrogen_source


@pytest.mark.parametrize("species", ["Hydrogen", "ParaHydrogen", "OrthoHydrogen"])
def test_departures_share_one_stable_eos(species):
    t = PropsSI("T", "P", 101325.0, "Q", 1, species)
    state = hydrogen_gas_departure(t, 101325.0, species)
    assert state.density == pytest.approx(PropsSI("D", "P", 101325., "Q", 1, species), rel=1e-8)
    assert state.compressibility == pytest.approx(state.ideal_density / state.density)
    assert state.enthalpy_departure == pytest.approx(state.enthalpy - state.ideal_enthalpy)
    assert -0.11 < state.relative_volume_departure < -0.08
    assert abs(hydrogen_gas_departure(295., 101325., species).relative_volume_departure) < 0.001


@pytest.mark.parametrize("t,p", [(18., 101325.), (13., 1000.), (math.nan, 1.), (20., -1.)])
def test_departure_rejects_invalid_or_metastable_gas(t, p):
    with pytest.raises(ValueError):
        hydrogen_gas_departure(t, p)


def _bound(species="Hydrogen"):
    return homogeneous_equilibrium_hydrogen_source(
        mass_flow=0.26527468341375925, orifice_diameter=0.012,
        upstream_temperature=20.186769311178438,
        upstream_pressure=277411.2291691, ambient_temperature=288.15,
        hydrogen_species=species,
    )


def _model(source, species="Hydrogen", **overrides):
    common = dict(
        ambient_temperature=288.15, ambient_pressure=101325.,
        ambient_density=PropsSI("D", "T", 288.15, "P", 101325., "Air"),
        fuel_molecular_weight=PropsSI("M", "Hydrogen"),
        ambient_molecular_weight=PropsSI("M", "Air"),
        fuel_heat_capacity=PropsSI("C", "T", 288.15, "P", 101325., "Hydrogen"),
        ambient_heat_capacity=PropsSI("C", "T", 288.15, "P", 101325., "Air"),
        equilibrium_air_condensation=True, temperature_dependent_phase_enthalpy=True,
        hydrogen_enthalpy_species=species, radial_points=41,
        conservative_establishment="entrained_mass",
    )
    common.update(overrides)
    return ConservedGaussianJet(source, **common)


# CoolProp provides the ortho EOS, but no ortho viscosity for the existing
# droplet-source correlation; do not substitute a different isomer's model.
@pytest.mark.parametrize("species", ["Hydrogen", "ParaHydrogen"])
@pytest.mark.parametrize("consistent", [False, True])
def test_pipe_phase_gaussian_energy_is_conserved_against_independent_upstream_state(species, consistent):
    bound = _bound(species)
    pipe, endpoint = bound.postflash, bound.phase_plane.endpoint
    model = _model(bound.source, species, consistent_phase_ambient=consistent)
    grid, tables = model._phase_enthalpy_tables
    h_upstream = PropsSI("H", "T", pipe.upstream_temperature,
                        "P", pipe.upstream_pressure, species)
    h_reference = np.interp(288.15, grid, tables["Hydrogen"])
    expected = pipe.mass_flow * (h_upstream - h_reference + 0.5 * pipe.upstream_velocity**2)
    assert endpoint.incoming_specific_kinetic_energy == pytest.approx(0.5 * pipe.upstream_velocity**2)
    assert endpoint.outgoing_specific_energy == pytest.approx(endpoint.incoming_specific_energy, rel=1e-12)
    assert model._established_target_fluxes(0.0)[3] == pytest.approx(expected, rel=1e-12)
    _, established = model.established_initial_state()
    assert model._fluxes(established)[4] == pytest.approx(expected, rel=1e-8)
    # Deliberately removing the phase ledger must expose the former defect.
    old_model = _model(replace(bound.source, enthalpy_boundary=None), species)
    old_energy = old_model._established_target_fluxes(0.0)[3]
    assert abs(old_energy - expected) / abs(expected) > 0.04


@pytest.mark.parametrize("change,match", [
    ({"ambient_temperature": 289.15}, "ambient reference"),
    ({"hydrogen_enthalpy_species": "ParaHydrogen"}, "hydrogen species"),
    ({"temperature_dependent_phase_enthalpy": False}, "component phase"),
    ({"conservative_establishment": "published"}, "conservative establishment"),
])
def test_phase_energy_boundary_rejects_incompatible_receiver(change, match):
    with pytest.raises(ValueError, match=match):
        _model(_bound().source, **change)


@pytest.mark.parametrize("ke", [-1., math.nan, math.inf])
def test_endpoint_rejects_invalid_incoming_kinetic_energy(ke):
    with pytest.raises(ValueError, match="kinetic energy"):
        multiphase_hydrogen_evaporation_endpoint(
            storage_temperature=20., incoming_specific_kinetic_energy=ke,
        )
