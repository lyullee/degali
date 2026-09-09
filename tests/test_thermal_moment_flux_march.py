"""Direct-flux inverse and RK4 conservation tests."""

from types import SimpleNamespace

import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_transverse_mixing import mixing
from degali.addons.reservoir_thermal import ReservoirShortSegment, encode_section
from degali.addons.thermal_moment_flux_march import (
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
)


def test_phase_partitioned_six_flux_inverse_recovers_a_perturbed_section(mixing):
    initial = encode_section(mixing.state, mixing.section.thermal_width_ratio)
    inverter = ThermalMomentFluxInverter(
        mixing.section.jetplume, mixing.section.thermodynamics,
        order=8, probes=257, quadrature_points=128,
    )
    target = ReservoirShortSegment.flux_values(mixing, order=8)
    guess = initial.copy()
    guess[[0, 2, 4, 7]] += [.01, -.01, .005, -.008]
    result = inverter.match(target, guess)
    assert np.max(np.abs(result.scaled_residuals)) < 1e-8
    assert result.fluxes == pytest.approx(target, rel=1e-8, abs=1e-9)
    assert result.parameters == pytest.approx(initial, rel=2e-6, abs=2e-6)


@pytest.mark.parametrize("target", [np.ones(5), np.full(6, np.nan), np.array([1., 0., 1., 0., 1., 1.])])
def test_flux_inverse_rejects_invalid_targets(mixing, target):
    inverter = ThermalMomentFluxInverter(mixing.section.jetplume, mixing.section.thermodynamics)
    with pytest.raises(ValueError, match="six finite fluxes"):
        inverter.match(target, encode_section(mixing.state, mixing.section.thermal_width_ratio))


class _IdentityInverter:
    def fluxes(self, parameters):
        return np.array([2., .2, 3., 0., 5., -4.])

    def match(self, target, initial_parameters, *, position=None):
        parameters = np.asarray(initial_parameters, float).copy()
        parameters[5:7] = position
        state = np.array([1., .1, 1., 0., 1., position[0], position[1]])
        return SimpleNamespace(parameters=parameters, state=state, fluxes=np.asarray(target),
            scaled_residuals=np.zeros(6), mixing=None)


def test_rk4_advances_six_primary_fluxes_and_hits_exact_x_target():
    source = np.array([.4, 0., .2, -.1, 3., -2.])
    def local(_match):
        return dict(ledger={"sources": source[:5]}, moment_rate=source[5],
            weak_budget_scaled_error=2e-10, minimum_chi_species=.2,
            minimum_chi_momentum=.3, edge_gradient_defects={"heat": .25}, valid=True)
    initial = np.zeros(8)
    march = FluxSpaceThermalMomentMarch(
        _IdentityInverter(), thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", local_evaluator=local,
    ).march(initial, .2, step=.03, target_x=.095)
    assert march.reached_target
    assert march.arc_length[-1] == pytest.approx(.095, abs=1e-12)
    assert march.fluxes[-1] == pytest.approx(march.fluxes[0]+.095*source, abs=1e-14)
    assert march.cumulative_sources[-1] == pytest.approx(.095*source, abs=1e-14)
    assert march.fluxes[-1]-march.fluxes[0] == pytest.approx(march.cumulative_sources[-1], abs=1e-14)
    assert march.maximum_weak_residual == pytest.approx(2e-10)
    assert march.minimum_sampled_diffusivity == pytest.approx(.2)
    assert march.maximum_edge_heat_defect == pytest.approx(.25)


def test_flux_march_rejects_invalid_physics_choices():
    with pytest.raises(ValueError):
        FluxSpaceThermalMomentMarch(_IdentityInverter(), thermal_species_ratio=0.,
                                    mechanical_work="reduced_buoyancy_work")
    with pytest.raises(ValueError):
        FluxSpaceThermalMomentMarch(_IdentityInverter(), thermal_species_ratio=1.,
                                    mechanical_work="none")
