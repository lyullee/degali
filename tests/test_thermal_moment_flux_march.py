"""Direct-flux inverse and RK4 conservation tests."""

from types import SimpleNamespace

import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_transverse_mixing import mixing
from degali.addons.reservoir_thermal import ReservoirShortSegment, encode_section
from degali.addons.thermal_moment_flux_march import (
    FluxInverseError,
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
    ThermalMomentFluxTrajectory,
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
    def __init__(self):
        self.tolerances = []

    def fluxes(self, parameters):
        return np.array([2., .2, 3., 0., 5., -4.])

    def match(self, target, initial_parameters, *, position=None, tolerance=None):
        self.tolerances.append(tolerance)
        parameters = np.asarray(initial_parameters, float).copy()
        parameters[5:7] = position
        state = np.array([1., .1, 1., 0., 1., position[0], position[1]])
        return SimpleNamespace(parameters=parameters, state=state, fluxes=np.asarray(target),
            scaled_residuals=np.zeros(6), mixing=None)


class _EndpointStepInverter(_IdentityInverter):
    def __init__(self):
        super().__init__()
        self.accepted_x = 0.

    def match(self, target, initial_parameters, *, position=None, tolerance=None):
        if tolerance is None and position[0]-self.accepted_x > .02:
            raise FluxInverseError("manufactured endpoint left the manifold")
        result = super().match(
            target, initial_parameters, position=position, tolerance=tolerance,
        )
        if tolerance is None:
            self.accepted_x = position[0]
        return result


class _EndpointWarmStartInverter(_IdentityInverter):
    def match(self, target, initial_parameters, *, position=None, tolerance=None):
        initial = np.asarray(initial_parameters, float)
        if tolerance is None and initial[5] == pytest.approx(position[0]):
            raise FluxInverseError("last-stage warm start selected the wrong branch")
        return super().match(
            target, initial_parameters, position=position, tolerance=tolerance,
        )


def test_rk4_advances_six_primary_fluxes_and_hits_exact_x_target():
    source = np.array([.4, 0., .2, -.1, 3., -2.])
    def local(_match):
        return dict(ledger={"sources": source[:5]}, moment_rate=source[5],
            weak_budget_scaled_error=2e-10, minimum_chi_species=.2,
            minimum_chi_momentum=.3, edge_gradient_defects={"heat": .25}, valid=True)
    initial = np.zeros(8)
    inverter = _IdentityInverter()
    march = FluxSpaceThermalMomentMarch(
        inverter, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", stage_inverse_tolerance=2e-5,
        local_evaluator=local,
    ).march(initial, .2, step=.03, target_x=.095)
    assert march.reached_target
    assert march.arc_length[-1] == pytest.approx(.095, abs=1e-12)
    assert march.fluxes[-1] == pytest.approx(march.fluxes[0]+.095*source, abs=1e-14)
    assert march.cumulative_sources[-1] == pytest.approx(.095*source, abs=1e-14)
    assert march.fluxes[-1]-march.fluxes[0] == pytest.approx(march.cumulative_sources[-1], abs=1e-14)
    assert march.maximum_weak_residual == pytest.approx(2e-10)
    assert march.minimum_sampled_diffusivity == pytest.approx(.2)
    assert march.maximum_edge_heat_defect == pytest.approx(.25)
    assert any(value == pytest.approx(2e-5) for value in inverter.tolerances)
    assert inverter.tolerances[-1] is None


def test_flux_march_rejects_invalid_physics_choices():
    with pytest.raises(ValueError):
        FluxSpaceThermalMomentMarch(_IdentityInverter(), thermal_species_ratio=0.,
                                    mechanical_work="reduced_buoyancy_work")
    with pytest.raises(ValueError):
        FluxSpaceThermalMomentMarch(_IdentityInverter(), thermal_species_ratio=1.,
                                    mechanical_work="none")
    with pytest.raises(ValueError, match="positivity domain"):
        FluxSpaceThermalMomentMarch(
            _IdentityInverter(), thermal_species_ratio=1.,
            mechanical_work="reduced_buoyancy_work", positivity_domain="corners",
        )


def test_only_inverse_failure_halves_and_retains_the_smaller_step():
    source = np.array([.4, 0., .2, -.1, 3., -2.])
    def local(_match):
        return dict(ledger={"sources": source[:5]}, moment_rate=source[5],
            weak_budget_scaled_error=0., minimum_chi_species=.2,
            minimum_chi_momentum=.3, edge_gradient_defects={"heat": .25}, valid=True)
    result = FluxSpaceThermalMomentMarch(
        _EndpointStepInverter(), thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", stage_inverse_tolerance=2e-5,
        local_evaluator=local,
    ).march(np.zeros(8), .05, step=.03, minimum_step=.005)
    assert result.stop_reason == "arc-length ceiling reached"
    assert result.rejected_steps == 1
    assert result.accepted_steps == 4
    assert result.minimum_accepted_step == pytest.approx(.005)
    assert result.fluxes[-1]-result.fluxes[0] == pytest.approx(.05*source, abs=1e-14)


def test_endpoint_retries_once_from_the_previous_accepted_section():
    source = np.array([.4, 0., .2, -.1, 3., -2.])
    def local(_match):
        return dict(ledger={"sources": source[:5]}, moment_rate=source[5],
            weak_budget_scaled_error=0., minimum_chi_species=.2,
            minimum_chi_momentum=.3, edge_gradient_defects={"heat": .25}, valid=True)
    result = FluxSpaceThermalMomentMarch(
        _EndpointWarmStartInverter(), thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", stage_inverse_tolerance=2e-5,
        local_evaluator=local,
    ).march(np.zeros(8), .01, step=.01, minimum_step=.001)
    assert result.accepted_steps == 1
    assert result.rejected_steps == 0
    assert result.endpoint_fallback_attempts == 1
    assert result.endpoint_fallback_successes == 1


def test_reduced_step_grows_only_after_eight_successes():
    source = np.array([.4, 0., .2, -.1, 3., -2.])
    def local(_match):
        return dict(ledger={"sources": source[:5]}, moment_rate=source[5],
            weak_budget_scaled_error=0., minimum_chi_species=.2,
            minimum_chi_momentum=.3, edge_gradient_defects={"heat": .25}, valid=True)
    result = FluxSpaceThermalMomentMarch(
        _EndpointStepInverter(), thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", stage_inverse_tolerance=2e-5,
        local_evaluator=local,
    ).march(np.zeros(8), .15, step=.03, minimum_step=.005)
    assert result.stop_reason == "arc-length ceiling reached"
    assert result.rejected_steps == 2
    assert result.step_growths == 1
    assert result.accepted_steps == 10
    assert result.minimum_accepted_step == pytest.approx(.015)


def test_local_physical_failure_is_not_retried_or_step_halved():
    def invalid(_match):
        return {"valid": False}
    result = FluxSpaceThermalMomentMarch(
        _IdentityInverter(), thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", stage_inverse_tolerance=2e-5,
        local_evaluator=invalid,
    ).march(np.zeros(8), .01, step=.01, minimum_step=.001)
    assert result.accepted_steps == 0
    assert result.rejected_steps == 0
    assert result.endpoint_fallback_attempts == 0
    assert "local thermal-moment closure" in result.stop_reason


def test_local_failure_reports_the_specific_physical_gate():
    def invalid(_match):
        return dict(
            ledger={"sources": np.zeros(5)}, moment_rate=0.,
            weak_budget_scaled_error=0., minimum_chi_species=-.2,
            minimum_chi_momentum=.3, edge_gradient_defects={"heat": 0.},
            family=SimpleNamespace(maximum_scaled_residual=0.),
            incoming=True, positive_diffusion=False,
            maximum_outward_mass=-.1, curvature_half_width=.01, valid=False,
        )
    result = FluxSpaceThermalMomentMarch(
        _IdentityInverter(), thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work", local_evaluator=invalid,
    ).march(np.zeros(8), .01, step=.01)
    assert "chi=(-2.000e-01,3.000e-01)" in result.stop_reason


def test_flux_trajectory_interpolates_parameters_and_queries_receptors(section):
    first = encode_section(
        np.array([1.13, .15, .02, .01, 25., 1., 1.5]), 1.02,
    )
    second = first.copy()
    second[[0, 1, 2, 4, 7]] += np.log([1.01, .99, 1.02, .98, 1.03])
    second[5:7] = [2., 1.6]
    trajectory = ThermalMomentFluxTrajectory(
        section.jetplume, section.thermodynamics, [second, first],
    )
    middle = trajectory.parameters_at(1.5)
    assert middle == pytest.approx(.5*(first+second))
    state, model = trajectory._section_at(1.5)
    assert trajectory.state_at(1.5) == pytest.approx(state)
    assert trajectory.temperature_at(1.5, 0., 1.5) == pytest.approx(
        model.point_temperature(state, 0., 1.5)
    )
    assert trajectory.concentration_at(1.5, 0., 1.5) == pytest.approx(
        100.*model.point_mole_fraction(state, 0., 1.5)
    )
    assert trajectory.parameters_at(2.+5e-9) == pytest.approx(second)
    assert trajectory.parameters_at(2.+2e-8) is None
    assert trajectory.state_at(.9) is None


def test_flux_trajectory_rejects_nonunique_downwind_coordinates(section):
    parameters = np.zeros((2, 8))
    parameters[:, 5] = 1.
    with pytest.raises(ValueError, match="increase strictly"):
        ThermalMomentFluxTrajectory(
            section.jetplume, section.thermodynamics, parameters,
        )
