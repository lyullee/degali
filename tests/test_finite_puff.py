import math

import numpy as np
import pytest

from degali.addons.finite_puff import (
    GaussianPuffConfig,
    integrate_finite_gaussian_puff,
)
from degali.addons.finite_release import FiniteReleasePuffHandoff
from degali.addons.transient_receptor import WindHistory


class _IdealThermodynamics:
    ambient_temperature = 300.0
    ambient_density = 1.2
    fuel_molecular_weight = 0.002016
    _humid_ambient_molecular_weight = 0.029
    _ambient_enthalpy = 3.0e5

    @staticmethod
    def _temperature_from_enthalpy(enthalpy, _fraction):
        return enthalpy / 1000.0

    @staticmethod
    def _density_from_temperature(temperature, fraction):
        inverse_mw = fraction / 0.002016 + (1.0 - fraction) / 0.029
        mixture_mw = 1.0 / inverse_mw
        return 1.2 * (300.0 / temperature) * (mixture_mw / 0.029)


def _handoff(*, temperature=280.0, height=4.0):
    mass = 10.0
    momentum = (20.0, 0.0, 0.0)
    energy = mass * 1000.0 * (temperature - 300.0) + 20.0
    return FiniteReleasePuffHandoff(
        status="transition_ready",
        source_duration_s=5.0,
        transition_time_s=5.0,
        transition_arc_length_m=10.0,
        centre_position_m=(10.0, 0.0, height),
        transverse_widths_m=(2.0, 1.0),
        longitudinal_half_length_m=5.0,
        bulk_velocity_m_s=(2.0, 0.0, 0.0),
        total_mass_kg=mass,
        hydrogen_mass_kg=0.01,
        momentum_kg_m_s=momentum,
        relative_energy_j=energy,
        hydrogen_flux_residual_kg_s=0.0,
        hydrogen_inventory_residual_kg=0.0,
        requires_transient_puff_continuation=True,
        qualification="test handoff",
    )


def test_native_puff_conserves_hydrogen_and_closes_entrained_mass():
    result = integrate_finite_gaussian_puff(
        _handoff(), _IdealThermodynamics(),
        GaussianPuffConfig(
            duration_s=1.0, time_step_s=0.01,
            wind_velocity_m_s=(4.0, 0.0, 0.0),
            shear_entrainment_coefficient=0.01,
            ambient_turbulence_coefficient=0.0,
            form_drag_coefficient=0.0,
        ),
    )
    assert result.maximum_relative_mass_residual < 1.0e-12
    assert result.maximum_relative_hydrogen_residual < 1.0e-12
    assert result.states[-1].total_mass_kg > result.states[0].total_mass_kg
    assert result.states[-1].bulk_hydrogen_mass_fraction \
        < result.states[0].bulk_hydrogen_mass_fraction
    assert result.states[-1].bulk_velocity_m_s[0] > 2.0


def test_dense_ground_puff_spreads_in_both_horizontal_directions():
    result = integrate_finite_gaussian_puff(
        _handoff(height=1.0), _IdealThermodynamics(),
        GaussianPuffConfig(
            duration_s=0.5, time_step_s=0.01,
            wind_velocity_m_s=(2.0, 0.0, 0.0),
            shear_entrainment_coefficient=0.0,
            ambient_turbulence_coefficient=0.0,
            form_drag_coefficient=0.0,
        ),
    )
    first, last = result.states[0], result.states[-1]
    assert first.ground_contact
    assert last.longitudinal_half_length_m > first.longitudinal_half_length_m
    assert last.lateral_half_width_m > first.lateral_half_width_m
    assert last.centre_position_m[2] >= last.vertical_half_width_m


def test_receptor_trace_is_queryable_and_rejects_extrapolation():
    result = integrate_finite_gaussian_puff(
        _handoff(), _IdealThermodynamics(),
        GaussianPuffConfig(
            duration_s=1.0, time_step_s=0.02,
            wind_velocity_m_s=(2.0, 0.0, 0.0),
            shear_entrainment_coefficient=0.0,
            ambient_turbulence_coefficient=0.0,
            gravity_spreading_coefficient=0.0,
            form_drag_coefficient=0.0,
        ),
    )
    trace = result.receptor_trace((11.0, 0.0, 4.0))
    assert np.all(np.isfinite(trace.mole_fraction))
    assert 0.0 < trace.peak_mole_fraction < 1.0
    assert 0.0 <= trace.peak_time_s <= 1.0
    with pytest.raises(ValueError, match="outside"):
        result.state_at(1.1)


def test_puff_requires_transition_ready_handoff():
    invalid = _handoff().__class__(
        **{**_handoff().__dict__, "status": "transition_beyond_trajectory"}
    )
    with pytest.raises(ValueError, match="transition-ready"):
        integrate_finite_gaussian_puff(
            invalid, _IdealThermodynamics(), GaussianPuffConfig(duration_s=1.0)
        )


def test_vector_wind_history_splits_boundaries_and_turns_puff_momentum():
    history = WindHistory(
        time_s=[0.0, 0.5, 1.0],
        speed_m_s=[4.0, 4.0, 4.0],
        direction_from_deg=[270.0, 180.0, 180.0],
    )
    result = integrate_finite_gaussian_puff(
        _handoff(), _IdealThermodynamics(),
        GaussianPuffConfig(
            duration_s=1.0, time_step_s=0.3,
            wind_history=history,
            shear_entrainment_coefficient=0.02,
            ambient_turbulence_coefficient=0.0,
            form_drag_coefficient=0.0,
        ),
    )
    times = np.array([state.elapsed_s for state in result.states])
    assert np.any(np.isclose(times, 0.5, atol=1.0e-12))
    assert result.states[-1].bulk_velocity_m_s[1] > 0.0
    assert result.states[-1].centre_position_m[1] > 0.0


def test_wind_history_must_cover_puff_window():
    history = WindHistory(
        time_s=[0.0, 0.5], speed_m_s=[2.0, 2.0],
        direction_from_deg=[270.0, 270.0],
    )
    with pytest.raises(ValueError, match="complete puff"):
        integrate_finite_gaussian_puff(
            _handoff(), _IdealThermodynamics(),
            GaussianPuffConfig(duration_s=1.0, wind_history=history),
        )
