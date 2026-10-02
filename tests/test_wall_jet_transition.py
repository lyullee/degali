import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.axisymmetric_jet import AxisymmetricJetSource, ConservedGaussianJet
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.wall_jet_transition import (
    FiniteWallJetCrosswind,
    WallJetTransitionConfig,
)
from degali.core.jetplume import JetCoefficients


def base_model(
    *, ground="geometry", wind=2.0, ground_heat_coefficient=0.0,
    ground_surface_temperature=None,
):
    jetplume = SimpleNamespace(
        th=None,
        k=JetCoefficients(sc=1.16**2),
        deltay=0.1,
        deltaz=0.08,
        betay=0.92,
        betaz=0.87,
        gammaz=-0.012,
        ustar=0.2,
        _split=lambda area, *_: (math.sqrt(area), math.sqrt(area)),
        _wind=lambda *_: wind,
        _wind_profile=lambda *_: (wind, 0.0),
    )
    source = AxisymmetricJetSource(
        diameter=0.01, velocity=100.0, density=0.6, temperature=50.0
    )
    thermodynamics = ConservedGaussianJet(
        source,
        ambient_temperature=295.0,
        ambient_pressure=101325.0,
        ambient_density=1.197,
        fuel_molecular_weight=0.00201588,
        ambient_molecular_weight=0.02896546,
        fuel_heat_capacity=14294.8,
        ambient_heat_capacity=1006.2,
        radial_points=41,
    )
    return IndependentEnergyCrosswind(
        jetplume,
        thermodynamics,
        quadrature_points=16,
        energy_transport="total",
        ground_interaction=ground,
        ground_heat_transfer_coefficient_w_m2_k=ground_heat_coefficient,
        ground_surface_temperature_k=ground_surface_temperature,
    )


def contact_state(*, density=1.1, attachment=1.0):
    return np.array([density, 0.02, 0.03, 0.08, 15.0, 2.0, 0.10, attachment])


def test_wall_jet_requires_explicit_ground_geometry():
    with pytest.raises(ValueError, match="requires geometry"):
        FiniteWallJetCrosswind(base_model(ground="free"))
    with pytest.raises(ValueError, match="positive"):
        WallJetTransitionConfig(response_depths=0.0)


def test_attachment_suppresses_lift_and_surface_shear_removes_momentum():
    wall = FiniteWallJetCrosswind(base_model())
    attached = contact_state(attachment=1.0)
    released = contact_state(attachment=0.0)
    attached_source = wall.source_terms(attached)
    released_source = wall.source_terms(released)
    diagnostic = wall.diagnostics(attached)

    assert 0.0 < diagnostic.geometric_contact_fraction < 1.0
    assert 0.0 < diagnostic.vertical_release_fraction < 1.0
    assert diagnostic.wall_shear_force_per_m > 0.0
    assert attached_source[2] < released_source[2]
    assert attached_source[3] < released_source[3]


def test_lift_richardson_releases_a_light_wall_jet_continuously():
    wall = FiniteWallJetCrosswind(base_model())
    dense = wall.diagnostics(contact_state(density=1.4, attachment=0.5))
    light = wall.diagnostics(contact_state(density=0.4, attachment=0.5))

    assert dense.lift_richardson == 0.0
    assert dense.equilibrium_attachment_fraction > 0.0
    assert light.lift_richardson > wall.config.critical_richardson
    assert light.equilibrium_attachment_fraction == pytest.approx(0.0)


def test_no_contact_recovers_free_vertical_and_streamwise_sources():
    geometry = base_model(ground="geometry")
    free = base_model(ground="free")
    wall = FiniteWallJetCrosswind(geometry)
    state = contact_state(attachment=0.0)
    state[6] = 5.0

    assert wall.diagnostics(state).geometric_contact_fraction == 0.0
    assert wall.source_terms(state) == pytest.approx(
        free.source_terms(state[:7]), rel=1.0e-13, abs=1.0e-13
    )


def test_ground_heat_boundary_is_opt_in_contact_only_and_one_directional():
    warm = base_model(
        ground_heat_coefficient=12.0,
        ground_surface_temperature=295.0,
    )
    contact = contact_state()[:7]
    adiabatic = warm.source_terms(contact)
    heated = warm.source_terms(contact, include_ground_heat_transfer=True)

    assert warm.ground_heat_source(contact) > 0.0
    assert heated[:4] == pytest.approx(adiabatic[:4])
    assert heated[4] - adiabatic[4] == pytest.approx(
        warm.ground_heat_source(contact)
    )

    airborne = contact.copy()
    airborne[6] = 5.0
    assert warm.ground_heat_source(airborne) == 0.0

    cold = base_model(
        ground_heat_coefficient=12.0,
        ground_surface_temperature=20.0,
    )
    assert cold.ground_heat_source(contact) == 0.0


def test_ground_heat_boundary_requires_physical_declared_inputs():
    with pytest.raises(ValueError, match="requires a surface temperature"):
        base_model(ground_heat_coefficient=1.0)
    with pytest.raises(ValueError, match="non-negative"):
        base_model(ground_heat_coefficient=-1.0)


def test_eight_state_observation_adapter_matches_the_base_section():
    base = base_model()
    wall = FiniteWallJetCrosswind(base)
    state = contact_state(attachment=0.4)

    assert wall.section_widths(state) == pytest.approx(
        base.section_widths(state[:7])
    )
    assert wall.point_mole_fraction(state, 0.2, 0.4) == pytest.approx(
        base.point_mole_fraction(state[:7], 0.2, 0.4)
    )
    assert wall.point_temperature(state, 0.2, 0.4) == pytest.approx(
        base.point_temperature(state[:7], 0.2, 0.4)
    )


def test_short_wall_jet_run_keeps_attachment_bounded_and_audits_fluxes():
    wall = FiniteWallJetCrosswind(
        base_model(), WallJetTransitionConfig(response_depths=2.0)
    )
    initial = contact_state(density=1.1, attachment=0.8)
    result = wall.solve(
        initial,
        maximum_distance=0.03,
        maximum_step=0.005,
        relative_tolerance=2.0e-6,
    )

    assert np.all((result.attachment_fraction >= 0.0) & (result.attachment_fraction <= 1.0))
    assert result.states.shape[1] == 8
    assert result.maximum_relative_balance_residual < 2.0e-3
