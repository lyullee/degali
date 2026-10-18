import math
from dataclasses import replace

import numpy as np
import pytest

from degali.addons.fast_obstacle_transport_3d import _source_weights_xyz
from degali.addons.high_pressure_source_handoff_3d import (
    ambient_pressure_jet_from_energy_nozzle,
    declared_ambient_pressure_hydrogen_jet,
    low_mach_handoff_from_energy_nozzle,
)
from degali.addons.notional import energy_conserving_notional_nozzle


@pytest.fixture(scope="module")
def hytunnel_nozzle():
    return energy_conserving_notional_nozzle(
        fluid="Hydrogen",
        storage_temperature=286.15,
        storage_pressure=11.515e6,
        mass_flow=0.005,
        orifice_diameter=0.001,
        ambient_pressure=98_500.0,
        allow_supercritical_gas=True,
    )


@pytest.fixture(scope="module")
def hytunnel_handoff(hytunnel_nozzle):
    return low_mach_handoff_from_energy_nozzle(
        hytunnel_nozzle,
        source_id="hytunnel-1mm-5g-3p5ms-coflow",
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_axial_velocity_m_s=3.5,
        maximum_mach_number=0.30,
    )


def test_hytunnel_handoff_closes_fluxes_before_resolved_injection(hytunnel_handoff):
    handoff = hytunnel_handoff
    assert handoff.numerically_admissible
    assert handoff.maximum_relative_conservation_residual < 1.0e-12
    assert handoff.mach_number < 0.30
    assert handoff.relative_density_departure == pytest.approx(-0.10)
    assert handoff.entrained_air_mass_flow_kg_s > handoff.source.mass_flow_kg_s
    assert handoff.plane_momentum_flux_n == pytest.approx(
        handoff.source_increment_momentum_flux_n
        + handoff.entrained_air_mass_flow_kg_s
        * handoff.ambient_axial_velocity_m_s
    )
    assert handoff.source_specific_momentum_m_s == pytest.approx(
        handoff.source.momentum_flux_n / handoff.source.mass_flow_kg_s
    )
    assert 0.5 < handoff.handoff_distance_m < 0.8
    assert 0.15 < handoff.equivalent_diameter_m < 0.25
    assert handoff.external_validation_required


def test_handoff_maps_to_a_resolved_cardinal_plane_without_losing_source_flux(
    hytunnel_handoff,
):
    fields = hytunnel_handoff.resolved_source_fields(
        nozzle_position_m=(0.0, 0.45, 0.45),
        grid_spacing_m=0.0064,
        duration_s=2.0,
    )
    kwargs = fields.as_config_kwargs()
    assert fields.spatially_resolved
    assert fields.cells_across_equivalent_diameter > 29.0
    assert fields.source_position_m[0] == pytest.approx(
        hytunnel_handoff.handoff_distance_m
    )
    assert fields.source_sigma_xyz_m[0] == pytest.approx(0.0064 / math.sqrt(12.0))
    assert fields.source_sigma_xyz_m[1:] == pytest.approx((
        hytunnel_handoff.moment_matched_radial_sigma_m,
        hytunnel_handoff.moment_matched_radial_sigma_m,
    ))
    assert fields.source_schedule.released_mass_kg == pytest.approx(0.010)
    assert kwargs["active_momentum_solver"] == "finite_volume_projection"
    assert kwargs["thermal_buoyancy_coupling"] is True
    assert kwargs["source_momentum_velocity_m_s"] == pytest.approx((
        hytunnel_handoff.source_specific_momentum_m_s, 0.0, 0.0,
    ))


def test_handoff_refuses_an_under_resolved_or_rotated_source(hytunnel_handoff):
    with pytest.raises(ValueError, match="under-resolved"):
        hytunnel_handoff.resolved_source_fields(
            nozzle_position_m=(0.0, 0.45, 0.45),
            grid_spacing_m=0.10,
            duration_s=1.0,
        )

    # The physical handoff supports arbitrary directions, while the exact
    # cell-integrated receiving quadrature deliberately accepts cardinal axes.
    from degali.addons.high_pressure_source_handoff_3d import (
        low_mach_handoff_from_ambient_pressure_jet,
    )
    rotated = low_mach_handoff_from_ambient_pressure_jet(
        hytunnel_handoff.source,
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        axis_unit_vector=(1.0, 1.0, 0.0),
    )
    with pytest.raises(ValueError, match="Cartesian-axis"):
        rotated.resolved_source_fields(
            nozzle_position_m=(0.0, 0.45, 0.45),
            grid_spacing_m=0.0064,
            duration_s=1.0,
        )


def test_published_source_requires_a_mass_consistent_complete_plane():
    with pytest.raises(ValueError, match="inconsistent with mass flow"):
        declared_ambient_pressure_hydrogen_jet(
            source_id="diameter-only-is-not-enough",
            source_model="published",
            mass_flow_kg_s=0.005,
            diameter_m=0.00725,
            velocity_m_s=2_000.0,
            density_kg_m3=0.167,
            temperature_k=143.0,
            ambient_pressure_pa=98_500.0,
        )


def test_anisotropic_cell_integrated_source_is_exactly_normalised():
    x = y = z = np.linspace(0.05, 0.95, 10)
    solid = np.zeros((10, 10, 10), dtype=bool)
    solid[4, 4, 4] = True
    weights = _source_weights_xyz(
        x, y, z, solid,
        (0.45, 0.55, 0.45),
        (0.02, 0.08, 0.12),
    )
    assert np.sum(weights) == pytest.approx(1.0, abs=1.0e-15)
    assert weights[4, 4, 4] == 0.0
    assert np.all(weights >= 0.0)


def test_energy_nozzle_conversion_rejects_an_inadmissible_state(hytunnel_nozzle):
    source = ambient_pressure_jet_from_energy_nozzle(hytunnel_nozzle)
    assert source.relative_mass_residual < 1.0e-9
    assert source.momentum_flux_n == pytest.approx(hytunnel_nozzle.momentum_flux)
    with pytest.raises(ValueError, match="not admissible"):
        ambient_pressure_jet_from_energy_nozzle(replace(
            hytunnel_nozzle, discharge_coefficient=1.1,
        ))
