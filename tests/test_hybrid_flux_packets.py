import numpy as np
import pytest

from degali.addons.hybrid_conservative_exchange import ConservativeInterfacePacket
from degali.addons.hybrid_flux_packets import (
    CartesianInterfacePlane,
    apply_incoming_momentum_characteristic_to_velocity,
    apply_incoming_momentum_to_velocity,
    apply_incoming_packet_to_fields,
    boundary_cell_mask,
    outgoing_packet_from_face_fluxes,
)


def test_outgoing_flux_packet_integrates_species_energy_and_momentum():
    total = [np.zeros((2, 2, 3)), np.zeros((2, 3, 2)), np.zeros((3, 2, 2))]
    species = [np.zeros_like(total[0]), np.zeros_like(total[1]), np.zeros_like(total[2])]
    enthalpy = [np.zeros_like(item) for item in total]
    momentum = [[np.zeros_like(item) for item in total] for _ in range(3)]
    total[0][:, :, -1] = 2.0
    species[0][:, :, -1] = 0.4
    enthalpy[0][:, :, -1] = 10.0
    momentum[0][0][:, :, -1] = 3.0
    plane = CartesianInterfacePlane(axis=0, face_index=-1, outward_sign=1.0, face_area_m2=0.5)

    packet = outgoing_packet_from_face_fluxes(
        total_mass_fluxes=total,
        species_mass_fluxes=species,
        plane=plane,
        time_step_s=0.2,
        sensible_enthalpy_fluxes=enthalpy,
        momentum_fluxes=momentum,
    )

    assert packet.total_mass_kg == pytest.approx(0.8)
    assert packet.species_mass_kg == pytest.approx(0.16)
    assert packet.sensible_enthalpy_j == pytest.approx(4.0)
    assert packet.momentum_kg_m_s == pytest.approx((1.2, 0.0, 0.0))


def test_packet_application_closes_field_integrals():
    field = np.zeros((2, 3, 4))
    enthalpy = np.zeros_like(field)
    mask = boundary_cell_mask(field.shape, axis=2, side="min", layers=1)
    packet = ConservativeInterfacePacket(.4, .2, 12.0, (1.0, 0.0, 0.0))
    updated, updated_enthalpy = apply_incoming_packet_to_fields(
        field, packet=packet, target_mask=mask, cell_volume_m3=0.5,
        sensible_enthalpy_j_m3=enthalpy,
    )
    assert np.sum(updated) * 0.5 == pytest.approx(.2)
    assert np.sum(updated_enthalpy) * 0.5 == pytest.approx(12.0)

    u, v, w = apply_incoming_momentum_to_velocity(
        np.zeros_like(field), np.zeros_like(field), np.zeros_like(field),
        packet=packet, target_mask=mask,
        density_kg_m3=np.ones_like(field), cell_volume_m3=0.5,
    )
    assert np.sum(u[mask]) * 0.5 == pytest.approx(1.0)
    assert np.all(v == 0.0) and np.all(w == 0.0)


def test_flux_packet_rejects_species_above_total():
    total = [np.zeros((1, 1, 2)), np.zeros((1, 2, 1)), np.zeros((2, 1, 1))]
    species = [np.zeros_like(item) for item in total]
    total[0][0, 0, -1] = 0.1
    species[0][0, 0, -1] = 0.2
    with pytest.raises(ValueError, match="exceeds outgoing total"):
        outgoing_packet_from_face_fluxes(
            total_mass_fluxes=total, species_mass_fluxes=species,
            plane=CartesianInterfacePlane(0, -1, 1.0, 1.0), time_step_s=1.0,
        )


def test_characteristic_momentum_maps_packet_velocity_without_accumulation():
    shape = (2, 2, 2)
    mask = boundary_cell_mask(shape, axis=2, side="min")
    packet = ConservativeInterfacePacket(
        total_mass_kg=2.0, species_mass_kg=0.5,
        momentum_kg_m_s=(0.4, 0.0, 0.0),
    )
    faces = (
        np.zeros((2, 2, 3)),
        np.zeros((2, 3, 2)),
        np.zeros((3, 2, 2)),
    )
    result = apply_incoming_momentum_characteristic_to_velocity(
        np.ones(shape) * 0.8,
        np.zeros(shape), np.zeros(shape),
        *faces,
        packet=packet, target_mask=mask,
        density_kg_m3=np.ones(shape), cell_volume_m3=1.0,
        plane=CartesianInterfacePlane(0, 0, -1.0, 1.0),
    )
    u, v, w, u_face, v_face, w_face = result
    np.testing.assert_allclose(u[mask], 0.2)
    np.testing.assert_allclose(u[~mask], 0.8)
    np.testing.assert_allclose(u_face[:, :, 0], 0.2)
    assert np.all(v == 0.0) and np.all(w == 0.0)
    assert np.all(v_face == 0.0) and np.all(w_face == 0.0)
