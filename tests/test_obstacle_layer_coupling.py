import numpy as np
import pytest

from degali.addons.obstacle_layer_coupling import (
    HeightResolvedBoundaryPacket,
    collapse_packet_series,
    export_height_resolved_packets,
)


def test_height_resolved_packet_collapses_without_inventory_change():
    total = np.array([[[2.0, 1.0], [0.5, 0.25]], [[1.0, 0.0], [0.0, 0.5]]])
    species = 0.4 * total
    momentum = np.zeros((3, *total.shape))
    momentum[0] = 3.0 * total
    packets = export_height_resolved_packets(
        time_s=[1.0, 1.5], dt_s=.5, face='z_max',
        total_rate_kg_s=total, species_rate_kg_s=species,
        momentum_rate_kg_m_s=momentum,
    )
    collapsed = collapse_packet_series(packets)
    assert collapsed.total_mass_kg == pytest.approx(total.sum() * .5)
    assert collapsed.species_mass_kg == pytest.approx(species.sum() * .5)
    assert collapsed.momentum_kg_m_s[0] == pytest.approx((3.0 * total).sum() * .5)
    assert collapsed.species_mass_kg <= collapsed.total_mass_kg


def test_packet_rejects_species_creation_and_shape_mismatch():
    with pytest.raises(ValueError, match='species rate cannot exceed'):
        HeightResolvedBoundaryPacket(
            time_s=0., dt_s=.1, face='x_max',
            total_rate_kg_s=np.ones((2, 2)),
            species_rate_kg_s=np.full((2, 2), 1.1),
            momentum_rate_kg_m_s=np.zeros((3, 2, 2)),
        )
    with pytest.raises(ValueError, match='shape'):
        export_height_resolved_packets(
            time_s=[0., .1], dt_s=.1, face='x_max',
            total_rate_kg_s=np.ones((1, 2, 2)),
            species_rate_kg_s=np.ones((1, 2, 2)),
        )
