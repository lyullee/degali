import numpy as np
import pytest
from degali.addons.interface_mass_application import apply_interface_mass_to_field

def test_interface_mass_application_preserves_requested_mass():
    field = np.zeros((2, 2, 2))
    updated = apply_interface_mass_to_field(
        field, target_mask=np.array([[[1, 0], [0, 1]], [[0, 0], [1, 0]]], dtype=bool),
        incoming_mass_kg=0.3, cell_volume_m3=0.5)
    assert np.sum(updated) * .5 == pytest.approx(.3)

def test_interface_mass_application_rejects_empty_target_for_positive_mass():
    with pytest.raises(ValueError, match="target cell"):
        apply_interface_mass_to_field(np.zeros((1, 1, 1)),
            target_mask=np.zeros((1, 1, 1), dtype=bool), incoming_mass_kg=1,
            cell_volume_m3=1)
