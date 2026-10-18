"""Apply conservative interface mass transfers to scalar fields."""
from __future__ import annotations
import math
import numpy as np

def apply_interface_mass_to_field(field_kg_m3: np.ndarray, *,
                                  target_mask: np.ndarray,
                                  incoming_mass_kg: float,
                                  cell_volume_m3: float) -> np.ndarray:
    """Return a field with incoming species mass distributed over target cells."""
    field = np.asarray(field_kg_m3, dtype=float).copy()
    mask = np.asarray(target_mask, dtype=bool)
    volume = float(cell_volume_m3)
    mass = float(incoming_mass_kg)
    if field.shape != mask.shape:
        raise ValueError("target_mask must match field shape")
    if not np.all(np.isfinite(field)) or np.any(field < 0):
        raise ValueError("field must be finite and non-negative")
    if not math.isfinite(volume) or volume <= 0 or not math.isfinite(mass) or mass < 0:
        raise ValueError("cell volume and incoming mass must be finite and non-negative")
    count = int(np.count_nonzero(mask))
    if count == 0 and mass > 0:
        raise ValueError("positive incoming mass requires at least one target cell")
    if count:
        field[mask] += mass / (count * volume)
    return field

__all__ = ["apply_interface_mass_to_field"]
