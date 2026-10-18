"""Concentration-blind atmospheric surface-layer mixing closure.

The neutral Monin--Obukhov limit gives the momentum eddy viscosity

    nu_t(z) = kappa * u_star * z,

because the neutral momentum similarity function is one.  This module only
implements that evidence-supported neutral limit.  It deliberately does not
infer a Monin--Obukhov length from a Pasquill class and does not select a
turbulent Schmidt number from concentration observations.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


ATMOSPHERIC_SURFACE_LAYER_MIXING_SCHEMA = (
    "degali.atmospheric-surface-layer-mixing.v1"
)


@dataclass(frozen=True)
class NeutralSurfaceLayerMixing:
    friction_velocity_m_s: float
    roughness_length_m: float
    source_id: str
    von_karman_constant: float = 0.4

    def __post_init__(self) -> None:
        for name in (
            "friction_velocity_m_s", "roughness_length_m",
            "von_karman_constant",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
            object.__setattr__(self, name, value)
        if not isinstance(self.source_id, str) or not self.source_id.strip():
            raise ValueError("source_id must be a non-empty string")

    def eddy_viscosity_m2_s(
        self, z_m: np.ndarray, solid_mask: np.ndarray,
    ) -> np.ndarray:
        """Return neutral surface-layer momentum eddy viscosity by height."""
        z = np.asarray(z_m, dtype=float)
        if z.ndim != 1 or z.size < 1 or np.any(~np.isfinite(z)) or np.any(z < 0.0):
            raise ValueError("z_m must be a finite non-negative one-dimensional array")
        solid = np.asarray(solid_mask)
        if solid.dtype != np.bool_ or solid.ndim != 3 or solid.shape[0] != z.size:
            raise ValueError(
                "solid_mask must be a three-dimensional boolean array matching z_m"
            )
        profile = (
            self.von_karman_constant
            * self.friction_velocity_m_s
            * z
        )
        viscosity = np.broadcast_to(
            profile[:, None, None], solid.shape,
        ).copy()
        viscosity[solid] = 0.0
        return viscosity

    def record(self) -> dict[str, object]:
        return {
            "schema": ATMOSPHERIC_SURFACE_LAYER_MIXING_SCHEMA,
            "model": "neutral_monin_obukhov_surface_layer",
            "formula": "nu_t = kappa * u_star * z",
            "friction_velocity_m_s": self.friction_velocity_m_s,
            "roughness_length_m": self.roughness_length_m,
            "von_karman_constant": self.von_karman_constant,
            "source_id": self.source_id,
            "stability_function_phi_m": 1.0,
            "limitations": [
                "neutral surface-layer reference only",
                "no obstacle-generated turbulence production beyond resolved momentum",
                "no concentration-derived coefficient selection",
            ],
        }


__all__ = [
    "ATMOSPHERIC_SURFACE_LAYER_MIXING_SCHEMA",
    "NeutralSurfaceLayerMixing",
]
