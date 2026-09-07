"""Opt-in derivative consistency for the existing atmospheric spread widths.

This does not validate the inherited atmospheric correlations below one metre.
It differentiates the widths already used by the research geometry for x>0.
"""
from __future__ import annotations

import math

from .energy_crosswind import IndependentEnergyCrosswind
from .enthalpy_profile import GaussianEnthalpyCrosswind


class _ContinuousAmbientSpreading:
    def _geometry(self, state):
        geometry = super()._geometry(state)
        x, theta = float(state[5]), float(state[3])
        if x <= 0:
            raise ValueError('continuous ambient-spreading option requires x>0')
        if x > 1.:
            return geometry
        jp = self.jetplume
        sya, sza = geometry[7:9]
        ct = math.cos(theta)
        dy = sya / x * jp.betay * ct
        dz = sza / x * (jp.betaz + 2. * jp.gammaz * math.log(x)) * ct
        return (*geometry[:9], dy, dz)


class ContinuousAmbientDensityCrosswind(_ContinuousAmbientSpreading, IndependentEnergyCrosswind):
    """Density-Gaussian research model with continuous background growth."""


class ContinuousAmbientEnthalpyCrosswind(_ContinuousAmbientSpreading, GaussianEnthalpyCrosswind):
    """Enthalpy-Gaussian variant; still requires its own passed interface."""


def with_continuous_ambient_spreading(original):
    """Build a separate model with the original thermodynamics and settings."""
    classes = {IndependentEnergyCrosswind: ContinuousAmbientDensityCrosswind,
               GaussianEnthalpyCrosswind: ContinuousAmbientEnthalpyCrosswind}
    if type(original) not in classes:
        raise ValueError('unsupported original model; do not silently remove other physics')
    candidate = classes[type(original)](original.jetplume, original.thermodynamics,
        quadrature_points=original.quadrature_points, energy_transport=original.energy_transport,
        houf_width_mapping=original.houf_width_mapping, ground_interaction=original.ground_interaction,
        thermal_relaxation_rate=original.thermal_relaxation_rate,
        phase_transition_lag_rate=original.phase_transition_lag_rate,
        turbulence_heat_exchange_rate=original.turbulence_heat_exchange_rate)
    candidate.velocity_shape_exponent = original.velocity_shape_exponent
    return candidate
