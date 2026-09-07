"""One formulation from dense ground cloud to risen plume.

The structural limit
--------------------
An LH₂ release passes continuously through neutral buoyancy: dense and
slumping near the source, buoyant and rising a few metres later.  Nothing in
the physics switches at that point.  But the models available do:

* :class:`~degali.core.closures.DegadisClosure` has no vertical momentum
  equation at all, so a cloud that becomes light is frozen on the ground;
* :class:`~degali.addons.buoyant.BuoyantClosure` adds the momentum but keeps
  DEGADIS's ground-layer entrainment, which is calibrated for air entering
  through the top of a thin layer and gives rise velocities above 20 m/s once
  the cloud has left the ground;
* :class:`~degali.addons.liftoff.LiftoffPlume` has the perimeter entrainment
  that bounds the rise, but no ground-layer entrainment at all, so it cannot
  describe the cloud before it lifts.

Using them in sequence needs a handover, and a handover needs a threshold --
here a height -- which is a free parameter standing in for physics.  Results
depend on it, and neither model is valid at it: the first is out of
calibration above the threshold and the second below.

What this does instead
----------------------
One closure, no threshold.  The cloud always has a cross-section with a
*ground-contact fraction*

.. math:: f = \\frac{L_{ground}}{L_{ground} + L_{free}}

which is 1 for a cloud lying flat on the ground and 0 once it has cleared it,
and which the geometry already computes.  Entrainment is then

.. math:: E = (1 - f)\\, E_{perimeter} + f\\, E_{layer}

so a ground-hugging cloud entrains through its top as DEGADIS says, a risen
plume entrains through its whole perimeter as the URAHFREP model says, and a
cloud in between does both in proportion to how much of it is where.  The
same vertical momentum equation runs throughout, with buoyancy scaled by the
same fraction so that a cloud resting on the ground cannot lift the part of
itself that is still on it.

Nothing here is fitted.  Both entrainment laws are the ones already validated
in their own regimes, and the blend is the geometry.

What it does not fix
--------------------
Three limits remain, and they are properties of integral models rather than of
this one:

* **The cross-section carries one thermodynamic state.**  Giannissi and
  co-workers find that condensation of atmospheric humidity dominates LH₂
  cloud buoyancy and happens in a *shell* where cold hydrogen meets moist air.
  :meth:`~degali.addons.liftoff.LiftoffPlume._section_density` integrates
  the mixing line over the concentration profile, which is closer than a
  single mean state, but every annulus is still assumed to be at adiabatic
  equilibrium.  A shell reaction is not that.
* **Air condensation is not modelled.**  Below 90 K oxygen condenses and below
  77 K nitrogen does, removing mass from the gas and releasing latent heat.
  On the LH₂ mixing line that is above 81 mole per cent -- near the source,
  which is exactly where the dense phase lives.
* **Equilibrium is assumed instantaneous.**  Across a 270 K temperature
  difference the condensation kinetics may not keep up with the mixing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core.closures import CloudState, DegadisClosure
from ..core.constants import GG, PI, VKC
from .liftoff import ALPHA, BETA, CD, GAMMA

#: State indices.
U_W, U_Z, U_R = range(3)


@dataclass
class UnifiedClosure:
    """Dense slumping and buoyant rise in one set of equations.

    Parameters
    ----------
    coefficient, alpha_wind, z0, delrmn
        The gravity-slumping constants, as
        :class:`~degali.core.closures.DegadisClosure` uses them.
    segment_length
        Width of the cloud's flat part, m.  For a cloud handed over from a
        source this is the width the source left it with.
    alpha, beta, gamma
        Perimeter-entrainment coefficients, from the URAHFREP report.
    drag
        Form drag on the rising cross-section.
    """

    coefficient: float
    alpha_wind: float
    z0: float
    segment_length: float
    delrmn: float = 0.0
    alpha: float = ALPHA
    beta: float = BETA
    gamma: float = GAMMA
    drag: float = CD

    n_states: int = 3
    state_names: tuple[str, ...] = ("w", "z", "radius")

    def __post_init__(self):
        self._dense = DegadisClosure(
            coefficient=self.coefficient, alpha=self.alpha_wind, z0=self.z0,
            delrmn=self.delrmn,
        )

    # -- geometry ----------------------------------------------------------

    def contact_fraction(self, radius: float, z: float) -> float:
        """How much of the cross-section perimeter is on the ground.

        One for a cloud lying flat, zero once the section has cleared, and a
        smooth transition between.  This is what replaces a handover height.
        """
        if radius <= 0.0:
            return 1.0
        if z >= radius:
            return 0.0
        if z <= -radius:
            return 1.0
        ratio = max(min(z / radius, 1.0), -1.0)
        ground = self.segment_length + 2.0 * radius * math.sqrt(
            max(1.0 - ratio * ratio, 0.0)
        )
        free = self.segment_length + 2.0 * radius * (PI - math.acos(ratio))
        total = ground + free
        return ground / total if total > 0.0 else 1.0

    # -- interface ---------------------------------------------------------

    def initial(self, state: CloudState) -> np.ndarray:
        """A cloud from a ground-level source starts flat, at rest."""
        return np.array([0.0, 0.0, max(state.heff, 0.1)])

    def _unpack(self, state: CloudState) -> tuple[float, float, float]:
        own = state.own if state.own is not None else np.zeros(3)
        return float(own[U_W]), float(own[U_Z]), max(float(own[U_R]), 1e-6)

    def lateral_rate(self, state: CloudState) -> float:
        """Gravity slumping while dense and in contact, nothing when clear.

        Both limits are already right: DEGADIS's slumping needs the cloud to
        be on the ground for gravity to spread it against, and a plume that
        has lifted spreads by entrainment rather than by slumping.
        """
        _w, z, radius = self._unpack(state)
        f = self.contact_fraction(radius, z)
        return f * self._dense.lateral_rate(state)

    def derivatives(self, state: CloudState) -> np.ndarray:
        w, z, radius = self._unpack(state)
        u = max(state.wind, 1e-6)
        f = self.contact_fraction(radius, z)
        d = np.zeros(3)

        # -- vertical momentum, blended --------------------------------------
        mass = max(state.mass_flux, 1e-12)
        buoyancy = GG * (state.rhoa - state.rho) * max(state.heff, 1e-6) / u
        # a cloud still resting on the ground can only lift the part that has
        # cleared, which is what the contact fraction measures
        buoyancy *= 1.0 - f
        width = self.segment_length + 2.0 * radius
        resistance = self.drag * state.rhoa / 2.0 * w * abs(w) * width
        d[U_W] = (buoyancy - resistance - w * state.mass_flux_rate) / mass

        dz = w / u
        if z <= 0.0 and dz < 0.0:
            dz = 0.0
            d[U_W] = max(d[U_W], 0.0)
        d[U_Z] = dz

        # -- the cross-section grows by entrainment ---------------------------
        # ground layer: air enters through the top of a thin layer, at the
        # velocity DEGADIS's Richardson-corrected law gives
        layer = VKC * state.ustar * (1.0 + self.alpha_wind)
        # free perimeter: shear, cross-flow and ambient turbulence
        theta = math.atan2(w, u)
        free = (
            self.alpha * abs(w)
            + self.beta * abs(u * math.sin(theta))
            + self.gamma * u
        )
        # one law, weighted by where the cloud actually is
        entrain = f * layer + (1.0 - f) * free
        d[U_R] = entrain / u
        return d

    def elevation(self, state: CloudState) -> float:
        _w, z, _r = self._unpack(state)
        return max(z, 0.0)

    def regime(self, state: CloudState) -> str:
        """Where the cloud is, in words."""
        _w, z, radius = self._unpack(state)
        f = self.contact_fraction(radius, z)
        if f > 0.95:
            return "on the ground"
        if f > 0.05:
            return "lifting"
        return "airborne"
