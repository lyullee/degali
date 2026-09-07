"""A ground-level closure that lets a buoyant cloud leave the ground.

DEGADIS's own closure is right for LNG and wrong for hydrogen, and the reason
is arithmetic rather than opinion.  On the adiabatic mixing line for liquid
hydrogen into ambient air at 15 °C, the cloud is denser than air only above
about 99.9 mole per cent; the minimum density ratio, 0.73, sits at 72 mole per
cent, and at the lower flammable limit of 4 per cent the cloud is still
lighter than air.  The whole flammable range is buoyant.

DEGADIS handles that case by doing nothing: below ``DELRMN`` the spreading
rate is set to zero and there is no vertical momentum equation, so the cloud
is frozen in place at ground level and disperses passively.  That keeps a
flammable layer at head height which in reality has risen away.

This closure replaces that behaviour with the vertical momentum balance the
original lacks, taken from the URAHFREP ground-truncated plume model (see
:mod:`degali.addons.liftoff` for the same equations in a standalone plume).
Two states are carried alongside the cloud:

===========  ===========================================================
``w``        rise velocity of the centroid, m/s
``z``        centroid height above the ground, m
===========  ===========================================================

with

.. math::
    \\frac{dw}{dx} = \\frac{g\\,(\\rho_a - \\rho)}{\\rho\\, u},
    \\qquad \\frac{dz}{dx} = \\frac{w}{u}

While the cloud is dense it behaves exactly as DEGADIS does -- it slumps and
stays down, and ``z`` is held at zero rather than allowed to go negative,
because the ground is there.  The transition needs no switch: the same
equation gives a negative ``w`` for a dense cloud and the floor at ``z = 0``
does the rest.

Lateral spreading while buoyant
-------------------------------
A cloud that is rising is not slumping, so gravity no longer drives lateral
growth.  But it is not passive either: rising into a wind shears it and mixes
it faster than ambient turbulence alone.  The URAHFREP report addresses this
by scaling the entrainment with the vertical velocity,

.. math:: u_{entrain} = u_{entrain}(\\text{passive})\\,(1 + \\lambda\\, w/u_a)

and reports :math:`w/u_a \\propto \\sqrt{Ri^*}` for a two-dimensional buoyant
ground plume.  The same factor is applied here to the lateral rate, which
keeps a rising cloud spreading somewhat faster than a passive one without
inventing a separate mechanism for it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core.closures import CloudState, DegadisClosure

#: Coefficient in the enhanced-mixing factor ``1 + lambda w/u_a``.  The
#: URAHFREP report gives the form but not a fitted value; 1.0 is the neutral
#: choice and makes the enhancement equal to the velocity ratio itself.
LAMBDA = 1.0

#: State indices within this closure's own vector.
B_W, B_Z = 0, 1


@dataclass
class BuoyantClosure:
    """Gravity slumping while dense, buoyant rise once light.

    Parameters
    ----------
    coefficient, alpha, z0, delrmn
        As :class:`~degali.core.closures.DegadisClosure`: while the cloud is
        dense this closure *is* that one.
    enhanced_mixing
        Apply the ``1 + lambda w/u_a`` factor to lateral spreading once the
        cloud is rising.  Off gives a cloud that rises without spreading any
        faster, which is the more conservative choice for a flammable
        envelope.
    handover_height
        Height at which this closure declares the cloud lifted off, m.
        Beyond it the closure stops integrating and :attr:`lifted_off` is set.

        This is a real limit rather than a convenience.  The entrainment the
        downwind model supplies is calibrated for a cloud hugging the ground:
        air enters through the top of a thin layer.  A cloud that has left the
        ground entrains through its whole perimeter, several times faster, and
        that is what bounds a real plume's rise.  Continuing to integrate with
        ground-layer entrainment gives rise velocities above 20 m/s, an order
        of magnitude faster than any buoyant plume goes, because the momentum
        is not being diluted at anything like the right rate.

        So the honest division of labour is: this closure decides *whether and
        where* the cloud leaves the ground, and
        :class:`~degali.addons.liftoff.LiftoffPlume` -- which has perimeter
        entrainment -- carries it after that.  Results past the handover
        height should be taken from the plume model, not from here.
    """

    coefficient: float
    alpha: float
    z0: float
    delrmn: float = 0.0
    enhanced_mixing: bool = True
    handover_height: float = 2.0
    lam: float = LAMBDA

    n_states: int = 2
    state_names: tuple[str, ...] = ("w", "z")

    #: Set once the cloud passes :attr:`handover_height`.  Its value is the
    #: downwind distance at which that happened.
    lifted_off: float = float("nan")

    def __post_init__(self):
        self._dense = DegadisClosure(
            coefficient=self.coefficient, alpha=self.alpha, z0=self.z0,
            delrmn=self.delrmn,
        )

    # -- interface ---------------------------------------------------------

    def initial(self, state: CloudState) -> np.ndarray:
        """A cloud handed over from a source starts on the ground, at rest."""
        return np.zeros(2)

    def _rise(self, state: CloudState) -> tuple[float, float]:
        own = state.own if state.own is not None else np.zeros(2)
        return float(own[B_W]), float(own[B_Z])

    def lateral_rate(self, state: CloudState) -> float:
        slump = self._dense.lateral_rate(state)
        if not state.buoyant:
            return slump
        if not self.enhanced_mixing:
            return 0.0
        # Rising, so gravity no longer spreads it; but the rise shears the
        # cloud against the wind and mixes it faster than ambient turbulence
        # alone. The report scales entrainment by (1 + lambda w/u_a); the same
        # factor is applied to the lateral growth a passive cloud would have,
        # taken as sigma_y growing with distance.
        w, _z = self._rise(state)
        passive = state.sz / max(state.dist, 1.0)
        return self.lam * abs(w) / max(state.wind, 1e-6) * passive

    def derivatives(self, state: CloudState) -> np.ndarray:
        w, z = self._rise(state)
        u = max(state.wind, 1e-6)
        d = np.zeros(2)
        if z >= self.handover_height:
            # beyond here the ground-layer entrainment is out of calibration;
            # the plume model takes over
            if not math.isfinite(self.lifted_off):
                self.lifted_off = state.dist
            return d
        # the equation DEGADIS has no place for, written as a momentum
        # balance so that entrainment dilutes the rise
        mass = max(state.mass_flux, 1e-12)
        buoyancy = 9.81 * (state.rhoa - state.rho) * max(state.heff, 1e-6) / u
        d[B_W] = (buoyancy - w * state.mass_flux_rate) / mass
        dz = w / u
        # the ground is a floor: a dense cloud cannot sink through it, and a
        # cloud resting on it has no downward velocity to accumulate
        if z <= 0.0 and dz < 0.0:
            dz = 0.0
            d[B_W] = max(d[B_W], 0.0)
        d[B_Z] = dz
        return d

    def elevation(self, state: CloudState) -> float:
        _w, z = self._rise(state)
        return max(z, 0.0)


def make_buoyant(downwind, **kw) -> BuoyantClosure:
    """Build a :class:`BuoyantClosure` matching a configured downwind model.

    Takes the slumping constants from the model itself so that the dense
    behaviour is identical to what DEGADIS would have done.
    """
    return BuoyantClosure(
        coefficient=downwind.consts.c_spread,
        alpha=downwind.alpha,
        z0=downwind.case.z0,
        delrmn=downwind.p.delrmn,
        **kw,
    )
