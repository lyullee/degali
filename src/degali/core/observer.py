"""Transient releases: the observer formulation.

Ports ``UIT.FOR``, ``TUPF.FOR`` and ``OB.FOR``.

Why observers
-------------
A steady release gives a steady plume, and ``DEG2S`` integrates it in
*distance*.  A release that starts, changes and stops does not, and DEGADIS
handles that by switching to a Lagrangian frame: it releases ``NOBS``
imaginary observers over the source at different times and follows each one
downwind, accumulating the contaminant it collects on the way.  Each observer
becomes its own quasi-steady plume, and the concentration at a fixed point and
time is read off whichever observer happens to be passing.

An observer is a vertical slab of the cloud, spanning the full crosswind
width, that travels at the cloud's own advection speed rather than the wind's.
That speed follows from the power-law wind profile and the way the cloud
deepens as it entrains, giving

.. math:: u_i(t) = (1+\\alpha)\\, \\aleph^{1+\\alpha} (t - t_0)^{\\alpha},
          \\qquad
          x_i(t) = \\left[\\aleph (t - t_0)\\right]^{1+\\alpha} - R_{max}

with :math:`\\aleph` (``ALEPH``) the wind-speed constant that ``DEG1`` computes
at the end of the source calculation, and the origin shifted so that an
observer released at :math:`t_0` starts at the upwind edge of the largest the
source ever gets.

What is integrated
------------------
Five quantities accumulate as the observer crosses the source, all *rates*
because the observer is a slab of finite width:

=============  =========================================================
``hwidth``     slab half-width, m
``mrate``      total mass collected
``crate``      contaminant mass collected
``bdarate``    dry-air mass collected
``hrate``      enthalpy collected
=============  =========================================================

integrated in *time*, from ``TUP`` (when the observer crosses the upwind edge
of the secondary source) to ``TDOWN`` (the downwind edge).  Both are found by
root-finding on :func:`edge`, the difference between the observer position and
the source radius at that instant, with :func:`~degali.core.numerics.limit`
narrowing the bracket first because the source radius is itself a function of
time.

An observer that is entirely outside the source at some instant collects
nothing, and ``OB`` returns zero derivatives while leaving its diagnostics
untouched -- ``OBOUT`` refuses to copy them forward in that case, so the last
*valid* layer state survives to be used as the initial condition for the
downwind integration.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .atmosphere import richardson_star
from .constants import VKC
from .entrainment import phi
from .numerics import limit, zbrent

#: State indices, from ``OB``'s ``DATA HWIDTH/1/,Mrate/2/,Crate/3/,
#: BDArate/4/,Hrate/5/``.
O_HWIDTH, O_MRATE, O_CRATE, O_BDARATE, O_HRATE = range(5)


@dataclass
class ObserverKinematics:
    """The Lagrangian frame: ``UIT``, ``XIT`` and ``T0OB``."""

    alpha: float
    aleph: float  #: wind-speed constant from ``DEG1``
    rmax: float  #: largest source radius reached, m

    @property
    def alpha1(self) -> float:
        return self.alpha + 1.0

    def velocity(self, t: float, t0: float) -> float:
        """``UIT``: observer speed, m/s."""
        return self.alpha1 * self.aleph**self.alpha1 * (t - t0) ** self.alpha

    def position(self, t: float, t0: float) -> float:
        """``XIT``: observer position, m, measured from the source centre.

        Before release the observer sits at the upwind edge; the original
        returns ``-RMAX`` rather than extrapolating backwards.
        """
        if t - t0 <= 0.0:
            return -self.rmax
        return (self.aleph * (t - t0)) ** self.alpha1 - self.rmax

    def release_time(self, x: float, t: float) -> float:
        """``T0OB``: the release time of the observer at ``x`` at time ``t``.

        The inverse of :meth:`position`, guarded at the upwind edge: when
        ``x`` is within 0.1 % of ``-RMAX`` the exponent would amplify
        round-off without limit, so the original returns ``t`` unchanged.
        """
        arg = 0.0
        check = abs(abs(x) - abs(self.rmax)) / (abs(x) + abs(self.rmax))
        if check > 0.001:
            arg = (x + self.rmax) ** (1.0 / self.alpha1) / self.aleph
        return t - arg


def edge(time: float, t0: float, direction: float, kin: ObserverKinematics, radg) -> float:
    """Port of ``EDGE``: observer position minus the source edge, m.

    ``direction`` is ``-1`` for the upwind edge and ``+1`` for the downwind
    one.  ``radg`` is the source-radius vector from ``/GEN3/``, interpolated
    in time.
    """
    return kin.position(time, t0) - direction * radg(time)


def crossing_upwind(
    t0: float, kin: ObserverKinematics, radg, *, tol: float = 5.0e-4
) -> float:
    """Port of ``TUPF``: when this observer crosses the upwind source edge.

    The initial bracket is the observer's release time and the time it reaches
    ``x = 0``; ``LIMIT`` then narrows it in twentieths before ``ZBRENT``
    finishes the job.
    """
    tmax = kin.rmax ** (1.0 / kin.alpha1) / kin.aleph + t0
    tmin = max(t0, 0.0)
    f = lambda t: edge(t, t0, -1.0, kin, radg)
    dt = (tmax - tmin) / 20.0
    hi, lo = limit(f, tmin, dt, tmax, tmin)
    return zbrent(f, lo, hi, tol)


def crossing_downwind(
    t0: float, kin: ObserverKinematics, radg, *, tol: float = 5.0e-4
) -> float:
    """Port of ``TDNF``: when this observer crosses the downwind source edge.

    Bracketed between the time it reaches ``x = 0`` and the time it reaches
    ``x = +RMAX``; ``LIMIT`` starts from the *upper* end here, not the lower.
    """
    tmin = max(kin.rmax ** (1.0 / kin.alpha1) / kin.aleph + t0, 0.0)
    tmax = (2.0 * kin.rmax) ** (1.0 / kin.alpha1) / kin.aleph + t0
    f = lambda t: edge(t, t0, +1.0, kin, radg)
    dt = (tmax - tmin) / 20.0
    hi, lo = limit(f, tmax, dt, tmax, tmin)
    return zbrent(f, lo, hi, tol)


@dataclass
class ObserverState:
    """Layer diagnostics, ``PRMT(8..12)`` in ``OB``."""

    cclay: float
    wclay: float
    walay: float
    enthlay: float
    rholay: float
    outside: bool  #: True when the observer sees no gas at this instant


class Observer:
    """The ``OB`` derivative routine for one observer."""

    def __init__(
        self, thermo, vectors, kin: ObserverKinematics, *, params,
        gammaf, ustar, rhoa, u0, z0, humid,
    ):
        self.th = thermo
        self.v = vectors
        self.kin = kin
        self.p = params
        self.gammaf = gammaf
        self.ustar = ustar
        self.rhoa = rhoa
        self.u0 = u0
        self.z0 = z0
        self.humid = humid
        #: Last *valid* layer state; ``OBOUT`` only advances this when the
        #: observer is over the source.
        self.last: ObserverState | None = None

    def derivatives(self, time, y, dery, t0, rlen) -> ObserverState:
        """Port of ``OB``."""
        kin = self.kin
        adiabatic = self.th.ambient.isofl == 1 or self.th.ambient.ihtfl == 0
        xi = kin.position(time, t0)
        rg = self.v.radg_at(time)

        dery[:] = 0.0
        if abs(xi) >= rg:
            # outside the source: collect nothing, and keep the last diagnostics
            state = ObserverState(0.0, 0.0, 0.0, 0.0, 0.0, outside=True)
            return state

        bipr = math.sqrt(rg * rg - xi * xi)
        ui = kin.velocity(time, t0)
        q = self.v.qstr_at(time)
        wc = self.v.srcwc_at(time)
        wa = self.v.srcwa_at(time)
        enth = self.v.srcen_at(time)

        wclay = y[O_CRATE] / y[O_MRATE]
        walay = y[O_BDARATE] / y[O_MRATE]
        enthl = y[O_HRATE] / y[O_MRATE] if not adiabatic else 0.0

        mix = self.th.properties(wclay, walay, enthl, ifl=1)
        rholay = mix.rho
        cclay = wclay * rholay

        cc = cclay * self.p.dellay
        # the centreline density is extrapolated from the layer with a
        # constant density-concentration slope, not looked up again
        rho = self.p.dellay * (rholay - self.rhoa) + self.rhoa

        szob = self.z0 * (
            y[O_CRATE] / y[O_HWIDTH] / cc / (self.u0 * self.z0 / kin.alpha1)
        ) ** (1.0 / kin.alpha1)
        heff = self.gammaf / kin.alpha1 * szob
        hlay = self.p.dellay * heff
        # the layer and the centreline give different Richardson numbers; the
        # larger wins, so entrainment is suppressed by whichever is more stable
        ristr = max(
            richardson_star(rho, self.rhoa, heff, self.ustar),
            richardson_star(rholay, self.rhoa, hlay, self.ustar),
            0.0,
        )
        welay = self.p.dellay * VKC * self.ustar * kin.alpha1 / phi(
            ristr, 0.0, self.p.iphifl
        )

        da = ui * bipr
        dery[O_HWIDTH] = da / rlen
        dery[O_CRATE] = q * da
        dery[O_MRATE] = (q / wc + self.rhoa * welay) * da
        dery[O_BDARATE] = (q * wa / wc + self.rhoa * welay / (1.0 + self.humid)) * da
        if not adiabatic:
            dery[O_HRATE] = q * enth / wc * da

        return ObserverState(cclay, wclay, walay, enthl, rholay, outside=False)

    def record(self, state: ObserverState) -> None:
        """Port of ``OBOUT``: advance the saved layer state, unless outside."""
        if not state.outside:
            self.last = state

    def initial_conditions(self, tup: float, xdown: float, xup: float, sz0er: float):
        """Port of the ``SSSUP`` block that seeds each observer.

        The slab starts as a sliver: half-width ``SZ0ER`` times half the
        source length it will traverse, with the material it would already
        hold at that width taken from the source vectors at ``TUP``.
        """
        rlen = xdown - xup
        rrr = rlen / 2.0
        y = np.zeros(5)
        y[O_HWIDTH] = sz0er * rrr
        crate = y[O_HWIDTH] * rrr * self.v.qstr_at(tup)
        y[O_MRATE] = crate / self.v.srcwc_at(tup)
        y[O_CRATE] = crate
        y[O_BDARATE] = y[O_MRATE] * self.v.srcwa_at(tup)
        y[O_HRATE] = y[O_MRATE] * self.v.srcen_at(tup)
        return y, rlen


def release_times(
    kin: ObserverKinematics, vectors, tend: float, nobs: int
) -> np.ndarray:
    """Port of the ``SSSUP`` block that spaces the observers in time.

    The earliest useful release time is the one that puts an observer at the
    upwind edge when the source first appears -- but for a low-wind case that
    grows a blanket, a *later* point on the radius history can imply an even
    earlier release, so every recorded radius is checked. The latest is the
    one that reaches the downwind edge exactly as the source ends. The
    observers are then spread evenly across that window, offset by half a
    spacing so none sits on an endpoint.
    """
    t01 = kin.release_time(vectors.radg_at(0.0), 0.0)
    for t, r in zip(vectors.time, vectors.radg):
        t01 = min(t01, kin.release_time(r, t))
    xend = vectors.radg_at(tend)
    t0f = kin.release_time(-xend, tend)
    dtob = (t0f - t01) / float(nobs)
    return t01 + dtob / 2.0 + dtob * np.arange(nobs)
