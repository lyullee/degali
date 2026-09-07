"""Turning observer profiles into cloud snapshots: a port of ``DEG3``.

Ports ``GETTS.FOR``, ``SORTS.FOR``, ``SORTS1.FOR`` and the ``DEG3`` driver.

The problem
-----------
:mod:`.transient` leaves ``NOBS`` profiles, each a curve of concentration
against distance for one observer, with an arrival time attached to every
point.  What a user wants instead is the *cloud at an instant*: at t = 60 s,
how far does the flammable region reach and how wide is it.

Every observer contributes exactly one point to each snapshot -- the place it
happens to be at that time.  So a snapshot is assembled by walking each
observer's profile, finding the pair of recorded points that straddle the
snapshot time, and interpolating between them.  With 30 observers a snapshot
has up to 30 points.

Two things then have to be fixed up:

**Order.** Observer 1 is released first and is therefore *furthest* downwind
at any given time, so the natural row order runs from far to near.  ``SORTS1``
reverses each column so distance increases down the snapshot.

**Along-wind dispersion.** Each observer's plume is quasi-steady and carries
no along-wind spreading of its own; the puff character comes only from the
observers being staggered. For a short release that under-predicts how much
the cloud smears out along the wind, so ``SORTS1`` convolves the snapshot with
a Gaussian of scale :math:`\\sigma_x = a\\,x^{p}`,

.. math::
    C^*(x) = \\frac{1}{\\sqrt{2\\pi}} \\int C(x')
             \\frac{1}{\\sigma_x}
             \\exp\\!\\left[-\\frac{(x-x')^2}{2\\sigma_x^2}\\right] dx'

integrated numerically over :math:`\\pm 4\\sigma_x`, and skipped entirely for
points closer than ``SIGXMD`` to the source, where the correction would be
meaningless. The coefficients come from the stability class.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import RT2, SQRTPI
from .downwind import series
from .numerics import gaminc
from .rkgst import Control, rkgst

#: Columns of a snapshot row.
SNAPSHOT_COLUMNS = (
    "dist", "dist0", "yc", "cc", "ccstr", "rho", "gamma", "temp", "sz", "sy", "b",
)


def sort_times(
    t0s: np.ndarray, *, rmax: float, aleph: float, alpha1: float,
    ert1: float = 0.0, erdt: float = 0.0, erntim: int = 0, check5: bool = False,
) -> np.ndarray:
    """Port of ``GETTS``: choose the instants to take snapshots at.

    With ``CHECK5`` the user's ``.ER3`` values are used directly.  Otherwise
    the window runs from when the *first* observer passes ``x = 2 R_max`` to
    when the *last* passes ``x = 6 R_max`` -- deliberately wide, because the
    right window depends on the wind speed and the original's stated intent is
    to "show the user where to look on the next go around".

    The interval is rounded to a whole number of seconds and floored at five,
    and at least four times are always produced.
    """
    if check5:
        t1 = ert1
        dt = erdt
        ntim = int(erntim)
    else:
        t1 = t0s[0] + (2.0 * rmax) ** (1.0 / alpha1) / aleph
        tf = t0s[-1] + (6.0 * rmax) ** (1.0 / alpha1) / aleph
        ntim = 20
        dt = 2.0 * (tf - t1) / float(ntim - 1)
        dt = float(round(dt))
        if dt < 5.0:
            dt = 5.0
            ntim = int((tf - t1) / dt) + 1
    t1 = float(round(t1))
    ntim = max(ntim, 4)
    return dt * np.arange(ntim) + t1


@dataclass
class Snapshot:
    """The cloud at one instant."""

    time: float
    rows: np.ndarray  #: (n, 11), see :data:`SNAPSHOT_COLUMNS`
    mass_above_ulc: float = 0.0  #: kg
    mass_above_llc: float = 0.0  #: kg

    def __len__(self) -> int:
        return len(self.rows)

    def column(self, name: str) -> np.ndarray:
        return self.rows[:, SNAPSHOT_COLUMNS.index(name)]

    @property
    def mass_between(self) -> float:
        return self.mass_above_llc - self.mass_above_ulc

    def distance_to(self, mole_fraction: float) -> float:
        """Furthest distance at which ``yc`` still exceeds ``mole_fraction``."""
        yc = self.column("yc")
        x = self.column("dist")
        above = yc >= mole_fraction
        if not above.any():
            return float("nan")
        i = int(np.max(np.flatnonzero(above)))
        if i == len(yc) - 1:
            return float(x[i])
        f = (yc[i] - mole_fraction) / (yc[i] - yc[i + 1])
        return float(x[i] + f * (x[i + 1] - x[i]))


class TimeSort:
    """The ``DEG3`` sort: observer profiles in, cloud snapshots out."""

    def __init__(
        self, transient, *, sigxco: float, sigxp: float, sigxmd: float,
        sigxfl: bool = True, alpha1: float, gammaf: float, table, gas, ambient,
    ):
        self.tr = transient
        self.sigxco = sigxco
        self.sigxp = sigxp
        self.sigxmd = sigxmd
        self.sigxfl = sigxfl
        self.alpha1 = alpha1
        self.gammaf = gammaf
        self.table = table
        self.gas = gas
        self.ambient = ambient

    # -- assembling a snapshot ---------------------------------------------

    def _interpolate(self, obs, t: float) -> np.ndarray | None:
        """One observer's contribution at time ``t``, or ``None``.

        The observer's recorded points carry arrival times; the pair that
        straddles ``t`` is interpolated linearly in every quantity.
        ``dist0`` is the observer's *first* recorded distance -- ``SORTS``
        sets ``TABLE(21)`` once, before the loop, and never updates it -- so
        ``x - x0`` is how far the plume has travelled from the source. That
        is what the along-wind correction scales with, not the spacing
        between output points.
        """
        rows = obs.rows
        if len(rows) < 2:
            return None
        times = rows[:, 1]
        if t < times[0] or t > times[-1]:
            return None
        j = int(np.searchsorted(times, t))
        j = min(max(j, 1), len(rows) - 1)
        t_prev, t_now = times[j - 1], times[j]
        frac = (t - t_prev) / (t_now - t_prev) if t_now != t_prev else 0.0
        prev, now = rows[j - 1], rows[j]
        interp = prev + (now - prev) * frac
        # dist, dist0, yc, cc, ccstr(filled later), rho, gamma, temp, sz, sy, b
        return np.array([
            interp[0], rows[0, 0], interp[2], interp[3], interp[3],
            interp[4], interp[5], interp[6], interp[7], interp[8], interp[9],
        ])

    def _assemble(self, t: float) -> np.ndarray:
        rows = [
            r for r in (self._interpolate(o, t) for o in self.tr.observers)
            if r is not None
        ]
        if not rows:
            return np.zeros((0, 11))
        arr = np.array(rows)
        # observer 1 is released first and so sits furthest downwind; reverse
        # so distance increases down the snapshot
        return arr[np.argsort(arr[:, 0])]

    # -- the along-wind correction -----------------------------------------

    def _correct_along_wind(self, rows: np.ndarray) -> None:
        """Port of the ``SORTS1``/``SORTSI`` convolution.  Modifies in place."""
        if not self.sigxfl or len(rows) < 2:
            rows[:, 4] = rows[:, 3]
            return

        dist = rows[:, 0]
        dist0 = rows[:, 1]
        cc = rows[:, 3]

        def concentration_at(x: float) -> tuple[float, float]:
            """Linear interpolation of ``cc`` and the dispersing distance."""
            i = int(np.searchsorted(dist, x))
            if i >= len(dist):
                return float(cc[-1]), float(x - dist0[-1])
            i = max(i, 1)
            im = i - 1
            denom = dist[im] - dist[i]
            frac = (x - dist[i]) / denom if denom != 0 else 0.0
            d0 = frac * (dist0[im] - dist0[i]) + dist0[i]
            return float(frac * (cc[im] - cc[i]) + cc[i]), float(x - d0)

        for i in range(len(rows)):
            spread = dist[i] - dist0[i]
            if spread < self.sigxmd:
                rows[i, 4] = cc[i]
                continue
            sx = self.sigxco * spread**self.sigxp
            half = 4.0 * sx
            lo = max(dist[0], dist[i] - half)
            hi = min(dist[-1], dist[i] + half)
            if hi <= lo:
                rows[i, 4] = cc[i]
                continue

            def integrand(x, y, dery, prmt):
                c, d = concentration_at(x)
                sxl = self.sigxco * d**self.sigxp if d > 0 else sx
                if sxl <= 0.0:
                    dery[0] = 0.0
                    return
                arg = (abs(dist[i] - x) / sxl) ** 2 / 2.0
                dery[0] = c / sxl if arg == 0.0 else (
                    c / sxl / math.exp(arg) if arg <= 80.0 else 0.0
                )

            def nothing(*_args):
                return None

            prmt = Control([lo, hi, half / 20.0, 0.003, 2.0 * sx] + [0.0] * 20)
            res = rkgst(integrand, nothing, prmt, [0.0], [1.0], ndim=1)
            rows[i, 4] = float(res.y[0]) / RT2 / SQRTPI

    # -- flammable mass ----------------------------------------------------

    def _flammable_mass(self, rows: np.ndarray) -> tuple[float, float]:
        """Port of the ``SORTS1`` block at label 430.

        Integrates the mass above each level of concern along the snapshot,
        using the same profile integrals as ``PSS`` but with the corrected
        centreline concentration.
        """
        from .downwind import _adiabat_m2

        adiabatic = self.ambient.isofl == 1 or self.ambient.ihtfl == 0
        a1 = self.alpha1
        hi = np.zeros(len(rows))
        low = np.zeros(len(rows))

        for j, r in enumerate(rows):
            cc, bb, sy, sz, gamma = r[4], r[10], r[9], r[8], r[6]
            if adiabatic:
                chi = self.table.from_mole_fraction(self.gas.ulc).cc
                clow = self.table.from_mole_fraction(self.gas.llc).cc
            else:
                chi = _adiabat_m2(self.table, self.gas, self.ambient, self.gas.ulc, gamma)
                clow = _adiabat_m2(self.table, self.gas, self.ambient, self.gas.llc, gamma)

            gamhi0 = 2.0 * cc * bb * sz / a1
            gammax = 2.0 * cc * sz * self.gammaf / a1 * (bb + SQRTPI / 2.0 * sy)
            if cc > clow > 0.0:
                wlow = math.log(cc / clow)
                low[j] = min(
                    gaminc(1.0 / a1, wlow) * gamhi0
                    + 2.0 * clow * sy * sz / a1 * series(wlow, a1),
                    gammax,
                )
            if cc > chi > 0.0:
                whi = math.log(cc / chi)
                hi[j] = min(
                    gaminc(1.0 / a1, whi) * gamhi0
                    + 2.0 * chi * sy * sz / a1 * series(whi, a1),
                    gammax,
                )

        dx = np.diff(rows[:, 0])
        mass_hi = float(np.sum((hi[1:] + hi[:-1]) / 2.0 * dx))
        mass_low = float(np.sum((low[1:] + low[:-1]) / 2.0 * dx))
        return mass_hi, mass_low

    # -- the driver ---------------------------------------------------------

    def run(self, times: np.ndarray, *, min_points: int = 3) -> list[Snapshot]:
        """Produce a snapshot for each requested time.

        Times at which fewer than three observers contribute are dropped:
        ``SORTS1`` calls one or two points "of little value" and trims the
        window from both ends rather than reporting a cloud built from them.
        """
        out: list[Snapshot] = []
        for t in times:
            rows = self._assemble(float(t))
            if len(rows) < min_points:
                continue
            self._correct_along_wind(rows)
            try:
                m_hi, m_low = self._flammable_mass(rows)
            except (OverflowError, ValueError):
                m_hi = m_low = 0.0
            out.append(
                Snapshot(time=float(t), rows=rows,
                         mass_above_ulc=m_hi, mass_above_llc=m_low)
            )
        return out
