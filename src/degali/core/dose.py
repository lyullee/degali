"""Concentration at a fixed point: a port of ``DEG4``.

``DEG3`` answers "where is the cloud at time t".  ``DEG4`` answers the
complementary question, "what does an observer standing at x see, and when" --
the concentration *time history* at a receptor, which is what a toxic-exposure
assessment integrates to get a dose.

Ports ``GETTD.FOR`` and ``DOSOUT.FOR``.

How it differs from ``DEG3``
---------------------------
The sorting machinery is identical; only the choice of times changes.  ``DEG3``
picks instants spanning the whole cloud lifetime.  ``GETTD`` instead brackets
the window during which the cloud is passing *this* receptor: from when the
first observer crosses it to when the last one does, split into at most
``MAXNT`` steps of whole seconds.

Then, rather than reporting a whole snapshot, ``DOSOUT`` interpolates each
snapshot to the receptor's own distance, giving one row per time.  The
centreline is reported directly; off-axis receptors are evaluated from the
profile shape,

.. math::
    C(y, z) = C_c \\exp\\!\\left[-\\left(\\frac{z}{S_z}\\right)^{1+\\alpha}
              - \\left(\\frac{y - B}{S_y}\\right)^{2}\\right]

with the crosswind term dropped inside the flat core, :math:`|y| \\le B`.  Up
to four such off-axis points can be requested per receptor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .timesort import SNAPSHOT_COLUMNS, Snapshot, TimeSort

#: Columns of a receptor time-history row.
HISTORY_COLUMNS = ("time", "yc", "cc", "rho", "gamma", "temp", "b", "sz", "sy")


@dataclass
class Receptor:
    """A fixed point to report concentration at.

    ``offsets`` are up to four ``(y, z)`` pairs in metres; the centreline
    history is always produced, these are extra.
    """

    x: float  #: downwind distance, m
    offsets: list[tuple[float, float]] = field(default_factory=list)


@dataclass
class DoseHistory:
    """Concentration against time at one receptor."""

    receptor: Receptor
    rows: np.ndarray  #: (n, 9), see :data:`HISTORY_COLUMNS`
    #: Mole fraction at each requested off-axis point, ``(n, len(offsets))``.
    offaxis: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))

    def column(self, name: str) -> np.ndarray:
        return self.rows[:, HISTORY_COLUMNS.index(name)]

    @property
    def peak(self) -> tuple[float, float]:
        """Peak centreline mole fraction and the time it occurs."""
        if len(self.rows) == 0:
            return 0.0, float("nan")
        i = int(np.argmax(self.rows[:, 1]))
        return float(self.rows[i, 1]), float(self.rows[i, 0])

    def dose(self, exponent: float = 1.0, column: int = 0) -> float:
        r"""Toxic load :math:`\int C^n \, dt`, in (mole fraction)^n seconds.

        ``exponent`` is the probit exponent ``n``; ``column`` selects the
        centreline (``0``) or an off-axis point (``1``-based into
        ``offsets``). Integrated by the trapezium rule over the reported
        times, which is what the original's tabulated output supports.
        """
        if len(self.rows) < 2:
            return 0.0
        t = self.rows[:, 0]
        c = self.rows[:, 1] if column == 0 else self.offaxis[:, column - 1]
        return float(np.trapezoid(np.power(np.maximum(c, 0.0), exponent), t))


def receptor_times(
    t0s: np.ndarray, x: float, *, oodist: float, rmax: float, aleph: float,
    alpha1: float, maxnt: int = 40,
) -> np.ndarray:
    """Port of ``GETTD``: the window during which the cloud passes ``x``.

    Bracketed by the arrival of the first and last observers at that distance,
    with the step rounded to whole seconds and floored at one, and at least
    four times produced.
    """

    def ts(t0: float) -> float:
        return t0 + (x - oodist + rmax) ** (1.0 / alpha1) / aleph

    t1 = ts(t0s[0])
    tf = ts(t0s[-1])
    ntim = maxnt
    dt = max(float(round((tf - t1) / float(ntim - 1))), 1.0)
    ntim = min(int((tf - t1) / dt) + 1, maxnt)
    ntim = max(ntim, 4)
    t1 = float(round(t1))
    return dt * np.arange(ntim) + t1


class DoseRun:
    """The ``DEG4`` driver."""

    def __init__(self, sorter: TimeSort, *, oodist: float = 0.0, alpha1: float):
        self.sorter = sorter
        self.oodist = oodist
        self.alpha1 = alpha1

    def _at_receptor(self, snap: Snapshot, x: float) -> np.ndarray | None:
        """Port of the ``DOSOUT`` interpolation to the receptor distance.

        Returns ``None`` when the whole snapshot lies downwind of the
        receptor -- the cloud has already passed, and the original reports the
        record as missing rather than extrapolating backwards.
        """
        dist = snap.column("dist") + self.oodist
        if dist[0] > x or dist[-1] < x:
            return None
        j = int(np.searchsorted(dist, x))
        j = min(max(j, 1), len(dist) - 1)
        denom = dist[j] - dist[j - 1]
        frac = (x - dist[j - 1]) / denom if denom != 0 else 0.0

        def lerp(name: str) -> float:
            col = snap.column(name)
            return float(col[j - 1] + frac * (col[j] - col[j - 1]))

        return np.array([
            snap.time, lerp("yc"), lerp("ccstr"), lerp("rho"), lerp("gamma"),
            lerp("temp"), lerp("b"), lerp("sz"), lerp("sy"),
        ])

    def _offaxis(self, row: np.ndarray, offsets) -> np.ndarray:
        """Mole fraction at each ``(y, z)`` offset, from the profile shape."""
        _t, yc, cc, _rho, gamma, _temp, b, sz, sy = row
        out = np.zeros(len(offsets))
        adiabatic = (
            self.sorter.ambient.isofl == 1 or self.sorter.ambient.ihtfl == 0
        )
        for k, (y, z) in enumerate(offsets):
            if y <= 0.0 and z <= 0.0:
                continue
            arg = (z / sz) ** self.alpha1 if sz > 0 else 0.0
            if y > b and sy > 0:
                arg += ((y - b) / sy) ** 2
            if arg >= 80.0:
                continue
            c_local = cc / math.exp(arg)
            if adiabatic:
                out[k] = self.sorter.table.from_concentration(c_local).yc
            else:
                # off the adiabatic line: invert the linear density relation
                out[k] = _mole_fraction_from_gamma(
                    self.sorter.table, self.sorter.gas, self.sorter.ambient,
                    c_local, gamma,
                )
        return out

    def run(self, receptors, t0s: np.ndarray, *, rmax: float, aleph: float,
            maxnt: int = 40) -> list[DoseHistory]:
        """Build a time history at each receptor."""
        out: list[DoseHistory] = []
        for r in receptors:
            times = receptor_times(
                t0s, r.x, oodist=self.oodist, rmax=rmax, aleph=aleph,
                alpha1=self.alpha1, maxnt=maxnt,
            )
            snaps = self.sorter.run(times)
            rows = []
            offs = []
            for s in snaps:
                row = self._at_receptor(s, r.x)
                if row is None:
                    continue
                rows.append(row)
                if r.offsets:
                    offs.append(self._offaxis(row, r.offsets))
            out.append(
                DoseHistory(
                    receptor=r,
                    rows=np.array(rows) if rows else np.zeros((0, 9)),
                    offaxis=np.array(offs) if offs else np.zeros((0, 0)),
                )
            )
        return out


def _mole_fraction_from_gamma(table, gas, ambient, cc: float, gamma: float) -> float:
    """Port of ``ADIABAT(ifl=-1)``: composition from concentration and slope.

    Once the ground has heated the cloud off the adiabatic mixing line the
    tabulated density no longer applies, but the linear relation
    :math:`\\rho = \\rho_a + \\gamma c` still does, so the mass fraction follows
    directly and the mole fraction from the mixing rule.
    """
    from .constants import WMA, WMW

    ccl = max(cc, 0.0)
    denom = table.rhoa + ccl * gamma
    if denom <= 0.0:
        return 0.0
    wc = ccl / denom
    wa = (1.0 - (1.0 + ambient.humsrc) * wc) / (1.0 + ambient.humid)
    wm = 1.0 / (wc / gas.mw + wa / WMA + (1.0 - wa - wc) / WMW)
    return wm / gas.mw * wc
