"""The ``DEG2S`` driver: carry the cloud downwind to the level of concern.

Ports the ``DEG2S`` main program, ``PSSOUTSS``, ``SSGOUTSS`` and ``SSOUT``.
:mod:`.downwind` holds the two sets of derivatives; this module runs them, and
in particular handles the handover between them, which is the delicate part.

The handover
------------
The dense phase ends when the flat core of the concentration profile is
consumed by its own Gaussian shoulders, :math:`B = B_{eff} - (\\sqrt\\pi/2) S_y
\\le 0`.  ``RKGST`` will have stepped past that point, so ``PSSOUT`` linearly
interpolates the state back to :math:`B = 0` using the fraction

.. math:: f = \\frac{B_{prev}}{B_{prev} - B_{now}}

and interpolates *distance*, concentration, added heat, layer density, mass
flux and the two flammable-mass integrals along with it.  The Gaussian stage
then starts from those interpolated values.

Continuity is imposed on concentration, not on :math:`S_y`.  The Gaussian
stage needs a :math:`\\sigma_y` consistent with the concentration it inherits,
so ``DEG2S`` inverts the Gaussian concentration formula for the
:math:`\\sigma_y` that reproduces the final dense-phase ``Cc``, then solves the
Pasquill-Gifford power law for the virtual origin ``XV`` that would have
produced it:

.. math::
    S_{y,T} = \\frac{E (1+\\alpha)}{u_0 S_z \\sqrt\\pi C_c}
              \\left(\\frac{z_0}{S_z}\\right)^{\\alpha}, \\qquad
    X_V = \\left(\\frac{S_{y,T}}{\\sqrt 2 \\delta_y}\\right)^{1/\\beta_y} - X_T

Output recording
----------------
Both controllers record a point when a monitored quantity has moved by more
than a relative criterion (``ODLP``, ``ODLG``) *or* when a fixed distance has
passed (``ODLLP``, ``ODLLG``, both 80 m), and both record the *previous* step
rather than the current one -- except immediately after a recorded point, when
they record the current one.  ``PSSOUT`` skips density, gamma, temperature and
:math:`S_y` from its error test; ``SSGOUT`` skips only :math:`B`, which is
identically zero there.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import RT2, SQRTPI
from .downwind import (
    G_DH,
    G_MHI,
    G_MLOW,
    G_RHOUH,
    P_BEFF,
    P_DH,
    P_MHI,
    P_MLOW,
    P_RHOUH,
    P_SY2,
    Downwind,
)
from .rkgst import Control, rkgst

#: Columns of a recorded profile point, in the order ``SSOUT`` writes them.
PROFILE_COLUMNS = (
    "dist", "yc", "cc", "rho", "gamma", "temp", "b", "sz", "sy",
)


@dataclass
class Profile:
    """The downwind concentration profile."""

    rows: np.ndarray  #: (n, 9), see :data:`PROFILE_COLUMNS`
    #: Distance at which the dense phase gave way to the Gaussian phase, m.
    transition: float = 0.0
    xv: float = 0.0  #: virtual origin, m
    mass_above_ufl: float = 0.0  #: kg
    mass_above_lfl: float = 0.0  #: kg
    n_dense: int = 0  #: how many rows came from the dense phase

    @property
    def mass_between(self) -> float:
        """Contaminant mass between the two levels of concern, kg."""
        return self.mass_above_lfl - self.mass_above_ufl

    def column(self, name: str) -> np.ndarray:
        return self.rows[:, PROFILE_COLUMNS.index(name)]

    def distance_to(self, mole_fraction: float) -> float:
        """Downwind distance at which ``yc`` falls to ``mole_fraction``.

        Interpolated in log-concentration against distance, which is how the
        profile actually behaves; returns ``nan`` if the cloud never reaches
        that level.
        """
        yc = self.column("yc")
        if mole_fraction >= yc[0] or mole_fraction <= yc[-1]:
            return float("nan")
        i = int(np.argmax(yc <= mole_fraction))
        x0, x1 = self.rows[i - 1, 0], self.rows[i, 0]
        y0, y1 = yc[i - 1], yc[i]
        f = (math.log(mole_fraction) - math.log(y0)) / (math.log(y1) - math.log(y0))
        return float(x0 + f * (x1 - x0))


class _Recorder:
    """Shared bookkeeping for ``PSSOUT`` and ``SSGOUT``.

    Both keep three vectors: the current step, the previous step, and the last
    point actually written.  The error test looks at the current step against
    *both* of the others, so a quantity that drifts slowly still triggers a
    record once the drift accumulates.
    """

    def __init__(self, width: int, oodist: float):
        self.width = width
        self.oodist = oodist
        self.curnt = np.zeros(width)
        self.bksp = np.zeros(width)
        self.out = np.zeros(width)
        self.nrec = 0
        self.rows: list[np.ndarray] = []

    def emit(self, vec: np.ndarray, thermo, gas, ambient, table, alpha1) -> None:
        """Port of ``SSOUT``: append a profile row.

        ``SSOUT`` also computes the half-widths to the two levels of concern
        at the elevation ``GASZZC``, which the listing prints but which are
        not fed back into the model; they are recomputed on demand instead of
        stored, so the row stays the nine physical quantities.
        """
        row = np.array(vec[:9], dtype=float)
        row[0] = row[0] + self.oodist
        self.rows.append(row)
        self.nrec += 1


class SteadyStateRun:
    """The ``DEG2S`` main program."""

    def __init__(self, downwind: Downwind, handoff, er2, *, oodist: float = 0.0):
        self.dw = downwind
        self.h = handoff
        self.er2 = er2
        self.oodist = oodist
        self.erate = handoff["ess"]
        self.yclow = handoff["yclow"]

    # -- initial conditions -------------------------------------------------

    def _initial(self):
        """Port of the ``DEG2S`` setup block.

        The material balance over the secondary source has to close before the
        integration starts: the take-up rate, the concentration and the
        geometry come from three different places in ``DEG1`` and need not be
        mutually consistent to round-off. ``DEG2S`` fixes it by adjusting
        whichever of :math:`\\sigma_{y0}` or :math:`\\sigma_{z0}` keeps the
        result physical.
        """
        h, dw = self.h, self.dw
        a1 = dw.alpha1
        rl, bb, sz0 = h["outl"], h["outb"], h["outsz"]
        cc = h["outcc"]
        if rl <= 0.0 or bb <= 0.0 or cc <= 0.0:
            raise ValueError(
                f"the source model produced a degenerate secondary source "
                f"(length {rl:g} m, half-width {bb:g} m, concentration "
                f"{cc:g} kg/m3); there is nothing for the downwind model to "
                f"start from"
            )
        qstr0 = self.erate / 2.0 / rl / bb
        sy0er = self.er2["sy0er"]

        ratio1 = (
            dw.case.u0 * dw.case.z0 / a1 / dw.case.z0**a1 * cc / bb / qstr0 / rl
        )
        ratio = ratio1 * sz0**a1 * (bb + SQRTPI / 2.0 * sy0er)
        if ratio < 1.0:
            sy0er = (1.0 / (ratio1 * sz0**a1) - bb) * 2.0 / SQRTPI
        else:
            sz0 = (1.0 / ((bb + SQRTPI / 2.0 * sy0er) * ratio1)) ** (1.0 / a1)

        lay = dw.th.table.from_concentration(cc / self.dw.p.dellay)
        ueff0 = dw.consts.ueff0
        # six cloud states, plus whatever the closure integrates alongside
        y = np.zeros(6 + dw.n_closure)
        y[P_RHOUH] = lay.rho * ueff0 * (sz0 / dw.case.z0) ** a1
        y[P_SY2] = sy0er * sy0er
        y[P_BEFF] = bb + SQRTPI / 2.0 * sy0er
        return y, rl, bb, sz0, cc, lay.rho

    # -- the dense phase ----------------------------------------------------

    def _run_dense(self, y, rl, sz0, recorder):
        dw = self.dw
        er2 = self.er2
        dw.sz_seed = sz0
        state: dict = {}
        stop: dict = {}

        def fct(x, yy, dd, prmt):
            state["st"] = dw.pss(x, yy, dd, self.erate)

        def outp(x, yy, dd, ihlf, ndim, prmt):
            st = state["st"]
            # PSSOUT copies PRMT(21) into PRMT(22) at the top: the seed for
            # the material-balance loop advances only at *output* points, not
            # at every derivative evaluation, so every trial step inside a
            # Runge-Kutta stage starts from the same sigma_z.
            dw.sz_seed = st.sz
            cur = np.array([
                x, st.yc, st.cc, st.rho, st.gamma, st.temp, st.bb, st.sz,
                math.sqrt(yy[P_SY2]), yy[P_DH], st.rholay, yy[P_RHOUH],
                yy[P_MHI], yy[P_MLOW],
            ])
            # PSSOUT increments LASTP on every call and resets it only on the
            # very first, so at the point of decision it is zero exactly when
            # the previous call recorded a point -- which is what selects
            # CURNT over BKSP below.
            recorder.lastp += 1
            if recorder.nrec == 0:
                recorder.curnt = cur.copy()
                recorder.lastp = 0
            recorder.bksp = recorder.curnt.copy()
            recorder.curnt = cur

            bksp, curnt = recorder.bksp, recorder.curnt

            if st.bb <= 0.0:
                # the flat core is gone: interpolate the state back to B = 0
                if bksp[0] != recorder.out[0]:
                    recorder.emit(bksp, None, None, None, None, None)
                frac = bksp[6] / (bksp[6] - curnt[6])
                stop["dist"] = bksp[0] - frac * (bksp[0] - curnt[0])
                stop["cc"] = bksp[2] - frac * (bksp[2] - curnt[2])
                yy[P_DH] = bksp[9] - frac * (bksp[9] - curnt[9])
                stop["rholay"] = bksp[10] - frac * (bksp[10] - curnt[10])
                yy[P_RHOUH] = bksp[11] - frac * (bksp[11] - curnt[11])
                yy[P_MHI] = bksp[12] - frac * (bksp[12] - curnt[12])
                yy[P_MLOW] = bksp[13] - frac * (bksp[13] - curnt[13])
                stop["reason"] = "B<=0"
                prmt.halt()
                return

            erm = 0.0
            if st.yc <= self.yclow:
                if recorder.nrec >= 5:
                    stop["dist"] = x
                    stop["again"] = True
                    if curnt[0] != recorder.out[0]:
                        recorder.emit(curnt, None, None, None, None, None)
                    prmt.halt()
                    return
                erm = er2["odlp"]  # force a record

            # error test: distance, mole fraction, concentration, then B and
            # sigma_z -- density, gamma, temperature and sigma_y are skipped
            for ii in (1, 2, 6, 7):
                denom = curnt[ii] + 1.0e-10
                erm = max(
                    erm,
                    abs((curnt[ii] - bksp[ii]) / denom),
                    abs((curnt[ii] - recorder.out[ii]) / denom),
                )

            dx = curnt[0] - recorder.out[0]
            if recorder.nrec != 0 and erm < er2["odlp"] and dx <= er2["odllp"]:
                return

            vec = curnt.copy() if recorder.lastp == 0 else bksp.copy()
            recorder.out = vec
            recorder.emit(vec, None, None, None, None, None)
            recorder.lastp = -1

        prmt = Control([rl / 2.0, 6.023e13, er2["stpp"], er2["errp"], er2["smxp"]] + [0.0] * 20)
        weights = [er2["wtszp"], er2["wtsyp"], er2["wtbep"], er2["wtdh"], 1.0, 1.0]
        # The closure's own states are integrated alongside but given zero
        # error weight, so they take no part in RKGST's step-size test. That
        # is what keeps the step sequence -- and therefore which points get
        # recorded, and therefore the output rows -- identical to DEGADIS's
        # whatever closure is fitted. Giving them unit weight instead lets a
        # state that starts at zero dominate the relative error test and
        # bisect the integration to a standstill.
        weights = weights + [0.0] * dw.n_closure
        recorder.lastp = 0
        # RKGST copies its initial vector, and PSSOUT writes the interpolated
        # B = 0 state back into that copy, so the caller must take the result.
        res = rkgst(fct, outp, prmt, y, weights, ndim=6 + dw.n_closure)
        return stop, res.y

    # -- the Gaussian phase -------------------------------------------------

    def _run_gaussian(self, y, stop, recorder):
        dw = self.dw
        er2 = self.er2
        a1 = dw.alpha1

        cc = stop["cc"]
        rholay = stop["rholay"]
        xt = stop["dist"]
        sz = (y[P_RHOUH] / rholay / dw.consts.ueff0) ** (1.0 / a1) * dw.case.z0
        syt = (
            self.erate * a1 * (dw.case.z0 / sz) ** dw.alpha
            / dw.case.u0 / sz / cc / SQRTPI
        )
        xv = (syt / RT2 / dw.deltay) ** (1.0 / dw.betay) - xt

        g = np.zeros(6)
        g[G_RHOUH] = y[P_RHOUH]
        g[G_DH] = y[P_DH]
        g[G_MHI] = y[P_MHI]
        g[G_MLOW] = y[P_MLOW]

        dw.sz_seed = sz
        state: dict = {}
        recorder.ri = 0
        recorder.rii = -100.0 / er2["stpg"]
        first = {"done": False}

        def fct(x, yy, dd, prmt):
            state["st"] = dw.ssg(x, yy, dd, self.erate, xv)

        def outp(x, yy, dd, ihlf, ndim, prmt):
            st = state["st"]
            dw.sz_seed = st.sz
            sy = RT2 * dw.deltay * (x + xv) ** dw.betay
            cur = np.array([x, st.yc, st.cc, st.rho, st.gamma, st.temp, 0.0, st.sz, sy])
            if not first["done"]:
                recorder.curnt = cur.copy()
                recorder.out = np.zeros(9)
                first["done"] = True
            recorder.ri += 1
            recorder.bksp = recorder.curnt.copy()
            recorder.curnt = cur

            if st.yc < self.yclow:
                recorder.emit(cur, None, None, None, None, None)
                prmt.halt()
                return

            erm = 0.0
            for ii in range(1, 9):
                if ii == 6:  # B is identically zero here
                    continue
                denom = recorder.curnt[ii] + 1.0e-10
                erm = max(
                    erm,
                    abs((recorder.curnt[ii] - recorder.bksp[ii]) / denom),
                    abs((recorder.curnt[ii] - recorder.out[ii]) / denom),
                )
            dx = recorder.curnt[0] - recorder.out[0]
            if recorder.ri != 1 and erm < er2["odlg"] and dx <= er2["odllg"]:
                return

            vec = (
                recorder.curnt.copy()
                if recorder.ri == recorder.rii + 1
                else recorder.bksp.copy()
            )
            recorder.out = vec
            recorder.emit(vec, None, None, None, None, None)
            recorder.ri = recorder.rii

        prmt = Control([xt, 6.023e23, er2["stpg"], er2["errg"], er2["smxg"]] + [0.0] * 20)
        weights = [er2["wtruh"], er2["wtdhg"], 1.0, 1.0]
        res = rkgst(fct, outp, prmt, g[:4], weights, ndim=4)
        return res.y, xt, xv

    # -- the driver ---------------------------------------------------------

    def run(self) -> Profile:
        """Integrate downwind until the concentration drops below ``YCLOW``."""
        y, rl, bb, sz0, cc, _rholay = self._initial()
        recorder = _Recorder(9, self.oodist)
        stop, y = self._run_dense(y, rl, sz0, recorder)
        n_dense = recorder.nrec

        if stop.get("again") or "cc" not in stop:
            # Either the cloud thinned below the level of concern while still
            # dense, or the dense phase ended without its flat core being
            # consumed -- a closure that lifts the cloud off the ground can
            # reach the concentration floor first. Both mean there is no
            # Gaussian stage to hand over to.
            return Profile(
                rows=np.array(recorder.rows), transition=stop.get("dist", 0.0),
                mass_above_ufl=float(y[P_MHI]), mass_above_lfl=float(y[P_MLOW]),
                n_dense=n_dense,
            )

        g, xt, xv = self._run_gaussian(y, stop, recorder)
        return Profile(
            rows=np.array(recorder.rows), transition=xt, xv=xv,
            mass_above_ufl=float(g[G_MHI]), mass_above_lfl=float(g[G_MLOW]),
            n_dense=n_dense,
        )
