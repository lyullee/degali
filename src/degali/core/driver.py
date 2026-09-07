"""The ``DEG1`` driver: run the source calculation to completion.

Ports the ``DEG1`` main program, ``SRC1O`` and ``NOBL.FOR``.  Where
:mod:`.blanket` holds the physics, this module holds the control flow, which
is where most of DEGADIS's remaining subtlety lives.

The decision the driver makes
-----------------------------
At every restart point the driver compares the primary source flux against the
atmospheric take-up ceiling over the same footprint:

* flux below the ceiling  ->  no blanket forms; :func:`no_blanket` walks the
  rest of the release on a fixed time grid.
* flux above the ceiling  ->  integrate the blanket until it either drains
  (``H <= SRCCUT``), reaches steady state, or the release ends.

These two can alternate.  A blanket that drains hands control to
:func:`no_blanket`, which watches the flux ratio and hands it straight back if
the source strengthens again.  The Fortran guards against ping-ponging with a
lookahead: if the blanket collapses before ``SRCSS`` seconds have elapsed, it
scans the *remaining* source table, and if the flux never again exceeds the
current take-up it latches ``REFLAG`` false so ``NOBL`` can no longer restart
the blanket.  That lookahead is reproduced in :meth:`SourceRun._may_restart`.

``SRC1O`` is not an output routine
----------------------------------
It is a stateful controller called from inside ``RKGST``, and it does four
things:

1. Records an output point when any of 13 monitored quantities has moved by
   more than ``SRCOER`` relative, measured against *both* the previous step
   and the last recorded point.
2. Declares steady state -- either immediately when the front has stopped
   moving, or after ``SRCSS`` seconds during which a reduced set of
   quantities stayed inside ``SRCOER``.  The reduced set excludes cloud
   depth, water fractions, Richardson number and enthalpy, which keep
   drifting long after the concentration field has settled.
3. Stops the integration when the blanket drains, and records the terminal
   point.
4. Tracks the peak take-up rate and the largest radius, which set the
   wind-speed constant ``ALEPH`` handed to the downwind model.

Its counters and buffers persist between calls (Fortran ``DATA`` plus
``/noauto``), so it is a class here rather than a function.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .atmosphere import friction_velocity
from .blanket import (
    I_E,
    I_M,
    I_MA,
    I_MC,
    I_P,
    I_R,
    Blanket,
    BlanketParameters,
    BlanketResult,
    BlanketState,
)
from .constants import PI, SQRTPI, VKC
from .entrainment import phi_hat
from .numerics import gamma
from .rkgst import Control, rkgst

#: The 13 quantities ``SRC1O`` monitors, in the order it writes them.
MONITORED = (
    "time", "radius", "height", "qstar", "sz",
    "yc", "ya", "rho", "ri", "wc", "wa", "enth", "temp",
)

#: Indices (0-based into :data:`MONITORED`) excluded from the steady-state
#: test.  The Fortran writes ``ii.ne.3 .and. ii.ne.9 .and. ii.ne.12 .and.
#: ii.ne.7 .and. ii.ne.11`` over 1-based indices: height, Richardson number,
#: enthalpy, air mole fraction and air mass fraction.  These drift long after
#: the concentration field has settled, so including them would prevent the
#: run from ever converging.
_SS_EXCLUDED = {2, 8, 11, 6, 10}


@dataclass
class DriverParameters(BlanketParameters):
    """Adds the ``DEG1``-level controls to the blanket closure constants."""

    stpin: float = 0.01  #: initial time step, s
    erbnd: float = 0.005  #: RKGST relative error bound
    #: Per-state error weights, ``DERY`` on entry to ``RKGST``.
    wtrg: float = 1.0
    wttm: float = 1.0
    wtyc: float = 1.0
    wtya: float = 1.0
    wteb: float = 1.0
    wtmb: float = 1.0
    ernobl: float = 1.0005  #: NOBL blanket-restart threshold on flux/qstar
    noblpt: int = 100  #: number of NOBL time steps
    stpmax: float = 6.0e3  #: maximum time step, s
    nobl_dtmin: float = 0.5  #: minimum NOBL time increment, s


class _OutputController:
    """Port of ``SRC1O``.  Stateful; one instance per ``RKGST`` call."""

    def __init__(self, run: "SourceRun"):
        self.run = run
        self.p = run.p
        self.i = 0
        self.iii = 0
        self.tlast = 0.0
        self.curnt = np.zeros(13)
        self.bksp = np.zeros(13)
        self.outp = np.zeros(13)
        self.records: list[np.ndarray] = []
        self.stop_reason: str | None = None

    def __call__(self, time, y, dery, ihlf, ndim, prmt) -> None:
        run, p = self.run, self.p
        st = run.last_state
        if st is None:
            return
        self.i += 1
        self.iii += 1

        run.ctx["vel"] = st.velocity
        run.ctx["ht_last"], run.ctx["hh_last"] = st.ht, st.hh

        hei = st.height
        if hei <= 0.0:
            return self._terminate(time, y, st, "drained", clamp=True)

        # peak take-up and largest radius: these set ALEPH downwind
        qsav = PI * y[I_R] * y[I_R] * st.qstar
        if qsav >= run.emax:
            run.emax = qsav
            run.rm = float(y[I_R])
            run.szm = st.sz
        run.rmax = max(run.rmax, float(y[I_R]))

        if hei <= p.srccut:
            return self._terminate(time, y, st, "srccut")
        if st.mixture.yc <= run.case.yclow and run.case.u0 == 0.0:
            return self._terminate(time, y, st, "yclow (no wind)")
        if time > run.tend + 1.0 and run.case.u0 == 0.0 and st.velocity == 0.0:
            return self._terminate(time, y, st, "stagnant (no wind)")

        cur = np.array(
            [time, y[I_R], hei, st.qstar, st.sz, st.mixture.yc, st.mixture.ya,
             st.mixture.rho, st.richardson, st.mixture.wc, st.mixture.wa,
             st.mixture.enthalpy, st.mixture.temp]
        )

        if self.i == 1:
            self.curnt = cur
            self.iii = 1
            return self._record(time)

        self.bksp = self.curnt.copy()
        self.curnt = cur

        erm = 0.0
        ermss = 0.0
        for ii in range(1, 13):
            div = self.curnt[ii] if self.curnt[ii] != 0.0 else p.srcoer
            er1 = abs((self.curnt[ii] - self.bksp[ii]) / div)
            er2 = abs((self.curnt[ii] - self.outp[ii]) / div)
            if ii not in _SS_EXCLUDED:
                ermss = max(er1, er2, ermss)
            erm = max(er1, er2, erm)

        if run.case.steady_state:
            if st.velocity == 0.0 and time > p.srcss:
                return self._steady(time, y, st)
            if ermss > p.srcoer:
                return self._record(time)
            if time - self.tlast > p.srcss:
                return self._steady(time, y, st)
            return None

        if erm < p.srcoer:
            return None
        return self._record(time)

    # -- outcomes ----------------------------------------------------------

    def _record(self, time: float) -> None:
        self.tlast = time
        if self.iii == 1:
            self.bksp = self.curnt.copy()
        self.outp = self.bksp.copy()
        self.iii = 0
        self.records.append(self.outp.copy())

    def _steady(self, time, y, st: BlanketState) -> None:
        """Label 122: steady state reached, close out the secondary source.

        Note the Fortran reuses the local ``Qstar`` here, recomputing it from
        ``PRMT(21)`` -- which holds ``ERTE``, the *primary* source rate, not
        the take-up -- and then falls through to the terminal write, so the
        last row of the listing carries the recomputed pair rather than the
        integrated one.
        """
        run = self.run
        run.result.outcc = st.mixture.wc * st.mixture.rho
        run.result.swcl = st.mixture.wc
        run.result.swal = st.mixture.wa
        run.result.srhl = st.mixture.rho
        run.result.senl = st.mixture.enthalpy
        run.result.outl = 2.0 * float(y[I_R])
        sz = st.sz
        if run.case.u0 != 0.0:
            alpha1 = run.alpha + 1.0
            qstar = run.last_erate / PI / y[I_R] ** 2
            sz = (
                alpha1 / run.case.u0 / run.case.z0
                * qstar * run.result.outl / run.result.outcc
            ) ** (1.0 / alpha1) * run.case.z0
        run.result.outsz = sz
        run.result.outb = PI * float(y[I_R]) ** 2 / run.result.outl / 2.0
        self._terminate(time, y, st, "steady state", qstar=qstar, sz=sz)

    def _terminate(
        self, time, y, st: BlanketState, why: str, clamp=False,
        qstar: float | None = None, sz: float | None = None,
    ) -> None:
        run = self.run
        self.i = -1
        if time >= run.tend:
            run.check3 = True
        if why in ("steady state",):
            run.check3 = True
        hei = st.height
        if clamp:
            hei = 0.0
            y[I_R] = min(run.rmax, y[I_R])
        run.tsc1 = time
        self.records.append(
            np.array([time, y[I_R], hei,
                      st.qstar if qstar is None else qstar,
                      st.sz if sz is None else sz, st.mixture.yc,
                      st.mixture.ya, st.mixture.rho, st.richardson,
                      st.mixture.wc, st.mixture.wa, st.mixture.enthalpy,
                      st.mixture.temp])
        )
        self.stop_reason = why
        run.prmt.halt()


class SourceRun:
    """The ``DEG1`` main program."""

    def __init__(self, case, thermo, params: DriverParameters | None = None):
        self.case = case
        self.th = thermo
        self.p = params or DriverParameters()
        self.ustar = friction_velocity(case.u0, case.z0, case.zr, case.rml)
        self.alpha = None
        self.gammaf = None
        self.rhoa = thermo.table.rhoa

        self.blanket = Blanket(case, thermo, self.p)
        self.blanket.ustar = self.ustar
        self.blanket.rhoa = self.rhoa

        self.result = BlanketResult()
        self.tend = case.tend
        self.tsc1 = case.tend
        self.emax = 0.0
        self.rm = 0.0
        self.szm = 0.0
        self.rmax = 0.0
        self.check3 = False
        self.reflag = True
        self.records: list[np.ndarray] = []
        self.last_state: BlanketState | None = None
        self.last_takeup = 0.0
        self.last_erate = 0.0  # PRMT(21): the primary source rate
        self.ctx: dict = {}
        self.prmt: list = []

    # -- helpers -----------------------------------------------------------

    def set_alpha(self, alpha: float) -> None:
        self.alpha = alpha
        self.gammaf = gamma(1.0 / (alpha + 1.0))
        self.blanket.alpha = alpha
        self.blanket.gammaf = self.gammaf

    def _takeup_ceiling(self, time: float) -> tuple[float, float]:
        """``(qstar, qstre)``: atmospheric ceiling and primary flux, kg/(m**2 s)."""
        src = self.case.source
        rl = 2.0 * src.radius_at(time)
        qstre = src.rate_at(time) / (PI * rl**2 / 4.0)
        rhop = src.rho_at(time)
        ccp = src.wc_at(time) * rhop
        qstar = 0.0
        if self.case.u0 != 0.0:
            qstar = (
                ccp * VKC * self.ustar * (self.alpha + 1.0)
                * self.p.dellay / (self.p.dellay - 1.0)
                / phi_hat(rhop, rl, rhoa=self.rhoa, alpha=self.alpha,
                          ustar=self.ustar, u0=self.case.u0, z0=self.case.z0,
                          gammaf=self.gammaf, dellay=self.p.dellay)
            )
        return qstar, qstre

    def _may_restart(self, time0: float, qstar: float, qstre: float) -> bool:
        """Port of the ``REFLAG`` lookahead (labels 190-300 of ``DEG1``).

        A blanket that collapses before ``SRCSS`` seconds is suspicious: the
        integration may be oscillating around the balance point.  Scan the
        remainder of the source table; if the flux never exceeds the current
        one again, disable blanket restarts for good.
        """
        if qstre <= qstar:
            return True
        if time0 >= self.p.srcss:
            return True
        src = self.case.source
        idx = np.searchsorted(src.time, time0, side="right")
        for ii in range(idx, len(src)):
            if src.rate[ii] == 0.0:
                self.reflag = False
                return True
            flux = src.rate[ii] / (PI * src.radius[ii] ** 2)
            if flux > qstre:
                return True
        return True

    # -- the blanket integration -------------------------------------------

    def _integrate_blanket(self, time0: float) -> _OutputController:
        p, case, src = self.p, self.case, self.case.source
        pwcp = src.wc_at(time0)
        rhop = src.rho_at(time0)

        y = np.zeros(6)
        y[I_R] = src.radius_at(time0)
        self.rmax = float(y[I_R])
        y[I_M] = max(case.gmass0 / pwcp, PI * y[I_R] ** 2 * 1.1 * p.srccut * rhop)
        htod = (case.gmass0 / pwcp) / rhop / 2.0 / PI / y[I_R] ** 3
        y[I_MC] = y[I_M] * pwcp
        pwap = (1.0 - pwcp) / (1.0 + case.ambient.humid)
        y[I_MA] = y[I_M] * pwap
        y[I_E] = y[I_M] * src.enthalpy_at(time0)
        y[I_P] = 0.0

        self.ctx = {
            "vuflag": htod > 0.1,
            "rmax": self.rmax,
            "vel": 0.0,
            "ht_last": htod * 2.0 * y[I_R],
            "hh_last": 0.0,
            "wc_last": pwcp,
        }
        self.last_state = None

        # PRMT(2) = 6.023e23 -- the Fortran leaves the end time effectively
        # unbounded and relies on SRC1O to stop the integration.
        self.prmt = Control(
            [time0, 6.023e23, p.stpin, p.erbnd, p.stpmax] + [0.0] * 20
        )
        weights = [p.wtrg, p.wttm, p.wtyc, p.wtya, p.wteb, p.wtmb]
        control = _OutputController(self)

        def fct(x, yy, dd, prmt):
            d, st = self.blanket.derivatives(x, yy, self.ctx)
            dd[:] = d
            self.last_state = st
            self.last_takeup = st.qstar * PI * yy[I_R] ** 2
            self.last_erate = st.rate

        rkgst(fct, control, self.prmt, y, weights, ndim=6)
        self.records.extend(control.records)
        self.result.formed = True
        return control

    # -- the blanket-free path ----------------------------------------------

    def no_blanket(self) -> float | None:
        """Port of ``NOBL``.

        Walks the release on a fixed grid of ``NOBLPT`` steps, floored at
        ``NOBL_DTMIN`` seconds.  At each step it checks whether the primary
        flux has climbed back above the take-up ceiling; if it has, and
        restarts are still permitted, it returns the time at which to restart
        the blanket.
        """
        p, case, src = self.p, self.case, self.case.source
        npts = p.noblpt
        deltat = (self.tend - self.tsc1) / float(npts)
        if deltat < p.nobl_dtmin:
            npts = int((self.tend - self.tsc1) / p.nobl_dtmin) + 1
            deltat = (self.tend - self.tsc1) / float(npts)
        t0 = self.tsc1
        alpha1 = self.alpha + 1.0

        for i in range(1, npts + 1):
            time = t0 + float(i) * deltat
            if i == npts:
                time = self.tend
            rl = 2.0 * src.radius_at(time)
            erate = src.rate_at(time)
            flux = erate / (PI * rl * rl / 4.0)
            pwcp = src.wc_at(time)
            pwap = (1.0 - pwcp) / (1.0 + case.ambient.humid)
            hprim = src.enthalpy_at(time)
            self.th.build_adiabatic_table(pwcp, pwap, hprim)
            rhop = src.rho_at(time)
            ccp = pwcp * rhop

            qstar = (
                ccp * VKC * self.ustar * alpha1 * p.dellay / (p.dellay - 1.0)
                / phi_hat(rhop, rl, rhoa=self.rhoa, alpha=self.alpha,
                          ustar=self.ustar, u0=case.u0, z0=case.z0,
                          gammaf=self.gammaf, dellay=p.dellay)
            )
            if abs(flux / qstar) > p.ernobl and self.reflag:
                self.check3 = True
                return time

            # close the material balance on the SZF solution
            twidth = PI * rl / 4.0
            qqqq = erate / case.u0 / case.z0 * alpha1 * case.z0**alpha1 / twidth
            cc = min(qqqq / self._szf(flux, rl, pwcp) ** alpha1, ccp)
            sz = (qqqq / cc) ** (1.0 / alpha1)
            mix = self.th.table.from_concentration(cc)

            if erate >= self.emax:
                self.emax = erate
                self.rm = src.radius_at(time)
                self.szm = sz
            self.rmax = max(self.rmax, src.radius_at(time))

            self.records.append(
                np.array([time, src.radius_at(time), 0.0, flux, sz,
                          mix.yc, mix.ya, mix.rho, 0.0, mix.wc, mix.wa,
                          mix.enthalpy, mix.temp])
            )

            if i == 3 and case.steady_state:
                rlist = src.radius_at(time)
                self.result.outcc = cc
                self.result.swcl = mix.wc
                self.result.swal = mix.wa
                self.result.senl = mix.enthalpy
                self.result.srhl = mix.rho
                self.result.outsz = sz
                self.result.outl = 2.0 * rlist
                self.result.outb = PI * rlist**2 / self.result.outl / 2.0
                self.result.formed = False
                return None
        return None

    def _szf(self, flux: float, rl: float, pwcp: float) -> float:
        """Placeholder for ``SZF``; the material-balance closure that follows
        in ``NOBL`` overwrites this value, so only its magnitude matters here.
        """
        alpha1 = self.alpha + 1.0
        cc = pwcp * self.th.table.rhoa
        uheff = flux * rl / max(cc, 1e-30)
        return max(
            (uheff * alpha1 / self.case.u0 / self.case.z0) ** (1.0 / alpha1)
            * self.case.z0,
            1e-6,
        )

    # -- the driver ---------------------------------------------------------

    def run(self) -> BlanketResult:
        """Execute the source calculation, alternating blanket and NOBL."""
        case = self.case
        if self.alpha is None:
            raise RuntimeError("call set_alpha() first")

        time0 = 0.0
        qstar, qstre = self._takeup_ceiling(time0)
        blanket_first = not (qstre < qstar and case.gmass0 == 0.0)
        if not blanket_first:
            self.tsc1 = 0.0

        guard = 0
        while True:
            guard += 1
            if guard > 100:
                raise RuntimeError("DEG1 blanket/NOBL restart loop did not settle")

            if blanket_first:
                self.check3 = False
                self._integrate_blanket(time0)
                if self.check3:
                    self.tend = self.tsc1
                    break
                time0 = self.tsc1
                qstar, qstre = self._takeup_ceiling(time0)
                self._may_restart(time0, qstar, qstre)
                if qstre > qstar and time0 >= self.p.srcss:
                    continue

            timeout = self.no_blanket()
            if self.check3 and timeout is not None:
                time0 = timeout
                blanket_first = True
                continue
            break

        self.rmax *= 1.01  # "GUARANTEE A GOOD VALUE"
        aleph = 0.0
        if case.u0 != 0.0:
            alpha1 = self.alpha + 1.0
            aleph = (
                case.u0 / self.gammaf * (self.szm / case.z0) ** self.alpha
                / (SQRTPI / 2.0 * self.rm + self.rmax) ** (self.alpha / alpha1)
                / alpha1
            )

        r = self.result
        r.history = [
            BlanketState(
                time=rec[0], radius=rec[1], height=rec[2], qstar=rec[3],
                sz=rec[4], mixture=None, richardson=rec[8], velocity=0.0,
            )
            for rec in self.records
        ]
        r.rm, r.szm, r.emax, r.rmax = self.rm, self.szm, self.emax, self.rmax
        r.aleph = aleph
        r.tend = self.tend
        r.ess = case.source.rate_at(0.0)
        self.raw_records = np.array(self.records) if self.records else np.zeros((0, 13))
        return r
