"""The transient supervisor: a port of ``SSSUP.FOR`` and the ``DEG2`` driver.

``DEG2S`` runs one plume.  ``SSSUP`` runs ``NOBS`` of them -- one per observer
-- and the concentration at a fixed point and time is later assembled from
whichever observers pass it.  For each observer the sequence is:

1. Find the times it crosses the upwind and downwind edges of the secondary
   source (:mod:`.observer`).
2. Integrate ``OB`` between them, accumulating the material the observer
   sweeps up.
3. Convert what it collected into a starting plume: half-width ``B``, source
   strength, layer composition and :math:`\\sigma_{z0}`.
4. Run ``PSS`` then ``SSG`` from the downwind edge, exactly as ``DEG2S`` does,
   until the concentration drops below ``YCLOW``.

Each recorded point carries a *time* as well as a distance -- the time that
observer reaches that distance, ``TS(t_0, x) = t_0 + (x + R_{max})^{1/(1+\\alpha)}
/ \\aleph``.  That is what makes a set of steady profiles into a transient
field: ``DEG3`` and ``DEG4`` later sort the points by distance and interpolate
in time across observers.

Differences from the steady-state path
--------------------------------------
* ``PSS`` integrates four states, not six: the flammable-mass integrals are
  only accumulated for steady-state runs, and ``SSG`` integrates two.
* Each observer rebuilds the adiabatic mixing table on *its own* layer
  composition, because different observers cross the source at different
  times and therefore collect different mixtures.
* The last table entry is replaced by an extrapolation to the observer's
  centreline concentration, since that concentration generally lies past the
  end of the mixing line the observer actually sampled.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import RT2, SQRTPI
from .downwind import (
    G_DH,
    G_RHOUH,
    P_BEFF,
    P_DH,
    P_RHOUH,
    P_SY2,
    Downwind,
)
from .observer import (
    O_CRATE,
    O_HWIDTH,
    Observer,
    ObserverKinematics,
    crossing_downwind,
    crossing_upwind,
    release_times,
)
from .rkgst import Control, rkgst

#: Columns of a transient profile point, matching the record ``PSSOUT`` and
#: ``SSGOUT`` write for ``DEG3``.
TRANSIENT_COLUMNS = (
    "dist", "time", "yc", "cc", "rho", "gamma", "temp", "sz", "sy", "b",
)


@dataclass
class ObserverResult:
    """One observer's contribution to the transient field."""

    index: int
    t0: float  #: release time, s
    tup: float  #: crosses the upwind source edge, s
    tdown: float  #: crosses the downwind source edge, s
    xup: float
    xdown: float
    erate: float  #: equivalent source strength, kg/s
    b: float  #: initial half-width, m
    sz0: float
    cc0: float  #: initial centreline concentration, kg/m**3
    rows: np.ndarray = field(default_factory=lambda: np.zeros((0, 10)))
    transition: float = 0.0  #: distance at which the dense phase ended, m
    xv: float = 0.0
    n_dense: int = 0


@dataclass
class TransientResult:
    """The whole transient run."""

    observers: list[ObserverResult]

    @property
    def rows(self) -> np.ndarray:
        """All recorded points, from every observer, stacked."""
        parts = [o.rows for o in self.observers if len(o.rows)]
        return np.vstack(parts) if parts else np.zeros((0, 10))

    def max_concentration_at(self, dist: float) -> tuple[float, float]:
        """Peak mole fraction at a downwind distance, and when it occurs.

        Interpolates each observer's profile to ``dist`` and takes the largest.
        This is the quantity a siting study reports; a full ``DEG3`` sort would
        also give the time history at that distance.
        """
        best_yc = 0.0
        best_t = float("nan")
        for o in self.observers:
            if len(o.rows) < 2:
                continue
            x = o.rows[:, 0]
            if dist < x.min() or dist > x.max():
                continue
            yc = float(np.interp(dist, x, o.rows[:, 2]))
            if yc > best_yc:
                best_yc = yc
                best_t = float(np.interp(dist, x, o.rows[:, 1]))
        return best_yc, best_t


class TransientRun:
    """The ``SSSUP`` supervisor."""

    def __init__(
        self, downwind: Downwind, handoff, er2, vectors, *, nobs: int | None = None,
        oodist: float = 0.0,
    ):
        self.dw = downwind
        self.h = handoff
        self.er2 = er2
        self.v = vectors
        self.nobs = int(nobs if nobs is not None else er2["nobs"])
        self.oodist = oodist
        self.yclow = handoff["yclow"]
        # RHOE lives in /parm/ and is read once from the handoff. It is *not*
        # the current mixing table's last row: SSSUP rebuilds that table for
        # every observer, so reading it from the table would let the clamp
        # drift observer by observer.
        self.rhoe = handoff["rhoe"]
        self.hmrte = 0.0  #: set by the caller from the source enthalpy
        self.kin = ObserverKinematics(
            alpha=handoff["alpha"], aleph=handoff["aleph"], rmax=handoff["rmax"]
        )

    # -- one observer -------------------------------------------------------

    def _crossings(self, t0: float) -> tuple[float, float, bool, bool]:
        """Port of the ``pup``/``pdn`` guards in ``SSSUP``.

        An observer released too late to reach the downwind edge before the
        source ends, or early enough that it has already passed the upwind
        edge when the spill begins, gets the endpoint rather than a root.
        Solving for a root that does not exist makes ``LIMIT`` abort.
        """
        kin, h = self.kin, self.h
        xend = self.v.radg_at(h["tend"])
        pdn = not (xend > kin.position(h["tend"], t0))
        r0 = self.v.radg_at(0.0)
        pup = not (t0 <= 0.0 and kin.position(0.0, t0) > -r0)
        tup = (
            crossing_upwind(t0, kin, self.v.radg_at, tol=self.er2["ertupf"])
            if pup else 0.0
        )
        tdown = (
            crossing_downwind(t0, kin, self.v.radg_at, tol=self.er2["ertdnf"])
            if pdn else h["tend"]
        )
        return tup, tdown, pup, pdn

    def _collect(self, t0: float, tup: float, tdown: float, xup: float, xdown: float):
        """Integrate ``OB`` across the source and return the layer it collected."""
        er2 = self.er2
        ob = Observer(
            self.dw.th, self.v, self.kin, params=self.dw.p,
            gammaf=self.dw.gammaf, ustar=self.dw.ustar, rhoa=self.dw.rhoa,
            u0=self.dw.case.u0, z0=self.dw.case.z0,
            humid=self.dw.th.ambient.humid,
        )
        y, rlen = ob.initial_conditions(tup, xdown, xup, er2["sz0er"])
        adiabatic = (
            self.dw.th.ambient.isofl == 1 or self.dw.th.ambient.ihtfl == 0
        )
        ndim = 4 if adiabatic else 5

        def fct(t, yy, dd, prmt):
            fct.state = ob.derivatives(t, yy, dd, t0, rlen)

        def outp(t, yy, dd, ihlf, ndim_, prmt):
            ob.record(fct.state)

        prmt = Control([
            tup, tdown, er2["stpo"], er2["erro"], max(1.0, (tdown - tup) / 50.0)
        ] + [0.0] * 20)
        weights = [er2["wtaio"], er2["wtqoo"], er2["wtszo"], 1.0, 1.0]
        res = rkgst(fct, outp, prmt, y, weights, ndim=ndim)
        return ob.last, res.y, rlen

    def _seed_plume(self, layer, y, rlen, tdown):
        """Port of the ``SSSUP`` block that turns a collected slab into a plume."""
        dw = self.dw
        a1 = dw.alpha1
        b = float(y[O_HWIDTH])
        area = b * rlen
        qstr0 = float(y[O_CRATE]) / area
        erate = 2.0 * qstr0 * rlen * b

        cc = min(layer.cclay * dw.p.dellay, self.rhoe)
        source = dw.case.source
        ccp = source.wc_at(tdown) * source.rho_at(tdown)
        cc = min(cc, ccp)

        sz0 = (qstr0 * rlen / cc * a1 / dw.case.u0 / dw.case.z0) ** (1.0 / a1) * dw.case.z0
        sy0er = 0.0
        ratio1 = (
            dw.case.u0 * dw.case.z0 / a1 / dw.case.z0**a1 * cc / b / qstr0 / rlen
        )
        ratio = ratio1 * sz0**a1 * (b + SQRTPI / 2.0 * sy0er)
        if ratio < 1.0:
            sy0er = (1.0 / (ratio1 * sz0**a1) - b) * 2.0 / SQRTPI
        else:
            sz0 = (1.0 / ((b + SQRTPI / 2.0 * sy0er) * ratio1)) ** (1.0 / a1)

        # each observer sampled its own mixture, so the mixing line is rebuilt
        dw.th.ambient.humsrc = (
            1.0 - layer.wclay - layer.walay * (1.0 + dw.th.ambient.humid)
        ) / layer.wclay
        dw.th.build_adiabatic_table(layer.wclay, layer.walay, layer.enthlay)
        if dw.th.ambient.isofl != 1:
            extend_table_to_centreline(
                dw.th.table, dw.case.gas, dw.th.ambient,
                cc=cc, cclay=layer.cclay, rhoa=dw.rhoa,
                hmrte=self.hmrte, gastem=dw.case.gas.temp,
            )
        return erate, b, sz0, sy0er, cc, layer.rholay

    def _profile(self, obs: ObserverResult, erate, b, sz0, sy0er, cc, rholay):
        """Run ``PSS`` then ``SSG`` for one observer, recording (x, t) points."""
        dw, er2, kin = self.dw, self.er2, self.kin
        a1 = dw.alpha1
        rows: list[np.ndarray] = []

        y = np.zeros(4)
        y[P_RHOUH] = rholay * dw.consts.ueff0 * (sz0 / dw.case.z0) ** a1
        y[P_SY2] = sy0er * sy0er
        y[P_BEFF] = b + SQRTPI / 2.0 * sy0er
        dw.sz_seed = sz0

        rec = {"n": 0, "lastp": 0, "curnt": None, "bksp": None, "out": np.zeros(12)}
        stop: dict = {}
        st_holder: dict = {}

        def ts(x: float) -> float:
            return obs.t0 + (x + kin.rmax) ** (1.0 / a1) / kin.aleph

        def emit(vec):
            rows.append(
                np.array([
                    vec[0] + self.oodist, ts(vec[0]), vec[1], vec[2], vec[3],
                    vec[4], vec[5], vec[6], vec[7], vec[8],
                ])
            )
            rec["n"] += 1

        def fct(x, yy, dd, prmt):
            st_holder["st"] = dw.pss(x, yy, dd, erate)

        def outp(x, yy, dd, ihlf, ndim, prmt):
            st = st_holder["st"]
            dw.sz_seed = st.sz
            rec["lastp"] += 1
            cur = np.array([
                x, st.yc, st.cc, st.rho, st.gamma, st.temp, st.sz,
                math.sqrt(yy[P_SY2]), st.bb, yy[P_DH], st.rholay, yy[P_RHOUH],
            ])
            if rec["n"] == 0:
                rec["curnt"] = cur.copy()
                rec["lastp"] = 0
            rec["bksp"] = rec["curnt"].copy()
            rec["curnt"] = cur
            bksp, curnt = rec["bksp"], rec["curnt"]

            if st.bb <= 0.0:
                frac = bksp[8] / (bksp[8] - curnt[8])
                stop["dist"] = bksp[0] - frac * (bksp[0] - curnt[0])
                stop["cc"] = bksp[2] - frac * (bksp[2] - curnt[2])
                yy[P_DH] = bksp[9] - frac * (bksp[9] - curnt[9])
                stop["rholay"] = bksp[10] - frac * (bksp[10] - curnt[10])
                yy[P_RHOUH] = bksp[11] - frac * (bksp[11] - curnt[11])
                prmt.halt()
                return

            erm = 0.0
            if st.yc <= self.yclow:
                stop["again"] = True
                if curnt[0] != rec["out"][0]:
                    emit(curnt)
                prmt.halt()
                return

            for ii in (1, 2, 6, 8):
                denom = curnt[ii] + 1.0e-10
                erm = max(
                    erm,
                    abs((curnt[ii] - bksp[ii]) / denom),
                    abs((curnt[ii] - rec["out"][ii]) / denom),
                )
            dx = curnt[0] - rec["out"][0]
            if rec["n"] != 0 and erm < er2["odlp"] and dx <= er2["odllp"]:
                return
            vec = curnt.copy() if rec["lastp"] == 0 else bksp.copy()
            rec["out"] = vec
            emit(vec)
            rec["lastp"] = -1

        prmt = Control([
            obs.xdown, 6.023e23, er2["stpp"], er2["errp"], er2["smxp"]
        ] + [0.0] * 20)
        res = rkgst(
            fct, outp, prmt, y,
            [er2["wtszp"], er2["wtsyp"], er2["wtbep"], er2["wtdh"]], ndim=4,
        )
        y = res.y
        obs.n_dense = rec["n"]

        if stop.get("again") or "dist" not in stop:
            obs.rows = np.array(rows) if rows else np.zeros((0, 10))
            return

        # -- the Gaussian stage ---------------------------------------------
        xt = stop["dist"]
        cc_t = stop["cc"]
        rholay_t = stop["rholay"]
        sz = (y[P_RHOUH] / rholay_t / dw.consts.ueff0) ** (1.0 / a1) * dw.case.z0
        syt = (
            erate * a1 * (dw.case.z0 / sz) ** dw.alpha
            / dw.case.u0 / sz / cc_t / SQRTPI
        )
        xv = (syt / RT2 / dw.deltay) ** (1.0 / dw.betay) - xt
        obs.transition, obs.xv = xt, xv

        g = np.array([y[P_RHOUH], y[P_DH]])
        dw.sz_seed = sz
        rec2 = {"ri": 0, "rii": -100.0 / er2["stpg"], "curnt": None,
                "bksp": None, "out": np.zeros(12), "first": False}

        def fct2(x, yy, dd, prmt):
            st_holder["st"] = dw.ssg(x, yy, dd, erate, xv)

        def outp2(x, yy, dd, ihlf, ndim, prmt):
            st = st_holder["st"]
            dw.sz_seed = st.sz
            sy = RT2 * dw.deltay * (x + xv) ** dw.betay
            cur = np.array([
                x, st.yc, st.cc, st.rho, st.gamma, st.temp, st.sz, sy, 0.0,
                yy[G_DH], st.rholay, yy[G_RHOUH],
            ])
            if not rec2["first"]:
                rec2["curnt"] = cur.copy()
                rec2["first"] = True
            rec2["ri"] += 1
            rec2["bksp"] = rec2["curnt"].copy()
            rec2["curnt"] = cur

            if st.yc < self.yclow:
                emit(cur)
                prmt.halt()
                return

            erm = 0.0
            for ii in range(1, 9):
                if ii == 8:  # B is identically zero here
                    continue
                denom = rec2["curnt"][ii] + 1.0e-10
                erm = max(
                    erm,
                    abs((rec2["curnt"][ii] - rec2["bksp"][ii]) / denom),
                    abs((rec2["curnt"][ii] - rec2["out"][ii]) / denom),
                )
            dx = rec2["curnt"][0] - rec2["out"][0]
            if rec2["ri"] != 1 and erm < er2["odlg"] and dx <= er2["odllg"]:
                return
            vec = (
                rec2["curnt"].copy()
                if rec2["ri"] == rec2["rii"] + 1
                else rec2["bksp"].copy()
            )
            rec2["out"] = vec
            emit(vec)
            rec2["ri"] = rec2["rii"]

        prmt = Control([xt, 6.023e23, er2["stpg"], er2["errg"], er2["smxg"]] + [0.0] * 20)
        rkgst(fct2, outp2, prmt, g, [er2["wtruh"], er2["wtdhg"]], ndim=2)
        obs.rows = np.array(rows) if rows else np.zeros((0, 10))

    # -- the supervisor -----------------------------------------------------

    def run(self) -> TransientResult:
        """Run every observer.  Failures are recorded, not raised.

        A single observer can fail to converge -- typically one released at
        the very edge of the window, whose slab is degenerate -- without
        invalidating the rest, and the original carries on in the same way.
        """
        t0s = release_times(self.kin, self.v, self.h["tend"], self.nobs)
        results: list[ObserverResult] = []
        for i, t0 in enumerate(t0s):
            tup, tdown, _pup, _pdn = self._crossings(t0)
            xup = self.kin.position(tup, t0)
            xdown = self.kin.position(tdown, t0)
            obs = ObserverResult(
                index=i, t0=t0, tup=tup, tdown=tdown, xup=xup, xdown=xdown,
                erate=0.0, b=0.0, sz0=0.0, cc0=0.0,
            )
            try:
                layer, y, rlen = self._collect(t0, tup, tdown, xup, xdown)
                if layer is None:
                    results.append(obs)
                    continue
                erate, b, sz0, sy0er, cc, rholay = self._seed_plume(
                    layer, y, rlen, tdown
                )
                obs.erate, obs.b, obs.sz0, obs.cc0 = erate, b, sz0, cc
                self._profile(obs, erate, b, sz0, sy0er, cc, rholay)
            except (RuntimeError, ValueError, OverflowError):
                pass
            results.append(obs)
        return TransientResult(observers=results)


def extend_table_to_centreline(
    table, gas, ambient, *, cc: float, cclay: float, rhoa: float,
    hmrte: float, gastem: float,
) -> None:
    """Port of the ``SSSUP`` block at labels 300-400.

    Each observer collects a mixture and ``SETDEN`` builds the adiabatic
    mixing line *for that mixture*, which by construction ends at the layer
    composition. But the quantity the downwind model starts from is the
    *centreline* concentration, ``dellay`` times the layer's, which lies past
    the end of that line. So the last table entry is replaced by an
    extrapolation to it.

    The density extrapolates linearly in concentration -- the mixing line is
    close to straight there -- and the enthalpy and temperature extrapolate
    against ``yc/rho``, on the reasoning that ratio is proportional to
    temperature. Both are then clamped: enthalpy between the release value and
    zero, temperature between the release temperature and ambient.

    If the extrapolation implies a mass fraction above one, the mixture is
    already pure contaminant and the entry becomes exactly that, with the
    temperature depending on which thermal model is active. Modifies ``table``
    in place.
    """
    import numpy as np

    n = table.n
    i = n - 1  # the entry to overwrite: the old pure-contaminant end point
    rholay_prev = float(table.rho[i - 1])
    rho = cc * (rholay_prev - rhoa) / cclay + rhoa
    wc = cc / rho

    if wc > 1.0:
        # already pure contaminant; the density *is* the concentration
        if ambient.isofl == 1:
            h_new, t_new = float(table.h[i - 1]), float(table.t[i - 1])
        elif ambient.ihtfl == 0:
            h_new, t_new = hmrte, gastem
        else:
            hhh = float(table.h[i - 1])
            if abs(hhh) > abs(hmrte):
                hhh = hmrte
            h_new = hhh
            t_new = min(gastem * gas.rho / cc, ambient.tamb)
        table.yc[i], table.cc[i], table.rho[i] = 1.0, cc, cc
        table.h[i], table.t[i] = h_new, t_new
        return

    from .constants import WMA, WMW

    wa = (1.0 - (1.0 + ambient.humsrc) * wc) / (1.0 + ambient.humid)
    wm = 1.0 / (wc / gas.mw + wa / WMA + (1.0 - wa - wc) / WMW)
    yc = min(max(wm / gas.mw * wc, 0.0), 1.0)

    denom = (
        table.yc[i - 2] / table.rho[i - 2] - table.yc[i - 1] / table.rho[i - 1]
    )
    if denom != 0.0:
        slope = (yc / rho - table.yc[i - 1] / table.rho[i - 1]) / denom
        h_new = slope * (table.h[i - 2] - table.h[i - 1]) + table.h[i - 1]
        h_new = min(max(hmrte, h_new), 0.0)
        t_new = slope * (table.t[i - 2] - table.t[i - 1]) + table.t[i - 1]
        t_new = min(max(gastem, t_new), ambient.tamb)
    else:
        h_new, t_new = float(table.h[i - 1]), float(table.t[i - 1])

    table.yc[i], table.cc[i], table.rho[i] = yc, cc, rho
    table.h[i], table.t[i] = h_new, t_new
