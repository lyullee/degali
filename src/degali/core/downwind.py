"""Downwind dispersion: ports of ``PSS.FOR``, ``SSG.FOR`` and ``SERIES.FOR``.

The cloud leaves the secondary source as a slab of half-width :math:`B` with
Gaussian edges, and DEGADIS carries it downwind in two stages.

Dense phase (``PSS``)
---------------------
While the cloud is still heavier than air it spreads laterally under its own
weight.  The concentration profile is a flat core of half-width :math:`B_{eff}`
with Gaussian shoulders of scale :math:`S_y`, and the vertical profile is a
stretched exponential set by :math:`S_z` and the wind-profile exponent.  Four
quantities are integrated in distance:

===============  ==========================================================
:math:`\\rho u H` mass flux per unit width through the layer
:math:`S_y^2`     lateral Gaussian variance
:math:`B_{eff}`   half-width including the shoulders
:math:`\\Delta h`  heat picked up from the ground, per unit mass flux
===============  ==========================================================

with two more, the mass above the upper and lower flammable limits, carried
along for the flammable-cloud report.

Gravity spreading gives
:math:`dB_{eff}/dx \\propto \\sqrt{\\Delta\\rho/\\rho_a}\\,(S_z/z_0)^{1/2-\\alpha}`
and stops when the density excess drops below ``DELRMN``.  Entrainment feeds
:math:`\\rho u H` at a rate suppressed by :math:`\\Phi`, exactly as in the source
model.

The flat core shrinks as the shoulders grow, and when :math:`B = B_{eff} -
(\\sqrt{\\pi}/2) S_y` reaches zero the dense phase ends.

Gaussian phase (``SSG``)
------------------------
Past that point the cloud is a plain Gaussian plume: :math:`S_y` follows the
Pasquill-Gifford power law directly, so only :math:`\\rho u H` and
:math:`\\Delta h` are still integrated.  A virtual origin ``XV`` and virtual
time ``T0`` are chosen so the two stages join smoothly.

The inner iteration
-------------------
Both stages face the same implicit problem: the concentration depends on
:math:`S_z`, and :math:`S_z` depends on the layer density, which depends on the
concentration.  The Fortran resolves it by fixed-point iteration seeded from
the previous step's :math:`S_z` (``PRMT(22)``), converging to ``rcrit`` --
2e-3 in ``PSS``, 2.5e-3 in ``SSG``.  Those are loose enough that the
converged value depends on the seed, so the seed has to be carried exactly as
the original carries it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .atmosphere import richardson_star
from .constants import PI, SQPIO2, SQRTPI, RT2, VKC
from .entrainment import phi, surface_exchange
from .closures import CloudState, Closure, DegadisClosure
from .numerics import gaminc

#: State indices in the dense phase, from ``PSS``'s
#: ``DATA irhouh/1/,iSY2/2/,iBEFF/3/,idh/4/,Mhi/5/,Mlow/6/``.
P_RHOUH, P_SY2, P_BEFF, P_DH, P_MHI, P_MLOW = range(6)

#: State indices in the Gaussian phase, from ``SSG``'s
#: ``DATA irhouh/1/,idh/2/, mhi/3/, mlow/4/``.
G_RHOUH, G_DH, G_MHI, G_MLOW = range(4)


def series(arg: float, alpha1: float, *, kmax: int = 100, error: float = 1.0e-4) -> float:
    r"""Port of ``SERIES.FOR``.

    Evaluates the expansion that appears when integrating the Gaussian
    shoulders of the concentration profile above a threshold:

    .. math::
        S(w) = \sum_{k \ge 0} \frac{2(k+1)}{2k+3}
               \frac{w^{k+3/2}}{p (p+1) \cdots (p+k)} \cdot \frac{2}{3}

    with :math:`p = 1/(1+\alpha)`.  The original stops when a term falls below
    ``1e-4`` relative and warns that it overflows past ``arg = 13.8``; both
    behaviours are kept, the overflow guard as an explicit error.
    """
    if arg > 13.8:
        raise OverflowError(
            f"SERIES overflows for arguments above 13.8, got {arg:.4g}"
        )
    pp = 1.0 / alpha1
    coeff = 2.0 / 3.0
    term = coeff / pp * arg**1.5
    total = term
    for kk in range(1, kmax + 1):
        term = 2.0 * (kk + 1) / (2 * kk + 3) / (pp + kk) * arg * term
        total += term
        if abs(term) / max(abs(total), error) <= error:
            return total
    raise RuntimeError(f"SERIES did not converge for argument {arg:.5g}")


@dataclass
class DownwindConstants:
    """The group constants ``DEG2S`` precomputes into ``PRMT``."""

    #: ``PRMT(9)``: coefficient of the gravity-spreading rate.
    c_spread: float
    #: ``PRMT(10)``: coefficient of the vertical growth rate (unused by PSS
    #: directly, but part of the original's setup).
    c_sz: float
    #: ``PRMT(18) = u0 z0 / (1 + alpha)``.
    ueff0: float
    #: ``PRMT(19) = rho_a kappa u* (1 + alpha)``.
    entrain: float

    @classmethod
    def build(cls, *, ce, gg, z0, alpha1, gammaf, u0, alpha, vkc, ustar, rhoa):
        return cls(
            c_spread=ce * math.sqrt(gg * z0 / alpha1 * gammaf) * gammaf / u0,
            c_sz=z0**alpha * vkc * ustar * alpha1 * alpha1 / u0,
            ueff0=u0 * z0 / alpha1,
            entrain=rhoa * vkc * ustar * alpha1,
        )


@dataclass
class DownwindState:
    """Diagnostics returned alongside the derivatives."""

    dist: float
    cc: float  #: centreline ground-level concentration, kg/m**3
    bb: float  #: half-width of the flat core, m
    sz: float
    sy: float
    yc: float
    rho: float
    temp: float
    gamma: float  #: ``(rho - rhoa)/cc``, the density-concentration slope
    rholay: float
    heff: float


class Downwind:
    """The steady-state downwind model: ``PSS`` then ``SSG``."""

    def __init__(
        self, case, thermo, params, *, alpha, gammaf, ustar,
        deltay=None, betay=None, rhoa=None, consts=None, closure=None,
    ):
        self.case = case
        self.th = thermo
        self.p = params
        self.alpha = alpha
        self.alpha1 = alpha + 1.0
        self.gammaf = gammaf
        self.ustar = ustar
        self.rhoa = thermo.table.rhoa if rhoa is None else rhoa
        # DEG2 reads deltay/betay from the .TR2 handoff rather than
        # recomputing them from the stability class, and the file carries only
        # seven significant figures. Recomputing instead shifts sigma_y by
        # ~4e-7 relative, which the concentration inherits directly, so the
        # handoff values are used when supplied.
        self.deltay = case.stability.deltay if deltay is None else deltay
        self.betay = case.stability.betay if betay is None else betay
        self.consts = consts or DownwindConstants.build(
            ce=params.ce, gg=9.81, z0=case.z0, alpha1=self.alpha1,
            gammaf=gammaf, u0=case.u0, alpha=alpha, vkc=VKC,
            ustar=ustar, rhoa=self.rhoa,
        )
        #: ``PRMT(22)``: the previous step's sigma_z, seeding the inner loop.
        self.sz_seed = 0.0
        #: What the cloud does with its density difference. The default is
        #: DEGADIS's own choice -- slump while dense, do nothing while buoyant
        #: -- and is what every parity test exercises. Substituting another
        #: changes only the spreading and any vertical motion; the mass and
        #: energy balances are common to all of them.
        self.closure: Closure = closure or DegadisClosure(
            coefficient=self.consts.c_spread, alpha=alpha,
            z0=case.z0, delrmn=params.delrmn,
        )
        #: Extra states the closure integrates, appended after the cloud's own.
        self.n_closure = self.closure.n_states

    # -- shared pieces -----------------------------------------------------

    def _close_material_balance(self, cc_of_sz, y_rhouh, dh, rcrit, maxiter, trap):
        """The fixed-point loop shared by ``PSS`` and ``SSG``.

        ``cc_of_sz`` maps a trial ``sz`` to the centreline concentration; the
        layer density that follows then fixes ``sz`` again.
        """
        sz0 = self.sz_seed
        sz = sz0
        for _ in range(maxiter + 1):
            cc = cc_of_sz(sz)
            cclay = cc / self.p.dellay
            lay = self.th.add_heat(cclay, dh)
            prod = max(y_rhouh / lay.rho / self.consts.ueff0, 1.0e-10)
            sz = prod ** (1.0 / self.alpha1) * self.case.z0
            if abs(sz - sz0) / (abs(sz) + abs(sz0) + 1.0e-10) <= rcrit:
                return sz, cc, cclay, lay
            sz0 = sz
        raise RuntimeError(f"downwind material balance did not converge (TRAP {trap})")

    def _thermal(self, cc, cclay, dh, heff, lay):
        """Centreline state, Richardson numbers and the entrainment factor."""
        table = self.th.table
        centre = table.from_concentration(cc)
        layer = table.from_concentration(cclay)
        temp = centre.temp
        rho = centre.rho
        cp = lay.cp
        rit = 0.0
        if self.th.ambient.isofl == 0 and self.th.ambient.ihtfl != 0:
            # ADDHEAT is called with the same RHO and TEMP variables the
            # preceding ADIABAT wrote, and overwrites both. Everything after
            # this point -- the Richardson number, the spreading rate, gamma
            # -- uses the *heated* centreline density, not the adiabatic one.
            hot = self.th.add_heat(cc, self.p.dellay * dh)
            temp, rho = hot.temp, hot.rho
            # ADDHEAT also overwrites CP, and this is the value SURFAC is
            # given further down -- the centreline one, not the layer one
            # computed inside the material-balance loop.
            cp = hot.cp
            wind = self.case.u0 * (heff / self.case.z0) ** self.alpha
            rit = max(
                9.81 * (self.th.ambient.tsurf - temp) / temp * heff / self.ustar / wind,
                0.0,
            )
        ristr = richardson_star(rho, self.rhoa, heff, self.ustar)
        return centre, layer, temp, rho, rit, phi(ristr, rit, self.p.iphifl), cp

    def _surface(self, lay, layer, heff, cp, layer_temp):
        """Ground exchange, with the argument mix the Fortran actually passes.

        ``ADDHEAT`` writes ``rholay``, ``temlay`` and ``cp``; the ``ADIABAT``
        call that follows puts its density in a *different* variable
        (``rholam``), so ``SURFAC`` always receives the heated layer density
        and the centreline heat capacity.

        The layer *temperature* differs between the two stages, and only by a
        letter: ``PSS`` writes the second ``ADIABAT``'s temperature into
        ``temlay``, overwriting the heated value, while ``SSG`` writes it into
        ``temlam`` and leaves ``temlay`` heated. So ``PSS`` computes the
        ground flux from the adiabatic layer temperature and ``SSG`` from the
        heated one -- and the same variable also drives the cut-off test that
        zeroes the flux. Whether that is intentional is not recoverable from
        the source; it is reproduced either way.
        """
        heigh = heff * self.p.dellay
        yw = min(max(1.0 - layer.yc - layer.ya, 0.0), 1.0)
        return surface_exchange(
            layer_temp, heigh, lay.rho, layer.wm, cp, yw,
            tsurf=self.th.ambient.tsurf, pamb=self.th.ambient.pamb,
            zr=self.case.zr, u0=self.case.u0, z0=self.case.z0,
            alpha=self.alpha, ustar=self.ustar,
            isofl=self.th.ambient.isofl, ihtfl=self.th.ambient.ihtfl,
            iwtfl=self.th.ambient.iwtfl, htco=self.th.ambient.htco,
            wtco=self.th.ambient.wtco,
            vapour_pressure=self.th.backend.water_vapour_pressure,
        )

    def _flammable_thresholds(self, gamma):
        """Concentrations at the upper and lower levels of concern.

        With heat transfer active the lookup uses the ``ifl = -2`` branch of
        ``ADIABAT``, which takes the local density-concentration slope
        ``gamma`` instead of the table -- the cloud is no longer on the
        adiabatic mixing line once the ground has warmed it.
        """
        gas = self.case.gas
        table = self.th.table
        if self.th.ambient.isofl == 1 or self.th.ambient.ihtfl == 0:
            return (
                table.from_mole_fraction(gas.ulc).cc,
                table.from_mole_fraction(gas.llc).cc,
            )
        return (
            _adiabat_m2(table, gas, self.th.ambient, gas.ulc, gamma),
            _adiabat_m2(table, gas, self.th.ambient, gas.llc, gamma),
        )

    # -- the dense phase ---------------------------------------------------

    def pss(self, dist, y, dery, erate):
        """Port of ``PSS``."""
        c = self.consts
        a1 = self.alpha1
        bb = y[P_BEFF] - SQRTPI / 2.0 * math.sqrt(y[P_SY2])

        def cc_of_sz(sz):
            return (
                erate * a1 / 2.0 / self.case.u0
                * (self.case.z0 / sz) ** self.alpha / sz / y[P_BEFF]
            )

        sz, cc, cclay, lay = self._close_material_balance(
            cc_of_sz, y[P_RHOUH], y[P_DH], 2.0e-3, 100, 32
        )
        self.sz_seed = sz
        heff = self.gammaf / a1 * sz
        centre, layer, temp, rho, rit, phi_v, cp = self._thermal(
            cc, cclay, y[P_DH], heff, lay
        )

        cloud = CloudState(
            dist=dist, beff=y[P_BEFF], sz=sz, heff=heff, rho=rho,
            rhoa=self.rhoa, cc=cc, temp=temp,
            wind=self.case.u0 * (max(heff, 1e-6) / self.case.z0) ** self.alpha,
            ustar=self.ustar, mass_flux=float(y[P_RHOUH]),
            own=y[6 : 6 + self.n_closure] if self.n_closure else None,
        )
        dery[P_BEFF] = self.closure.lateral_rate(cloud)

        if y[P_BEFF] > 0.0:
            dery[P_SY2] = (
                8.0 * self.betay / PI * y[P_BEFF] ** 2
                * (self.deltay * SQPIO2 / y[P_BEFF])
                ** (1.0 / self.betay)
            )
        else:
            dery[P_SY2] = 0.0

        _watrte, qrte = self._surface(lay, layer, heff, cp, layer.temp)
        if temp >= self.th.ambient.tsurf or layer.temp >= self.th.ambient.tamb:
            qrte = 0.0

        ruhb = lay.rho * c.ueff0 * (sz / self.case.z0) ** a1 * y[P_BEFF]
        druhb = c.entrain * y[P_BEFF] / phi_v
        dery[P_DH] = (qrte * y[P_BEFF] / self.p.dellay - y[P_DH] * druhb) / ruhb
        dery[P_RHOUH] = (druhb - y[P_RHOUH] * dery[P_BEFF]) / y[P_BEFF]
        if self.n_closure:
            cloud.mass_flux_rate = float(dery[P_RHOUH])
            dery[6 : 6 + self.n_closure] = self.closure.derivatives(cloud)

        gamma = (rho - self.rhoa) / cc
        # transient runs integrate only four states: the flammable-mass
        # integrals are a steady-state report, and SSSUP passes NDIM = 4
        if len(dery) > P_MLOW:
            dery[P_MHI] = dery[P_MLOW] = 0.0
        if self.case.steady_state and len(dery) > P_MLOW:
            chi, clow = self._flammable_thresholds(gamma)
            gamhi0 = 2.0 * cc * bb * sz / a1
            gammax = (
                2.0 * cc * sz * self.gammaf / a1
                * (bb + SQRTPI / 2.0 * math.sqrt(y[P_SY2]))
            )
            if cc > clow:
                wlow = math.log(cc / clow)
                dery[P_MLOW] = min(
                    gaminc(1.0 / a1, wlow) * gamhi0
                    + 2.0 * clow * math.sqrt(y[P_SY2]) * sz / a1 * series(wlow, a1),
                    gammax,
                )
            if cc > chi:
                whi = math.log(cc / chi)
                dery[P_MHI] = min(
                    gaminc(1.0 / a1, whi) * gamhi0
                    + 2.0 * chi * math.sqrt(y[P_SY2]) * sz / a1 * series(whi, a1),
                    gammax,
                )

        return DownwindState(
            dist=dist, cc=cc, bb=bb, sz=sz, sy=math.sqrt(y[P_SY2]),
            yc=centre.yc, rho=rho, temp=temp, gamma=gamma,
            rholay=lay.rho, heff=heff,
        )

    # -- the Gaussian phase -------------------------------------------------

    def ssg(self, dist, y, dery, erate, xv):
        """Port of ``SSG``."""
        c = self.consts
        a1 = self.alpha1
        sy = RT2 * self.deltay * (dist + xv) ** self.betay

        def cc_of_sz(sz):
            return (
                erate * a1 * (self.case.z0 / sz) ** self.alpha
                / self.case.u0 / sz / SQRTPI / sy
            )

        sz, cc, cclay, lay = self._close_material_balance(
            cc_of_sz, y[G_RHOUH], y[G_DH], 2.5e-3, 20, 33
        )
        self.sz_seed = sz
        heff = self.gammaf * sz / a1
        centre, layer, temp, rho, rit, phi_v, cp = self._thermal(
            cc, cclay, y[G_DH], heff, lay
        )

        dery[G_RHOUH] = c.entrain / phi_v
        _watrte, qrte = self._surface(lay, layer, heff, cp, lay.temp)
        if temp >= self.th.ambient.tsurf or lay.temp >= self.th.ambient.tamb:
            qrte = 0.0
        dery[G_DH] = (qrte / self.p.dellay - y[G_DH] * dery[G_RHOUH]) / y[G_RHOUH]

        gamma = (rho - self.rhoa) / cc
        if len(dery) > G_MLOW:
            dery[G_MHI] = dery[G_MLOW] = 0.0
        if self.case.steady_state and len(dery) > G_MLOW:
            chi, clow = self._flammable_thresholds(gamma)
            gammax = SQRTPI * cc * sz * sy * self.gammaf / a1
            if cc > clow:
                wlow = math.log(cc / clow)
                dery[G_MLOW] = min(
                    2.0 * clow * sy * sz / a1 * series(wlow, a1), gammax
                )
            if cc > chi:
                whi = math.log(cc / chi)
                dery[G_MHI] = min(2.0 * chi * sy * sz / a1 * series(whi, a1), gammax)

        return DownwindState(
            dist=dist, cc=cc, bb=0.0, sz=sz, sy=sy, yc=centre.yc, rho=rho,
            temp=temp, gamma=gamma, rholay=lay.rho, heff=heff,
        )


def _adiabat_m2(table, gas, ambient, yc: float, gamma: float) -> float:
    """Port of ``ADIABAT(ifl=-2)``: concentration at a mole fraction, given
    the local density slope ``gamma`` rather than the adiabatic table.

    Used once the ground has heated the cloud off the adiabatic mixing line,
    where the tabulated density no longer applies but the linear relation
    :math:`\\rho = \\rho_a + \\gamma c` still does.
    """
    from .constants import WMA, WMW

    ycl = max(yc, 0.0)
    ya = (1.0 - (1.0 + gas.mw * ambient.humsrc / WMW) * ycl) / (
        1.0 + ambient.humid * WMA / WMW
    )
    yw = 1.0 - ya - ycl
    wm = ycl * gas.mw + ya * WMA + yw * WMW
    wc = gas.mw / wm * ycl
    return wc * table.rhoa / (1.0 - gamma * wc)
