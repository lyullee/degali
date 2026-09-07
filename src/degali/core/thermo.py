"""Mixture thermodynamics: contaminant + dry air + water (vapour and liquid).

Ports ``TPROP.FOR`` in full:

=====================  ==================================================
Fortran                Python
=====================  ==================================================
``WATVP``              :meth:`ThermoBackend.water_vapour_pressure`
``CPC``                :meth:`ThermoBackend.cp_contaminant`
``ENTHAL``             :meth:`Thermo.enthalpy`
``TPROP``              :meth:`Thermo.properties`
``SETENT``             :meth:`Thermo.reference_enthalpies`
``SETDEN``             :meth:`Thermo.build_adiabatic_table`
``ADIABAT``            :class:`AdiabaticTable`
``ADDHEAT``            :meth:`Thermo.add_heat`
=====================  ==================================================

Physical model
--------------
The mixture is treated as three species: contaminant (mole fraction ``yc``),
dry air (``ya``) and water (``yw``).  Air and water vapour are ideal gases;
water in excess of saturation condenses to liquid at a fixed density; the
contaminant contributes a partial volume ``wc * T / T_e / rho_e``, which is
the ideal-gas scaling of its user-supplied release-state density and is what
lets the same expression describe an aerosol.  Enthalpy is referenced to the
ambient temperature.

Backends
--------
``backend="legacy"``
    Bit-faithful to 1989.  Water vapour pressure is
    :math:`\\exp(14.683943 - 5407/T)` atm, the gas constant is ``0.08205``,
    air and water-vapour heat capacities are constants, and the contaminant
    heat capacity follows the two-parameter ``CPK``/``CPP`` correlation.
``backend="coolprop"``
    Water properties (vapour pressure, vapour ``cp``, liquid density) and, when
    a CoolProp fluid name is supplied, the contaminant ``cp`` and density come
    from CoolProp's equations of state.  Everything structural -- the mixing
    rules, the condensation criterion, the enthalpy reference -- is unchanged,
    so results remain comparable term by term.

Use ``legacy`` to prove the port, ``coolprop`` to improve it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .constants import (
    ATM_TO_PA,
    CPA,
    CPW,
    DHFUS,
    DHVAP,
    IGEN,
    RGAS,
    RHOWL,
    WMA,
    WMW,
)
from .numerics import RootBracketError, brentq, limit, zbrent

Backend = Literal["legacy", "coolprop"]


# ==========================================================================
# Property backends
# ==========================================================================


class ThermoBackend:
    """Interface for pure-component property evaluation."""

    name: str = "abstract"

    def water_vapour_pressure(self, temp: float) -> float:
        """Water vapour pressure, atm."""
        raise NotImplementedError

    def cp_air(self, temp: float) -> float:
        return CPA

    def cp_water_vapour(self, temp: float) -> float:
        return CPW

    def rho_liquid_water(self, temp: float) -> float:
        return RHOWL

    def latent_heat(self, temp: float) -> float:
        """Combined heat of vaporisation and (below 273.15 K) fusion, J/kg."""
        frac = 0.0
        if temp < 273.15:
            frac = min((273.15 - temp) / 10.0, 1.0)  # ENTHAL: deltaf = 10 K
        return DHVAP + DHFUS * frac


class LegacyBackend(ThermoBackend):
    """The 1989 correlations, unchanged."""

    name = "legacy"

    def water_vapour_pressure(self, temp: float) -> float:
        # WATVP.  The commented-out alternative in the Fortran,
        #   6.0298e-3 * exp(5407*(1/273.15 - 1/T)),
        # is algebraically the same expression.
        return math.exp(14.683943 - 5407.0 / temp)


class CoolPropBackend(ThermoBackend):
    """Water and (optionally) contaminant properties from CoolProp EOS.

    Every property here is a function of temperature alone at a fixed ambient
    pressure, and the root finders that drive the model evaluate them
    thousands of times over a narrow range.  A direct ``PropsSI`` call costs
    tens of microseconds, which makes an uncached backend about a hundred
    times slower than the correlations it replaces -- slow enough that a
    transient run stops being practical.

    So each property is tabulated once on a fine grid and interpolated
    thereafter.  Interpolation keeps the result a smooth function of
    temperature, which matters because the integrators effectively
    differentiate these properties: caching on a rounded temperature would be
    faster still, but it would put small steps into them and the adaptive step
    control would chase those steps.

    Parameters
    ----------
    contaminant
        A CoolProp fluid name (e.g. ``"Methane"``, ``"Ammonia"``), usually
        supplied by :func:`degali.core.fluids.resolve`.  Without one, water
        and air still come from the EOS and the contaminant keeps the 1989
        correlations.
    """

    name = "coolprop"

    #: Points in each interpolation grid.
    POINTS = 4001
    #: Default grid span, K.  Wide enough for any release DEGADIS models --
    #: cryogenic liquefied gases at one end, heated ground at the other -- so
    #: that it is built once and never rebuilt.  Growing a grid outward as
    #: queries arrive instead costs hundreds of rebuilds, because the model
    #: walks the temperature range gradually rather than jumping to its
    #: extremes.
    SPAN = (80.0, 400.0)
    #: Margin added when a query does fall outside, K.
    MARGIN = 50.0
    #: Absolute floor on any grid, K.  Liquid hydrogen boils at 20 K, well
    #: below :attr:`SPAN`, and simply widening by :attr:`MARGIN` puts the
    #: lower end of the grid below absolute zero.  That is not merely untidy:
    #: a sizeable part of a 4001-point grid is then spent on temperatures that
    #: do not exist, and the interpolation that remains is coarser for it.
    FLOOR = 4.0

    def __init__(
        self, contaminant: str | None = None, span=None, *,
        force_contaminant_gas: bool = False,
    ):
        from CoolProp.CoolProp import PropsSI  # imported lazily

        self._props = PropsSI
        self.contaminant = contaminant
        # A flashing-jet source is explicitly a gas plume with any remaining
        # droplets carried as a separate mass.  At cryogenic temperatures a
        # plain T,P query can instead select bulk liquid at the *total*
        # ambient pressure, even though the contaminant in the mixture is at
        # its partial pressure.  Opt in for those gas-plume paths rather than
        # silently changing the historical dense-gas/pool backend everywhere.
        self.force_contaminant_gas = force_contaminant_gas
        self._tmin_water = 273.16  # CoolProp water triple point
        self._grids: dict[str, tuple] = {}
        #: Properties for which an equation-of-state lookup failed and a
        #: fallback was used, with a count.
        #:
        #: These paths used to fail silently: a caller asking for the
        #: ``coolprop`` backend could get the 1989 correlation for some
        #: property over some temperature range and never know. That is the
        #: same class of defect as a misplaced bare ``raise`` -- the answer
        #: changes and the error does not surface -- and it cost this project
        #: a day once already, when a two-phase gap was being filled with the
        #: nearest valid value.
        #:
        #: Nothing is silenced now; it is counted, and
        #: :meth:`fallback_report` says what happened.
        self.fallbacks: dict[str, int] = {}
        if span is None and contaminant is not None:
            # start the grid at the contaminant's own boiling point rather
            # than at a value chosen for LNG: hydrogen boils at 20 K and
            # would otherwise force a rebuild on the first cold lookup
            try:
                boil = float(PropsSI("T", "P", ATM_TO_PA, "Q", 0, contaminant))
                span = (max(min(self.SPAN[0], boil - 20.0), self.FLOOR),
                        self.SPAN[1])
            except Exception:
                span = None
        self.SPAN = span or self.SPAN

    def _interp(self, key: str, temp: float, compute, pressure: float = 0.0) -> float:
        """Evaluate a property from a precomputed grid, extending it if needed.

        The grid spans whatever range the model has asked for so far, plus a
        margin, and is rebuilt when a query falls outside it. Linear
        interpolation on 4001 points is smooth, deterministic and about two
        orders of magnitude cheaper than the underlying call -- rounding the
        temperature to a cache key instead would be faster still but would put
        small steps into functions the integrators differentiate.
        """
        grid = self._grids.get(key)
        if grid is None or not (grid[0] <= temp <= grid[1]):
            lo, hi = (grid[0], grid[1]) if grid else self.SPAN
            if temp < lo:
                lo = max(temp - self.MARGIN, self.FLOOR)
            if temp > hi:
                hi = temp + self.MARGIN
            xs = np.linspace(lo, hi, self.POINTS)
            ys = np.empty(self.POINTS)
            ok = np.zeros(self.POINTS, dtype=bool)
            for i, x in enumerate(xs):
                try:
                    ys[i] = compute(float(x), pressure)
                    ok[i] = True
                except Exception:
                    ok[i] = False
            if not ok.any():
                raise ValueError(f"no usable {key!r} data over {lo}-{hi} K")
            if not ok.all():
                # Part of the span is outside the fluid's single-phase region
                # -- a cryogenic release sits below its own saturation line at
                # ambient pressure. Hold the nearest valid value across the
                # gap rather than abandoning the whole grid, which would
                # silently drop back to the 1989 correlation for the entire
                # run.
                ys = np.interp(xs, xs[ok], ys[ok])
            grid = (lo, hi, xs, ys)
            self._grids[key] = grid
        return float(np.interp(temp, grid[2], grid[3]))

    def _vapour_pressure(self, t: float, _p: float = 0.0) -> float:
        # Below the triple point CoolProp's saturation line is undefined for
        # 'Water'; extrapolate the ice sublimation curve with the legacy form,
        # scaled so the two agree at the triple point.
        if t < self._tmin_water:
            p_tp = (
                self._props("P", "T", self._tmin_water, "Q", 0, "Water") / ATM_TO_PA
            )
            legacy_tp = math.exp(14.683943 - 5407.0 / self._tmin_water)
            return math.exp(14.683943 - 5407.0 / t) * (p_tp / legacy_tp)
        return self._props("P", "T", t, "Q", 0, "Water") / ATM_TO_PA

    def water_vapour_pressure(self, temp: float) -> float:
        return self._interp("pw", temp, self._vapour_pressure)

    def cp_air(self, temp: float) -> float:
        # Air is two-phase below about 79 K at one atmosphere and CoolProp
        # declines to give a pseudo-pure heat capacity there. The grid can
        # reach that far when the release is cryogenic, so the evaluation is
        # clamped; at those temperatures the mixture is almost pure
        # contaminant and the air term carries no weight.
        return self._interp(
            "cpa", temp,
            lambda t, _p: float(
                self._props("C", "T", max(t, 100.0), "P", ATM_TO_PA, "Air")
            ),
        )

    def cp_water_vapour(self, temp: float) -> float:
        return self._interp(
            "cpw", temp,
            lambda t, _p: float(
                self._props("C", "T", max(t, self._tmin_water + 0.1), "Q", 1, "Water")
            ),
        )

    def rho_liquid_water(self, temp: float) -> float:
        return self._interp(
            "rhow", temp,
            lambda t, _p: float(
                self._props(
                    "D", "T", min(max(t, self._tmin_water), 640.0), "Q", 0, "Water"
                )
            ),
        )

    def _latent(self, t: float, _p: float = 0.0) -> float:
        tt = max(t, self._tmin_water + 0.1)
        try:
            dh = float(
                self._props("H", "T", tt, "Q", 1, "Water")
                - self._props("H", "T", tt, "Q", 0, "Water")
            )
        except Exception:
            self.fallbacks["latent heat"] = self.fallbacks.get("latent heat", 0) + 1
            dh = DHVAP
        if t < 273.15:
            dh += DHFUS * min((273.15 - t) / 10.0, 1.0)
        return dh

    def latent_heat(self, temp: float) -> float:
        return self._interp("dhv", temp, self._latent)

    def cp_contaminant_eos(self, temp: float, pamb_atm: float) -> float | None:
        if self.contaminant is None:
            return None
        try:
            return self._interp(
                "cpc", temp,
                lambda t, p: float(
                    self._props(
                        "C", "T|gas" if self.force_contaminant_gas else "T",
                        t, "P", p * ATM_TO_PA, self.contaminant,
                    )
                ),
                pamb_atm,
            )
        except Exception:
            self.fallbacks["cp of contaminant"] = (
                self.fallbacks.get("cp of contaminant", 0) + 1
            )
            return None

    def rho_contaminant_eos(self, temp: float, pamb_atm: float) -> float | None:
        if self.contaminant is None:
            return None
        try:
            return self._interp(
                "rhoc", temp,
                lambda t, p: float(
                    self._props(
                        "D", "T|gas" if self.force_contaminant_gas else "T",
                        t, "P", p * ATM_TO_PA, self.contaminant,
                    )
                ),
                pamb_atm,
            )
        except Exception:
            self.fallbacks["density of contaminant"] = (
                self.fallbacks.get("density of contaminant", 0) + 1
            )
            return None


    def fallback_report(self) -> str:
        """What the equation of state could not supply, and how often.

        An empty string means every property came from the equation of state,
        which is what asking for this backend is meant to get.
        """
        if not self.fallbacks:
            return ""
        return "; ".join(
            f"{name}: {count} lookups fell back"
            for name, count in sorted(self.fallbacks.items())
        )


def make_backend(backend: Backend = "legacy", *, contaminant: str | None = None):
    if backend == "legacy":
        return LegacyBackend()
    if backend == "coolprop":
        return CoolPropBackend(contaminant=contaminant)
    raise ValueError(f"unknown backend {backend!r}")


# ==========================================================================
# State containers
# ==========================================================================


@dataclass
class GasProperties:
    """Contaminant description -- ``/cgprop/`` in the Fortran."""

    name: str = "GAS"
    mw: float = WMA  #: GASMW, kg/kmol
    temp: float = 298.0  #: GASTEM, release temperature, K
    rho: float = 1.0  #: GASRHO, density at release T and ambient P, kg/m**3
    cpk: float = 0.0  #: GASCPK, q1 in the heat-capacity correlation
    cpp: float = 1.0  #: GASCPP, p1 in the heat-capacity correlation
    ulc: float = 1.0  #: GASULC, upper level of concern, mole fraction
    llc: float = 1.0e-12  #: GASLLC, lower level of concern, mole fraction
    zzc: float = 0.0  #: GASZZC, elevation for contour calculations, m
    coolprop_name: str | None = None  #: EOS fluid name, if available


@dataclass
class AmbientConditions:
    """Atmosphere and surface -- ``/comatm/`` in the Fortran."""

    tamb: float = 298.0  #: ambient temperature, K
    pamb: float = 1.0  #: ambient pressure, atm
    humid: float = 0.0  #: absolute humidity, kg water / kg bone-dry air
    tsurf: float = 298.0  #: surface temperature, K
    isofl: int = 0  #: 1 = "isothermal" run, density read from a table
    ihtfl: int = 0  #: heat-transfer flag (0 = adiabatic mixing only)
    iwtfl: int = 0  #: ground-to-cloud water transfer flag
    htco: float = 0.0  #: fixed heat-transfer coefficient
    wtco: float = 0.0  #: fixed water-transfer coefficient
    humsrc: float = 0.0  #: water carried by the source, kg water / kg contam


@dataclass
class MixtureState:
    """Result of a property evaluation."""

    wc: float  #: contaminant mass fraction
    wa: float  #: dry-air mass fraction
    yc: float  #: contaminant mole fraction
    ya: float  #: dry-air mole fraction
    wm: float  #: mixture molecular weight, kg/kmol
    temp: float  #: K
    rho: float  #: kg/m**3
    cp: float  #: J/(kg K)
    enthalpy: float  #: J/kg, referenced to ``tamb``

    @property
    def ww(self) -> float:
        return 1.0 - self.wc - self.wa

    @property
    def yw(self) -> float:
        return 1.0 - self.yc - self.ya

    @property
    def cc(self) -> float:
        """Contaminant concentration, kg contaminant / m**3 mixture."""
        return self.wc * self.rho


# ==========================================================================
# Adiabatic mixing table
# ==========================================================================


@dataclass
class AdiabaticTable:
    """Port of ``/GEN2/ DEN(5, IGEN)`` and the ``ADIABAT`` lookups.

    Columns, in the Fortran's order:

    ==========  ==========================================
    ``yc``      contaminant mole fraction
    ``cc``      contaminant concentration, kg/m**3
    ``rho``     mixture density, kg/m**3
    ``h``       mixture enthalpy, J/kg
    ``t``       mixture temperature, K
    ==========  ==========================================

    Row 0 is pure ambient air, the last row is pure contaminant.  All lookups
    are piecewise linear and, crucially, *interpolate in different variables
    depending on the entry point* -- ``ADIABAT`` interpolates density against
    concentration for ``ifl=0`` but against mass fraction for ``ifl=1``.  That
    asymmetry is deliberate in the original and is preserved here.
    """

    yc: np.ndarray
    cc: np.ndarray
    rho: np.ndarray
    h: np.ndarray
    t: np.ndarray
    humid: float = 0.0
    humsrc: float = 0.0
    gasmw: float = WMA

    @property
    def rhoa(self) -> float:
        """Pure ambient air density, ``DEN(3,1)``."""
        return float(self.rho[0])

    @property
    def rhoe(self) -> float:
        """Pure contaminant density, ``DEN(3,NP)``."""
        return float(self.rho[-1])

    @property
    def n(self) -> int:
        return len(self.yc)

    def as_array(self) -> np.ndarray:
        """The table as a ``(n, 5)`` array in ``DEN`` column order."""
        return np.column_stack([self.yc, self.cc, self.rho, self.h, self.t])

    # -- the four ADIABAT entry points -------------------------------------

    def _tail_interp(self, i: int, wcl: float) -> tuple[float, float]:
        """Label 8000 of ``ADIABAT``: enthalpy and temperature by mass fraction."""
        w1 = self.cc[i - 1] / self.rho[i - 1]
        w2 = self.cc[i] / self.rho[i]
        slope_h = (self.h[i] - self.h[i - 1]) / (w2 - w1)
        enthalpy = (wcl - w1) * slope_h + self.h[i - 1]
        slope_t = (self.t[i] - self.t[i - 1]) / (w2 - w1)
        temp = (wcl - w1) * slope_t + self.t[i - 1]
        return float(enthalpy), float(temp)

    def from_concentration(self, cc: float) -> MixtureState:
        """``ADIABAT(ifl=0)``: everything from contaminant concentration."""
        ccl = max(cc, 0.0)
        i = 1
        while i < self.n and ccl > self.cc[i]:
            i += 1
        if i >= self.n:  # past the end of the table: clamp, as the Fortran does
            i = self.n - 1
            ccl = min(ccl, float(self.cc[i]))
        slope = (self.rho[i] - self.rho[i - 1]) / (self.cc[i] - self.cc[i - 1])
        rho = (ccl - self.cc[i - 1]) * slope + self.rho[i - 1]
        wc = ccl / rho
        wa = (1.0 - (1.0 + self.humsrc) * wc) / (1.0 + self.humid)
        ww = 1.0 - wa - wc
        wm = 1.0 / (wc / self.gasmw + wa / WMA + ww / WMW)
        enthalpy, temp = self._tail_interp(i, wc)
        return MixtureState(
            wc=wc,
            wa=wa,
            yc=wm / self.gasmw * wc,
            ya=wm / WMA * wa,
            wm=wm,
            temp=temp,
            rho=float(rho),
            cp=CPA,
            enthalpy=enthalpy,
        )

    def from_mass_fraction(self, wc: float, wa: float | None = None) -> MixtureState:
        """``ADIABAT(ifl=1)``: everything from contaminant mass fraction.

        ``wa`` is genuinely an *input* on this path.  In the Fortran it is a
        dummy argument that ``ADIABAT`` writes only in the out-of-range
        branches (``wc < 0`` or ``wc > 1``); for a normal call it reads
        whatever the caller passed and uses it to form the molecular weight,
        hence the mole fraction, hence which panel of the table is
        interpolated on.

        That makes the result depend on a value most callers never set.
        ``TPROP`` passes its own consistent ``wa``, but ``SZLOCAL`` passes an
        uninitialised local, which under ``/noauto`` is static and therefore
        zero -- and stays zero, because nothing ever writes it.  So ``SZF``
        effectively runs its lookups as though the mixture contained no dry
        air at all, selecting panels one or two positions off and shifting the
        density by up to 0.5 per cent.

        This is a bug in DEGADIS, not in the port.  It is reproduced because
        the whole point of the port is to reproduce DEGADIS; callers that want
        the intended behaviour pass ``wa`` explicitly.  ``None`` keeps the
        thermodynamically consistent value.
        """
        wcl = wc
        if wa is None:
            wa = (1.0 - (1.0 + self.humsrc) * wc) / (1.0 + self.humid)
        if wc < 0.0:
            wcl = 0.0
            wa = 1.0 / (1.0 + self.humid)
        if wc > 1.0:
            wcl = 1.0
            wa = 0.0
        ww = 1.0 - wa - wcl
        wm = 1.0 / (wcl / self.gasmw + wa / WMA + ww / WMW)
        yc = wm / self.gasmw * wcl
        ya = wm / WMA * wa

        i = 1
        while i < self.n and yc > self.yc[i]:
            i += 1
        if i >= self.n:
            i = self.n - 1  # extrapolate on the last panel
        slope = (self.rho[i] - self.rho[i - 1]) / (self.cc[i] - self.cc[i - 1])
        rho = (self.rho[i - 1] - self.cc[i - 1] * slope) / (1.0 - slope * wcl)

        j = 1
        while j < self.n and wcl > self.cc[j] / self.rho[j]:
            j += 1
        if j >= self.n:
            j = self.n - 1
        enthalpy, temp = self._tail_interp(j, wcl)
        return MixtureState(
            wc=wcl,
            wa=wa,
            yc=yc,
            ya=ya,
            wm=wm,
            temp=temp,
            rho=float(rho),
            cp=CPA,
            enthalpy=enthalpy,
        )

    def from_mole_fraction(self, yc: float) -> MixtureState:
        """``ADIABAT(ifl=2)``: everything from contaminant mole fraction."""
        i = 1
        while i < self.n and yc > self.yc[i]:
            i += 1
        if i >= self.n:
            i = self.n - 1
        f = (yc - self.yc[i - 1]) / (self.yc[i] - self.yc[i - 1])
        cc = float(self.cc[i - 1] + f * (self.cc[i] - self.cc[i - 1]))
        return self.from_concentration(cc)


# ==========================================================================
# The main thermodynamic object
# ==========================================================================


@dataclass
class Thermo:
    """Mixture property evaluator -- the Python face of ``TPROP``."""

    gas: GasProperties
    ambient: AmbientConditions
    backend: ThermoBackend = field(default_factory=LegacyBackend)
    table: AdiabaticTable | None = None

    #: Enthalpy of the material streams entering the cloud, J/kg (``SETENT``).
    hmrte: float = 0.0  #: contaminant from the primary source
    harte: float = 0.0  #: entrained ambient air (zero by construction)
    hwrte: float = 0.0  #: water evaporated from the surface

    #: When true, reproduce DEGADIS's own numerics: invert enthalpy for
    #: temperature with the ported ``ZBRENT`` at its 1e-3 K tolerance, and
    #: form heat capacities as quotients of differences.  When false, solve
    #: the inversion to machine precision and take the heat capacity as an
    #: analytic derivative.  See :meth:`add_heat` for why that matters.
    legacy_numerics: bool = True

    # -- pure-component pieces ---------------------------------------------

    def cp_contaminant(self, temp: float) -> float:
        """Port of ``CPC``.

        The Fortran models the *molar* heat capacity as
        ``3.33e4 + q1 * p1 * T**(p1-1)`` and, when evaluating a mean over a
        temperature interval, replaces the derivative term by the divided
        difference ``q1 * (T**p1 - Te**p1)/(T - Te)``.  The result is divided
        by the molecular weight to give J/(kg K).

        With ``backend="coolprop"`` and a known fluid this is replaced by the
        EOS value averaged the same way.
        """
        te = self.gas.temp
        if isinstance(self.backend, CoolPropBackend):
            cp = self.backend.cp_contaminant_eos(temp, self.ambient.pamb)
            if cp is not None:
                return cp
        con = 3.33e4
        if temp != te:
            cpc = con + self.gas.cpk * (temp**self.gas.cpp - te**self.gas.cpp) / (
                temp - te
            )
        else:
            cpc = con + self.gas.cpk * self.gas.cpp * te ** (self.gas.cpp - 1.0)
        return cpc / self.gas.mw

    def rho_contaminant(self, temp: float) -> float:
        """Contaminant density at ``temp`` and ambient pressure, kg/m**3.

        Legacy: ideal scaling of the release-state density,
        ``rho_e * T_e / T``.  CoolProp: the EOS value when available.
        """
        if isinstance(self.backend, CoolPropBackend):
            rho = self.backend.rho_contaminant_eos(temp, self.ambient.pamb)
            if rho is not None:
                return rho
        return self.gas.rho * self.gas.temp / temp

    # -- enthalpy ----------------------------------------------------------

    def _condensate(self, wc: float, wa: float, wm: float, temp: float) -> float:
        """Liquid water mass fraction, ``conden`` in ``ENTHAL``/``TPROP``."""
        ww = 1.0 - wc - wa
        ya = wa * wm / WMA
        yc = wc * wm / self.gas.mw
        ywsat = self.backend.water_vapour_pressure(temp) / self.ambient.pamb
        if ywsat >= 1.0:
            return 0.0
        wwsat = WMW / wm * ywsat * (ya + yc) / (1.0 - ywsat)
        return max(0.0, ww - wwsat)

    def enthalpy(self, wc: float, wa: float, temp: float) -> float:
        """Port of ``ENTHAL``: mixture enthalpy referenced to ``tamb``, J/kg."""
        amb = self.ambient
        ww = 1.0 - wa - wc
        wm = 1.0 / (wc / self.gas.mw + wa / WMA + ww / WMW)
        conden = self._condensate(wc, wa, wm, temp)
        dh = self.backend.latent_heat(temp)
        return (
            wc * self.cp_contaminant(temp) * (temp - amb.tamb)
            - conden * dh
            + ww * self.backend.cp_water_vapour(temp) * (temp - amb.tamb)
            + wa * self.backend.cp_air(temp) * (temp - amb.tamb)
        )

    def reference_enthalpies(self) -> tuple[float, float, float]:
        """Port of ``SETENT``.  Sets and returns ``(hmrte, harte, hwrte)``."""
        amb = self.ambient
        if amb.isofl == 1:
            self.hmrte = self.harte = self.hwrte = 0.0
            return 0.0, 0.0, 0.0
        self.hmrte = self.cp_contaminant(self.gas.temp) * (self.gas.temp - amb.tamb)
        self.harte = 0.0  # entrained air is at tamb, the reference state
        self.hwrte = (
            self.backend.cp_water_vapour(amb.tsurf) * (amb.tsurf - amb.tamb)
            if amb.iwtfl != 0
            else 0.0
        )
        return self.hmrte, self.harte, self.hwrte

    def heat_capacity(
        self, wc: float, wa: float, temp: float, *, half_width: float = 5.0
    ) -> float:
        """Mean mixture heat capacity about ``temp``, J/(kg K).

        A *secant* over a finite interval, not a derivative.  Mixture enthalpy
        has a slope discontinuity where water begins to condense, so the true
        derivative jumps there; a secant over a finite interval does not,
        because enthalpy itself is continuous.  Narrowing the interval towards
        a derivative recovers the jump and makes the downwind integrator
        bisect to a standstill -- which is how this width was chosen.

        ``TPROP`` uses a one-sided 10 K difference; this is the same width
        centred, which removes the bias.
        """
        hi = self.enthalpy(wc, wa, temp + half_width)
        lo = self.enthalpy(wc, wa, temp - half_width)
        return (hi - lo) / (2.0 * half_width)

    def _invert_enthalpy(
        self, wc: float, wa: float, enth: float, tmin: float, tmax: float,
        tol: float,
    ) -> float:
        """Temperature at which the mixture holds ``enth``.

        ``legacy_numerics`` selects between the ported ``ZBRENT`` at the
        original's tolerance and an accurate solve.
        """
        residual = lambda t: enth - self.enthalpy(wc, wa, t)
        if self.legacy_numerics:
            return zbrent(residual, tmin, tmax, tol)
        return brentq(residual, tmin, tmax, xtol=1e-12, rtol=1e-14)

    # -- density -----------------------------------------------------------

    def density(self, wc: float, wa: float, wm: float, temp: float) -> float:
        """Partial-volume mixing rule, label 400 of ``TPROP``.

        .. math::
            \\frac{1}{\\rho} = \\frac{w_a}{\\rho_a(T)}
                             + \\frac{w_w - w_{cond}}{\\rho_w(T)}
                             + \\frac{w_{cond}}{\\rho_{w,liq}}
                             + \\frac{w_c}{\\rho_c(T)}
        """
        amb = self.ambient
        ww = 1.0 - wc - wa
        conden = self._condensate(wc, wa, wm, temp)
        rho_air = amb.pamb * WMA / RGAS / temp
        rho_wv = amb.pamb * WMW / RGAS / temp
        v = (
            wa / rho_air
            + (ww - conden) / rho_wv
            + conden / self.backend.rho_liquid_water(temp)
            + wc / self.rho_contaminant(temp)
        )
        return 1.0 / v

    # -- the TPROP entry point ---------------------------------------------

    def properties(
        self,
        wc: float,
        wa: float,
        enth: float,
        *,
        ifl: int = 1,
        temp: float | None = None,
        tol: float = 1.0e-3,
    ) -> MixtureState:
        """Port of ``TPROP``.

        Parameters
        ----------
        wc, wa
            Contaminant and dry-air mass fractions.
        enth
            Mixture enthalpy, J/kg.  Ignored (and recomputed) for ``ifl == 0``.
        temp
            Required for ``ifl == -1`` only.  In the Fortran the temperature is
            a separate argument that ``TPROP`` leaves untouched on that path,
            so the caller supplies it and only the density is computed.
        ifl
            The Fortran's dispatch flag:

            * ``-1`` -- density only, at the supplied ``temp``.
            * ``0``  -- compute the enthalpy of adiabatic mixing first.
            * ``1``  -- from enthalpy, but short-circuit through the adiabatic
              table when heat transfer is switched off.
            * ``2``  -- from enthalpy, always by inversion (used by ``SETDEN``).
        """
        amb = self.ambient
        gas = self.gas
        ww = 1.0 - wc - wa
        wm = 1.0 / (wc / gas.mw + wa / WMA + ww / WMW)
        yc = wm / gas.mw * wc
        ya = wm / WMA * wa

        # "isothermal" simulation: the density relation is a supplied table
        if amb.isofl == 1:
            if self.table is None:
                raise RuntimeError("isofl=1 requires an adiabatic/density table")
            return self.table.from_mass_fraction(wc)

        if ifl == 0:
            enth = wc * self.cp_contaminant(gas.temp) * (gas.temp - amb.tamb) + (
                ww - wa * amb.humid
            ) * self.backend.cp_water_vapour(amb.tsurf) * (amb.tsurf - amb.tamb)

        # heat transfer off => adiabatic mixing is exact, use the table
        if ifl == 1 and amb.ihtfl == 0:
            if self.table is None:
                raise RuntimeError("ihtfl=0 requires the adiabatic table")
            return self.table.from_mass_fraction(wc)

        if ifl == -1:
            if temp is None:
                raise ValueError("ifl=-1 requires an explicit temperature")
        else:
            tmin = min(gas.temp, amb.tsurf, amb.tamb)
            tmax = max(gas.temp, amb.tsurf, amb.tamb)
            elow = self.enthalpy(wc, wa, tmin)
            ehigh = self.enthalpy(wc, wa, tmax)
            if enth < elow:
                temp, enth = tmin, elow
            elif enth > ehigh:
                temp, enth = tmax, ehigh
            else:
                temp = self._invert_enthalpy(wc, wa, enth, tmin, tmax, tol)

        rho = self.density(wc, wa, wm, temp)

        if self.legacy_numerics:
            # TPROP's mean cp over a 10 K interval above temp, floored at cp_air
            t2 = temp + 10.0
            h2 = self.enthalpy(wc, wa, t2)
            cp = (enth - h2) / (temp - t2)
            cp = max(cp, self.backend.cp_air(temp))
        else:
            cp = max(
                self.heat_capacity(wc, wa, temp), self.backend.cp_air(temp)
            )

        return MixtureState(
            wc=wc, wa=wa, yc=yc, ya=ya, wm=wm, temp=temp, rho=rho, cp=cp, enthalpy=enth
        )

    # -- adiabatic mixing table --------------------------------------------

    def build_adiabatic_table(
        self,
        wc: float = 1.0,
        wa: float = 0.0,
        enthalpy: float | None = None,
        *,
        npts: int = 200,
        tcrit: float = 0.002,
        iback: int = 25,
        maxrows: int | None = None,
        exact_grid: bool = False,
    ) -> AdiabaticTable:
        """Port of ``SETDEN``.

        Walks the adiabatic mixing line from pure air to pure contaminant on a
        uniform grid of ``npts`` dilution steps and records only those points
        needed so that linear interpolation between recorded points reproduces
        every skipped point to within ``tcrit`` relative error, looking back at
        most ``iback`` points.  This adaptive thinning is why the printed table
        has irregular spacing, and reproducing it exactly is required to match
        reference output -- every later lookup interpolates on *these* nodes.
        """
        amb = self.ambient
        gas = self.gas
        if enthalpy is None:
            enthalpy = self.hmrte
        if amb.isofl == 1:
            raise RuntimeError("SETDEN is not used for isofl=1 runs")

        # IGEN = 42 is a Fortran array dimension, not a physical limit, and it
        # is adequate only for the substances DEGADIS was built for. The
        # adaptive thinning keeps whatever nodes linear interpolation needs,
        # so a strongly curved mixing line needs more of them: liquid
        # hydrogen spans 20 to 289 K, a factor of fourteen, against LNG's
        # factor of under three, and overflows 42. The bound is raised rather
        # than the tolerance loosened, because loosening it would silently
        # degrade every lookup.
        if maxrows is None:
            maxrows = max(IGEN, npts)

        zero = 1.0e-20
        rows: list[tuple[float, float, float, float, float]] = [
            (
                0.0,
                0.0,
                amb.pamb * (1.0 + amb.humid) * WMW / (RGAS * (WMW / WMA + amb.humid))
                / amb.tamb,
                0.0,
                amb.tamb,
            )
        ]

        ils = npts - 1
        backsp: list[tuple[float, float, float, float, float]] = []

        for i in range(ils, 0, -1):
            # NOTE: the Fortran writes ``(float(i)/float(iils))`` and FLOAT is
            # single precision, so the dilution grid carries ~1e-7 relative
            # error.  Reproduced here because it propagates into every table
            # node; ``exact_grid=True`` uses the intended double-precision grid.
            frac = float(i) / float(npts)
            if not exact_grid:
                frac = float(np.float32(i) / np.float32(npts))
            zbda = frac / (1.0 + amb.humid)
            zw = zbda * amb.humid
            zg = 1.0 - zbda - zw
            enmix = zg * enthalpy
            zbda = zbda + zg * wa
            zg = zg * wc
            st = self.properties(zg, zbda, enmix, ifl=2)
            curnt = (st.yc, zg * st.rho, st.rho, enmix, st.temp)

            if i == ils:
                backsp = [curnt]
                continue

            last = rows[-1]
            err = 0.0
            for cand in backsp:
                slope = (last[1] - curnt[1]) / (last[0] - curnt[0])
                ccint = (cand[0] - curnt[0]) * slope + curnt[1]
                err = max(err, 2.0 * abs(cand[1] - ccint) / (abs(cand[1] + ccint) + zero))
                slope = (last[2] - curnt[2]) / (last[0] - curnt[0])
                rhoint = (cand[0] - curnt[0]) * slope + curnt[2]
                err = max(
                    err, 2.0 * abs(cand[2] - rhoint) / (abs(cand[2] + rhoint) + zero)
                )
                wccal = cand[1] / rhoint
                w1 = curnt[1] / curnt[2]
                w2 = last[1] / last[2]
                slope = (last[3] - curnt[3]) / (w2 - w1)
                entint = (wccal - w1) * slope + curnt[3]
                err = max(
                    err, 2.0 * abs(cand[3] - entint) / (abs(cand[3] + entint) + zero)
                )
                slope = (last[4] - curnt[4]) / (w2 - w1)
                temint = (wccal - w1) * slope + curnt[4]
                err = max(
                    err, 2.0 * abs(cand[4] - temint) / (abs(cand[4] + temint) + zero)
                )

            if err <= tcrit and len(backsp) < iback:
                backsp.append(curnt)
                continue

            rows.append(backsp[-1])
            if len(rows) >= maxrows:
                raise RuntimeError("adiabatic table overflow (TRAP 28)")
            backsp = [curnt]

        # final row: pure contaminant
        if wc == 1.0:
            rows.append((1.0, gas.rho, gas.rho, enthalpy, gas.temp))
        else:
            st = self.properties(wc, wa, enthalpy, ifl=2)
            rows.append((st.yc, wc * st.rho, st.rho, enthalpy, st.temp))

        arr = np.array(rows, dtype=float)
        self.table = AdiabaticTable(
            yc=arr[:, 0],
            cc=arr[:, 1],
            rho=arr[:, 2],
            h=arr[:, 3],
            t=arr[:, 4],
            humid=amb.humid,
            humsrc=amb.humsrc,
            gasmw=gas.mw,
        )
        return self.table

    # -- ADDHEAT -----------------------------------------------------------

    def add_heat(self, cc: float, dh: float, *, tol: float = 1.0e-3) -> MixtureState:
        """Port of ``ADDHEAT``: adiabatic mixing plus an enthalpy increment.

        ``dh`` is the extra enthalpy per unit mass acquired from the ground.
        Negative or zero ``dh`` returns the adiabatic state unchanged (the
        correlations are not valid for a colder surface).
        """
        if self.table is None:
            raise RuntimeError("add_heat requires the adiabatic table")
        amb = self.ambient
        st = self.table.from_concentration(cc)
        amt = st.temp
        if amb.isofl == 1 or amb.ihtfl == 0 or dh <= 0.0:
            return st

        enth = st.enthalpy + dh
        if enth > 0.0:
            temp = amb.tamb
        else:
            tmax = max(self.gas.temp, amb.tsurf, amb.tamb)
            ehi = self.enthalpy(st.wc, st.wa, tmax)
            if enth > ehi:
                temp = tmax
            else:
                residual = lambda t: enth - self.enthalpy(st.wc, st.wa, t)
                tmin = amt
                try:
                    temp = self._invert_enthalpy(
                        st.wc, st.wa, enth, tmin, tmax, tol
                    )
                except RootBracketError:
                    # The adiabatic mixing temperature stops bracketing the
                    # root once enough heat has been added, so ADDHEAT falls
                    # back on LIMIT to widen the interval before retrying.
                    tinc = (tmin + tmax) / 100.0
                    hi, lo = limit(residual, tmin, tinc, 2.0 * tmax, tmin / 2.0)
                    temp = self._invert_enthalpy(
                        st.wc, st.wa, enth, min(lo, hi), max(lo, hi), tol
                    )

        rho = self.density(st.wc, st.wa, st.wm, temp)
        if self.legacy_numerics:
            # ADDHEAT forms the heat capacity as dh/(T - T_adiabatic). Early in
            # the dense phase that denominator is a few hundredths of a kelvin
            # while ZBRENT resolves T to only 1e-3 K, so the quotient carries
            # percent-level noise -- which then feeds the ground heat flux and
            # integrates back into dh. It is the one place the port loses
            # accuracy to the original's own numerics.
            cp = max(dh / (temp - amt), CPA) if temp != amt else CPA
        elif temp != amt:
            # The same secant, but over an accurately located interval. With
            # the temperature resolved to 1e-12 K instead of 1e-3 K the
            # quotient is well conditioned however small the interval is, and
            # it is exactly the mean heat capacity over the heating the cloud
            # actually received.
            cp = max(dh / (temp - amt), CPA)
        else:
            cp = max(self.heat_capacity(st.wc, st.wa, amt), CPA)
        return MixtureState(
            wc=st.wc,
            wa=st.wa,
            yc=st.yc,
            ya=st.ya,
            wm=st.wm,
            temp=temp,
            rho=rho,
            cp=cp,
            enthalpy=enth,
        )
