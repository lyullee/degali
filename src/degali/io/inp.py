"""Reading and writing DEGADIS input decks.

The ``RUN_NAME.INP`` file is the deck ``DEGINP`` (or ``DEGBRIDG``, for jet
releases that reach the ground) hands to ``DEG1``.  It is free-format: values
are whitespace- or comma-separated and the *order* carries all the meaning.
``IO.FOR`` reads it with a fixed sequence of list-directed ``READ`` statements,
so a deck is only interpretable by walking that sequence.

This module does the walk once and produces a :class:`Case`, after which
everything is named.  :func:`write_inp` goes the other way, so a case built in
Python can be handed to the original Fortran for comparison.

Deck layout (from ``IO.FOR`` and the DEGADIS readme)::

    TITLE1 .. TITLE4                 4 lines, 80 characters each
    U0, Z0, ZR                       wind speed (m/s) at height (m); roughness (m)
    ISTAB                            Pasquill-Gifford class, 1=A .. 6=F
    OODIST, AVTIME                   jet touchdown distance (m); averaging time (s)
    INDVEL, RML                      2 = use the supplied Monin-Obukhov length
    TAMB, PAMB, HUMID, RELHUM        K; atm; kg/kg BDA; per cent
    ISOFL, TSURF                     1 = density supplied as a table; surface T (K)
    IHTFL, HTCO                      heat-transfer prescription; fixed coefficient
    IWTFL, WTCO                      water-transfer prescription; fixed coefficient
    GASNAM                           3-character label
    GASMW, GASTEM, GASRHO            kg/kmol; K; kg/m**3
    CPK, CPP                         heat-capacity correlation constants
    GASULC, GASLLC, GASZZC           upper and lower contour mole fractions; height (m)
    [ NP, then NP rows of 5 ]        density table, only when ISOFL != 0
    YCLOW                            mole fraction at which to stop
    GMASS0                           initial contaminant mass over the source (kg)
    NT                               number of source description rows
    NT rows of:
        PTIME, ET, R1T, PWC, PTEMP, PFRACV
    CHECK4                           T for a steady-state simulation
    TINP                             timestamp written by DEGINP

Quirks reproduced
-----------------
* ``GASULC``/``GASLLC`` are floored at 2e-12 and 1e-12.
* ``CPP == 0`` is a shorthand for a constant heat capacity: the reader then
  sets ``CPP = 1`` and ``CPK = GASMW*CPK - 3.33e4``, inverting the ``CPC``
  correlation so that a user who types "1800 J/kg/K" gets it.
* The last two source rows are overwritten with the properties of row
  ``NT-2``; only their times and (zero) rates are honoured.  ``TEND`` is
  ``PTIME(NT-2)``, not the last time in the table.
* For a transient run (``CHECK4`` false) ``YCLOW`` is divided by 5, so the
  integration is carried a factor of five below the requested cut-off.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core.atmosphere import StabilityDefaults, absolute_humidity, stability_defaults
from ..core.constants import PI, POUND
from ..core.thermo import (
    AmbientConditions,
    CoolPropBackend,
    GasProperties,
    LegacyBackend,
    Thermo,
    ThermoBackend,
)


@dataclass
class SourceTable:
    """The primary source description, ``/GEN1/`` in the Fortran.

    Every column is a function of ``time`` and is interpolated linearly by
    ``AFGEN2``.  ``enthalpy`` and ``rho`` are not read from the deck; they are
    derived from ``wc`` and ``temp`` when the case is loaded.
    """

    time: np.ndarray  #: PTIME, s
    rate: np.ndarray  #: ET, contaminant mass rate, kg/s
    radius: np.ndarray  #: R1T, source radius, m
    wc: np.ndarray  #: PWC, contaminant mass fraction
    temp: np.ndarray  #: PTEMP, K
    fracv: np.ndarray  #: PFRACV, liquid mass fraction
    enthalpy: np.ndarray = field(default_factory=lambda: np.zeros(0))  #: PENTH, J/kg
    rho: np.ndarray = field(default_factory=lambda: np.zeros(0))  #: PRHO, kg/m**3

    def __len__(self) -> int:
        return len(self.time)

    @property
    def tend(self) -> float:
        """``TEND``: the last time at which the source is still running."""
        return float(self.time[-3])

    def _at(self, column: np.ndarray, t: float) -> float:
        return float(np.interp(t, self.time, column))

    def rate_at(self, t: float) -> float:
        return self._at(self.rate, t)

    def radius_at(self, t: float) -> float:
        return self._at(self.radius, t)

    def wc_at(self, t: float) -> float:
        return self._at(self.wc, t)

    def enthalpy_at(self, t: float) -> float:
        return self._at(self.enthalpy, t)

    def rho_at(self, t: float) -> float:
        return self._at(self.rho, t)


@dataclass
class Case:
    """A fully specified DEGADIS simulation."""

    titles: list[str]
    u0: float  #: wind speed, m/s
    z0: float  #: reference height for u0, m
    zr: float  #: surface roughness length, m
    istab: int  #: Pasquill-Gifford class, 1..6
    oodist: float  #: jet touchdown distance, m (0 for ground-level releases)
    avtime: float  #: averaging time, s
    indvel: int  #: 2 = RML supplied by the user
    rml: float  #: Monin-Obukhov length, m (0 = infinite)
    gas: GasProperties
    ambient: AmbientConditions
    relhum: float  #: relative humidity, per cent
    yclow: float  #: mole fraction at which to stop the integration
    gmass0: float  #: initial contaminant mass over the source, kg
    source: SourceTable
    steady_state: bool  #: CHECK4
    density_table: np.ndarray | None = None  #: DEN, only when ISOFL != 0
    timestamp: str = ""
    #: True for an "HSE type" spill: no initial rate but a finite initial mass.
    instantaneous: bool = False
    stability: StabilityDefaults | None = None

    @property
    def tend(self) -> float:
        return self.source.tend

    def make_thermo(
        self,
        backend: ThermoBackend | str | None = None,
        *,
        legacy_numerics: bool | None = None,
    ) -> Thermo:
        """Build the thermodynamic evaluator and its adiabatic mixing table.

        ``backend`` may be an object, or the string ``"legacy"`` or
        ``"coolprop"``. Choosing ``"coolprop"`` also resolves the deck's
        contaminant to a CoolProp fluid by label and molecular weight, so that
        its heat capacity and density come from an equation of state rather
        than the deck's two fitted constants; when no fluid matches, water and
        air still come from the EOS and the contaminant keeps its
        correlations.

        ``legacy_numerics`` defaults to ``True`` for the legacy backend and
        ``False`` for CoolProp: someone who has asked for real properties has
        not asked for the 1989 root-finder tolerances as well.
        """
        if isinstance(backend, str):
            from ..core.fluids import resolve
            from ..core.thermo import make_backend

            if backend == "coolprop":
                match = resolve(self.gas.name, self.gas.mw)
                self.gas.coolprop_name = match.fluid
                backend = make_backend("coolprop", contaminant=match.fluid)
            else:
                backend = make_backend(backend)
        if legacy_numerics is None:
            legacy_numerics = not isinstance(backend, CoolPropBackend)
        th = Thermo(
            gas=self.gas,
            ambient=self.ambient,
            backend=backend or LegacyBackend(),
            legacy_numerics=legacy_numerics,
        )
        th.reference_enthalpies()
        if self.ambient.isofl == 1:
            if self.density_table is None:
                raise ValueError("ISOFL=1 requires a density table in the deck")
            from ..core.thermo import AdiabaticTable

            d = self.density_table
            th.table = AdiabaticTable(
                yc=d[:, 0], cc=d[:, 1], rho=d[:, 2], h=d[:, 3], t=d[:, 4],
                humid=self.ambient.humid, humsrc=self.ambient.humsrc,
                gasmw=self.gas.mw,
            )
        else:
            th.build_adiabatic_table(1.0, 0.0, th.hmrte)
        return th


class _Tokens:
    """List-directed read over a free-format deck.

    Fortran's list-directed ``READ`` treats commas and whitespace alike and
    keeps consuming lines until the requested items are satisfied, so the
    reader cannot be line-oriented.
    """

    def __init__(self, lines: list[str], start: int):
        self._lines = lines
        self._i = start
        self._buf: list[str] = []

    def take(self, n: int) -> list[str]:
        out: list[str] = []
        while len(out) < n:
            if not self._buf:
                if self._i >= len(self._lines):
                    raise EOFError("unexpected end of INP deck")
                self._buf = self._lines[self._i].replace(",", " ").split()
                self._i += 1
                continue
            out.append(self._buf.pop(0))
        return out

    def floats(self, n: int) -> list[float]:
        return [float(t) for t in self.take(n)]

    def ints(self, n: int) -> list[int]:
        return [int(float(t)) for t in self.take(n)]

    def line(self) -> str:
        self._buf = []
        if self._i >= len(self._lines):
            raise EOFError("unexpected end of INP deck")
        text = self._lines[self._i]
        self._i += 1
        return text


def read_inp(path: str | Path, *, backend: ThermoBackend | None = None) -> Case:
    """Read a ``RUN_NAME.INP`` deck.

    Port of ``IO.FOR``, including its derived quantities: the stability
    defaults from ``ATMDEF``, the humidity reconciliation, the ``CPP == 0``
    shorthand, and the source enthalpy and density columns.
    """
    raw = Path(path).read_text(errors="replace").replace("\r\n", "\n")
    lines = raw.split("\n")
    be = backend or LegacyBackend()

    titles = [lines[i].rstrip() for i in range(4)]
    tok = _Tokens(lines, 4)

    u0, z0, zr = tok.floats(3)
    (istab,) = tok.ints(1)
    oodist, avtime = tok.floats(2)
    stab = stability_defaults(zr, istab, avtime)

    indvel_s, rml_s = tok.take(2)
    indvel = int(float(indvel_s))
    rml = float(rml_s) if indvel == 2 else stab.rml

    tamb, pamb, humid, relhum = tok.floats(4)
    humid, relhum = absolute_humidity(
        tamb, pamb, be.water_vapour_pressure, humid=humid, relhum=relhum
    )

    isofl, tsurf = tok.take(2)
    isofl = int(float(isofl))
    tsurf = float(tsurf)
    ihtfl_s, htco_s = tok.take(2)
    iwtfl_s, wtco_s = tok.take(2)

    gasnam = tok.line().strip()

    gasmw, gastem, gasrho = tok.floats(3)
    cpk, cpp = tok.floats(2)
    if cpp == 0.0:
        # shorthand: CPK is a constant mass heat capacity in J/(kg K)
        cpp = 1.0
        cpk = gasmw * cpk - 3.33e4

    gasulc, gasllc, gaszzc = tok.floats(3)
    gasulc = max(gasulc, 2.0e-12)
    gasllc = max(gasllc, 1.0e-12)

    density_table = None
    if isofl != 0:
        (np_rows,) = tok.ints(1)
        density_table = np.array(
            [tok.floats(5) for _ in range(np_rows)], dtype=float
        )

    (yclow,) = tok.floats(1)
    yclow = min(max(yclow, 1.0e-12), gasllc)
    (gmass0,) = tok.floats(1)

    (nt,) = tok.ints(1)
    rows = np.array([tok.floats(6) for _ in range(nt)], dtype=float)

    check4 = tok.take(1)[0].strip().upper().startswith("T")
    try:
        timestamp = tok.line().strip()
    except EOFError:
        timestamp = ""

    if not check4:
        yclow = yclow / 5.0  # transients are carried a factor of 5 further

    gas = GasProperties(
        name=gasnam, mw=gasmw, temp=gastem, rho=gasrho, cpk=cpk, cpp=cpp,
        ulc=gasulc, llc=gasllc, zzc=gaszzc,
    )
    ambient = AmbientConditions(
        tamb=tamb, pamb=pamb, humid=humid, tsurf=tsurf, isofl=isofl,
        ihtfl=int(float(ihtfl_s)), iwtfl=int(float(iwtfl_s)),
        htco=float(htco_s), wtco=float(wtco_s), humsrc=0.0,
    )

    source = SourceTable(
        time=rows[:, 0], rate=rows[:, 1], radius=rows[:, 2],
        wc=rows[:, 3], temp=rows[:, 4], fracv=rows[:, 5],
        enthalpy=np.zeros(nt), rho=np.zeros(nt),
    )

    # Derived source columns (IO.FOR labels 111 and 112).  Only the first
    # NT-2 rows are evaluated; the last two inherit row NT-3's properties,
    # because their rates are zero and only their times matter.
    th = Thermo(gas=gas, ambient=ambient, backend=be)
    if isofl == 1 and density_table is not None:
        # TPROP routes every property call through the supplied table when
        # ISOFL is set, so it has to be installed before the source columns
        # are derived
        from ..core.thermo import AdiabaticTable

        d = density_table
        th.table = AdiabaticTable(
            yc=d[:, 0], cc=d[:, 1], rho=d[:, 2], h=d[:, 3], t=d[:, 4],
            humid=humid, humsrc=0.0, gasmw=gasmw,
        )
    for i in range(nt - 2):
        wa = (1.0 - source.wc[i]) / (1.0 + humid)
        h = th.enthalpy(source.wc[i], wa, source.temp[i])
        source.enthalpy[i] = h
        source.rho[i] = th.properties(
            source.wc[i], wa, h, ifl=-1, temp=source.temp[i]
        ).rho
    for i in (nt - 2, nt - 1):
        source.wc[i] = source.wc[nt - 3]
        source.temp[i] = source.temp[nt - 3]
        source.fracv[i] = source.fracv[nt - 3]
        source.enthalpy[i] = source.enthalpy[nt - 3]
        source.rho[i] = source.rho[nt - 3]

    return Case(
        titles=titles, u0=u0, z0=z0, zr=zr, istab=istab, oodist=oodist,
        avtime=avtime, indvel=indvel, rml=rml, gas=gas, ambient=ambient,
        relhum=relhum, yclow=yclow, gmass0=gmass0, source=source,
        steady_state=check4, density_table=density_table, timestamp=timestamp,
        instantaneous=bool(source.rate[0] == 0.0 and gmass0 != 0.0),
        stability=stab,
    )


def write_inp(case: Case, path: str | Path) -> None:
    """Write a deck the original Fortran can read.

    Round-tripping through this and :func:`read_inp` is how a case constructed
    in Python is handed to the reference implementation for comparison.
    """
    g, a, s = case.gas, case.ambient, case.source
    out: list[str] = [f"{t:<80.80}" for t in (case.titles + [""] * 4)[:4]]
    e = lambda *v: out.append(" ".join(f"{x:.15E}" if isinstance(x, float) else str(x) for x in v))

    e(case.u0, case.z0, case.zr)
    out.append(f"{case.istab:12d}")
    e(case.oodist, case.avtime)
    out.append(f"{case.indvel:12d} {case.rml:.15E}")
    e(a.tamb, a.pamb, a.humid, case.relhum)
    out.append(f"{a.isofl:12d} {a.tsurf:.15E}")
    out.append(f"{a.ihtfl:12d} {a.htco:.15E}")
    out.append(f"{a.iwtfl:12d} {a.wtco:.15E}")
    out.append(g.name)
    e(g.mw, g.temp, g.rho)
    e(g.cpk, g.cpp)
    e(g.ulc, g.llc, g.zzc)
    if a.isofl != 0 and case.density_table is not None:
        out.append(f"{len(case.density_table):12d}")
        for row in case.density_table:
            e(*[float(x) for x in row])
    # YCLOW is divided by 5 on read for transients; undo that on write
    e(case.yclow * (1.0 if case.steady_state else 5.0))
    e(case.gmass0)
    out.append(f"{len(s):12d}")
    for i in range(len(s)):
        e(s.time[i], s.rate[i], s.radius[i], s.wc[i], s.temp[i], s.fracv[i])
    out.append(" T" if case.steady_state else " F")
    out.append(f" {case.timestamp}")
    Path(path).write_text("\n".join(out) + "\n")
