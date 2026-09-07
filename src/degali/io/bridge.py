"""Handing a landed jet to the ground-level model: a port of ``DEGBRIDG``.

``JETPLU`` stops at touchdown and reports three numbers: how far downwind the
plume landed, the ground-level centreline concentration there, and the
half-width of the landed footprint.  ``DEGBRIDG`` turns those into a
ground-level DEGADIS case -- an equivalent area source that ``DEG1`` can start
from.

The equivalence
---------------
The landed plume is replaced by a circular source of radius equal to the
plume's half-width, releasing the same contaminant mass rate but *already
diluted* to the concentration the plume had when it landed.  So the source
composition is not the pure release material: it is the mixture at ``cc``,
read off the mixing line, with the corresponding temperature.  That is the
whole point of the bridge -- the ground-level model would otherwise restart
from pure contaminant and over-predict everything downwind.

The touchdown distance becomes ``OODIST``, an offset added to every reported
distance so the two stages share one coordinate.

Density model
-------------
The jet deck's ``NDEN`` flag selects between three treatments, and the bridge
carries the choice across:

``NDEN < 0``
    Two-point table: linear between ambient and release density, isothermal.
``NDEN == 0``
    No table; the ground-level model computes the adiabatic mixing line
    itself, with heat transfer if ``INDHT`` says so.
``NDEN > 0``
    The user's own density table, treated as isothermal.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core.constants import PI, RGAS, WMA, WMW
from ..core.thermo import (
    AdiabaticTable,
    AmbientConditions,
    GasProperties,
    LegacyBackend,
    Thermo,
)
from .inp import Case, SourceTable
from .jetdeck import Touchdown

#: ``DEGBRIDG`` writes a source that runs for this long when the jet deck
#: declares a continuous release: ``if(check4) tend = 60230.D0``.
CONTINUOUS_TEND = 60230.0


@dataclass
class JetInput:
    """A ``.IN`` jet deck -- what the user writes and ``DEGBRIDG`` reads.

    Distinct from the ``.INO`` deck :func:`~degali.io.jetdeck.read_ino`
    handles, which ``JETPLUIN`` derives from this one. ``DEGBRIDG`` goes back
    to the original because it needs fields ``JETPLUIN`` consumed rather than
    passed on -- ``INDVEL``, ``INDHT``, the raw ``NDEN`` flag.
    """

    titles: list[str]
    u0: float
    z0: float
    zr: float
    indvel: int
    istab: int
    rml: float
    tamb: float
    pamb: float
    relhum: float
    tsurf: float
    gasnam: str
    gasmw: float
    avtime: float
    temjet: float  #: jet release temperature, K
    gasulc: float
    gasllc: float
    gaszzc: float
    indht: int  #: heat-transfer prescription for the ground-level stage
    cpk: float
    cpp: float
    nden: int
    den: np.ndarray  #: (n, 3) user density table, when ``nden > 0``
    erate: float
    elejet: float
    diajet: float
    tend: float
    distmx: float


def read_in(path: str | Path) -> JetInput:
    """Read a ``.IN`` jet deck, following ``DEGBRIDG``'s read sequence.

    The file is human-written and carries trailing comments on every line, so
    unlike the machine-written decks it must be read line by line: only the
    leading numbers on each line are data.
    """
    text = Path(path).read_text(errors="replace").replace("\r\n", "\n")
    lines = text.split("\n")
    titles = [lines[i].rstrip() for i in range(4)]
    idx = 4

    def row(n: int) -> list[float]:
        nonlocal idx
        parts = lines[idx].split()
        idx += 1
        return [float(x) for x in parts[:n]]

    def text_row() -> str:
        nonlocal idx
        value = lines[idx].split()[0] if lines[idx].split() else ""
        idx += 1
        return value

    u0, z0 = row(2)
    (zr,) = row(1)
    indvel, istab, rml = row(3)
    tamb, pamb, relhum = row(3)
    (tsurf,) = row(1)
    gasnam = text_row()
    (gasmw,) = row(1)
    (avtime,) = row(1)
    (temjet,) = row(1)
    gasulc, gasllc, gaszzc = row(3)
    indht, cpk, cpp = row(3)
    (nden,) = row(1)
    nden = int(nden)
    den = np.array([row(3) for _ in range(nden)]) if nden > 0 else np.zeros((0, 3))
    (erate,) = row(1)
    elejet, diajet = row(2)
    (tend,) = row(1)
    (distmx,) = row(1)

    return JetInput(
        titles=titles, u0=u0, z0=z0, zr=zr, indvel=int(indvel),
        istab=int(istab), rml=rml, tamb=tamb, pamb=pamb, relhum=relhum,
        tsurf=tsurf, gasnam=gasnam, gasmw=gasmw, avtime=avtime,
        temjet=temjet, gasulc=gasulc, gasllc=gasllc, gaszzc=gaszzc,
        indht=int(indht), cpk=cpk, cpp=cpp, nden=nden, den=den,
        erate=erate, elejet=elejet, diajet=diajet, tend=tend, distmx=distmx,
    )


def bridge(touchdown: Touchdown, deck: JetInput) -> Case:
    """Port of ``DEGBRIDG``: build a ground-level case from a jet touchdown.

    Raises
    ------
    ValueError
        If the plume never landed. ``JETPLU`` signals that with a touchdown
        distance of zero, and there is no ground-level stage to set up.
    """
    if not touchdown.lands:
        raise ValueError(
            "the plume never reached the ground at a concentration of "
            "interest; there is no ground-level stage to bridge to"
        )

    backend = LegacyBackend()
    vp = backend.water_vapour_pressure(deck.tamb)
    sat = WMW / WMA * vp / (deck.pamb - vp)
    humid = deck.relhum / 100.0 * sat

    rhoa = (
        deck.pamb * (1.0 + humid) * WMW / (RGAS * (WMW / WMA + humid)) / deck.tamb
    )
    rhoe = deck.pamb * deck.gasmw / RGAS / deck.temjet

    # -- density model -----------------------------------------------------
    nden = max(deck.nden, -1)
    if nden < 0:
        isofl = 1
        table_rows = np.array([
            [0.0, 0.0, rhoa, 0.0, deck.tamb],
            [1.0, rhoe, rhoe, 0.0, deck.tamb],
        ])
    elif nden == 0:
        isofl = 0
        table_rows = None
    else:
        isofl = 1
        table_rows = np.column_stack([
            deck.den[:, 0], deck.den[:, 1], deck.den[:, 2],
            np.zeros(nden), np.full(nden, deck.tamb),
        ])

    gas = GasProperties(
        name=deck.gasnam, mw=deck.gasmw, temp=deck.temjet, rho=rhoe,
        cpk=deck.cpk, cpp=deck.cpp, ulc=deck.gasulc, llc=deck.gasllc,
        zzc=deck.gaszzc,
    )
    ambient = AmbientConditions(
        tamb=deck.tamb, pamb=deck.pamb, humid=humid, tsurf=deck.tsurf,
        isofl=isofl, ihtfl=deck.indht, iwtfl=0, htco=0.0, wtco=0.0,
    )

    # -- the diluted source composition -------------------------------------
    thermo = Thermo(gas=gas, ambient=ambient, backend=backend)
    thermo.reference_enthalpies()
    if table_rows is not None:
        thermo.table = AdiabaticTable(
            yc=table_rows[:, 0], cc=table_rows[:, 1], rho=table_rows[:, 2],
            h=table_rows[:, 3], t=table_rows[:, 4], humid=humid, humsrc=0.0,
            gasmw=deck.gasmw,
        )
        # the last table row is the release density, which is also what
        # DEGBRIDG writes as GASRHO
        gas.rho = float(table_rows[-1, 2])
    else:
        thermo.build_adiabatic_table(1.0, 0.0, thermo.hmrte)
    landed = thermo.table.from_concentration(touchdown.concentration)

    # -- the equivalent area source ------------------------------------------
    continuous = deck.tend <= 0.0
    tend = CONTINUOUS_TEND if continuous else deck.tend
    radius = touchdown.halfwidth
    times = np.array([0.0, tend, tend + 1.0, tend + 2.0])
    rates = np.array([deck.erate, deck.erate, 0.0, 0.0])
    radii = np.array([radius, radius, 0.0, 0.0])
    wc = np.full(4, landed.wc)
    temp = np.full(4, landed.temp)
    fracv = np.ones(4)

    source = SourceTable(
        time=times, rate=rates, radius=radii, wc=wc, temp=temp, fracv=fracv,
        enthalpy=np.zeros(4), rho=np.zeros(4),
    )
    for i in range(2):
        wa = (1.0 - source.wc[i]) / (1.0 + humid)
        h = thermo.enthalpy(source.wc[i], wa, source.temp[i])
        source.enthalpy[i] = h
        source.rho[i] = thermo.properties(
            source.wc[i], wa, h, ifl=-1, temp=source.temp[i]
        ).rho
    source.enthalpy[2:] = source.enthalpy[1]
    source.rho[2:] = source.rho[1]

    from ..core.atmosphere import stability_defaults

    stability = stability_defaults(deck.zr, deck.istab, deck.avtime)
    rml = deck.rml if deck.indvel == 2 else stability.rml

    return Case(
        titles=deck.titles,
        u0=deck.u0, z0=deck.z0, zr=deck.zr, istab=deck.istab,
        oodist=touchdown.distance, avtime=deck.avtime,
        indvel=deck.indvel, rml=rml,
        gas=gas, ambient=ambient, relhum=deck.relhum,
        # DEGBRIDG divides the level of concern by 1.02 rather than the 5 the
        # ground-level reader would apply, because the jet stage has already
        # carried the cloud most of the way down
        yclow=deck.gasllc / 1.02,
        gmass0=0.0, source=source, steady_state=True,
        density_table=table_rows, timestamp="", instantaneous=False,
        stability=stability,
    )
