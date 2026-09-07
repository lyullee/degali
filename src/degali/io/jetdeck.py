"""Reading jet/plume input decks (``.INO``) and touchdown files (``.IND``).

``JETPLUIN`` writes the ``.INO`` deck that ``JETPLU`` reads; ``JETPLU`` writes
the ``.IND`` file that ``DEGBRIDG`` reads to build a ground-level DEGADIS
deck. Both are free-format and positional, like the other DEGADIS files, so
they can only be read by walking the original's ``READ`` sequence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class JetDeck:
    """A ``.INO`` deck, in the order ``JETPLU`` consumes it."""

    titles: list[str]
    u0: float
    z0: float
    zr: float
    istab: int
    rml: float
    ustar: float
    tamb: float
    pamb: float
    humid: float
    relhum: float
    tsurf: float
    avtime: float
    deltay: float
    betay: float
    deltaz: float
    betaz: float
    gammaz: float
    #: The .INO deck carries no contaminant label -- JETPLUIN consumes it and
    #: does not pass it on -- so a CoolProp fluid can only be resolved from
    #: the molecular weight here. The .IN deck read by
    #: :func:`degali.io.bridge.read_in` does carry one.
    gasmw: float
    gastem: float
    gasulc: float
    gasllc: float
    gaszzc: float
    gascpk: float
    gascpp: float
    nden: int
    den: np.ndarray  #: (n, 5) density table, empty when ``nden0 == 0``
    isofl: int
    erate: float  #: contaminant mass release rate, kg/s
    elejet: float  #: orifice elevation, m
    diajet: float  #: orifice diameter, m
    alfa1: float
    alfa2: float
    distmx: float  #: maximum integration step, m
    tend: float  #: release duration, s (0 for a continuous release)

    @property
    def yclow(self) -> float:
        """Mole fraction at which the jet integration gives up.

        A quarter of the lower level of concern for a finite release, the
        level itself for a continuous one.
        """
        return self.gasllc / 4.0 if self.tend > 0.0 else self.gasllc


def read_ino(path: str | Path) -> JetDeck:
    """Read a ``.INO`` jet deck."""
    lines = Path(path).read_text(errors="replace").replace("\r\n", "\n").split("\n")
    titles = [lines[i].rstrip() for i in range(4)]
    toks: list[str] = []
    for line in lines[4:]:
        toks.extend(line.split())
    pos = 0

    def take(n: int) -> list[float]:
        nonlocal pos
        out = [float(t) for t in toks[pos:pos + n]]
        pos += n
        return out

    u0, z0 = take(2)
    (zr,) = take(1)
    istab, rml, ustar = take(3)
    tamb, pamb, humid, relhum, tsurf = take(5)
    avtime, deltay, betay = take(3)
    deltaz, betaz, gammaz = take(3)
    (gasmw,) = take(1)
    (gastem,) = take(1)
    gasulc, gasllc, gaszzc = take(3)
    gascpk, gascpp = take(2)
    nden, nden0 = (int(x) for x in take(2))
    den = np.array([take(5) for _ in range(nden0)]) if nden0 > 0 else np.zeros((0, 5))
    (isofl,) = take(1)
    (erate,) = take(1)
    elejet, diajet = take(2)
    alfa1, alfa2 = take(2)
    distmx, tend = take(2)

    return JetDeck(
        titles=titles, u0=u0, z0=z0, zr=zr, istab=int(istab), rml=rml,
        ustar=ustar, tamb=tamb, pamb=pamb, humid=humid, relhum=relhum,
        tsurf=tsurf, avtime=avtime, deltay=deltay, betay=betay,
        deltaz=deltaz, betaz=betaz, gammaz=gammaz, gasmw=gasmw,
        gastem=gastem, gasulc=gasulc, gasllc=gasllc, gaszzc=gaszzc,
        gascpk=gascpk, gascpp=gascpp, nden=int(nden), den=den,
        isofl=int(isofl), erate=erate, elejet=elejet, diajet=diajet,
        alfa1=alfa1, alfa2=alfa2, distmx=distmx, tend=tend,
    )


@dataclass
class Touchdown:
    """A ``.IND`` file: what ``JETPLU`` hands to ``DEGBRIDG``.

    ``distance == 0`` is the signal that the plume never reached the ground at
    a concentration of interest, and that no ground-level run should follow.
    """

    distance: float  #: downwind distance at touchdown, m
    concentration: float  #: ground-level centreline concentration, kg/m**3
    halfwidth: float  #: ``delta * sigma_y`` at touchdown, m

    @property
    def lands(self) -> bool:
        return self.distance > 0.0


def read_ind(path: str | Path) -> Touchdown:
    """Read a ``.IND`` touchdown file."""
    values = [float(x) for x in Path(path).read_text().split()[:3]]
    return Touchdown(*values)


def write_ind(td: Touchdown, path: str | Path) -> None:
    """Write a ``.IND`` file the original ``DEGBRIDG`` can read."""
    Path(path).write_text(
        f" {td.distance!r} {td.concentration!r} {td.halfwidth!r}\n"
    )
