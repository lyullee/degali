"""Reading the ``.TR2`` handoff file written by ``DEG1`` and read by ``DEG2``.

``TRANS`` writes it and ``STRT2`` reads it back.  It is the entire interface
between the source calculation and the downwind calculation: the primary
source table, the adiabatic mixing table, the ``/GEN3/`` source vectors, and
every scalar the second stage needs.  There are no labels -- the file is a
sequence of unformatted numbers whose meaning is fixed by position, so it can
only be read by walking ``TRANS`` in order.

Reading it lets ``degali`` start a downwind calculation from a source
calculation the *Fortran* performed, which is how the two halves of the port
are validated independently of each other.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Handoff:
    """Everything ``DEG1`` passes to ``DEG2``."""

    titles: list[str]
    source: np.ndarray  #: (n, 8) PTIME, ET, R1T, PWC, PTEMP, PFRACV, PENTH, PRHO
    den: np.ndarray  #: (m, 5) the adiabatic mixing table
    gen3: np.ndarray  #: (k, 7) time, radg, qstr, srcden, srcwc, srcwa, srcen
    scalars: dict = field(default_factory=dict)

    def __getitem__(self, key):
        return self.scalars[key]


def _tokens(lines, start):
    i = start
    buf: list[str] = []

    def take(n):
        nonlocal i, buf
        out = []
        while len(out) < n:
            if not buf:
                buf = lines[i].replace(",", " ").split()
                i += 1
                continue
            out.append(buf.pop(0))
        return out

    def line():
        nonlocal i, buf
        buf = []
        i += 1
        return lines[i - 1]

    def pos():
        return i

    return take, line, pos


def read_tr2(path: str | Path) -> Handoff:
    """Read a ``.TR2`` file, following ``TRANS``/``STRT2`` field for field."""
    lines = Path(path).read_text(errors="replace").replace("\r\n", "\n").split("\n")
    titles = [lines[k].rstrip() for k in range(4)]
    take, line, pos = _tokens(lines, 4)

    n = int(take(1)[0])
    source = np.array([[float(x) for x in take(8)] for _ in range(n)])
    m = int(take(1)[0])
    den = np.array([[float(x) for x in take(5)] for _ in range(m)])
    k = int(take(1)[0])
    gen3 = np.array([[float(x) for x in take(7)] for _ in range(k)])

    tinp_tsrc = line()
    tobs_tsrt = line()
    s: dict = {"tinp": tinp_tsrc[:24].strip(), "tsrc": tinp_tsrc[24:].strip()}

    s["oodist"], s["avtime"] = (float(x) for x in take(2))
    for name, value in zip(
        ("u0", "z0", "zr", "rml", "ustar"), (float(x) for x in take(5))
    ):
        s[name] = value
    for name, value in zip(
        ("vkc", "gg", "rhoe", "rhoa", "deltay"), (float(x) for x in take(5))
    ):
        s[name] = value
    s["betay"], s["gammaf"], s["yclow"] = (float(x) for x in take(3))
    for name, value in zip(
        ("rm", "szm", "emax", "rmax", "tsc1"), (float(x) for x in take(5))
    ):
        s[name] = value
    s["aleph"], s["tend"] = (float(x) for x in take(2))

    flags = take(6)
    for name, value in zip(
        ("check1", "check2", "again", "check3", "check4", "check5"), flags
    ):
        s[name] = value.strip().upper().startswith("T")

    s["alpha"] = float(take(1)[0])
    s["gasnam"] = line().strip()
    s["gasmw"], s["gastem"], s["gasrho"] = (float(x) for x in take(3))
    s["gascpk"], s["gascpp"] = (float(x) for x in take(2))
    s["gasulc"], s["gasllc"], s["gaszzc"] = (float(x) for x in take(3))
    s["istab"] = int(take(1)[0])
    s["tamb"], s["pamb"], s["humid"] = (float(x) for x in take(3))
    s["isofl"], s["tsurf"] = int(take(1)[0]), float(take(1)[0])
    s["ihtfl"], s["htco"] = int(take(1)[0]), float(take(1)[0])
    s["iwtfl"], s["wtco"] = int(take(1)[0]), float(take(1)[0])
    s["sigxco"], s["sigxp"], s["sigxmd"] = (float(x) for x in take(3))

    if s["check4"]:
        s["ess"], s["slen"], s["swid"] = (float(x) for x in take(3))
        s["outcc"], s["outsz"], s["outb"], s["outl"] = (float(x) for x in take(4))
        s["swcl"], s["swal"], s["senl"], s["srhl"] = (float(x) for x in take(4))

    return Handoff(titles=titles, source=source, den=den, gen3=gen3, scalars=s)
