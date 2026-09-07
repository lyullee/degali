"""Building the ``/GEN3/`` source vectors: a port of ``CRFG.FOR``.

``SRC1O`` already thins the blanket integration down to the points where
something changed by more than ``SRCOER``.  ``CRFG`` thins it *again*, and for
a different reason: what the downwind model consumes is not a list of states
but six functions of time,

===========  ========================================================
``radg``     blanket radius
``qstr``     take-up flux
``srcden``   density
``srcwc``    contaminant mass fraction
``srcwa``    dry-air mass fraction
``srcen``    enthalpy
===========  ========================================================

which it interpolates linearly.  So the question is not "did anything change"
but "can the points I keep reproduce the points I drop".  ``CRFG`` answers it
with a lookahead: it holds up to ``NTAB`` records in a buffer, and for each
candidate end point checks that a straight line from the last kept point to
the candidate reproduces *every* buffered record in between to within
``CRFGER``.  The moment one fails, the previous candidate is kept and the
buffer restarts.

This is why the table printed in the ``.scl`` listing has fewer rows than the
integration recorded: the listing is the thinned set, written as ``CRFG``
selects it. Comparing against ``.scl`` therefore tests both stages at once.

Two details are load-bearing:

* The first kept point has its time forced to zero, whatever the integration
  reported, so the vectors always start at the release instant.
* The relative errors are taken against the *buffered* value with a
  ``1e-20`` guard in the denominator for the five quantities that can be
  zero, but with no guard for radius and density, which never are.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MAXL, POUND

#: Number of fields in a record written by ``SRC1O``.
IOUSRC = 13

#: Buffer depth: ``ntab = ntab0/iousrc`` with ``ntab0 = 910`` in ``DEG1.for``.
NTAB = 910 // IOUSRC

#: Column indices used by ``CRFG``, 0-based (the Fortran's ``iti``, ``irg``,
#: ``iqs``, ``idn``, ``iwc``, ``iwa``, ``ien`` minus one).
I_TIME, I_RADG, I_QSTR, I_DEN, I_WC, I_WA, I_EN = 0, 1, 3, 7, 9, 10, 11

_ZERO = 1.0e-20


@dataclass
class SourceVectors:
    """``/GEN3/``: the secondary source as functions of time."""

    time: np.ndarray
    radg: np.ndarray  #: blanket radius, m
    qstr: np.ndarray  #: take-up flux, kg/(m**2 s)
    srcden: np.ndarray  #: density, kg/m**3
    srcwc: np.ndarray  #: contaminant mass fraction
    srcwa: np.ndarray  #: dry-air mass fraction
    srcen: np.ndarray  #: enthalpy, J/kg
    #: The rows as ``CRFG`` writes them to the listing, for comparison against
    #: the ``.scl`` file: time, radius, height, qstar, sz, yc, rho, T, Ri.
    listing: np.ndarray = None

    def __len__(self) -> int:
        return len(self.time)

    def at(self, t: float) -> dict:
        return {
            name: float(np.interp(t, self.time, getattr(self, name)))
            for name in ("radg", "qstr", "srcden", "srcwc", "srcwa", "srcen")
        }

    def _interp(self, column: np.ndarray, t: float) -> float:
        return float(np.interp(t, self.time, column))

    def radg_at(self, t: float) -> float:
        return self._interp(self.radg, t)

    def qstr_at(self, t: float) -> float:
        return self._interp(self.qstr, t)

    def srcden_at(self, t: float) -> float:
        return self._interp(self.srcden, t)

    def srcwc_at(self, t: float) -> float:
        return self._interp(self.srcwc, t)

    def srcwa_at(self, t: float) -> float:
        return self._interp(self.srcwa, t)

    def srcen_at(self, t: float) -> float:
        return self._interp(self.srcen, t)

    @classmethod
    def from_handoff(cls, gen3: np.ndarray) -> "SourceVectors":
        """Build from the ``(k, 7)`` block in a ``.TR2`` handoff file."""
        a = np.asarray(gen3, dtype=float)
        return cls(
            time=a[:, 0], radg=a[:, 1], qstr=a[:, 2], srcden=a[:, 3],
            srcwc=a[:, 4], srcwa=a[:, 5], srcen=a[:, 6],
        )


def build_source_vectors(
    records: np.ndarray, *, crfger: float = 0.008, ntab: int = NTAB
) -> SourceVectors:
    """Port of ``CRFG``.

    Parameters
    ----------
    records
        The raw ``SRC1O`` output, ``(n, 13)``.
    crfger
        ``RER``: the relative error a kept point must reproduce the dropped
        ones to.
    ntab
        Buffer depth.  Running out of buffer is not an error -- the Fortran
        logs a note and keeps the oldest candidate -- because a larger buffer
        would only have let it drop more points, not change the answer.
    """
    if len(records) < 2:
        raise ValueError("CRFG needs at least two records")

    rows = [np.asarray(r, dtype=float) for r in records]
    kept: list[np.ndarray] = []
    listing: list[np.ndarray] = []

    # The first vector point is the initial record with its time forced to 0.
    first = rows[0].copy()
    kept.append(first)
    listing.append(first)
    last = first

    cur = rows[1]
    idx = 2
    buffer: list[np.ndarray] = []

    def _err(candidate: np.ndarray, stored: np.ndarray, ratio: float) -> float:
        """Relative error of the linear interpolant at the stored point."""
        err = 0.0
        for col, guarded in (
            (I_RADG, False), (I_QSTR, True), (I_DEN, False),
            (I_WC, True), (I_WA, True), (I_EN, True),
        ):
            interp = (candidate[col] - last[col]) * ratio + last[col]
            denom = stored[col] + (_ZERO if guarded else 0.0)
            err = max(err, abs(interp - stored[col]) / denom)
        return err

    while True:
        buffer = []
        overflow = True
        for _ in range(2, ntab + 1):
            buffer.append(cur)
            if idx >= len(rows):
                cur_next = None
            else:
                cur_next = rows[idx]
                idx += 1
            if cur_next is None:
                # EOF (label 900): keep the final record and stop
                kept.append(cur)
                listing.append(cur)
                return _finish(kept, listing)

            timeslot = last[I_TIME]
            span = cur_next[I_TIME] - timeslot
            failed = False
            for stored in buffer:
                ratio = (stored[I_TIME] - timeslot) / span if span else 0.0
                if _err(cur_next, stored, ratio) > crfger:
                    failed = True
                    break
            cur = cur_next
            if failed:
                overflow = False
                break

        # label 150 (or the NTAB overflow note): keep the last record that
        # satisfied the criterion, which is the newest buffered one.
        keep = buffer[-1]
        kept.append(keep)
        listing.append(keep)
        last = keep
        if overflow:
            # buffer exhausted without a failure; the Fortran writes a note
            # and carries on from the same point
            continue


def _finish(kept: list[np.ndarray], listing: list[np.ndarray]) -> SourceVectors:
    arr = np.array(kept)
    if len(arr) + 1 > MAXL:
        raise RuntimeError("CRFG overflowed the GEN3 vectors (TRAP 5)")
    listing_arr = np.array(listing)
    return SourceVectors(
        time=arr[:, I_TIME],
        radg=arr[:, I_RADG],
        qstr=arr[:, I_QSTR],
        srcden=arr[:, I_DEN],
        srcwc=arr[:, I_WC],
        srcwa=arr[:, I_WA],
        srcen=arr[:, I_EN],
        listing=listing_arr[:, [0, 1, 2, 3, 4, 5, 7, 12, 8]],
    )
