"""Exact port of ``RKGST.FOR`` -- the integrator behind every DEGADIS module.

``RKGST`` is a fourth-order Runge-Kutta-Gill scheme with Richardson step
doubling.  DEGADIS uses it for the source blanket (``DEG1``), the downwind
integrations (``DEG2``, ``DEG2S``), the observer time histories (``DEG3``,
``DEG4``), the jet/plume trajectory (``JETPLU``) and two internal quadratures
(``ALPH``, ``SZF``).

Why keep it
-----------
:func:`scipy.integrate.solve_ivp` is better in every way that matters for new
work -- dense output, event handling, stiff solvers, tighter error control.
But three properties of ``RKGST`` are load-bearing for *reproducing* DEGADIS:

1. The error test is **relative and per-component with weights**,
   ``max_i w_i * (2/15) * |y2_i - y1_i| / |y2_i + y1_i|``, not the L2 norm of
   an absolute-plus-relative tolerance that SciPy uses.
2. Output points are emitted from inside the stepper via the ``OUTP``
   callback, and the callback is what decides when to *stop* (by setting
   ``PRMT(5)``) and what to record.  Several DEGADIS routines do real work in
   ``OUTP``, including mutating ``PRMT`` entries the derivative function then
   reads back.
3. The step-size controller mixes the original halving/doubling ladder with a
   Press-style ``0.9 * h * (tol/err)**0.25`` expansion, clipped at
   ``PRMT(5)``.  The sequence of steps it takes -- and therefore the exact
   points at which output is recorded -- is not reproducible any other way.

The port is deliberately literal, including the ``AUX`` scratch layout and the
``ISTEP``/``IHLF`` bookkeeping, so it can be read side by side with the
Fortran.  Modern integration lives in the ``solve_ivp``-based paths; this is
the reference oracle they are checked against.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from collections.abc import MutableSequence
from typing import Callable, Sequence

import numpy as np

#: ``IHLF`` return codes, as documented in the Fortran header.
IHLF_OK = 0
IHLF_TOO_MANY_BISECTIONS = 11
IHLF_ZERO_STEP = 12
IHLF_WRONG_SIGN = 13

# Runge-Kutta-Gill coefficients (label 20 of RKGST).
_A = (0.5, 1.0 - math.sqrt(0.5), 1.0 + math.sqrt(0.5), 1.0 / 6.0)
_B = (2.0, 1.0, 1.0, 2.0)
_C = (0.5, _A[1], _A[2], 0.5)

ERRSET = 1.0


class StopIntegration(Exception):
    """Raised by an output callback to set ``PRMT(5)`` and stop, as the
    Fortran does by assigning a non-zero value."""


@dataclass
class RKGSTResult:
    x: float
    y: np.ndarray
    dery: np.ndarray
    ihlf: int
    nfev: int
    nout: int


class Control(MutableSequence):
    """Named access to ``PRMT``, the integrator's control array.

    ``RKGST`` takes a single array carrying, in the first five slots, the
    integration bounds, the initial step, the error tolerance and a stop flag
    -- and, beyond them, whatever the caller's derivative routine wants to
    stash.  Reading ``prmt[9]`` at a call site tells nobody anything.

    This wraps the same list.  It *is* a list as far as the integrator and the
    ported derivative routines are concerned, so the numerics are untouched
    and every index the Fortran uses still works; the names are for the
    reader.

    ==========  =====  ====================================================
    ``lower``     0    start of the integration range
    ``upper``     1    end of it
    ``step``      2    initial step; ``RKGST`` doubles it before the first
    ``tolerance`` 3    error bound per step
    ``stop``      4    set non-zero by ``OUTP`` to terminate; also carries
                       the maximum step on entry
    ==========  =====  ====================================================
    """

    __slots__ = ("_values",)

    def __init__(self, values):
        self._values = list(values)

    # -- the named five ----------------------------------------------------

    @property
    def lower(self) -> float:
        return self._values[0]

    @property
    def upper(self) -> float:
        return self._values[1]

    @property
    def step(self) -> float:
        return self._values[2]

    @property
    def tolerance(self) -> float:
        return self._values[3]

    @property
    def stop(self) -> float:
        return self._values[4]

    @stop.setter
    def stop(self, value: float) -> None:
        self._values[4] = value

    def halt(self) -> None:
        """Terminate at the end of this step, as ``PRMT(5) = 1`` does."""
        self._values[4] = 1.0

    # -- and the list it still is -----------------------------------------

    def __getitem__(self, i):
        return self._values[i]

    def __setitem__(self, i, v):
        self._values[i] = v

    def __delitem__(self, i):
        del self._values[i]

    def __len__(self) -> int:
        return len(self._values)

    def insert(self, i, v) -> None:
        self._values.insert(i, v)

    def __repr__(self) -> str:
        return (
            f"Control(lower={self.lower:g}, upper={self.upper:g}, "
            f"step={self.step:g}, tolerance={self.tolerance:g}, "
            f"stop={self.stop:g}, +{len(self._values) - 5} caller slots)"
        )


def rkgst(
    fct: Callable[[float, np.ndarray, np.ndarray, list], None],
    outp: Callable[[float, np.ndarray, np.ndarray, int, int, list], None],
    prmt: list,
    y: Sequence[float],
    dery: Sequence[float],
    ndim: int | None = None,
) -> RKGSTResult:
    """Integrate ``dy/dx = f(x, y)`` exactly as ``RKGST`` does.

    Parameters
    ----------
    fct
        Derivative routine ``fct(x, y, dery, prmt)``.  It writes into ``dery``
        in place and may read/write ``prmt`` beyond index 4, exactly like the
        Fortran ``FCT``.
    outp
        Output routine ``outp(x, y, dery, ihlf, ndim, prmt)``.  Set
        ``prmt[4]`` non-zero (or raise :class:`StopIntegration`) to terminate.
    prmt
        Mutable list, 1-based in the Fortran and 0-based here:

        * ``prmt[0]`` lower bound of the interval
        * ``prmt[1]`` upper bound
        * ``prmt[2]`` initial increment
        * ``prmt[3]`` upper error bound
        * ``prmt[4]`` maximum step size on input, termination flag on output
        * ``prmt[5:]`` free for communication between ``fct`` and ``outp``
    y
        Initial values; returned updated.
    dery
        Per-component error *weights* on input (centred on one), derivatives
        on output.  ``error_i = prmt[3] / weight_i``.
    """
    y = np.array(y, dtype=float)
    dery = np.array(dery, dtype=float)
    if ndim is None:
        ndim = len(y)

    aux = np.zeros((8, ndim))
    nfev = 0
    nout = 0

    def _fct(xx, yy, dd):
        nonlocal nfev
        nfev += 1
        fct(xx, yy, dd, prmt)

    def _outp(xx, yy, dd, ih):
        nonlocal nout
        nout += 1
        try:
            outp(xx, yy, dd, ih, ndim, prmt)
        except StopIntegration:
            prmt[4] = 1.0

    aux[7, :] = 2.0 / 15.0 * dery[:ndim]
    x = prmt[0]
    xend = prmt[1]
    h = prmt[2]
    stpmin = abs(h / 1024.0)
    stpmax = abs(prmt[4])
    prmt[4] = 0.0
    _fct(x, y, dery)

    # label "ERROR TEST": IF(H*(XEND-X)) 380, 370, 20
    test = h * (xend - x)
    if test == 0.0:
        ihlf = IHLF_ZERO_STEP
        _fct(x, y, dery)
        _outp(x, y, dery, ihlf)
        return RKGSTResult(x, y, dery, ihlf, nfev, nout)
    if test < 0.0:
        ihlf = IHLF_WRONG_SIGN
        _fct(x, y, dery)
        _outp(x, y, dery, ihlf, )
        return RKGSTResult(x, y, dery, ihlf, nfev, nout)

    aux[0, :] = y[:ndim]
    aux[1, :] = dery[:ndim]
    aux[2, :] = 0.0
    aux[5, :] = 0.0
    irec = 0
    h = h + h  # the outer step is twice the specified step
    ihlf = -1
    istep = 0
    iend = 0
    delt = 0.0

    def _finish(code: int) -> RKGSTResult:
        _fct(x, y, dery)
        _outp(x, y, dery, code)
        return RKGSTResult(x, y, dery, code, nfev, nout)

    while True:  # label 40 -- start of a (double) Runge-Kutta step
        probe = (x + h - xend) * h
        if probe > 0.0:
            h = xend - x
            iend = 1
        elif probe == 0.0:
            iend = 1

        # label 70 -- record the initial values of this step
        _fct(x, y, dery)
        _outp(x, y, dery, irec)
        if prmt[4] != 0.0:
            return RKGSTResult(x, y, dery, IHLF_OK, nfev, nout)

        itest = 0
        while True:  # label 90
            istep += 1

            # ---- innermost Runge-Kutta-Gill loop (labels 100-140) --------
            for j in range(4):
                aj, bj, cj = _A[j], _B[j], _C[j]
                r1 = h * dery[:ndim]
                r2 = aj * (r1 - bj * aux[5, :])
                y[:ndim] += r2
                r2 = r2 + r2 + r2
                aux[5, :] += r2 - cj * r1
                if j < 3:
                    if j != 1:  # IF(J-3) 130,140,130 -> skip the x bump at j==2
                        x = x + h / 2.0
                    _fct(x, y, dery)

            # ---- accuracy test (label 150) -------------------------------
            if itest == 0:
                # the step just taken was the *double* step
                aux[3, :] = y[:ndim]
                itest = 1
                istep = istep + istep - 2
                # fall through to label 180 (halve and redo)
                ihlf += 1
                x = x - h
                h = h / 2.0
                y[:ndim] = aux[0, :]
                dery[:ndim] = aux[1, :]
                aux[5, :] = aux[2, :]
                continue

            imod = istep // 2
            if istep - imod - imod != 0:  # odd: only the first half is done
                _fct(x, y, dery)
                aux[4, :] = y[:ndim]
                aux[6, :] = dery[:ndim]
                continue

            # label 230 -- relative, weighted, per-component error
            delt = 0.0
            for i in range(ndim):
                arg = abs(aux[3, i] + y[i])
                if arg == 0.0:
                    arg = 0.25 * abs(aux[3, i])
                if arg == 0.0:
                    arg = ERRSET
                rer = aux[7, i] * abs(aux[3, i] - y[i]) / arg
                delt = max(delt, rer)

            if delt > prmt[3]:  # label 250 -- error too great
                if abs(h) < stpmin:
                    return _finish(IHLF_TOO_MANY_BISECTIONS)
                aux[3, :] = aux[4, :]
                istep = istep + istep - 4
                x = x - h
                iend = 0
                ihlf += 1
                x = x - h
                h = h / 2.0
                y[:ndim] = aux[0, :]
                dery[:ndim] = aux[1, :]
                aux[5, :] = aux[2, :]
                continue

            break  # label 280 -- result values are good

        # ---- accept the step (labels 280-310) ----------------------------
        _fct(x, y, dery)
        aux[0, :] = y[:ndim]
        aux[1, :] = dery[:ndim]
        aux[2, :] = aux[5, :]
        y[:ndim] = aux[4, :]
        dery[:ndim] = aux[6, :]
        _fct(x - h, y, dery)
        _outp(x - h, y, dery, ihlf)
        if prmt[4] != 0.0:
            return RKGSTResult(x - h, y, dery, IHLF_OK, nfev, nout)

        y[:ndim] = aux[0, :]
        dery[:ndim] = aux[1, :]
        irec = ihlf
        if iend > 0:  # label 390 -- interval complete
            _fct(x, y, dery)
            _outp(x, y, dery, ihlf)
            return RKGSTResult(x, y, dery, IHLF_OK, nfev, nout)

        # ---- label 320 -- double the increment ---------------------------
        ihlf -= 1
        istep = istep // 2
        h = h + h
        if abs(h) >= stpmax:
            continue

        # Press-style expansion, clipped at stpmax; resets the halving ladder
        if delt != 0.0:
            trial = 0.9 * abs(h) * (abs(prmt[3] / delt)) ** 0.25
        else:
            trial = 10.0 * abs(h)
        h = math.copysign(min(trial, stpmax), h)
        ihlf = -1
        istep = 0


def quadrature(
    integrand: Callable[[float], float],
    a: float,
    b: float,
    step: float,
    errbnd: float,
    *,
    stpmax: float | None = None,
) -> float:
    """One-dimensional quadrature the way ``ALPH`` and ``SZF`` do it.

    Both set up a single-equation system ``dy/dx = integrand(x)`` with
    ``y(a) = 0`` and read the result off ``y(b)``.  Reproducing this -- rather
    than calling :func:`scipy.integrate.quad` -- matters because the ``RKGST``
    error bound is loose (``ERBNDZ = 0.005`` for ``ALPH``) and the resulting
    quadrature error feeds straight into the fitted wind-profile exponent.
    """
    if stpmax is None:
        stpmax = abs(b - a)
    prmt = [a, b, step, errbnd, stpmax] + [0.0] * 20

    def fct(x, y, dery, p):
        dery[0] = integrand(x)

    def outp(x, y, dery, ihlf, ndim, p):
        return None

    res = rkgst(fct, outp, prmt, [0.0], [1.0], ndim=1)
    if res.ihlf >= 10:
        raise RuntimeError(f"RKGST quadrature failed, IHLF={res.ihlf}")
    return float(res.y[0])
