"""Numerical helpers.

DEGADIS 2.1 shipped its own root finder (``ZBRENT``), ODE integrator
(``RKGST``, Runge-Kutta-Gill with step halving), special functions
(``GAMMA``, ``INCGAMMA``, ``ERF``, ``SERIES``), table interpolators
(``AFGEN``, ``AFGEN2``) and a Gauss elimination routine (``SIMUL``).

Everything here has a modern equivalent in SciPy/NumPy.  Two things are still
worth keeping:

1. An exact transcription of ``ZBRENT``.  The Fortran calls it with loose
   tolerances (``acrit = 0.001`` K on temperature inversion), so the root it
   returns is a specific point *near* the true root, not the true root.  When
   reproducing reference output bit for bit, that specific point matters.
2. ``series`` and ``gaminc``, which are DEGADIS-specific combinations of
   incomplete gamma functions used in the flammable-mass integrals; the
   original series expansions converge to a fixed relative criterion.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy.optimize import brentq as _brentq
from scipy.special import gamma as _gamma
from scipy.special import gammainc as _gammainc


class RootBracketError(ValueError):
    """The supplied interval does not bracket a root (``ZBRENT`` ``ierr = 2``)."""


class RootConvergenceError(RuntimeError):
    """``ZBRENT`` exhausted ``itmax`` iterations (``ierr = 1``)."""


def zbrent(
    func: Callable[[float], float],
    x1: float,
    x2: float,
    tol: float,
    *,
    itmax: int = 100,
    eps: float = 3.0e-8,
) -> float:
    """Exact transcription of ``ZBRENT.FOR`` (Press et al., 1st ed., p. 253).

    Kept verbatim rather than delegating to :func:`scipy.optimize.brentq`
    because the convergence test ``tol1 = 2*eps*|b| + tol/2`` and the
    ``b += sign(tol1, xm)`` fallback step place the returned root at a
    reproducible point that SciPy's slightly different test would not hit.

    ``tol`` is an *absolute* tolerance on the independent variable.
    """
    a, b = float(x1), float(x2)
    fa, fb = func(a), func(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if math.copysign(1.0, fa) * math.copysign(1.0, fb) > 0.0:
        raise RootBracketError(f"f({a})={fa} and f({b})={fb} do not bracket a root")

    c, fc = b, fb
    d = e = 0.0
    for _ in range(itmax):
        if math.copysign(1.0, fb) * math.copysign(1.0, fc) > 0.0:
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb

        tol1 = 2.0 * eps * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        if abs(xm) <= tol1 or fb == 0.0:
            return b

        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q)):
                e, d = d, p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d

        a, fa = b, fb
        b = b + d if abs(d) > tol1 else b + math.copysign(tol1, xm)
        fb = func(b)

    raise RootConvergenceError(f"zbrent failed to converge in {itmax} iterations")


def brentq(func, a, b, *, xtol=2e-12, rtol=8.9e-16, maxiter=100):
    """:func:`scipy.optimize.brentq`, raising this module's own exceptions.

    Callers switch between this and :func:`zbrent` on a flag, so the two must
    fail the same way -- otherwise a bracket-widening fallback written for one
    of them silently stops working for the other.
    """
    try:
        return _brentq(func, a, b, xtol=xtol, rtol=rtol, maxiter=maxiter)
    except ValueError as exc:
        raise RootBracketError(str(exc)) from exc
    except RuntimeError as exc:
        raise RootConvergenceError(str(exc)) from exc


def gaminc(a: float, x: float) -> float:
    r"""Port of ``GAMINC`` in ``INCGAMMA.FOR``: the *unregularised* lower
    incomplete gamma function.

    .. math:: \gamma(a, x) = \int_0^x t^{a-1} e^{-t} \, dt

    The routine computes the Numerical Recipes ``gammp`` -- which *is*
    regularised -- and then undoes the normalisation on the last line,
    ``gaminc = exp(log(aa) + gln)``. SciPy's :func:`scipy.special.gammainc` is
    the regularised form, so it has to be multiplied back by
    :math:`\Gamma(a)`.

    The distinction is easy to miss and expensive: in ``PSS`` the result is
    compared against a cap built from :math:`\Gamma(1/(1+\alpha))`, so using
    the regularised form makes the flammable-mass derivative low by exactly
    that factor -- 6.5 per cent for a typical wind profile -- while leaving
    every other quantity in the routine correct.
    """
    if x <= 0.0:
        return 0.0
    return float(_gammainc(a, x)) * float(_gamma(a))


def gamma(x: float) -> float:
    """Port of ``GAMMA.FOR`` (Lanczos series in the original)."""
    return float(_gamma(x))


def afgen(
    xtab: np.ndarray, ytab: np.ndarray, x: float, *, sentinel: float = -1.0e-20
) -> float:
    """Port of ``AFGEN``/``AFGEN2``: piecewise-linear table lookup.

    The Fortran tables are fixed-length arrays terminated by the ``POUND``
    sentinel (``-1e-20``); values outside the table are clamped to the end
    points.  ``sentinel`` selects where the valid data stops.
    """
    x_arr = np.asarray(xtab, dtype=float)
    y_arr = np.asarray(ytab, dtype=float)
    stop = np.argmax(x_arr == sentinel)
    if x_arr[stop] == sentinel:
        x_arr, y_arr = x_arr[:stop], y_arr[:stop]
    return float(np.interp(x, x_arr, y_arr))


def series(w: float, alpha1: float, *, crit: float = 1.0e-7) -> float:
    """Port of ``SERIES.FOR``.

    Evaluates :math:`\\sum_{n} w^{n/\\alpha_1}` style expansion used in the
    UFL/LFL mass integrals.  The concrete form is filled in when
    ``downwind.py`` lands; the signature is fixed here so callers are stable.
    """
    raise NotImplementedError("ported alongside downwind.py")


def limit(
    func: Callable[[float], float],
    x0: float,
    xinc: float,
    xhigh: float,
    xlow: float,
    *,
    rmax: float = 400.0,
) -> tuple[float, float]:
    """Port of ``LIMIT.FOR``: widen a bracket until it contains a sign change.

    ``ZBRENT`` needs a bracketing interval and DEGADIS cannot always supply
    one -- ``ADDHEAT`` in particular derives its lower bound from the adiabatic
    mixing temperature, which stops bracketing the root once the ground has
    added enough heat. ``LIMIT`` then steps outward from ``x0`` in both
    directions by growing multiples of ``xinc`` until either side changes sign,
    and returns the interval straddling that crossing.

    The original comment calls it "expensive"; the search is linear in the
    number of increments, capped at 400.

    Returns
    -------
    (xhigh, xlow)
        The narrowed bracket, in the Fortran's argument order.
    """
    rrr = 1.0
    aflag = bflag = True
    fc = func(x0)
    fa = fb = 0.0
    aaa = bbb = x0

    while True:
        if aflag:
            aaa = x0 + xinc * rrr
            if aaa >= xhigh:
                aaa = xhigh
                aflag = False
            fa = func(aaa)
            if math.copysign(1.0, fa) * math.copysign(1.0, fc) < 0.0:
                return aaa, aaa - xinc * 1.01

        if bflag:
            bbb = x0 - xinc * rrr
            if bbb <= xlow:
                bbb = xlow
                bflag = False
            fb = func(bbb)
            if math.copysign(1.0, fb) * math.copysign(1.0, fc) < 0.0:
                return bbb + xinc * 1.01, bbb

        if math.copysign(1.0, fa) * fb > 0.0:
            rrr += 1.0
            if rrr > rmax:
                raise RootBracketError("LIMIT failed to find a bracket")
            continue

        return aaa, bbb
