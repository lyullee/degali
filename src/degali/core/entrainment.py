"""Air entrainment closure and surface exchange.

Ports ``RIPHIF.FOR`` and ``SURFAC.FOR``:

=====================  ==================================================
Fortran                Python
=====================  ==================================================
``RIF``                :func:`richardson_star` (in :mod:`.atmosphere`)
``RIFT``               :func:`richardson_thermal` (in :mod:`.atmosphere`)
``PHIF``               :func:`phi`
``PHIHAT``             :func:`phi_hat`
``GSERIES``            :func:`hypergeometric_series`
``SURFAC``             :func:`surface_exchange`
=====================  ==================================================

The physics
-----------
Vertical entrainment of air into a dense cloud is suppressed by stable
stratification.  DEGADIS expresses this through :math:`\\Phi(Ri^*)`, a
correction to the neutral entrainment velocity:

.. math::
    w_e = \\frac{\\delta_{lay}\\, \\kappa\\, u_* (1+\\alpha)}{\\Phi}

so a larger :math:`\\Phi` means less entrainment.  Five prescriptions ship with
the model, selected by ``IPHIFL`` in the ``ER1`` file; the test cases use
``IPHIFL = 3``.  All of them are of the form "constant, plus a power of
:math:`Ri^*` when the layer is stably stratified", optionally divided by a
thermal-stratification correction built from :math:`Ri_t`.

:func:`phi_hat` is different in kind.  It is the *fetch-averaged* value of
:math:`\\Phi` needed to close the maximum atmospheric take-up rate
:math:`q^*_{max}` over the secondary source, and the original evaluates it in
closed form as a Gauss hypergeometric function :math:`{}_2F_1`, summed
directly.  Two branches are needed because the series only converges for
:math:`|z| < 1`; the second uses the standard :math:`z \\to 1/z` connection
formula.  Here :func:`scipy.special.hyp2f1` handles both branches at once, but
the original series is retained under ``legacy=True`` because its convergence
criterion (``crit = 7e-5``, and a hard clamp of ``|z|`` at 0.999) makes it
return a specific approximation rather than the true value.
"""

from __future__ import annotations

import math
from typing import Callable

from scipy.special import gamma as _gamma
from scipy.special import hyp2f1 as _hyp2f1

from .constants import GG, VKC

#: Constant appearing in the ``PHIHAT`` closure, ``data phic/3.1D0/``.
PHIC = 3.1

#: The Fortran writes the cube-root exponent as the literal ``0.333333333333``
#: rather than ``1/3``.  The 1e-12 truncation is visible in the heat-transfer
#: coefficient, so the literal is reproduced.
_ONE_THIRD = 0.333333333333


# ==========================================================================
# PHIF -- the entrainment suppression function
# ==========================================================================


def phi(ri: float, rit: float = 0.0, iphifl: int = 3) -> float:
    """Port of ``PHIF``.

    Parameters
    ----------
    ri
        Bulk Richardson number :math:`Ri^*` based on the density excess.
    rit
        Thermal Richardson number :math:`Ri_t`; used only by prescriptions 3
        and 4, which divide by :math:`\\sqrt{1 + 0.25\\,Ri_t^{2/3}}` to account
        for convective mixing driven by a warm ground.
    iphifl
        Prescription selector, 1-5:

        1. Neutral value 0.74; the original Havens and Spicer (1985) fit.
        2. Neutral value 0.88 with a steeper stable branch.
        3. As 2, with the thermal correction applied (**the model default**,
           and what all five EPA test cases use).
        4. As 3, but the unstable branch is flattened to the neutral value.
        5. Constant 0.88 -- entrainment independent of stratification.
    """
    if iphifl == 1:
        if ri < 0.0:
            return 0.74 / (1.0 + 0.65 * abs(ri) ** 0.6)
        if ri == 0.0:
            return 0.74
        return 0.74 + 0.25 * ri**0.7 + 1.2e-7 * ri**3

    if iphifl == 2:
        if ri < 0.0:
            return 0.88 / (1.0 + 0.65 * abs(ri) ** 0.6)
        if ri == 0.0:
            return 0.88
        return 0.88 + 9.9e-2 * ri**1.04 + 1.4e-25 * ri**5.7

    if iphifl in (3, 4):
        corr1 = 0.25 * rit**0.666666 + 1.0
        corr = math.sqrt(corr1)
        riw = ri / corr1
        if ri > 0.0:
            return (0.88 + 9.9e-2 * riw**1.04 + 1.4e-25 * riw**5.7) / corr
        if iphifl == 4 or ri == 0.0:
            return 0.88 / corr
        return 0.88 / (1.0 + 0.65 * abs(riw) ** 0.6) / corr

    if iphifl == 5:
        return 0.88

    raise ValueError(f"IPHIFL must be 1..5, got {iphifl!r} (TRAP 29)")


def entrainment_velocity(
    ustar: float, alpha1: float, dellay: float, phi_value: float
) -> float:
    """Neutral entrainment velocity divided by :math:`\\Phi`.

    This is the ``wel = dellay * vkc*ustar*alpha1/phi`` line that appears in
    ``SZLOCAL`` and, in equivalent form, in ``PSS`` and ``SSG``.
    """
    return dellay * VKC * ustar * alpha1 / phi_value


# ==========================================================================
# PHIHAT -- fetch-averaged Phi for the secondary-source take-up rate
# ==========================================================================


def hypergeometric_series(
    aaa: float, bbb: float, ccc: float, zzz: float, *, crit: float = 7.0e-5
) -> float:
    """Port of ``GSERIES``.

    Sums :math:`{}_2F_1(a, b; c; z)` term by term, stopping when a term falls
    below ``crit`` in relative size, and returns the mean of the last two
    partial sums.  ``|z|`` is clamped at 0.999 -- the original does this
    explicitly "to avoid excessive execution times", and the clamp changes the
    answer, so it is reproduced rather than fixed.
    """
    z = math.copysign(0.999, zzz) if abs(zzz) > 0.999 else zzz

    term = aaa * bbb / ccc * z
    sumo = 1.0 + term
    for k in range(1, 100001):
        rk = float(k)
        term = term * (aaa + rk) * (bbb + rk) / (ccc + rk) * z / (rk + 1.0)
        total = sumo + term
        if abs(term / total) <= crit:
            return (total + sumo) / 2.0
        sumo = total
    raise RuntimeError("GSERIES did not converge in RIPHIF")


def phi_hat(
    rho: float,
    fetch: float,
    *,
    rhoa: float,
    alpha: float,
    ustar: float,
    u0: float,
    z0: float,
    gammaf: float,
    dellay: float,
    legacy: bool = True,
) -> float:
    r"""Port of ``PHIHAT``.

    The maximum atmospheric take-up rate over a secondary source of streamwise
    extent ``fetch`` is

    .. math::
        q^*_{max} = C_c\, \kappa\, u_*\, (1+\alpha)
                    \frac{\delta_{lay}}{\delta_{lay}-1}
                    \frac{1}{\hat\Phi}

    where :math:`\hat\Phi` is the value of :math:`\Phi` averaged along the
    fetch.  Because :math:`\Phi` follows a power law in :math:`Ri^*` and
    :math:`Ri^*` itself grows as a power of fetch, the average has a closed
    form in terms of :math:`{}_2F_1`.

    A cloud no denser than ambient entrains at the neutral rate, so
    :math:`\hat\Phi = 0.88`.
    """
    if rho <= rhoa:
        return 0.88

    alpha1 = alpha + 1.0
    pow_ = 1.0 / alpha1
    p1 = 1.04 / alpha1
    p1i = 1.0 / p1
    p2 = 1.0 + p1i
    p3 = (alpha - 0.04) / 1.04
    p4 = (1.08 - alpha) / 1.04

    ci = GG * (rho - rhoa) / rhoa * z0 / ustar**2 * gammaf / alpha1
    ci *= (VKC * ustar * alpha1**2 / u0 / z0 / PHIC * dellay / (dellay - 1.0)) ** pow_

    ri = ci * fetch**pow_
    zzz = -0.099 * ri**1.04 / 0.88

    series: Callable[..., float]
    series = hypergeometric_series if legacy else (lambda a, b, c, z: float(_hyp2f1(a, b, c, z)))

    if abs(zzz) < 1.0:
        return 0.88 / series(1.0, p1i, p2, zzz)

    zinv = 1.0 / zzz
    return 0.88 / (
        -zinv / (1.0 - p1) * series(1.0, -p3, p4, zinv)
        + float(_gamma(p1i)) * float(_gamma(p4)) / (p1 - 1.0) * (-zinv) ** p1i
    )


# ==========================================================================
# SURFAC -- ground-to-cloud heat and water transfer
# ==========================================================================


def surface_exchange(
    temp: float,
    height: float,
    rho: float,
    wm: float,
    cp: float,
    yw: float,
    *,
    tsurf: float,
    pamb: float,
    zr: float,
    u0: float,
    z0: float,
    alpha: float,
    ustar: float,
    isofl: int,
    ihtfl: int,
    iwtfl: int,
    htco: float,
    wtco: float,
    vapour_pressure: Callable[[float], float],
) -> tuple[float, float]:
    """Port of ``SURFAC``: heat and water fluxes from the ground into the cloud.

    Returns
    -------
    (watrte, qrte)
        Water mass flux, kg/(m**2 s), and heat flux, W/m**2.  Both are zero
        when the surface is colder than the cloud, because the correlations
        are one-directional; the caller is expected not to let the cloud cool.

    Heat transfer prescriptions (``IHTFL``)
    ---------------------------------------
    ``-1``
        Fixed coefficient ``HTCO``.
    ``0``
        No heat transfer; the cloud follows the adiabatic mixing line.
    ``1``
        Local correlation: the larger of a natural-convection term
        ``18 * ((rho/M)**2 |dT|)**(1/3)`` and a forced term
        ``1.22 rho c_p u_*^2 / u(h/2)``, evaluated at half the cloud depth.
    ``2``
        LLNL form, ``h = HTCO * rho * c_p``.
    ``3``
        Colenbrander's method: natural term uses the film temperature and the
        forced term is referenced to the 10 m wind.

    Water transfer (``IWTFL``) is analogous: a fixed coefficient for
    ``IWTFL < 0``, otherwise the larger of a natural-convection mass-transfer
    coefficient and one derived from the heat-transfer coefficient by the
    Chilton-Colburn analogy (the ``20.7 * h / (c_p M)`` term).
    """
    if isofl == 1 or ihtfl == 0 or height <= 0.0:
        return 0.0, 0.0

    deltem = tsurf - temp
    if deltem < 0.0:
        return 0.0, 0.0

    hhh = max(height / 2.0, zr)
    topvel = u0 * (hhh / z0) ** alpha
    pnat = ((rho / wm) ** 2 * abs(deltem)) ** _ONE_THIRD

    if ihtfl == 1:
        hn = 18.0 * pnat
        hf = 1.22 * rho * cp * ustar**2 / topvel if u0 != 0.0 else 0.0
        ho = max(hn, hf)
    elif ihtfl == 2:
        ho = htco * rho * cp
    elif ihtfl == 3:
        avtemp = (tsurf + temp) / 2.0
        hn = 89.0 * (deltem / avtemp**2) ** _ONE_THIRD
        u10 = u0 * (10.0 / z0) ** alpha
        hf = 1.22 * rho * cp * ustar**2 / u10 if u0 != 0.0 else 0.0
        ho = max(hn, hf)
    else:  # ihtfl == -1
        ho = htco

    qrte = max(ho * deltem, 0.0)

    if iwtfl == 0:
        return 0.0, qrte

    if iwtfl > 0:
        fn = 9.9e-3 * pnat
        ff = 20.7 * ho / cp / wm
        fo = max(fn, ff)
    else:
        fo = wtco

    pw_cloud = min(vapour_pressure(temp), yw * pamb)
    watrte = fo * (vapour_pressure(tsurf) - pw_cloud) / pamb
    return max(watrte, 0.0), qrte
