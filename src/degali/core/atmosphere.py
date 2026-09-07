"""Atmospheric boundary layer description.

Ports the following DEGADIS 2.1 units:

===================  =====================================================
Fortran              Python
===================  =====================================================
``ATMDEF.FOR``       :func:`stability_defaults`
``PSIF.FOR``         :func:`psi`
``ALPH.FOR``         :func:`friction_velocity`, :func:`fit_alpha`
``IO.FOR`` (part)    :func:`absolute_humidity`, :func:`ambient_density`
===================  =====================================================

Modernisation notes
-------------------
* ``ALPH`` drove ``ZBRENT`` around ``ALPHI``, which in turn drove ``RKGST``
  (Runge-Kutta-Gill) over a one-dimensional quadrature.  Here the quadrature is
  :func:`scipy.integrate.quad` and the root find is
  :func:`scipy.optimize.brentq`.  Both are far more accurate than the original
  fixed-tolerance schemes, so the fitted ``alpha`` can differ in the 5th
  decimal; :func:`fit_alpha` accepts ``legacy_rkgst=True`` to reproduce the
  original integrator exactly when bit-comparison is wanted.
* The stability tables are pure data and are reproduced verbatim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from scipy.integrate import quad
from scipy.optimize import brentq

from .constants import GG, RGAS, VKC, WMA, WMW
from .numerics import zbrent
from .rkgst import quadrature

StabilityClass = Literal["A", "B", "C", "D", "E", "F"]

#: Pasquill-Gifford class -> 1-based ``ISTAB`` index used by the Fortran.
STABILITY_INDEX = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6}
STABILITY_NAME = {v: k for k, v in STABILITY_INDEX.items()}


@dataclass(frozen=True)
class StabilityDefaults:
    """Output of ``ATMDEF``.

    Attributes
    ----------
    rml
        Monin-Obukhov length, m.  ``0.0`` is the sentinel the Fortran uses for
        *infinite* (neutral, class D); :func:`psi` treats it that way.
    deltay, betay
        Lateral dispersion parameters, :math:`\\sigma_y` grows as
        ``sqrt(2) * deltay * (x + x_v)**betay``.
    deltaz, betaz, gammaz
        Vertical (passive, Gaussian-phase) dispersion parameters.
    sigxco, sigxp, sigxmd
        Along-wind dispersion correction coefficients used by DEG3/DEG4.
    """

    rml: float
    deltay: float
    betay: float
    deltaz: float
    betaz: float
    gammaz: float
    sigxco: float
    sigxp: float
    sigxmd: float


# ``ATMDEF.FOR``: (min averaging time, deltay coefficient, deltaz, betaz,
#                  gammaz, rml coefficient, rml exponent,
#                  sigxco, sigxp, sigxmd)
_ATMDEF_TABLE = {
    1: (18.4, 0.423, 107.66, -1.7172, 0.2770, -11.43, 0.103, 0.02, 1.22, 130.0),
    2: (18.4, 0.313, 0.1355, 0.8752, 0.0136, -25.98, 0.171, 0.02, 1.22, 130.0),
    3: (18.4, 0.210, 0.09623, 0.9477, -0.0020, -123.4, 0.304, 0.02, 1.22, 130.0),
    4: (18.3, 0.136, 0.04134, 1.1737, -0.0316, 0.0, 0.0, 0.04, 1.14, 100.0),
    5: (11.4, 0.102, 0.02275, 1.3010, -0.0450, 123.4, 0.304, 0.17, 0.97, 50.0),
    6: (4.6, 0.0674, 0.01122, 1.4024, -0.0540, 25.98, 0.171, 0.17, 0.97, 50.0),
}


def stability_defaults(
    zr: float, istab: int | StabilityClass, avtime: float
) -> StabilityDefaults:
    """Port of ``ATMDEF``.

    Parameters
    ----------
    zr
        Surface roughness length, m.
    istab
        Pasquill-Gifford class, either the letter ``"A"``..``"F"`` or the
        1-based integer the Fortran used.
    avtime
        Averaging time, s.  Clipped from below by a class-dependent minimum.
    """
    if isinstance(istab, str):
        istab = STABILITY_INDEX[istab.strip().upper()]
    if istab not in _ATMDEF_TABLE:
        raise ValueError(f"stability index must be 1..6, got {istab!r}")

    tmin, dycoef, deltaz, betaz, gammaz, rmlc, rmlp, sxc, sxp, sxm = _ATMDEF_TABLE[istab]
    timeav = max(avtime, tmin)
    deltay = dycoef * (timeav / 600.0) ** 0.2
    # class D uses rml = 0. as the "infinite" sentinel; rmlc is 0. there so the
    # expression below reproduces it without a special case.
    rml = rmlc * zr**rmlp if rmlc != 0.0 else 0.0

    return StabilityDefaults(
        rml=rml,
        deltay=deltay,
        betay=0.9,
        deltaz=deltaz,
        betaz=betaz,
        gammaz=gammaz,
        sigxco=sxc,
        sigxp=sxp,
        sigxmd=sxm,
    )


def psi(z: float, rml: float) -> float:
    """Port of ``PSIF``: the Businger stability correction to the log profile.

    ``rml == 0`` means an infinite Monin-Obukhov length (neutral), for which
    the correction vanishes.
    """
    if rml < 0.0:
        a = (1.0 - 15.0 * z / rml) ** 0.25
        return (
            2.0 * math.log((1.0 + a) / 2.0)
            + math.log((1.0 + a * a) / 2.0)
            - 2.0 * math.atan(a)
            + math.pi / 2.0
        )
    if rml == 0.0:
        return 0.0
    return -4.7 * z / rml


def wind_log(z: float, ustar: float, zr: float, rml: float) -> float:
    """Monin-Obukhov log-law wind speed at height ``z`` (m/s)."""
    return ustar / VKC * (math.log((z + zr) / zr) - psi(z, rml))


def wind_power(z: float, u0: float, z0: float, alpha: float) -> float:
    """Power-law wind speed at height ``z`` (m/s).

    This is the profile DEGADIS actually integrates against;
    :func:`fit_alpha` chooses ``alpha`` so that it best matches
    :func:`wind_log`.
    """
    return u0 * (z / z0) ** alpha


def friction_velocity(u0: float, z0: float, zr: float, rml: float) -> float:
    """Port of the ``USTAR`` calculation at the top of ``ALPH``.

    Note DEGADIS uses a von Karman constant of 0.35, not the more common 0.40.
    """
    if u0 == 0.0:
        return 0.0
    return u0 * VKC / (math.log((z0 + zr) / zr) - psi(z0, rml))


def _alpha_integrand(
    alpha: float,
    u0: float,
    z0: float,
    zr: float,
    rml: float,
    ustar: float,
    ialpfl: int,
):
    """Port of ``ARG``: the integrand of the weighted profile-fit residual."""

    def f(z: float) -> float:
        w = 1.0 if ialpfl == 2 else 1.0 / (1.0 + z)
        ubst = wind_log(z, ustar, zr, rml)
        ualp = u0 * (z / z0) ** alpha
        return w * (ubst - ualp) * math.log(z / z0) * ualp

    return f


def _alpha_residual(
    alpha: float,
    u0: float,
    z0: float,
    zr: float,
    rml: float,
    ustar: float,
    zlow: float,
    ialpfl: int,
    *,
    legacy: bool = True,
    stpinz: float = -0.1,
    erbndz: float = 0.005,
) -> float:
    r"""Port of ``ALPHI``.

    The residual whose root defines ``alpha``:

    .. math::
        \int_{z_0}^{z_{low}} W(z)\,[u_{log}(z) - u_{pow}(z)]
        \, \ln(z/z_0) \, u_{pow}(z) \; dz = 0

    with :math:`W = 1/(1+z)` for ``ialpfl == 1`` and :math:`W = 1` for
    ``ialpfl == 2``.  The integration runs *downward*, from ``z0`` to
    ``zlow``; the original passes ``STPINZ = -0.1`` for exactly that reason.

    ``legacy=True`` integrates with the ported ``RKGST`` at the original's
    ``ERBNDZ = 0.005``.  That bound is loose enough that the quadrature error
    shifts the root by about 2e-5 relative, so reproducing reference output
    means reproducing the quadrature rather than improving it.
    ``legacy=False`` uses :func:`scipy.integrate.quad`.
    """
    integrand = _alpha_integrand(alpha, u0, z0, zr, rml, ustar, ialpfl)
    if legacy:
        return quadrature(integrand, z0, zlow, stpinz, erbndz, stpmax=abs(zlow - z0))
    # accurate path: split near the lower limit, where the slope has a log kink
    mid = min(z0, zlow * 10.0) if zlow > 0 else z0
    val = 0.0
    if mid > zlow:
        val += quad(integrand, zlow, mid, limit=200)[0]
    if z0 > mid:
        val += quad(integrand, mid, z0, limit=200)[0]
    return -val


def fit_alpha(
    u0: float,
    z0: float,
    zr: float,
    rml: float,
    *,
    ustar: float | None = None,
    zlow: float = 0.01,
    xli: float = 0.05,
    xri: float = 0.50,
    eps: float = 1.0e-5,
    stpinz: float = -0.1,
    erbndz: float = 0.005,
    ialpfl: int = 1,
    alpco: float = 0.2,
    legacy: bool = True,
) -> float:
    """Port of ``ALPH``: fit the power-law exponent to the log wind profile.

    Parameters follow the ``ER1`` parameter file: ``XLI``/``XRI`` bracket the
    search, ``EPS`` is the convergence criterion, ``ZLOW`` the lower limit of
    the fit, ``STPINZ``/``ERBNDZ`` control the quadrature, and ``IALPFL``
    selects the weight function (``0`` means "do not fit, use ``ALPCO``").

    With ``legacy=True`` (the default) the original ``RKGST`` quadrature and
    ``ZBRENT`` root find are used, reproducing DEGADIS to round-off.  With
    ``legacy=False`` the fit is solved accurately; the two typically differ by
    about 2e-5 relative, which is the original's own convergence tolerance.

    Returns
    -------
    float
        The exponent ``alpha``; ``0.0`` for a windless simulation, matching
        the Fortran.
    """
    if u0 == 0.0:
        return 0.0
    if ialpfl == 0:
        return alpco
    if ustar is None:
        ustar = friction_velocity(u0, z0, zr, rml)

    lo = max(zlow, zr)  # ALPHI: "to take care of large zr"

    def f(a: float) -> float:
        return _alpha_residual(
            a,
            u0,
            z0,
            zr,
            rml,
            ustar,
            lo,
            ialpfl,
            legacy=legacy,
            stpinz=stpinz,
            erbndz=erbndz,
        )

    if legacy:
        return zbrent(f, xli, xri, eps)
    return brentq(f, xli, xri, xtol=eps, rtol=1e-15)


def saturation_humidity(tamb: float, pamb: float, vapour_pressure) -> float:
    """Saturation absolute humidity, kg water / kg bone-dry air.

    ``vapour_pressure`` is a callable ``T -> p`` in the same units as ``pamb``
    (the legacy backend supplies ``watvp`` in atm).
    """
    vp = vapour_pressure(tamb)
    return WMW / WMA * vp / (pamb - vp)


def absolute_humidity(
    tamb: float,
    pamb: float,
    vapour_pressure,
    *,
    humid: float = 0.0,
    relhum: float = 0.0,
) -> tuple[float, float]:
    """Port of the humidity block of ``IO.FOR``.

    ``relhum`` (per cent) wins if non-zero; otherwise ``humid`` (kg/kg BDA) is
    used to back out the relative humidity.  Returns ``(humid, relhum)``.
    """
    sat = saturation_humidity(tamb, pamb, vapour_pressure)
    if relhum > 0.0:
        return relhum / 100.0 * sat, relhum
    if humid > 0.0:
        return humid, 100.0 * humid / sat
    return 0.0, 0.0


def ambient_density(tamb: float, pamb: float, humid: float) -> float:
    """Port of the ``RHOA`` expression standardised by MCB#6 (1991).

    .. math::
        \\rho_a = \\frac{P (1 + w) M_w}{R T \\, (M_w/M_a + w)}

    with ``w`` the absolute humidity in kg water per kg bone-dry air.
    """
    return pamb * (1.0 + humid) * WMW / (RGAS * (WMW / WMA + humid)) / tamb


def richardson_star(rho: float, rhoa: float, heff: float, ustar: float) -> float:
    """Port of ``RIF``: the bulk Richardson number ``Ri*``."""
    return GG * (rho - rhoa) / rhoa * heff / ustar / ustar


def richardson_thermal(
    tsurf: float, temp: float, heff: float, ustar: float, wind: float
) -> float:
    """Port of ``RIFT``: the thermal Richardson number ``Ri_t``, clipped at 0."""
    return max(GG * (tsurf - temp) / temp * heff / ustar / wind, 0.0)
