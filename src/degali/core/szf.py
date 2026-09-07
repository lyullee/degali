"""Vertical dispersion over a source with no blanket: a port of ``SZF.FOR``.

When the atmosphere can carry the release away as fast as it is produced,
there is no dense pool to integrate in time.  What there is instead is a layer
that thickens along the fetch as it picks up contaminant from the ground and
air from above.  ``SZF`` integrates that growth *in space*, from the upwind
edge of the source to the downwind edge, and reports the layer state at the
end.

The single integrated quantity is

.. math:: Y = \\rho \\, \\delta_{lay} \\, u_{eff} \\, H_{eff}

the mass flux per unit width through the layer, which grows by entrainment
from above and by contaminant from below:

.. math:: \\frac{dY}{dx} = w_e \\rho_a + \\frac{Q}{w_{cp}}

Everything else is recovered algebraically at each step.  The layer-average
contaminant mass fraction follows from a material balance over the fetch so
far, :math:`w_{c,lay} = Qx/Y`; the density comes from the adiabatic mixing
line; the centreline concentration is ``dellay`` times the layer average; and
:math:`\\sigma_z` is inverted from the power-law wind profile.  Entrainment is
then suppressed by :math:`\\Phi(Ri^*)` evaluated on the *centreline* density.

Note that the thermal Richardson number is passed as zero here, so the
ground-heating correction that :func:`~degali.core.entrainment.phi` applies
elsewhere is switched off along the source.
"""

from __future__ import annotations

from dataclasses import dataclass

from .atmosphere import richardson_star
from .constants import VKC
from .entrainment import phi
from .rkgst import rkgst


@dataclass
class LayerState:
    """The layer at the downwind edge of the source."""

    sz: float  #: vertical dispersion parameter, m
    cclay: float  #: layer-average contaminant concentration, kg/m**3
    wclay: float  #: layer-average contaminant mass fraction
    rholay: float  #: layer-average density, kg/m**3
    cc: float  #: centreline concentration, kg/m**3


def sigma_z_over_source(
    q: float,
    length: float,
    wcp: float,
    *,
    table,
    alpha: float,
    gammaf: float,
    u0: float,
    z0: float,
    ustar: float,
    rhoa: float,
    dellay: float,
    iphifl: int = 3,
    szstp0: float = 0.01,
    szerr: float = 0.001,
) -> LayerState:
    """Port of ``SZF``.

    Parameters
    ----------
    q
        Contaminant flux from the source, kg/(m**2 s).
    length
        Streamwise extent of the source, m.
    wcp
        Contaminant mass fraction of the primary source material.
    table
        The adiabatic mixing table, used for both lookups.

    Returns
    -------
    LayerState
        The layer at ``x = length``.  ``SZLOCO`` copies the diagnostics at
        every accepted step, so what comes back is the last *valid* state
        rather than whatever the final derivative call happened to leave
        behind -- a distinction that matters because ``RKGST`` calls the
        derivative routine at trial points it later rejects.
    """
    alpha1 = alpha + 1.0
    diag: dict = {}

    def fct(x, y, dery, prmt):
        if y[0] <= 0.0:
            wclay = 0.0
            rholay = rhoa
            cclay = 0.0
            cc = 0.0
            rho = rhoa
        else:
            # contaminant material balance over the fetch so far
            wclay = q * x / y[0]
            # SZLOCAL passes an uninitialised WALAY to ADIABAT, which reads
            # it. Under /noauto that local is static and zero, and nothing
            # ever writes it, so the lookup runs with wa = 0 for the whole
            # integration. See AdiabaticTable.from_mass_fraction.
            lay = table.from_mass_fraction(wclay, wa=0.0)
            cclay, rholay = lay.cc, lay.rho
            cc = cclay * dellay
            rho = table.from_concentration(cc).rho  # centreline

        uheff = y[0] / rholay / dellay
        sz = (uheff / u0 / z0 * alpha1) ** (1.0 / alpha1) * z0
        heff = gammaf * sz / alpha1
        ristar = richardson_star(rho, rhoa, heff, ustar)
        wel = dellay * VKC * ustar * alpha1 / phi(ristar, 0.0, iphifl)
        dery[0] = wel * rhoa + q / wcp

        diag["trial"] = (sz, cclay, wclay, rholay, cc)

    def outp(x, y, dery, ihlf, ndim, prmt):
        # SZLOCO: promote the trial diagnostics to the accepted ones
        diag["last"] = diag["trial"]

    prmt = [0.0, length, szstp0, szerr, length] + [0.0] * 12
    res = rkgst(fct, outp, prmt, [0.0], [1.0], ndim=1)
    if res.ihlf >= 10:
        raise RuntimeError(f"SZF failed to converge, IHLF={res.ihlf} (TRAP 3)")

    sz, cclay, wclay, rholay, cc = diag["last"]
    return LayerState(sz=sz, cclay=cclay, wclay=wclay, rholay=rholay, cc=cc)
