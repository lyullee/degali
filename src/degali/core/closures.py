"""Swappable closures for the ground-level cloud.

DEGADIS makes one particular choice about what a ground-level cloud does with
its density difference: it spreads sideways under gravity while it is denser
than air, and does nothing at all once it is lighter.  That choice is right
for LNG, where the cloud stays dense across the whole flammable range, and
wrong for hydrogen, where it is buoyant everywhere below 99.9 mole per cent.

So the choice is made pluggable rather than fixed.  A closure supplies two
things to the downwind model:

* the lateral spreading rate, ``dB_eff/dx``;
* optionally, extra state of its own -- a rise velocity, an elevation --
  integrated alongside the cloud.

:class:`DegadisClosure` reproduces the original exactly and is the default, so
nothing changes unless a caller asks for something else.  Alternatives live in
:mod:`degali.addons`.

The interface is deliberately narrow.  A closure sees the local cloud state
and returns rates; it does not get to rewrite the mass or energy balances,
which are common to every variant and are validated against the Fortran.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, Sequence

import numpy as np


@dataclass
class CloudState:
    """What a closure is told about the cloud at one downwind station."""

    dist: float  #: downwind distance, m
    beff: float  #: half-width including the Gaussian shoulders, m
    sz: float  #: vertical dispersion parameter, m
    heff: float  #: effective depth, m
    rho: float  #: centreline density, kg/m**3
    rhoa: float  #: ambient density, kg/m**3
    cc: float  #: centreline concentration, kg/m**3
    temp: float  #: K
    wind: float  #: ambient speed at the cloud, m/s
    ustar: float  #: friction velocity, m/s
    #: Mass flux per unit width, ``rho u H``, kg/(m s), and its downwind
    #: rate.  A closure that carries momentum needs both: momentum per unit
    #: mass is diluted by the air the cloud entrains, and leaving that term
    #: out lets a buoyant cloud accelerate without limit.
    mass_flux: float = 0.0
    mass_flux_rate: float = 0.0
    #: The closure's own state vector, as it left it.
    own: np.ndarray = None

    @property
    def delrho(self) -> float:
        """Density excess over ambient. Negative for a buoyant cloud."""
        return self.rho - self.rhoa

    @property
    def buoyant(self) -> bool:
        return self.rho < self.rhoa


class Closure(Protocol):
    """What the downwind model requires of a spreading closure."""

    #: Extra states this closure integrates alongside the cloud.
    n_states: int

    #: Column labels for those states, for reporting.
    state_names: Sequence[str]

    def initial(self, state: CloudState) -> np.ndarray:
        """Initial values for the closure's own states."""

    def lateral_rate(self, state: CloudState) -> float:
        """``dB_eff/dx``, m/m."""

    def derivatives(self, state: CloudState) -> np.ndarray:
        """``d/dx`` of the closure's own states."""

    def elevation(self, state: CloudState) -> float:
        """Height of the cloud centroid above the ground, m.

        Zero for a closure that keeps the cloud on the ground. A non-zero
        value means the concentration profile has to be evaluated about that
        height rather than about the surface.
        """


class DegadisClosure:
    """DEGADIS 2.1 exactly: gravity slumping, and nothing when buoyant.

    ``PSS`` computes

    .. code-block:: fortran

        DERY(iBEFF) = 0.D0
        IF(delrho .GT. delrmn) DERY(iBEFF) = PRMT(9)*sqrt(delrho/rhoa)

    where ``PRMT(9)`` collects the constants of
    :math:`C_E \\sqrt{g z_0 \\Gamma / (1+\\alpha)} \\, \\Gamma / u_0`, and the
    rate carries a further :math:`(S_z/z_0)^{1/2-\\alpha}`.

    The guard is what matters here. Below ``DELRMN`` the cloud stops
    spreading; there is no vertical momentum equation, so it also cannot rise.
    A cloud that becomes lighter than air is frozen in place and disperses
    passively for ever.
    """

    n_states = 0
    state_names: tuple[str, ...] = ()

    def __init__(self, coefficient: float, alpha: float, z0: float,
                 delrmn: float = 0.0):
        self.coefficient = coefficient  # PRMT(9)
        self.alpha = alpha
        self.z0 = z0
        self.delrmn = delrmn

    def initial(self, state: CloudState) -> np.ndarray:
        return np.zeros(0)

    def lateral_rate(self, state: CloudState) -> float:
        if state.delrho <= self.delrmn:
            return 0.0
        return (
            self.coefficient
            * math.sqrt(state.delrho / state.rhoa)
            * (state.sz / self.z0) ** (0.5 - self.alpha)
        )

    def derivatives(self, state: CloudState) -> np.ndarray:
        return np.zeros(0)

    def elevation(self, state: CloudState) -> float:
        return 0.0
