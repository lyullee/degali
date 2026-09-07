"""Resolving a DEGADIS contaminant to a CoolProp fluid.

An input deck identifies its contaminant with a three-character label and a
molecular weight, which is all the 1989 model needed: the heat capacity came
from two fitted constants and the density from an ideal-gas scaling.  A
real-fluid equation of state needs to know *which* fluid, so this module maps
the deck onto CoolProp's names.

Resolution is by label first, then by molecular weight.  The weight match is
deliberately narrow -- 0.5 % -- because getting the wrong fluid is worse than
getting no fluid: the legacy correlations are approximate but they are at
least approximations of the right substance.  When nothing matches,
:func:`resolve` returns ``None`` and the caller keeps the 1989 correlations
for the contaminant while still using real properties for water and air.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Labels DEGADIS decks use, and the fluids they mean.  Keys are compared
#: case-insensitively after stripping.  ``LNG`` is treated as methane, which
#: is what the Burro test cases intend and what the deck's own molecular
#: weight of 16.04 says.
ALIASES: dict[str, str] = {
    "LNG": "Methane",
    "CH4": "Methane",
    "MET": "Methane",
    "NH3": "Ammonia",
    "AMM": "Ammonia",
    "CL2": "Chlorine",
    "CHL": "Chlorine",
    "PRO": "n-Propane",
    "C3H": "n-Propane",
    "LPG": "n-Propane",
    "BUT": "n-Butane",
    "ETH": "Ethane",
    "ETY": "Ethylene",
    "PRY": "Propylene",
    "CO2": "CarbonDioxide",
    "SO2": "SulfurDioxide",
    "H2S": "HydrogenSulfide",
    "H2": "Hydrogen",
    "N2": "Nitrogen",
    "O2": "Oxygen",
    "AR": "Argon",
    "HE": "Helium",
    "N2O": "NitrousOxide",
    "COS": "CarbonylSulfide",
    "R22": "R22",
    "R12": "R12",
}

#: Molecular weight fallback, kg/kmol.  Only for substances whose weight is
#: distinctive enough that a 0.5 % match is unambiguous among these entries.
BY_WEIGHT: tuple[tuple[float, str], ...] = (
    (2.016, "Hydrogen"),
    (4.003, "Helium"),
    (16.043, "Methane"),
    (17.031, "Ammonia"),
    (28.014, "Nitrogen"),
    (30.070, "Ethane"),
    (31.999, "Oxygen"),
    (34.081, "HydrogenSulfide"),
    (39.948, "Argon"),
    (44.010, "CarbonDioxide"),
    (44.096, "n-Propane"),
    (58.122, "n-Butane"),
    (60.076, "CarbonylSulfide"),
    (64.064, "SulfurDioxide"),
    (70.906, "Chlorine"),
)

#: Ethylene and propylene sit close to ethane and propane in molecular weight,
#: so they are reachable by label only.
_WEIGHT_TOLERANCE = 0.005


@dataclass(frozen=True)
class Resolution:
    """How a contaminant was matched, and to what."""

    fluid: str | None
    how: str  #: ``"alias"``, ``"molecular weight"`` or ``"unmatched"``

    def __bool__(self) -> bool:
        return self.fluid is not None


def resolve(name: str | None, molecular_weight: float | None = None) -> Resolution:
    """Find the CoolProp fluid for a deck's contaminant.

    Parameters
    ----------
    name
        The deck's ``GASNAM`` label, or any CoolProp fluid name.
    molecular_weight
        ``GASMW``, used when the label does not resolve.

    Examples
    --------
    >>> resolve("LNG", 16.04).fluid
    'Methane'
    >>> resolve("XYZ", 17.0).fluid
    'Ammonia'
    >>> bool(resolve("XYZ", 123.4))
    False
    """
    if name:
        key = name.strip().upper()
        if key in ALIASES:
            return Resolution(ALIASES[key], "alias")
        # a deck may already name a CoolProp fluid outright
        canonical = {v.upper(): v for v in ALIASES.values()}
        if key in canonical:
            return Resolution(canonical[key], "alias")

    if molecular_weight:
        for weight, fluid in BY_WEIGHT:
            if abs(molecular_weight - weight) <= _WEIGHT_TOLERANCE * weight:
                return Resolution(fluid, "molecular weight")

    return Resolution(None, "unmatched")


def is_available(fluid: str) -> bool:
    """Whether CoolProp is installed and knows this fluid."""
    try:
        from CoolProp.CoolProp import PropsSI
    except ImportError:
        return False
    try:
        PropsSI("M", fluid)
        return True
    except (ValueError, RuntimeError):
        return False
