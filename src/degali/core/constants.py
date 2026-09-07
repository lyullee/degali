"""Physical and numerical constants.

Ported from the DEGADIS 2.1 include files ``DEGIN.PRM``, ``DEG1.PRM``,
``DEG2.PRM``, ``DEG3.PRM`` and ``DEG4.PRM``, plus the ``DATA`` statements in
``DEG1.for`` that were effectively compile-time constants.

The values are reproduced exactly.  Where the original carried a rounded
literal (for example ``rgas = 0.08205`` rather than the CODATA value) the
literal is kept under its original name and the modern value is offered
separately, so that ``backend="legacy"`` can reproduce the Fortran bit for bit
while ``backend="coolprop"`` can use better physics.
"""

from __future__ import annotations

import math

# --------------------------------------------------------------------------
# DEGIN.PRM -- shared by every DEGADIS program
# --------------------------------------------------------------------------

WMA = 28.96  #: molecular weight of dry air, kg/kmol
WMW = 18.02  #: molecular weight of water, kg/kmol
CPA = 1006.3  #: heat capacity of dry air, J/(kg K)
CPW = 1865.0  #: heat capacity of water vapour, J/(kg K)
RHOWL = 1000.0  #: liquid water density, kg/m**3
DHVAP = 2.5023e6  #: latent heat of vaporisation of water, J/kg
DHFUS = 0.33e6  #: latent heat of fusion of water, J/kg
RGAS = 0.08205  #: gas constant, atm m**3 / (kmol K)

#: The original hard-codes ``pi`` to 12 digits in DEGIN.PRM; ``math.pi`` is
#: identical to double precision, so the two are interchangeable.
PI = 3.14159265358

# DEG1.PRM / DEG2.PRM derived constants (kept for readability of ported code)
SQRTPI = 1.772453851  #: sqrt(pi)
RT2 = 1.414213562  #: sqrt(2)
SQPIO2 = 1.253314137  #: sqrt(pi/2)

# Array bounds.  These were fixed dimensions in Fortran; in Python they only
# survive as the limits the *input format* implies (e.g. a maximum of 42
# source-description rows in an INP file written by DEGINP).
IGEN = 42  #: max entries in the source table / adiabatic density table
MAXL = 650  #: length of the /GEN3/ secondary-source output vectors
MAXNOB = 50  #: max number of observers (transient runs)
MAXNT = 40  #: max number of sort times
NDOS = 10  #: max dosage receptors

# --------------------------------------------------------------------------
# DATA statements in DEG1.for
# --------------------------------------------------------------------------

GG = 9.81  #: gravitational acceleration, m/s**2
VKC = 0.35  #: von Karman constant as used by DEGADIS (note: not 0.40)
POUND = -1.0e-20  #: end-of-table sentinel used throughout the Fortran

# --------------------------------------------------------------------------
# Unit conversions
# --------------------------------------------------------------------------

ATM_TO_PA = 101325.0

# --------------------------------------------------------------------------
# Modern reference values, used only by the non-legacy backends
# --------------------------------------------------------------------------

R_UNIVERSAL = 8.31446261815324  #: J/(mol K), CODATA 2018 (exact by definition)
RGAS_MODERN = R_UNIVERSAL / ATM_TO_PA * 1000.0  #: atm m**3/(kmol K) = 0.0820574
WMA_MODERN = 28.9647  #: dry air, kg/kmol (Picard et al. 2008)
WMW_MODERN = 18.01528  #: water, kg/kmol

__all__ = [n for n in dir() if n.isupper()] + ["math"]
