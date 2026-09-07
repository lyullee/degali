"""DEGALI -- Dense Gas Dispersion for Liquid Hydrogen.

DEGADIS (DEnse GAs DISpersion) models the atmospheric dispersion of clouds
heavier than air.  It was written by Tom Spicer and Jerry Havens at the
University of Arkansas for the US Coast Guard and the Gas Research Institute,
released by EPA in 1989, and is still named in 49 CFR 193.2059 as an
acceptable means of determining LNG vapour dispersion exclusion zones.

``degali`` uses a verified Python reimplementation of DEGADIS 2.1 as its
foundation.  The goal is not a transliteration:
the numerics are modernised (SciPy integrators and root finders in place of
the bundled Runge-Kutta-Gill and Brent routines, CoolProp equations of state
in place of the 1989 correlations), while every step is validated against the
original Fortran running in the same repository.

Two backends are available throughout:

``legacy``
    Bit-faithful to DEGADIS 2.1, including its approximations and its
    round-off quirks.  Use this to demonstrate that the port is correct.
``coolprop``
    Real-fluid properties and accurate quadrature.  Use this for new work.

Getting started
---------------

.. code-block:: python

    from degali import run_steady

    profile, source = run_steady("B9.INP")
    print(profile.distance_to(0.05))   # metres to the LNG lower flammable limit
    print(profile.mass_above_lfl)      # kg

or from a shell::

    degali steady B9.INP
    degali transient B9T.INP --snapshot 60
    degali jet EX2.INO --bridge EX2.IN
"""

__version__ = "0.1.0"

from .run import (
    Receptor,
    SourceResult,
    TransientOutput,
    run_jet,
    run_jet_to_ground,
    run_source,
    run_steady,
    run_transient,
)
from .lh2 import (
    LH2CoupledResearchResult,
    LH2CrosswindHandoff,
    LH2ExpandedSource,
    LH2NearFieldResearchResult,
    handoff_lh2_near_field_to_crosswind,
    lh2_source_from_measured_throat,
    run_lh2_near_field_research,
    run_lh2_crosswind_research,
)

__all__ = [
    "__version__",
    "Receptor", "SourceResult", "TransientOutput",
    "run_source", "run_steady", "run_transient", "run_jet",
    "run_jet_to_ground",
    "LH2ExpandedSource", "LH2NearFieldResearchResult", "LH2CrosswindHandoff",
    "LH2CoupledResearchResult",
    "handoff_lh2_near_field_to_crosswind",
    "lh2_source_from_measured_throat", "run_lh2_near_field_research",
    "run_lh2_crosswind_research",
]
