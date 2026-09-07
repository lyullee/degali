"""Validating degali: against the original Fortran, and against measurement.

Two distinct jobs live here.

:mod:`~degali.validation.reference` drives the original DEGADIS 2.1 Fortran,
built from source in this repository, and extracts its internal state at full
double precision.  That is how the port was shown to reproduce the model.

The rest compares the model against field trials.
:mod:`~degali.validation.rediphem` reads the Risoe REDIPHEM database of
heavy-gas releases, :mod:`~degali.validation.trialcase` turns a trial into a
DEGADIS deck while recording every assumption it had to make,
:mod:`~degali.validation.compare` pairs predictions with measurements, and
:mod:`~degali.validation.statistics` reduces the pairs to the measures the
published evaluations report.

The two are independent. Reproducing DEGADIS says nothing about whether
DEGADIS is right; measuring DEGADIS against trials says nothing about whether
this is DEGADIS. Both are needed.
"""

from .compare import Comparison, SeriesResult, compare, compare_series
from .rediphem import Trial, load
from .reference import PORTABILITY_PATCHES, TEST_CASES, Reference, ReferenceRun
from .statistics import Statistics, statistics, table
from .flashing import FlashResult, equivalent_source, plume_half_width
from .trialcase import (
    SMEDIS_SOURCES,
    SUBSTANCES,
    EquivalentSource,
    Substance,
    TrialCase,
    computed_source,
    jet_to_case,
    puff_to_case,
    to_case,
)

__all__ = [
    "Reference", "ReferenceRun", "TEST_CASES", "PORTABILITY_PATCHES",
    "Trial", "load",
    "Substance", "SUBSTANCES", "TrialCase", "to_case",
    "EquivalentSource", "SMEDIS_SOURCES", "computed_source", "jet_to_case",
    "puff_to_case",
    "FlashResult", "equivalent_source", "plume_half_width",
    "Comparison", "SeriesResult", "compare", "compare_series",
    "Statistics", "statistics", "table",
]
