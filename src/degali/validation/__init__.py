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
from .model_comparison import (
    ComparisonCase,
    ModelComparisonReport,
    ModelPrediction,
    compare_models,
)
from .screening_gate import ScreeningDecision, evaluate_screening
from .ffi_source_state import (
    FFI_SOURCE_STATE_ENVELOPE_SCHEMA,
    FfiResidualRow,
    FfiSourceState,
    FfiSourceStateEnvelope,
    FfiSourceStateEnvelopeCase,
    ffi_reference_provenance,
    ffi_source_state_envelope_report,
    run_ffi_source_state_envelope,
)
from .flashing import FlashResult, equivalent_source, plume_half_width
from .integration import (
    BranchDecision,
    ValidationEvidence,
    assess_validation_portfolio,
    validation_portfolio,
    validation_portfolio_dict,
)
from .sandia_pool import (
    PoolContourLowerBound,
    PoolContourScreen,
    read_pool_contour_lower_bounds,
    screen_contour_lower_bound,
    screen_pool_contours,
    screen_report,
)
from .preslhy_e31 import (
    E31AxisTemperatureHistory,
    E31ReleaseLineTemperatureHistory,
    E31AxisTemperatureEnvelope,
    E31AxisTemperatureComparison,
    E31BlowdownPressureComparison,
    E31VesselInventoryChange,
    E31InventoryReleaseTiming,
    E31HighPressureRun,
    E31PressureRun,
    E31IdealChokedSourceBound,
    E31NozzlePressureRise,
    E31PressureHistory,
    E31ValveInterval,
    axis_temperature_envelope,
    compare_axis_temperature_envelope,
    compare_blowdown_pressure_stations,
    vessel_inventory_change,
    vessel_inventory_release_timing,
    homogeneous_equilibrium_source_bound,
    ideal_choked_hydrogen_source_bound,
    nozzle_pressure_rise,
    read_e31_high_pressure_run,
    read_e31_pressure_run,
    valve_open_interval,
)
from .preslhy_e34 import (
    E34EvaporationWindow,
    E34PoolHistory,
    evaporation_window,
    read_e34_pool_history,
)
from .open_channel_h2 import (
    OpenChannelHydrogenRun,
    ThresholdWindow,
    read_open_channel_h2_csv,
    threshold_window,
)
from .elvhys import (
    ElvhysTcsTest,
    ElvhysTimeSeries,
    read_elvhys_stream,
    read_elvhys_test,
)
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
    "ValidationEvidence", "BranchDecision", "assess_validation_portfolio",
    "validation_portfolio", "validation_portfolio_dict",
    "PoolContourLowerBound", "PoolContourScreen", "read_pool_contour_lower_bounds",
    "screen_contour_lower_bound", "screen_pool_contours", "screen_report",
    "E31PressureHistory", "E31AxisTemperatureHistory", "E31ReleaseLineTemperatureHistory", "E31HighPressureRun", "E31PressureRun",
    "E31ValveInterval", "E31NozzlePressureRise", "E31AxisTemperatureEnvelope",
    "E31AxisTemperatureComparison",
    "E31BlowdownPressureComparison",
    "E31VesselInventoryChange",
    "E31InventoryReleaseTiming",
    "E31IdealChokedSourceBound",
    "read_e31_high_pressure_run", "read_e31_pressure_run", "valve_open_interval", "nozzle_pressure_rise",
    "axis_temperature_envelope",
    "compare_axis_temperature_envelope",
    "compare_blowdown_pressure_stations",
    "vessel_inventory_change",
    "vessel_inventory_release_timing",
    "homogeneous_equilibrium_source_bound",
    "ideal_choked_hydrogen_source_bound",
    "E34PoolHistory", "E34EvaporationWindow", "read_e34_pool_history",
    "evaporation_window",
    "OpenChannelHydrogenRun", "ThresholdWindow", "read_open_channel_h2_csv",
    "threshold_window",
    "ElvhysTcsTest", "ElvhysTimeSeries", "read_elvhys_stream", "read_elvhys_test",
    "Comparison", "SeriesResult", "compare", "compare_series",
    "Statistics", "statistics", "table",
    "ComparisonCase", "ModelPrediction", "ModelComparisonReport", "compare_models",
    "ScreeningDecision", "evaluate_screening",
    "FFI_SOURCE_STATE_ENVELOPE_SCHEMA", "FfiSourceState", "FfiResidualRow",
    "FfiSourceStateEnvelopeCase", "FfiSourceStateEnvelope",
    "ffi_reference_provenance", "run_ffi_source_state_envelope",
    "ffi_source_state_envelope_report",
]
