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
    E34DynamicWindHydrogenHistory,
    E34EvaporationWindow,
    E34FillPrehistoryCandidate,
    E34FillPrehistoryPolicy,
    E34PoolHistory,
    E34SourceWindowCandidate,
    E34SourceWindowPolicy,
    E34WindHistory,
    evaporation_window,
    read_e34_dynamic_wind_hydrogen_history,
    read_e34_pool_history,
    read_e34_wind_history,
    source_only_evaporation_candidates,
    source_only_fill_prehistory_candidate,
)
from .open_channel_h2 import (
    OpenChannelHydrogenRun,
    ThresholdWindow,
    read_open_channel_h2_csv,
    threshold_window,
)
from .field_manifest_runtime import verify_runtime_manifest
from .open_channel_h2_obstacle import (
    ChannelSensorReceptor,
    FFIObstacleCase,
    OpenChannelH2ObstacleBenchmark,
    conditional_ffi_o1_crossbeam_assembly,
    conditional_ffi_o1_thin_top_crossbeam_assembly,
    load_open_channel_h2_obstacle_benchmark,
    load_registered_ffi_obstacle_geometry,
)
from .elvhys import (
    ElvhysTcsTest,
    ElvhysTimeSeries,
    read_elvhys_stream,
    read_elvhys_test,
)
from .usn_carpark import (
    UsnCarparkRun,
    UsnCarparkSensor,
    UsnThresholdMetrics,
    UsnWellMixedComparison,
    compare_well_mixed_baseline,
    integrated_positive_source_mass_g,
    read_usn_carpark_mat,
    threshold_metrics as usn_threshold_metrics,
    well_mixed_concentration_history,
)
from .h2safe import (
    CONCENTRATION_UNIT_STATUS as H2SAFE_CONCENTRATION_UNIT_STATUS,
    DATASET_DOI as H2SAFE_DATASET_DOI,
    DATASET_URL as H2SAFE_DATASET_URL,
    LAB1_VENT_BOUNDARIES as H2SAFE_LAB1_VENT_BOUNDARIES,
    LAB2_RECIRCULATION_FANS as H2SAFE_LAB2_RECIRCULATION_FANS,
    LAB2_SUPPLY_BOUNDARIES as H2SAFE_LAB2_SUPPLY_BOUNDARIES,
    PUBLISHED_VERTICAL_AXIS as H2SAFE_PUBLISHED_VERTICAL_AXIS,
    RELEASES as H2SAFE_RELEASES,
    SENSOR_CALIBRATION_STATUS as H2SAFE_SENSOR_CALIBRATION_STATUS,
    SOURCE_CLOCK_ALIGNMENT_STATUS as H2SAFE_SOURCE_CLOCK_ALIGNMENT_STATUS,
    H2SafeReleaseMetadata,
    H2SafeRecirculationFan,
    H2SafeRun,
    H2SafeSensor,
    H2SafeStlBounds,
    H2SafeVentBoundary,
    h2safe_raw_xyz_to_solver_xyz,
    read_h2safe_ascii_stl_bounds,
    read_h2safe_run,
    read_h2safe_step_product_names,
)
from .exjet import (
    sample_exjet_spatial_profile_predictions,
    score_exjet_profile_points,
)
from .hytunnel_ps import (
    HyTunnelConcentrationDataset,
    HyTunnelConcentrationPoint,
    HyTunnelConfiguration,
    HyTunnelFlowDataset,
    HyTunnelFlowPoint,
    read_hytunnel_concentration_workbook,
    read_hytunnel_flow_workbook,
)
from .thorney_ti021 import (
    THORNEY_TI021_EVIDENCE_SCHEMA,
    ThorneyTi021Dataset,
    Ti021ArcObservation,
    Ti021DoseSensor,
    read_thorney_ti021_workbook,
    thorney_ti021_evidence_record,
)
from .obstacle_flow_canonical import (
    POISEUILLE_CHANNEL_VERIFICATION_SCHEMA,
    D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA,
    PoiseuilleChannelVerificationConfig,
    PoiseuilleChannelVerificationResult,
    D3Q19PoiseuilleVerificationConfig,
    D3Q19PoiseuilleVerificationResult,
    verify_poiseuille_channel,
    verify_d3q19_poiseuille_channel,
)
from .obstacle_flow_convergence import (
    RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA,
    ObstacleReceptorConvergenceMetric,
    ResolvedObstacleGridCase,
    ResolvedObstacleGridConvergenceResult,
    resolved_obstacle_grid_convergence_record,
    run_resolved_obstacle_grid_convergence,
)
from .aij_case_h import (
    AIJ_CASE_H_DATASET_SCHEMA,
    AIJ_CASE_H_FLOW_COMPARISON_SCHEMA,
    AIJ_CASE_H_RECORD_DOI,
    AIJ_CASE_H_DECLARED_README_DOI,
    AijCaseHApproachPoint,
    AijCaseHObservation,
    AijCaseHDataset,
    AijCaseHFlowComparisonPoint,
    AijCaseHFlowComparison,
    read_aij_case_h,
    compare_aij_case_h_centerplane_flow,
    compare_aij_case_h_centerplane_flow_3d,
    aij_case_h_flow_comparison_record,
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
    "E34PoolHistory", "E34EvaporationWindow", "E34DynamicWindHydrogenHistory",
    "E34WindHistory", "E34SourceWindowPolicy", "E34SourceWindowCandidate",
    "E34FillPrehistoryPolicy", "E34FillPrehistoryCandidate",
    "read_e34_pool_history", "read_e34_dynamic_wind_hydrogen_history",
    "read_e34_wind_history", "evaporation_window",
    "source_only_evaporation_candidates", "source_only_fill_prehistory_candidate",
    "OpenChannelHydrogenRun", "ThresholdWindow", "read_open_channel_h2_csv",
    "threshold_window",
    "ElvhysTcsTest", "ElvhysTimeSeries", "read_elvhys_stream", "read_elvhys_test",
    "UsnCarparkRun", "UsnCarparkSensor", "UsnThresholdMetrics",
    "UsnWellMixedComparison", "read_usn_carpark_mat", "usn_threshold_metrics",
    "integrated_positive_source_mass_g", "well_mixed_concentration_history",
    "compare_well_mixed_baseline",
    "H2SAFE_DATASET_DOI", "H2SAFE_DATASET_URL", "H2SAFE_RELEASES",
    "H2SAFE_LAB1_VENT_BOUNDARIES", "H2SAFE_LAB2_SUPPLY_BOUNDARIES",
    "H2SAFE_LAB2_RECIRCULATION_FANS",
    "H2SAFE_PUBLISHED_VERTICAL_AXIS", "h2safe_raw_xyz_to_solver_xyz",
    "H2SAFE_CONCENTRATION_UNIT_STATUS", "H2SAFE_SOURCE_CLOCK_ALIGNMENT_STATUS",
    "H2SAFE_SENSOR_CALIBRATION_STATUS", "H2SafeReleaseMetadata",
    "H2SafeVentBoundary", "H2SafeRecirculationFan",
    "H2SafeSensor", "H2SafeRun", "H2SafeStlBounds", "read_h2safe_run",
    "read_h2safe_ascii_stl_bounds", "read_h2safe_step_product_names",
    "HyTunnelConfiguration", "HyTunnelConcentrationPoint",
    "HyTunnelConcentrationDataset", "HyTunnelFlowPoint", "HyTunnelFlowDataset",
    "read_hytunnel_concentration_workbook", "read_hytunnel_flow_workbook",
    "THORNEY_TI021_EVIDENCE_SCHEMA", "Ti021ArcObservation",
    "Ti021DoseSensor", "ThorneyTi021Dataset", "read_thorney_ti021_workbook",
    "thorney_ti021_evidence_record",
    "POISEUILLE_CHANNEL_VERIFICATION_SCHEMA",
    "D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA",
    "PoiseuilleChannelVerificationConfig",
    "PoiseuilleChannelVerificationResult", "verify_poiseuille_channel",
    "D3Q19PoiseuilleVerificationConfig", "D3Q19PoiseuilleVerificationResult",
    "verify_d3q19_poiseuille_channel",
    "RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA",
    "ObstacleReceptorConvergenceMetric", "ResolvedObstacleGridCase",
    "ResolvedObstacleGridConvergenceResult",
    "run_resolved_obstacle_grid_convergence",
    "resolved_obstacle_grid_convergence_record",
    "AIJ_CASE_H_DATASET_SCHEMA", "AIJ_CASE_H_FLOW_COMPARISON_SCHEMA",
    "AIJ_CASE_H_RECORD_DOI", "AIJ_CASE_H_DECLARED_README_DOI",
    "AijCaseHApproachPoint", "AijCaseHObservation", "AijCaseHDataset",
    "AijCaseHFlowComparisonPoint", "AijCaseHFlowComparison",
    "read_aij_case_h", "compare_aij_case_h_centerplane_flow",
    "compare_aij_case_h_centerplane_flow_3d",
    "aij_case_h_flow_comparison_record",
    "Comparison", "SeriesResult", "compare", "compare_series",
    "Statistics", "statistics", "table",
    "ComparisonCase", "ModelPrediction", "ModelComparisonReport", "compare_models",
    "ScreeningDecision", "evaluate_screening",
    "FFI_SOURCE_STATE_ENVELOPE_SCHEMA", "FfiSourceState", "FfiResidualRow",
    "FfiSourceStateEnvelopeCase", "FfiSourceStateEnvelope",
    "ffi_reference_provenance", "run_ffi_source_state_envelope",
    "ffi_source_state_envelope_report",
    "score_exjet_profile_points",
    "sample_exjet_spatial_profile_predictions",
]
