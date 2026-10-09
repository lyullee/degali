"""Failure-safe field screening workflow around the reduced-order operator."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from itertools import product
import math
from typing import TYPE_CHECKING

from .field_contracts import (
    BoundedValue,
    FieldApplicability,
    FieldCoordinateReference,
    FieldScenario,
    FieldValidationEvidence,
    SensorModel,
    assess_field_applicability,
)
from .field_distributed_source import FieldDistributedVapourSource
from .field_source_io import FieldAtmosphericSourceSchedule
from .field_history import MeasuredHistoryQualityAssessment
from .field_geometry import (
    FieldObstacleGeometryUncertainty,
    WindPlaneObstacleProjection,
    project_cuboid_to_wind_plane,
)
from .field_jet import FieldJetScalarHandoff, jet_scalar_handoff_location
from .field_meteorology import (
    FieldStabilityAlternatives,
    FieldWindHistory,
    StabilityScalarMixingClosure,
)
from .field_lh2 import (
    FieldLH2SourcePreparation,
    build_lh2_saturation_table_for_release,
    prepare_field_lh2_flash,
)
from .field_observation import FieldSensorTrace, apply_sensor_model
from .semi_fv_obstacle import (
    DistributedScalarSource,
    RectangularObstacle2D,
    SemiFVConfig,
    SemiFVReceptor,
    SemiFVRefinementStudy,
    SemiFVResult,
    SourceRateSchedule,
    run_semi_fv_refinement_study as run_transport_refinement_study,
    solve_semi_fv_obstacle,
)
from .site_geometry import AxisAlignedCuboid, OrientedCuboid, WindFrame
from .finite_release import SteadyWindApplicability

if TYPE_CHECKING:
    from .lh2_property_table import LH2SaturationTable


_FIELD_MASS_RESIDUAL_REL_TOL = 1.0e-8
_FIELD_MASS_RESIDUAL_ABS_TOL_KG = 1.0e-12


@dataclass(frozen=True)
class FieldSensorDeployment:
    """One labelled field detector evaluated on the shared local transport run."""

    label: str
    sensor: SensorModel

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("field sensor deployment label must be non-empty")
        if not isinstance(self.sensor, SensorModel):
            raise TypeError("field sensor deployment sensor must be a SensorModel")


@dataclass(frozen=True)
class FieldSensorDeploymentResult:
    """A sensor trace or an explicit reason it was withheld from the 2-D plane."""

    label: str
    sensor: SensorModel
    trace: FieldSensorTrace | None
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("field sensor result label must be non-empty")
        if not isinstance(self.sensor, SensorModel):
            raise TypeError("field sensor result sensor must be a SensorModel")
        if self.trace is not None and not isinstance(self.trace, FieldSensorTrace):
            raise TypeError("field sensor result trace must be a FieldSensorTrace or None")
        if not isinstance(self.warnings, tuple):
            raise TypeError("field sensor result warnings must be a tuple")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("field sensor result warnings must be non-empty strings")
        if self.trace is None and not self.warnings:
            raise ValueError("withheld field sensor result requires at least one warning")


@dataclass(frozen=True)
class FieldSemiFVRequest:
    """All explicit inputs for one local direct-vapour field screening run.

    ``direct_vapour_schedule`` is a post-flash atmospheric hydrogen-vapour
    source record. It is never inferred from the upstream ``ReleaseSource``;
    when supplied it replaces the nominal flash vapour *rate* in local
    transport while retaining the nominal flash only as an auditable source
    contract check.

    ``lh2_validation_available`` is only a compatibility flag.  A true value
    is accepted by construction for legacy Python callers but is fail-closed
    at screening time unless ``validation_evidence`` fingerprints the matched
    external dataset.

    The three ambient scalar values are nominal environmental boundaries.
    Their optional ``*_uncertainty`` companions are explicit deterministic
    bounds; the field envelope varies them through flash and sensor conversion
    rather than treating them as a probability distribution.
    """

    scenario: FieldScenario
    transport: SemiFVConfig = field(default_factory=SemiFVConfig)
    coordinate_reference: FieldCoordinateReference | None = None
    obstacle: AxisAlignedCuboid | OrientedCuboid | None = None
    obstacles: tuple[AxisAlignedCuboid | OrientedCuboid, ...] = ()
    obstacle_geometry_uncertainty: tuple[FieldObstacleGeometryUncertainty, ...] = ()
    ambient_temperature_k: float = 295.0
    ambient_pressure_pa: float = 101325.0
    ambient_air_density_kg_m3: float = 1.2
    lateral_sensor_tolerance_m: float | None = None
    lateral_source_tolerance_m: float | None = None
    lh2_validation_available: bool = False
    property_table: "LH2SaturationTable | None" = None
    direct_vapour_schedule: SourceRateSchedule | None = None
    direct_vapour_warnings: tuple[str, ...] = ()
    direct_vapour_effective_area_m2: float | None = None
    stability_mixing_closure: StabilityScalarMixingClosure | None = None
    stability_alternatives: FieldStabilityAlternatives | None = None
    stability_uncertainty_resolved: bool = False
    phase_routing_uncertainty_resolved: bool = False
    wind_history: FieldWindHistory | None = None
    jet_scalar_handoff: FieldJetScalarHandoff | None = None
    post_release_duration_s: float = 0.0
    distributed_vapour_sources: tuple[FieldDistributedVapourSource, ...] = ()
    measured_history_quality: MeasuredHistoryQualityAssessment | None = None
    measured_history_provenance: tuple[tuple[str, str], ...] = ()
    measured_history_source_uncertainty_resolved: bool = False
    sensor_deployments: tuple[FieldSensorDeployment, ...] = ()
    validation_evidence: FieldValidationEvidence | None = None
    ambient_temperature_uncertainty_k: BoundedValue | None = None
    ambient_pressure_uncertainty_pa: BoundedValue | None = None
    ambient_air_density_uncertainty_kg_m3: BoundedValue | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.scenario, FieldScenario):
            raise TypeError("scenario must be a FieldScenario")
        if not isinstance(self.transport, SemiFVConfig):
            raise TypeError("transport must be a SemiFVConfig")
        for name, value in (
            ("obstacles", self.obstacles),
            ("obstacle_geometry_uncertainty", self.obstacle_geometry_uncertainty),
            ("direct_vapour_warnings", self.direct_vapour_warnings),
            ("distributed_vapour_sources", self.distributed_vapour_sources),
            ("sensor_deployments", self.sensor_deployments),
        ):
            if not isinstance(value, tuple):
                raise TypeError(f"{name} must be a tuple")
        if self.coordinate_reference is not None and not isinstance(
            self.coordinate_reference, FieldCoordinateReference
        ):
            raise TypeError("coordinate_reference must be a FieldCoordinateReference or None")
        if self.validation_evidence is not None and not isinstance(
            self.validation_evidence, FieldValidationEvidence
        ):
            raise TypeError("validation_evidence must be a FieldValidationEvidence or None")
        for name, value in {
            "ambient_temperature_k": self.ambient_temperature_k,
            "ambient_pressure_pa": self.ambient_pressure_pa,
            "ambient_air_density_kg_m3": self.ambient_air_density_kg_m3,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        ambient_bounds = {
            "ambient_temperature_uncertainty_k": (
                self.ambient_temperature_uncertainty_k,
                self.ambient_temperature_k,
                "ambient temperature",
                "K",
            ),
            "ambient_pressure_uncertainty_pa": (
                self.ambient_pressure_uncertainty_pa,
                self.ambient_pressure_pa,
                "ambient pressure",
                "Pa",
            ),
            "ambient_air_density_uncertainty_kg_m3": (
                self.ambient_air_density_uncertainty_kg_m3,
                self.ambient_air_density_kg_m3,
                "ambient air density",
                "kg/m3",
            ),
        }
        for name, (bound, nominal, label, unit) in ambient_bounds.items():
            if bound is None:
                continue
            if not isinstance(bound, BoundedValue):
                raise TypeError(f"{name} must be a BoundedValue or None")
            if bound.unit != unit:
                raise ValueError(f"{name} must use unit {unit!r}")
            if bound.lower <= 0.0:
                raise ValueError(f"{label} uncertainty lower bound must be positive")
            if not math.isclose(
                bound.nominal, float(nominal), rel_tol=1.0e-9, abs_tol=1.0e-12,
            ):
                raise ValueError(f"{name}.nominal must match {label} nominal value")
        if self.lateral_sensor_tolerance_m is not None and (
            isinstance(self.lateral_sensor_tolerance_m, bool)
            or not math.isfinite(float(self.lateral_sensor_tolerance_m))
            or self.lateral_sensor_tolerance_m < 0.0
        ):
            raise ValueError("lateral_sensor_tolerance_m must be finite and non-negative")
        if self.lateral_source_tolerance_m is not None and (
            isinstance(self.lateral_source_tolerance_m, bool)
            or not math.isfinite(float(self.lateral_source_tolerance_m))
            or self.lateral_source_tolerance_m < 0.0
        ):
            raise ValueError("lateral_source_tolerance_m must be finite and non-negative")
        if (
            isinstance(self.post_release_duration_s, bool)
            or not math.isfinite(float(self.post_release_duration_s))
            or self.post_release_duration_s < 0.0
        ):
            raise ValueError("post_release_duration_s must be finite and non-negative")
        if self.post_release_duration_s > 0.0 and self.scenario.temporal_mode != "transient":
            raise ValueError("post_release_duration_s requires a transient field scenario")
        if self.obstacle is not None and self.obstacles:
            raise ValueError("declare either obstacle or obstacles, not both")
        cuboid_types = (AxisAlignedCuboid, OrientedCuboid)
        if self.obstacle is not None and not isinstance(self.obstacle, cuboid_types):
            raise TypeError("obstacle must be an AxisAlignedCuboid or OrientedCuboid")
        if any(not isinstance(item, cuboid_types) for item in self.obstacles):
            raise TypeError("obstacles must contain only AxisAlignedCuboid or OrientedCuboid values")
        labels = tuple(item.label for item in self.obstacles)
        if len(set(labels)) != len(labels):
            raise ValueError("declared field obstacle labels must be unique")
        declared_obstacles = ((self.obstacle,) if self.obstacle is not None else self.obstacles)
        if any(
            not isinstance(item, FieldObstacleGeometryUncertainty)
            for item in self.obstacle_geometry_uncertainty
        ):
            raise TypeError(
                "obstacle_geometry_uncertainty must contain only "
                "FieldObstacleGeometryUncertainty values"
            )
        uncertainty_labels = tuple(
            item.obstacle.label for item in self.obstacle_geometry_uncertainty
        )
        if len(set(uncertainty_labels)) != len(uncertainty_labels):
            raise ValueError("obstacle geometry uncertainty labels must be unique")
        declared_labels = tuple(item.label for item in declared_obstacles)
        if any(label not in declared_labels for label in uncertainty_labels):
            raise ValueError(
                "obstacle geometry uncertainty must refer to a declared obstacle"
            )
        if any(
            item.obstacle != next(
                obstacle for obstacle in declared_obstacles
                if obstacle.label == item.obstacle.label
            )
            for item in self.obstacle_geometry_uncertainty
        ):
            raise ValueError(
                "obstacle geometry uncertainty nominal geometry must match the declared obstacle"
            )
        if self.direct_vapour_schedule is not None:
            if not isinstance(self.direct_vapour_schedule, SourceRateSchedule):
                raise TypeError("direct_vapour_schedule must be a SourceRateSchedule")
            if self.direct_vapour_schedule.source_id.strip().lower() == "declared":
                raise ValueError(
                    "direct_vapour_schedule requires an explicit source_id provenance"
                )
            if self.scenario.temporal_mode != "transient" or self.scenario.source.duration_s is None:
                raise ValueError("direct_vapour_schedule requires a transient field scenario")
            if (
                self.scenario.source.duration_uncertainty is not None
                and not self.scenario.source.duration_uncertainty.is_exact
            ):
                raise ValueError(
                    "direct_vapour_schedule cannot be paired with unresolved "
                    "duration_uncertainty; use a time-aligned source envelope"
                )
            if not math.isclose(
                self.direct_vapour_schedule.duration_s,
                self.scenario.source.duration_s,
                rel_tol=1.0e-9, abs_tol=1.0e-12,
            ):
                raise ValueError("direct_vapour_schedule duration must match source.duration_s")
            self.direct_vapour_schedule.require_zero_endpoint(
                "direct_vapour_schedule"
            )
        if any(not isinstance(warning, str) or not warning.strip() for warning in self.direct_vapour_warnings):
            raise ValueError("direct_vapour_warnings must contain non-empty strings")
        if self.direct_vapour_effective_area_m2 is not None and (
            isinstance(self.direct_vapour_effective_area_m2, bool)
            or not math.isfinite(float(self.direct_vapour_effective_area_m2))
            or self.direct_vapour_effective_area_m2 <= 0.0
        ):
            raise ValueError("direct_vapour_effective_area_m2 must be positive and finite when supplied")
        if self.stability_mixing_closure is not None:
            if not isinstance(self.stability_mixing_closure, StabilityScalarMixingClosure):
                raise TypeError("stability_mixing_closure must be a StabilityScalarMixingClosure")
            # Resolve this when the request is constructed, rather than after
            # source preparation has already completed.  There is deliberately
            # no default or nearest-class fallback in the closure.
            self.stability_mixing_closure.selected_diffusivity_m2_s(
                self.scenario.weather.stability
            )
        if self.stability_alternatives is not None:
            if not isinstance(self.stability_alternatives, FieldStabilityAlternatives):
                raise TypeError("stability_alternatives must be a FieldStabilityAlternatives")
            if self.stability_mixing_closure is None:
                raise ValueError(
                    "stability alternatives require an explicit stability_mixing_closure"
                )
            for stability in self.stability_alternatives.classes_for(
                self.scenario.weather.stability
            ):
                self.stability_mixing_closure.selected_diffusivity_m2_s(stability)
        if not isinstance(self.stability_uncertainty_resolved, bool):
            raise TypeError("stability_uncertainty_resolved must be boolean")
        if self.stability_uncertainty_resolved and self.stability_alternatives is None:
            raise ValueError(
                "stability_uncertainty_resolved requires declared stability alternatives"
            )
        if not isinstance(self.phase_routing_uncertainty_resolved, bool):
            raise TypeError("phase_routing_uncertainty_resolved must be boolean")
        if self.wind_history is not None and not isinstance(self.wind_history, FieldWindHistory):
            raise TypeError("wind_history must be a FieldWindHistory")
        if self.measured_history_quality is not None:
            if not isinstance(self.measured_history_quality, MeasuredHistoryQualityAssessment):
                raise TypeError("measured_history_quality must be MeasuredHistoryQualityAssessment or None")
            if not self.measured_history_quality.approved:
                raise ValueError("measured_history_quality must be approved when attached to a field request")
            if self.direct_vapour_schedule is None:
                raise ValueError("measured_history_quality requires a direct_vapour_schedule")
        if not isinstance(self.measured_history_source_uncertainty_resolved, bool):
            raise TypeError("measured_history_source_uncertainty_resolved must be boolean")
        if (
            self.measured_history_source_uncertainty_resolved
            and self.measured_history_quality is None
        ):
            raise ValueError(
                "measured_history_source_uncertainty_resolved requires measured_history_quality"
            )
        provenance = tuple(self.measured_history_provenance)
        if any(
            not isinstance(key, str) or not key.strip()
            or not isinstance(value, str) or not value.strip()
            for key, value in provenance
        ):
            raise ValueError("measured_history_provenance must contain non-empty string key/value pairs")
        if len({key for key, _value in provenance}) != len(provenance):
            raise ValueError("measured_history_provenance keys must be unique")
        if provenance and self.measured_history_quality is None:
            raise ValueError("measured_history_provenance requires measured_history_quality")
        object.__setattr__(self, "measured_history_provenance", provenance)
        if any(not isinstance(item, FieldSensorDeployment) for item in self.sensor_deployments):
            raise TypeError("sensor_deployments must contain only FieldSensorDeployment values")
        sensor_labels = tuple(item.label for item in self.sensor_deployments)
        if len(set(sensor_labels)) != len(sensor_labels):
            raise ValueError("field sensor deployment labels must be unique")
        if "field_sensor" in sensor_labels:
            raise ValueError("field_sensor is reserved for FieldScenario.sensor")
        if self.jet_scalar_handoff is not None:
            if not isinstance(self.jet_scalar_handoff, FieldJetScalarHandoff):
                raise TypeError("jet_scalar_handoff must be a FieldJetScalarHandoff")
            if self.direct_vapour_schedule is not None:
                raise ValueError(
                    "jet_scalar_handoff cannot be combined with a direct_vapour_schedule"
                )
        if any(
            not isinstance(item, FieldDistributedVapourSource)
            for item in self.distributed_vapour_sources
        ):
            raise TypeError(
                "distributed_vapour_sources must contain only FieldDistributedVapourSource values"
            )
        if self.phase_routing_uncertainty_resolved and not any(
            source.source_kind == "pool_vapour"
            for source in self.distributed_vapour_sources
        ):
            raise ValueError(
                "phase_routing_uncertainty_resolved requires a pool_vapour distributed source"
            )
        labels = tuple(item.label for item in self.distributed_vapour_sources)
        if len(set(labels)) != len(labels) or "primary" in labels:
            raise ValueError("field distributed vapour source labels must be unique and not 'primary'")
        source_duration = (
            self.scenario.source.duration_s
            if self.scenario.temporal_mode == "transient"
            else self.transport.duration_s
        )
        assert source_duration is not None
        total_duration = source_duration + self.post_release_duration_s
        if any(item.schedule.duration_s > total_duration + 1.0e-12
               for item in self.distributed_vapour_sources):
            raise ValueError("field distributed vapour source schedule duration exceeds transport duration")

    def ambient_uncertainty_fields(self) -> dict[str, BoundedValue]:
        """Return ambient boundary bounds that must be propagated as field corners."""
        values = {
            "ambient_temperature_k": self.ambient_temperature_uncertainty_k,
            "ambient_pressure_pa": self.ambient_pressure_uncertainty_pa,
            "ambient_air_density_kg_m3": self.ambient_air_density_uncertainty_kg_m3,
        }
        return {
            name: bound for name, bound in values.items() if bound is not None
        }


@dataclass(frozen=True)
class FieldSemiFVScreeningResult:
    """Structured result; a blocked run contains no transport field."""

    request: FieldSemiFVRequest
    applicability: FieldApplicability
    source_preparation: FieldLH2SourcePreparation | None
    wind_frame: WindFrame | None
    obstacle_projection: WindPlaneObstacleProjection | None
    transport: SemiFVResult | None
    sensor_trace: FieldSensorTrace | None
    wind_history_assessment: SteadyWindApplicability | None = None
    obstacle_projections: tuple[WindPlaneObstacleProjection, ...] = ()
    sensor_results: tuple[FieldSensorDeploymentResult, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("field screening result request must be a FieldSemiFVRequest")
        if not isinstance(self.applicability, FieldApplicability):
            raise TypeError("field screening result applicability must be a FieldApplicability")
        if self.source_preparation is not None and not isinstance(
            self.source_preparation, FieldLH2SourcePreparation
        ):
            raise TypeError(
                "field screening result source_preparation must be a "
                "FieldLH2SourcePreparation or None"
            )
        if self.wind_frame is not None and not isinstance(self.wind_frame, WindFrame):
            raise TypeError("field screening result wind_frame must be a WindFrame or None")
        if self.obstacle_projection is not None and not isinstance(
            self.obstacle_projection, WindPlaneObstacleProjection
        ):
            raise TypeError(
                "field screening result obstacle_projection must be a "
                "WindPlaneObstacleProjection or None"
            )
        if self.transport is not None and not isinstance(self.transport, SemiFVResult):
            raise TypeError("field screening result transport must be a SemiFVResult or None")
        if self.sensor_trace is not None and not isinstance(self.sensor_trace, FieldSensorTrace):
            raise TypeError(
                "field screening result sensor_trace must be a FieldSensorTrace or None"
            )
        if self.wind_history_assessment is not None and not isinstance(
            self.wind_history_assessment, SteadyWindApplicability
        ):
            raise TypeError(
                "field screening result wind_history_assessment must be a "
                "SteadyWindApplicability or None"
            )
        if not isinstance(self.obstacle_projections, tuple):
            raise TypeError("field screening result obstacle_projections must be a tuple")
        if not isinstance(self.sensor_results, tuple):
            raise TypeError("field screening result sensor_results must be a tuple")
        if any(
            not isinstance(item, WindPlaneObstacleProjection)
            for item in self.obstacle_projections
        ):
            raise TypeError(
                "field screening result obstacle_projections must contain only "
                "WindPlaneObstacleProjection values"
            )
        if any(
            not isinstance(item, FieldSensorDeploymentResult)
            for item in self.sensor_results
        ):
            raise TypeError(
                "field screening result sensor_results must contain only "
                "FieldSensorDeploymentResult values"
            )
        labels = tuple(item.label for item in self.sensor_results)
        if len(set(labels)) != len(labels):
            raise ValueError("field screening result sensor labels must be unique")
        if self.transport is None and (
            self.sensor_trace is not None
            or any(item.trace is not None for item in self.sensor_results)
        ):
            raise ValueError(
                "field screening result sensor traces require a transport result"
            )

    @property
    def completed(self) -> bool:
        return self.transport is not None and self.applicability.status != "blocked"


@dataclass(frozen=True)
class FieldSemiFVEnvelopeCase:
    """One deterministic field-input corner and its screening result."""

    values: tuple[tuple[str, float | str], ...]
    result: FieldSemiFVScreeningResult

    def __post_init__(self) -> None:
        if not isinstance(self.values, tuple):
            raise TypeError("field envelope corner values must be a tuple")
        if not isinstance(self.result, FieldSemiFVScreeningResult):
            raise TypeError(
                "field envelope case result must be a FieldSemiFVScreeningResult"
            )
        seen: set[str] = set()
        for key, value in self.values:
            if not isinstance(key, str) or not key.strip():
                raise ValueError("field envelope corner keys must be non-empty strings")
            if key in seen:
                raise ValueError("field envelope corner keys must be unique")
            seen.add(key)
            if isinstance(value, bool):
                raise TypeError("field envelope corner values cannot be boolean")
            if isinstance(value, (int, float)):
                if not math.isfinite(float(value)):
                    raise ValueError("field envelope corner numeric values must be finite")
            elif not isinstance(value, str) or not value.strip():
                raise TypeError(
                    "field envelope corner values must be finite numbers or non-empty strings"
                )

    @property
    def corner(self) -> dict[str, float | str]:
        return dict(self.values)


@dataclass(frozen=True)
class FieldSemiFVEnvelope:
    """Untruncated deterministic envelope over every declared input corner."""

    request: FieldSemiFVRequest
    cases: tuple[FieldSemiFVEnvelopeCase, ...]
    property_table_used: bool
    warnings: tuple[str, ...] = ()
    atmospheric_source_schedule: FieldAtmosphericSourceSchedule | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("field envelope request must be a FieldSemiFVRequest")
        if not isinstance(self.cases, tuple):
            raise TypeError("field uncertainty envelope cases must be a tuple")
        if not self.cases:
            raise ValueError("field uncertainty envelope needs at least one case")
        if any(not isinstance(item, FieldSemiFVEnvelopeCase) for item in self.cases):
            raise TypeError(
                "field uncertainty envelope cases must contain only "
                "FieldSemiFVEnvelopeCase values"
            )
        selections = tuple(tuple(case.values) for case in self.cases)
        if len(set(selections)) != len(selections):
            raise ValueError(
                "field uncertainty envelope corner selections must be unique"
            )
        if not isinstance(self.property_table_used, bool):
            raise TypeError("field uncertainty envelope property_table_used must be boolean")
        if not isinstance(self.warnings, tuple):
            raise TypeError("field envelope warnings must be a tuple")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("field uncertainty envelope warnings must be non-empty strings")
        if self.atmospheric_source_schedule is not None and not isinstance(
            self.atmospheric_source_schedule, FieldAtmosphericSourceSchedule
        ):
            raise TypeError(
                "field uncertainty envelope atmospheric_source_schedule must be a "
                "FieldAtmosphericSourceSchedule or None"
            )

    @property
    def completed_case_count(self) -> int:
        return sum(case.result.completed for case in self.cases)

    @property
    def blocked_case_count(self) -> int:
        return sum(case.result.applicability.status == "blocked" for case in self.cases)

    @property
    def selections(self) -> tuple[tuple[tuple[str, float | str], ...], ...]:
        """Return every actual selection, including categorical stability cases."""
        return tuple(case.values for case in self.cases)


@dataclass(frozen=True)
class FieldSensorArrayUncertaintyCase:
    """One calibration corner post-processed from a shared true scalar field."""

    values: tuple[tuple[str, float], ...]
    sensor_results: tuple[FieldSensorDeploymentResult, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.values, tuple):
            raise TypeError("sensor-array corner values must be a tuple of key/value pairs")
        seen: set[str] = set()
        for item in self.values:
            if not isinstance(item, tuple) or len(item) != 2:
                raise TypeError("sensor-array corner entries must be 2-item tuples")
            key, value = item
            if not isinstance(key, str) or not key.strip():
                raise ValueError("sensor-array corner keys must be non-empty strings")
            if key in seen:
                raise ValueError("sensor-array corner keys must be unique")
            seen.add(key)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError("sensor-array corner values must be numeric")
            if not math.isfinite(float(value)):
                raise ValueError("sensor-array corner values must be finite")
        if not self.sensor_results:
            raise ValueError("sensor-array uncertainty case needs at least one sensor result")
        if any(not isinstance(item, FieldSensorDeploymentResult) for item in self.sensor_results):
            raise TypeError(
                "sensor-array uncertainty case results must contain only "
                "FieldSensorDeploymentResult values"
            )
        labels = tuple(item.label for item in self.sensor_results)
        if len(set(labels)) != len(labels):
            raise ValueError("sensor-array uncertainty case sensor labels must be unique")

    @property
    def corner(self) -> dict[str, float]:
        return dict(self.values)


@dataclass(frozen=True)
class FieldSensorArrayUncertaintyEnvelope:
    """Calibration uncertainty envelope without re-solving unchanged transport."""

    screening: FieldSemiFVScreeningResult
    cases: tuple[FieldSensorArrayUncertaintyCase, ...]
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError(
                "sensor-array uncertainty envelope screening must be a "
                "FieldSemiFVScreeningResult"
            )
        if any(
            not isinstance(item, FieldSensorArrayUncertaintyCase)
            for item in self.cases
        ):
            raise TypeError(
                "sensor-array uncertainty envelope cases must contain only "
                "FieldSensorArrayUncertaintyCase values"
            )
        if self.screening.completed and not self.cases:
            raise ValueError(
                "completed sensor-array uncertainty envelope requires at least one case"
            )
        corner_keys = tuple(tuple(item.values) for item in self.cases)
        if len(set(corner_keys)) != len(corner_keys):
            raise ValueError("sensor-array uncertainty envelope corner selections must be unique")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("sensor-array uncertainty envelope warnings must be non-empty strings")

    @property
    def completed_case_count(self) -> int:
        return len(self.cases) if self.screening.completed else 0


@dataclass(frozen=True)
class FieldSemiFVRefinementResult:
    """Numerical-refinement audit tied to one declared field screening case."""

    screening: FieldSemiFVScreeningResult
    study: SemiFVRefinementStudy | None
    applicability: FieldApplicability

    def __post_init__(self) -> None:
        if not isinstance(self.screening, FieldSemiFVScreeningResult):
            raise TypeError(
                "field refinement screening must be a FieldSemiFVScreeningResult"
            )
        if self.study is not None and not isinstance(self.study, SemiFVRefinementStudy):
            raise TypeError(
                "field refinement study must be a SemiFVRefinementStudy or None"
            )
        if not isinstance(self.applicability, FieldApplicability):
            raise TypeError("field refinement applicability must be a FieldApplicability")

    @property
    def completed(self) -> bool:
        return self.study is not None and self.applicability.status != "blocked"


def _merge_applicability(
    *items: FieldApplicability,
    warnings: tuple[str, ...] = (),
) -> FieldApplicability:
    reasons = tuple(reason for item in items for reason in item.reasons)
    all_warnings = tuple(warning for item in items for warning in item.warnings) + warnings
    if reasons:
        return FieldApplicability("blocked", reasons, all_warnings, False)
    status = "conditional" if (
        any(item.status == "conditional" for item in items) or all_warnings
    ) else "accepted"
    return FieldApplicability(
        status, (), all_warnings,
        all(item.uncertainty_complete for item in items),
    )


def _projection_applicability(
    projections: tuple[WindPlaneObstacleProjection, ...],
) -> FieldApplicability:
    if not projections:
        return FieldApplicability("accepted", uncertainty_complete=True)
    reasons = tuple(reason for projection in projections for reason in projection.reasons)
    warnings = tuple(warning for projection in projections for warning in projection.warnings)
    if reasons:
        return FieldApplicability("blocked", reasons=reasons, warnings=warnings)
    status = "conditional" if any(
        projection.status == "conditional" for projection in projections
    ) else "accepted"
    return FieldApplicability(
        status, warnings=warnings, uncertainty_complete=True,
    )


def _transport_conservation_applicability(
    transport: SemiFVResult,
) -> FieldApplicability:
    """Fail closed when local scalar transport loses too much mass.

    The semi-FV solver reports its worst-step inventory residual. A field
    result must not silently pass that diagnostic through to an operational
    screen, especially when obstacle routing or multiple source schedules are
    active. This numerical gate is not a physical uncertainty interval; the
    full residual remains in the transport diagnostics for audit.
    """
    diagnostics = transport.diagnostics
    residual = float(diagnostics.maximum_mass_residual_kg)
    final_residual = float(diagnostics.final_mass_residual_kg)
    source_ledger_residual = float(diagnostics.source_mass_ledger_residual_kg)
    source_schedule_residual = float(
        diagnostics.maximum_source_mass_schedule_residual_kg
    )
    scale = max(
        abs(float(diagnostics.mass_injected_kg)),
        abs(float(diagnostics.mass_domain_kg)) + abs(float(diagnostics.mass_outflow_kg)),
        _FIELD_MASS_RESIDUAL_ABS_TOL_KG,
    )
    tolerance = max(
        _FIELD_MASS_RESIDUAL_ABS_TOL_KG,
        _FIELD_MASS_RESIDUAL_REL_TOL * scale,
    )
    ledger_scale = max(
        abs(float(diagnostics.mass_injected_kg)),
        abs(float(diagnostics.mass_domain_kg)) + abs(float(diagnostics.mass_outflow_kg)),
        _FIELD_MASS_RESIDUAL_ABS_TOL_KG,
    )
    ledger_tolerance = max(
        _FIELD_MASS_RESIDUAL_ABS_TOL_KG,
        _FIELD_MASS_RESIDUAL_REL_TOL * ledger_scale,
    )
    if not math.isfinite(residual) or residual > tolerance:
        return FieldApplicability(
            "blocked",
            reasons=(
                "semi-FV transport mass-conservation residual "
                f"{residual:.6g} kg exceeds the numerical gate {tolerance:.6g} kg",
            ),
            warnings=(
                "transport diagnostics are retained, but the field result is withheld until "
                "the source/obstacle discretisation closes its inventory",
            ),
            uncertainty_complete=False,
        )
    if not math.isfinite(final_residual) or abs(final_residual) > tolerance:
        return FieldApplicability(
            "blocked",
            reasons=(
                "semi-FV final mass-conservation residual "
                f"{final_residual:.6g} kg exceeds the numerical gate {tolerance:.6g} kg",
            ),
            warnings=(
                "transport final inventory does not close against injected mass and outflow; "
                "the field result is withheld until the transport ledger closes",
            ),
            uncertainty_complete=False,
        )
    if not math.isfinite(source_ledger_residual) or abs(source_ledger_residual) > ledger_tolerance:
        return FieldApplicability(
            "blocked",
            reasons=(
                "semi-FV source-mass ledger residual "
                f"{source_ledger_residual:.6g} kg exceeds the numerical gate "
                f"{ledger_tolerance:.6g} kg",
            ),
            warnings=(
                "transport diagnostics are retained, but the field result is withheld until "
                "the per-source injection ledger closes against total injected mass",
            ),
            uncertainty_complete=False,
        )
    if not math.isfinite(source_schedule_residual) or abs(source_schedule_residual) > ledger_tolerance:
        return FieldApplicability(
            "blocked",
            reasons=(
                "semi-FV source-mass schedule residual "
                f"{source_schedule_residual:.6g} kg exceeds the numerical gate "
                f"{ledger_tolerance:.6g} kg",
            ),
            warnings=(
                "transport diagnostics are retained, but the field result is withheld until "
                "each source schedule closes against the injected-mass ledger",
            ),
            uncertainty_complete=False,
        )
    return FieldApplicability("accepted", uncertainty_complete=True)


def _declared_global_obstacles(
    request: FieldSemiFVRequest,
) -> tuple[AxisAlignedCuboid | OrientedCuboid, ...]:
    if request.obstacle is not None:
        return (request.obstacle,)
    return request.obstacles


def _local_obstacles_from_projections(
    projections: tuple[WindPlaneObstacleProjection, ...],
) -> RectangularObstacle2D | tuple[RectangularObstacle2D, ...] | None:
    local = tuple(
        projection.obstacle for projection in projections
        if projection.obstacle is not None
    )
    if not local:
        return None
    return local[0] if len(local) == 1 else local


def _wind_history_applicability(
    request: FieldSemiFVRequest,
) -> tuple[FieldApplicability, SteadyWindApplicability | None]:
    """Gate a fixed-wind-plane calculation against a measured wind record."""
    history = request.wind_history
    if history is None:
        return FieldApplicability("accepted", uncertainty_complete=True), None
    duration = (
        request.scenario.source.duration_s
        if request.scenario.temporal_mode == "transient"
        else request.transport.duration_s
    )
    assert duration is not None
    try:
        assessment = history.assess(duration + request.post_release_duration_s)
    except ValueError as error:
        return FieldApplicability(
            "blocked",
            reasons=(f"measured wind-history gate cannot be evaluated: {error}",),
        ), None
    if not assessment.applicable:
        return FieldApplicability(
            "blocked",
            reasons=tuple(
                f"measured wind-history gate failed: {reason}"
                for reason in assessment.reasons
                if reason != "wind_window_within_declared_steady_limits"
            ),
        ), assessment
    try:
        _record_speed, _record_direction, speed_difference, direction_difference = (
            history.nominal_weather_difference(request.scenario.weather)
        )
    except ValueError as error:
        return FieldApplicability(
            "blocked",
            reasons=(
                "measured wind-history nominal consistency cannot be evaluated: "
                f"{error}",
            ),
        ), assessment
    mismatch_reasons = []
    if speed_difference > history.maximum_nominal_speed_relative_difference:
        mismatch_reasons.append("wind_speed_nominal_difference_exceeds_declared_limit")
    if direction_difference > history.maximum_nominal_direction_difference_deg:
        mismatch_reasons.append("wind_direction_nominal_difference_exceeds_declared_limit")
    if mismatch_reasons:
        return FieldApplicability(
            "blocked",
            reasons=tuple(
                f"measured wind-history gate failed: {reason}"
                for reason in mismatch_reasons
            ),
        ), assessment
    return FieldApplicability(
        "conditional",
        warnings=(
            "measured wind history passed the declared steady-wind gate; "
            "the local fixed-wind plane uses a nominal wind consistent with the measured record",
        ),
        uncertainty_complete=True,
    ), assessment


def _receptor_for_sensor(
    request: FieldSemiFVRequest,
    frame: WindFrame,
    transport: SemiFVConfig,
    sensor: SensorModel,
    label: str,
) -> tuple[SemiFVReceptor | None, tuple[str, ...]]:
    downwind, crosswind = frame.local(sensor.position_m[0], sensor.position_m[1])
    tolerance = request.lateral_sensor_tolerance_m
    if tolerance is None:
        tolerance = 0.5 * transport.length_m / transport.nx
    if abs(crosswind) > tolerance:
        return None, (
            f"sensor {label!r} is outside the local wind plane; no off-centre concentration was inferred",
        )
    if downwind < 0.0 or downwind > transport.length_m or not 0.0 <= sensor.position_m[2] <= transport.height_m:
        return None, (
            f"sensor {label!r} lies outside the local semi-FV domain; no sensor trace was produced",
        )
    return SemiFVReceptor(label, downwind, sensor.position_m[2]), ()


def _declared_sensor_deployments(
    request: FieldSemiFVRequest,
) -> tuple[FieldSensorDeployment, ...]:
    """Keep the legacy scenario sensor while allowing named additional detectors."""
    primary = () if request.scenario.sensor is None else (
        FieldSensorDeployment("field_sensor", request.scenario.sensor),
    )
    return primary + request.sensor_deployments


def _transport_config_for_screening(
    request: FieldSemiFVRequest,
    vapour_mass_flow_kg_s: float,
    local_obstacles: RectangularObstacle2D | tuple[RectangularObstacle2D, ...] | None,
    distributed_sources: tuple[DistributedScalarSource, ...] = (),
) -> SemiFVConfig:
    source = request.scenario.source
    source_location = _transport_source_location(request)
    source_duration = (
        source.duration_s if request.scenario.temporal_mode == "transient"
        else request.transport.duration_s
    )
    assert source_duration is not None
    duration = source_duration + request.post_release_duration_s
    source_schedule = request.direct_vapour_schedule
    source_rate = vapour_mass_flow_kg_s
    if source_schedule is not None:
        source_rate = 0.0
    elif request.post_release_duration_s > 0.0:
        # A nominal transient source is converted to a transparent finite
        # rate record so post-release transport cannot keep feeding it.
        source_schedule = SourceRateSchedule(
            (0.0, source_duration, duration),
            (vapour_mass_flow_kg_s, 0.0, 0.0),
            source_id="field:finite-nominal-direct-vapour",
        )
        source_rate = 0.0
    diffusivity = (
        request.transport.diffusivity_m2_s
        if request.stability_mixing_closure is None
        else request.stability_mixing_closure.selected_diffusivity_m2_s(
            request.scenario.weather.stability
        )
    )
    return replace(
        request.transport,
        duration_s=duration,
        wind_speed_m_s=request.scenario.weather.speed_m_s.nominal,
        source_rate_kg_s=source_rate,
        source_schedule=source_schedule,
        source_height_m=source_location[2],
        source_sigma_m=(
            request.transport.source_sigma_m
            if request.jet_scalar_handoff is None
            else request.jet_scalar_handoff.vertical_sigma_m
        ),
        diffusivity_m2_s=diffusivity,
        obstacle=local_obstacles if local_obstacles is not None else request.transport.obstacle,
        distributed_sources=distributed_sources,
    )


def _transport_source_location(
    request: FieldSemiFVRequest,
) -> tuple[float, float, float]:
    """Return the original or conservation-screened jet source location."""
    handoff = request.jet_scalar_handoff
    if handoff is None:
        return tuple(float(value) for value in request.scenario.source.location_m)
    return jet_scalar_handoff_location(
        request.scenario.source, request.scenario.weather, handoff,
    )


def _local_distributed_sources(
    request: FieldSemiFVRequest,
    frame: WindFrame,
) -> tuple[DistributedScalarSource, ...]:
    """Map declared global vapour sources into the one supported wind plane."""
    tolerance = request.lateral_source_tolerance_m
    if tolerance is None:
        tolerance = 0.5 * request.transport.length_m / request.transport.nx
    converted = []
    for source in request.distributed_vapour_sources:
        downwind, crosswind = frame.local(source.position_m[0], source.position_m[1])
        if abs(crosswind) > tolerance:
            raise ValueError(
                f"distributed vapour source {source.label!r} is outside the local wind plane"
            )
        if downwind < 0.0 or downwind > request.transport.length_m:
            raise ValueError(
                f"distributed vapour source {source.label!r} lies outside the local downwind domain"
            )
        if source.position_m[2] > request.transport.height_m:
            raise ValueError(
                f"distributed vapour source {source.label!r} lies above the local domain"
            )
        converted.append(DistributedScalarSource(
            source.label, downwind, source.position_m[2], source.vertical_sigma_m,
            source.schedule,
        ))
    return tuple(converted)


def run_field_semi_fv_screening(request: FieldSemiFVRequest) -> FieldSemiFVScreeningResult:
    """Run the direct-flash-vapour reduced-order field screening path.

    The calculation preserves the existing free-field models: it is an opt-in
    local scalar transport screen, not a replacement for the validated
    DEGADIS-compatible path.  It never adds post-flash liquid, rainout, or
    pool vapour to the direct nozzle-vapour source without an explicit phase
    routing calculation.
    """
    scenario = request.scenario
    if scenario.model_family != "degali":
        applicability = FieldApplicability(
            "blocked",
            reasons=("local semi-FV screening is a DEGALI path, not a SLABx solver",),
        )
        return FieldSemiFVScreeningResult(
            request, applicability, None, None, None, None, None
        )

    source = scenario.source
    try:
        transport_source_location = _transport_source_location(request)
    except ValueError as error:
        applicability = FieldApplicability("blocked", reasons=(str(error),))
        return FieldSemiFVScreeningResult(
            request, applicability, None, None, None, None, None,
        )
    wind_to_earth_rad = scenario.weather.wind_to_math_radians()
    wind_to_local_rad = (
        wind_to_earth_rad if request.coordinate_reference is None
        else request.coordinate_reference.earth_to_local_math_radians(wind_to_earth_rad)
    )
    frame = WindFrame(
        origin_x_m=transport_source_location[0],
        origin_y_m=transport_source_location[1],
        direction_rad=wind_to_local_rad,
    )
    try:
        local_distributed_sources = _local_distributed_sources(request, frame)
    except ValueError as error:
        applicability = FieldApplicability("blocked", reasons=(str(error),))
        return FieldSemiFVScreeningResult(
            request, applicability, None, frame, None, None, None,
        )
    projections = tuple(
        project_cuboid_to_wind_plane(cuboid, frame)
        for cuboid in _declared_global_obstacles(request)
    )
    # Preserve the original singular field for callers that submitted one
    # obstacle, while exposing every projection in the plural audit record.
    projection = projections[0] if len(projections) == 1 else None
    local_obstacles = _local_obstacles_from_projections(projections)
    wind_history_applicability, wind_history_assessment = _wind_history_applicability(
        request
    )
    field_gate = assess_field_applicability(
        scenario,
        obstacle_present=local_obstacles is not None,
        obstacle_supported=not any(item.status == "blocked" for item in projections),
        lh2_validation_available=request.lh2_validation_available,
        validation_evidence=request.validation_evidence,
    )
    applicability = _merge_applicability(
        field_gate,
        _projection_applicability(projections),
        wind_history_applicability,
    )
    if applicability.status == "blocked":
        return FieldSemiFVScreeningResult(
            request, applicability, None, frame, projection, None, None,
            wind_history_assessment, obstacle_projections=projections,
        )
    source_preparation = prepare_field_lh2_flash(
        scenario,
        ambient_temperature_k=request.ambient_temperature_k,
        ambient_pressure_pa=request.ambient_pressure_pa,
        lh2_validation_available=request.lh2_validation_available,
        validation_evidence=request.validation_evidence,
        property_table=request.property_table,
    )
    applicability = _merge_applicability(
        applicability,
        source_preparation.applicability,
    )
    if applicability.status == "blocked" or source_preparation.flash_result is None:
        return FieldSemiFVScreeningResult(
            request, applicability, source_preparation, frame, projection, None, None,
            wind_history_assessment, obstacle_projections=projections,
        )

    flash = source_preparation.flash_result.flash
    schedule = request.direct_vapour_schedule
    has_distributed_source = bool(local_distributed_sources)
    if schedule is None and flash.vapour_mass_flow <= 0.0 and not has_distributed_source:
        blocked = _merge_applicability(
            applicability,
            FieldApplicability(
                "blocked", reasons=("flash produces no direct atmospheric vapour source",)
            ),
        )
        return FieldSemiFVScreeningResult(
            request, blocked, source_preparation, frame, projection, None, None,
            wind_history_assessment, obstacle_projections=projections,
        )
    handoff = request.jet_scalar_handoff
    if handoff is not None and not math.isclose(
        flash.vapour_mass_flow,
        handoff.hydrogen_mass_flow_kg_s,
        rel_tol=handoff.maximum_flash_rate_relative_residual,
        abs_tol=1.0e-12,
    ):
        relative_difference = abs(
            flash.vapour_mass_flow - handoff.hydrogen_mass_flow_kg_s
        ) / max(handoff.hydrogen_mass_flow_kg_s, 1.0e-30)
        blocked = _merge_applicability(
            applicability,
            FieldApplicability(
                "blocked",
                reasons=(
                    "jet scalar handoff H2 rate does not match the freshly prepared "
                    f"direct-flash rate (relative difference {relative_difference:.6g})",
                ),
            ),
        )
        return FieldSemiFVScreeningResult(
            request, blocked, source_preparation, frame, projection, None, None,
            wind_history_assessment, obstacle_projections=projections,
        )
    if schedule is not None and schedule.released_mass_kg <= 0.0 and not has_distributed_source:
        blocked = _merge_applicability(
            applicability,
            FieldApplicability(
                "blocked", reasons=("direct-vapour schedule releases no atmospheric hydrogen mass",)
            ),
        )
        return FieldSemiFVScreeningResult(
            request, blocked, source_preparation, frame, projection, None, None,
            wind_history_assessment, obstacle_projections=projections,
        )

    transport_config = _transport_config_for_screening(
        request, flash.vapour_mass_flow, local_obstacles, local_distributed_sources,
    )
    receptor_entries = tuple(
        (
            deployment,
            *_receptor_for_sensor(
                request, frame, transport_config, deployment.sensor, deployment.label,
            ),
        )
        for deployment in _declared_sensor_deployments(request)
    )
    extra_warnings = [
        warning for _deployment, _receptor, warnings in receptor_entries
        for warning in warnings
    ]
    extra_warnings.append(
        "semi-FV transports direct flash vapour as a scalar; post-flash liquid requires explicit rainout/pool routing"
    )
    if schedule is not None:
        history_kind = dict(request.measured_history_provenance).get("history_kind")
        if history_kind == "pressure_driven_orifice":
            extra_warnings.append(
                "declared post-flash direct-vapour schedule replaces the nominal flash rate; "
                "upstream pressure/temperature history was re-flashed through the explicit pressure-driven throat closure"
            )
        else:
            extra_warnings.append(
                "declared post-flash direct-vapour schedule replaces the nominal flash rate; "
                "upstream pressure/temperature history was not re-flashed"
            )
        extra_warnings.extend(request.direct_vapour_warnings)
    if handoff is not None:
        extra_warnings.extend(handoff.warnings)
    for distributed in request.distributed_vapour_sources:
        extra_warnings.append(
            f"distributed vapour source {distributed.label!r} uses declared "
            f"{distributed.source_kind!r} handoff with evidence {distributed.evidence_id!r}"
        )
        extra_warnings.extend(distributed.warnings)
    for deployment in request.sensor_deployments:
        if any(not value.is_exact for value in deployment.sensor.uncertainty_fields().values()):
            extra_warnings.append(
                f"sensor deployment {deployment.label!r} uses nominal calibration values in this single screen; run a declared sensor-array uncertainty study before an operational decision"
            )
    if request.post_release_duration_s > 0.0:
        extra_warnings.append(
            "semi-FV continues the finite scalar inventory after source shutoff; "
            "this is not a three-dimensional Gaussian-puff or gravity-spreading closure"
        )
    if request.stability_mixing_closure is None and scenario.weather.stability != "neutral":
        extra_warnings.append(
            "stability class is recorded but not converted into a turbulence closure; set diffusivity explicitly"
        )
    elif request.stability_mixing_closure is not None:
        closure = request.stability_mixing_closure
        extra_warnings.append(
            "stability class uses caller-declared scalar mixing closure "
            f"{closure.model_id!r} with evidence {closure.evidence_id!r}; "
            "this is not a turbulence-resolving closure"
        )
    if scenario.surface.heat_transfer_w_m2_k.nominal > 0.0:
        extra_warnings.append(
            "surface heat-transfer input is retained for phase routing but is not used by the direct-vapour screen"
        )
    try:
        transport = solve_semi_fv_obstacle(
            transport_config,
            receptors=tuple(
                receptor for _deployment, receptor, _warnings in receptor_entries
                if receptor is not None
            ),
        )
    except ValueError as error:
        blocked = _merge_applicability(
            applicability,
            FieldApplicability("blocked", reasons=(str(error),)),
            warnings=tuple(extra_warnings),
        )
        return FieldSemiFVScreeningResult(
            request, blocked, source_preparation, frame, projection, None, None,
            wind_history_assessment, obstacle_projections=projections,
        )

    conservation_applicability = _transport_conservation_applicability(transport)
    traces_by_label = {trace.receptor.label: trace for trace in transport.receptor_traces}
    sensor_results = []
    sensor_trace = None
    for deployment, receptor, warnings in receptor_entries:
        true_trace = None if receptor is None else traces_by_label.get(deployment.label)
        trace = None if true_trace is None else apply_sensor_model(
            true_trace.time_s,
            true_trace.concentration_kg_m3,
            deployment.sensor,
            ambient_air_density_kg_m3=request.ambient_air_density_kg_m3,
        )
        result_warnings = warnings
        if receptor is not None and trace is None:
            result_warnings = result_warnings + (
                f"sensor {deployment.label!r} receptor trace was unavailable after transport",
            )
        sensor_results.append(FieldSensorDeploymentResult(
            deployment.label, deployment.sensor, trace, result_warnings,
        ))
        if deployment.label == "field_sensor":
            sensor_trace = trace
    final_applicability = _merge_applicability(
        applicability,
        FieldApplicability(
            transport.diagnostics.applicability,
            warnings=transport.diagnostics.warnings,
            uncertainty_complete=True,
        ),
        conservation_applicability,
        warnings=tuple(extra_warnings),
    )
    return FieldSemiFVScreeningResult(
        request, final_applicability, source_preparation, frame, projection,
        transport, sensor_trace, wind_history_assessment,
        obstacle_projections=projections,
        sensor_results=tuple(sensor_results),
    )


def _sensor_at_calibration_corner(
    sensor: SensorModel,
    label: str,
    values: dict[str, float],
) -> SensorModel:
    """Fix a detector's bounded calibration fields at one named corner."""
    fields = sensor.uncertainty_fields()

    def exact(name: str) -> BoundedValue:
        original = fields[name]
        key = f"{label}.{name}"
        return BoundedValue(
            values[key], unit=original.unit, source=f"{original.source}; sensor-array-corner",
        )

    return replace(
        sensor,
        response_time_s=exact("sensor_response_time_s"),
        gain=exact("sensor_gain"),
        bias_mole_fraction=exact("sensor_bias_mole_fraction"),
    )


def run_field_sensor_array_uncertainty_envelope(
    request: FieldSemiFVRequest,
    *,
    max_cases: int = 128,
) -> FieldSensorArrayUncertaintyEnvelope:
    """Vary every detector calibration while preserving one true transport field.

    Sensor response, gain, bias and averaging are observational operators;
    they do not alter the conservative scalar transport. This function first
    runs the field screen once, then applies every declared lower/upper
    calibration corner to its exact receptor concentration histories. It
    deliberately rejects withheld detectors because an off-plane sensor has
    no true local trace on which a calibration envelope could be based.
    """
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(max_cases, int) or max_cases < 1:
        raise ValueError("max_cases must be a positive integer")
    deployments = _declared_sensor_deployments(request)
    if not deployments:
        raise ValueError("sensor-array uncertainty envelope requires at least one declared sensor")
    screening = run_field_semi_fv_screening(request)
    if not screening.completed or screening.transport is None:
        return FieldSensorArrayUncertaintyEnvelope(
            screening, (),
            ("base field screen did not complete; no sensor calibration envelope was produced",),
        )
    nominal_results = {item.label: item for item in screening.sensor_results}
    withheld = [
        deployment.label for deployment in deployments
        if deployment.label not in nominal_results or nominal_results[deployment.label].trace is None
    ]
    if withheld:
        raise ValueError(
            "sensor-array uncertainty envelope cannot evaluate withheld detectors: "
            + ", ".join(withheld)
        )
    fields = {
        f"{deployment.label}.{name}": value
        for deployment in deployments
        for name, value in deployment.sensor.uncertainty_fields().items()
    }
    choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in fields.items()
    }
    count = math.prod(len(values) for values in choices.values())
    if count > max_cases:
        raise ValueError(f"{count} sensor-array calibration corners exceed max_cases={max_cases}")
    traces_by_label = {trace.receptor.label: trace for trace in screening.transport.receptor_traces}
    names = tuple(choices)
    cases = []
    for selected in product(*(choices[name] for name in names)):
        values = {name: float(value) for name, value in zip(names, selected)}
        sensor_results = []
        for deployment in deployments:
            sensor = _sensor_at_calibration_corner(deployment.sensor, deployment.label, values)
            true_trace = traces_by_label[deployment.label]
            trace = apply_sensor_model(
                true_trace.time_s,
                true_trace.concentration_kg_m3,
                sensor,
                ambient_air_density_kg_m3=request.ambient_air_density_kg_m3,
            )
            sensor_results.append(FieldSensorDeploymentResult(
                deployment.label, sensor, trace,
            ))
        cases.append(FieldSensorArrayUncertaintyCase(
            tuple((name, values[name]) for name in names), tuple(sensor_results),
        ))
    return FieldSensorArrayUncertaintyEnvelope(
        screening=screening,
        cases=tuple(cases),
        warnings=(
            "sensor-array calibration corners reuse one completed conservative scalar transport field; only response, gain and bias vary",
            "this is deterministic lower/upper sensitivity, not a probability interval or a sensor-field calibration validation",
        ),
    )


def run_field_semi_fv_refinement_study(
    request: FieldSemiFVRequest,
    *,
    refinement_factors: tuple[int, ...] = (1, 2),
    relative_tolerance: float = 0.05,
    absolute_floor_kg_m3: float = 1.0e-12,
    max_cell_steps: int = 20_000_000,
) -> FieldSemiFVRefinementResult:
    """Audit a completed field screen at its declared in-plane sensor.

    This is deliberately opt-in because it recomputes the local transport at
    multiple grid/time resolutions. It reports numerical resolution only; it
    does not remove the field screen's source, geometry or physical-validation
    limitations.
    """
    screening = run_field_semi_fv_screening(request)
    if not screening.completed:
        return FieldSemiFVRefinementResult(
            screening, None, screening.applicability,
        )
    if request.scenario.sensor is None:
        applicability = _merge_applicability(
            screening.applicability,
            FieldApplicability(
                "blocked", reasons=("field refinement study requires a declared sensor",),
            ),
        )
        return FieldSemiFVRefinementResult(screening, None, applicability)
    if (
        screening.source_preparation is None
        or screening.source_preparation.flash_result is None
        or screening.wind_frame is None
    ):
        applicability = _merge_applicability(
            screening.applicability,
            FieldApplicability("blocked", reasons=("screening source preparation is unavailable",)),
        )
        return FieldSemiFVRefinementResult(screening, None, applicability)
    flash = screening.source_preparation.flash_result.flash
    projections = screening.obstacle_projections or (
        () if screening.obstacle_projection is None else (screening.obstacle_projection,)
    )
    transport_config = _transport_config_for_screening(
        request,
        flash.vapour_mass_flow,
        _local_obstacles_from_projections(projections),
        _local_distributed_sources(request, screening.wind_frame),
    )
    receptor, receptor_warnings = _receptor_for_sensor(
        request, screening.wind_frame, transport_config,
        request.scenario.sensor, "field_sensor",
    )
    if receptor is None:
        applicability = _merge_applicability(
            screening.applicability,
            FieldApplicability(
                "blocked",
                reasons=("field sensor cannot be represented in the local refinement plane",),
            ),
            warnings=receptor_warnings,
        )
        return FieldSemiFVRefinementResult(screening, None, applicability)
    try:
        study = run_transport_refinement_study(
            transport_config,
            receptors=(receptor,),
            refinement_factors=refinement_factors,
            relative_tolerance=relative_tolerance,
            absolute_floor_kg_m3=absolute_floor_kg_m3,
            max_cell_steps=max_cell_steps,
        )
    except ValueError as error:
        applicability = _merge_applicability(
            screening.applicability,
            FieldApplicability("blocked", reasons=(str(error),)),
            warnings=receptor_warnings,
        )
        return FieldSemiFVRefinementResult(screening, None, applicability)
    numerical = FieldApplicability(
        "accepted" if study.converged else "conditional",
        warnings=study.warnings + receptor_warnings,
        uncertainty_complete=study.converged,
    )
    return FieldSemiFVRefinementResult(
        screening, study, _merge_applicability(screening.applicability, numerical),
    )


def _distributed_source_corner_cases(
    sources: tuple[FieldDistributedVapourSource, ...],
    *,
    max_cases: int,
) -> tuple[
    tuple[tuple[FieldDistributedVapourSource, ...], tuple[tuple[str, float | str], ...]],
    ...,
]:
    """Enumerate explicit geometry/width/rate corners for supplemental sources."""
    if not sources:
        return (((), ()),)
    fields: dict[str, BoundedValue] = {}
    for source in sources:
        fields.update(source.uncertainty_fields())
    choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in fields.items()
    }
    schedule_choices = tuple(source.schedule_corners() for source in sources)
    count = math.prod(len(values) for values in choices.values()) * math.prod(
        len(values) for values in schedule_choices
    )
    if count > max_cases:
        raise ValueError(
            f"{count} distributed-source uncertainty corners exceed max_cases={max_cases}"
        )
    names = tuple(choices)
    cases = []
    geometry_products = product(*(choices[name] for name in names))
    schedule_products = tuple(product(*schedule_choices))
    for selected in geometry_products:
        values = {name: float(value) for name, value in zip(names, selected)}
        geometry_selection = tuple((name, values[name]) for name in names)
        for selected_schedules in schedule_products:
            schedule_by_label = {
                source.label: (label, schedule)
                for source, (label, schedule) in zip(sources, selected_schedules)
            }
            fixed_sources = tuple(
                (
                    source.at_corner({
                        name: values[name] for name in source.uncertainty_fields()
                    })
                    if source.uncertainty_fields() else source
                ).at_schedule_corner(schedule_by_label[source.label][0])
                for source in sources
            )
            schedule_selection = tuple(
                (
                    f"distributed_source.{source.label}.schedule",
                    schedule_by_label[source.label][0],
                )
                for source in sources
                if source.has_schedule_uncertainty
            )
            cases.append((fixed_sources, geometry_selection + schedule_selection))
    return tuple(cases)


def _obstacle_geometry_corner_cases(
    request: FieldSemiFVRequest,
    *,
    max_cases: int,
) -> tuple[tuple[tuple[AxisAlignedCuboid | OrientedCuboid, ...], tuple[tuple[str, float], ...]], ...]:
    """Enumerate explicit global-obstacle geometry corners."""
    declared = _declared_global_obstacles(request)
    uncertainties = request.obstacle_geometry_uncertainty
    if not uncertainties:
        return ((declared, ()),)
    fields: dict[str, BoundedValue] = {}
    for uncertainty in uncertainties:
        fields.update(uncertainty.uncertainty_fields())
    choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in fields.items()
    }
    count = math.prod(len(values) for values in choices.values())
    if count > max_cases:
        raise ValueError(
            f"{count} obstacle-geometry uncertainty corners exceed max_cases={max_cases}"
        )
    names = tuple(choices)
    by_label = {item.obstacle.label: item for item in uncertainties}
    cases = []
    for selected in product(*(choices[name] for name in names)):
        values = {name: float(value) for name, value in zip(names, selected)}
        fixed = tuple(
            by_label[obstacle.label].at_corner({
                name: values[name]
                for name in by_label[obstacle.label].uncertainty_fields()
            })
            if obstacle.label in by_label else obstacle
            for obstacle in declared
        )
        cases.append((fixed, tuple((name, values[name]) for name in names)))
    return tuple(cases)


def run_field_semi_fv_envelope(
    request: FieldSemiFVRequest,
    *,
    max_cases: int = 128,
    table_nodes: int = 161,
    atmospheric_source_schedule: FieldAtmosphericSourceSchedule | None = None,
) -> FieldSemiFVEnvelope:
    """Run every declared uncertainty corner without probabilistic relabelling.

    The optional table is built once across the release-temperature interval
    and reused in every corner. Explicit geometry/vertical-width bounds on
    distributed vapour sources and declared obstacle geometry bounds are
    expanded alongside the primary scenario corners; each source shape and
    obstacle mask remains a separate deterministic case. If the table cannot
    be built, the reference direct path remains available and a warning
    records why no table was used.
    """
    if atmospheric_source_schedule is not None and not isinstance(
        atmospheric_source_schedule, FieldAtmosphericSourceSchedule
    ):
        raise TypeError(
            "atmospheric_source_schedule must be a FieldAtmosphericSourceSchedule or None"
        )
    if (
        atmospheric_source_schedule is not None
        and atmospheric_source_schedule.evidence.source_kind
        not in {
            "declared_atmospheric_vapour",
            "post_flash_atmospheric_vapour",
        }
    ):
        raise ValueError(
            "atmospheric_source_schedule can replace only a declared or "
            "post_flash atmospheric-vapour source boundary; pool_vapour and "
            "droplet_evaporation require an explicit "
            "distributed/phase-routing handoff with source location and scalar width"
        )
    if request.direct_vapour_schedule is not None and atmospheric_source_schedule is not None:
        raise ValueError(
            "request direct_vapour_schedule and atmospheric_source_schedule are mutually exclusive"
        )
    if request.direct_vapour_schedule is not None:
        raise ValueError(
            "a direct-vapour schedule needs a time-aligned source uncertainty envelope; "
            "scalar FieldScenario corners would not vary its history"
        )
    source_schedule_corners: tuple[tuple[str, SourceRateSchedule], ...] = ()
    envelope_request = request
    if atmospheric_source_schedule is not None:
        if request.scenario.temporal_mode != "transient" or request.scenario.source.duration_s is None:
            raise ValueError("atmospheric_source_schedule requires a transient scenario duration")
        if not math.isclose(
            atmospheric_source_schedule.schedule.duration_s,
            request.scenario.source.duration_s,
            rel_tol=1.0e-9, abs_tol=1.0e-12,
        ):
            raise ValueError("atmospheric source schedule duration must match source.duration_s")
        if (
            request.scenario.source.duration_uncertainty is not None
            and not request.scenario.source.duration_uncertainty.is_exact
        ):
            raise ValueError(
                "atmospheric_source_schedule cannot represent non-exact release-duration "
                "uncertainty; provide time-aligned source schedule corners"
            )
        all_fields = request.scenario.uncertainty_fields()
        nominal_values = {name: value.nominal for name, value in all_fields.items()}
        nominal_source_scenario = request.scenario.at_corner(nominal_values)
        # A supplied atmospheric history replaces the source-rate and upstream
        # thermodynamic boundary.  Source geometry is still a field input: do
        # not collapse declared location/direction bounds merely because the
        # atmospheric rate schedule is fixed.
        source_geometry = replace(
            nominal_source_scenario.source,
            location_m=request.scenario.source.location_m,
            location_uncertainty_m=request.scenario.source.location_uncertainty_m,
            direction_m=request.scenario.source.direction_m,
            direction_uncertainty_m=request.scenario.source.direction_uncertainty_m,
        )
        envelope_request = replace(
            request,
            scenario=replace(
                nominal_source_scenario,
                source=source_geometry,
                weather=request.scenario.weather,
                surface=request.scenario.surface,
                sensor=request.scenario.sensor,
            ),
        )
        source_schedule_corners = atmospheric_source_schedule.corner_schedules()
    if any(
        not value.is_exact
        for deployment in request.sensor_deployments
        for value in deployment.sensor.uncertainty_fields().values()
    ):
        raise ValueError(
            "additional sensor deployment uncertainty is not represented by FieldScenario corners; "
            "run a declared sensor-array uncertainty study instead"
        )
    # A distributed atmospheric schedule is a separate source branch, but it
    # still has to fit the actual local transport window of every primary
    # release-duration corner.  Check this before constructing corner
    # requests: otherwise ``dataclasses.replace`` would fail halfway through
    # an envelope with a generic duration error and leave the time-alignment
    # requirement implicit.
    duration_uncertainty = envelope_request.scenario.source.duration_uncertainty
    if duration_uncertainty is None or duration_uncertainty.is_exact:
        duration_values = (
            envelope_request.scenario.source.duration_s
            if envelope_request.scenario.source.duration_s is not None
            else envelope_request.transport.duration_s,
        )
    else:
        duration_values = duration_uncertainty.corners()
    for source_duration in duration_values:
        total_duration = float(source_duration) + envelope_request.post_release_duration_s
        for distributed in envelope_request.distributed_vapour_sources:
            for schedule_label, schedule in distributed.schedule_corners():
                if schedule.duration_s > total_duration + 1.0e-12:
                    raise ValueError(
                        f"distributed vapour source {distributed.label!r} {schedule_label} "
                        f"schedule duration {schedule.duration_s:g} s exceeds the primary "
                        f"release-duration uncertainty corner {float(source_duration):g} s; "
                        "provide a time-aligned distributed-source envelope or explicit "
                        "post-release continuation"
                    )
    corners = envelope_request.scenario.corner_cases(max_cases=max_cases)
    ambient_fields = envelope_request.ambient_uncertainty_fields()
    ambient_choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in ambient_fields.items()
    }
    ambient_names = tuple(ambient_choices)
    ambient_corner_cases = tuple(
        (
            {name: float(value) for name, value in zip(ambient_names, selected)},
            tuple((name, float(value)) for name, value in zip(ambient_names, selected)),
        )
        for selected in product(*(ambient_choices[name] for name in ambient_names))
    ) if ambient_names else (({}, ()),)
    stability_alternatives = envelope_request.stability_alternatives
    stability_classes = (
        (request.scenario.weather.stability,)
        if stability_alternatives is None
        else stability_alternatives.classes_for(request.scenario.weather.stability)
    )
    total_cases = len(corners) * len(ambient_corner_cases) * len(stability_classes)
    distributed_source_cases = _distributed_source_corner_cases(
        envelope_request.distributed_vapour_sources, max_cases=max_cases,
    )
    obstacle_geometry_cases = _obstacle_geometry_corner_cases(
        envelope_request, max_cases=max_cases,
    )
    if total_cases * len(distributed_source_cases) * len(obstacle_geometry_cases) > max_cases:
        raise ValueError(
            f"{total_cases * len(distributed_source_cases) * len(obstacle_geometry_cases)} "
            f"field uncertainty corners including ambient, distributed-source and obstacle "
            f"geometry plus stability alternatives "
            f"exceed max_cases={max_cases}"
        )
    table = envelope_request.property_table
    warnings: list[str] = []
    if table is None:
        try:
            table_kwargs = {
                "ambient_pressure_pa": envelope_request.ambient_pressure_pa,
                "nodes": table_nodes,
            }
            ambient_pressure_bound = envelope_request.ambient_pressure_uncertainty_pa
            if ambient_pressure_bound is not None and not ambient_pressure_bound.is_exact:
                table_kwargs["ambient_pressure_bounds_pa"] = (
                    ambient_pressure_bound.lower,
                    ambient_pressure_bound.upper,
                )
            table = build_lh2_saturation_table_for_release(
                envelope_request.scenario.source,
                **table_kwargs,
            )
        except (ValueError, RuntimeError) as error:
            warnings.append(
                "LH2 saturation table was not constructed; direct property path used: "
                f"{error}"
            )
    cases = []
    schedule_multiplier = max(1, len(source_schedule_corners))
    if (
        len(corners) * len(ambient_corner_cases) * len(stability_classes)
        * len(distributed_source_cases) * len(obstacle_geometry_cases)
        * schedule_multiplier > max_cases
    ):
        raise ValueError(
            f"{len(corners) * len(ambient_corner_cases) * len(stability_classes) * len(distributed_source_cases) * len(obstacle_geometry_cases) * schedule_multiplier} "
            f"field uncertainty corners including ambient, distributed-source, obstacle-geometry and source schedule corners "
            f"exceed max_cases={max_cases}"
        )
    for values in corners:
        for ambient_values, ambient_selection in ambient_corner_cases:
            for stability in stability_classes:
                scenario = envelope_request.scenario.at_corner(values)
                if stability_alternatives is not None:
                    scenario = replace(
                        scenario, weather=replace(scenario.weather, stability=stability),
                    )
                for distributed_sources, distributed_selection in distributed_source_cases:
                    for obstacle_geometry, obstacle_selection in obstacle_geometry_cases:
                        for schedule_label, schedule in (source_schedule_corners or (("nominal", None),)):
                            corner_request = replace(
                                envelope_request,
                                scenario=scenario,
                                property_table=table,
                                direct_vapour_schedule=schedule,
                                distributed_vapour_sources=distributed_sources,
                                ambient_temperature_k=ambient_values.get(
                                    "ambient_temperature_k", envelope_request.ambient_temperature_k,
                                ),
                                ambient_pressure_pa=ambient_values.get(
                                    "ambient_pressure_pa", envelope_request.ambient_pressure_pa,
                                ),
                                ambient_air_density_kg_m3=ambient_values.get(
                                    "ambient_air_density_kg_m3", envelope_request.ambient_air_density_kg_m3,
                                ),
                                ambient_temperature_uncertainty_k=None,
                                ambient_pressure_uncertainty_pa=None,
                                ambient_air_density_uncertainty_kg_m3=None,
                                obstacle=(
                                    obstacle_geometry[0]
                                    if envelope_request.obstacle is not None
                                    else None
                                ),
                                obstacles=(
                                    () if envelope_request.obstacle is not None
                                    else obstacle_geometry
                                ),
                                obstacle_geometry_uncertainty=(),
                            )
                            selection: tuple[tuple[str, float | str], ...] = tuple(
                                (name, float(value)) for name, value in values.items()
                            )
                            selection += ambient_selection
                            if stability_alternatives is not None:
                                selection += (("weather_stability", stability),)
                            selection += distributed_selection + obstacle_selection
                            if atmospheric_source_schedule is not None:
                                selection += (("source_schedule", schedule_label),)
                            cases.append(FieldSemiFVEnvelopeCase(
                                selection, run_field_semi_fv_screening(corner_request),
                            ))
    return FieldSemiFVEnvelope(
        request, tuple(cases), property_table_used=table is not None,
        warnings=tuple(warnings), atmospheric_source_schedule=atmospheric_source_schedule,
    )


__all__ = [
    "FieldSensorDeployment", "FieldSensorDeploymentResult",
    "FieldSemiFVRequest", "FieldSemiFVScreeningResult", "FieldSemiFVEnvelopeCase",
    "FieldSemiFVEnvelope", "FieldSensorArrayUncertaintyCase",
    "FieldSensorArrayUncertaintyEnvelope", "FieldSemiFVRefinementResult",
    "run_field_semi_fv_screening", "run_field_semi_fv_envelope",
    "run_field_semi_fv_refinement_study", "run_field_sensor_array_uncertainty_envelope",
]
