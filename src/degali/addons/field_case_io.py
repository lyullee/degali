"""Strict JSON contracts for direct-flash and phase-routed field screening.

The base schema covers a declared direct-flash field case, explicit bounded
obstacle geometry, already-atmospheric distributed source handoffs, and its
separately quality-gated historian attachment.  A source may either declare
its mass-flow bound directly or explicitly opt in to the pressure-driven LH2
throat adapter; omitting both remains invalid.  The dedicated phase-routing
schema adds only the existing conservative rainout/dynamic-pool transport
adapter.  It does not discover plant tags, infer phase from a level trace,
convert normal volume flow, or manufacture atmospheric launch or wake
coefficients.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping
from dataclasses import dataclass

from .field_contracts import (
    BoundedValue,
    CircularBoundedValue,
    FieldCoordinateReference,
    FieldValidationEvidence,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
)
from .field_json import strict_json_loads
from .droplet_rainout import DropletClass
from .field_meteorology import FieldStabilityAlternatives, StabilityScalarMixingClosure
from .field_decision import FieldConditionalReviewAuthorization
from .field_phase_routing import FieldPhaseRoutingConfig, FieldPhaseRoutingUncertainty
from .field_pool_launch import PoolVapourLaunchBoundary
from .field_distributed_source import FieldDistributedVapourSource
from .field_geometry import FieldObstacleGeometryUncertainty
from .field_historian_io import (
    HistorianCsvChannel,
    ImportedMeasuredReleaseHistory,
    ImportedPressureDrivenMeasuredHistory,
    MeasuredHistoryCsvMap,
    PressureDrivenHistoryCsvMap,
    direct_vapour_schedule_from_imported_history,
    direct_vapour_schedule_from_imported_pressure_driven_history,
    read_measured_history_csv,
    read_pressure_driven_history_csv,
)
from .field_history import (
    FieldMeasuredFlashSchedule,
    MeasuredHistoryQualityCriteria,
    request_with_measured_flash_schedule,
    request_with_pressure_driven_history_schedule,
)
from .field_lh2 import release_with_pressure_driven_lh2_mass_flow
from .field_workflow import FieldSensorDeployment, FieldSemiFVRequest
from .semi_fv_obstacle import SemiFVConfig, SourceRateSchedule
from .site_geometry import AxisAlignedCuboid, OrientedCuboid


FIELD_SCREENING_INPUT_SCHEMA = "degali.field-screening-input.v1"
FIELD_PHASE_ROUTING_SCREENING_INPUT_SCHEMA = (
    "degali.field-phase-routing-screening-input.v1"
)


@dataclass(frozen=True)
class FieldPhaseRoutingTransportInput:
    """Strict user inputs for the conservative phase/pool transport adapter."""

    config: FieldPhaseRoutingConfig
    pool_launch: PoolVapourLaunchBoundary
    pool_vertical_sigma_m: float
    droplet_population: tuple[DropletClass, ...]
    droplet_population_evidence_id: str
    phase_routing_evidence_id: str
    table_nodes: int = 161
    phase_uncertainty: FieldPhaseRoutingUncertainty | None = None
    pool_vertical_sigma_uncertainty: BoundedValue | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.config, FieldPhaseRoutingConfig):
            raise TypeError("phase-routing config must be FieldPhaseRoutingConfig")
        if self.phase_uncertainty is not None and not isinstance(
            self.phase_uncertainty, FieldPhaseRoutingUncertainty
        ):
            raise TypeError(
                "phase_uncertainty must be FieldPhaseRoutingUncertainty or None"
            )
        if self.phase_uncertainty is not None:
            # Keep direct Python construction subject to the same nominal
            # anchor/simplex checks as the strict JSON parser.
            self.phase_uncertainty.config_corners(self.config)
        if self.pool_vertical_sigma_uncertainty is not None and not isinstance(
            self.pool_vertical_sigma_uncertainty, BoundedValue
        ):
            raise TypeError(
                "pool_vertical_sigma_uncertainty must be BoundedValue or None"
            )
        if self.pool_vertical_sigma_uncertainty is not None:
            bound = self.pool_vertical_sigma_uncertainty
            if (
                not isinstance(bound.source, str)
                or not bound.source.strip()
                or bound.source.lower() == "unspecified"
            ):
                raise ValueError(
                    "pool_vertical_sigma_uncertainty requires an explicit source"
                )
            if bound.lower <= 0.0 or bound.upper <= 0.0:
                raise ValueError(
                    "pool_vertical_sigma_uncertainty bounds must be positive"
                )
            if not math.isclose(
                float(self.pool_vertical_sigma_m), bound.nominal,
                rel_tol=1.0e-12, abs_tol=1.0e-15,
            ):
                raise ValueError(
                    "pool_vertical_sigma_uncertainty.nominal must match "
                    "pool_vertical_sigma_m"
                )
        if not isinstance(self.pool_launch, PoolVapourLaunchBoundary):
            raise TypeError("pool_launch must be PoolVapourLaunchBoundary")
        if (
            not math.isfinite(float(self.pool_vertical_sigma_m))
            or self.pool_vertical_sigma_m <= 0.0
        ):
            raise ValueError("pool_vertical_sigma_m must be positive and finite")
        if not self.droplet_population or not all(
            isinstance(item, DropletClass) for item in self.droplet_population
        ):
            raise TypeError("droplet_population must contain declared DropletClass values")
        configured_population = tuple(
            self.config.phase_model_options.get("droplet_classes", ())
        )
        if configured_population != self.droplet_population:
            raise ValueError(
                "phase-routing config droplet classes must match the declared population"
            )
        if (
            not isinstance(self.droplet_population_evidence_id, str)
            or not self.droplet_population_evidence_id.strip()
            or self.droplet_population_evidence_id == "unspecified"
        ):
            raise ValueError("droplet_population_evidence_id must be explicitly declared")
        if (
            not isinstance(self.phase_routing_evidence_id, str)
            or not self.phase_routing_evidence_id.strip()
            or self.phase_routing_evidence_id == "unspecified"
        ):
            raise ValueError("phase_routing_evidence_id must be explicitly declared")
        if isinstance(self.table_nodes, bool) or not isinstance(self.table_nodes, int):
            raise TypeError("table_nodes must be an integer")
        if self.table_nodes < 4:
            raise ValueError("table_nodes must be at least 4")


@dataclass(frozen=True)
class FieldScreeningCase:
    """Parsed case plus optional promoted measured-source provenance."""

    request: FieldSemiFVRequest
    source_request: FieldSemiFVRequest | None = None
    imported_history: ImportedMeasuredReleaseHistory | ImportedPressureDrivenMeasuredHistory | None = None
    measured_schedule: FieldMeasuredFlashSchedule | None = None
    phase_routing_transport: FieldPhaseRoutingTransportInput | None = None
    conditional_review: FieldConditionalReviewAuthorization | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("field case request must be FieldSemiFVRequest")
        if self.source_request is not None and not isinstance(
            self.source_request, FieldSemiFVRequest
        ):
            raise TypeError("source_request must be FieldSemiFVRequest or None")
        if self.imported_history is None and self.measured_schedule is not None:
            raise ValueError("measured schedule requires an imported history")
        if self.imported_history is not None and not isinstance(
            self.imported_history,
            (ImportedMeasuredReleaseHistory, ImportedPressureDrivenMeasuredHistory),
        ):
            raise TypeError(
                "imported_history must be an imported measured or pressure-driven history or None"
            )
        if self.measured_schedule is not None and not isinstance(
            self.measured_schedule, FieldMeasuredFlashSchedule
        ):
            raise TypeError("measured_schedule must be FieldMeasuredFlashSchedule or None")
        if self.phase_routing_transport is not None and not isinstance(
            self.phase_routing_transport, FieldPhaseRoutingTransportInput
        ):
            raise TypeError(
                "phase_routing_transport must be FieldPhaseRoutingTransportInput or None"
            )
        if self.imported_history is not None and self.phase_routing_transport is not None:
            raise ValueError("measured history and phase-routing transport are mutually exclusive")
        if self.conditional_review is not None and not isinstance(
            self.conditional_review, FieldConditionalReviewAuthorization
        ):
            raise TypeError(
                "conditional_review must be FieldConditionalReviewAuthorization or None"
            )


def _object(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _keys(
    value: Mapping[str, Any],
    name: str,
    *,
    required: set[str],
    optional: set[str] | None = None,
) -> None:
    present = set(value)
    optional = set() if optional is None else optional
    missing = sorted(required - present)
    unknown = sorted(present - required - optional)
    if missing or unknown:
        detail = []
        if missing:
            detail.append("missing=" + ", ".join(missing))
        if unknown:
            detail.append("unknown=" + ", ".join(unknown))
        raise ValueError(f"{name} keys are invalid: " + "; ".join(detail))


def _string(value: object, name: str, *, declared: bool = False) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    result = value.strip()
    if declared and result.lower() == "unspecified":
        raise ValueError(f"{name} must be explicitly declared")
    return result


def _number(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number, not boolean")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _ambient_bound(
    data: Mapping[str, Any],
    *,
    nominal_key: str,
    uncertainty_key: str,
    default: float,
    unit: str,
) -> tuple[float, BoundedValue | None]:
    """Parse one ambient nominal plus an optional explicit bounded interval."""
    bound = None
    if uncertainty_key in data and data[uncertainty_key] is not None:
        bound = _bounded(data[uncertainty_key], uncertainty_key, unit=unit)
    nominal = (
        default
        if nominal_key not in data or data[nominal_key] is None
        else _number(data[nominal_key], nominal_key)
    )
    if bound is not None:
        if nominal_key not in data or data[nominal_key] is None:
            nominal = bound.nominal
        elif not math.isclose(
            nominal, bound.nominal, rel_tol=1.0e-9, abs_tol=1.0e-12,
        ):
            raise ValueError(
                f"{uncertainty_key}.nominal must match {nominal_key}"
            )
    return nominal, bound


def _integer(value: object, name: str, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer of at least {minimum}")
    return value


def _boolean(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be boolean")
    return value


def _point3(value: object, name: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must be a three-value coordinate array")
    return tuple(_number(item, f"{name}[{index}]") for index, item in enumerate(value))


def _source_vector(
    value: object, name: str, *, unit: str,
) -> tuple[tuple[float, float, float], tuple[BoundedValue, BoundedValue, BoundedValue] | None]:
    """Parse a legacy numeric vector or three explicit bounded components."""
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must be a three-value vector array")
    mapping_flags = tuple(isinstance(item, Mapping) for item in value)
    if any(mapping_flags) and not all(mapping_flags):
        raise ValueError(
            f"{name} must use either three numbers or three bounded objects"
        )
    if not any(mapping_flags):
        return _point3(value, name), None
    bounds = tuple(
        _bounded(item, f"{name}[{index}]", unit=unit)
        for index, item in enumerate(value)
    )
    nominal = tuple(float(bound.nominal) for bound in bounds)
    return nominal, bounds  # type: ignore[return-value]


def _source_location(
    value: object,
) -> tuple[tuple[float, float, float], tuple[BoundedValue, BoundedValue, BoundedValue] | None]:
    """Parse a source location in metres."""
    return _source_vector(value, "scenario.source.location_m", unit="m")


def _bounded(value: object, name: str, *, unit: str) -> BoundedValue:
    data = _object(value, name)
    _keys(
        data, name, required={"nominal", "unit", "source"},
        optional={"lower", "upper"},
    )
    if data["unit"] != unit:
        raise ValueError(f"{name}.unit must be {unit!r}")
    return BoundedValue(
        _number(data["nominal"], f"{name}.nominal"),
        None if "lower" not in data else _number(data["lower"], f"{name}.lower"),
        None if "upper" not in data else _number(data["upper"], f"{name}.upper"),
        unit=unit,
        source=_string(data["source"], f"{name}.source", declared=True),
    )


def _circular_direction(value: object, name: str) -> CircularBoundedValue:
    """Parse a declared degree interval, including one which crosses north."""
    data = _object(value, name)
    _keys(
        data, name, required={"nominal", "unit", "source"},
        optional={"lower", "upper"},
    )
    if data["unit"] != "deg":
        raise ValueError(f"{name}.unit must be 'deg'")
    return CircularBoundedValue(
        _number(data["nominal"], f"{name}.nominal"),
        None if "lower" not in data else _number(data["lower"], f"{name}.lower"),
        None if "upper" not in data else _number(data["upper"], f"{name}.upper"),
        unit="deg",
        source=_string(data["source"], f"{name}.source", declared=True),
    )


def _coordinate_reference(value: object) -> FieldCoordinateReference:
    data = _object(value, "coordinate_reference")
    _keys(
        data, "coordinate_reference",
        required={
            "coordinate_system_id", "origin_id", "x_axis_bearing_math_to_deg",
            "vertical_datum_id", "evidence_id",
        },
    )
    return FieldCoordinateReference(
        coordinate_system_id=_string(
            data["coordinate_system_id"], "coordinate_reference.coordinate_system_id", declared=True,
        ),
        origin_id=_string(data["origin_id"], "coordinate_reference.origin_id", declared=True),
        x_axis_bearing_math_to_deg=_number(
            data["x_axis_bearing_math_to_deg"],
            "coordinate_reference.x_axis_bearing_math_to_deg",
        ),
        vertical_datum_id=_string(
            data["vertical_datum_id"], "coordinate_reference.vertical_datum_id", declared=True,
        ),
        evidence_id=_string(data["evidence_id"], "coordinate_reference.evidence_id", declared=True),
    )


def _validation_evidence(value: object) -> FieldValidationEvidence:
    data = _object(value, "validation_evidence")
    _keys(
        data,
        "validation_evidence",
        required={
            "dataset_id", "path", "sha256", "row_count", "source_boundary_id",
            "weather_id", "obstacle_geometry_id", "receptor_geometry_id",
            "temporal_operator_id", "common_clock_id", "scope",
        },
    )
    return FieldValidationEvidence(
        dataset_id=_string(data["dataset_id"], "validation_evidence.dataset_id", declared=True),
        path=_string(data["path"], "validation_evidence.path", declared=True),
        sha256=_string(data["sha256"], "validation_evidence.sha256", declared=True),
        row_count=_integer(data["row_count"], "validation_evidence.row_count", minimum=1),
        source_boundary_id=_string(
            data["source_boundary_id"], "validation_evidence.source_boundary_id", declared=True,
        ),
        weather_id=_string(data["weather_id"], "validation_evidence.weather_id", declared=True),
        obstacle_geometry_id=_string(
            data["obstacle_geometry_id"], "validation_evidence.obstacle_geometry_id", declared=True,
        ),
        receptor_geometry_id=_string(
            data["receptor_geometry_id"], "validation_evidence.receptor_geometry_id", declared=True,
        ),
        temporal_operator_id=_string(
            data["temporal_operator_id"], "validation_evidence.temporal_operator_id", declared=True,
        ),
        common_clock_id=_string(
            data["common_clock_id"], "validation_evidence.common_clock_id", declared=True,
        ),
        scope=_string(data["scope"], "validation_evidence.scope", declared=True),
    )


def _conditional_review(value: object) -> FieldConditionalReviewAuthorization:
    data = _object(value, "conditional_review")
    _keys(
        data,
        "conditional_review",
        required={
            "review_id", "reviewer_id", "reviewer_role", "reviewed_at_utc",
            "evidence_id", "allowed_scope",
        },
    )
    return FieldConditionalReviewAuthorization(
        review_id=_string(data["review_id"], "conditional_review.review_id", declared=True),
        reviewer_id=_string(
            data["reviewer_id"], "conditional_review.reviewer_id", declared=True,
        ),
        reviewer_role=_string(
            data["reviewer_role"], "conditional_review.reviewer_role", declared=True,
        ),
        reviewed_at_utc=_string(
            data["reviewed_at_utc"], "conditional_review.reviewed_at_utc", declared=True,
        ),
        evidence_id=_string(
            data["evidence_id"], "conditional_review.evidence_id", declared=True,
        ),
        allowed_scope=_string(
            data["allowed_scope"], "conditional_review.allowed_scope", declared=True,
        ),
    )


def _source(value: object) -> ReleaseSource:
    data = _object(value, "scenario.source")
    _keys(
        data, "scenario.source",
        required={
            "fluid", "location_m", "direction_m", "upstream_pressure",
            "upstream_temperature", "pressure_reference",
            "opening_area_m2", "discharge_coefficient", "liquid_fraction",
            "liquid_fraction_basis", "flash_model", "duration_s",
        },
        optional={
            "duration_uncertainty", "metadata", "mass_flow_kg_s",
            "pressure_driven_mass_flow",
        },
    )
    has_declared_rate = "mass_flow_kg_s" in data
    has_pressure_derivation = "pressure_driven_mass_flow" in data
    if has_declared_rate == has_pressure_derivation:
        raise ValueError(
            "scenario.source must declare exactly one of mass_flow_kg_s or "
            "pressure_driven_mass_flow"
        )
    metadata = data.get("metadata", {})
    if not isinstance(metadata, Mapping) or any(
        not isinstance(key, str) or not isinstance(item, str)
        for key, item in metadata.items()
    ):
        raise ValueError("scenario.source.metadata must map strings to strings")
    duration = data["duration_s"]
    location_m, location_uncertainty_m = _source_location(data["location_m"])
    direction_m, direction_uncertainty_m = _source_vector(
        data["direction_m"], "scenario.source.direction_m", unit="1",
    )
    duration_uncertainty = None
    if isinstance(duration, Mapping):
        if "duration_uncertainty" in data:
            raise ValueError(
                "scenario.source cannot declare duration_s as a bounded object "
                "and duration_uncertainty together"
            )
        duration_uncertainty = _bounded(
            duration, "scenario.source.duration_s", unit="s",
        )
        duration_value: float | None = duration_uncertainty.nominal
    elif duration is None:
        if "duration_uncertainty" in data:
            raise ValueError(
                "scenario.source.duration_uncertainty requires duration_s"
            )
        duration_value = None
    else:
        duration_value = _number(duration, "scenario.source.duration_s")
        if "duration_uncertainty" in data:
            duration_uncertainty = _bounded(
                data["duration_uncertainty"],
                "scenario.source.duration_uncertainty",
                unit="s",
            )
    pressure_derivation = None
    if has_pressure_derivation:
        pressure_data = _object(
            data["pressure_driven_mass_flow"],
            "scenario.source.pressure_driven_mass_flow",
        )
        _keys(
            pressure_data,
            "scenario.source.pressure_driven_mass_flow",
            required={"ambient_pressure_pa", "source_id"},
            optional={"allow_supercritical_gas"},
        )
        ambient_pressure = _bounded(
            pressure_data["ambient_pressure_pa"],
            "scenario.source.pressure_driven_mass_flow.ambient_pressure_pa",
            unit="Pa",
        )
        pressure_derivation = (
            ambient_pressure,
            _string(
                pressure_data["source_id"],
                "scenario.source.pressure_driven_mass_flow.source_id",
                declared=True,
            ),
            (
                False
                if "allow_supercritical_gas" not in pressure_data
                else _boolean(
                    pressure_data["allow_supercritical_gas"],
                    "scenario.source.pressure_driven_mass_flow.allow_supercritical_gas",
                )
            ),
        )
    source = ReleaseSource(
        fluid=_string(data["fluid"], "scenario.source.fluid"),
        location_m=location_m,
        direction_m=direction_m,
        upstream_pressure=_bounded(data["upstream_pressure"], "scenario.source.upstream_pressure", unit="Pa"),
        upstream_temperature=_bounded(data["upstream_temperature"], "scenario.source.upstream_temperature", unit="K"),
        pressure_reference=_string(data["pressure_reference"], "scenario.source.pressure_reference"),
        mass_flow_kg_s=(
            _bounded(data["mass_flow_kg_s"], "scenario.source.mass_flow_kg_s", unit="kg/s")
            if has_declared_rate
            else BoundedValue(0.0, unit="kg/s", source="pressure-driven-adapter")
        ),
        opening_area_m2=_bounded(data["opening_area_m2"], "scenario.source.opening_area_m2", unit="m2"),
        discharge_coefficient=_bounded(data["discharge_coefficient"], "scenario.source.discharge_coefficient", unit="1"),
        liquid_fraction=_bounded(data["liquid_fraction"], "scenario.source.liquid_fraction", unit="1"),
        liquid_fraction_basis=_string(data["liquid_fraction_basis"], "scenario.source.liquid_fraction_basis"),
        flash_model=_string(data["flash_model"], "scenario.source.flash_model", declared=True),
        duration_s=duration_value,
        duration_uncertainty=duration_uncertainty,
        metadata={str(key): str(item) for key, item in metadata.items()},
        location_uncertainty_m=location_uncertainty_m,
        direction_uncertainty_m=direction_uncertainty_m,
    )
    if pressure_derivation is None:
        return source
    ambient_pressure, source_id, allow_supercritical_gas = pressure_derivation
    return release_with_pressure_driven_lh2_mass_flow(
        source,
        ambient_pressure_pa=ambient_pressure.nominal,
        ambient_pressure_bounds_pa=(ambient_pressure.lower, ambient_pressure.upper),
        source_id=source_id,
        allow_supercritical_gas=allow_supercritical_gas,
    )


def _weather(value: object) -> WeatherState:
    data = _object(value, "scenario.weather")
    _keys(
        data, "scenario.weather",
        required={"speed_m_s", "direction_deg", "direction_convention", "stability", "reference_height_m"},
    )
    return WeatherState(
        speed_m_s=_bounded(data["speed_m_s"], "scenario.weather.speed_m_s", unit="m/s"),
        direction_deg=_circular_direction(data["direction_deg"], "scenario.weather.direction_deg"),
        direction_convention=_string(data["direction_convention"], "scenario.weather.direction_convention"),
        stability=_string(data["stability"], "scenario.weather.stability"),
        reference_height_m=_number(data["reference_height_m"], "scenario.weather.reference_height_m"),
    )


def _surface(value: object) -> SurfaceBoundary:
    data = _object(value, "scenario.surface")
    _keys(
        data, "scenario.surface",
        required={"heat_transfer_w_m2_k", "surface_temperature_k", "substrate", "evidence_id"},
    )
    return SurfaceBoundary(
        heat_transfer_w_m2_k=_bounded(
            data["heat_transfer_w_m2_k"],
            "scenario.surface.heat_transfer_w_m2_k", unit="W/m2/K",
        ),
        surface_temperature_k=_bounded(
            data["surface_temperature_k"],
            "scenario.surface.surface_temperature_k", unit="K",
        ),
        substrate=_string(data["substrate"], "scenario.surface.substrate", declared=True),
        evidence_id=_string(data["evidence_id"], "scenario.surface.evidence_id", declared=True),
    )


def _sensor(value: object, name: str) -> SensorModel:
    data = _object(value, name)
    _keys(
        data, name,
        required={"position_m", "response_time_s", "gain", "bias_mole_fraction", "averaging_time_s"},
    )
    return SensorModel(
        position_m=_point3(data["position_m"], f"{name}.position_m"),
        response_time_s=_bounded(data["response_time_s"], f"{name}.response_time_s", unit="s"),
        gain=_bounded(data["gain"], f"{name}.gain", unit="1"),
        bias_mole_fraction=_bounded(
            data["bias_mole_fraction"], f"{name}.bias_mole_fraction", unit="mole_fraction",
        ),
        averaging_time_s=_number(data["averaging_time_s"], f"{name}.averaging_time_s"),
    )


def _distributed_source(value: object, index: int) -> FieldDistributedVapourSource:
    """Parse one explicit already-atmospheric global source for strict cases."""
    name = f"distributed_vapour_sources[{index}]"
    data = _object(value, name)
    _keys(
        data,
        name,
        required={"label", "position_m", "vertical_sigma_m", "schedule", "evidence_id"},
        optional={"source_kind", "warnings"},
    )
    position, position_uncertainty = _source_vector(
        data["position_m"], f"{name}.position_m", unit="m",
    )
    sigma = _bounded(data["vertical_sigma_m"], f"{name}.vertical_sigma_m", unit="m")
    schedule_data = _object(data["schedule"], f"{name}.schedule")
    _keys(
        schedule_data,
        f"{name}.schedule",
        required={"time_s", "rate_kg_s", "source_id"},
        optional={"rate_lower_kg_s", "rate_upper_kg_s", "rate_operator"},
    )
    times = schedule_data["time_s"]
    rates = schedule_data["rate_kg_s"]
    if not isinstance(times, list) or not isinstance(rates, list):
        raise ValueError(f"{name}.schedule time_s and rate_kg_s must be arrays")
    parsed_times = tuple(
        _number(item, f"{name}.schedule.time_s[{item_index}]")
        for item_index, item in enumerate(times)
    )
    parsed_rates = tuple(
        _number(item, f"{name}.schedule.rate_kg_s[{item_index}]")
        for item_index, item in enumerate(rates)
    )
    if len(parsed_rates) < 2 or parsed_rates[-1] != 0.0:
        raise ValueError(
            f"{name}.schedule must have a zero final endpoint rate"
        )
    schedule = SourceRateSchedule(
        parsed_times,
        parsed_rates,
        source_id=_string(
            schedule_data["source_id"],
            f"{name}.schedule.source_id",
            declared=True,
        ),
        rate_operator=_string(
            schedule_data.get("rate_operator", "piecewise_constant"),
            f"{name}.schedule.rate_operator",
        ),
    )
    if schedule.released_mass_kg <= 0.0:
        raise ValueError(f"{name}.schedule must release positive mass")
    lower_schedule = upper_schedule = None
    has_lower = "rate_lower_kg_s" in schedule_data
    has_upper = "rate_upper_kg_s" in schedule_data
    if has_lower != has_upper:
        raise ValueError(
            f"{name}.schedule rate_lower_kg_s and rate_upper_kg_s must be supplied together"
        )
    if has_lower and has_upper:
        lower_rates = schedule_data["rate_lower_kg_s"]
        upper_rates = schedule_data["rate_upper_kg_s"]
        if not isinstance(lower_rates, list) or not isinstance(upper_rates, list):
            raise ValueError(
                f"{name}.schedule rate_lower_kg_s and rate_upper_kg_s must be arrays"
            )
        if len(lower_rates) != len(parsed_rates) or len(upper_rates) != len(parsed_rates):
            raise ValueError(
                f"{name}.schedule rate bounds must match the nominal time axis"
            )
        parsed_lower = tuple(
            _number(item, f"{name}.schedule.rate_lower_kg_s[{item_index}]")
            for item_index, item in enumerate(lower_rates)
        )
        parsed_upper = tuple(
            _number(item, f"{name}.schedule.rate_upper_kg_s[{item_index}]")
            for item_index, item in enumerate(upper_rates)
        )
        if any(value < 0.0 for value in (*parsed_lower, *parsed_upper)):
            raise ValueError(f"{name}.schedule rate bounds must be non-negative")
        if any(
            not lower <= nominal <= upper
            for lower, nominal, upper in zip(parsed_lower, parsed_rates, parsed_upper)
        ):
            raise ValueError(
                f"{name}.schedule rate bounds must contain nominal rates"
            )
        if parsed_lower[-1] != 0.0 or parsed_upper[-1] != 0.0:
            raise ValueError(f"{name}.schedule rate bounds must end at zero")
        lower_schedule = SourceRateSchedule(
            parsed_times, parsed_lower,
            source_id=schedule.source_id,
            rate_operator=schedule.rate_operator,
        )
        upper_schedule = SourceRateSchedule(
            parsed_times, parsed_upper,
            source_id=schedule.source_id,
            rate_operator=schedule.rate_operator,
        )
    warnings_value = data.get("warnings", [])
    if not isinstance(warnings_value, list) or any(
        not isinstance(item, str) or not item.strip() for item in warnings_value
    ):
        raise ValueError(f"{name}.warnings must be an array of non-empty strings")
    return FieldDistributedVapourSource(
        label=_string(data["label"], f"{name}.label", declared=True),
        position_m=position,
        vertical_sigma_m=sigma.nominal,
        schedule=schedule,
        evidence_id=_string(data["evidence_id"], f"{name}.evidence_id", declared=True),
        source_kind=_string(
            data.get("source_kind", "declared_atmospheric_vapour"),
            f"{name}.source_kind",
            declared=True,
        ),
        warnings=tuple(warnings_value),
        position_uncertainty_m=position_uncertainty,
        vertical_sigma_uncertainty=sigma,
        lower_schedule=lower_schedule,
        upper_schedule=upper_schedule,
    )


def _scenario(value: object) -> FieldScenario:
    data = _object(value, "scenario")
    _keys(
        data, "scenario",
        required={"source", "weather", "surface", "model_family", "temporal_mode"},
        optional={"sensor"},
    )
    sensor = None if "sensor" not in data or data["sensor"] is None else _sensor(data["sensor"], "scenario.sensor")
    return FieldScenario(
        source=_source(data["source"]), weather=_weather(data["weather"]),
        surface=_surface(data["surface"]), sensor=sensor,
        model_family=_string(data["model_family"], "scenario.model_family"),
        temporal_mode=_string(data["temporal_mode"], "scenario.temporal_mode"),
    )


def _transport(value: object) -> SemiFVConfig:
    data = _object(value, "transport")
    _keys(
        data, "transport",
        required={
            "length_m", "height_m", "nx", "nz", "time_step_s", "duration_s",
            "diffusivity_m2_s", "gravitational_settling_m_s", "source_sigma_m",
        },
    )
    return SemiFVConfig(
        length_m=_number(data["length_m"], "transport.length_m"),
        height_m=_number(data["height_m"], "transport.height_m"),
        nx=_integer(data["nx"], "transport.nx", minimum=4),
        nz=_integer(data["nz"], "transport.nz", minimum=4),
        time_step_s=_number(data["time_step_s"], "transport.time_step_s"),
        duration_s=_number(data["duration_s"], "transport.duration_s"),
        diffusivity_m2_s=_number(data["diffusivity_m2_s"], "transport.diffusivity_m2_s"),
        gravitational_settling_m_s=_number(
            data["gravitational_settling_m_s"], "transport.gravitational_settling_m_s",
        ),
        source_sigma_m=_number(data["source_sigma_m"], "transport.source_sigma_m"),
    )


def _point2(value: object, name: str) -> tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{name} must be a two-value coordinate array")
    return tuple(_number(item, f"{name}[{index}]") for index, item in enumerate(value))


def _obstacle_scalar(
    value: object,
    name: str,
    *,
    unit: str,
    angular: bool = False,
) -> tuple[float, BoundedValue | CircularBoundedValue | None]:
    if isinstance(value, Mapping):
        bound = (
            _circular_direction(value, name)
            if angular else _bounded(value, name, unit=unit)
        )
        return float(bound.nominal), bound
    return _number(value, name), None


def _obstacles(
    value: object,
) -> tuple[
    tuple[AxisAlignedCuboid | OrientedCuboid, ...],
    tuple[FieldObstacleGeometryUncertainty, ...],
]:
    if not isinstance(value, list):
        raise ValueError("obstacles must be a JSON array")
    obstacles = []
    uncertainties = []
    axis_keys = {"x_min_m", "x_max_m", "y_min_m", "y_max_m", "z_min_m", "z_max_m", "label"}
    axis_uncertainty_keys = axis_keys | {"uncertainty_evidence_id"}
    oriented_keys = {
        "center_m", "length_m", "width_m", "z_min_m", "z_max_m", "long_axis_bearing_deg",
        "long_axis_bearing_convention", "label",
    }
    oriented_uncertainty_keys = oriented_keys | {"uncertainty_evidence_id"}
    for index, item in enumerate(value):
        name = f"obstacles[{index}]"
        data = _object(item, name)
        keyset = set(data)
        if keyset == axis_keys or keyset == axis_uncertainty_keys:
            uncertain = keyset == axis_uncertainty_keys
            raw_names = ("x_min_m", "x_max_m", "y_min_m", "y_max_m", "z_min_m", "z_max_m")
            parsed = tuple(
                _obstacle_scalar(data[field], f"{name}.{field}", unit="m")
                for field in raw_names
            )
            flags = tuple(bound is not None for _nominal, bound in parsed)
            if uncertain and not all(flags):
                raise ValueError(
                    f"{name} bounded obstacle geometry must bound every linear dimension"
                )
            if not uncertain and any(flags):
                raise ValueError(f"{name} has an unexpected bounded geometry value")
            obstacle = AxisAlignedCuboid(
                *(item[0] for item in parsed),
                _string(data["label"], f"{name}.label", declared=True),
            )
            obstacles.append(obstacle)
            if uncertain:
                uncertainties.append(FieldObstacleGeometryUncertainty(
                    obstacle,
                    tuple(zip(raw_names, (item[1] for item in parsed))),
                    _string(
                        data["uncertainty_evidence_id"],
                        f"{name}.uncertainty_evidence_id", declared=True,
                    ),
                ))
        elif keyset == oriented_keys or keyset == oriented_uncertainty_keys:
            uncertain = keyset == oriented_uncertainty_keys
            convention = _string(
                data["long_axis_bearing_convention"],
                f"{name}.long_axis_bearing_convention",
            )
            if convention != "math_to":
                raise ValueError(
                    f"{name}.long_axis_bearing_convention must be 'math_to'"
                )
            center = data["center_m"]
            if not isinstance(center, (list, tuple)) or len(center) != 2:
                raise ValueError(f"{name}.center_m must be a two-value coordinate array")
            center_parsed = tuple(
                _obstacle_scalar(item, f"{name}.center_m[{idx}]", unit="m")
                for idx, item in enumerate(center)
            )
            parsed = center_parsed + tuple(
                _obstacle_scalar(data[field], f"{name}.{field}", unit="m")
                for field in ("length_m", "width_m", "z_min_m", "z_max_m")
            ) + (
                _obstacle_scalar(
                    data["long_axis_bearing_deg"],
                    f"{name}.long_axis_bearing_deg",
                    unit="deg", angular=True,
                ),
            )
            flags = tuple(bound is not None for _nominal, bound in parsed)
            if uncertain and not all(flags):
                raise ValueError(
                    f"{name} bounded obstacle geometry must bound every dimension and bearing"
                )
            if not uncertain and any(flags):
                raise ValueError(f"{name} has an unexpected bounded geometry value")
            obstacle = OrientedCuboid(
                (parsed[0][0], parsed[1][0]), parsed[2][0], parsed[3][0],
                parsed[4][0], parsed[5][0], parsed[6][0],
                _string(data["label"], f"{name}.label", declared=True),
            )
            obstacles.append(obstacle)
            if uncertain:
                field_names = (
                    "center_x_m", "center_y_m", "length_m", "width_m",
                    "z_min_m", "z_max_m", "long_axis_bearing_deg",
                )
                uncertainties.append(FieldObstacleGeometryUncertainty(
                    obstacle,
                    tuple(zip(field_names, (item[1] for item in parsed))),
                    _string(
                        data["uncertainty_evidence_id"],
                        f"{name}.uncertainty_evidence_id", declared=True,
                    ),
                ))
        else:
            raise ValueError(
                f"{name} must use exactly the axis-aligned or oriented-cuboid obstacle keys"
            )
    return tuple(obstacles), tuple(uncertainties)


def _sensor_deployments(value: object) -> tuple[FieldSensorDeployment, ...]:
    if not isinstance(value, list):
        raise ValueError("sensor_deployments must be a JSON array")
    deployments = []
    for index, item in enumerate(value):
        name = f"sensor_deployments[{index}]"
        data = _object(item, name)
        _keys(data, name, required={"label", "sensor"})
        deployments.append(FieldSensorDeployment(
            _string(data["label"], f"{name}.label", declared=True),
            _sensor(data["sensor"], f"{name}.sensor"),
        ))
    return tuple(deployments)


def _mixing_closure(value: object) -> StabilityScalarMixingClosure:
    data = _object(value, "stability_mixing_closure")
    _keys(data, "stability_mixing_closure", required={"diffusivity_m2_s", "evidence_id"}, optional={"model_id"})
    diffusivity = _object(data["diffusivity_m2_s"], "stability_mixing_closure.diffusivity_m2_s")
    return StabilityScalarMixingClosure(
        {str(key): _number(item, f"stability_mixing_closure.diffusivity_m2_s.{key}") for key, item in diffusivity.items()},
        evidence_id=_string(data["evidence_id"], "stability_mixing_closure.evidence_id", declared=True),
        model_id=_string(data.get("model_id", "declared_stability_scalar_diffusivity"), "stability_mixing_closure.model_id"),
    )


def _stability_alternatives(value: object) -> FieldStabilityAlternatives:
    data = _object(value, "stability_alternatives")
    _keys(
        data, "stability_alternatives", required={"alternatives", "evidence_id"},
    )
    raw_alternatives = data["alternatives"]
    if not isinstance(raw_alternatives, list):
        raise ValueError("stability_alternatives.alternatives must be a JSON array")
    return FieldStabilityAlternatives(
        tuple(
            _string(item, f"stability_alternatives.alternatives[{index}]", declared=True)
            for index, item in enumerate(raw_alternatives)
        ),
        evidence_id=_string(
            data["evidence_id"], "stability_alternatives.evidence_id", declared=True,
        ),
    )


def _phase_gas_numerics(value: object) -> dict[str, object]:
    data = _object(value, "phase_routing_transport.phase_routing.gas_numerics")
    _keys(
        data,
        "phase_routing_transport.phase_routing.gas_numerics",
        required={
            "maximum_nearfield_distance_m", "radial_points",
            "nearfield_maximum_step_m", "nearfield_relative_tolerance",
            "crosswind_maximum_distance_m", "crosswind_maximum_step_m",
            "puff_time_step_s",
        },
    )
    positive = {
        "maximum_nearfield_distance": _number(
            data["maximum_nearfield_distance_m"],
            "phase_routing_transport.phase_routing.gas_numerics.maximum_nearfield_distance_m",
        ),
        "nearfield_maximum_step": _number(
            data["nearfield_maximum_step_m"],
            "phase_routing_transport.phase_routing.gas_numerics.nearfield_maximum_step_m",
        ),
        "nearfield_relative_tolerance": _number(
            data["nearfield_relative_tolerance"],
            "phase_routing_transport.phase_routing.gas_numerics.nearfield_relative_tolerance",
        ),
        "crosswind_maximum_distance": _number(
            data["crosswind_maximum_distance_m"],
            "phase_routing_transport.phase_routing.gas_numerics.crosswind_maximum_distance_m",
        ),
        "crosswind_maximum_step": _number(
            data["crosswind_maximum_step_m"],
            "phase_routing_transport.phase_routing.gas_numerics.crosswind_maximum_step_m",
        ),
        "puff_time_step": _number(
            data["puff_time_step_s"],
            "phase_routing_transport.phase_routing.gas_numerics.puff_time_step_s",
        ),
    }
    for name, item in positive.items():
        if item <= 0.0:
            raise ValueError(
                "phase_routing_transport.phase_routing.gas_numerics values "
                "must be positive"
            )
    positive["radial_points"] = _integer(
        data["radial_points"],
        "phase_routing_transport.phase_routing.gas_numerics.radial_points",
        minimum=41,
    )
    return positive


def _phase_liquid_numerics(value: object) -> dict[str, object]:
    data = _object(value, "phase_routing_transport.phase_routing.liquid_numerics")
    _keys(
        data,
        "phase_routing_transport.phase_routing.liquid_numerics",
        required={"maximum_droplet_time_s"},
    )
    maximum_time = _number(
        data["maximum_droplet_time_s"],
        "phase_routing_transport.phase_routing.liquid_numerics.maximum_droplet_time_s",
    )
    if maximum_time <= 0.0:
        raise ValueError(
            "phase_routing_transport.phase_routing.liquid_numerics.maximum_droplet_time_s "
            "must be positive"
        )
    return {"maximum_droplet_time_s": maximum_time}


def _phase_droplet_population(
    value: object,
) -> tuple[tuple[DropletClass, ...], str]:
    data = _object(value, "phase_routing_transport.phase_routing.droplet_population")
    _keys(
        data,
        "phase_routing_transport.phase_routing.droplet_population",
        required={"classes", "evidence_id"},
    )
    classes = _phase_droplet_classes(
        data["classes"],
        "phase_routing_transport.phase_routing.droplet_population.classes",
    )
    evidence_id = _string(
        data["evidence_id"],
        "phase_routing_transport.phase_routing.droplet_population.evidence_id",
        declared=True,
    )
    return classes, evidence_id


def _phase_droplet_classes(
    value: object,
    name: str,
) -> tuple[DropletClass, ...]:
    raw_classes = value
    if not isinstance(raw_classes, list) or not raw_classes:
        raise ValueError(
            f"{name} "
            "must be a non-empty JSON array"
        )
    classes = []
    base_name = name
    for index, value in enumerate(raw_classes):
        class_name = f"{base_name}[{index}]"
        item = _object(value, class_name)
        _keys(item, class_name, required={"diameter_m", "mass_fraction"})
        diameter = _number(item["diameter_m"], f"{class_name}.diameter_m")
        fraction = _number(item["mass_fraction"], f"{class_name}.mass_fraction")
        if diameter <= 0.0 or fraction <= 0.0:
            raise ValueError(f"{class_name} values must be positive")
        classes.append(DropletClass(diameter, fraction))
    if not math.isclose(
        sum(item.mass_fraction for item in classes),
        1.0,
        rel_tol=1.0e-10,
        abs_tol=1.0e-12,
    ):
        raise ValueError(f"{base_name} mass fractions must sum to one")
    return tuple(classes)


def _phase_routing_config(
    value: object,
) -> tuple[
    FieldPhaseRoutingConfig,
    tuple[DropletClass, ...],
    str,
    str,
    FieldPhaseRoutingUncertainty | None,
]:
    data = _object(value, "phase_routing_transport.phase_routing")
    _keys(
        data,
        "phase_routing_transport.phase_routing",
        required={
            "post_release_duration_s", "puff_duration_s", "pool_area_m2",
            "pool_time_step_s", "evaporation_coefficient_m2_s",
            "gas_numerics", "liquid_numerics", "droplet_population",
            "evidence_id",
        },
        optional={
            "relative_humidity_pct", "roughness_m", "averaging_time_s",
            "wind_reference_height_m", "pool_model", "uncertainty",
        },
    )
    droplet_population, droplet_evidence_id = _phase_droplet_population(
        data["droplet_population"]
    )
    liquid_options = _phase_liquid_numerics(data["liquid_numerics"])
    liquid_options["droplet_classes"] = droplet_population
    evidence_id = _string(
        data["evidence_id"],
        "phase_routing_transport.phase_routing.evidence_id",
        declared=True,
    )
    config = FieldPhaseRoutingConfig(
        post_release_duration_s=_number(
            data["post_release_duration_s"],
            "phase_routing_transport.phase_routing.post_release_duration_s",
        ),
        puff_duration_s=_number(
            data["puff_duration_s"],
            "phase_routing_transport.phase_routing.puff_duration_s",
        ),
        pool_area_m2=_number(
            data["pool_area_m2"],
            "phase_routing_transport.phase_routing.pool_area_m2",
        ),
        pool_time_step_s=_number(
            data["pool_time_step_s"],
            "phase_routing_transport.phase_routing.pool_time_step_s",
        ),
        evaporation_coefficient_m2_s=_number(
            data["evaporation_coefficient_m2_s"],
            "phase_routing_transport.phase_routing.evaporation_coefficient_m2_s",
        ),
        relative_humidity_pct=_number(
            data.get("relative_humidity_pct", 0.0),
            "phase_routing_transport.phase_routing.relative_humidity_pct",
        ),
        roughness_m=_number(
            data.get("roughness_m", 0.001),
            "phase_routing_transport.phase_routing.roughness_m",
        ),
        averaging_time_s=_number(
            data.get("averaging_time_s", 60.0),
            "phase_routing_transport.phase_routing.averaging_time_s",
        ),
        wind_reference_height_m=_number(
            data.get("wind_reference_height_m", 10.0),
            "phase_routing_transport.phase_routing.wind_reference_height_m",
        ),
        pool_model=_string(
            data.get("pool_model", "dynamic"),
            "phase_routing_transport.phase_routing.pool_model",
        ),
        gas_model_options=_phase_gas_numerics(data["gas_numerics"]),
        phase_model_options=liquid_options,
    )
    phase_uncertainty = None
    if "uncertainty" in data:
        uncertainty_data = _object(
            data["uncertainty"],
            "phase_routing_transport.phase_routing.uncertainty",
        )
        uncertainty_name = "phase_routing_transport.phase_routing.uncertainty"
        if not uncertainty_data:
            raise ValueError(f"{uncertainty_name} must declare at least one bound")
        _keys(
            uncertainty_data,
            uncertainty_name,
            required=set(),
            optional={
                "post_release_duration_s", "puff_duration_s", "pool_area_m2",
                "pool_time_step_s", "evaporation_coefficient_m2_s",
                "droplet_population",
            },
        )
        units = {
            "post_release_duration_s": "s",
            "puff_duration_s": "s",
            "pool_area_m2": "m2",
            "pool_time_step_s": "s",
            "evaporation_coefficient_m2_s": "m2/s",
        }
        parsed: dict[str, BoundedValue] = {}
        for name, unit in units.items():
            if name in uncertainty_data:
                parsed[name] = _bounded(
                    uncertainty_data[name], f"{uncertainty_name}.{name}", unit=unit,
                )
        population_corners = None
        population_evidence_id = None
        if "droplet_population" in uncertainty_data:
            population_name = f"{uncertainty_name}.droplet_population"
            population_data = _object(uncertainty_data["droplet_population"], population_name)
            _keys(
                population_data,
                population_name,
                required={"corners", "evidence_id"},
            )
            raw_corners = population_data["corners"]
            if not isinstance(raw_corners, list) or not raw_corners:
                raise ValueError(f"{population_name}.corners must be a non-empty JSON array")
            population_corners = []
            for index, corner in enumerate(raw_corners):
                corner_name = f"{population_name}.corners[{index}]"
                corner_data = _object(corner, corner_name)
                _keys(corner_data, corner_name, required={"classes"})
                population_corners.append(
                    _phase_droplet_classes(
                        corner_data["classes"], f"{corner_name}.classes",
                    )
                )
            population_evidence_id = _string(
                population_data["evidence_id"],
                f"{population_name}.evidence_id",
                declared=True,
            )
        phase_uncertainty = FieldPhaseRoutingUncertainty(
            **parsed,
            droplet_population_corners=(
                None if population_corners is None else tuple(population_corners)
            ),
            droplet_population_uncertainty_evidence_id=population_evidence_id,
        )
        # Validate the interval anchor and positivity at input time, not only
        # when a later envelope happens to be executed.
        phase_uncertainty.config_corners(config)
    return config, droplet_population, droplet_evidence_id, evidence_id, phase_uncertainty


def _pool_launch(value: object) -> PoolVapourLaunchBoundary:
    data = _object(value, "phase_routing_transport.pool_launch")
    _keys(
        data,
        "phase_routing_transport.pool_launch",
        required={"source_height_m", "closure_id", "evidence_id"},
    )
    return PoolVapourLaunchBoundary(
        source_height_m=_number(
            data["source_height_m"],
            "phase_routing_transport.pool_launch.source_height_m",
        ),
        closure_id=_string(
            data["closure_id"],
            "phase_routing_transport.pool_launch.closure_id",
            declared=True,
        ),
        evidence_id=_string(
            data["evidence_id"],
            "phase_routing_transport.pool_launch.evidence_id",
            declared=True,
        ),
    )


def _phase_routing_transport(value: object) -> FieldPhaseRoutingTransportInput:
    data = _object(value, "phase_routing_transport")
    _keys(
        data,
        "phase_routing_transport",
        required={"phase_routing", "pool_launch", "pool_vertical_sigma_m"},
        optional={"table_nodes", "pool_vertical_sigma_uncertainty"},
    )
    (
        config,
        droplet_population,
        droplet_evidence_id,
        phase_routing_evidence_id,
        phase_uncertainty,
    ) = _phase_routing_config(data["phase_routing"])
    raw_sigma = data["pool_vertical_sigma_m"]
    sigma_uncertainty = None
    if isinstance(raw_sigma, Mapping):
        sigma_uncertainty = _bounded(
            raw_sigma,
            "phase_routing_transport.pool_vertical_sigma_m",
            unit="m",
        )
        sigma_nominal = sigma_uncertainty.nominal
        if "pool_vertical_sigma_uncertainty" in data:
            raise ValueError(
                "phase_routing_transport cannot declare pool_vertical_sigma_m "
                "and pool_vertical_sigma_uncertainty together"
            )
    else:
        sigma_nominal = _number(
            raw_sigma,
            "phase_routing_transport.pool_vertical_sigma_m",
        )
        if "pool_vertical_sigma_uncertainty" in data:
            sigma_uncertainty = _bounded(
                data["pool_vertical_sigma_uncertainty"],
                "phase_routing_transport.pool_vertical_sigma_uncertainty",
                unit="m",
            )
    return FieldPhaseRoutingTransportInput(
        config=config,
        pool_launch=_pool_launch(data["pool_launch"]),
        pool_vertical_sigma_m=sigma_nominal,
        droplet_population=droplet_population,
        droplet_population_evidence_id=droplet_evidence_id,
        phase_routing_evidence_id=phase_routing_evidence_id,
        table_nodes=_integer(
            data.get("table_nodes", 161),
            "phase_routing_transport.table_nodes",
            minimum=4,
        ),
        phase_uncertainty=phase_uncertainty,
        pool_vertical_sigma_uncertainty=sigma_uncertainty,
    )


def _historian_channel(value: object, name: str, *, unit: str) -> HistorianCsvChannel:
    data = _object(value, name)
    _keys(
        data, name,
        required={"value_column", "unit", "source_id", "calibration_evidence_id"},
        optional={
            "response_time_s", "time_offset_s", "lower_column", "upper_column",
            "absolute_half_width", "relative_half_width",
        },
    )
    if data["unit"] != unit:
        raise ValueError(f"{name}.unit must be {unit!r}")
    optional_column = lambda key: (
        None if key not in data or data[key] is None else _string(data[key], f"{name}.{key}")
    )
    optional_number = lambda key: (
        None if key not in data or data[key] is None else _number(data[key], f"{name}.{key}")
    )
    return HistorianCsvChannel(
        value_column=_string(data["value_column"], f"{name}.value_column"),
        unit=unit,
        source_id=_string(data["source_id"], f"{name}.source_id", declared=True),
        calibration_evidence_id=_string(
            data["calibration_evidence_id"], f"{name}.calibration_evidence_id", declared=True,
        ),
        response_time_s=0.0 if "response_time_s" not in data else _number(data["response_time_s"], f"{name}.response_time_s"),
        time_offset_s=0.0 if "time_offset_s" not in data else _number(data["time_offset_s"], f"{name}.time_offset_s"),
        lower_column=optional_column("lower_column"),
        upper_column=optional_column("upper_column"),
        absolute_half_width=optional_number("absolute_half_width"),
        relative_half_width=optional_number("relative_half_width"),
    )


def _measured_history_map(value: object) -> MeasuredHistoryCsvMap:
    data = _object(value, "measured_history.map")
    _keys(
        data, "measured_history.map",
        required={
            "time_s_column", "event_id", "event_evidence_id", "phase_evidence_id",
            "pressure_pa", "temperature_k", "mass_flow_kg_s",
        },
        optional={"liquid_fraction"},
    )
    return MeasuredHistoryCsvMap(
        time_s_column=_string(data["time_s_column"], "measured_history.map.time_s_column"),
        event_id=_string(data["event_id"], "measured_history.map.event_id", declared=True),
        event_evidence_id=_string(
            data["event_evidence_id"], "measured_history.map.event_evidence_id", declared=True,
        ),
        phase_evidence_id=_string(
            data["phase_evidence_id"], "measured_history.map.phase_evidence_id", declared=True,
        ),
        pressure_pa=_historian_channel(data["pressure_pa"], "measured_history.map.pressure_pa", unit="Pa"),
        temperature_k=_historian_channel(data["temperature_k"], "measured_history.map.temperature_k", unit="K"),
        mass_flow_kg_s=_historian_channel(data["mass_flow_kg_s"], "measured_history.map.mass_flow_kg_s", unit="kg/s"),
        liquid_fraction=(
            None if "liquid_fraction" not in data or data["liquid_fraction"] is None
            else _historian_channel(data["liquid_fraction"], "measured_history.map.liquid_fraction", unit="1")
        ),
    )


def _pressure_driven_history_map(value: object) -> PressureDrivenHistoryCsvMap:
    data = _object(value, "pressure_driven_history.map")
    _keys(
        data, "pressure_driven_history.map",
        required={
            "time_s_column", "event_id", "event_evidence_id", "phase_evidence_id",
            "pressure_pa", "temperature_k",
        },
        optional={"liquid_fraction"},
    )
    return PressureDrivenHistoryCsvMap(
        time_s_column=_string(data["time_s_column"], "pressure_driven_history.map.time_s_column"),
        event_id=_string(data["event_id"], "pressure_driven_history.map.event_id", declared=True),
        event_evidence_id=_string(
            data["event_evidence_id"], "pressure_driven_history.map.event_evidence_id", declared=True,
        ),
        phase_evidence_id=_string(
            data["phase_evidence_id"], "pressure_driven_history.map.phase_evidence_id", declared=True,
        ),
        pressure_pa=_historian_channel(
            data["pressure_pa"], "pressure_driven_history.map.pressure_pa", unit="Pa",
        ),
        temperature_k=_historian_channel(
            data["temperature_k"], "pressure_driven_history.map.temperature_k", unit="K",
        ),
        liquid_fraction=(
            None if "liquid_fraction" not in data or data["liquid_fraction"] is None
            else _historian_channel(
                data["liquid_fraction"], "pressure_driven_history.map.liquid_fraction", unit="1",
            )
        ),
    )


def _measured_history_quality(
    value: object,
    *,
    context: str = "measured_history",
) -> MeasuredHistoryQualityCriteria:
    data = _object(value, f"{context}.quality_criteria")
    _keys(
        data, f"{context}.quality_criteria",
        required={
            "maximum_sample_interval_s", "maximum_response_time_s",
            "maximum_absolute_time_offset_s", "maximum_relative_half_width", "evidence_id",
        },
    )
    return MeasuredHistoryQualityCriteria(
        maximum_sample_interval_s=_number(
            data["maximum_sample_interval_s"], f"{context}.quality_criteria.maximum_sample_interval_s",
        ),
        maximum_response_time_s=_number(
            data["maximum_response_time_s"], f"{context}.quality_criteria.maximum_response_time_s",
        ),
        maximum_absolute_time_offset_s=_number(
            data["maximum_absolute_time_offset_s"], f"{context}.quality_criteria.maximum_absolute_time_offset_s",
        ),
        maximum_relative_half_width=_number(
            data["maximum_relative_half_width"], f"{context}.quality_criteria.maximum_relative_half_width",
        ),
        evidence_id=_string(
            data["evidence_id"], f"{context}.quality_criteria.evidence_id", declared=True,
        ),
    )


def _attach_measured_history(
    request: FieldSemiFVRequest,
    value: object,
    *,
    base_directory: Path | None,
) -> tuple[FieldSemiFVRequest, ImportedMeasuredReleaseHistory, FieldMeasuredFlashSchedule]:
    if base_directory is None:
        raise ValueError("measured_history requires a case file path so its CSV path can be resolved")
    data = _object(value, "measured_history")
    _keys(data, "measured_history", required={"csv_path", "map", "quality_criteria"})
    csv_path = Path(_string(data["csv_path"], "measured_history.csv_path"))
    if not csv_path.is_absolute():
        csv_path = base_directory / csv_path
    imported = read_measured_history_csv(csv_path, _measured_history_map(data["map"]))
    quality = _measured_history_quality(data["quality_criteria"])
    schedule = direct_vapour_schedule_from_imported_history(
        request.scenario.source, imported, quality_criteria=quality,
    )
    return request_with_measured_flash_schedule(request, schedule), imported, schedule


def _attach_pressure_driven_history(
    request: FieldSemiFVRequest,
    value: object,
    *,
    base_directory: Path | None,
) -> tuple[FieldSemiFVRequest, ImportedPressureDrivenMeasuredHistory, FieldMeasuredFlashSchedule]:
    if base_directory is None:
        raise ValueError(
            "pressure_driven_history requires a case file path so its CSV path can be resolved"
        )
    data = _object(value, "pressure_driven_history")
    _keys(
        data, "pressure_driven_history", required={"csv_path", "map", "quality_criteria"},
    )
    csv_path = Path(_string(data["csv_path"], "pressure_driven_history.csv_path"))
    if not csv_path.is_absolute():
        csv_path = base_directory / csv_path
    imported = read_pressure_driven_history_csv(
        csv_path, _pressure_driven_history_map(data["map"]),
    )
    quality = _measured_history_quality(
        data["quality_criteria"], context="pressure_driven_history",
    )
    schedule = direct_vapour_schedule_from_imported_pressure_driven_history(
        request.scenario.source, imported, quality_criteria=quality,
    )
    attached = request_with_pressure_driven_history_schedule(
        request, imported.history, schedule,
    )
    return attached, imported, schedule


def field_screening_case_from_mapping(
    value: object,
    *,
    base_directory: str | Path | None = None,
) -> FieldScreeningCase:
    """Create one field case from a strict JSON-compatible mapping.

    The base schema re-flashes the declared nominal source, with a controlled
    ``measured_history`` or ``pressure_driven_history`` object as its sole schedule exception.
    The dedicated
    phase-routing schema can instead invoke the conservative rainout/pool
    ledger plus an explicit ground-level scalar launch. Arbitrary post-flash
    schedules, imported pool ledgers and jet handoffs remain outside both
    formats so their provenance gates cannot be bypassed.
    """
    data = _object(value, "field screening input")
    schema = data.get("schema")
    supported_schemas = {
        FIELD_SCREENING_INPUT_SCHEMA,
        FIELD_PHASE_ROUTING_SCREENING_INPUT_SCHEMA,
    }
    if schema not in supported_schemas:
        raise ValueError(
            "field screening input schema must be one of "
            + ", ".join(repr(item) for item in sorted(supported_schemas))
        )
    required = {
        "schema", "scenario", "coordinate_reference", "transport",
        "lh2_validation_available",
    }
    optional = {
        "obstacles", "ambient_temperature_k", "ambient_pressure_pa",
        "ambient_air_density_kg_m3", "ambient_temperature_uncertainty_k",
        "ambient_pressure_uncertainty_pa", "ambient_air_density_uncertainty_kg_m3",
        "lateral_sensor_tolerance_m",
        "lateral_source_tolerance_m", "post_release_duration_s",
        "sensor_deployments", "stability_mixing_closure", "stability_alternatives",
        "conditional_review",
        "validation_evidence",
    }
    if schema == FIELD_SCREENING_INPUT_SCHEMA:
        optional.update({
            "measured_history", "pressure_driven_history", "distributed_vapour_sources",
        })
    else:
        required.add("phase_routing_transport")
    _keys(data, "field screening input", required=required, optional=optional)
    validation_flag = _boolean(data["lh2_validation_available"], "lh2_validation_available")
    validation_evidence = (
        None if "validation_evidence" not in data
        else _validation_evidence(data["validation_evidence"])
    )
    if validation_flag and validation_evidence is None:
        raise ValueError(
            "lh2_validation_available=True requires validation_evidence"
        )
    optional_number = lambda key, default: (
        default if key not in data or data[key] is None else _number(data[key], key)
    )
    scenario = _scenario(data["scenario"])
    transport = _transport(data["transport"])
    conditional_review = (
        None if "conditional_review" not in data
        else _conditional_review(data["conditional_review"])
    )
    phase_routing_transport = (
        None if schema == FIELD_SCREENING_INPUT_SCHEMA
        else _phase_routing_transport(data["phase_routing_transport"])
    )
    distributed_sources = ()
    if "distributed_vapour_sources" in data:
        raw_sources = data["distributed_vapour_sources"]
        if not isinstance(raw_sources, list):
            raise ValueError("distributed_vapour_sources must be an array")
        distributed_sources = tuple(
            _distributed_source(item, index)
            for index, item in enumerate(raw_sources)
        )
    if phase_routing_transport is None:
        post_release_duration_s = optional_number("post_release_duration_s", 0.0)
    else:
        post_release_duration_s = phase_routing_transport.config.post_release_duration_s
        if "post_release_duration_s" in data and not math.isclose(
            _number(data["post_release_duration_s"], "post_release_duration_s"),
            post_release_duration_s,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        ):
            raise ValueError(
                "post_release_duration_s must match "
                "phase_routing_transport.phase_routing.post_release_duration_s"
            )
    if scenario.temporal_mode == "transient":
        required_post_release_duration_s = post_release_duration_s
        if phase_routing_transport is not None and phase_routing_transport.phase_uncertainty is not None:
            required_post_release_duration_s = phase_routing_transport.phase_uncertainty.max_value(
                phase_routing_transport.config, "post_release_duration_s",
            )
        source_duration_s = scenario.source.duration_s
        source_duration_label = "source.duration_s"
        if scenario.source.duration_uncertainty is not None:
            source_duration_s = scenario.source.duration_uncertainty.upper
            source_duration_label = "maximum declared source duration_s"
        expected_duration_s = source_duration_s + required_post_release_duration_s
        if not math.isclose(
            transport.duration_s, expected_duration_s, rel_tol=1.0e-9, abs_tol=1.0e-12,
        ):
            duration_label = (
                "maximum declared post_release_duration_s"
                if phase_routing_transport is not None
                and phase_routing_transport.phase_uncertainty is not None
                else "post_release_duration_s"
            )
            raise ValueError(
                "transport.duration_s must equal "
                f"{source_duration_label} plus {duration_label} "
                "for a transient field case"
            )
    parsed_obstacles = ((), ()) if "obstacles" not in data else _obstacles(data["obstacles"])
    ambient_temperature_k, ambient_temperature_uncertainty_k = _ambient_bound(
        data,
        nominal_key="ambient_temperature_k",
        uncertainty_key="ambient_temperature_uncertainty_k",
        default=295.0,
        unit="K",
    )
    ambient_pressure_pa, ambient_pressure_uncertainty_pa = _ambient_bound(
        data,
        nominal_key="ambient_pressure_pa",
        uncertainty_key="ambient_pressure_uncertainty_pa",
        default=101325.0,
        unit="Pa",
    )
    ambient_air_density_kg_m3, ambient_air_density_uncertainty_kg_m3 = _ambient_bound(
        data,
        nominal_key="ambient_air_density_kg_m3",
        uncertainty_key="ambient_air_density_uncertainty_kg_m3",
        default=1.2,
        unit="kg/m3",
    )
    request = FieldSemiFVRequest(
        scenario=scenario,
        coordinate_reference=_coordinate_reference(data["coordinate_reference"]),
        transport=transport,
        obstacles=parsed_obstacles[0],
        obstacle_geometry_uncertainty=parsed_obstacles[1],
        ambient_temperature_k=ambient_temperature_k,
        ambient_pressure_pa=ambient_pressure_pa,
        ambient_air_density_kg_m3=ambient_air_density_kg_m3,
        ambient_temperature_uncertainty_k=ambient_temperature_uncertainty_k,
        ambient_pressure_uncertainty_pa=ambient_pressure_uncertainty_pa,
        ambient_air_density_uncertainty_kg_m3=ambient_air_density_uncertainty_kg_m3,
        lateral_sensor_tolerance_m=optional_number("lateral_sensor_tolerance_m", None),
        lateral_source_tolerance_m=optional_number("lateral_source_tolerance_m", None),
        lh2_validation_available=validation_flag,
        validation_evidence=validation_evidence,
        post_release_duration_s=post_release_duration_s,
        distributed_vapour_sources=distributed_sources,
        sensor_deployments=() if "sensor_deployments" not in data else _sensor_deployments(data["sensor_deployments"]),
        stability_mixing_closure=(
            None if "stability_mixing_closure" not in data
            else _mixing_closure(data["stability_mixing_closure"])
        ),
        stability_alternatives=(
            None if "stability_alternatives" not in data
            else _stability_alternatives(data["stability_alternatives"])
        ),
    )
    if phase_routing_transport is not None:
        return FieldScreeningCase(
            request,
            source_request=request,
            phase_routing_transport=phase_routing_transport,
            conditional_review=conditional_review,
        )
    if "measured_history" in data and "pressure_driven_history" in data:
        raise ValueError(
            "field screening input cannot declare both measured_history and pressure_driven_history"
        )
    if "measured_history" not in data and "pressure_driven_history" not in data:
        return FieldScreeningCase(
            request,
            source_request=request,
            conditional_review=conditional_review,
        )
    if "pressure_driven_history" in data:
        attached, imported, schedule = _attach_pressure_driven_history(
            request, data["pressure_driven_history"],
            base_directory=None if base_directory is None else Path(base_directory),
        )
    else:
        attached, imported, schedule = _attach_measured_history(
            request, data["measured_history"],
            base_directory=None if base_directory is None else Path(base_directory),
        )
    return FieldScreeningCase(
        attached, source_request=request, imported_history=imported,
        measured_schedule=schedule, conditional_review=conditional_review,
    )


def field_screening_request_from_mapping(value: object) -> FieldSemiFVRequest:
    """Create the base field request from an in-memory non-historian case."""
    return field_screening_case_from_mapping(value).request


def read_field_screening_case_json(path: str | Path) -> FieldScreeningCase:
    """Read one UTF-8 case file; parse errors retain the input file location."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    try:
        value = strict_json_loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid field-screening JSON at {source}: {error.msg}") from error
    return field_screening_case_from_mapping(value, base_directory=source.parent)


def read_field_screening_request_json(path: str | Path) -> FieldSemiFVRequest:
    """Read one case and return its prepared base field request."""
    return read_field_screening_case_json(path).request


__all__ = [
    "FIELD_SCREENING_INPUT_SCHEMA", "FIELD_PHASE_ROUTING_SCREENING_INPUT_SCHEMA",
    "FieldPhaseRoutingTransportInput", "FieldPhaseRoutingUncertainty",
    "FieldConditionalReviewAuthorization",
    "FieldScreeningCase",
    "field_screening_case_from_mapping", "field_screening_request_from_mapping",
    "read_field_screening_case_json", "read_field_screening_request_json",
]
