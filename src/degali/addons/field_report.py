"""Portable, failure-safe reports for field semi-FV screening results.

The report is an audit record, not a claim that every conditional result is
fit for an operational clearance decision.  In particular, a blocked result
is serialised with its reasons and without fabricated field or sensor values.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .field_contracts import FieldApplicability, FieldScenario, SensorModel
from .field_distributed_source import FieldDistributedVapourSource
from .field_droplet_launch import FieldDropletHandoffRefinementStudy
from .field_history import (
    FieldJointMeasuredHistoryScreeningEnvelope,
    FieldMeasuredHistoryScreeningEnvelope,
)
from .field_superposition import FieldSensorSuperposition
from .field_workflow import (
    FieldSemiFVEnvelope,
    FieldSemiFVRefinementResult,
    FieldSemiFVScreeningResult,
    FieldSensorArrayUncertaintyEnvelope,
)


FIELD_SCREENING_REPORT_SCHEMA = "degali.field-screening.v1"
FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA = "degali.field-sensor-superposition.v1"
FIELD_DROPLET_HANDOFF_REFINEMENT_REPORT_SCHEMA = "degali.field-droplet-handoff-refinement.v1"
FIELD_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA = "degali.field-measured-history-envelope.v1"
FIELD_JOINT_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA = "degali.field-joint-measured-history-envelope.v1"
FIELD_SENSOR_ARRAY_UNCERTAINTY_ENVELOPE_REPORT_SCHEMA = "degali.field-sensor-array-uncertainty-envelope.v1"
FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA = "degali.field-semi-fv-envelope.v1"


def _json_native(value: Any) -> Any:
    """Convert NumPy/dataclass residue to the exact JSON value shape."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_json_native(item) for item in value.tolist()]
    if is_dataclass(value):
        return _json_native(asdict(value))
    if isinstance(value, dict):
        return {str(key): _json_native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_native(item) for item in value]
    return value


def _applicability_record(value: FieldApplicability) -> dict[str, object]:
    return {
        "status": value.status,
        "reasons": list(value.reasons),
        "warnings": list(value.warnings),
        "uncertainty_complete": value.uncertainty_complete,
    }


def _sensor_record(sensor: SensorModel | None) -> dict[str, object] | None:
    if sensor is None:
        return None
    return {
        "position_m": list(sensor.position_m),
        "averaging_time_s": sensor.averaging_time_s,
        "uncertainty": {
            name: value.as_dict() for name, value in sensor.uncertainty_fields().items()
        },
    }


def _scenario_record(scenario: FieldScenario) -> dict[str, object]:
    record = scenario.as_dict()
    record["sensor"] = _sensor_record(scenario.sensor)
    return record


def _distributed_source_record(
    source: FieldDistributedVapourSource,
) -> dict[str, object]:
    """Serialize a supplemental source, including retained geometry bounds."""
    record: dict[str, object] = {
        "label": source.label,
        "source_kind": source.source_kind,
        "position_m": list(source.position_m),
        "vertical_sigma_m": source.vertical_sigma_m,
        "schedule_source_id": source.schedule.source_id,
        "schedule_rate_operator": source.schedule.rate_operator,
        "schedule_duration_s": source.schedule.duration_s,
        "released_mass_kg": source.schedule.released_mass_kg,
        "evidence_id": source.evidence_id,
        "warnings": list(source.warnings),
    }
    if source.position_uncertainty_m is not None:
        record["position_uncertainty_m"] = [
            value.as_dict() for value in source.position_uncertainty_m
        ]
    if source.vertical_sigma_uncertainty is not None:
        record["vertical_sigma_uncertainty"] = (
            source.vertical_sigma_uncertainty.as_dict()
        )
    if source.has_schedule_uncertainty:
        record["schedule_uncertainty"] = {
            "lower_rate_kg_s": list(source.lower_schedule.rate_kg_s),
            "upper_rate_kg_s": list(source.upper_schedule.rate_kg_s),
            "interpretation": "deterministic source-rate corners on the declared common time axis; not a probability interval",
        }
    return record


def _transport_input_record(result: FieldSemiFVScreeningResult) -> dict[str, object]:
    config = result.request.transport
    schedule = result.request.direct_vapour_schedule
    closure = result.request.stability_mixing_closure
    alternatives = result.request.stability_alternatives
    history_quality = result.request.measured_history_quality
    return {
        "length_m": config.length_m,
        "height_m": config.height_m,
        "nx": config.nx,
        "nz": config.nz,
        "requested_time_step_s": config.time_step_s,
        "requested_duration_s": config.duration_s,
        "post_release_duration_s": result.request.post_release_duration_s,
        "ambient_boundary": {
            "temperature_k": result.request.ambient_temperature_k,
            "pressure_pa": result.request.ambient_pressure_pa,
            "air_density_kg_m3": result.request.ambient_air_density_kg_m3,
            "uncertainty": {
                name: value.as_dict()
                for name, value in result.request.ambient_uncertainty_fields().items()
            },
        },
        "diffusivity_m2_s": config.diffusivity_m2_s,
        "stability_mixing_closure": (
            None if closure is None else {
                "model_id": closure.model_id,
                "evidence_id": closure.evidence_id,
                "selected_stability": result.request.scenario.weather.stability,
                "selected_diffusivity_m2_s": closure.selected_diffusivity_m2_s(
                    result.request.scenario.weather.stability
                ),
            }
        ),
        "stability_alternatives": (
            None if alternatives is None else {
                "alternatives": list(alternatives.alternatives),
                "classes_evaluated_with_nominal": list(
                    alternatives.classes_for(result.request.scenario.weather.stability)
                ),
                "evidence_id": alternatives.evidence_id,
                "uncertainty_resolved_for_this_case": (
                    result.request.stability_uncertainty_resolved
                ),
            }
        ),
        "gravitational_settling_m_s": config.gravitational_settling_m_s,
        "source_sigma_m": config.source_sigma_m,
        "additional_sensor_deployments": [
            {"label": deployment.label, "sensor": _sensor_record(deployment.sensor)}
            for deployment in result.request.sensor_deployments
        ],
        "declared_direct_vapour_schedule": (
            None if schedule is None else {
                "source_id": schedule.source_id,
                "rate_operator": schedule.rate_operator,
                "duration_s": schedule.duration_s,
                "released_mass_kg": schedule.released_mass_kg,
                "effective_area_m2": result.request.direct_vapour_effective_area_m2,
                "warnings": list(result.request.direct_vapour_warnings),
            }
        ),
        "measured_history_quality": (
            None if history_quality is None else {
                "approved": history_quality.approved,
                "criteria": asdict(history_quality.criteria),
                "channel_metrics": [
                    {
                        "channel": metric[0],
                        "maximum_sample_interval_s": metric[1],
                        "response_time_s": metric[2],
                        "absolute_time_offset_s": metric[3],
                        "maximum_relative_half_width": metric[4],
                    }
                    for metric in history_quality.channel_metrics
                ],
                "reasons": list(history_quality.reasons),
            }
        ),
        "measured_history_provenance": dict(result.request.measured_history_provenance),
        "measured_history_source_uncertainty_resolved": (
            result.request.measured_history_source_uncertainty_resolved
        ),
        "phase_routing_uncertainty_resolved": (
            result.request.phase_routing_uncertainty_resolved
        ),
        "distributed_vapour_sources": [
            _distributed_source_record(source)
            for source in result.request.distributed_vapour_sources
        ],
    }


def _wind_history_record(result: FieldSemiFVScreeningResult) -> dict[str, object] | None:
    declared = result.request.wind_history
    if declared is None:
        return None
    time, _speed, _direction = declared.history.arrays()
    assessment = result.wind_history_assessment
    record_speed, record_direction, speed_difference, direction_difference = (
        declared.nominal_weather_difference(result.request.scenario.weather)
    )
    return {
        "source_id": declared.source_id,
        "evidence_id": declared.evidence_id,
        "sample_count": len(time),
        "time_coverage_s": [float(time[0]), float(time[-1])],
        "maximum_direction_span_deg": declared.maximum_direction_span_deg,
        "maximum_speed_range_fraction": declared.maximum_speed_range_fraction,
        "maximum_nominal_speed_relative_difference": (
            declared.maximum_nominal_speed_relative_difference
        ),
        "maximum_nominal_direction_difference_deg": (
            declared.maximum_nominal_direction_difference_deg
        ),
        "record_mean_speed_m_s": record_speed,
        "record_mean_wind_to_direction_deg": record_direction,
        "nominal_speed_relative_difference": speed_difference,
        "nominal_direction_difference_deg": direction_difference,
        "steady_wind_assessment": (
            None if assessment is None else asdict(assessment)
        ),
        "use": (
            "static-wind applicability gate only; transport retains the "
            "declared scenario nominal wind after a measured-record consistency check"
        ),
    }


def _jet_handoff_record(result: FieldSemiFVScreeningResult) -> dict[str, object] | None:
    handoff = result.request.jet_scalar_handoff
    if handoff is None:
        return None
    location = None
    if result.wind_frame is not None:
        location = [
            result.wind_frame.origin_x_m,
            result.wind_frame.origin_y_m,
            result.request.scenario.source.location_m[2] + handoff.elevation_offset_m,
        ]
    return {
        "source_id": handoff.source_id,
        "evidence_id": handoff.evidence_id,
        "model_id": handoff.model_id,
        "streamline_distance_m": handoff.streamline_distance_m,
        "downwind_offset_m": handoff.downwind_offset_m,
        "elevation_offset_m": handoff.elevation_offset_m,
        "resolved_location_m": location,
        "vertical_sigma_m": handoff.vertical_sigma_m,
        "hydrogen_mass_flow_kg_s": handoff.hydrogen_mass_flow_kg_s,
        "species_flux_relative_residual": handoff.species_flux_relative_residual,
        "maximum_flash_rate_relative_residual": handoff.maximum_flash_rate_relative_residual,
        "maximum_wind_misalignment_deg": handoff.maximum_wind_misalignment_deg,
        "warnings": list(handoff.warnings),
    }


def _obstacle_projection_record(
    projection: object,
) -> dict[str, object]:
    """Serialise one wind-plane obstacle projection without inferred wakes."""
    # Kept structural rather than importing the geometry class: reports also
    # serialise an already-completed result from a caller's older process.
    return {
        "status": projection.status,
        "reasons": list(projection.reasons),
        "warnings": list(projection.warnings),
        "downwind_bounds_m": list(projection.downwind_bounds_m),
        "crosswind_bounds_m": list(projection.crosswind_bounds_m),
        "local_obstacle": (
            None if projection.obstacle is None else asdict(projection.obstacle)
        ),
    }


def _declared_obstacle_record(obstacle: object) -> dict[str, object]:
    """Preserve pre-projection global obstacle geometry in every report."""
    if not is_dataclass(obstacle):
        raise TypeError("declared field obstacle must be a dataclass")
    return {
        "geometry_type": type(obstacle).__name__,
        "geometry": asdict(obstacle),
    }


def field_screening_report(result: FieldSemiFVScreeningResult) -> dict[str, object]:
    """Return a JSON-serialisable audit report for one field screen."""
    if not isinstance(result, FieldSemiFVScreeningResult):
        raise TypeError("result must be a FieldSemiFVScreeningResult")
    report: dict[str, Any] = {
        "schema": FIELD_SCREENING_REPORT_SCHEMA,
        "scenario": _scenario_record(result.request.scenario),
        "coordinate_reference": (
            None if result.request.coordinate_reference is None
            else asdict(result.request.coordinate_reference)
        ),
        "validation_evidence": (
            None if result.request.validation_evidence is None
            else result.request.validation_evidence.as_record()
        ),
        "applicability": _applicability_record(result.applicability),
        "transport_input": _transport_input_record(result),
        "wind_history": _wind_history_record(result),
        "jet_scalar_handoff": _jet_handoff_record(result),
        "wind_frame": (
            None if result.wind_frame is None else {
                "origin_m": [result.wind_frame.origin_x_m, result.wind_frame.origin_y_m],
                "wind_to_math_radians": result.wind_frame.direction_rad,
                "wind_to_local_math_radians": result.wind_frame.direction_rad,
                "wind_to_earth_math_radians": result.request.scenario.weather.wind_to_math_radians(),
            }
        ),
    }
    preparation = result.source_preparation
    report["source_preparation"] = None if preparation is None else {
        "applicability": _applicability_record(preparation.applicability),
        "flash_closure_diagnostics": dict(preparation.flash_closure_diagnostics),
        "flash": (
            None if preparation.flash_result is None else {
                "model": preparation.flash_result.model,
                "effective_area_m2": preparation.flash_result.effective_area_m2,
                "effective_diameter_m": preparation.flash_result.effective_diameter_m,
                "direct_vapour_mass_flow_kg_s": preparation.flash_result.flash.vapour_mass_flow,
                "postflash_liquid_mass_flow_kg_s": preparation.flash_result.flash.liquid_mass_flow,
                "mass_residual_kg_s": preparation.flash_result.mass_residual_kg_s,
                "phase_mass_residual": preparation.flash_result.phase_mass_residual,
                "momentum_residual": preparation.flash_result.momentum_residual,
                "energy_residual": preparation.flash_result.energy_residual,
                "conservative": preparation.flash_result.conservative,
                "closure_tolerances": preparation.flash_result.closure_tolerances,
                "closure_warnings": list(preparation.flash_result.closure_warnings),
                "warnings": list(preparation.flash_result.warnings),
            }
        ),
    }
    projection = result.obstacle_projection
    declared = (
        (result.request.obstacle,) if result.request.obstacle is not None
        else result.request.obstacles
    )
    report["declared_global_obstacles"] = [
        _declared_obstacle_record(item) for item in declared
    ]
    report["declared_obstacle_geometry_uncertainty"] = [
        item.as_record() for item in result.request.obstacle_geometry_uncertainty
    ]
    report["obstacle_projection"] = (
        None if projection is None else _obstacle_projection_record(projection)
    )
    projections = result.obstacle_projections or (
        () if projection is None else (projection,)
    )
    report["obstacle_projections"] = [
        _obstacle_projection_record(item) for item in projections
    ]
    transport = result.transport
    report["transport_result"] = None if transport is None else {
        "diagnostics": asdict(transport.diagnostics),
        "maximum_concentration_kg_m3": transport.maximum_concentration_kg_m3,
        "receptor_count": len(transport.receptor_traces),
    }
    trace = result.sensor_trace
    report["sensor_result"] = None if trace is None else {
        "time_start_s": float(trace.time_s[0]),
        "time_end_s": float(trace.time_s[-1]),
        "maximum_true_mole_fraction": float(np.max(trace.true_mole_fraction)),
        "maximum_indicated_mole_fraction": float(np.max(trace.indicated_mole_fraction)),
        "final_indicated_mole_fraction": float(trace.indicated_mole_fraction[-1]),
        "ambient_air_density_kg_m3": trace.ambient_air_density_kg_m3,
    }
    report["sensor_results"] = [
        {
            "label": sensor_result.label,
            "position_m": list(sensor_result.sensor.position_m),
            "withheld": sensor_result.trace is None,
            "warnings": list(sensor_result.warnings),
            "result": None if sensor_result.trace is None else {
                "time_start_s": float(sensor_result.trace.time_s[0]),
                "time_end_s": float(sensor_result.trace.time_s[-1]),
                "maximum_true_mole_fraction": float(np.max(sensor_result.trace.true_mole_fraction)),
                "maximum_indicated_mole_fraction": float(np.max(sensor_result.trace.indicated_mole_fraction)),
                "final_indicated_mole_fraction": float(sensor_result.trace.indicated_mole_fraction[-1]),
            },
        }
        for sensor_result in result.sensor_results
    ]
    # ``allow_nan=False`` in the writer is the final guard, but this avoids a
    # vague JSON failure if a future result type adds an invalid scalar.
    for section in (report["transport_result"], report["sensor_result"]):
        if section is not None:
            for value in section.values():
                if isinstance(value, float) and not math.isfinite(value):
                    raise ValueError("field report contains a non-finite numeric result")
    return _json_native(report)


def write_field_screening_report(
    result: FieldSemiFVScreeningResult,
    path: str | Path,
) -> Path:
    """Write one new report explicitly chosen by the caller and return its path.

    Report artifacts are evidence records, so an existing path is never
    replaced implicitly.  Exclusive creation also keeps this boundary
    fail-safe if two writers race for the same destination.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(field_screening_report(result), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n"
    )
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(payload)
    return destination


def field_droplet_handoff_refinement_report(
    study: FieldDropletHandoffRefinementStudy,
) -> dict[str, object]:
    """Return the auditable trajectory-to-field-source refinement record."""
    if not isinstance(study, FieldDropletHandoffRefinementStudy):
        raise TypeError("study must be a FieldDropletHandoffRefinementStudy")
    return _json_native({
        "schema": FIELD_DROPLET_HANDOFF_REFINEMENT_REPORT_SCHEMA,
        "droplet_boundary": asdict(study.droplet_boundary),
        "launch": asdict(study.launch),
        "release_duration_s": study.release_duration_s,
        "maximum_droplet_time_s": study.maximum_droplet_time_s,
        "reference_trajectory_segment_duration_s": (
            study.reference_trajectory_segment_duration_s
        ),
        "relative_tolerance": study.relative_tolerance,
        "converged": study.converged,
        "warnings": list(study.warnings),
        "cases": [
            {
                "trajectory_segment_duration_s": case.trajectory_segment_duration_s,
                "source_count": case.source_count,
                "in_flight_vapour_mass_kg": case.in_flight_vapour_mass_kg,
                "sensor_true_peak_mole_fraction": case.sensor_true_peak_mole_fraction,
                "sensor_true_dose_mole_fraction_s": case.sensor_true_dose_mole_fraction_s,
                "applicability": _applicability_record(case.screening.applicability),
                "field_screening": field_screening_report(case.screening),
            }
            for case in study.cases
        ],
        "peak_relative_changes": [
            {"trajectory_segment_duration_s": duration, "relative_change": change}
            for duration, change in study.peak_relative_changes
        ],
        "dose_relative_changes": [
            {"trajectory_segment_duration_s": duration, "relative_change": change}
            for duration, change in study.dose_relative_changes
        ],
    })


def field_measured_history_envelope_report(
    envelope: FieldMeasuredHistoryScreeningEnvelope,
) -> dict[str, object]:
    """Serialise source-history uncertainty cases with their field outcomes."""
    if not isinstance(envelope, FieldMeasuredHistoryScreeningEnvelope):
        raise TypeError("envelope must be a FieldMeasuredHistoryScreeningEnvelope")
    quality = envelope.quality_assessment
    return _json_native({
        "schema": FIELD_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA,
        "quality_assessment": {
            "approved": quality.approved,
            "criteria": asdict(quality.criteria),
            "channel_metrics": [
                {
                    "channel": metric[0],
                    "maximum_sample_interval_s": metric[1],
                    "response_time_s": metric[2],
                    "absolute_time_offset_s": metric[3],
                    "maximum_relative_half_width": metric[4],
                }
                for metric in quality.channel_metrics
            ],
            "reasons": list(quality.reasons),
        },
        "completed_case_count": envelope.completed_case_count,
        "blocked_case_count": envelope.blocked_case_count,
        "property_table_used": envelope.property_table_used,
        "warnings": list(envelope.warnings),
        "cases": [
            {
                "selection": dict(case.selection),
                "direct_vapour_mass_kg": case.schedule.direct_vapour_mass_kg,
                "unrouted_postflash_liquid_mass_kg": (
                    case.schedule.unrouted_postflash_liquid_mass_kg
                ),
                "field_screening": field_screening_report(case.result),
            }
            for case in envelope.cases
        ],
    })


def field_semi_fv_envelope_report(envelope: FieldSemiFVEnvelope) -> dict[str, object]:
    """Serialise deterministic nominal-source/weather/sensor corner results.

    This report is intentionally not an operational decision: every corner
    needs its own numerical-refinement and scope review before a caller can
    make a bounded operational use of the envelope.
    """
    if not isinstance(envelope, FieldSemiFVEnvelope):
        raise TypeError("envelope must be a FieldSemiFVEnvelope")
    return _json_native({
        "schema": FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA,
        "property_table_used": envelope.property_table_used,
        "completed_case_count": envelope.completed_case_count,
        "blocked_case_count": envelope.blocked_case_count,
        "atmospheric_source_schedule": (
            None if envelope.atmospheric_source_schedule is None
            else envelope.atmospheric_source_schedule.as_record()
        ),
        "declared_distributed_vapour_sources": [
            _distributed_source_record(source)
            for source in envelope.request.distributed_vapour_sources
        ],
        "declared_obstacle_geometry_uncertainty": [
            item.as_record()
            for item in envelope.request.obstacle_geometry_uncertainty
        ],
        "warnings": list(envelope.warnings) + [
            "deterministic field corners are sensitivity cases, not a probability interval or an operational decision",
        ],
        "cases": [
            {
                "selection": dict(case.values),
                "field_screening": field_screening_report(case.result),
            }
            for case in envelope.cases
        ],
    })


def field_joint_measured_history_envelope_report(
    envelope: FieldJointMeasuredHistoryScreeningEnvelope,
) -> dict[str, object]:
    """Serialise source-history and independent field-input corner outcomes."""
    if not isinstance(envelope, FieldJointMeasuredHistoryScreeningEnvelope):
        raise TypeError("envelope must be a FieldJointMeasuredHistoryScreeningEnvelope")
    quality = envelope.quality_assessment
    return _json_native({
        "schema": FIELD_JOINT_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA,
        "quality_assessment": {
            "approved": quality.approved,
            "criteria": asdict(quality.criteria),
            "channel_metrics": [
                {
                    "channel": metric[0],
                    "maximum_sample_interval_s": metric[1],
                    "response_time_s": metric[2],
                    "absolute_time_offset_s": metric[3],
                    "maximum_relative_half_width": metric[4],
                }
                for metric in quality.channel_metrics
            ],
            "reasons": list(quality.reasons),
        },
        "completed_case_count": envelope.completed_case_count,
        "blocked_case_count": envelope.blocked_case_count,
        "property_table_used": envelope.property_table_used,
        "warnings": list(envelope.warnings),
        "cases": [
            {
                "source_selection": dict(case.source_selection),
                "field_selection": dict(case.field_values),
                "direct_vapour_mass_kg": case.schedule.direct_vapour_mass_kg,
                "unrouted_postflash_liquid_mass_kg": (
                    case.schedule.unrouted_postflash_liquid_mass_kg
                ),
                "field_screening": field_screening_report(case.result),
            }
            for case in envelope.cases
        ],
    })


def field_sensor_array_uncertainty_envelope_report(
    envelope: FieldSensorArrayUncertaintyEnvelope,
) -> dict[str, object]:
    """Serialise detector calibration corners sharing one true transport field."""
    if not isinstance(envelope, FieldSensorArrayUncertaintyEnvelope):
        raise TypeError("envelope must be a FieldSensorArrayUncertaintyEnvelope")
    return _json_native({
        "schema": FIELD_SENSOR_ARRAY_UNCERTAINTY_ENVELOPE_REPORT_SCHEMA,
        "base_field_screening": field_screening_report(envelope.screening),
        "completed_case_count": envelope.completed_case_count,
        "warnings": list(envelope.warnings),
        "cases": [
            {
                "selection": dict(case.values),
                "sensor_results": [
                    {
                        "label": sensor_result.label,
                        "position_m": list(sensor_result.sensor.position_m),
                        "maximum_true_mole_fraction": float(np.max(sensor_result.trace.true_mole_fraction)),
                        "maximum_indicated_mole_fraction": float(np.max(sensor_result.trace.indicated_mole_fraction)),
                        "final_indicated_mole_fraction": float(sensor_result.trace.indicated_mole_fraction[-1]),
                    }
                    for sensor_result in case.sensor_results
                    if sensor_result.trace is not None
                ],
            }
            for case in envelope.cases
        ],
    })


def field_refinement_report(result: FieldSemiFVRefinementResult) -> dict[str, object]:
    """Return a screening report extended with the optional refinement audit."""
    if not isinstance(result, FieldSemiFVRefinementResult):
        raise TypeError("result must be a FieldSemiFVRefinementResult")
    report = field_screening_report(result.screening)
    report["refinement_applicability"] = _applicability_record(result.applicability)
    study = result.study
    report["numerical_refinement"] = None if study is None else {
        "relative_tolerance": study.relative_tolerance,
        "converged": study.converged,
        "estimated_cell_steps": study.estimated_cell_steps,
        "warnings": list(study.warnings),
        "cases": [
            {
                "refinement_factor": factor,
                "diagnostics": asdict(case.diagnostics),
            }
            for factor, case in study.cases
        ],
        "receptor_changes": [asdict(change) for change in study.receptor_changes],
    }
    return _json_native(report)


def write_field_refinement_report(
    result: FieldSemiFVRefinementResult,
    path: str | Path,
) -> Path:
    """Write a new screening-plus-refinement report explicitly chosen by the caller.

    Existing evidence artifacts are never overwritten implicitly; exclusive
    creation preserves that guarantee even when writers race.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(field_refinement_report(result), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n"
    )
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(payload)
    return destination


def field_sensor_superposition_report(result: FieldSensorSuperposition) -> dict[str, object]:
    """Return a compact audit record for a common-sensor source sum."""
    if not isinstance(result, FieldSensorSuperposition):
        raise TypeError("result must be a FieldSensorSuperposition")
    trace = result.sensor_trace
    return _json_native({
        "schema": FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA,
        "applicability": _applicability_record(result.applicability),
        "branches": [
            {
                "label": branch.label,
                "status": branch.result.applicability.status,
                "source_location_m": list(branch.result.request.scenario.source.location_m),
                "post_release_duration_s": branch.result.request.post_release_duration_s,
            }
            for branch in result.branches
        ],
        "time_start_s": float(result.time_s[0]),
        "time_end_s": float(result.time_s[-1]),
        "maximum_scalar_concentration_kg_m3": float(np.max(result.concentration_kg_m3)),
        "maximum_true_mole_fraction": float(np.max(trace.true_mole_fraction)),
        "maximum_indicated_mole_fraction": float(np.max(trace.indicated_mole_fraction)),
        "warnings": list(result.warnings),
    })


__all__ = [
    "FIELD_SCREENING_REPORT_SCHEMA", "FIELD_SENSOR_SUPERPOSITION_REPORT_SCHEMA",
    "FIELD_DROPLET_HANDOFF_REFINEMENT_REPORT_SCHEMA",
    "FIELD_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA",
    "FIELD_JOINT_MEASURED_HISTORY_ENVELOPE_REPORT_SCHEMA",
    "FIELD_SENSOR_ARRAY_UNCERTAINTY_ENVELOPE_REPORT_SCHEMA",
    "FIELD_SEMI_FV_ENVELOPE_REPORT_SCHEMA",
    "field_screening_report", "field_droplet_handoff_refinement_report",
    "field_measured_history_envelope_report",
    "field_joint_measured_history_envelope_report",
    "field_semi_fv_envelope_report",
    "field_sensor_array_uncertainty_envelope_report",
    "field_sensor_superposition_report",
    "write_field_screening_report", "field_refinement_report",
    "write_field_refinement_report",
]
