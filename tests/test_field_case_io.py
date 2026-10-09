import json
import hashlib
import math
from types import SimpleNamespace

import pytest

import degali.addons.field_operational_phase_transport as field_operational_phase_transport
import degali.addons.field_execution as field_execution
from degali.cli import main
from degali.addons.field_case_io import (
    FIELD_PHASE_ROUTING_SCREENING_INPUT_SCHEMA,
    FIELD_SCREENING_INPUT_SCHEMA,
    field_screening_case_from_mapping,
    field_screening_request_from_mapping,
    read_field_screening_case_json,
    read_field_screening_request_json,
)
from degali.addons.field_decision import FieldOperationalScreeningDecision
from degali.addons.field_historian_io import ImportedPressureDrivenMeasuredHistory
from degali.addons.field_workflow import run_field_semi_fv_screening
from degali.addons.site_geometry import OrientedCuboid


def _bounded(value, unit, source, *, relative=0.0):
    return {
        "nominal": value,
        "lower": value * (1.0 - relative),
        "upper": value * (1.0 + relative),
        "unit": unit,
        "source": source,
    }


def _sensor(position, source):
    return {
        "position_m": position,
        "response_time_s": _bounded(0.25, "s", source + "-t90"),
        "gain": _bounded(1.0, "1", source + "-gain"),
        "bias_mole_fraction": _bounded(0.0, "mole_fraction", source + "-bias"),
        "averaging_time_s": 0.1,
    }


def _case():
    return {
        "schema": FIELD_SCREENING_INPUT_SCHEMA,
        "lh2_validation_available": False,
        "coordinate_reference": {
            "coordinate_system_id": "plant-grid-rev-A",
            "origin_id": "transfer-skid-origin",
            "x_axis_bearing_math_to_deg": 0.0,
            "vertical_datum_id": "site-grade",
            "evidence_id": "layout-drawing-A",
        },
        "conditional_review": {
            "review_id": "field-review-A",
            "reviewer_id": "process-safety-reviewer-01",
            "reviewer_role": "process_safety_engineer",
            "reviewed_at_utc": "2026-10-05T00:00:00Z",
            "evidence_id": "field-screening-review-record-A",
            "allowed_scope": "conditional_screening_only",
        },
        "scenario": {
            "model_family": "degali",
            "temporal_mode": "transient",
            "source": {
                "fluid": "lh2",
                "location_m": [0.0, 0.0, 0.5],
                "direction_m": [1.0, 0.0, 0.0],
                "upstream_pressure": _bounded(0.4e6, "Pa", "PT-source"),
                "upstream_temperature": _bounded(26.084, "K", "TT-source"),
                "pressure_reference": "absolute",
                "mass_flow_kg_s": _bounded(0.265, "kg/s", "FT-source", relative=0.05),
                "opening_area_m2": _bounded(math.pi * 0.012**2 / 4.0 / 0.8, "m2", "drawing-A"),
                "discharge_coefficient": _bounded(0.8, "1", "Cd-test"),
                "liquid_fraction": _bounded(0.922, "1", "phase-review-A"),
                "liquid_fraction_basis": "upstream",
                "flash_model": "homogeneous_equilibrium",
                "duration_s": 1.0,
                "metadata": {"release_boundary_evidence_id": "scenario-review-A"},
            },
            "weather": {
                "speed_m_s": _bounded(2.0, "m/s", "met-A", relative=0.1),
                "direction_deg": _bounded(270.0, "deg", "met-A"),
                "direction_convention": "meteorological_from",
                "stability": "neutral",
                "reference_height_m": 10.0,
            },
            "surface": {
                "heat_transfer_w_m2_k": _bounded(10.0, "W/m2/K", "surface-A"),
                "surface_temperature_k": _bounded(293.15, "K", "surface-A"),
                "substrate": "concrete",
                "evidence_id": "surface-review-A",
            },
            "sensor": _sensor([1.0, 0.0, 0.5], "H2-01"),
        },
        "transport": {
            "length_m": 4.0,
            "height_m": 2.0,
            "nx": 20,
            "nz": 10,
            "time_step_s": 0.005,
            "duration_s": 1.0,
            "diffusivity_m2_s": 0.5,
            "gravitational_settling_m_s": 0.0,
            "source_sigma_m": 0.2,
        },
        "obstacles": [{
            "x_min_m": 2.0,
            "x_max_m": 2.5,
            "y_min_m": -0.25,
            "y_max_m": 0.25,
            "z_min_m": 0.0,
            "z_max_m": 0.6,
            "label": "valve-skid",
        }],
        "sensor_deployments": [{
            "label": "H2-02",
            "sensor": _sensor([1.5, 0.0, 0.5], "H2-02"),
        }],
    }


def _measured_history_case(tmp_path):
    historian = tmp_path / "vent-event.csv"
    historian.write_text(
        "time_s,pressure_pa,temperature_k,mass_flow_kg_s,liquid_fraction\n"
        "0,400000,26.084,0.265,0.922\n"
        "1,401000,26.200,0,0.922\n",
        encoding="utf-8",
    )
    case = _case()
    case["measured_history"] = {
        "csv_path": historian.name,
        "map": {
            "time_s_column": "time_s",
            "event_id": "vent-event-A",
            "event_evidence_id": "event-window-review-A",
            "phase_evidence_id": "phase-state-review-A",
            "pressure_pa": {
                "value_column": "pressure_pa", "unit": "Pa", "source_id": "PT-01",
                "calibration_evidence_id": "PT-01-cal", "absolute_half_width": 500.0,
            },
            "temperature_k": {
                "value_column": "temperature_k", "unit": "K", "source_id": "TT-01",
                "calibration_evidence_id": "TT-01-cal", "absolute_half_width": 0.25,
            },
            "mass_flow_kg_s": {
                "value_column": "mass_flow_kg_s", "unit": "kg/s", "source_id": "FT-01",
                "calibration_evidence_id": "FT-01-cal", "relative_half_width": 0.05,
            },
            "liquid_fraction": {
                "value_column": "liquid_fraction", "unit": "1", "source_id": "LT-01",
                "calibration_evidence_id": "LT-01-basis", "absolute_half_width": 0.0,
            },
        },
        "quality_criteria": {
            "maximum_sample_interval_s": 1.0,
            "maximum_response_time_s": 1.0,
            "maximum_absolute_time_offset_s": 0.1,
            "maximum_relative_half_width": 0.1,
            "evidence_id": "historian-quality-procedure-A",
        },
    }
    return case


def _pressure_driven_history_case(tmp_path):
    historian = tmp_path / "pressure-driven-event.csv"
    historian.write_text(
        "time_s,pressure_pa,temperature_k\n"
        "0,400000,26.076\n"
        "1,390000,25.949\n",
        encoding="utf-8",
    )
    case = _case()
    source = case["scenario"]["source"]
    source["mass_flow_kg_s"] = _bounded(0.0, "kg/s", "pressure-history-placeholder")
    case["pressure_driven_history"] = {
        "csv_path": historian.name,
        "map": {
            "time_s_column": "time_s",
            "event_id": "pressure-event-A",
            "event_evidence_id": "pressure-event-window-review-A",
            "phase_evidence_id": "pressure-phase-review-A",
            "pressure_pa": {
                "value_column": "pressure_pa", "unit": "Pa", "source_id": "PT-ORIFICE-01",
                "calibration_evidence_id": "PT-ORIFICE-01-cal", "absolute_half_width": 500.0,
            },
            "temperature_k": {
                "value_column": "temperature_k", "unit": "K", "source_id": "TT-ORIFICE-01",
                "calibration_evidence_id": "TT-ORIFICE-01-cal", "absolute_half_width": 0.05,
            },
        },
        "quality_criteria": {
            "maximum_sample_interval_s": 1.0,
            "maximum_response_time_s": 1.0,
            "maximum_absolute_time_offset_s": 0.1,
            "maximum_relative_half_width": 0.1,
            "evidence_id": "pressure-history-quality-A",
        },
    }
    return case


def _phase_routing_case():
    case = _case()
    case["schema"] = FIELD_PHASE_ROUTING_SCREENING_INPUT_SCHEMA
    case["transport"]["duration_s"] = 2.0
    case["phase_routing_transport"] = {
        "phase_routing": {
            "post_release_duration_s": 1.0,
            "puff_duration_s": 1.0,
            "pool_area_m2": 0.5,
            "pool_time_step_s": 0.1,
            "evaporation_coefficient_m2_s": 1.0e-7,
            "pool_model": "dynamic",
            "evidence_id": "phase-routing-assumptions-A",
            "gas_numerics": {
                "maximum_nearfield_distance_m": 0.2,
                "radial_points": 41,
                "nearfield_maximum_step_m": 0.001,
                "nearfield_relative_tolerance": 2.0e-6,
                "crosswind_maximum_distance_m": 10.0,
                "crosswind_maximum_step_m": 0.05,
                "puff_time_step_s": 0.02,
            },
            "liquid_numerics": {"maximum_droplet_time_s": 1.0},
            "droplet_population": {
                "classes": [
                    {"diameter_m": 1.0e-4, "mass_fraction": 0.25},
                    {"diameter_m": 1.0e-3, "mass_fraction": 0.75},
                ],
                "evidence_id": "spray-population-review-A",
            },
        },
        "pool_launch": {
            "source_height_m": 0.0,
            "closure_id": "ground_level_scalar",
            "evidence_id": "pool-launch-review-A",
        },
        "pool_vertical_sigma_m": 0.2,
        "table_nodes": 81,
    }
    return case


def test_json_case_contract_builds_a_reproducible_local_obstacle_screen(tmp_path):
    source = tmp_path / "field-case.json"
    source.write_text(json.dumps(_case()), encoding="utf-8")

    request = read_field_screening_request_json(source)
    assert request.scenario.source.pressure_reference == "absolute"
    assert request.coordinate_reference is not None
    assert request.coordinate_reference.coordinate_system_id == "plant-grid-rev-A"
    assert request.scenario.sensor is not None
    assert request.obstacles[0].label == "valve-skid"
    assert request.sensor_deployments[0].label == "H2-02"

    parsed = read_field_screening_case_json(source)
    assert parsed.conditional_review is not None
    assert parsed.conditional_review.review_id == "field-review-A"

    result = run_field_semi_fv_screening(request)
    assert result.completed
    assert result.transport is not None
    assert result.transport.diagnostics.obstacle_count == 1
    assert {item.label for item in result.sensor_results} == {"field_sensor", "H2-02"}


def test_json_case_supports_explicit_pressure_driven_mass_flow_adapter():
    case = _case()
    source = case["scenario"]["source"]
    source.pop("mass_flow_kg_s")
    source["pressure_driven_mass_flow"] = {
        "ambient_pressure_pa": _bounded(
            101325.0, "Pa", "ambient-pressure-review-A", relative=0.01,
        ),
        "source_id": "orifice-json-review-A",
    }

    request = field_screening_request_from_mapping(case)
    derived = request.scenario.source.mass_flow_kg_s
    assert derived.unit == "kg/s"
    assert derived.source == "orifice-json-review-A"
    # Ambient-pressure uncertainty need not widen a choked-flow result.  The
    # deterministic envelope must contain the nominal value, but equality is
    # physically valid when every declared corner has the same mass flux.
    assert 0.0 < derived.lower <= derived.nominal <= derived.upper
    assert request.scenario.source.metadata["mass_flow_derivation"] == (
        "pressure_driven_homogeneous_equilibrium_throat"
    )
    assert request.scenario.source.metadata["mass_flow_derivation_source_id"] == (
        "orifice-json-review-A"
    )
    assert request.scenario.source.metadata[
        "mass_flow_derivation_ambient_pressure_bounds_pa"
    ]
    source_corners = request.scenario.corner_cases(max_cases=32)
    assert all("mass_flow_kg_s" not in values for values in source_corners)
    ambient_corners = {
        values["pressure_driven_ambient_pressure_pa"]
        for values in source_corners
    }
    assert sorted(ambient_corners) == pytest.approx(
        sorted((
            request.scenario.source.pressure_driven_mass_flow.ambient_pressure_pa.lower,
            request.scenario.source.pressure_driven_mass_flow.ambient_pressure_pa.upper,
        ))
    )
    derived_corner_rates = {
        request.scenario.at_corner(values).source.mass_flow_kg_s.nominal
        for values in source_corners
    }
    assert sorted(derived_corner_rates) == pytest.approx(
        sorted((derived.lower, derived.upper))
    )
    result = run_field_semi_fv_screening(request)
    assert result.completed

    malformed = _case()
    malformed["scenario"]["source"]["pressure_driven_mass_flow"] = {
        "ambient_pressure_pa": _bounded(
            101325.0, "Pa", "ambient-pressure-review-A",
        ),
        "source_id": "orifice-json-review-A",
    }
    with pytest.raises(ValueError, match="exactly one of mass_flow_kg_s"):
        field_screening_request_from_mapping(malformed)

    missing = _case()
    missing["scenario"]["source"].pop("mass_flow_kg_s")
    with pytest.raises(ValueError, match="exactly one of mass_flow_kg_s"):
        field_screening_request_from_mapping(missing)

    unsupported_corner = _case()
    unsupported_corner["scenario"]["source"].pop("mass_flow_kg_s")
    unsupported_corner["scenario"]["source"]["pressure_driven_mass_flow"] = {
        "ambient_pressure_pa": {
            "nominal": 101325.0,
            "lower": 101325.0,
            "upper": 500000.0,
            "unit": "Pa",
            "source": "ambient-pressure-review-A",
        },
        "source_id": "orifice-json-review-A",
    }
    with pytest.raises(ValueError, match="corner cannot be evaluated"):
        field_screening_request_from_mapping(unsupported_corner)


def test_json_case_supports_pressure_driven_history_without_fabricating_measured_flow(tmp_path):
    source = tmp_path / "pressure-history-case.json"
    source.write_text(json.dumps(_pressure_driven_history_case(tmp_path)), encoding="utf-8")

    parsed = read_field_screening_case_json(source)

    assert isinstance(parsed.imported_history, ImportedPressureDrivenMeasuredHistory)
    assert not hasattr(parsed.imported_history.history, "mass_flow_kg_s")
    assert parsed.request.scenario.source.mass_flow_kg_s.nominal > 0.0
    assert dict(parsed.request.measured_history_provenance)["history_kind"] == (
        "pressure_driven_orifice"
    )
    result = run_field_semi_fv_screening(parsed.request)
    assert result.completed
    assert result.transport is not None


def test_json_case_accepts_a_provenanced_distributed_source_with_shape_bounds():
    case = _case()
    case["obstacles"] = []
    case["sensor_deployments"] = []
    case["distributed_vapour_sources"] = [{
        "label": "pool-vapour-json-A",
        "source_kind": "pool_vapour",
        "position_m": [
            _bounded(1.0, "m", "pool-layout-A", relative=0.1),
            _bounded(0.0, "m", "pool-layout-A"),
            _bounded(0.1, "m", "pool-layout-A"),
        ],
        "vertical_sigma_m": _bounded(0.15, "m", "pool-width-A", relative=1.0 / 3.0),
        "schedule": {
            "time_s": [0.0, 0.5, 1.0],
            "rate_kg_s": [0.02, 0.01, 0.0],
            "rate_lower_kg_s": [0.01, 0.005, 0.0],
            "rate_upper_kg_s": [0.03, 0.015, 0.0],
            "source_id": "pool-ledger-json-A",
            "rate_operator": "linear",
        },
        "evidence_id": "pool-evidence-json-A",
    }]

    parsed = field_screening_case_from_mapping(case)
    source = parsed.request.distributed_vapour_sources[0]

    assert source.source_kind == "pool_vapour"
    assert source.position_uncertainty_m is not None
    assert source.vertical_sigma_uncertainty is not None
    assert source.uncertainty_fields()[
        "distributed_source.pool-vapour-json-A.position_x_m"
    ].corners() == pytest.approx((0.9, 1.1))
    assert source.has_schedule_uncertainty
    assert source.schedule.rate_operator == "linear"
    assert source.schedule_corners()[0][1].rate_kg_s == (0.01, 0.005, 0.0)
    assert run_field_semi_fv_screening(parsed.request).completed


def test_json_case_rejects_a_distributed_source_without_a_zero_endpoint():
    case = _case()
    case["distributed_vapour_sources"] = [{
        "label": "bad-pool-json-A",
        "position_m": [1.0, 0.0, 0.1],
        "vertical_sigma_m": _bounded(0.15, "m", "pool-width-A"),
        "schedule": {
            "time_s": [0.0, 0.5, 1.0],
            "rate_kg_s": [0.02, 0.01, 0.01],
            "source_id": "pool-ledger-json-A",
        },
        "evidence_id": "pool-evidence-json-A",
    }]

    with pytest.raises(ValueError, match="zero final endpoint"):
        field_screening_case_from_mapping(case)


def test_json_case_accepts_bounded_release_duration_and_propagates_it_to_corners():
    case = _case()
    case["scenario"]["source"]["duration_s"] = _bounded(
        1.0, "s", "release-timing-review-A", relative=0.25,
    )
    case["transport"]["duration_s"] = 1.25
    parsed = field_screening_case_from_mapping(case)

    assert parsed.request.scenario.source.duration_s == pytest.approx(1.0)
    assert parsed.request.scenario.source.duration_uncertainty is not None
    assert {corner["duration_s"] for corner in parsed.request.scenario.corner_cases()} == {
        0.75, 1.25,
    }


def test_json_case_accepts_bounded_ambient_boundary_values():
    case = _case()
    case.update({
        "ambient_temperature_uncertainty_k": _bounded(
            295.0, "K", "ambient-mast-A", relative=0.02,
        ),
        "ambient_pressure_uncertainty_pa": _bounded(
            101325.0, "Pa", "ambient-mast-A", relative=0.01,
        ),
        "ambient_air_density_uncertainty_kg_m3": _bounded(
            1.2, "kg/m3", "ambient-mast-A", relative=0.05,
        ),
    })

    parsed = field_screening_case_from_mapping(case)

    assert parsed.request.ambient_temperature_k == pytest.approx(295.0)
    assert parsed.request.ambient_temperature_uncertainty_k is not None
    assert parsed.request.ambient_pressure_uncertainty_pa is not None
    assert parsed.request.ambient_air_density_uncertainty_kg_m3 is not None
    assert parsed.request.ambient_uncertainty_fields()[
        "ambient_temperature_k"
    ].corners() == pytest.approx((289.1, 300.9))


def test_json_case_rejects_ambient_bound_nominal_mismatch():
    case = _case()
    case.update({
        "ambient_temperature_k": 294.0,
        "ambient_temperature_uncertainty_k": _bounded(
            295.0, "K", "ambient-mast-A", relative=0.02,
        ),
    })

    with pytest.raises(ValueError, match="nominal must match ambient_temperature_k"):
        field_screening_case_from_mapping(case)


def test_json_case_accepts_numeric_duration_with_separate_bound():
    case = _case()
    case["scenario"]["source"]["duration_uncertainty"] = _bounded(
        1.0, "s", "release-timing-review-A", relative=0.25,
    )
    case["transport"]["duration_s"] = 1.25
    parsed = field_screening_case_from_mapping(case)

    assert parsed.request.scenario.source.duration_s == pytest.approx(1.0)
    assert parsed.request.scenario.source.duration_uncertainty is not None


def test_json_case_rejects_transport_domain_shorter_than_duration_upper_corner():
    case = _case()
    case["scenario"]["source"]["duration_uncertainty"] = _bounded(
        1.0, "s", "release-timing-review-A", relative=0.25,
    )
    with pytest.raises(ValueError, match="maximum declared source duration_s"):
        field_screening_case_from_mapping(case)


def test_json_case_accepts_bounded_source_location_coordinates():
    case = _case()
    case["scenario"]["source"]["location_m"] = [
        _bounded(0.0, "m", "layout-x-review", relative=0.1),
        _bounded(0.0, "m", "layout-y-review", relative=0.1),
        _bounded(0.5, "m", "layout-z-review", relative=0.1),
    ]
    parsed = field_screening_case_from_mapping(case)

    assert parsed.request.scenario.source.location_m == pytest.approx((0.0, 0.0, 0.5))
    assert parsed.request.scenario.source.location_uncertainty_m is not None
    corners = parsed.request.scenario.corner_cases(max_cases=8)
    assert len(corners) == 8
    assert {corner["source_location_z_m"] for corner in corners} == {0.45, 0.55}


def test_json_case_rejects_mixed_numeric_and_bounded_source_location_coordinates():
    case = _case()
    case["scenario"]["source"]["location_m"] = [
        _bounded(0.0, "m", "layout-x-review"), 0.0, 0.5,
    ]
    with pytest.raises(ValueError, match="either three numbers or three bounded objects"):
        field_screening_case_from_mapping(case)


def test_json_case_accepts_bounded_source_direction_components():
    case = _case()
    case["scenario"]["source"]["direction_m"] = [
        _bounded(1.0, "1", "orientation-review-A", relative=0.1),
        _bounded(0.0, "1", "orientation-review-A"),
        _bounded(0.0, "1", "orientation-review-A"),
    ]
    parsed = field_screening_case_from_mapping(case)

    assert parsed.request.scenario.source.direction_m == pytest.approx((1.0, 0.0, 0.0))
    assert parsed.request.scenario.source.direction_uncertainty_m is not None
    corners = parsed.request.scenario.corner_cases(max_cases=8)
    assert {corner["source_direction_x"] for corner in corners} == {0.9, 1.1}


def test_json_case_rejects_zero_vector_source_direction_corner():
    case = _case()
    case["scenario"]["source"]["direction_m"] = [
        _bounded(1.0, "1", "orientation-review-A", relative=1.0),
        _bounded(0.0, "1", "orientation-review-A"),
        _bounded(0.0, "1", "orientation-review-A"),
    ]
    with pytest.raises(ValueError, match="zero-vector corner"):
        field_screening_case_from_mapping(case)


def test_json_case_rejects_unbacked_validation_flag():
    case = _case()
    case["lh2_validation_available"] = True
    with pytest.raises(ValueError, match="requires validation_evidence"):
        field_screening_case_from_mapping(case)


def test_json_case_preserves_fingerprinted_validation_evidence():
    case = _case()
    case["lh2_validation_available"] = True
    case["validation_evidence"] = {
        "dataset_id": "lh2-obstacle-trial-A",
        "path": "evidence/trial-A.csv",
        "sha256": "0" * 64,
        "row_count": 12,
        "source_boundary_id": "source-A",
        "weather_id": "met-A",
        "obstacle_geometry_id": "obstacle-A",
        "receptor_geometry_id": "detectors-A",
        "temporal_operator_id": "t90-average-A",
        "common_clock_id": "clock-A",
        "scope": "lh2_obstacle_transport",
    }
    parsed = field_screening_case_from_mapping(case)

    assert parsed.request.validation_evidence is not None
    assert parsed.request.validation_evidence.dataset_id == "lh2-obstacle-trial-A"
    assert parsed.request.validation_evidence.supports_obstacle_transport


def test_json_case_rejects_duplicate_duration_uncertainty_declarations():
    case = _case()
    case["scenario"]["source"]["duration_s"] = _bounded(
        1.0, "s", "release-timing-review-A", relative=0.25,
    )
    case["scenario"]["source"]["duration_uncertainty"] = _bounded(
        1.0, "s", "release-timing-review-A", relative=0.25,
    )
    with pytest.raises(ValueError, match="bounded object and duration_uncertainty together"):
        field_screening_case_from_mapping(case)


def test_json_case_reader_rejects_normalized_duplicate_keys(tmp_path):
    source = tmp_path / "duplicate-keys.json"
    source.write_text(
        '{"schema": "degali.field-screening-input.v1", '
        '" SCHEMA ": "degali.field-screening-input.v1"}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate or normalization-colliding"):
        read_field_screening_case_json(source)


def test_json_case_contract_refuses_typos_and_unit_substitution():
    malformed = _case()
    malformed["scenario"]["source"]["mass_flow_kg_s"]["unit"] = "Nm3/h"
    with pytest.raises(ValueError, match="mass_flow_kg_s.unit"):
        field_screening_request_from_mapping(malformed)

    malformed = _case()
    malformed["transport"]["implicit_wake_factor"] = 0.5
    with pytest.raises(ValueError, match="unknown=implicit_wake_factor"):
        field_screening_request_from_mapping(malformed)

    malformed = _case()
    malformed["transport"]["duration_s"] = 2.0
    with pytest.raises(ValueError, match="transport.duration_s must equal"):
        field_screening_request_from_mapping(malformed)

    malformed = _case()
    del malformed["coordinate_reference"]
    with pytest.raises(ValueError, match="missing=coordinate_reference"):
        field_screening_request_from_mapping(malformed)

    malformed = _case()
    malformed["conditional_review"]["allowed_scope"] = "design_basis"
    with pytest.raises(ValueError, match="conditional_screening_only"):
        field_screening_request_from_mapping(malformed)

    malformed = _case()
    malformed["conditional_review"]["reviewed_at_utc"] = "2026-10-05"
    with pytest.raises(ValueError, match="ISO-8601 UTC timestamp"):
        field_screening_request_from_mapping(malformed)


def test_json_case_accepts_a_wind_direction_interval_that_crosses_north():
    case = _case()
    case["scenario"]["weather"]["direction_deg"] = {
        "nominal": 0.0, "lower": 350.0, "upper": 10.0,
        "unit": "deg", "source": "met-A",
    }

    request = field_screening_request_from_mapping(case)

    direction = request.scenario.weather.direction_deg
    assert direction.corners() == pytest.approx((350.0, 10.0))
    assert {corner["wind_direction_deg"] for corner in request.scenario.corner_cases()} == {350.0, 10.0}


def test_json_case_accepts_an_oriented_cuboid_with_an_explicit_global_bearing():
    case = _case()
    case["obstacles"] = [{
        "center_m": [2.25, 0.0], "length_m": 1.0, "width_m": 0.5,
        "z_min_m": 0.0, "z_max_m": 0.6, "long_axis_bearing_deg": 30.0,
        "long_axis_bearing_convention": "math_to", "label": "rotated-valve-skid",
    }]

    request = field_screening_request_from_mapping(case)

    assert isinstance(request.obstacles[0], OrientedCuboid)
    assert request.obstacles[0].long_axis_bearing_deg == pytest.approx(30.0)

    case["obstacles"][0]["long_axis_bearing_convention"] = "meteorological_from"
    with pytest.raises(ValueError, match="long_axis_bearing_convention"):
        field_screening_request_from_mapping(case)


def test_json_case_accepts_bounded_axis_aligned_obstacle_geometry():
    case = _case()
    def bound(value, lower=None, upper=None):
        return {
            "nominal": value,
            "lower": value if lower is None else lower,
            "upper": value if upper is None else upper,
            "unit": "m",
            "source": "layout-revision-A",
        }

    case["obstacles"] = [{
        "x_min_m": bound(2.0, 1.8, 2.2),
        "x_max_m": bound(2.5),
        "y_min_m": bound(-0.25),
        "y_max_m": bound(0.25),
        "z_min_m": bound(0.0),
        "z_max_m": bound(0.6),
        "label": "valve-skid",
        "uncertainty_evidence_id": "layout-obstacle-review-A",
    }]

    request = field_screening_request_from_mapping(case)

    assert request.obstacles[0].label == "valve-skid"
    assert len(request.obstacle_geometry_uncertainty) == 1
    uncertainty = request.obstacle_geometry_uncertainty[0]
    assert uncertainty.evidence_id == "layout-obstacle-review-A"
    assert uncertainty.uncertainty_fields()["obstacle.valve-skid.x_min_m"].lower == pytest.approx(1.8)


def test_json_case_accepts_bounded_oriented_obstacle_geometry():
    case = _case()
    linear = lambda value: {
        "nominal": value, "lower": value, "upper": value,
        "unit": "m", "source": "layout-revision-B",
    }
    case["obstacles"] = [{
        "center_m": [linear(2.25), linear(0.0)],
        "length_m": linear(1.0), "width_m": linear(0.5),
        "z_min_m": linear(0.0), "z_max_m": linear(0.6),
        "long_axis_bearing_deg": {
            "nominal": 30.0, "lower": 25.0, "upper": 35.0,
            "unit": "deg", "source": "layout-revision-B",
        },
        "long_axis_bearing_convention": "math_to",
        "label": "rotated-valve-skid",
        "uncertainty_evidence_id": "layout-obstacle-review-B",
    }]

    request = field_screening_request_from_mapping(case)

    assert isinstance(request.obstacles[0], OrientedCuboid)
    assert len(request.obstacle_geometry_uncertainty) == 1
    assert request.obstacle_geometry_uncertainty[0].uncertainty_fields()[
        "obstacle.rotated-valve-skid.long_axis_bearing_deg"
    ].corners() == pytest.approx((25.0, 35.0))


def test_json_case_contract_requires_stability_evidence_and_per_class_mixing_closure():
    case = _case()
    case["stability_alternatives"] = {
        "alternatives": ["stable"], "evidence_id": "met-stability-review-A",
    }
    with pytest.raises(ValueError, match="require an explicit stability_mixing_closure"):
        field_screening_request_from_mapping(case)

    case["stability_mixing_closure"] = {
        "diffusivity_m2_s": {"neutral": 0.5, "stable": 0.2},
        "evidence_id": "met-mixing-A",
    }
    request = field_screening_request_from_mapping(case)
    assert request.stability_alternatives is not None
    assert request.stability_alternatives.classes_for("neutral") == ("neutral", "stable")


def test_json_case_contract_does_not_accept_a_postflash_schedule_without_its_dedicated_gate():
    malformed = _case()
    malformed["direct_vapour_schedule"] = {
        "time_s": [0.0, 1.0], "rate_kg_s": [0.1, 0.0],
    }
    with pytest.raises(ValueError, match="unknown=direct_vapour_schedule"):
        field_screening_request_from_mapping(malformed)


def test_phase_routing_schema_builds_only_the_dedicated_conservative_pool_handoff():
    parsed = field_screening_case_from_mapping(_phase_routing_case())

    phase = parsed.phase_routing_transport
    assert phase is not None
    assert phase.config.post_release_duration_s == pytest.approx(1.0)
    assert phase.config.pool_model == "dynamic"
    assert phase.config.gas_model_options["crosswind_maximum_distance"] == pytest.approx(10.0)
    assert phase.config.phase_model_options["maximum_droplet_time_s"] == pytest.approx(1.0)
    assert len(phase.droplet_population) == 2
    assert phase.droplet_population_evidence_id == "spray-population-review-A"
    assert phase.phase_routing_evidence_id == "phase-routing-assumptions-A"
    assert phase.pool_launch.evidence_id == "pool-launch-review-A"
    assert phase.pool_vertical_sigma_m == pytest.approx(0.2)
    assert phase.table_nodes == 81
    assert parsed.request.post_release_duration_s == pytest.approx(1.0)

    legacy = _case()
    legacy["phase_routing_transport"] = _phase_routing_case()["phase_routing_transport"]
    with pytest.raises(ValueError, match="unknown=phase_routing_transport"):
        field_screening_case_from_mapping(legacy)


def test_phase_routing_schema_retains_explicit_phase_scalar_uncertainty():
    case = _phase_routing_case()
    case["transport"]["duration_s"] = 2.0
    case["phase_routing_transport"]["phase_routing"]["uncertainty"] = {
        "pool_area_m2": {
            "nominal": 0.5, "lower": 0.4, "upper": 0.6,
            "unit": "m2", "source": "pool-footprint-review-A",
        },
        "evaporation_coefficient_m2_s": {
            "nominal": 1.0e-7, "lower": 5.0e-8, "upper": 2.0e-7,
            "unit": "m2/s", "source": "evaporation-review-A",
        },
    }

    parsed = field_screening_case_from_mapping(case)
    uncertainty = parsed.phase_routing_transport.phase_uncertainty
    assert uncertainty is not None
    assert uncertainty.pool_area_m2 is not None
    assert uncertainty.pool_area_m2.corners() == pytest.approx((0.4, 0.6))
    assert uncertainty.evaporation_coefficient_m2_s is not None
    assert uncertainty.evaporation_coefficient_m2_s.corners() == pytest.approx(
        (5.0e-8, 2.0e-7)
    )

    case["phase_routing_transport"]["phase_routing"]["uncertainty"]["pool_area_m2"]["unit"] = "m"
    with pytest.raises(ValueError, match="pool_area_m2.unit"):
        field_screening_case_from_mapping(case)


def test_phase_routing_schema_retains_pool_launch_width_uncertainty():
    case = _phase_routing_case()
    case["phase_routing_transport"]["pool_vertical_sigma_uncertainty"] = {
        "nominal": 0.2, "lower": 0.1, "upper": 0.3,
        "unit": "m", "source": "pool-width-review-A",
    }
    parsed = field_screening_case_from_mapping(case)
    phase = parsed.phase_routing_transport
    assert phase is not None
    assert phase.pool_vertical_sigma_m == pytest.approx(0.2)
    assert phase.pool_vertical_sigma_uncertainty is not None
    assert phase.pool_vertical_sigma_uncertainty.corners() == pytest.approx((0.1, 0.3))

    case["phase_routing_transport"]["pool_vertical_sigma_m"] = {
        "nominal": 0.2, "lower": 0.1, "upper": 0.3,
        "unit": "m", "source": "pool-width-review-B",
    }
    del case["phase_routing_transport"]["pool_vertical_sigma_uncertainty"]
    parsed_inline = field_screening_case_from_mapping(case)
    assert parsed_inline.phase_routing_transport is not None
    assert parsed_inline.phase_routing_transport.pool_vertical_sigma_uncertainty is not None

    case["phase_routing_transport"]["pool_vertical_sigma_m"]["unit"] = "m2"
    with pytest.raises(ValueError, match="pool_vertical_sigma_m.unit"):
        field_screening_case_from_mapping(case)


def test_phase_routing_schema_retains_complete_droplet_population_corners():
    case = _phase_routing_case()
    case["phase_routing_transport"]["phase_routing"]["uncertainty"] = {
        "droplet_population": {
            "evidence_id": "spray-population-envelope-A",
            "corners": [
                {
                    "classes": [
                        {"diameter_m": 1.0e-4, "mass_fraction": 0.25},
                        {"diameter_m": 1.0e-3, "mass_fraction": 0.75},
                    ],
                },
                {
                    "classes": [
                        {"diameter_m": 8.0e-5, "mass_fraction": 0.5},
                        {"diameter_m": 1.2e-3, "mass_fraction": 0.5},
                    ],
                },
            ],
        },
    }

    parsed = field_screening_case_from_mapping(case)
    uncertainty = parsed.phase_routing_transport.phase_uncertainty
    assert uncertainty is not None
    assert uncertainty.droplet_population_corners is not None
    assert len(uncertainty.droplet_population_corners) == 2
    assert uncertainty.droplet_population_uncertainty_evidence_id == (
        "spray-population-envelope-A"
    )

    case["phase_routing_transport"]["phase_routing"]["uncertainty"][
        "droplet_population"
    ]["corners"][1]["classes"][1]["mass_fraction"] = 0.6
    with pytest.raises(ValueError, match="sum to one"):
        field_screening_case_from_mapping(case)


def test_phase_routing_schema_rejects_ambiguous_or_uncontrolled_source_paths():
    case = _phase_routing_case()
    case["post_release_duration_s"] = 2.0
    with pytest.raises(ValueError, match="must match"):
        field_screening_case_from_mapping(case)

    case = _phase_routing_case()
    case["measured_history"] = {}
    with pytest.raises(ValueError, match="unknown=measured_history"):
        field_screening_case_from_mapping(case)

    case = _phase_routing_case()
    case["phase_routing_transport"]["phase_routing"]["wake_factor"] = 0.5
    with pytest.raises(ValueError, match="unknown=wake_factor"):
        field_screening_case_from_mapping(case)

    case = _phase_routing_case()
    droplet_population = (
        case["phase_routing_transport"]["phase_routing"]["droplet_population"]
    )
    droplet_population["classes"][1]["mass_fraction"] = 0.5
    with pytest.raises(ValueError, match="mass fractions must sum to one"):
        field_screening_case_from_mapping(case)


def test_json_case_can_attach_only_a_quality_gated_and_provenanced_measured_history(tmp_path):
    source = tmp_path / "measured-field-case.json"
    source.write_text(json.dumps(_measured_history_case(tmp_path)), encoding="utf-8")

    case = read_field_screening_case_json(source)
    assert case.imported_history is not None
    assert case.measured_schedule is not None
    assert case.request.direct_vapour_schedule is not None
    assert case.request.measured_history_quality is not None
    assert case.request.measured_history_quality.approved
    assert not case.request.measured_history_source_uncertainty_resolved
    assert dict(case.request.measured_history_provenance)["event_id"] == "vent-event-A"

    with pytest.raises(ValueError, match="requires a case file path"):
        field_screening_request_from_mapping(_measured_history_case(tmp_path))


def test_field_screen_cli_emits_an_auditable_withheld_record_for_diagnostic_only_run(
    tmp_path, capsys,
):
    source = tmp_path / "field-case.json"
    source.write_text(json.dumps(_case()), encoding="utf-8")

    code = main(["field-screen", str(source), "--no-refinement", "--require-screening"])
    payload = json.loads(capsys.readouterr().out)

    assert code == 2
    assert payload["schema"] == "degali.field-screening-execution.v1"
    assert len(payload["input"]["sha256"]) == 64
    assert payload["input"]["refinement_requested"] is False
    assert payload["operational_screening"]["status"] == "withheld"
    assert payload["operational_screening"]["design_basis_allowed"] is False


def test_field_screen_cli_writes_a_new_conditional_screening_record(tmp_path, capsys):
    source = tmp_path / "field-case.json"
    output = tmp_path / "field-execution.json"
    exact_case = _case()
    exact_case["scenario"]["source"]["mass_flow_kg_s"]["lower"] = 0.265
    exact_case["scenario"]["source"]["mass_flow_kg_s"]["upper"] = 0.265
    exact_case["scenario"]["weather"]["speed_m_s"]["lower"] = 2.0
    exact_case["scenario"]["weather"]["speed_m_s"]["upper"] = 2.0
    source.write_text(json.dumps(exact_case), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--output", str(output),
        "--allow-conditional", "--relative-tolerance", "1.0",
    ])
    captured = capsys.readouterr()
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert code == 0
    assert "conditional_allowed" in captured.out
    assert payload["input"]["refinement_requested"] is True
    assert payload["input"]["conditional_review"]["review_id"] == "field-review-A"
    assert payload["input"]["conditional_review"]["applied_to_this_execution"] is True
    assert payload["operational_screening"]["status"] == "conditional_allowed"
    assert payload["operational_screening"]["approval_allowed"] is False

    assert main(["field-screen", str(source), "--output", str(output)]) == 1
    assert "refusing to overwrite" in capsys.readouterr().err


def test_field_screen_cli_refuses_conditional_opt_in_without_a_review_record(
    tmp_path, capsys,
):
    source = tmp_path / "unreviewed-field-case.json"
    case = _case()
    del case["conditional_review"]
    source.write_text(json.dumps(case), encoding="utf-8")

    assert main([
        "field-screen", str(source), "--allow-conditional",
    ]) == 1
    assert "requires a strict conditional_review record" in capsys.readouterr().err


def test_field_screen_cli_propagates_all_declared_case_bounds_before_allowing_screening(
    tmp_path, capsys,
):
    source = tmp_path / "uncertain-field-case.json"
    output = tmp_path / "uncertain-field-execution.json"
    source.write_text(json.dumps(_case()), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--uncertainty-envelope",
        "--allow-conditional", "--relative-tolerance", "1.0",
        "--output", str(output), "--require-screening",
    ])
    captured = capsys.readouterr()
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert code == 0
    assert "conditional_allowed" in captured.out
    assert payload["input"]["uncertainty_envelope_requested"] is True
    assert payload["field_report"]["schema"] == "degali.field-operational-uncertainty-envelope.v1"
    assert payload["operational_screening"]["status"] == "conditional_allowed"
    assert len(payload["field_report"]["cases"]) == 4
    assert main(["field-verify", str(output)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True


def test_field_screen_cli_runs_sensor_calibration_operational_envelope(
    tmp_path, capsys,
):
    source = tmp_path / "sensor-envelope-case.json"
    output = tmp_path / "sensor-envelope-execution.json"
    case = _case()
    # Keep the source/weather and geometry exact so this test exercises the
    # detector-calibration envelope rather than a separate physical envelope.
    case["scenario"]["source"]["mass_flow_kg_s"] = _bounded(
        0.265, "kg/s", "FT-source",
    )
    case["scenario"]["weather"]["speed_m_s"] = _bounded(
        2.0, "m/s", "met-A",
    )
    case["obstacles"] = []
    for sensor in (
        case["scenario"]["sensor"],
        case["sensor_deployments"][0]["sensor"],
    ):
        sensor["response_time_s"] = _bounded(
            0.25, "s", "calibration-t90", relative=0.2,
        )
        sensor["gain"] = _bounded(
            1.0, "1", "calibration-gain", relative=0.1,
        )
        sensor["bias_mole_fraction"] = {
            "nominal": 0.001,
            "lower": 0.0,
            "upper": 0.002,
            "unit": "mole_fraction",
            "source": "calibration-bias",
        }
    source.write_text(json.dumps(case), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--sensor-array-envelope",
        "--allow-conditional", "--relative-tolerance", "1.0",
        "--max-uncertainty-cases", "64", "--output", str(output),
        "--require-screening",
    ])
    payload = json.loads(output.read_text(encoding="utf-8"))
    capsys.readouterr()

    assert code == 0
    assert payload["input"]["sensor_array_envelope_requested"] is True
    assert payload["input"]["uncertainty_envelope_kind"] == "sensor_calibration"
    assert payload["field_report"]["schema"] == (
        "degali.field-operational-sensor-array-envelope.v1"
    )
    assert payload["operational_screening"]["status"] == "conditional_allowed"
    assert len(payload["field_report"]["cases"]) == 64
    assert main(["field-verify", str(output)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True
    tampered = json.loads(output.read_text(encoding="utf-8"))
    tampered["field_report"]["cases"][0]["field_result"][
        "sensor_result"
    ]["maximum_true_mole_fraction"] += 1.0e-6
    output.write_text(json.dumps(tampered, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(output)]) == 1
    assert "does not match a fresh deterministic recomputation" in capsys.readouterr().err


def test_field_screen_cli_runs_joint_source_sensor_operational_envelope(
    tmp_path, capsys,
):
    source = tmp_path / "joint-source-sensor-case.json"
    output = tmp_path / "joint-source-sensor-execution.json"
    case = _case()
    case["scenario"]["source"]["mass_flow_kg_s"] = _bounded(
        0.265, "kg/s", "FT-source",
    )
    case["scenario"]["weather"]["speed_m_s"] = _bounded(
        2.0, "m/s", "met-A",
    )
    case["sensor_deployments"][0]["sensor"]["response_time_s"] = _bounded(
        0.25, "s", "H2-02-t90", relative=0.2,
    )
    source.write_text(json.dumps(case), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--joint-source-sensor-envelope",
        "--allow-conditional", "--relative-tolerance", "1.0",
        "--max-uncertainty-cases", "4", "--output", str(output),
        "--require-screening",
    ])
    payload = json.loads(output.read_text(encoding="utf-8"))
    capsys.readouterr()

    assert code == 0
    assert payload["input"]["joint_source_sensor_envelope_requested"] is True
    assert payload["input"]["uncertainty_envelope_kind"] == "source_sensor"
    assert payload["field_report"]["schema"] == (
        "degali.field-operational-source-sensor-envelope.v1"
    )
    assert payload["operational_screening"]["status"] == "conditional_allowed"
    assert len(payload["field_report"]["cases"]) == 2
    assert main(["field-verify", str(output)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True


def test_field_screen_cli_runs_a_fingerprinted_atmospheric_source_envelope(
    tmp_path, capsys,
):
    case_path = tmp_path / "field-case.json"
    source_csv = tmp_path / "atmospheric-schedule.csv"
    source_case_path = tmp_path / "atmospheric-source-case.json"
    output = tmp_path / "field-source-execution.json"
    case_path.write_text(json.dumps(_case()), encoding="utf-8")
    source_csv.write_text(
        "time_s,rate_kg_s,rate_lower_kg_s,rate_upper_kg_s\n"
        "0,0.01,0.005,0.02\n"
        "0.5,0.02,0.01,0.04\n"
        "1.0,0,0,0\n",
        encoding="utf-8",
    )
    source_case_path.write_text(json.dumps({
        "schema": "degali.field-atmospheric-schedule-input.v1",
        "csv_path": source_csv.name,
        "evidence": {
            "dataset_id": "source-case-A",
            "path": source_csv.name,
            "sha256": hashlib.sha256(source_csv.read_bytes()).hexdigest(),
            "row_count": 3,
            "source_boundary_id": "source-boundary-A",
            "common_clock_id": "clock-A",
            "source_kind": "post_flash_atmospheric_vapour",
        },
    }), encoding="utf-8")

    code = main([
        "field-screen", str(case_path),
        "--atmospheric-source-case", str(source_case_path),
        "--uncertainty-envelope", "--allow-conditional",
        "--relative-tolerance", "1.0", "--output", str(output),
        "--require-screening",
    ])
    payload = json.loads(output.read_text(encoding="utf-8"))
    capsys.readouterr()

    assert code == 0
    assert payload["input"]["uncertainty_envelope_kind"] == "atmospheric_source"
    assert payload["input"]["atmospheric_source_schedule_case"]["schedule"][
        "evidence"
    ]["dataset_id"] == "source-case-A"
    assert payload["field_report"]["field_uncertainty_envelope"][
        "atmospheric_source_schedule"
    ]["evidence"]["common_clock_id"] == "clock-A"
    assert {
        item["selection"]["source_schedule"]
        for item in payload["field_report"]["cases"]
    } == {"lower", "nominal", "upper"}


def test_field_screen_cli_requires_an_envelope_for_atmospheric_source_case(
    tmp_path, capsys,
):
    case_path = tmp_path / "field-case.json"
    case_path.write_text(json.dumps(_case()), encoding="utf-8")
    source_case_path = tmp_path / "atmospheric-source-case.json"
    source_case_path.write_text("{}", encoding="utf-8")

    assert main([
        "field-screen", str(case_path),
        "--atmospheric-source-case", str(source_case_path),
    ]) == 1
    assert "requires --uncertainty-envelope" in capsys.readouterr().err


def test_field_screen_cli_records_release_duration_uncertainty_corners(tmp_path, capsys):
    source = tmp_path / "duration-envelope-case.json"
    case = _case()
    case["scenario"]["source"]["mass_flow_kg_s"]["lower"] = 0.265
    case["scenario"]["source"]["mass_flow_kg_s"]["upper"] = 0.265
    case["scenario"]["weather"]["speed_m_s"]["lower"] = 2.0
    case["scenario"]["weather"]["speed_m_s"]["upper"] = 2.0
    case["scenario"]["source"]["duration_s"] = _bounded(
        1.0, "s", "release-timing-review-A", relative=0.25,
    )
    case["transport"]["duration_s"] = 1.25
    source.write_text(json.dumps(case), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--uncertainty-envelope", "--no-refinement",
        "--require-screening",
    ])
    payload = json.loads(capsys.readouterr().out)

    assert code == 2
    assert len(payload["field_report"]["cases"]) == 2
    assert {
        item["selection"]["duration_s"] for item in payload["field_report"]["cases"]
    } == {0.75, 1.25}


def test_field_screen_cli_propagates_declared_stability_alternatives(tmp_path, capsys):
    source = tmp_path / "stability-field-case.json"
    output = tmp_path / "stability-field-execution.json"
    case = _case()
    case["scenario"]["source"]["mass_flow_kg_s"]["lower"] = 0.265
    case["scenario"]["source"]["mass_flow_kg_s"]["upper"] = 0.265
    case["scenario"]["weather"]["speed_m_s"]["lower"] = 2.0
    case["scenario"]["weather"]["speed_m_s"]["upper"] = 2.0
    case["stability_mixing_closure"] = {
        "diffusivity_m2_s": {"neutral": 0.5, "stable": 0.2},
        "evidence_id": "met-mixing-A",
    }
    case["stability_alternatives"] = {
        "alternatives": ["stable"], "evidence_id": "met-stability-review-A",
    }
    source.write_text(json.dumps(case), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--uncertainty-envelope",
        "--max-uncertainty-cases", "4", "--allow-conditional",
        "--relative-tolerance", "1.0", "--output", str(output),
        "--require-screening",
    ])
    payload = json.loads(output.read_text(encoding="utf-8"))
    capsys.readouterr()

    assert code == 0
    cases = payload["field_report"]["cases"]
    assert {item["selection"]["weather_stability"] for item in cases} == {
        "neutral", "stable",
    }
    assert all(
        item["field_result"]["transport_input"]["stability_alternatives"]
        ["uncertainty_resolved_for_this_case"]
        for item in cases
    )


def test_field_screen_cli_keeps_a_provenanced_measured_source_withheld_until_enveloped(
    tmp_path, capsys,
):
    source = tmp_path / "measured-field-case.json"
    source.write_text(json.dumps(_measured_history_case(tmp_path)), encoding="utf-8")

    code = main(["field-screen", str(source), "--no-refinement", "--require-screening"])
    payload = json.loads(capsys.readouterr().out)

    assert code == 2
    assert payload["input"]["measured_history_imported"] is True
    assert payload["input"]["measured_history_event_id"] == "vent-event-A"
    transport_input = payload["field_report"]["transport_input"]
    assert transport_input["measured_history_provenance"]["event_id"] == "vent-event-A"
    assert transport_input["measured_history_source_uncertainty_resolved"] is False
    assert payload["operational_screening"]["status"] == "withheld"


def test_field_verify_rejects_changed_measured_history_csv_provenance(
    tmp_path, capsys,
):
    source = tmp_path / "measured-field-case.json"
    output = tmp_path / "measured-field-execution.json"
    source.write_text(json.dumps(_measured_history_case(tmp_path)), encoding="utf-8")

    assert main([
        "field-screen", str(source), "--no-refinement", "--output", str(output),
    ]) == 0
    capsys.readouterr()
    assert main(["field-verify", str(output)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True

    historian = tmp_path / "vent-event.csv"
    historian.write_text(
        historian.read_text(encoding="utf-8").replace("401000", "402000"),
        encoding="utf-8",
    )
    assert main(["field-verify", str(output)]) == 1
    assert "measured_history_provenance" in capsys.readouterr().err


def test_field_screen_cli_refines_and_aggregates_every_measured_history_corner(
    tmp_path, capsys,
):
    source = tmp_path / "measured-field-envelope-case.json"
    output = tmp_path / "measured-field-envelope-execution.json"
    case = _measured_history_case(tmp_path)
    case["scenario"]["weather"]["speed_m_s"]["lower"] = 2.0
    case["scenario"]["weather"]["speed_m_s"]["upper"] = 2.0
    case["scenario"]["source"]["location_m"] = [
        {
            "nominal": 0.0,
            "lower": -0.1,
            "upper": 0.1,
            "unit": "m",
            "source": "layout-review-measured-cli",
        },
        {
            "nominal": 0.0,
            "lower": 0.0,
            "upper": 0.0,
            "unit": "m",
            "source": "layout-review-measured-cli",
        },
        {
            "nominal": 0.5,
            "lower": 0.5,
            "upper": 0.5,
            "unit": "m",
            "source": "layout-review-measured-cli",
        },
    ]
    case["measured_history"]["map"]["pressure_pa"]["absolute_half_width"] = 0.0
    case["measured_history"]["map"]["temperature_k"]["absolute_half_width"] = 0.0
    source.write_text(json.dumps(case), encoding="utf-8")

    code = main([
        "field-screen", str(source), "--uncertainty-envelope",
        "--max-uncertainty-cases", "4", "--allow-conditional",
        "--relative-tolerance", "1.0", "--output", str(output),
        "--require-screening",
    ])
    captured = capsys.readouterr()
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert code == 0
    assert "conditional_allowed" in captured.out
    assert payload["input"]["uncertainty_envelope_kind"] == "joint_measured_history"
    assert payload["field_report"]["schema"] == "degali.field-operational-measured-history-envelope.v1"
    assert payload["operational_screening"]["status"] == "conditional_allowed"
    assert len(payload["field_report"]["cases"]) == 4
    assert {
        item["field_selection"]["source_location_x_m"]
        for item in payload["field_report"]["cases"]
    } == {-0.1, 0.1}

    assert main(["field-verify", str(output)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True

    tampered = json.loads(output.read_text(encoding="utf-8"))
    field_result = tampered["field_report"]["cases"][0]["field_result"]
    field_result["transport_result"]["maximum_concentration_kg_m3"] += 1.0e-6
    output.write_text(json.dumps(tampered, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(output)]) == 1
    assert "does not match a fresh deterministic recomputation" in capsys.readouterr().err


def test_field_screen_cli_refines_pressure_driven_history_and_replays_it(
    tmp_path, capsys,
):
    source = tmp_path / "pressure-driven-field-envelope-case.json"
    output = tmp_path / "pressure-driven-field-envelope-execution.json"
    case = _pressure_driven_history_case(tmp_path)
    case["pressure_driven_history"]["map"]["pressure_pa"]["absolute_half_width"] = 0.0
    case["pressure_driven_history"]["map"]["temperature_k"]["absolute_half_width"] = 0.0
    case["scenario"]["source"]["opening_area_m2"]["lower"] *= 0.9
    case["scenario"]["source"]["opening_area_m2"]["upper"] *= 1.1
    case["scenario"]["source"]["discharge_coefficient"]["lower"] = 0.72
    case["scenario"]["source"]["discharge_coefficient"]["upper"] = 0.88
    case["scenario"]["weather"]["speed_m_s"]["lower"] = 2.0
    case["scenario"]["weather"]["speed_m_s"]["upper"] = 2.0
    source.write_text(json.dumps(case), encoding="utf-8")

    nominal_case = read_field_screening_case_json(source)
    assert nominal_case.request.measured_history_source_uncertainty_resolved is False

    code = main([
        "field-screen", str(source), "--uncertainty-envelope",
        "--max-uncertainty-cases", "4", "--allow-conditional",
        "--relative-tolerance", "1.0", "--output", str(output),
        "--require-screening",
    ])
    captured = capsys.readouterr()
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert code == 0
    assert "conditional_allowed" in captured.out
    assert payload["input"]["uncertainty_envelope_kind"] == (
        "joint_pressure_driven_history"
    )
    assert payload["field_report"]["operational_screening"]["status"] == (
        "conditional_allowed"
    )
    assert payload["field_report"]["cases"][0]["field_result"]["transport_input"][
        "measured_history_provenance"
    ]["history_kind"] == "pressure_driven_orifice"
    assert main(["field-verify", str(output)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True
    historian = tmp_path / "pressure-driven-event.csv"
    historian.write_text(
        historian.read_text(encoding="utf-8").replace("390000", "391000"),
        encoding="utf-8",
    )
    assert main(["field-verify", str(output)]) == 1
    assert "measured_history_provenance" in capsys.readouterr().err


def test_field_screen_cli_routes_the_strict_phase_pool_case_to_the_complete_envelope(
    tmp_path, capsys, monkeypatch,
):
    source = tmp_path / "phase-pool-field-case.json"
    source.write_text(json.dumps(_phase_routing_case()), encoding="utf-8")
    decision = FieldOperationalScreeningDecision(
        "withheld", False, False, ("test holdback",), (),
        refinement_required=True, refinement_available=False,
        uncertainty_resolution_required=True, uncertainty_resolved=False,
    )
    captured_call = {}

    def fake_run(request, config, launch, **kwargs):
        captured_call.update(
            request=request, config=config, launch=launch, kwargs=kwargs,
        )
        return SimpleNamespace(operational_decision=decision)

    monkeypatch.setattr(
        field_operational_phase_transport,
        "run_field_operational_phase_routing_transport_envelope",
        fake_run,
    )
    monkeypatch.setattr(
        field_operational_phase_transport,
        "field_operational_phase_routing_transport_envelope_report",
        lambda _result: {
            "schema": "degali.field-operational-phase-routing-transport-envelope.v1",
            "cases": [],
        },
    )

    code = main([
        "field-screen", str(source), "--uncertainty-envelope",
        "--no-refinement", "--max-uncertainty-cases", "8", "--require-screening",
    ])
    payload = json.loads(capsys.readouterr().out)

    assert code == 2
    assert payload["input"]["uncertainty_envelope_kind"] == "phase_routing_transport"
    assert payload["input"]["phase_routing_transport_declared"] is True
    assert (
        payload["input"]["phase_routing_transport"]["phase_routing"]
        ["gas_model_options"]["crosswind_maximum_distance"]
        == pytest.approx(10.0)
    )
    assert (
        payload["input"]["phase_routing_transport"]["pool_launch"]["evidence_id"]
        == "pool-launch-review-A"
    )
    assert (
        payload["input"]["phase_routing_transport"]["phase_routing"]
        ["droplet_population"]["evidence_id"]
        == "spray-population-review-A"
    )
    assert payload["operational_screening"]["status"] == "withheld"
    assert captured_call["config"].pool_area_m2 == pytest.approx(0.5)
    assert captured_call["launch"].evidence_id == "pool-launch-review-A"
    assert captured_call["kwargs"]["pool_vertical_sigma_m"] == pytest.approx(0.2)
    assert captured_call["kwargs"]["table_nodes"] == 81
    assert captured_call["kwargs"]["max_cases"] == 8
    assert captured_call["kwargs"]["include_refinement"] is False

    execution = tmp_path / "phase-pool-execution.json"
    assert main([
        "field-screen", str(source), "--uncertainty-envelope", "--no-refinement",
        "--max-uncertainty-cases", "8", "--output", str(execution),
    ]) == 0
    capsys.readouterr()
    monkeypatch.setattr(
        field_execution,
        "run_field_operational_phase_routing_transport_envelope",
        fake_run,
    )
    monkeypatch.setattr(
        field_execution,
        "field_operational_phase_routing_transport_envelope_report",
        lambda _result: {
            "schema": "degali.field-operational-phase-routing-transport-envelope.v1",
            "cases": [],
        },
    )
    assert main(["field-verify", str(execution)]) == 0
    verification = json.loads(capsys.readouterr().out)
    assert verification["report_recomputed"] is True

    tampered = json.loads(execution.read_text(encoding="utf-8"))
    tampered["input"]["uncertainty_envelope_kind"] = "nominal_field"
    execution.write_text(json.dumps(tampered, indent=2) + "\n", encoding="utf-8")
    assert main(["field-verify", str(execution)]) == 1
    assert "complete phase envelope" in capsys.readouterr().err


def test_field_screen_cli_requires_the_complete_envelope_for_a_phase_pool_case(
    tmp_path, capsys,
):
    source = tmp_path / "phase-pool-field-case.json"
    source.write_text(json.dumps(_phase_routing_case()), encoding="utf-8")

    assert main(["field-screen", str(source)]) == 1
    assert "require --uncertainty-envelope" in capsys.readouterr().err
