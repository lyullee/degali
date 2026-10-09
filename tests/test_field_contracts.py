import pytest
import math

from degali.addons.field_contracts import (
    BoundedValue,
    CircularBoundedValue,
    FieldCoordinateReference,
    FieldApplicability,
    FieldScenario,
    FieldValidationEvidence,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
    assess_field_applicability,
    bounded,
)


def test_bounded_value_keeps_interval_without_probability_claim():
    value = bounded(10.0, relative=0.1, unit="bar", source="calibration")
    assert value.corners() == pytest.approx((9.0, 11.0))
    assert not value.is_exact
    assert value.as_dict()["source"] == "calibration"


def test_bounded_numeric_boundaries_reject_boolean_coercion():
    with pytest.raises(ValueError, match="not boolean"):
        BoundedValue(True)
    with pytest.raises(ValueError, match="not boolean"):
        CircularBoundedValue(True)


def test_field_boundary_records_reject_untyped_nested_values():
    with pytest.raises(TypeError, match="upstream_pressure"):
        ReleaseSource(upstream_pressure=1.0)
    with pytest.raises(TypeError, match="speed_m_s"):
        WeatherState(speed_m_s=2.0)
    with pytest.raises(TypeError, match="response_time_s"):
        SensorModel((0.0, 0.0, 0.0), response_time_s=1.0)
    with pytest.raises(TypeError, match="heat_transfer_w_m2_k"):
        SurfaceBoundary(heat_transfer_w_m2_k=10.0)
    with pytest.raises(TypeError, match="source must be"):
        FieldScenario(source="not-a-release")


@pytest.mark.parametrize(
    "factory",
    [
        lambda: ReleaseSource(upstream_pressure=BoundedValue(101325.0, unit="K")),
        lambda: ReleaseSource(upstream_temperature=BoundedValue(293.15, unit="Pa")),
        lambda: ReleaseSource(mass_flow_kg_s=BoundedValue(1.0, unit="Pa")),
        lambda: ReleaseSource(opening_area_m2=BoundedValue(1.0e-6, unit="Pa")),
        lambda: ReleaseSource(discharge_coefficient=BoundedValue(0.8, unit="kg/s")),
        lambda: ReleaseSource(liquid_fraction=BoundedValue(1.0, unit="kg/s")),
        lambda: ReleaseSource(
            duration_s=1.0,
            duration_uncertainty=BoundedValue(1.0, 0.5, 1.5, unit="kg/s", source="timing"),
        ),
        lambda: WeatherState(speed_m_s=BoundedValue(2.0, unit="Pa")),
        lambda: WeatherState(direction_deg=BoundedValue(180.0, unit="m/s")),
        lambda: WeatherState(
            direction_deg=CircularBoundedValue(180.0, unit="m/s")
        ),
        lambda: SensorModel(
            (0.0, 0.0, 0.0), response_time_s=BoundedValue(1.0, unit="m/s")
        ),
        lambda: SensorModel((0.0, 0.0, 0.0), gain=BoundedValue(1.0, unit="s")),
        lambda: SensorModel(
            (0.0, 0.0, 0.0), bias_mole_fraction=BoundedValue(0.0, unit="s")
        ),
        lambda: SurfaceBoundary(heat_transfer_w_m2_k=BoundedValue(10.0, unit="K")),
        lambda: SurfaceBoundary(surface_temperature_k=BoundedValue(293.15, unit="Pa")),
    ],
)
def test_field_boundary_records_reject_explicitly_wrong_units(factory):
    with pytest.raises(ValueError, match="must use unit"):
        factory()


def test_release_source_validates_direction_and_exposes_fields():
    source = ReleaseSource(
        mass_flow_kg_s=BoundedValue(1.0, 0.8, 1.2, "kg/s", "historian"),
        liquid_fraction=BoundedValue(0.7, 0.5, 0.9),
        duration_s=5.0,
    )
    assert source.is_transient
    assert source.direction_unit == pytest.approx((1.0, 0.0, 0.0))
    assert set(source.uncertainty_fields()) == {
        "upstream_pressure", "upstream_temperature", "mass_flow_kg_s",
        "opening_area_m2", "discharge_coefficient", "liquid_fraction",
    }


def test_release_source_metadata_is_an_immutable_snapshot():
    metadata = {"release_boundary_evidence_id": "review-A"}
    source = ReleaseSource(metadata=metadata)
    metadata["release_boundary_evidence_id"] = "tampered"

    assert source.metadata["release_boundary_evidence_id"] == "review-A"
    with pytest.raises(TypeError):
        source.metadata["new"] = "mutation"
    with pytest.raises(ValueError, match="strings to strings"):
        ReleaseSource(metadata={"bad": 1})


def test_release_duration_uncertainty_is_propagated_and_fixed_at_each_corner():
    scenario = FieldScenario(
        source=ReleaseSource(
            mass_flow_kg_s=BoundedValue(1.0),
            duration_s=10.0,
            duration_uncertainty=BoundedValue(10.0, 8.0, 12.0, "s", "timing-review"),
        ),
        temporal_mode="transient",
    )

    assert scenario.uncertainty_fields()["duration_s"].corners() == pytest.approx((8.0, 12.0))
    corners = scenario.corner_cases()
    assert {corner["duration_s"] for corner in corners} == {8.0, 12.0}
    fixed = scenario.at_corner(next(corner for corner in corners if corner["duration_s"] == 8.0))
    assert fixed.source.duration_s == pytest.approx(8.0)
    assert fixed.source.duration_uncertainty is None
    assert "duration_s" not in fixed.uncertainty_fields()


def test_release_duration_uncertainty_requires_matching_positive_nominal():
    with pytest.raises(ValueError, match="nominal must match duration_s"):
        ReleaseSource(
            duration_s=10.0,
            duration_uncertainty=BoundedValue(9.0, 8.0, 10.0, "s", "timing-review"),
        )


def test_release_location_uncertainty_is_propagated_to_geometry_corners():
    scenario = FieldScenario(
        source=ReleaseSource(
            location_m=(1.0, 2.0, 0.5),
            location_uncertainty_m=(
                BoundedValue(1.0, 0.9, 1.1, "m", "layout-review"),
                BoundedValue(2.0, 1.8, 2.2, "m", "layout-review"),
                BoundedValue(0.5, 0.45, 0.55, "m", "elevation-review"),
            ),
            discharge_coefficient=BoundedValue(0.8),
        ),
        surface=SurfaceBoundary(heat_transfer_w_m2_k=BoundedValue(10.0)),
    )

    corners = scenario.corner_cases(max_cases=8)
    assert len(corners) == 8
    assert {corner["source_location_x_m"] for corner in corners} == {0.9, 1.1}
    fixed = scenario.at_corner(corners[0])
    assert fixed.source.location_uncertainty_m is None
    assert all(value.is_exact for value in fixed.source.uncertainty_fields().values())


def test_release_location_uncertainty_requires_matching_nominal_and_explicit_source():
    with pytest.raises(ValueError, match=r"location_uncertainty_m\[0\].nominal"):
        ReleaseSource(
            location_m=(1.0, 2.0, 0.5),
            location_uncertainty_m=(
                BoundedValue(1.1, 0.9, 1.1, "m", "layout-review"),
                BoundedValue(2.0, 2.0, 2.0, "m", "layout-review"),
                BoundedValue(0.5, 0.5, 0.5, "m", "layout-review"),
            ),
        )


def test_release_geometry_uncertainty_requires_physical_units():
    with pytest.raises(ValueError, match=r"location_uncertainty_m\[0\].*unit 'm'"):
        ReleaseSource(
            location_uncertainty_m=(
                BoundedValue(0.0, -0.1, 0.1, "kg/s", "layout-review"),
                BoundedValue(0.0, 0.0, 0.0, "m", "layout-review"),
                BoundedValue(0.0, 0.0, 0.0, "m", "layout-review"),
            ),
        )
    with pytest.raises(ValueError, match=r"direction_uncertainty_m\[0\].*unit '1'"):
        ReleaseSource(
            direction_uncertainty_m=(
                BoundedValue(1.0, 1.0, 1.0, "m", "orientation-review"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review"),
            ),
        )


def test_release_direction_uncertainty_propagates_nonzero_vector_corners():
    scenario = FieldScenario(
        source=ReleaseSource(
            direction_m=(1.0, 0.0, 0.0),
            direction_uncertainty_m=(
                BoundedValue(1.0, 0.8, 1.2, "1", "orientation-review"),
                BoundedValue(0.0, -0.1, 0.1, "1", "orientation-review"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review"),
            ),
            discharge_coefficient=BoundedValue(0.8),
        ),
        surface=SurfaceBoundary(heat_transfer_w_m2_k=BoundedValue(10.0)),
    )

    corners = scenario.corner_cases(max_cases=4)
    assert len(corners) == 4
    fixed = scenario.at_corner(corners[0])
    assert fixed.source.direction_uncertainty_m is None
    assert fixed.source.direction_unit != (0.0, 0.0, 0.0)


def test_release_direction_uncertainty_rejects_a_zero_vector_corner():
    with pytest.raises(ValueError, match="zero-vector corner"):
        ReleaseSource(
            direction_m=(1.0, 0.0, 0.0),
            direction_uncertainty_m=(
                BoundedValue(1.0, 0.0, 1.0, "1", "orientation-review"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review"),
            ),
        )


def test_field_status_is_fail_closed():
    with pytest.raises(ValueError):
        FieldApplicability("blocked")
    assert FieldApplicability("conditional", warnings=("model form",)).status == "conditional"
    with pytest.raises(TypeError, match="uncertainty_complete"):
        FieldApplicability("accepted", uncertainty_complete=1)
    with pytest.raises(ValueError, match="warnings"):
        FieldApplicability("conditional", warnings=("",))
    with pytest.raises(ValueError, match="complete uncertainty"):
        FieldApplicability("blocked", reasons=("blocked",), uncertainty_complete=True)
    with pytest.raises(ValueError, match="blocking reasons"):
        FieldApplicability("conditional", reasons=("not resolved",))
    assert FieldApplicability(
        "accepted", warnings=("scope limit",), uncertainty_complete=True
    ).warnings == ("scope limit",)


def test_sensor_and_weather_boundaries_are_checked():
    assert WeatherState().stability == "neutral"
    assert WeatherState(direction_deg=BoundedValue(270.0)).wind_to_math_radians() == pytest.approx(0.0)
    assert WeatherState(
        direction_deg=BoundedValue(90.0), direction_convention="math_to"
    ).wind_to_math_radians() == pytest.approx(math.pi / 2.0)
    with pytest.raises(ValueError):
        SensorModel((0.0, 0.0, -0.1))


def test_release_and_weather_boundaries_reject_boolean_geometry_and_height():
    with pytest.raises(ValueError, match="location_m"):
        ReleaseSource(location_m=(0.0, True, 0.0))
    with pytest.raises(ValueError, match="duration_s"):
        ReleaseSource(duration_s=True)
    with pytest.raises(ValueError, match="reference_height_m"):
        WeatherState(reference_height_m=True)


def test_sensor_contract_rejects_nonphysical_numeric_coercions():
    with pytest.raises(ValueError, match="position_m"):
        SensorModel((0.0, True, 0.0))
    with pytest.raises(ValueError, match="averaging time"):
        SensorModel((0.0, 0.0, 0.0), averaging_time_s=float("nan"))
    with pytest.raises(ValueError, match="averaging time"):
        SensorModel((0.0, 0.0, 0.0), averaging_time_s=True)
    with pytest.raises(ValueError, match="gain"):
        SensorModel(
            (0.0, 0.0, 0.0),
            gain=BoundedValue(0.5, -0.1, 0.5, source="bad-calibration"),
        )


def test_coordinate_reference_rotates_true_wind_bearings_into_the_site_grid():
    reference = FieldCoordinateReference(
        coordinate_system_id="plant-grid-rev-A", origin_id="transfer-skid-origin",
        x_axis_bearing_math_to_deg=90.0, vertical_datum_id="site-grade",
        evidence_id="layout-drawing-A",
    )

    assert reference.earth_to_local_math_radians(0.0) == pytest.approx(1.5 * math.pi)
    assert reference.x_axis_bearing_math_to_deg == pytest.approx(90.0)
    with pytest.raises(ValueError, match="coordinate_system_id"):
        FieldCoordinateReference("", "origin", 0.0, "datum", "evidence")


def test_circular_wind_direction_interval_crosses_north_without_reversing_the_sector():
    direction = CircularBoundedValue(0.0, 350.0, 10.0, "deg", "met-mast")
    weather = WeatherState(direction_deg=direction)
    scenario = FieldScenario(
        source=ReleaseSource(duration_s=1.0), weather=weather, temporal_mode="transient",
    )

    assert direction.span_deg == pytest.approx(20.0)
    assert direction.corners() == pytest.approx((350.0, 10.0))
    assert weather.wind_to_math_radians() == pytest.approx(1.5 * math.pi)
    assert {corner["wind_direction_deg"] for corner in scenario.corner_cases()} == {350.0, 10.0}
    assert scenario.as_dict()["weather"]["uncertainty"]["wind_direction_deg"]["crosses_zero_deg"]

    with pytest.raises(ValueError, match="inside the declared circular interval"):
        CircularBoundedValue(90.0, 350.0, 10.0, "deg", "met-mast")


def test_field_scenario_exposes_deterministic_uncertainty_corners_and_gate():
    scenario = FieldScenario(
        source=ReleaseSource(
            mass_flow_kg_s=BoundedValue(1.0, 0.9, 1.1),
            duration_s=10.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0, 1.0, 3.0)),
        surface=SurfaceBoundary(),
        temporal_mode="transient",
    )
    corners = scenario.corner_cases(max_cases=1024)
    assert len(corners) == 2 ** sum(
        not value.is_exact for value in scenario.uncertainty_fields().values()
    )
    status = assess_field_applicability(scenario, obstacle_present=True)
    assert status.status == "conditional"
    assert not status.reasons
    assert status.uncertainty_complete
    fixed = scenario.at_corner(corners[0])
    assert all(value.is_exact for value in fixed.uncertainty_fields().values())


def test_field_scenario_blocks_unresolved_flash():
    scenario = FieldScenario(source=ReleaseSource(flash_model="unresolved"))
    status = assess_field_applicability(scenario)
    assert status.status == "blocked"
    assert "flash" in status.reasons[0]


def test_field_gate_requires_absolute_pressure_and_declared_phase_basis():
    scenario = FieldScenario(source=ReleaseSource(pressure_reference="gauge"))
    status = assess_field_applicability(scenario)
    assert status.status == "blocked"
    assert "absolute" in status.reasons[0]

    scenario = FieldScenario(source=ReleaseSource(liquid_fraction_basis="unspecified"))
    status = assess_field_applicability(scenario)
    assert status.status == "blocked"
    assert "basis" in status.reasons[0]


def test_validation_flag_requires_fingerprinted_matched_evidence():
    scenario = FieldScenario(source=ReleaseSource())
    blocked = assess_field_applicability(
        scenario, lh2_validation_available=True,
    )
    assert blocked.status == "blocked"
    assert "FieldValidationEvidence" in blocked.reasons[0]

    evidence = FieldValidationEvidence(
        dataset_id="lh2-obstacle-trial-A",
        path="evidence/trial-A.csv",
        sha256="0" * 64,
        row_count=12,
        source_boundary_id="source-A",
        weather_id="met-A",
        obstacle_geometry_id="obstacle-A",
        receptor_geometry_id="detectors-A",
        temporal_operator_id="t90-average-A",
        common_clock_id="clock-A",
    )
    accepted = assess_field_applicability(
        scenario,
        obstacle_present=True,
        lh2_validation_available=True,
        validation_evidence=evidence,
    )
    assert accepted.status == "conditional"
    assert not any("no site-specific" in warning for warning in accepted.warnings)


def test_free_field_validation_evidence_does_not_qualify_an_obstacle():
    evidence = FieldValidationEvidence(
        dataset_id="lh2-free-field-A",
        path="evidence/free-field-A.csv",
        sha256="1" * 64,
        row_count=5,
        source_boundary_id="source-A",
        weather_id="met-A",
        obstacle_geometry_id="none",
        receptor_geometry_id="detectors-A",
        temporal_operator_id="t90-average-A",
        common_clock_id="clock-A",
        scope="lh2_free_field",
    )
    status = assess_field_applicability(
        FieldScenario(source=ReleaseSource()),
        obstacle_present=True,
        validation_evidence=evidence,
    )
    assert status.status == "conditional"
    assert any("free-field only" in warning for warning in status.warnings)


def test_fingerprinted_validation_evidence_without_a_score_is_conditional():
    evidence = FieldValidationEvidence(
        dataset_id="lh2-free-field-score-pending",
        path="evidence/free-field-score-pending.csv",
        sha256="2" * 64,
        row_count=5,
        source_boundary_id="source-A",
        weather_id="met-A",
        obstacle_geometry_id="none",
        receptor_geometry_id="detectors-A",
        temporal_operator_id="t90-average-A",
        common_clock_id="clock-A",
        scope="lh2_free_field",
    )

    status = assess_field_applicability(
        FieldScenario(source=ReleaseSource()),
        validation_evidence=evidence,
    )

    assert status.status == "conditional"
    assert any("matched field-validation score" in warning for warning in status.warnings)
