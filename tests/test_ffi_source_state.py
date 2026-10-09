from dataclasses import replace
import json

import pytest

from degali.addons.field_contracts import BoundedValue, CircularBoundedValue
from degali.validation.ffi_source_state import (
    FFI_SOURCE_STATE_ENVELOPE_SCHEMA,
    FfiSourceState,
    ffi_source_state_envelope_report,
    run_ffi_source_state_envelope,
)
from degali.validation.spadeadam import load
from degali.cli import main


def _trial(test: int = 6):
    return next(item for item in load() if item.test == test)


def _state(trial):
    exact = FfiSourceState.from_trial(trial)
    return replace(
        exact,
        rate_kg_s=BoundedValue(
            trial.rate, 0.95 * trial.rate, 1.05 * trial.rate,
            unit="kg/s", source="test source-rate tolerance",
        ),
        wind_m_s=BoundedValue(
            trial.wind_low, trial.wind_low - 0.2, trial.wind_low + 0.2,
            unit="m/s", source="test low-mast wind tolerance",
        ),
        ambient_temperature_k=BoundedValue(
            277.15, 276.15, 278.15,
            unit="K", source="test ambient-temperature tolerance",
        ),
        wind_direction_from_deg=CircularBoundedValue(
            trial.wind_direction, trial.wind_direction - 2.0,
            trial.wind_direction + 2.0,
            unit="deg", source="test wind-direction tolerance",
        ),
    )


def test_ffi_source_state_envelope_keeps_sensor_residuals_as_deterministic_corners():
    trial = _trial()
    readings = tuple(
        reading for reading in trial.readings
        if abs(reading.radius - 30.0) < 0.5 and not reading.over_range
    )
    result = run_ffi_source_state_envelope(
        trial, _state(trial), readings=readings, max_cases=16,
    )

    assert result.status == "complete"
    assert len(result.cases) == 16
    assert all(case.metrics["available_rows"] == len(readings) for case in result.cases)
    report = ffi_source_state_envelope_report(result)
    assert report["schema"] == FFI_SOURCE_STATE_ENVELOPE_SCHEMA
    assert report["scope"].startswith("deterministic source/weather/ambient")
    assert report["disposition"]["operational_screening_allowed"] is False
    assert report["disposition"]["validation_qualified"] is False
    sensor = next(iter(report["sensor_envelope"].values()))
    assert sensor["observed_is_lower_bound"] is True
    assert sensor["predicted_lower_vol_pct"] <= sensor["predicted_upper_vol_pct"]
    assert sensor["residual_lower_vol_pct_point"] <= sensor["residual_upper_vol_pct_point"]
    assert sensor["worst_abs_residual_vol_pct_point"] >= 0.0
    assert set(sensor["worst_residual_selection"]) == set(result.source_state.uncertainty_fields())
    assert 0.0 <= sensor["lower_bound_satisfaction_fraction"] <= 1.0
    assert sensor["worst_lower_bound_deficit_vol_pct_point"] >= 0.0
    assert sensor["lower_bound_constraint_status"] in {"satisfied", "violated"}
    assert any(
        value["lower_bound_constraint_status"] == "violated"
        for value in report["sensor_envelope"].values()
    )
    assert any(
        value["worst_lower_bound_deficit_vol_pct_point"] > 1.0
        for value in report["sensor_envelope"].values()
    )
    assert all(
        case.metrics["lower_bound_rows"] == len(readings)
        for case in result.cases
    )
    assert all(row["observed_is_lower_bound"] for row in report["cases"][0]["rows"])
    operator = report["observation_operator"]
    assert len(operator["arc_max_table"]) == len(result.cases)
    assert len(operator["sensor_height_table"]) == 3 * len(result.cases)
    assert all(row["case_status"] == "complete" for row in operator["arc_max_table"])
    assert all(row["available_sensor_count"] == 15 for row in operator["arc_max_table"])
    assert all(row["withheld_sensor_count"] == 0 for row in operator["sensor_height_table"])
    assert all(
        row["lower_bound_constraint_status"] in {"satisfied", "violated"}
        for row in operator["arc_max_table"]
    )
    assert {row["height_m"] for row in operator["sensor_height_table"]} == {0.1, 1.0, 1.8}


def test_ffi_source_state_envelope_rejects_non_horizontal_and_unbounded_case_count():
    trial = _trial(1)
    with pytest.raises(ValueError, match="outdoor horizontal"):
        run_ffi_source_state_envelope(trial, FfiSourceState.from_trial(trial))

    horizontal = _trial()
    with pytest.raises(ValueError, match="exceed max_cases"):
        run_ffi_source_state_envelope(horizontal, _state(horizontal), max_cases=2)


def test_ffi_source_state_requires_evidence_for_every_boundary():
    trial = _trial()
    exact = FfiSourceState.from_trial(trial)
    with pytest.raises(ValueError, match="explicit evidence source"):
        replace(
            exact,
            ambient_temperature_k=BoundedValue(
                277.15, unit="K", source="unspecified",
            ),
        )


def test_ffi_source_state_requires_direction_and_scalar_bound_types():
    trial = _trial()
    exact = FfiSourceState.from_trial(trial)
    with pytest.raises(TypeError, match="wind_direction_from_deg"):
        replace(
            exact,
            wind_direction_from_deg=BoundedValue(
                trial.wind_direction, unit="deg", source="test",
            ),
        )
    with pytest.raises(TypeError, match="rate_kg_s"):
        replace(
            exact,
            rate_kg_s=CircularBoundedValue(
                trial.rate, unit="kg/s", source="test",
            ),
        )


def test_ffi_source_state_withholds_an_all_upwind_corner_instead_of_zero_filling():
    trial = _trial()
    exact = FfiSourceState.from_trial(trial)
    state = replace(
        exact,
        wind_direction_from_deg=CircularBoundedValue(
            90.0, 90.0, 90.0, unit="deg", source="test direction corner",
        ),
    )
    readings = tuple(
        reading for reading in trial.readings
        if abs(reading.radius - 30.0) < 0.5 and not reading.over_range
    )
    result = run_ffi_source_state_envelope(trial, state, readings=readings)

    assert result.status == "blocked"
    assert result.cases[0].status == "blocked"
    assert all(row.predicted_vol_pct is None for row in result.cases[0].rows)
    assert all(row.withheld_reason for row in result.cases[0].rows)
    operator = ffi_source_state_envelope_report(result)["observation_operator"]
    assert all(row["operator_status"] == "withheld" for row in operator["arc_max_table"])
    assert all(row["projected_sensor_max_vol_pct"] is None for row in operator["sensor_height_table"])
    assert all(
        row["lower_bound_constraint_status"] == "not_applicable"
        for row in operator["sensor_height_table"]
    )
    assert all(
        row["withheld_sensor_count"] == row["sensor_count"]
        for row in operator["sensor_height_table"]
    )


def test_ffi_source_state_cli_writes_versioned_execution_record(tmp_path):
    output = tmp_path / "ffi-source-state-execution.json"
    code = main([
        "ffi-source-state", "--test", "6", "--radius", "30",
        "--rate-bounds", "0.82,0.833,0.85",
        "--wind-bounds", "2.2,2.3,2.4", "--max-cases", "4",
        "--output", str(output), "--require-complete",
    ])

    assert code == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["schema"] == "degali.ffi-source-state-execution.v1"
    provenance = payload["input"]["reference_provenance"]
    assert set(provenance["files"]) == {"conditions.csv", "sensors.csv"}
    assert all(len(item["sha256"]) == 64 for item in provenance["files"].values())
    assert payload["source_state_envelope"]["schema"] == FFI_SOURCE_STATE_ENVELOPE_SCHEMA
    assert payload["source_state_envelope"]["status"] == "complete"
    assert len(payload["source_state_envelope"]["cases"]) == 4
