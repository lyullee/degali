from dataclasses import replace
import math

import pytest

from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    WeatherState,
)
from degali.addons.field_historian_io import (
    HistorianCsvChannel,
    MeasuredHistoryCsvMap,
    PressureDrivenHistoryCsvMap,
    direct_vapour_schedule_from_imported_history,
    direct_vapour_schedule_from_imported_pressure_driven_history,
    read_measured_history_csv,
    read_pressure_driven_history_csv,
)
from degali.addons.field_history import (
    MeasuredHistoryQualityCriteria,
    request_with_measured_flash_schedule,
)
from degali.addons.field_report import field_screening_report
from degali.addons.field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening
from degali.addons.semi_fv_obstacle import SemiFVConfig


def _mapping():
    return MeasuredHistoryCsvMap(
        time_s_column="time_s", event_id="vent-event-2026-10-05-01",
        event_evidence_id="event-window-review-01",
        phase_evidence_id="vent-state-review-01",
        pressure_pa=HistorianCsvChannel(
            "pressure_pa", "Pa", "PT-1105A", "PT-1105A-cal-2026", absolute_half_width=500.0,
        ),
        temperature_k=HistorianCsvChannel(
            "temperature_k", "K", "TT-1106", "TT-1106-cal-2026", absolute_half_width=0.25,
        ),
        mass_flow_kg_s=HistorianCsvChannel(
            "mass_flow_kg_s", "kg/s", "FT-1101", "FT-1101-cal-2026", relative_half_width=0.05,
        ),
        liquid_fraction=HistorianCsvChannel(
            "liquid_fraction", "1", "LT-1101", "LT-1101-basis-2026", absolute_half_width=0.0,
        ),
    )


def _pressure_mapping():
    return PressureDrivenHistoryCsvMap(
        time_s_column="time_s", event_id="pressure-event-2026-10-05-01",
        event_evidence_id="pressure-event-window-review-01",
        phase_evidence_id="pressure-phase-review-01",
        pressure_pa=HistorianCsvChannel(
            "pressure_pa", "Pa", "PT-ORIFICE-01", "PT-ORIFICE-cal-2026", absolute_half_width=500.0,
        ),
        temperature_k=HistorianCsvChannel(
            "temperature_k", "K", "TT-ORIFICE-01", "TT-ORIFICE-cal-2026", absolute_half_width=0.05,
        ),
    )


def _release():
    return ReleaseSource(
        fluid="lh2", upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.2, unit="kg/s"),
        opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"),
        discharge_coefficient=BoundedValue(0.8), liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium", duration_s=1.0,
    )


def _quality():
    return MeasuredHistoryQualityCriteria(
        maximum_sample_interval_s=1.0, maximum_response_time_s=1.0,
        maximum_absolute_time_offset_s=0.1, maximum_relative_half_width=0.1,
        evidence_id="historian-quality-procedure-01",
    )


def test_csv_history_import_requires_explicit_si_mapping_and_preserves_provenance(tmp_path):
    source = tmp_path / "event.csv"
    source.write_text(
        "time_s,pressure_pa,temperature_k,mass_flow_kg_s,liquid_fraction\n"
        "0,400000,26.084,0.2,0.922\n"
        "1,401000,26.2,0,0.922\n",
        encoding="utf-8",
    )
    imported = read_measured_history_csv(source, _mapping())

    assert imported.history.release_time_s == (0.0, 1.0)
    assert imported.history.mass_flow_kg_s.lower == pytest.approx((0.19, 0.0))
    assert imported.history.pressure_pa.lower == pytest.approx((399500.0, 400500.0))
    assert imported.provenance.row_count == 2
    assert imported.provenance.event_id == "vent-event-2026-10-05-01"
    assert len(imported.provenance.sha256) == 64

    schedule = direct_vapour_schedule_from_imported_history(
        _release(), imported, quality_criteria=_quality(),
    )
    assert schedule.quality_assessment is not None
    assert schedule.quality_assessment.approved
    assert any("sha256=" in warning for warning in schedule.warnings)
    accelerated = direct_vapour_schedule_from_imported_history(
        _release(), imported, quality_criteria=_quality(),
        auto_property_table_minimum_intervals=1,
    )
    assert dict(accelerated.provenance)["property_table_used"] == "true"
    assert dict(accelerated.provenance)["property_table_auto_built"] == "true"
    assert any("saturation table" in warning for warning in accelerated.warnings)

    request = FieldSemiFVRequest(
        FieldScenario(
            source=_release(),
            weather=WeatherState(
                speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0),
            ),
            sensor=SensorModel((1.0, 0.0, 0.5), averaging_time_s=0.1),
            temporal_mode="transient",
        ),
        transport=SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            time_step_s=0.005, duration_s=1.0, source_sigma_m=0.2,
        ),
    )
    attached = request_with_measured_flash_schedule(request, schedule)
    report = field_screening_report(run_field_semi_fv_screening(attached))
    provenance = report["transport_input"]["measured_history_provenance"]
    assert provenance["import_format"] == "degali.measured-history-csv.v1"
    assert provenance["source_sha256"] == imported.provenance.sha256
    assert provenance["phase_evidence_id"] == "vent-state-review-01"
    assert provenance["channel.mass_flow_kg_s.calibration_evidence_id"] == "FT-1101-cal-2026"


def test_csv_history_import_rejects_a_missing_required_source_channel(tmp_path):
    source = tmp_path / "pressure-only.csv"
    source.write_text(
        "time_s,pressure_pa,temperature_k\n0,400000,26\n1,401000,26.2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="mass_flow_kg_s"):
        read_measured_history_csv(source, _mapping())


def test_csv_history_import_requires_explicit_phase_evidence():
    with pytest.raises(ValueError, match="phase_evidence_id"):
        replace(_mapping(), phase_evidence_id="unspecified")


def test_pressure_driven_csv_import_keeps_pressure_as_boundary_and_fingerprints_file(tmp_path):
    source = tmp_path / "pressure-driven-event.csv"
    source.write_text(
        "time_s,pressure_pa,temperature_k\n"
        "0,400000,26.076\n"
        "1,390000,25.949\n",
        encoding="utf-8",
    )
    imported = read_pressure_driven_history_csv(source, _pressure_mapping())

    assert imported.history.duration_s == pytest.approx(1.0)
    assert not hasattr(imported.history, "mass_flow_kg_s")
    schedule = direct_vapour_schedule_from_imported_pressure_driven_history(
        replace(_release(), mass_flow_kg_s=BoundedValue(0.0, unit="kg/s")),
        imported,
        quality_criteria=_quality(),
    )
    assert schedule.total_measured_mass_kg > 0.0
    assert schedule.quality_assessment is not None
    assert schedule.quality_assessment.approved
    provenance = dict(schedule.provenance)
    assert provenance["history_kind"] == "pressure_driven_orifice"
    assert provenance["import_format"] == "degali.pressure-driven-history-csv.v1"
    assert provenance["source_sha256"] == imported.provenance.sha256


def test_pressure_driven_csv_nominal_schedule_keeps_orifice_bounds_unresolved(tmp_path):
    source = tmp_path / "pressure-driven-bounded-event.csv"
    source.write_text(
        "time_s,pressure_pa,temperature_k\n"
        "0,400000,26.076\n"
        "1,390000,25.949\n",
        encoding="utf-8",
    )
    imported = read_pressure_driven_history_csv(source, _pressure_mapping())
    base = replace(_release(), mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"))
    area = base.opening_area_m2.nominal
    release = replace(
        base,
        opening_area_m2=BoundedValue(area, 0.9 * area, 1.1 * area, "m2", "area-review"),
        discharge_coefficient=BoundedValue(0.8, 0.72, 0.88, "1", "cd-review"),
    )

    schedule = direct_vapour_schedule_from_imported_pressure_driven_history(
        release, imported, quality_criteria=_quality(),
    )

    assert schedule.total_measured_mass_kg > 0.0
    assert schedule.source_uncertainty_resolved is False
    assert any("bounds are unresolved" in warning for warning in schedule.warnings)


def test_csv_history_import_refuses_nonphysical_calibrated_bounds(tmp_path):
    source = tmp_path / "bad-pressure.csv"
    source.write_text(
        "time_s,pressure_pa,temperature_k,mass_flow_kg_s,liquid_fraction\n"
        "0,400,26,0.2,0.922\n1,401,26.2,0,0.922\n",
        encoding="utf-8",
    )
    mapping = MeasuredHistoryCsvMap(
        time_s_column="time_s", event_id="event-a", event_evidence_id="event-a-review",
        phase_evidence_id="event-a-phase-review",
        pressure_pa=HistorianCsvChannel(
            "pressure_pa", "Pa", "PT", "PT-cal", absolute_half_width=500.0,
        ),
        temperature_k=_mapping().temperature_k,
        mass_flow_kg_s=_mapping().mass_flow_kg_s,
        liquid_fraction=_mapping().liquid_fraction,
    )

    with pytest.raises(ValueError, match="pressure_pa lower bound"):
        read_measured_history_csv(source, mapping)


@pytest.mark.parametrize(
    ("pressure_values", "message"),
    (
        ((400500, 401000, 400000), "bounds must satisfy lower <= upper"),
        ((400500, 400000, 400100), "nominal value lies outside declared bounds"),
    ),
)
def test_csv_history_import_rejects_inconsistent_paired_calibration_bounds(
    tmp_path, pressure_values, message,
):
    source = tmp_path / "inconsistent-pressure-bounds.csv"
    nominal, lower, upper = pressure_values
    source.write_text(
        "time_s,pressure_pa,pressure_lower_pa,pressure_upper_pa,temperature_k,mass_flow_kg_s\n"
        f"0,{nominal},{lower},{upper},26,0.2\n"
        "1,401000,400500,401500,26.2,0\n",
        encoding="utf-8",
    )
    mapping = replace(
        _mapping(),
        pressure_pa=HistorianCsvChannel(
            "pressure_pa", "Pa", "PT", "PT-cal",
            lower_column="pressure_lower_pa", upper_column="pressure_upper_pa",
        ),
        liquid_fraction=None,
    )

    with pytest.raises(ValueError, match=message):
        read_measured_history_csv(source, mapping)


@pytest.mark.parametrize(
    ("name", "contents", "message"),
    (
        (
            "duplicate-header.csv",
            "time_s,pressure_pa,temperature_k,mass_flow_kg_s,mass_flow_kg_s\n"
            "0,400000,26,0.2,0.2\n1,401000,26.2,0,0\n",
            "duplicate columns",
        ),
        (
            "empty-header.csv",
            "time_s,pressure_pa,,mass_flow_kg_s\n0,400000,26,0.2\n1,401000,26.2,0\n",
            "empty column",
        ),
        (
            "extra-value.csv",
            "time_s,pressure_pa,temperature_k,mass_flow_kg_s\n"
            "0,400000,26,0.2,unexpected\n1,401000,26.2,0\n",
            "different number of fields",
        ),
    ),
)
def test_csv_history_import_rejects_ambiguous_csv_structure(
    tmp_path, name, contents, message,
):
    source = tmp_path / name
    source.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        read_measured_history_csv(source, _mapping())
