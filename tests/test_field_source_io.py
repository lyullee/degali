import hashlib
import json

import pytest

from degali.cli import main
from degali.addons.field_source_io import (
    FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA,
    FieldAtmosphericSourceSchedule,
    FieldSourceScheduleEvidence,
    field_source_schedule_case_from_mapping,
    field_source_schedule_from_csv,
    read_field_source_schedule_json,
    write_field_source_schedule_json,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule


def _csv(tmp_path, *, final_rate=0.0):
    path = tmp_path / "schedule.csv"
    path.write_text(
        "time_s,rate_kg_s\n"
        f"0,0.2\n1,0.1\n3,{final_rate}\n",
        encoding="utf-8",
    )
    return path


def _evidence(path, *, source_kind="post_flash_atmospheric_vapour"):
    return FieldSourceScheduleEvidence(
        dataset_id="schedule-a",
        path=str(path.resolve()),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        row_count=3,
        source_boundary_id="source-a",
        common_clock_id="clock-a",
        source_kind=source_kind,
    )


def test_atmospheric_schedule_csv_is_fingerprinted_without_phase_inference(tmp_path):
    path = _csv(tmp_path)
    imported = field_source_schedule_from_csv(path, evidence=_evidence(path))

    assert imported.schedule.source_id == "source-a"
    assert imported.schedule.released_mass_kg == pytest.approx(0.4)
    assert imported.evidence.common_clock_id == "clock-a"
    assert "no liquid/phase inference" in imported.as_record()["scope"]


def test_linear_atmospheric_schedule_operator_is_fingerprinted_and_mass_consistent(tmp_path):
    path = _csv(tmp_path)
    imported = field_source_schedule_from_csv(
        path, evidence=_evidence(path), rate_operator="linear",
    )

    assert imported.schedule.rate_operator == "linear"
    assert imported.schedule.released_mass_kg == pytest.approx(0.25)
    assert imported.as_record()["schedule"]["rate_operator"] == "linear"


def test_strict_atmospheric_schedule_case_preserves_linear_operator(tmp_path):
    path = _csv(tmp_path)
    evidence = _evidence(path)
    case = {
        "schema": FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA,
        "csv_path": path.name,
        "rate_operator": "linear",
        "evidence": {
            "dataset_id": evidence.dataset_id,
            "path": path.name,
            "sha256": evidence.sha256,
            "row_count": evidence.row_count,
            "source_boundary_id": evidence.source_boundary_id,
            "common_clock_id": evidence.common_clock_id,
            "source_kind": evidence.source_kind,
        },
    }

    parsed = field_source_schedule_case_from_mapping(case, base_directory=tmp_path)
    assert parsed.schedule.rate_operator == "linear"
    assert parsed.schedule.released_mass_kg == pytest.approx(0.25)


def test_declared_atmospheric_source_kind_is_preserved(tmp_path):
    path = _csv(tmp_path)
    imported = field_source_schedule_from_csv(
        path,
        evidence=_evidence(path, source_kind="declared_atmospheric_vapour"),
    )

    assert imported.evidence.source_kind == "declared_atmospheric_vapour"
    assert imported.as_record()["evidence"]["source_kind"] == "declared_atmospheric_vapour"


def test_typed_source_schedule_writer_pins_and_round_trips_case(tmp_path):
    path = _csv(tmp_path)
    parsed = field_source_schedule_from_csv(path, evidence=_evidence(path))
    output = tmp_path / "exported" / "source-case.json"

    exported = write_field_source_schedule_json(parsed, output)

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["csv_path"] == "../schedule.csv"
    assert payload["evidence"]["path"] == "schedule.csv"
    assert len(payload["evidence"]["sha256"]) == 64
    assert exported.schedule == parsed.schedule
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_field_source_schedule_json(parsed, output)


def test_rate_bounds_are_preserved_as_deterministic_source_corners(tmp_path):
    path = tmp_path / "bounded-schedule.csv"
    path.write_text(
        "time_s,rate_kg_s,rate_lower_kg_s,rate_upper_kg_s\n"
        "0,0.2,0.1,0.3\n1,0.1,0.05,0.2\n3,0,0,0\n",
        encoding="utf-8",
    )
    imported = field_source_schedule_from_csv(path, evidence=_evidence(path))

    assert [label for label, _schedule in imported.corner_schedules()] == [
        "lower", "nominal", "upper",
    ]
    assert imported.corner_schedules()[0][1].released_mass_kg == pytest.approx(0.2)
    assert imported.as_record()["uncertainty"]["interpretation"].startswith("deterministic")


def test_schedule_requires_a_zero_final_endpoint(tmp_path):
    path = _csv(tmp_path, final_rate=0.1)
    with pytest.raises(ValueError, match="final endpoint rate"):
        field_source_schedule_from_csv(path, evidence=_evidence(path))


def test_schedule_rejects_a_zero_mass_atmospheric_source(tmp_path):
    path = tmp_path / "zero-schedule.csv"
    path.write_text(
        "time_s,rate_kg_s\n0,0\n1,0\n3,0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="positive mass"):
        field_source_schedule_from_csv(path, evidence=_evidence(path))


def test_schedule_rejects_missing_values_and_normalized_header_collisions(tmp_path):
    missing = tmp_path / "missing-value.csv"
    missing.write_text(
        "time_s,rate_kg_s\n0,0.2\n1\n3,0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="fewer values than header"):
        field_source_schedule_from_csv(missing, evidence=_evidence(missing))

    collision = tmp_path / "normalized-header.csv"
    collision.write_text(
        "time_s, rate_kg_s,rate_kg_s\n0,0.2,0.2\n1,0.1,0.1\n3,0,0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate columns"):
        field_source_schedule_from_csv(collision, evidence=_evidence(collision))


def test_schedule_rejects_header_only_csv(tmp_path):
    path = tmp_path / "header-only.csv"
    path.write_text("time_s,rate_kg_s\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no data rows"):
        field_source_schedule_from_csv(path, evidence=_evidence(path))


def test_schedule_rejects_changed_csv_digest(tmp_path):
    path = _csv(tmp_path)
    evidence = _evidence(path)
    path.write_text(path.read_text(encoding="utf-8") + "", encoding="utf-8")
    path.write_text("time_s,rate_kg_s\n0,0.3\n1,0.1\n3,0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256"):
        field_source_schedule_from_csv(path, evidence=evidence)


def test_schedule_rejects_ambiguous_csv_structure_and_column_roles(tmp_path):
    duplicate = tmp_path / "duplicate-header.csv"
    duplicate.write_text(
        "time_s,time_s,rate_kg_s\n0,0,0.2\n1,1,0.1\n3,3,0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate columns"):
        field_source_schedule_from_csv(duplicate, evidence=_evidence(duplicate))

    extra = tmp_path / "extra-value.csv"
    extra.write_text(
        "time_s,rate_kg_s\n0,0.2,unexpected\n1,0.1,unexpected\n3,0,unexpected\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="more values than header"):
        field_source_schedule_from_csv(extra, evidence=_evidence(extra))

    path = _csv(tmp_path)
    with pytest.raises(ValueError, match="column roles must be distinct"):
        field_source_schedule_from_csv(
            path, evidence=_evidence(path), time_column="rate_kg_s",
        )


def test_strict_schedule_case_resolves_relative_csv_and_cli(tmp_path, capsys):
    path = _csv(tmp_path)
    evidence = _evidence(path, source_kind="pool_vapour")
    case = {
        "schema": FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA,
        "csv_path": "schedule.csv",
        "evidence": {
            "dataset_id": evidence.dataset_id,
            "path": "schedule.csv",
            "sha256": evidence.sha256,
            "row_count": evidence.row_count,
            "source_boundary_id": evidence.source_boundary_id,
            "common_clock_id": evidence.common_clock_id,
            "source_kind": evidence.source_kind,
        },
    }
    case_path = tmp_path / "source-case.json"
    case_path.write_text(json.dumps(case), encoding="utf-8")
    parsed = field_source_schedule_case_from_mapping(case, base_directory=tmp_path)
    assert parsed.schedule.duration_s == 3.0
    assert main(["field-source", str(case_path)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["atmospheric_source_schedule"]["evidence"]["source_kind"] == "pool_vapour"


def test_strict_schedule_json_rejects_normalized_duplicate_keys(tmp_path):
    path = tmp_path / "duplicate-keys.json"
    path.write_text(
        '{"schema": "degali.field-atmospheric-schedule-input.v1", '
        '" schema ": "degali.field-atmospheric-schedule-input.v1"}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate or normalization-colliding"):
        read_field_source_schedule_json(path)


def test_unknown_source_kind_is_rejected(tmp_path):
    path = _csv(tmp_path)
    with pytest.raises(ValueError, match="source_kind"):
        _evidence(path, source_kind="liquid_guess")


def test_direct_rate_bound_schedule_requires_zero_endpoints_and_column_names(tmp_path):
    path = _csv(tmp_path)
    evidence = _evidence(path)
    nonzero_nominal = SourceRateSchedule((0.0, 1.0, 3.0), (0.2, 0.1, 0.01), source_id="source-a")
    with pytest.raises(ValueError, match="schedule final endpoint"):
        FieldAtmosphericSourceSchedule(
            nonzero_nominal, evidence, "time_s", "rate_kg_s",
        )
    nominal = SourceRateSchedule((0.0, 1.0, 3.0), (0.2, 0.1, 0.0), source_id="source-a")
    lower = SourceRateSchedule((0.0, 1.0, 3.0), (0.1, 0.05, 0.01), source_id="source-a")
    upper = SourceRateSchedule((0.0, 1.0, 3.0), (0.3, 0.2, 0.0), source_id="source-a")
    with pytest.raises(ValueError, match="lower_schedule final endpoint"):
        FieldAtmosphericSourceSchedule(
            nominal, evidence, "time_s", "rate_kg_s", lower, upper,
            "lower", "upper",
        )
    valid_lower = SourceRateSchedule((0.0, 1.0, 3.0), (0.1, 0.05, 0.0), source_id="source-a")
    with pytest.raises(ValueError, match="lower_rate_column"):
        FieldAtmosphericSourceSchedule(
            nominal, evidence, "time_s", "rate_kg_s", valid_lower, upper,
            "", "upper",
        )
    with pytest.raises(ValueError, match="supplied together"):
        FieldAtmosphericSourceSchedule(
            nominal, evidence, "time_s", "rate_kg_s", None, None,
            "lower", "upper",
        )
