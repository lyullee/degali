"""Strict import boundary for already-atmospheric field source schedules.

This adapter is intentionally downstream of phase physics.  It accepts a
piecewise-constant hydrogen-vapour rate that has already been established by a
source/flash or pool ledger, fingerprints the CSV, and preserves its common
clock and source-boundary identifiers.  It never converts liquid flow,
pressure, level, or a filename into an atmospheric leak rate.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import math
import json
import os
from collections import Counter
from pathlib import Path
from typing import Literal, Mapping

from .semi_fv_obstacle import SourceRateSchedule
from .field_json import strict_json_loads


FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA = "degali.field-atmospheric-schedule-input.v1"
FIELD_SOURCE_SCHEDULE_EXECUTION_SCHEMA = "degali.field-atmospheric-schedule-execution.v1"
_SOURCE_KINDS = {
    "declared_atmospheric_vapour",
    "post_flash_atmospheric_vapour",
    "pool_vapour",
    "droplet_evaporation",
}


def _string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


@dataclass(frozen=True)
class FieldSourceScheduleEvidence:
    """Evidence identity for one already-atmospheric source schedule CSV."""

    dataset_id: str
    path: str
    sha256: str
    row_count: int
    source_boundary_id: str
    common_clock_id: str
    source_kind: Literal[
        "declared_atmospheric_vapour",
        "post_flash_atmospheric_vapour", "pool_vapour", "droplet_evaporation"
    ]

    def __post_init__(self) -> None:
        for name, value in {
            "dataset_id": self.dataset_id,
            "path": self.path,
            "source_boundary_id": self.source_boundary_id,
            "common_clock_id": self.common_clock_id,
        }.items():
            if not isinstance(value, str) or not value.strip() or value.strip().lower() == "unspecified":
                raise ValueError(f"{name} must be explicitly declared")
        if (
            not isinstance(self.sha256, str)
            or len(self.sha256) != 64
            or any(char not in "0123456789abcdef" for char in self.sha256)
        ):
            raise ValueError("sha256 must be a lowercase SHA-256 digest")
        if not isinstance(self.row_count, int) or isinstance(self.row_count, bool) or self.row_count < 2:
            raise ValueError("row_count must be an integer of at least two")
        if self.source_kind not in _SOURCE_KINDS:
            raise ValueError(
                "source_kind must be one of " + ", ".join(sorted(_SOURCE_KINDS))
            )

    def as_record(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "path": self.path,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "source_boundary_id": self.source_boundary_id,
            "common_clock_id": self.common_clock_id,
            "source_kind": self.source_kind,
        }


@dataclass(frozen=True)
class FieldAtmosphericSourceSchedule:
    """A fingerprinted atmospheric source schedule ready for field transport."""

    schedule: SourceRateSchedule
    evidence: FieldSourceScheduleEvidence
    time_column: str
    rate_column: str
    lower_schedule: SourceRateSchedule | None = None
    upper_schedule: SourceRateSchedule | None = None
    lower_rate_column: str | None = None
    upper_rate_column: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.schedule, SourceRateSchedule):
            raise TypeError("schedule must be a SourceRateSchedule")
        if not isinstance(self.evidence, FieldSourceScheduleEvidence):
            raise TypeError("evidence must be FieldSourceScheduleEvidence")
        if self.schedule.source_id != self.evidence.source_boundary_id:
            raise ValueError("schedule source_id must match evidence source_boundary_id")
        if self.schedule.released_mass_kg <= 0.0:
            raise ValueError("atmospheric source schedule must release positive mass")
        self.schedule.require_zero_endpoint("schedule")
        if (self.lower_schedule is None) != (self.upper_schedule is None):
            raise ValueError("lower_schedule and upper_schedule must be supplied together")
        if (self.lower_rate_column is None) != (self.upper_rate_column is None):
            raise ValueError("lower_rate_column and upper_rate_column must be supplied together")
        if (self.lower_schedule is None) != (self.lower_rate_column is None):
            raise ValueError(
                "lower/upper schedule objects and rate columns must be supplied together"
            )
        if self.lower_schedule is not None and self.lower_rate_column is None:
            raise ValueError("rate-bound schedules require their CSV column names")
        if self.lower_schedule is not None and self.upper_schedule is not None:
            for name, corner in (("lower_schedule", self.lower_schedule), ("upper_schedule", self.upper_schedule)):
                if not isinstance(corner, SourceRateSchedule):
                    raise TypeError(f"{name} must be a SourceRateSchedule")
                if (
                    corner.source_id != self.schedule.source_id
                    or corner.time_s != self.schedule.time_s
                    or corner.rate_operator != self.schedule.rate_operator
                ):
                    raise ValueError(
                        f"{name} must use the nominal schedule source, time axis, and rate operator"
                    )
                corner.require_zero_endpoint(name)
            for lower, nominal, upper in zip(
                self.lower_schedule.rate_kg_s,
                self.schedule.rate_kg_s,
                self.upper_schedule.rate_kg_s,
            ):
                if not lower <= nominal <= upper:
                    raise ValueError("source schedule rate bounds must contain nominal rates")
        column_names = [self.time_column, self.rate_column]
        column_names.extend(
            value for value in (self.lower_rate_column, self.upper_rate_column)
            if value is not None
        )
        for name, value in {
            "time_column": self.time_column, "rate_column": self.rate_column,
            "lower_rate_column": self.lower_rate_column,
            "upper_rate_column": self.upper_rate_column,
        }.items():
            if value is not None and (
                not isinstance(value, str)
                or not value.strip()
                or any(character in value for character in "\r\n")
            ):
                raise ValueError(f"{name} must be a non-empty single-line string")
        duplicates = sorted(
            name for name, count in Counter(column_names).items() if count > 1
        )
        if duplicates:
            raise ValueError(
                "source schedule column roles must be distinct: "
                + ", ".join(duplicates)
            )

    def as_record(self) -> dict[str, object]:
        uncertainty = None
        if self.lower_schedule is not None and self.upper_schedule is not None:
            uncertainty = {
                "lower_rate_kg_s": list(self.lower_schedule.rate_kg_s),
                "upper_rate_kg_s": list(self.upper_schedule.rate_kg_s),
                "interpretation": "deterministic source-rate corners; not a probability interval",
            }
        return {
            "schema": "degali.field-atmospheric-schedule.v1",
            "evidence": self.evidence.as_record(),
            "time_column": self.time_column,
            "rate_column": self.rate_column,
            "lower_rate_column": self.lower_rate_column,
            "upper_rate_column": self.upper_rate_column,
            "schedule": {
                "source_id": self.schedule.source_id,
                "rate_operator": self.schedule.rate_operator,
                "time_s": list(self.schedule.time_s),
                "rate_kg_s": list(self.schedule.rate_kg_s),
                "duration_s": self.schedule.duration_s,
                "released_mass_kg": self.schedule.released_mass_kg,
            },
            "uncertainty": uncertainty,
            "scope": (
                "already-atmospheric source schedule only; no liquid/phase inference, "
                "hazard-distance interpolation, or design-basis claim"
            ),
        }

    def corner_schedules(self) -> tuple[tuple[str, SourceRateSchedule], ...]:
        """Return deterministic lower/nominal/upper source-rate corners."""
        if self.lower_schedule is None or self.upper_schedule is None:
            return (("nominal", self.schedule),)
        return (
            ("lower", self.lower_schedule),
            ("nominal", self.schedule),
            ("upper", self.upper_schedule),
        )


def field_source_schedule_from_csv(
    path: str | Path,
    *,
    evidence: FieldSourceScheduleEvidence,
    time_column: str = "time_s",
    rate_column: str = "rate_kg_s",
    lower_rate_column: str = "rate_lower_kg_s",
    upper_rate_column: str = "rate_upper_kg_s",
    rate_operator: str = "piecewise_constant",
) -> FieldAtmosphericSourceSchedule:
    """Read and fingerprint a monotone atmospheric source-rate schedule.

    Each row is an interval-start/end-point record for the default operator;
    with ``rate_operator='linear'`` the rows are time-node values instead. The
    final row is required to carry a zero rate so that the schedule cannot
    silently continue beyond its declared duration.
    """
    if not isinstance(evidence, FieldSourceScheduleEvidence):
        raise TypeError("evidence must be FieldSourceScheduleEvidence")
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    if Path(evidence.path).resolve() != source.resolve():
        raise ValueError("source schedule evidence path does not match the CSV path")
    column_roles = [time_column, rate_column, lower_rate_column, upper_rate_column]
    if any(
        not isinstance(name, str) or not name.strip()
        or any(character in name for character in "\r\n")
        for name in column_roles
    ):
        raise ValueError("source schedule column names must be non-empty and single-line")
    duplicates = sorted(
        name for name, count in Counter(column_roles).items() if count > 1
    )
    if duplicates:
        raise ValueError(
            "source schedule column roles must be distinct: "
            + ", ".join(duplicates)
        )
    with source.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        raw_fieldnames = next(reader, None)
        fieldnames = tuple(raw_fieldnames or ())
        if not fieldnames:
            raise ValueError("source schedule CSV must contain a header row")
        if any(not name.strip() for name in fieldnames):
            raise ValueError("source schedule CSV header contains an empty column name")
        normalized = tuple(name.strip() for name in fieldnames)
        header_duplicates = sorted(
            name for name, count in Counter(normalized).items() if count > 1
        )
        if header_duplicates:
            raise ValueError(
                "source schedule CSV header contains duplicate columns: "
                + ", ".join(header_duplicates)
            )
        required = {time_column, rate_column}
        missing = sorted(required - set(fieldnames))
        if missing:
            raise ValueError("source schedule CSV is missing required columns: " + ", ".join(missing))
        rows: list[dict[str, str]] = []
        for row_number, values in enumerate(reader, start=2):
            if len(values) != len(fieldnames):
                detail = (
                    "more values than header columns"
                    if len(values) > len(fieldnames)
                    else "fewer values than header columns"
                )
                raise ValueError(
                    f"invalid source schedule CSV row {row_number}: {detail}; "
                    f"expected {len(fieldnames)} values, got {len(values)}"
                )
            rows.append(dict(zip(fieldnames, values)))
    if not rows:
        raise ValueError("source schedule CSV contains no data rows")
    if len(rows) != evidence.row_count:
        raise ValueError(
            f"source schedule row_count {len(rows)} does not match evidence {evidence.row_count}"
        )
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != evidence.sha256:
        raise ValueError("source schedule CSV SHA-256 does not match evidence")
    has_lower = lower_rate_column in fieldnames
    has_upper = upper_rate_column in fieldnames
    if has_lower != has_upper:
        raise ValueError("source schedule lower/upper rate columns must be supplied together")
    times: list[float] = []
    rates: list[float] = []
    lower_rates: list[float] = []
    upper_rates: list[float] = []
    for row_number, row in enumerate(rows, start=2):
        if None in row:
            raise ValueError(
                f"invalid source schedule CSV row {row_number}: more values than header columns"
            )
        try:
            time = float(row[time_column])
            rate = float(row[rate_column])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"invalid source schedule CSV row {row_number}") from error
        if not math.isfinite(time) or not math.isfinite(rate):
            raise ValueError(f"source schedule row {row_number} contains non-finite values")
        if rate < 0.0:
            raise ValueError(f"source schedule row {row_number} has a negative rate")
        times.append(time)
        rates.append(rate)
        if has_lower and has_upper:
            try:
                lower = float(row[lower_rate_column])
                upper = float(row[upper_rate_column])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"invalid source schedule bounds at row {row_number}") from error
            if not all(math.isfinite(value) for value in (lower, upper)) or lower < 0.0 or upper < 0.0:
                raise ValueError(f"source schedule bounds at row {row_number} are invalid")
            if not lower <= rate <= upper:
                raise ValueError(f"source schedule bounds at row {row_number} do not contain nominal rate")
            lower_rates.append(lower)
            upper_rates.append(upper)
    if times[0] != 0.0:
        raise ValueError("source schedule must start at time zero")
    if any(later <= earlier for earlier, later in zip(times, times[1:])):
        raise ValueError("source schedule times must be strictly increasing")
    if rates[-1] != 0.0:
        raise ValueError("source schedule final endpoint rate must be zero")
    schedule = SourceRateSchedule(
        tuple(times), tuple(rates), source_id=evidence.source_boundary_id,
        rate_operator=rate_operator,
    )
    lower_schedule = upper_schedule = None
    if has_lower and has_upper:
        if lower_rates[-1] != 0.0 or upper_rates[-1] != 0.0:
            raise ValueError("source schedule lower/upper final endpoint rates must be zero")
        lower_schedule = SourceRateSchedule(
            tuple(times), tuple(lower_rates), source_id=evidence.source_boundary_id,
            rate_operator=rate_operator,
        )
        upper_schedule = SourceRateSchedule(
            tuple(times), tuple(upper_rates), source_id=evidence.source_boundary_id,
            rate_operator=rate_operator,
        )
    return FieldAtmosphericSourceSchedule(
        schedule, evidence, time_column, rate_column, lower_schedule, upper_schedule,
        lower_rate_column if has_lower else None,
        upper_rate_column if has_upper else None,
    )


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def field_source_schedule_case_from_mapping(
    value: object,
    *,
    base_directory: str | Path | None = None,
) -> FieldAtmosphericSourceSchedule:
    """Parse the strict JSON wrapper for one atmospheric schedule CSV."""
    data = _mapping(value, "field atmospheric schedule input")
    allowed = {
        "schema", "csv_path", "time_column", "rate_column",
        "lower_rate_column", "upper_rate_column", "rate_operator", "evidence",
    }
    unknown = set(data) - allowed
    missing = {"schema", "csv_path", "evidence"} - set(data)
    if unknown:
        raise ValueError("field atmospheric schedule input has unknown keys: " + ", ".join(sorted(unknown)))
    if missing:
        raise ValueError("field atmospheric schedule input is missing: " + ", ".join(sorted(missing)))
    if data["schema"] != FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA:
        raise ValueError(f"schema must be {FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA!r}")
    csv_path = Path(_string(data["csv_path"], "csv_path"))
    if not csv_path.is_absolute():
        if base_directory is None:
            raise ValueError("relative csv_path requires base_directory")
        csv_path = Path(base_directory) / csv_path
    evidence_data = _mapping(data["evidence"], "evidence")
    required = {
        "dataset_id", "path", "sha256", "row_count", "source_boundary_id",
        "common_clock_id", "source_kind",
    }
    unknown_evidence = set(evidence_data) - required
    missing_evidence = required - set(evidence_data)
    if unknown_evidence:
        raise ValueError("evidence has unknown keys: " + ", ".join(sorted(unknown_evidence)))
    if missing_evidence:
        raise ValueError("evidence is missing: " + ", ".join(sorted(missing_evidence)))
    evidence_path = Path(_string(evidence_data["path"], "evidence.path"))
    if not evidence_path.is_absolute():
        evidence_path = csv_path.parent / evidence_path
    evidence = FieldSourceScheduleEvidence(
        dataset_id=_string(evidence_data["dataset_id"], "evidence.dataset_id"),
        path=str(evidence_path.resolve()),
        sha256=_string(evidence_data["sha256"], "evidence.sha256"),
        row_count=evidence_data["row_count"] if isinstance(evidence_data["row_count"], int) else -1,
        source_boundary_id=_string(evidence_data["source_boundary_id"], "evidence.source_boundary_id"),
        common_clock_id=_string(evidence_data["common_clock_id"], "evidence.common_clock_id"),
        source_kind=_string(evidence_data["source_kind"], "evidence.source_kind"),
    )
    time_column = _string(data.get("time_column", "time_s"), "time_column")
    rate_column = _string(data.get("rate_column", "rate_kg_s"), "rate_column")
    lower_rate_column = _string(
        data.get("lower_rate_column", "rate_lower_kg_s"), "lower_rate_column",
    )
    upper_rate_column = _string(
        data.get("upper_rate_column", "rate_upper_kg_s"), "upper_rate_column",
    )
    rate_operator = _string(
        data.get("rate_operator", "piecewise_constant"), "rate_operator",
    )
    return field_source_schedule_from_csv(
        csv_path, evidence=evidence, time_column=time_column, rate_column=rate_column,
        lower_rate_column=lower_rate_column, upper_rate_column=upper_rate_column,
        rate_operator=rate_operator,
    )


def read_field_source_schedule_json(path: str | Path) -> FieldAtmosphericSourceSchedule:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    try:
        value = strict_json_loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid atmospheric schedule JSON at {source}: {error.msg}") from error
    return field_source_schedule_case_from_mapping(value, base_directory=source.parent)


def write_field_source_schedule_json(
    value: FieldAtmosphericSourceSchedule,
    path: str | Path,
) -> FieldAtmosphericSourceSchedule:
    """Write and re-read a strict source-schedule case with fresh provenance."""
    if not isinstance(value, FieldAtmosphericSourceSchedule):
        raise TypeError("value must be a FieldAtmosphericSourceSchedule")
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(
            f"refusing to overwrite existing atmospheric schedule case: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(value.evidence.path).resolve()
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)
    observed_digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if observed_digest != value.evidence.sha256:
        raise ValueError("source schedule CSV changed since its provenance was recorded")
    lower_column = value.lower_rate_column or "rate_lower_kg_s"
    upper_column = value.upper_rate_column or "rate_upper_kg_s"
    verified = field_source_schedule_from_csv(
        csv_path,
        evidence=value.evidence,
        time_column=value.time_column,
        rate_column=value.rate_column,
        lower_rate_column=lower_column,
        upper_rate_column=upper_column,
        rate_operator=value.schedule.rate_operator,
    )
    if verified != value:
        raise ValueError("source schedule CSV does not match the typed schedule")
    record = {
        "schema": FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA,
        "csv_path": Path(os.path.relpath(csv_path, destination.parent.resolve())).as_posix(),
        "time_column": value.time_column,
        "rate_column": value.rate_column,
        "lower_rate_column": lower_column,
        "upper_rate_column": upper_column,
        "rate_operator": value.schedule.rate_operator,
        "evidence": {
            "dataset_id": value.evidence.dataset_id,
            "path": csv_path.name,
            "sha256": observed_digest,
            "row_count": value.evidence.row_count,
            "source_boundary_id": value.evidence.source_boundary_id,
            "common_clock_id": value.evidence.common_clock_id,
            "source_kind": value.evidence.source_kind,
        },
    }
    payload = json.dumps(
        record, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False,
    ) + "\n"
    try:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing atmospheric schedule case: {destination}"
        ) from error
    return read_field_source_schedule_json(destination)


def field_source_schedule_record(value: FieldAtmosphericSourceSchedule) -> dict[str, object]:
    if not isinstance(value, FieldAtmosphericSourceSchedule):
        raise TypeError("value must be a FieldAtmosphericSourceSchedule")
    return value.as_record()


__all__ = [
    "FIELD_SOURCE_SCHEDULE_INPUT_SCHEMA", "FIELD_SOURCE_SCHEDULE_EXECUTION_SCHEMA",
    "FieldSourceScheduleEvidence", "FieldAtmosphericSourceSchedule",
    "field_source_schedule_from_csv", "field_source_schedule_case_from_mapping",
    "read_field_source_schedule_json", "write_field_source_schedule_json",
    "field_source_schedule_record",
]
