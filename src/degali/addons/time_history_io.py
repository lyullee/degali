"""Strict CSV import and provenance for source and wind histories."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
from typing import Sequence

from .transient_receptor import SourceHistory, WindHistory


TIME_HISTORY_CSV_INPUT_SCHEMA = "degali.source-weather-time-history-csv.v1"


def _string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _headers(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        raw = reader.fieldnames
        if raw is None or not raw:
            raise ValueError("history CSV must declare a header")
        if any(name is None or not str(name).strip() for name in raw):
            raise ValueError("history CSV headers must be non-empty")
        normalized = [str(name).strip() for name in raw]
        if len(set(normalized)) != len(normalized):
            raise ValueError("history CSV headers must be unique")
        rows: list[dict[str, str]] = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"history CSV row {row_number} has more values than headers")
            if any(value is None for value in row.values()):
                raise ValueError(f"history CSV row {row_number} has missing values")
            rows.append({str(key).strip(): str(value).strip() for key, value in row.items()})
    if len(rows) < 2:
        raise ValueError("history CSV must contain at least two data rows")
    return normalized, rows


def _required_float(row: dict[str, str], column: str, row_number: int) -> float:
    text = row.get(column, "")
    if not text:
        raise ValueError(f"history CSV row {row_number} column {column!r} is empty")
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError(
            f"history CSV row {row_number} column {column!r} is not numeric"
        ) from error
    if not math.isfinite(value):
        raise ValueError(f"history CSV row {row_number} column {column!r} is not finite")
    return value


def _require_columns(headers: Sequence[str], columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in headers]
    if missing:
        raise ValueError("history CSV is missing columns: " + ", ".join(missing))


@dataclass(frozen=True)
class SourceHistoryCsvMap:
    """Explicit SI column map for an atmospheric source history."""

    time_s_column: str
    mass_rate_kg_s_column: str
    event_id: str
    evidence_id: str
    pressure_pa_column: str | None = None
    temperature_k_column: str | None = None
    source_direction_to_deg_column: str | None = None
    rate_operator: str = "piecewise_constant"

    def __post_init__(self) -> None:
        for name, value in {
            "time_s_column": self.time_s_column,
            "mass_rate_kg_s_column": self.mass_rate_kg_s_column,
            "event_id": self.event_id,
            "evidence_id": self.evidence_id,
        }.items():
            _string(value, name)
        optional = (
            self.pressure_pa_column,
            self.temperature_k_column,
            self.source_direction_to_deg_column,
        )
        if (self.pressure_pa_column is None) != (self.temperature_k_column is None):
            raise ValueError("pressure and temperature columns must be supplied together")
        if any(value is not None and not str(value).strip() for value in optional):
            raise ValueError("optional source history columns must be non-empty")
        if self.rate_operator not in {"piecewise_constant", "linear"}:
            raise ValueError("rate_operator must be piecewise_constant or linear")


@dataclass(frozen=True)
class WindHistoryCsvMap:
    """Explicit meteorological wind column map."""

    time_s_column: str
    speed_m_s_column: str
    direction_from_deg_column: str
    event_id: str
    evidence_id: str

    def __post_init__(self) -> None:
        for name, value in {
            "time_s_column": self.time_s_column,
            "speed_m_s_column": self.speed_m_s_column,
            "direction_from_deg_column": self.direction_from_deg_column,
            "event_id": self.event_id,
            "evidence_id": self.evidence_id,
        }.items():
            _string(value, name)


@dataclass(frozen=True)
class ImportedSourceHistory:
    history: SourceHistory
    path: str
    sha256: str
    row_count: int
    event_id: str
    evidence_id: str

    def as_record(self) -> dict[str, object]:
        return {
            "kind": "source_history",
            "path": self.path,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "event_id": self.event_id,
            "evidence_id": self.evidence_id,
            "released_mass_kg": self.history.released_mass_kg(),
        }


@dataclass(frozen=True)
class ImportedWindHistory:
    history: WindHistory
    path: str
    sha256: str
    row_count: int
    event_id: str
    evidence_id: str

    def as_record(self) -> dict[str, object]:
        time, speed, direction = self.history.arrays()
        return {
            "kind": "wind_history",
            "path": self.path,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "event_id": self.event_id,
            "evidence_id": self.evidence_id,
            "time_range_s": [float(time[0]), float(time[-1])],
            "speed_range_m_s": [float(speed.min()), float(speed.max())],
            "direction_range_from_deg": [float(direction.min()), float(direction.max())],
        }


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_source_history_csv(
    path: str | Path,
    mapping: SourceHistoryCsvMap,
) -> ImportedSourceHistory:
    if not isinstance(mapping, SourceHistoryCsvMap):
        raise TypeError("mapping must be a SourceHistoryCsvMap")
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    headers, rows = _headers(source)
    optional = [
        column for column in (
            mapping.pressure_pa_column,
            mapping.temperature_k_column,
            mapping.source_direction_to_deg_column,
        ) if column is not None
    ]
    _require_columns(headers, [mapping.time_s_column, mapping.mass_rate_kg_s_column, *optional])
    time: list[float] = []
    rate: list[float] = []
    pressure: list[float] | None = [] if mapping.pressure_pa_column else None
    temperature: list[float] | None = [] if mapping.temperature_k_column else None
    direction: list[float] | None = [] if mapping.source_direction_to_deg_column else None
    for row_number, row in enumerate(rows, start=2):
        time.append(_required_float(row, mapping.time_s_column, row_number))
        rate.append(_required_float(row, mapping.mass_rate_kg_s_column, row_number))
        if pressure is not None:
            pressure.append(_required_float(row, mapping.pressure_pa_column, row_number))
            temperature.append(_required_float(row, mapping.temperature_k_column, row_number))
        if direction is not None:
            direction.append(_required_float(row, mapping.source_direction_to_deg_column, row_number))
    history = SourceHistory(
        time_s=tuple(time), mass_rate_kg_s=tuple(rate),
        pressure_pa=None if pressure is None else tuple(pressure),
        temperature_k=None if temperature is None else tuple(temperature),
        source_direction_to_deg=None if direction is None else tuple(direction),
        rate_operator=mapping.rate_operator,
    )
    history.arrays()
    return ImportedSourceHistory(
        history=history, path=str(source), sha256=_digest(source), row_count=len(rows),
        event_id=mapping.event_id, evidence_id=mapping.evidence_id,
    )


def read_wind_history_csv(
    path: str | Path,
    mapping: WindHistoryCsvMap,
) -> ImportedWindHistory:
    if not isinstance(mapping, WindHistoryCsvMap):
        raise TypeError("mapping must be a WindHistoryCsvMap")
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    headers, rows = _headers(source)
    _require_columns(headers, [
        mapping.time_s_column, mapping.speed_m_s_column,
        mapping.direction_from_deg_column,
    ])
    time = []
    speed = []
    direction = []
    for row_number, row in enumerate(rows, start=2):
        time.append(_required_float(row, mapping.time_s_column, row_number))
        speed.append(_required_float(row, mapping.speed_m_s_column, row_number))
        direction.append(_required_float(row, mapping.direction_from_deg_column, row_number))
    history = WindHistory(tuple(time), tuple(speed), tuple(direction))
    history.arrays()
    return ImportedWindHistory(
        history=history, path=str(source), sha256=_digest(source), row_count=len(rows),
        event_id=mapping.event_id, evidence_id=mapping.evidence_id,
    )


__all__ = [
    "TIME_HISTORY_CSV_INPUT_SCHEMA",
    "SourceHistoryCsvMap", "WindHistoryCsvMap",
    "ImportedSourceHistory", "ImportedWindHistory",
    "read_source_history_csv", "read_wind_history_csv",
]
