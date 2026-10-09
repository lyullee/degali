"""Fail-safe CSV import for a selected measured LH2 source-history event.

The importer accepts only an already selected, SI-unit CSV interval. It does
not infer a release period from valve signals, convert a normal-volume flow to
mass flow, estimate a phase fraction from a level trace, or use an operating
pressure trend as a discharge rate. Those are engineering source-boundary
decisions which need their own evidence before this module is called.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, replace
import hashlib
import math
from pathlib import Path
from typing import TYPE_CHECKING

from .field_contracts import BoundedValue, ReleaseSource
from .field_history import (
    FieldMeasuredFlashSchedule,
    PressureDrivenMeasuredHistory,
    MeasuredHistoryQualityCriteria,
    MeasuredReleaseHistory,
    MeasuredTimeSeries,
    direct_vapour_schedule_from_pressure_driven_history,
    direct_vapour_schedule_from_measured_history,
)

if TYPE_CHECKING:
    from .lh2_property_table import LH2SaturationTable
@dataclass(frozen=True)
class HistorianCsvChannel:
    """One calibrated SI channel in a selected historian CSV export.

    Either paired lower/upper columns or one/two declared half-width bounds
    are required. The maximum of declared absolute and relative half-width is
    used so combining certificate limits cannot accidentally narrow a bound.
    """

    value_column: str
    unit: str
    source_id: str
    calibration_evidence_id: str
    response_time_s: float = 0.0
    time_offset_s: float = 0.0
    lower_column: str | None = None
    upper_column: str | None = None
    absolute_half_width: float | None = None
    relative_half_width: float | None = None

    def __post_init__(self) -> None:
        for name, value in {
            "value_column": self.value_column,
            "source_id": self.source_id,
            "calibration_evidence_id": self.calibration_evidence_id,
        }.items():
            if not isinstance(value, str) or not value.strip() or value == "unspecified":
                raise ValueError(f"{name} must be a declared non-empty string")
        if not isinstance(self.unit, str) or not self.unit.strip():
            raise ValueError("unit must be a non-empty string")
        paired_bounds = self.lower_column is not None or self.upper_column is not None
        if paired_bounds and (not self.lower_column or not self.upper_column):
            raise ValueError("lower_column and upper_column must be declared together")
        declared_width = (
            self.absolute_half_width is not None or self.relative_half_width is not None
        )
        if paired_bounds and declared_width:
            raise ValueError("use either paired bound columns or declared half-width bounds, not both")
        if not paired_bounds and not declared_width:
            raise ValueError("channel uncertainty must be declared as bound columns or half-width")
        for name, value in {
            "response_time_s": self.response_time_s,
            "time_offset_s": self.time_offset_s,
        }.items():
            if not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
        if self.response_time_s < 0.0:
            raise ValueError("response_time_s must be non-negative")
        for name, value in {
            "absolute_half_width": self.absolute_half_width,
            "relative_half_width": self.relative_half_width,
        }.items():
            if value is not None and (
                not math.isfinite(float(value)) or float(value) < 0.0
            ):
                raise ValueError(f"{name} must be non-negative and finite when supplied")

    def required_columns(self) -> tuple[str, ...]:
        return tuple(
            item for item in (self.value_column, self.lower_column, self.upper_column)
            if item is not None
        )

    def bounds(self, nominal: float, row: dict[str, str], *, row_number: int) -> tuple[float, float]:
        if self.lower_column is not None:
            lower = _required_float(row, self.lower_column, row_number)
            upper = _required_float(row, self.upper_column, row_number)
            if lower > upper:
                raise ValueError(
                    f"CSV row {row_number} {self.value_column!r} bounds must satisfy lower <= upper"
                )
            if not lower <= nominal <= upper:
                raise ValueError(
                    f"CSV row {row_number} {self.value_column!r} nominal value lies outside declared bounds"
                )
            return lower, upper
        half_widths = [
            0.0 if self.absolute_half_width is None else float(self.absolute_half_width),
            0.0 if self.relative_half_width is None else abs(nominal) * float(self.relative_half_width),
        ]
        half_width = max(half_widths)
        return nominal - half_width, nominal + half_width


@dataclass(frozen=True)
class MeasuredHistoryCsvMap:
    """Explicit schema and provenance for one selected event CSV.

    ``time_s_column`` is a monotonic numeric *relative-seconds* column. An
    absolute timestamp or spreadsheet serial must be normalised outside this
    reader with its clock-origin provenance retained in ``event_evidence_id``.
    This avoids locale-specific timestamp parsing and hidden time-zone shifts.
    """

    time_s_column: str
    event_id: str
    event_evidence_id: str
    phase_evidence_id: str
    pressure_pa: HistorianCsvChannel
    temperature_k: HistorianCsvChannel
    mass_flow_kg_s: HistorianCsvChannel
    liquid_fraction: HistorianCsvChannel | None = None

    def __post_init__(self) -> None:
        for name, value in {
            "time_s_column": self.time_s_column,
            "event_id": self.event_id,
            "event_evidence_id": self.event_evidence_id,
            "phase_evidence_id": self.phase_evidence_id,
        }.items():
            if not isinstance(value, str) or not value.strip() or value == "unspecified":
                raise ValueError(f"{name} must be a declared non-empty string")
        expected = {
            "pressure_pa": "Pa", "temperature_k": "K", "mass_flow_kg_s": "kg/s",
        }
        for name, unit in expected.items():
            channel = getattr(self, name)
            if not isinstance(channel, HistorianCsvChannel):
                raise TypeError(f"{name} must be a HistorianCsvChannel")
            if channel.unit != unit:
                raise ValueError(f"{name} must declare SI unit {unit!r}")
        if self.liquid_fraction is not None:
            if not isinstance(self.liquid_fraction, HistorianCsvChannel):
                raise TypeError("liquid_fraction must be a HistorianCsvChannel or None")
            if self.liquid_fraction.unit != "1":
                raise ValueError("liquid_fraction must declare SI unit '1'")

    def channels(self) -> dict[str, HistorianCsvChannel]:
        result = {
            "pressure_pa": self.pressure_pa,
            "temperature_k": self.temperature_k,
            "mass_flow_kg_s": self.mass_flow_kg_s,
        }
        if self.liquid_fraction is not None:
            result["liquid_fraction"] = self.liquid_fraction
        return result


@dataclass(frozen=True)
class PressureDrivenHistoryCsvMap:
    """Explicit schema for a pressure/temperature-only orifice history CSV."""

    time_s_column: str
    event_id: str
    event_evidence_id: str
    phase_evidence_id: str
    pressure_pa: HistorianCsvChannel
    temperature_k: HistorianCsvChannel
    liquid_fraction: HistorianCsvChannel | None = None

    def __post_init__(self) -> None:
        for name, value in {
            "time_s_column": self.time_s_column,
            "event_id": self.event_id,
            "event_evidence_id": self.event_evidence_id,
            "phase_evidence_id": self.phase_evidence_id,
        }.items():
            if not isinstance(value, str) or not value.strip() or value == "unspecified":
                raise ValueError(f"{name} must be a declared non-empty string")
        for name, unit in {"pressure_pa": "Pa", "temperature_k": "K"}.items():
            channel = getattr(self, name)
            if not isinstance(channel, HistorianCsvChannel):
                raise TypeError(f"{name} must be a HistorianCsvChannel")
            if channel.unit != unit:
                raise ValueError(f"{name} must declare SI unit {unit!r}")
        if self.liquid_fraction is not None:
            if not isinstance(self.liquid_fraction, HistorianCsvChannel):
                raise TypeError("liquid_fraction must be a HistorianCsvChannel or None")
            if self.liquid_fraction.unit != "1":
                raise ValueError("liquid_fraction must declare SI unit '1'")

    def channels(self) -> dict[str, HistorianCsvChannel]:
        result = {
            "pressure_pa": self.pressure_pa,
            "temperature_k": self.temperature_k,
        }
        if self.liquid_fraction is not None:
            result["liquid_fraction"] = self.liquid_fraction
        return result


@dataclass(frozen=True)
class MeasuredHistoryCsvProvenance:
    """Fingerprint and stated mapping for an imported source history."""

    path: str
    sha256: str
    row_count: int
    event_id: str
    event_evidence_id: str
    phase_evidence_id: str
    channel_evidence: tuple[tuple[str, str], ...]

    def __post_init__(self) -> None:
        for name, value in {
            "path": self.path,
            "event_id": self.event_id,
            "event_evidence_id": self.event_evidence_id,
            "phase_evidence_id": self.phase_evidence_id,
        }.items():
            if not isinstance(value, str) or not value.strip() or value == "unspecified":
                raise ValueError(f"CSV provenance {name} must be a declared non-empty string")
        if len(self.sha256) != 64 or any(value not in "0123456789abcdef" for value in self.sha256):
            raise ValueError("CSV provenance requires a lowercase SHA-256 digest")
        if not isinstance(self.row_count, int) or self.row_count < 2:
            raise ValueError("CSV provenance row_count must be at least two")
        evidence = tuple(self.channel_evidence)
        if any(
            not isinstance(name, str) or not name.strip()
            or not isinstance(value, str) or not value.strip() or value == "unspecified"
            for name, value in evidence
        ):
            raise ValueError("CSV channel evidence must contain declared non-empty string pairs")
        if len({name for name, _value in evidence}) != len(evidence):
            raise ValueError("CSV channel evidence names must be unique")
        object.__setattr__(self, "channel_evidence", evidence)


@dataclass(frozen=True)
class ImportedMeasuredReleaseHistory:
    """Measured history plus the artefact fingerprint required to audit it."""

    history: MeasuredReleaseHistory
    mapping: MeasuredHistoryCsvMap
    provenance: MeasuredHistoryCsvProvenance

    def __post_init__(self) -> None:
        if not isinstance(self.history, MeasuredReleaseHistory):
            raise TypeError("imported history must be a MeasuredReleaseHistory")
        if not isinstance(self.mapping, MeasuredHistoryCsvMap):
            raise TypeError("imported mapping must be a MeasuredHistoryCsvMap")
        if not isinstance(self.provenance, MeasuredHistoryCsvProvenance):
            raise TypeError("imported provenance must be a MeasuredHistoryCsvProvenance")


@dataclass(frozen=True)
class ImportedPressureDrivenMeasuredHistory:
    """Pressure-driven history plus the fingerprinted CSV mapping."""

    history: PressureDrivenMeasuredHistory
    mapping: PressureDrivenHistoryCsvMap
    provenance: MeasuredHistoryCsvProvenance

    def __post_init__(self) -> None:
        if not isinstance(self.history, PressureDrivenMeasuredHistory):
            raise TypeError(
                "imported pressure-driven history must be a PressureDrivenMeasuredHistory"
            )
        if not isinstance(self.mapping, PressureDrivenHistoryCsvMap):
            raise TypeError("pressure-driven mapping must be a PressureDrivenHistoryCsvMap")
        if not isinstance(self.provenance, MeasuredHistoryCsvProvenance):
            raise TypeError("imported provenance must be a MeasuredHistoryCsvProvenance")


def _import_provenance_pairs(
    imported: ImportedMeasuredReleaseHistory | ImportedPressureDrivenMeasuredHistory,
) -> tuple[tuple[str, str], ...]:
    """Flatten the controlled CSV mapping for the portable screening report."""
    provenance = imported.provenance
    pairs: list[tuple[str, str]] = [
        (
            "import_format",
            "degali.pressure-driven-history-csv.v1"
            if isinstance(imported, ImportedPressureDrivenMeasuredHistory)
            else "degali.measured-history-csv.v1",
        ),
        (
            "history_kind",
            "pressure_driven_orifice"
            if isinstance(imported, ImportedPressureDrivenMeasuredHistory)
            else "measured_flow",
        ),
        ("event_id", provenance.event_id),
        ("event_evidence_id", provenance.event_evidence_id),
        ("phase_evidence_id", provenance.phase_evidence_id),
        ("source_path", provenance.path),
        ("source_sha256", provenance.sha256),
        ("source_row_count", str(provenance.row_count)),
    ]
    for name, channel in imported.mapping.channels().items():
        prefix = f"channel.{name}"
        pairs.extend((
            (f"{prefix}.unit", channel.unit),
            (f"{prefix}.source_id", channel.source_id),
            (f"{prefix}.calibration_evidence_id", channel.calibration_evidence_id),
            (f"{prefix}.response_time_s", repr(float(channel.response_time_s))),
            (f"{prefix}.time_offset_s", repr(float(channel.time_offset_s))),
        ))
        if channel.lower_column is not None:
            pairs.extend((
                (f"{prefix}.lower_column", channel.lower_column),
                (f"{prefix}.upper_column", channel.upper_column or ""),
            ))
        else:
            pairs.extend((
                (f"{prefix}.absolute_half_width", repr(channel.absolute_half_width)),
                (f"{prefix}.relative_half_width", repr(channel.relative_half_width)),
            ))
    return tuple(pairs)


def _required_float(row: dict[str, str], column: str | None, row_number: int) -> float:
    assert column is not None
    raw = row.get(column)
    if raw is None or not raw.strip():
        raise ValueError(f"CSV row {row_number} has no value for required column {column!r}")
    try:
        value = float(raw)
    except ValueError as error:
        raise ValueError(f"CSV row {row_number} has invalid numeric value in {column!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"CSV row {row_number} has non-finite value in {column!r}")
    return value


def _series_from_rows(
    name: str,
    channel: HistorianCsvChannel,
    times: tuple[float, ...],
    rows: list[dict[str, str]],
) -> MeasuredTimeSeries:
    nominal = tuple(
        _required_float(row, channel.value_column, row_number)
        for row_number, row in enumerate(rows, start=2)
    )
    bounds = tuple(
        channel.bounds(value, row, row_number=row_number)
        for row_number, (value, row) in enumerate(zip(nominal, rows), start=2)
    )
    lower, upper = tuple(item[0] for item in bounds), tuple(item[1] for item in bounds)
    if name in {"pressure_pa", "temperature_k"} and min(lower) <= 0.0:
        raise ValueError(f"CSV {name} lower bound must stay positive")
    if name == "mass_flow_kg_s" and min(lower) < 0.0:
        raise ValueError("CSV mass_flow_kg_s lower bound must stay non-negative")
    if name == "liquid_fraction" and (min(lower) < 0.0 or max(upper) > 1.0):
        raise ValueError("CSV liquid_fraction bounds must lie in [0, 1]")
    return MeasuredTimeSeries(
        time_s=times, nominal=nominal, lower=lower, upper=upper,
        unit=channel.unit, source_id=channel.source_id,
        response_time_s=channel.response_time_s, time_offset_s=channel.time_offset_s,
    )


def _read_history_csv(
    path: str | Path,
    mapping: MeasuredHistoryCsvMap | PressureDrivenHistoryCsvMap,
) -> ImportedMeasuredReleaseHistory | ImportedPressureDrivenMeasuredHistory:
    """Read one selected SI history CSV without interpolation."""
    if not isinstance(mapping, (MeasuredHistoryCsvMap, PressureDrivenHistoryCsvMap)):
        raise TypeError(
            "mapping must be a MeasuredHistoryCsvMap or PressureDrivenHistoryCsvMap"
        )
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    with source.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        raw_headers = next(reader, None)
        headers = tuple(raw_headers or ())
        if not headers:
            raise ValueError("history CSV requires a non-empty header")
        if any(not header.strip() for header in headers):
            raise ValueError("history CSV header contains an empty column")
        duplicates = sorted(
            header for header, count in Counter(headers).items() if count > 1
        )
        if duplicates:
            raise ValueError(
                "history CSV header contains duplicate columns: "
                + ", ".join(duplicates)
            )
        rows: list[dict[str, str]] = []
        for row in reader:
            if not row:
                continue
            if len(row) != len(headers):
                raise ValueError(
                    "history CSV row has a different number of fields than its header"
                )
            rows.append(dict(zip(headers, row)))
        required = {mapping.time_s_column}
        required.update(
            column for channel in mapping.channels().values()
            for column in channel.required_columns()
        )
        missing = sorted(required - set(headers))
        if missing:
            raise ValueError("history CSV is missing required columns: " + ", ".join(missing))
    if len(rows) < 2:
        raise ValueError("history CSV requires at least two data rows")
    times = tuple(
        _required_float(row, mapping.time_s_column, row_number)
        for row_number, row in enumerate(rows, start=2)
    )
    if any(later <= earlier for earlier, later in zip(times, times[1:])):
        raise ValueError("history CSV time_s column must be strictly increasing")
    series = {
        name: _series_from_rows(name, channel, times, rows)
        for name, channel in mapping.channels().items()
    }
    if isinstance(mapping, PressureDrivenHistoryCsvMap):
        history: MeasuredReleaseHistory | PressureDrivenMeasuredHistory = PressureDrivenMeasuredHistory(
            pressure_pa=series["pressure_pa"], temperature_k=series["temperature_k"],
            liquid_fraction=series.get("liquid_fraction"),
            event_id=mapping.event_id, phase_evidence_id=mapping.phase_evidence_id,
        )
    else:
        history = MeasuredReleaseHistory(
            pressure_pa=series["pressure_pa"], temperature_k=series["temperature_k"],
            mass_flow_kg_s=series["mass_flow_kg_s"],
            liquid_fraction=series.get("liquid_fraction"),
        )
    provenance = MeasuredHistoryCsvProvenance(
        path=str(source.resolve()), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        row_count=len(rows), event_id=mapping.event_id,
        event_evidence_id=mapping.event_evidence_id,
        phase_evidence_id=mapping.phase_evidence_id,
        channel_evidence=tuple(
            (name, channel.calibration_evidence_id)
            for name, channel in mapping.channels().items()
        ),
    )
    if isinstance(mapping, PressureDrivenHistoryCsvMap):
        return ImportedPressureDrivenMeasuredHistory(history, mapping, provenance)
    return ImportedMeasuredReleaseHistory(history, mapping, provenance)


def read_measured_history_csv(
    path: str | Path,
    mapping: MeasuredHistoryCsvMap,
) -> ImportedMeasuredReleaseHistory:
    """Read one selected P/T/mass-flow SI history CSV without interpolation."""
    if not isinstance(mapping, MeasuredHistoryCsvMap):
        raise TypeError("mapping must be a MeasuredHistoryCsvMap")
    imported = _read_history_csv(path, mapping)
    assert isinstance(imported, ImportedMeasuredReleaseHistory)
    return imported


def read_pressure_driven_history_csv(
    path: str | Path,
    mapping: PressureDrivenHistoryCsvMap,
) -> ImportedPressureDrivenMeasuredHistory:
    """Read one selected P/T(/liquid-fraction) orifice history CSV."""
    if not isinstance(mapping, PressureDrivenHistoryCsvMap):
        raise TypeError("mapping must be a PressureDrivenHistoryCsvMap")
    imported = _read_history_csv(path, mapping)
    assert isinstance(imported, ImportedPressureDrivenMeasuredHistory)
    return imported


def direct_vapour_schedule_from_imported_history(
    release: ReleaseSource,
    imported: ImportedMeasuredReleaseHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria | None = None,
    property_table: "LH2SaturationTable | None" = None,
    auto_property_table: bool = True,
    auto_property_table_minimum_intervals: int = 128,
    table_nodes: int = 161,
    **flash_kwargs: object,
) -> FieldMeasuredFlashSchedule:
    """Flash an imported history and retain its file fingerprint in warnings.

    A supplied table is used as declared. For long selected historian windows,
    the default automatic path builds one bounded table over every declared
    temperature bound and reuses it across source intervals. Small windows
    retain the direct reference path because constructing a table would cost
    more than the saved property calls. Neither path extrapolates a table.
    """
    if not isinstance(imported, ImportedMeasuredReleaseHistory):
        raise TypeError("imported must be an ImportedMeasuredReleaseHistory")
    if not isinstance(auto_property_table, bool):
        raise TypeError("auto_property_table must be boolean")
    if (
        not isinstance(auto_property_table_minimum_intervals, int)
        or auto_property_table_minimum_intervals < 1
    ):
        raise ValueError("auto_property_table_minimum_intervals must be a positive integer")
    if not isinstance(table_nodes, int) or table_nodes < 4:
        raise ValueError("table_nodes must be an integer of at least four")
    table = property_table
    table_warning: tuple[str, ...] = ()
    auto_built = False
    interval_count = len(imported.history.release_time_s) - 1
    if table is None and auto_property_table and interval_count >= auto_property_table_minimum_intervals:
        try:
            from .field_lh2 import build_lh2_saturation_table_for_temperature_bounds

            temperatures = imported.history.temperature_k
            table = build_lh2_saturation_table_for_temperature_bounds(
                min(temperatures.lower), max(temperatures.upper),
                ambient_pressure_pa=float(flash_kwargs.get("ambient_pressure_pa", 101325.0)),
                nodes=table_nodes,
            )
            auto_built = True
            table_warning = (
                "long historian window uses one bounded LH2 saturation table across all declared temperature bounds",
            )
        except (ValueError, RuntimeError) as error:
            table_warning = (
                "LH2 saturation table was not constructed for the imported historian; "
                f"direct property path used: {error}",
            )
    schedule = direct_vapour_schedule_from_measured_history(
        release, imported.history, quality_criteria=quality_criteria,
        property_table=table, **flash_kwargs,
    )
    provenance = imported.provenance
    return FieldMeasuredFlashSchedule(
        schedule=schedule.schedule,
        history_duration_s=schedule.history_duration_s,
        selection=schedule.selection,
        total_measured_mass_kg=schedule.total_measured_mass_kg,
        direct_vapour_mass_kg=schedule.direct_vapour_mass_kg,
        unrouted_postflash_liquid_mass_kg=schedule.unrouted_postflash_liquid_mass_kg,
        alignment_method=schedule.alignment_method,
        warnings=schedule.warnings + table_warning + (
            "measured-history CSV import "
            f"event_id={provenance.event_id!r}, sha256={provenance.sha256}, "
            f"event_evidence_id={provenance.event_evidence_id!r}",
        ),
        quality_assessment=schedule.quality_assessment,
        provenance=_import_provenance_pairs(imported) + (
            ("property_table_used", str(table is not None).lower()),
            ("property_table_auto_built", str(auto_built).lower()),
            ("property_table_interval_count", str(interval_count)),
        ),
        source_uncertainty_resolved=schedule.source_uncertainty_resolved,
    )


def direct_vapour_schedule_from_imported_pressure_driven_history(
    release: ReleaseSource,
    imported: ImportedPressureDrivenMeasuredHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria | None = None,
    property_table: "LH2SaturationTable | None" = None,
    auto_property_table: bool = True,
    auto_property_table_minimum_intervals: int = 128,
    table_nodes: int = 161,
    **flash_kwargs: object,
) -> FieldMeasuredFlashSchedule:
    """Derive a pressure-driven schedule while retaining CSV provenance."""
    if not isinstance(imported, ImportedPressureDrivenMeasuredHistory):
        raise TypeError(
            "imported must be an ImportedPressureDrivenMeasuredHistory"
        )
    if not isinstance(auto_property_table, bool):
        raise TypeError("auto_property_table must be boolean")
    if (
        not isinstance(auto_property_table_minimum_intervals, int)
        or auto_property_table_minimum_intervals < 1
    ):
        raise ValueError("auto_property_table_minimum_intervals must be a positive integer")
    if not isinstance(table_nodes, int) or table_nodes < 4:
        raise ValueError("table_nodes must be an integer of at least four")
    table = property_table
    table_warning: tuple[str, ...] = ()
    auto_built = False
    interval_count = len(imported.history.release_time_s) - 1
    if table is None and auto_property_table and interval_count >= auto_property_table_minimum_intervals:
        try:
            from .field_lh2 import build_lh2_saturation_table_for_temperature_bounds

            temperatures = imported.history.temperature_k
            table = build_lh2_saturation_table_for_temperature_bounds(
                min(temperatures.lower), max(temperatures.upper),
                ambient_pressure_pa=float(flash_kwargs.get("ambient_pressure_pa", 101325.0)),
                nodes=table_nodes,
            )
            auto_built = True
            table_warning = (
                "long pressure-driven historian window uses one bounded LH2 saturation table across all declared temperature bounds",
            )
        except (ValueError, RuntimeError) as error:
            table_warning = (
                "LH2 saturation table was not constructed for the pressure-driven historian; "
                f"direct property path used: {error}",
            )
    source_for_nominal_schedule = release
    unresolved_orifice_bounds = (
        not release.opening_area_m2.is_exact
        or not release.discharge_coefficient.is_exact
    )
    if unresolved_orifice_bounds:
        source_for_nominal_schedule = replace(
            release,
            opening_area_m2=BoundedValue(
                release.opening_area_m2.nominal,
                unit=release.opening_area_m2.unit,
                source=f"{release.opening_area_m2.source}; nominal-case",
            ),
            discharge_coefficient=BoundedValue(
                release.discharge_coefficient.nominal,
                unit=release.discharge_coefficient.unit,
                source=f"{release.discharge_coefficient.source}; nominal-case",
            ),
        )
    schedule = direct_vapour_schedule_from_pressure_driven_history(
        source_for_nominal_schedule,
        imported.history,
        quality_criteria=quality_criteria,
        property_table=table,
        **flash_kwargs,
    )
    if unresolved_orifice_bounds:
        schedule = replace(
            schedule,
            warnings=schedule.warnings + (
                "pressure-driven opening-area/discharge-coefficient bounds are unresolved in the nominal case; use a joint envelope for coherent throat recomputation",
            ),
            provenance=schedule.provenance + (
                ("opening_area_selection", "nominal"),
                ("discharge_coefficient_selection", "nominal"),
            ),
            source_uncertainty_resolved=False,
        )
    provenance = imported.provenance
    return FieldMeasuredFlashSchedule(
        schedule=schedule.schedule,
        history_duration_s=schedule.history_duration_s,
        selection=schedule.selection,
        total_measured_mass_kg=schedule.total_measured_mass_kg,
        direct_vapour_mass_kg=schedule.direct_vapour_mass_kg,
        unrouted_postflash_liquid_mass_kg=schedule.unrouted_postflash_liquid_mass_kg,
        alignment_method=schedule.alignment_method,
        warnings=schedule.warnings + table_warning + (
            "pressure-driven historian CSV import "
            f"event_id={provenance.event_id!r}, sha256={provenance.sha256}, "
            f"event_evidence_id={provenance.event_evidence_id!r}",
        ),
        quality_assessment=schedule.quality_assessment,
        provenance=_import_provenance_pairs(imported) + (
            ("property_table_used", str(table is not None).lower()),
            ("property_table_auto_built", str(auto_built).lower()),
            ("property_table_interval_count", str(interval_count)),
        ),
        source_uncertainty_resolved=schedule.source_uncertainty_resolved,
    )


__all__ = [
    "HistorianCsvChannel", "MeasuredHistoryCsvMap", "PressureDrivenHistoryCsvMap",
    "MeasuredHistoryCsvProvenance", "ImportedMeasuredReleaseHistory",
    "ImportedPressureDrivenMeasuredHistory", "read_measured_history_csv",
    "read_pressure_driven_history_csv",
    "direct_vapour_schedule_from_imported_history",
    "direct_vapour_schedule_from_imported_pressure_driven_history",
]
