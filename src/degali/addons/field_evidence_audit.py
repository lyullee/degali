"""Read-only readiness audit for site validation evidence.

The audit deliberately stops before constructing :class:`FieldValidationEvidence`.
It inventories structured CSV/JSON/XLSX files and reports whether the five evidence
channels needed for an obstacle-transport comparison are visible:
source boundary, weather, obstacle geometry, fixed-receptor observations and a
common clock.  A matching filename or a pool-radius result is never enough to
promote a dataset.  The result is therefore a triage report, not validation.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Literal, Mapping
import zipfile


FIELD_EVIDENCE_AUDIT_SCHEMA = "degali.field-evidence-audit.v1"
FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA = "degali.field-evidence-audit-execution.v1"
FIELD_EVIDENCE_AUDIT_GATE_CODES = frozenset({
    "scan_incomplete",
    "channels_missing",
    "channels_partial",
    "coherent_candidate_unqualified",
    "cross_file_candidate_unqualified",
    "promotion_not_allowed",
})

_CHANNELS = (
    "source_boundary",
    "weather",
    "obstacle_geometry",
    "receptor_observations",
    "common_clock",
)
# Historical artifacts used the narrower aliases below.  Keep those groups
# available for strict read-back while allowing the current inventory to
# recognise common exports that explicitly label a source boundary, reference
# wind speed, or a measured mean concentration.
_LEGACY_SOURCE_BOUNDARY_FIELD_GROUPS = (
    ("source_boundary_id", "source_id", "release_id"),
    (
        "mass_flow", "mass_flow_kg_s", "release_rate",
        "release_rate_kg_s", "leak_rate", "source_rate",
        "source_rate_kg_s", "mass_rate_kg_s", "mass_rate_reported_kg_s",
        "rate_kg_s", "inflow_kg_s", "erate",
    ),
)
_LEGACY_WEATHER_FIELD_GROUPS = (
    ("weather_id", "wind_id", "meteorology_id"),
    (
        "wind_speed", "wind_speed_m_s", "wind_m_s",
        "mean_wind_speed_m_s", "wind_speed_high_m_s",
        "wind_speed_low_m_s", "wind_u", "wind_v",
    ),
)
_LEGACY_RECEPTOR_FIELD_GROUPS = (
    ("sensor", "sensor_id", "receptor", "receptor_id"),
    ("x_downwind_m", "x_m", "x"),
    ("y_crosswind_m", "y_m", "y"),
    ("height_m", "z_m", "z", "elevation_m"),
    (
        "observed_mole_fraction", "mole_fraction", "h2_mole_fraction",
        "concentration", "h2_concentration", "volume_percent",
        "observed_volume_percent", "observed_vol_pct",
        "observed_mean_vol_pct", "observed_time_mean_vol_pct",
    ),
)
_LEGACY_COMMON_CLOCK_FIELD_GROUPS = (
    ("common_clock_id", "clock_id", "synchronized_clock_id"),
)
_LEGACY_COLLECTION_FIELD_GROUPS = {
    "source_boundary": _LEGACY_SOURCE_BOUNDARY_FIELD_GROUPS,
    "weather": _LEGACY_WEATHER_FIELD_GROUPS,
    "receptor_observations": _LEGACY_RECEPTOR_FIELD_GROUPS,
    "common_clock": _LEGACY_COMMON_CLOCK_FIELD_GROUPS,
}
_COMMON_CLOCK_FIELD_GROUPS = (
    (
        "common_clock_id", "clock_id", "synchronized_clock_id",
        "synchronization_id", "sync_id", "time_sync_id", "clock_sync_id",
    ),
)
# Keep the collection contract machine-readable without weakening the
# conservative channel detector below.  Each inner tuple is an OR-group of
# accepted normalized field aliases; groups are AND-ed within one channel.
_CHANNEL_COLLECTION_SPECS = {
    "source_boundary": {
        "required_field_groups": (
            (
                "source_boundary_id", "source_id", "release_id",
                "source_boundary", "source_bound",
            ),
            (
                "mass_flow", "mass_flow_kg_s", "release_rate",
                "release_rate_kg_s", "leak_rate", "source_rate",
                "source_rate_kg_s", "mass_rate_kg_s", "mass_rate_reported_kg_s",
                "rate_kg_s", "inflow_kg_s", "erate",
            ),
        ),
        "qualification": (
            "Join a non-negative physical source rate to the same event and "
            "phase/source-boundary identity used by the receptor record."
        ),
        "missing_action": (
            "Export the event's atmospheric source-boundary ID and non-negative "
            "mass-rate history with units and phase/source provenance."
        ),
    },
    "weather": {
        "required_field_groups": (
            ("weather_id", "wind_id", "meteorology_id"),
            (
                "wind_speed", "wind_speed_m_s", "wind_m_s",
                "mean_wind_speed_m_s", "wind_speed_high_m_s",
                "wind_speed_low_m_s", "wind_u", "wind_v",
                "wind_ref_m_s", "u_ref_m_s", "wind_10m_m_s",
            ),
        ),
        "qualification": (
            "Match wind speed/direction (and any stability record) to the same "
            "event interval and common clock as the source and receptors."
        ),
        "missing_action": (
            "Export event-linked wind speed/direction, stability metadata and "
            "the measurement interval with an explicit weather ID."
        ),
    },
    "obstacle_geometry": {
        "required_field_groups": (
            ("obstacle_geometry_id", "obstacle_id", "geometry_id", "geometry", "obstacle", "barrier", "building", "wall"),
            ("height", "width", "length", "polygon", "vertices", "position", "geometry"),
        ),
        "qualification": (
            "Reconcile obstacle geometry, coordinate reference and the model's "
            "wind-plane representation for the same event."
        ),
        "missing_action": (
            "Export the event obstacle geometry or geometry ID with dimensions, "
            "coordinates, vertical datum and coordinate-reference evidence."
        ),
    },
    "receptor_observations": {
        "required_field_groups": (
            ("sensor", "sensor_id", "receptor", "receptor_id"),
            ("x_downwind_m", "x_m", "x"),
            ("y_crosswind_m", "y_m", "y"),
            ("height_m", "z_m", "z", "elevation_m"),
            (
                "observed_mole_fraction", "mole_fraction", "h2_mole_fraction",
                "concentration", "h2_concentration", "volume_percent",
                "observed_volume_percent", "observed_vol_pct",
                "observed_mean_vol_pct", "observed_time_mean_vol_pct",
                "mean_c_pct", "observed_mean_c_pct", "mean_vol_pct",
            ),
        ),
        "qualification": (
            "Verify sensor identity, calibration/response metadata, units and "
            "the declared temporal averaging operator for quantitative H2 data."
        ),
        "missing_action": (
            "Export fixed-receptor H2 concentration time series with sensor ID, "
            "coordinates, height, units, calibration and averaging metadata."
        ),
    },
    "common_clock": {
        "required_field_groups": _COMMON_CLOCK_FIELD_GROUPS,
        "qualification": (
            "Prove that source, weather, obstacle and receptor timestamps share "
            "the declared clock or an explicit synchronization transform."
        ),
        "missing_action": (
            "Export a common-clock/synchronization ID and timestamp alignment "
            "metadata covering every selected evidence channel."
        ),
    },
}
_CONCENTRATION_FIELD_NAMES = (
    "observed_mole_fraction",
    "mole_fraction",
    "h2_mole_fraction",
    "concentration",
    "h2_concentration",
    "volume_percent",
    "observed_volume_percent",
    "observed_vol_pct",
    "observed_mean_vol_pct",
    "observed_time_mean_vol_pct",
    "mean_c_pct",
    "observed_mean_c_pct",
    "mean_vol_pct",
)
_CLOCK_ID_FIELD_NAMES = (
    "common_clock_id", "clock_id", "synchronized_clock_id",
    "synchronization_id", "sync_id", "time_sync_id", "clock_sync_id",
)
_SYNCHRONIZATION_METADATA_FIELD_NAMES = (
    "synchron", "synchronization", "synchronization_method",
    "time_alignment", "time_alignment_id", "common_time",
    "common_time_id", "clock_offset_s",
)
_SKIP_DIRECTORIES = {
    ".git", ".venv", "venv", "__pycache__", "build", "dist", ".mypy_cache",
    ".pytest_cache", "site-packages", "node_modules", "__pypackages__",
}
# Keep Unicode letters/digits in inventory keys.  External workbooks often use
# Korean (or another local language) labels alongside English telemetry names;
# treating those labels as empty strings makes an otherwise readable workbook
# look corrupt before channel detection has a chance to inspect it.  Channel
# promotion still depends on the explicit ASCII aliases below, so preserving a
# local-language key cannot infer a validation channel by itself.
_TOKEN_RE = re.compile(r"[^\w]+", re.UNICODE)
_PROCESS_TAG_RE = re.compile(r"^(?:pt|tt|lt|dpt)\d+[a-z0-9_]*$")
_MAX_CSV_SAMPLE_ROWS = 200
_MISSING_VALUE_TOKENS = frozenset({
    "", "na", "n_a", "none", "null", "unknown", "unspecified",
    "not_recorded", "not_available", "missing", "-",
})
_DIAGNOSTIC_NOTE_PATTERNS = (
    (
        "fixed_receptor_without_concentration",
        "fixed-receptor geometry or threshold metadata found without a quantitative concentration column;",
    ),
    (
        "threshold_summary_only",
        "threshold crossing summary is not an observed concentration time series",
    ),
    (
        "weather_without_identity",
        "meteorology fields found without an explicit weather identifier and matched interval",
    ),
    (
        "release_rate_without_source_identity",
        "release-rate field found without an explicit source-boundary identifier",
    ),
    (
        "operating_history_without_source_boundary",
        "operating pressure/temperature/level history found without an explicit atmospheric source-rate boundary",
    ),
    (
        "instrument_history_without_source_boundary",
        "instrument-tag operating history found without an explicit atmospheric source-rate boundary",
    ),
    (
        "sensor_time_without_clock",
        "sensor time axis found without a common-clock or synchronization identifier",
    ),
    (
        "unreadable_structured_file",
        "unreadable structured file:",
    ),
    (
        "sampling_limited",
        "sampled from at most",
    ),
    (
        "empty_structured_file",
        "structured file contains no data rows; header/keys were not promoted",
    ),
)
_DIAGNOSTIC_NOTE_KEYS = tuple(key for key, _prefix in _DIAGNOSTIC_NOTE_PATTERNS)


def _normalise_key(value: object) -> str:
    return _TOKEN_RE.sub("_", str(value).strip().lower()).strip("_")


def _json_object_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON keys before Python's parser can overwrite them."""
    seen_raw: set[str] = set()
    seen_normalized: dict[str, str] = {}
    output: dict[str, object] = {}
    for key, value in pairs:
        if key in seen_raw:
            raise ValueError(f"JSON object contains duplicate key: {key!r}")
        normalized = _normalise_key(key)
        if normalized and normalized in seen_normalized:
            previous = seen_normalized[normalized]
            raise ValueError(
                "JSON object contains keys with the same normalized name: "
                f"{previous!r} and {key!r}"
            )
        seen_raw.add(key)
        if normalized:
            seen_normalized[normalized] = key
        output[key] = value
    return output


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _diagnostic_counts(
    candidates: tuple["FieldEvidenceCandidate", ...],
) -> dict[str, int]:
    """Count candidate files by conservative near-miss diagnostic category."""
    counts = Counter({key: 0 for key in _DIAGNOSTIC_NOTE_KEYS})
    for candidate in candidates:
        for key, prefix in _DIAGNOSTIC_NOTE_PATTERNS:
            if any(prefix in note for note in candidate.notes):
                counts[key] += 1
    return dict(counts)


def _collection_requirements(
    channel_paths: tuple[tuple[str, tuple[str, ...]], ...],
    *,
    common_clock_field_groups: tuple[tuple[str, ...], ...] | None = None,
    field_groups_by_channel: Mapping[str, tuple[tuple[str, ...], ...]] | None = None,
) -> dict[str, dict[str, object]]:
    """Build a deterministic, actionable collection plan from channel paths.

    This is a triage aid only.  A detected path is still an unqualified
    candidate and therefore receives the same manual event/clock/geometry
    qualification note as a missing path receives a collection action.
    """
    paths = dict(channel_paths)
    if set(paths) != set(_CHANNELS):
        raise ValueError("collection requirements need every required channel")
    requirements: dict[str, dict[str, object]] = {}
    for channel in _CHANNELS:
        spec = _CHANNEL_COLLECTION_SPECS[channel]
        field_groups = spec["required_field_groups"]
        if field_groups_by_channel is not None and channel in field_groups_by_channel:
            field_groups = field_groups_by_channel[channel]
        if channel == "common_clock" and common_clock_field_groups is not None:
            field_groups = common_clock_field_groups
        candidate_paths = paths[channel]
        requirements[channel] = {
            "status": "detected" if candidate_paths else "missing",
            "candidate_paths": list(candidate_paths),
            "required_field_groups": [
                list(group) for group in field_groups
            ],
            "qualification": spec["qualification"],
            "next_action": (
                spec["qualification"] if candidate_paths
                else spec["missing_action"]
            ),
        }
    return requirements


def _flatten_json_keys(value: object, prefix: str = "") -> set[str]:
    keys: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            normal = _normalise_key(key)
            if normal:
                keys.add(normal)
                keys.update(_flatten_json_keys(child, f"{prefix}_{normal}".strip("_")))
    elif isinstance(value, list):
        for child in value[:200]:
            keys.update(_flatten_json_keys(child, prefix))
    return keys


def _json_value_has_data(value: object) -> bool:
    """Return whether a JSON value contains at least one non-blank leaf."""
    if value is None:
        return False
    if isinstance(value, str):
        return _normalise_key(value) not in _MISSING_VALUE_TOKENS
    if isinstance(value, dict):
        return any(_json_value_has_data(child) for child in value.values())
    if isinstance(value, list):
        return any(_json_value_has_data(child) for child in value)
    return True


def _flatten_json_values(
    value: object, output: dict[str, list[object]] | None = None,
) -> dict[str, list[object]]:
    """Collect bounded JSON field values for value-level channel checks."""
    if output is None:
        output = {}
    if isinstance(value, dict):
        for key, child in value.items():
            normal = _normalise_key(key)
            if normal:
                output.setdefault(normal, []).append(child)
            _flatten_json_values(child, output)
    elif isinstance(value, list):
        for child in value[:200]:
            _flatten_json_values(child, output)
    return output


def _flatten_json_record(
    value: object, output: dict[str, object] | None = None,
) -> dict[str, object]:
    """Flatten one JSON object while retaining its field co-occurrence."""
    if output is None:
        output = {}
    if isinstance(value, dict):
        for key, child in value.items():
            normal = _normalise_key(key)
            if isinstance(child, (dict, list)):
                _flatten_json_record(child, output)
            elif normal:
                output[normal] = child
    elif isinstance(value, list):
        for child in value:
            _flatten_json_record(child, output)
    return output


def _json_record_maps(value: object) -> tuple[dict[str, object], ...]:
    """Return bounded JSON record maps for same-record channel checks."""
    if isinstance(value, list):
        return tuple(
            _flatten_json_record(item)
            for item in value[:200]
            if isinstance(item, dict)
        )
    if isinstance(value, dict):
        for key in ("rows", "observations", "records", "data"):
            child = value.get(key)
            if isinstance(child, list):
                metadata: dict[str, object] = {}
                for metadata_key, metadata_value in value.items():
                    if metadata_key not in {"rows", "observations", "records", "data"}:
                        _flatten_json_record(
                            {metadata_key: metadata_value}, metadata,
                        )
                return tuple(
                    {**metadata, **_flatten_json_record(item)}
                    for item in child[:200]
                    if isinstance(item, dict)
                )
        return (_flatten_json_record(value),)
    return ()


def _channels_for_keys(
    keys: set[str], values: dict[str, list[object]] | None = None,
    records: tuple[dict[str, object], ...] | None = None,
) -> tuple[str, ...]:
    """Classify structured fields and, when available, their values."""
    compact = {key.replace("_", "") for key in keys}
    has = lambda *names: any(name in keys or name.replace("_", "") in compact for name in names)

    def matching_values(names: tuple[str, ...]) -> tuple[object, ...]:
        if values is None:
            return ()
        wanted = {name.replace("_", "") for name in names}
        return tuple(
            item
            for key, items in values.items()
            if key in names or key.replace("_", "") in wanted
            for item in items
        )

    def has_value(*names: str) -> bool:
        if values is None:
            return has(*names)
        return any(_json_value_has_data(item) for item in matching_values(names))

    def has_numeric(*names: str, minimum: float | None = None) -> bool:
        if values is None:
            return has(*names)
        for item in matching_values(names):
            if isinstance(item, (dict, list, tuple, bool)):
                continue
            try:
                number = float(item)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(number):
                continue
            if minimum is not None and number < minimum:
                continue
            return True
        return False

    def record_has_value(record: dict[str, object], *names: str) -> bool:
        return any(
            key in record and _json_value_has_data(record[key])
            for key in names
        )

    def record_has_numeric(
        record: dict[str, object], *names: str, minimum: float | None = None,
    ) -> bool:
        for name in names:
            if name not in record or isinstance(record[name], bool):
                continue
            try:
                number = float(record[name])
            except (TypeError, ValueError):
                continue
            if math.isfinite(number) and (
                minimum is None or number >= minimum
            ):
                return True
        return False

    def any_record(predicate) -> bool:
        if records is None:
            return False
        return any(predicate(record) for record in records)

    source_rate = has(
        "mass_flow", "mass_flow_kg_s", "release_rate", "release_rate_kg_s",
        "leak_rate", "source_rate", "source_rate_kg_s", "mass_rate_kg_s",
        "mass_rate_reported_kg_s", "rate_kg_s", "inflow_kg_s", "erate",
    )
    source_id = has(
        "source_boundary_id", "source_id", "release_id",
        "source_boundary", "source_bound",
    )
    weather_id = has("weather_id", "wind_id", "meteorology_id")
    wind = has(
        "wind_speed", "wind_speed_m_s", "wind_direction", "wind_direction_deg",
        "wind_m_s", "mean_wind_speed_m_s", "wind_speed_high_m_s",
        "wind_speed_low_m_s", "wind_from_deg", "mean_wind_from_deg",
        "wind_u", "wind_v", "anemometer", "wind_ref_m_s", "u_ref_m_s",
        "wind_10m_m_s",
    )
    obstacle_id = has("obstacle_geometry_id", "obstacle_id", "geometry_id")
    obstacle_shape = (
        has("obstacle", "barrier", "building", "wall")
        and has("geometry", "height", "width", "length", "polygon", "vertices", "position")
    )
    sensor = has("sensor", "sensor_id", "receptor", "receptor_id")
    concentration = has(*_CONCENTRATION_FIELD_NAMES)
    coordinates = (
        has("x_downwind_m", "x_m", "x")
        and has("y_crosswind_m", "y_m", "y")
        and has("height_m", "z_m", "z", "elevation_m")
    )
    clock = has(*_CLOCK_ID_FIELD_NAMES)
    timestamp = has("timestamp", "time_s", "observation_time_s", "datetime", "utc_time")

    detected: list[str] = []
    # A physical rate without an explicit boundary identifier is only a
    # literature/source inventory, not a site source boundary that can be
    # joined to a weather and receptor record.
    source_valid = (
        any_record(lambda record: (
            record_has_value(
                record, "source_boundary_id", "source_id", "release_id",
                "source_boundary", "source_bound",
            )
            and record_has_numeric(
                record,
                "mass_flow", "mass_flow_kg_s", "release_rate", "release_rate_kg_s",
                "leak_rate", "source_rate", "source_rate_kg_s", "mass_rate_kg_s",
                "mass_rate_reported_kg_s", "rate_kg_s", "inflow_kg_s", "erate",
                minimum=0.0,
            )
        ))
        if records is not None else (
            has_value(
                "source_boundary_id", "source_id", "release_id",
                "source_boundary", "source_bound",
            )
            and has_numeric(
                "mass_flow", "mass_flow_kg_s", "release_rate", "release_rate_kg_s",
                "leak_rate", "source_rate", "source_rate_kg_s", "mass_rate_kg_s",
                "mass_rate_reported_kg_s", "rate_kg_s", "inflow_kg_s", "erate",
                minimum=0.0,
            )
        )
    )
    if source_rate and source_id and source_valid:
        detected.append("source_boundary")
    # A wind/temperature pair without an explicit weather identity cannot be
    # joined safely to a source, obstacle or receptor record. Keep it as a
    # diagnostic note below rather than promoting it to a validation channel.
    weather_valid = (
        any_record(lambda record: (
            record_has_value(record, "weather_id", "wind_id", "meteorology_id")
            and record_has_numeric(
                record,
                "wind_speed", "wind_speed_m_s", "wind_m_s", "mean_wind_speed_m_s",
                "wind_speed_high_m_s", "wind_speed_low_m_s", "wind_u", "wind_v",
                "wind_ref_m_s", "u_ref_m_s", "wind_10m_m_s",
            )
        ))
        if records is not None else (
            has_value("weather_id", "wind_id", "meteorology_id")
            and has_numeric(
                "wind_speed", "wind_speed_m_s", "wind_m_s", "mean_wind_speed_m_s",
            "wind_speed_high_m_s", "wind_speed_low_m_s", "wind_u", "wind_v",
            "wind_ref_m_s", "u_ref_m_s", "wind_10m_m_s",
            )
        )
    )
    if weather_id and wind and weather_valid:
        detected.append("weather")
    obstacle_valid = (
        any_record(lambda record: (
            record_has_value(record, "obstacle_geometry_id", "obstacle_id", "geometry_id")
            or record_has_value(
                record,
                "geometry", "height", "width", "length", "polygon", "vertices", "position",
            )
        ))
        if records is not None else (
            has_value("obstacle_geometry_id", "obstacle_id", "geometry_id")
            or has_value(
                "geometry", "height", "width", "length", "polygon", "vertices", "position",
            )
        )
    )
    if (obstacle_id or obstacle_shape) and obstacle_valid:
        detected.append("obstacle_geometry")
    receptor_valid = (
        any_record(lambda record: (
            record_has_value(record, "sensor", "sensor_id", "receptor", "receptor_id")
            and record_has_numeric(
                record,
                *_CONCENTRATION_FIELD_NAMES, minimum=0.0,
            )
            and record_has_numeric(record, "x_downwind_m", "x_m", "x")
            and record_has_numeric(record, "y_crosswind_m", "y_m", "y")
            and record_has_numeric(record, "height_m", "z_m", "z", "elevation_m")
        ))
        if records is not None else (
            has_value("sensor", "sensor_id", "receptor", "receptor_id")
            and has_numeric(*_CONCENTRATION_FIELD_NAMES, minimum=0.0)
            and has_numeric("x_downwind_m", "x_m", "x")
            and has_numeric("y_crosswind_m", "y_m", "y")
            and has_numeric("height_m", "z_m", "z", "elevation_m")
        )
    )
    if sensor and concentration and coordinates and receptor_valid:
        detected.append("receptor_observations")
    clock_valid = (
        any_record(lambda record: (
            record_has_value(record, *_CLOCK_ID_FIELD_NAMES)
            or (
                record_has_value(record, "timestamp", "time_s", "observation_time_s", "datetime", "utc_time")
                and record_has_value(record, *_SYNCHRONIZATION_METADATA_FIELD_NAMES)
            )
        ))
        if records is not None else (
            (clock and has_value(*_CLOCK_ID_FIELD_NAMES))
            or (
                timestamp
                and has(*_SYNCHRONIZATION_METADATA_FIELD_NAMES)
                and has_value("timestamp", "time_s", "observation_time_s", "datetime", "utc_time")
                and has_value(*_SYNCHRONIZATION_METADATA_FIELD_NAMES)
            )
        )
    )
    if (clock or (timestamp and has(*_SYNCHRONIZATION_METADATA_FIELD_NAMES))) and clock_valid:
        detected.append("common_clock")
    return tuple(detected)


@dataclass(frozen=True)
class FieldEvidenceCandidate:
    """One structured file with conservatively detected evidence channels."""

    path: str
    file_type: Literal["csv", "json", "xlsx"]
    row_count: int | None
    fields: tuple[str, ...]
    detected_channels: tuple[str, ...]
    notes: tuple[str, ...] = ()
    sha256: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.path, str) or not self.path.strip():
            raise ValueError("audit candidate path must be non-empty")
        if self.file_type not in {"csv", "json", "xlsx"}:
            raise ValueError("audit candidate file_type is unsupported")
        if self.row_count is not None and (
            isinstance(self.row_count, bool)
            or not isinstance(self.row_count, int)
            or self.row_count < 0
        ):
            raise ValueError("audit candidate row_count must be non-negative")
        for name, values in {
            "fields": self.fields,
            "detected_channels": self.detected_channels,
            "notes": self.notes,
        }.items():
            if not isinstance(values, tuple):
                raise TypeError(f"audit candidate {name} must be a tuple")
            if any(not isinstance(value, str) or not value.strip() for value in values):
                raise ValueError(f"audit candidate {name} must contain non-empty strings")
        if self.sha256 is not None and (
            not isinstance(self.sha256, str)
            or len(self.sha256) != 64
            or any(character not in "0123456789abcdef" for character in self.sha256)
        ):
            raise ValueError("audit candidate sha256 must be a lowercase SHA-256 digest")
        unknown = set(self.detected_channels) - set(_CHANNELS)
        if unknown:
            raise ValueError(f"audit candidate has unknown channels: {sorted(unknown)}")
        if len(set(self.detected_channels)) != len(self.detected_channels):
            raise ValueError("audit candidate detected_channels must be unique")

    def as_record(self) -> dict[str, object]:
        return {
            "path": self.path,
            "file_type": self.file_type,
            "row_count": self.row_count,
            "fields": list(self.fields),
            "detected_channels": list(self.detected_channels),
            "notes": list(self.notes),
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class FieldEvidenceAudit:
    """Fail-safe inventory of validation-evidence readiness."""

    root: str
    scanned_files: int
    candidates: tuple[FieldEvidenceCandidate, ...]
    status: Literal["withheld", "partial", "candidate_complete"]
    channel_paths: tuple[tuple[str, tuple[str, ...]], ...]
    reasons: tuple[str, ...]
    scan_complete: bool = True
    coherent_candidate_paths: tuple[str, ...] = ()
    max_files: int = 20000
    gate_codes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.root, str) or not self.root.strip():
            raise ValueError("audit root must be non-empty")
        if isinstance(self.scanned_files, bool) or not isinstance(self.scanned_files, int) or self.scanned_files < 0:
            raise ValueError("audit scanned_files must be non-negative")
        if isinstance(self.max_files, bool) or not isinstance(self.max_files, int) or self.max_files <= 0:
            raise ValueError("audit max_files must be a positive integer")
        if self.status not in {"withheld", "partial", "candidate_complete"}:
            raise ValueError("audit status is unsupported")
        if not isinstance(self.scan_complete, bool):
            raise TypeError("audit scan_complete must be boolean")
        if not isinstance(self.coherent_candidate_paths, tuple):
            raise TypeError("audit coherent_candidate_paths must be a tuple")
        if any(
            not isinstance(path, str) or not path.strip()
            for path in self.coherent_candidate_paths
        ):
            raise ValueError(
                "audit coherent_candidate_paths must contain non-empty strings"
            )
        if len(set(self.coherent_candidate_paths)) != len(self.coherent_candidate_paths):
            raise ValueError("audit coherent_candidate_paths must be unique")
        if not isinstance(self.candidates, tuple):
            raise TypeError("audit candidates must be a tuple")
        if any(not isinstance(candidate, FieldEvidenceCandidate) for candidate in self.candidates):
            raise TypeError("audit candidates must contain FieldEvidenceCandidate values")
        candidate_paths = tuple(candidate.path for candidate in self.candidates)
        if len(set(candidate_paths)) != len(candidate_paths):
            raise ValueError("audit candidates must have unique paths")
        if self.scanned_files < len(candidate_paths):
            raise ValueError(
                "audit scanned_files cannot be smaller than the candidate count"
            )
        candidate_by_path = {candidate.path: candidate for candidate in self.candidates}
        if set(self.coherent_candidate_paths) - set(candidate_by_path):
            raise ValueError(
                "audit coherent_candidate_paths must refer to candidates"
            )
        if any(
            set(candidate_by_path[path].detected_channels) != set(_CHANNELS)
            for path in self.coherent_candidate_paths
        ):
            raise ValueError(
                "audit coherent candidate must expose every required channel"
            )
        if not isinstance(self.channel_paths, tuple):
            raise TypeError("audit channel_paths must be a tuple")
        if any(
            not isinstance(item, tuple) or len(item) != 2
            for item in self.channel_paths
        ):
            raise TypeError("audit channel_paths entries must be (channel, paths) tuples")
        paths = dict(self.channel_paths)
        if set(paths) != set(_CHANNELS):
            raise ValueError("audit channel_paths must contain every required channel")
        if len(paths) != len(self.channel_paths):
            raise ValueError("audit channel_paths must contain each channel once")
        for channel, channel_values in self.channel_paths:
            if not isinstance(channel, str) or channel not in _CHANNELS:
                raise ValueError("audit channel_paths contains an unsupported channel")
            if not isinstance(channel_values, tuple):
                raise TypeError("audit channel paths must be tuples")
            if any(not isinstance(path, str) or not path.strip() for path in channel_values):
                raise ValueError("audit channel paths must contain non-empty strings")
            if len(set(channel_values)) != len(channel_values):
                raise ValueError("audit channel paths must be unique")
        expected_channel_paths = {
            channel: tuple(
                candidate.path
                for candidate in self.candidates
                if channel in candidate.detected_channels
            )
            for channel in _CHANNELS
        }
        if paths != expected_channel_paths:
            raise ValueError(
                "audit channel_paths must match candidate detected_channels"
            )
        expected_coherent_paths = tuple(
            candidate.path
            for candidate in self.candidates
            if set(candidate.detected_channels) == set(_CHANNELS)
        )
        if self.coherent_candidate_paths != expected_coherent_paths:
            raise ValueError(
                "audit coherent_candidate_paths must match candidate detected_channels"
            )
        if not isinstance(self.reasons, tuple):
            raise TypeError("audit reasons must be a tuple")
        if any(not isinstance(reason, str) or not reason.strip() for reason in self.reasons):
            raise ValueError("audit reasons must contain non-empty strings")
        if not isinstance(self.gate_codes, tuple):
            raise TypeError("audit gate_codes must be a tuple")
        if any(code not in FIELD_EVIDENCE_AUDIT_GATE_CODES for code in self.gate_codes):
            raise ValueError("audit gate_codes contain an unsupported code")
        if len(set(self.gate_codes)) != len(self.gate_codes):
            raise ValueError("audit gate_codes must be unique")
        if self.status == "withheld" and any(paths.values()):
            raise ValueError("withheld audit cannot contain detected channels")
        if self.status == "partial" and not any(paths.values()):
            raise ValueError("partial audit must contain at least one detected channel")
        if self.status == "partial" and self.scan_complete and all(paths.values()):
            raise ValueError(
                "complete audit with every channel detected must be candidate_complete"
            )
        if self.status == "candidate_complete" and not all(paths.values()):
            raise ValueError("candidate_complete audit requires every channel")
        if self.status == "candidate_complete" and not self.scan_complete:
            raise ValueError("incomplete audit cannot be candidate_complete")
        expected_status = (
            "withheld" if not any(paths.values()) else "partial"
        ) if not self.scan_complete else (
            "withheld" if not any(paths.values()) else
            "partial" if not all(paths.values()) else "candidate_complete"
        )
        if self.status != expected_status:
            raise ValueError(
                "audit status does not match scan completeness and channel paths"
            )

    @property
    def promotion_allowed(self) -> bool:
        """Always false: an inventory cannot prove identity, clock, or geometry."""
        return False

    @property
    def diagnostic_counts(self) -> dict[str, int]:
        """Return reproducible counts of candidate near-miss diagnostics.

        Counts are candidate files, not rows or occurrences.  A file is counted
        at most once per category, so a large historian export cannot dominate
        the readiness summary merely by containing many repeated records.
        """
        return _diagnostic_counts(self.candidates)

    @property
    def collection_requirements(self) -> dict[str, dict[str, object]]:
        """Return the next evidence-collection action for every channel.

        Detected candidates are deliberately labelled ``detected`` rather than
        ``qualified``: the audit still cannot establish event identity, clock,
        calibration or geometry equivalence.
        """
        return _collection_requirements(self.channel_paths)

    def as_record(self) -> dict[str, object]:
        return {
            "schema": FIELD_EVIDENCE_AUDIT_SCHEMA,
            "root": self.root,
            "scanned_files": self.scanned_files,
            "max_files": self.max_files,
            "scan_complete": self.scan_complete,
            "status": self.status,
            "promotion_allowed": False,
            "channel_paths": {channel: list(paths) for channel, paths in self.channel_paths},
            "coherent_candidate_paths": list(self.coherent_candidate_paths),
            "candidates": [candidate.as_record() for candidate in self.candidates],
            "diagnostic_counts": self.diagnostic_counts,
            "collection_requirements": self.collection_requirements,
            "reasons": list(self.reasons),
            "gate_codes": list(self.gate_codes),
        }


def _xlsx_sample(
    path: Path,
) -> tuple[
    set[str], tuple[str, ...], bool, dict[str, list[object]], tuple[dict[str, object], ...],
]:
    """Read bounded worksheet samples without treating workbook metadata as evidence."""
    try:
        import openpyxl
    except ImportError as error:  # pragma: no cover - optional dependency
        raise RuntimeError("XLSX audit requires the optional openpyxl package") from error
    fields: set[str] = set()
    has_data_row = False
    values: dict[str, list[object]] = {}
    records: list[dict[str, object]] = []
    notes: list[str] = [
        "XLSX fields were sampled from at most the first 40 rows and 200 cells per sheet; "
        "row_count is not a dataset count",
    ]
    book = None
    usable_sheet_count = 0
    malformed_sheet_count = 0
    try:
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        for sheet in book.worksheets:
            sampled_rows: list[tuple[object, ...]] = []
            for row_index, row in enumerate(sheet.iter_rows(values_only=True), start=1):
                if row_index > 40:
                    break
                sampled_rows.append(tuple(row[:200]))
            header_index: int | None = None
            for index, row_values in enumerate(sampled_rows):
                nonempty = any(
                    value is not None and str(value).strip()
                    for value in row_values
                )
                if nonempty:
                    header_index = index
                    break
            if header_index is None:
                continue
            header_row = sampled_rows[header_index]
            headers = tuple(
                _normalise_key(value) if value is not None else ""
                for value in header_row
            )
            try:
                if any(not header for header in headers):
                    raise ValueError("XLSX header contains an empty field")
                duplicate_headers = sorted(
                    header for header, count in Counter(headers).items()
                    if count > 1
                )
                if duplicate_headers:
                    raise ValueError(
                        "XLSX header contains duplicate fields: "
                        + ", ".join(duplicate_headers)
                    )
            except ValueError as error:
                malformed_sheet_count += 1
                notes.append(
                    f"worksheet {sheet.title!r} skipped: {error}"
                )
                continue
            usable_sheet_count += 1
            fields.update(headers)
            for row_index, row_values in enumerate(
                sampled_rows[header_index + 1:], start=header_index + 1
            ):
                nonempty = any(
                    value is not None and str(value).strip()
                    for value in row_values
                )
                if nonempty:
                    has_data_row = True
                    record: dict[str, object] = {}
                    for header, value in zip(headers, row_values):
                        if header and value is not None:
                            values.setdefault(header, []).append(value)
                            record[header] = value
                    records.append(record)
            for row_values in sampled_rows:
                for value in row_values:
                    if isinstance(value, str) and value.strip():
                        normal = _normalise_key(value)
                        if normal:
                            fields.add(normal)
        if usable_sheet_count == 0 and malformed_sheet_count:
            raise ValueError("all XLSX worksheets have ambiguous headers")
    finally:
        if book is not None:
            book.close()
    return fields, notes, has_data_row, values, tuple(records)


def _candidate_from_file(path: Path) -> FieldEvidenceCandidate | None:
    suffix = path.suffix.lower()
    if suffix not in {".csv", ".json", ".xlsx"}:
        return None
    notes: list[str] = []
    row_count: int | None = None
    fields: set[str] = set()
    has_data_records = True
    field_values: dict[str, list[object]] = {}
    records: tuple[dict[str, object], ...] = ()
    sha256: str | None = None
    try:
        sha256 = _sha256_file(path)
        if suffix == ".csv":
            with path.open("r", newline="", encoding="utf-8-sig") as handle:
                reader = csv.reader(handle)
                raw_headers = next(reader, None)
                if not raw_headers:
                    fields = set()
                    rows: list[tuple[str, ...]] = []
                    row_count = 0
                else:
                    headers = tuple(_normalise_key(value) for value in raw_headers)
                    if any(not header for header in headers):
                        raise ValueError("CSV header contains an empty field")
                    duplicate_headers = sorted(
                        header for header, count in Counter(headers).items()
                        if count > 1
                    )
                    if duplicate_headers:
                        raise ValueError(
                            "CSV header contains duplicate fields: "
                            + ", ".join(duplicate_headers)
                        )
                    rows = []
                    row_count = 0
                    has_data_records = False
                    for row in reader:
                        # Blank physical lines are not observations and are
                        # ignored, matching DictReader's empty-row behavior.
                        if not row:
                            continue
                        if len(row) != len(headers):
                            raise ValueError(
                                "CSV row has a different number of fields than its header"
                            )
                        row_count += 1
                        has_data_records = has_data_records or any(
                            value is not None and str(value).strip()
                            for value in row
                        )
                        if len(rows) < _MAX_CSV_SAMPLE_ROWS:
                            rows.append(tuple(row))
                    fields = set(headers)
                    if row_count > _MAX_CSV_SAMPLE_ROWS:
                        notes.append(
                            "CSV fields and channel values were sampled from at most "
                            f"the first {_MAX_CSV_SAMPLE_ROWS} data rows; row_count is "
                            "the complete dataset count"
                        )
                csv_records: list[dict[str, object]] = []
                for row in rows:
                    record = dict(zip(headers, row)) if raw_headers else {}
                    csv_records.append(record)
                    for key, value in record.items():
                        field_values.setdefault(key, []).append(value)
                records = tuple(csv_records)
            if not fields:
                notes.append("CSV has no header fields")
        elif suffix == ".json":
            payload = json.loads(
                path.read_text(encoding="utf-8-sig"),
                object_pairs_hook=_json_object_pairs,
            )
            fields = _flatten_json_keys(payload)
            field_values = _flatten_json_values(payload)
            records = _json_record_maps(payload)
            if isinstance(payload, list):
                row_count = len(payload)
                has_data_records = any(_json_value_has_data(item) for item in payload)
            elif isinstance(payload, dict):
                has_data_records = _json_value_has_data(payload)
                for key in ("rows", "observations", "records", "data"):
                    value = payload.get(key)
                    if isinstance(value, list):
                        row_count = len(value)
                        has_data_records = any(
                            _json_value_has_data(item) for item in value
                        )
                        break
        else:
            fields, xlsx_notes, has_data_records, field_values, records = _xlsx_sample(path)
            notes.extend(xlsx_notes)
    except (
        OSError, UnicodeError, csv.Error, json.JSONDecodeError, ValueError,
        RuntimeError, zipfile.BadZipFile,
    ) as error:
        file_type = "csv" if suffix == ".csv" else "json" if suffix == ".json" else "xlsx"
        return FieldEvidenceCandidate(
            path=str(path.resolve()), file_type=file_type,
            row_count=None, fields=tuple(), detected_channels=tuple(),
            notes=(f"unreadable structured file: {type(error).__name__}",),
            sha256=sha256,
        )
    # A header-only CSV/JSON object is not a measured source, weather record,
    # or fixed-receptor observation.  Keep it as a diagnostic candidate so the
    # audit explains what was found, but do not promote its field names into a
    # readiness channel.  JSON geometry metadata without a row-list keeps
    # ``row_count=None`` and remains eligible because it can be a real object
    # rather than a tabular export.
    if suffix != ".xlsx":
        has_data_records = has_data_records and row_count != 0
    channels = (
        _channels_for_keys(fields, field_values, records)
        if has_data_records else ()
    )
    if not has_data_records:
        notes.append(
            "structured file contains no data rows; header/keys were not promoted "
            "to an evidence channel"
        )
    compact = {key.replace("_", "") for key in fields}
    has = lambda *names: any(
        name in fields or name.replace("_", "") in compact for name in names
    )
    sensor_geometry = (
        has("sensor", "sensor_id", "receptor", "receptor_id")
        and (
            (
                has("x_downwind_m", "x_m", "x")
                and has("y_crosswind_m", "y_m", "y")
                and has("height_m", "z_m", "z", "elevation_m")
            )
            or (has("radius_m") and has("angle_deg", "bearing_deg") and has("height_m"))
        )
    )
    if sensor_geometry and "receptor_observations" not in channels:
        notes.append(
            "fixed-receptor geometry or threshold metadata found without a quantitative "
            "concentration column; not a validation observation"
        )
    if has(
        "first_above_4pct_s", "final_below_4pct_s", "threshold_crossing_s",
        "published_first_to_last_span_s",
    ):
        notes.append(
            "threshold crossing summary is not an observed concentration time series"
        )
    if (
        has(
            "wind_speed", "wind_speed_m_s", "wind_m_s", "mean_wind_speed_m_s",
            "wind_speed_high_m_s", "wind_speed_low_m_s", "wind_direction",
            "wind_direction_deg", "wind_from_deg", "mean_wind_from_deg",
            "wind_ref_m_s", "u_ref_m_s", "wind_10m_m_s",
        )
        and "weather" not in channels
    ):
        notes.append(
            "meteorology fields found without an explicit weather identifier and matched interval"
        )
    if not has(
        "source_boundary_id", "source_id", "release_id",
        "source_boundary", "source_bound",
    ) and has(
        "mass_flow", "mass_flow_kg_s", "release_rate", "release_rate_kg_s",
        "leak_rate", "source_rate", "source_rate_kg_s", "mass_rate_kg_s",
        "mass_rate_reported_kg_s", "rate_kg_s", "inflow_kg_s", "erate",
    ):
        notes.append(
            "release-rate field found without an explicit source-boundary identifier"
        )
    if (
        has(
            "pressure", "pressure_pa", "upstream_pressure", "tank_pressure",
            "temperature", "temperature_k", "tank_temperature", "level",
            "liquid_level", "fill_level",
        )
        and has("timestamp", "time_s", "observation_time_s", "datetime", "utc_time")
        and "source_boundary" not in channels
    ):
        notes.append(
            "operating pressure/temperature/level history found without an explicit "
            "atmospheric source-rate boundary; not a dispersion source"
        )
    elif (
        has("timestamp", "time_s", "observation_time_s", "datetime", "utc_time")
        and any(_PROCESS_TAG_RE.fullmatch(key) for key in fields)
        and "source_boundary" not in channels
    ):
        notes.append(
            "instrument-tag operating history found without an explicit atmospheric "
            "source-rate boundary; not a dispersion source"
        )
    if (
        has("timestamp", "time_s", "observation_time_s", "datetime", "utc_time")
        and not has(*_CLOCK_ID_FIELD_NAMES)
        and not has(*_SYNCHRONIZATION_METADATA_FIELD_NAMES)
        and sensor_geometry
    ):
        notes.append(
            "sensor time axis found without a common-clock or synchronization identifier"
        )
    if not channels and not notes:
        return None
    return FieldEvidenceCandidate(
        path=str(path.resolve()),
        file_type="csv" if suffix == ".csv" else "json" if suffix == ".json" else "xlsx",
        row_count=row_count, fields=tuple(sorted(fields)), detected_channels=channels,
        notes=tuple(notes), sha256=sha256,
    )


def audit_field_evidence(root: str | Path, *, max_files: int = 20000) -> FieldEvidenceAudit:
    """Inventory structured evidence channels below ``root`` without importing data.

    ``max_files`` is a safety bound.  Excluded virtual-environment/build trees
    keep a repository scan useful while preserving explicit external data paths.
    """
    base = Path(root)
    if not base.is_dir():
        raise NotADirectoryError(base)
    if not isinstance(max_files, int) or max_files <= 0:
        raise ValueError("max_files must be a positive integer")
    files: list[Path] = []
    scan_complete = True
    for path in base.rglob("*"):
        if not path.is_file() or any(
            part in _SKIP_DIRECTORIES
            or part.startswith(".venv")
            or (part.startswith("venv") and len(part) > 4)
            for part in path.parts
        ):
            continue
        if path.suffix.lower() not in {".csv", ".json", ".xlsx"}:
            continue
        if len(files) >= max_files:
            scan_complete = False
            break
        files.append(path)
    candidates = tuple(
        candidate for path in sorted(files) if (candidate := _candidate_from_file(path)) is not None
    )
    channel_paths = tuple(
        (channel, tuple(candidate.path for candidate in candidates if channel in candidate.detected_channels))
        for channel in _CHANNELS
    )
    coherent_candidate_paths = tuple(
        candidate.path
        for candidate in candidates
        if set(candidate.detected_channels) == set(_CHANNELS)
    )
    detected = [paths for _, paths in channel_paths]
    gate_codes: list[str] = ["promotion_not_allowed"]
    if not scan_complete:
        gate_codes.append("scan_incomplete")
        missing = [channel for channel, paths in channel_paths if not paths]
        status: Literal["withheld", "partial", "candidate_complete"] = (
            "withheld" if not any(detected) else "partial"
        )
        gate_codes.append("channels_missing" if not any(detected) else "channels_partial")
        reasons = (
            f"scan stopped at max_files={max_files}; the directory was not fully inspected",
            "required channels are missing in the scanned prefix: " + ", ".join(missing)
            if missing else
            "all channels appeared in the scanned prefix, but candidate_complete is withheld "
            "because the directory was not fully inspected",
        )
    elif not any(detected):
        status: Literal["withheld", "partial", "candidate_complete"] = "withheld"
        gate_codes.append("channels_missing")
        reasons = ("no structured file exposed any required field-validation channel",)
    elif not all(detected):
        status = "partial"
        gate_codes.append("channels_partial")
        missing = [channel for channel, paths in channel_paths if not paths]
        reasons = ("required channels are missing: " + ", ".join(missing),)
    else:
        status = "candidate_complete"
        if coherent_candidate_paths:
            gate_codes.append("coherent_candidate_unqualified")
            reasons = (
                "all channels were detected, including at least one same-file coherent "
                "candidate, but identity, common clock, geometry and calibration still "
                "require manual qualification",
            )
        else:
            gate_codes.append("cross_file_candidate_unqualified")
            reasons = (
                "all channels were detected across candidates, but no single file "
                "contained every channel; cross-file identity, common clock, obstacle "
                "geometry and calibration still require manual qualification",
            )
    reasons += (
        "this audit never creates FieldValidationEvidence and cannot qualify a model or design basis",
    )
    return FieldEvidenceAudit(
        root=str(base.resolve()), scanned_files=len(files), candidates=candidates,
        status=status, channel_paths=channel_paths, reasons=reasons,
        scan_complete=scan_complete,
        coherent_candidate_paths=coherent_candidate_paths,
        max_files=max_files,
        gate_codes=tuple(dict.fromkeys(gate_codes)),
    )


def field_evidence_audit_record(audit: FieldEvidenceAudit) -> dict[str, object]:
    """Return the stable JSON record for a read-only readiness audit."""
    if not isinstance(audit, FieldEvidenceAudit):
        raise TypeError("audit must be a FieldEvidenceAudit")
    return audit.as_record()


def _audit_mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _audit_string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _audit_strings(value: object, name: str) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be a JSON array")
    return tuple(_audit_string(item, f"{name}[]") for item in value)


def _audit_diagnostic_counts(value: object, name: str) -> dict[str, int]:
    data = _audit_mapping(value, name)
    unknown = sorted(set(data) - set(_DIAGNOSTIC_NOTE_KEYS))
    if unknown:
        raise ValueError(f"{name} keys are invalid: unknown=" + ", ".join(unknown))
    counts: dict[str, int] = {}
    for key in _DIAGNOSTIC_NOTE_KEYS:
        # Older v1 artifacts may predate a newly added diagnostic category.
        # Missing categories are therefore interpreted as zero, while any
        # category present in the artifact is still checked exactly.
        count = data.get(key, 0)
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f"{name}.{key} must be a non-negative integer")
        counts[key] = count
    return counts


def _audit_collection_requirements(
    value: object,
    name: str,
    audit: FieldEvidenceAudit,
) -> None:
    """Validate the optional machine-actionable collection plan exactly."""
    data = _audit_mapping(value, name)
    if set(data) != set(_CHANNELS):
        raise ValueError(f"{name} must contain exactly the required channels")
    normalized: dict[str, dict[str, object]] = {}
    for channel in _CHANNELS:
        item = _audit_mapping(data[channel], f"{name}.{channel}")
        _audit_keys(
            item,
            f"{name}.{channel}",
            required={
                "status", "candidate_paths", "required_field_groups",
                "qualification", "next_action",
            },
        )
        status = _audit_string(item["status"], f"{name}.{channel}.status")
        if status not in {"missing", "detected"}:
            raise ValueError(f"{name}.{channel}.status is unsupported")
        candidate_paths = _audit_strings(
            item["candidate_paths"], f"{name}.{channel}.candidate_paths",
        )
        groups_value = item["required_field_groups"]
        if not isinstance(groups_value, list):
            raise ValueError(
                f"{name}.{channel}.required_field_groups must be a JSON array"
            )
        groups = [
            list(_audit_strings(group, f"{name}.{channel}.required_field_groups[]"))
            for group in groups_value
        ]
        normalized[channel] = {
            "status": status,
            "candidate_paths": list(candidate_paths),
            "required_field_groups": groups,
            "qualification": _audit_string(
                item["qualification"], f"{name}.{channel}.qualification",
            ),
            "next_action": _audit_string(
                item["next_action"], f"{name}.{channel}.next_action",
            ),
        }
    expected = audit.collection_requirements
    if normalized != expected:
        # Keep already-published v1 artifacts readable when their optional
        # collection plan predates the expanded source/weather/receptor and
        # synchronization aliases. New audits emit the expanded contract;
        # only this exact historical representation is accepted as equivalent.
        legacy_expected = _collection_requirements(
            audit.channel_paths,
            field_groups_by_channel=_LEGACY_COLLECTION_FIELD_GROUPS,
        )
        legacy_clock_only = _collection_requirements(
            audit.channel_paths,
            common_clock_field_groups=_LEGACY_COMMON_CLOCK_FIELD_GROUPS,
        )
        legacy_non_clock = _collection_requirements(
            audit.channel_paths,
            common_clock_field_groups=_COMMON_CLOCK_FIELD_GROUPS,
            field_groups_by_channel={
                "source_boundary": _LEGACY_SOURCE_BOUNDARY_FIELD_GROUPS,
                "weather": _LEGACY_WEATHER_FIELD_GROUPS,
                "receptor_observations": _LEGACY_RECEPTOR_FIELD_GROUPS,
            },
        )
        if (
            normalized == legacy_expected
            or normalized == legacy_clock_only
            or normalized == legacy_non_clock
        ):
            return
        raise ValueError(f"{name} must match channel paths and collection contract")


def _audit_keys(
    value: Mapping[str, object], name: str, *, required: set[str],
    optional: set[str] | None = None,
) -> None:
    optional = set() if optional is None else optional
    missing = sorted(required - set(value))
    unknown = sorted(set(value) - required - optional)
    if missing or unknown:
        parts: list[str] = []
        if missing:
            parts.append("missing=" + ", ".join(missing))
        if unknown:
            parts.append("unknown=" + ", ".join(unknown))
        raise ValueError(f"{name} keys are invalid: " + "; ".join(parts))


def field_evidence_audit_from_mapping(value: object) -> FieldEvidenceAudit:
    """Parse one strict :func:`field_evidence_audit_record` mapping."""
    data = _audit_mapping(value, "field evidence audit")
    _audit_keys(
        data,
        "field evidence audit",
        required={
            "schema", "root", "scanned_files", "scan_complete", "status",
            "promotion_allowed", "channel_paths", "coherent_candidate_paths",
            "candidates", "reasons",
        },
        optional={
            "max_files", "gate_codes", "diagnostic_counts",
            "collection_requirements",
        },
    )
    if data["schema"] != FIELD_EVIDENCE_AUDIT_SCHEMA:
        raise ValueError("schema must be " + repr(FIELD_EVIDENCE_AUDIT_SCHEMA))
    if data["promotion_allowed"] is not False:
        raise ValueError("promotion_allowed must be false")
    scanned_files = data["scanned_files"]
    if isinstance(scanned_files, bool) or not isinstance(scanned_files, int) or scanned_files < 0:
        raise ValueError("scanned_files must be a non-negative integer")
    max_files = data.get("max_files", 20000)
    if isinstance(max_files, bool) or not isinstance(max_files, int) or max_files <= 0:
        raise ValueError("max_files must be a positive integer")
    scan_complete = data["scan_complete"]
    if not isinstance(scan_complete, bool):
        raise ValueError("scan_complete must be boolean")
    candidates_value = data["candidates"]
    if not isinstance(candidates_value, list):
        raise ValueError("candidates must be a JSON array")
    candidates: list[FieldEvidenceCandidate] = []
    for index, candidate_value in enumerate(candidates_value):
        candidate_data = _audit_mapping(candidate_value, f"candidates[{index}]")
        _audit_keys(
            candidate_data,
            f"candidates[{index}]",
            required={
                "path", "file_type", "row_count", "fields", "detected_channels", "notes",
            },
            optional={"sha256"},
        )
        row_count = candidate_data["row_count"]
        if row_count is not None and (
            isinstance(row_count, bool) or not isinstance(row_count, int) or row_count < 0
        ):
            raise ValueError(f"candidates[{index}].row_count must be non-negative or null")
        sha256 = candidate_data.get("sha256")
        if sha256 is not None:
            sha256 = _audit_string(sha256, f"candidates[{index}].sha256").lower()
        fields = _audit_strings(candidate_data["fields"], f"candidates[{index}].fields")
        detected_channels = _audit_strings(
            candidate_data["detected_channels"], f"candidates[{index}].detected_channels",
        )
        notes = _audit_strings(candidate_data["notes"], f"candidates[{index}].notes")
        candidates.append(
            FieldEvidenceCandidate(
                path=_audit_string(candidate_data["path"], f"candidates[{index}].path"),
                file_type=_audit_string(
                    candidate_data["file_type"], f"candidates[{index}].file_type",
                ),
                row_count=row_count,
                fields=fields,
                detected_channels=detected_channels,
                notes=notes,
                sha256=sha256,
            )
        )
    channel_data = _audit_mapping(data["channel_paths"], "channel_paths")
    if set(channel_data) != set(_CHANNELS):
        raise ValueError("channel_paths must contain exactly the required channels")
    channel_paths = tuple(
        (
            channel,
            _audit_strings(channel_data[channel], f"channel_paths.{channel}"),
        )
        for channel in _CHANNELS
    )
    audit = FieldEvidenceAudit(
        root=_audit_string(data["root"], "root"),
        scanned_files=scanned_files,
        candidates=tuple(candidates),
        status=_audit_string(data["status"], "status"),
        channel_paths=channel_paths,
        reasons=_audit_strings(data["reasons"], "reasons"),
        scan_complete=scan_complete,
        coherent_candidate_paths=_audit_strings(
            data["coherent_candidate_paths"], "coherent_candidate_paths",
        ),
        max_files=max_files,
        gate_codes=(
            _audit_strings(data["gate_codes"], "gate_codes")
            if "gate_codes" in data else ()
        ),
    )
    if "diagnostic_counts" in data:
        supplied_counts = _audit_diagnostic_counts(
            data["diagnostic_counts"], "diagnostic_counts",
        )
        if supplied_counts != audit.diagnostic_counts:
            raise ValueError("diagnostic_counts must match candidate notes")
    if "collection_requirements" in data:
        _audit_collection_requirements(
            data["collection_requirements"], "collection_requirements", audit,
        )
    return audit


def _verify_audit_candidate_files(
    audit: FieldEvidenceAudit, *, require_digest: bool = False,
) -> None:
    root = Path(audit.root).resolve()
    for candidate in audit.candidates:
        path = Path(candidate.path).resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise ValueError(
                f"audit candidate path escapes the audit root: {path}"
            ) from error
        if not path.is_file():
            raise FileNotFoundError(f"audit candidate file is missing: {path}")
        if require_digest and candidate.sha256 is None:
            raise ValueError(
                f"audit candidate has no SHA-256 provenance: {path}"
            )
        if candidate.sha256 is not None:
            current = _sha256_file(path)
            if current != candidate.sha256:
                raise ValueError(
                    f"audit candidate changed since its provenance was recorded: {path}"
                )


def write_field_evidence_audit_json(
    audit: FieldEvidenceAudit, path: str | Path,
) -> FieldEvidenceAudit:
    """Write and re-read a strict, file-pinned field-evidence audit artifact."""
    if not isinstance(audit, FieldEvidenceAudit):
        raise TypeError("audit must be a FieldEvidenceAudit")
    output = Path(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing field evidence audit: {output}")
    _verify_audit_candidate_files(audit, require_digest=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA,
        "input": {"root": audit.root, "max_files": audit.max_files},
        "field_evidence_audit": field_evidence_audit_record(audit),
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    try:
        with output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(rendered)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing field evidence audit: {output}"
        ) from error
    return read_field_evidence_audit_json(output)


def read_field_evidence_audit_json(path: str | Path) -> FieldEvidenceAudit:
    """Read a strict audit artifact and verify every pinned candidate file."""
    case_path = Path(path)
    payload = json.loads(
        case_path.read_text(encoding="utf-8-sig"),
        object_pairs_hook=_json_object_pairs,
    )
    wrapper = _audit_mapping(payload, "field evidence audit execution")
    _audit_keys(
        wrapper,
        "field evidence audit execution",
        required={"schema", "input", "field_evidence_audit"},
    )
    if wrapper["schema"] != FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA:
        raise ValueError("schema must be " + repr(FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA))
    input_data = _audit_mapping(wrapper["input"], "input")
    _audit_keys(input_data, "input", required={"root"}, optional={"max_files"})
    audit = field_evidence_audit_from_mapping(wrapper["field_evidence_audit"])
    input_root = Path(_audit_string(input_data["root"], "input.root")).resolve()
    if input_root != Path(audit.root).resolve():
        raise ValueError("input.root must match field_evidence_audit.root")
    if "max_files" in input_data:
        input_max_files = input_data["max_files"]
        if (
            isinstance(input_max_files, bool)
            or not isinstance(input_max_files, int)
            or input_max_files <= 0
        ):
            raise ValueError("input.max_files must be a positive integer")
        if input_max_files != audit.max_files:
            raise ValueError("input.max_files must match field_evidence_audit.max_files")
    _verify_audit_candidate_files(audit, require_digest=True)
    return audit


__all__ = [
    "FIELD_EVIDENCE_AUDIT_SCHEMA", "FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA",
    "FieldEvidenceCandidate", "FieldEvidenceAudit", "audit_field_evidence",
    "field_evidence_audit_record", "field_evidence_audit_from_mapping",
    "write_field_evidence_audit_json", "read_field_evidence_audit_json",
]
