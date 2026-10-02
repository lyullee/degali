"""Local intake for the public ELVHYS WP4.2 confined cryogenic-H2 tests.

The HSE campaign at DOI 10.18710/JXJP0H supplies independently useful
concentration, temperature, pressure and ventilation histories for a 1 m3
transfer connection space. It does *not* supply an independently measured
hydrogen mass-flow boundary. This module makes the observations auditable
while preventing their use as an outdoor LH2-pool plume score.

No ELVHYS file or derived time series is distributed with DEGALI.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np


_STREAMS = ("CONC", "TEMP", "PRES", "FLMT", "MISC")


@dataclass(frozen=True)
class ElvhysTimeSeries:
    """One local ELVHYS stream, retaining its published clock and channels."""

    source_path: Path
    time_s: np.ndarray
    channels: dict[str, np.ndarray]
    stream: str

    @property
    def channel_names(self) -> tuple[str, ...]:
        return tuple(self.channels)


@dataclass(frozen=True)
class ElvhysTcsTest:
    """The available time series for one declared ELVHYS test number."""

    test_number: int
    streams: dict[str, ElvhysTimeSeries]
    dataset_doi: str = "10.18710/JXJP0H"
    enclosure_volume_m3: float = 1.0
    source_geometry: str = "confined_transfer_connection_space"

    @property
    def hydrogen_mass_flow_measured(self) -> bool:
        """The published FLMT stream is ventilation flow, not H2 mass flow."""
        return False

    @property
    def quantitative_outdoor_lh2_pool_validation_allowed(self) -> bool:
        """The confined, no-flow-boundary data cannot score an outdoor pool."""
        return False

    @property
    def validation_scope(self) -> str:
        return "confined_cryogenic_h2_observation_only"

    def stream(self, name: str) -> ElvhysTimeSeries:
        try:
            return self.streams[name.upper()]
        except KeyError as error:
            raise ValueError(f"ELVHYS test {self.test_number} lacks stream {name!r}") from error


def _float(value: str | None, *, column: str, path: Path) -> float:
    if value is None or not value.strip():
        return math.nan
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(f"invalid ELVHYS value in {column!r} of {path.name}") from error


def read_elvhys_stream(path: str | Path, *, stream: str) -> ElvhysTimeSeries:
    """Read one user-supplied public ELVHYS CSV stream.

    Only the explicit ``Time`` column defines time. Missing channel readings
    remain NaN; no baseline subtraction, delay correction or resampling is
    performed in this intake layer.
    """
    source = Path(path)
    declared = stream.upper()
    if declared not in _STREAMS:
        raise ValueError(f"unknown ELVHYS stream {stream!r}")
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        if "Time" not in fields:
            raise ValueError(f"ELVHYS {source.name} is missing the 'Time' column")
        channel_names = tuple(name for name in fields if name != "Time")
        if not channel_names:
            raise ValueError(f"ELVHYS {source.name} has no measurement channels")
        rows = list(reader)

    time_s = np.asarray(
        [_float(row.get("Time"), column="Time", path=source) for row in rows],
        dtype=float,
    )
    if time_s.size < 2 or not np.all(np.isfinite(time_s)):
        raise ValueError(f"ELVHYS {source.name} needs at least two finite timestamps")
    if np.any(np.diff(time_s) < 0.0):
        raise ValueError(f"ELVHYS {source.name} timestamps must be nondecreasing")
    channels = {
        name: np.asarray(
            [_float(row.get(name), column=name, path=source) for row in rows], dtype=float,
        )
        for name in channel_names
    }
    return ElvhysTimeSeries(
        source_path=source, time_s=time_s, channels=channels, stream=declared,
    )


def _find_stream_file(root: Path, test_number: int, stream: str) -> Path | None:
    prefix = f"ELE402HSE{test_number:03d}"
    files = sorted(
        path for path in root.rglob(f"{prefix}*{stream}*.csv") if path.is_file()
    )
    if len(files) > 1:
        raise ValueError(
            f"ELVHYS test {test_number} has multiple {stream} files in {root}"
        )
    return files[0] if files else None


def read_elvhys_test(
    directory: str | Path,
    test_number: int,
    *,
    required_streams: tuple[str, ...] = ("CONC", "TEMP", "PRES", "FLMT"),
) -> ElvhysTcsTest:
    """Read a declared ELVHYS test from a local public-data directory.

    The test number is not inferred from timing or signal shape. The caller
    declares it, and each selected stream must have the official filename
    prefix. ``MISC`` is optional because some analyses need no weather data.
    """
    root = Path(directory)
    if not isinstance(test_number, int) or test_number <= 0:
        raise ValueError("ELVHYS test_number must be a positive integer")
    normalized = tuple(name.upper() for name in required_streams)
    unknown = set(normalized) - set(_STREAMS)
    if unknown:
        raise ValueError("unknown required ELVHYS streams: " + ", ".join(sorted(unknown)))
    streams: dict[str, ElvhysTimeSeries] = {}
    for name in _STREAMS:
        file = _find_stream_file(root, test_number, name)
        if file is None:
            if name in normalized:
                raise ValueError(f"ELVHYS test {test_number} lacks required {name} CSV")
            continue
        streams[name] = read_elvhys_stream(file, stream=name)
    return ElvhysTcsTest(test_number=test_number, streams=streams)


__all__ = [
    "ElvhysTcsTest", "ElvhysTimeSeries", "read_elvhys_stream", "read_elvhys_test",
]
