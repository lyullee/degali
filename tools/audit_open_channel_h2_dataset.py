"""Audit the public open-channel H2 source/sensor timing dataset.

The dataset is useful evidence for hydrogen source and detector timing, but it
does not contain the weather/geometry/common-clock package needed for an LH2
dispersion qualification.  This tool deliberately reports that boundary and
never promotes the runs to atmospheric-LH2 validation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from degali.validation.open_channel_h2 import OpenChannelHydrogenRun, read_open_channel_h2_csv, threshold_window


SCHEMA = "degali.open-channel-h2-audit.v1"
DATASET_DOI = "10.23642/usn.26117989.v2"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_stats(values: np.ndarray) -> tuple[float, float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (math.nan, math.nan, math.nan)
    return float(np.min(finite)), float(np.mean(finite)), float(np.max(finite))


def _integrate_trapezoid(values: np.ndarray, times: np.ndarray) -> float:
    """Integrate without depending on NumPy's removed ``np.trapz`` alias."""

    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    if values.ndim != 1 or times.shape != values.shape or values.size < 2:
        raise ValueError("trapezoid inputs must be matching one-dimensional histories")
    return float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(times)))


@dataclass(frozen=True)
class OpenChannelRunSummary:
    """Machine-readable, non-promoting summary of one local CSV."""

    file: str
    sha256: str
    dataset_doi: str
    flow_column: str
    sensor_count: int
    source_clock_start_s: float
    source_clock_end_s: float
    sensor_clock_start_s: float
    sensor_clock_end_s: float
    source_clock_duration_s: float
    sensor_clock_duration_s: float
    source_sensor_start_delta_s: float
    source_rate_min_g_s: float
    source_rate_mean_g_s: float
    source_rate_max_g_s: float
    negative_source_rate_count: int
    integrated_source_mass_g: float
    threshold_percent: float
    threshold_sensor_count: int
    threshold_arrival_min_s: float | None
    threshold_arrival_median_s: float | None
    threshold_arrival_max_s: float | None
    threshold_peak_min_percent: float | None
    threshold_peak_median_percent: float | None
    threshold_peak_max_percent: float | None
    validation_scope: str
    quantitative_lh2_pool_validation_allowed: bool


def _summary(path: Path, *, threshold_percent: float) -> OpenChannelRunSummary:
    run: OpenChannelHydrogenRun = read_open_channel_h2_csv(path)
    flow_min, flow_mean, flow_max = _finite_stats(run.mass_flow_g_s)
    integrated_mass = _integrate_trapezoid(run.mass_flow_g_s, run.flow_time_s)
    windows = [
        threshold_window(run, sensor, threshold_percent=threshold_percent)
        for sensor in run.sensor_names
    ]
    windows = [window for window in windows if window is not None]
    arrivals = np.asarray([window.arrival_s for window in windows], dtype=float)
    peaks = np.asarray([window.peak_percent for window in windows], dtype=float)

    def optional_stats(values: np.ndarray) -> tuple[float | None, float | None, float | None]:
        if values.size == 0:
            return None, None, None
        return float(np.min(values)), float(np.median(values)), float(np.max(values))

    arrival_min, arrival_median, arrival_max = optional_stats(arrivals)
    peak_min, peak_median, peak_max = optional_stats(peaks)
    return OpenChannelRunSummary(
        file=str(path.resolve()),
        sha256=_sha256(path),
        dataset_doi=run.dataset_doi,
        flow_column=run.flow_column,
        sensor_count=len(run.sensor_names),
        source_clock_start_s=float(run.flow_time_s[0]),
        source_clock_end_s=float(run.flow_time_s[-1]),
        sensor_clock_start_s=float(run.sensor_time_s[0]),
        sensor_clock_end_s=float(run.sensor_time_s[-1]),
        source_clock_duration_s=float(run.flow_time_s[-1] - run.flow_time_s[0]),
        sensor_clock_duration_s=float(run.sensor_time_s[-1] - run.sensor_time_s[0]),
        source_sensor_start_delta_s=float(run.sensor_time_s[0] - run.flow_time_s[0]),
        source_rate_min_g_s=flow_min,
        source_rate_mean_g_s=flow_mean,
        source_rate_max_g_s=flow_max,
        negative_source_rate_count=int(np.count_nonzero(run.mass_flow_g_s < 0.0)),
        integrated_source_mass_g=integrated_mass,
        threshold_percent=float(threshold_percent),
        threshold_sensor_count=len(windows),
        threshold_arrival_min_s=arrival_min,
        threshold_arrival_median_s=arrival_median,
        threshold_arrival_max_s=arrival_max,
        threshold_peak_min_percent=peak_min,
        threshold_peak_median_percent=peak_median,
        threshold_peak_max_percent=peak_max,
        validation_scope=run.validation_scope,
        quantitative_lh2_pool_validation_allowed=run.quantitative_lh2_pool_validation_allowed,
    )


def audit_open_channel_h2_dataset(
    root: str | Path,
    *,
    threshold_percent: float = 4.0,
) -> dict[str, Any]:
    """Audit extracted CSVs without promoting them to field validation."""

    if not math.isfinite(threshold_percent) or threshold_percent <= 0.0:
        raise ValueError("threshold_percent must be a finite positive value")
    directory = Path(root)
    if not directory.is_dir():
        return {
            "schema": SCHEMA,
            "status": "not_available",
            "dataset_doi": DATASET_DOI,
            "root": str(directory.resolve()),
            "files": [],
            "errors": [],
            "promotion_allowed": False,
            "gate_codes": ["dataset_not_available", "not_lh2_pool_validation"],
        }
    paths = sorted(directory.rglob("*.csv"))
    summaries: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    for path in paths:
        try:
            summaries.append(asdict(_summary(path, threshold_percent=threshold_percent)))
        except (OSError, ValueError, TypeError) as error:
            errors.append({"file": str(path.resolve()), "error": str(error)})
    status = "boundary_only" if paths and not errors else "withheld"
    if not paths:
        status = "not_available"
    return {
        "schema": SCHEMA,
        "status": status,
        "dataset_doi": DATASET_DOI,
        "root": str(directory.resolve()),
        "file_count": len(paths),
        "parsed_file_count": len(summaries),
        "files": summaries,
        "errors": errors,
        "promotion_allowed": False,
        "quantitative_lh2_pool_validation_allowed": False,
        "scope": "hydrogen_source_and_sensor_timing_in_channel_only",
        "gate_codes": [
            "source_sensor_timing_boundary",
            "not_lh2_pool_validation",
            "weather_geometry_common_clock_missing",
        ],
    }


def write_audit_json(record: dict[str, Any], path: str | Path) -> None:
    """Write an audit record without replacing an existing artifact."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(record, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="directory containing extracted public CSV files")
    parser.add_argument("--threshold-percent", type=float, default=4.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    record = audit_open_channel_h2_dataset(args.root, threshold_percent=args.threshold_percent)
    rendered = json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False)
    if args.output is None:
        print(rendered)
    else:
        write_audit_json(record, args.output)
        print(f"wrote {args.output} status={record['status']} files={record.get('file_count', 0)}")
    return 0 if record["status"] in {"boundary_only", "not_available"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
