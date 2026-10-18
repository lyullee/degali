"""Evidence-gated intake for the public USN/FFI obstructed H2 channel cases.

This adapter connects the article's fixed sensor coordinates to the raw FFI
source/sensor histories while preserving their independent clocks. It does not
turn the confined channel into an outdoor LH2 case and withholds prediction
scoring until obstacle geometry and the event clock are resolved.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .open_channel_h2 import OpenChannelHydrogenRun, read_open_channel_h2_csv
from ..addons.site_geometry import AxisAlignedCuboid


SCHEMA = "degali.ffi-open-channel-obstacle-benchmark.v1"


@dataclass(frozen=True)
class ChannelSensorReceptor:
    """One article-registered concentration sensor in channel coordinates."""

    sensor_id: str
    position_m: tuple[float, float, float]
    placement: str
    along_centerline: bool


@dataclass(frozen=True)
class FFIObstacleCase:
    """One article case and all raw experimental repeats mapped to it."""

    article_case: str
    trials: tuple[str, ...]
    obstacle: bool
    article_les_source_rate_g_s: float
    runs: tuple[OpenChannelHydrogenRun, ...]


@dataclass(frozen=True)
class OpenChannelH2ObstacleBenchmark:
    """Hash-verified observations plus explicit promotion boundaries."""

    sensor_receptors: tuple[ChannelSensorReceptor, ...]
    cases: tuple[FFIObstacleCase, ...]
    geometry_quantitative_allowed: bool
    obstacle_h2_validation_allowed: bool
    degali_score_allowed: bool
    gate_codes: tuple[str, ...]
    registry_path: Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_number(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number")
    return number


def load_open_channel_h2_obstacle_benchmark(
    registry_path: str | Path,
    extracted_data_root: str | Path,
) -> OpenChannelH2ObstacleBenchmark:
    """Load the registered FFI sensor map and verify every mapped raw CSV.

    The loader is intentionally strict: ambiguous filenames, changed file
    hashes, sensor-ID mismatches, invalid coordinates, or malformed case
    crosswalks stop intake before a model comparison can be assembled.
    """

    registry_file = Path(registry_path)
    data_root = Path(extracted_data_root)
    try:
        record = json.loads(registry_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"unable to read FFI obstacle registry: {error}") from error
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError(f"FFI obstacle registry schema must be {SCHEMA!r}")

    facility = record.get("facility")
    if not isinstance(facility, dict):
        raise ValueError("FFI obstacle registry is missing facility geometry")
    dimensions = facility.get("channel_dimensions_m")
    if not isinstance(dimensions, dict):
        raise ValueError("FFI obstacle registry is missing channel dimensions")
    bounds = (
        _finite_number(dimensions.get("x_width"), "channel x width"),
        _finite_number(dimensions.get("y_length"), "channel y length"),
        _finite_number(dimensions.get("z_height"), "channel z height"),
    )
    if any(value <= 0.0 for value in bounds):
        raise ValueError("FFI channel dimensions must be positive")

    sensor_rows = record.get("sensors")
    if not isinstance(sensor_rows, list) or len(sensor_rows) != 29:
        raise ValueError("FFI obstacle registry must contain all 29 article sensors")
    sensors: list[ChannelSensorReceptor] = []
    for row in sensor_rows:
        if not isinstance(row, dict):
            raise ValueError("FFI sensor entries must be objects")
        sensor_id = row.get("id")
        placement = row.get("placement")
        centerline = row.get("centerline")
        if not isinstance(sensor_id, str) or not sensor_id.strip():
            raise ValueError("FFI sensor id must be non-empty text")
        if not isinstance(placement, str) or not placement.strip():
            raise ValueError(f"FFI sensor {sensor_id!r} has no placement class")
        if not isinstance(centerline, bool):
            raise ValueError(f"FFI sensor {sensor_id!r} centerline flag must be boolean")
        position = tuple(
            _finite_number(row.get(key), f"{sensor_id} {key}")
            for key in ("x_m", "y_m", "z_m")
        )
        if any(value < 0.0 or value > bound for value, bound in zip(position, bounds)):
            raise ValueError(f"FFI sensor {sensor_id!r} lies outside the channel")
        sensors.append(ChannelSensorReceptor(sensor_id, position, placement, centerline))
    sensor_ids = {sensor.sensor_id for sensor in sensors}
    if len(sensor_ids) != len(sensors):
        raise ValueError("FFI sensor registry contains duplicate sensor IDs")
    expected_ids = {f"sensor {number}" for number in range(1, 31) if number != 16}
    if sensor_ids != expected_ids:
        raise ValueError("FFI sensor registry IDs do not match the published 29-sensor layout")

    cases_raw = record.get("case_crosswalk")
    if not isinstance(cases_raw, list) or not cases_raw:
        raise ValueError("FFI obstacle registry has no article-case crosswalk")
    cases: list[FFIObstacleCase] = []
    seen_paths: set[Path] = set()
    seen_cases: set[str] = set()
    for case in cases_raw:
        if not isinstance(case, dict):
            raise ValueError("FFI case crosswalk entries must be objects")
        article_case = case.get("article_case")
        trials = case.get("article_trials")
        filenames = case.get("dataset_files")
        hashes = case.get("dataset_csv_sha256")
        obstacle = case.get("obstacle")
        rates = _finite_number(
            case.get("article_les_source_rate_g_s"),
            f"{article_case} LES source rate",
        )
        if not isinstance(article_case, str) or not article_case.strip():
            raise ValueError("FFI article-case ID must be non-empty text")
        if article_case in seen_cases:
            raise ValueError(f"duplicate FFI article case {article_case!r}")
        seen_cases.add(article_case)
        if (
            not isinstance(trials, list) or not trials
            or not all(isinstance(item, str) and item.strip() for item in trials)
        ):
            raise ValueError(f"FFI case {article_case!r} has invalid article trials")
        if not isinstance(filenames, list) or len(filenames) != len(trials):
            raise ValueError(f"FFI case {article_case!r} trial/file counts differ")
        if not isinstance(hashes, list) or len(hashes) != len(filenames):
            raise ValueError(f"FFI case {article_case!r} file/hash counts differ")
        if not isinstance(obstacle, bool):
            raise ValueError(f"FFI case {article_case!r} obstacle flag must be boolean")
        if rates <= 0.0:
            raise ValueError(f"FFI case {article_case!r} LES source rate must be positive")

        runs: list[OpenChannelHydrogenRun] = []
        for filename, expected_hash in zip(filenames, hashes):
            if not isinstance(filename, str) or Path(filename).name != filename:
                raise ValueError(f"FFI dataset filename must be a basename: {filename!r}")
            if (
                not isinstance(expected_hash, str)
                or len(expected_hash) != 64
                or any(char not in "0123456789abcdef" for char in expected_hash)
            ):
                raise ValueError(f"FFI CSV SHA-256 is malformed for {filename!r}")
            matches = tuple(data_root.rglob(filename))
            if len(matches) != 1:
                raise ValueError(
                    f"expected exactly one FFI CSV named {filename!r}; found {len(matches)}"
                )
            source_path = matches[0].resolve()
            if source_path in seen_paths:
                raise ValueError(f"FFI raw CSV is mapped more than once: {filename!r}")
            seen_paths.add(source_path)
            actual_hash = _sha256(source_path)
            if actual_hash != expected_hash:
                raise ValueError(f"FFI CSV SHA-256 mismatch for {filename!r}")
            run = read_open_channel_h2_csv(source_path)
            if set(run.sensor_names) != sensor_ids:
                raise ValueError(
                    f"FFI CSV {filename!r} sensor IDs differ from the registered layout"
                )
            runs.append(run)
        cases.append(FFIObstacleCase(
            article_case=article_case,
            trials=tuple(trials),
            obstacle=obstacle,
            article_les_source_rate_g_s=rates,
            runs=tuple(runs),
        ))

    obstacle = record.get("obstacle")
    if not isinstance(obstacle, dict):
        raise ValueError("FFI obstacle registry is missing obstacle geometry evidence")
    geometry_allowed = obstacle.get("quantitative_geometry_allowed") is True
    promotion = record.get("promotion")
    if not isinstance(promotion, dict):
        raise ValueError("FFI obstacle registry is missing promotion gates")
    gate_codes = promotion.get("gate_codes")
    if not isinstance(gate_codes, list) or not all(
        isinstance(item, str) and item for item in gate_codes
    ):
        raise ValueError("FFI promotion gate codes must be non-empty text entries")
    obstacle_validation_allowed = (
        promotion.get("obstacle_h2_validation_allowed") is True
        and geometry_allowed
        and not gate_codes
    )
    score_allowed = (
        promotion.get("degali_score_available") is True
        and obstacle_validation_allowed
    )
    if not geometry_allowed and not gate_codes:
        raise ValueError("withheld FFI geometry requires explicit gate codes")

    return OpenChannelH2ObstacleBenchmark(
        sensor_receptors=tuple(sensors),
        cases=tuple(cases),
        geometry_quantitative_allowed=geometry_allowed,
        obstacle_h2_validation_allowed=obstacle_validation_allowed,
        degali_score_allowed=score_allowed,
        gate_codes=tuple(gate_codes),
        registry_path=registry_file.resolve(),
    )


def conditional_ffi_o1_crossbeam_assembly(
    *,
    channel_width_m: float = 0.9,
    channel_height_m: float = 0.8,
    y_start_m: float = 2.1,
    assembly_length_m: float = 1.0,
    reported_blockage_ratio: float = 0.368,
    crossbeam_count: int = 2,
    edge_margin_m: float = 0.05,
) -> tuple[AxisAlignedCuboid, ...]:
    """Construct an explicit *conditional* O1 two-top-crossbeam surrogate.

    The article identifies O1 as two top cross beams and reports a 36.8%
    frontal blockage, but does not publish beam section dimensions or the
    leading-edge datum.  This constructor preserves the reported blockage in
    the channel cross-section and spreads the two beams over the reported
    one-metre assembly length.  It is deliberately named conditional so a
    caller cannot mistake the geometry for machine-registered facility CAD.
    """
    values = (
        channel_width_m, channel_height_m, y_start_m, assembly_length_m,
        reported_blockage_ratio, edge_margin_m,
    )
    if any(not math.isfinite(float(value)) or float(value) <= 0.0 for value in values):
        raise ValueError("crossbeam dimensions and blockage must be positive and finite")
    if not isinstance(crossbeam_count, int) or isinstance(crossbeam_count, bool):
        raise TypeError("crossbeam_count must be an integer")
    if crossbeam_count < 1:
        raise ValueError("crossbeam_count must be positive")
    width = channel_width_m - 2.0 * edge_margin_m
    if width <= 0.0 or edge_margin_m >= channel_width_m / 2.0:
        raise ValueError("edge margin leaves no transverse beam width")
    if reported_blockage_ratio >= 1.0:
        raise ValueError("reported blockage ratio must be below one")
    beam_height = (
        reported_blockage_ratio * channel_width_m * channel_height_m
        / (crossbeam_count * width)
    )
    if beam_height >= channel_height_m:
        raise ValueError("reported blockage cannot be represented by top beams")
    segment = assembly_length_m / crossbeam_count
    # Keep a sub-cell clearance at the top because the solver's embedded
    # obstacle contract requires strict interior bounds.  The reported
    # blockage is preserved by retaining the nominal channel cross-section in
    # the beam-height calculation.
    z_top = channel_height_m - edge_margin_m
    z_min = z_top - beam_height
    beams: list[AxisAlignedCuboid] = []
    for index in range(crossbeam_count):
        lower_y = y_start_m + index * segment
        upper_y = y_start_m + (index + 1) * segment
        beams.append(AxisAlignedCuboid(
            edge_margin_m, channel_width_m - edge_margin_m,
            lower_y, upper_y, z_min, z_top,
            f"FFI O1 conditional top crossbeam {index + 1}/{crossbeam_count}",
        ))
    return tuple(beams)


def load_registered_ffi_obstacle_geometry(
    path: str | Path,
    *,
    article_case: str = "O1",
    channel_dimensions_m: tuple[float, float, float] = (0.9, 5.8, 0.8),
) -> tuple[AxisAlignedCuboid, ...]:
    """Load a machine-registered FFI obstacle assembly without guessing.

    The input is intentionally narrower than a CAD reader: an upstream
    geometry-registration step must export an exact axis-aligned cuboid
    decomposition in the article frame.  Conditional article dimensions,
    inferred blockage or a non-registered status are rejected.  This keeps
    the resolved solver from silently substituting the two-beam surrogate
    when a real facility geometry is requested.

    Expected schema::

        {
          "schema": "degali.ffi-obstacle-geometry.v1",
          "registration_status": "machine_registered_exact",
          "article_case": "O1",
          "coordinate_frame": "ffi-open-channel-v1",
          "obstacles": [{"shape": "axis_aligned_cuboid", ...bounds...}]
        }
    """
    geometry_path = Path(path).resolve()
    if not geometry_path.is_file():
        raise FileNotFoundError(geometry_path)
    try:
        record = json.loads(geometry_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("registered FFI geometry must be valid UTF-8 JSON") from error
    if not isinstance(record, dict):
        raise ValueError("registered FFI geometry must be a JSON object")
    if record.get("schema") != "degali.ffi-obstacle-geometry.v1":
        raise ValueError("registered FFI geometry schema is unsupported")
    if record.get("registration_status") != "machine_registered_exact":
        raise ValueError(
            "quantitative FFI geometry requires registration_status="
            "machine_registered_exact"
        )
    if record.get("article_case") != article_case:
        raise ValueError(
            f"registered FFI geometry article_case must be {article_case!r}"
        )
    if record.get("coordinate_frame") != "ffi-open-channel-v1":
        raise ValueError("registered FFI geometry coordinate_frame is unsupported")
    dimensions = tuple(_finite_number(value, "channel dimension") for value in channel_dimensions_m)
    if len(dimensions) != 3 or any(value <= 0.0 for value in dimensions):
        raise ValueError("channel_dimensions_m must be three positive values")
    obstacles = record.get("obstacles")
    if not isinstance(obstacles, list) or not obstacles:
        raise ValueError("registered FFI geometry requires a non-empty obstacles list")
    labels: set[str] = set()
    loaded: list[AxisAlignedCuboid] = []
    for index, item in enumerate(obstacles):
        if not isinstance(item, dict) or item.get("shape") != "axis_aligned_cuboid":
            raise ValueError(
                f"registered FFI obstacle {index} must be an axis_aligned_cuboid"
            )
        label = item.get("label", f"registered obstacle {index + 1}")
        if not isinstance(label, str) or not label.strip() or label in labels:
            raise ValueError("registered FFI obstacle labels must be unique text")
        labels.add(label)
        bounds = tuple(
            _finite_number(item.get(name), f"obstacle {index} {name}")
            for name in (
                "x_min_m", "x_max_m", "y_min_m", "y_max_m", "z_min_m", "z_max_m"
            )
        )
        obstacle = AxisAlignedCuboid(*bounds, label=label)
        width, length, height = dimensions
        tolerance = 1.0e-10
        if not (
            -tolerance <= obstacle.x_min_m < obstacle.x_max_m <= width + tolerance
            and -tolerance <= obstacle.y_min_m < obstacle.y_max_m <= length + tolerance
            and -tolerance <= obstacle.z_min_m < obstacle.z_max_m <= height + tolerance
        ):
            raise ValueError(
                f"registered FFI obstacle {label!r} lies outside channel bounds"
            )
        loaded.append(obstacle)
    return tuple(loaded)


def conditional_ffi_o1_thin_top_crossbeam_assembly(
    *,
    channel_width_m: float = 0.9,
    channel_height_m: float = 0.8,
    y_start_m: float = 2.1,
    assembly_length_m: float = 1.0,
    reported_blockage_ratio: float = 0.368,
    crossbeam_count: int = 2,
    beam_height_m: float = 0.15,
    top_clearance_m: float = 0.01,
) -> tuple[AxisAlignedCuboid, ...]:
    """Construct a thin upper-beam O1 sensitivity surrogate.

    The article reports two *top* cross beams and a 0.368 blockage ratio but
    does not publish their section or lateral edge coordinates.  This
    sensitivity geometry assumes the beams span nearly the full channel width
    and uses the reported blockage to infer their vertical section.  With the
    default 0.15 m height, the implied span is 0.8832 m and the two beams sit
    immediately below the ceiling.  It is a documented shape hypothesis, not
    machine-registered facility CAD.
    """
    values = (
        channel_width_m, channel_height_m, y_start_m, assembly_length_m,
        reported_blockage_ratio, beam_height_m, top_clearance_m,
    )
    if any(not math.isfinite(float(value)) or float(value) <= 0.0 for value in values):
        raise ValueError("thin crossbeam dimensions and blockage must be positive and finite")
    if not isinstance(crossbeam_count, int) or isinstance(crossbeam_count, bool):
        raise TypeError("crossbeam_count must be an integer")
    if crossbeam_count < 1:
        raise ValueError("crossbeam_count must be positive")
    if reported_blockage_ratio >= 1.0:
        raise ValueError("reported blockage ratio must be below one")
    if beam_height_m + top_clearance_m >= channel_height_m:
        raise ValueError("thin top beams must leave positive channel clearance")
    projected_area = reported_blockage_ratio * channel_width_m * channel_height_m
    span_width = projected_area / (crossbeam_count * beam_height_m)
    if not 0.0 < span_width < channel_width_m:
        raise ValueError("beam height does not permit an interior span for the reported blockage")
    edge_margin = 0.5 * (channel_width_m - span_width)
    segment = assembly_length_m / crossbeam_count
    z_top = channel_height_m - top_clearance_m
    z_min = z_top - beam_height_m
    beams: list[AxisAlignedCuboid] = []
    for index in range(crossbeam_count):
        lower_y = y_start_m + index * segment
        upper_y = y_start_m + (index + 1) * segment
        beams.append(AxisAlignedCuboid(
            edge_margin, channel_width_m - edge_margin,
            lower_y, upper_y, z_min, z_top,
            f"FFI O1 conditional thin top crossbeam {index + 1}/{crossbeam_count}",
        ))
    return tuple(beams)


__all__ = [
    "ChannelSensorReceptor", "FFIObstacleCase", "OpenChannelH2ObstacleBenchmark",
    "conditional_ffi_o1_crossbeam_assembly",
    "conditional_ffi_o1_thin_top_crossbeam_assembly",
    "load_registered_ffi_obstacle_geometry",
    "load_open_channel_h2_obstacle_benchmark",
]
