"""Read the public H2SAFE indoor helium-surrogate experiments.

The H2SAFE record (DOI ``10.7799/17118570``) contains five controlled
helium releases in two full-scale laboratories.  The downloaded package
provides a one-second relative clock, sensor coordinates, concentration
signals, STEP/STL geometry and a separate HVAC/release document.

The package does *not* declare the concentration unit in the CSV headers,
does not connect the relative CSV clock to the documented release clock, and
does not include event calibration certificates.  This reader therefore
preserves the native values and missing markers without converting them to
vol.%.  H2SAFE is component evidence for indoor geometry, forced ventilation
and transient sensor response; it is not direct hydrogen or LH2 quantitative
validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import math
from pathlib import Path
import re
from typing import Literal, Mapping

import numpy as np


DATASET_DOI = "10.7799/17118570"
DATASET_URL = "https://data.nlr.gov/submissions/330"
SPECIES = "helium"
CONCENTRATION_UNIT_STATUS = "not_declared_in_downloaded_source"
SOURCE_CLOCK_ALIGNMENT_STATUS = "not_declared"
SENSOR_CALIBRATION_STATUS = "missing"
PUBLISHED_VERTICAL_AXIS = "Y"
_STEP_PRODUCT = re.compile(r"PRODUCT\s*\(\s*'((?:''|[^'])*)'", re.IGNORECASE)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class H2SafeReleaseMetadata:
    """Release information transcribed from the hash-pinned HVAC document."""

    lab_id: str
    test_id: str
    test_date: str | None
    nozzle_diameter_mm: float
    orientation: str
    pressure_psi: float
    release_rate_lower: float
    release_rate_upper: float
    release_rate_unit: str
    release_duration_s: float
    release_duration_basis: str
    documented_clock_window: str | None
    release_position_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if self.lab_id not in {"Lab1_Small", "Lab2_Large"}:
            raise ValueError("H2SAFE lab_id is unsupported")
        if not self.test_id.startswith("test"):
            raise ValueError("H2SAFE test_id must start with test")
        finite = (
            self.nozzle_diameter_mm,
            self.pressure_psi,
            self.release_rate_lower,
            self.release_rate_upper,
            self.release_duration_s,
            *self.release_position_m,
        )
        if not all(math.isfinite(float(value)) for value in finite):
            raise ValueError("H2SAFE release metadata must be finite")
        if self.nozzle_diameter_mm <= 0.0 or self.pressure_psi <= 0.0:
            raise ValueError("H2SAFE nozzle diameter and pressure must be positive")
        if not 0.0 < self.release_rate_lower <= self.release_rate_upper:
            raise ValueError("H2SAFE release-rate bounds are invalid")
        if self.release_duration_s <= 0.0:
            raise ValueError("H2SAFE release duration must be positive")

    @property
    def solver_position_m(self) -> tuple[float, float, float]:
        return h2safe_raw_xyz_to_solver_xyz(self.release_position_m)


@dataclass(frozen=True)
class H2SafeVentBoundary:
    vent_id: str
    role: Literal["inlet", "outlet"]
    area_m2: float
    speed_m_s: float
    angle_text: str
    location_text: str | None = None

    def __post_init__(self) -> None:
        if not self.vent_id.strip() or self.role not in {"inlet", "outlet"}:
            raise ValueError("H2SAFE vent identity/role is invalid")
        if not math.isfinite(self.area_m2) or self.area_m2 <= 0.0:
            raise ValueError("H2SAFE vent area must be finite and positive")
        if not math.isfinite(self.speed_m_s) or self.speed_m_s < 0.0:
            raise ValueError("H2SAFE vent speed must be finite and non-negative")

    @property
    def volumetric_flow_m3_s(self) -> float:
        return self.area_m2 * self.speed_m_s


@dataclass(frozen=True)
class H2SafeRecirculationFan:
    fan_id: str
    speed_m_s: float
    pressure_drop_pa: float
    inlet_area_m2: float
    outlet_area_m2: float

    def __post_init__(self) -> None:
        values = (
            self.speed_m_s, self.pressure_drop_pa,
            self.inlet_area_m2, self.outlet_area_m2,
        )
        if not self.fan_id.strip() or not all(math.isfinite(value) for value in values):
            raise ValueError("H2SAFE recirculation fan metadata is invalid")
        if any(value <= 0.0 for value in values):
            raise ValueError("H2SAFE recirculation fan values must be positive")


_LAB1_SOURCE = (-5.186701676, 1.580460097, -4.322127)
_LAB2_SOURCE = (-11.63093726, 1.566073152, 15.86309373)

RELEASES: Mapping[tuple[str, str], H2SafeReleaseMetadata] = {
    ("Lab1_Small", "test1"): H2SafeReleaseMetadata(
        lab_id="Lab1_Small", test_id="test1", test_date=None,
        nozzle_diameter_mm=2.0, orientation="vertical", pressure_psi=40.0,
        release_rate_lower=100.0, release_rate_upper=100.0,
        release_rate_unit="SLM", release_duration_s=1389.0,
        release_duration_basis="declared_seconds",
        documented_clock_window=None, release_position_m=_LAB1_SOURCE,
    ),
    ("Lab2_Large", "test1"): H2SafeReleaseMetadata(
        lab_id="Lab2_Large", test_id="test1", test_date="2025-05-09",
        nozzle_diameter_mm=0.5, orientation="horizontal", pressure_psi=220.0,
        release_rate_lower=92.0, release_rate_upper=100.0,
        release_rate_unit="L/min", release_duration_s=600.0,
        release_duration_basis="nominal_from_minute_resolution_clock_text",
        documented_clock_window="12:31-12:41 PM", release_position_m=_LAB2_SOURCE,
    ),
    ("Lab2_Large", "test2"): H2SafeReleaseMetadata(
        lab_id="Lab2_Large", test_id="test2", test_date="2025-05-09",
        nozzle_diameter_mm=0.5, orientation="vertical", pressure_psi=220.0,
        release_rate_lower=92.0, release_rate_upper=100.0,
        release_rate_unit="L/min", release_duration_s=1740.0,
        release_duration_basis="nominal_from_minute_resolution_clock_text",
        documented_clock_window="2:11-2:40 PM", release_position_m=_LAB2_SOURCE,
    ),
    ("Lab2_Large", "test3"): H2SafeReleaseMetadata(
        lab_id="Lab2_Large", test_id="test3", test_date="2025-05-13",
        nozzle_diameter_mm=1.0, orientation="vertical", pressure_psi=45.0,
        release_rate_lower=100.0, release_rate_upper=100.0,
        release_rate_unit="SLM", release_duration_s=2576.0,
        release_duration_basis="nominal_from_mixed_second_minute_clock_text",
        documented_clock_window="10:54:04 AM-11:37 AM",
        release_position_m=_LAB2_SOURCE,
    ),
    ("Lab2_Large", "test4"): H2SafeReleaseMetadata(
        lab_id="Lab2_Large", test_id="test4", test_date="2025-05-13",
        nozzle_diameter_mm=1.0, orientation="vertical", pressure_psi=430.0,
        release_rate_lower=680.0, release_rate_upper=700.0,
        release_rate_unit="SLM", release_duration_s=480.0,
        release_duration_basis="nominal_from_minute_resolution_clock_text",
        documented_clock_window="3:29-3:37 PM", release_position_m=_LAB2_SOURCE,
    ),
}

_IN2_TO_M2 = 0.00064516
LAB1_VENT_BOUNDARIES = (
    H2SafeVentBoundary("snorkel_outlet", "outlet", 0.007126, 2.59, "normal"),
    H2SafeVentBoundary("outlet_vent", "outlet", 0.264194, 6.05, "normal"),
    H2SafeVentBoundary("metal_duct_1_left", "inlet", 0.104839, 2.33, "normal"),
    H2SafeVentBoundary("metal_duct_1_right", "inlet", 0.104839, 1.86, "normal"),
    H2SafeVentBoundary("metal_duct_2_left", "inlet", 0.104839, 2.41, "18 deg up"),
    H2SafeVentBoundary("metal_duct_2_right", "inlet", 0.104839, 2.43, "18 deg up"),
    H2SafeVentBoundary("metal_duct_3_left", "inlet", 0.104839, 1.04, "normal"),
    H2SafeVentBoundary("metal_duct_3_right", "inlet", 0.104839, 2.16, "normal"),
    H2SafeVentBoundary("metal_duct_4_left", "inlet", 0.104839, 1.20, "normal"),
    H2SafeVentBoundary("metal_duct_4_right", "inlet", 0.104839, 2.00, "normal"),
)
LAB2_SUPPLY_BOUNDARIES = tuple(
    H2SafeVentBoundary(
        f"supply_duct_{number}", "inlet", 288.0 * _IN2_TO_M2, speed,
        "not_declared", location,
    )
    for number, speed, location in (
        (1, 1.1684, "door-office"),
        (2, 1.09728, None),
        (3, 1.5494, None),
        (4, 0.0, None),
        (5, 0.0, None),
        (6, 2.19964, None),
        (7, 0.0, None),
        (8, 1.6764, "garage"),
    )
)
LAB2_RECIRCULATION_FANS = (
    H2SafeRecirculationFan(
        "recirculation_fan_1", 11.8364, 85.8,
        80.0 * _IN2_TO_M2, 144.0 * _IN2_TO_M2,
    ),
    H2SafeRecirculationFan(
        "recirculation_fan_2", 14.478, 128.0,
        80.0 * _IN2_TO_M2, 144.0 * _IN2_TO_M2,
    ),
)


@dataclass(frozen=True)
class H2SafeSensor:
    sensor_id: str
    position_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not self.sensor_id.strip() or self.sensor_id != self.sensor_id.lower():
            raise ValueError("H2SAFE sensor IDs must be normalized lowercase")
        if len(self.position_m) != 3 or not all(
            math.isfinite(value) for value in self.position_m
        ):
            raise ValueError("H2SAFE sensor position must contain three finite values")

    @property
    def solver_position_m(self) -> tuple[float, float, float]:
        return h2safe_raw_xyz_to_solver_xyz(self.position_m)


def h2safe_raw_xyz_to_solver_xyz(
    position_m: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Map published X/Y(vertical)/Z coordinates to solver x/y/z(vertical)."""

    if len(position_m) != 3 or not all(math.isfinite(value) for value in position_m):
        raise ValueError("H2SAFE coordinate transform requires three finite values")
    raw_x, raw_vertical_y, raw_z = position_m
    return float(raw_x), float(raw_z), float(raw_vertical_y)


@dataclass(frozen=True)
class H2SafeRun:
    """One H2SAFE event in its native, explicitly uncalibrated representation."""

    source_path: Path
    source_sha256: str
    registry_path: Path
    registry_sha256: str
    release: H2SafeReleaseMetadata
    registry: tuple[H2SafeSensor, ...]
    sensor_ids: tuple[str, ...]
    time_s: np.ndarray
    concentration_raw: np.ndarray

    def __post_init__(self) -> None:
        if self.source_path.name.lower() != f"{self.release.test_id}.csv":
            raise ValueError("H2SAFE source filename and release metadata differ")
        if len(self.source_sha256) != 64 or len(self.registry_sha256) != 64:
            raise ValueError("H2SAFE SHA-256 digest is malformed")
        if self.time_s.ndim != 1 or self.time_s.size < 2:
            raise ValueError("H2SAFE time must be a one-dimensional series")
        if not np.all(np.isfinite(self.time_s)) or np.any(np.diff(self.time_s) <= 0.0):
            raise ValueError("H2SAFE time must be finite and strictly increasing")
        if self.concentration_raw.shape != (self.time_s.size, len(self.sensor_ids)):
            raise ValueError("H2SAFE concentration matrix shape is inconsistent")
        if len(set(self.sensor_ids)) != len(self.sensor_ids):
            raise ValueError("H2SAFE trace sensor IDs must be unique")
        registry_ids = {sensor.sensor_id for sensor in self.registry}
        missing = sorted(set(self.sensor_ids) - registry_ids)
        if missing:
            raise ValueError(
                "H2SAFE trace contains sensors absent from the registry: "
                + ", ".join(missing)
            )
        finite = self.concentration_raw[np.isfinite(self.concentration_raw)]
        if finite.size == 0 or np.any(finite < 0.0):
            raise ValueError("H2SAFE native concentrations must include non-negative values")

    @property
    def event_id(self) -> str:
        return f"H2SAFE-{self.release.lab_id}-{self.release.test_id}"

    @property
    def sample_interval_s(self) -> float | None:
        intervals = np.diff(self.time_s)
        first = float(intervals[0])
        return first if np.allclose(intervals, first, rtol=0.0, atol=1.0e-12) else None

    @property
    def fully_missing_sensor_ids(self) -> tuple[str, ...]:
        return tuple(
            sensor_id for index, sensor_id in enumerate(self.sensor_ids)
            if np.all(~np.isfinite(self.concentration_raw[:, index]))
        )

    @property
    def finite_sensor_ids(self) -> tuple[str, ...]:
        missing = set(self.fully_missing_sensor_ids)
        return tuple(sensor_id for sensor_id in self.sensor_ids if sensor_id not in missing)

    @property
    def native_maximum(self) -> float:
        return float(np.nanmax(self.concentration_raw))

    @property
    def validation_scope(self) -> str:
        return "helium_indoor_transport_component_only"

    @property
    def quantitative_hydrogen_validation_allowed(self) -> bool:
        return False

    @property
    def vol_percent_conversion_allowed(self) -> bool:
        return False


@dataclass(frozen=True)
class H2SafeStlBounds:
    """Axis-aligned bounds read directly from a published ASCII STL."""

    path: Path
    sha256: str
    vertex_record_count: int
    minimum_m: tuple[float, float, float]
    maximum_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if len(self.sha256) != 64 or self.vertex_record_count <= 0:
            raise ValueError("H2SAFE STL bounds provenance is invalid")
        if not all(math.isfinite(value) for value in (*self.minimum_m, *self.maximum_m)):
            raise ValueError("H2SAFE STL bounds must be finite")
        if any(lower > upper for lower, upper in zip(self.minimum_m, self.maximum_m)):
            raise ValueError("H2SAFE STL minimum exceeds maximum")

    def contains(
        self, point_m: tuple[float, float, float], *, tolerance_m: float = 1.0e-5,
    ) -> bool:
        if not math.isfinite(tolerance_m) or tolerance_m < 0.0:
            raise ValueError("H2SAFE STL containment tolerance must be non-negative")
        return all(
            lower - tolerance_m <= value <= upper + tolerance_m
            for value, lower, upper in zip(point_m, self.minimum_m, self.maximum_m)
        )


def _read_registry(path: Path) -> tuple[H2SafeSensor, ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["id", "X", "Y", "Z"]:
            raise ValueError("H2SAFE sensor registry header must be id,X,Y,Z")
        sensors: list[H2SafeSensor] = []
        for row in reader:
            sensor_id = str(row["id"]).strip().lower()
            try:
                position = tuple(float(row[key]) for key in ("X", "Y", "Z"))
            except (TypeError, ValueError) as error:
                raise ValueError(f"H2SAFE sensor {sensor_id!r} has invalid coordinates") from error
            sensors.append(H2SafeSensor(sensor_id, position))  # type: ignore[arg-type]
    if not sensors or len({item.sensor_id for item in sensors}) != len(sensors):
        raise ValueError("H2SAFE sensor registry is empty or contains duplicate IDs")
    return tuple(sensors)


def _native_value(value: str, *, row: int, sensor_id: str) -> float:
    text = value.strip()
    if text in {"", "#N/A", "#NA", "NA", "N/A"}:
        return math.nan
    try:
        result = float(text)
    except ValueError as error:
        raise ValueError(
            f"H2SAFE {sensor_id} row {row} is not numeric or a missing marker"
        ) from error
    if not math.isfinite(result):
        raise ValueError(f"H2SAFE {sensor_id} row {row} must be finite or missing")
    return result


def read_h2safe_run(
    extracted_root: str | Path, *, lab_id: str, test_id: str,
) -> H2SafeRun:
    """Read one unmodified CSV event and its full laboratory sensor registry."""

    release = RELEASES.get((lab_id, test_id))
    if release is None:
        raise ValueError(f"unsupported H2SAFE event {lab_id}/{test_id}")
    root = Path(extracted_root)
    lab = root / lab_id
    registry_path = lab / "sensors.csv"
    source_path = lab / f"{test_id}.csv"
    registry = _read_registry(registry_path)
    with source_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("H2SAFE event CSV is empty") from error
        if len(header) < 2 or header[0].strip().lower() != "time":
            raise ValueError("H2SAFE event CSV first column must be time")
        sensor_ids = tuple(item.strip().lower() for item in header[1:])
        if any(not item for item in sensor_ids):
            raise ValueError("H2SAFE event CSV contains an empty sensor ID")
        times: list[float] = []
        values: list[list[float]] = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"H2SAFE event row {row_number} has {len(row)} fields; "
                    f"expected {len(header)}"
                )
            try:
                time_s = float(row[0])
            except ValueError as error:
                raise ValueError(f"H2SAFE row {row_number} time is not numeric") from error
            times.append(time_s)
            values.append([
                _native_value(value, row=row_number, sensor_id=sensor_id)
                for sensor_id, value in zip(sensor_ids, row[1:], strict=True)
            ])
    return H2SafeRun(
        source_path=source_path.resolve(),
        source_sha256=_sha256(source_path),
        registry_path=registry_path.resolve(),
        registry_sha256=_sha256(registry_path),
        release=release,
        registry=registry,
        sensor_ids=sensor_ids,
        time_s=np.asarray(times, dtype=float),
        concentration_raw=np.asarray(values, dtype=float),
    )


def read_h2safe_ascii_stl_bounds(path: str | Path) -> H2SafeStlBounds:
    """Read native metre-coordinate bounds without adding a mesh dependency."""

    source = Path(path)
    minimum = np.full(3, np.inf)
    maximum = np.full(3, -np.inf)
    count = 0
    try:
        with source.open("r", encoding="ascii") as stream:
            for line_number, line in enumerate(stream, start=1):
                fields = line.split()
                if not fields or fields[0].lower() != "vertex":
                    continue
                if len(fields) != 4:
                    raise ValueError(
                        f"H2SAFE STL vertex line {line_number} is malformed"
                    )
                try:
                    point = np.asarray(fields[1:4], dtype=float)
                except ValueError as error:
                    raise ValueError(
                        f"H2SAFE STL vertex line {line_number} is not numeric"
                    ) from error
                if not np.all(np.isfinite(point)):
                    raise ValueError(
                        f"H2SAFE STL vertex line {line_number} is not finite"
                    )
                minimum = np.minimum(minimum, point)
                maximum = np.maximum(maximum, point)
                count += 1
    except UnicodeDecodeError as error:
        raise ValueError(
            "H2SAFE STL must be the published ASCII representation"
        ) from error
    if count == 0:
        raise ValueError("H2SAFE STL contains no ASCII vertex records")
    return H2SafeStlBounds(
        path=source.resolve(),
        sha256=_sha256(source),
        vertex_record_count=count,
        minimum_m=tuple(float(value) for value in minimum),
        maximum_m=tuple(float(value) for value in maximum),
    )


def read_h2safe_step_product_names(path: str | Path) -> tuple[str, ...]:
    """Return unique assembly product names in their published STEP order."""

    source = Path(path)
    try:
        text = source.read_text(encoding="ascii")
    except UnicodeDecodeError as error:
        raise ValueError("H2SAFE STEP must be the published ASCII file") from error
    names: list[str] = []
    seen: set[str] = set()
    for match in _STEP_PRODUCT.finditer(text):
        name = match.group(1).replace("''", "'").strip()
        if name and name not in seen:
            names.append(name)
            seen.add(name)
    if not names:
        raise ValueError("H2SAFE STEP contains no named PRODUCT entities")
    return tuple(names)


__all__ = [
    "DATASET_DOI", "DATASET_URL", "SPECIES", "CONCENTRATION_UNIT_STATUS",
    "SOURCE_CLOCK_ALIGNMENT_STATUS", "SENSOR_CALIBRATION_STATUS", "RELEASES",
    "PUBLISHED_VERTICAL_AXIS", "h2safe_raw_xyz_to_solver_xyz",
    "LAB1_VENT_BOUNDARIES", "LAB2_SUPPLY_BOUNDARIES",
    "LAB2_RECIRCULATION_FANS", "H2SafeVentBoundary", "H2SafeRecirculationFan",
    "H2SafeReleaseMetadata", "H2SafeSensor", "H2SafeRun", "H2SafeStlBounds",
    "read_h2safe_run", "read_h2safe_ascii_stl_bounds",
    "read_h2safe_step_product_names",
]
