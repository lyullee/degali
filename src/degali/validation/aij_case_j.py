"""Strict reader for the independent AIJ Case J non-isothermal dataset.

Case J contains building-wake wind-tunnel statistics with ethylene tracer. It
can qualify component checks for thermal/wake/concentration transport, but it
is neither transient data nor hydrogen-release validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import math
from pathlib import Path


AIJ_CASE_J_DATASET_SCHEMA = "degali.aij-case-j-dataset.v1"
AIJ_CASE_J_RECORD_DOI = "10.5281/zenodo.15447888"

_EXPECTED_MD5 = {
    "readme_caseJ.md": "55bd48aae19fc7c913ba329b57e4a62f",
    "AF_caseJ.csv": "200ac848d1bbdf0db3fc4679d9d8d5f3",
    "BD-ST_caseJ.csv": "a4b0c21abdb9c2ef7e397ce870241c68",
    "RS_caseJ.csv": "7e71dbea4b1b8ce02830a589375fe6a7",
}

_AF_COLUMNS = {
    "x(m)", "y(m)", "z(m)", "U(m/s)", "V(m/s)", "W(m/s)",
    "u_rms(m/s)", "v_rms(m/s)", "w_rms(m/s)", "uv(m2/s2)",
    "uw(m2/s2)", "k(m2/s2)", "T(K)", "t_rms(K)", "ut(m/s K)",
    "vt(m/s K)", "wt(m/s K)",
}
_SURFACE_COLUMNS = {"x(m)", "y(m)", "z(m)", "Building_surface_temperature(K)"}
_RS_COLUMNS = {
    "No.", "x(m)", "y(m)", "z(m)", "U(m/s)", "V(m/s)", "W(m/s)",
    "u_rms(m/s)", "v_rms(m/s)", "w_rms(m/s)", "k(m2/s2)", "T(K)",
    "t_rms(K)", "ut(m/s K)", "vt(m/s K)", "wt(m/s K)", "C(ppm)",
    "c_rms(ppm)",
}


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_rows(path: Path, expected_columns: set[str], label: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or ()) != expected_columns:
            raise ValueError(f"AIJ Case J {label} CSV header is unsupported")
        rows = list(reader)
    if not rows:
        raise ValueError(f"AIJ Case J {label} CSV is empty")
    return rows


def _number(row: dict[str, str], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"AIJ Case J {key} must be numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"AIJ Case J {key} must be finite")
    return value


@dataclass(frozen=True)
class AijCaseJApproachPoint:
    x_m: float
    y_m: float
    z_m: float
    u_m_s: float
    v_m_s: float
    w_m_s: float
    temperature_k: float
    turbulent_kinetic_energy_m2_s2: float


@dataclass(frozen=True)
class AijCaseJSurfacePoint:
    x_m: float
    y_m: float
    z_m: float
    temperature_k: float


@dataclass(frozen=True)
class AijCaseJObservation:
    point_id: int
    x_m: float
    y_m: float
    z_m: float
    u_m_s: float
    v_m_s: float
    w_m_s: float
    temperature_k: float
    turbulent_kinetic_energy_m2_s2: float
    concentration_ppm: float
    concentration_vol_percent: float


@dataclass(frozen=True)
class AijCaseJDataset:
    schema: str
    record_doi: str
    readme_md5: str
    approach_sha256: str
    surface_temperature_sha256: str
    observations_sha256: str
    approach_points: tuple[AijCaseJApproachPoint, ...]
    surface_points: tuple[AijCaseJSurfacePoint, ...]
    observations: tuple[AijCaseJObservation, ...]
    building_width_m: float = 0.08
    building_depth_m: float = 0.08
    building_height_m: float = 0.16
    reference_speed_m_s: float = 1.46
    floor_temperature_k: float = 318.95
    approach_temperature_k: float = 285.32
    reference_reynolds_number: float = 16_000.0
    bulk_richardson_number: float = -0.086
    source_mixture_volume_flow_m3_s: float = 9.17e-6
    source_tracer_mole_fraction: float = 0.05
    source_temperature_k: float = 303.55
    source_hole_diameter_m: float = 0.005
    source_downstream_offset_building_heights: float = 0.25
    tracer: str = "ethylene"
    warnings: tuple[str, ...] = (
        "non-isothermal unstable-stratification wind-tunnel component benchmark",
        "ethylene tracer, not hydrogen",
        "published mean/RMS statistics, not a transient time series",
        "README DOI differs from the Zenodo record DOI",
        "README reports Re≈43,000; the 2026 AIJ benchmark review reports Re≈16,000",
    )

    @property
    def inlet_velocity_profile(self) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Return the measured AF z/U profile with an explicit ground no-slip point."""
        points = sorted(self.approach_points, key=lambda item: item.z_m)
        return (
            (0.0,) + tuple(item.z_m for item in points),
            (0.0,) + tuple(item.u_m_s for item in points),
        )

    @property
    def inlet_temperature_profile(self) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Return the measured AF z/T profile; solver interpolation holds endpoints."""
        points = sorted(self.approach_points, key=lambda item: item.z_m)
        return (
            tuple(item.z_m for item in points),
            tuple(item.temperature_k for item in points),
        )

    @property
    def mean_building_surface_temperature_k(self) -> float:
        """Arithmetic mean of measured surface points, not area-weighted."""
        if not self.surface_points:
            raise ValueError("AIJ Case J surface temperature data are empty")
        return math.fsum(point.temperature_k for point in self.surface_points) / len(self.surface_points)

    @property
    def building_surface_temperature_range_k(self) -> tuple[float, float]:
        """Measured pointwise min/max; the solver currently accepts a uniform wall value."""
        if not self.surface_points:
            raise ValueError("AIJ Case J surface temperature data are empty")
        temperatures = tuple(point.temperature_k for point in self.surface_points)
        return min(temperatures), max(temperatures)


def read_aij_case_j(directory: str | Path, *, verify_official_md5: bool = True) -> AijCaseJDataset:
    """Load the Zenodo CSV package; ppm is converted to vol% as ppm / 10,000."""
    root = Path(directory).resolve()
    paths = {name: root / name for name in _EXPECTED_MD5}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing AIJ Case J files: " + ", ".join(missing))
    digests = {name: _md5(path) for name, path in paths.items()}
    if verify_official_md5:
        mismatches = [
            f"{name} expected {_EXPECTED_MD5[name]}, found {digests[name]}"
            for name in paths if digests[name] != _EXPECTED_MD5[name]
        ]
        if mismatches:
            raise ValueError("AIJ Case J Zenodo MD5 mismatch: " + "; ".join(mismatches))

    approach_rows = _read_rows(paths["AF_caseJ.csv"], _AF_COLUMNS, "approach")
    surface_rows = _read_rows(paths["BD-ST_caseJ.csv"], _SURFACE_COLUMNS, "surface temperature")
    observation_rows = _read_rows(paths["RS_caseJ.csv"], _RS_COLUMNS, "observation")
    if (len(approach_rows), len(surface_rows), len(observation_rows)) != (20, 121, 444):
        raise ValueError(
            "AIJ Case J row counts differ from Zenodo v1: "
            f"approach={len(approach_rows)}, surface={len(surface_rows)}, "
            f"observations={len(observation_rows)}"
        )

    approach = tuple(AijCaseJApproachPoint(
        *(_number(row, key) for key in (
            "x(m)", "y(m)", "z(m)", "U(m/s)", "V(m/s)", "W(m/s)",
            "T(K)", "k(m2/s2)",
        ))
    ) for row in approach_rows)
    surface = tuple(AijCaseJSurfacePoint(
        *(_number(row, key) for key in (
            "x(m)", "y(m)", "z(m)", "Building_surface_temperature(K)",
        ))
    ) for row in surface_rows)
    observations = tuple(AijCaseJObservation(
        int(_number(row, "No.")),
        _number(row, "x(m)"), _number(row, "y(m)"), _number(row, "z(m)"),
        _number(row, "U(m/s)"), _number(row, "V(m/s)"), _number(row, "W(m/s)"),
        _number(row, "T(K)"), _number(row, "k(m2/s2)"),
        _number(row, "C(ppm)"), _number(row, "C(ppm)") / 10_000.0,
    ) for row in observation_rows)
    ids = tuple(item.point_id for item in observations)
    if len(set(ids)) != len(ids):
        raise ValueError("AIJ Case J observation point IDs must be unique")
    if any(item.concentration_ppm < 0.0 for item in observations):
        raise ValueError("AIJ Case J mean concentration must not be negative")

    return AijCaseJDataset(
        schema=AIJ_CASE_J_DATASET_SCHEMA,
        record_doi=AIJ_CASE_J_RECORD_DOI,
        readme_md5=digests["readme_caseJ.md"],
        approach_sha256=hashlib.sha256(paths["AF_caseJ.csv"].read_bytes()).hexdigest(),
        surface_temperature_sha256=hashlib.sha256(paths["BD-ST_caseJ.csv"].read_bytes()).hexdigest(),
        observations_sha256=hashlib.sha256(paths["RS_caseJ.csv"].read_bytes()).hexdigest(),
        approach_points=approach,
        surface_points=surface,
        observations=observations,
    )
