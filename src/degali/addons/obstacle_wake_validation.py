"""Local intake for a quantitative, neutral-gas obstacle-wake benchmark.

The package deliberately does not bundle third-party measurements and this
module does not turn a neutral-tracer wake into an LH2 correction.  It makes
the public AIJ Case-H measurements auditable when a user supplies their local
CSV copy, before any future wake parameterisation is proposed.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import math
from pathlib import Path

from degali.validation.smedis import Sensor, SmedisTrial, read_trial


_REQUIRED_COLUMNS = frozenset({
    "No.", "x (m)", "y (m)", "z (m)", "U (m/s)", "V (m/s)", "W (m/s)",
    "u_rms (m/s)", "v_rms (m/s)", "w_rsm (m/s)", "k (m2/s2)",
    "<c> (ppm)", "c_rms (ppm)",
})
_MISSING_VALUE = -999.9


@dataclass(frozen=True)
class NeutralObstacleWakePoint:
    """One AIJ Case-H wind-tunnel observation in its published coordinates."""

    identifier: int
    x_m: float
    y_m: float
    z_m: float
    velocity_m_s: tuple[float, float, float]
    velocity_rms_m_s: tuple[float, float, float]
    turbulent_kinetic_energy_m2_s2: float
    concentration_ppm: float
    concentration_rms_ppm: float

    def __post_init__(self) -> None:
        """Reject corrupt rows without turning missing sentinels into data.

        AIJ Case-H uses ``-999.9`` for missing values and the reader maps those
        entries to ``nan``.  Missing velocity/concentration values therefore
        remain admissible, while infinities and physically impossible negative
        variances/concentrations are rejected at the typed boundary.
        """
        if isinstance(self.identifier, bool) or not isinstance(self.identifier, int):
            raise TypeError("Case-H point identifier must be an integer")
        if self.identifier <= 0:
            raise ValueError("Case-H point identifier must be positive")
        for name, value in {
            "x_m": self.x_m, "y_m": self.y_m, "z_m": self.z_m,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"Case-H {name} must be finite")
        for name, values in {
            "velocity_m_s": self.velocity_m_s,
            "velocity_rms_m_s": self.velocity_rms_m_s,
        }.items():
            if not isinstance(values, tuple) or len(values) != 3:
                raise TypeError(f"Case-H {name} must contain three values")
            if any(isinstance(value, bool) or math.isinf(float(value)) for value in values):
                raise ValueError(f"Case-H {name} cannot contain infinite values")
        for name, value in {
            "turbulent_kinetic_energy_m2_s2": self.turbulent_kinetic_energy_m2_s2,
            "concentration_ppm": self.concentration_ppm,
            "concentration_rms_ppm": self.concentration_rms_ppm,
        }.items():
            if isinstance(value, bool) or math.isinf(float(value)):
                raise ValueError(f"Case-H {name} cannot be infinite")
        if math.isfinite(float(self.turbulent_kinetic_energy_m2_s2)) and self.turbulent_kinetic_energy_m2_s2 < 0.0:
            raise ValueError("Case-H turbulent kinetic energy cannot be negative")
        # The published ``<c>`` channel is a background-subtracted mean and
        # legitimately contains small signed values near the detection floor.
        # Preserve those observations instead of clipping them or treating
        # them as corrupt.  The RMS channel, in contrast, is a magnitude and
        # must remain non-negative.
        if (
            math.isfinite(float(self.concentration_rms_ppm))
            and self.concentration_rms_ppm < 0.0
        ):
            raise ValueError("Case-H concentration_rms_ppm cannot be negative")
        if any(
            math.isfinite(float(value)) and value < 0.0
            for value in self.velocity_rms_m_s
        ):
            raise ValueError("Case-H velocity RMS values cannot be negative")

    @property
    def concentration_measured(self) -> bool:
        return math.isfinite(self.concentration_ppm)

    @property
    def velocity_measured(self) -> bool:
        return all(math.isfinite(value) for value in self.velocity_m_s)


@dataclass(frozen=True)
class NeutralObstacleWakeBenchmark:
    """A local Case-H intake with a deliberately narrow physical scope."""

    source_path: Path
    points: tuple[NeutralObstacleWakePoint, ...]
    building_width_m: float = 0.10
    building_depth_m: float = 0.10
    building_height_m: float = 0.20
    source_diameter_m: float = 0.004
    source_gas: str = "ethylene"
    neutrally_stratified: bool = True

    @property
    def quantitative_wake_prediction_allowed(self) -> bool:
        """Whether this benchmark alone permits an LH2 wake prediction.

        It intentionally returns ``False``: Case H resolves an excellent
        neutral-gas cuboid wake, but neither cryogenic density evolution nor
        a liquid/condensed source boundary.
        """
        return False

    @property
    def validation_scope(self) -> str:
        return "neutral_gas_single_cuboid_only"


@dataclass(frozen=True)
class DenseGasFenceObservation:
    """A colocated SMEDIS control/fence concentration observation.

    Values retain the published volume-percent convention.  They describe a
    matched dense-gas field trial, not a transferable obstacle multiplier.
    """

    position_m: tuple[float, float, float]
    control_percent: float
    fence_percent: float
    control_std_percent: float
    fence_std_percent: float

    def __post_init__(self) -> None:
        if not isinstance(self.position_m, tuple) or len(self.position_m) != 3:
            raise TypeError("SMEDIS fence position must contain three values")
        if any(
            isinstance(value, bool) or not math.isfinite(float(value))
            for value in self.position_m
        ):
            raise ValueError("SMEDIS fence position must contain finite values")
        for name, value in {
            "control_percent": self.control_percent,
            "fence_percent": self.fence_percent,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"SMEDIS {name} must be finite")
            if value < 0.0:
                raise ValueError(f"SMEDIS {name} cannot be negative")
        for name, value in {
            "control_std_percent": self.control_std_percent,
            "fence_std_percent": self.fence_std_percent,
        }.items():
            if isinstance(value, bool) or math.isinf(float(value)):
                raise ValueError(f"SMEDIS {name} cannot be infinite")
            if math.isfinite(float(value)) and value < 0.0:
                raise ValueError(f"SMEDIS {name} cannot be negative")

    @property
    def fence_to_control_ratio(self) -> float:
        """Return the observed ratio only where a positive control exists."""
        if self.control_percent <= 0.0:
            return math.nan
        return self.fence_percent / self.control_percent


@dataclass(frozen=True)
class DenseGasFenceBenchmark:
    """A no-fit, paired SMEDIS dense-gas fence observation benchmark."""

    control_path: Path
    fence_path: Path
    observations: tuple[DenseGasFenceObservation, ...]
    substance: str
    control_fence_count: int
    fence_count: int

    @property
    def quantitative_lh2_prediction_allowed(self) -> bool:
        """Whether this propane field pair licenses an LH2 wake closure.

        It intentionally returns ``False``.  The pair can falsify a proposed
        scalar-fence model, but it does not provide an LH2 cryogenic source or
        an independently validated density/heat coupling.
        """
        return False

    @property
    def validation_scope(self) -> str:
        return "paired_dense_gas_fence_observation_only"

    @property
    def positive_control_observations(self) -> tuple[DenseGasFenceObservation, ...]:
        return tuple(item for item in self.observations if item.control_percent > 0.0)


def _number(row: dict[str, str], column: str) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid Case-H value in {column!r}") from error
    return math.nan if value == _MISSING_VALUE else value


def read_aij_case_h(path: str | Path) -> NeutralObstacleWakeBenchmark:
    """Read a user-supplied AIJ Case-H ``RS_caseH.csv`` file.

    The file is publicly available under CC BY 4.0 from its Zenodo record,
    but is intentionally not a DEGALI package asset.  The reader rejects an
    altered schema rather than guessing units or a missing velocity component.
    """
    source = Path(path)
    with source.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = frozenset(reader.fieldnames or ())
        missing = _REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(
                "AIJ Case-H CSV is missing required columns: "
                + ", ".join(sorted(missing))
            )
        points = []
        for row in reader:
            identifier_value = _number(row, "No.")
            if not math.isfinite(identifier_value) or not identifier_value.is_integer():
                raise ValueError("AIJ Case-H point identifier must be an integer")
            points.append(NeutralObstacleWakePoint(
                identifier=int(identifier_value),
                x_m=_number(row, "x (m)"),
                y_m=_number(row, "y (m)"),
                z_m=_number(row, "z (m)"),
                velocity_m_s=(
                    _number(row, "U (m/s)"), _number(row, "V (m/s)"),
                    _number(row, "W (m/s)"),
                ),
                velocity_rms_m_s=(
                    _number(row, "u_rms (m/s)"), _number(row, "v_rms (m/s)"),
                    _number(row, "w_rsm (m/s)"),
                ),
                turbulent_kinetic_energy_m2_s2=_number(row, "k (m2/s2)"),
                concentration_ppm=_number(row, "<c> (ppm)"),
                concentration_rms_ppm=_number(row, "c_rms (ppm)"),
            ))
    if not points:
        raise ValueError("AIJ Case-H CSV contains no observations")
    if len({point.identifier for point in points}) != len(points):
        raise ValueError("AIJ Case-H point identifiers must be unique")
    return NeutralObstacleWakeBenchmark(source_path=source, points=tuple(points))


_MATCHED_FENCE_CONDITIONS = {
    # SMEDIS spreadsheets use these documented variants across series.  The
    # reader lower-cases labels but otherwise deliberately leaves them intact.
    "release rate": ("release rate", "release rate (kg/s)"),
    "release duration": ("release duration", "release duration  (s)"),
    "release point x": ("release point x",),
    "release point y": ("release point y", "y"),
    "release point z": ("release point z", "z"),
    "site average windspeed at zref": ("site average windspeed at zref",),
    "ideal wind direction": (
        "ideal wind direction", "ideal wind direction (deg from n)",
    ),
    "surface roughness": ("surface roughness", "surface roughness  (m)"),
}


def _sensor_groups(trial: SmedisTrial) -> dict[tuple[float, float, float], list[Sensor]]:
    """Group repeated colocated channels without discarding their records."""
    groups: dict[tuple[float, float, float], list[Sensor]] = {}
    for sensor in trial.sensors:
        groups.setdefault((sensor.x, sensor.y, sensor.z), []).append(sensor)
    return groups


def _required_condition(
    trial: SmedisTrial, label: str, aliases: tuple[str, ...] = (),
) -> float:
    for key in (label, *aliases):
        value = trial.conditions.get(key, math.nan)
        if math.isfinite(value):
            return value
    raise ValueError(f"SMEDIS trial {trial.path.name} lacks {label!r}")


def pair_smedis_fence_trials(
    control: SmedisTrial, fence: SmedisTrial,
) -> DenseGasFenceBenchmark:
    """Pair a declared no-fence control with a single-fence SMEDIS trial.

    This operation refuses unmatched release, meteorology or sensor layouts.
    Repeated coordinates are retained in spreadsheet order because SMEDIS
    publishes them as separate channels at the same location.
    """
    control_count = _required_condition(control, "number of fences")
    fence_count = _required_condition(fence, "number of fences")
    if control_count != 0.0 or fence_count != 1.0:
        raise ValueError("SMEDIS pair must be a zero-fence control and one-fence trial")
    if control.substance.strip().lower() != fence.substance.strip().lower():
        raise ValueError("SMEDIS fence pair must release the same substance")
    for label, aliases in _MATCHED_FENCE_CONDITIONS.items():
        left = _required_condition(control, label, aliases)
        right = _required_condition(fence, label, aliases)
        if not math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"SMEDIS fence pair differs in {label!r}")

    left_groups = _sensor_groups(control)
    right_groups = _sensor_groups(fence)
    if left_groups.keys() != right_groups.keys():
        raise ValueError("SMEDIS fence pair has different sensor positions")
    observations = []
    for position in sorted(left_groups):
        left_sensors = left_groups[position]
        right_sensors = right_groups[position]
        if len(left_sensors) != len(right_sensors):
            raise ValueError("SMEDIS fence pair has different repeated sensor counts")
        for left, right in zip(left_sensors, right_sensors, strict=True):
            observations.append(DenseGasFenceObservation(
                position_m=position,
                control_percent=left.mean,
                fence_percent=right.mean,
                control_std_percent=left.std,
                fence_std_percent=right.std,
            ))
    if not observations:
        raise ValueError("SMEDIS fence pair contains no concentration observations")
    return DenseGasFenceBenchmark(
        control_path=control.path,
        fence_path=fence.path,
        observations=tuple(observations),
        substance=control.substance,
        control_fence_count=int(control_count),
        fence_count=int(fence_count),
    )


def read_smedis_fence_pair(
    control_path: str | Path, fence_path: str | Path,
) -> DenseGasFenceBenchmark:
    """Read user-held public SMEDIS spreadsheets as one declared pair."""
    return pair_smedis_fence_trials(read_trial(control_path), read_trial(fence_path))


__all__ = [
    "DenseGasFenceBenchmark", "DenseGasFenceObservation",
    "NeutralObstacleWakeBenchmark", "NeutralObstacleWakePoint",
    "pair_smedis_fence_trials", "read_aij_case_h", "read_smedis_fence_pair",
]
