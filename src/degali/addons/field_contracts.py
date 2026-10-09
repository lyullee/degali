"""Portable field-input contracts for the deployable DEGALI model.

The contracts in this module are deliberately independent of a plant historian
or a particular thermodynamic backend.  They make the quantities that drive a
release explicit and attach a bounded uncertainty to them.  A bounded interval
is used when a defensible probability distribution is unavailable; the model
must not silently turn an engineering tolerance into a probability.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from itertools import product
import math
from types import MappingProxyType
from typing import Mapping, Sequence


def _finite(value: float, name: str) -> float:
    # ``bool`` is an ``int`` subclass in Python.  Accepting it here would
    # silently turn a malformed Python boundary such as ``True`` into the
    # physical value ``1.0`` and could alter a source or weather corner.
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite and not boolean")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class BoundedValue:
    """A nominal value with an explicit bounded engineering interval."""

    nominal: float
    lower: float | None = None
    upper: float | None = None
    unit: str = ""
    source: str = "unspecified"

    def __post_init__(self) -> None:
        nominal = _finite(self.nominal, "nominal")
        lower = nominal if self.lower is None else _finite(self.lower, "lower")
        upper = nominal if self.upper is None else _finite(self.upper, "upper")
        if lower > upper:
            raise ValueError("lower must not exceed upper")
        if nominal < lower or nominal > upper:
            raise ValueError("nominal must lie inside [lower, upper]")
        object.__setattr__(self, "nominal", nominal)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def is_exact(self) -> bool:
        return self.lower == self.upper

    def corners(self) -> tuple[float, float]:
        return float(self.lower), float(self.upper)

    def as_dict(self) -> dict[str, object]:
        return {
            "nominal": self.nominal,
            "lower": self.lower,
            "upper": self.upper,
            "unit": self.unit,
            "source": self.source,
        }


@dataclass(frozen=True)
class CircularBoundedValue:
    """A bounded angular quantity on a 0–360 degree circle.

    Unlike :class:`BoundedValue`, a directional interval may cross north.  For
    example, a 350°–10° meteorological record is the short 20° interval, not
    an invalid lower/upper pair or a 340° sector through south.  The object
    retains both declared endpoints for deterministic sensitivity corners; it
    does not assign a probability density or interpolate an angular history.
    """

    nominal: float
    lower: float | None = None
    upper: float | None = None
    unit: str = "deg"
    source: str = "unspecified"

    def __post_init__(self) -> None:
        nominal_raw = _finite(self.nominal, "nominal")
        lower_raw = nominal_raw if self.lower is None else _finite(self.lower, "lower")
        upper_raw = nominal_raw if self.upper is None else _finite(self.upper, "upper")
        raw_span = abs(upper_raw - lower_raw)
        if raw_span >= 360.0 and not math.isclose(raw_span, 0.0, abs_tol=1.0e-12):
            raise ValueError("circular direction interval must span less than 360 degrees")

        nominal = nominal_raw % 360.0
        lower = lower_raw % 360.0
        upper = upper_raw % 360.0
        span = (upper - lower) % 360.0
        offset = (nominal - lower) % 360.0
        if span == 0.0:
            if not math.isclose(offset, 0.0, abs_tol=1.0e-12):
                raise ValueError("nominal direction must lie inside the declared circular interval")
        elif offset > span + 1.0e-12:
            raise ValueError("nominal direction must lie inside the declared circular interval")

        object.__setattr__(self, "nominal", nominal)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def is_exact(self) -> bool:
        return self.lower == self.upper

    @property
    def span_deg(self) -> float:
        """Positive span from ``lower`` to ``upper`` in increasing degrees."""
        return (float(self.upper) - float(self.lower)) % 360.0

    def corners(self) -> tuple[float, float]:
        return float(self.lower), float(self.upper)

    def as_dict(self) -> dict[str, object]:
        return {
            "nominal": self.nominal,
            "lower": self.lower,
            "upper": self.upper,
            "unit": self.unit,
            "source": self.source,
            "circular": True,
            "span_deg": self.span_deg,
            "crosses_zero_deg": not self.is_exact and self.lower > self.upper,
        }


@dataclass(frozen=True)
class FieldCoordinateReference:
    """Declared site-grid relation to true east/north and the vertical datum.

    Field locations are expressed in a plant drawing's local x/y/z grid while
    meteorological bearings are conventionally relative to true east/north.
    ``x_axis_bearing_math_to_deg`` is the mathematical direction of the site
    grid's positive x axis, measured counter-clockwise from true east. It lets
    the wind boundary be rotated into the same grid before source, detector and
    obstacle coordinates enter the local wind-plane calculation.
    """

    coordinate_system_id: str
    origin_id: str
    x_axis_bearing_math_to_deg: float
    vertical_datum_id: str
    evidence_id: str

    def __post_init__(self) -> None:
        for name, value in {
            "coordinate_system_id": self.coordinate_system_id,
            "origin_id": self.origin_id,
            "vertical_datum_id": self.vertical_datum_id,
            "evidence_id": self.evidence_id,
        }.items():
            if not isinstance(value, str) or not value.strip() or value.strip().lower() == "unspecified":
                raise ValueError(f"field coordinate reference requires a declared {name}")
        if not math.isfinite(float(self.x_axis_bearing_math_to_deg)):
            raise ValueError("x_axis_bearing_math_to_deg must be finite")
        object.__setattr__(
            self,
            "x_axis_bearing_math_to_deg",
            float(self.x_axis_bearing_math_to_deg) % 360.0,
        )

    def earth_to_local_math_radians(self, earth_direction_rad: float) -> float:
        """Rotate a true-east/north mathematical bearing into the site grid."""
        if not math.isfinite(float(earth_direction_rad)):
            raise ValueError("earth_direction_rad must be finite")
        return (float(earth_direction_rad) - math.radians(self.x_axis_bearing_math_to_deg)) % (2.0 * math.pi)


@dataclass(frozen=True)
class FieldValidationEvidence:
    """Fingerprint of a matched LH2 concentration-validation dataset.

    The field solver cannot infer validation from a boolean flag.  A usable
    evidence record therefore names the common-clock source, weather,
    obstacle, receptor and temporal operators and fingerprints the external
    artefact.  This is a qualification boundary, not a claim that the
    reduced-order model has been validated merely because the record exists.
    """

    dataset_id: str
    path: str
    sha256: str
    row_count: int
    source_boundary_id: str
    weather_id: str
    obstacle_geometry_id: str
    receptor_geometry_id: str
    temporal_operator_id: str
    common_clock_id: str
    scope: str = "lh2_obstacle_transport"

    def __post_init__(self) -> None:
        for name, value in {
            "dataset_id": self.dataset_id,
            "path": self.path,
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "obstacle_geometry_id": self.obstacle_geometry_id,
            "receptor_geometry_id": self.receptor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
            "common_clock_id": self.common_clock_id,
        }.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip().lower() == "unspecified"
            ):
                raise ValueError(f"validation evidence {name} must be explicitly declared")
        if (
            not isinstance(self.sha256, str)
            or len(self.sha256) != 64
            or any(char not in "0123456789abcdef" for char in self.sha256)
        ):
            raise ValueError("validation evidence sha256 must be a lowercase SHA-256 digest")
        if isinstance(self.row_count, bool) or not isinstance(self.row_count, int) or self.row_count <= 0:
            raise ValueError("validation evidence row_count must be a positive integer")
        if self.scope not in {"lh2_free_field", "lh2_obstacle_transport"}:
            raise ValueError(
                "validation evidence scope must be lh2_free_field or lh2_obstacle_transport"
            )
        if self.scope == "lh2_obstacle_transport" and self.obstacle_geometry_id.strip().lower() in {
            "none", "not_applicable", "n/a",
        }:
            raise ValueError(
                "obstacle-transport validation evidence requires an obstacle geometry identifier"
            )

    @property
    def supports_obstacle_transport(self) -> bool:
        return self.scope == "lh2_obstacle_transport"

    def as_record(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "path": self.path,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "obstacle_geometry_id": self.obstacle_geometry_id,
            "receptor_geometry_id": self.receptor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
            "common_clock_id": self.common_clock_id,
            "scope": self.scope,
        }


def bounded(value: float, *, relative: float = 0.0, unit: str = "",
            source: str = "unspecified") -> BoundedValue:
    """Create a symmetric bounded value without inventing a distribution."""
    value = _finite(value, "value")
    relative = _finite(relative, "relative")
    if relative < 0.0:
        raise ValueError("relative must be non-negative")
    span = abs(value) * relative
    return BoundedValue(value, value - span, value + span, unit, source)


def _validate_optional_unit(value: BoundedValue, expected: str, name: str) -> None:
    """Reject an explicitly wrong unit while retaining legacy blank labels.

    Direct Python callers historically used ``BoundedValue`` without a unit
    label.  The strict case-file boundary already requires exact labels, so
    the direct contract treats an empty label as legacy/unspecified but never
    lets a non-empty, incorrect label pass into a field calculation.
    """
    unit = value.unit
    if not isinstance(unit, str) or (unit and unit != expected):
        raise ValueError(f"{name} must use unit '{expected}' when a unit is declared")


@dataclass(frozen=True)
class PressureDrivenMassFlowBoundary:
    """Typed provenance for a rate derived from the declared throat closure.

    The ambient-pressure interval is a source-boundary input, not an
    independent post-hoc uncertainty on the resulting mass flow.  Keeping it
    typed lets :class:`FieldScenario` recompute the derived rate coherently at
    each deterministic source corner instead of crossing a derived rate
    independently with the pressure/orifice inputs that produced it.
    """

    source_id: str
    ambient_pressure_pa: BoundedValue
    allow_supercritical_gas: bool = False

    def __post_init__(self) -> None:
        if (
            not isinstance(self.source_id, str)
            or not self.source_id.strip()
            or self.source_id.strip().lower() == "unspecified"
        ):
            raise ValueError("pressure-driven mass-flow source_id must be explicitly declared")
        if not isinstance(self.ambient_pressure_pa, BoundedValue):
            raise TypeError("pressure-driven ambient pressure must be a BoundedValue")
        _validate_optional_unit(
            self.ambient_pressure_pa, "Pa", "pressure-driven ambient pressure"
        )
        if (
            not isinstance(self.ambient_pressure_pa.source, str)
            or not self.ambient_pressure_pa.source.strip()
            or self.ambient_pressure_pa.source.strip().lower() == "unspecified"
        ):
            raise ValueError(
                "pressure-driven ambient pressure requires an explicit source"
            )
        if self.ambient_pressure_pa.lower <= 0.0:
            raise ValueError("pressure-driven ambient pressure must be positive")
        if not isinstance(self.allow_supercritical_gas, bool):
            raise TypeError("allow_supercritical_gas must be boolean")
        object.__setattr__(self, "source_id", self.source_id.strip())


@dataclass(frozen=True)
class ReleaseSource:
    """Plant-independent release boundary passed to a dispersion solver.

    ``mass_flow_kg_s`` is the nominal upstream mass-flow boundary used for a
    nominal flash calculation. A transient ``duration_s`` by itself therefore
    means a constant nominal source, not an inferred plant historian. A
    separately declared post-flash atmospheric schedule is required before a
    local transport calculation may use a time-varying source. Pressure
    reference and phase-fraction basis are explicit because neither can be
    inferred safely from a historian tag.
    """

    fluid: str = "hydrogen"
    location_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    direction_m: tuple[float, float, float] = (1.0, 0.0, 0.0)
    upstream_pressure: BoundedValue = field(
        default_factory=lambda: BoundedValue(101325.0, unit="Pa", source="default")
    )
    upstream_temperature: BoundedValue = field(
        default_factory=lambda: BoundedValue(293.15, unit="K", source="default")
    )
    pressure_reference: str = "absolute"
    mass_flow_kg_s: BoundedValue = field(
        default_factory=lambda: BoundedValue(0.0, unit="kg/s", source="default")
    )
    opening_area_m2: BoundedValue = field(
        default_factory=lambda: BoundedValue(1.0e-6, unit="m2", source="default")
    )
    discharge_coefficient: BoundedValue = field(
        default_factory=lambda: BoundedValue(0.8, 0.5, 1.0, "1", "engineering")
    )
    liquid_fraction: BoundedValue = field(
        default_factory=lambda: BoundedValue(1.0, 1.0, 1.0, "1", "declared")
    )
    liquid_fraction_basis: str = "upstream"
    flash_model: str = "declared"
    duration_s: float | None = None
    duration_uncertainty: BoundedValue | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    pressure_driven_mass_flow: PressureDrivenMassFlowBoundary | None = None
    location_uncertainty_m: tuple[BoundedValue, BoundedValue, BoundedValue] | None = None
    direction_uncertainty_m: tuple[BoundedValue, BoundedValue, BoundedValue] | None = None

    def __post_init__(self) -> None:
        for name, value in {
            "upstream_pressure": self.upstream_pressure,
            "upstream_temperature": self.upstream_temperature,
            "mass_flow_kg_s": self.mass_flow_kg_s,
            "opening_area_m2": self.opening_area_m2,
            "discharge_coefficient": self.discharge_coefficient,
            "liquid_fraction": self.liquid_fraction,
        }.items():
            if not isinstance(value, BoundedValue):
                raise TypeError(f"{name} must be a BoundedValue")
        for name, value, expected in (
            ("upstream_pressure", self.upstream_pressure, "Pa"),
            ("upstream_temperature", self.upstream_temperature, "K"),
            ("mass_flow_kg_s", self.mass_flow_kg_s, "kg/s"),
            ("opening_area_m2", self.opening_area_m2, "m2"),
            ("discharge_coefficient", self.discharge_coefficient, "1"),
            ("liquid_fraction", self.liquid_fraction, "1"),
        ):
            _validate_optional_unit(value, expected, name)
        if not isinstance(self.fluid, str) or not self.fluid.strip():
            raise ValueError("fluid must be a non-empty string")
        for name, value in (("location_m", self.location_m), ("direction_m", self.direction_m)):
            if len(value) != 3 or any(isinstance(v, bool) for v in value) or not all(math.isfinite(float(v)) for v in value):
                raise ValueError(f"{name} must contain three finite values")
        if self.location_uncertainty_m is not None:
            if len(self.location_uncertainty_m) != 3 or not all(
                isinstance(value, BoundedValue) for value in self.location_uncertainty_m
            ):
                raise TypeError(
                    "location_uncertainty_m must contain three BoundedValue values"
                )
            for index, (nominal, bound) in enumerate(
                zip(self.location_m, self.location_uncertainty_m)
            ):
                if bound.unit != "m":
                    raise ValueError(
                        f"location_uncertainty_m[{index}] must use unit 'm'"
                    )
                if (
                    not isinstance(bound.source, str)
                    or not bound.source.strip()
                    or bound.source.strip().lower() == "unspecified"
                ):
                    raise ValueError(
                        f"location_uncertainty_m[{index}] requires an explicit source"
                    )
                if not math.isclose(
                    bound.nominal, float(nominal), rel_tol=1.0e-9, abs_tol=1.0e-12
                ):
                    raise ValueError(
                        f"location_uncertainty_m[{index}].nominal must match location_m"
                    )
        if self.direction_uncertainty_m is not None:
            if len(self.direction_uncertainty_m) != 3 or not all(
                isinstance(value, BoundedValue) for value in self.direction_uncertainty_m
            ):
                raise TypeError(
                    "direction_uncertainty_m must contain three BoundedValue values"
                )
            for index, (nominal, bound) in enumerate(
                zip(self.direction_m, self.direction_uncertainty_m)
            ):
                if bound.unit != "1":
                    raise ValueError(
                        f"direction_uncertainty_m[{index}] must use unit '1'"
                    )
                if (
                    not isinstance(bound.source, str)
                    or not bound.source.strip()
                    or bound.source.strip().lower() == "unspecified"
                ):
                    raise ValueError(
                        f"direction_uncertainty_m[{index}] requires an explicit source"
                    )
                if not math.isclose(
                    bound.nominal, float(nominal), rel_tol=1.0e-9, abs_tol=1.0e-12
                ):
                    raise ValueError(
                        f"direction_uncertainty_m[{index}].nominal must match direction_m"
                    )
            for combination in product(
                *(bound.corners() for bound in self.direction_uncertainty_m)
            ):
                if math.sqrt(sum(float(value) ** 2 for value in combination)) <= 0.0:
                    raise ValueError(
                        "direction_uncertainty_m includes a zero-vector corner"
                    )
        norm = math.sqrt(sum(float(v) ** 2 for v in self.direction_m))
        if norm <= 0.0:
            raise ValueError("direction_m must be non-zero")
        if self.duration_s is not None and (
            isinstance(self.duration_s, bool)
            or not math.isfinite(float(self.duration_s))
            or float(self.duration_s) <= 0.0
        ):
            raise ValueError("duration_s must be positive and finite when provided")
        if self.duration_uncertainty is not None:
            if not isinstance(self.duration_uncertainty, BoundedValue):
                raise TypeError("duration_uncertainty must be a BoundedValue or None")
            _validate_optional_unit(self.duration_uncertainty, "s", "duration_uncertainty")
            if (
                not isinstance(self.duration_uncertainty.source, str)
                or not self.duration_uncertainty.source.strip()
                or self.duration_uncertainty.source.strip().lower() == "unspecified"
            ):
                raise ValueError(
                    "duration_uncertainty requires an explicit source"
                )
            if self.duration_s is None:
                raise ValueError("duration_uncertainty requires duration_s")
            if self.duration_uncertainty.lower <= 0.0:
                raise ValueError("duration uncertainty must be positive")
            if not math.isclose(
                self.duration_uncertainty.nominal,
                float(self.duration_s),
                rel_tol=1.0e-9,
                abs_tol=1.0e-12,
            ):
                raise ValueError(
                    "duration_uncertainty nominal must match duration_s"
                )
        if self.upstream_pressure.lower <= 0.0:
            raise ValueError("upstream pressure must be positive")
        if self.upstream_temperature.lower <= 0.0:
            raise ValueError("upstream temperature must be positive")
        if self.mass_flow_kg_s.lower < 0.0:
            raise ValueError("mass flow cannot be negative")
        if self.opening_area_m2.lower <= 0.0:
            raise ValueError("opening area must be positive")
        if self.discharge_coefficient.lower <= 0.0 or self.discharge_coefficient.upper > 1.0:
            raise ValueError("discharge coefficient must lie in (0, 1]")
        if self.liquid_fraction.lower < 0.0 or self.liquid_fraction.upper > 1.0:
            raise ValueError("liquid fraction must lie in [0, 1]")
        if self.pressure_reference not in {"absolute", "gauge", "unspecified"}:
            raise ValueError("pressure_reference must be absolute, gauge, or unspecified")
        if self.liquid_fraction_basis not in {"upstream", "post_flash", "unspecified"}:
            raise ValueError(
                "liquid_fraction_basis must be upstream, post_flash, or unspecified"
            )
        if not isinstance(self.flash_model, str) or not self.flash_model.strip():
            raise ValueError("flash_model must be a non-empty string")
        if not isinstance(self.metadata, Mapping):
            raise TypeError("metadata must be a mapping of strings to strings")
        if any(
            not isinstance(key, str) or not isinstance(value, str)
            for key, value in self.metadata.items()
        ):
            raise ValueError("metadata must map strings to strings")
        if self.pressure_driven_mass_flow is not None and not isinstance(
            self.pressure_driven_mass_flow, PressureDrivenMassFlowBoundary
        ):
            raise TypeError(
                "pressure_driven_mass_flow must be a PressureDrivenMassFlowBoundary or None"
            )
        if self.pressure_driven_mass_flow is not None and (
            self.mass_flow_kg_s.lower == 0.0
            or self.mass_flow_kg_s.upper == 0.0
        ):
            raise ValueError(
                "pressure-driven mass-flow boundary must carry a positive derived rate"
            )
        if self.pressure_driven_mass_flow is not None and (
            not isinstance(self.mass_flow_kg_s.source, str)
            or self.mass_flow_kg_s.source.strip()
            != self.pressure_driven_mass_flow.source_id
        ):
            raise ValueError(
                "pressure-driven mass-flow rate source must match its typed source_id"
            )
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def direction_unit(self) -> tuple[float, float, float]:
        norm = math.sqrt(sum(float(v) ** 2 for v in self.direction_m))
        return tuple(float(v) / norm for v in self.direction_m)

    @property
    def is_transient(self) -> bool:
        return self.duration_s is not None

    def uncertainty_fields(self) -> dict[str, BoundedValue]:
        values = {
            "upstream_pressure": self.upstream_pressure,
            "upstream_temperature": self.upstream_temperature,
            "opening_area_m2": self.opening_area_m2,
            "discharge_coefficient": self.discharge_coefficient,
            "liquid_fraction": self.liquid_fraction,
        }
        if self.pressure_driven_mass_flow is None:
            values["mass_flow_kg_s"] = self.mass_flow_kg_s
        else:
            values[
                "pressure_driven_ambient_pressure_pa"
            ] = self.pressure_driven_mass_flow.ambient_pressure_pa
        if self.duration_uncertainty is not None:
            values["duration_s"] = self.duration_uncertainty
        if self.location_uncertainty_m is not None:
            values.update({
                name: bound
                for name, bound in zip(
                    ("source_location_x_m", "source_location_y_m", "source_location_z_m"),
                    self.location_uncertainty_m,
                )
            })
        if self.direction_uncertainty_m is not None:
            values.update({
                name: bound
                for name, bound in zip(
                    (
                        "source_direction_x", "source_direction_y",
                        "source_direction_z",
                    ),
                    self.direction_uncertainty_m,
                )
            })
        return values

    def as_dict(self) -> dict[str, object]:
        return {
            "fluid": self.fluid,
            "location_m": tuple(self.location_m),
            "direction_m": tuple(self.direction_m),
            "duration_s": self.duration_s,
            "flash_model": self.flash_model,
            "pressure_reference": self.pressure_reference,
            "liquid_fraction_basis": self.liquid_fraction_basis,
            "is_transient": self.is_transient,
            "uncertainty": {k: v.as_dict() for k, v in self.uncertainty_fields().items()},
            "metadata": dict(self.metadata),
            "pressure_driven_mass_flow": (
                None
                if self.pressure_driven_mass_flow is None
                else {
                    "source_id": self.pressure_driven_mass_flow.source_id,
                    "ambient_pressure_pa": self.pressure_driven_mass_flow.ambient_pressure_pa.as_dict(),
                    "allow_supercritical_gas": self.pressure_driven_mass_flow.allow_supercritical_gas,
                }
            ),
        }


@dataclass(frozen=True)
class WeatherState:
    """Meteorological boundary with bounded speed, direction, and stability."""

    speed_m_s: BoundedValue = field(
        default_factory=lambda: BoundedValue(2.0, 2.0, 2.0, "m/s", "default")
    )
    direction_deg: BoundedValue | CircularBoundedValue = field(
        default_factory=lambda: CircularBoundedValue(0.0, 0.0, 0.0, "deg", "default")
    )
    direction_convention: str = "meteorological_from"
    stability: str = "neutral"
    reference_height_m: float = 10.0

    def __post_init__(self) -> None:
        if not isinstance(self.speed_m_s, BoundedValue):
            raise TypeError("speed_m_s must be a BoundedValue")
        _validate_optional_unit(self.speed_m_s, "m/s", "speed_m_s")
        if self.speed_m_s.lower < 0.0:
            raise ValueError("wind speed cannot be negative")
        if not isinstance(self.direction_deg, (BoundedValue, CircularBoundedValue)):
            raise TypeError("direction_deg must be a BoundedValue or CircularBoundedValue")
        if isinstance(self.direction_deg, BoundedValue):
            _validate_optional_unit(self.direction_deg, "deg", "direction_deg")
        elif not isinstance(self.direction_deg.unit, str) or (
            self.direction_deg.unit and self.direction_deg.unit != "deg"
        ):
            raise ValueError("direction_deg must use unit 'deg' when a unit is declared")
        if isinstance(self.reference_height_m, bool) or not math.isfinite(float(self.reference_height_m)) or self.reference_height_m <= 0.0:
            raise ValueError("reference_height_m must be positive and finite")
        if self.stability not in {"very_unstable", "unstable", "neutral", "stable", "very_stable"}:
            raise ValueError("unsupported stability class")
        if self.direction_convention not in {"meteorological_from", "math_to"}:
            raise ValueError("direction_convention must be meteorological_from or math_to")

    def uncertainty_fields(self) -> dict[str, BoundedValue | CircularBoundedValue]:
        return {
            "wind_speed_m_s": self.speed_m_s,
            "wind_direction_deg": self.direction_deg,
        }

    def wind_to_math_radians(self) -> float:
        """Return the nominal wind-*to* direction from global +x, CCW."""
        direction = self.direction_deg.nominal
        if self.direction_convention == "meteorological_from":
            return math.radians((270.0 - direction) % 360.0)
        return math.radians(direction % 360.0)


@dataclass(frozen=True)
class SensorModel:
    """Observation operator after the physical concentration field.

    ``position_m`` uses the same global x/y/z Cartesian frame as
    :attr:`ReleaseSource.location_m`; field workflows rotate it with the
    explicitly declared weather direction before sampling a wind-plane model.
    ``response_time_s`` is the sensor ``t90`` used by the first-order response
    operator.
    """

    position_m: tuple[float, float, float]
    response_time_s: BoundedValue = field(
        default_factory=lambda: BoundedValue(1.0, 1.0, 1.0, "s", "declared")
    )
    gain: BoundedValue = field(
        default_factory=lambda: BoundedValue(1.0, 1.0, 1.0, "1", "declared")
    )
    bias_mole_fraction: BoundedValue = field(
        default_factory=lambda: BoundedValue(0.0, 0.0, 0.0, "mole_fraction", "declared")
    )
    averaging_time_s: float = 0.0

    def __post_init__(self) -> None:
        for name, value in {
            "response_time_s": self.response_time_s,
            "gain": self.gain,
            "bias_mole_fraction": self.bias_mole_fraction,
        }.items():
            if not isinstance(value, BoundedValue):
                raise TypeError(f"{name} must be a BoundedValue")
        for name, value, expected in (
            ("response_time_s", self.response_time_s, "s"),
            ("gain", self.gain, "1"),
            ("bias_mole_fraction", self.bias_mole_fraction, "mole_fraction"),
        ):
            _validate_optional_unit(value, expected, name)
        if len(self.position_m) != 3 or any(isinstance(v, bool) for v in self.position_m) or not all(math.isfinite(float(v)) for v in self.position_m):
            raise ValueError("position_m must contain three finite values")
        if self.position_m[2] < 0.0:
            raise ValueError("sensor height cannot be negative")
        if self.gain.lower <= 0.0:
            raise ValueError("sensor gain must be positive")
        if isinstance(self.averaging_time_s, bool) or not math.isfinite(float(self.averaging_time_s)):
            raise ValueError("sensor averaging time must be finite and not boolean")
        if self.response_time_s.lower <= 0.0 or self.averaging_time_s < 0.0:
            raise ValueError("sensor response and averaging times must be non-negative")

    def uncertainty_fields(self) -> dict[str, BoundedValue | CircularBoundedValue]:
        return {
            "sensor_response_time_s": self.response_time_s,
            "sensor_gain": self.gain,
            "sensor_bias_mole_fraction": self.bias_mole_fraction,
        }


@dataclass(frozen=True)
class SurfaceBoundary:
    """Surface heat-transfer closure used by cryogenic pool/rainout paths."""

    heat_transfer_w_m2_k: BoundedValue = field(
        default_factory=lambda: BoundedValue(10.0, 5.0, 20.0, "W/m2/K", "engineering")
    )
    surface_temperature_k: BoundedValue = field(
        default_factory=lambda: BoundedValue(293.15, 293.15, 293.15, "K", "declared")
    )
    substrate: str = "unspecified"
    evidence_id: str = "unspecified"

    def __post_init__(self) -> None:
        if not isinstance(self.heat_transfer_w_m2_k, BoundedValue):
            raise TypeError("heat_transfer_w_m2_k must be a BoundedValue")
        if not isinstance(self.surface_temperature_k, BoundedValue):
            raise TypeError("surface_temperature_k must be a BoundedValue")
        _validate_optional_unit(
            self.heat_transfer_w_m2_k, "W/m2/K", "heat_transfer_w_m2_k"
        )
        _validate_optional_unit(self.surface_temperature_k, "K", "surface_temperature_k")
        if self.heat_transfer_w_m2_k.lower < 0.0:
            raise ValueError("surface heat transfer cannot be negative")
        if self.surface_temperature_k.lower <= 0.0:
            raise ValueError("surface temperature must be positive")
        if not isinstance(self.substrate, str) or not self.substrate.strip():
            raise ValueError("substrate must be a non-empty string")
        if not isinstance(self.evidence_id, str) or not self.evidence_id.strip():
            raise ValueError("surface evidence_id must be a non-empty string")

    def uncertainty_fields(self) -> dict[str, BoundedValue]:
        return {
            "surface_heat_transfer_w_m2_k": self.heat_transfer_w_m2_k,
            "surface_temperature_k": self.surface_temperature_k,
        }


@dataclass(frozen=True)
class FieldScenario:
    """Complete, auditable input envelope for a field calculation.

    ``model_family`` and ``temporal_mode`` are explicit so a comparison with
    SLABx or a steady/transient replay cannot be hidden in a solver default.
    ``corner_cases`` below is a bounded deterministic sensitivity envelope; it
    intentionally makes no probabilistic claim.
    """

    source: ReleaseSource
    weather: WeatherState = field(default_factory=WeatherState)
    surface: SurfaceBoundary = field(default_factory=SurfaceBoundary)
    sensor: SensorModel | None = None
    model_family: str = "degali"
    temporal_mode: str = "steady"

    def __post_init__(self) -> None:
        if not isinstance(self.source, ReleaseSource):
            raise TypeError("source must be a ReleaseSource")
        if not isinstance(self.weather, WeatherState):
            raise TypeError("weather must be a WeatherState")
        if not isinstance(self.surface, SurfaceBoundary):
            raise TypeError("surface must be a SurfaceBoundary")
        if self.sensor is not None and not isinstance(self.sensor, SensorModel):
            raise TypeError("sensor must be a SensorModel or None")
        if self.model_family not in {"degali", "slabx"}:
            raise ValueError("model_family must be degali or slabx")
        if self.temporal_mode not in {"steady", "transient"}:
            raise ValueError("temporal_mode must be steady or transient")
        if self.temporal_mode == "transient" and not self.source.is_transient:
            raise ValueError("transient scenario requires source.duration_s")

    def uncertainty_fields(self) -> dict[str, BoundedValue]:
        values = dict(self.source.uncertainty_fields())
        values.update(self.weather.uncertainty_fields())
        values.update(self.surface.uncertainty_fields())
        if self.sensor is not None:
            values.update(self.sensor.uncertainty_fields())
        return values

    def corner_cases(self, *, max_cases: int = 1024) -> tuple[dict[str, float], ...]:
        """Return deterministic lower/upper combinations for sensitivity runs."""
        if max_cases < 1:
            raise ValueError("max_cases must be positive")
        fields = self.uncertainty_fields()
        choices = {
            name: ((value.nominal,) if value.is_exact else value.corners())
            for name, value in fields.items()
        }
        count = math.prod(len(values) for values in choices.values())
        if count > max_cases:
            raise ValueError(f"{count} uncertainty corners exceed max_cases={max_cases}")
        names = tuple(choices)
        return tuple(
            {name: float(value) for name, value in zip(names, combination)}
            for combination in product(*(choices[name] for name in names))
        )

    def at_corner(self, values: Mapping[str, float]) -> "FieldScenario":
        """Fix every uncertain input at one fully specified corner.

        The method accepts dictionaries returned by :meth:`corner_cases`.
        Requiring the complete key set prevents a partial sensitivity run from
        silently retaining another uncertain quantity at its nominal value.
        """
        fields = self.uncertainty_fields()
        provided, expected = set(values), set(fields)
        if provided != expected:
            missing = sorted(expected - provided)
            extra = sorted(provided - expected)
            raise ValueError(
                "corner keys must match uncertainty fields; "
                f"missing={missing}, extra={extra}"
            )

        def exact(name: str) -> BoundedValue:
            original = fields[name]
            return BoundedValue(
                values[name], unit=original.unit, source=f"{original.source}; corner"
            )

        source = replace(
            self.source,
            upstream_pressure=exact("upstream_pressure"),
            upstream_temperature=exact("upstream_temperature"),
            opening_area_m2=exact("opening_area_m2"),
            discharge_coefficient=exact("discharge_coefficient"),
            liquid_fraction=exact("liquid_fraction"),
        )
        if self.source.pressure_driven_mass_flow is None:
            source = replace(source, mass_flow_kg_s=exact("mass_flow_kg_s"))
        else:
            # The mass flow is a function of the source/orifice and the typed
            # ambient-pressure boundary.  Recompute it at the same corner
            # rather than crossing the already-derived rate independently.
            from .field_lh2 import release_with_pressure_driven_lh2_mass_flow

            boundary = self.source.pressure_driven_mass_flow
            ambient = exact("pressure_driven_ambient_pressure_pa")
            derivation_metadata = {
                key: value
                for key, value in self.source.metadata.items()
                if not key.startswith("mass_flow_derivation_")
            }
            source = replace(
                source,
                mass_flow_kg_s=BoundedValue(0.0, unit="kg/s", source="pressure-driven-corner"),
                metadata=derivation_metadata,
                pressure_driven_mass_flow=None,
            )
            source = release_with_pressure_driven_lh2_mass_flow(
                source,
                ambient_pressure_pa=ambient.nominal,
                source_id=boundary.source_id,
                allow_supercritical_gas=boundary.allow_supercritical_gas,
            )
        location_names = (
            "source_location_x_m", "source_location_y_m", "source_location_z_m",
        )
        if all(name in fields for name in location_names):
            source = replace(
                source,
                location_m=tuple(float(values[name]) for name in location_names),
                location_uncertainty_m=None,
            )
        direction_names = (
            "source_direction_x", "source_direction_y", "source_direction_z",
        )
        if all(name in fields for name in direction_names):
            source = replace(
                source,
                direction_m=tuple(float(values[name]) for name in direction_names),
                direction_uncertainty_m=None,
            )
        if "duration_s" in fields:
            source = replace(
                source,
                duration_s=float(values["duration_s"]),
                duration_uncertainty=None,
            )
        weather = replace(
            self.weather,
            speed_m_s=exact("wind_speed_m_s"),
            direction_deg=exact("wind_direction_deg"),
        )
        surface = replace(
            self.surface,
            heat_transfer_w_m2_k=exact("surface_heat_transfer_w_m2_k"),
            surface_temperature_k=exact("surface_temperature_k"),
        )
        sensor = self.sensor
        if sensor is not None:
            sensor = replace(
                sensor,
                response_time_s=exact("sensor_response_time_s"),
                gain=exact("sensor_gain"),
                bias_mole_fraction=exact("sensor_bias_mole_fraction"),
            )
        return replace(self, source=source, weather=weather, surface=surface, sensor=sensor)

    def as_dict(self) -> dict[str, object]:
        return {
            "source": self.source.as_dict(),
            "weather": {
                "stability": self.weather.stability,
                "reference_height_m": self.weather.reference_height_m,
                "direction_convention": self.weather.direction_convention,
                "uncertainty": {k: v.as_dict() for k, v in self.weather.uncertainty_fields().items()},
            },
            "surface": {
                "substrate": self.surface.substrate,
                "evidence_id": self.surface.evidence_id,
                "uncertainty": {k: v.as_dict() for k, v in self.surface.uncertainty_fields().items()},
            },
            "model_family": self.model_family,
            "temporal_mode": self.temporal_mode,
            "sensor_present": self.sensor is not None,
        }


@dataclass(frozen=True)
class FieldApplicability:
    """Fail-closed status attached to every deployable calculation."""

    status: str
    reasons: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    uncertainty_complete: bool = False

    def __post_init__(self) -> None:
        if self.status not in {"accepted", "conditional", "blocked"}:
            raise ValueError("status must be accepted, conditional, or blocked")
        if not isinstance(self.uncertainty_complete, bool):
            raise TypeError("uncertainty_complete must be boolean")
        if any(not isinstance(value, str) or not value.strip() for value in self.reasons):
            raise ValueError("applicability reasons must be non-empty strings")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("applicability warnings must be non-empty strings")
        if self.status == "blocked" and not self.reasons:
            raise ValueError("blocked status requires at least one reason")
        if self.status == "blocked" and self.uncertainty_complete:
            raise ValueError("blocked applicability cannot claim complete uncertainty")
        if self.status != "blocked" and self.reasons:
            raise ValueError("non-blocked applicability cannot contain blocking reasons")


def assess_field_applicability(
    scenario: FieldScenario,
    *,
    obstacle_present: bool = False,
    obstacle_supported: bool = True,
    lh2_validation_available: bool = False,
    validation_evidence: FieldValidationEvidence | None = None,
) -> FieldApplicability:
    """Fail-closed gate for a deployable screening calculation.

    Missing validation does not silently become an accepted result.  It is
    reported as ``conditional`` so the caller can still use the result for
    screening while preserving the limitation in an audit trail.  A legacy
    ``lh2_validation_available=True`` flag without a fingerprinted evidence
    record is stronger than missing validation and is therefore blocked.
    """
    reasons: list[str] = []
    warnings: list[str] = []
    if validation_evidence is not None and not isinstance(
        validation_evidence, FieldValidationEvidence
    ):
        raise TypeError("validation_evidence must be a FieldValidationEvidence or None")
    if scenario.source.fluid.strip().lower() not in {"hydrogen", "h2", "lh2"}:
        reasons.append("field contract currently supports hydrogen releases only")
    if scenario.source.pressure_reference != "absolute":
        reasons.append("upstream pressure must be declared as absolute before thermodynamic use")
    if scenario.source.liquid_fraction_basis == "unspecified":
        reasons.append("liquid fraction basis is unresolved")
    if obstacle_present and not obstacle_supported:
        reasons.append("obstacle geometry is outside the supported reduced-order closure")
    if scenario.source.flash_model.strip().lower() in {"", "unknown", "unresolved"}:
        reasons.append("flash state is unresolved")
    if scenario.model_family == "slabx":
        warnings.append("SLABx comparison is a model-form sensitivity, not a replacement calibration")
    if obstacle_present:
        warnings.append("obstacle result uses semi-finite-volume routing and is not obstacle-resolved CFD")
    if validation_evidence is not None:
        warnings.append(
            "fingerprinted LH2 validation evidence is present, but a separate matched "
            "field-validation score is still required before model qualification",
        )
        if obstacle_present and not validation_evidence.supports_obstacle_transport:
            warnings.append(
                "validation evidence is free-field only and does not qualify the declared obstacle transport"
            )
    elif lh2_validation_available:
        reasons.append(
            "lh2_validation_available=True requires an explicit FieldValidationEvidence record"
        )
    else:
        warnings.append("no site-specific LH2 obstacle validation dataset was supplied")
    if reasons:
        return FieldApplicability("blocked", tuple(reasons), tuple(warnings), False)
    status = "conditional" if warnings else "accepted"
    return FieldApplicability(status, (), tuple(warnings), True)


__all__ = [
    "BoundedValue", "CircularBoundedValue", "FieldCoordinateReference", "FieldValidationEvidence", "bounded", "PressureDrivenMassFlowBoundary", "ReleaseSource", "WeatherState", "SensorModel",
    "SurfaceBoundary", "FieldScenario", "FieldApplicability",
    "assess_field_applicability",
]
