"""Explicit stability-to-scalar-mixing closure for local field screening."""

from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .finite_release import SteadyWindApplicability, assess_steady_wind_applicability
from .field_contracts import WeatherState
from .transient_receptor import WindHistory


_STABILITIES = {"very_unstable", "unstable", "neutral", "stable", "very_stable"}


@dataclass(frozen=True)
class FieldStabilityAlternatives:
    """Evidence-backed categorical stability alternatives for a field envelope.

    The nominal stability remains on :class:`WeatherState`; this contract adds
    alternative classes that must be evaluated alongside it.  It deliberately
    carries no probability weights and does not select a turbulence correlation.
    Each class requires a caller-declared diffusivity through
    :class:`StabilityScalarMixingClosure` before it can enter a semi-FV
    envelope.
    """

    alternatives: tuple[str, ...]
    evidence_id: str

    def __post_init__(self) -> None:
        alternatives = tuple(self.alternatives)
        if not alternatives:
            raise ValueError("stability alternatives require at least one declared class")
        if any(not isinstance(item, str) or item not in _STABILITIES for item in alternatives):
            raise ValueError("stability alternatives contain an unsupported class")
        if not isinstance(self.evidence_id, str) or not self.evidence_id.strip() or self.evidence_id == "unspecified":
            raise ValueError("stability alternatives require a declared evidence_id")
        object.__setattr__(self, "alternatives", alternatives)

    def classes_for(self, nominal: str) -> tuple[str, ...]:
        """Return nominal plus alternatives in stable order without duplicates."""
        if nominal not in _STABILITIES:
            raise ValueError(f"unsupported stability class: {nominal!r}")
        return tuple(dict.fromkeys((nominal, *self.alternatives)))


@dataclass(frozen=True)
class StabilityScalarMixingClosure:
    """Caller-declared effective diffusivity by atmospheric stability class.

    This is an auditable scalar-transport closure, not a turbulence-resolving
    RANS/LES model. No stability correlation is selected automatically; the
    caller must provide a positive value for every stability class they use.
    """

    diffusivity_m2_s: Mapping[str, float]
    evidence_id: str
    model_id: str = "declared_stability_scalar_diffusivity"

    def __post_init__(self) -> None:
        if not isinstance(self.evidence_id, str) or not self.evidence_id.strip() or self.evidence_id == "unspecified":
            raise ValueError("stability mixing closure requires a declared evidence_id")
        if not isinstance(self.model_id, str) or not self.model_id.strip():
            raise ValueError("stability mixing closure model_id must be non-empty")
        values = dict(self.diffusivity_m2_s)
        if not values:
            raise ValueError("stability mixing closure requires at least one declared class")
        unknown = set(values) - _STABILITIES
        if unknown:
            raise ValueError(f"unsupported stability classes in mixing closure: {sorted(unknown)}")
        for stability, value in values.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"mixing diffusivity for {stability!r} must be positive and finite")
        object.__setattr__(
            self,
            "diffusivity_m2_s",
            MappingProxyType({name: float(value) for name, value in values.items()}),
        )

    def selected_diffusivity_m2_s(self, stability: str) -> float:
        """Return the explicitly declared value, never a fallback correlation."""
        if stability not in _STABILITIES:
            raise ValueError(f"unsupported stability class: {stability!r}")
        try:
            return self.diffusivity_m2_s[stability]
        except KeyError as error:
            raise ValueError(
                f"mixing closure has no diffusivity for scenario stability {stability!r}"
            ) from error


@dataclass(frozen=True)
class FieldWindHistory:
    """Provenanced measured wind record used to gate a static field screen.

    The local semi-FV operator has one fixed wind plane.  This wrapper does
    not hide direction changes behind a scalar mean: it assesses whether a
    measured record remains inside caller-declared steady-wind limits for the
    release window.  The accepted calculation still uses the explicit
    :class:`~degali.addons.field_contracts.WeatherState` boundary.
    """

    history: WindHistory
    source_id: str
    evidence_id: str
    maximum_direction_span_deg: float = 20.0
    maximum_speed_range_fraction: float = 0.25
    maximum_nominal_speed_relative_difference: float = 0.25
    maximum_nominal_direction_difference_deg: float = 20.0

    def __post_init__(self) -> None:
        if not isinstance(self.history, WindHistory):
            raise TypeError("history must be a WindHistory")
        for name, value in {
            "source_id": self.source_id,
            "evidence_id": self.evidence_id,
        }.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip().lower() == "unspecified"
            ):
                raise ValueError(f"field wind history requires a declared {name}")
        for name, value in {
            "maximum_direction_span_deg": self.maximum_direction_span_deg,
            "maximum_speed_range_fraction": self.maximum_speed_range_fraction,
            "maximum_nominal_speed_relative_difference": (
                self.maximum_nominal_speed_relative_difference
            ),
            "maximum_nominal_direction_difference_deg": (
                self.maximum_nominal_direction_difference_deg
            ),
        }.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        # Validate now and retain immutable values so a mutable caller list
        # cannot alter a scenario after it has been audited.
        time, speed, direction = self.history.arrays()
        object.__setattr__(
            self,
            "history",
            WindHistory(tuple(float(value) for value in time),
                        tuple(float(value) for value in speed),
                        tuple(float(value) for value in direction)),
        )

    def assess(self, duration_s: float) -> SteadyWindApplicability:
        """Apply the declared static-wind gate over the release window."""
        if not math.isfinite(float(duration_s)) or duration_s <= 0.0:
            raise ValueError("wind-history assessment duration must be positive and finite")
        return assess_steady_wind_applicability(
            self.history,
            start_s=0.0,
            end_s=float(duration_s),
            maximum_direction_span_deg=self.maximum_direction_span_deg,
            maximum_speed_range_fraction=self.maximum_speed_range_fraction,
        )

    def nominal_weather_difference(
        self,
        weather: WeatherState,
    ) -> tuple[float, float, float, float]:
        """Compare declared nominal weather to the record's mean wind vector.

        Returns ``(record_speed, record_to_direction_deg, speed_difference,
        direction_difference_deg)``. Direction is vector-averaged, rather
        than arithmetically averaged across the 0/360-degree boundary.
        """
        if not isinstance(weather, WeatherState):
            raise TypeError("weather must be a WeatherState")
        _time, speed, direction_from = self.history.arrays()
        directions_to = np.deg2rad((270.0 - direction_from) % 360.0)
        mean_x = float(np.mean(speed * np.cos(directions_to)))
        mean_y = float(np.mean(speed * np.sin(directions_to)))
        record_speed = math.hypot(mean_x, mean_y)
        if record_speed <= 1.0e-12:
            raise ValueError("wind-history mean vector is indeterminate")
        record_direction = math.degrees(math.atan2(mean_y, mean_x)) % 360.0
        nominal_speed = weather.speed_m_s.nominal
        speed_difference = abs(nominal_speed - record_speed) / record_speed
        nominal_direction = math.degrees(weather.wind_to_math_radians()) % 360.0
        direction_difference = abs(
            (nominal_direction - record_direction + 180.0) % 360.0 - 180.0
        )
        return record_speed, record_direction, speed_difference, direction_difference


__all__ = [
    "FieldStabilityAlternatives", "StabilityScalarMixingClosure", "FieldWindHistory",
]
