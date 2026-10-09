"""Explicit source-plane state ledger for model-to-model handoff.

The ledger is deliberately upstream of dispersion.  It records the state a
source calculation claims to hand to an atmospheric model, but it never
turns pressure, liquid fraction, or a static source-state row into an
atmospheric leak history.  A liquid-bearing ledger must be resolved by the
flash/rainout/pool route before :meth:`SourceStateLedger.to_source_history`
can be used.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal, Mapping

from .semi_fv_obstacle import SourceRateSchedule
from .transient_receptor import SourceHistory


SOURCE_STATE_LEDGER_SCHEMA = "degali.source-state-ledger.v1"


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class SourceState:
    """One atmospheric source-plane state on the common event clock."""

    time_s: float
    h2_rate_kg_s: float
    h2_mass_fraction: float
    temperature_k: float
    density_kg_m3: float
    area_m2: float
    velocity_m_s: float = 0.0
    liquid_fraction: float = 0.0
    height_m: float = 0.0
    bearing_to_deg: float = 0.0

    def __post_init__(self) -> None:
        for name in (
            "time_s", "h2_rate_kg_s", "h2_mass_fraction", "temperature_k",
            "density_kg_m3", "area_m2", "velocity_m_s", "liquid_fraction",
            "height_m", "bearing_to_deg",
        ):
            value = _finite(getattr(self, name), name)
            object.__setattr__(self, name, value)
        if self.time_s < 0.0:
            raise ValueError("time_s must be non-negative")
        if self.h2_rate_kg_s < 0.0:
            raise ValueError("h2_rate_kg_s must be non-negative")
        if not 0.0 < self.h2_mass_fraction <= 1.0:
            raise ValueError("h2_mass_fraction must lie in (0, 1]")
        if self.temperature_k <= 0.0 or self.density_kg_m3 <= 0.0 or self.area_m2 <= 0.0:
            raise ValueError("temperature, density, and area must be positive")
        if self.velocity_m_s < 0.0 or self.height_m < 0.0:
            raise ValueError("velocity and height must be non-negative")
        if not 0.0 <= self.liquid_fraction <= 1.0:
            raise ValueError("liquid_fraction must lie in [0, 1]")
        object.__setattr__(self, "bearing_to_deg", self.bearing_to_deg % 360.0)

    @property
    def carrier_rate_kg_s(self) -> float:
        """Total H2+carrier rate implied by the declared H2 fraction."""

        return self.h2_rate_kg_s / self.h2_mass_fraction

    @property
    def momentum_flux_n(self) -> float:
        """One-dimensional carrier momentum flux at the source plane [N]."""

        return self.carrier_rate_kg_s * self.velocity_m_s


@dataclass(frozen=True)
class SourceStateLedger:
    """Hashable-in-content source handoff with an explicit temporal operator."""

    substance: str
    stage: str
    states: tuple[SourceState, ...]
    duration_s: float
    rate_operator: Literal["piecewise_constant", "linear"] = "piecewise_constant"
    observation_operator: str = "not specified"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("substance", "stage", "observation_operator"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
        if not isinstance(self.states, tuple) or len(self.states) < 2:
            raise ValueError("source state ledger requires at least two time states")
        if any(not isinstance(state, SourceState) for state in self.states):
            raise TypeError("states must contain SourceState values")
        duration = _finite(self.duration_s, "duration_s")
        if duration <= 0.0:
            raise ValueError("duration_s must be positive")
        object.__setattr__(self, "duration_s", duration)
        if self.rate_operator not in {"piecewise_constant", "linear"}:
            raise ValueError("rate_operator must be piecewise_constant or linear")
        times = tuple(state.time_s for state in self.states)
        if abs(times[0]) > 1.0e-12:
            raise ValueError("source state ledger must start at event time zero")
        if any(right <= left for left, right in zip(times, times[1:])):
            raise ValueError("source state times must be strictly increasing")
        if not math.isclose(times[-1], duration, rel_tol=1.0e-12, abs_tol=1.0e-12):
            raise ValueError("duration_s must equal the final source-state time")
        if any(state.time_s > duration for state in self.states):
            raise ValueError("source state time exceeds duration_s")
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def has_zero_endpoint(self) -> bool:
        return self.states[-1].h2_rate_kg_s == 0.0

    @property
    def released_h2_mass_kg(self) -> float:
        if self.rate_operator == "linear":
            return math.fsum(
                0.5 * (left.h2_rate_kg_s + right.h2_rate_kg_s) * (right.time_s - left.time_s)
                for left, right in zip(self.states, self.states[1:])
            )
        return math.fsum(
            left.h2_rate_kg_s * (right.time_s - left.time_s)
            for left, right in zip(self.states, self.states[1:])
        )

    @property
    def released_carrier_mass_kg(self) -> float:
        if self.rate_operator == "linear":
            return math.fsum(
                0.5 * (left.carrier_rate_kg_s + right.carrier_rate_kg_s)
                * (right.time_s - left.time_s)
                for left, right in zip(self.states, self.states[1:])
            )
        return math.fsum(
            left.carrier_rate_kg_s * (right.time_s - left.time_s)
            for left, right in zip(self.states, self.states[1:])
        )

    def require_zero_endpoint(self) -> None:
        """Require an explicit shutoff before atmospheric transport."""

        if not self.has_zero_endpoint:
            raise ValueError("source state ledger must end with an explicit zero H2 rate")

    def to_source_history(self, *, allow_liquid: bool = False) -> SourceHistory:
        """Convert a resolved, gas-only ledger into the causal source history."""

        self.require_zero_endpoint()
        if not allow_liquid and any(state.liquid_fraction > 1.0e-12 for state in self.states):
            raise ValueError(
                "liquid-bearing source state requires an explicit flash/rainout/pool handoff"
            )
        return SourceHistory(
            time_s=tuple(state.time_s for state in self.states),
            mass_rate_kg_s=tuple(state.h2_rate_kg_s for state in self.states),
            source_direction_to_deg=tuple(state.bearing_to_deg for state in self.states),
            rate_operator=self.rate_operator,
        )

    def to_source_rate_schedule(
        self,
        *,
        source_id: str | None = None,
        allow_liquid: bool = False,
    ) -> SourceRateSchedule:
        """Expose the resolved H2 rate to conservative obstacle transport.

        This is intentionally a narrow handoff: the obstacle operator receives
        only the already-atmospheric hydrogen schedule.  Liquid/flash state is
        not inferred here and must be resolved by the phase-routing path first.
        """
        self.require_zero_endpoint()
        if not allow_liquid and any(
            state.liquid_fraction > 1.0e-12 for state in self.states
        ):
            raise ValueError(
                "liquid-bearing source state requires an explicit flash/rainout/pool handoff"
            )
        resolved_source_id = source_id or self.metadata.get("source_id", "ledger")
        if not isinstance(resolved_source_id, str) or not resolved_source_id.strip():
            raise ValueError("source_id must be a non-empty string")
        return SourceRateSchedule(
            time_s=tuple(state.time_s for state in self.states),
            rate_kg_s=tuple(state.h2_rate_kg_s for state in self.states),
            source_id=resolved_source_id.strip(),
            rate_operator=self.rate_operator,
        )

    def as_record(self) -> dict[str, Any]:
        return {
            "schema": SOURCE_STATE_LEDGER_SCHEMA,
            "substance": self.substance,
            "stage": self.stage,
            "duration_s": self.duration_s,
            "rate_operator": self.rate_operator,
            "observation_operator": self.observation_operator,
            "metadata": dict(self.metadata),
            "released_h2_mass_kg": self.released_h2_mass_kg,
            "released_carrier_mass_kg": self.released_carrier_mass_kg,
            "states": [asdict(state) for state in self.states],
            "transport_ready": self.has_zero_endpoint and all(
                state.liquid_fraction <= 1.0e-12 for state in self.states
            ),
        }

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "SourceStateLedger":
        if data.get("schema", SOURCE_STATE_LEDGER_SCHEMA) != SOURCE_STATE_LEDGER_SCHEMA:
            raise ValueError("source state ledger schema is not supported")
        raw_states = data.get("states")
        if not isinstance(raw_states, list):
            raise ValueError("source state ledger states must be a JSON array")
        return cls(
            substance=str(data["substance"]),
            stage=str(data["stage"]),
            states=tuple(SourceState(**dict(row)) for row in raw_states),
            duration_s=float(data["duration_s"]),
            rate_operator=str(data.get("rate_operator", "piecewise_constant")),
            observation_operator=str(data.get("observation_operator", "not specified")),
            metadata=dict(data.get("metadata", {})),
        )

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.as_record(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def read_json(cls, path: str | Path) -> "SourceStateLedger":
        return cls.from_mapping(json.loads(Path(path).read_text(encoding="utf-8")))


__all__ = ["SOURCE_STATE_LEDGER_SCHEMA", "SourceState", "SourceStateLedger"]
