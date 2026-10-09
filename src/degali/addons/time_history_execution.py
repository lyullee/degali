"""Manifest-bound source/weather/receptor time-history execution.

This module joins the packet transport operator to the existing hash-pinned
``FieldEvidenceManifest`` without promoting a model.  A valid file package
with missing sensor/operator metadata remains ``conditional``; digest drift,
ID mismatch, or an uncovered clock is ``withheld``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal, Sequence

from .field_evidence_manifest import FieldEvidenceManifest
from .transient_receptor import (
    FixedReceptor,
    PacketSampler,
    ReceptorTrace,
    SourceHistory,
    WindHistory,
    replay_fixed_receptors_with_source_history,
)


TIME_HISTORY_EXECUTION_SCHEMA = "degali.source-weather-receptor-time-history-execution.v1"
TIME_HISTORY_EXECUTION_GATE_CODES = frozenset({
    "manifest_files_verified",
    "manifest_files_withheld",
    "manifest_metadata_complete",
    "manifest_metadata_conditional",
    "event_id_matched",
    "event_id_mismatch",
    "common_clock_matched",
    "common_clock_mismatch",
    "receptor_geometry_matched",
    "receptor_geometry_mismatch",
    "temporal_operator_matched",
    "temporal_operator_mismatch",
    "time_history_valid",
    "time_history_withheld",
    "execution_conditional",
    "execution_withheld",
})


def _non_empty(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


@dataclass(frozen=True)
class TimeHistoryExecutionGate:
    """Decision gate for one manifest-bound time-history run."""

    status: Literal["accepted", "conditional", "withheld"]
    manifest_id: str
    event_id: str
    common_clock_id: str
    source_boundary_id: str
    weather_id: str
    receptor_geometry_id: str
    temporal_operator_id: str
    reasons: tuple[str, ...]
    required_actions: tuple[str, ...]
    gate_codes: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.status not in {"accepted", "conditional", "withheld"}:
            raise ValueError("time-history gate status is invalid")
        for name, value in {
            "manifest_id": self.manifest_id,
            "event_id": self.event_id,
            "common_clock_id": self.common_clock_id,
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "receptor_geometry_id": self.receptor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
        }.items():
            _non_empty(value, name)
        if not isinstance(self.reasons, tuple) or any(
            not isinstance(item, str) or not item.strip() for item in self.reasons
        ):
            raise ValueError("time-history gate reasons must be non-empty strings")
        if not isinstance(self.required_actions, tuple) or any(
            not isinstance(item, str) or not item.strip()
            for item in self.required_actions
        ):
            raise ValueError("time-history gate required_actions must be non-empty strings")
        if not isinstance(self.gate_codes, tuple) or any(
            code not in TIME_HISTORY_EXECUTION_GATE_CODES for code in self.gate_codes
        ) or len(set(self.gate_codes)) != len(self.gate_codes):
            raise ValueError("time-history gate_codes are invalid")
        if self.status == "withheld" and not self.reasons:
            raise ValueError("withheld time-history gate requires a reason")

    @property
    def executable(self) -> bool:
        return self.status != "withheld"

    def as_record(self) -> dict[str, object]:
        return {
            "status": self.status,
            "manifest_id": self.manifest_id,
            "event_id": self.event_id,
            "common_clock_id": self.common_clock_id,
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "receptor_geometry_id": self.receptor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
            "reasons": list(self.reasons),
            "required_actions": list(self.required_actions),
            "gate_codes": list(self.gate_codes),
            "execution_allowed": self.executable,
        }


@dataclass(frozen=True)
class TimeHistoryExecutionResult:
    """Trace summary and gate for one manifest-bound execution."""

    gate: TimeHistoryExecutionGate
    packet_count: int
    released_mass_kg: float | None
    traces: tuple[ReceptorTrace, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.gate, TimeHistoryExecutionGate):
            raise TypeError("gate must be a TimeHistoryExecutionGate")
        if isinstance(self.packet_count, bool) or self.packet_count < 0:
            raise ValueError("packet_count must be a non-negative integer")
        if self.released_mass_kg is not None and (
            not math.isfinite(float(self.released_mass_kg))
            or self.released_mass_kg <= 0.0
        ):
            raise ValueError("released_mass_kg must be positive and finite when supplied")
        if not isinstance(self.traces, tuple) or any(
            not isinstance(item, ReceptorTrace) for item in self.traces
        ):
            raise TypeError("traces must contain ReceptorTrace values")
        if not self.gate.executable and self.traces:
            raise ValueError("withheld time-history execution cannot carry traces")

    def as_record(self) -> dict[str, object]:
        return {
            "schema": TIME_HISTORY_EXECUTION_SCHEMA,
            "gate": self.gate.as_record(),
            "packet_count": self.packet_count,
            "released_mass_kg": self.released_mass_kg,
            "receptors": [
                {
                    "sensor": trace.receptor.name,
                    "sample_count": int(trace.time_s.size),
                    "true_maximum_mole_fraction": float(trace.true_mole_fraction.max()),
                    "indicated_maximum_mole_fraction": float(trace.indicated_mole_fraction.max()),
                }
                for trace in self.traces
            ],
        }


def _mismatch(
    actual: str,
    expected: str,
    label: str,
) -> tuple[str, str]:
    return (
        f"{label} mismatch: manifest={expected!r}, execution={actual!r}",
        f"bind the {label} to the same FieldEvidenceManifest before rerunning",
    )


def run_manifest_bound_time_history(
    manifest: FieldEvidenceManifest,
    *,
    source_event_id: str,
    weather_event_id: str,
    source_common_clock_id: str,
    weather_common_clock_id: str,
    receptor_geometry_id: str,
    temporal_operator_id: str,
    source_history: SourceHistory,
    wind_history: WindHistory,
    receptors: Sequence[FixedReceptor],
    packet_sampler: PacketSampler,
    observation_time_s: Sequence[float] | None = None,
    packet_subdivisions_per_interval: int = 1,
) -> TimeHistoryExecutionResult:
    """Verify a manifest join and run the causal packet operator if permitted."""
    if not isinstance(manifest, FieldEvidenceManifest):
        raise TypeError("manifest must be a FieldEvidenceManifest")
    values = {
        "source_event_id": source_event_id,
        "weather_event_id": weather_event_id,
        "source_common_clock_id": source_common_clock_id,
        "weather_common_clock_id": weather_common_clock_id,
        "receptor_geometry_id": receptor_geometry_id,
        "temporal_operator_id": temporal_operator_id,
    }
    for name, value in values.items():
        values[name] = _non_empty(value, name)
    reasons: list[str] = []
    actions: list[str] = []
    codes: list[str] = []
    try:
        manifest.verify_files()
    except (OSError, ValueError) as error:
        gate = TimeHistoryExecutionGate(
            status="withheld",
            manifest_id=manifest.manifest_id,
            event_id=manifest.event_id,
            common_clock_id=manifest.common_clock_id,
            source_boundary_id=manifest.source_boundary_id,
            weather_id=manifest.weather_id,
            receptor_geometry_id=manifest.receptor_geometry_id,
            temporal_operator_id=manifest.temporal_operator_id,
            reasons=(f"manifest files could not be verified: {error}",),
            required_actions=("restore the hash-pinned manifest files and rerun verification",),
            gate_codes=("manifest_files_withheld", "execution_withheld"),
        )
        return TimeHistoryExecutionResult(gate, 0, None, ())
    codes.append("manifest_files_verified")

    readiness = manifest.readiness_record()
    if readiness["status"] == "conditional":
        codes.extend(("manifest_metadata_conditional", "execution_conditional"))
        reasons.extend(
            f"manifest metadata is conditional: {item}"
            for item in readiness["missing_requirements"]
        )
        actions.append("supply the missing sensor registry/operator metadata before promotion")
    else:
        codes.append("manifest_metadata_complete")

    comparisons = (
        ("event_id", str(source_event_id), manifest.event_id, "event_id_matched", "event_id_mismatch"),
        ("event_id", str(weather_event_id), manifest.event_id, "event_id_matched", "event_id_mismatch"),
        ("common_clock_id", str(source_common_clock_id), manifest.common_clock_id, "common_clock_matched", "common_clock_mismatch"),
        ("common_clock_id", str(weather_common_clock_id), manifest.common_clock_id, "common_clock_matched", "common_clock_mismatch"),
        ("receptor_geometry_id", str(receptor_geometry_id), manifest.receptor_geometry_id, "receptor_geometry_matched", "receptor_geometry_mismatch"),
        ("temporal_operator_id", str(temporal_operator_id), manifest.temporal_operator_id, "temporal_operator_matched", "temporal_operator_mismatch"),
    )
    matched: set[str] = set()
    for label, actual, expected, match_code, mismatch_code in comparisons:
        if actual != expected:
            reason, action = _mismatch(actual, expected, label)
            reasons.append(reason)
            actions.append(action)
            if mismatch_code not in codes:
                codes.append(mismatch_code)
        else:
            matched.add(match_code)
    codes.extend(code for code in sorted(matched) if code not in codes)
    if any(code.endswith("_mismatch") for code in codes):
        codes.append("execution_withheld")
        gate = TimeHistoryExecutionGate(
            status="withheld", manifest_id=manifest.manifest_id,
            event_id=manifest.event_id, common_clock_id=manifest.common_clock_id,
            source_boundary_id=manifest.source_boundary_id, weather_id=manifest.weather_id,
            receptor_geometry_id=manifest.receptor_geometry_id,
            temporal_operator_id=manifest.temporal_operator_id,
            reasons=tuple(dict.fromkeys(reasons)),
            required_actions=tuple(dict.fromkeys(actions)),
            gate_codes=tuple(dict.fromkeys(codes)),
        )
        return TimeHistoryExecutionResult(gate, 0, None, ())

    try:
        packets = source_history.to_packets(
            subdivisions_per_interval=packet_subdivisions_per_interval
        )
        traces = replay_fixed_receptors_with_source_history(
            source_history, wind_history, receptors, packet_sampler,
            observation_time_s=observation_time_s,
            packet_subdivisions_per_interval=packet_subdivisions_per_interval,
        )
    except (TypeError, ValueError, RuntimeError) as error:
        gate = TimeHistoryExecutionGate(
            status="withheld", manifest_id=manifest.manifest_id,
            event_id=manifest.event_id, common_clock_id=manifest.common_clock_id,
            source_boundary_id=manifest.source_boundary_id, weather_id=manifest.weather_id,
            receptor_geometry_id=manifest.receptor_geometry_id,
            temporal_operator_id=manifest.temporal_operator_id,
            reasons=(f"source/weather/receptor time history is not executable: {error}",),
            required_actions=("supply a common covered clock and valid source/weather history",),
            gate_codes=tuple(dict.fromkeys((*codes, "time_history_withheld", "execution_withheld"))),
        )
        return TimeHistoryExecutionResult(gate, 0, None, ())
    codes.append("time_history_valid")
    gate = TimeHistoryExecutionGate(
        status="conditional" if "execution_conditional" in codes else "accepted",
        manifest_id=manifest.manifest_id, event_id=manifest.event_id,
        common_clock_id=manifest.common_clock_id,
        source_boundary_id=manifest.source_boundary_id, weather_id=manifest.weather_id,
        receptor_geometry_id=manifest.receptor_geometry_id,
        temporal_operator_id=manifest.temporal_operator_id,
        reasons=tuple(dict.fromkeys(reasons)),
        required_actions=tuple(dict.fromkeys(actions)),
        gate_codes=tuple(dict.fromkeys(codes)),
    )
    return TimeHistoryExecutionResult(
        gate, len(packets), source_history.released_mass_kg(), traces,
    )


__all__ = [
    "TIME_HISTORY_EXECUTION_SCHEMA",
    "TIME_HISTORY_EXECUTION_GATE_CODES",
    "TimeHistoryExecutionGate",
    "TimeHistoryExecutionResult",
    "run_manifest_bound_time_history",
]
