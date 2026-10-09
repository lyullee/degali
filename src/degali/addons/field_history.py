"""Time-aligned measured LH2 release histories for a fail-safe flash boundary.

This module deliberately makes alignment choices visible. It linearly
interpolates declared measurements onto the flow-clock interval starts, but
does not deconvolve a slow instrument or invent data outside a channel's
recorded coverage. Each output is a direct post-flash vapour schedule; liquid
from the flash remains outside the scalar atmospheric source.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
import math
from typing import Literal, Mapping, TYPE_CHECKING

import numpy as np

from .field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
)
from .field_lh2 import (
    lh2_flash_source_from_release,
    release_with_pressure_driven_lh2_mass_flow,
)
from .semi_fv_obstacle import SourceRateSchedule

if TYPE_CHECKING:
    from .lh2_property_table import LH2SaturationTable


_Bound = Literal["nominal", "lower", "upper"]


def _validate_history_selection(
    selection: object,
    name: str,
    *,
    allow_field_values: bool = False,
) -> None:
    """Validate deterministic historian/field corner key-value pairs."""
    if not isinstance(selection, tuple):
        raise TypeError(f"{name} must be a tuple of key/value pairs")
    seen: set[str] = set()
    for item in selection:
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError(f"{name} entries must be 2-item tuples")
        key, value = item
        if not isinstance(key, str) or not key.strip():
            raise ValueError(f"{name} keys must be non-empty strings")
        if key in seen:
            raise ValueError(f"{name} keys must be unique")
        seen.add(key)
        if allow_field_values:
            if isinstance(value, bool):
                raise TypeError(f"{name} values cannot be boolean")
            if isinstance(value, (int, float)):
                if not math.isfinite(float(value)):
                    raise ValueError(f"{name} numeric values must be finite")
            elif not isinstance(value, str) or not value.strip():
                raise TypeError(
                    f"{name} values must be finite numbers or non-empty strings"
                )
        elif value not in {"nominal", "lower", "upper"}:
            raise ValueError(
                f"{name} values must be nominal, lower, or upper"
            )


@dataclass(frozen=True)
class MeasuredTimeSeries:
    """One calibrated measurement channel with per-sample bounded values.

    ``physical_time_s = recorded time_s + time_offset_s``. A positive offset
    therefore moves the recorded trace later on the release clock. The
    response time is audit metadata only: no deconvolution is performed.
    """

    time_s: tuple[float, ...]
    nominal: tuple[float, ...]
    lower: tuple[float, ...] | None = None
    upper: tuple[float, ...] | None = None
    unit: str = ""
    source_id: str = "unspecified"
    time_offset_s: float = 0.0
    response_time_s: float = 0.0

    def __post_init__(self) -> None:
        if len(self.time_s) < 2 or len(self.time_s) != len(self.nominal):
            raise ValueError("measured series needs equal time and nominal sequences of at least two values")
        lower = self.nominal if self.lower is None else self.lower
        upper = self.nominal if self.upper is None else self.upper
        if len(lower) != len(self.time_s) or len(upper) != len(self.time_s):
            raise ValueError("measured series lower and upper sequences must match time length")
        time = tuple(float(value) for value in self.time_s)
        nominal = tuple(float(value) for value in self.nominal)
        lower_values = tuple(float(value) for value in lower)
        upper_values = tuple(float(value) for value in upper)
        if not all(math.isfinite(value) for value in (*time, *nominal, *lower_values, *upper_values)):
            raise ValueError("measured series values must be finite")
        if any(second <= first for first, second in zip(time, time[1:])):
            raise ValueError("measured series times must be strictly increasing")
        if any(not low <= value <= high for low, value, high in zip(lower_values, nominal, upper_values)):
            raise ValueError("measured series nominal must lie inside each [lower, upper] interval")
        if not isinstance(self.source_id, str) or not self.source_id.strip():
            raise ValueError("measured series source_id must be non-empty")
        if not all(math.isfinite(float(value)) for value in (self.time_offset_s, self.response_time_s)):
            raise ValueError("measured series time offset and response time must be finite")
        if self.response_time_s < 0.0:
            raise ValueError("measured series response_time_s must be non-negative")
        object.__setattr__(self, "time_s", time)
        object.__setattr__(self, "nominal", nominal)
        object.__setattr__(self, "lower", lower_values)
        object.__setattr__(self, "upper", upper_values)

    @property
    def physical_time_s(self) -> tuple[float, ...]:
        return tuple(value + self.time_offset_s for value in self.time_s)

    @property
    def has_uncertainty(self) -> bool:
        return any(low != high for low, high in zip(self.lower, self.upper))

    def values_at(self, physical_time_s: tuple[float, ...], *, bound: _Bound = "nominal") -> tuple[float, ...]:
        """Linearly sample one declared bound without extrapolation."""
        if bound not in {"nominal", "lower", "upper"}:
            raise ValueError("bound must be nominal, lower, or upper")
        target = np.asarray(physical_time_s, dtype=float)
        if target.ndim != 1 or target.size == 0 or not np.all(np.isfinite(target)):
            raise ValueError("target physical times must be a non-empty finite sequence")
        axis = np.asarray(self.physical_time_s, dtype=float)
        tolerance = 1.0e-12 * max(1.0, abs(axis[0]), abs(axis[-1]))
        if np.any(target < axis[0] - tolerance) or np.any(target > axis[-1] + tolerance):
            raise ValueError(
                f"measurement {self.source_id!r} does not cover the requested aligned time range"
            )
        source = np.asarray(getattr(self, bound), dtype=float)
        return tuple(float(value) for value in np.interp(target, axis, source))


@dataclass(frozen=True)
class MeasuredReleaseHistory:
    """Pressure, temperature and flow channels needed for time-resolved flash.

    The flow channel supplies the release clock. Other channels are sampled at
    each flow interval start after their own declared time offsets.
    """

    pressure_pa: MeasuredTimeSeries
    temperature_k: MeasuredTimeSeries
    mass_flow_kg_s: MeasuredTimeSeries
    liquid_fraction: MeasuredTimeSeries | None = None

    def __post_init__(self) -> None:
        expected = {
            "pressure_pa": "Pa",
            "temperature_k": "K",
            "mass_flow_kg_s": "kg/s",
        }
        for name, unit in expected.items():
            value = getattr(self, name)
            if not isinstance(value, MeasuredTimeSeries):
                raise TypeError(f"{name} must be a MeasuredTimeSeries")
            if value.unit != unit:
                raise ValueError(f"{name} must declare unit={unit!r}")
        if self.liquid_fraction is not None:
            if not isinstance(self.liquid_fraction, MeasuredTimeSeries):
                raise TypeError("liquid_fraction must be a MeasuredTimeSeries")
            if self.liquid_fraction.unit != "1":
                raise ValueError("liquid_fraction must declare unit='1'")

    @property
    def release_time_s(self) -> tuple[float, ...]:
        axis = self.mass_flow_kg_s.physical_time_s
        origin = axis[0]
        return tuple(value - origin for value in axis)

    @property
    def duration_s(self) -> float:
        return self.release_time_s[-1]

    def channels(self) -> dict[str, MeasuredTimeSeries]:
        result = {
            "pressure_pa": self.pressure_pa,
            "temperature_k": self.temperature_k,
            "mass_flow_kg_s": self.mass_flow_kg_s,
        }
        if self.liquid_fraction is not None:
            result["liquid_fraction"] = self.liquid_fraction
        return result


@dataclass(frozen=True)
class PressureDrivenMeasuredHistory:
    """Explicit pressure/temperature history for a declared orifice source.

    This is intentionally separate from :class:`MeasuredReleaseHistory`:
    pressure and temperature do not become a leak rate unless the caller also
    supplies an explicit absolute-pressure source boundary (opening area,
    discharge coefficient and source identity).  The adapter below recomputes
    the throat rate at each selected pressure/temperature corner and never
    treats a pressure trend as a measured mass-flow channel.
    """

    pressure_pa: MeasuredTimeSeries
    temperature_k: MeasuredTimeSeries
    liquid_fraction: MeasuredTimeSeries | None = None
    event_id: str = "unspecified"
    phase_evidence_id: str = "unspecified"

    def __post_init__(self) -> None:
        for name, value, unit in (
            ("pressure_pa", self.pressure_pa, "Pa"),
            ("temperature_k", self.temperature_k, "K"),
        ):
            if not isinstance(value, MeasuredTimeSeries):
                raise TypeError(f"{name} must be a MeasuredTimeSeries")
            if value.unit != unit:
                raise ValueError(f"{name} must declare unit={unit!r}")
        if self.liquid_fraction is not None:
            if not isinstance(self.liquid_fraction, MeasuredTimeSeries):
                raise TypeError("liquid_fraction must be a MeasuredTimeSeries")
            if self.liquid_fraction.unit != "1":
                raise ValueError("liquid_fraction must declare unit='1'")
        for name, value in {
            "event_id": self.event_id,
            "phase_evidence_id": self.phase_evidence_id,
        }.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip().lower() == "unspecified"
            ):
                raise ValueError(f"pressure-driven history requires a declared {name}")
        if self.pressure_pa.physical_time_s[-1] <= self.pressure_pa.physical_time_s[0]:
            raise ValueError("pressure-driven history must span a positive duration")

    @property
    def release_time_s(self) -> tuple[float, ...]:
        axis = self.pressure_pa.physical_time_s
        origin = axis[0]
        return tuple(value - origin for value in axis)

    @property
    def duration_s(self) -> float:
        return self.release_time_s[-1]

    def channels(self) -> dict[str, MeasuredTimeSeries]:
        result = {
            "pressure_pa": self.pressure_pa,
            "temperature_k": self.temperature_k,
        }
        if self.liquid_fraction is not None:
            result["liquid_fraction"] = self.liquid_fraction
        return result


@dataclass(frozen=True)
class MeasuredHistoryQualityCriteria:
    """Declared promotion limits for a measured transient source boundary.

    These are operational acceptance limits, not a measurement-error model.
    A history may still be flashed for diagnosis without passing them, but it
    cannot be attached to a field transient screening request until it does.
    """

    maximum_sample_interval_s: float
    maximum_response_time_s: float
    maximum_absolute_time_offset_s: float
    maximum_relative_half_width: float
    evidence_id: str

    def __post_init__(self) -> None:
        positive = {"maximum_sample_interval_s": self.maximum_sample_interval_s}
        non_negative = {
            "maximum_response_time_s": self.maximum_response_time_s,
            "maximum_absolute_time_offset_s": self.maximum_absolute_time_offset_s,
            "maximum_relative_half_width": self.maximum_relative_half_width,
        }
        for name, value in positive.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        for name, value in non_negative.items():
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be non-negative and finite")
        if not isinstance(self.evidence_id, str) or not self.evidence_id.strip() or self.evidence_id == "unspecified":
            raise ValueError("measured-history quality criteria require a declared evidence_id")


@dataclass(frozen=True)
class MeasuredHistoryQualityAssessment:
    """Auditable result of checking historian suitability for transient use."""

    criteria: MeasuredHistoryQualityCriteria
    approved: bool
    channel_metrics: tuple[tuple[str, float, float, float, float], ...]
    reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.criteria, MeasuredHistoryQualityCriteria):
            raise TypeError("quality assessment criteria must be MeasuredHistoryQualityCriteria")
        if not isinstance(self.approved, bool):
            raise TypeError("quality assessment approved must be boolean")
        if self.approved != (not self.reasons):
            raise ValueError("quality assessment approval must match its reasons")
        for metric in self.channel_metrics:
            if len(metric) != 5 or not isinstance(metric[0], str) or not metric[0]:
                raise ValueError("quality channel metrics need a channel name and four values")
            if not all(math.isfinite(float(value)) and float(value) >= 0.0 for value in metric[1:]):
                raise ValueError("quality channel metrics must be finite and non-negative")


def assess_measured_history_quality(
    history: MeasuredReleaseHistory,
    criteria: MeasuredHistoryQualityCriteria,
) -> MeasuredHistoryQualityAssessment:
    """Check declared historian timing and interval uncertainty limits.

    The report evaluates each channel's largest sample gap, stated t90-like
    response time, clock offset and interval half-width. A nonzero interval
    around a zero signal is treated as unbounded relative uncertainty and
    fails the declared criterion instead of being divided by an arbitrary
    unit-dependent floor.
    """
    if not isinstance(history, MeasuredReleaseHistory):
        raise TypeError("history must be a MeasuredReleaseHistory")
    if not isinstance(criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("criteria must be MeasuredHistoryQualityCriteria")
    metrics = []
    reasons = []
    for name, channel in history.channels().items():
        maximum_gap = max(np.diff(np.asarray(channel.physical_time_s, dtype=float)))
        maximum_half_width = 0.0
        for nominal, lower, upper in zip(channel.nominal, channel.lower, channel.upper):
            half_width = 0.5 * (upper - lower)
            if nominal == 0.0:
                # Keep the audit metric JSON-safe while guaranteeing failure
                # of any finite declared relative-uncertainty limit.
                relative = (
                    0.0 if half_width == 0.0
                    else max(criteria.maximum_relative_half_width + 1.0, 1.0)
                )
            else:
                relative = half_width / abs(nominal)
            maximum_half_width = max(maximum_half_width, relative)
        metrics.append((
            name, float(maximum_gap), float(channel.response_time_s),
            float(abs(channel.time_offset_s)), float(maximum_half_width),
        ))
        if maximum_gap > criteria.maximum_sample_interval_s:
            reasons.append(f"{name}:sample_interval_exceeds_declared_limit")
        if channel.response_time_s > criteria.maximum_response_time_s:
            reasons.append(f"{name}:response_time_exceeds_declared_limit")
        if abs(channel.time_offset_s) > criteria.maximum_absolute_time_offset_s:
            reasons.append(f"{name}:time_offset_exceeds_declared_limit")
        if maximum_half_width > criteria.maximum_relative_half_width:
            reasons.append(f"{name}:uncertainty_half_width_exceeds_declared_limit")
    return MeasuredHistoryQualityAssessment(
        criteria=criteria,
        approved=not reasons,
        channel_metrics=tuple(metrics),
        reasons=tuple(reasons),
    )


@dataclass(frozen=True)
class FieldMeasuredFlashSchedule:
    """A time-aligned direct-vapour schedule plus explicit phase ledger."""

    schedule: SourceRateSchedule
    history_duration_s: float
    selection: tuple[tuple[str, _Bound], ...]
    total_measured_mass_kg: float
    direct_vapour_mass_kg: float
    unrouted_postflash_liquid_mass_kg: float
    alignment_method: str
    warnings: tuple[str, ...]
    quality_assessment: MeasuredHistoryQualityAssessment | None = None
    provenance: tuple[tuple[str, str], ...] = ()
    source_uncertainty_resolved: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.schedule, SourceRateSchedule):
            raise TypeError("measured-flash schedule must be a SourceRateSchedule")
        self.schedule.require_zero_endpoint("measured-flash schedule")
        _validate_history_selection(self.selection, "measured-flash selection")
        for name, value in {
            "history_duration_s": self.history_duration_s,
            "total_measured_mass_kg": self.total_measured_mass_kg,
            "direct_vapour_mass_kg": self.direct_vapour_mass_kg,
            "unrouted_postflash_liquid_mass_kg": self.unrouted_postflash_liquid_mass_kg,
        }.items():
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
        if not math.isclose(self.schedule.duration_s, self.history_duration_s, rel_tol=1.0e-9, abs_tol=1.0e-12):
            raise ValueError("history duration must match direct-vapour schedule duration")
        if self.direct_vapour_mass_kg > self.total_measured_mass_kg + 1.0e-9:
            raise ValueError("direct vapour mass cannot exceed measured mass")
        if not math.isclose(
            self.schedule.released_mass_kg,
            self.direct_vapour_mass_kg,
            rel_tol=1.0e-8,
            abs_tol=1.0e-10,
        ):
            raise ValueError(
                "direct vapour mass must match the measured-flash schedule mass"
            )
        ledger_tolerance = 1.0e-8 * max(self.total_measured_mass_kg, 1.0)
        if abs(
            self.direct_vapour_mass_kg
            + self.unrouted_postflash_liquid_mass_kg
            - self.total_measured_mass_kg
        ) > ledger_tolerance:
            raise ValueError(
                "direct vapour and unrouted liquid masses must close the measured-mass ledger"
            )
        if self.quality_assessment is not None and not isinstance(
            self.quality_assessment, MeasuredHistoryQualityAssessment
        ):
            raise TypeError("quality_assessment must be MeasuredHistoryQualityAssessment or None")
        if not isinstance(self.source_uncertainty_resolved, bool):
            raise TypeError("source_uncertainty_resolved must be boolean")
        if not isinstance(self.alignment_method, str) or not self.alignment_method.strip():
            raise ValueError("measured-flash alignment_method must be non-empty")
        if not isinstance(self.warnings, tuple):
            raise TypeError("measured-flash warnings must be a tuple")
        if any(not isinstance(item, str) or not item.strip() for item in self.warnings):
            raise ValueError("measured-flash warnings must contain non-empty strings")
        provenance = tuple(self.provenance)
        if any(
            not isinstance(key, str) or not key.strip()
            or not isinstance(value, str) or not value.strip()
            for key, value in provenance
        ):
            raise ValueError("measured-flash provenance must contain non-empty string key/value pairs")
        if len({key for key, _value in provenance}) != len(provenance):
            raise ValueError("measured-flash provenance keys must be unique")
        object.__setattr__(self, "provenance", provenance)


@dataclass(frozen=True)
class FieldMeasuredFlashEnvelopeCase:
    selection: tuple[tuple[str, _Bound], ...]
    schedule: FieldMeasuredFlashSchedule

    def __post_init__(self) -> None:
        _validate_history_selection(self.selection, "measured-flash envelope selection")
        if not isinstance(self.schedule, FieldMeasuredFlashSchedule):
            raise TypeError("measured-flash envelope schedule must be FieldMeasuredFlashSchedule")


@dataclass(frozen=True)
class FieldMeasuredFlashEnvelope:
    """Deterministic global measurement-bound corners, not a probability band."""

    cases: tuple[FieldMeasuredFlashEnvelopeCase, ...]
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.cases, tuple) or not self.cases:
            raise ValueError("measured-flash envelope needs at least one case")
        if any(not isinstance(item, FieldMeasuredFlashEnvelopeCase) for item in self.cases):
            raise TypeError(
                "measured-flash envelope cases must contain only "
                "FieldMeasuredFlashEnvelopeCase values"
            )
        selections = tuple(tuple(case.selection) for case in self.cases)
        if len(set(selections)) != len(selections):
            raise ValueError("measured-flash envelope selections must be unique")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("measured-flash envelope warnings must be non-empty strings")


@dataclass(frozen=True)
class FieldMeasuredHistoryScreeningCase:
    """One measured-source uncertainty selection propagated to field transport."""

    selection: tuple[tuple[str, _Bound | float], ...]
    schedule: FieldMeasuredFlashSchedule
    result: "FieldSemiFVScreeningResult"

    def __post_init__(self) -> None:
        _validate_history_selection(
            self.selection, "measured-history screening selection",
            allow_field_values=True,
        )
        if not isinstance(self.schedule, FieldMeasuredFlashSchedule):
            raise TypeError("measured-history screening schedule must be FieldMeasuredFlashSchedule")
        from .field_workflow import FieldSemiFVScreeningResult
        if not isinstance(self.result, FieldSemiFVScreeningResult):
            raise TypeError(
                "measured-history screening result must be FieldSemiFVScreeningResult"
            )


@dataclass(frozen=True)
class FieldMeasuredHistoryScreeningEnvelope:
    """Deterministic source-history envelope at one fixed field configuration."""

    request: "FieldSemiFVRequest"
    quality_assessment: MeasuredHistoryQualityAssessment
    cases: tuple[FieldMeasuredHistoryScreeningCase, ...]
    warnings: tuple[str, ...]
    property_table_used: bool = False

    def __post_init__(self) -> None:
        from .field_workflow import FieldSemiFVRequest
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("measured-history screening request must be FieldSemiFVRequest")
        if not isinstance(self.quality_assessment, MeasuredHistoryQualityAssessment):
            raise TypeError(
                "measured-history screening quality_assessment must be "
                "MeasuredHistoryQualityAssessment"
            )
        if not isinstance(self.cases, tuple) or not self.cases:
            raise ValueError("measured-history screening envelope needs at least one case")
        if any(not isinstance(item, FieldMeasuredHistoryScreeningCase) for item in self.cases):
            raise TypeError(
                "measured-history screening cases must contain only "
                "FieldMeasuredHistoryScreeningCase values"
            )
        selections = tuple(tuple(case.selection) for case in self.cases)
        if len(set(selections)) != len(selections):
            raise ValueError("measured-history screening selections must be unique")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("measured-history screening warnings must be non-empty strings")
        if not isinstance(self.property_table_used, bool):
            raise TypeError("measured-history screening property_table_used must be boolean")

    @property
    def completed_case_count(self) -> int:
        return sum(case.result.completed for case in self.cases)

    @property
    def blocked_case_count(self) -> int:
        return sum(case.result.applicability.status == "blocked" for case in self.cases)


@dataclass(frozen=True)
class FieldJointMeasuredHistoryScreeningCase:
    """One source-history and independent field-input corner.

    ``source_selection`` only contains channels which have been re-flashed
    from the historian.  ``field_values`` contains ambient boundaries,
    source placement/direction, meteorology and the detector observation
    operator.  Keeping them separate makes it impossible for a report reader
    to mistake a static source assumption for a second, independent
    measurement variation.
    """

    source_selection: tuple[tuple[str, _Bound], ...]
    field_values: tuple[tuple[str, float | str], ...]
    schedule: FieldMeasuredFlashSchedule
    result: "FieldSemiFVScreeningResult"

    def __post_init__(self) -> None:
        _validate_history_selection(self.source_selection, "joint history source selection")
        _validate_history_selection(
            self.field_values, "joint history field selection", allow_field_values=True,
        )
        if not isinstance(self.schedule, FieldMeasuredFlashSchedule):
            raise TypeError("joint history schedule must be FieldMeasuredFlashSchedule")
        from .field_workflow import FieldSemiFVScreeningResult
        if not isinstance(self.result, FieldSemiFVScreeningResult):
            raise TypeError("joint history result must be FieldSemiFVScreeningResult")


@dataclass(frozen=True)
class FieldJointMeasuredHistoryScreeningEnvelope:
    """Deterministic joint source-history, weather, and detector envelope.

    This is a bounded sensitivity calculation, not a confidence or
    probability interval.  Surface parameters are deliberately outside this
    direct-vapour transport envelope: after-flash liquid is explicitly not
    routed through a pool/evaporation source in this workflow.
    """

    request: "FieldSemiFVRequest"
    quality_assessment: MeasuredHistoryQualityAssessment
    cases: tuple[FieldJointMeasuredHistoryScreeningCase, ...]
    warnings: tuple[str, ...]
    property_table_used: bool = False

    def __post_init__(self) -> None:
        from .field_workflow import FieldSemiFVRequest
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("joint history screening request must be FieldSemiFVRequest")
        if not isinstance(self.quality_assessment, MeasuredHistoryQualityAssessment):
            raise TypeError(
                "joint history screening quality_assessment must be "
                "MeasuredHistoryQualityAssessment"
            )
        if not isinstance(self.cases, tuple) or not self.cases:
            raise ValueError("joint history screening envelope needs at least one case")
        if any(not isinstance(item, FieldJointMeasuredHistoryScreeningCase) for item in self.cases):
            raise TypeError(
                "joint history screening cases must contain only "
                "FieldJointMeasuredHistoryScreeningCase values"
            )
        selections = tuple(
            (tuple(case.source_selection), tuple(case.field_values))
            for case in self.cases
        )
        if len(set(selections)) != len(selections):
            raise ValueError("joint history screening selections must be unique")
        if any(not isinstance(value, str) or not value.strip() for value in self.warnings):
            raise ValueError("joint history screening warnings must be non-empty strings")
        if not isinstance(self.property_table_used, bool):
            raise TypeError("joint history screening property_table_used must be boolean")

    @property
    def completed_case_count(self) -> int:
        return sum(case.result.completed for case in self.cases)

    @property
    def blocked_case_count(self) -> int:
        return sum(case.result.applicability.status == "blocked" for case in self.cases)


def _selection(
    history: MeasuredReleaseHistory,
    requested: Mapping[str, _Bound] | None,
) -> dict[str, _Bound]:
    channels = history.channels()
    selected = dict(requested or {})
    unknown = set(selected) - set(channels)
    if unknown:
        raise ValueError(f"unknown measured-history selections: {sorted(unknown)}")
    for name in channels:
        selected.setdefault(name, "nominal")
        if selected[name] not in {"nominal", "lower", "upper"}:
            raise ValueError("measurement selection values must be nominal, lower or upper")
    return selected


def direct_vapour_schedule_from_measured_history(
    release: ReleaseSource,
    history: MeasuredReleaseHistory,
    *,
    selection: Mapping[str, _Bound] | None = None,
    ambient_temperature_k: float = 295.0,
    ambient_pressure_pa: float = 101325.0,
    droplet_size_coefficient: float = 15.0,
    property_table: "LH2SaturationTable | None" = None,
    source_id: str = "measured-history:postflash-direct-vapour",
    quality_criteria: MeasuredHistoryQualityCriteria | None = None,
) -> FieldMeasuredFlashSchedule:
    """Flash aligned measured source states into a direct-vapour schedule.

    The release duration and the aligned flow duration must match exactly;
    this prevents an incomplete historian window from being silently treated
    as a complete transient release. Static opening, discharge coefficient
    and any undeclared liquid fraction remain the caller's explicit release
    assumptions. Measurement uncertainty is selected by named deterministic
    bounds, never converted to a probability distribution.
    """
    if not isinstance(release, ReleaseSource):
        raise TypeError("release must be a ReleaseSource")
    if not isinstance(history, MeasuredReleaseHistory):
        raise TypeError("history must be a MeasuredReleaseHistory")
    if release.duration_s is None:
        raise ValueError("measured source history requires release.duration_s")
    if (
        release.duration_uncertainty is not None
        and not release.duration_uncertainty.is_exact
    ):
        raise ValueError(
            "measured source history fixes the time axis; unresolved "
            "duration_uncertainty requires a time-aligned source envelope"
        )
    if not math.isclose(release.duration_s, history.duration_s, rel_tol=1.0e-9, abs_tol=1.0e-12):
        raise ValueError("release.duration_s must match the aligned flow-history duration")
    if not isinstance(source_id, str) or not source_id.strip():
        raise ValueError("source_id must be non-empty")
    if quality_criteria is not None and not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria or None")
    quality_assessment = (
        None if quality_criteria is None
        else assess_measured_history_quality(history, quality_criteria)
    )
    selected = _selection(history, selection)
    flow_clock = history.mass_flow_kg_s.physical_time_s
    sampled = {
        name: channel.values_at(flow_clock, bound=selected[name])
        for name, channel in history.channels().items()
    }
    if "liquid_fraction" not in sampled:
        sampled["liquid_fraction"] = (release.liquid_fraction.nominal,) * len(flow_clock)
    flash_rates: list[float] = []
    liquid_rates: list[float] = []
    for index in range(len(flow_clock) - 1):
        measured_flow = sampled["mass_flow_kg_s"][index]
        if measured_flow < 0.0:
            raise ValueError("aligned measured mass flow cannot be negative")
        if measured_flow == 0.0:
            flash_rates.append(0.0)
            liquid_rates.append(0.0)
            continue
        state = replace(
            release,
            upstream_pressure=BoundedValue(
                sampled["pressure_pa"][index], unit="Pa",
                source=f"{history.pressure_pa.source_id}; aligned-{selected['pressure_pa']}",
            ),
            upstream_temperature=BoundedValue(
                sampled["temperature_k"][index], unit="K",
                source=f"{history.temperature_k.source_id}; aligned-{selected['temperature_k']}",
            ),
            mass_flow_kg_s=BoundedValue(
                measured_flow, unit="kg/s",
                source=f"{history.mass_flow_kg_s.source_id}; aligned-{selected['mass_flow_kg_s']}",
            ),
            liquid_fraction=BoundedValue(
                sampled["liquid_fraction"][index], unit="1",
                source=(
                    f"{history.liquid_fraction.source_id}; aligned-{selected['liquid_fraction']}"
                    if history.liquid_fraction is not None else release.liquid_fraction.source
                ),
            ),
        )
        flash_result = lh2_flash_source_from_release(
            state,
            ambient_temperature_k=ambient_temperature_k,
            ambient_pressure_pa=ambient_pressure_pa,
            droplet_size_coefficient=droplet_size_coefficient,
            property_table=property_table,
        )
        if not flash_result.conservative:
            raise ValueError(
                "measured-history flash interval is not conservatively closed: "
                + "; ".join(flash_result.closure_warnings)
            )
        flash = flash_result.flash
        flash_rates.append(float(flash.vapour_mass_flow))
        liquid_rates.append(float(flash.liquid_mass_flow))
    relative_times = history.release_time_s
    schedule = SourceRateSchedule(
        relative_times, tuple(flash_rates) + (0.0,), source_id=source_id,
    )
    intervals = np.diff(np.asarray(relative_times, dtype=float))
    measured_mass = float(np.dot(np.asarray(sampled["mass_flow_kg_s"][:-1]), intervals))
    liquid_mass = float(np.dot(np.asarray(liquid_rates), intervals))
    warnings = [
        "pressure, temperature and flow are linearly interpolated onto declared flow interval starts; no extrapolation was used",
        "direct post-flash vapour only is supplied to atmospheric scalar transport; post-flash liquid remains unrouted",
    ]
    delayed = [
        name for name, channel in history.channels().items()
        if channel.response_time_s > 0.0
    ]
    if delayed:
        warnings.append(
            "measurement response times are recorded but not deconvolved for: " + ", ".join(delayed)
        )
    if history.liquid_fraction is None:
        warnings.append("liquid fraction is not time-resolved; ReleaseSource nominal fraction was used")
    if quality_assessment is None:
        warnings.append(
            "measured history has not passed declared operational quality criteria; it cannot be promoted to a field transient source"
        )
    elif quality_assessment.approved:
        warnings.append(
            "measured history passed declared operational quality criteria "
            f"with evidence_id={quality_assessment.criteria.evidence_id!r}"
        )
    else:
        warnings.append(
            "measured history failed declared operational quality criteria: "
            + ", ".join(quality_assessment.reasons)
        )
    return FieldMeasuredFlashSchedule(
        schedule=schedule,
        history_duration_s=history.duration_s,
        selection=tuple(sorted(selected.items())),
        total_measured_mass_kg=measured_mass,
        direct_vapour_mass_kg=schedule.released_mass_kg,
        unrouted_postflash_liquid_mass_kg=liquid_mass,
        alignment_method="linear-to-flow-interval-start",
        warnings=tuple(warnings),
        quality_assessment=quality_assessment,
        source_uncertainty_resolved=not any(
            channel.has_uncertainty for channel in history.channels().values()
        ),
    )


def direct_vapour_schedule_envelope_from_measured_history(
    release: ReleaseSource,
    history: MeasuredReleaseHistory,
    *,
    max_cases: int = 64,
    **kwargs: object,
) -> FieldMeasuredFlashEnvelope:
    """Evaluate every global lower/upper measurement-bound corner.

    A corner applies one selected bound to every sample in that channel. These
    are transparent sensitivity cases, not independent per-sample random
    draws and not confidence intervals.
    """
    if not isinstance(max_cases, int) or max_cases < 1:
        raise ValueError("max_cases must be a positive integer")
    channels = history.channels()
    options = {
        name: (("lower", "upper") if channel.has_uncertainty else ("nominal",))
        for name, channel in channels.items()
    }
    names = tuple(options)
    count = math.prod(len(options[name]) for name in names)
    if count > max_cases:
        raise ValueError(f"{count} measured-history corners exceed max_cases={max_cases}")
    cases = []
    for values in product(*(options[name] for name in names)):
        selected = dict(zip(names, values))
        schedule = direct_vapour_schedule_from_measured_history(
            release, history, selection=selected, **kwargs,
        )
        cases.append(FieldMeasuredFlashEnvelopeCase(tuple(sorted(selected.items())), schedule))
    return FieldMeasuredFlashEnvelope(
        tuple(cases),
        warnings=(
            "history envelope uses global lower/upper channel selections; it is deterministic sensitivity, not a probability interval",
        ),
    )


def _pressure_driven_selection(
    history: PressureDrivenMeasuredHistory,
    requested: Mapping[str, _Bound] | None,
) -> dict[str, _Bound]:
    channels = history.channels()
    selected = dict(requested or {})
    unknown = set(selected) - set(channels)
    if unknown:
        raise ValueError(f"unknown pressure-driven history selections: {sorted(unknown)}")
    for name in channels:
        selected.setdefault(name, "nominal")
        if selected[name] not in {"nominal", "lower", "upper"}:
            raise ValueError(
                "pressure-driven history selection values must be nominal, lower or upper"
            )
    return selected


def direct_vapour_schedule_from_pressure_driven_history(
    release: ReleaseSource,
    history: PressureDrivenMeasuredHistory,
    *,
    selection: Mapping[str, _Bound] | None = None,
    ambient_temperature_k: float = 295.0,
    ambient_pressure_pa: float = 101325.0,
    droplet_size_coefficient: float = 15.0,
    property_table: "LH2SaturationTable | None" = None,
    source_id: str = "pressure-driven-history:postflash-direct-vapour",
    quality_criteria: MeasuredHistoryQualityCriteria | None = None,
    allow_supercritical_gas: bool = False,
) -> FieldMeasuredFlashSchedule:
    """Derive a coherent time-varying source from an explicit orifice history.

    Pressure and temperature are only boundary states.  Each interval is
    recomputed through the declared absolute-pressure throat closure with the
    caller's fixed opening area, discharge coefficient, ambient pressure and
    source identity.  A pressure trend alone can therefore never become a
    leak rate.  The envelope helper below re-flashes lower/nominal/upper
    pressure, temperature and liquid-fraction choices together.
    """
    if not isinstance(release, ReleaseSource):
        raise TypeError("release must be a ReleaseSource")
    if not isinstance(history, PressureDrivenMeasuredHistory):
        raise TypeError("history must be a PressureDrivenMeasuredHistory")
    if release.duration_s is None:
        raise ValueError("pressure-driven history requires release.duration_s")
    if release.duration_uncertainty is not None and not release.duration_uncertainty.is_exact:
        raise ValueError(
            "pressure-driven history fixes the time axis; unresolved duration_uncertainty requires a separate time-aligned source envelope"
        )
    if not math.isclose(
        release.duration_s, history.duration_s, rel_tol=1.0e-9, abs_tol=1.0e-12,
    ):
        raise ValueError("release.duration_s must match the pressure-history duration")
    if release.pressure_reference != "absolute":
        raise ValueError("pressure-driven history requires absolute upstream pressure")
    if release.mass_flow_kg_s.lower != 0.0 or release.mass_flow_kg_s.upper != 0.0:
        raise ValueError(
            "pressure-driven history will not overwrite a non-zero measured mass-flow boundary"
        )
    if release.pressure_driven_mass_flow is not None:
        raise ValueError("pressure-driven history requires an unattached source")
    if not release.opening_area_m2.is_exact or not release.discharge_coefficient.is_exact:
        raise ValueError(
            "pressure-driven history requires exact opening_area_m2 and discharge_coefficient; evaluate their uncertainty as a separate source envelope"
        )
    if history.liquid_fraction is None and not release.liquid_fraction.is_exact:
        raise ValueError(
            "pressure-driven history requires a time-resolved liquid fraction or an exact ReleaseSource liquid_fraction"
        )
    if not isinstance(allow_supercritical_gas, bool):
        raise TypeError("allow_supercritical_gas must be boolean")
    if not isinstance(source_id, str) or not source_id.strip() or source_id.strip().lower() == "unspecified":
        raise ValueError("pressure-driven history source_id must be explicitly declared")
    for name, value in {
        "ambient_temperature_k": ambient_temperature_k,
        "ambient_pressure_pa": ambient_pressure_pa,
        "droplet_size_coefficient": droplet_size_coefficient,
    }.items():
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if quality_criteria is not None and not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria or None")

    selected = _pressure_driven_selection(history, selection)
    clock = history.pressure_pa.physical_time_s
    pressure_values = history.pressure_pa.values_at(clock, bound=selected["pressure_pa"])
    temperature_values = history.temperature_k.values_at(clock, bound=selected["temperature_k"])
    liquid_values = (
        history.liquid_fraction.values_at(clock, bound=selected["liquid_fraction"])
        if history.liquid_fraction is not None
        else (release.liquid_fraction.nominal,) * len(clock)
    )
    flash_rates: list[float] = []
    liquid_rates: list[float] = []
    derived_rates: list[float] = []
    for index in range(len(clock) - 1):
        state = replace(
            release,
            upstream_pressure=BoundedValue(
                pressure_values[index], unit="Pa",
                source=f"{history.pressure_pa.source_id}; aligned-{selected['pressure_pa']}",
            ),
            upstream_temperature=BoundedValue(
                temperature_values[index], unit="K",
                source=f"{history.temperature_k.source_id}; aligned-{selected['temperature_k']}",
            ),
            mass_flow_kg_s=BoundedValue(
                0.0, unit="kg/s", source="pressure-driven-history-placeholder",
            ),
            liquid_fraction=BoundedValue(
                liquid_values[index], unit="1",
                source=(
                    f"{history.liquid_fraction.source_id}; aligned-{selected['liquid_fraction']}"
                    if history.liquid_fraction is not None else release.liquid_fraction.source
                ),
            ),
            pressure_driven_mass_flow=None,
            duration_s=None,
            duration_uncertainty=None,
        )
        derived = release_with_pressure_driven_lh2_mass_flow(
            state, ambient_pressure_pa=ambient_pressure_pa,
            source_id=source_id, allow_supercritical_gas=allow_supercritical_gas,
        )
        derived_rates.append(float(derived.mass_flow_kg_s.nominal))
        flash_result = lh2_flash_source_from_release(
            derived, ambient_temperature_k=ambient_temperature_k,
            ambient_pressure_pa=ambient_pressure_pa,
            droplet_size_coefficient=droplet_size_coefficient,
            property_table=property_table,
        )
        if not flash_result.conservative:
            raise ValueError(
                "pressure-driven history flash interval is not conservatively closed: "
                + "; ".join(flash_result.closure_warnings)
            )
        flash_rates.append(float(flash_result.flash.vapour_mass_flow))
        liquid_rates.append(float(flash_result.flash.liquid_mass_flow))

    relative_times = history.release_time_s
    schedule = SourceRateSchedule(
        relative_times, tuple(flash_rates) + (0.0,), source_id=source_id,
    )
    intervals = np.diff(np.asarray(clock, dtype=float))
    total_mass = float(np.dot(np.asarray(derived_rates), intervals))
    liquid_mass = float(np.dot(np.asarray(liquid_rates), intervals))
    derived_flow = MeasuredTimeSeries(
        time_s=history.pressure_pa.time_s,
        nominal=tuple(derived_rates) + (0.0,),
        unit="kg/s", source_id=f"{source_id}:derived-mass-flow",
        time_offset_s=history.pressure_pa.time_offset_s,
        response_time_s=max(history.pressure_pa.response_time_s, history.temperature_k.response_time_s),
    )
    quality_history = MeasuredReleaseHistory(
        pressure_pa=history.pressure_pa, temperature_k=history.temperature_k,
        mass_flow_kg_s=derived_flow, liquid_fraction=history.liquid_fraction,
    )
    quality_assessment = (
        None if quality_criteria is None
        else assess_measured_history_quality(quality_history, quality_criteria)
    )
    warnings = [
        "mass flow was recomputed interval-by-interval from explicit absolute pressure, opening area, discharge coefficient, ambient pressure and source identity; pressure trend alone was not treated as a leak rate",
        "direct post-flash vapour only is supplied to atmospheric scalar transport; post-flash liquid remains unrouted",
    ]
    if quality_assessment is None:
        warnings.append(
            "pressure-driven history has no declared operational quality criteria; it remains diagnostic and cannot be promoted to a field transient source"
        )
    elif quality_assessment.approved:
        warnings.append(
            "pressure-driven history passed declared operational quality criteria "
            f"with evidence_id={quality_assessment.criteria.evidence_id!r}"
        )
    else:
        warnings.append(
            "pressure-driven history failed declared operational quality criteria: "
            + ", ".join(quality_assessment.reasons)
        )
    provenance = (
        ("history_kind", "pressure_driven_orifice"),
        ("event_id", history.event_id),
        ("phase_evidence_id", history.phase_evidence_id),
        ("source_boundary_id", source_id.strip()),
        ("ambient_pressure_pa", repr(float(ambient_pressure_pa))),
        ("selection", repr(tuple(sorted(selected.items())))),
    )
    return FieldMeasuredFlashSchedule(
        schedule=schedule, history_duration_s=history.duration_s,
        selection=tuple(sorted(selected.items())),
        total_measured_mass_kg=total_mass,
        direct_vapour_mass_kg=schedule.released_mass_kg,
        unrouted_postflash_liquid_mass_kg=liquid_mass,
        alignment_method="pressure-clock-to-interval-start",
        warnings=tuple(warnings), quality_assessment=quality_assessment,
        provenance=provenance,
        source_uncertainty_resolved=not any(
            channel.has_uncertainty for channel in history.channels().values()
        ),
    )


def direct_vapour_schedule_envelope_from_pressure_driven_history(
    release: ReleaseSource,
    history: PressureDrivenMeasuredHistory,
    *,
    max_cases: int = 64,
    **kwargs: object,
) -> FieldMeasuredFlashEnvelope:
    """Recompute every coherent historian and orifice-source bound corner."""
    if not isinstance(max_cases, int) or max_cases < 1:
        raise ValueError("max_cases must be a positive integer")
    channels = dict(history.channels())
    # Unlike a measured-flow historian, this route still predicts the throat
    # rate. Opening area and discharge coefficient are therefore source
    # boundary inputs rather than post-hoc multipliers and must be crossed
    # coherently with each historian state selection.
    source_fields = {
        "opening_area_m2": release.opening_area_m2,
        "discharge_coefficient": release.discharge_coefficient,
    }
    if any(not isinstance(value, BoundedValue) for value in source_fields.values()):
        raise TypeError(
            "pressure-driven source area and discharge coefficient must be BoundedValue values"
        )
    options = {
        name: (("lower", "upper") if channel.has_uncertainty else ("nominal",))
        for name, channel in channels.items()
    }
    options.update({
        name: (("lower", "upper") if not value.is_exact else ("nominal",))
        for name, value in source_fields.items()
    })
    names = tuple(options)
    count = math.prod(len(options[name]) for name in names)
    if count > max_cases:
        raise ValueError(f"{count} pressure-driven history corners exceed max_cases={max_cases}")
    cases = []
    for combination in product(*(options[name] for name in names)):
        selected = dict(zip(names, combination))
        corner_source = release
        for name, value in source_fields.items():
            choice = selected[name]
            selected_value = (
                value.nominal if choice == "nominal"
                else value.lower if choice == "lower"
                else value.upper
            )
            corner_source = replace(
                corner_source,
                **{
                    name: BoundedValue(
                        selected_value,
                        unit=value.unit,
                        source=f"{value.source}; selected-{choice}",
                    )
                },
            )
        schedule = direct_vapour_schedule_from_pressure_driven_history(
            corner_source,
            history,
            selection={
                name: choice
                for name, choice in selected.items()
                if name in channels
            },
            **kwargs,
        )
        cases.append(FieldMeasuredFlashEnvelopeCase(tuple(sorted(selected.items())), schedule))
    return FieldMeasuredFlashEnvelope(
        tuple(cases),
        warnings=(
            "pressure-driven history envelope recomputes throat flow and flash for each coherent global pressure/temperature/liquid-fraction/opening-area/discharge-coefficient corner; it is deterministic sensitivity, not a probability interval",
        ),
    )


def request_with_measured_flash_schedule(
    request: "FieldSemiFVRequest",
    measured_schedule: FieldMeasuredFlashSchedule,
) -> "FieldSemiFVRequest":
    """Attach a measured-flash schedule without dropping its alignment limits.

    The helper is intentionally explicit: it refuses to replace an existing
    atmospheric source schedule, and the request constructor rechecks duration
    consistency against its `ReleaseSource`.
    """
    from .field_workflow import FieldSemiFVRequest

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(measured_schedule, FieldMeasuredFlashSchedule):
        raise TypeError("measured_schedule must be a FieldMeasuredFlashSchedule")
    assessment = measured_schedule.quality_assessment
    if assessment is None or not assessment.approved:
        detail = (
            "no declared quality assessment"
            if assessment is None else ", ".join(assessment.reasons)
        )
        raise ValueError(
            "measured-history schedule cannot be promoted to transient field transport without an approved quality assessment: "
            + detail
        )
    if request.direct_vapour_schedule is not None:
        raise ValueError("request already declares a direct-vapour schedule")
    provenance_warning = (
        "direct-vapour schedule was generated by measured-history flash with "
        f"selection={dict(measured_schedule.selection)!r}",
    )
    return replace(
        request,
        direct_vapour_schedule=measured_schedule.schedule,
        direct_vapour_warnings=tuple(request.direct_vapour_warnings)
        + provenance_warning + measured_schedule.warnings,
        measured_history_quality=assessment,
        measured_history_provenance=measured_schedule.provenance,
        measured_history_source_uncertainty_resolved=(
            measured_schedule.source_uncertainty_resolved
        ),
    )


def _history_property_table(
    request: "FieldSemiFVRequest",
    history: MeasuredReleaseHistory,
    *,
    table_nodes: int,
) -> tuple["LH2SaturationTable | None", tuple[str, ...]]:
    """Build one bounded table for all declared historian-temperature bounds.

    The table is built once per envelope rather than once per P/T/flow corner
    or per source interval. A supplied request table remains authoritative;
    no expanded or silent replacement table is constructed around it.
    """
    if not isinstance(table_nodes, int) or table_nodes < 4:
        raise ValueError("table_nodes must be an integer of at least four")
    if request.property_table is not None:
        return request.property_table, ()
    if request.scenario.source.fluid.strip().lower() not in {"hydrogen", "h2", "lh2"}:
        return None, (
            "LH2 saturation table was not constructed for the measured-history envelope: unsupported source fluid",
        )
    temperatures = history.temperature_k
    try:
        from .field_lh2 import build_lh2_saturation_table_for_temperature_bounds

        table_kwargs: dict[str, object] = {
            "ambient_pressure_pa": request.ambient_pressure_pa,
            "nodes": table_nodes,
        }
        pressure_bound = request.ambient_pressure_uncertainty_pa
        if pressure_bound is not None and not pressure_bound.is_exact:
            table_kwargs["ambient_pressure_bounds_pa"] = (
                pressure_bound.lower, pressure_bound.upper,
            )
        table = build_lh2_saturation_table_for_temperature_bounds(
            min(temperatures.lower), max(temperatures.upper), **table_kwargs,
        )
    except (ValueError, RuntimeError) as error:
        return None, (
            "LH2 saturation table was not constructed for the measured-history envelope; "
            f"direct property path used: {error}",
        )
    return table, (
        "one bounded LH2 saturation table covers all declared historian temperature bounds and is reused across source-history corners",
    )


def run_field_measured_history_envelope(
    request: "FieldSemiFVRequest",
    history: MeasuredReleaseHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria,
    max_cases: int = 64,
    source_id: str = "measured-history:field-envelope-direct-vapour",
    table_nodes: int = 161,
) -> FieldMeasuredHistoryScreeningEnvelope:
    """Propagate global P/T/flow bound selections into field sensor results.

    This complements :func:`run_field_semi_fv_envelope`: it keeps the declared
    weather, sensor, obstacle and transport configuration fixed while every
    uncertain measurement channel is selected globally at its lower or upper
    bound and re-flashed on the aligned flow clock.  Declared source placement
    and direction bounds are also crossed once before wind-plane projection;
    they are field geometry, not another historian source-state selection.
    It is deterministic source sensitivity, not a confidence interval and not
    a Cartesian product with the separate field-input envelope.
    """
    from .field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(history, MeasuredReleaseHistory):
        raise TypeError("history must be a MeasuredReleaseHistory")
    if not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria")
    if request.direct_vapour_schedule is not None:
        raise ValueError("measured-history field envelope requires no existing direct-vapour schedule")
    quality = assess_measured_history_quality(history, quality_criteria)
    if not quality.approved:
        raise ValueError(
            "measured-history field envelope requires approved quality criteria: "
            + ", ".join(quality.reasons)
        )
    property_table, table_warnings = _history_property_table(
        request, history, table_nodes=table_nodes,
    )
    source_envelope = direct_vapour_schedule_envelope_from_measured_history(
        request.scenario.source,
        history,
        max_cases=max_cases,
        ambient_temperature_k=request.ambient_temperature_k,
        ambient_pressure_pa=request.ambient_pressure_pa,
        property_table=property_table,
        source_id=source_id,
        quality_criteria=quality_criteria,
    )
    geometry_corner_cases = _history_source_geometry_corner_cases(request)
    total_cases = len(source_envelope.cases) * len(geometry_corner_cases)
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} measured-history/geometry corners exceed max_cases={max_cases}"
        )
    cases = []
    for source_case in source_envelope.cases:
        for geometry_values, geometry_selection in geometry_corner_cases:
            geometry_request = replace(
                request,
                scenario=_scenario_at_history_source_geometry_corner(
                    request.scenario, geometry_values,
                ),
            )
            field_request = request_with_measured_flash_schedule(
                geometry_request, source_case.schedule,
            )
            cases.append(FieldMeasuredHistoryScreeningCase(
                selection=source_case.selection + geometry_selection,
                schedule=source_case.schedule,
                result=run_field_semi_fv_screening(field_request),
            ))
    return FieldMeasuredHistoryScreeningEnvelope(
        request=request,
        quality_assessment=quality,
        cases=tuple(cases),
        warnings=source_envelope.warnings + (
            "measured-history field envelope varies global P/T/flow/liquid-fraction source selections and source location/direction corners; weather, sensor, obstacle and other field-input uncertainty must be evaluated separately",
        ) + table_warnings,
        property_table_used=property_table is not None,
    )


def _pressure_driven_field_request_with_schedule(
    request: "FieldSemiFVRequest",
    history: PressureDrivenMeasuredHistory,
    source_case: FieldMeasuredFlashEnvelopeCase,
) -> "FieldSemiFVRequest":
    """Attach one pressure-driven schedule with a positive typed bridge source.

    ``FieldSemiFVRequest`` still asks the LH2 applicability bridge to prepare one
    nominal flash before the scheduled atmospheric source is transported.  A
    pressure-only historian has no measured flow to populate that bridge, so use
    the corner's integrated throat-derived mean and its first selected P/T state.
    The attached schedule remains authoritative for time-varying injection.
    """
    from .field_workflow import FieldSemiFVRequest

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(history, PressureDrivenMeasuredHistory):
        raise TypeError("history must be a PressureDrivenMeasuredHistory")
    if not isinstance(source_case, FieldMeasuredFlashEnvelopeCase):
        raise TypeError("source_case must be FieldMeasuredFlashEnvelopeCase")
    representative_flow = source_case.schedule.total_measured_mass_kg / history.duration_s
    if not math.isfinite(representative_flow) or representative_flow <= 0.0:
        raise ValueError(
            "pressure-driven history produced a non-positive representative source flow"
        )
    selected = dict(source_case.selection)
    source_fields = {
        "opening_area_m2": request.scenario.source.opening_area_m2,
        "discharge_coefficient": request.scenario.source.discharge_coefficient,
    }

    def selected_source_value(name: str) -> BoundedValue:
        original = source_fields[name]
        choice = selected.get(name, "nominal")
        if choice not in {"nominal", "lower", "upper"}:
            raise ValueError(
                f"pressure-driven history source selection {name!r} is unsupported"
            )
        value = (
            original.nominal if choice == "nominal"
            else original.lower if choice == "lower"
            else original.upper
        )
        return BoundedValue(
            value,
            unit=original.unit,
            source=f"{original.source}; selected-{choice}",
        )
    pressure0 = history.pressure_pa.values_at(
        (history.pressure_pa.physical_time_s[0],),
        bound=selected["pressure_pa"],
    )[0]
    temperature0 = history.temperature_k.values_at(
        (history.pressure_pa.physical_time_s[0],),
        bound=selected["temperature_k"],
    )[0]
    liquid0 = (
        history.liquid_fraction.values_at(
            (history.pressure_pa.physical_time_s[0],),
            bound=selected["liquid_fraction"],
        )[0]
        if history.liquid_fraction is not None
        else request.scenario.source.liquid_fraction.nominal
    )
    field_source = replace(
        request.scenario.source,
        upstream_pressure=BoundedValue(
            pressure0, unit="Pa", source="pressure-driven-history field bridge"
        ),
        upstream_temperature=BoundedValue(
            temperature0, unit="K", source="pressure-driven-history field bridge"
        ),
        mass_flow_kg_s=BoundedValue(
            representative_flow,
            unit="kg/s",
            source=f"{source_case.schedule.schedule.source_id}; interval-integrated-history-mean",
        ),
        liquid_fraction=BoundedValue(
            liquid0, unit="1", source="pressure-driven-history field bridge"
        ),
        opening_area_m2=selected_source_value("opening_area_m2"),
        discharge_coefficient=selected_source_value("discharge_coefficient"),
    )
    bridge_request = replace(
        request,
        scenario=replace(request.scenario, source=field_source),
        direct_vapour_warnings=tuple(request.direct_vapour_warnings) + (
            "semi-FV applicability uses the interval-integrated pressure-driven history mean as a positive typed upstream bridge; the attached time-varying schedule remains authoritative for atmospheric injection",
        ),
    )
    return request_with_measured_flash_schedule(
        bridge_request, source_case.schedule,
    )


def _history_source_geometry_corner_cases(
    request: "FieldSemiFVRequest",
) -> tuple[
    tuple[dict[str, float], tuple[tuple[str, float], ...]],
    ...,
]:
    """Enumerate history source location/direction corners.

    The history schedule is independent of source placement, but the
    field transport and receptor/obstacle applicability are not.  Keep these
    geometry choices outside the throat/flash source envelope and fix them on
    the request before the schedule is attached.  This mirrors the generic
    field envelope's source-geometry treatment without crossing the history's
    upstream historian values a second time.
    """
    source = request.scenario.source
    fields: dict[str, BoundedValue] = {}
    if source.location_uncertainty_m is not None:
        fields.update({
            name: value
            for name, value in zip(
                (
                    "source_location_x_m", "source_location_y_m",
                    "source_location_z_m",
                ),
                source.location_uncertainty_m,
            )
        })
    if source.direction_uncertainty_m is not None:
        fields.update({
            name: value
            for name, value in zip(
                (
                    "source_direction_x", "source_direction_y",
                    "source_direction_z",
                ),
                source.direction_uncertainty_m,
            )
        })
    if not fields:
        return (({}, ()),)
    choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in fields.items()
    }
    names = tuple(choices)
    return tuple(
        (
            {name: float(value) for name, value in zip(names, selected)},
            tuple(
                (name, float(value))
                for name, value in zip(names, selected)
            ),
        )
        for selected in product(*(choices[name] for name in names))
    )


def _scenario_at_history_source_geometry_corner(
    scenario: FieldScenario,
    values: Mapping[str, float],
) -> FieldScenario:
    """Fix only history source placement/direction at one corner."""
    source = scenario.source
    location_names = (
        "source_location_x_m", "source_location_y_m", "source_location_z_m",
    )
    direction_names = (
        "source_direction_x", "source_direction_y", "source_direction_z",
    )
    if any(name in values for name in location_names):
        if not all(name in values for name in location_names):
            raise ValueError(
                "history source location corner must specify all three coordinates"
            )
        source = replace(
            source,
            location_m=tuple(float(values[name]) for name in location_names),
            location_uncertainty_m=None,
        )
    if any(name in values for name in direction_names):
        if not all(name in values for name in direction_names):
            raise ValueError(
                "history source direction corner must specify all three components"
            )
        source = replace(
            source,
            direction_m=tuple(float(values[name]) for name in direction_names),
            direction_uncertainty_m=None,
        )
    return replace(scenario, source=source)


def request_with_pressure_driven_history_schedule(
    request: "FieldSemiFVRequest",
    history: PressureDrivenMeasuredHistory,
    measured_schedule: FieldMeasuredFlashSchedule,
) -> "FieldSemiFVRequest":
    """Attach one approved pressure-driven schedule to a field request."""
    if not isinstance(measured_schedule, FieldMeasuredFlashSchedule):
        raise TypeError("measured_schedule must be a FieldMeasuredFlashSchedule")
    source_case = FieldMeasuredFlashEnvelopeCase(
        selection=measured_schedule.selection,
        schedule=measured_schedule,
    )
    return _pressure_driven_field_request_with_schedule(
        request, history, source_case,
    )


def run_field_pressure_driven_history_envelope(
    request: "FieldSemiFVRequest",
    history: PressureDrivenMeasuredHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria,
    max_cases: int = 64,
    source_id: str = "pressure-driven-history:field-envelope-direct-vapour",
    table_nodes: int = 161,
    allow_supercritical_gas: bool = False,
) -> FieldMeasuredHistoryScreeningEnvelope:
    """Propagate a pressure-driven orifice historian envelope into field transport.

    The pressure and temperature channels are treated as upstream boundary states,
    not as a proxy leak-rate signal.  Every source corner first recomputes the
    throat flow from the explicit opening area, discharge coefficient, ambient
    pressure and source identity, then flashes that interval into the direct-vapour
    schedule used by the semi-FV screening.  This source-only helper requires one
    fixed ambient boundary; use the joint pressure-driven history envelope below
    when ambient temperature/pressure corners must be combined with weather and
    detector uncertainty.
    """
    from .field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(history, PressureDrivenMeasuredHistory):
        raise TypeError("history must be a PressureDrivenMeasuredHistory")
    if not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria")
    if request.direct_vapour_schedule is not None:
        raise ValueError(
            "pressure-driven history field envelope requires no existing direct-vapour schedule"
        )
    if request.ambient_pressure_uncertainty_pa is not None and not request.ambient_pressure_uncertainty_pa.is_exact:
        raise ValueError(
            "pressure-driven history field envelope requires one fixed ambient pressure; "
            "use the joint pressure-driven history envelope to propagate ambient corners"
        )
    if not isinstance(allow_supercritical_gas, bool):
        raise TypeError("allow_supercritical_gas must be boolean")

    # The pressure-driven schedule derives its own mass-flow channel.  Build a
    # bounded LH2 table from the historian temperature limits exactly once, just
    # as the measured-flow path does.
    property_table, table_warnings = _history_property_table(
        request, history, table_nodes=table_nodes,
    )
    source_envelope = direct_vapour_schedule_envelope_from_pressure_driven_history(
        request.scenario.source,
        history,
        max_cases=max_cases,
        ambient_temperature_k=request.ambient_temperature_k,
        ambient_pressure_pa=request.ambient_pressure_pa,
        property_table=property_table,
        source_id=source_id,
        quality_criteria=quality_criteria,
        allow_supercritical_gas=allow_supercritical_gas,
    )
    if not source_envelope.cases:
        raise ValueError("pressure-driven history envelope produced no source cases")
    first_assessment = source_envelope.cases[0].schedule.quality_assessment
    if first_assessment is None:
        raise ValueError(
            "pressure-driven history field envelope requires an attached quality assessment"
        )
    if not first_assessment.approved:
        raise ValueError(
            "pressure-driven history field envelope requires approved quality criteria: "
            + ", ".join(first_assessment.reasons)
        )
    if any(
        case.schedule.quality_assessment != first_assessment
        for case in source_envelope.cases
    ):
        raise ValueError(
            "pressure-driven history source cases produced inconsistent quality assessments"
        )

    geometry_corner_cases = _history_source_geometry_corner_cases(request)
    total_cases = len(source_envelope.cases) * len(geometry_corner_cases)
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} pressure-driven history/geometry corners exceed max_cases={max_cases}"
        )
    cases = []
    for source_case in source_envelope.cases:
        for geometry_values, geometry_selection in geometry_corner_cases:
            geometry_request = replace(
                request,
                scenario=_scenario_at_history_source_geometry_corner(
                    request.scenario, geometry_values,
                ),
            )
            field_request = _pressure_driven_field_request_with_schedule(
                geometry_request, history, source_case,
            )
            cases.append(FieldMeasuredHistoryScreeningCase(
                selection=source_case.selection + geometry_selection,
                schedule=source_case.schedule,
                result=run_field_semi_fv_screening(field_request),
            ))
    return FieldMeasuredHistoryScreeningEnvelope(
        request=request,
        quality_assessment=first_assessment,
        cases=tuple(cases),
        warnings=source_envelope.warnings + (
            "pressure-driven history field envelope coherently varies global pressure/temperature/liquid-fraction bounds and recomputes throat flow before field transport",
            "pressure-driven opening area and discharge coefficient are included as coherent source-boundary corners before throat flow and flash recomputation",
            "pressure-driven source location and direction corners are fixed before wind-plane projection, receptor applicability and obstacle transport",
        ) + table_warnings,
        property_table_used=property_table is not None,
    )


def run_field_joint_pressure_driven_history_envelope(
    request: "FieldSemiFVRequest",
    history: PressureDrivenMeasuredHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria,
    max_cases: int = 64,
    source_id: str = "pressure-driven-history:joint-field-envelope-direct-vapour",
    table_nodes: int = 161,
    allow_supercritical_gas: bool = False,
) -> FieldJointMeasuredHistoryScreeningEnvelope:
    """Propagate pressure-driven source and independent field corners jointly.

    Unlike the source-only helper, this function includes ambient
    temperature/pressure, wind, declared stability alternatives and detector
    calibration corners.  Ambient pressure is therefore a real source-boundary
    corner: each value is passed into the throat closure before flashing, while
    the corresponding field request has that corner fixed and carries no
    unresolved ambient bound.
    """
    from .field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(history, PressureDrivenMeasuredHistory):
        raise TypeError("history must be a PressureDrivenMeasuredHistory")
    if not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria")
    if not isinstance(max_cases, int) or max_cases < 1:
        raise ValueError("max_cases must be a positive integer")
    if not isinstance(allow_supercritical_gas, bool):
        raise TypeError("allow_supercritical_gas must be boolean")
    if request.direct_vapour_schedule is not None:
        raise ValueError(
            "joint pressure-driven history envelope requires no existing direct-vapour schedule"
        )
    property_table, table_warnings = _history_property_table(
        request, history, table_nodes=table_nodes,
    )

    ambient_fields = request.ambient_uncertainty_fields()
    ambient_names = tuple(ambient_fields)
    ambient_choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in ambient_fields.items()
    }
    ambient_corner_cases = tuple(
        (
            {name: float(value) for name, value in zip(ambient_names, selected)},
            tuple((name, float(value)) for name, value in zip(ambient_names, selected)),
        )
        for selected in product(*(ambient_choices[name] for name in ambient_names))
    ) if ambient_names else (({}, ()),)

    source_envelopes = []
    for ambient_values, ambient_selection in ambient_corner_cases:
        source_envelope = direct_vapour_schedule_envelope_from_pressure_driven_history(
            request.scenario.source,
            history,
            max_cases=max_cases,
            ambient_temperature_k=ambient_values.get(
                "ambient_temperature_k", request.ambient_temperature_k,
            ),
            ambient_pressure_pa=ambient_values.get(
                "ambient_pressure_pa", request.ambient_pressure_pa,
            ),
            property_table=property_table,
            source_id=source_id,
            quality_criteria=quality_criteria,
            allow_supercritical_gas=allow_supercritical_gas,
        )
        source_envelopes.append((ambient_values, ambient_selection, source_envelope))
    first_assessment = source_envelopes[0][2].cases[0].schedule.quality_assessment
    if first_assessment is None:
        raise ValueError(
            "joint pressure-driven history envelope requires an attached quality assessment"
        )
    if not first_assessment.approved:
        raise ValueError(
            "joint pressure-driven history envelope requires approved quality criteria: "
            + ", ".join(first_assessment.reasons)
        )
    if any(
        case.schedule.quality_assessment != first_assessment
        for _ambient_values, _ambient_selection, source_envelope in source_envelopes
        for case in source_envelope.cases
    ):
        raise ValueError(
            "joint pressure-driven history source cases produced inconsistent quality assessments"
        )
    source_warnings = tuple(dict.fromkeys(
        warning
        for _ambient_values, _ambient_selection, source_envelope in source_envelopes
        for warning in source_envelope.warnings
    ))

    weather_choices = _joint_history_weather_choices(request.scenario)
    weather_names = tuple(weather_choices)
    weather_combinations = tuple(
        product(*(weather_choices[name] for name in weather_names))
    )
    stability_classes = _joint_history_stability_classes(request)
    sensor_choices = _joint_history_sensor_calibration_choices(request)
    sensor_names = tuple(sensor_choices)
    sensor_combinations = tuple(
        product(*(sensor_choices[name] for name in sensor_names))
    )
    geometry_corner_cases = _history_source_geometry_corner_cases(request)
    total_source_cases = sum(
        len(source_envelope.cases)
        for _ambient_values, _ambient_selection, source_envelope in source_envelopes
    )
    total_cases = (
        total_source_cases
        * len(geometry_corner_cases)
        * len(weather_combinations)
        * len(stability_classes)
        * len(sensor_combinations)
    )
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} joint pressure-driven-history/field corners exceed max_cases={max_cases}"
        )

    cases = []
    for ambient_values, ambient_selection, source_envelope in source_envelopes:
        for source_case in source_envelope.cases:
            for geometry_values, geometry_selection in geometry_corner_cases:
                for weather_combination in weather_combinations:
                    weather_values = {
                        name: float(value)
                        for name, value in zip(weather_names, weather_combination)
                    }
                    weather_scenario = _scenario_at_joint_history_weather_corner(
                        request.scenario, weather_values,
                    )
                    weather_scenario = _scenario_at_history_source_geometry_corner(
                        weather_scenario, geometry_values,
                    )
                    for stability in stability_classes:
                        corner_scenario = weather_scenario
                        if request.stability_alternatives is not None:
                            corner_scenario = replace(
                                weather_scenario,
                                weather=replace(weather_scenario.weather, stability=stability),
                            )
                        corner_request = replace(
                            request,
                            scenario=corner_scenario,
                            ambient_temperature_k=ambient_values.get(
                                "ambient_temperature_k", request.ambient_temperature_k,
                            ),
                            ambient_pressure_pa=ambient_values.get(
                                "ambient_pressure_pa", request.ambient_pressure_pa,
                            ),
                            ambient_air_density_kg_m3=ambient_values.get(
                                "ambient_air_density_kg_m3", request.ambient_air_density_kg_m3,
                            ),
                            ambient_temperature_uncertainty_k=None,
                            ambient_pressure_uncertainty_pa=None,
                            ambient_air_density_uncertainty_kg_m3=None,
                        )
                        field_request = _pressure_driven_field_request_with_schedule(
                            corner_request, history, source_case,
                        )
                        physical_result = run_field_semi_fv_screening(field_request)
                        for sensor_combination in sensor_combinations:
                            sensor_values = {
                                name: float(value)
                                for name, value in zip(sensor_names, sensor_combination)
                            }
                            values: dict[str, float | str] = (
                                weather_values
                                | dict(ambient_values)
                                | dict(geometry_values)
                                | sensor_values
                            )
                            if request.stability_alternatives is not None:
                                values["weather_stability"] = stability
                            field_names = (
                                *weather_names, *ambient_names,
                                *tuple(name for name, _value in geometry_selection),
                                *sensor_names,
                            )
                            if request.stability_alternatives is not None:
                                field_names += ("weather_stability",)
                            cases.append(FieldJointMeasuredHistoryScreeningCase(
                                source_selection=source_case.selection,
                                field_values=tuple((name, values[name]) for name in field_names),
                                schedule=source_case.schedule,
                                result=_joint_history_calibrated_result(physical_result, sensor_values),
                            ))
    ambient_warning = (
        "ambient temperature, pressure and air-density boundary corners are propagated through pressure-driven throat recomputation, flash and detector conversion",
    ) if ambient_names else ()
    return FieldJointMeasuredHistoryScreeningEnvelope(
        request=request,
        quality_assessment=first_assessment,
        cases=tuple(cases),
        warnings=source_warnings + (
            "joint pressure-driven envelope varies coherent historian P/T/liquid-fraction corners, explicit throat flow and independent wind/stability/detector inputs; pressure trend alone was never treated as a leak rate",
            "pressure-driven opening area and discharge coefficient are crossed as source-boundary corners before throat flow and flash recomputation",
            "pressure-driven source location and direction corners are fixed before wind-plane projection, receptor applicability and obstacle transport",
            "detector calibration corners reuse each exact pressure-driven source/weather transport field; joint envelope is deterministic sensitivity, not a probability interval or validation result",
        ) + ambient_warning + table_warnings,
        property_table_used=property_table is not None,
    )


def _joint_history_weather_choices(scenario: FieldScenario) -> dict[str, tuple[float, ...]]:
    """Return physical field inputs independent of a measured source history.

    Pressure, temperature, flow and liquid fraction are intentionally absent:
    each is replaced by the aligned historian channel before the atmospheric
    vapour schedule is built.  For measured-flow histories, static opening
    area and discharge coefficient are likewise not varied because measured
    flow replaces the orifice prediction.  The pressure-driven joint path
    crosses those two source fields in its own throat-derived source envelope
    before calling this independent-weather helper.  Surface values are not
    included in direct-vapour transport; see the public function warning.
    """
    fields = dict(scenario.weather.uncertainty_fields())
    return {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in fields.items()
    }


def _joint_history_stability_classes(request: "FieldSemiFVRequest") -> tuple[str, ...]:
    """Return declared stability alternatives without inventing a class weight."""
    alternatives = request.stability_alternatives
    if alternatives is None:
        return (request.scenario.weather.stability,)
    return alternatives.classes_for(request.scenario.weather.stability)


def _scenario_at_joint_history_weather_corner(
    scenario: FieldScenario,
    values: Mapping[str, float],
) -> FieldScenario:
    """Fix only independent weather values at one physical corner."""
    expected = set(_joint_history_weather_choices(scenario))
    if set(values) != expected:
        raise ValueError(
            "joint measured-history weather-corner keys must match independent weather inputs; "
            f"missing={sorted(expected - set(values))}, extra={sorted(set(values) - expected)}"
        )

    def exact(original: BoundedValue, name: str) -> BoundedValue:
        return BoundedValue(
            float(values[name]), unit=original.unit, source=f"{original.source}; joint-history corner",
        )

    weather = replace(
        scenario.weather,
        speed_m_s=exact(scenario.weather.speed_m_s, "wind_speed_m_s"),
        direction_deg=exact(scenario.weather.direction_deg, "wind_direction_deg"),
    )
    return replace(scenario, weather=weather)


def _joint_history_sensor_deployments(
    request: "FieldSemiFVRequest",
) -> tuple["FieldSensorDeployment", ...]:
    """Return all observation operators in a stable, namespaced order."""
    from .field_workflow import FieldSensorDeployment

    primary = () if request.scenario.sensor is None else (
        FieldSensorDeployment("field_sensor", request.scenario.sensor),
    )
    return primary + request.sensor_deployments


def _joint_history_sensor_calibration_choices(
    request: "FieldSemiFVRequest",
) -> dict[str, tuple[float, ...]]:
    """Return every observation-only calibration corner with a detector prefix."""
    return {
        f"{deployment.label}.{name}": (
            (value.nominal,) if value.is_exact else value.corners()
        )
        for deployment in _joint_history_sensor_deployments(request)
        for name, value in deployment.sensor.uncertainty_fields().items()
    }


def _sensor_at_joint_history_calibration_corner(
    sensor: SensorModel,
    label: str,
    values: Mapping[str, float],
) -> SensorModel:
    """Fix a detector calibration without changing the true scalar field."""
    fields = sensor.uncertainty_fields()

    def exact(name: str) -> BoundedValue:
        original = fields[name]
        return BoundedValue(
            float(values[f"{label}.{name}"]), unit=original.unit,
            source=f"{original.source}; joint-history sensor corner",
        )

    return replace(
        sensor,
        response_time_s=exact("sensor_response_time_s"),
        gain=exact("sensor_gain"),
        bias_mole_fraction=exact("sensor_bias_mole_fraction"),
    )


def _joint_history_calibrated_result(
    screening: "FieldSemiFVScreeningResult",
    calibration_values: Mapping[str, float],
) -> "FieldSemiFVScreeningResult":
    """Reapply detector corners to one completed physical transport result.

    The true concentration traces belong to the conserved semi-FV solve and
    are intentionally reused.  A detector which has no in-plane true trace is
    rejected instead of inheriting a nominal or an interpolated value.
    """
    from .field_observation import apply_sensor_model
    from .field_workflow import FieldSensorDeployment, FieldSensorDeploymentResult

    request = screening.request
    deployments = _joint_history_sensor_deployments(request)
    calibrated = tuple(
        FieldSensorDeployment(
            deployment.label,
            _sensor_at_joint_history_calibration_corner(
                deployment.sensor, deployment.label, calibration_values,
            ),
        )
        for deployment in deployments
    )
    selected = {deployment.label: deployment.sensor for deployment in calibrated}
    scenario = replace(request.scenario, sensor=selected.get("field_sensor"))
    calibrated_request = replace(
        request,
        scenario=scenario,
        sensor_deployments=tuple(
            deployment for deployment in calibrated if deployment.label != "field_sensor"
        ),
    )
    if not screening.completed or screening.transport is None:
        return replace(screening, request=calibrated_request)

    traces_by_label = {
        trace.receptor.label: trace for trace in screening.transport.receptor_traces
    }
    missing = [
        deployment.label for deployment in calibrated
        if deployment.label not in traces_by_label
    ]
    if missing:
        raise ValueError(
            "joint measured-history calibration envelope cannot evaluate withheld detectors: "
            + ", ".join(missing)
        )
    sensor_results = tuple(
        FieldSensorDeploymentResult(
            deployment.label,
            deployment.sensor,
            apply_sensor_model(
                traces_by_label[deployment.label].time_s,
                traces_by_label[deployment.label].concentration_kg_m3,
                deployment.sensor,
                ambient_air_density_kg_m3=calibrated_request.ambient_air_density_kg_m3,
            ),
        )
        for deployment in calibrated
    )
    warnings = tuple(
        warning for warning in screening.applicability.warnings
        if "uses nominal calibration values in this single screen" not in warning
    )
    return replace(
        screening,
        request=calibrated_request,
        applicability=replace(screening.applicability, warnings=warnings),
        sensor_trace=next(
            (item.trace for item in sensor_results if item.label == "field_sensor"), None
        ),
        sensor_results=sensor_results,
    )


def run_field_joint_measured_history_envelope(
    request: "FieldSemiFVRequest",
    history: MeasuredReleaseHistory,
    *,
    quality_criteria: MeasuredHistoryQualityCriteria,
    max_cases: int = 64,
    source_id: str = "measured-history:joint-field-envelope-direct-vapour",
    table_nodes: int = 161,
) -> FieldJointMeasuredHistoryScreeningEnvelope:
    """Jointly propagate historian and independent field-input corners.

    This is the operationally useful counterpart to the separate source and
    field envelopes.  It re-flashes global historian P/T/flow/liquid-fraction
    bounds and any declared ambient temperature/pressure corners, then combines
    them with source placement/direction, ambient air-density, wind speed/direction
    and every declared in-plane detector response/gain/bias bound.  It explicitly
    does *not* re-vary static source P/T/flow/opening/Cd/liquid-fraction values,
    which would double count information already substituted by the measured
    history.

    Physical transport is solved once per source-history/ambient/weather corner.
    Detector calibration and ambient-density corners are then applied to the
    exact true receptor traces without changing the source history. A detector
    withheld from the local wind plane is rejected rather than being assigned a
    nominal calibration or fabricated concentration.
    """
    from .field_workflow import FieldSemiFVRequest, run_field_semi_fv_screening

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(history, MeasuredReleaseHistory):
        raise TypeError("history must be a MeasuredReleaseHistory")
    if not isinstance(quality_criteria, MeasuredHistoryQualityCriteria):
        raise TypeError("quality_criteria must be MeasuredHistoryQualityCriteria")
    if not isinstance(max_cases, int) or max_cases < 1:
        raise ValueError("max_cases must be a positive integer")
    if request.direct_vapour_schedule is not None:
        raise ValueError("joint measured-history envelope requires no existing direct-vapour schedule")
    quality = assess_measured_history_quality(history, quality_criteria)
    if not quality.approved:
        raise ValueError(
            "joint measured-history envelope requires approved quality criteria: "
            + ", ".join(quality.reasons)
        )
    property_table, table_warnings = _history_property_table(
        request, history, table_nodes=table_nodes,
    )

    ambient_fields = request.ambient_uncertainty_fields()
    ambient_names = tuple(ambient_fields)
    ambient_choices = {
        name: ((value.nominal,) if value.is_exact else value.corners())
        for name, value in ambient_fields.items()
    }
    ambient_corner_cases = tuple(
        (
            {name: float(value) for name, value in zip(ambient_names, selected)},
            tuple((name, float(value)) for name, value in zip(ambient_names, selected)),
        )
        for selected in product(*(ambient_choices[name] for name in ambient_names))
    ) if ambient_names else (({}, ()),)
    source_envelopes = []
    for ambient_values, ambient_selection in ambient_corner_cases:
        source_envelope = direct_vapour_schedule_envelope_from_measured_history(
            request.scenario.source,
            history,
            max_cases=max_cases,
            ambient_temperature_k=ambient_values.get(
                "ambient_temperature_k", request.ambient_temperature_k,
            ),
            ambient_pressure_pa=ambient_values.get(
                "ambient_pressure_pa", request.ambient_pressure_pa,
            ),
            property_table=property_table,
            source_id=source_id,
            quality_criteria=quality_criteria,
        )
        source_envelopes.append((ambient_values, ambient_selection, source_envelope))
    source_warnings = tuple(dict.fromkeys(
        warning
        for _ambient_values, _ambient_selection, source_envelope in source_envelopes
        for warning in source_envelope.warnings
    ))
    weather_choices = _joint_history_weather_choices(request.scenario)
    weather_names = tuple(weather_choices)
    weather_combinations = tuple(
        product(*(weather_choices[name] for name in weather_names))
    )
    stability_classes = _joint_history_stability_classes(request)
    sensor_choices = _joint_history_sensor_calibration_choices(request)
    sensor_names = tuple(sensor_choices)
    sensor_combinations = tuple(
        product(*(sensor_choices[name] for name in sensor_names))
    )
    geometry_corner_cases = _history_source_geometry_corner_cases(request)
    total_source_cases = sum(
        len(source_envelope.cases)
        for _ambient_values, _ambient_selection, source_envelope in source_envelopes
    )
    total_cases = (
        total_source_cases
        * len(geometry_corner_cases)
        * len(weather_combinations)
        * len(stability_classes)
        * len(sensor_combinations)
    )
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} joint measured-history/field corners exceed max_cases={max_cases}"
        )

    cases = []
    for ambient_values, ambient_selection, source_envelope in source_envelopes:
        for source_case in source_envelope.cases:
            for geometry_values, geometry_selection in geometry_corner_cases:
                for weather_combination in weather_combinations:
                    weather_values = {
                        name: float(value)
                        for name, value in zip(weather_names, weather_combination)
                    }
                    weather_scenario = _scenario_at_joint_history_weather_corner(
                        request.scenario, weather_values,
                    )
                    weather_scenario = _scenario_at_history_source_geometry_corner(
                        weather_scenario, geometry_values,
                    )
                    for stability in stability_classes:
                        corner_scenario = weather_scenario
                        if request.stability_alternatives is not None:
                            corner_scenario = replace(
                                weather_scenario,
                                weather=replace(weather_scenario.weather, stability=stability),
                            )
                        corner_request = replace(
                            request,
                            scenario=corner_scenario,
                            ambient_temperature_k=ambient_values.get(
                                "ambient_temperature_k", request.ambient_temperature_k,
                            ),
                            ambient_pressure_pa=ambient_values.get(
                                "ambient_pressure_pa", request.ambient_pressure_pa,
                            ),
                            ambient_air_density_kg_m3=ambient_values.get(
                                "ambient_air_density_kg_m3", request.ambient_air_density_kg_m3,
                            ),
                            ambient_temperature_uncertainty_k=None,
                            ambient_pressure_uncertainty_pa=None,
                            ambient_air_density_uncertainty_kg_m3=None,
                        )
                        field_request = request_with_measured_flash_schedule(
                            corner_request, source_case.schedule,
                        )
                        physical_result = run_field_semi_fv_screening(field_request)
                        for sensor_combination in sensor_combinations:
                            sensor_values = {
                                name: float(value)
                                for name, value in zip(sensor_names, sensor_combination)
                            }
                            values: dict[str, float | str] = (
                                weather_values
                                | dict(ambient_values)
                                | dict(geometry_values)
                                | sensor_values
                            )
                            if request.stability_alternatives is not None:
                                values["weather_stability"] = stability
                            field_names = (
                                *weather_names, *ambient_names,
                                *tuple(name for name, _value in geometry_selection),
                                *sensor_names,
                            )
                            if request.stability_alternatives is not None:
                                field_names += ("weather_stability",)
                            cases.append(FieldJointMeasuredHistoryScreeningCase(
                                source_selection=source_case.selection,
                                field_values=tuple((name, values[name]) for name in field_names),
                                schedule=source_case.schedule,
                                result=_joint_history_calibrated_result(physical_result, sensor_values),
                            ))
    ambient_warning = (
        "ambient temperature, pressure and air-density boundary corners are propagated "
        "through measured-history re-flash and detector conversion",
    ) if ambient_names else ()
    return FieldJointMeasuredHistoryScreeningEnvelope(
        request=request,
        quality_assessment=quality,
        cases=tuple(cases),
        warnings=source_warnings + (
            "joint envelope varies history-reflashed P/T/flow/liquid-fraction channels, source location/direction, independent wind, declared stability alternatives and detector calibration inputs; source opening area and discharge coefficient are retained because measured flow replaces the orifice-flow prediction",
            "surface heat-transfer and temperature are retained as declared but do not affect a direct-vapour schedule with post-flash liquid unrouted; evaluate a phase-routing/pool source envelope before using surface uncertainty in an operational decision",
            "detector calibration corners reuse each exact source-history/weather transport field; joint envelope is deterministic sensitivity, not a probability interval or a sensor-field calibration validation",
        ) + ambient_warning + table_warnings,
        property_table_used=property_table is not None,
    )


__all__ = [
    "MeasuredTimeSeries", "MeasuredReleaseHistory", "PressureDrivenMeasuredHistory",
    "MeasuredHistoryQualityCriteria",
    "MeasuredHistoryQualityAssessment", "assess_measured_history_quality", "FieldMeasuredFlashSchedule",
    "FieldMeasuredFlashEnvelopeCase", "FieldMeasuredFlashEnvelope",
    "FieldMeasuredHistoryScreeningCase", "FieldMeasuredHistoryScreeningEnvelope",
    "FieldJointMeasuredHistoryScreeningCase", "FieldJointMeasuredHistoryScreeningEnvelope",
    "direct_vapour_schedule_from_measured_history",
    "direct_vapour_schedule_envelope_from_measured_history",
    "direct_vapour_schedule_from_pressure_driven_history",
    "direct_vapour_schedule_envelope_from_pressure_driven_history",
    "request_with_measured_flash_schedule", "request_with_pressure_driven_history_schedule",
    "run_field_measured_history_envelope",
    "run_field_joint_pressure_driven_history_envelope",
    "run_field_joint_measured_history_envelope",
]
