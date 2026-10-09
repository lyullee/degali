"""Event-balanced metrics for fixed-receptor validation.

Sensor rows are repeated observations within an experiment.  This module
keeps pooled sensor metrics for traceability, but also reports a macro score
that gives each release event equal weight.  One-sided lower-bound readings
are excluded from symmetric error metrics and are evaluated separately.

The functions here do not select a source, fit a residual, or promote a
comparison to an engineering decision.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal, Sequence


EVENT_BALANCED_METRICS_SCHEMA = "degali.event-balanced-metrics.v1"
ObservationKind = Literal["exact", "lower_bound"]


@dataclass(frozen=True)
class EventMetricObservation:
    """One observed/modelled concentration pair at one receptor."""

    event_id: str
    sensor_id: str
    observed_mole_fraction: float
    predicted_mole_fraction: float
    observation_kind: ObservationKind = "exact"

    def __post_init__(self) -> None:
        for name, value in (("event_id", self.event_id), ("sensor_id", self.sensor_id)):
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
        for name, value in (
            ("observed_mole_fraction", self.observed_mole_fraction),
            ("predicted_mole_fraction", self.predicted_mole_fraction),
        ):
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite and numeric")
            if not 0.0 <= float(value) <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")
        if self.observation_kind not in {"exact", "lower_bound"}:
            raise ValueError("observation_kind must be exact or lower_bound")


@dataclass(frozen=True)
class EventMetricSummary:
    """Metrics for one release event."""

    event_id: str
    row_count: int
    exact_count: int
    lower_bound_count: int
    mean_absolute_error_mole_fraction: float | None
    root_mean_square_error_mole_fraction: float | None
    mean_bias_model_minus_observed: float | None
    maximum_absolute_error_mole_fraction: float | None
    threshold_counts: dict[str, int]
    lower_bound_satisfied_count: int | None
    lower_bound_satisfaction_fraction: float | None


@dataclass(frozen=True)
class EventBalancedMetricScore:
    """Pooled and equal-event-weighted scores."""

    threshold_mole_fraction: float
    event_count: int
    row_count: int
    exact_row_count: int
    lower_bound_row_count: int
    pooled_mean_absolute_error_mole_fraction: float | None
    pooled_root_mean_square_error_mole_fraction: float | None
    pooled_mean_bias_model_minus_observed: float | None
    pooled_maximum_absolute_error_mole_fraction: float | None
    macro_mean_absolute_error_mole_fraction: float | None
    macro_root_mean_square_error_mole_fraction: float | None
    macro_mean_bias_model_minus_observed: float | None
    threshold_counts: dict[str, int]
    lower_bound_satisfied_count: int | None
    lower_bound_satisfaction_fraction: float | None
    events: tuple[EventMetricSummary, ...]

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-safe, schema-tagged record."""
        return {
            "schema": EVENT_BALANCED_METRICS_SCHEMA,
            "threshold_mole_fraction": self.threshold_mole_fraction,
            "counts": {
                "events": self.event_count,
                "rows": self.row_count,
                "exact_rows": self.exact_row_count,
                "lower_bound_rows": self.lower_bound_row_count,
            },
            "pooled": {
                "mean_absolute_error_mole_fraction": self.pooled_mean_absolute_error_mole_fraction,
                "root_mean_square_error_mole_fraction": self.pooled_root_mean_square_error_mole_fraction,
                "mean_bias_model_minus_observed": self.pooled_mean_bias_model_minus_observed,
                "maximum_absolute_error_mole_fraction": self.pooled_maximum_absolute_error_mole_fraction,
                "threshold_counts": dict(self.threshold_counts),
            },
            "macro_equal_event_weight": {
                "mean_absolute_error_mole_fraction": self.macro_mean_absolute_error_mole_fraction,
                "root_mean_square_error_mole_fraction": self.macro_root_mean_square_error_mole_fraction,
                "mean_bias_model_minus_observed": self.macro_mean_bias_model_minus_observed,
            },
            "lower_bound_constraint": {
                "satisfied_count": self.lower_bound_satisfied_count,
                "satisfaction_fraction": self.lower_bound_satisfaction_fraction,
            },
            "events": [
                {
                    "event_id": item.event_id,
                    "row_count": item.row_count,
                    "exact_count": item.exact_count,
                    "lower_bound_count": item.lower_bound_count,
                    "mean_absolute_error_mole_fraction": item.mean_absolute_error_mole_fraction,
                    "root_mean_square_error_mole_fraction": item.root_mean_square_error_mole_fraction,
                    "mean_bias_model_minus_observed": item.mean_bias_model_minus_observed,
                    "maximum_absolute_error_mole_fraction": item.maximum_absolute_error_mole_fraction,
                    "threshold_counts": dict(item.threshold_counts),
                    "lower_bound_satisfied_count": item.lower_bound_satisfied_count,
                    "lower_bound_satisfaction_fraction": item.lower_bound_satisfaction_fraction,
                }
                for item in self.events
            ],
        }


def _symmetric_metrics(rows: Sequence[EventMetricObservation]) -> tuple[float | None, float | None, float | None, float | None]:
    if not rows:
        return None, None, None, None
    errors = [float(row.predicted_mole_fraction - row.observed_mole_fraction) for row in rows]
    return (
        math.fsum(abs(value) for value in errors) / len(errors),
        math.sqrt(math.fsum(value * value for value in errors) / len(errors)),
        math.fsum(errors) / len(errors),
        max(abs(value) for value in errors),
    )


def _threshold_counts(rows: Sequence[EventMetricObservation], threshold: float) -> dict[str, int]:
    exact = [row for row in rows if row.observation_kind == "exact"]
    observed_positive = [row.observed_mole_fraction >= threshold for row in exact]
    predicted_positive = [row.predicted_mole_fraction >= threshold for row in exact]
    return {
        "tp": sum(actual and predicted for actual, predicted in zip(observed_positive, predicted_positive)),
        "fn": sum(actual and not predicted for actual, predicted in zip(observed_positive, predicted_positive)),
        "fp": sum(not actual and predicted for actual, predicted in zip(observed_positive, predicted_positive)),
        "tn": sum(not actual and not predicted for actual, predicted in zip(observed_positive, predicted_positive)),
        "scored_exact_rows": len(exact),
    }


def _lower_bound_metrics(rows: Sequence[EventMetricObservation]) -> tuple[int | None, float | None]:
    lower = [row for row in rows if row.observation_kind == "lower_bound"]
    if not lower:
        return None, None
    satisfied = sum(row.predicted_mole_fraction >= row.observed_mole_fraction - 1.0e-12 for row in lower)
    return satisfied, satisfied / len(lower)


def score_event_balanced_metrics(
    observations: Sequence[EventMetricObservation],
    *,
    threshold_mole_fraction: float = 0.04,
) -> EventBalancedMetricScore:
    """Score fixed-receptor pairs with pooled and equal-event weighting.

    Every ``(event_id, sensor_id)`` pair must be unique.  Exact observations
    contribute to symmetric error and threshold metrics.  Lower-bound rows are
    scored only against their one-sided constraint.
    """
    if isinstance(threshold_mole_fraction, bool) or not math.isfinite(float(threshold_mole_fraction)):
        raise ValueError("threshold_mole_fraction must be finite and numeric")
    if not 0.0 <= float(threshold_mole_fraction) <= 1.0:
        raise ValueError("threshold_mole_fraction must lie in [0, 1]")
    rows = tuple(observations)
    if not rows:
        raise ValueError("at least one metric observation is required")
    if any(not isinstance(row, EventMetricObservation) for row in rows):
        raise TypeError("observations must contain EventMetricObservation values")
    keys = [(row.event_id, row.sensor_id) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("event_id and sensor_id pairs must be unique")
    events: list[EventMetricSummary] = []
    for event_id in sorted({row.event_id for row in rows}):
        event_rows = tuple(row for row in rows if row.event_id == event_id)
        exact = tuple(row for row in event_rows if row.observation_kind == "exact")
        lower = tuple(row for row in event_rows if row.observation_kind == "lower_bound")
        mae, rmse, bias, maximum = _symmetric_metrics(exact)
        satisfied, fraction = _lower_bound_metrics(event_rows)
        events.append(EventMetricSummary(
            event_id=event_id,
            row_count=len(event_rows),
            exact_count=len(exact),
            lower_bound_count=len(lower),
            mean_absolute_error_mole_fraction=mae,
            root_mean_square_error_mole_fraction=rmse,
            mean_bias_model_minus_observed=bias,
            maximum_absolute_error_mole_fraction=maximum,
            threshold_counts=_threshold_counts(event_rows, float(threshold_mole_fraction)),
            lower_bound_satisfied_count=satisfied,
            lower_bound_satisfaction_fraction=fraction,
        ))
    exact_rows = tuple(row for row in rows if row.observation_kind == "exact")
    lower_rows = tuple(row for row in rows if row.observation_kind == "lower_bound")
    pooled = _symmetric_metrics(exact_rows)
    lower_satisfied, lower_fraction = _lower_bound_metrics(rows)
    macro_rows = [item for item in events if item.mean_absolute_error_mole_fraction is not None]
    macro_mae = (math.fsum(float(item.mean_absolute_error_mole_fraction) for item in macro_rows) / len(macro_rows)) if macro_rows else None
    macro_rmse = (math.fsum(float(item.root_mean_square_error_mole_fraction) for item in macro_rows) / len(macro_rows)) if macro_rows else None
    macro_bias = (math.fsum(float(item.mean_bias_model_minus_observed) for item in macro_rows) / len(macro_rows)) if macro_rows else None
    pooled_threshold = _threshold_counts(rows, float(threshold_mole_fraction))
    return EventBalancedMetricScore(
        threshold_mole_fraction=float(threshold_mole_fraction),
        event_count=len(events),
        row_count=len(rows),
        exact_row_count=len(exact_rows),
        lower_bound_row_count=len(lower_rows),
        pooled_mean_absolute_error_mole_fraction=pooled[0],
        pooled_root_mean_square_error_mole_fraction=pooled[1],
        pooled_mean_bias_model_minus_observed=pooled[2],
        pooled_maximum_absolute_error_mole_fraction=pooled[3],
        macro_mean_absolute_error_mole_fraction=macro_mae,
        macro_root_mean_square_error_mole_fraction=macro_rmse,
        macro_mean_bias_model_minus_observed=macro_bias,
        threshold_counts=pooled_threshold,
        lower_bound_satisfied_count=lower_satisfied,
        lower_bound_satisfaction_fraction=lower_fraction,
        events=tuple(events),
    )


__all__ = [
    "EVENT_BALANCED_METRICS_SCHEMA",
    "EventMetricObservation",
    "EventMetricSummary",
    "EventBalancedMetricScore",
    "score_event_balanced_metrics",
]
