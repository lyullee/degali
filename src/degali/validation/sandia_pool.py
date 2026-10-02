"""Bounded validation of steady cross-wind LH2 pool contours.

The Sandia 2025 report publishes contour *figures*, not the underlying sensor
time series.  A contour is consequently retained as a one-sided geometric
constraint: when an ``x_H2`` contour is visibly present at a down-wind plane,
the real centreline must reach at least that plane.  A model whose centreline
falls below that fraction earlier is falsified; a model reaching farther is
not thereby validated.

No experiment values are embedded here.  Callers load a local reduction from
an ignored path, record the figure/panel and its one-grid resolution, and
keep the original report outside the public package.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PoolContourLowerBound:
    """A visibly resolved contour gives a lower bound on down-wind reach."""

    identifier: str
    figure: str
    panel: str
    contour_mole_fraction: float
    minimum_reach_m: float
    spatial_resolution_m: float
    all_inflow_vaporization_supported: bool
    source_note: str

    def __post_init__(self):
        if not self.identifier or not self.figure or not self.panel or not self.source_note:
            raise ValueError("identifier, figure, panel and source note are required")
        if not 0.0 < self.contour_mole_fraction < 1.0:
            raise ValueError("contour mole fraction must lie strictly between zero and one")
        if self.minimum_reach_m <= 0.0 or self.spatial_resolution_m <= 0.0:
            raise ValueError("reach and spatial resolution must be positive")


@dataclass(frozen=True)
class PoolContourScreen:
    """One-sided comparison; it is never an accuracy score or promotion gate."""

    identifier: str
    model_lfl_reach_m: float
    observed_minimum_reach_m: float
    spatial_resolution_m: float
    outcome: str
    source_assumption_supported: bool
    field_accuracy_claim_allowed: bool = False
    default_promotion_allowed: bool = False


def screen_contour_lower_bound(
    observation: PoolContourLowerBound, model_lfl_reach_m: float,
) -> PoolContourScreen:
    """Reject only a model that ends before the resolved contour plane.

    ``model_lfl_reach_m`` must be the model centreline distance to the same
    mole fraction.  The observed contour can lie off centreline and is based
    on finite thermocouple interpolation, so equality is not a target.  One
    grid spacing is allowed only for the drawing/readout resolution; it is
    not a fitted tolerance.
    """
    value = float(model_lfl_reach_m)
    if value < 0.0 or value != value or value == float("inf"):
        raise ValueError("model LFL reach must be finite and nonnegative")
    threshold = observation.minimum_reach_m - observation.spatial_resolution_m
    outcome = "falsified_below_resolved_reach" if value < threshold else "not_falsified"
    return PoolContourScreen(
        identifier=observation.identifier,
        model_lfl_reach_m=value,
        observed_minimum_reach_m=observation.minimum_reach_m,
        spatial_resolution_m=observation.spatial_resolution_m,
        outcome=outcome,
        source_assumption_supported=observation.all_inflow_vaporization_supported,
    )


def read_pool_contour_lower_bounds(path: str | Path) -> tuple[PoolContourLowerBound, ...]:
    """Read a local, user-controlled contour reduction; do not ship it."""
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(raw, list) or not raw:
        raise ValueError("contour reduction must be a non-empty JSON list")
    rows = tuple(PoolContourLowerBound(**row) for row in raw)
    identifiers = [row.identifier for row in rows]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("contour-reduction identifiers must be unique")
    return rows


def screen_pool_contours(
    observations: Iterable[PoolContourLowerBound], predictions: dict[str, float],
) -> tuple[PoolContourScreen, ...]:
    """Require exactly one declared model reach for every local observation."""
    rows = tuple(observations)
    if set(predictions) != {row.identifier for row in rows}:
        raise ValueError("predictions must match the contour-reduction identifiers exactly")
    return tuple(screen_contour_lower_bound(row, predictions[row.identifier]) for row in rows)


def screen_report(screens: Iterable[PoolContourScreen]) -> dict:
    """Return a JSON-safe no-score report for the local audit tool."""
    rows = tuple(screens)
    return {
        "third_party_observations_distributed": False,
        "accuracy_score_calculated": False,
        "default_promotion_allowed": False,
        "all_source_assumptions_supported": all(row.source_assumption_supported for row in rows),
        "all_lower_bounds_not_falsified": all(row.outcome == "not_falsified" for row in rows),
        "screens": [asdict(row) for row in rows],
    }


__all__ = [
    "PoolContourLowerBound", "PoolContourScreen", "screen_contour_lower_bound",
    "read_pool_contour_lower_bounds", "screen_pool_contours", "screen_report",
]
