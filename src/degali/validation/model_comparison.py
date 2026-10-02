"""Apples-to-apples comparison of dispersion models.

Model-to-model claims are easy to overstate when each program receives a
different source term, wind direction, averaging window, or receptor
operator.  This module makes those choices explicit and refuses to rank
predictions whose case definitions are not compatible.  It deliberately
does not ship HyRAM, PHAST, EFFECTS, or any third-party observations; a
caller supplies predictions produced by the programs it is licensed to run.

The report can therefore be generated in a public repository without
redistributing controlled software or experimental data.  A result is still
useful when only one external prediction is available: the paired statistics
are retained, while the absent ranking is reported honestly.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Mapping, Sequence

import numpy as np

from .statistics import Statistics, statistics


@dataclass(frozen=True)
class ComparisonCase:
    """Frozen metadata that must match before model ranking is allowed.

    ``wind_direction_rad`` is a *to* direction in global coordinates.  The
    receptor operator describes how a model value was reduced at the sensor
    (for example ``centreline_1s_peak`` versus ``arc_max_60s_mean``); it is
    part of the physical comparison, not presentation metadata.
    """

    source_mode: str
    source_rate_kg_s: float
    wind_speed_m_s: float
    wind_direction_rad: float
    release_height_m: float
    receptor_operator: str
    averaging_time_s: float
    phase_closure: str = "declared"
    geometry: str = "open_horizontal"

    def __post_init__(self) -> None:
        if not isinstance(self.source_mode, str) or not self.source_mode.strip():
            raise ValueError("source_mode must be a non-empty string")
        if not isinstance(self.receptor_operator, str) or not self.receptor_operator.strip():
            raise ValueError("receptor_operator must be a non-empty string")
        if not isinstance(self.phase_closure, str) or not self.phase_closure.strip():
            raise ValueError("phase_closure must be a non-empty string")
        if not isinstance(self.geometry, str) or not self.geometry.strip():
            raise ValueError("geometry must be a non-empty string")
        values = (
            self.source_rate_kg_s,
            self.wind_speed_m_s,
            self.wind_direction_rad,
            self.release_height_m,
            self.averaging_time_s,
        )
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("case numerical fields must be finite")
        if self.source_rate_kg_s <= 0.0:
            raise ValueError("source_rate_kg_s must be positive")
        if self.wind_speed_m_s < 0.0:
            raise ValueError("wind_speed_m_s must be non-negative")
        if self.release_height_m < 0.0:
            raise ValueError("release_height_m must be non-negative")
        if self.averaging_time_s < 0.0:
            raise ValueError("averaging_time_s must be non-negative")


@dataclass(frozen=True)
class ModelPrediction:
    """One model output paired with the exact case that produced it."""

    name: str
    values: tuple[float, ...] | Sequence[float]
    case: ComparisonCase
    calibrated_to_observations: bool = False
    provenance: str = "unspecified"

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("model name must be a non-empty string")
        values = np.asarray(self.values, dtype=float)
        if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
            raise ValueError("model values must be a finite, non-empty 1-D sequence")
        object.__setattr__(self, "values", tuple(float(value) for value in values))


@dataclass(frozen=True)
class ModelComparisonReport:
    """Paired metrics plus an auditable ranking boundary."""

    reference_case: ComparisonCase
    results: dict[str, Statistics]
    eligible_models: tuple[str, ...]
    excluded_models: dict[str, str]
    warnings: tuple[str, ...]

    @property
    def ranking_allowed(self) -> bool:
        """Whether at least two compatible, uncalibrated models can be ranked."""
        return len(self.eligible_models) >= 2

    @property
    def best_fac2_model(self) -> str | None:
        """Return the best eligible FAC2 model, or ``None`` if ranking is barred."""
        if not self.ranking_allowed:
            return None
        return max(self.eligible_models, key=lambda name: self.results[name].fac2)

    def as_dict(self) -> dict:
        """Return JSON-ready metadata and metrics without raw observations."""
        return {
            "reference_case": asdict(self.reference_case),
            "results": {name: asdict(item) for name, item in self.results.items()},
            "eligible_models": list(self.eligible_models),
            "excluded_models": dict(self.excluded_models),
            "warnings": list(self.warnings),
            "ranking_allowed": self.ranking_allowed,
            "best_fac2_model": self.best_fac2_model,
        }


def _case_mismatches(reference: ComparisonCase, candidate: ComparisonCase) -> list[str]:
    mismatches: list[str] = []
    text_fields = ("source_mode", "receptor_operator", "phase_closure", "geometry")
    for field in text_fields:
        if getattr(reference, field) != getattr(candidate, field):
            mismatches.append(field)
    numeric_fields = (
        "source_rate_kg_s", "wind_speed_m_s", "wind_direction_rad",
        "release_height_m", "averaging_time_s",
    )
    for field in numeric_fields:
        left, right = float(getattr(reference, field)), float(getattr(candidate, field))
        if not math.isclose(left, right, rel_tol=1.0e-9, abs_tol=1.0e-12):
            mismatches.append(field)
    return mismatches


def compare_models(
    observed: Sequence[float],
    predictions: Sequence[ModelPrediction],
    *,
    case: ComparisonCase,
    floor: float | None = None,
) -> ModelComparisonReport:
    """Compare predictions and rank only compatible, uncalibrated outputs.

    All paired statistics are returned, including excluded models, so a
    reviewer can see exactly what was computed.  A prediction is excluded
    from ranking when its case metadata differs, it was calibrated to the
    supplied observations, or its values cannot be paired.  This is a
    comparison protocol, not a license to claim universal superiority.
    """
    obs = np.asarray(observed, dtype=float)
    if obs.ndim != 1 or obs.size == 0 or not np.all(np.isfinite(obs)):
        raise ValueError("observed values must be a finite, non-empty 1-D sequence")
    if floor is not None and (not math.isfinite(float(floor)) or floor < 0.0):
        raise ValueError("floor must be finite and non-negative")
    if not predictions:
        raise ValueError("at least one model prediction is required")
    names = [prediction.name for prediction in predictions]
    if len(set(names)) != len(names):
        raise ValueError("model names must be unique")

    results: dict[str, Statistics] = {}
    eligible: list[str] = []
    excluded: dict[str, str] = {}
    warnings: list[str] = []
    for prediction in predictions:
        values = np.asarray(prediction.values, dtype=float)
        reason: list[str] = []
        mismatches = _case_mismatches(case, prediction.case)
        if mismatches:
            reason.append("case mismatch: " + ", ".join(mismatches))
        if prediction.calibrated_to_observations:
            reason.append("prediction was calibrated to these observations")
        if values.shape != obs.shape:
            reason.append(f"length mismatch: observed={obs.size}, prediction={values.size}")
        if values.shape == obs.shape:
            try:
                results[prediction.name] = statistics(obs, values, floor=floor)
            except ValueError as exc:
                reason.append(f"statistics unavailable: {exc}")
        if reason:
            excluded[prediction.name] = "; ".join(reason)
            warnings.append(f"{prediction.name} excluded from ranking: {excluded[prediction.name]}")
        else:
            eligible.append(prediction.name)

    if len(eligible) < 2:
        warnings.append(
            "fewer than two compatible uncalibrated models are available; "
            "no superiority ranking is permitted"
        )
    return ModelComparisonReport(
        reference_case=case,
        results=results,
        eligible_models=tuple(eligible),
        excluded_models=excluded,
        warnings=tuple(warnings),
    )


__all__ = [
    "ComparisonCase", "ModelPrediction", "ModelComparisonReport", "compare_models",
]
