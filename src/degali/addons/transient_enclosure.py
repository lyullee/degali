"""Bounded transient hydrogen balance for a ventilated enclosure.

The free-field plume solvers are not an enclosure flow solver.  This module
therefore provides a deliberately separate reduced-order path for measured or
declared enclosure forcing.  It advances an H2 volume fraction in an
effective mixed volume with three conservative exchanges:

* a pure-H2 source volume displaces mixture from the fixed control volume;
* measured mechanical ventilation removes the mixed concentration;
* a buoyancy-driven opening exchange scales with the square root of the
  instantaneous mixture-density deficit.

The analytic constant-forcing update keeps the volume fraction in [0, 1] and
closes the H2 mass balance at every step.  It predicts an enclosure-average
history, not a spatial sensor map and not a CFD velocity field.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np


TRANSIENT_ENCLOSURE_SCHEMA = "degali.transient-enclosure.v1"


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _vector(values: Sequence[float], name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size < 2 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite one-dimensional vector")
    return result


@dataclass(frozen=True)
class EnclosureForcing:
    """Common-clock atmospheric source and mechanical-ventilation histories."""

    time_s: Sequence[float]
    h2_mass_rate_kg_s: Sequence[float]
    ventilation_flow_m3_s: Sequence[float]
    source_id: str = "declared"

    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        time = _vector(self.time_s, "forcing time_s")
        rate = _vector(self.h2_mass_rate_kg_s, "forcing h2_mass_rate_kg_s")
        ventilation = _vector(
            self.ventilation_flow_m3_s, "forcing ventilation_flow_m3_s"
        )
        if rate.shape != time.shape or ventilation.shape != time.shape:
            raise ValueError("enclosure forcing arrays must share one clock")
        if np.any(np.diff(time) <= 0.0):
            raise ValueError("enclosure forcing time must be strictly increasing")
        if np.any(rate < 0.0) or np.any(ventilation < 0.0):
            raise ValueError("enclosure source and ventilation must be non-negative")
        if not isinstance(self.source_id, str) or not self.source_id.strip():
            raise ValueError("enclosure forcing source_id must be non-empty")
        return time, rate, ventilation


@dataclass(frozen=True)
class TransientEnclosureConfig:
    """Physical geometry and declared reduced-order closure coefficients."""

    enclosure_volume_m3: float
    effective_mixed_volume_fraction: float
    effective_opening_area_m2: float
    transport_delay_s: float
    opening_height_m: float
    hydrogen_density_kg_m3: float = 0.08375
    ambient_density_kg_m3: float = 1.225
    gravity_m_s2: float = 9.81

    def __post_init__(self) -> None:
        for name in (
            "enclosure_volume_m3",
            "effective_mixed_volume_fraction",
            "effective_opening_area_m2",
            "transport_delay_s",
            "opening_height_m",
            "hydrogen_density_kg_m3",
            "ambient_density_kg_m3",
            "gravity_m_s2",
        ):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.enclosure_volume_m3 <= 0.0 or self.opening_height_m <= 0.0:
            raise ValueError("enclosure volume and opening height must be positive")
        if not 0.0 < self.effective_mixed_volume_fraction <= 1.0:
            raise ValueError("effective mixed-volume fraction must lie in (0, 1]")
        if self.effective_opening_area_m2 < 0.0 or self.transport_delay_s < 0.0:
            raise ValueError("effective opening area and transport delay cannot be negative")
        if self.hydrogen_density_kg_m3 <= 0.0 or self.ambient_density_kg_m3 <= 0.0:
            raise ValueError("hydrogen and ambient densities must be positive")
        if self.hydrogen_density_kg_m3 >= self.ambient_density_kg_m3:
            raise ValueError("this buoyant-H2 closure requires hydrogen lighter than ambient")
        if self.gravity_m_s2 <= 0.0:
            raise ValueError("gravity must be positive")

    @property
    def effective_mixed_volume_m3(self) -> float:
        return self.enclosure_volume_m3 * self.effective_mixed_volume_fraction


@dataclass(frozen=True)
class TransientEnclosureResult:
    """Bounded enclosure-average history and conservative mass audit."""

    schema: str
    time_s: np.ndarray
    hydrogen_volume_fraction: np.ndarray
    source_volume_flow_m3_s: np.ndarray
    ventilation_flow_m3_s: np.ndarray
    buoyant_opening_exchange_m3_s: np.ndarray
    inventory_kg: np.ndarray
    cumulative_injected_kg: np.ndarray
    cumulative_removed_kg: np.ndarray
    mass_residual_kg: np.ndarray
    diagnostics: dict[str, object]

    @property
    def concentration_vol_percent(self) -> np.ndarray:
        return 100.0 * self.hydrogen_volume_fraction

    @property
    def maximum_mass_residual_kg(self) -> float:
        return float(np.max(np.abs(self.mass_residual_kg)))


def solve_transient_enclosure(
    forcing: EnclosureForcing,
    evaluation_time_s: Sequence[float],
    config: TransientEnclosureConfig,
) -> TransientEnclosureResult:
    """Advance the bounded enclosure balance on a requested event clock."""

    if not isinstance(forcing, EnclosureForcing):
        raise TypeError("forcing must be EnclosureForcing")
    if not isinstance(config, TransientEnclosureConfig):
        raise TypeError("config must be TransientEnclosureConfig")
    forcing_time, mass_rate, ventilation = forcing.arrays()
    time = _vector(evaluation_time_s, "evaluation_time_s")
    if np.any(np.diff(time) <= 0.0) or time[0] < 0.0:
        raise ValueError("evaluation time must be non-negative and strictly increasing")
    if time[-1] > forcing_time[-1] + 1.0e-12:
        raise ValueError("forcing history must cover the complete evaluation time")
    if time[0] < max(0.0, forcing_time[0]) - 1.0e-12:
        raise ValueError("forcing history must cover the evaluation start")

    count = time.size
    alpha = np.zeros(count, dtype=float)
    source_volume = np.zeros(count, dtype=float)
    mechanical = np.zeros(count, dtype=float)
    buoyant = np.zeros(count, dtype=float)
    inventory = np.zeros(count, dtype=float)
    injected = np.zeros(count, dtype=float)
    removed = np.zeros(count, dtype=float)
    residual = np.zeros(count, dtype=float)
    volume = config.effective_mixed_volume_m3
    density_contrast = (
        config.ambient_density_kg_m3 - config.hydrogen_density_kg_m3
    ) / config.ambient_density_kg_m3

    for index in range(1, count):
        dt = float(time[index] - time[index - 1])
        midpoint = 0.5 * float(time[index] + time[index - 1])
        delayed_source_time = midpoint - config.transport_delay_s
        source_mass_rate = 0.0
        if delayed_source_time >= 0.0:
            source_mass_rate = float(np.interp(
                delayed_source_time, forcing_time, mass_rate
            ))
        source_q = source_mass_rate / config.hydrogen_density_kg_m3
        ventilation_q = float(np.interp(midpoint, forcing_time, ventilation))
        opening_q = config.effective_opening_area_m2 * math.sqrt(
            config.gravity_m_s2
            * config.opening_height_m
            * density_contrast
            * max(alpha[index - 1], 0.0)
        )
        total_outflow_q = source_q + ventilation_q + opening_q
        decay = total_outflow_q / volume
        equilibrium = 0.0 if total_outflow_q <= 0.0 else source_q / total_outflow_q
        if decay <= 0.0:
            next_alpha = alpha[index - 1]
            integral_alpha_dt = alpha[index - 1] * dt
        else:
            exponential = math.exp(-decay * dt)
            next_alpha = equilibrium + (alpha[index - 1] - equilibrium) * exponential
            integral_alpha_dt = (
                equilibrium * dt
                + (alpha[index - 1] - equilibrium) * (1.0 - exponential) / decay
            )
        alpha[index] = min(1.0, max(0.0, next_alpha))
        source_volume[index] = source_q
        mechanical[index] = ventilation_q
        buoyant[index] = opening_q
        injected[index] = (
            injected[index - 1]
            + config.hydrogen_density_kg_m3 * source_q * dt
        )
        removed[index] = (
            removed[index - 1]
            + config.hydrogen_density_kg_m3 * total_outflow_q * integral_alpha_dt
        )
        inventory[index] = config.hydrogen_density_kg_m3 * volume * alpha[index]
        residual[index] = inventory[index] + removed[index] - injected[index]

    return TransientEnclosureResult(
        schema=TRANSIENT_ENCLOSURE_SCHEMA,
        time_s=time,
        hydrogen_volume_fraction=alpha,
        source_volume_flow_m3_s=source_volume,
        ventilation_flow_m3_s=mechanical,
        buoyant_opening_exchange_m3_s=buoyant,
        inventory_kg=inventory,
        cumulative_injected_kg=injected,
        cumulative_removed_kg=removed,
        mass_residual_kg=residual,
        diagnostics={
            "source_id": forcing.source_id,
            "effective_mixed_volume_m3": volume,
            "bounded_volume_fraction": True,
            "source_displacement_included": True,
            "opening_exchange_closure": "A_eff*sqrt(g*H*mixture_density_deficit)",
            "prediction_scope": "enclosure-average transient concentration; not a spatial sensor map or CFD field",
        },
    )


__all__ = [
    "TRANSIENT_ENCLOSURE_SCHEMA",
    "EnclosureForcing",
    "TransientEnclosureConfig",
    "TransientEnclosureResult",
    "solve_transient_enclosure",
]
