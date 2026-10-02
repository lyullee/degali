"""Explicit substrate-to-pool heat balance for a quiescent LH2 pool.

The atmospheric pool path accepts a vapour mass rate.  This module supplies
that rate from a declared, one-dimensional substrate boundary instead of
silently treating all released liquid as immediate vapour.  It is deliberately
small: it resolves conduction in a homogeneous solid beneath a pool held at
the hydrogen saturation temperature.  Pool spread, droplet trajectories,
water freezing, boiling nucleation and wind-driven surface transfer are *not*
inferred here.  A caller may, however, supply an audited liquid inflow rate.

The model is most useful when a known liquid inventory or rainout rate covers
a known, approximately fixed footprint.  The E3.4 reader in
``degali.validation.preslhy_e34`` is an independent observation operator for
that regime; neither module contains campaign data or calibrated material
constants.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class SolidSubstrate:
    """Declared homogeneous solid properties for the pool heat boundary.

    ``depth_m`` is the numerical depth, whose far face remains at
    ``initial_temperature_k``.  It therefore represents a finite
    approximation to a semi-infinite substrate, and must be made deeper than
    the penetration depth ``sqrt(alpha * duration)`` before treating a run as
    insensitive to that boundary.
    """

    conductivity_w_m_k: float
    density_kg_m3: float
    heat_capacity_j_kg_k: float
    initial_temperature_k: float
    depth_m: float
    cells: int = 48

    def __post_init__(self) -> None:
        for name in (
            "conductivity_w_m_k", "density_kg_m3", "heat_capacity_j_kg_k",
            "initial_temperature_k", "depth_m",
        ):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.cells < 3:
            raise ValueError("cells must be at least three")

    @property
    def diffusivity_m2_s(self) -> float:
        return self.conductivity_w_m_k / (
            self.density_kg_m3 * self.heat_capacity_j_kg_k
        )


@dataclass(frozen=True)
class PoolEvaporationStep:
    """One post-rainout pool state, evaluated at the end of its timestep."""

    elapsed_s: float
    heat_flux_w_m2: float
    vapour_rate_kg_s: float
    evaporated_mass_kg: float
    remaining_liquid_kg: float | None
    surface_adjacent_temperature_k: float
    liquid_inflow_rate_kg_s: float = 0.0
    cumulative_liquid_inflow_kg: float = 0.0


@dataclass(frozen=True)
class PoolEvaporationResult:
    """Time-resolved substrate heat balance and its physical audit."""

    steps: tuple[PoolEvaporationStep, ...]
    area_m2: float
    saturation_temperature_k: float
    latent_heat_j_kg: float
    substrate: SolidSubstrate
    finite_depth_ratio: float

    @property
    def total_evaporated_mass_kg(self) -> float:
        return self.steps[-1].evaporated_mass_kg if self.steps else 0.0

    @property
    def mean_vapour_rate_kg_s(self) -> float:
        if not self.steps or self.steps[-1].elapsed_s <= 0.0:
            return 0.0
        return self.total_evaporated_mass_kg / self.steps[-1].elapsed_s

    @property
    def total_liquid_inflow_kg(self) -> float:
        return self.steps[-1].cumulative_liquid_inflow_kg if self.steps else 0.0

    @property
    def depth_is_effectively_semi_infinite(self) -> bool:
        """Whether the declared numerical depth exceeds five penetration depths."""
        return self.finite_depth_ratio >= 5.0


def semi_infinite_heat_flux(
    substrate: SolidSubstrate,
    *,
    elapsed_s: float,
    saturation_temperature_k: float = 20.27,
) -> float:
    """Analytic solid-to-pool flux after a step to a constant surface temperature.

    This is provided as a mesh-independent check on early-time numerical
    solutions.  It diverges at exactly zero time, so callers must give a
    positive elapsed time and use a time average for a finite first step.
    """

    if elapsed_s <= 0.0:
        raise ValueError("elapsed_s must be positive")
    if saturation_temperature_k <= 0.0:
        raise ValueError("saturation_temperature_k must be positive")
    delta_t = substrate.initial_temperature_k - saturation_temperature_k
    return (
        substrate.conductivity_w_m_k * delta_t
        / math.sqrt(math.pi * substrate.diffusivity_m2_s * elapsed_s)
    )


def _tridiagonal_solve(lower: np.ndarray, diagonal: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve a tridiagonal system without turning each thermal step into a dense solve."""

    diagonal = diagonal.copy()
    rhs = rhs.copy()
    for index in range(1, len(diagonal)):
        factor = lower[index - 1] / diagonal[index - 1]
        diagonal[index] -= factor * upper[index - 1]
        rhs[index] -= factor * rhs[index - 1]
    out = np.empty_like(rhs)
    out[-1] = rhs[-1] / diagonal[-1]
    for index in range(len(diagonal) - 2, -1, -1):
        out[index] = (rhs[index] - upper[index] * out[index + 1]) / diagonal[index]
    return out


def substrate_conduction_evaporation(
    substrate: SolidSubstrate,
    *,
    area_m2: float,
    duration_s: float,
    time_step_s: float,
    saturation_temperature_k: float = 20.27,
    latent_heat_j_kg: float = 4.46e5,
    initial_liquid_mass_kg: float | None = None,
    liquid_inflow_rate_kg_s: float = 0.0,
    inflow_duration_s: float = 0.0,
) -> PoolEvaporationResult:
    """Integrate a declared post-rainout LH2 pool over a solid substrate.

    The substrate equation is solved by backward Euler on a cell-centred
    one-dimensional grid.  The surface is held at the supplied saturation
    temperature while liquid remains; the deep face stays at the declared
    initial temperature.  Consequently the result is a conduction-only
    lower-complexity boundary, not a fitted total evaporation correlation.

    If ``initial_liquid_mass_kg`` is supplied, evaporation stops on depletion
    rather than continuing to create vapour from an absent pool.  A constant
    declared rainout/feed rate may be active for ``inflow_duration_s``; this
    permits concurrent fixed-footprint pool formation and evaporation without
    inventing a spreading law.
    """

    if area_m2 <= 0.0 or duration_s <= 0.0 or time_step_s <= 0.0:
        raise ValueError("area_m2, duration_s and time_step_s must be positive")
    if saturation_temperature_k <= 0.0 or latent_heat_j_kg <= 0.0:
        raise ValueError("saturation_temperature_k and latent_heat_j_kg must be positive")
    if initial_liquid_mass_kg is not None and initial_liquid_mass_kg < 0.0:
        raise ValueError("initial_liquid_mass_kg cannot be negative")
    if not math.isfinite(liquid_inflow_rate_kg_s) or liquid_inflow_rate_kg_s < 0.0:
        raise ValueError("liquid_inflow_rate_kg_s must be finite and non-negative")
    if not math.isfinite(inflow_duration_s) or inflow_duration_s < 0.0:
        raise ValueError("inflow_duration_s must be finite and non-negative")
    if liquid_inflow_rate_kg_s > 0.0 and inflow_duration_s <= 0.0:
        raise ValueError("positive liquid inflow requires a positive duration")
    if liquid_inflow_rate_kg_s > 0.0 and initial_liquid_mass_kg is None:
        raise ValueError("liquid inflow requires a finite initial pool inventory")

    n = substrate.cells
    dz = substrate.depth_m / n
    alpha = substrate.diffusivity_m2_s
    temperatures = np.full(n, substrate.initial_temperature_k, dtype=float)
    elapsed = 0.0
    evaporated = 0.0
    remaining = initial_liquid_mass_kg
    cumulative_inflow = 0.0
    out: list[PoolEvaporationStep] = []

    while elapsed < duration_s - 1.0e-12:
        dt = min(time_step_s, duration_s - elapsed)
        inflow_interval = max(
            0.0, min(elapsed + dt, inflow_duration_s) - elapsed
        )
        inflow_mass = liquid_inflow_rate_kg_s * inflow_interval
        cumulative_inflow += inflow_mass
        if remaining is not None:
            remaining += inflow_mass
        liquid_present = remaining is None or remaining > 0.0
        if liquid_present:
            r = alpha * dt / (dz * dz)
            diagonal = np.full(n, 1.0 + 2.0 * r)
            lower = np.full(n - 1, -r)
            upper = np.full(n - 1, -r)
            # Surface-to-first-centre distance is dz/2; the far face is kept
            # at the initial state.  Both give the 3r/2r end stencil.
            diagonal[0] = 1.0 + 3.0 * r
            diagonal[-1] = 1.0 + 3.0 * r
            rhs = temperatures.copy()
            rhs[0] += 2.0 * r * saturation_temperature_k
            rhs[-1] += 2.0 * r * substrate.initial_temperature_k
            next_temperatures = _tridiagonal_solve(lower, diagonal, upper, rhs)
            flux = max(
                2.0 * substrate.conductivity_w_m_k
                * (next_temperatures[0] - saturation_temperature_k) / dz,
                0.0,
            )
            potential = flux * area_m2 * dt / latent_heat_j_kg
            if remaining is not None and potential > remaining:
                potential = remaining
                flux = potential * latent_heat_j_kg / (area_m2 * dt)
            temperatures = next_temperatures
        else:
            flux = 0.0
            potential = 0.0

        evaporated += potential
        if remaining is not None:
            remaining = max(remaining - potential, 0.0)
        elapsed += dt
        out.append(PoolEvaporationStep(
            elapsed_s=float(elapsed), heat_flux_w_m2=float(flux),
            vapour_rate_kg_s=float(potential / dt),
            evaporated_mass_kg=float(evaporated), remaining_liquid_kg=remaining,
            surface_adjacent_temperature_k=float(temperatures[0]),
            liquid_inflow_rate_kg_s=float(inflow_mass / dt),
            cumulative_liquid_inflow_kg=float(cumulative_inflow),
        ))

    penetration = math.sqrt(alpha * duration_s)
    return PoolEvaporationResult(
        steps=tuple(out), area_m2=float(area_m2),
        saturation_temperature_k=float(saturation_temperature_k),
        latent_heat_j_kg=float(latent_heat_j_kg), substrate=substrate,
        finite_depth_ratio=float(substrate.depth_m / max(penetration, 1.0e-30)),
    )


__all__ = [
    "SolidSubstrate", "PoolEvaporationStep", "PoolEvaporationResult",
    "semi_infinite_heat_flux", "substrate_conduction_evaporation",
]
