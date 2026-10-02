"""Conservative finite-release handoff and steady-wind applicability checks.

DEGALI's jet/crosswind solvers are steady spatial models.  A finite release
must stop using an indefinitely fed plume once the material travel clock
reaches the source duration.  This module locates that point and converts the
full-cross-section fluxes into a finite inventory without importing SLAB's
quarter-cloud bookkeeping.

The result is deliberately a *puff handoff*, not a hidden puff solver.  It
prevents steady-plume extrapolation beyond its valid source clock and exposes
all mass, momentum and energy totals required by a transient continuation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from .transient_receptor import WindHistory


@dataclass(frozen=True)
class FiniteReleasePuffHandoff:
    """Full-cross-section inventory at the finite plume-to-puff transition."""

    status: str
    source_duration_s: float
    transition_time_s: float | None
    transition_arc_length_m: float | None
    centre_position_m: tuple[float, float, float] | None
    transverse_widths_m: tuple[float, float] | None
    longitudinal_half_length_m: float | None
    bulk_velocity_m_s: tuple[float, float, float] | None
    total_mass_kg: float | None
    hydrogen_mass_kg: float | None
    momentum_kg_m_s: tuple[float, float, float] | None
    relative_energy_j: float | None
    hydrogen_flux_residual_kg_s: float | None
    hydrogen_inventory_residual_kg: float | None
    requires_transient_puff_continuation: bool
    qualification: str


@dataclass(frozen=True)
class SteadyWindApplicability:
    """Evidence-based gate for applying one steady/yawed plume to a time window."""

    applicable: bool
    start_s: float
    end_s: float
    mean_speed_m_s: float
    speed_range_fraction: float
    direction_span_deg: float
    maximum_speed_range_fraction: float
    maximum_direction_span_deg: float
    reasons: tuple[str, ...]


def _result_array(result: Any, name: str) -> np.ndarray:
    value = result[name] if isinstance(result, dict) else getattr(result, name)
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"result {name} contains non-finite values")
    return array


def _transition_layout(
    model: Any, states: np.ndarray, fluxes: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return positions, momenta, energies and energy-column index."""
    if states.ndim != 2 or fluxes.ndim != 2 or len(states) != len(fluxes):
        raise ValueError("states and fluxes must be aligned two-dimensional arrays")
    if states.shape[1] == 7 and fluxes.shape[1] == 5:
        positions = np.column_stack((states[:, 5], np.zeros(len(states)), states[:, 6]))
        momenta = np.column_stack((fluxes[:, 2], np.zeros(len(states)), fluxes[:, 3]))
        return positions, momenta, fluxes[:, 4], 4
    if states.shape[1] == 10 and fluxes.shape[1] == 6:
        positions = states[:, 7:10]
        momenta = fluxes[:, 2:5]
        return positions, momenta, fluxes[:, 5], 5
    raise ValueError(
        "expected a seven-state independent-energy or ten-state yawed trajectory"
    )


def _interpolate_rows(values: np.ndarray, index: int, fraction: float) -> np.ndarray:
    return values[index] + fraction * (values[index + 1] - values[index])


def finite_release_puff_handoff(
    model: Any,
    result: Any,
    *,
    source_duration_s: float,
    source_hydrogen_mass_flow_kg_s: float,
    upstream_material_time_s: float = 0.0,
    hydrogen_flux_relative_tolerance: float = 1.0e-6,
) -> FiniteReleasePuffHandoff:
    """Locate the source-clock transition and create a conserved puff inventory.

    Material time is integrated as ``ds / (|momentum flux| / mass flux)``.
    The switch occurs at the declared source duration, matching the finite
    release clock used by SLAB/SLABx.  The local full-section fluxes are then
    multiplied by that duration.  A mismatch between the declared H2 source
    rate and transported H2 flux is rejected; this catches the carrier-vs-H2
    clock error that otherwise moves the transition distance silently.

    The returned state requires a separate transient puff continuation.  It
    must not be sampled downstream with the steady plume after ``status`` is
    ``"transition_ready"``.
    """

    for name, value in {
        "source_duration_s": source_duration_s,
        "source_hydrogen_mass_flow_kg_s": source_hydrogen_mass_flow_kg_s,
        "hydrogen_flux_relative_tolerance": hydrogen_flux_relative_tolerance,
    }.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(upstream_material_time_s) or upstream_material_time_s < 0.0:
        raise ValueError("upstream_material_time_s must be finite and non-negative")
    if upstream_material_time_s >= source_duration_s:
        raise ValueError(
            "the plume-to-puff transition lies upstream of this trajectory"
        )

    arc = _result_array(result, "arc_length")
    states = _result_array(result, "states")
    fluxes = _result_array(result, "fluxes")
    if arc.ndim != 1 or len(arc) != len(states) or len(arc) < 2:
        raise ValueError("trajectory requires at least two aligned arc stations")
    if np.any(np.diff(arc) <= 0.0):
        raise ValueError("arc length must be strictly increasing")
    positions, momenta, energies, _energy_column = _transition_layout(
        model, states, fluxes
    )
    mass_flux = fluxes[:, 0]
    if np.any(mass_flux <= 0.0):
        raise ValueError("total mass flux must remain positive")
    bulk_speed = np.linalg.norm(momenta, axis=1) / mass_flux
    if np.any(bulk_speed <= 0.0):
        raise ValueError("bulk material speed must remain positive")
    increments = np.diff(arc) * 0.5 * (
        1.0 / bulk_speed[:-1] + 1.0 / bulk_speed[1:]
    )
    material_time = upstream_material_time_s + np.r_[0.0, np.cumsum(increments)]
    target = source_duration_s
    if material_time[-1] < target:
        return FiniteReleasePuffHandoff(
            status="transition_beyond_trajectory",
            source_duration_s=float(source_duration_s),
            transition_time_s=None,
            transition_arc_length_m=None,
            centre_position_m=None,
            transverse_widths_m=None,
            longitudinal_half_length_m=None,
            bulk_velocity_m_s=None,
            total_mass_kg=None,
            hydrogen_mass_kg=None,
            momentum_kg_m_s=None,
            relative_energy_j=None,
            hydrogen_flux_residual_kg_s=None,
            hydrogen_inventory_residual_kg=None,
            requires_transient_puff_continuation=False,
            qualification=(
                "the supplied steady trajectory ends before the finite-release "
                "material clock reaches the source duration"
            ),
        )

    upper = int(np.searchsorted(material_time, target, side="left"))
    if upper == 0:
        lower = 0
        fraction = 0.0
    else:
        lower = upper - 1
        fraction = (
            (target - material_time[lower])
            / (material_time[upper] - material_time[lower])
        )
    state = _interpolate_rows(states, lower, fraction)
    flux = _interpolate_rows(fluxes, lower, fraction)
    position = _interpolate_rows(positions, lower, fraction)
    momentum_flux = _interpolate_rows(momenta, lower, fraction)
    energy_flux = float(
        energies[lower] + fraction * (energies[lower + 1] - energies[lower])
    )
    arc_at_transition = float(
        arc[lower] + fraction * (arc[lower + 1] - arc[lower])
    )
    hydrogen_residual = float(flux[1] - source_hydrogen_mass_flow_kg_s)
    hydrogen_scale = max(abs(source_hydrogen_mass_flow_kg_s), 1.0e-30)
    if abs(hydrogen_residual) / hydrogen_scale > hydrogen_flux_relative_tolerance:
        raise ValueError(
            "transported H2 flux does not match the declared H2 source rate; "
            "do not use carrier-mixture mass for the finite-release clock"
        )

    if states.shape[1] == 7:
        widths = model.section_widths(state)
    else:
        widths = model.base.section_widths(model.proxy(state))
    widths = tuple(float(value) for value in widths)
    if len(widths) != 2 or min(widths) <= 0.0:
        raise ValueError("transition cross-section returned invalid widths")

    total_mass = float(flux[0] * source_duration_s)
    hydrogen_mass = float(source_hydrogen_mass_flow_kg_s * source_duration_s)
    momentum = momentum_flux * source_duration_s
    velocity = momentum_flux / flux[0]
    speed = float(np.linalg.norm(velocity))
    return FiniteReleasePuffHandoff(
        status="transition_ready",
        source_duration_s=float(source_duration_s),
        transition_time_s=float(target),
        transition_arc_length_m=arc_at_transition,
        centre_position_m=tuple(float(value) for value in position),
        transverse_widths_m=widths,
        longitudinal_half_length_m=float(0.5 * speed * source_duration_s),
        bulk_velocity_m_s=tuple(float(value) for value in velocity),
        total_mass_kg=total_mass,
        hydrogen_mass_kg=hydrogen_mass,
        momentum_kg_m_s=tuple(float(value) for value in momentum),
        relative_energy_j=float(energy_flux * source_duration_s),
        hydrogen_flux_residual_kg_s=hydrogen_residual,
        hydrogen_inventory_residual_kg=float(
            flux[1] * source_duration_s - hydrogen_mass
        ),
        requires_transient_puff_continuation=True,
        qualification=(
            "conservative full-cross-section plume-to-puff handoff; downstream "
            "puff entrainment, gravity spreading and thermodynamics are not "
            "integrated by this interface"
        ),
    )


def assess_steady_wind_applicability(
    history: WindHistory,
    *,
    start_s: float | None = None,
    end_s: float | None = None,
    maximum_direction_span_deg: float = 20.0,
    maximum_speed_range_fraction: float = 0.25,
) -> SteadyWindApplicability:
    """Reject one-steady-plume use when a declared wind window is too variable.

    Thresholds are applicability choices, not fitted physics.  Direction is
    unwrapped before its range is taken, so a 359-to-1 degree record has a
    two-degree rather than a 358-degree span.
    """

    if maximum_direction_span_deg <= 0.0 or maximum_speed_range_fraction <= 0.0:
        raise ValueError("steady-wind thresholds must be positive")
    time, speed, direction = history.arrays()
    start = float(time[0] if start_s is None else start_s)
    end = float(time[-1] if end_s is None else end_s)
    if not math.isfinite(start) or not math.isfinite(end) or end <= start:
        raise ValueError("wind-assessment window must be finite and increasing")
    if start < time[0] or end > time[-1]:
        raise ValueError("wind-assessment window lies outside the measured record")
    selected = (time >= start) & (time <= end)
    sample_time = np.r_[start, time[selected], end]
    sample_time = np.unique(sample_time)
    sample_speed = np.interp(sample_time, time, speed)
    angle = np.unwrap(np.deg2rad(direction))
    sample_angle = np.interp(sample_time, time, angle)
    mean_speed = float(np.mean(sample_speed))
    speed_fraction = float(np.ptp(sample_speed) / mean_speed)
    direction_span = float(np.rad2deg(np.ptp(sample_angle)))
    reasons: list[str] = []
    if speed_fraction > maximum_speed_range_fraction:
        reasons.append("wind_speed_range_exceeds_steady_limit")
    if direction_span > maximum_direction_span_deg:
        reasons.append("wind_direction_span_exceeds_steady_limit")
    if not reasons:
        reasons.append("wind_window_within_declared_steady_limits")
    return SteadyWindApplicability(
        applicable=len(reasons) == 1 and reasons[0].startswith("wind_window_within"),
        start_s=start,
        end_s=end,
        mean_speed_m_s=mean_speed,
        speed_range_fraction=speed_fraction,
        direction_span_deg=direction_span,
        maximum_speed_range_fraction=float(maximum_speed_range_fraction),
        maximum_direction_span_deg=float(maximum_direction_span_deg),
        reasons=tuple(reasons),
    )


__all__ = [
    "FiniteReleasePuffHandoff", "SteadyWindApplicability",
    "finite_release_puff_handoff", "assess_steady_wind_applicability",
]
