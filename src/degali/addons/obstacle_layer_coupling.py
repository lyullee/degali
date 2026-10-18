"""Conservative height-resolved obstacle-to-layer interface.

The resolved obstacle solver and the dense-layer solver use different state
representations.  This module is deliberately narrow: a boundary crossing is
exported as a time packet containing height-resolved total/species/momentum
rates, and can be collapsed once (without renormalisation) for a layer
boundary condition.  It does not claim that the two solvers have the same
transport physics; it prevents the interface from silently creating or
destroying inventory.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

import numpy as np

from .hybrid_conservative_exchange import ConservativeInterfacePacket


_FACES = {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"}


def _as_rate(value: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite 2-D rate array")
    if np.any(array < 0.0):
        raise ValueError(f"{name} must be non-negative")
    return array.copy()


@dataclass(frozen=True)
class HeightResolvedBoundaryPacket:
    """One obstacle-face packet before vertical collapse.

    Arrays are indexed ``(height, transverse)`` and contain outward positive
    rates in kg/s or kg m/s.  Momentum is stored as three matching arrays.
    ``face`` is retained so a receiver can choose the correct normal boundary
    rather than treating all exported mass as an x-outlet.
    """

    time_s: float
    dt_s: float
    face: str
    total_rate_kg_s: np.ndarray
    species_rate_kg_s: np.ndarray
    momentum_rate_kg_m_s: np.ndarray

    def __post_init__(self) -> None:
        time_s, dt_s = float(self.time_s), float(self.dt_s)
        if not math.isfinite(time_s) or not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("packet time and positive dt must be finite")
        if self.face not in _FACES:
            raise ValueError(f"unsupported boundary face: {self.face!r}")
        total = _as_rate(self.total_rate_kg_s, "total_rate_kg_s")
        species = _as_rate(self.species_rate_kg_s, "species_rate_kg_s")
        momentum = np.asarray(self.momentum_rate_kg_m_s, dtype=float)
        if momentum.shape != (3, *total.shape) or not np.all(np.isfinite(momentum)):
            raise ValueError("momentum_rate_kg_m_s must have shape (3, nz, ntransverse)")
        if np.any(species > total + 1.0e-12 * np.maximum(total, 1.0)):
            raise ValueError("species rate cannot exceed total rate")
        object.__setattr__(self, "time_s", time_s)
        object.__setattr__(self, "dt_s", dt_s)
        object.__setattr__(self, "total_rate_kg_s", total)
        object.__setattr__(self, "species_rate_kg_s", species)
        object.__setattr__(self, "momentum_rate_kg_m_s", momentum.copy())

    @property
    def total_mass_kg(self) -> float:
        return float(self.total_rate_kg_s.sum() * self.dt_s)

    @property
    def species_mass_kg(self) -> float:
        return float(self.species_rate_kg_s.sum() * self.dt_s)

    @property
    def momentum_kg_m_s(self) -> tuple[float, float, float]:
        return tuple(float(value * self.dt_s) for value in self.momentum_rate_kg_m_s.sum(axis=(1, 2)))

    def collapse_to_layer_packet(self) -> ConservativeInterfacePacket:
        """Collapse height/transverse bins exactly into one layer packet."""

        return ConservativeInterfacePacket(
            total_mass_kg=self.total_mass_kg,
            species_mass_kg=self.species_mass_kg,
            momentum_kg_m_s=self.momentum_kg_m_s,
        )


def export_height_resolved_packets(
    *,
    time_s: Sequence[float],
    dt_s: float,
    face: str,
    total_rate_kg_s: np.ndarray,
    species_rate_kg_s: np.ndarray,
    momentum_rate_kg_m_s: np.ndarray | None = None,
) -> tuple[HeightResolvedBoundaryPacket, ...]:
    """Export one packet per solver interval from a boundary-rate history.

    ``time_s`` has length ``n_steps`` and all rate histories have shape
    ``(n_steps, nz, ntransverse)``.  If no total or momentum history is
    available, callers must explicitly pass species as total; this prevents an
    accidental interpretation of a passive-species flux as a mixture flux.
    """

    times = np.asarray(time_s, dtype=float)
    total = np.asarray(total_rate_kg_s, dtype=float)
    species = np.asarray(species_rate_kg_s, dtype=float)
    if times.ndim != 1 or not np.all(np.isfinite(times)):
        raise ValueError("time_s must be a finite 1-D array")
    if total.ndim != 3 or species.shape != total.shape or total.shape[0] != len(times):
        raise ValueError("rate histories must have shape (n_steps, nz, ntransverse)")
    if momentum_rate_kg_m_s is None:
        momentum = np.zeros((3, *total.shape), dtype=float)
    else:
        momentum = np.asarray(momentum_rate_kg_m_s, dtype=float)
        if momentum.shape != (3, *total.shape):
            raise ValueError("momentum history must have shape (3, n_steps, nz, ntransverse)")
    return tuple(
        HeightResolvedBoundaryPacket(
            time_s=float(times[index]), dt_s=dt_s, face=face,
            total_rate_kg_s=total[index], species_rate_kg_s=species[index],
            momentum_rate_kg_m_s=momentum[:, index],
        )
        for index in range(len(times))
    )


def collapse_packet_series(
    packets: Iterable[HeightResolvedBoundaryPacket],
) -> ConservativeInterfacePacket:
    """Return the exact integrated layer packet for a packet sequence."""

    total = species = enthalpy = 0.0
    momentum = np.zeros(3, dtype=float)
    for packet in packets:
        total += packet.total_mass_kg
        species += packet.species_mass_kg
        momentum += np.asarray(packet.momentum_kg_m_s)
    return ConservativeInterfacePacket(total, species, enthalpy, tuple(momentum))


__all__ = [
    "HeightResolvedBoundaryPacket",
    "export_height_resolved_packets",
    "collapse_packet_series",
]
