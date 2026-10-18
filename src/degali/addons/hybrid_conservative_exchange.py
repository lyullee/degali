"""Conservative scalar/energy/momentum exchange for a coarse/patch pair.

This module is the narrow interface between the fast outer transport and a
resolved obstacle patch.  It does not claim that either solver is accurate;
it guarantees that a coupled execution cannot silently create H2, sensible
enthalpy, or momentum at the interface.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


def _finite(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _momentum(value: Iterable[float]) -> tuple[float, float, float]:
    array = np.asarray(tuple(value), dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError("momentum must contain three finite components")
    return tuple(float(item) for item in array)


@dataclass(frozen=True)
class ConservativeInterfacePacket:
    """One directed packet crossing a coarse/patch interface."""

    total_mass_kg: float
    species_mass_kg: float
    sensible_enthalpy_j: float = 0.0
    momentum_kg_m_s: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self) -> None:
        total = _finite(self.total_mass_kg, "total_mass_kg")
        species = _finite(self.species_mass_kg, "species_mass_kg")
        _finite(self.sensible_enthalpy_j, "sensible_enthalpy_j")
        _momentum(self.momentum_kg_m_s)
        if total < 0.0 or species < 0.0 or species > total + 1.0e-12:
            raise ValueError("packet species mass must lie between zero and total mass")

    @classmethod
    def zero(cls) -> "ConservativeInterfacePacket":
        return cls(0.0, 0.0)


@dataclass(frozen=True)
class ConservativeControlVolumeState:
    """Aggregated state used by the interface ledger.

    The production solvers may map these quantities to boundary faces or
    control-volume strips.  The exchange kernel intentionally operates on
    integrated quantities so the equal-and-opposite update is unambiguous.
    """

    total_mass_kg: float
    species_mass_kg: float
    sensible_enthalpy_j: float = 0.0
    momentum_kg_m_s: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self) -> None:
        total = _finite(self.total_mass_kg, "total_mass_kg")
        species = _finite(self.species_mass_kg, "species_mass_kg")
        _finite(self.sensible_enthalpy_j, "sensible_enthalpy_j")
        _momentum(self.momentum_kg_m_s)
        if total < 0.0 or species < 0.0 or species > total + 1.0e-12:
            raise ValueError("state species mass must lie between zero and total mass")

    def receive(self, packet: ConservativeInterfacePacket) -> "ConservativeControlVolumeState":
        return ConservativeControlVolumeState(
            self.total_mass_kg + packet.total_mass_kg,
            self.species_mass_kg + packet.species_mass_kg,
            self.sensible_enthalpy_j + packet.sensible_enthalpy_j,
            tuple(
                a + b for a, b in zip(self.momentum_kg_m_s, packet.momentum_kg_m_s)
            ),
        )

    def send(self, packet: ConservativeInterfacePacket) -> "ConservativeControlVolumeState":
        tolerance = 1.0e-12
        if packet.total_mass_kg > self.total_mass_kg + tolerance:
            raise ValueError("outgoing packet exceeds donor total mass")
        if packet.species_mass_kg > self.species_mass_kg + tolerance:
            raise ValueError("outgoing packet exceeds donor species mass")
        return ConservativeControlVolumeState(
            max(0.0, self.total_mass_kg - packet.total_mass_kg),
            max(0.0, self.species_mass_kg - packet.species_mass_kg),
            self.sensible_enthalpy_j - packet.sensible_enthalpy_j,
            tuple(
                a - b for a, b in zip(self.momentum_kg_m_s, packet.momentum_kg_m_s)
            ),
        )


def exchange_conservative_pair(
    coarse: ConservativeControlVolumeState,
    patch: ConservativeControlVolumeState,
    *,
    coarse_to_patch: ConservativeInterfacePacket,
    patch_to_coarse: ConservativeInterfacePacket,
) -> tuple[ConservativeControlVolumeState, ConservativeControlVolumeState]:
    """Apply equal-and-opposite packets and return updated states."""

    coarse_after_send = coarse.send(coarse_to_patch)
    patch_after_send = patch.send(patch_to_coarse)
    return (
        coarse_after_send.receive(patch_to_coarse),
        patch_after_send.receive(coarse_to_patch),
    )


def conserved_totals(
    coarse: ConservativeControlVolumeState,
    patch: ConservativeControlVolumeState,
) -> ConservativeInterfacePacket:
    """Return the pair's conserved total for a post-run audit."""

    return ConservativeInterfacePacket(
        coarse.total_mass_kg + patch.total_mass_kg,
        coarse.species_mass_kg + patch.species_mass_kg,
        coarse.sensible_enthalpy_j + patch.sensible_enthalpy_j,
        tuple(
            a + b for a, b in zip(coarse.momentum_kg_m_s, patch.momentum_kg_m_s)
        ),
    )


__all__ = [
    "ConservativeInterfacePacket",
    "ConservativeControlVolumeState",
    "exchange_conservative_pair",
    "conserved_totals",
]
