"""Conservative coupling primitives for coarse fields and obstacle patches.

The coupling contract is expressed in species mass, rather than concentration
values.  This keeps a local obstacle patch from creating or destroying H2 when
its mesh or time step differs from the surrounding transport model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math


@dataclass(frozen=True)
class InterfaceFluxLedger:
    """Accumulated H2 mass crossing one coarse/patch interface."""

    incoming_kg: float = 0.0
    outgoing_kg: float = 0.0

    def __post_init__(self) -> None:
        for name in ("incoming_kg", "outgoing_kg"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")

    @property
    def residual_kg(self) -> float:
        return self.incoming_kg - self.outgoing_kg

    def add(self, *, incoming_kg: float = 0.0, outgoing_kg: float = 0.0) -> "InterfaceFluxLedger":
        incoming = float(incoming_kg)
        outgoing = float(outgoing_kg)
        if not math.isfinite(incoming) or not math.isfinite(outgoing):
            raise ValueError("interface flux must be finite")
        if incoming < 0.0 or outgoing < 0.0:
            raise ValueError("interface flux must be non-negative")
        return InterfaceFluxLedger(
            self.incoming_kg + incoming,
            self.outgoing_kg + outgoing,
        )


@dataclass
class HybridInterfaceCoupler:
    """Stateful conservative coarse/patch exchange for one time step stream."""

    ledger: InterfaceFluxLedger = field(default_factory=InterfaceFluxLedger)

    def transfer_to_patch(self, donor_mass_kg: float, requested_mass_kg: float) -> float:
        transferred, _ = conservative_interface_transfer(
            donor_mass_kg=donor_mass_kg,
            requested_mass_kg=requested_mass_kg,
        )
        self.ledger = self.ledger.add(outgoing_kg=transferred)
        return transferred

    def transfer_to_coarse(self, donor_mass_kg: float, requested_mass_kg: float) -> float:
        transferred, _ = conservative_interface_transfer(
            donor_mass_kg=donor_mass_kg,
            requested_mass_kg=requested_mass_kg,
        )
        self.ledger = self.ledger.add(incoming_kg=transferred)
        return transferred

    def require_balanced(self, tolerance_kg: float = 1.0e-12) -> None:
        tolerance = float(tolerance_kg)
        if not math.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("tolerance_kg must be finite and non-negative")
        if abs(self.ledger.residual_kg) > tolerance:
            raise RuntimeError(
                "hybrid interface mass residual exceeds tolerance: "
                f"{self.ledger.residual_kg:.6g} kg"
            )

    def exchange_step(
        self,
        *,
        coarse_mass_kg: float,
        patch_mass_kg: float,
        coarse_to_patch_kg: float = 0.0,
        patch_to_coarse_kg: float = 0.0,
    ) -> tuple[float, float]:
        """Apply one conservative exchange and return updated domain masses."""
        to_patch = self.transfer_to_patch(coarse_mass_kg, coarse_to_patch_kg)
        to_coarse = self.transfer_to_coarse(patch_mass_kg, patch_to_coarse_kg)
        return coarse_mass_kg - to_patch + to_coarse, patch_mass_kg + to_patch - to_coarse


def conservative_interface_transfer(
    *,
    donor_mass_kg: float,
    requested_mass_kg: float,
) -> tuple[float, float]:
    """Return ``(transferred, donor_remainder)`` without creating mass.

    A requested patch transfer is clipped to the available donor mass.  The
    caller applies the same transferred amount with opposite signs to the two
    domains and records it in :class:`InterfaceFluxLedger`.
    """
    donor = float(donor_mass_kg)
    requested = float(requested_mass_kg)
    if not math.isfinite(donor) or not math.isfinite(requested):
        raise ValueError("mass values must be finite")
    if donor < 0.0 or requested < 0.0:
        raise ValueError("mass values must be non-negative")
    transferred = min(donor, requested)
    return transferred, donor - transferred


__all__ = [
    "InterfaceFluxLedger", "HybridInterfaceCoupler",
    "conservative_interface_transfer",
]
