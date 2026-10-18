"""Execution-plan metadata for coarse transport plus a local obstacle patch."""

from __future__ import annotations

from dataclasses import dataclass
from .hybrid_amr import ObstaclePatchBounds


@dataclass(frozen=True)
class HybridAMRRunPlan:
    """A reproducible coarse/local-patch execution plan.

    This is deliberately a plan object: the solver is required to consume it
    and record the same geometry and clocks in its output before promotion.
    """

    domain_m: tuple[float, float, float]
    coarse_cell_size_m: float
    patch_cell_size_m: float
    patch: ObstaclePatchBounds
    duration_s: float
    coarse_time_step_s: float
    patch_time_step_s: float

    @property
    def patch_volume_fraction(self) -> float:
        domain_volume = self.domain_m[0] * self.domain_m[1] * self.domain_m[2]
        return self.patch.volume_m3 / domain_volume

    @property
    def substeps_per_coarse_step(self) -> int:
        ratio = self.coarse_time_step_s / self.patch_time_step_s
        rounded = round(ratio)
        if rounded < 1 or abs(ratio - rounded) > 1.0e-12:
            raise ValueError("patch time step must divide coarse time step")
        return rounded

    @property
    def partition_weights(self) -> tuple[float, float]:
        """Return ``(coarse, patch)`` geometric control-volume weights."""
        fraction = self.patch_volume_fraction
        if not 0.0 < fraction < 1.0:
            raise ValueError("patch must occupy a strict subset of the domain")
        return 1.0 - fraction, fraction


__all__ = ["HybridAMRRunPlan"]
