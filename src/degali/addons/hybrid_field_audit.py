"""Audits for a coarse/local result pair before field stitching.

The audit is deliberately fail-closed: it does not average incompatible
fields.  It only accepts results with a common time axis and explicit mass
arrays, then reports conservation and temporal alignment evidence.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

def threshold_region_metrics(field_vol_percent: np.ndarray, *, threshold_vol_percent: float = 4.0,
                             cell_volume_m3: float, coordinates_m: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None) -> dict[str, object]:
    """Return reproducible 4-vol% region metrics for a scalar field.

    ``field_vol_percent`` is indexed ``(z, y, x)`` and may contain NaNs on
    solid cells.  The function deliberately reports geometry only; it does
    not infer a safety classification or replace a time-history exceedance
    analysis.
    """
    field = np.asarray(field_vol_percent, dtype=float)
    volume = float(cell_volume_m3)
    threshold = float(threshold_vol_percent)
    if field.ndim != 3 or not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("field must be 3-D and cell_volume_m3 must be positive")
    if not np.isfinite(threshold) or threshold < 0.0:
        raise ValueError("threshold_vol_percent must be finite and non-negative")
    mask = np.isfinite(field) & (field >= threshold)
    values = field[mask]
    result: dict[str, object] = {
        "threshold_vol_percent": threshold,
        "cell_count": int(np.count_nonzero(mask)),
        "volume_m3": float(np.count_nonzero(mask) * volume),
        "maximum_vol_percent": (None if values.size == 0 else float(np.max(values))),
        "centroid_m": None,
    }
    if coordinates_m is not None:
        if len(coordinates_m) != 3:
            raise ValueError("coordinates_m must contain z, y, and x arrays")
        z, y, x = (np.asarray(axis, dtype=float) for axis in coordinates_m)
        if any(axis.ndim != 1 for axis in (z, y, x)) or field.shape != (len(z), len(y), len(x)):
            raise ValueError("coordinate axes do not match field shape")
        if values.size:
            kk, jj, ii = np.nonzero(mask)
            result["centroid_m"] = [
                float(np.mean(z[kk])), float(np.mean(y[jj])), float(np.mean(x[ii]))
            ]
    return result

def infer_net_internal_flux_rate(result: object, *, time_tolerance_s: float = 1e-12) -> np.ndarray:
    """Infer net internal mass flux required by a solver's balance history.

    Positive values mean the domain needs net mass entering through an internal
    interface.  This is a diagnostic closure signal, not a replacement for
    face-resolved fluxes; it is useful for constructing the next coupled run.
    """
    required = ("time_s", "domain_mass_kg", "injected_mass_kg", "outflow_mass_kg")
    missing = [name for name in required if not hasattr(result, name)]
    if missing:
        raise TypeError(f"result missing balance histories: {missing}")
    t = np.asarray(result.time_s, dtype=float)
    mass = np.asarray(result.domain_mass_kg, dtype=float)
    injected = np.asarray(result.injected_mass_kg, dtype=float)
    outflow = np.asarray(result.outflow_mass_kg, dtype=float)
    if any(a.ndim != 1 or a.shape != t.shape for a in (mass, injected, outflow)):
        raise ValueError("solver balance histories must have matching one-dimensional shapes")
    if len(t) < 2 or np.any(~np.isfinite(np.r_[t, mass, injected, outflow])):
        raise ValueError("solver balance histories must be finite and non-empty")
    dt = np.diff(t)
    if np.any(dt <= float(time_tolerance_s)):
        raise ValueError("solver time history must be strictly increasing")
    return np.diff(mass) / dt - (np.diff(injected) / dt - np.diff(outflow) / dt)

def aggregate_substeps_to_coarse(*, fine_time_s, fine_mass_kg, coarse_time_s,
                                 time_tolerance_s: float = 1e-10):
    """Sample fine-grid endpoint states at an exact coarse clock.

    No interpolation or mass creation is performed.  Every coarse timestamp
    must exist in the fine history; otherwise the operation fails closed.
    """
    fine_t = np.asarray(fine_time_s, dtype=float)
    fine_m = np.asarray(fine_mass_kg, dtype=float)
    coarse_t = np.asarray(coarse_time_s, dtype=float)
    if fine_t.ndim != 1 or fine_m.shape != fine_t.shape or coarse_t.ndim != 1:
        raise ValueError("fine time/mass and coarse time histories must be one-dimensional")
    if len(fine_t) == 0 or not np.all(np.isfinite(np.r_[fine_t, fine_m, coarse_t])):
        raise ValueError("histories must be non-empty and finite")
    if np.any(np.diff(fine_t) <= 0) or np.any(np.diff(coarse_t) <= 0):
        raise ValueError("time histories must be strictly increasing")
    indices = []
    for target in coarse_t:
        matches = np.flatnonzero(np.isclose(fine_t, target, rtol=0, atol=time_tolerance_s))
        if len(matches) == 0:
            raise ValueError("coarse timestamp is not an exact fine-grid endpoint")
        indices.append(int(matches[-1]))
    return fine_m[np.asarray(indices, dtype=int)].copy()

@dataclass(frozen=True)
class HybridFieldAudit:
    common_time_s: np.ndarray
    coarse_mass_kg: np.ndarray
    patch_mass_kg: np.ndarray
    total_mass_kg: np.ndarray
    initial_total_mass_kg: float
    maximum_total_mass_drift_kg: float
    time_aligned: bool
    quantitative_ready: bool

def audit_hybrid_mass_histories(*, coarse_time_s, coarse_mass_kg,
                                patch_time_s, patch_mass_kg,
                                time_tolerance_s: float = 1e-10,
                                mass_tolerance_kg: float = 1e-10) -> HybridFieldAudit:
    """Combine two *non-overlapping* mass histories on an identical clock."""
    t0 = np.asarray(coarse_time_s, dtype=float)
    t1 = np.asarray(patch_time_s, dtype=float)
    m0 = np.asarray(coarse_mass_kg, dtype=float)
    m1 = np.asarray(patch_mass_kg, dtype=float)
    if t0.ndim != 1 or t1.ndim != 1 or m0.shape != t0.shape or m1.shape != t1.shape:
        raise ValueError("time and mass histories must be one-dimensional and paired")
    if len(t0) == 0 or not np.all(np.isfinite(np.r_[t0, t1, m0, m1])):
        raise ValueError("histories must be non-empty and finite")
    if np.any(m0 < 0) or np.any(m1 < 0):
        raise ValueError("mass histories must be non-negative")
    if t0.shape != t1.shape or not np.allclose(t0, t1, rtol=0, atol=time_tolerance_s):
        raise ValueError("coarse and patch histories require a common time axis")
    total = m0 + m1
    drift = total - total[0]
    maximum = float(np.max(np.abs(drift)))
    return HybridFieldAudit(
        common_time_s=t0.copy(), coarse_mass_kg=m0.copy(), patch_mass_kg=m1.copy(),
        total_mass_kg=total, initial_total_mass_kg=float(total[0]),
        maximum_total_mass_drift_kg=maximum, time_aligned=True,
        quantitative_ready=bool(maximum <= float(mass_tolerance_kg)),
    )

def audit_solver_results(coarse_result: object, patch_result: object, *,
                         mass_tolerance_kg: float = 1e-10) -> HybridFieldAudit:
    """Audit two solver result objects using their public histories.

    Both DEGALI obstacle result types expose ``time_s`` and
    ``domain_mass_kg``.  Missing or malformed histories are rejected rather
    than silently reconstructed from concentration fields.
    """
    required = ("time_s", "domain_mass_kg")
    for name, result in (("coarse", coarse_result), ("patch", patch_result)):
        missing = [attribute for attribute in required if not hasattr(result, attribute)]
        if missing:
            raise TypeError(f"{name}_result missing required histories: {missing}")
    return audit_hybrid_mass_histories(
        coarse_time_s=getattr(coarse_result, "time_s"),
        coarse_mass_kg=getattr(coarse_result, "domain_mass_kg"),
        patch_time_s=getattr(patch_result, "time_s"),
        patch_mass_kg=getattr(patch_result, "domain_mass_kg"),
        mass_tolerance_kg=mass_tolerance_kg,
    )

def audit_solver_results_with_patch_substeps(coarse_result: object,
                                             patch_result: object, *,
                                             mass_tolerance_kg: float = 1e-10,
                                             time_tolerance_s: float = 1e-10) -> HybridFieldAudit:
    """Audit solver results when patch output contains finer time steps."""
    for name, result in (("coarse", coarse_result), ("patch", patch_result)):
        if not all(hasattr(result, field) for field in ("time_s", "domain_mass_kg")):
            raise TypeError(f"{name}_result must expose time_s and domain_mass_kg")
    coarse_t = np.asarray(coarse_result.time_s, dtype=float)
    patch_mass = aggregate_substeps_to_coarse(
        fine_time_s=patch_result.time_s,
        fine_mass_kg=patch_result.domain_mass_kg,
        coarse_time_s=coarse_t,
        time_tolerance_s=time_tolerance_s,
    )
    return audit_hybrid_mass_histories(
        coarse_time_s=coarse_t,
        coarse_mass_kg=coarse_result.domain_mass_kg,
        patch_time_s=coarse_t,
        patch_mass_kg=patch_mass,
        time_tolerance_s=time_tolerance_s,
        mass_tolerance_kg=mass_tolerance_kg,
    )

def audit_partitioned_solver_results(coarse_result: object, patch_result: object,
                                     *, coarse_weight: float, patch_weight: float,
                                     mass_tolerance_kg: float = 1e-10,
                                     time_tolerance_s: float = 1e-10) -> HybridFieldAudit:
    """Audit results when coarse/patch masses are partitioned, not additive.

    The weights represent non-overlapping control-volume fractions after an
    explicit overlap remap.  This prevents double-counting the patch mass that
    caused the independent-run diagnostic to drift.
    """
    wc, wp = float(coarse_weight), float(patch_weight)
    if not np.isfinite(wc) or not np.isfinite(wp) or wc < 0 or wp < 0:
        raise ValueError("partition weights must be finite and non-negative")
    if abs(wc + wp - 1.0) > 1e-12:
        raise ValueError("partition weights must sum to one")
    for name, result in (("coarse", coarse_result), ("patch", patch_result)):
        if not all(hasattr(result, field) for field in ("time_s", "domain_mass_kg")):
            raise TypeError(f"{name}_result must expose time_s and domain_mass_kg")
    coarse_t = np.asarray(coarse_result.time_s, dtype=float)
    patch_m = aggregate_substeps_to_coarse(
        fine_time_s=patch_result.time_s, fine_mass_kg=patch_result.domain_mass_kg,
        coarse_time_s=coarse_t, time_tolerance_s=time_tolerance_s)
    coarse_m = np.asarray(coarse_result.domain_mass_kg, dtype=float)
    return audit_hybrid_mass_histories(
        coarse_time_s=coarse_t, coarse_mass_kg=wc * coarse_m,
        patch_time_s=coarse_t, patch_mass_kg=wp * patch_m,
        time_tolerance_s=time_tolerance_s, mass_tolerance_kg=mass_tolerance_kg)

def audit_partitioned_solver_results_with_plan(coarse_result: object,
                                                patch_result: object, plan: object,
                                                *, mass_tolerance_kg: float = 1e-10,
                                                time_tolerance_s: float = 1e-10) -> HybridFieldAudit:
    """Use geometric partition weights recorded by ``HybridAMRRunPlan``."""
    if not hasattr(plan, "partition_weights"):
        raise TypeError("plan must expose partition_weights")
    coarse_weight, patch_weight = plan.partition_weights
    return audit_partitioned_solver_results(
        coarse_result, patch_result, coarse_weight=coarse_weight,
        patch_weight=patch_weight, mass_tolerance_kg=mass_tolerance_kg,
        time_tolerance_s=time_tolerance_s)

__all__ = ["HybridFieldAudit", "aggregate_substeps_to_coarse",
           "threshold_region_metrics",
           "infer_net_internal_flux_rate",
           "audit_hybrid_mass_histories", "audit_solver_results",
           "audit_solver_results_with_patch_substeps",
           "audit_partitioned_solver_results",
           "audit_partitioned_solver_results_with_plan"]
