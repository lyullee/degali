"""Shared-clock conservative executor for the coarse/local prototype."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
from .hybrid_amr_plan import HybridAMRRunPlan
from .hybrid_interface import HybridInterfaceCoupler
from .hybrid_amr import PatchFluxAccumulator, cartesian_face_coordinates
from .hybrid_conservative_exchange import (
    ConservativeControlVolumeState,
    ConservativeInterfacePacket,
    exchange_conservative_pair,
)

MassRequest = Callable[[float, float, float], tuple[float, float]]

@dataclass(frozen=True)
class HybridAMRExecutionResult:
    """Mass-only execution evidence; not a field-validation result."""
    coarse_mass_kg: float
    patch_mass_kg: float
    elapsed_s: float
    coarse_steps: int
    patch_steps: int
    interface_residual_kg: float

class HybridAMRExecutor:
    """Run a shared coarse/patch clock with conservative scalar exchange."""
    def __init__(self, plan: HybridAMRRunPlan, *, tolerance_kg: float = 1e-12):
        self.plan = plan
        self.tolerance_kg = float(tolerance_kg)
        self.coupler = HybridInterfaceCoupler()

    def run(self, *, coarse_mass_kg: float, patch_mass_kg: float,
            request: MassRequest) -> HybridAMRExecutionResult:
        coarse, patch = float(coarse_mass_kg), float(patch_mass_kg)
        if coarse < 0 or patch < 0:
            raise ValueError("initial masses must be non-negative")
        initial_total = coarse + patch
        substeps = self.plan.substeps_per_coarse_step
        n_coarse = round(self.plan.duration_s / self.plan.coarse_time_step_s)
        if abs(n_coarse * self.plan.coarse_time_step_s - self.plan.duration_s) > 1e-12:
            raise ValueError("duration must be an integer number of coarse steps")
        n_patch = n_coarse * substeps
        dt = self.plan.patch_time_step_s
        for index in range(n_patch):
            requested = request(index * dt, coarse, patch)
            if len(requested) != 2:
                raise ValueError("request must return two mass transfers")
            coarse, patch = self.coupler.exchange_step(
                coarse_mass_kg=coarse, patch_mass_kg=patch,
                coarse_to_patch_kg=requested[0], patch_to_coarse_kg=requested[1])
        conservation_residual = (coarse + patch) - initial_total
        if abs(conservation_residual) > self.tolerance_kg:
            raise RuntimeError("hybrid execution mass conservation residual exceeds tolerance")
        return HybridAMRExecutionResult(
            coarse_mass_kg=coarse, patch_mass_kg=patch,
            elapsed_s=self.plan.duration_s, coarse_steps=n_coarse,
            patch_steps=n_patch,
            interface_residual_kg=conservation_residual)

    def run_with_flux_rates(self, *, coarse_mass_kg: float, patch_mass_kg: float,
                            time_s, coarse_to_patch_kg_s, patch_to_coarse_kg_s,
                            time_tolerance_s: float = 1e-10) -> HybridAMRExecutionResult:
        """Run using measured/calculated interface mass-flux rates.

        Rates are sampled only at exact patch-clock endpoints.  The requested
        mass for each substep is rate multiplied by the patch timestep; no
        interpolation is silently introduced.
        """
        import numpy as np
        times = np.asarray(time_s, dtype=float)
        forward = np.asarray(coarse_to_patch_kg_s, dtype=float)
        reverse = np.asarray(patch_to_coarse_kg_s, dtype=float)
        if times.ndim != 1 or forward.shape != times.shape or reverse.shape != times.shape:
            raise ValueError("flux rates and time_s must have matching one-dimensional shapes")
        if len(times) == 0 or np.any(~np.isfinite(np.r_[times, forward, reverse])):
            raise ValueError("flux schedule must be finite and non-empty")
        if np.any(forward < 0) or np.any(reverse < 0) or np.any(np.diff(times) <= 0):
            raise ValueError("flux rates must be non-negative and time must increase")
        subdt = self.plan.patch_time_step_s
        def request(elapsed, _coarse, _patch):
            matches = np.flatnonzero(np.isclose(times, elapsed, rtol=0, atol=time_tolerance_s))
            if len(matches) == 0:
                raise ValueError("flux schedule lacks an exact patch-clock endpoint")
            index = int(matches[-1])
            return float(forward[index] * subdt), float(reverse[index] * subdt)
        return self.run(coarse_mass_kg=coarse_mass_kg, patch_mass_kg=patch_mass_kg,
                        request=request)

    def run_bidirectional(self, *, coarse_state, patch_state,
                          coarse_step, patch_step,
                          coarse_mass_kg: float = 0.0,
                          patch_mass_kg: float = 0.0):
        """Advance two state callbacks with equal-and-opposite interface flux.

        Each callback receives ``(state, elapsed_s, dt_s, incoming_mass_kg)``
        and returns ``(new_state, outgoing_mass_kg)``.  The returned outgoing
        mass is limited by the donor inventory and applied to the other state
        on the same substep.  This is the orchestration contract for plugging
        real field solvers into the conservative exchange kernel.
        """
        import math
        if not callable(coarse_step) or not callable(patch_step):
            raise TypeError("coarse_step and patch_step must be callable")
        n = round(self.plan.duration_s / self.plan.patch_time_step_s)
        if n < 1 or abs(n * self.plan.patch_time_step_s - self.plan.duration_s) > 1e-12:
            raise ValueError("duration must be an integer number of patch steps")
        coarse_mass = float(coarse_mass_kg)
        patch_mass = float(patch_mass_kg)
        if coarse_mass < 0.0 or patch_mass < 0.0:
            raise ValueError("initial masses must be non-negative")
        t = 0.0
        history = []
        incoming_coarse = incoming_patch = 0.0
        for _ in range(n):
            coarse_state, coarse_out = coarse_step(coarse_state, t, self.plan.patch_time_step_s, incoming_coarse)
            patch_state, patch_out = patch_step(patch_state, t, self.plan.patch_time_step_s, incoming_patch)
            coarse_out, patch_out = float(coarse_out), float(patch_out)
            if not math.isfinite(coarse_out) or not math.isfinite(patch_out) or coarse_out < 0 or patch_out < 0:
                raise ValueError("step callbacks must return finite non-negative outgoing masses")
            coarse_to_patch = min(coarse_mass, coarse_out)
            patch_to_coarse = min(patch_mass, patch_out)
            coarse_mass += patch_to_coarse - coarse_to_patch
            patch_mass += coarse_to_patch - patch_to_coarse
            incoming_coarse, incoming_patch = patch_to_coarse, coarse_to_patch
            history.append((t, coarse_to_patch, patch_to_coarse, coarse_mass, patch_mass))
            t += self.plan.patch_time_step_s
        return coarse_state, patch_state, tuple(history)

    def run_bidirectional_transactional(self, *, coarse_state, patch_state,
                                        coarse_step, patch_step,
                                        coarse_mass_kg: float = 0.0,
                                        patch_mass_kg: float = 0.0):
        """Two-phase-safe variant passing donor availability to each callback.

        Callbacks receive ``(state, t, dt, incoming_mass, available_mass)`` and
        must return ``(new_state, accepted_outgoing_mass)``.  Because the
        available amount is supplied before state mutation, field updates and
        executor bookkeeping cannot diverge through post-hoc clipping.
        """
        import math
        cm, pm = float(coarse_mass_kg), float(patch_mass_kg)
        if cm < 0 or pm < 0:
            raise ValueError("initial masses must be non-negative")
        n = round(self.plan.duration_s / self.plan.patch_time_step_s)
        if n < 1 or abs(n * self.plan.patch_time_step_s - self.plan.duration_s) > 1e-12:
            raise ValueError("duration must be an integer number of patch steps")
        ci = pi = 0.0
        history = []
        for index in range(n):
            t = index * self.plan.patch_time_step_s
            coarse_state, co = coarse_step(coarse_state, t, self.plan.patch_time_step_s, pi, cm)
            patch_state, po = patch_step(patch_state, t, self.plan.patch_time_step_s, ci, pm)
            co, po = float(co), float(po)
            if (not math.isfinite(co) or not math.isfinite(po) or co < 0 or po < 0
                    or co > cm + 1e-12 or po > pm + 1e-12):
                raise ValueError("callback returned outgoing mass above available donor mass")
            cm -= co; pm -= po
            ci, pi = co, po
            cm += pi; pm += ci
            history.append((t, co, po, cm, pm))
        return coarse_state, patch_state, tuple(history)

    def run_bidirectional_conservative(
        self,
        *,
        coarse_state: ConservativeControlVolumeState,
        patch_state: ConservativeControlVolumeState,
        coarse_step,
        patch_step,
    ):
        """Advance a shared clock with species/energy/momentum exchange.

        Each callback receives ``(state, t, dt, incoming_packet)`` and returns
        ``(updated_state, outgoing_packet)``.  The callback may represent a
        fast outer step or a resolved patch step, but the interface update is
        always applied by this executor after both callbacks return.  A packet
        that exceeds the donor total/species mass is rejected; no post-hoc
        clipping is allowed because it would break the field solver ledger.
        """
        import math
        if not isinstance(coarse_state, ConservativeControlVolumeState):
            raise TypeError("coarse_state must be ConservativeControlVolumeState")
        if not isinstance(patch_state, ConservativeControlVolumeState):
            raise TypeError("patch_state must be ConservativeControlVolumeState")
        if not callable(coarse_step) or not callable(patch_step):
            raise TypeError("coarse_step and patch_step must be callable")
        n = round(self.plan.duration_s / self.plan.patch_time_step_s)
        if n < 1 or abs(n * self.plan.patch_time_step_s - self.plan.duration_s) > 1e-12:
            raise ValueError("duration must be an integer number of patch steps")
        zero = ConservativeInterfacePacket.zero()
        incoming_coarse = zero
        incoming_patch = zero
        history = []
        initial = (
            coarse_state.total_mass_kg + patch_state.total_mass_kg,
            coarse_state.species_mass_kg + patch_state.species_mass_kg,
            coarse_state.sensible_enthalpy_j + patch_state.sensible_enthalpy_j,
            tuple(
                a + b for a, b in zip(
                    coarse_state.momentum_kg_m_s, patch_state.momentum_kg_m_s,
                )
            ),
        )

        def ledger_matches(total) -> bool:
            scalar_ok = all(
                math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-10)
                for a, b in zip(total[:3], initial[:3])
            )
            momentum_ok = all(
                math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-10)
                for a, b in zip(total[3], initial[3])
            )
            return scalar_ok and momentum_ok

        for index in range(n):
            t = index * self.plan.patch_time_step_s
            # The previous exchange has already applied the packets to the
            # receiver state.  They are passed here as an audit-visible event,
            # not applied a second time.
            coarse_state, coarse_out = coarse_step(
                coarse_state, t, self.plan.patch_time_step_s, incoming_coarse,
            )
            patch_state, patch_out = patch_step(
                patch_state, t, self.plan.patch_time_step_s, incoming_patch,
            )
            if not isinstance(coarse_state, ConservativeControlVolumeState):
                raise TypeError("coarse_step must return ConservativeControlVolumeState")
            if not isinstance(patch_state, ConservativeControlVolumeState):
                raise TypeError("patch_step must return ConservativeControlVolumeState")
            if not isinstance(coarse_out, ConservativeInterfacePacket):
                raise TypeError("coarse_step must return ConservativeInterfacePacket")
            if not isinstance(patch_out, ConservativeInterfacePacket):
                raise TypeError("patch_step must return ConservativeInterfacePacket")
            coarse_state, patch_state = exchange_conservative_pair(
                coarse_state, patch_state,
                coarse_to_patch=coarse_out,
                patch_to_coarse=patch_out,
            )
            incoming_coarse = patch_out
            incoming_patch = coarse_out
            total = (
                coarse_state.total_mass_kg + patch_state.total_mass_kg,
                coarse_state.species_mass_kg + patch_state.species_mass_kg,
                coarse_state.sensible_enthalpy_j + patch_state.sensible_enthalpy_j,
                tuple(
                    a + b for a, b in zip(
                        coarse_state.momentum_kg_m_s, patch_state.momentum_kg_m_s,
                    )
                ),
            )
            if not ledger_matches(total):
                raise RuntimeError("conservative hybrid field ledger drifted")
            history.append((t, coarse_out, patch_out, coarse_state, patch_state))
        return coarse_state, patch_state, tuple(history)

__all__ = ["HybridAMRExecutionResult", "HybridAMRExecutor", "MassRequest"]

def run_solver_with_patch_flux_collector(solver, config, collector, **kwargs):
    """Invoke a solver through the explicit face-flux coupling contract."""
    if not callable(solver) or not callable(collector):
        raise TypeError("solver and collector must be callable")
    try:
        return solver(config, face_flux_callback=collector, **kwargs)
    except TypeError as exc:
        if "face_flux_callback" in str(exc):
            raise TypeError("solver does not implement face_flux_callback contract") from exc
        raise

__all__.append("run_solver_with_patch_flux_collector")

def prepare_patch_flux_collector(config: object, patch: object, *, tolerance_m: float = 1e-10):
    """Build a collector only when patch bounds coincide with solver faces."""
    if not hasattr(config, "domain_m") or not hasattr(config, "cells"):
        raise TypeError("config must expose domain_m and cells")
    x_face, y_face, z_face = cartesian_face_coordinates(config.domain_m, config.cells)
    tol = float(tolerance_m)
    if tol < 0:
        raise ValueError("tolerance_m must be non-negative")
    for values, lo, hi in ((x_face, patch.x_min, patch.x_max),
                           (y_face, patch.y_min, patch.y_max),
                           (z_face, patch.z_min, patch.z_max)):
        if not any(abs(float(v) - float(lo)) <= tol for v in values):
            raise ValueError("patch lower bound is not a solver face")
        if not any(abs(float(v) - float(hi)) <= tol for v in values):
            raise ValueError("patch upper bound is not a solver face")
    return PatchFluxAccumulator(x_face, y_face, z_face, patch)

__all__.append("prepare_patch_flux_collector")

def run_solver_with_prepared_patch_flux(solver, config, patch, **kwargs):
    """Prepare an aligned collector, run the solver, and return both objects."""
    collector = prepare_patch_flux_collector(config, patch)
    result = run_solver_with_patch_flux_collector(solver, config, collector, **kwargs)
    return result, collector

__all__.append("run_solver_with_prepared_patch_flux")
