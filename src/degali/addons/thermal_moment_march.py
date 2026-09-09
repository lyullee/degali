"""Guarded adaptive march for the opt-in conservative thermal-moment state.

This module supplies numerical continuation and an audit trail around an
already selected local closure.  It does not select the thermal/species
diffusivity ratio, repair an invalid local state, or extrapolate past failure.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class ThermalMomentMarchResult:
    arc_length: np.ndarray
    states: np.ndarray
    rhs_calls: int
    rejected_steps: int
    reached_target: bool
    stop_reason: str
    attempted_step: float
    maximum_local_error: float
    maximum_weak_residual: float
    minimum_sampled_diffusivity: float
    maximum_edge_heat_defect: float


class GuardedThermalMomentMarch:
    """Adaptive midpoint continuation with strict local-closure rejection."""

    def __init__(self, evaluator):
        if not callable(evaluator):
            raise TypeError("a callable thermal-moment evaluator is required")
        self.evaluator = evaluator

    @staticmethod
    def _validate_options(
        length, initial_step, maximum_step, minimum_step, tolerance, maximum_steps
    ):
        values = (length, initial_step, maximum_step, minimum_step, tolerance)
        if not all(math.isfinite(float(value)) and float(value) > 0.0 for value in values):
            raise ValueError("march lengths, steps and tolerance must be positive and finite")
        if not minimum_step <= initial_step <= maximum_step:
            raise ValueError("steps must satisfy minimum <= initial <= maximum")
        if isinstance(maximum_steps, bool) or int(maximum_steps) != maximum_steps or maximum_steps < 1:
            raise ValueError("maximum_steps must be a positive integer")

    def _rhs(self, state, diagnostics):
        out = self.evaluator(np.asarray(state[:8], dtype=float))
        required = (
            "rates", "ledger", "moment_rate", "weak_budget_scaled_error",
            "minimum_chi_species", "minimum_chi_momentum",
            "edge_gradient_defects", "valid",
        )
        if any(key not in out for key in required):
            raise ValueError("thermal-moment evaluator returned an incomplete result")
        if not out["valid"]:
            raise ValueError("local thermal-moment closure is invalid")
        rates = np.asarray(out["rates"], dtype=float)
        sources = np.asarray(out["ledger"]["sources"], dtype=float)
        rhs = np.r_[rates, sources, float(out["moment_rate"])]
        if rates.shape != (8,) or sources.shape != (5,) or rhs.shape != (14,):
            raise ValueError("thermal-moment evaluator returned inconsistent dimensions")
        if not np.all(np.isfinite(rhs)):
            raise ValueError("thermal-moment right-hand side is nonfinite")
        diagnostics["rhs_calls"] += 1
        diagnostics["maximum_weak_residual"] = max(
            diagnostics["maximum_weak_residual"],
            float(out["weak_budget_scaled_error"]),
        )
        diagnostics["minimum_sampled_diffusivity"] = min(
            diagnostics["minimum_sampled_diffusivity"],
            float(out["minimum_chi_species"]),
            float(out["minimum_chi_momentum"]),
        )
        diagnostics["maximum_edge_heat_defect"] = max(
            diagnostics["maximum_edge_heat_defect"],
            float(out["edge_gradient_defects"]["heat"]),
        )
        return rhs

    def _midpoint(self, state, step, diagnostics):
        first = self._rhs(state, diagnostics)
        return state + step * self._rhs(state + 0.5 * step * first, diagnostics)

    def march(
        self,
        initial_parameters,
        length,
        *,
        initial_step=0.002,
        maximum_step=0.005,
        minimum_step=1.0e-6,
        tolerance=2.0e-5,
        maximum_steps=4000,
        target=None,
    ):
        """March until arc length or an optional predicate is reached.

        The 14-component state contains eight physical/log parameters followed
        by five cumulative sources and the cumulative thermal second moment.
        The last accepted state is always returned when a local gate fails.
        """
        self._validate_options(
            length, initial_step, maximum_step, minimum_step, tolerance, maximum_steps
        )
        parameters = np.asarray(initial_parameters, dtype=float)
        if parameters.shape != (8,) or not np.all(np.isfinite(parameters)):
            raise ValueError("initial_parameters must contain eight finite values")
        if target is not None and not callable(target):
            raise TypeError("target must be callable")
        state = np.r_[parameters, np.zeros(6)]
        arc = [0.0]
        states = [state.copy()]
        step = min(float(initial_step), float(length))
        diagnostics = {
            "rhs_calls": 0,
            "maximum_weak_residual": 0.0,
            "minimum_sampled_diffusivity": math.inf,
            "maximum_edge_heat_defect": 0.0,
        }
        rejected = 0
        maximum_error = 0.0
        attempted = step
        reason = "arc-length ceiling reached"
        reached = bool(target is not None and target(state[:8]))
        if reached:
            reason = "target reached"

        while not reached and arc[-1] < length and len(arc) - 1 < maximum_steps:
            step = min(step, float(length) - arc[-1])
            attempted = step
            try:
                full = self._midpoint(state, step, diagnostics)
                half = self._midpoint(state, 0.5 * step, diagnostics)
                refined = self._midpoint(half, 0.5 * step, diagnostics)
                scale = np.maximum(np.maximum(np.abs(full[:8]), np.abs(refined[:8])), 1.0)
                error = float(np.max(np.abs(refined[:8] - full[:8]) / scale))
                if not math.isfinite(error) or error > tolerance:
                    raise ArithmeticError(f"local step error {error:.3e} exceeds tolerance")
            except (ArithmeticError, FloatingPointError, RuntimeError, ValueError, np.linalg.LinAlgError) as exc:
                rejected += 1
                if step <= minimum_step * (1.0 + 1.0e-12):
                    reason = f"minimum step rejected: {exc}"
                    break
                step = max(0.5 * step, minimum_step)
                continue

            state = refined
            arc.append(arc[-1] + step)
            states.append(state.copy())
            maximum_error = max(maximum_error, error)
            reached = bool(target is not None and target(state[:8]))
            if reached:
                reason = "target reached"
                break
            if error < tolerance / 8.0:
                step = min(1.5 * step, maximum_step)

        if not reached and len(arc) - 1 >= maximum_steps:
            reason = "accepted-step budget exhausted"
        minimum_diffusivity = diagnostics["minimum_sampled_diffusivity"]
        if math.isinf(minimum_diffusivity):
            minimum_diffusivity = math.nan
        return ThermalMomentMarchResult(
            arc_length=np.asarray(arc),
            states=np.asarray(states),
            rhs_calls=diagnostics["rhs_calls"],
            rejected_steps=rejected,
            reached_target=reached,
            stop_reason=reason,
            attempted_step=float(attempted),
            maximum_local_error=maximum_error,
            maximum_weak_residual=diagnostics["maximum_weak_residual"],
            minimum_sampled_diffusivity=float(minimum_diffusivity),
            maximum_edge_heat_defect=diagnostics["maximum_edge_heat_defect"],
        )


__all__ = ["GuardedThermalMomentMarch", "ThermalMomentMarchResult"]
