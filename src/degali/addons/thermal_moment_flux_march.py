"""Direct conservative-flux continuation for the opt-in thermal model.

The six transported quantities are primary numerical state variables.  A
physical Gaussian section is reconstructed at every Runge--Kutta stage; an
invalid inverse or local closure stops the march instead of being clipped.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import least_squares

from .buoyancy_profile import BuoyancyConstrainedEnthalpySection
from .reservoir_thermal import ReservoirThermalMoments, ReservoirShortSegment, decode_section, encode_section
from .thermal_moments import EnthalpyMomentOperators
from .transverse_mixing import ConservativeTransverseMixing


_FLUX_FLOORS = np.array([1e-12, 1e-12, 1., 1., 1., 1.])


@dataclass(frozen=True)
class ThermalMomentFluxMatch:
    parameters: np.ndarray
    state: np.ndarray
    thermal_width_ratio: float
    fluxes: np.ndarray
    scaled_residuals: np.ndarray
    solver_evaluations: int
    mixing: object


@dataclass(frozen=True)
class ThermalMomentFluxMarchResult:
    arc_length: np.ndarray
    fluxes: np.ndarray
    parameters: np.ndarray
    cumulative_sources: np.ndarray
    rhs_calls: int
    accepted_steps: int
    reached_target: bool
    stop_reason: str
    maximum_inverse_residual: float
    maximum_weak_residual: float
    minimum_sampled_diffusivity: float
    maximum_edge_heat_defect: float


class ThermalMomentFluxInverter:
    """Invert six section fluxes without introducing a fitted closure."""

    def __init__(self, jetplume, thermodynamics, *, order=8, probes=257,
                 quadrature_points=128, tolerance=1e-8):
        if isinstance(order, bool) or int(order) != order or order < 4:
            raise ValueError("flux quadrature order must be an integer of at least four")
        if isinstance(probes, bool) or int(probes) != probes or probes < 33:
            raise ValueError("phase probes must be an integer of at least 33")
        if not math.isfinite(tolerance) or tolerance <= 0.:
            raise ValueError("inverse tolerance must be positive and finite")
        self.jp, self.th = jetplume, thermodynamics
        self.order, self.probes = int(order), int(probes)
        self.quadrature_points, self.tolerance = int(quadrature_points), float(tolerance)

    @staticmethod
    def scales(values):
        values = np.asarray(values, float)
        return np.maximum(np.abs(values), _FLUX_FLOORS)

    @staticmethod
    def _validate_target(target):
        target = np.asarray(target, float)
        if (target.shape != (6,) or not np.all(np.isfinite(target))
                or target[0] <= 0. or target[1] <= 0.
                or math.hypot(target[2], target[3]) <= 0.):
            raise ValueError("target must contain six finite fluxes with positive mass, hydrogen and momentum")
        return target

    def _decode(self, reduced, target, position):
        rho, area, uc, beta = np.exp(reduced)
        theta = math.atan2(target[3], target[2])
        section = BuoyancyConstrainedEnthalpySection(
            self.jp, self.th, thermal_width_ratio=beta,
            quadrature_points=self.quadrature_points,
        )
        state = np.array([rho, .1, area, theta, uc, position[0], position[1]])
        i1 = section.k.profile_integral(1.)
        iu = section.k.profile_integral(1.+section.velocity_shape_exponent)
        denominator = area*(section._wind(state)*math.cos(theta)*i1+uc*iu)
        if not math.isfinite(denominator) or denominator <= 0.:
            raise ValueError("hydrogen-flux elimination has a nonpositive denominator")
        centre_hydrogen_density = target[1]/denominator
        state[1] = centre_hydrogen_density/rho
        section._physical(state)
        mixing = ConservativeTransverseMixing(section, state, probes=self.probes)
        return section, state, mixing

    def _values_and_jacobian(self, reduced, target, position):
        section, state, mixing = self._decode(reduced, target, position)
        values = ReservoirShortSegment.flux_values(mixing, order=self.order)
        full = mixing.flux_jacobian(self.order)
        i1 = section.k.profile_integral(1.)
        iu = section.k.profile_integral(1.+section.velocity_shape_exponent)
        theta, uc = state[3], state[4]
        wind = section._wind(state)*math.cos(theta)
        plus, minus = state.copy(), state.copy()
        step = 1e-5
        plus[2] *= math.exp(step)
        minus[2] *= math.exp(-step)
        wind_area = ((section._wind(plus)-section._wind(minus))*math.cos(theta)/(2*step))
        denominator = wind*i1+uc*iu
        dlog_centre_hydrogen = np.array([
            0., -(1.+wind_area*i1/denominator), -uc*iu/denominator, 0.,
        ])
        columns = np.array([0, 2, 4, 7])
        reduced_five = full[:, columns]+full[:, 1, None]*dlog_centre_hydrogen
        axial = math.cos(theta)*reduced_five[2]+math.sin(theta)*reduced_five[3]
        moment2 = EnthalpyMomentOperators(section).reduced_moment_jacobian(
            state, quadrature_points=self.quadrature_points,
        )
        jacobian = np.vstack([reduced_five[0], axial, reduced_five[4], moment2])
        return values, jacobian, state, mixing

    def fluxes(self, parameters):
        state, beta = decode_section(parameters)
        section = BuoyancyConstrainedEnthalpySection(
            self.jp, self.th, thermal_width_ratio=beta,
            quadrature_points=self.quadrature_points,
        )
        mixing = ConservativeTransverseMixing(section, state, probes=self.probes)
        return ReservoirShortSegment.flux_values(mixing, order=self.order)

    def match(self, target, initial_parameters, *, position=None):
        target = self._validate_target(target)
        initial = np.asarray(initial_parameters, float)
        state0, _ = decode_section(initial)
        if position is None:
            position = state0[5:7]
        position = np.asarray(position, float)
        if position.shape != (2,) or not np.all(np.isfinite(position)):
            raise ValueError("inverse position must contain two finite coordinates")
        theta = math.atan2(target[3], target[2])
        axial_target = math.hypot(target[2], target[3])
        reduced_target = np.array([target[0], axial_target, target[4], target[5]])
        reduced_scales = np.maximum(np.abs(reduced_target), [1e-12, 1., 1., 1.])
        lower = np.log([1e-3, 1e-14, 1e-8, .5])
        upper = np.log([100., 100., 1e4, 2.])
        start = np.clip(initial[[0, 2, 4, 7]], lower+1e-10, upper-1e-10)
        cache = {}

        def evaluate(parameters):
            key = tuple(np.asarray(parameters, float))
            if key not in cache:
                cache.clear()
                cache[key] = self._values_and_jacobian(parameters, target, position)
            return cache[key]

        def residual(parameters):
            try:
                values = evaluate(parameters)[0]
                reduced_values = np.array([
                    values[0], math.hypot(values[2], values[3]), values[4], values[5],
                ])
                return (reduced_values-reduced_target)/reduced_scales
            except (ArithmeticError, FloatingPointError, RuntimeError, ValueError, np.linalg.LinAlgError):
                return np.full(4, 1e6)

        def jacobian(parameters):
            try:
                return evaluate(parameters)[1]/reduced_scales[:, None]
            except (ArithmeticError, FloatingPointError, RuntimeError, ValueError, np.linalg.LinAlgError):
                return np.zeros((4, 4))

        result = least_squares(
            residual, start, jac=jacobian, bounds=(lower, upper), x_scale="jac",
            xtol=1e-11, ftol=1e-11, gtol=1e-11, max_nfev=80,
        )
        values, _, state, mixing = evaluate(result.x)
        scaled = (values-target)/self.scales(target)
        beta = float(np.exp(result.x[3]))
        if (not result.success or np.max(np.abs(scaled)) > self.tolerance
                or beta <= .5+1e-8 or beta >= 2.-1e-8):
            raise ValueError(
                f"six-flux inverse failed: success={result.success}, "
                f"residual={np.max(np.abs(scaled)):.3e}, beta={beta:.9g}"
            )
        parameters = encode_section(state, beta)
        parameters[3], parameters[5:7] = theta, position
        return ThermalMomentFluxMatch(
            parameters=parameters, state=state, thermal_width_ratio=beta,
            fluxes=values, scaled_residuals=scaled,
            solver_evaluations=int(result.nfev), mixing=mixing,
        )


class FluxSpaceThermalMomentMarch:
    """Classical RK4 march whose primary state is six physical fluxes."""

    def __init__(self, inverter, *, thermal_species_ratio, mechanical_work,
                 local_evaluator=None):
        if not hasattr(inverter, "match") or not hasattr(inverter, "fluxes"):
            raise TypeError("a six-flux inverter is required")
        if mechanical_work != "reduced_buoyancy_work":
            raise ValueError("explicitly select reduced_buoyancy_work")
        if not math.isfinite(thermal_species_ratio) or thermal_species_ratio <= 0.:
            raise ValueError("thermal/species diffusivity ratio must be positive and finite")
        self.inverter = inverter
        self.ratio, self.work = float(thermal_species_ratio), mechanical_work
        self.local_evaluator = local_evaluator

    def _local(self, match):
        if self.local_evaluator is not None:
            return self.local_evaluator(match)
        return ReservoirThermalMoments(
            match.mixing, thermal_species_ratio=self.ratio,
            mechanical_work=self.work,
        ).evaluate(order=self.inverter.order, probes=self.inverter.probes)

    def _rhs(self, primary, guess, diagnostics):
        match = self.inverter.match(primary[:6], guess, position=primary[6:8])
        out = self._local(match)
        required = (
            "ledger", "moment_rate", "weak_budget_scaled_error",
            "minimum_chi_species", "minimum_chi_momentum",
            "edge_gradient_defects", "valid",
        )
        if any(key not in out for key in required) or not out["valid"]:
            raise ValueError("local thermal-moment closure is invalid or incomplete")
        sources = np.r_[np.asarray(out["ledger"]["sources"], float), float(out["moment_rate"])]
        if sources.shape != (6,) or not np.all(np.isfinite(sources)):
            raise ValueError("local six-flux source vector is invalid")
        theta = match.state[3]
        rhs = np.r_[sources, math.cos(theta), math.sin(theta)]
        diagnostics["rhs_calls"] += 1
        diagnostics["maximum_inverse_residual"] = max(
            diagnostics["maximum_inverse_residual"],
            float(np.max(np.abs(match.scaled_residuals))),
        )
        diagnostics["maximum_weak_residual"] = max(
            diagnostics["maximum_weak_residual"], float(out["weak_budget_scaled_error"]),
        )
        diagnostics["minimum_sampled_diffusivity"] = min(
            diagnostics["minimum_sampled_diffusivity"],
            float(out["minimum_chi_species"]), float(out["minimum_chi_momentum"]),
        )
        diagnostics["maximum_edge_heat_defect"] = max(
            diagnostics["maximum_edge_heat_defect"],
            float(out["edge_gradient_defects"]["heat"]),
        )
        return rhs, match.parameters

    def _rk4(self, primary, step, guess, diagnostics):
        k1, p1 = self._rhs(primary, guess, diagnostics)
        k2, p2 = self._rhs(primary+.5*step*k1, p1, diagnostics)
        k3, p3 = self._rhs(primary+.5*step*k2, p2, diagnostics)
        k4, p4 = self._rhs(primary+step*k3, p3, diagnostics)
        weighted = step*(k1+2*k2+2*k3+k4)/6.
        return primary+weighted, p4, weighted[:6]

    def march(self, initial_parameters, length, *, step, target_x=None,
              maximum_steps=4000, target_tolerance=1e-8, progress=None):
        values = (length, step, target_tolerance)
        if not all(math.isfinite(float(v)) and float(v) > 0. for v in values):
            raise ValueError("length, step and target tolerance must be positive and finite")
        if isinstance(maximum_steps, bool) or int(maximum_steps) != maximum_steps or maximum_steps < 1:
            raise ValueError("maximum_steps must be a positive integer")
        parameters = np.asarray(initial_parameters, float)
        state, _ = decode_section(parameters)
        if target_x is not None and not math.isfinite(float(target_x)):
            raise ValueError("target_x must be finite")
        if progress is not None and not callable(progress):
            raise TypeError("progress must be callable")
        initial_fluxes = np.asarray(self.inverter.fluxes(parameters), float)
        primary = np.r_[initial_fluxes, state[5:7]]
        arc, flux_history, parameter_history = [0.], [initial_fluxes.copy()], [parameters.copy()]
        cumulative = [np.zeros(6)]
        diagnostics = dict(rhs_calls=0, maximum_inverse_residual=0., maximum_weak_residual=0.,
                           minimum_sampled_diffusivity=math.inf, maximum_edge_heat_defect=0.)
        reason, reached = "arc-length ceiling reached", False

        for _ in range(int(maximum_steps)):
            if arc[-1] >= length:
                break
            if target_x is not None and abs(primary[6]-target_x) <= target_tolerance:
                reached, reason = True, "target x reached"
                break
            ds = min(float(step), float(length)-arc[-1])
            if target_x is not None:
                direction = math.cos(decode_section(parameters)[0][3])
                remaining = target_x-primary[6]
                if remaining*direction <= 0.:
                    reason = "target x cannot be approached in the current direction"
                    break
                ds = min(ds, abs(remaining/direction))
            try:
                next_primary, next_guess, source_increment = self._rk4(
                    primary, ds, parameters, diagnostics,
                )
                if target_x is not None and (next_primary[6]-target_x)*(primary[6]-target_x) <= 0.:
                    # Correct the final arc step using the RK-computed x increment.
                    for _ in range(4):
                        error = next_primary[6]-target_x
                        if abs(error) <= target_tolerance:
                            break
                        slope = (next_primary[6]-primary[6])/ds
                        corrected = ds-(next_primary[6]-target_x)/slope
                        if not 0. < corrected <= float(step)*(1.+1e-12):
                            raise ValueError("target-x step correction left its bracket")
                        ds = corrected
                        next_primary, next_guess, source_increment = self._rk4(
                            primary, ds, parameters, diagnostics,
                        )
                endpoint = self.inverter.match(
                    next_primary[:6], next_guess, position=next_primary[6:8],
                )
            except (ArithmeticError, FloatingPointError, RuntimeError, ValueError, np.linalg.LinAlgError) as exc:
                reason = f"step rejected: {exc}"
                break
            primary, parameters = next_primary, endpoint.parameters
            arc.append(arc[-1]+ds)
            flux_history.append(primary[:6].copy())
            parameter_history.append(parameters.copy())
            cumulative.append(cumulative[-1]+source_increment)
            diagnostics["maximum_inverse_residual"] = max(
                diagnostics["maximum_inverse_residual"],
                float(np.max(np.abs(endpoint.scaled_residuals))),
            )
            if progress is not None:
                progress(len(arc)-1, arc[-1], primary.copy())
            if target_x is not None and abs(primary[6]-target_x) <= target_tolerance:
                reached, reason = True, "target x reached"
                break
        else:
            reason = "accepted-step budget exhausted"

        minimum = diagnostics["minimum_sampled_diffusivity"]
        if math.isinf(minimum):
            minimum = math.nan
        return ThermalMomentFluxMarchResult(
            arc_length=np.asarray(arc), fluxes=np.asarray(flux_history),
            parameters=np.asarray(parameter_history), cumulative_sources=np.asarray(cumulative),
            rhs_calls=diagnostics["rhs_calls"], accepted_steps=len(arc)-1,
            reached_target=reached, stop_reason=reason,
            maximum_inverse_residual=diagnostics["maximum_inverse_residual"],
            maximum_weak_residual=diagnostics["maximum_weak_residual"],
            minimum_sampled_diffusivity=float(minimum),
            maximum_edge_heat_defect=diagnostics["maximum_edge_heat_defect"],
        )


__all__ = [
    "ThermalMomentFluxInverter", "ThermalMomentFluxMatch",
    "FluxSpaceThermalMomentMarch", "ThermalMomentFluxMarchResult",
]
