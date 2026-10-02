"""Opt-in finite wall-jet attachment and lift-off for the LH2 research path.

The legacy JETPLU state has no memory of ground contact.  It can either treat
the section as free or remove the buried perimeter every time the section
overlaps the ground.  This module adds one bounded state: the fraction of the
section that remains dynamically attached to the surface.  The state relaxes
towards a local Richardson-number lift-off criterion over a local plume-depth
scale.  While attached, surface shear removes streamwise momentum and the
positive vertical force is reduced continuously.

This is a research closure, not a DEGADIS parity path or a validated wall-jet
correlation.  It deliberately wraps :class:`IndependentEnergyCrosswind` so
the established default calculations remain unchanged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

from .energy_crosswind import (
    E_THETA,
    IndependentEnergyCrosswind,
)


W_ATTACHMENT = 7


@dataclass(frozen=True)
class WallJetTransitionConfig:
    """Physics controls for the finite wall-contact state.

    ``critical_richardson`` uses the same bulk lift-off definition documented
    for :class:`~degali.core.jetplume.JetCoefficients`: ``g H Delta-rho /
    (rho_a u_*^2)``.  Thirty is the published Briggs-scale diagnostic already
    retained by DEGALI.  ``response_depths`` sets a numerical/physical memory
    length in local above-ground plume depths; it is exposed and must not be
    calibrated against the Test 6 concentration.  ``wall_shear_multiplier``
    multiplies the directly available surface stress ``rho_a u_*^2``.
    """

    critical_richardson: float = 30.0
    response_depths: float = 1.0
    wall_shear_multiplier: float = 1.0

    def __post_init__(self) -> None:
        values = (
            self.critical_richardson,
            self.response_depths,
            self.wall_shear_multiplier,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("wall-jet controls must be finite")
        if self.critical_richardson <= 0.0 or self.response_depths <= 0.0:
            raise ValueError(
                "wall-jet Richardson threshold and response depth must be positive"
            )
        if self.wall_shear_multiplier < 0.0:
            raise ValueError("wall-jet shear multiplier must be non-negative")


@dataclass(frozen=True)
class WallJetDiagnostics:
    """Local attachment, lift-off and wall-stress terms."""

    geometric_contact_fraction: float
    attachment_fraction: float
    equilibrium_attachment_fraction: float
    lift_richardson: float
    vertical_release_fraction: float
    response_length_m: float
    wall_shear_force_per_m: float


@dataclass(frozen=True)
class FiniteWallJetResult:
    """Eight-state wall-jet trajectory and conservative-flux audit."""

    arc_length: np.ndarray
    states: np.ndarray
    fluxes: np.ndarray
    sources: np.ndarray
    cumulative_sources: np.ndarray
    balance_residuals: np.ndarray
    temperatures: np.ndarray
    centre_mole_fractions: np.ndarray
    maximum_relative_balance_residual: float

    @property
    def attachment_fraction(self) -> np.ndarray:
        return self.states[:, W_ATTACHMENT]


class FiniteWallJetCrosswind:
    """Add finite surface attachment to an independent-energy plume.

    The wrapped base must use ``ground_interaction='geometry'`` or
    ``'surface_layer'`` so its ellipse supplies the exposed perimeter,
    ground-contact chord and cleared fraction.  No base object is mutated.
    """

    def __init__(
        self,
        base: IndependentEnergyCrosswind,
        config: WallJetTransitionConfig | None = None,
    ) -> None:
        if not isinstance(base, IndependentEnergyCrosswind):
            raise TypeError("finite wall jet requires IndependentEnergyCrosswind")
        if base.ground_interaction not in {"geometry", "surface_layer"}:
            raise ValueError(
                "finite wall jet requires geometry or surface_layer ground interaction"
            )
        self.base = base
        self.config = config or WallJetTransitionConfig()

    @staticmethod
    def _smoothstep(value: float) -> float:
        value = min(max(float(value), 0.0), 1.0)
        return value * value * (3.0 - 2.0 * value)

    @staticmethod
    def _logit(value: float) -> float:
        bounded = min(max(float(value), 1.0e-9), 1.0 - 1.0e-9)
        return math.log(bounded / (1.0 - bounded))

    @staticmethod
    def _logistic(value: float) -> float:
        clipped = float(np.clip(value, -35.0, 35.0))
        return 1.0 / (1.0 + math.exp(-clipped))

    def _split_state(self, state: np.ndarray) -> tuple[np.ndarray, float]:
        state = np.asarray(state, dtype=float)
        if state.shape != (8,) or not np.all(np.isfinite(state)):
            raise ValueError("finite wall-jet state must contain eight finite values")
        self.base._physical(state[:7])
        attachment = float(state[W_ATTACHMENT])
        if not 0.0 <= attachment <= 1.0:
            raise ValueError("wall attachment fraction must lie in [0, 1]")
        return state[:7], attachment

    def diagnostics(self, state: np.ndarray) -> WallJetDiagnostics:
        """Evaluate the local wall-contact state without changing the model."""
        section, attachment = self._split_state(state)
        rho, _fraction, sysz, theta, _uc, _x, z = self.base._physical(section)
        (
            _sy, sz, _ua, _perimeter, ground_width, cleared_fraction,
            _alpha_wind, _sya, _sza, _dsya, _dsza,
        ) = self.base._geometry(section)
        ct = max(math.cos(theta), 1.0e-6)
        half_depth = self.base.k.delta * sz * ct
        layer_depth = max(z + half_depth, math.sqrt(sysz), 1.0e-6)
        contact = min(max(1.0 - cleared_fraction, 0.0), 1.0)
        density_deficit = max(self.base.rhoa - rho, 0.0)
        lift_richardson = (
            9.81 * layer_depth * density_deficit
            / (
                self.base.rhoa
                * max(float(self.base.jetplume.ustar), 1.0e-6) ** 2
            )
        )
        released = self._smoothstep(
            lift_richardson / self.config.critical_richardson
        )
        equilibrium = contact * (1.0 - released)
        vertical_release = 1.0 - attachment * contact
        response_length = self.config.response_depths * layer_depth
        wall_shear = (
            attachment
            * self.config.wall_shear_multiplier
            * self.base.rhoa
            * max(float(self.base.jetplume.ustar), 0.0) ** 2
            * ground_width
        )
        return WallJetDiagnostics(
            geometric_contact_fraction=float(contact),
            attachment_fraction=float(attachment),
            equilibrium_attachment_fraction=float(equilibrium),
            lift_richardson=float(lift_richardson),
            vertical_release_fraction=float(vertical_release),
            response_length_m=float(response_length),
            wall_shear_force_per_m=float(wall_shear),
        )

    def initial_state(
        self,
        base_state: np.ndarray,
        *,
        attachment_fraction: float | None = None,
    ) -> np.ndarray:
        """Lift a seven-state interface into the finite wall-jet state."""
        base_state = np.asarray(base_state, dtype=float)
        self.base._physical(base_state)
        trial = np.concatenate((base_state, [0.0]))
        if attachment_fraction is None:
            attachment_fraction = self.diagnostics(
                trial
            ).equilibrium_attachment_fraction
        attachment_fraction = float(attachment_fraction)
        if not 0.0 <= attachment_fraction <= 1.0:
            raise ValueError("initial wall attachment must lie in [0, 1]")
        return np.concatenate((base_state, [attachment_fraction]))

    def fluxes(self, state: np.ndarray) -> np.ndarray:
        """Return the five transported fluxes; attachment stores no flux."""
        section, _attachment = self._split_state(state)
        return self.base._as_array(self.base.integral_fluxes(section))

    def section_widths(self, state: np.ndarray) -> tuple[float, float]:
        """Return Gaussian widths for the eight-state trajectory adapter."""
        section, _attachment = self._split_state(state)
        return self.base.section_widths(section)

    def point_mole_fraction(
        self, state: np.ndarray, lateral: float, height: float
    ) -> float:
        """Evaluate an exact receptor using the wrapped thermodynamics."""
        section, _attachment = self._split_state(state)
        return self.base.point_mole_fraction(section, lateral, height)

    def point_temperature(
        self, state: np.ndarray, lateral: float, height: float
    ) -> float:
        """Evaluate receptor temperature using the wrapped thermodynamics."""
        section, _attachment = self._split_state(state)
        return self.base.point_temperature(section, lateral, height)

    @property
    def thermodynamics(self):
        """Expose the wrapped closure to existing trajectory observers."""
        return self.base.thermodynamics

    def source_terms(self, state: np.ndarray, **kwargs) -> np.ndarray:
        """Return conservative sources including wall shear and lift release."""
        section, _attachment = self._split_state(state)
        source = self.base.source_terms(section, **kwargs).copy()
        diagnostic = self.diagnostics(state)
        _rho, _fraction, _sysz, theta, _uc, _x, _z = self.base._physical(section)
        (
            _sy, _sz, ua, perimeter, _ground_width, _cleared,
            _alpha_wind, _sya, _sza, _dsya, _dsza,
        ) = self.base._geometry(section)
        st, ct = math.sin(theta), math.cos(theta)
        drag = self.base.k.cd * perimeter * self.base.rhoa / 2.0 * (ua * st) ** 2
        raw_vertical = (
            self.base.buoyancy_force(section)
            - math.copysign(1.0, theta) * drag * ct
        )
        if raw_vertical > 0.0:
            source[3] = raw_vertical * diagnostic.vertical_release_fraction
        else:
            source[3] = raw_vertical
        source[2] -= diagnostic.wall_shear_force_per_m * max(ct, 0.0)
        if not np.all(np.isfinite(source)):
            raise RuntimeError("finite wall-jet source became non-finite")
        return source

    def _attachment_rate(self, state: np.ndarray) -> tuple[float, float]:
        diagnostic = self.diagnostics(state)
        attachment = diagnostic.attachment_fraction
        target = diagnostic.equilibrium_attachment_fraction
        eta_rate = (
            self._logit(target) - self._logit(attachment)
        ) / diagnostic.response_length_m
        physical_rate = attachment * (1.0 - attachment) * eta_rate
        return float(physical_rate), float(eta_rate)

    def derivatives(self, _s: float, state: np.ndarray, **kwargs) -> np.ndarray:
        """Solve five flux balances, two trajectory equations and attachment."""
        section, _attachment = self._split_state(state)
        fluxes = self.fluxes(state)
        jacobian = self.base.flux_jacobian(section)
        sources = self.source_terms(state, **kwargs)
        column_scales = np.maximum(
            np.abs(section[:5]), np.array([1.0, 0.1, 1.0e-3, 1.0, 1.0])
        )
        row_scales = np.maximum(
            np.abs(fluxes), np.array([1.0e-3, 1.0e-4, 1.0, 1.0, 1.0])
        )
        scaled = jacobian * column_scales[np.newaxis, :] / row_scales[:, np.newaxis]
        condition = float(np.linalg.cond(scaled))
        if not math.isfinite(condition) or condition > 1.0e12:
            raise RuntimeError(
                f"finite wall-jet flux Jacobian is ill-conditioned ({condition:.3e})"
            )
        physical = column_scales * np.linalg.solve(scaled, sources / row_scales)
        attachment_rate, _eta_rate = self._attachment_rate(state)
        return np.concatenate((
            physical,
            [math.cos(section[E_THETA]), math.sin(section[E_THETA])],
            [attachment_rate],
        ))

    def _to_internal(self, state: np.ndarray) -> np.ndarray:
        section, attachment = self._split_state(state)
        return np.concatenate((self.base._to_internal(section), [self._logit(attachment)]))

    def _from_internal(self, state: np.ndarray) -> np.ndarray:
        state = np.asarray(state, dtype=float)
        if state.shape != (8,) or not np.all(np.isfinite(state)):
            raise ValueError("finite wall-jet internal state must contain eight values")
        return np.concatenate((
            self.base._from_internal(state[:7]),
            [self._logistic(state[W_ATTACHMENT])],
        ))

    def solve(
        self,
        initial_state: np.ndarray,
        *,
        maximum_distance: float,
        maximum_step: float = 0.05,
        relative_tolerance: float = 2.0e-6,
        method: str = "flux_RK4",
        **source_options,
    ) -> FiniteWallJetResult:
        """Integrate the research wall-jet state and audit all five balances."""
        initial_state = np.asarray(initial_state, dtype=float)
        self._split_state(initial_state)
        if maximum_distance <= 0.0 or maximum_step <= 0.0:
            raise ValueError("integration distance and maximum step must be positive")

        def right_hand_side(s: float, internal: np.ndarray) -> np.ndarray:
            state = self._from_internal(internal)
            physical = self.derivatives(s, state, **source_options)
            section = state[:7]
            fraction = section[1]
            _attachment_rate, eta_rate = self._attachment_rate(state)
            return np.array([
                physical[0] / section[0],
                physical[1] / max(fraction * (1.0 - fraction), 1.0e-14),
                physical[2] / section[2],
                physical[3],
                physical[4] / section[4],
                physical[5],
                physical[6],
                eta_rate,
            ])

        cumulative = None
        if method == "flux_RK4":
            steps = int(math.ceil(maximum_distance / maximum_step))
            integration_t = np.linspace(0.0, maximum_distance, steps + 1)
            states = np.empty((steps + 1, 8))
            states[0] = initial_state
            cumulative = np.zeros((steps + 1, 5))

            def stage_state(target_fluxes, guess, position, eta):
                section = guess[:7].copy()
                section[E_THETA] = math.atan2(
                    target_fluxes[3], target_fluxes[2]
                )
                section[5:7] = position
                matched = self.base._match_flux_array(target_fluxes, section)
                return np.concatenate((matched, [self._logistic(eta)]))

            for index in range(steps):
                step = float(integration_t[index + 1] - integration_t[index])
                current = states[index]
                current_fluxes = self.fluxes(current)
                eta = self._logit(current[W_ATTACHMENT])

                source1 = self.source_terms(current, **source_options)
                path1 = np.array([
                    math.cos(current[E_THETA]), math.sin(current[E_THETA])
                ])
                eta1 = self._attachment_rate(current)[1]
                state2 = stage_state(
                    current_fluxes + 0.5 * step * source1,
                    current,
                    current[5:7] + 0.5 * step * path1,
                    eta + 0.5 * step * eta1,
                )

                source2 = self.source_terms(state2, **source_options)
                path2 = np.array([
                    math.cos(state2[E_THETA]), math.sin(state2[E_THETA])
                ])
                eta2 = self._attachment_rate(state2)[1]
                state3 = stage_state(
                    current_fluxes + 0.5 * step * source2,
                    state2,
                    current[5:7] + 0.5 * step * path2,
                    eta + 0.5 * step * eta2,
                )

                source3 = self.source_terms(state3, **source_options)
                path3 = np.array([
                    math.cos(state3[E_THETA]), math.sin(state3[E_THETA])
                ])
                eta3 = self._attachment_rate(state3)[1]
                state4 = stage_state(
                    current_fluxes + step * source3,
                    state3,
                    current[5:7] + step * path3,
                    eta + step * eta3,
                )

                source4 = self.source_terms(state4, **source_options)
                path4 = np.array([
                    math.cos(state4[E_THETA]), math.sin(state4[E_THETA])
                ])
                eta4 = self._attachment_rate(state4)[1]
                flux_increment = step / 6.0 * (
                    source1 + 2.0 * source2 + 2.0 * source3 + source4
                )
                path_increment = step / 6.0 * (
                    path1 + 2.0 * path2 + 2.0 * path3 + path4
                )
                eta_increment = step / 6.0 * (
                    eta1 + 2.0 * eta2 + 2.0 * eta3 + eta4
                )
                states[index + 1] = stage_state(
                    current_fluxes + flux_increment,
                    state4,
                    current[5:7] + path_increment,
                    eta + eta_increment,
                )
                cumulative[index + 1] = cumulative[index] + flux_increment
        else:
            integration = solve_ivp(
                right_hand_side,
                (0.0, float(maximum_distance)),
                self._to_internal(initial_state),
                method=method,
                rtol=relative_tolerance,
                atol=relative_tolerance * 1.0e-3,
                max_step=maximum_step,
            )
            if not integration.success:
                raise RuntimeError(
                    f"finite wall-jet integration failed: {integration.message}"
                )
            integration_t = integration.t
            states = np.array([
                self._from_internal(column) for column in integration.y.T
            ])
        fluxes = np.array([self.fluxes(state) for state in states])
        sources = np.array([
            self.source_terms(state, **source_options) for state in states
        ])
        if cumulative is None:
            cumulative = cumulative_trapezoid(
                sources, integration_t, axis=0, initial=0.0
            )
        balances = fluxes - fluxes[0] - cumulative
        scales = np.maximum(
            np.max(np.abs(fluxes), axis=0),
            np.array([1.0e-3, 1.0e-4, 1.0, 1.0, 1.0]),
        )
        maximum_balance = float(np.max(np.abs(balances) / scales))
        temperatures = np.array([
            self.base.centre_temperature(state[:7]) for state in states
        ])
        mole_fractions = np.array([
            self.base.centre_mole_fraction(state[:7]) for state in states
        ])
        return FiniteWallJetResult(
            arc_length=integration_t,
            states=states,
            fluxes=fluxes,
            sources=sources,
            cumulative_sources=cumulative,
            balance_residuals=balances,
            temperatures=temperatures,
            centre_mole_fractions=mole_fractions,
            maximum_relative_balance_residual=maximum_balance,
        )
