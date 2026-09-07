"""Independent-energy Gaussian state for the LH2 crosswind research path.

The Fortran-compatible :class:`~degali.core.jetplume.JetPlume` carries one
thermodynamic scalar.  Density, composition, temperature and enthalpy are all
read from one mixing line, leaving four physical states for four balances.
This module supplies the five-state section representation needed before a
fifth, total-energy balance can be integrated.  It is deliberately separate
from :mod:`degali.core`: no legacy JETPLU result is changed.

The boundary projection and downstream five-balance ODE follow the frozen
sequence in ``docs/prereg-crosswind-independent-energy-state.md``. Existing
JETPLU remains available unchanged as the legacy and atmospheric baseline.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import least_squares

from ..core.constants import R_UNIVERSAL, VKC
from ..core.entrainment import phi
from ..core.jetplume import JetIntegralFluxes, ellipse

if TYPE_CHECKING:
    from .axisymmetric_jet import ConservedGaussianJet
    from ..core.jetplume import JetPlume


# Centre density, centre H2 mass fraction, width product, trajectory angle,
# excess centre velocity, downwind distance and elevation.
E_RHO, E_Y, E_SYSZ, E_THETA, E_UC, E_X, E_Z = range(7)


@dataclass(frozen=True)
class IndependentEnergyProjection:
    """Five-flux projection into the independent-energy crosswind state."""

    state: np.ndarray
    fluxes: JetIntegralFluxes
    relative_residuals: dict[str, float]
    quadrature_residual: float
    quadrature_points: int
    success: bool
    message: str


@dataclass(frozen=True)
class IndependentEnergyJetResult:
    """Integrated seven-state plume and direct conservation audit."""

    arc_length: np.ndarray
    states: np.ndarray
    fluxes: np.ndarray
    sources: np.ndarray
    cumulative_sources: np.ndarray
    balance_residuals: np.ndarray
    temperatures: np.ndarray
    centre_mole_fractions: np.ndarray
    maximum_relative_balance_residual: float


class IndependentEnergyCrosswind:
    """Five-state Gaussian cross-section with independent density and H2.

    The velocity profile and finite Gaussian domain are exactly JETPLU's.  The
    density and H2 mass-density profiles instead follow the unignited HyRAM+
    form independently.  Temperature and enthalpy are evaluated by the
    conserved near-field phase model from local density and composition.
    """

    def __init__(
        self,
        jetplume: "JetPlume",
        thermodynamics: "ConservedGaussianJet",
        *,
        quadrature_points: int = 32,
        energy_transport: str = "total",
        houf_width_mapping: str = "velocity",
        ground_interaction: str = "free",
        thermal_relaxation_rate: float = 0.0,
        phase_transition_lag_rate: float = 0.0,
        turbulence_heat_exchange_rate: float = 0.0,
    ):
        if quadrature_points < 8:
            raise ValueError("crosswind quadrature requires at least eight points")
        if energy_transport not in {"total", "enthalpy"}:
            raise ValueError("energy transport must be 'total' or 'enthalpy'")
        if houf_width_mapping not in {"velocity", "scalar"}:
            raise ValueError("Houf width mapping must be 'velocity' or 'scalar'")
        if ground_interaction not in {"free", "geometry", "surface_layer"}:
            raise ValueError(
                "ground interaction must be 'free', 'geometry' or "
                "'surface_layer'"
            )
        if not np.isfinite(thermal_relaxation_rate) or thermal_relaxation_rate < 0.0:
            raise ValueError(
                "thermal-relaxation rate must be a finite non-negative number"
            )
        if not np.isfinite(phase_transition_lag_rate) or phase_transition_lag_rate < 0.0:
            raise ValueError(
                "phase-transition lag rate must be a finite non-negative number"
            )
        if not np.isfinite(turbulence_heat_exchange_rate) or turbulence_heat_exchange_rate < 0.0:
            raise ValueError(
                "turbulence heat exchange rate must be a finite non-negative "
                "number"
            )
        self.jetplume = jetplume
        self.thermodynamics = thermodynamics
        self.th = jetplume.th
        self.k = jetplume.k
        self.rhoa = float(thermodynamics.ambient_density)
        self.velocity_shape_exponent = thermodynamics.spreading_ratio**2
        self.quadrature_points = int(quadrature_points)
        self.energy_transport = energy_transport
        self.houf_width_mapping = houf_width_mapping
        self.ground_interaction = ground_interaction
        self.thermal_relaxation_rate = float(thermal_relaxation_rate)
        self.phase_transition_lag_rate = float(phase_transition_lag_rate)
        self.turbulence_heat_exchange_rate = float(
            turbulence_heat_exchange_rate
        )
        self._quadrature_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def _quadrature(self, points: int) -> tuple[np.ndarray, np.ndarray]:
        """Return flattened Gaussian exponent and unit-width area weights."""
        if points not in self._quadrature_cache:
            nodes, weights = np.polynomial.legendre.leggauss(points)
            coordinate = self.k.delta * nodes
            grid_y, grid_z = np.meshgrid(
                coordinate, coordinate, indexing="ij"
            )
            weight_y, weight_z = np.meshgrid(weights, weights, indexing="ij")
            exponent = math.pi / 8.0 * (grid_y**2 + grid_z**2)
            unit_weights = (
                math.pi / 4.0 * self.k.delta**2 * weight_y * weight_z
            )
            # Profiles depend only on this exponent. Merge EXACTLY equal
            # nodes (no rounding/coarsening) using reflection/transpose
            # symmetry, retaining the original tensor-rule weights.
            unique, inverse = np.unique(exponent.ravel(), return_inverse=True)
            self._quadrature_cache[points] = (
                unique, np.bincount(inverse, weights=unit_weights.ravel())
            )
        return self._quadrature_cache[points]

    def _mean_section_quantities(
        self,
        state: np.ndarray,
        *,
        use_condensed_air: bool = False,
    ) -> tuple[float, float, float, float]:
        """Return section-mean (rel. enthalpy, temperature, H2 fraction, cp).

        The mean values are mass-flux weighted and use the current
        thermodynamic branch: either condensed-equilibrium or frozen ideal-gas
        closure.
        """
        rho_c, fraction_c, sysz, theta, uc, _x, _z = self._physical(state)
        exponent, unit_weights = self._quadrature(self.quadrature_points)
        scalar_shape = np.exp(-exponent)
        velocity_shape = np.exp(-self.velocity_shape_exponent * exponent)
        density = self.rhoa + (rho_c - self.rhoa) * scalar_shape
        fuel_mass_density = rho_c * fraction_c * scalar_shape
        fraction = np.minimum(
            np.maximum(fuel_mass_density / density, 1.0e-14), 1.0
        )
        velocity = self._wind(state) * math.cos(theta) + uc * velocity_shape
        if use_condensed_air and self.thermodynamics.equilibrium_air_condensation:
            temperature, rho_h = self.thermodynamics._condensed_air_state(
                density, fraction
            )
        else:
            temperature = self.thermodynamics._temperature_from_density(
                density, fraction
            )
            rho_h = density * self.thermodynamics._mixture_enthalpy(
                temperature, fraction
            )
        weights = sysz * unit_weights
        mass_flux = float(np.sum(density * velocity * weights))
        if mass_flux <= 0.0:
            return 0.0, float(self.thermodynamics.ambient_temperature), fraction_c, 0.0
        relative_enthalpy_flux = float(np.sum(
            (rho_h - self.thermodynamics._ambient_enthalpy * density)
            * velocity * weights
        ))
        temperature_flux = float(np.sum(
            temperature * density * velocity * weights
        ))
        fraction_flux = float(np.sum(
            fraction * density * velocity * weights
        ))
        mean_relative_enthalpy = relative_enthalpy_flux / mass_flux
        mean_temperature = temperature_flux / mass_flux
        mean_fraction = fraction_flux / mass_flux
        mean_cp = self._mixture_heat_capacity(mean_temperature, mean_fraction)
        return (
            mean_relative_enthalpy,
            mean_temperature,
            mean_fraction,
            mean_cp,
        )

    def _mixture_heat_capacity(self, temperature: float, fraction: float) -> float:
        """Return mixture heat capacity by local enthalpy slope."""
        try:
            delta = max(
                1.0e-3,
                1.0e-4 * max(abs(float(temperature)), 1.0),
            )
            above = float(self.thermodynamics._mixture_enthalpy(
                float(temperature) + delta, float(fraction)
            ))
            below = float(self.thermodynamics._mixture_enthalpy(
                float(temperature) - delta, float(fraction)
            ))
            heat_capacity = (above - below) / (2.0 * delta)
        except (FloatingPointError, OverflowError):
            heat_capacity = 0.0
        if not np.isfinite(heat_capacity) or heat_capacity <= 0.0:
            return 1000.0
        return float(heat_capacity)

    @staticmethod
    def _physical(state: np.ndarray) -> tuple[float, ...]:
        if np.shape(state) != (7,) or not np.all(np.isfinite(state)):
            raise ValueError("independent-energy state must contain seven finite values")
        rho_c, fraction_c, sysz, theta, uc, x, z = map(float, state)
        if min(rho_c, sysz, uc) <= 0.0:
            raise ValueError("density, width product and excess velocity must be positive")
        if not 0.0 < fraction_c <= 1.0:
            raise ValueError("centre H2 mass fraction must lie in (0, 1]")
        if not -math.pi / 2.0 < theta < math.pi / 2.0:
            raise ValueError("crosswind trajectory angle must lie within +/- pi/2")
        return rho_c, fraction_c, sysz, theta, uc, x, z

    def _wind(self, state: np.ndarray) -> float:
        """JETPLU cross-section-averaged wind for one seven-state section."""
        _rho_c, _fraction_c, sysz, theta, _uc, x, z = self._physical(state)
        sya = self.jetplume.deltay * max(x, 0.0) ** self.jetplume.betay
        if x > 0.0:
            sza = (
                self.jetplume.deltaz * x**self.jetplume.betaz
                * math.exp(self.jetplume.gammaz * math.log(x) ** 2)
            )
        else:
            sza = 0.0
        _sy, sz = self.jetplume._split(sysz, sya, sza)
        return float(self.jetplume._wind(z, sz, math.cos(theta)))

    def _geometry(
        self, state: np.ndarray
    ) -> tuple[float, ...]:
        """Return widths, wind, contact geometry and ambient spreading."""
        _rho_c, _fraction_c, sysz, theta, _uc, x, z = self._physical(state)
        ct = math.cos(theta)
        sya = self.jetplume.deltay * max(x, 0.0) ** self.jetplume.betay
        if x > 0.0:
            sza = (
                self.jetplume.deltaz * x**self.jetplume.betaz
                * math.exp(self.jetplume.gammaz * math.log(x) ** 2)
            )
        else:
            sza = 0.0
        dsya = dsza = 0.0
        if x > 1.0:
            dsya = sya / x * self.jetplume.betay * ct
            dsza = (
                sza / x
                * (
                    self.jetplume.betaz
                    + 2.0 * self.jetplume.gammaz * math.log(x)
                )
                * ct
            )
        sy, sz = self.jetplume._split(sysz, sya, sza)
        half_depth = self.k.delta * sz * ct
        ground_width = 0.0
        cleared_fraction = 1.0
        if self.ground_interaction != "free" and z < half_depth:
            segment = max(
                min(z / max(ct, 1.0e-6), self.k.delta * sz),
                -self.k.delta * sz,
            )
            perimeter, _area = ellipse(
                sy * self.k.delta, sz * self.k.delta, segment
            )
            ratio = segment / max(self.k.delta * sz, 1.0e-30)
            ground_width = (
                2.0 * self.k.delta * sy
                * math.sqrt(max(1.0 - ratio * ratio, 0.0))
            )
            cleared_fraction = min(
                max(z / max(half_depth, 1.0e-9), 0.0), 1.0
            )
        else:
            perimeter, _area = ellipse(
                sy * self.k.delta, sz * self.k.delta, sz * self.k.delta
            )
        ua, alpha_wind = self.jetplume._wind_profile(z, sz, ct)
        return (
            sy, sz, float(ua), perimeter, ground_width, cleared_fraction,
            float(alpha_wind), sya, sza, dsya, dsza,
        )

    def profiles(
        self, state: np.ndarray, *, quadrature_points: int | None = None
    ) -> tuple[np.ndarray, ...]:
        """Return velocity, density, H2 fraction, temperature and ``rho*h``."""
        rho_c, fraction_c, _sysz, theta, uc, _x, _z = self._physical(state)
        points = self.quadrature_points if quadrature_points is None else int(
            quadrature_points
        )
        exponent, _weights = self._quadrature(points)
        scalar_shape = np.exp(-exponent)
        velocity_shape = np.exp(-self.velocity_shape_exponent * exponent)
        density = self.rhoa + (rho_c - self.rhoa) * scalar_shape
        fuel_mass_density = rho_c * fraction_c * scalar_shape
        fraction = fuel_mass_density / density
        if np.any(fraction <= 0.0) or np.any(fraction > 1.0 + 1.0e-10):
            raise ValueError("independent Gaussian profiles produced invalid composition")
        fraction = np.minimum(fraction, 1.0)
        velocity = self._wind(state) * math.cos(theta) + uc * velocity_shape
        if np.any(velocity <= 0.0):
            raise ValueError("independent Gaussian profile has reverse axial flow")
        if self.thermodynamics.equilibrium_air_condensation:
            temperature, rho_h = self.thermodynamics._condensed_air_state(
                density, fraction
            )
        else:
            temperature = self.thermodynamics._temperature_from_density(
                density, fraction
            )
            rho_h = density * self.thermodynamics._mixture_enthalpy(
                temperature, fraction
            )
        return velocity, density, fraction, temperature, rho_h

    def integral_fluxes(
        self, state: np.ndarray, *, quadrature_points: int | None = None
    ) -> JetIntegralFluxes:
        """Integrate all five conserved transports represented by ``state``."""
        _rho_c, _fraction_c, sysz, theta, _uc, _x, _z = self._physical(state)
        points = self.quadrature_points if quadrature_points is None else int(
            quadrature_points
        )
        _exponent, unit_weights = self._quadrature(points)
        velocity, density, fraction, _temperature, rho_h = self.profiles(
            state, quadrature_points=points
        )
        weights = sysz * unit_weights
        mass = float(np.sum(density * velocity * weights))
        species = float(np.sum(density * fraction * velocity * weights))
        axial_momentum = float(np.sum(density * velocity**2 * weights))
        energy_density = (
            rho_h - self.thermodynamics._ambient_enthalpy * density
        ) * velocity
        if self.energy_transport == "total":
            energy_density = energy_density + 0.5 * density * velocity**3
        energy = float(np.sum(energy_density * weights))
        return JetIntegralFluxes(
            total_mass=mass,
            contaminant_mass=species,
            momentum_x=axial_momentum * math.cos(theta),
            momentum_z=axial_momentum * math.sin(theta),
            energy=energy,
        )

    def centre_temperature(self, state: np.ndarray) -> float:
        """Temperature selected by the phase closure at the centre state."""
        rho_c, fraction_c, *_rest = self._physical(state)
        if self.thermodynamics.equilibrium_air_condensation:
            temperature, _rho_h = self.thermodynamics._condensed_air_state(
                np.array([rho_c]), np.array([fraction_c])
            )
            return float(temperature[0])
        return float(
            self.thermodynamics._temperature_from_density(rho_c, fraction_c)
        )

    def centre_mole_fraction(self, state: np.ndarray) -> float:
        """Convert the independent centre mass fraction to mole fraction."""
        _rho_c, fraction_c, *_rest = self._physical(state)
        numerator = fraction_c / self.thermodynamics.fuel_molecular_weight
        denominator = numerator + (
            (1.0 - fraction_c)
            / self.thermodynamics._humid_ambient_molecular_weight
        )
        return float(numerator / denominator)

    def buoyancy_force(self, state: np.ndarray) -> float:
        """Return the signed section-integrated vertical buoyancy force."""
        rho_c, _fraction_c, sysz, *_rest = self._physical(state)
        return float(
            9.81 * (self.rhoa - rho_c) * sysz
            * self.k.profile_integral(1.0)
        )

    def _buoyancy_work(self, state: np.ndarray, buoyancy: float) -> float:
        """Return a conservative reduced buoyancy-work contribution."""
        rho_c, _fraction_c, sysz, theta, uc, _x, _z = self._physical(state)
        if buoyancy == 0.0:
            return 0.0
        exponent, unit_weights = self._quadrature(self.quadrature_points)
        scalar_shape = np.exp(-exponent)
        velocity_shape = np.exp(-self.velocity_shape_exponent * exponent)
        density = self.rhoa + (rho_c - self.rhoa) * scalar_shape
        velocity = self._wind(state) * math.cos(theta) + uc * velocity_shape
        return float(
            9.81 * math.sin(theta)
            * np.sum((self.rhoa - density) * velocity * (sysz * unit_weights))
        )

    def _mean_relative_specific_enthalpy(self, state: np.ndarray) -> float:
        """Return section-mean specific enthalpy relative to ambient."""
        mean_relative_enthalpy, *_rest = self._mean_section_quantities(
            state, use_condensed_air=self.thermodynamics.equilibrium_air_condensation
        )
        return float(mean_relative_enthalpy)

    def _mean_relative_specific_enthalpy_frozen(self, state: np.ndarray) -> float:
        """Return section-mean specific enthalpy as if no condensation occurs."""
        mean_relative_enthalpy, *_rest = self._mean_section_quantities(
            state, use_condensed_air=False
        )
        return float(mean_relative_enthalpy)

    def _mean_temperature(self, state: np.ndarray) -> float:
        """Return section-mean temperature."""
        _mean_relative_enthalpy, mean_temperature, _mean_fraction, _mean_cp = (
            self._mean_section_quantities(
                state,
                use_condensed_air=self.thermodynamics.equilibrium_air_condensation,
            )
        )
        return float(mean_temperature)

    def _thermal_relaxation_source(self, state: np.ndarray, mass_source: float) -> float:
        """Return ambient-relaxation energy contribution."""
        if self.thermal_relaxation_rate <= 0.0 or mass_source <= 0.0:
            return 0.0
        return -self.thermal_relaxation_rate * mass_source * (
            self._mean_relative_specific_enthalpy(state)
        )

    def _phase_transition_relaxation_source(
        self, state: np.ndarray, mass_source: float
    ) -> float:
        """Return finite-rate phase-transition relaxation correction."""
        if (
            self.phase_transition_lag_rate <= 0.0
            or mass_source <= 0.0
            or not self.thermodynamics.equilibrium_air_condensation
        ):
            return 0.0
        mean_frozen, frozen_temperature, _, _ = self._mean_section_quantities(
            state, use_condensed_air=False
        )
        mean_equilibrium, equilibrium_temperature, *_ = self._mean_section_quantities(
            state, use_condensed_air=self.thermodynamics.equilibrium_air_condensation
        )
        phase_gap = mean_frozen - mean_equilibrium
        # Condensation strength is linked to physical subcooling, while this
        # rate remains a user-selected relaxation time scale.
        subcooling = max(
            0.0,
            (self.thermodynamics.ambient_temperature
             - equilibrium_temperature)
        )
        activity = subcooling / (2.0 * self.thermodynamics.ambient_temperature)
        activity = np.clip(activity, 0.0, 1.0)
        return (
            self.phase_transition_lag_rate
            * mass_source * phase_gap * activity
        )

    def _turbulent_heat_exchange_source(
        self, state: np.ndarray, mass_source: float
    ) -> float:
        """Return turbulent heat exchange across the plume perimeter."""
        if (
            self.turbulence_heat_exchange_rate <= 0.0
            or mass_source <= 0.0
        ):
            return 0.0
        _, sysz, _ua, perimeter, *_rest = self._geometry(state)
        if perimeter <= 0.0 or sysz <= 0.0:
            return 0.0
        mean_relative_enthalpy, mean_temperature, _, mean_cp = (
            self._mean_section_quantities(
                state,
                use_condensed_air=self.thermodynamics.equilibrium_air_condensation,
            )
        )
        # Perimeter-to-width ratio is the smallest geometry-invariant measure of
        # turbulent transport area in a quasi-elliptic section.
        transfer_scale = perimeter / max(4.0 * math.sqrt(sysz), 1.0e-6)
        thermal_deficit = mean_relative_enthalpy
        if mean_cp > 0.0:
            thermal_deficit = mean_cp * (
                mean_temperature - self.thermodynamics.ambient_temperature
            )
        return -self.turbulence_heat_exchange_rate * mass_source * (
            thermal_deficit * transfer_scale
        )

    def houf_velocity_width(self, state: np.ndarray) -> float:
        """Return HyRAM's velocity e-folding width for an elliptic section.

        ``sysz`` stores the product of scalar Gaussian standard deviations,
        while HyRAM's ``B`` is the velocity-profile e-folding width. The
        velocity profile is narrower by the spreading ratio ``lambda``.
        ``scalar`` retains the superseded 2026-09-05 mapping only so that the
        original research result remains reproducible.
        """
        _rho_c, _fraction_c, sysz, *_rest = self._physical(state)
        width = math.sqrt(max(2.0 * sysz, 1.0e-30))
        if self.houf_width_mapping == "velocity":
            width /= self.thermodynamics.spreading_ratio
        return float(width)

    def section_widths(self, state: np.ndarray) -> tuple[float, float]:
        """Return JETPLU-compatible lateral and vertical Gaussian sigmas."""
        sy, sz, *_rest = self._geometry(state)
        return float(sy), float(sz)

    def point_mole_fraction(
        self, state: np.ndarray, lateral: float, height: float
    ) -> float:
        """H2 mole fraction at a receptor, including the ground image."""
        rho_c, fraction_c, _sysz, _theta, _uc, _x, centre_z = self._physical(
            state
        )
        sy, sz = self.section_widths(state)
        lateral_shape = math.exp(-0.5 * (float(lateral) / sy) ** 2)
        direct = math.exp(-0.5 * ((float(height) - centre_z) / sz) ** 2)
        image = math.exp(-0.5 * ((float(height) + centre_z) / sz) ** 2)
        shape = lateral_shape * (direct + image)
        density = self.rhoa + (rho_c - self.rhoa) * shape
        fuel_density = rho_c * fraction_c * shape
        if density <= 0.0 or fuel_density < 0.0 or fuel_density > density:
            raise RuntimeError("ground image produced a non-physical local mixture")
        fraction = fuel_density / density
        numerator = fraction / self.thermodynamics.fuel_molecular_weight
        denominator = numerator + (
            (1.0 - fraction)
            / self.thermodynamics._humid_ambient_molecular_weight
        )
        return float(numerator / denominator)

    def point_temperature(
        self, state: np.ndarray, lateral: float, height: float
    ) -> float:
        """Local temperature at a receptor, including the ground image."""
        rho_c, fraction_c, _sysz, _theta, _uc, _x, centre_z = self._physical(
            state
        )
        sy, sz = self.section_widths(state)
        lateral_shape = math.exp(-0.5 * (float(lateral) / sy) ** 2)
        direct = math.exp(-0.5 * ((float(height) - centre_z) / sz) ** 2)
        image = math.exp(-0.5 * ((float(height) + centre_z) / sz) ** 2)
        shape = lateral_shape * (direct + image)
        density = self.rhoa + (rho_c - self.rhoa) * shape
        fuel_density = rho_c * fraction_c * shape
        if density <= 0.0 or fuel_density < 0.0 or fuel_density > density:
            raise RuntimeError("ground image produced a non-physical local mixture")
        fraction = fuel_density / density
        if self.thermodynamics.equilibrium_air_condensation:
            temperature, _rho_h = self.thermodynamics._condensed_air_state(
                np.array([density]), np.array([fraction])
            )
            return float(temperature[0])
        return float(
            self.thermodynamics._temperature_from_density(density, fraction)
        )

    @staticmethod
    def _as_array(fluxes: JetIntegralFluxes) -> np.ndarray:
        return np.array([
            fluxes.total_mass,
            fluxes.contaminant_mass,
            fluxes.momentum_x,
            fluxes.momentum_z,
            fluxes.energy,
        ])

    def project(
        self,
        target: JetIntegralFluxes,
        initial_state: np.ndarray,
        *,
        relative_tolerance: float = 1.0e-8,
        quadrature_tolerance: float = 1.0e-5,
        maximum_quadrature_points: int = 512,
    ) -> IndependentEnergyProjection:
        """Solve the five centre states against five target fluxes."""
        initial = np.asarray(initial_state, dtype=float)
        rho0, fraction0, sysz0, theta0, uc0, x, z = self._physical(initial)
        target_values = self._as_array(target)
        scales = np.maximum(np.abs(target_values), np.array([
            1.0e-12, 1.0e-12, 1.0, 1.0, 1.0,
        ]))

        def encode(state: np.ndarray) -> np.ndarray:
            rho_c, fraction_c, sysz, theta, uc = state[:5]
            bounded_fraction = min(max(fraction_c, 1.0e-12), 1.0 - 1.0e-12)
            return np.array([
                math.log(rho_c),
                math.log(bounded_fraction / (1.0 - bounded_fraction)),
                math.log(sysz),
                theta,
                math.log(uc),
            ])

        def decode(parameters: np.ndarray) -> np.ndarray:
            logistic = float(np.clip(parameters[1], -35.0, 35.0))
            fraction_c = 1.0 / (1.0 + math.exp(-logistic))
            return np.array([
                math.exp(float(parameters[0])),
                fraction_c,
                math.exp(float(parameters[2])),
                float(parameters[3]),
                math.exp(float(parameters[4])),
                x,
                z,
            ])

        lower = np.array([
            math.log(1.0e-3), -25.0, math.log(1.0e-14),
            -math.pi / 2.0 + 1.0e-8, math.log(1.0e-8),
        ])
        upper = np.array([
            math.log(100.0), 25.0, math.log(100.0),
            math.pi / 2.0 - 1.0e-8, math.log(1.0e4),
        ])
        parameters = encode(np.array([
            rho0, fraction0, sysz0, theta0, uc0, x, z,
        ]))
        points = self.quadrature_points
        quadrature_residual = math.inf
        fit = None
        state = initial.copy()
        final_refit = False

        while True:
            def residual(values: np.ndarray) -> np.ndarray:
                try:
                    represented = self._as_array(self.integral_fluxes(
                        decode(values), quadrature_points=points
                    ))
                except (RuntimeError, ValueError, FloatingPointError):
                    return np.full(5, 1.0e6)
                if not np.all(np.isfinite(represented)):
                    return np.full(5, 1.0e6)
                return (represented - target_values) / scales

            parameters = np.minimum(
                np.maximum(parameters, lower + 1.0e-12), upper - 1.0e-12
            )
            fit = least_squares(
                residual, parameters, bounds=(lower, upper),
                xtol=1.0e-12, ftol=1.0e-12, gtol=1.0e-12,
                max_nfev=1000,
            )
            parameters = fit.x
            state = decode(parameters)
            coarse = self.integral_fluxes(state, quadrature_points=points)
            if final_refit:
                break
            if points >= maximum_quadrature_points:
                break
            fine_points = min(2 * points, maximum_quadrature_points)
            fine = self.integral_fluxes(state, quadrature_points=fine_points)
            quadrature_residual = self.quadrature_error(coarse, fine)
            points = fine_points
            if quadrature_residual <= quadrature_tolerance:
                # Refit at the accepted finer order so the reported physical
                # flux residuals, rather than only the quadrature comparison,
                # meet the frozen gate.
                final_refit = True
                continue
            if points >= maximum_quadrature_points:
                break
            continue

        fluxes = self.integral_fluxes(state, quadrature_points=points)
        residual_values = (self._as_array(fluxes) - target_values) / scales
        names = ("total_mass", "hydrogen", "momentum_x", "momentum_z", "energy")
        relative_residuals = dict(zip(names, map(float, np.abs(residual_values))))
        success = bool(
            fit is not None
            and fit.success
            and max(relative_residuals.values()) <= relative_tolerance
            and quadrature_residual <= quadrature_tolerance
        )
        message = "projection converged" if success else (
            "five-flux projection or quadrature convergence failed"
            if fit is None else str(fit.message)
        )
        return IndependentEnergyProjection(
            state=state,
            fluxes=fluxes,
            relative_residuals=relative_residuals,
            quadrature_residual=float(quadrature_residual),
            quadrature_points=points,
            success=success,
            message=message,
        )

    def quadrature_error(self, coarse: JetIntegralFluxes, fine: JetIntegralFluxes) -> float:
        """Legacy density-Gaussian mass/momentum are analytic; audit energy."""
        return abs(fine.energy-coarse.energy)/max(abs(fine.energy), 1.)

    def source_terms(
        self,
        state: np.ndarray,
        *,
        include_work: bool = False,
        include_thermal_relaxation: bool = False,
        include_phase_transition: bool = False,
        include_turbulent_heat_exchange: bool = False,
    ) -> np.ndarray:
        """Mass, species, vector-momentum and energy sources per metre.

        `include_work` enables a reduced buoyancy-work closure for research
        calculations.  It is disabled by default to preserve the original
        default behavior.  `include_thermal_relaxation` adds an ambient
        relaxation term for thermal equilibration.
        `include_phase_transition` adds finite-rate phase-equilibrium relaxation.
        `include_turbulent_heat_exchange` adds a boundary-like turbulent heat
        exchange term.
        """
        rho_c, _fraction_c, sysz, theta, uc, _x, _z = self._physical(state)
        (
            _sy, sz, ua, perimeter, ground_width, cleared_fraction,
            alpha_wind, sya, sza, dsya, dsza,
        ) = self._geometry(state)
        st, ct = math.sin(theta), math.cos(theta)
        uentr = max(uc, 0.0)

        richardson = (
            9.81 * abs(rho_c - self.rhoa) * math.sqrt(sysz)
            / max(self.rhoa * max(uentr, 1.0e-8) ** 2, 1.0e-12)
        )
        alpha = self.k.alfa1
        if self.k.plume_transition:
            ratio = min(abs(richardson) / self.k.RI_PLUME, 1.0)
            alpha = self.k.ALPHA_JET + (
                self.k.ALPHA_PLUME - self.k.ALPHA_JET
            ) * ratio**2
        # The seven-trial mechanism-isolation test selects Ricou--Spalding
        # local density scaling after 10D. Apply it independently of a failed
        # single-scalar projection having restored the legacy switch.
        alpha /= math.sqrt(max(rho_c / self.rhoa, 1.0e-12))
        shear_entrainment = alpha * uentr * perimeter

        width = self.houf_velocity_width(state)
        centre_velocity = max(uc + ua * ct, 1.0e-9)
        contrast = abs(self.rhoa - rho_c)
        buoyant_entrainment = 0.0
        if contrast > 1.0e-12 and st > 0.0:
            local_froude = (
                centre_velocity**2 * max(rho_c, 1.0e-12)
                / (9.81 * width * contrast)
            )
            buoyant_entrainment = (
                self.thermodynamics._buoyancy_coefficient
                / max(local_froude, 1.0e-30)
                * 2.0 * math.pi * centre_velocity * width * st
            )
        cap = (
            self.thermodynamics.plume_entrainment_limit
            * 2.0 * math.pi * width * centre_velocity
        )
        local_entrainment = min(
            shear_entrainment + buoyant_entrainment, cap
        )
        crossflow_entrainment = (
            self.k.alfa2 * ua * ct * abs(st) * perimeter
        )
        ambient_spread_entrainment = (
            self.k.profile_integral(0.0)
            * ua * (sza * dsya + sya * dsza)
        )
        ground_entrainment = 0.0
        if self.ground_interaction == "surface_layer" and ground_width > 0.0:
            half_depth = self.k.delta * sz * ct
            layer_depth = max(state[E_Z] + half_depth, 1.0e-9)
            ri_layer = (
                9.81 * max(rho_c - self.rhoa, 0.0) * layer_depth
                / (
                    self.rhoa
                    * max(self.jetplume.ustar, 1.0e-6) ** 2
                    * self.k.delta
                )
            )
            layer_velocity = (
                VKC * self.jetplume.ustar * (1.0 + alpha_wind)
                / phi(ri_layer, 0.0, 3)
            )
            ground_entrainment = layer_velocity * ground_width
        entrainment = (
            local_entrainment
            + crossflow_entrainment
            + ambient_spread_entrainment
            + ground_entrainment
        )
        mass_source = self.rhoa * entrainment

        drag = self.k.cd * perimeter * self.rhoa / 2.0 * (ua * st) ** 2
        buoyancy = self.buoyancy_force(state)
        momentum_x_source = ua * mass_source + drag * abs(st)
        momentum_z_source = (
            buoyancy - math.copysign(1.0, theta) * drag * ct
        )
        if (
            self.ground_interaction != "free"
            and cleared_fraction < 1.0
            and momentum_z_source > 0.0
        ):
            momentum_z_source *= cleared_fraction
        energy_source = (
            0.5 * mass_source * ua**2
            if self.energy_transport == "total" else 0.0
        )
        if include_work and self.energy_transport == "total":
            energy_source += self._buoyancy_work(state, buoyancy)
        if include_thermal_relaxation:
            energy_source += self._thermal_relaxation_source(state, mass_source)
        if include_phase_transition:
            energy_source += self._phase_transition_relaxation_source(
                state, mass_source
            )
        if include_turbulent_heat_exchange:
            energy_source += self._turbulent_heat_exchange_source(
                state, mass_source
            )
        return np.array([
            mass_source,
            0.0,
            momentum_x_source,
            momentum_z_source,
            energy_source,
        ])

    def _match_flux_array(
        self,
        target_values: np.ndarray,
        initial_state: np.ndarray,
        *,
        relative_tolerance: float = 1.0e-9,
    ) -> np.ndarray:
        """Recover five physical states from nearby conserved fluxes."""
        initial_state = np.asarray(initial_state, dtype=float)
        _rho0, _fraction0, sysz0, _theta0, uc0, x, z = self._physical(
            initial_state
        )
        target_values = np.asarray(target_values, dtype=float)
        target_mass = float(target_values[0])
        target_species = float(target_values[1])
        momentum_target = math.hypot(target_values[2], target_values[3])
        theta_target = math.atan2(target_values[3], target_values[2])
        target_energy = float(target_values[4])
        if min(target_mass, target_species, momentum_target) <= 0.0:
            raise RuntimeError("independent-energy target fluxes must be positive")

        i0 = self.k.profile_integral(0.0)
        i1 = self.k.profile_integral(1.0)
        iu = self.k.profile_integral(self.velocity_shape_exponent)
        i1u = self.k.profile_integral(1.0 + self.velocity_shape_exponent)
        i2u = self.k.profile_integral(2.0 * self.velocity_shape_exponent)
        i1_2u = self.k.profile_integral(
            1.0 + 2.0 * self.velocity_shape_exponent
        )

        def decode(parameters: np.ndarray) -> np.ndarray:
            sysz = math.exp(float(parameters[0]))
            uc = math.exp(float(parameters[1]))
            wind_state = np.array([
                self.rhoa, 0.1, sysz, theta_target, uc, x, z,
            ])
            ambient_axial = self._wind(wind_state) * math.cos(theta_target)
            mass_ambient = self.rhoa * (ambient_axial * i0 + uc * iu)
            mass_density_shape = ambient_axial * i1 + uc * i1u
            density_excess = (
                target_mass / sysz - mass_ambient
            ) / mass_density_shape
            rho_c = self.rhoa + density_excess
            fraction_c = target_species / (
                sysz * rho_c * mass_density_shape
            )
            return np.array([
                rho_c, fraction_c, sysz, theta_target, uc, x, z,
            ])

        parameters = np.log([sysz0, uc0])
        lower = np.array([
            math.log(1.0e-14), math.log(1.0e-8),
        ])
        upper = np.array([
            math.log(100.0), math.log(1.0e4),
        ])

        def residual(values: np.ndarray) -> np.ndarray:
            try:
                state = decode(values)
                rho_c, fraction_c, sysz, _theta, uc, _x, _z = (
                    self._physical(state)
                )
                ambient_axial = self._wind(state) * math.cos(theta_target)
                density_excess = rho_c - self.rhoa
                momentum_ambient = self.rhoa * (
                    ambient_axial**2 * i0
                    + 2.0 * ambient_axial * uc * iu
                    + uc**2 * i2u
                )
                momentum_density = density_excess * (
                    ambient_axial**2 * i1
                    + 2.0 * ambient_axial * uc * i1u
                    + uc**2 * i1_2u
                )
                represented_momentum = sysz * (
                    momentum_ambient + momentum_density
                )
                represented_energy = self.integral_fluxes(state).energy
            except (RuntimeError, ValueError, FloatingPointError):
                return np.full(2, 1.0e6)
            return np.array([
                (represented_momentum - momentum_target)
                / max(abs(momentum_target), 1.0),
                (represented_energy - target_energy)
                / max(abs(target_energy), 1.0),
            ])

        scales = np.maximum(np.abs(target_values), np.array([
            1.0e-12, 1.0e-12, 1.0, 1.0, 1.0,
        ]))
        clipped = np.minimum(
            np.maximum(parameters, lower + 1.0e-12), upper - 1.0e-12
        )

        def fit_from(start: np.ndarray):
            fit = least_squares(
                residual,
                np.minimum(
                    np.maximum(start, lower + 1.0e-12), upper - 1.0e-12
                ),
                bounds=(lower, upper),
                xtol=1.0e-11,
                ftol=1.0e-11,
                gtol=1.0e-11,
                max_nfev=200,
            )
            try:
                state = decode(fit.x)
                represented = self._as_array(self.integral_fluxes(state))
                error = (represented - target_values) / scales
                score = float(np.max(np.abs(error)))
            except (RuntimeError, ValueError, FloatingPointError):
                state = np.full(7, np.nan)
                error = np.full(5, np.inf)
                score = math.inf
            return fit, state, error, score

        fit, state, error, score = fit_from(clipped)
        # A nearby RK flux target can occasionally leave the two-variable
        # inversion on a flat local branch.  Retry only after that strict
        # solve fails, using symmetric factor-of-two starts in physical
        # width-area and excess velocity.  This changes neither equations nor
        # acceptance tolerance and the lowest full five-flux residual wins.
        if not fit.success or score > relative_tolerance:
            offset = math.log(2.0)
            for delta_width in (-offset, 0.0, offset):
                for delta_velocity in (-offset, 0.0, offset):
                    if delta_width == 0.0 and delta_velocity == 0.0:
                        continue
                    candidate = fit_from(
                        clipped + np.array([delta_width, delta_velocity])
                    )
                    candidate_fit, _state, _error, candidate_score = candidate
                    if (
                        candidate_score < score
                        or (
                            candidate_score == score
                            and candidate_fit.success and not fit.success
                        )
                    ):
                        fit, state, error, score = candidate
        if (
            not fit.success
            or not np.all(np.isfinite(state))
            or score > relative_tolerance
        ):
            raise RuntimeError(
                "independent-energy flux-state inversion failed: "
                f"{fit.message}; maximum residual {score:.3e}"
            )
        return state

    def flux_jacobian(self, state: np.ndarray) -> np.ndarray:
        """Finite-difference Jacobian of five fluxes by five centre states."""
        state = np.asarray(state, dtype=float)
        self._physical(state)
        jacobian = np.empty((5, 5))
        floors = np.array([1.0, 0.1, 1.0e-3, 1.0, 1.0])
        base = self._as_array(self.integral_fluxes(state))
        for column in range(5):
            step = 2.0e-6 * max(abs(state[column]), floors[column])
            plus, minus = state.copy(), state.copy()
            plus[column] += step
            minus[column] -= step
            one_sided = (
                column in (E_RHO, E_Y, E_SYSZ, E_UC)
                and minus[column] <= 0.0
            )
            if column == E_Y and plus[column] >= 1.0:
                jacobian[:, column] = (
                    base - self._as_array(self.integral_fluxes(minus))
                ) / step
            elif one_sided:
                jacobian[:, column] = (
                    self._as_array(self.integral_fluxes(plus)) - base
                ) / step
            else:
                jacobian[:, column] = (
                    self._as_array(self.integral_fluxes(plus))
                    - self._as_array(self.integral_fluxes(minus))
                ) / (2.0 * step)
        return jacobian

    def derivatives(
        self,
        _s: float,
        state: np.ndarray,
        *,
        include_work: bool = False,
        include_thermal_relaxation: bool = False,
        include_phase_transition: bool = False,
        include_turbulent_heat_exchange: bool = False,
    ) -> np.ndarray:
        """Solve the five conservative balances and two trajectory equations."""
        state = np.asarray(state, dtype=float)
        self._physical(state)
        fluxes = self._as_array(self.integral_fluxes(state))
        jacobian = self.flux_jacobian(state)
        sources = self.source_terms(
            state,
            include_work=include_work,
            include_thermal_relaxation=include_thermal_relaxation,
            include_phase_transition=include_phase_transition,
            include_turbulent_heat_exchange=include_turbulent_heat_exchange,
        )
        column_scales = np.maximum(
            np.abs(state[:5]), np.array([1.0, 0.1, 1.0e-3, 1.0, 1.0])
        )
        row_scales = np.maximum(
            np.abs(fluxes), np.array([1.0e-3, 1.0e-4, 1.0, 1.0, 1.0])
        )
        scaled_jacobian = (
            jacobian * column_scales[np.newaxis, :]
            / row_scales[:, np.newaxis]
        )
        condition = float(np.linalg.cond(scaled_jacobian))
        if not math.isfinite(condition) or condition > 1.0e12:
            raise RuntimeError(
                f"independent-energy flux Jacobian is ill-conditioned "
                f"({condition:.3e})"
            )
        scaled_derivatives = np.linalg.solve(
            scaled_jacobian, sources / row_scales
        )
        physical = column_scales * scaled_derivatives
        return np.concatenate((
            physical,
            [math.cos(state[E_THETA]), math.sin(state[E_THETA])],
        ))

    @staticmethod
    def _to_internal(state: np.ndarray) -> np.ndarray:
        rho_c, fraction_c, sysz, theta, uc, x, z = map(float, state)
        fraction_c = min(max(fraction_c, 1.0e-14), 1.0 - 1.0e-14)
        return np.array([
            math.log(rho_c),
            math.log(fraction_c / (1.0 - fraction_c)),
            math.log(sysz),
            theta,
            math.log(uc),
            x,
            z,
        ])

    @staticmethod
    def _from_internal(state: np.ndarray) -> np.ndarray:
        logistic = float(np.clip(state[E_Y], -35.0, 35.0))
        return np.array([
            math.exp(float(state[E_RHO])),
            1.0 / (1.0 + math.exp(-logistic)),
            math.exp(float(state[E_SYSZ])),
            float(state[E_THETA]),
            math.exp(float(state[E_UC])),
            float(state[E_X]),
            float(state[E_Z]),
        ])

    def solve(
        self,
        initial_state: np.ndarray,
        *,
        maximum_distance: float,
        maximum_step: float = 0.05,
        relative_tolerance: float = 2.0e-6,
        method: str = "flux_RK4",
        include_work: bool = False,
        include_thermal_relaxation: bool = False,
        include_phase_transition: bool = False,
        include_turbulent_heat_exchange: bool = False,
    ) -> IndependentEnergyJetResult:
        """Integrate the seven-state plume over additional arc length.

        By default, `include_work=False` for backward compatibility with the
        pre-registered baseline behavior.
        """
        initial_state = np.asarray(initial_state, dtype=float)
        self._physical(initial_state)
        if maximum_distance <= 0.0 or maximum_step <= 0.0:
            raise ValueError("integration distance and maximum step must be positive")

        def right_hand_side(s: float, internal: np.ndarray) -> np.ndarray:
            state = self._from_internal(internal)
            physical = self.derivatives(
                s,
                state,
                include_work=include_work,
                include_thermal_relaxation=include_thermal_relaxation,
                include_phase_transition=include_phase_transition,
                include_turbulent_heat_exchange=include_turbulent_heat_exchange,
            )
            fraction = state[E_Y]
            return np.array([
                physical[E_RHO] / state[E_RHO],
                physical[E_Y] / max(fraction * (1.0 - fraction), 1.0e-14),
                physical[E_SYSZ] / state[E_SYSZ],
                physical[E_THETA],
                physical[E_UC] / state[E_UC],
                physical[E_X],
                physical[E_Z],
            ])

        cumulative = None
        if method == "flux_RK4":
            steps = int(math.ceil(maximum_distance / maximum_step))
            arc_length = np.linspace(0.0, maximum_distance, steps + 1)
            states = np.empty((steps + 1, 7))
            states[0] = initial_state
            cumulative = np.zeros((steps + 1, 5))

            def stage_state(
                target_fluxes: np.ndarray,
                guess: np.ndarray,
                position: np.ndarray,
            ) -> np.ndarray:
                guess = guess.copy()
                guess[E_THETA] = math.atan2(
                    target_fluxes[3], target_fluxes[2]
                )
                guess[E_X:E_Z + 1] = position
                return self._match_flux_array(target_fluxes, guess)

            for index in range(steps):
                step = float(arc_length[index + 1] - arc_length[index])
                current = states[index]
                current_fluxes = self._as_array(self.integral_fluxes(current))
                source1 = self.source_terms(
                    current,
                    include_work=include_work,
                    include_thermal_relaxation=include_thermal_relaxation,
                    include_phase_transition=include_phase_transition,
                    include_turbulent_heat_exchange=include_turbulent_heat_exchange,
                )
                path1 = np.array([
                    math.cos(current[E_THETA]),
                    math.sin(current[E_THETA]),
                ])
                state2 = stage_state(
                    current_fluxes + 0.5 * step * source1,
                    current,
                    current[E_X:E_Z + 1] + 0.5 * step * path1,
                )
                source2 = self.source_terms(
                    state2,
                    include_work=include_work,
                    include_thermal_relaxation=include_thermal_relaxation,
                    include_phase_transition=include_phase_transition,
                    include_turbulent_heat_exchange=include_turbulent_heat_exchange,
                )
                path2 = np.array([
                    math.cos(state2[E_THETA]),
                    math.sin(state2[E_THETA]),
                ])
                state3 = stage_state(
                    current_fluxes + 0.5 * step * source2,
                    state2,
                    current[E_X:E_Z + 1] + 0.5 * step * path2,
                )
                source3 = self.source_terms(
                    state3,
                    include_work=include_work,
                    include_thermal_relaxation=include_thermal_relaxation,
                    include_phase_transition=include_phase_transition,
                    include_turbulent_heat_exchange=include_turbulent_heat_exchange,
                )
                path3 = np.array([
                    math.cos(state3[E_THETA]),
                    math.sin(state3[E_THETA]),
                ])
                state4 = stage_state(
                    current_fluxes + step * source3,
                    state3,
                    current[E_X:E_Z + 1] + step * path3,
                )
                source4 = self.source_terms(
                    state4,
                    include_work=include_work,
                    include_thermal_relaxation=include_thermal_relaxation,
                    include_phase_transition=include_phase_transition,
                    include_turbulent_heat_exchange=include_turbulent_heat_exchange,
                )
                path4 = np.array([
                    math.cos(state4[E_THETA]),
                    math.sin(state4[E_THETA]),
                ])
                flux_increment = step / 6.0 * (
                    source1 + 2.0 * source2 + 2.0 * source3 + source4
                )
                path_increment = step / 6.0 * (
                    path1 + 2.0 * path2 + 2.0 * path3 + path4
                )
                states[index + 1] = stage_state(
                    current_fluxes + flux_increment,
                    state4,
                    current[E_X:E_Z + 1] + path_increment,
                )
                cumulative[index + 1] = cumulative[index] + flux_increment
            integration_t = arc_length
        elif method == "state_RK4":
            steps = int(math.ceil(maximum_distance / maximum_step))
            arc_length = np.linspace(0.0, maximum_distance, steps + 1)
            internal = np.empty((steps + 1, 7))
            internal[0] = self._to_internal(initial_state)
            for index in range(steps):
                step = float(arc_length[index + 1] - arc_length[index])
                current = internal[index]
                k1 = right_hand_side(arc_length[index], current)
                k2 = right_hand_side(
                    arc_length[index] + 0.5 * step,
                    current + 0.5 * step * k1,
                )
                k3 = right_hand_side(
                    arc_length[index] + 0.5 * step,
                    current + 0.5 * step * k2,
                )
                k4 = right_hand_side(
                    arc_length[index] + step, current + step * k3
                )
                internal[index + 1] = current + step / 6.0 * (
                    k1 + 2.0 * k2 + 2.0 * k3 + k4
                )
            integration_t = arc_length
            integration_y = internal
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
                    f"independent-energy crosswind integration failed: "
                    f"{integration.message}"
                )
            integration_t = integration.t
            integration_y = integration.y.T
        if method != "flux_RK4":
            states = np.array([
                self._from_internal(column) for column in integration_y
            ])
        fluxes = np.array([
            self._as_array(self.integral_fluxes(state)) for state in states
        ])
        sources = np.array([
            self.source_terms(
                state,
                include_work=include_work,
                include_thermal_relaxation=include_thermal_relaxation,
                include_phase_transition=include_phase_transition,
                include_turbulent_heat_exchange=include_turbulent_heat_exchange,
            ) for state in states
        ])
        if cumulative is None:
            cumulative = cumulative_trapezoid(
                sources, integration_t, axis=0, initial=0.0
            )
        balances = fluxes - fluxes[0] - cumulative
        flux_scales = np.maximum(
            np.max(np.abs(fluxes), axis=0),
            np.array([1.0e-3, 1.0e-4, 1.0, 1.0, 1.0]),
        )
        maximum_balance = float(np.max(np.abs(balances) / flux_scales))
        temperatures = np.array([
            self.centre_temperature(state) for state in states
        ])
        mole_fractions = np.array([
            self.centre_mole_fraction(state) for state in states
        ])
        return IndependentEnergyJetResult(
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
