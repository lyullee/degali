"""Constant-pressure thermodynamics for a conservative binary gas mixture.

Unlike :mod:`hydrogen_air_thermodynamics`, this closure does not hold a fixed
carrier-air inventory in each Eulerian cell.  The primary species variable is
its conserved mass concentration and the total mixture density is recovered
from a binary ideal-gas equation of state.  This is the thermodynamic contract
needed by dense initial-cloud experiments such as Thorney Island Trial 021.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ..core.constants import R_UNIVERSAL


@dataclass(frozen=True)
class BinaryIdealGasMixtureState:
    """Primitive and conserved fields for a binary ideal-gas mixture."""

    temperature_k: np.ndarray
    species_mole_fraction: np.ndarray
    species_mass_fraction: np.ndarray
    mixture_molar_mass_kg_mol: np.ndarray
    mixture_heat_capacity_j_kg_k: np.ndarray
    mixture_density_kg_m3: np.ndarray
    species_mass_concentration_kg_m3: np.ndarray
    carrier_mass_concentration_kg_m3: np.ndarray
    sensible_enthalpy_density_j_m3: np.ndarray
    buoyancy_acceleration_m_s2: np.ndarray
    relative_density_departure: np.ndarray


@dataclass(frozen=True)
class BinaryIdealGasMixtureClosure:
    """Calorically perfect binary-mixture closure at fixed pressure."""

    species_name: str
    species_molar_mass_kg_mol: float
    species_heat_capacity_j_kg_k: float
    carrier_name: str = "dry_air"
    carrier_molar_mass_kg_mol: float = 28.96546e-3
    carrier_heat_capacity_j_kg_k: float = 1005.0
    ambient_temperature_k: float = 291.15
    ambient_pressure_pa: float = 100_500.0
    gravity_m_s2: float = 9.80665
    minimum_temperature_k: float = 150.0
    maximum_temperature_k: float = 1_000.0

    def __post_init__(self) -> None:
        for name in ("species_name", "carrier_name"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
        for name in (
            "species_molar_mass_kg_mol", "species_heat_capacity_j_kg_k",
            "carrier_molar_mass_kg_mol", "carrier_heat_capacity_j_kg_k",
            "ambient_temperature_k", "ambient_pressure_pa",
            "minimum_temperature_k", "maximum_temperature_k",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
            if float(value) <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.minimum_temperature_k >= self.ambient_temperature_k:
            raise ValueError("minimum_temperature_k must be below ambient temperature")
        if self.maximum_temperature_k <= self.ambient_temperature_k:
            raise ValueError("maximum_temperature_k must exceed ambient temperature")
        if not math.isfinite(float(self.gravity_m_s2)) or self.gravity_m_s2 < 0.0:
            raise ValueError("gravity_m_s2 must be finite and non-negative")

    @property
    def ambient_density_kg_m3(self) -> float:
        return (
            self.ambient_pressure_pa * self.carrier_molar_mass_kg_mol
            / (R_UNIVERSAL * self.ambient_temperature_k)
        )

    def mass_fraction_from_mole_fraction(
        self, species_mole_fraction: np.ndarray | float,
    ) -> np.ndarray:
        x = np.asarray(species_mole_fraction, dtype=float)
        if np.any(~np.isfinite(x)) or np.any((x < 0.0) | (x > 1.0)):
            raise ValueError("species_mole_fraction must lie in [0, 1]")
        numerator = x * self.species_molar_mass_kg_mol
        denominator = numerator + (
            (1.0 - x) * self.carrier_molar_mass_kg_mol
        )
        return numerator / denominator

    def mole_fraction_from_mass_fraction(
        self, species_mass_fraction: np.ndarray | float,
    ) -> np.ndarray:
        y = np.asarray(species_mass_fraction, dtype=float)
        if np.any(~np.isfinite(y)) or np.any((y < 0.0) | (y > 1.0)):
            raise ValueError("species_mass_fraction must lie in [0, 1]")
        species_moles = y / self.species_molar_mass_kg_mol
        carrier_moles = (1.0 - y) / self.carrier_molar_mass_kg_mol
        return species_moles / (species_moles + carrier_moles)

    def _primitive_from_mass_fraction_and_enthalpy(
        self,
        species_mass_fraction: np.ndarray,
        sensible_enthalpy_density_j_m3: np.ndarray,
        *,
        validate_temperature: bool,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        y, enthalpy = np.broadcast_arrays(
            np.asarray(species_mass_fraction, dtype=float),
            np.asarray(sensible_enthalpy_density_j_m3, dtype=float),
        )
        reciprocal_molar_mass = (
            y / self.species_molar_mass_kg_mol
            + (1.0 - y) / self.carrier_molar_mass_kg_mol
        )
        molar_mass = 1.0 / reciprocal_molar_mass
        heat_capacity = (
            y * self.species_heat_capacity_j_kg_k
            + (1.0 - y) * self.carrier_heat_capacity_j_kg_k
        )
        nondimensional_enthalpy = (
            enthalpy * R_UNIVERSAL
            / (self.ambient_pressure_pa * molar_mass * heat_capacity)
        )
        denominator = 1.0 - nondimensional_enthalpy
        temperature = self.ambient_temperature_k / denominator
        density = (
            self.ambient_pressure_pa * molar_mass
            / (R_UNIVERSAL * temperature)
        )
        if validate_temperature and (
            np.any(~np.isfinite(temperature))
            or np.any(denominator <= 0.0)
            or np.any(
                (temperature < self.minimum_temperature_k)
                | (temperature > self.maximum_temperature_k)
            )
        ):
            raise ValueError(
                "reconstructed binary-mixture temperature lies outside bounds"
            )
        return temperature, density, molar_mass, heat_capacity

    def state_from_mass_fraction_and_sensible_enthalpy(
        self,
        species_mass_fraction: np.ndarray | float,
        sensible_enthalpy_density_j_m3: np.ndarray | float,
    ) -> BinaryIdealGasMixtureState:
        y, enthalpy = np.broadcast_arrays(
            np.asarray(species_mass_fraction, dtype=float),
            np.asarray(sensible_enthalpy_density_j_m3, dtype=float),
        )
        if np.any(~np.isfinite(y)) or np.any((y < 0.0) | (y > 1.0)):
            raise ValueError("species_mass_fraction must lie in [0, 1]")
        if np.any(~np.isfinite(enthalpy)):
            raise ValueError("sensible enthalpy density must be finite")
        temperature, density, molar_mass, heat_capacity = (
            self._primitive_from_mass_fraction_and_enthalpy(
                y, enthalpy, validate_temperature=True,
            )
        )
        x = self.mole_fraction_from_mass_fraction(y)
        relative_departure = density / self.ambient_density_kg_m3 - 1.0
        return BinaryIdealGasMixtureState(
            temperature_k=temperature,
            species_mole_fraction=x,
            species_mass_fraction=y,
            mixture_molar_mass_kg_mol=molar_mass,
            mixture_heat_capacity_j_kg_k=heat_capacity,
            mixture_density_kg_m3=density,
            species_mass_concentration_kg_m3=density * y,
            carrier_mass_concentration_kg_m3=density * (1.0 - y),
            sensible_enthalpy_density_j_m3=enthalpy,
            buoyancy_acceleration_m_s2=-self.gravity_m_s2 * relative_departure,
            relative_density_departure=relative_departure,
        )

    def state_from_mole_fraction_and_temperature(
        self,
        species_mole_fraction: np.ndarray | float,
        temperature_k: np.ndarray | float,
    ) -> BinaryIdealGasMixtureState:
        x, temperature = np.broadcast_arrays(
            np.asarray(species_mole_fraction, dtype=float),
            np.asarray(temperature_k, dtype=float),
        )
        if np.any(~np.isfinite(x)) or np.any((x < 0.0) | (x > 1.0)):
            raise ValueError("species_mole_fraction must lie in [0, 1]")
        if np.any(~np.isfinite(temperature)) or np.any(
            (temperature < self.minimum_temperature_k)
            | (temperature > self.maximum_temperature_k)
        ):
            raise ValueError("temperature_k lies outside closure bounds")
        y = self.mass_fraction_from_mole_fraction(x)
        molar_mass = 1.0 / (
            y / self.species_molar_mass_kg_mol
            + (1.0 - y) / self.carrier_molar_mass_kg_mol
        )
        heat_capacity = (
            y * self.species_heat_capacity_j_kg_k
            + (1.0 - y) * self.carrier_heat_capacity_j_kg_k
        )
        density = (
            self.ambient_pressure_pa * molar_mass
            / (R_UNIVERSAL * temperature)
        )
        enthalpy = density * heat_capacity * (
            temperature - self.ambient_temperature_k
        )
        return self.state_from_mass_fraction_and_sensible_enthalpy(y, enthalpy)

    def state_from_mass_fraction_and_specific_sensible_enthalpy(
        self,
        species_mass_fraction: np.ndarray | float,
        specific_sensible_enthalpy_j_kg: np.ndarray | float,
    ) -> BinaryIdealGasMixtureState:
        """Recover a fixed-pressure inlet state from conserved flux ratios.

        At a unidirectional boundary, species-flux/total-mass-flux gives the
        mass fraction and enthalpy-flux/total-mass-flux gives specific sensible
        enthalpy.  For the calorically perfect closure this determines
        temperature directly before the ideal-gas density is evaluated.
        """

        y, specific_enthalpy = np.broadcast_arrays(
            np.asarray(species_mass_fraction, dtype=float),
            np.asarray(specific_sensible_enthalpy_j_kg, dtype=float),
        )
        if np.any(~np.isfinite(y)) or np.any((y < 0.0) | (y > 1.0)):
            raise ValueError("species_mass_fraction must lie in [0, 1]")
        if np.any(~np.isfinite(specific_enthalpy)):
            raise ValueError("specific sensible enthalpy must be finite")
        heat_capacity = (
            y * self.species_heat_capacity_j_kg_k
            + (1.0 - y) * self.carrier_heat_capacity_j_kg_k
        )
        temperature = self.ambient_temperature_k + (
            specific_enthalpy / heat_capacity
        )
        if np.any(
            (temperature < self.minimum_temperature_k)
            | (temperature > self.maximum_temperature_k)
        ):
            raise ValueError("inlet temperature lies outside closure bounds")
        molar_mass = 1.0 / (
            y / self.species_molar_mass_kg_mol
            + (1.0 - y) / self.carrier_molar_mass_kg_mol
        )
        density = (
            self.ambient_pressure_pa * molar_mass
            / (R_UNIVERSAL * temperature)
        )
        return self.state_from_mass_fraction_and_sensible_enthalpy(
            y, density * specific_enthalpy,
        )

    def state_from_species_mass_concentration_and_sensible_enthalpy(
        self,
        species_mass_concentration_kg_m3: np.ndarray | float,
        sensible_enthalpy_density_j_m3: np.ndarray | float,
    ) -> BinaryIdealGasMixtureState:
        """Invert ``rho*Y`` and sensible enthalpy without a fixed air inventory."""
        concentration, enthalpy = np.broadcast_arrays(
            np.asarray(species_mass_concentration_kg_m3, dtype=float),
            np.asarray(sensible_enthalpy_density_j_m3, dtype=float),
        )
        if np.any(~np.isfinite(concentration)) or np.any(concentration < 0.0):
            raise ValueError(
                "species mass concentration must be finite and non-negative"
            )
        if np.any(~np.isfinite(enthalpy)):
            raise ValueError("sensible enthalpy density must be finite")

        # Let c be the known species mass concentration and a the unknown
        # carrier concentration.  Combining p/(R*T) = c/Ms + a/Ma with
        # H = (cp_s*c + cp_a*a) * (T - T0) gives a quadratic in a.  Solving
        # it directly avoids 64 full-array EOS evaluations per corrector call.
        cp_s = self.species_heat_capacity_j_kg_k
        cp_a = self.carrier_heat_capacity_j_kg_k
        molar_s = self.species_molar_mass_kg_mol
        molar_a = self.carrier_molar_mass_kg_mol
        pressure_over_r = self.ambient_pressure_pa / R_UNIVERSAL
        _, pure_density, _, _ = self._primitive_from_mass_fraction_and_enthalpy(
            np.ones_like(concentration), enthalpy, validate_temperature=False,
        )
        if np.any(~np.isfinite(pure_density)) or np.any(pure_density <= 0.0):
            raise ValueError("enthalpy does not admit a physical pure-species state")
        if np.any(concentration > pure_density * (1.0 + 1.0e-12)):
            raise ValueError(
                "species concentration exceeds the fixed-pressure pure-species bound"
            )
        d0 = cp_s * concentration
        e0 = self.ambient_temperature_k * d0 + enthalpy
        f0 = concentration / molar_s
        f1 = 1.0 / molar_a
        coefficient_a = self.ambient_temperature_k * cp_a * f1
        coefficient_b = e0 * f1 + self.ambient_temperature_k * cp_a * f0 - pressure_over_r * cp_a
        coefficient_c = e0 * f0 - pressure_over_r * d0
        discriminant = coefficient_b * coefficient_b - 4.0 * coefficient_a * coefficient_c
        disc_scale = np.maximum(
            coefficient_b * coefficient_b
            + np.abs(4.0 * coefficient_a * coefficient_c),
            1.0,
        )
        disc_tolerance = 32.0 * np.finfo(float).eps * disc_scale
        if np.any(discriminant < -disc_tolerance):
            raise ValueError("conserved state has no physical binary-mixture root")
        discriminant = np.maximum(discriminant, 0.0)

        # The q-form of the quadratic formula avoids cancellation when one
        # root is much smaller than the other.  Select the non-negative root;
        # the zero-species ambient-air state correctly selects the larger root.
        root_disc = np.sqrt(discriminant)
        q = -0.5 * (coefficient_b + np.where(coefficient_b >= 0.0, root_disc, -root_disc))
        root_one = q / coefficient_a
        root_two = np.divide(
            coefficient_c, q, out=np.zeros_like(coefficient_c), where=q != 0.0,
        )
        carrier_concentration = np.where(root_one >= 0.0, root_one, root_two)
        pure_limit = concentration >= pure_density * (1.0 - 1.0e-12)
        carrier_concentration = np.where(
            pure_limit & (concentration > 0.0), 0.0, carrier_concentration,
        )
        if np.any(~np.isfinite(carrier_concentration)) or np.any(
            carrier_concentration < 0.0
        ):
            raise ValueError("conserved state has no non-negative carrier concentration")
        total_density = concentration + carrier_concentration
        y = np.divide(
            concentration, total_density, out=np.zeros_like(concentration),
            where=total_density > 0.0,
        )
        temperature = self.ambient_temperature_k + np.divide(
            enthalpy,
            cp_s * concentration + cp_a * carrier_concentration,
            out=np.zeros_like(enthalpy),
            where=(cp_s * concentration + cp_a * carrier_concentration) > 0.0,
        )
        if np.any(~np.isfinite(temperature)) or np.any(
            (temperature < self.minimum_temperature_k)
            | (temperature > self.maximum_temperature_k)
        ):
            raise ValueError("reconstructed binary-mixture temperature lies outside bounds")
        state = self.state_from_mass_fraction_and_sensible_enthalpy(y, enthalpy)
        error = np.abs(
            state.species_mass_concentration_kg_m3 - concentration
        )
        tolerance = 2.0e-12 * np.maximum(1.0, concentration)
        if np.any(error > tolerance):
            raise RuntimeError("binary-mixture conserved-state inversion failed")
        return state


__all__ = [
    "BinaryIdealGasMixtureClosure",
    "BinaryIdealGasMixtureState",
]
