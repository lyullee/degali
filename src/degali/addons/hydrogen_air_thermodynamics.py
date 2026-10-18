"""Thermodynamic closure for an active H2/air thermal scalar.

The resolved obstacle prototype conserves hydrogen mass concentration and a
signed sensible-enthalpy density.  This module maps those two conserved
quantities to temperature, hydrogen mole fraction, ideal-gas mixture density
and vertical buoyancy acceleration.

The closure deliberately keeps the ambient-air carrier inventory fixed per
Eulerian control volume.  That is the standard active-scalar approximation,
not a replacement for a variable-density low-Mach continuity/EOS solve.  A
relative-density-departure field is therefore returned so callers can fail
closed when the approximation is used outside its declared range.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ..core.constants import R_UNIVERSAL
from .field_observation import (
    DRY_AIR_MOLAR_MASS_KG_MOL,
    HYDROGEN_MOLAR_MASS_KG_MOL,
)


@dataclass(frozen=True)
class HydrogenAirThermodynamicState:
    """Primitive fields reconstructed from the two conserved scalars."""

    temperature_k: np.ndarray
    hydrogen_mole_fraction: np.ndarray
    hydrogen_mass_fraction: np.ndarray
    mixture_molar_mass_kg_mol: np.ndarray
    mixture_density_kg_m3: np.ndarray
    buoyancy_acceleration_m_s2: np.ndarray
    relative_density_departure: np.ndarray
    active_scalar_applicability_mask: np.ndarray


@dataclass(frozen=True)
class HydrogenAirActiveScalarClosure:
    """Constant-cp ideal-gas closure with an explicit applicability bound."""

    ambient_temperature_k: float = 293.15
    ambient_pressure_pa: float = 101_325.0
    ambient_air_density_kg_m3: float = 1.204
    dry_air_heat_capacity_j_kg_k: float = 1005.0
    hydrogen_heat_capacity_j_kg_k: float = 14_300.0
    gravity_m_s2: float = 9.80665
    maximum_relative_density_departure: float = 0.10
    minimum_temperature_k: float = 14.0
    maximum_temperature_k: float = 1_000.0

    def __post_init__(self) -> None:
        for name in (
            "ambient_temperature_k", "ambient_pressure_pa",
            "ambient_air_density_kg_m3", "dry_air_heat_capacity_j_kg_k",
            "hydrogen_heat_capacity_j_kg_k",
            "maximum_relative_density_departure", "minimum_temperature_k",
            "maximum_temperature_k",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
            if float(value) <= 0.0:
                raise ValueError(f"{name} must be positive")
        gravity = self.gravity_m_s2
        if isinstance(gravity, bool) or not math.isfinite(float(gravity)):
            raise ValueError("gravity_m_s2 must be finite")
        if float(gravity) < 0.0:
            raise ValueError("gravity_m_s2 must be non-negative")
        if self.minimum_temperature_k >= self.ambient_temperature_k:
            raise ValueError("minimum_temperature_k must be below ambient temperature")
        if self.maximum_temperature_k <= self.ambient_temperature_k:
            raise ValueError("maximum_temperature_k must exceed ambient temperature")

    def source_specific_sensible_enthalpy_j_kg(
        self, source_temperature_k: float,
    ) -> float:
        """Return source H2 sensible enthalpy relative to ambient H2 gas."""
        temperature = float(source_temperature_k)
        if not math.isfinite(temperature) or not (
            self.minimum_temperature_k <= temperature <= self.maximum_temperature_k
        ):
            raise ValueError("source_temperature_k lies outside closure bounds")
        return self.hydrogen_heat_capacity_j_kg_k * (
            temperature - self.ambient_temperature_k
        )

    def neutral_buoyancy_temperature_k(
        self, hydrogen_mole_fraction: np.ndarray | float,
    ) -> np.ndarray:
        """Temperature at which an ideal H2/air mixture matches ambient density."""
        fraction = np.asarray(hydrogen_mole_fraction, dtype=float)
        if np.any(~np.isfinite(fraction)) or np.any((fraction < 0.0) | (fraction > 1.0)):
            raise ValueError("hydrogen_mole_fraction must lie in [0, 1]")
        mixture_molar_mass = (
            fraction * HYDROGEN_MOLAR_MASS_KG_MOL
            + (1.0 - fraction) * DRY_AIR_MOLAR_MASS_KG_MOL
        )
        return (
            self.ambient_temperature_k * mixture_molar_mass
            / DRY_AIR_MOLAR_MASS_KG_MOL
        )

    def state_from_conserved(
        self,
        hydrogen_mass_concentration_kg_m3: np.ndarray | float,
        sensible_enthalpy_density_j_m3: np.ndarray | float,
    ) -> HydrogenAirThermodynamicState:
        """Recover primitive fields without altering either conserved scalar."""
        concentration, enthalpy = np.broadcast_arrays(
            np.asarray(hydrogen_mass_concentration_kg_m3, dtype=float),
            np.asarray(sensible_enthalpy_density_j_m3, dtype=float),
        )
        if np.any(~np.isfinite(concentration)) or np.any(concentration < 0.0):
            raise ValueError("hydrogen mass concentration must be finite and non-negative")
        if np.any(~np.isfinite(enthalpy)):
            raise ValueError("sensible enthalpy density must be finite")

        volumetric_heat_capacity = (
            self.ambient_air_density_kg_m3 * self.dry_air_heat_capacity_j_kg_k
            + concentration * self.hydrogen_heat_capacity_j_kg_k
        )
        temperature = self.ambient_temperature_k + enthalpy / volumetric_heat_capacity
        if np.any(~np.isfinite(temperature)) or np.any(
            (temperature < self.minimum_temperature_k)
            | (temperature > self.maximum_temperature_k)
        ):
            minimum = float(np.min(temperature))
            maximum = float(np.max(temperature))
            raise ValueError(
                "reconstructed temperature lies outside closure bounds: "
                f"[{minimum:.6g}, {maximum:.6g}] K"
            )

        hydrogen_moles = concentration / HYDROGEN_MOLAR_MASS_KG_MOL
        carrier_air_moles = (
            self.ambient_air_density_kg_m3 / DRY_AIR_MOLAR_MASS_KG_MOL
        )
        mole_fraction = hydrogen_moles / (hydrogen_moles + carrier_air_moles)
        mass_fraction = concentration / (
            concentration + self.ambient_air_density_kg_m3
        )
        mixture_molar_mass = (
            mole_fraction * HYDROGEN_MOLAR_MASS_KG_MOL
            + (1.0 - mole_fraction) * DRY_AIR_MOLAR_MASS_KG_MOL
        )
        mixture_density = (
            self.ambient_pressure_pa * mixture_molar_mass
            / (R_UNIVERSAL * temperature)
        )
        buoyancy_reference_density = (
            self.ambient_pressure_pa * DRY_AIR_MOLAR_MASS_KG_MOL
            / (R_UNIVERSAL * self.ambient_temperature_k)
        )
        relative_departure = (
            (mixture_density - buoyancy_reference_density)
            / buoyancy_reference_density
        )
        buoyancy = -self.gravity_m_s2 * relative_departure
        applicability = (
            np.abs(relative_departure) <= self.maximum_relative_density_departure
        )
        return HydrogenAirThermodynamicState(
            temperature_k=temperature,
            hydrogen_mole_fraction=mole_fraction,
            hydrogen_mass_fraction=mass_fraction,
            mixture_molar_mass_kg_mol=mixture_molar_mass,
            mixture_density_kg_m3=mixture_density,
            buoyancy_acceleration_m_s2=buoyancy,
            relative_density_departure=relative_departure,
            active_scalar_applicability_mask=applicability,
        )


__all__ = [
    "HydrogenAirActiveScalarClosure",
    "HydrogenAirThermodynamicState",
]
