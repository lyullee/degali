"""Bounded tabular saturation properties for repeated LH2 flash calculations.

The table is deliberately restricted to the saturated-hydrogen states used by
the field equilibrium-flash adapter.  An out-of-range state is an error, not a
silent fallback to a different thermodynamic path.  This makes a speed-up
auditable and keeps direct CoolProp evaluation available as the reference.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class LH2SaturationTable:
    """Linear saturated-H2 property table with explicit valid temperature span."""

    temperature_k: np.ndarray
    pressure_pa: np.ndarray
    liquid_density_kg_m3: np.ndarray
    vapour_density_kg_m3: np.ndarray
    liquid_enthalpy_j_kg: np.ndarray
    vapour_enthalpy_j_kg: np.ndarray
    liquid_viscosity_pa_s: np.ndarray
    surface_tension_n_m: np.ndarray
    hydrogen_species: str = "Hydrogen"

    def __post_init__(self) -> None:
        arrays = (
            self.temperature_k,
            self.pressure_pa,
            self.liquid_density_kg_m3,
            self.vapour_density_kg_m3,
            self.liquid_enthalpy_j_kg,
            self.vapour_enthalpy_j_kg,
            self.liquid_viscosity_pa_s,
            self.surface_tension_n_m,
        )
        size = len(self.temperature_k)
        if size < 4 or any(np.asarray(value).shape != (size,) for value in arrays):
            raise ValueError("LH2 saturation table requires equally sized one-dimensional arrays")
        if any(not np.all(np.isfinite(value)) for value in arrays):
            raise ValueError("LH2 saturation table values must be finite")
        if np.any(np.diff(self.temperature_k) <= 0.0) or np.any(np.diff(self.pressure_pa) <= 0.0):
            raise ValueError("LH2 saturation table temperature and pressure must increase")
        if self.hydrogen_species not in {"Hydrogen", "ParaHydrogen", "OrthoHydrogen"}:
            raise ValueError("unsupported hydrogen species")

    @classmethod
    def build(
        cls,
        minimum_temperature_k: float,
        maximum_temperature_k: float,
        *,
        nodes: int = 161,
        hydrogen_species: str = "Hydrogen",
    ) -> "LH2SaturationTable":
        """Build a bounded reference table directly from CoolProp once."""
        from CoolProp.CoolProp import PropsSI

        if not all(math.isfinite(float(value)) for value in (
            minimum_temperature_k, maximum_temperature_k,
        )):
            raise ValueError("table temperatures must be finite")
        if nodes < 4:
            raise ValueError("LH2 saturation table requires at least four nodes")
        triple = float(PropsSI("Ttriple", hydrogen_species))
        critical = float(PropsSI("Tcrit", hydrogen_species))
        if not triple < minimum_temperature_k < maximum_temperature_k < critical:
            raise ValueError(
                "table bounds must lie strictly between hydrogen triple and critical temperatures"
            )
        temperature = np.linspace(minimum_temperature_k, maximum_temperature_k, nodes)

        def values(output: str, quality: int) -> np.ndarray:
            return np.array([
                float(PropsSI(output, "T", temp, "Q", quality, hydrogen_species))
                for temp in temperature
            ])

        return cls(
            temperature_k=temperature,
            pressure_pa=values("P", 0),
            liquid_density_kg_m3=values("D", 0),
            vapour_density_kg_m3=values("D", 1),
            liquid_enthalpy_j_kg=values("H", 0),
            vapour_enthalpy_j_kg=values("H", 1),
            liquid_viscosity_pa_s=values("V", 0),
            surface_tension_n_m=values("I", 0),
            hydrogen_species=hydrogen_species,
        )

    @property
    def minimum_temperature_k(self) -> float:
        return float(self.temperature_k[0])

    @property
    def maximum_temperature_k(self) -> float:
        return float(self.temperature_k[-1])

    @property
    def minimum_pressure_pa(self) -> float:
        return float(self.pressure_pa[0])

    @property
    def maximum_pressure_pa(self) -> float:
        return float(self.pressure_pa[-1])

    def _by_temperature(self, values: np.ndarray, temperature_k: float) -> float:
        if not self.minimum_temperature_k <= temperature_k <= self.maximum_temperature_k:
            raise ValueError(
                "temperature lies outside the LH2 saturation-table domain "
                f"[{self.minimum_temperature_k:g}, {self.maximum_temperature_k:g}] K"
            )
        return float(np.interp(temperature_k, self.temperature_k, values))

    def _by_pressure(self, values: np.ndarray, pressure_pa: float) -> float:
        if not self.minimum_pressure_pa <= pressure_pa <= self.maximum_pressure_pa:
            raise ValueError(
                "pressure lies outside the LH2 saturation-table domain "
                f"[{self.minimum_pressure_pa:g}, {self.maximum_pressure_pa:g}] Pa"
            )
        return float(np.interp(pressure_pa, self.pressure_pa, values))

    @staticmethod
    def _mixture(liquid: float, vapour: float, quality: float, *, density: bool) -> float:
        if not 0.0 <= quality <= 1.0:
            raise ValueError("vapour quality must lie in [0, 1]")
        if density:
            return 1.0 / ((1.0 - quality) / liquid + quality / vapour)
        return liquid + quality * (vapour - liquid)

    def props(
        self,
        output: str,
        key1: str,
        value1: float,
        key2: str,
        value2: float,
        fluid: str,
    ) -> float | None:
        """Return a supported PropsSI-equivalent value, else ``None``.

        Returning ``None`` only means that this narrow saturation table does
        not own the requested property (for example air density).  A supported
        H2 saturation query outside the table raises ``ValueError``.
        """
        if fluid != self.hydrogen_species:
            return None
        inputs = {key1: float(value1), key2: float(value2)}
        if set(inputs) == {"T", "Q"}:
            temperature, quality = inputs["T"], inputs["Q"]
            pressure = self._by_temperature(self.pressure_pa, temperature)
            liquid_density = self._by_temperature(self.liquid_density_kg_m3, temperature)
            vapour_density = self._by_temperature(self.vapour_density_kg_m3, temperature)
            liquid_enthalpy = self._by_temperature(self.liquid_enthalpy_j_kg, temperature)
            vapour_enthalpy = self._by_temperature(self.vapour_enthalpy_j_kg, temperature)
            if output == "P" and quality == 0.0:
                return pressure
            if output == "D":
                return self._mixture(liquid_density, vapour_density, quality, density=True)
            if output == "H":
                return self._mixture(liquid_enthalpy, vapour_enthalpy, quality, density=False)
            if output == "V" and quality == 0.0:
                return self._by_temperature(self.liquid_viscosity_pa_s, temperature)
            if output == "I" and quality == 0.0:
                return self._by_temperature(self.surface_tension_n_m, temperature)
            return None
        if set(inputs) == {"P", "Q"}:
            pressure, quality = inputs["P"], inputs["Q"]
            temperature = self._by_pressure(self.temperature_k, pressure)
            liquid_density = self._by_pressure(self.liquid_density_kg_m3, pressure)
            vapour_density = self._by_pressure(self.vapour_density_kg_m3, pressure)
            liquid_enthalpy = self._by_pressure(self.liquid_enthalpy_j_kg, pressure)
            vapour_enthalpy = self._by_pressure(self.vapour_enthalpy_j_kg, pressure)
            if output == "T":
                if not 0.0 <= quality <= 1.0:
                    raise ValueError("vapour quality must lie in [0, 1]")
                return temperature
            if output == "D":
                return self._mixture(liquid_density, vapour_density, quality, density=True)
            if output == "H":
                return self._mixture(liquid_enthalpy, vapour_enthalpy, quality, density=False)
            if output == "V" and quality == 0.0:
                return self._by_pressure(self.liquid_viscosity_pa_s, pressure)
            if output == "I" and quality == 0.0:
                return self._by_pressure(self.surface_tension_n_m, pressure)
            return None
        return None


__all__ = ["LH2SaturationTable"]
