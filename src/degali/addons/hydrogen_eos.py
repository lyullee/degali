"""Pure-H2 HEOS departure screening, not an H2/air mixture EOS.

CoolProp implements Leachman et al. (2009), doi:10.1063/1.3160306.
Volume and caloric departures are reported together. Applying only Z to a
multicomponent, condensing jet is not a thermodynamically closed correction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class HydrogenGasDeparture:
    hydrogen_species: str
    temperature: float
    pressure: float
    density: float
    ideal_density: float
    compressibility: float
    enthalpy: float
    ideal_enthalpy: float
    enthalpy_departure: float

    @property
    def relative_volume_departure(self) -> float:
        """(v_real - v_ideal) / v_ideal at fixed T and P."""
        return self.compressibility - 1.0


def hydrogen_gas_departure(
    temperature: float, pressure: float, hydrogen_species: str = "Hydrogen",
) -> HydrogenGasDeparture:
    """Return the stable pure-gas volume and enthalpy departures in SI.

    ``pressure`` is the pressure of the pure-component comparison, not an
    asserted Dalton partial pressure in a nonideal mixture. Sub-saturation
    temperature states are rejected instead of imposing a metastable gas.
    Saturated vapour itself is allowed, within a 1e-8 relative pressure
    tolerance for the endpoint root solver.
    """
    import CoolProp as CP
    from CoolProp.CoolProp import AbstractState, PropsSI

    if hydrogen_species not in {"Hydrogen", "ParaHydrogen", "OrthoHydrogen"}:
        raise ValueError("invalid hydrogen spin species")
    if not all(math.isfinite(x) and x > 0 for x in (temperature, pressure)):
        raise ValueError("temperature and pressure must be finite and positive")
    triple = float(PropsSI("Ttriple", hydrogen_species))
    critical = float(PropsSI("Tcrit", hydrogen_species))
    if temperature < triple:
        raise ValueError("hydrogen temperature lies below the triple point")
    if temperature < critical:
        saturation = float(PropsSI("P", "T", temperature, "Q", 1, hydrogen_species))
        if pressure > saturation * (1.0 + 1.0e-8):
            raise ValueError("hydrogen is not a stable single gas at this T and P")
    state = AbstractState("HEOS", hydrogen_species)
    state.specify_phase(CP.iphase_gas)
    state.update(CP.PT_INPUTS, pressure, temperature)
    density = state.rhomass()
    ideal_density = pressure * state.molar_mass() / (state.gas_constant() * temperature)
    enthalpy, ideal_enthalpy = state.hmass(), state.hmass_idealgas()
    return HydrogenGasDeparture(
        hydrogen_species, temperature, pressure, density, ideal_density,
        ideal_density / density, enthalpy, ideal_enthalpy,
        enthalpy - ideal_enthalpy,
    )
