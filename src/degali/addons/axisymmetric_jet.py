"""Conserved-energy axisymmetric gas jet/plume.

This is an independent implementation of the published integral balances in
SAND2025-04942, equations 65--92.  It is deliberately separate from the
Fortran-compatible :mod:`degali.core.jetplume` path: the latter must remain
bit-for-bit available, while this model supplies the missing near-nozzle
physics for a quiescent cryogenic jet.

The five conservation laws are evaluated as flux integrals and differentiated
numerically.  This avoids transcribing another implementation's expanded
algebra and makes every term auditable against its physical balance:
continuity, horizontal and vertical momentum, species, and total energy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp, trapezoid
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq, least_squares

from ..core.constants import CPW, DHFUS, DHVAP, RHOWL, WMW_MODERN
from .hydrogen_eos import hydrogen_gas_departure


R_UNIVERSAL = 8.31446261815324  # J/(mol K); molecular weights are kg/mol
GRAVITY = 9.80665
_CONDENSED_LOOKUP_CACHE: dict[tuple[object, ...], tuple[np.ndarray, ...]] = {}


def phase_dry_air_composition(include_argon: bool = False) -> tuple[float, tuple[float, ...]]:
    """Return dry molar mass (kg/mol) and N2/O2/Ar mass fractions.

    These are the explicit species of the research phase-volume closure,
    not the pseudo-pure CoolProp ``Air`` composition.
    """
    mole_n, mole_o, mole_ar = 0.78084, 0.20946, 0.00934 if include_argon else 0.0
    mole_scale = mole_n + mole_o + mole_ar
    x_n, x_o, x_ar = mole_n / mole_scale, mole_o / mole_scale, mole_ar / mole_scale
    mw_n, mw_o, mw_ar = 0.0280134, 0.0319988, 0.039948
    mw_dry = x_n * mw_n + x_o * mw_o + x_ar * mw_ar
    return mw_dry, (x_n * mw_n / mw_dry, x_o * mw_o / mw_dry, x_ar * mw_ar / mw_dry)


def phase_ambient_from_rh(
    temperature: float, pressure: float, relative_humidity: float,
    *, include_argon: bool = False,
) -> tuple[float, float, float]:
    """Return dry molar mass, kg-water/kg-dry-air, and ideal phase density.

    Use the same saturation table, species and gas-volume equation as the
    phase inversion. This opt-in consistency helper is NOT a real-air EOS.
    """
    if not all(math.isfinite(v) for v in (temperature, pressure, relative_humidity)):
        raise ValueError("phase ambient inputs must be finite")
    if not 14.5 <= temperature <= 300.0 or pressure <= 0.0:
        raise ValueError("phase ambient requires 14.5--300 K and positive pressure")
    if not 0.0 <= relative_humidity <= 100.0:
        raise ValueError("relative humidity must lie between 0 and 100 percent")
    table = _air_phase_property_table()
    pv = relative_humidity / 100.0 * float(np.interp(
        temperature, table["temperature"], table["water_saturation"]
    ))
    if pv >= pressure:
        raise ValueError("water vapour pressure must be below ambient pressure")
    mw_dry, _ = phase_dry_air_composition(include_argon)
    mw_water = WMW_MODERN / 1000.0
    humidity = mw_water / mw_dry * pv / (pressure - pv)
    mw_humid = (1.0 + humidity) / (1.0 / mw_dry + humidity / mw_water)
    return mw_dry, humidity, pressure * mw_humid / (R_UNIVERSAL * temperature)


@lru_cache(maxsize=8)
def _ideal_component_enthalpy_table(
    species: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a low-density Helmholtz ideal-gas enthalpy table."""
    import CoolProp as CP
    from CoolProp.CoolProp import AbstractState

    temperature = np.linspace(14.1, 400.0, 1931)
    state = AbstractState("HEOS", species)
    values = []
    for value in temperature:
        try:
            state.update(CP.DmolarT_INPUTS, 1.0, float(value))
            values.append(state.hmass_idealgas())
        except ValueError:
            values.append(np.nan)
    values = np.asarray(values)
    valid = np.flatnonzero(np.isfinite(values))
    if len(valid) < 2:
        raise RuntimeError(f"insufficient ideal-gas enthalpy data for {species}")
    first = int(valid[0])
    if first:
        slope = (
            (values[valid[1]] - values[first])
            / (temperature[valid[1]] - temperature[first])
        )
        values[:first] = values[first] + slope * (
            temperature[:first] - temperature[first]
        )
    last = int(valid[-1])
    interior_missing = (
        ~np.isfinite(values)
        & (temperature >= temperature[first])
        & (temperature <= temperature[last])
    )
    values[interior_missing] = np.interp(
        temperature[interior_missing], temperature[valid], values[valid]
    )
    if last < len(values) - 1:
        slope = (
            (values[last] - values[valid[-2]])
            / (temperature[last] - temperature[valid[-2]])
        )
        values[last + 1:] = values[last] + slope * (
            temperature[last + 1:] - temperature[last]
        )
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"non-finite ideal-gas enthalpy data for {species}")
    return temperature, values


@lru_cache(maxsize=4)
def _ideal_gas_enthalpy_table(
    fuel: str, ambient: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return low-density CoolProp ideal-gas enthalpy tables."""
    import CoolProp as CP
    from CoolProp.CoolProp import AbstractState

    temperature = np.linspace(14.1, 400.0, 1931)
    enthalpies = []
    for fluid in (fuel, ambient):
        state = AbstractState("HEOS", fluid)
        values = []
        for value in temperature:
            # This exposes the Helmholtz ideal term without asking for a
            # metastable condensed-air pressure/temperature flash.
            try:
                state.update(CP.DmolarT_INPUTS, 1.0, float(value))
                values.append(state.hmass_idealgas())
            except ValueError:
                values.append(np.nan)
        values = np.asarray(values)
        valid = np.flatnonzero(np.isfinite(values))
        if len(valid) < 2:
            raise RuntimeError(f"insufficient ideal-gas enthalpy data for {fluid}")
        first = int(valid[0])
        if first:
            slope = (
                (values[valid[1]] - values[first])
                / (temperature[valid[1]] - temperature[first])
            )
            values[:first] = values[first] + slope * (
                temperature[:first] - temperature[first]
            )
        enthalpies.append(values)
    return temperature, enthalpies[0], enthalpies[1]


@lru_cache(maxsize=1)
def _air_phase_property_table() -> dict[str, np.ndarray]:
    """Stable N2/O2/Ar/H2O saturation and condensed-property table."""
    from CoolProp.CoolProp import PropsSI

    from .cryogenic_air import (
        _CRITICAL_TEMPERATURE,
        _SOLID_DENSITY,
        _SOLID_SUBLIMATION_ENTHALPY,
        _TRIPLE_TEMPERATURE,
        air_saturation_pressure,
    )

    temperature = np.linspace(14.1, 300.0, 721)
    out: dict[str, np.ndarray] = {"temperature": temperature}
    for species in ("Nitrogen", "Oxygen"):
        saturation, latent, condensed_density = [], [], []
        triple = _TRIPLE_TEMPERATURE[species]
        critical = _CRITICAL_TEMPERATURE[species]
        for value in temperature:
            if value >= critical:
                saturation.append(1.0e30)
                latent.append(0.0)
                condensed_density.append(_SOLID_DENSITY[species])
            elif value < triple:
                saturation.append(air_saturation_pressure(species, float(value)))
                latent.append(_SOLID_SUBLIMATION_ENTHALPY[species])
                condensed_density.append(_SOLID_DENSITY[species])
            else:
                saturation.append(float(PropsSI("P", "T", value, "Q", 1, species)))
                latent.append(float(
                    PropsSI("H", "T", value, "Q", 1, species)
                    - PropsSI("H", "T", value, "Q", 0, species)
                ))
                condensed_density.append(float(
                    PropsSI("D", "T", value, "Q", 0, species)
                ))
        key = species.lower()
        out[f"{key}_saturation"] = np.array(saturation)
        out[f"{key}_latent"] = np.array(latent)
        out[f"{key}_density"] = np.array(condensed_density)

    # Argon is 0.934 mol % of dry atmospheric air. CoolProp starts at the
    # fluid triple point, so below it use a constant-latent
    # Clausius-Clapeyron continuation anchored to that same point. The
    # 7.79 kJ/mol sublimation enthalpy is the NIST thin-film measurement used
    # in the pre-registered completeness screen; it is not fitted here.
    argon_mw = float(PropsSI("M", "Argon"))
    argon_triple = float(PropsSI("Ttriple", "Argon"))
    argon_critical = float(PropsSI("Tcrit", "Argon"))
    argon_triple_pressure = float(PropsSI("ptriple", "Argon"))
    argon_sublimation = 7.79e3 / argon_mw
    argon_solid_density = float(PropsSI(
        "D", "T", argon_triple, "Q", 0, "Argon"
    ))
    argon_saturation, argon_latent, argon_density = [], [], []
    for value in temperature:
        if value >= argon_critical:
            argon_saturation.append(1.0e30)
            argon_latent.append(0.0)
            argon_density.append(argon_solid_density)
        elif value < argon_triple:
            argon_saturation.append(
                argon_triple_pressure * math.exp(
                    -(argon_sublimation * argon_mw) / R_UNIVERSAL
                    * (1.0 / float(value) - 1.0 / argon_triple)
                )
            )
            argon_latent.append(argon_sublimation)
            argon_density.append(argon_solid_density)
        else:
            argon_saturation.append(float(
                PropsSI("P", "T", value, "Q", 1, "Argon")
            ))
            argon_latent.append(float(
                PropsSI("H", "T", value, "Q", 1, "Argon")
                - PropsSI("H", "T", value, "Q", 0, "Argon")
            ))
            argon_density.append(float(
                PropsSI("D", "T", value, "Q", 0, "Argon")
            ))
    out["argon_saturation"] = np.array(argon_saturation)
    out["argon_latent"] = np.array(argon_latent)
    out["argon_density"] = np.array(argon_density)

    water_triple = 273.16
    water_triple_pressure = float(
        PropsSI("P", "T", water_triple, "Q", 0, "Water")
    )
    legacy_triple_pressure = math.exp(
        14.683943 - 5407.0 / water_triple
    ) * 101325.0
    water_saturation, water_latent, water_density = [], [], []
    for value in temperature:
        if value < water_triple:
            water_saturation.append(
                math.exp(14.683943 - 5407.0 / float(value))
                * 101325.0 * water_triple_pressure
                / legacy_triple_pressure
            )
            water_latent.append(
                DHVAP + DHFUS * min((273.15 - float(value)) / 10.0, 1.0)
            )
            water_density.append(RHOWL)
        else:
            water_saturation.append(float(
                PropsSI("P", "T", value, "Q", 0, "Water")
            ))
            water_latent.append(float(
                PropsSI("H", "T", value, "Q", 1, "Water")
                - PropsSI("H", "T", value, "Q", 0, "Water")
            ))
            water_density.append(float(
                PropsSI("D", "T", value, "Q", 0, "Water")
            ))
    out["water_saturation"] = np.array(water_saturation)
    out["water_latent"] = np.array(water_latent)
    out["water_density"] = np.array(water_density)
    return out


@dataclass(frozen=True)
class SourceEnthalpyBoundary:
    """Explicit phase-source enthalpy, J/kg of the complete mixture.

    Air species use their ambient-temperature gas as zero. The absolute
    H2 reference is recorded so the receiving component ideal-gas tables
    can change reference without discarding residual or phase enthalpy.
    Kinetic energy is excluded and is added by the flux balance exactly once.
    This is an energy boundary, not a replacement mixture EOS.
    """

    specific_enthalpy: float
    ambient_temperature: float
    hydrogen_species: str
    hydrogen_reference_enthalpy: float


@dataclass(frozen=True)
class AxisymmetricJetSource:
    """One internally consistent ambient-pressure plug-flow source plane."""

    diameter: float
    velocity: float
    density: float
    temperature: float
    mass_fraction: float = 1.0
    theta: float = math.pi / 2.0
    x: float = 0.0
    y: float = 0.0
    enthalpy_boundary: SourceEnthalpyBoundary | None = None

    @property
    def area(self) -> float:
        return math.pi * self.diameter**2 / 4.0

    @property
    def mass_flow(self) -> float:
        return self.density * self.velocity * self.area

    @property
    def fuel_mass_flow(self) -> float:
        return self.mass_fraction * self.mass_flow


@dataclass(frozen=True)
class InitialEntrainmentResult:
    """Conservative plug state after HyRAM+ initial heating equations 58--64."""

    source: AxisymmetricJetSource
    initial_mass_flow: float
    entrained_air_mass_flow: float
    length: float
    fuel_mass_residual: float
    momentum_residual: float
    energy_residual: float
    geometry_mass_residual: float


def entrain_and_heat_initial_plug(
    source: AxisymmetricJetSource,
    *,
    minimum_temperature: float,
    ambient_temperature: float,
    ambient_pressure: float,
    ambient_density: float,
    fuel_molecular_weight: float,
    ambient_molecular_weight: float,
    fuel_heat_capacity: float,
    ambient_heat_capacity: float,
    momentum_entrainment_beta: float = 0.28,
) -> InitialEntrainmentResult:
    """Apply published plug-flow initial entrainment and heating balances."""
    if not math.isclose(source.mass_fraction, 1.0, abs_tol=1.0e-12):
        raise ValueError("initial entrainment currently requires a pure-fuel plug")
    if min(
        minimum_temperature, ambient_temperature, ambient_pressure,
        ambient_density, fuel_molecular_weight, ambient_molecular_weight,
        fuel_heat_capacity, ambient_heat_capacity, momentum_entrainment_beta,
    ) <= 0.0:
        raise ValueError("initial entrainment inputs must be positive")
    if minimum_temperature >= ambient_temperature:
        raise ValueError("minimum source temperature must be below ambient")
    if source.temperature >= minimum_temperature:
        return InitialEntrainmentResult(
            source=source,
            initial_mass_flow=source.mass_flow,
            entrained_air_mass_flow=0.0,
            length=0.0,
            fuel_mass_residual=0.0,
            momentum_residual=0.0,
            energy_residual=0.0,
            geometry_mass_residual=0.0,
        )

    h_air_min = ambient_heat_capacity * minimum_temperature
    h_fuel_min = fuel_heat_capacity * minimum_temperature
    h_air_ambient = ambient_heat_capacity * ambient_temperature
    h_in_total = (
        fuel_heat_capacity * source.temperature + 0.5 * source.velocity**2
    )

    def energy_balance(fuel_fraction: float) -> float:
        # Equation 61 multiplied by Y/m_dot_in avoids cancellation at Y -> 0.
        return (
            (1.0 - fuel_fraction) * (h_air_min - h_air_ambient)
            + fuel_fraction * (h_fuel_min - h_in_total)
            + 0.5 * source.velocity**2 * fuel_fraction**2
        )

    fuel_fraction = float(brentq(
        energy_balance, 1.0e-12, 1.0, xtol=1.0e-14, rtol=1.0e-14
    ))
    mass_out = source.mass_flow / fuel_fraction
    air_mass = mass_out - source.mass_flow
    velocity_out = source.velocity * fuel_fraction
    rho_air_min = (
        ambient_pressure * ambient_molecular_weight
        / (R_UNIVERSAL * minimum_temperature)
    )
    rho_fuel_min = (
        ambient_pressure * fuel_molecular_weight
        / (R_UNIVERSAL * minimum_temperature)
    )
    density_out = 1.0 / (
        (1.0 - fuel_fraction) / rho_air_min
        + fuel_fraction / rho_fuel_min
    )
    diameter_out = math.sqrt(
        4.0 * mass_out / (math.pi * density_out * velocity_out)
    )
    momentum_entrainment = momentum_entrainment_beta * math.sqrt(
        source.area * source.density * source.velocity**2 / ambient_density
    )
    length = (
        (1.0 - fuel_fraction) * air_mass
        / (ambient_density * momentum_entrainment)
    )
    heated = AxisymmetricJetSource(
        diameter=diameter_out,
        velocity=velocity_out,
        density=density_out,
        temperature=minimum_temperature,
        mass_fraction=fuel_fraction,
        theta=source.theta,
        x=source.x + length * math.cos(source.theta),
        y=source.y + length * math.sin(source.theta),
    )
    output_energy = mass_out * (
        fuel_fraction * h_fuel_min
        + (1.0 - fuel_fraction) * h_air_min
        + 0.5 * velocity_out**2
    )
    input_energy = (
        source.mass_flow * h_in_total + air_mass * h_air_ambient
    )
    scale_energy = max(abs(input_energy), abs(output_energy), 1.0)
    return InitialEntrainmentResult(
        source=heated,
        initial_mass_flow=source.mass_flow,
        entrained_air_mass_flow=air_mass,
        length=length,
        fuel_mass_residual=abs(
            heated.fuel_mass_flow - source.mass_flow
        ) / source.mass_flow,
        momentum_residual=abs(
            mass_out * velocity_out - source.mass_flow * source.velocity
        ) / max(abs(source.mass_flow * source.velocity), 1.0e-30),
        energy_residual=abs(output_energy - input_energy) / scale_energy,
        geometry_mass_residual=abs(
            heated.mass_flow - mass_out
        ) / mass_out,
    )


@dataclass(frozen=True)
class AxisymmetricJetResult:
    """Integrated centreline solution in increasing streamline distance."""

    S: np.ndarray
    velocity: np.ndarray
    width: np.ndarray
    density: np.ndarray
    mass_fraction: np.ndarray
    theta: np.ndarray
    x: np.ndarray
    y: np.ndarray
    temperature: np.ndarray
    mass_flux: np.ndarray
    species_flux: np.ndarray
    momentum_flux: np.ndarray
    energy_flux: np.ndarray
    radiative_heat_added: np.ndarray

    def state_at_s(self, distance: float) -> np.ndarray:
        """Interpolate the seven ODE variables at streamline distance ``S``."""
        if not self.S[0] <= distance <= self.S[-1]:
            raise ValueError(
                f"streamline distance {distance:g} m lies outside the "
                f"computed range {self.S[0]:g}--{self.S[-1]:g} m"
            )
        values = (
            self.velocity, self.width, self.density, self.mass_fraction,
            self.theta, self.x, self.y,
        )
        return np.array([
            np.interp(distance, self.S, value) for value in values
        ])

    def state_at_y(self, distance: float) -> np.ndarray:
        """Interpolate the seven ODE variables at a vertical coordinate."""
        values = (
            self.velocity, self.width, self.density, self.mass_fraction,
            self.theta, self.x, self.y,
        )
        return np.array([
            np.interp(distance, self.y, value) for value in values
        ])


class ConservedGaussianJet:
    r"""Axisymmetric Gaussian jet with mass, momentum and energy conservation.

    Velocity, density and fuel mass density follow

    .. math::

       v=v_c e^{-r^2/B^2},\quad
       \rho=\rho_a+(\rho_c-\rho_a)e^{-r^2/(\lambda_h B)^2},\quad
       \rho Y=\rho_cY_c e^{-r^2/(\lambda_Y B)^2}.

    Entrainment is the published source-momentum plus local-buoyancy closure.
    No coefficient in this class is inferred from the Hecht--Panda data.
    """

    def __init__(
        self,
        source: AxisymmetricJetSource,
        *,
        ambient_temperature: float,
        ambient_pressure: float,
        ambient_density: float,
        fuel_molecular_weight: float,
        ambient_molecular_weight: float,
        fuel_heat_capacity: float,
        ambient_heat_capacity: float,
        spreading_ratio: float = 1.16,
        thermodynamic_spreading_ratio: float | None = None,
        thermodynamic_profile: str = "density",
        ideal_gas_enthalpy_fluids: tuple[str, str] | None = None,
        momentum_entrainment_beta: float = 0.28,
        plume_entrainment_limit: float = 0.082,
        radial_limit: float = 5.0,
        radial_points: int = 241,
        conservative_establishment: bool | str = False,
        equilibrium_air_condensation: bool = False,
        equilibrium_dry_air_condensation: bool = True,
        ambient_absolute_humidity: float = 0.0,
        ambient_coflow_velocity: float = 0.0,
        temperature_dependent_phase_enthalpy: bool = False,
        hydrogen_enthalpy_species: str = "Hydrogen",
        equilibrium_argon_condensation: bool = False,
        consistent_phase_ambient: bool = False,
        radiative_absorptivity: float = 0.0,
        energy_transport: str = "total",
        hydrogen_nonideal_volume_correction: bool = False,
    ):
        if min(
            source.diameter, source.velocity, source.density,
            source.temperature, ambient_temperature, ambient_pressure,
            ambient_density, fuel_molecular_weight,
            ambient_molecular_weight, fuel_heat_capacity,
            ambient_heat_capacity, spreading_ratio,
        ) <= 0.0:
            raise ValueError("jet source and thermodynamic inputs must be positive")
        if not 0.0 < source.mass_fraction <= 1.0:
            raise ValueError("source mass fraction must lie in (0, 1]")
        if radial_limit < 4.0 or radial_points < 41:
            raise ValueError("radial quadrature must cover at least 4B with 41 points")

        self.source = source
        self.ambient_temperature = float(ambient_temperature)
        self.ambient_pressure = float(ambient_pressure)
        self.ambient_density = float(ambient_density)
        self.fuel_molecular_weight = float(fuel_molecular_weight)
        self.ambient_molecular_weight = float(ambient_molecular_weight)
        self.ambient_absolute_humidity = float(ambient_absolute_humidity)
        if self.ambient_absolute_humidity < 0.0:
            raise ValueError("ambient absolute humidity cannot be negative")
        self.ambient_coflow_velocity = float(ambient_coflow_velocity)
        if self.ambient_coflow_velocity < 0.0:
            raise ValueError("ambient co-flow velocity cannot be negative")
        self._ambient_velocity_x = (
            self.ambient_coflow_velocity * math.cos(source.theta)
        )
        self._ambient_velocity_y = (
            self.ambient_coflow_velocity * math.sin(source.theta)
        )
        self.water_molecular_weight = WMW_MODERN / 1000.0
        self._humid_ambient_molecular_weight = (
            (1.0 + self.ambient_absolute_humidity)
            / (
                1.0 / self.ambient_molecular_weight
                + self.ambient_absolute_humidity / self.water_molecular_weight
            )
        )
        self.fuel_heat_capacity = float(fuel_heat_capacity)
        self.ambient_heat_capacity = float(ambient_heat_capacity)
        self._humid_ambient_heat_capacity = (
            self.ambient_heat_capacity + self.ambient_absolute_humidity * CPW
        ) / (1.0 + self.ambient_absolute_humidity)
        self.equilibrium_argon_condensation = bool(
            equilibrium_argon_condensation
        )
        self.equilibrium_dry_air_condensation = bool(
            equilibrium_dry_air_condensation
        )
        if (
            self.equilibrium_argon_condensation
            and not self.equilibrium_dry_air_condensation
        ):
            raise ValueError(
                "argon condensation requires dry-air phase equilibrium"
            )
        self.radiative_absorptivity = float(radiative_absorptivity)
        if not 0.0 <= self.radiative_absorptivity <= 1.0:
            raise ValueError("radiative absorptivity must lie between 0 and 1")
        if energy_transport not in {"total", "enthalpy"}:
            raise ValueError("energy transport must be 'total' or 'enthalpy'")
        self.energy_transport = energy_transport
        self.hydrogen_nonideal_volume_correction = bool(
            hydrogen_nonideal_volume_correction
        )
        mw_dry, dry_fractions = phase_dry_air_composition(self.equilibrium_argon_condensation)
        (self._dry_nitrogen_mass_fraction, self._dry_oxygen_mass_fraction,
         self._dry_argon_mass_fraction) = dry_fractions
        self.consistent_phase_ambient = bool(consistent_phase_ambient)
        if self.consistent_phase_ambient:
            if not (equilibrium_air_condensation and temperature_dependent_phase_enthalpy):
                raise ValueError("consistent phase ambient requires component phase thermodynamics")
            if not 14.5 <= self.ambient_temperature <= 300.0:
                raise ValueError("consistent phase ambient requires 14.5--300 K")
            if not math.isfinite(self.ambient_absolute_humidity):
                raise ValueError("ambient absolute humidity must be finite")
            # Override caller-supplied Air EOS properties explicitly: the
            # far-field endpoint must obey the same EOS as local phase gas.
            self.ambient_molecular_weight = mw_dry
            self._humid_ambient_molecular_weight = (
                (1.0 + self.ambient_absolute_humidity)
                / (1.0 / mw_dry + self.ambient_absolute_humidity / self.water_molecular_weight)
            )
            self.ambient_density = (
                self.ambient_pressure * self._humid_ambient_molecular_weight
                / (R_UNIVERSAL * self.ambient_temperature)
            )
            phase_table = _air_phase_property_table()
            dry_mass = 1.0 / (1.0 + self.ambient_absolute_humidity)
            for species, mass, mw in zip(
                ("nitrogen", "oxygen", "argon", "water"),
                (*[dry_mass * w for w in dry_fractions],
                 dry_mass * self.ambient_absolute_humidity),
                (0.0280134, 0.0319988, 0.039948, self.water_molecular_weight),
            ):
                partial = self.ambient_pressure * mass / mw * self._humid_ambient_molecular_weight
                saturation = float(np.interp(
                    self.ambient_temperature, phase_table["temperature"],
                    phase_table[f"{species}_saturation"],
                ))
                if partial > saturation * (1.0 + 1e-10):
                    raise ValueError("consistent phase ambient must be an uncondensed gas state")
        self.spreading_ratio = float(spreading_ratio)
        self.thermodynamic_spreading_ratio = float(
            spreading_ratio if thermodynamic_spreading_ratio is None
            else thermodynamic_spreading_ratio
        )
        if self.thermodynamic_spreading_ratio <= 0.0:
            raise ValueError("thermodynamic spreading ratio must be positive")
        self.thermodynamic_profile = str(thermodynamic_profile)
        if self.thermodynamic_profile not in {
            "density", "temperature_mass_fraction"
        }:
            raise ValueError(
                "thermodynamic profile must be 'density' or "
                "'temperature_mass_fraction'"
            )
        self.momentum_entrainment_beta = float(momentum_entrainment_beta)
        self.plume_entrainment_limit = float(plume_entrainment_limit)
        self.radial_limit = float(radial_limit)
        self.radial_points = int(radial_points)
        if conservative_establishment is True:
            establishment = "entrained_mass"
        elif conservative_establishment is False:
            establishment = "published"
        else:
            establishment = str(conservative_establishment)
        if establishment not in {
            "published", "entrained_mass", "source_flux", "scalar_peak",
        }:
            raise ValueError(
                "establishment must be 'published', 'entrained_mass', "
                "'source_flux' or 'scalar_peak'"
            )
        if (
            establishment == "published" and (
                self.thermodynamic_profile != "density"
                or not math.isclose(
                    self.thermodynamic_spreading_ratio,
                    self.spreading_ratio,
                    rel_tol=0.0,
                    abs_tol=1.0e-14,
                )
            )
        ):
            raise ValueError(
                "differential thermodynamic spreading requires a conservative "
                "establishment"
            )
        self.establishment = establishment
        self.equilibrium_air_condensation = bool(equilibrium_air_condensation)
        self.temperature_dependent_phase_enthalpy = bool(
            temperature_dependent_phase_enthalpy
        )
        self.hydrogen_enthalpy_species = str(hydrogen_enthalpy_species)
        if self.hydrogen_enthalpy_species not in {
            "Hydrogen", "ParaHydrogen", "OrthoHydrogen"
        }:
            raise ValueError(
                "hydrogen enthalpy species must be 'Hydrogen', "
                "'ParaHydrogen' or 'OrthoHydrogen'"
            )
        if (
            self.hydrogen_enthalpy_species != "Hydrogen"
            and not self.temperature_dependent_phase_enthalpy
        ):
            raise ValueError(
                "a hydrogen spin-isomer enthalpy requires "
                "temperature-dependent phase enthalpy"
            )
        if (
            self.temperature_dependent_phase_enthalpy
            and not self.equilibrium_air_condensation
        ):
            raise ValueError(
                "temperature-dependent phase enthalpy requires equilibrium "
                "air condensation"
            )
        if (
            self.equilibrium_argon_condensation
            and not self.equilibrium_air_condensation
        ):
            raise ValueError(
                "argon condensation requires equilibrium air condensation"
            )
        if (
            self.equilibrium_air_condensation
            and self.thermodynamic_profile != "density"
            and self.establishment != "entrained_mass"
        ):
            raise ValueError(
                "equilibrium air condensation with an independent thermal "
                "profile requires four-flux conservative establishment"
            )
        self.ideal_gas_enthalpy_fluids = ideal_gas_enthalpy_fluids
        if (
            self.ideal_gas_enthalpy_fluids is not None
            and len(self.ideal_gas_enthalpy_fluids) != 2
        ):
            raise ValueError(
                "ideal-gas enthalpy fluids must contain fuel and ambient names"
            )
        if (
            self.equilibrium_air_condensation
            and self.ideal_gas_enthalpy_fluids is not None
        ):
            raise ValueError(
                "ideal-gas enthalpy tables cannot be combined with condensed "
                "air"
            )
        if (
            self.ambient_absolute_humidity > 0.0
            and not self.equilibrium_air_condensation
        ):
            raise ValueError(
                "humid ambient air currently requires equilibrium condensation"
            )
        self._ideal_enthalpy_table = (
            None if self.ideal_gas_enthalpy_fluids is None
            else _ideal_gas_enthalpy_table(*self.ideal_gas_enthalpy_fluids)
        )
        self._phase_enthalpy_tables = None
        if self.temperature_dependent_phase_enthalpy:
            component_tables = {
                species: _ideal_component_enthalpy_table(species)[1]
                for species in ("Nitrogen", "Oxygen", "Argon", "Water")
            }
            component_tables["Hydrogen"] = _ideal_component_enthalpy_table(
                self.hydrogen_enthalpy_species
            )[1]
            self._phase_enthalpy_tables = (
                _ideal_component_enthalpy_table(
                    self.hydrogen_enthalpy_species
                )[0],
                component_tables,
            )
        self._source_boundary_enthalpy = None
        if source.enthalpy_boundary is not None:
            boundary = source.enthalpy_boundary
            if self.establishment == "published":
                raise ValueError("phase source enthalpy requires conservative establishment")
            if self._phase_enthalpy_tables is None or not self.equilibrium_air_condensation:
                raise ValueError("phase source enthalpy requires component phase thermodynamics")
            if not all(math.isfinite(value) for value in (
                boundary.specific_enthalpy, boundary.ambient_temperature,
                boundary.hydrogen_reference_enthalpy,
            )):
                raise ValueError("source enthalpy boundary must be finite")
            if not math.isclose(boundary.ambient_temperature, self.ambient_temperature,
                                rel_tol=0.0, abs_tol=1e-8):
                raise ValueError("source enthalpy ambient reference does not match the jet")
            if boundary.hydrogen_species != self.hydrogen_enthalpy_species:
                raise ValueError("source enthalpy hydrogen species does not match the jet")
            h_grid, h_tables = self._phase_enthalpy_tables
            self._source_boundary_enthalpy = boundary.specific_enthalpy + source.mass_fraction * (
                boundary.hydrogen_reference_enthalpy
                - float(np.interp(self.ambient_temperature, h_grid, h_tables["Hydrogen"]))
            )
        self._condensed_lookup = None
        if self.equilibrium_air_condensation:
            # Build once per process; subsequent sources share the immutable
            # component-property grid.
            _air_phase_property_table()
            key = (
                self.ambient_temperature, self.ambient_pressure,
                self.ambient_density, self.fuel_molecular_weight,
                self.fuel_heat_capacity, self.ambient_heat_capacity,
                self.ambient_absolute_humidity,
                float(self.temperature_dependent_phase_enthalpy),
                float(self.equilibrium_argon_condensation),
                float(self.equilibrium_dry_air_condensation),
                self.hydrogen_enthalpy_species,
                self.consistent_phase_ambient,
            )
            if key not in _CONDENSED_LOOKUP_CACHE:
                fraction_grid = np.linspace(0.0, 1.0, 161)
                ideal_temperature_grid = np.linspace(14.5, 300.0, 201)
                if self.consistent_phase_ambient:
                    # The exactly known ambient endpoint must be a node,
                    # not an interpolation across two different cp values.
                    ideal_temperature_grid = np.unique(np.r_[
                        ideal_temperature_grid, self.ambient_temperature,
                    ])
                ideal_temperature_mesh, fraction_mesh = np.meshgrid(
                    ideal_temperature_grid, fraction_grid, indexing="ij"
                )
                mw_mesh = 1.0 / (
                    fraction_mesh / self.fuel_molecular_weight
                    + (1.0 - fraction_mesh)
                    / self._humid_ambient_molecular_weight
                )
                rho_mesh = (
                    self.ambient_pressure * mw_mesh
                    / (R_UNIVERSAL * ideal_temperature_mesh)
                )
                temperature, rho_h = self._condensed_air_state_exact(
                    rho_mesh.ravel(), fraction_mesh.ravel()
                )
                _CONDENSED_LOOKUP_CACHE[key] = (
                    ideal_temperature_grid, fraction_grid,
                    temperature.reshape(rho_mesh.shape),
                    rho_h.reshape(rho_mesh.shape),
                )
            ideal_temperature_grid, fraction_grid, temperature, rho_h = (
                _CONDENSED_LOOKUP_CACHE[key]
            )
            self._condensed_lookup = (
                ideal_temperature_grid,
                fraction_grid,
                RegularGridInterpolator(
                    (ideal_temperature_grid, fraction_grid), temperature,
                    bounds_error=False,
                ),
                RegularGridInterpolator(
                    (ideal_temperature_grid, fraction_grid), rho_h,
                    bounds_error=False,
                ),
            )
        self.establishment_residuals: dict[str, float] | None = None

        self._eta = np.linspace(0.0, self.radial_limit, self.radial_points)
        self._ambient_enthalpy = float(self._mixture_enthalpy(
            self.ambient_temperature, 0.0
        ))
        source_relative_velocity = abs(
            source.velocity - self.ambient_coflow_velocity
        )
        self._source_froude = self._densimetric_froude(
            source_relative_velocity, source.diameter, source.density
        )
        self._buoyancy_coefficient = self._source_buoyancy_coefficient(
            self._source_froude
        )
        self._momentum_entrainment = (
            self.momentum_entrainment_beta
            * math.sqrt(
                source.area * source.density * source_relative_velocity**2
                / self.ambient_density
            )
        )

    @staticmethod
    def _source_buoyancy_coefficient(froude: float) -> float:
        if froude < 268.0:
            return 17.313 - 0.11665 * froude + 2.0771e-4 * froude**2
        return 0.97

    def _densimetric_froude(
        self, velocity: float, width: float, density: float
    ) -> float:
        contrast = abs(self.ambient_density - density)
        if contrast <= 1.0e-15:
            return math.inf
        return velocity / math.sqrt(
            GRAVITY * width * contrast / max(density, 1.0e-30)
        )

    def _mixture_molar_weight(self, fuel_fraction: float | np.ndarray) -> float | np.ndarray:
        """Return mixture molar mass from mass fraction for local closure equations."""
        fraction = np.asarray(fuel_fraction)
        return 1.0 / (
            fraction / self.fuel_molecular_weight
            + (1.0 - fraction) / self._humid_ambient_molecular_weight
        )

    def _fuel_mole_fraction(self, fuel_fraction: float | np.ndarray) -> float | np.ndarray:
        """Return hydrogen mole fraction from local hydrogen mass fraction."""
        fraction = np.asarray(fuel_fraction)
        ambient_mass = 1.0 - fraction
        denominator = (
            fraction / self.fuel_molecular_weight
            + ambient_mass / self._humid_ambient_molecular_weight
        )
        return np.where(
            np.abs(denominator) > 0.0,
            fraction / self.fuel_molecular_weight / denominator,
            0.0,
        )

    def _hydrogen_gas_compressibility(
        self, temperature: float | np.ndarray, fuel_fraction: float | np.ndarray
    ) -> float | np.ndarray:
        """Return local hydrogen compressibility from CoolProp mixture partial pressure."""
        if not self.hydrogen_nonideal_volume_correction:
            return np.ones_like(np.asarray(temperature, dtype=float))

        temperature = np.asarray(temperature, dtype=float)
        mole_fraction = np.asarray(self._fuel_mole_fraction(fuel_fraction), dtype=float)
        temperature, mole_fraction = np.broadcast_arrays(temperature, mole_fraction)
        flat_temperature = temperature.ravel()
        flat_mole_fraction = mole_fraction.ravel()

        compressibility = np.ones_like(flat_temperature)
        for index, (value, fraction) in enumerate(zip(
            flat_temperature, flat_mole_fraction
        )):
            if not np.isfinite(value) or fraction <= 1.0e-14:
                continue
            pressure = self.ambient_pressure * fraction
            try:
                result = hydrogen_gas_departure(
                    float(value), float(pressure),
                    self.hydrogen_enthalpy_species,
                )
            except (ValueError, RuntimeError):
                continue
            if result.compressibility > 0.0:
                compressibility[index] = result.compressibility
        return compressibility.reshape(temperature.shape)

    def _density_from_temperature(
        self, temperature: float | np.ndarray, fuel_fraction: float | np.ndarray
    ) -> float | np.ndarray:
        """Return gas density from local closure, with optional non-ideal correction."""
        mole_weight = np.asarray(self._mixture_molar_weight(fuel_fraction), dtype=float)
        temperature = np.asarray(temperature, dtype=float)
        mole_weight, temperature = np.broadcast_arrays(mole_weight, temperature)
        compressibility = np.asarray(
            self._hydrogen_gas_compressibility(temperature, fuel_fraction),
            dtype=float,
        )
        return (
            self.ambient_pressure * mole_weight
            / (R_UNIVERSAL * temperature * compressibility)
        ).astype(float)

    def _temperature_from_density(
        self, density: float | np.ndarray, fuel_fraction: float | np.ndarray
    ) -> float | np.ndarray:
        """Invert closure relation for temperature from a local density."""
        density = np.asarray(density, dtype=float)
        if np.any(density <= 0.0):
            raise ValueError("density must be positive for temperature inversion")

        mole_weight = np.asarray(self._mixture_molar_weight(fuel_fraction), dtype=float)
        density = np.asarray(density, dtype=float)
        mole_weight, density = np.broadcast_arrays(mole_weight, density)
        ideal_temperature = (
            self.ambient_pressure * mole_weight / (R_UNIVERSAL * density)
        )
        if not self.hydrogen_nonideal_volume_correction:
            return ideal_temperature

        temperature = np.clip(ideal_temperature, 14.0, 400.0)
        for _ in range(32):
            compressibility = np.asarray(
                self._hydrogen_gas_compressibility(temperature, fuel_fraction)
            )
            target = ideal_temperature / compressibility
            if np.allclose(temperature, target, rtol=1.0e-10, atol=0.0):
                break
            temperature = 0.5 * temperature + 0.5 * target
            temperature = np.clip(temperature, 14.0, 400.0)
        return float(temperature) if temperature.ndim == 0 else temperature

    def _mixture_enthalpy(
        self, temperature: float | np.ndarray, fuel_fraction: float | np.ndarray
    ) -> float | np.ndarray:
        """Return consistently referenced mixture specific enthalpy."""
        if self._phase_enthalpy_tables is not None:
            grid, tables = self._phase_enthalpy_tables
            values = np.asarray(temperature)

            def component(species: str) -> np.ndarray:
                table = tables[species]
                return np.interp(values, grid, table) - np.interp(
                    self.ambient_temperature, grid, table
                )

            fraction = np.asarray(fuel_fraction)
            ambient_mass = 1.0 - fraction
            dry_mass = ambient_mass / (1.0 + self.ambient_absolute_humidity)
            result = (
                fraction * component("Hydrogen")
                + dry_mass * self._dry_nitrogen_mass_fraction
                * component("Nitrogen")
                + dry_mass * self._dry_oxygen_mass_fraction
                * component("Oxygen")
                + dry_mass * self._dry_argon_mass_fraction
                * component("Argon")
                + dry_mass * self.ambient_absolute_humidity
                * component("Water")
            )
        elif self._ideal_enthalpy_table is None:
            cp = self._humid_ambient_heat_capacity + np.asarray(fuel_fraction) * (
                self.fuel_heat_capacity - self._humid_ambient_heat_capacity
            )
            result = cp * np.asarray(temperature)
        else:
            grid, fuel_h, ambient_h = self._ideal_enthalpy_table
            values = np.asarray(temperature)

            def interpolate(table: np.ndarray) -> np.ndarray:
                result = np.interp(values, grid, table)
                low_slope = (table[1] - table[0]) / (grid[1] - grid[0])
                high_slope = (table[-1] - table[-2]) / (grid[-1] - grid[-2])
                result = np.where(
                    values < grid[0],
                    table[0] + low_slope * (values - grid[0]),
                    result,
                )
                return np.where(
                    values > grid[-1],
                    table[-1] + high_slope * (values - grid[-1]),
                    result,
                )

            h_fuel = interpolate(fuel_h)
            h_ambient = interpolate(ambient_h)
            # Independent per-species reference offsets cancel exactly from
            # the balances because fuel and total mass fluxes are conserved.
            # Removing them at ambient T prevents catastrophic cancellation
            # between two unrelated CoolProp reference conventions.
            h_fuel = h_fuel - np.interp(
                self.ambient_temperature, grid, fuel_h
            )
            h_ambient = h_ambient - np.interp(
                self.ambient_temperature, grid, ambient_h
            )
            fraction = np.asarray(fuel_fraction)
            result = fraction * h_fuel + (1.0 - fraction) * h_ambient
        return float(result) if np.ndim(result) == 0 else result

    def _temperature_from_enthalpy(
        self, enthalpy: float, fuel_fraction: float
    ) -> float:
        """Invert the selected monotone caloric relation for an initial guess."""
        if self._phase_enthalpy_tables is not None:
            grid = self._phase_enthalpy_tables[0]
            mixture_h = self._mixture_enthalpy(grid, fuel_fraction)
            if enthalpy < mixture_h[0] or enthalpy > mixture_h[-1]:
                raise RuntimeError(
                    "establishment enthalpy lies outside phase property table"
                )
            return float(np.interp(enthalpy, mixture_h, grid))
        if self._ideal_enthalpy_table is None:
            cp = self._humid_ambient_heat_capacity + fuel_fraction * (
                self.fuel_heat_capacity - self._humid_ambient_heat_capacity
            )
            return enthalpy / cp
        grid, fuel_h, ambient_h = self._ideal_enthalpy_table
        fuel_h = fuel_h - np.interp(self.ambient_temperature, grid, fuel_h)
        ambient_h = ambient_h - np.interp(
            self.ambient_temperature, grid, ambient_h
        )
        mixture_h = fuel_fraction * fuel_h + (1.0 - fuel_fraction) * ambient_h
        if enthalpy < mixture_h[0] or enthalpy > mixture_h[-1]:
            raise RuntimeError("establishment enthalpy lies outside property table")
        return float(np.interp(enthalpy, mixture_h, grid))

    def source_specific_enthalpy(self) -> float:
        """Source enthalpy in this model's reference, retaining a phase ledger."""
        if self._source_boundary_enthalpy is not None:
            return self._source_boundary_enthalpy
        return float(self._mixture_enthalpy(
            self.source.temperature, self.source.mass_fraction
        ))

    def _published_initial_state(self) -> tuple[float, np.ndarray]:
        """Return the algebraic Winters/HyRAM established-flow state."""
        src = self.source
        lam2 = self.spreading_ratio**2
        factor = (lam2 + 1.0) / (2.0 * lam2)
        y_c = factor * src.mass_fraction
        # The published all-gas algebra supplies only a numerical initial
        # guess for a phase source. Its inverse cannot represent condensed
        # latent enthalpy. The conservative solve below uses the explicit
        # source ledger in _established_target_fluxes, without clipping it.
        h_plug = self._mixture_enthalpy(src.temperature, src.mass_fraction)
        h_c = self._ambient_enthalpy + factor * (
            h_plug - self._ambient_enthalpy
        )
        temperature_c = self._temperature_from_enthalpy(h_c, y_c)
        rho_c = self._density_from_temperature(temperature_c, y_c)
        width = src.diameter / math.sqrt(
            2.0 * (2.0 * lam2 + 1.0)
            / (
                lam2 * src.density / self.ambient_density
                + lam2 + 1.0
            )
        )

        froude_squared = self._source_froude**2
        if froude_squared >= 40.0:
            lengths = 6.2
        elif froude_squared >= 5.0:
            lengths = 3.9 + 0.057 * froude_squared
        elif froude_squared >= 1.0:
            lengths = 2.075 + 0.425 * froude_squared
        else:
            lengths = 0.0
        S = lengths * src.diameter
        state = np.array([
            src.velocity, width, rho_c, y_c, src.theta,
            src.x + S * math.cos(src.theta),
            src.y + S * math.sin(src.theta),
        ])
        return S, state

    def _established_target_fluxes(self, distance: float) -> np.ndarray:
        """Mass, species, momentum and energy targets after development."""
        src = self.source
        volume_flux = src.density * src.velocity * src.area
        energy = volume_flux * (
            self.source_specific_enthalpy()
            - self._ambient_enthalpy
        )
        if self.energy_transport == "total":
            energy += (
                0.5 * volume_flux * src.velocity**2
                + 0.5 * self.ambient_density * self._momentum_entrainment
                * distance * self.ambient_coflow_velocity**2
            )
        return np.array([
            src.mass_flow
                + self.ambient_density * self._momentum_entrainment * distance,
            src.mass_fraction * volume_flux,
            src.density * src.velocity**2 * src.area
                + self.ambient_density * self._momentum_entrainment
                * distance * self.ambient_coflow_velocity,
            energy,
        ])

    def established_initial_state(self) -> tuple[float, np.ndarray]:
        """Convert the source plug profile to an established Gaussian state.

        The optional conservative conversion retains the published centreline
        development length and obtains velocity, width, density and fuel
        fraction by matching total mass, species, momentum and total energy.
        ``source_flux`` uses the Station-3 source-plane fluxes for those four
        targets, so the Zone-IV profile conversion itself adds no ambient mass.
        """
        distance, published = self._published_initial_state()
        target_distance = 0.0 if self.establishment == "source_flux" else distance
        targets = self._established_target_fluxes(target_distance)

        if self.establishment == "published":
            fluxes = self._fluxes(published)
            actual = np.array([
                fluxes[0], fluxes[3],
                math.hypot(fluxes[1], fluxes[2]), fluxes[4],
            ])
            residual = (actual - targets) / np.maximum(np.abs(targets), 1.0e-30)
            self.establishment_residuals = dict(zip(
                ("mass", "species", "momentum", "energy"), map(float, residual)
            ))
            return distance, published

        if self.establishment == "scalar_peak":
            invariant_targets = targets[1:]

            def decode_scalar(values: np.ndarray) -> np.ndarray:
                candidate = published.copy()
                candidate[0] = math.exp(float(np.clip(values[0], -50.0, 20.0)))
                candidate[1] = math.exp(float(np.clip(values[1], -50.0, 20.0)))
                candidate[2] = math.exp(float(np.clip(values[2], -50.0, 20.0)))
                return candidate

            initial_scalar = np.log(published[:3])

            def residual_scalar(values: np.ndarray) -> np.ndarray:
                candidate_fluxes = self._fluxes(decode_scalar(values))
                actual = np.array([
                    candidate_fluxes[3],
                    math.hypot(candidate_fluxes[1], candidate_fluxes[2]),
                    candidate_fluxes[4],
                ])
                return (
                    (actual - invariant_targets)
                    / np.maximum(np.abs(invariant_targets), 1.0e-30)
                )

            solution = least_squares(
                residual_scalar, initial_scalar, xtol=1.0e-13,
                ftol=1.0e-13, gtol=1.0e-13, max_nfev=300,
            )
            state = decode_scalar(solution.x)
            errors = residual_scalar(solution.x)
            self.establishment_residuals = dict(zip(
                ("species", "momentum", "energy"), map(float, errors)
            ))
            if not solution.success or np.max(np.abs(errors)) > 1.0e-8:
                raise RuntimeError(
                    "scalar-constrained Gaussian establishment did not close: "
                    f"{self.establishment_residuals}"
                )
            return distance, state

        def decode(values: np.ndarray) -> np.ndarray:
            velocity = math.exp(float(np.clip(values[0], -50.0, 20.0)))
            width = math.exp(float(np.clip(values[1], -50.0, 20.0)))
            density = math.exp(float(np.clip(values[2], -50.0, 20.0)))
            logistic = float(np.clip(values[3], -35.0, 35.0))
            fraction = 1.0 / (1.0 + math.exp(-logistic))
            candidate = published.copy()
            candidate[:4] = (velocity, width, density, fraction)
            return candidate

        fraction0 = min(max(published[3], 1.0e-12), 1.0 - 1.0e-12)
        initial = np.array([
            math.log(published[0]),
            math.log(published[1]),
            math.log(published[2]),
            math.log(fraction0 / (1.0 - fraction0)),
        ])

        def residual(values: np.ndarray) -> np.ndarray:
            candidate_fluxes = self._fluxes(decode(values))
            actual = np.array([
                candidate_fluxes[0],
                candidate_fluxes[3],
                math.hypot(candidate_fluxes[1], candidate_fluxes[2]),
                candidate_fluxes[4],
            ])
            return (actual - targets) / np.maximum(np.abs(targets), 1.0e-30)

        solve_options = {}
        if self._ideal_enthalpy_table is not None:
            # Keep numerical trial states finite without constraining any
            # physically plausible established source solution.
            solve_options["bounds"] = (
                np.array([
                    math.log(1.0e-3),
                    math.log(self.source.diameter * 1.0e-3),
                    math.log(1.0e-3),
                    -20.0,
                ]),
                np.array([
                    math.log(5.0e3),
                    math.log(self.source.diameter * 1.0e2),
                    math.log(1.0e2),
                    20.0,
                ]),
            )
        solution = least_squares(
            residual, initial, xtol=1.0e-13, ftol=1.0e-13,
            gtol=1.0e-13, max_nfev=300, **solve_options,
        )
        state = decode(solution.x)
        errors = residual(solution.x)
        self.establishment_residuals = dict(zip(
            ("mass", "species", "momentum", "energy"), map(float, errors)
        ))
        if not solution.success or np.max(np.abs(errors)) > 1.0e-8:
            raise RuntimeError(
                "conservative Gaussian establishment did not close: "
                f"{self.establishment_residuals}"
            )
        return distance, state

    def _profiles_at_eta(
        self, state: np.ndarray, eta: np.ndarray
    ) -> tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
    ]:
        velocity_c, width, rho_c, y_c = state[:4]
        velocity = velocity_c * np.exp(-eta**2)
        species_shape = np.exp(-(eta / self.spreading_ratio) ** 2)
        thermodynamic_shape = np.exp(
            -(eta / self.thermodynamic_spreading_ratio) ** 2
        )
        if self.thermodynamic_profile == "density":
            fuel_mass_density = rho_c * y_c * species_shape
            density = self.ambient_density + (
                rho_c - self.ambient_density
            ) * thermodynamic_shape
            fuel_fraction = fuel_mass_density / density
            molecular_weight = self._mixture_molar_weight(fuel_fraction)
            temperature = self._temperature_from_density(
                density, fuel_fraction
            )
        else:
            molecular_weight_c = self._mixture_molar_weight(y_c)
            temperature_c = self._temperature_from_density(rho_c, y_c)
            temperature = self.ambient_temperature + (
                temperature_c - self.ambient_temperature
            ) * thermodynamic_shape
            fuel_fraction = y_c * species_shape
            molecular_weight = self._mixture_molar_weight(fuel_fraction)
            density = self._density_from_temperature(temperature, fuel_fraction)
        if self.equilibrium_air_condensation:
            temperature, rho_h = self._condensed_air_state(
                density, fuel_fraction
            )
        else:
            rho_h = density * self._mixture_enthalpy(
                temperature, fuel_fraction
            )
        return (
            velocity, density, fuel_fraction, molecular_weight, temperature,
            rho_h,
        )

    def _profiles(
        self, state: np.ndarray
    ) -> tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
    ]:
        return self._profiles_at_eta(state, self._eta)

    def _condensed_air_state(
        self, density: np.ndarray, fuel_fraction: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Interpolate equilibrium N2/O2 temperature and enthalpy."""
        if self._condensed_lookup is None:
            return self._condensed_air_state_exact(density, fuel_fraction)
        ideal_temperature_grid, fraction_grid, temperature_table, enthalpy_table = (
            self._condensed_lookup
        )
        shape = np.shape(density)
        molecular_weight = 1.0 / (
            fuel_fraction / self.fuel_molecular_weight
            + (1.0 - fuel_fraction)
            / self._humid_ambient_molecular_weight
        )
        ideal_temperature = (
            self.ambient_pressure * molecular_weight
            / (R_UNIVERSAL * density)
        )
        points = np.column_stack((
            np.clip(
                np.ravel(ideal_temperature),
                ideal_temperature_grid[0], ideal_temperature_grid[-1],
            ),
            np.clip(
                np.ravel(fuel_fraction), fraction_grid[0], fraction_grid[-1]
            ),
        ))
        return (
            temperature_table(points).reshape(shape),
            enthalpy_table(points).reshape(shape),
        )

    def _condensed_air_state_exact(
        self, density: np.ndarray, fuel_fraction: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve equilibrium N2/O2/Ar/H2O temperature and enthalpy."""
        table = _air_phase_property_table()
        grid = table["temperature"]
        mw_h, mw_n, mw_o, mw_ar = (
            self.fuel_molecular_weight, 0.0280134, 0.0319988, 0.039948
        )
        mw_w = self.water_molecular_weight
        w_n = self._dry_nitrogen_mass_fraction
        w_o = self._dry_oxygen_mass_fraction
        w_ar = self._dry_argon_mass_fraction

        def phase_at(temperature: np.ndarray):
            p_n = np.interp(
                temperature, grid, table["nitrogen_saturation"]
            )
            p_o = np.interp(
                temperature, grid, table["oxygen_saturation"]
            )
            p_ar = np.interp(temperature, grid, table["argon_saturation"])
            p_w = np.interp(temperature, grid, table["water_saturation"])
            mass_h = np.maximum(fuel_fraction, 1.0e-14)
            ambient_mass = np.maximum(1.0 - fuel_fraction, 0.0)
            dry_mass = ambient_mass / (1.0 + self.ambient_absolute_humidity)
            mass_n = dry_mass * w_n
            mass_o = dry_mass * w_o
            mass_ar = dry_mass * w_ar
            mass_w = dry_mass * self.ambient_absolute_humidity
            mole_h = mass_h / mw_h
            mole_n = mass_n / mw_n
            mole_o = mass_o / mw_o
            mole_ar = mass_ar / mw_ar
            mole_w = mass_w / mw_w

            assigned = np.zeros_like(temperature, dtype=bool)
            gas_n = np.zeros_like(temperature)
            gas_o = np.zeros_like(temperature)
            gas_ar = np.zeros_like(temperature)
            gas_w = np.zeros_like(temperature)
            total_moles = np.zeros_like(temperature)
            for mask in range(16):
                if not self.equilibrium_dry_air_condensation and mask & 7:
                    continue
                cond_n, cond_o, cond_ar, cond_w = (
                    bool(mask & 1), bool(mask & 2), bool(mask & 4),
                    bool(mask & 8),
                )
                q_n = p_n / self.ambient_pressure if cond_n else 0.0
                q_o = p_o / self.ambient_pressure if cond_o else 0.0
                q_ar = p_ar / self.ambient_pressure if cond_ar else 0.0
                q_w = p_w / self.ambient_pressure if cond_w else 0.0
                denominator = 1.0 - q_n - q_o - q_ar - q_w
                valid = denominator > 1.0e-12
                numerator = mole_h
                if not cond_n:
                    numerator = numerator + mole_n
                if not cond_o:
                    numerator = numerator + mole_o
                if not cond_ar:
                    numerator = numerator + mole_ar
                if not cond_w:
                    numerator = numerator + mole_w
                candidate_total = numerator / np.maximum(denominator, 1.0e-30)
                candidate_n = q_n * candidate_total if cond_n else mole_n
                candidate_o = q_o * candidate_total if cond_o else mole_o
                candidate_ar = (
                    q_ar * candidate_total if cond_ar else mole_ar
                )
                candidate_w = q_w * candidate_total if cond_w else mole_w
                tolerance = 1.0e-10 * np.maximum(candidate_total, 1.0)
                if self.equilibrium_dry_air_condensation:
                    if cond_n:
                        valid &= candidate_n <= mole_n + tolerance
                    else:
                        valid &= (
                            candidate_n
                            / np.maximum(candidate_total, 1.0e-30)
                            * self.ambient_pressure
                            <= p_n * (1.0 + 1.0e-10)
                        )
                    if cond_o:
                        valid &= candidate_o <= mole_o + tolerance
                    else:
                        valid &= (
                            candidate_o
                            / np.maximum(candidate_total, 1.0e-30)
                            * self.ambient_pressure
                            <= p_o * (1.0 + 1.0e-10)
                        )
                    if cond_ar:
                        valid &= candidate_ar <= mole_ar + tolerance
                    else:
                        valid &= (
                            candidate_ar
                            / np.maximum(candidate_total, 1.0e-30)
                            * self.ambient_pressure
                            <= p_ar * (1.0 + 1.0e-10)
                        )
                if cond_w:
                    valid &= candidate_w <= mole_w + tolerance
                else:
                    valid &= (
                        candidate_w / np.maximum(candidate_total, 1.0e-30)
                        * self.ambient_pressure <= p_w * (1.0 + 1.0e-10)
                    )
                take = valid & ~assigned
                gas_n[take] = candidate_n[take]
                gas_o[take] = candidate_o[take]
                gas_ar[take] = candidate_ar[take]
                gas_w[take] = candidate_w[take]
                total_moles[take] = candidate_total[take]
                assigned |= take
            if not np.all(assigned):
                raise RuntimeError(
                    "no physical N2/O2/Ar/H2O equilibrium phase split"
                )
            return (
                mass_h, mass_n, mass_o, mass_ar, mass_w,
                gas_n * mw_n, gas_o * mw_o, gas_ar * mw_ar,
                gas_w * mw_w, total_moles,
            )

        def density_at(temperature: np.ndarray):
            (
                mass_h, mass_n, mass_o, mass_ar, mass_w,
                gas_n, gas_o, gas_ar, gas_w, gas_moles,
            ) = phase_at(
                temperature
            )
            cond_n, cond_o, cond_ar, cond_w = (
                mass_n - gas_n, mass_o - gas_o, mass_ar - gas_ar,
                mass_w - gas_w,
            )
            rho_n = np.interp(temperature, grid, table["nitrogen_density"])
            rho_o = np.interp(temperature, grid, table["oxygen_density"])
            rho_ar = np.interp(temperature, grid, table["argon_density"])
            rho_w = np.interp(temperature, grid, table["water_density"])
            specific_volume = (
                gas_moles * R_UNIVERSAL * temperature / self.ambient_pressure
                + cond_n / rho_n + cond_o / rho_o + cond_ar / rho_ar
                + cond_w / rho_w
            )
            return 1.0 / specific_volume

        low = np.full_like(density, 14.1)
        high = np.full_like(density, max(self.ambient_temperature, 300.0))
        low_density, high_density = density_at(low), density_at(high)
        # Matched gas EOS reaches the 300 K boundary exactly; two equivalent
        # molar-mass sums can differ by a few floating-point ulps there.
        roundoff = 8 * np.finfo(float).eps if self.consistent_phase_ambient else 0.0
        if (np.any(density > low_density * (1.0 + roundoff))
                or np.any(high_density > density * (1.0 + roundoff))):
            raise RuntimeError(
                "bulk density lies outside the air-phase temperature bracket "
                f"(relative excess: low={np.max(density / low_density - 1):.3g}, "
                f"high={np.max(high_density / density - 1):.3g})"
            )
        # 32 bisections resolve the 286 K bracket below 7e-8 K, far beneath
        # both property-table resolution and experimental temperature error.
        for _ in range(32):
            middle = 0.5 * (low + high)
            too_dense = density_at(middle) > density
            low = np.where(too_dense, middle, low)
            high = np.where(too_dense, high, middle)
        temperature = 0.5 * (low + high)
        (
            mass_h, mass_n, mass_o, mass_ar, mass_w,
            gas_n, gas_o, gas_ar, gas_w, gas_moles,
        ) = phase_at(temperature)
        cond_n, cond_o, cond_ar, cond_w = (
            mass_n - gas_n, mass_o - gas_o, mass_ar - gas_ar,
            mass_w - gas_w,
        )
        latent_n = np.interp(temperature, grid, table["nitrogen_latent"])
        latent_o = np.interp(temperature, grid, table["oxygen_latent"])
        latent_ar = np.interp(temperature, grid, table["argon_latent"])
        latent_w = np.interp(temperature, grid, table["water_latent"])

        if self._phase_enthalpy_tables is None:
            raw_cp_n = 3.5 * R_UNIVERSAL / mw_n
            raw_cp_o = 3.5 * R_UNIVERSAL / mw_o
            raw_cp_ar = 2.5 * R_UNIVERSAL / mw_ar
            air_scale = self.ambient_heat_capacity / (
                w_n * raw_cp_n + w_o * raw_cp_o + w_ar * raw_cp_ar
            )
            cp_n, cp_o, cp_ar = (
                raw_cp_n * air_scale, raw_cp_o * air_scale,
                raw_cp_ar * air_scale,
            )
            specific_enthalpy = (
                mass_h * self.fuel_heat_capacity * temperature
                + mass_n * cp_n * temperature
                + mass_o * cp_o * temperature
                + mass_ar * cp_ar * temperature
                + mass_w * CPW * temperature
                - cond_n * latent_n - cond_o * latent_o
                - cond_ar * latent_ar - cond_w * latent_w
            )
        else:
            grid_h, tables = self._phase_enthalpy_tables

            def component(species: str) -> np.ndarray:
                table_h = tables[species]
                return np.interp(temperature, grid_h, table_h) - np.interp(
                    self.ambient_temperature, grid_h, table_h
                )

            specific_enthalpy = (
                mass_h * component("Hydrogen")
                + mass_n * component("Nitrogen")
                + mass_o * component("Oxygen")
                + mass_ar * component("Argon")
                + mass_w * component("Water")
                - cond_n * latent_n - cond_o * latent_o
                - cond_ar * latent_ar - cond_w * latent_w
            )
        rho_h = density * specific_enthalpy
        if self.consistent_phase_ambient:
            ambient = (fuel_fraction == 0.0) & np.isclose(
                density, self.ambient_density, rtol=8 * np.finfo(float).eps, atol=0.0,
            )
            # Analytical all-gas endpoint, validated in the constructor.
            # This removes bisection roundoff only; no finite-Y offset is fitted.
            temperature = np.where(ambient, self.ambient_temperature, temperature)
            rho_h = np.where(ambient, 0.0, rho_h)
        return temperature, rho_h

    def _fluxes(self, state: np.ndarray) -> np.ndarray:
        velocity_c, width, _rho_c, _y_c, theta = state[:5]
        if velocity_c <= 0.0 or width <= 0.0:
            raise ValueError("non-positive Gaussian jet velocity or width")
        velocity, density, fuel_fraction, _molecular_weight, _temperature, rho_h = (
            self._profiles(state)
        )
        area_weight = 2.0 * math.pi * width**2 * self._eta
        mass = trapezoid(density * velocity * area_weight, self._eta)
        momentum = trapezoid(
            density * velocity**2 * area_weight, self._eta
        )
        species = trapezoid(
            density * velocity * fuel_fraction * area_weight, self._eta
        )
        energy_density = (
            rho_h - self._ambient_enthalpy * density
        ) * velocity
        if self.energy_transport == "total":
            energy_density = energy_density + 0.5 * density * velocity**3
        energy = trapezoid(energy_density * area_weight, self._eta)
        return np.array([
            mass,
            momentum * math.cos(theta),
            momentum * math.sin(theta),
            species,
            energy,
        ])

    def enthalpy_flux(self, state: np.ndarray) -> float:
        """Return ambient-relative mixture-enthalpy flux for one section.

        This is the fifth balance used by Li et al. (2026), equation 35. The
        default conserved jet continues to use total enthalpy plus resolved
        kinetic energy in :meth:`_fluxes`; exposing this separate integral
        permits a controlled downstream energy-partition comparison.
        """
        velocity, density, _fraction, _mw, _temperature, rho_h = self._profiles(
            state
        )
        area_weight = 2.0 * math.pi * state[1] ** 2 * self._eta
        return float(trapezoid(
            (
                rho_h - self._ambient_enthalpy * density
            ) * velocity * area_weight,
            self._eta,
        ))

    def buoyancy_force(self, state: np.ndarray) -> float:
        """Return the signed section-integrated vertical buoyancy force."""
        _velocity, density, _fraction, _mw, _temperature, _rho_h = (
            self._profiles(state)
        )
        area_weight = 2.0 * math.pi * state[1] ** 2 * self._eta
        return float(trapezoid(
            GRAVITY * (self.ambient_density - density) * area_weight,
            self._eta,
        ))

    def _flux_jacobian(self, state: np.ndarray) -> np.ndarray:
        """Central finite-difference Jacobian of the five conserved fluxes."""
        jacobian = np.empty((5, 5))
        scales = (1.0, 1.0e-3, 1.0, 0.1, 1.0)
        for column in range(5):
            step = 2.0e-6 * max(abs(state[column]), scales[column])
            plus, minus = state.copy(), state.copy()
            plus[column] += step
            minus[column] -= step
            if column in (0, 1, 2, 3) and minus[column] <= 0.0:
                base = self._fluxes(state)
                jacobian[:, column] = (
                    self._fluxes(plus) - base
                ) / step
            else:
                jacobian[:, column] = (
                    self._fluxes(plus) - self._fluxes(minus)
                ) / (2.0 * step)
        return jacobian

    def _entrainment(self, state: np.ndarray) -> float:
        velocity, width, density, _fraction, theta = state[:5]
        relative_velocity = math.sqrt(max(
            velocity**2 + self.ambient_coflow_velocity**2
            - 2.0 * velocity * self.ambient_coflow_velocity
            * math.cos(theta - self.source.theta),
            0.0,
        ))
        contrast = abs(self.ambient_density - density)
        buoyancy = 0.0
        if (
            contrast > 1.0e-15
            and relative_velocity > 1.0e-15
            and math.sin(theta) > 0.0
        ):
            local_froude = (
                relative_velocity**2 * density
                / (GRAVITY * width * contrast)
            )
            buoyancy = (
                self._buoyancy_coefficient / local_froude
                * 2.0 * math.pi * relative_velocity * width * math.sin(theta)
            )
        entrainment = self._momentum_entrainment + buoyancy
        cap = (
            self.plume_entrainment_limit
            * 2.0 * math.pi * width * relative_velocity
        )
        return min(entrainment, cap)

    def derivatives(self, _S: float, state: np.ndarray) -> np.ndarray:
        velocity, width, density, _fraction, theta = state[:5]
        if min(velocity, width, density) <= 0.0:
            raise RuntimeError("axisymmetric jet reached a non-physical state")
        entrainment = self._entrainment(state)
        sources = np.array([
            self.ambient_density * entrainment,
            self.ambient_density * entrainment * self._ambient_velocity_x,
            self.ambient_density * entrainment * self._ambient_velocity_y
                + math.pi * self.spreading_ratio**2 * width**2
                * GRAVITY * (self.ambient_density - density),
            0.0,
            (
                0.5 * self.ambient_density * entrainment
                * self.ambient_coflow_velocity**2
                if self.energy_transport == "total" else 0.0
            ) + self._radiative_power_per_length(state),
        ])
        jacobian = self._flux_jacobian(state)
        physical = np.linalg.solve(jacobian, sources)
        return np.concatenate((physical, [math.cos(theta), math.sin(theta)]))

    def _radiative_power_per_length(self, state: np.ndarray) -> float:
        """Grey absorbed room radiation, W per metre of centreline.

        The full ``radial_limit * B`` envelope is deliberately used as the
        absorbing cylinder. ``radiative_absorptivity=1`` is therefore the
        pre-registered perfect-black upper bound, not a gas emissivity model.
        """
        if self.radiative_absorptivity == 0.0:
            return 0.0
        sigma = 5.670374419e-8
        width = float(state[1])
        temperature = self.centreline_temperature(state)
        return float(
            self.radiative_absorptivity
            * 2.0 * math.pi * self.radial_limit * width * sigma
            * max(self.ambient_temperature**4 - temperature**4, 0.0)
        )

    def centreline_temperature(self, state: np.ndarray) -> float:
        """Ideal-gas mixture temperature corresponding to a centreline state."""
        return float(self._profiles(state)[4][0])

    def radial_state(
        self, state: np.ndarray, radius: float
    ) -> tuple[float, float, float]:
        """Return ``(mass_fraction, density, temperature)`` at one radius."""
        width = float(state[1])
        eta = abs(float(radius)) / width
        _v, rho, fraction, _mw, temperature, _rho_h = self._profiles_at_eta(
            state, np.array([eta])
        )
        return float(fraction[0]), float(rho[0]), float(temperature[0])

    def half_width(self, state: np.ndarray, variable: str) -> float:
        """Radius where fuel or temperature excursion is half its centre value."""
        centre_y, _centre_rho, centre_t = self.radial_state(state, 0.0)
        if variable == "mass":
            target = 0.5 * centre_y
            index = 0
        elif variable == "temperature":
            target = 0.5 * (centre_t + self.ambient_temperature)
            index = 2
        else:
            raise ValueError("variable must be 'mass' or 'temperature'")

        def residual(radius: float) -> float:
            return self.radial_state(state, radius)[index] - target

        return float(brentq(residual, 0.0, 20.0 * state[1]))

    def solve(
        self,
        *,
        maximum_distance: float,
        minimum_mass_fraction: float = 7.0e-4,
        maximum_step: float | None = None,
        relative_tolerance: float = 2.0e-7,
        method: str = "RK45",
    ) -> AxisymmetricJetResult:
        """Integrate the established-flow equations along the centreline."""
        S0, state0 = self.established_initial_state()
        if maximum_distance <= S0:
            raise ValueError("maximum distance must exceed the development zone")
        if maximum_step is None:
            maximum_step = min(0.5 * self.source.diameter, maximum_distance / 200.0)

        def dilute(_S: float, state: np.ndarray) -> float:
            return state[3] - minimum_mass_fraction

        dilute.terminal = True
        dilute.direction = -1.0
        solution = solve_ivp(
            self.derivatives,
            (S0, maximum_distance),
            state0,
            method=method,
            rtol=relative_tolerance,
            atol=np.array([
                1.0e-7, 1.0e-11, 1.0e-10, 1.0e-10,
                1.0e-10, 1.0e-10, 1.0e-10,
            ]),
            max_step=maximum_step,
            events=dilute,
        )
        if not solution.success:
            raise RuntimeError(f"axisymmetric jet integration failed: {solution.message}")

        states = solution.y.T
        fluxes = np.array([self._fluxes(state) for state in states])
        temperatures = np.array([
            self.centreline_temperature(state) for state in states
        ])
        momentum = np.hypot(fluxes[:, 1], fluxes[:, 2])
        radiative_power = np.array([
            self._radiative_power_per_length(state) for state in states
        ])
        radiative_heat = cumulative_trapezoid(
            radiative_power, solution.t, initial=0.0
        )
        return AxisymmetricJetResult(
            S=solution.t,
            velocity=states[:, 0],
            width=states[:, 1],
            density=states[:, 2],
            mass_fraction=states[:, 3],
            theta=states[:, 4],
            x=states[:, 5],
            y=states[:, 6],
            temperature=temperatures,
            mass_flux=fluxes[:, 0],
            momentum_flux=momentum,
            species_flux=fluxes[:, 3],
            energy_flux=fluxes[:, 4],
            radiative_heat_added=radiative_heat,
        )


__all__ = [
    "AxisymmetricJetSource", "AxisymmetricJetResult", "ConservedGaussianJet",
    "InitialEntrainmentResult", "entrain_and_heat_initial_plug",
]
