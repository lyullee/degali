"""The notional nozzle: where a flashing jet has finished expanding.

A liquefied gas discharged through an orifice does not leave at the orifice
diameter.  It leaves choked, at a pressure above ambient, and then expands --
flashing further as it does -- to a plane where the pressure has equalised.
That plane is much wider and slower than the orifice, and it is the correct
starting point for a plume model.

Skipping it is not a small error.  Starting the Ooms integration at the
orifice diameter puts a 25 mm release three centimetres wide at 0.35 m
downstream and predicts 100 mole per cent where the PRESLHY near-field array
measured 75.  The jet's own entrainment cannot recover from a source that
narrow.

What is conserved
-----------------
Between the orifice plane (1) and the notional nozzle (2), mass and momentum:

.. math::
    \\dot m = \\rho_1 u_1 A_1 = \\rho_2 u_2 A_2, \\qquad
    \\dot m u_2 = \\dot m u_1 + (P_1 - P_a) A_1

so the pressure the jet is still carrying at the orifice is converted into
velocity, and the area follows from continuity once the density at ambient
pressure is known.  This is Birch and co-workers' 1987 construction, which
conserves momentum as well as mass; their 1984 version and the mass-only
form of Ewan and Moodie are also in use and give effective diameters
differing by tens of per cent.

For a flashing liquid the density at the notional nozzle is the homogeneous
two-phase density after an isenthalpic flash to ambient, which is what
:attr:`~degali.validation.flashing.FlashResult.orifice_density` already
computes.

``isentropic_throat`` and ``energy_conserving_notional_nozzle`` implement the
HyRAM+ 6.0 section-3.2 sequence.  They supersede the older convenience
function below for research-source work: the throat pressure maximises the
homogeneous-equilibrium mass flux, and the atmospheric plane closes mass,
momentum *and* total energy.  The measured mass flow is retained by inferring
the discharge coefficient rather than replacing the experiment with the
ideal-nozzle flow.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..core.constants import ATM_TO_PA, PI


@dataclass(frozen=True)
class NotionalNozzle:
    """The expanded jet, at ambient pressure."""

    diameter: float  #: effective diameter, m
    velocity: float  #: m/s
    density: float  #: kg/m**3
    orifice_velocity: float  #: m/s, for comparison
    orifice_pressure: float  #: Pa, absolute
    orifice_diameter: float  #: m, as supplied

    @property
    def expansion_ratio(self) -> float:
        """Effective diameter over orifice diameter."""
        return self.diameter / self.orifice_diameter


@dataclass(frozen=True)
class IsentropicThroat:
    """Maximum-mass-flux state on the storage isentrope."""

    storage_temperature: float
    storage_pressure: float
    storage_enthalpy: float
    storage_entropy: float
    pressure: float
    temperature: float
    density: float
    enthalpy: float
    velocity: float
    mass_flux: float
    choked: bool
    relative_energy_residual: float

    @property
    def pressure_ratio(self) -> float:
        """Throat pressure divided by storage pressure."""
        return self.pressure / self.storage_pressure


@dataclass(frozen=True)
class EnergyConservingNotionalNozzle:
    """Yuceil--Otugen atmospheric plane based on a measured mass flow."""

    throat: IsentropicThroat
    mass_flow: float
    orifice_diameter: float
    discharge_coefficient: float
    pressure_thrust: float
    advective_momentum: float
    momentum_flux: float
    velocity: float
    enthalpy: float
    temperature: float | None
    density: float | None
    diameter: float | None
    relative_mass_residual: float | None
    relative_momentum_residual: float
    relative_energy_residual: float | None

    @property
    def specific_momentum(self) -> float:
        """Axial momentum flux per unit measured H2 mass flow, m/s."""
        return self.momentum_flux / self.mass_flow

    @property
    def pressure_to_advective_momentum(self) -> float:
        """Pressure thrust divided by the throat advective term."""
        return self.pressure_thrust / max(self.advective_momentum, 1.0e-30)

    @property
    def admissible(self) -> bool:
        """Whether measured flow and atmospheric energy state are physical."""
        return (
            0.0 < self.discharge_coefficient <= 1.0
            and self.temperature is not None
            and self.density is not None
            and self.diameter is not None
        )


@dataclass(frozen=True)
class MeasuredThroatExpansion:
    """Atmospheric source obtained from a supplied choked-gas state.

    Unlike :class:`EnergyConservingNotionalNozzle`, this result does not infer
    the throat from a saturated storage state.  It is for experiments that
    report the throat pressure, temperature, density and velocity directly.
    """

    mass_flow: float
    throat_diameter: float
    throat_pressure: float
    throat_temperature: float
    throat_density: float
    throat_velocity: float
    pressure_thrust: float
    momentum_flux: float
    total_specific_energy: float
    velocity: float
    enthalpy: float
    temperature: float
    density: float
    diameter: float
    relative_input_mass_residual: float
    relative_mass_residual: float
    relative_momentum_residual: float
    relative_energy_residual: float


def expand_measured_throat_to_ambient(
    *, fluid: str, mass_flow: float, throat_diameter: float,
    throat_pressure: float, throat_temperature: float,
    throat_density: float, throat_velocity: float,
    ambient_pressure: float = ATM_TO_PA,
) -> MeasuredThroatExpansion:
    """Apply HyRAM+ equations 47--51 to a supplied throat plane.

    Mass fixes the atmospheric area, while throat advective momentum plus
    positive pressure thrust fixes velocity.  Static enthalpy at the
    atmospheric plane follows from conservation of total specific energy,
    ``h + u**2/2``.  No storage-state or ideal choking calculation is added.
    """
    from CoolProp.CoolProp import PropsSI

    values = (
        mass_flow, throat_diameter, throat_pressure, throat_temperature,
        throat_density, throat_velocity, ambient_pressure,
    )
    if any(value <= 0.0 for value in values):
        raise ValueError("throat expansion inputs must be positive")

    throat_area = PI * throat_diameter**2 / 4.0
    input_mass_flow = throat_density * throat_velocity * throat_area
    input_mass_residual = abs(input_mass_flow - mass_flow) / mass_flow
    if input_mass_residual > 0.05:
        raise ValueError(
            "reported throat state is inconsistent with mass flow "
            f"({100.0 * input_mass_residual:.1f}% residual)"
        )

    pressure_thrust = max(throat_pressure - ambient_pressure, 0.0) * throat_area
    advective_momentum = mass_flow * throat_velocity
    momentum_flux = advective_momentum + pressure_thrust
    velocity = momentum_flux / mass_flow

    throat_enthalpy = float(PropsSI(
        "H", "T", throat_temperature, "P", throat_pressure, fluid
    ))
    total_specific_energy = throat_enthalpy + 0.5 * throat_velocity**2
    enthalpy = total_specific_energy - 0.5 * velocity**2
    temperature = float(PropsSI(
        "T", "P", ambient_pressure, "H", enthalpy, fluid
    ))
    density = float(PropsSI(
        "D", "P", ambient_pressure, "H", enthalpy, fluid
    ))
    area = mass_flow / (density * velocity)
    diameter = math.sqrt(4.0 * area / PI)

    mass_residual = abs(density * velocity * area - mass_flow) / mass_flow
    momentum_residual = abs(
        mass_flow * velocity - momentum_flux
    ) / max(abs(momentum_flux), 1.0)
    energy_residual = abs(
        enthalpy + 0.5 * velocity**2 - total_specific_energy
    ) / max(abs(total_specific_energy), 1.0)
    return MeasuredThroatExpansion(
        mass_flow=mass_flow,
        throat_diameter=throat_diameter,
        throat_pressure=throat_pressure,
        throat_temperature=throat_temperature,
        throat_density=throat_density,
        throat_velocity=throat_velocity,
        pressure_thrust=pressure_thrust,
        momentum_flux=momentum_flux,
        total_specific_energy=total_specific_energy,
        velocity=velocity,
        enthalpy=enthalpy,
        temperature=temperature,
        density=density,
        diameter=diameter,
        relative_input_mass_residual=input_mass_residual,
        relative_mass_residual=mass_residual,
        relative_momentum_residual=momentum_residual,
        relative_energy_residual=energy_residual,
    )


def isentropic_throat(
    *, fluid: str, storage_temperature: float,
    ambient_pressure: float = ATM_TO_PA,
    storage_pressure: float | None = None,
) -> IsentropicThroat:
    """Find the homogeneous-equilibrium critical state.

    By default the storage state is saturated liquid.  ``storage_pressure``
    supplies an independently measured absolute pressure for a compressed or
    subcooled liquid; this is essential when caloric temperature and driving
    pressure do not lie on the saturation curve.  At every pressure on the
    resulting isentrope, velocity follows ``u**2/2 + h = h0``.  The returned
    pressure maximises ``rho*u`` between ambient and storage pressure,
    exactly as in HyRAM+ equations 44--45.  Optimising log-pressure avoids
    loss of resolution near the lower bound.
    """
    from CoolProp.CoolProp import PropsSI
    from scipy.optimize import minimize_scalar

    if storage_temperature <= 0.0 or ambient_pressure <= 0.0:
        raise ValueError("temperature and pressure must be positive")
    critical_temperature = float(PropsSI("Tcrit", fluid))
    if storage_temperature >= critical_temperature:
        raise ValueError(
            f"{fluid} is supercritical at {storage_temperature:.1f} K"
        )
    if storage_pressure is None:
        storage_pressure = float(
            PropsSI("P", "T", storage_temperature, "Q", 0, fluid)
        )
        h0 = float(PropsSI("H", "T", storage_temperature, "Q", 0, fluid))
        s0 = float(PropsSI("S", "T", storage_temperature, "Q", 0, fluid))
    else:
        storage_pressure = float(storage_pressure)
        if storage_pressure <= 0.0:
            raise ValueError("storage pressure must be positive")
        h0 = float(
            PropsSI("H", "T", storage_temperature, "P", storage_pressure, fluid)
        )
        s0 = float(
            PropsSI("S", "T", storage_temperature, "P", storage_pressure, fluid)
        )
    if storage_pressure <= ambient_pressure:
        raise ValueError("storage pressure must exceed ambient pressure")

    def state(pressure: float) -> tuple[float, float, float, float, float]:
        enthalpy = float(PropsSI("H", "P", pressure, "S", s0, fluid))
        density = float(PropsSI("D", "P", pressure, "S", s0, fluid))
        temperature = float(PropsSI("T", "P", pressure, "S", s0, fluid))
        velocity = math.sqrt(max(2.0 * (h0 - enthalpy), 0.0))
        return density * velocity, enthalpy, density, temperature, velocity

    low, high = math.log(ambient_pressure), math.log(storage_pressure)
    optimum = minimize_scalar(
        lambda log_pressure: -state(math.exp(log_pressure))[0],
        bounds=(low, high), method="bounded", options={"xatol": 1.0e-12},
    )
    candidates = (low, float(optimum.x), high)
    log_pressure = max(candidates, key=lambda value: state(math.exp(value))[0])
    pressure = math.exp(log_pressure)
    mass_flux, enthalpy, density, temperature, velocity = state(pressure)
    energy_residual = abs(h0 - enthalpy - 0.5 * velocity**2) / max(
        abs(h0), abs(enthalpy), 1.0
    )
    return IsentropicThroat(
        storage_temperature=storage_temperature,
        storage_pressure=storage_pressure,
        storage_enthalpy=h0,
        storage_entropy=s0,
        pressure=pressure,
        temperature=temperature,
        density=density,
        enthalpy=enthalpy,
        velocity=velocity,
        mass_flux=mass_flux,
        choked=pressure > ambient_pressure + max(0.01, 1.0e-8 * ambient_pressure),
        relative_energy_residual=energy_residual,
    )


def energy_conserving_notional_nozzle(
    *, fluid: str, storage_temperature: float, mass_flow: float,
    orifice_diameter: float, ambient_pressure: float = ATM_TO_PA,
    storage_pressure: float | None = None,
) -> EnergyConservingNotionalNozzle:
    """Expand the HEM throat to ambient pressure without creating energy.

    This is the Yuceil--Otugen option in HyRAM+.  With the HyRAM discharge-
    coefficient convention, the throat advective momentum is
    ``m_dot*u_t*Cd`` and pressure thrust acts on the geometric orifice area.
    The ambient-pressure thermodynamic state is solved from total enthalpy.
    An impossible state is returned with ``None`` state fields so callers can
    report the source incompatibility rather than silently clip it.
    """
    from CoolProp.CoolProp import PropsSI

    if mass_flow <= 0.0 or orifice_diameter <= 0.0:
        raise ValueError("mass flow and orifice diameter must be positive")
    throat = isentropic_throat(
        fluid=fluid, storage_temperature=storage_temperature,
        ambient_pressure=ambient_pressure,
        storage_pressure=storage_pressure,
    )
    area = PI * orifice_diameter**2 / 4.0
    discharge_coefficient = mass_flow / max(throat.mass_flux * area, 1.0e-30)
    advective = mass_flow * throat.velocity * discharge_coefficient
    pressure_thrust = max(throat.pressure - ambient_pressure, 0.0) * area
    momentum = advective + pressure_thrust
    velocity = momentum / mass_flow
    enthalpy = throat.storage_enthalpy - 0.5 * velocity**2
    momentum_residual = abs(
        mass_flow * velocity - advective - pressure_thrust
    ) / max(abs(momentum), 1.0)

    temperature = density = diameter = None
    mass_residual = energy_residual = None
    try:
        temperature = float(
            PropsSI("T", "P", ambient_pressure, "H", enthalpy, fluid)
        )
        density = float(
            PropsSI("D", "P", ambient_pressure, "H", enthalpy, fluid)
        )
        effective_area = mass_flow / (density * velocity)
        diameter = math.sqrt(4.0 * effective_area / PI)
        mass_residual = abs(
            density * velocity * effective_area - mass_flow
        ) / mass_flow
        energy_residual = abs(
            enthalpy + 0.5 * velocity**2 - throat.storage_enthalpy
        ) / max(abs(throat.storage_enthalpy), 1.0)
    except (ValueError, OverflowError):
        # Some measured-flow/geometric-area combinations make pressure work
        # demand more energy than the storage stream carries.  That is a
        # source-model incompatibility, not a state to clip into existence.
        pass

    return EnergyConservingNotionalNozzle(
        throat=throat,
        mass_flow=mass_flow,
        orifice_diameter=orifice_diameter,
        discharge_coefficient=discharge_coefficient,
        pressure_thrust=pressure_thrust,
        advective_momentum=advective,
        momentum_flux=momentum,
        velocity=velocity,
        enthalpy=enthalpy,
        temperature=temperature,
        density=density,
        diameter=diameter,
        relative_mass_residual=mass_residual,
        relative_momentum_residual=momentum_residual,
        relative_energy_residual=energy_residual,
    )


def notional_nozzle(
    *,
    mass_flow: float,
    orifice_diameter: float,
    orifice_density: float,
    expanded_density: float,
    orifice_pressure: float,
    ambient_pressure: float = ATM_TO_PA,
) -> NotionalNozzle:
    """Expand a choked release to ambient pressure.

    Parameters
    ----------
    mass_flow
        kg/s of released material, excluding entrained air.
    orifice_diameter
        m.
    orifice_density
        Density of the fluid *at* the orifice, kg/m**3.  For a flashing liquid
        that is close to the saturated liquid density: little has flashed yet
        at the orifice pressure.
    expanded_density
        Density once expanded to ambient, kg/m**3 -- for a flashing release,
        the homogeneous two-phase density after the flash.
    orifice_pressure
        Absolute pressure at the orifice plane, Pa.  Values at or below
        ambient mean the jet is not under-expanded and the notional nozzle
        reduces to the orifice.

    Returns
    -------
    NotionalNozzle
        With ``diameter`` at least the orifice diameter: an expansion cannot
        make the jet narrower.
    """
    area = PI * orifice_diameter**2 / 4.0
    u1 = mass_flow / max(orifice_density * area, 1e-30)
    excess = max(orifice_pressure - ambient_pressure, 0.0)
    u2 = u1 + excess * area / max(mass_flow, 1e-30)
    a2 = mass_flow / max(expanded_density * u2, 1e-30)
    d2 = max(math.sqrt(4.0 * a2 / PI), orifice_diameter)

    return NotionalNozzle(
        diameter=d2, velocity=u2, density=expanded_density,
        orifice_velocity=u1, orifice_pressure=orifice_pressure,
        orifice_diameter=orifice_diameter,
    )


def choking_pressure(fluid: str, storage_temperature: float) -> float:
    """Critical pressure for a saturated liquid through a sharp orifice, Pa.

    The homogeneous equilibrium model puts it close to the saturation pressure
    of the upstream fluid, which is what is returned.  Above the critical
    temperature there is no saturation pressure and the fluid is not a
    liquefied gas at these conditions.
    """
    from CoolProp.CoolProp import PropsSI

    if storage_temperature >= float(PropsSI("Tcrit", fluid)):
        raise ValueError(
            f"{fluid} is supercritical at {storage_temperature:.1f} K"
        )
    return float(PropsSI("P", "T", storage_temperature, "Q", 0, fluid))
