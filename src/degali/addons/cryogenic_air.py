"""Cryogenic-air phase checks and the Li et al. (2026) Zone-III model.

This module deliberately separates two things that the paper combines:

* :func:`li2026_zone3` reproduces equations (13)--(24), including the
  paper's nitrogen-only and energy-bookkeeping assumptions.
* :func:`air_saturation_pressure` supplies a phase-correct solid-vapour
  pressure below the triple point.  CoolProp's liquid saturation ancillary
  can return a number there, but that number is a metastable extrapolation,
  not a solid--vapour equilibrium pressure.

The separation matters for liquid-hydrogen releases.  The four cold tests in
Li et al. have Station-2 temperatures of 51--55 K, already below nitrogen's
63.15 K triple point.  PRESLHY saturated-liquid releases are colder still.
Consequently the paper model is useful as a reproducible comparison, but its
``liquid nitrogen`` closure is not silently treated as valid LH2 physics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from scipy.optimize import brentq


# Stable solid-vapour data, stored as (K, Pa).  N2 is Table 7-11 of NBS
# Circular 564.  O2 is Aoyama and Kanda's solid data as transcribed in Georgia
# Tech report A-593; the visibly non-monotone 41.58 K point is omitted in
# favour of the neighbouring 41.91 K value.  Interpolation is linear in
# ln(P) against 1/T, the Clausius-Clapeyron coordinates.  Below the lowest
# datum the lowest-temperature segment is extrapolated: a constant-latent-
# heat approximation which, unlike a liquid ancillary, preserves P -> 0.
_MMHG_TO_PA = 133.322368
_MOLECULAR_WEIGHT = {
    "Hydrogen": 0.00201588,
    "Nitrogen": 0.0280134,
    "Oxygen": 0.0319988,
}
_HYDROGEN_SPECIES = frozenset({
    "Hydrogen", "ParaHydrogen", "OrthoHydrogen",
})
_TRIPLE_TEMPERATURE = {
    "Hydrogen": 13.957,
    "Nitrogen": 63.151,
    "Oxygen": 54.361,
}
_CRITICAL_TEMPERATURE = {"Nitrogen": 126.192, "Oxygen": 154.581}
_SOLID_VAPOUR_TABLES = {
    "Nitrogen": tuple((t, p * _MMHG_TO_PA) for t, p in (
        (52.0, 5.7), (54.0, 10.2), (56.0, 17.6), (58.0, 29.4),
        (60.0, 47.2), (63.156, 94.0),
    )),
    "Oxygen": tuple((t, p * _MMHG_TO_PA) for t, p in (
        (36.0, 0.00006), (37.62, 0.00020), (39.41, 0.00092),
        (41.91, 0.0034), (43.09, 0.0101), (44.11, 0.022),
        (46.09, 0.045), (49.38, 0.178), (50.74, 0.291),
        (53.02, 0.708), (54.36, 1.20),
    )),
}


def _solid_vapour_pressure(species: str, temperature: float) -> float:
    table = _SOLID_VAPOUR_TABLES[species]
    upper = 1
    while upper < len(table) and temperature > table[upper][0]:
        upper += 1
    if upper == len(table):
        lower, upper = len(table) - 2, len(table) - 1
    elif upper == 0:  # unreachable with upper initialised to one
        lower, upper = 0, 1
    else:
        lower = upper - 1
    t0, p0 = table[lower]
    t1, p1 = table[upper]
    inverse_temperature = 1.0 / temperature
    weight = (inverse_temperature - 1.0 / t0) / (1.0 / t1 - 1.0 / t0)
    return float(math.exp(math.log(p0) + weight * (math.log(p1) - math.log(p0))))


@lru_cache(maxsize=8192)
def air_saturation_pressure(species: str, temperature: float) -> float:
    """Return the stable condensed-phase vapour pressure in Pa.

    Above the triple point this is the ordinary liquid-vapour saturation
    pressure from CoolProp.  Below it, use the NBS solid-vapour correlation
    instead of extrapolating a liquid correlation into the solid domain.

    Only nitrogen and oxygen are provided because those are the bulk air
    components relevant to the current LH2 near-field work.
    """
    from CoolProp.CoolProp import PropsSI

    canonical = {
        "n2": "Nitrogen", "nitrogen": "Nitrogen",
        "o2": "Oxygen", "oxygen": "Oxygen",
    }.get(species.strip().lower())
    if canonical is None:
        raise ValueError(f"no cryogenic-air phase correlation for {species!r}")
    if temperature <= 0.0:
        raise ValueError("temperature must be positive")

    critical = _CRITICAL_TEMPERATURE[canonical]
    if temperature >= critical:
        return float("inf")
    triple = _TRIPLE_TEMPERATURE[canonical]
    if temperature < triple:
        return _solid_vapour_pressure(canonical, temperature)
    return float(PropsSI("P", "T", temperature, "Q", 1, canonical))


@dataclass(frozen=True)
class AirPhaseEquilibrium:
    """Ideal-gas/condensed split of N2 and O2 at fixed T and P."""

    temperature: float
    pressure: float
    hydrogen_flow: float
    nitrogen_total_flow: float
    oxygen_total_flow: float
    nitrogen_gas_flow: float
    oxygen_gas_flow: float
    nitrogen_condensed_flow: float
    oxygen_condensed_flow: float
    hydrogen_partial_pressure: float
    nitrogen_partial_pressure: float
    oxygen_partial_pressure: float

    @property
    def condensed_flow(self) -> float:
        return self.nitrogen_condensed_flow + self.oxygen_condensed_flow

    @property
    def gas_flow(self) -> float:
        return (
            self.hydrogen_flow
            + self.nitrogen_gas_flow
            + self.oxygen_gas_flow
        )


@dataclass(frozen=True)
class HydrogenEvaporationEndpoint:
    """First state with no liquid H2, allowing solid N2 and O2."""

    hydrogen_species: str
    temperature: float
    air_ratio: float
    hydrogen_mass_fraction: float
    nitrogen_gas_ratio: float
    oxygen_gas_ratio: float
    nitrogen_solid_ratio: float
    oxygen_solid_ratio: float
    hydrogen_partial_pressure: float
    nitrogen_partial_pressure: float
    oxygen_partial_pressure: float
    specific_momentum: float
    particle_velocity_fraction: float
    specific_kinetic_energy: float
    relative_energy_residual: float
    incoming_specific_kinetic_energy: float = 0.0
    incoming_specific_energy: float = 0.0
    outgoing_specific_energy: float = 0.0


@dataclass(frozen=True)
class MultiphaseHydrogenSourcePlane:
    """Momentum- and volume-consistent plane at the H2 evaporation point."""

    endpoint: HydrogenEvaporationEndpoint
    hydrogen_flow: float
    nitrogen_flow: float
    oxygen_flow: float
    total_flow: float
    momentum_flux: float
    kinetic_energy_flow: float
    gas_velocity: float
    particle_velocity: float
    velocity: float
    density: float
    diameter: float
    formation_distance: float


def equilibrium_air_phase_split(
    *, temperature: float, pressure: float, hydrogen_flow: float,
    nitrogen_flow: float, oxygen_flow: float,
) -> AirPhaseEquilibrium:
    """Split entrained N2/O2 between gas and stable condensed phases.

    The calculation is an active-set solution of component partial-pressure
    equilibrium.  For a condensed component ``i``, ``p_i = p_sat,i``; for an
    uncondensed component its whole supplied mass remains gaseous and its
    partial pressure must not exceed saturation.  This avoids the pure-air
    shortcut in Li et al. equation (17) and permits N2 and O2 to change phase
    simultaneously.

    Flows are kg/s.  Hydrogen is assumed gaseous; this helper is intended for
    the post-flash near field, not the liquid-H2 nozzle interior.
    """
    if temperature <= 0.0 or pressure <= 0.0:
        raise ValueError("temperature and pressure must be positive")
    if min(hydrogen_flow, nitrogen_flow, oxygen_flow) < 0.0:
        raise ValueError("component flows cannot be negative")
    if hydrogen_flow <= 0.0:
        raise ValueError("a positive gaseous-hydrogen flow is required")

    species = ("Nitrogen", "Oxygen")
    totals = {
        "Nitrogen": nitrogen_flow / _MOLECULAR_WEIGHT["Nitrogen"],
        "Oxygen": oxygen_flow / _MOLECULAR_WEIGHT["Oxygen"],
    }
    n_h2 = hydrogen_flow / _MOLECULAR_WEIGHT["Hydrogen"]
    saturation = {
        name: air_saturation_pressure(name, temperature) for name in species
    }

    solution = None
    # Two species make exhaustive active-set enumeration clearer and safer
    # than a convergence loop, particularly close to a phase boundary.
    for mask in range(4):
        condensed = {
            name for index, name in enumerate(species) if mask & (1 << index)
        }
        q_sum = sum(saturation[name] / pressure for name in condensed)
        if q_sum >= 1.0:
            continue
        n_total = (
            n_h2 + sum(totals[name] for name in species if name not in condensed)
        ) / (1.0 - q_sum)
        gas_moles = {
            name: (
                saturation[name] / pressure * n_total
                if name in condensed else totals[name]
            )
            for name in species
        }
        tolerance = 1.0e-10 * max(n_total, 1.0)
        if any(
            gas_moles[name] > totals[name] + tolerance for name in condensed
        ):
            continue
        if any(
            gas_moles[name] / n_total * pressure
            > saturation[name] * (1.0 + 1.0e-10)
            for name in species if name not in condensed
        ):
            continue
        solution = n_total, gas_moles
        break
    if solution is None:
        raise ValueError("no consistent N2/O2 phase split at this state")

    n_total, gas_moles = solution
    mw_n2 = _MOLECULAR_WEIGHT["Nitrogen"]
    mw_o2 = _MOLECULAR_WEIGHT["Oxygen"]
    n2_gas = min(gas_moles["Nitrogen"] * mw_n2, nitrogen_flow)
    o2_gas = min(gas_moles["Oxygen"] * mw_o2, oxygen_flow)
    n2_condensed = max(nitrogen_flow - n2_gas, 0.0)
    o2_condensed = max(oxygen_flow - o2_gas, 0.0)
    if n2_condensed <= 1.0e-12 * max(nitrogen_flow, 1.0):
        n2_condensed = 0.0
        n2_gas = nitrogen_flow
    if o2_condensed <= 1.0e-12 * max(oxygen_flow, 1.0):
        o2_condensed = 0.0
        o2_gas = oxygen_flow
    return AirPhaseEquilibrium(
        temperature=temperature,
        pressure=pressure,
        hydrogen_flow=hydrogen_flow,
        nitrogen_total_flow=nitrogen_flow,
        oxygen_total_flow=oxygen_flow,
        nitrogen_gas_flow=n2_gas,
        oxygen_gas_flow=o2_gas,
        nitrogen_condensed_flow=n2_condensed,
        oxygen_condensed_flow=o2_condensed,
        hydrogen_partial_pressure=n_h2 / n_total * pressure,
        nitrogen_partial_pressure=gas_moles["Nitrogen"] / n_total * pressure,
        oxygen_partial_pressure=gas_moles["Oxygen"] / n_total * pressure,
    )


def multiphase_hydrogen_evaporation_endpoint(
    *, storage_temperature: float, ambient_temperature: float = 298.15,
    ambient_pressure: float = 101325.0,
    nitrogen_mole_fraction: float = 0.78084,
    oxygen_mole_fraction: float = 0.20946,
    specific_momentum: float = 0.0,
    particle_velocity_fraction: float = 1.0,
    storage_pressure: float | None = None,
    hydrogen_species: str = "Hydrogen",
    incoming_specific_kinetic_energy: float = 0.0,
) -> HydrogenEvaporationEndpoint:
    """Solve the post-flash LH2 endpoint with stable condensed air.

    One kilogram per second of liquid H2 is used as the basis.  It is
    saturated by default; ``storage_pressure`` instead supplies an absolute
    pressure for an independently measured compressed/subcooled state.  At
    the endpoint H2 has just fully evaporated, so its partial pressure equals
    its saturation pressure.  N2 and O2 are simultaneously limited by their
    stable solid-vapour pressures.  Their supplied mass then follows directly
    from a single species-consistent enthalpy balance.

    ``incoming_specific_kinetic_energy`` is J/kg of incoming H2, not per
    kilogram of the air-entrained mixture. It defaults to zero for a resting
    reservoir. A measured moving pipe plane must supply ``u_pipe**2/2``;
    pressure work is already included in enthalpy and is not added again.

    This replaces the impossible 20 K *gaseous-air* endpoint without imposing
    a 68 or 77 K handoff.  The result is a hydrogen gas carrying condensed
    air, which still requires the transported source zone before a
    single-phase dispersion model may use it.
    """
    from CoolProp.CoolProp import PropsSI

    if hydrogen_species not in _HYDROGEN_SPECIES:
        raise ValueError(
            "hydrogen species must be 'Hydrogen', 'ParaHydrogen' or "
            "'OrthoHydrogen'"
        )
    if min(storage_temperature, ambient_temperature, ambient_pressure) <= 0.0:
        raise ValueError("temperatures and pressure must be positive")
    if storage_temperature >= float(PropsSI("Tcrit", hydrogen_species)):
        raise ValueError("storage state is not saturated liquid hydrogen")
    if ambient_temperature <= storage_temperature:
        raise ValueError("ambient must be warmer than stored liquid hydrogen")
    if specific_momentum < 0.0:
        raise ValueError("specific momentum cannot be negative")
    if (
        not math.isfinite(incoming_specific_kinetic_energy)
        or incoming_specific_kinetic_energy < 0.0
    ):
        raise ValueError("incoming specific kinetic energy must be finite and nonnegative")
    if not 0.0 <= particle_velocity_fraction <= 1.0:
        raise ValueError("particle velocity fraction must lie from zero to one")

    mole_sum = nitrogen_mole_fraction + oxygen_mole_fraction
    x_n2, x_o2 = (
        nitrogen_mole_fraction / mole_sum,
        oxygen_mole_fraction / mole_sum,
    )
    mw_n2 = _MOLECULAR_WEIGHT["Nitrogen"]
    mw_o2 = _MOLECULAR_WEIGHT["Oxygen"]
    mw_air = x_n2 * mw_n2 + x_o2 * mw_o2
    w_n2, w_o2 = x_n2 * mw_n2 / mw_air, x_o2 * mw_o2 / mw_air

    def pressure_residual(temp: float) -> float:
        return (
            float(PropsSI("P", "T", temp, "Q", 1, hydrogen_species))
            + air_saturation_pressure("Nitrogen", temp)
            + air_saturation_pressure("Oxygen", temp)
            - ambient_pressure
        )

    temperature = float(brentq(
        pressure_residual,
        float(PropsSI("Ttriple", hydrogen_species)) + 0.05,
        float(PropsSI("Tcrit", hydrogen_species)) - 0.05,
        xtol=1.0e-10,
    ))
    p_h2 = float(PropsSI(
        "P", "T", temperature, "Q", 1, hydrogen_species
    ))
    p_n2 = air_saturation_pressure("Nitrogen", temperature)
    p_o2 = air_saturation_pressure("Oxygen", temperature)
    n_h2 = 1.0 / float(PropsSI("M", hydrogen_species))
    n_total = n_h2 * ambient_pressure / p_h2
    n2_gas = n_total * p_n2 / ambient_pressure * mw_n2
    o2_gas = n_total * p_o2 / ambient_pressure * mw_o2

    if storage_pressure is None:
        absolute_h_store = PropsSI(
            "H", "T", storage_temperature, "Q", 0, hydrogen_species
        )
    else:
        if storage_pressure <= 0.0:
            raise ValueError("storage pressure must be positive")
        absolute_h_store = PropsSI(
            "H", "T", storage_temperature, "P", storage_pressure,
            hydrogen_species,
        )
    h_store = float(
        absolute_h_store
        - PropsSI(
            "H", "T|gas", ambient_temperature, "P", ambient_pressure,
            hydrogen_species,
        )
    )
    incoming_energy = h_store + incoming_specific_kinetic_energy
    h_h2 = _relative_gas_enthalpy(
        hydrogen_species, temperature, ambient_temperature, p_h2,
        ambient_pressure,
    )
    h_n2_gas = _relative_gas_enthalpy(
        "Nitrogen", temperature, ambient_temperature, ambient_pressure
    )
    h_o2_gas = _relative_gas_enthalpy(
        "Oxygen", temperature, ambient_temperature, ambient_pressure
    )
    h_n2_solid = _relative_solid_enthalpy(
        "Nitrogen", temperature, ambient_temperature, ambient_pressure
    )
    h_o2_solid = _relative_solid_enthalpy(
        "Oxygen", temperature, ambient_temperature, ambient_pressure
    )
    fixed = (
        h_h2
        + n2_gas * (h_n2_gas - h_n2_solid)
        + o2_gas * (h_o2_gas - h_o2_solid)
    )
    condensed_air_enthalpy = w_n2 * h_n2_solid + w_o2 * h_o2_solid
    gas_air_ratio = n2_gas + o2_gas

    def kinetic_energy(ratio: float) -> float:
        gas_ratio = 1.0 + gas_air_ratio
        condensed_ratio = max(ratio - gas_air_ratio, 0.0)
        momentum_mass = (
            gas_ratio + particle_velocity_fraction * condensed_ratio
        )
        velocity = specific_momentum / momentum_mass
        return 0.5 * velocity**2 * (
            gas_ratio
            + particle_velocity_fraction**2 * condensed_ratio
        )

    # Per unit H2 flow, q=1 reduces to
    # 0.5*specific_momentum**2/(1 + air_ratio).  q=0 assigns all momentum to
    # the gas, matching Li et al.'s stationary-condensate limiting case.
    # Closing the kinetic term here prevents the downstream transport ledger
    # from creating energy after storage enthalpy has been spent on phases.
    if specific_momentum == 0.0:
        air_ratio = (incoming_energy - fixed) / condensed_air_enthalpy
    else:
        def total_energy_residual(ratio: float) -> float:
            return (
                fixed + ratio * condensed_air_enthalpy
                + kinetic_energy(ratio)
                - incoming_energy
            )

        minimum_ratio = max(
            n2_gas / max(w_n2, 1.0e-30),
            o2_gas / max(w_o2, 1.0e-30),
            0.0,
        )
        if total_energy_residual(minimum_ratio) < 0.0:
            raise ValueError("no positive condensed-air total-energy balance")
        high = max(
            1.0,
            (incoming_energy - fixed) / condensed_air_enthalpy * 2.0,
            minimum_ratio * 2.0,
        )
        while total_energy_residual(high) > 0.0:
            high *= 2.0
            if high > 1.0e6:
                raise ValueError("condensed-air energy root is not bounded")
        air_ratio = float(brentq(
            total_energy_residual, minimum_ratio, high, xtol=1.0e-12
        ))
    if air_ratio <= 0.0:
        raise ValueError("no positive condensed-air balance at H2 evaporation")
    n2_total, o2_total = air_ratio * w_n2, air_ratio * w_o2
    if n2_total < n2_gas or o2_total < o2_gas:
        raise ValueError("air supply is insufficient for the saturated gas split")
    n2_solid, o2_solid = n2_total - n2_gas, o2_total - o2_gas
    phase_outgoing = (
        h_h2 + n2_gas * h_n2_gas + o2_gas * h_o2_gas
        + n2_solid * h_n2_solid + o2_solid * h_o2_solid
    )
    specific_kinetic_energy = kinetic_energy(air_ratio)
    outgoing = phase_outgoing + specific_kinetic_energy
    residual = abs(outgoing - incoming_energy) / max(abs(incoming_energy), 1.0)
    return HydrogenEvaporationEndpoint(
        hydrogen_species=hydrogen_species,
        temperature=temperature,
        air_ratio=air_ratio,
        hydrogen_mass_fraction=1.0 / (1.0 + air_ratio),
        nitrogen_gas_ratio=n2_gas,
        oxygen_gas_ratio=o2_gas,
        nitrogen_solid_ratio=n2_solid,
        oxygen_solid_ratio=o2_solid,
        hydrogen_partial_pressure=p_h2,
        nitrogen_partial_pressure=p_n2,
        oxygen_partial_pressure=p_o2,
        specific_momentum=specific_momentum,
        particle_velocity_fraction=particle_velocity_fraction,
        specific_kinetic_energy=specific_kinetic_energy,
        relative_energy_residual=residual,
        incoming_specific_kinetic_energy=incoming_specific_kinetic_energy,
        incoming_specific_energy=incoming_energy,
        outgoing_specific_energy=outgoing,
    )


def multiphase_hydrogen_source_plane(
    *, hydrogen_flow: float, orifice_diameter: float,
    orifice_density: float, storage_temperature: float,
    ambient_temperature: float = 298.15,
    ambient_pressure: float = 101325.0,
    nitrogen_mole_fraction: float = 0.78084,
    oxygen_mole_fraction: float = 0.20946,
    specific_momentum: float | None = None,
    include_kinetic_energy: bool = False,
    beta_a: float = 0.281,
    particle_velocity_fraction: float = 1.0,
    storage_pressure: float | None = None,
    hydrogen_species: str = "Hydrogen",
    incoming_specific_kinetic_energy: float = 0.0,
) -> MultiphaseHydrogenSourcePlane:
    """Build the first multiphase source plane after liquid H2 evaporates.

    Ambient air entering the evaporation control volume has zero streamwise
    momentum.  The common mixture velocity therefore follows from
    ``m_H2 u_orifice = m_total u_plane``.  Volume is the saturated-H2 vapour
    volume plus the condensed N2/O2 particle volume; the equilibrium gas
    fractions of air at this 20 K endpoint are negligible but are included.

    This plane is still multiphase and must be marched with
    :func:`transported_condensed_air_source` before it can be supplied to the
    single-phase JetPlume equations.
    """
    from CoolProp.CoolProp import PropsSI

    if min(hydrogen_flow, orifice_diameter, orifice_density, beta_a) <= 0.0:
        raise ValueError(
            "flow, orifice diameter, density and beta_a must be positive"
        )
    if incoming_specific_kinetic_energy != 0.0 and not include_kinetic_energy:
        raise ValueError("incoming kinetic energy requires outgoing kinetic energy")

    orifice_area = math.pi * orifice_diameter**2 / 4.0
    orifice_velocity = hydrogen_flow / (orifice_density * orifice_area)
    if specific_momentum is None:
        specific_momentum = orifice_velocity
    if specific_momentum <= 0.0:
        raise ValueError("specific momentum must be positive")

    endpoint = multiphase_hydrogen_evaporation_endpoint(
        storage_temperature=storage_temperature,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        nitrogen_mole_fraction=nitrogen_mole_fraction,
        oxygen_mole_fraction=oxygen_mole_fraction,
        specific_momentum=specific_momentum if include_kinetic_energy else 0.0,
        particle_velocity_fraction=particle_velocity_fraction,
        storage_pressure=storage_pressure,
        hydrogen_species=hydrogen_species,
        incoming_specific_kinetic_energy=incoming_specific_kinetic_energy,
    )
    nitrogen_flow = hydrogen_flow * (
        endpoint.nitrogen_gas_ratio + endpoint.nitrogen_solid_ratio
    )
    oxygen_flow = hydrogen_flow * (
        endpoint.oxygen_gas_ratio + endpoint.oxygen_solid_ratio
    )
    total_flow = hydrogen_flow + nitrogen_flow + oxygen_flow
    momentum_flux = hydrogen_flow * specific_momentum
    gas_flow = hydrogen_flow * (
        1.0 + endpoint.nitrogen_gas_ratio + endpoint.oxygen_gas_ratio
    )
    condensed_flow = total_flow - gas_flow
    gas_velocity = momentum_flux / (
        gas_flow + particle_velocity_fraction * condensed_flow
    )
    particle_velocity = particle_velocity_fraction * gas_velocity
    kinetic_energy_flow = hydrogen_flow * endpoint.specific_kinetic_energy
    velocity = momentum_flux / total_flow

    # The endpoint lies effectively on the pure-H2 saturation boundary:
    # N2/O2 vapour pressures are below 2e-8 Pa.  Use CoolProp for the
    # non-ideal saturated H2 volume and ideal partial-gas volumes for the
    # vanishing air-vapour correction.
    rho_h2 = _hydrogen_gas_density(
        endpoint.temperature,
        endpoint.hydrogen_partial_pressure,
        hydrogen_species,
    )
    gas_volume = hydrogen_flow / rho_h2
    condensed_volume = hydrogen_flow * (
        endpoint.nitrogen_solid_ratio / _SOLID_DENSITY["Nitrogen"]
        + endpoint.oxygen_solid_ratio / _SOLID_DENSITY["Oxygen"]
    )
    volume_flow = gas_volume + condensed_volume
    density = total_flow / volume_flow
    diameter = math.sqrt(
        4.0 * total_flow / (math.pi * density * velocity)
    )
    rho_ambient = float(PropsSI(
        "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
    ))
    entrainment_scale = beta_a * math.sqrt(momentum_flux / rho_ambient)
    formation_distance = (
        (1.0 - endpoint.hydrogen_mass_fraction)
        * endpoint.air_ratio * hydrogen_flow
        / (entrainment_scale * rho_ambient)
    )
    return MultiphaseHydrogenSourcePlane(
        endpoint=endpoint,
        hydrogen_flow=hydrogen_flow,
        nitrogen_flow=nitrogen_flow,
        oxygen_flow=oxygen_flow,
        total_flow=total_flow,
        momentum_flux=momentum_flux,
        kinetic_energy_flow=kinetic_energy_flow,
        gas_velocity=gas_velocity,
        particle_velocity=particle_velocity,
        velocity=velocity,
        density=density,
        diameter=diameter,
        formation_distance=formation_distance,
    )


def particle_relaxation_time(
    *, diameter: float, particle_density: float, gas_viscosity: float,
) -> float:
    """Low-Re Stokes velocity-relaxation time, seconds."""
    if min(diameter, particle_density, gas_viscosity) <= 0.0:
        raise ValueError("diameter, density and viscosity must be positive")
    return particle_density * diameter**2 / (18.0 * gas_viscosity)


def ranz_marshall_transfer_number(
    *, particle_reynolds: float, transport_number: float,
) -> float:
    """Uncorrected sphere Nusselt/Sherwood number.

    This is the form used by Sandia's thin-skin spray model,
    ``2 * (1 + Re**0.5 * Pr**(1/3) / 3)``.  ``transport_number`` is Prandtl
    for heat transfer and Schmidt for mass transfer.  The function deliberately
    excludes a Stefan-flow correction: it is the fastest-transfer side of a
    diagnostic bound, not a complete finite-rate phase-change model.
    """
    if particle_reynolds < 0.0:
        raise ValueError("particle Reynolds number cannot be negative")
    if transport_number <= 0.0:
        raise ValueError("transport number must be positive")
    return 2.0 * (
        1.0
        + math.sqrt(particle_reynolds) * transport_number ** (1.0 / 3.0)
        / 3.0
    )


def minimum_heat_limited_sublimation_time(
    *, species: str, diameter: float, gas_temperature: float,
    particle_temperature: float, gas_thermal_conductivity: float,
    gas_prandtl: float, particle_reynolds: float = 0.0,
) -> float:
    """Shortest heat-limited lifetime of one condensed-air sphere, seconds.

    All convective heat is assigned to latent heat, the particle is held at
    its initial temperature, and Stefan-flow resistance is omitted.  These
    choices maximise phase-change rate.  From a spherical energy balance,

    ``d(diameter**2)/dt = -4 Nu k (T_g - T_p) / (rho_p L)``.

    The result is therefore a *lower bound* on survival time.  It is useful
    for rejecting instantaneous-equilibrium assumptions, but it cannot choose
    a particle diameter or replace a nucleation/growth model.
    """
    if species not in ("Nitrogen", "Oxygen"):
        raise ValueError("species must be Nitrogen or Oxygen")
    for name, value in {
        "diameter": diameter,
        "gas_temperature": gas_temperature,
        "particle_temperature": particle_temperature,
        "gas_thermal_conductivity": gas_thermal_conductivity,
        "gas_prandtl": gas_prandtl,
    }.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if particle_reynolds < 0.0:
        raise ValueError("particle Reynolds number cannot be negative")
    if gas_temperature <= particle_temperature:
        return math.inf

    density = _particle_density(species, particle_temperature)
    if particle_temperature < _TRIPLE_TEMPERATURE[species]:
        latent = _solid_sublimation_enthalpy(species, particle_temperature)
    else:
        from CoolProp.CoolProp import PropsSI

        latent = float(
            PropsSI("H", "T", particle_temperature, "Q", 1, species)
            - PropsSI("H", "T", particle_temperature, "Q", 0, species)
        )
    nusselt = ranz_marshall_transfer_number(
        particle_reynolds=particle_reynolds,
        transport_number=gas_prandtl,
    )
    return (
        density * latent * diameter**2
        / (
            4.0 * nusselt * gas_thermal_conductivity
            * (gas_temperature - particle_temperature)
        )
    )


def particle_terminal_velocity(
    *, diameter: float, particle_density: float, gas_density: float,
    gas_viscosity: float, gravity: float = 9.80665,
) -> float:
    """Settling speed from Schiller--Naumann drag, m/s.

    Stokes drag is recovered as particle Reynolds number tends to zero.  The
    nonlinear correction prevents a 100-micrometre stress case from being
    assigned the unrealistically large Stokes terminal speed.
    """
    if min(diameter, particle_density, gas_density, gas_viscosity, gravity) <= 0:
        raise ValueError("particle and gas properties must be positive")
    if particle_density <= gas_density:
        return 0.0

    def drag_coefficient(speed: float) -> float:
        reynolds = gas_density * speed * diameter / gas_viscosity
        if reynolds <= 1.0e-15:
            return float("inf")
        if reynolds < 1000.0:
            return 24.0 / reynolds * (1.0 + 0.15 * reynolds**0.687)
        return 0.44

    def residual(speed: float) -> float:
        if speed <= 0.0:
            return -(
                (particle_density - gas_density)
                * math.pi * diameter**3 / 6.0 * gravity
            )
        cd = drag_coefficient(speed)
        drag = (
            0.5 * cd * gas_density * math.pi * diameter**2 / 4.0 * speed**2
        )
        weight = (
            (particle_density - gas_density) * math.pi * diameter**3 / 6.0
            * gravity
        )
        return drag - weight

    high = max(
        particle_density * diameter**2 * gravity / (18.0 * gas_viscosity),
        1.0e-6,
    )
    while residual(high) < 0.0:
        high *= 2.0
    return float(brentq(residual, 0.0, high, xtol=1.0e-13))


@dataclass(frozen=True)
class CondensedAirStep:
    """One station in the finite-slip source-zone march."""

    distance: float
    temperature: float
    velocity: float
    diameter: float
    bulk_density: float
    hydrogen_mass_fraction: float
    nitrogen_condensed_flow: float
    oxygen_condensed_flow: float
    nitrogen_dropped_flow: float
    oxygen_dropped_flow: float
    mass_residual: float
    momentum_residual: float
    energy_residual: float


@dataclass(frozen=True)
class CondensedAirSource:
    """Result of marching the cryogenic multiphase source zone."""

    temperature: float
    velocity: float
    diameter: float
    density: float
    hydrogen_mass_fraction: float
    nitrogen_gas_flow: float
    oxygen_gas_flow: float
    nitrogen_condensed_flow: float
    oxygen_condensed_flow: float
    nitrogen_dropped_flow: float
    oxygen_dropped_flow: float
    entrained_air_flow: float
    distance: float
    handoff_reached: bool
    maximum_mass_residual: float
    maximum_momentum_residual: float
    maximum_energy_residual: float
    rows: tuple[CondensedAirStep, ...]

    @property
    def retained_condensed_flow(self) -> float:
        return self.nitrogen_condensed_flow + self.oxygen_condensed_flow


_R_UNIVERSAL = 8.31446261815324  # J/(mol K)
_IDEAL_CP = {
    "Nitrogen": 3.5 * _R_UNIVERSAL / 0.0280134,
    "Oxygen": 3.5 * _R_UNIVERSAL / 0.0319988,
}
_SOLID_DENSITY = {"Nitrogen": 1026.0, "Oxygen": 1300.0}
_SOLID_SUBLIMATION_ENTHALPY = {
    # Clausius-Clapeyron slopes fitted over the stable solid tables.  A
    # single slope keeps enthalpy continuous when the pressure lookup crosses
    # an interpolation knot; the tabulated pressure itself remains piecewise.
    "Nitrogen": _R_UNIVERSAL * 825.81417481 / 0.0280134,
    "Oxygen": _R_UNIVERSAL * 1047.42224397 / 0.0319988,
}


@lru_cache(maxsize=16384)
def _relative_gas_enthalpy(
    species: str, temperature: float, ambient_temperature: float,
    pressure: float, reference_pressure: float | None = None,
) -> float:
    """Species enthalpy relative to ambient gas, J/kg."""
    if reference_pressure is None:
        reference_pressure = pressure
    if species in _HYDROGEN_SPECIES:
        from CoolProp.CoolProp import PropsSI

        return float(
            PropsSI("H", "T|gas", temperature, "P", pressure, species)
            - PropsSI(
                "H", "T|gas", ambient_temperature, "P", reference_pressure,
                species,
            )
        )
    return _IDEAL_CP[species] * (temperature - ambient_temperature)


@lru_cache(maxsize=16384)
def _hydrogen_gas_density(
    temperature: float,
    partial_pressure: float,
    species: str = "Hydrogen",
) -> float:
    """H2 density on the gas branch at its component partial pressure."""
    from CoolProp.CoolProp import PropsSI

    if species not in _HYDROGEN_SPECIES:
        raise ValueError("invalid hydrogen spin species")
    return float(PropsSI(
        "D", "T|gas", temperature, "P", max(partial_pressure, 1.0e-9),
        species,
    ))


@lru_cache(maxsize=8192)
def _solid_sublimation_enthalpy(species: str, temperature: float) -> float:
    """Table-fit Clausius-Clapeyron sublimation enthalpy, J/kg."""
    del temperature  # constant-latent approximation over the solid range
    return _SOLID_SUBLIMATION_ENTHALPY[species]


@lru_cache(maxsize=16384)
def _relative_condensed_enthalpy(
    species: str, temperature: float, ambient_temperature: float,
    pressure: float,
) -> float:
    """Stable liquid/solid enthalpy relative to ambient gas, J/kg."""
    from CoolProp.CoolProp import PropsSI

    gas = _relative_gas_enthalpy(
        species, temperature, ambient_temperature, pressure
    )
    triple = _TRIPLE_TEMPERATURE[species]
    if temperature < triple:
        return gas - _solid_sublimation_enthalpy(species, temperature)
    return _relative_liquid_enthalpy(
        species, temperature, ambient_temperature, pressure
    )


def _relative_liquid_enthalpy(
    species: str, temperature: float, ambient_temperature: float,
    pressure: float,
) -> float:
    """Saturated-liquid enthalpy relative to ambient gas, J/kg."""
    from CoolProp.CoolProp import PropsSI

    gas = _relative_gas_enthalpy(
        species, temperature, ambient_temperature, pressure
    )
    latent = float(
        PropsSI("H", "T", temperature, "Q", 1, species)
        - PropsSI("H", "T", temperature, "Q", 0, species)
    )
    return gas - latent


def _relative_solid_enthalpy(
    species: str, temperature: float, ambient_temperature: float,
    pressure: float,
) -> float:
    """Solid enthalpy relative to ambient gas, J/kg."""
    return _relative_gas_enthalpy(
        species, temperature, ambient_temperature, pressure
    ) - _solid_sublimation_enthalpy(species, temperature)


@lru_cache(maxsize=8192)
def _particle_density(species: str, temperature: float) -> float:
    from CoolProp.CoolProp import PropsSI

    if temperature < _TRIPLE_TEMPERATURE[species]:
        return _SOLID_DENSITY[species]
    return float(PropsSI("D", "T", temperature, "Q", 0, species))


@lru_cache(maxsize=256)
def transported_condensed_air_source(
    *, hydrogen_flow: float, station2_temperature: float,
    station2_velocity: float, station2_density: float,
    station2_diameter: float, particle_diameter: float,
    initial_nitrogen_flow: float = 0.0, initial_oxygen_flow: float = 0.0,
    ambient_temperature: float = 298.15, ambient_pressure: float = 101325.0,
    nitrogen_mole_fraction: float = 0.78084,
    oxygen_mole_fraction: float = 0.20946,
    beta_a: float = 0.281, maximum_distance: float = 2.0,
    steps: int = 400, handoff_fraction: float = 0.01,
    maximum_air_increment: float = 0.02,
    station2_kinetic_energy: float | None = None,
    stationary_condensate: bool = False,
) -> CondensedAirSource:
    """March a phase-equilibrium, settling condensed-air source zone.

    This is an off-line source model, not yet a production JetPlume option.
    Ambient air is entrained into a uniform plug flow.  Both N2 and O2 may
    condense or freeze; normally retained particles share the axial plug
    velocity and settle out over their Schiller--Naumann terminal velocity.
    Dropped mass removes its own enthalpy, kinetic energy and axial momentum.
    Retained condensed material may re-evaporate at later, warmer stations.

    ``stationary_condensate`` is Li et al.'s opposite kinematic bound.  Every
    equilibrium condensed increment leaves with zero axial velocity, so it
    removes enthalpy but no axial momentum and particle diameter has no role.
    ``station2_kinetic_energy`` lets that two-velocity endpoint carry its
    separately closed kinetic-energy flux into this ledger.

    The common axial velocity is a controlled first implementation.  The
    particle relaxation-time diagnostic quantifies when this approximation
    fails; a two-velocity extension must precede production use of the 100 um
    stress case.
    """
    from CoolProp.CoolProp import PropsSI

    values = {
        "hydrogen_flow": hydrogen_flow,
        "station2_temperature": station2_temperature,
        "station2_velocity": station2_velocity,
        "station2_density": station2_density,
        "station2_diameter": station2_diameter,
        "particle_diameter": particle_diameter,
        "ambient_temperature": ambient_temperature,
        "ambient_pressure": ambient_pressure,
        "beta_a": beta_a,
        "maximum_distance": maximum_distance,
        "handoff_fraction": handoff_fraction,
        "maximum_air_increment": maximum_air_increment,
    }
    for name, value in values.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if steps < 2:
        raise ValueError("steps must be at least two")
    if min(initial_nitrogen_flow, initial_oxygen_flow) < 0.0:
        raise ValueError("initial component flows cannot be negative")
    if station2_temperature >= ambient_temperature:
        raise ValueError("Station-2 temperature must be below ambient")
    if not 0.0 < handoff_fraction < 1.0:
        raise ValueError("handoff_fraction must lie between zero and one")

    # Normalise the bulk-air fractions after omitting Ar and traces.  This
    # keeps mass exactly closed inside the two-component phase model.
    mole_sum = nitrogen_mole_fraction + oxygen_mole_fraction
    x_n2, x_o2 = (
        nitrogen_mole_fraction / mole_sum,
        oxygen_mole_fraction / mole_sum,
    )
    mw_n2, mw_o2 = 0.0280134, 0.0319988
    mw_air = x_n2 * mw_n2 + x_o2 * mw_o2
    w_n2, w_o2 = x_n2 * mw_n2 / mw_air, x_o2 * mw_o2 / mw_air

    area2 = math.pi * station2_diameter**2 / 4.0
    initial_retained_flow = (
        hydrogen_flow + initial_nitrogen_flow + initial_oxygen_flow
    )
    continuity_flow = station2_density * station2_velocity * area2
    if abs(continuity_flow / initial_retained_flow - 1.0) > 0.05:
        raise ValueError(
            "Station-2 mass flow, density, velocity and diameter disagree "
            "by more than 5%"
        )

    initial_split = equilibrium_air_phase_split(
        temperature=station2_temperature, pressure=ambient_pressure,
        hydrogen_flow=hydrogen_flow, nitrogen_flow=initial_nitrogen_flow,
        oxygen_flow=initial_oxygen_flow,
    )
    h_h2_initial = _relative_gas_enthalpy(
        "Hydrogen", station2_temperature, ambient_temperature,
        initial_split.hydrogen_partial_pressure, ambient_pressure,
    )
    initial_phase_enthalpy = (
        hydrogen_flow * h_h2_initial
        + initial_split.nitrogen_gas_flow * _relative_gas_enthalpy(
            "Nitrogen", station2_temperature, ambient_temperature,
            ambient_pressure,
        )
        + initial_split.oxygen_gas_flow * _relative_gas_enthalpy(
            "Oxygen", station2_temperature, ambient_temperature,
            ambient_pressure,
        )
        + initial_split.nitrogen_condensed_flow * _relative_condensed_enthalpy(
            "Nitrogen", station2_temperature, ambient_temperature,
            ambient_pressure,
        )
        + initial_split.oxygen_condensed_flow * _relative_condensed_enthalpy(
            "Oxygen", station2_temperature, ambient_temperature,
            ambient_pressure,
        )
    )
    momentum = initial_retained_flow * station2_velocity
    if station2_kinetic_energy is not None and station2_kinetic_energy < 0.0:
        raise ValueError("Station-2 kinetic-energy flow cannot be negative")
    if stationary_condensate and station2_kinetic_energy is None:
        raise ValueError(
            "stationary condensate requires its two-velocity kinetic energy"
        )
    initial_kinetic_energy = (
        0.5 * initial_retained_flow * station2_velocity**2
        if station2_kinetic_energy is None else station2_kinetic_energy
    )
    energy = initial_phase_enthalpy + initial_kinetic_energy
    initial_momentum = momentum
    n2_total, o2_total = initial_nitrogen_flow, initial_oxygen_flow
    n2_dropped = o2_dropped = 0.0
    dropped_momentum = 0.0
    entrained_air = initial_nitrogen_flow + initial_oxygen_flow
    temperature = station2_temperature
    peak_condensed = initial_split.condensed_flow
    max_mass_residual = max_momentum_residual = max_energy_residual = 0.0
    rho_ambient = float(PropsSI(
        "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
    ))
    nominal_dx = maximum_distance / steps
    distance = 0.0
    rows: list[CondensedAirStep] = []

    def phase_enthalpy(
        temp: float, split: AirPhaseEquilibrium,
        condensed_enthalpy: dict[str, float] | None = None,
    ) -> float:
        condensed_enthalpy = condensed_enthalpy or {}
        enthalpy = (
            hydrogen_flow * _relative_gas_enthalpy(
                "Hydrogen", temp, ambient_temperature,
                split.hydrogen_partial_pressure, ambient_pressure,
            )
            + split.nitrogen_gas_flow * _relative_gas_enthalpy(
                "Nitrogen", temp, ambient_temperature, ambient_pressure
            )
            + split.oxygen_gas_flow * _relative_gas_enthalpy(
                "Oxygen", temp, ambient_temperature, ambient_pressure
            )
        )
        if split.nitrogen_condensed_flow > 0.0:
            h_nitrogen = condensed_enthalpy.get("Nitrogen")
            if h_nitrogen is None:
                h_nitrogen = _relative_condensed_enthalpy(
                    "Nitrogen", temp, ambient_temperature, ambient_pressure
                )
            enthalpy += (
                split.nitrogen_condensed_flow * h_nitrogen
            )
        if split.oxygen_condensed_flow > 0.0:
            h_oxygen = condensed_enthalpy.get("Oxygen")
            if h_oxygen is None:
                h_oxygen = _relative_condensed_enthalpy(
                    "Oxygen", temp, ambient_temperature, ambient_pressure
                )
            enthalpy += (
                split.oxygen_condensed_flow * h_oxygen
            )
        return enthalpy

    def solve_temperature(lower: float):
        retained_flow = hydrogen_flow + n2_total + o2_total

        def residual(
            temp: float, condensed_enthalpy: dict[str, float] | None = None,
        ):
            split = equilibrium_air_phase_split(
                temperature=temp, pressure=ambient_pressure,
                hydrogen_flow=hydrogen_flow, nitrogen_flow=n2_total,
                oxygen_flow=o2_total,
            )
            kinetic_mass_flow = (
                split.gas_flow if stationary_condensate else retained_flow
            )
            kinetic = momentum**2 / (2.0 * kinetic_mass_flow)
            return phase_enthalpy(temp, split, condensed_enthalpy) + kinetic - energy

        low = max(15.0, min(lower, ambient_temperature - 1.0e-6))
        boundaries = [low]
        boundaries.extend(
            triple for triple in (
                _TRIPLE_TEMPERATURE["Oxygen"],
                _TRIPLE_TEMPERATURE["Nitrogen"],
            ) if low < triple < ambient_temperature
        )
        boundaries.append(ambient_temperature)
        epsilon = 1.0e-5

        def ordinary_state(temp: float):
            split = equilibrium_air_phase_split(
                temperature=temp, pressure=ambient_pressure,
                hydrogen_flow=hydrogen_flow, nitrogen_flow=n2_total,
                oxygen_flow=o2_total,
            )
            condensed_h = {
                species: _relative_condensed_enthalpy(
                    species, temp, ambient_temperature, ambient_pressure
                )
                for species, flow in (
                    ("Nitrogen", split.nitrogen_condensed_flow),
                    ("Oxygen", split.oxygen_condensed_flow),
                ) if flow > 0.0
            }
            particle_rho = {
                species: _particle_density(species, temp)
                for species, flow in (
                    ("Nitrogen", split.nitrogen_condensed_flow),
                    ("Oxygen", split.oxygen_condensed_flow),
                ) if flow > 0.0
            }
            return split, condensed_h, particle_rho

        for start, end in zip(boundaries[:-1], boundaries[1:]):
            a = start if start == low else start + epsilon
            b = end if end == ambient_temperature else end - epsilon
            if a <= b:
                fa, fb = residual(a), residual(b)
                if fa == 0.0 or fa * fb <= 0.0:
                    temp = a if fa == 0.0 else float(
                        brentq(residual, a, b, xtol=1.0e-8)
                    )
                    split, condensed_h, particle_rho = ordinary_state(temp)
                    scale = max(abs(energy), hydrogen_flow * 1.0e5, 1.0)
                    return (
                        temp, split, abs(residual(temp)) / scale,
                        condensed_h, particle_rho,
                    )

            # A pure species melts isothermally at its triple point.  If the
            # energy lies inside that latent-heat jump, solve the solid/liquid
            # fraction instead of asking a temperature root to cross a
            # discontinuity.
            if end == ambient_temperature:
                continue
            triple = end
            species = (
                "Oxygen" if math.isclose(
                    triple, _TRIPLE_TEMPERATURE["Oxygen"]
                ) else "Nitrogen"
            )
            split = equilibrium_air_phase_split(
                temperature=triple, pressure=ambient_pressure,
                hydrogen_flow=hydrogen_flow, nitrogen_flow=n2_total,
                oxygen_flow=o2_total,
            )
            condensed_flow = (
                split.oxygen_condensed_flow if species == "Oxygen"
                else split.nitrogen_condensed_flow
            )
            if condensed_flow <= 0.0:
                continue
            solid_h = _relative_solid_enthalpy(
                species, triple, ambient_temperature, ambient_pressure
            )
            liquid_h = _relative_liquid_enthalpy(
                species, triple, ambient_temperature, ambient_pressure
            )
            base_h = {
                other: _relative_condensed_enthalpy(
                    other, triple, ambient_temperature, ambient_pressure
                )
                for other, flow in (
                    ("Nitrogen", split.nitrogen_condensed_flow),
                    ("Oxygen", split.oxygen_condensed_flow),
                ) if other != species and flow > 0.0
            }
            f_solid = residual(triple, {**base_h, species: solid_h})
            f_liquid = residual(triple, {**base_h, species: liquid_h})
            if f_solid * f_liquid <= 0.0 and f_liquid != f_solid:
                liquid_fraction = -f_solid / (f_liquid - f_solid)
                effective_h = solid_h + liquid_fraction * (liquid_h - solid_h)
                condensed_h = {**base_h, species: effective_h}
                solid_rho = _SOLID_DENSITY[species]
                liquid_rho = float(PropsSI(
                    "D", "T", triple, "Q", 0, species
                ))
                effective_rho = 1.0 / (
                    (1.0 - liquid_fraction) / solid_rho
                    + liquid_fraction / liquid_rho
                )
                particle_rho = {
                    other: _particle_density(other, triple)
                    for other, flow in (
                        ("Nitrogen", split.nitrogen_condensed_flow),
                        ("Oxygen", split.oxygen_condensed_flow),
                    ) if other != species and flow > 0.0
                }
                particle_rho[species] = effective_rho
                return triple, split, 0.0, condensed_h, particle_rho

        if low > 15.0:
            return solve_temperature(15.0)
        raise ValueError("no source-zone temperature closes the energy balance")

    split = initial_split
    handoff_reached = False
    final_density = station2_density
    final_diameter = station2_diameter
    final_velocity = station2_velocity

    for _index in range(1, 100001):
        if distance >= maximum_distance:
            break
        entrainment_area_rate = beta_a * math.sqrt(
            max(momentum, 0.0) / rho_ambient
        )
        dx = min(
            nominal_dx,
            maximum_distance - distance,
            maximum_air_increment * hydrogen_flow
            / max(rho_ambient * entrainment_area_rate, 1.0e-30),
        )
        distance += dx
        added_air = rho_ambient * entrainment_area_rate * dx
        entrained_air += added_air
        n2_total += added_air * w_n2
        o2_total += added_air * w_o2

        (
            temperature, split, e_residual, condensed_enthalpy,
            particle_density,
        ) = solve_temperature(temperature)
        max_energy_residual = max(max_energy_residual, e_residual)
        retained_flow = hydrogen_flow + n2_total + o2_total
        velocity_mass_flow = (
            split.gas_flow if stationary_condensate else retained_flow
        )
        velocity = momentum / velocity_mass_flow

        # With Dalton partial pressures, m_H2/rho_H2(T,p_H2) is the common
        # mixture volume.  It retains H2 non-ideality (Z=0.91 at the 20 K
        # endpoint) and tends to the ideal total-mole volume as Z tends to 1.
        gas_volume_flow = hydrogen_flow / _hydrogen_gas_density(
            temperature, split.hydrogen_partial_pressure
        )
        n2_particle_density = (
            particle_density["Nitrogen"]
            if split.nitrogen_condensed_flow > 0.0 else 1.0
        )
        o2_particle_density = (
            particle_density["Oxygen"]
            if split.oxygen_condensed_flow > 0.0 else 1.0
        )
        condensed_volume_flow = 0.0 if stationary_condensate else (
            split.nitrogen_condensed_flow / n2_particle_density
            + split.oxygen_condensed_flow / o2_particle_density
        )
        volume_flow = gas_volume_flow + condensed_volume_flow
        area = volume_flow / max(velocity, 1.0e-12)
        diameter = 2.0 * math.sqrt(area / math.pi)
        radius = diameter / 2.0
        bulk_density = retained_flow / volume_flow
        gas_density = split.gas_flow / gas_volume_flow
        gas_viscosity = float(PropsSI(
            "V", "T|gas", temperature, "P",
            max(split.hydrogen_partial_pressure, 1.0e-9), "Hydrogen"
        ))
        residence = dx / max(velocity, 1.0e-12)

        drops = {}
        for species, condensed_flow, particle_density in (
            ("Nitrogen", split.nitrogen_condensed_flow, n2_particle_density),
            ("Oxygen", split.oxygen_condensed_flow, o2_particle_density),
        ):
            if stationary_condensate:
                fraction = 1.0
            else:
                terminal = particle_terminal_velocity(
                    diameter=particle_diameter,
                    particle_density=particle_density,
                    gas_density=gas_density,
                    gas_viscosity=gas_viscosity,
                )
                fraction = 1.0 - math.exp(
                    -terminal * residence / max(radius, 1.0e-12)
                )
            drops[species] = min(max(fraction, 0.0), 1.0) * condensed_flow

        n2_drop, o2_drop = drops["Nitrogen"], drops["Oxygen"]
        dropped = n2_drop + o2_drop
        if dropped > 0.0:
            drop_energy = (
                n2_drop * condensed_enthalpy.get("Nitrogen", 0.0)
                + o2_drop * condensed_enthalpy.get("Oxygen", 0.0)
                + (0.0 if stationary_condensate else 0.5 * dropped * velocity**2)
            )
            energy -= drop_energy
            momentum_loss = 0.0 if stationary_condensate else dropped * velocity
            momentum -= momentum_loss
            dropped_momentum += momentum_loss
            n2_total -= n2_drop
            o2_total -= o2_drop
            n2_dropped += n2_drop
            o2_dropped += o2_drop
            (
                temperature, split, e_residual, condensed_enthalpy,
                particle_density,
            ) = solve_temperature(temperature)
            max_energy_residual = max(max_energy_residual, e_residual)
            retained_flow = hydrogen_flow + n2_total + o2_total
            velocity_mass_flow = (
                split.gas_flow if stationary_condensate else retained_flow
            )
            velocity = momentum / velocity_mass_flow

            gas_volume_flow = hydrogen_flow / _hydrogen_gas_density(
                temperature, split.hydrogen_partial_pressure
            )
            condensed_volume_flow = 0.0
            if not stationary_condensate and split.nitrogen_condensed_flow > 0.0:
                condensed_volume_flow += (
                    split.nitrogen_condensed_flow
                    / particle_density["Nitrogen"]
                )
            if not stationary_condensate and split.oxygen_condensed_flow > 0.0:
                condensed_volume_flow += (
                    split.oxygen_condensed_flow
                    / particle_density["Oxygen"]
                )
            volume_flow = gas_volume_flow + condensed_volume_flow
            area = volume_flow / max(velocity, 1.0e-12)
            diameter = 2.0 * math.sqrt(area / math.pi)
            bulk_density = retained_flow / volume_flow

        condensed = split.condensed_flow
        peak_condensed = max(peak_condensed, condensed)
        mass_residual = abs(
            hydrogen_flow + entrained_air
            - (hydrogen_flow + n2_total + o2_total + n2_dropped + o2_dropped)
        ) / max(hydrogen_flow + entrained_air, 1.0e-30)
        max_mass_residual = max(max_mass_residual, mass_residual)
        momentum_residual = abs(
            initial_momentum - dropped_momentum - momentum
        ) / max(initial_momentum, 1.0e-30)
        max_momentum_residual = max(max_momentum_residual, momentum_residual)

        rows.append(CondensedAirStep(
            distance=distance, temperature=temperature, velocity=velocity,
            diameter=diameter, bulk_density=bulk_density,
            hydrogen_mass_fraction=hydrogen_flow / retained_flow,
            nitrogen_condensed_flow=split.nitrogen_condensed_flow,
            oxygen_condensed_flow=split.oxygen_condensed_flow,
            nitrogen_dropped_flow=n2_dropped,
            oxygen_dropped_flow=o2_dropped,
            mass_residual=mass_residual,
            momentum_residual=momentum_residual,
            energy_residual=e_residual,
        ))
        final_density, final_diameter, final_velocity = (
            bulk_density, diameter, velocity
        )
        if peak_condensed > 0.0 and condensed <= handoff_fraction * peak_condensed:
            handoff_reached = True
            break

    retained_flow = hydrogen_flow + n2_total + o2_total
    return CondensedAirSource(
        temperature=temperature, velocity=final_velocity,
        diameter=final_diameter, density=final_density,
        hydrogen_mass_fraction=hydrogen_flow / retained_flow,
        nitrogen_gas_flow=split.nitrogen_gas_flow,
        oxygen_gas_flow=split.oxygen_gas_flow,
        nitrogen_condensed_flow=split.nitrogen_condensed_flow,
        oxygen_condensed_flow=split.oxygen_condensed_flow,
        nitrogen_dropped_flow=n2_dropped,
        oxygen_dropped_flow=o2_dropped,
        entrained_air_flow=entrained_air,
        distance=rows[-1].distance,
        handoff_reached=handoff_reached,
        maximum_mass_residual=max_mass_residual,
        maximum_momentum_residual=max_momentum_residual,
        maximum_energy_residual=max_energy_residual,
        rows=tuple(rows),
    )


@dataclass(frozen=True)
class LiZone3:
    """Plug-flow source at Station 3 of Li et al. (2026)."""

    temperature: float
    hydrogen_mass_fraction: float
    nitrogen_mass_fraction: float
    oxygen_mass_fraction: float
    entrained_air_flow: float
    condensed_nitrogen_flow: float
    gas_flow: float
    velocity: float
    density: float
    diameter: float
    zone_length: float
    nitrogen_condensed_fraction: float
    nitrogen_saturation_pressure: float
    used_subtriple_liquid_extrapolation: bool

    @property
    def gas_mass_fraction_sum(self) -> float:
        return (
            self.hydrogen_mass_fraction
            + self.nitrogen_mass_fraction
            + self.oxygen_mass_fraction
        )


def _gas_property(props, output: str, fluid: str, temperature: float,
                  pressure: float) -> float:
    """Evaluate the gas branch used by the paper's component mixture."""
    return float(props(output, "T|gas", temperature, "P", pressure, fluid))


def li2026_zone3(
    *,
    hydrogen_flow: float,
    station2_temperature: float,
    station2_velocity: float,
    station2_density: float,
    station2_diameter: float,
    ambient_temperature: float = 298.15,
    ambient_pressure: float = 101325.0,
    station3_temperature: float = 77.35,
    nitrogen_mole_fraction: float = 0.78,
    oxygen_mole_fraction: float = 0.22,
    beta_a: float = 0.281,
    allow_subtriple_liquid_extrapolation: bool = False,
) -> LiZone3:
    """Reproduce the paper's nitrogen-dropout Zone III.

    The implementation reduces the paper's three nonlinear equations to one
    scalar root in entrained-air flow.  It retains the exact published
    assumptions: only N2 condenses, condensed N2 has zero axial velocity, all
    component gas properties are evaluated at total ambient pressure, and
    ``gamma + h_LN2 = h_N2`` in equations (15)--(16).

    By default, a sub-triple-point Station 2 is rejected.  Set
    ``allow_subtriple_liquid_extrapolation=True`` only to reproduce the
    published 51--55 K calculations; it asks CoolProp for the same metastable
    liquid saturation extrapolation used in equation (17), not stable solid
    nitrogen physics.
    """
    from CoolProp.CoolProp import PropsSI as props

    positive = {
        "hydrogen_flow": hydrogen_flow,
        "station2_temperature": station2_temperature,
        "station2_velocity": station2_velocity,
        "station2_density": station2_density,
        "station2_diameter": station2_diameter,
        "ambient_temperature": ambient_temperature,
        "ambient_pressure": ambient_pressure,
        "station3_temperature": station3_temperature,
        "beta_a": beta_a,
    }
    for name, value in positive.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if not math.isclose(
        nitrogen_mole_fraction + oxygen_mole_fraction, 1.0,
        rel_tol=0.0, abs_tol=1.0e-12,
    ):
        raise ValueError("nitrogen and oxygen mole fractions must sum to one")

    if station2_temperature >= station3_temperature:
        return LiZone3(
            temperature=station2_temperature,
            hydrogen_mass_fraction=1.0,
            nitrogen_mass_fraction=0.0,
            oxygen_mass_fraction=0.0,
            entrained_air_flow=0.0,
            condensed_nitrogen_flow=0.0,
            gas_flow=hydrogen_flow,
            velocity=station2_velocity,
            density=station2_density,
            diameter=station2_diameter,
            zone_length=0.0,
            nitrogen_condensed_fraction=0.0,
            nitrogen_saturation_pressure=float("nan"),
            used_subtriple_liquid_extrapolation=False,
        )

    n2_triple = float(props("Ttriple", "Nitrogen"))
    below_triple = station2_temperature < n2_triple
    if below_triple and not allow_subtriple_liquid_extrapolation:
        raise ValueError(
            "Li et al. equation (17) requests liquid-N2 saturation at "
            f"{station2_temperature:.2f} K, below the {n2_triple:.2f} K "
            "nitrogen triple point; use a solid-air model, or explicitly "
            "enable the metastable paper reproduction"
        )

    p_sat_n2 = float(
        props("P", "T", station2_temperature, "Q", 1, "Nitrogen")
    )
    sigma = min(max(
        1.0 - p_sat_n2 / (nitrogen_mole_fraction * ambient_pressure),
        0.0,
    ), 1.0)

    mw_n2 = float(props("M", "Nitrogen"))
    mw_o2 = float(props("M", "Oxygen"))
    mw_air = nitrogen_mole_fraction * mw_n2 + oxygen_mole_fraction * mw_o2
    w_n2 = nitrogen_mole_fraction * mw_n2 / mw_air
    w_o2 = oxygen_mole_fraction * mw_o2 / mw_air

    h2_in = float(props(
        "H", "T", station2_temperature, "P", ambient_pressure, "Hydrogen"
    ))
    h2_out = float(props(
        "H", "T", station3_temperature, "P", ambient_pressure, "Hydrogen"
    ))
    h_n2_out = _gas_property(
        props, "H", "Nitrogen", station3_temperature, ambient_pressure
    )
    h_o2_out = _gas_property(
        props, "H", "Oxygen", station3_temperature, ambient_pressure
    )
    # A CoolProp `Air` enthalpy has a different arbitrary reference from the
    # pure-component calls in Eq. (20).  Use the same species references at
    # both boundaries; this is also what reproduces Fig. 4.
    h_air = (
        w_n2 * float(props(
            "H", "T", ambient_temperature, "P", ambient_pressure, "Nitrogen"
        ))
        + w_o2 * float(props(
            "H", "T", ambient_temperature, "P", ambient_pressure, "Oxygen"
        ))
    )

    def state(air_flow: float) -> tuple[float, float, float]:
        condensed = air_flow * w_n2 * sigma
        gas_flow = hydrogen_flow + air_flow - condensed
        velocity = hydrogen_flow * station2_velocity / gas_flow
        return condensed, gas_flow, velocity

    def energy_residual(air_flow: float) -> float:
        condensed, gas_flow, velocity = state(air_flow)
        gas_n2 = air_flow * w_n2 * (1.0 - sigma)
        gas_o2 = air_flow * w_o2
        incoming = hydrogen_flow * (
            h2_in + 0.5 * station2_velocity**2
        ) + air_flow * h_air
        # Eqs. (15)--(16): h_LN2 + gamma is h_N2.  Writing it this way makes
        # the cancellation visible instead of pretending a latent term
        # independently heats the gas.
        outgoing = (
            hydrogen_flow * h2_out
            + (gas_n2 + condensed) * h_n2_out
            + gas_o2 * h_o2_out
            + 0.5 * gas_flow * velocity**2
        )
        return incoming - outgoing

    low, high = 0.0, max(hydrogen_flow, 1.0e-12)
    f_low, f_high = energy_residual(low), energy_residual(high)
    while f_low * f_high > 0.0 and high < 1.0e5 * hydrogen_flow:
        high *= 2.0
        f_high = energy_residual(high)
    if f_low * f_high > 0.0:
        raise ValueError("Li et al. Zone-III energy balance has no positive root")
    air_flow = float(brentq(energy_residual, low, high, xtol=1.0e-13))
    condensed, gas_flow, velocity = state(air_flow)

    y_h2 = hydrogen_flow / gas_flow
    y_n2 = air_flow * w_n2 * (1.0 - sigma) / gas_flow
    y_o2 = air_flow * w_o2 / gas_flow
    rho_h2 = _gas_property(
        props, "D", "Hydrogen", station3_temperature, ambient_pressure
    )
    rho_n2 = _gas_property(
        props, "D", "Nitrogen", station3_temperature, ambient_pressure
    )
    rho_o2 = _gas_property(
        props, "D", "Oxygen", station3_temperature, ambient_pressure
    )
    density = 1.0 / (y_h2 / rho_h2 + y_n2 / rho_n2 + y_o2 / rho_o2)
    diameter = 2.0 * math.sqrt(
        gas_flow / (density * velocity * math.pi)
    )
    rho_ambient = float(props(
        "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
    ))
    entrainment = beta_a * math.sqrt(
        math.pi * station2_diameter**2 / 4.0
        * station2_density * station2_velocity**2 / rho_ambient
    )
    zone_length = (1.0 - y_h2) * air_flow / (entrainment * rho_ambient)

    return LiZone3(
        temperature=station3_temperature,
        hydrogen_mass_fraction=y_h2,
        nitrogen_mass_fraction=y_n2,
        oxygen_mass_fraction=y_o2,
        entrained_air_flow=air_flow,
        condensed_nitrogen_flow=condensed,
        gas_flow=gas_flow,
        velocity=velocity,
        density=density,
        diameter=diameter,
        zone_length=zone_length,
        nitrogen_condensed_fraction=sigma,
        nitrogen_saturation_pressure=p_sat_n2,
        used_subtriple_liquid_extrapolation=below_triple,
    )
