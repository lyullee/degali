"""Atmospheric flash and initial LH2 droplet source reconstruction.

The ordinary LH2 source collapses the release to the point where all liquid
hydrogen has evaporated.  This module retains the state immediately after
pressure expansion: vapour quality, residual liquid flow and the initial
monodisperse droplet scale.  It is a source for a subsequent two-phase march,
not a single-phase dispersion boundary.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .axisymmetric_jet import AxisymmetricJetSource
    from .cryogenic_air import MultiphaseHydrogenSourcePlane


@dataclass(frozen=True)
class FlashingHydrogenDropletSource:
    """Mass-, momentum- and energy-consistent post-flash LH2 spray plane."""

    hydrogen_species: str
    mass_flow: float
    orifice_diameter: float
    upstream_temperature: float
    upstream_pressure: float
    upstream_quality: float | None
    upstream_density: float
    upstream_velocity: float
    ambient_pressure: float
    postflash_temperature: float
    postflash_quality: float
    postflash_density: float
    postflash_velocity: float
    postflash_diameter: float
    liquid_mass_flow: float
    vapour_mass_flow: float
    liquid_density: float
    liquid_kinematic_viscosity: float
    surface_tension: float
    jet_reynolds: float
    jet_weber: float
    shattered: bool
    droplet_size_coefficient: float
    droplet_diameter: float | None
    droplet_number_flux: float
    mass_residual: float
    momentum_residual: float
    energy_residual: float

    @property
    def liquid_mass_fraction(self) -> float:
        return 1.0 - self.postflash_quality


@dataclass(frozen=True)
class HomogeneousEquilibriumHydrogenSource:
    """First all-H2-vapour plane reached by collective air entrainment.

    This is the fast common-velocity bound.  It retains the measured
    post-flash liquid inventory in :attr:`postflash`, spends the entrained
    air enthalpy on that inventory, and starts a gas/condensed-air jet only at
    the point where liquid H2 reaches zero.
    """

    postflash: FlashingHydrogenDropletSource
    phase_plane: "MultiphaseHydrogenSourcePlane"
    source: "AxisymmetricJetSource"
    hydrogen_mass_residual: float
    total_mass_residual: float
    momentum_residual: float
    energy_residual: float

    @property
    def formation_distance(self) -> float:
        return self.phase_plane.formation_distance


def flashing_hydrogen_droplet_source(
    *,
    mass_flow: float,
    orifice_diameter: float,
    upstream_temperature: float,
    upstream_pressure: float,
    ambient_temperature: float = 298.15,
    ambient_pressure: float = 101325.0,
    upstream_quality: float | None = None,
    droplet_size_coefficient: float = 15.0,
    hydrogen_species: str = "Hydrogen",
) -> FlashingHydrogenDropletSource:
    """Expand a measured LH2 pipe state and retain its residual liquid.

    The atmospheric velocity follows the Yellow Book mass/momentum source
    balance, ``u_f = u_e + (P_e-P_a) A_e / m_dot``.  The flash quality then
    follows total specific-energy conservation.  The representative droplet
    diameter uses its Appleton/Wheatley correlation: the shattered-jet branch
    has the recommended ``C_ds=15`` and the documented 10--20 uncertainty
    range.  This coefficient is a published atomisation input and is never
    inferred from a dispersion observation.
    """
    from CoolProp.CoolProp import PropsSI

    if hydrogen_species not in {
        "Hydrogen", "ParaHydrogen", "OrthoHydrogen",
    }:
        raise ValueError(
            "hydrogen species must be 'Hydrogen', 'ParaHydrogen' or "
            "'OrthoHydrogen'"
        )
    for name, value in {
        "mass_flow": mass_flow,
        "orifice_diameter": orifice_diameter,
        "upstream_temperature": upstream_temperature,
        "upstream_pressure": upstream_pressure,
        "ambient_temperature": ambient_temperature,
        "ambient_pressure": ambient_pressure,
        "droplet_size_coefficient": droplet_size_coefficient,
    }.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if upstream_pressure <= ambient_pressure:
        raise ValueError("upstream pressure must exceed ambient pressure")
    if not 10.0 <= droplet_size_coefficient <= 20.0:
        raise ValueError(
            "droplet size coefficient must lie in the published 10--20 range"
        )

    if upstream_quality is None:
        upstream_quality_value = None
        upstream_density = float(PropsSI(
            "D", "T", upstream_temperature, "P", upstream_pressure,
            hydrogen_species,
        ))
        upstream_enthalpy = float(PropsSI(
            "H", "T", upstream_temperature, "P", upstream_pressure,
            hydrogen_species,
        ))
    else:
        upstream_quality_value = float(upstream_quality)
        if not 0.0 <= upstream_quality_value <= 1.0:
            raise ValueError("upstream quality must lie between zero and one")
        saturation_pressure = float(PropsSI(
            "P", "T", upstream_temperature, "Q", 0, hydrogen_species
        ))
        if not math.isclose(
            upstream_pressure, saturation_pressure, rel_tol=0.05
        ):
            raise ValueError(
                "a supplied upstream quality requires a saturation-consistent "
                "temperature and pressure"
            )
        upstream_density = float(PropsSI(
            "D", "T", upstream_temperature, "Q", upstream_quality_value,
            hydrogen_species,
        ))
        upstream_enthalpy = float(PropsSI(
            "H", "T", upstream_temperature, "Q", upstream_quality_value,
            hydrogen_species,
        ))

    area = math.pi * orifice_diameter**2 / 4.0
    upstream_velocity = mass_flow / (upstream_density * area)
    momentum_flux = (
        mass_flow * upstream_velocity
        + (upstream_pressure - ambient_pressure) * area
    )
    postflash_velocity = momentum_flux / mass_flow
    total_specific_energy = upstream_enthalpy + 0.5 * upstream_velocity**2
    postflash_enthalpy = total_specific_energy - 0.5 * postflash_velocity**2

    liquid_enthalpy = float(PropsSI(
        "H", "P", ambient_pressure, "Q", 0, hydrogen_species
    ))
    vapour_enthalpy = float(PropsSI(
        "H", "P", ambient_pressure, "Q", 1, hydrogen_species
    ))
    raw_quality = (
        (postflash_enthalpy - liquid_enthalpy)
        / (vapour_enthalpy - liquid_enthalpy)
    )
    postflash_quality = min(max(raw_quality, 0.0), 1.0)
    if 0.0 < raw_quality < 1.0:
        postflash_temperature = float(PropsSI(
            "T", "P", ambient_pressure, "Q", 0, hydrogen_species
        ))
        liquid_density = float(PropsSI(
            "D", "P", ambient_pressure, "Q", 0, hydrogen_species
        ))
        vapour_density = float(PropsSI(
            "D", "P", ambient_pressure, "Q", 1, hydrogen_species
        ))
        postflash_density = 1.0 / (
            (1.0 - postflash_quality) / liquid_density
            + postflash_quality / vapour_density
        )
    else:
        postflash_temperature = float(PropsSI(
            "T", "P", ambient_pressure, "H", postflash_enthalpy,
            hydrogen_species,
        ))
        postflash_density = float(PropsSI(
            "D", "P", ambient_pressure, "H", postflash_enthalpy,
            hydrogen_species,
        ))
        liquid_density = float(PropsSI(
            "D", "P", ambient_pressure, "Q", 0, hydrogen_species
        ))

    postflash_area = mass_flow / (postflash_density * postflash_velocity)
    postflash_diameter = math.sqrt(4.0 * postflash_area / math.pi)
    liquid_mass_flow = mass_flow * (1.0 - postflash_quality)
    vapour_mass_flow = mass_flow * postflash_quality

    liquid_viscosity = float(PropsSI(
        "V", "P", ambient_pressure, "Q", 0, hydrogen_species
    ))
    liquid_kinematic_viscosity = liquid_viscosity / liquid_density
    surface_tension = float(PropsSI(
        "I", "P", ambient_pressure, "Q", 0, hydrogen_species
    ))
    radius = 0.5 * postflash_diameter
    jet_reynolds = (
        2.0 * radius * postflash_velocity / liquid_kinematic_viscosity
    )
    jet_weber = (
        2.0 * radius * postflash_velocity**2 * liquid_density
        / surface_tension
    )
    normal_boiling_temperature = float(PropsSI(
        "T", "P", ambient_pressure, "Q", 0, hydrogen_species
    ))
    non_shattered = (
        jet_weber < jet_reynolds ** (-0.45) * 1.0e6
        and upstream_temperature < 1.11 * normal_boiling_temperature
    )
    shattered = not non_shattered

    droplet_diameter = None
    droplet_number_flux = 0.0
    if liquid_mass_flow > 0.0:
        if non_shattered:
            droplet_diameter = 3.78 * radius * (
                1.0 + 3.0 * math.sqrt(jet_weber) / jet_reynolds
            )
        else:
            ambient_density = float(PropsSI(
                "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
            ))
            droplet_diameter = (
                droplet_size_coefficient * surface_tension
                / (postflash_velocity**2 * ambient_density)
            )
        droplet_mass = (
            liquid_density * math.pi * droplet_diameter**3 / 6.0
        )
        droplet_number_flux = liquid_mass_flow / droplet_mass

    mass_residual = abs(
        postflash_density * postflash_velocity * postflash_area - mass_flow
    ) / mass_flow
    momentum_residual = abs(
        mass_flow * postflash_velocity - momentum_flux
    ) / max(abs(momentum_flux), 1.0)
    energy_residual = abs(
        postflash_enthalpy + 0.5 * postflash_velocity**2
        - total_specific_energy
    ) / max(abs(total_specific_energy), 1.0)
    return FlashingHydrogenDropletSource(
        hydrogen_species=hydrogen_species,
        mass_flow=mass_flow,
        orifice_diameter=orifice_diameter,
        upstream_temperature=upstream_temperature,
        upstream_pressure=upstream_pressure,
        upstream_quality=upstream_quality_value,
        upstream_density=upstream_density,
        upstream_velocity=upstream_velocity,
        ambient_pressure=ambient_pressure,
        postflash_temperature=postflash_temperature,
        postflash_quality=postflash_quality,
        postflash_density=postflash_density,
        postflash_velocity=postflash_velocity,
        postflash_diameter=postflash_diameter,
        liquid_mass_flow=liquid_mass_flow,
        vapour_mass_flow=vapour_mass_flow,
        liquid_density=liquid_density,
        liquid_kinematic_viscosity=liquid_kinematic_viscosity,
        surface_tension=surface_tension,
        jet_reynolds=jet_reynolds,
        jet_weber=jet_weber,
        shattered=shattered,
        droplet_size_coefficient=droplet_size_coefficient,
        droplet_diameter=droplet_diameter,
        droplet_number_flux=droplet_number_flux,
        mass_residual=mass_residual,
        momentum_residual=momentum_residual,
        energy_residual=energy_residual,
    )


def homogeneous_equilibrium_hydrogen_source(
    *,
    mass_flow: float,
    orifice_diameter: float,
    upstream_temperature: float,
    upstream_pressure: float,
    ambient_temperature: float = 298.15,
    ambient_pressure: float = 101325.0,
    upstream_quality: float | None = None,
    droplet_size_coefficient: float = 15.0,
    hydrogen_species: str = "Hydrogen",
    theta: float = 0.0,
    x: float = 0.0,
    y: float = 0.0,
) -> HomogeneousEquilibriumHydrogenSource:
    """Build a phase-safe source after the measured LH2 flash.

    The post-flash plane is reconstructed first, including its residual
    liquid.  The existing species-resolved N2/O2 energy balance then finds
    the amount of ambient air required to finish H2 evaporation while total
    momentum and kinetic energy are retained.  The reported formation length
    is the published constant-entrainment integral estimate; no receptor data
    or fitted evaporation distance enters this construction.

    Individual droplet heat-transfer resistance is assumed negligible here.
    Consequently this is a fastest homogeneous-equilibrium bound, not the
    later two-temperature finite-rate model.
    """
    from CoolProp.CoolProp import PropsSI

    from .axisymmetric_jet import AxisymmetricJetSource, SourceEnthalpyBoundary
    from .cryogenic_air import multiphase_hydrogen_source_plane

    postflash = flashing_hydrogen_droplet_source(
        mass_flow=mass_flow,
        orifice_diameter=orifice_diameter,
        upstream_temperature=upstream_temperature,
        upstream_pressure=upstream_pressure,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        upstream_quality=upstream_quality,
        droplet_size_coefficient=droplet_size_coefficient,
        hydrogen_species=hydrogen_species,
    )
    plane = multiphase_hydrogen_source_plane(
        hydrogen_flow=mass_flow,
        orifice_diameter=orifice_diameter,
        orifice_density=postflash.upstream_density,
        storage_temperature=upstream_temperature,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        specific_momentum=postflash.postflash_velocity,
        include_kinetic_energy=True,
        incoming_specific_kinetic_energy=0.5 * postflash.upstream_velocity**2,
        storage_pressure=upstream_pressure,
        hydrogen_species=hydrogen_species,
    )
    total_flow = plane.total_flow
    endpoint = plane.endpoint
    boundary = SourceEnthalpyBoundary(
        specific_enthalpy=(
            endpoint.outgoing_specific_energy - endpoint.specific_kinetic_energy
        ) * mass_flow / total_flow,
        ambient_temperature=ambient_temperature,
        hydrogen_species=hydrogen_species,
        hydrogen_reference_enthalpy=float(PropsSI(
            "H", "T|gas", ambient_temperature, "P", ambient_pressure,
            hydrogen_species,
        )),
    )
    source = AxisymmetricJetSource(
        diameter=plane.diameter,
        velocity=plane.velocity,
        density=plane.density,
        temperature=plane.endpoint.temperature,
        mass_fraction=mass_flow / total_flow,
        theta=theta,
        x=x + plane.formation_distance * math.cos(theta),
        y=y + plane.formation_distance * math.sin(theta),
        enthalpy_boundary=boundary,
    )
    source_momentum = source.mass_flow * source.velocity
    return HomogeneousEquilibriumHydrogenSource(
        postflash=postflash,
        phase_plane=plane,
        source=source,
        hydrogen_mass_residual=abs(source.fuel_mass_flow - mass_flow)
        / mass_flow,
        total_mass_residual=abs(source.mass_flow - total_flow) / total_flow,
        momentum_residual=abs(source_momentum - plane.momentum_flux)
        / max(abs(plane.momentum_flux), 1.0),
        energy_residual=plane.endpoint.relative_energy_residual,
    )


def minimum_heat_limited_hydrogen_evaporation_time(
    *,
    diameter: float,
    gas_temperature: float,
    droplet_temperature: float,
    gas_thermal_conductivity: float,
    gas_prandtl: float,
    particle_reynolds: float = 0.0,
    pressure: float = 101325.0,
    hydrogen_species: str = "Hydrogen",
) -> float:
    """Fastest convective lifetime of one LH2 sphere, seconds.

    The droplet remains at its initial temperature, all convective heat goes
    into sensible heating plus vaporisation, and Stefan-flow resistance is
    omitted.  These choices maximise evaporation.  The returned value is a
    lower bound on survival, useful only as a diagnostic: a droplet inside a
    cold jet sees the local gas temperature and shares the finite heat
    capacity of that jet rather than seeing the ambient directly.
    """
    from CoolProp.CoolProp import PropsSI

    if hydrogen_species not in {
        "Hydrogen", "ParaHydrogen", "OrthoHydrogen",
    }:
        raise ValueError("invalid hydrogen spin species")
    for name, value in {
        "diameter": diameter,
        "gas_temperature": gas_temperature,
        "droplet_temperature": droplet_temperature,
        "gas_thermal_conductivity": gas_thermal_conductivity,
        "gas_prandtl": gas_prandtl,
        "pressure": pressure,
    }.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if particle_reynolds < 0.0:
        raise ValueError("particle Reynolds number cannot be negative")
    if gas_temperature <= droplet_temperature:
        return math.inf
    saturation_temperature = float(PropsSI(
        "T", "P", pressure, "Q", 0, hydrogen_species
    ))
    if droplet_temperature > saturation_temperature + 1.0e-6:
        raise ValueError("droplet is not liquid at the supplied pressure")

    liquid_density = float(PropsSI(
        "D", "T|liquid", droplet_temperature, "P", pressure,
        hydrogen_species,
    ))
    liquid_enthalpy = float(PropsSI(
        "H", "T|liquid", droplet_temperature, "P", pressure,
        hydrogen_species,
    ))
    vapour_enthalpy = float(PropsSI(
        "H", "P", pressure, "Q", 1, hydrogen_species
    ))
    phase_change_enthalpy = vapour_enthalpy - liquid_enthalpy
    if phase_change_enthalpy <= 0.0:
        raise ValueError("non-positive hydrogen evaporation enthalpy")

    # Keep this module independent of the condensed-air implementation while
    # using the same uncorrected Ranz--Marshall diagnostic correlation.
    nusselt = 2.0 * (
        1.0
        + math.sqrt(particle_reynolds) * gas_prandtl ** (1.0 / 3.0) / 3.0
    )
    return (
        liquid_density * phase_change_enthalpy * diameter**2
        / (
            4.0 * nusselt * gas_thermal_conductivity
            * (gas_temperature - droplet_temperature)
        )
    )


def gasflow_phase_relaxation_coefficient(
    *,
    droplet_diameter: float,
    diffusion_coefficient: float,
    schmidt_number: float,
    particle_reynolds: float = 0.0,
) -> float:
    """GASFLOW-MPI droplet phase-relaxation coefficient, 1/s.

    PRESLHY D3.1 equation 45 writes ``c = 6 Sh D / d**2`` with
    ``Sh = 2 + 0.6 Re**0.5 Sc**(1/3)``.  ``D`` is the effective vapour-to-
    droplet mass diffusivity; callers must state whether it contains only
    molecular diffusion or also a turbulent contribution.  No empirical Lee
    relaxation constant is introduced here.
    """
    for name, value in {
        "droplet_diameter": droplet_diameter,
        "diffusion_coefficient": diffusion_coefficient,
        "schmidt_number": schmidt_number,
    }.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if particle_reynolds < 0.0:
        raise ValueError("particle Reynolds number cannot be negative")
    sherwood = 2.0 + 0.6 * math.sqrt(particle_reynolds) * (
        schmidt_number ** (1.0 / 3.0)
    )
    return 6.0 * sherwood * diffusion_coefficient / droplet_diameter**2


def critical_diffusivity_for_phase_delay(
    *,
    droplet_diameter: float,
    delay_time: float,
    schmidt_number: float,
    particle_reynolds: float = 0.0,
) -> float:
    """Effective diffusivity that would make ``1/c == delay_time``.

    This inversion is useful as an applicability screen: if the required
    diffusivity is orders of magnitude below any defensible gas diffusivity,
    finite-rate mass transfer cannot sustain droplets over that residence
    time within the GASFLOW closure.
    """
    for name, value in {
        "droplet_diameter": droplet_diameter,
        "delay_time": delay_time,
        "schmidt_number": schmidt_number,
    }.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if particle_reynolds < 0.0:
        raise ValueError("particle Reynolds number cannot be negative")
    sherwood = 2.0 + 0.6 * math.sqrt(particle_reynolds) * (
        schmidt_number ** (1.0 / 3.0)
    )
    return droplet_diameter**2 / (6.0 * sherwood * delay_time)


def critical_droplet_diameter_for_phase_delay(
    *,
    diffusion_coefficient: float,
    delay_time: float,
    schmidt_number: float,
    particle_reynolds: float = 0.0,
) -> float:
    """Droplet diameter that would make ``1/c == delay_time``, in metres.

    This is the diameter inversion of the PRESLHY D3.1/GASFLOW-MPI closure
    used by :func:`gasflow_phase_relaxation_coefficient`.  It provides a
    direct scale comparison with a measured or correlated droplet diameter;
    it is not a fitted residence-time parameter.
    """
    for name, value in {
        "diffusion_coefficient": diffusion_coefficient,
        "delay_time": delay_time,
        "schmidt_number": schmidt_number,
    }.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if particle_reynolds < 0.0:
        raise ValueError("particle Reynolds number cannot be negative")
    sherwood = 2.0 + 0.6 * math.sqrt(particle_reynolds) * (
        schmidt_number ** (1.0 / 3.0)
    )
    return math.sqrt(
        6.0 * sherwood * diffusion_coefficient * delay_time
    )
