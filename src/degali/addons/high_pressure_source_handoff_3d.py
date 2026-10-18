"""Conservative high-pressure H2 source handoff for the resolved 3-D tier.

The obstacle solver is a low-Mach model and must not receive a choked or
supersonic notional-nozzle velocity directly.  This module mixes an
ambient-pressure H2 source with explicitly accounted ambient air until a
declared physical Mach limit is reached.  Hydrogen mass, plane momentum and
total energy are closed before the state is regularised as a finite-width
volumetric source.

The entrainment distance uses the same source-momentum coefficient as the
existing conserved axisymmetric jet.  It is a modelled location, not a new
measurement or a concentration-fitted parameter.  The receiving resolved
solver injects H2 mass and the *source increment* of momentum/enthalpy; the
entrained ambient mass and its background momentum already belong to the
Eulerian carrier field.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING

from ..core.constants import R_UNIVERSAL
from .semi_fv_obstacle import SourceRateSchedule

if TYPE_CHECKING:
    from .notional import EnergyConservingNotionalNozzle


_H2_MOLAR_MASS_KG_MOL = 0.00201588
_AIR_MOLAR_MASS_KG_MOL = 0.02896546


def _positive_finite(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _unit_vector(values: tuple[float, float, float]) -> tuple[float, float, float]:
    if not isinstance(values, tuple) or len(values) != 3:
        raise TypeError("axis_unit_vector must be a three-value tuple")
    vector = tuple(float(value) for value in values)
    if any(not math.isfinite(value) for value in vector):
        raise ValueError("axis_unit_vector must be finite")
    magnitude = math.sqrt(math.fsum(value * value for value in vector))
    if magnitude <= 0.0:
        raise ValueError("axis_unit_vector must be non-zero")
    return tuple(value / magnitude for value in vector)


@dataclass(frozen=True)
class AmbientPressureHydrogenJet:
    """Complete pure-H2 source plane at ambient static pressure."""

    source_id: str
    source_model: str
    mass_flow_kg_s: float
    diameter_m: float
    velocity_m_s: float
    density_kg_m3: float
    temperature_k: float
    momentum_flux_n: float
    total_specific_energy_j_kg: float
    relative_mass_residual: float

    def __post_init__(self) -> None:
        if not isinstance(self.source_id, str) or not self.source_id.strip():
            raise ValueError("source_id must be non-empty")
        if not isinstance(self.source_model, str) or not self.source_model.strip():
            raise ValueError("source_model must be non-empty")
        for name in (
            "mass_flow_kg_s", "diameter_m", "velocity_m_s",
            "density_kg_m3", "temperature_k", "momentum_flux_n",
        ):
            object.__setattr__(self, name, _positive_finite(getattr(self, name), name))
        energy = float(self.total_specific_energy_j_kg)
        if not math.isfinite(energy):
            raise ValueError("total_specific_energy_j_kg must be finite")
        object.__setattr__(self, "total_specific_energy_j_kg", energy)
        residual = float(self.relative_mass_residual)
        if not math.isfinite(residual) or residual < 0.0:
            raise ValueError("relative_mass_residual must be finite and non-negative")
        object.__setattr__(self, "relative_mass_residual", residual)


@dataclass(frozen=True)
class ResolvedSourceConfigFields3D:
    """Fields that connect one accepted handoff to the resolved solver."""

    source_schedule: SourceRateSchedule
    source_position_m: tuple[float, float, float]
    source_sigma_m: float
    source_sigma_xyz_m: tuple[float, float, float]
    source_momentum_velocity_m_s: tuple[float, float, float]
    source_temperature_k: float
    cells_across_equivalent_diameter: float
    spatially_resolved: bool

    def as_config_kwargs(self) -> dict[str, object]:
        """Return only :class:`ResolvedObstacleFlow3DConfig` constructor keys."""
        return {
            "source_schedule": self.source_schedule,
            "source_position_m": self.source_position_m,
            "source_sigma_m": self.source_sigma_m,
            "source_sigma_xyz_m": self.source_sigma_xyz_m,
            "source_momentum_velocity_m_s": self.source_momentum_velocity_m_s,
            "source_temperature_k": self.source_temperature_k,
            "thermal_buoyancy_coupling": True,
            "active_momentum_solver": "finite_volume_projection",
        }


@dataclass(frozen=True)
class LowMachHydrogenHandoff3D:
    """Auditable mixed plane and its conservative resolved-source increment."""

    source: AmbientPressureHydrogenJet
    axis_unit_vector: tuple[float, float, float]
    ambient_temperature_k: float
    ambient_pressure_pa: float
    ambient_density_kg_m3: float
    ambient_axial_velocity_m_s: float
    maximum_mach_number: float
    maximum_relative_density_departure: float
    entrainment_beta: float
    handoff_distance_m: float
    entrained_air_mass_flow_kg_s: float
    total_mass_flow_kg_s: float
    hydrogen_mass_fraction: float
    mixture_velocity_m_s: float
    mixture_temperature_k: float
    mixture_density_kg_m3: float
    relative_density_departure: float
    mixture_molar_mass_kg_mol: float
    mixture_heat_capacity_j_kg_k: float
    mixture_sound_speed_m_s: float
    mach_number: float
    equivalent_diameter_m: float
    moment_matched_radial_sigma_m: float
    plane_momentum_flux_n: float
    source_increment_momentum_flux_n: float
    plane_sensible_enthalpy_flux_w: float
    source_increment_total_energy_flux_w: float
    equivalent_hydrogen_temperature_k: float
    relative_hydrogen_mass_residual: float
    relative_total_mass_residual: float
    relative_momentum_residual: float
    relative_total_energy_residual: float
    source_mass_consistency_residual: float
    external_validation_required: bool = True

    @property
    def maximum_relative_conservation_residual(self) -> float:
        return max(
            self.relative_hydrogen_mass_residual,
            self.relative_total_mass_residual,
            self.relative_momentum_residual,
            self.relative_total_energy_residual,
        )

    @property
    def numerically_admissible(self) -> bool:
        return (
            self.source_mass_consistency_residual <= 0.05
            and self.maximum_relative_conservation_residual <= 1.0e-10
            and self.mach_number <= self.maximum_mach_number * (1.0 + 1.0e-10)
            and abs(self.relative_density_departure)
            <= self.maximum_relative_density_departure * (1.0 + 1.0e-10)
        )

    @property
    def source_specific_momentum_m_s(self) -> float:
        """Momentum increment per kg of injected H2 used by the FV source."""
        return self.source_increment_momentum_flux_n / self.source.mass_flow_kg_s

    def resolved_source_fields(
        self,
        *,
        nozzle_position_m: tuple[float, float, float],
        grid_spacing_m: float,
        duration_s: float,
        minimum_cells_across_diameter: float = 4.0,
        require_spatial_resolution: bool = True,
    ) -> ResolvedSourceConfigFields3D:
        """Regularise the plane as a moment-matched, grid-thick source slab.

        The radial Gaussian matches the second radial moment of a uniform disk,
        giving ``sigma = radius/2 = diameter/4``.  The axial Gaussian matches
        the variance of one uniform grid-cell slab, ``dx/sqrt(12)``.  Arbitrary
        rotations are deliberately rejected because the current source
        quadrature is exactly cell-integrated only when the plane normal is a
        Cartesian axis.
        """
        if not self.numerically_admissible:
            raise RuntimeError("low-Mach handoff did not pass its conservation screen")
        if not isinstance(nozzle_position_m, tuple) or len(nozzle_position_m) != 3:
            raise TypeError("nozzle_position_m must be a three-value tuple")
        origin = tuple(float(value) for value in nozzle_position_m)
        if any(not math.isfinite(value) for value in origin):
            raise ValueError("nozzle_position_m must be finite")
        spacing = _positive_finite(grid_spacing_m, "grid_spacing_m")
        duration = _positive_finite(duration_s, "duration_s")
        minimum_cells = _positive_finite(
            minimum_cells_across_diameter, "minimum_cells_across_diameter",
        )

        axis_index = max(range(3), key=lambda index: abs(self.axis_unit_vector[index]))
        if not math.isclose(
            abs(self.axis_unit_vector[axis_index]), 1.0,
            rel_tol=0.0, abs_tol=1.0e-12,
        ):
            raise ValueError(
                "resolved source quadrature currently requires a Cartesian-axis jet"
            )
        position = tuple(
            origin[index] + self.handoff_distance_m * self.axis_unit_vector[index]
            for index in range(3)
        )
        axial_sigma = spacing / math.sqrt(12.0)
        sigmas = [self.moment_matched_radial_sigma_m] * 3
        sigmas[axis_index] = axial_sigma
        cells_across = self.equivalent_diameter_m / spacing
        spatially_resolved = cells_across >= minimum_cells
        if require_spatial_resolution and not spatially_resolved:
            raise ValueError(
                "low-Mach source plane is under-resolved: "
                f"{cells_across:.3g} cells across diameter < {minimum_cells:.3g}"
            )
        momentum_velocity = tuple(
            self.source_specific_momentum_m_s * value
            for value in self.axis_unit_vector
        )
        schedule = SourceRateSchedule(
            time_s=(0.0, duration),
            rate_kg_s=(self.source.mass_flow_kg_s, 0.0),
            source_id=(
                f"{self.source.source_id}:low-mach-handoff-M"
                f"{self.maximum_mach_number:g}"
            ),
        )
        return ResolvedSourceConfigFields3D(
            source_schedule=schedule,
            source_position_m=position,
            source_sigma_m=self.moment_matched_radial_sigma_m,
            source_sigma_xyz_m=tuple(sigmas),
            source_momentum_velocity_m_s=momentum_velocity,
            source_temperature_k=self.equivalent_hydrogen_temperature_k,
            cells_across_equivalent_diameter=cells_across,
            spatially_resolved=spatially_resolved,
        )


def ambient_pressure_jet_from_energy_nozzle(
    nozzle: "EnergyConservingNotionalNozzle",
    *,
    source_id: str = "energy-conserving-notional-nozzle",
) -> AmbientPressureHydrogenJet:
    """Convert an admissible HyRAM/Yuceil--Otugen result without re-solving it."""
    if not nozzle.admissible:
        raise ValueError("energy-conserving notional nozzle is not admissible")
    assert nozzle.diameter is not None
    assert nozzle.density is not None
    assert nozzle.temperature is not None
    area = math.pi * nozzle.diameter**2 / 4.0
    reconstructed = nozzle.density * nozzle.velocity * area
    residual = abs(reconstructed - nozzle.mass_flow) / nozzle.mass_flow
    return AmbientPressureHydrogenJet(
        source_id=source_id,
        source_model="HyRAM+ Yuceil-Otugen energy-conserving notional nozzle",
        mass_flow_kg_s=nozzle.mass_flow,
        diameter_m=nozzle.diameter,
        velocity_m_s=nozzle.velocity,
        density_kg_m3=nozzle.density,
        temperature_k=nozzle.temperature,
        momentum_flux_n=nozzle.momentum_flux,
        total_specific_energy_j_kg=nozzle.throat.storage_enthalpy,
        relative_mass_residual=max(
            residual, float(nozzle.relative_mass_residual or 0.0),
        ),
    )


def declared_ambient_pressure_hydrogen_jet(
    *,
    source_id: str,
    source_model: str,
    mass_flow_kg_s: float,
    diameter_m: float,
    velocity_m_s: float,
    density_kg_m3: float,
    temperature_k: float,
    ambient_pressure_pa: float,
    momentum_flux_n: float | None = None,
    total_specific_energy_j_kg: float | None = None,
) -> AmbientPressureHydrogenJet:
    """Build a complete published ambient-pressure source with audit closure.

    A published diameter alone is intentionally insufficient.  Density,
    velocity and temperature (or their traceable equivalents) are all required
    so that mass, momentum and energy cannot be silently invented.
    """
    from CoolProp.CoolProp import PropsSI

    mass_flow = _positive_finite(mass_flow_kg_s, "mass_flow_kg_s")
    diameter = _positive_finite(diameter_m, "diameter_m")
    velocity = _positive_finite(velocity_m_s, "velocity_m_s")
    density = _positive_finite(density_kg_m3, "density_kg_m3")
    temperature = _positive_finite(temperature_k, "temperature_k")
    pressure = _positive_finite(ambient_pressure_pa, "ambient_pressure_pa")
    area = math.pi * diameter**2 / 4.0
    residual = abs(density * velocity * area - mass_flow) / mass_flow
    if residual > 0.05:
        raise ValueError(
            "declared ambient-pressure source is inconsistent with mass flow "
            f"({100.0 * residual:.3g}% residual)"
        )
    if momentum_flux_n is None:
        momentum_flux_n = mass_flow * velocity
    if total_specific_energy_j_kg is None:
        static_enthalpy = float(PropsSI(
            "H", "T", temperature, "P", pressure, "Hydrogen",
        ))
        total_specific_energy_j_kg = static_enthalpy + 0.5 * velocity**2
    return AmbientPressureHydrogenJet(
        source_id=source_id,
        source_model=source_model,
        mass_flow_kg_s=mass_flow,
        diameter_m=diameter,
        velocity_m_s=velocity,
        density_kg_m3=density,
        temperature_k=temperature,
        momentum_flux_n=float(momentum_flux_n),
        total_specific_energy_j_kg=float(total_specific_energy_j_kg),
        relative_mass_residual=residual,
    )


def low_mach_handoff_from_ambient_pressure_jet(
    source: AmbientPressureHydrogenJet,
    *,
    ambient_temperature_k: float,
    ambient_pressure_pa: float,
    ambient_density_kg_m3: float | None = None,
    ambient_axial_velocity_m_s: float = 0.0,
    axis_unit_vector: tuple[float, float, float] = (1.0, 0.0, 0.0),
    maximum_mach_number: float = 0.30,
    maximum_relative_density_departure: float = 0.10,
    entrainment_beta: float = 0.28,
    hydrogen_heat_capacity_j_kg_k: float | None = None,
    air_heat_capacity_j_kg_k: float | None = None,
) -> LowMachHydrogenHandoff3D:
    """Find the first top-hat mixed plane satisfying the low-Mach gate."""
    from CoolProp.CoolProp import PropsSI

    ambient_temperature = _positive_finite(
        ambient_temperature_k, "ambient_temperature_k",
    )
    ambient_pressure = _positive_finite(ambient_pressure_pa, "ambient_pressure_pa")
    ambient_velocity = float(ambient_axial_velocity_m_s)
    if not math.isfinite(ambient_velocity) or ambient_velocity < 0.0:
        raise ValueError("ambient_axial_velocity_m_s must be finite and non-negative")
    if ambient_velocity >= source.velocity_m_s:
        raise ValueError("ambient axial velocity must be below source velocity")
    maximum_mach = _positive_finite(maximum_mach_number, "maximum_mach_number")
    if maximum_mach > 0.30:
        raise ValueError("low-Mach handoff limit cannot exceed 0.30")
    maximum_density_departure = _positive_finite(
        maximum_relative_density_departure,
        "maximum_relative_density_departure",
    )
    if maximum_density_departure > 1.0:
        raise ValueError("maximum_relative_density_departure cannot exceed 1")
    beta = _positive_finite(entrainment_beta, "entrainment_beta")
    axis = _unit_vector(axis_unit_vector)
    rho_air = (
        float(PropsSI("D", "T", ambient_temperature, "P", ambient_pressure, "Air"))
        if ambient_density_kg_m3 is None
        else _positive_finite(ambient_density_kg_m3, "ambient_density_kg_m3")
    )
    cp_h2 = (
        float(PropsSI("C", "T", ambient_temperature, "P", ambient_pressure, "Hydrogen"))
        if hydrogen_heat_capacity_j_kg_k is None
        else _positive_finite(
            hydrogen_heat_capacity_j_kg_k, "hydrogen_heat_capacity_j_kg_k",
        )
    )
    cp_air = (
        float(PropsSI("C", "T", ambient_temperature, "P", ambient_pressure, "Air"))
        if air_heat_capacity_j_kg_k is None
        else _positive_finite(air_heat_capacity_j_kg_k, "air_heat_capacity_j_kg_k")
    )
    h2_ambient_enthalpy = float(PropsSI(
        "H", "T", ambient_temperature, "P", ambient_pressure, "Hydrogen",
    ))

    h2_flow = source.mass_flow_kg_s
    source_energy_increment = h2_flow * (
        source.total_specific_energy_j_kg - h2_ambient_enthalpy
    )
    source_momentum_increment = source.momentum_flux_n
    reference_density = (
        ambient_pressure * _AIR_MOLAR_MASS_KG_MOL
        / (R_UNIVERSAL * ambient_temperature)
    )

    def mixed_state(total_flow: float) -> tuple[float, ...]:
        air_flow = total_flow - h2_flow
        plane_momentum = source_momentum_increment + air_flow * ambient_velocity
        velocity = plane_momentum / total_flow
        plane_energy = (
            source_energy_increment + 0.5 * air_flow * ambient_velocity**2
        )
        kinetic = 0.5 * total_flow * velocity**2
        heat_capacity_flow = h2_flow * cp_h2 + air_flow * cp_air
        temperature = ambient_temperature + (plane_energy - kinetic) / heat_capacity_flow
        if not math.isfinite(temperature) or temperature <= 0.0:
            return (math.inf,) * 12
        h2_fraction = h2_flow / total_flow
        molar_mass = 1.0 / (
            h2_fraction / _H2_MOLAR_MASS_KG_MOL
            + (1.0 - h2_fraction) / _AIR_MOLAR_MASS_KG_MOL
        )
        cp_mixture = heat_capacity_flow / total_flow
        gas_constant = R_UNIVERSAL / molar_mass
        cv_mixture = cp_mixture - gas_constant
        if cv_mixture <= 0.0:
            return (math.inf,) * 12
        gamma = cp_mixture / cv_mixture
        sound_speed = math.sqrt(gamma * gas_constant * temperature)
        mach = velocity / sound_speed
        density = ambient_pressure * molar_mass / (R_UNIVERSAL * temperature)
        density_departure = (density - reference_density) / reference_density
        sensible = plane_energy - kinetic
        return (
            mach, air_flow, plane_momentum, velocity, temperature, molar_mass,
            cp_mixture, sound_speed, density, density_departure, sensible,
            plane_energy,
        )

    def outside_receiving_limits(state: tuple[float, ...]) -> bool:
        return (
            state[0] > maximum_mach
            or abs(state[9]) > maximum_density_departure
        )

    lower = h2_flow
    lower_state = mixed_state(lower)
    if not outside_receiving_limits(lower_state):
        total_flow = lower
    else:
        upper = 2.0 * lower
        for _ in range(80):
            if not outside_receiving_limits(mixed_state(upper)):
                break
            upper *= 2.0
        else:
            raise RuntimeError("could not bracket a physical low-Mach mixed plane")
        for _ in range(100):
            midpoint = 0.5 * (lower + upper)
            if outside_receiving_limits(mixed_state(midpoint)):
                lower = midpoint
            else:
                upper = midpoint
        total_flow = upper

    (
        mach, air_flow, plane_momentum, velocity, temperature, molar_mass,
        cp_mixture, sound_speed, density, density_departure, sensible,
        plane_energy,
    ) = mixed_state(total_flow)
    area = total_flow / (density * velocity)
    diameter = math.sqrt(4.0 * area / math.pi)
    relative_velocity = source.velocity_m_s - ambient_velocity
    entrainment_per_length = beta * math.sqrt(
        rho_air * source.density_kg_m3
        * (math.pi * source.diameter_m**2 / 4.0)
        * relative_velocity**2
    )
    distance = air_flow / entrainment_per_length
    equivalent_h2_temperature = ambient_temperature + sensible / (h2_flow * cp_h2)
    if not 14.0 <= equivalent_h2_temperature <= 1_000.0:
        raise RuntimeError(
            "equivalent resolved-source temperature lies outside the active-scalar bounds"
        )

    hydrogen_residual = abs(h2_flow - source.mass_flow_kg_s) / h2_flow
    total_mass_residual = abs(total_flow - h2_flow - air_flow) / total_flow
    momentum_residual = abs(
        total_flow * velocity - plane_momentum
    ) / max(abs(plane_momentum), 1.0)
    energy_residual = abs(
        sensible + 0.5 * total_flow * velocity**2 - plane_energy
    ) / max(abs(plane_energy), 1.0)
    return LowMachHydrogenHandoff3D(
        source=source,
        axis_unit_vector=axis,
        ambient_temperature_k=ambient_temperature,
        ambient_pressure_pa=ambient_pressure,
        ambient_density_kg_m3=rho_air,
        ambient_axial_velocity_m_s=ambient_velocity,
        maximum_mach_number=maximum_mach,
        maximum_relative_density_departure=maximum_density_departure,
        entrainment_beta=beta,
        handoff_distance_m=distance,
        entrained_air_mass_flow_kg_s=air_flow,
        total_mass_flow_kg_s=total_flow,
        hydrogen_mass_fraction=h2_flow / total_flow,
        mixture_velocity_m_s=velocity,
        mixture_temperature_k=temperature,
        mixture_density_kg_m3=density,
        relative_density_departure=density_departure,
        mixture_molar_mass_kg_mol=molar_mass,
        mixture_heat_capacity_j_kg_k=cp_mixture,
        mixture_sound_speed_m_s=sound_speed,
        mach_number=mach,
        equivalent_diameter_m=diameter,
        moment_matched_radial_sigma_m=diameter / 4.0,
        plane_momentum_flux_n=plane_momentum,
        source_increment_momentum_flux_n=source_momentum_increment,
        plane_sensible_enthalpy_flux_w=sensible,
        source_increment_total_energy_flux_w=source_energy_increment,
        equivalent_hydrogen_temperature_k=equivalent_h2_temperature,
        relative_hydrogen_mass_residual=hydrogen_residual,
        relative_total_mass_residual=total_mass_residual,
        relative_momentum_residual=momentum_residual,
        relative_total_energy_residual=energy_residual,
        source_mass_consistency_residual=source.relative_mass_residual,
    )


def low_mach_handoff_from_energy_nozzle(
    nozzle: "EnergyConservingNotionalNozzle",
    **kwargs: object,
) -> LowMachHydrogenHandoff3D:
    """Convenience route from the existing energy-conserving nozzle model."""
    source_id = str(kwargs.pop("source_id", "energy-conserving-notional-nozzle"))
    source = ambient_pressure_jet_from_energy_nozzle(nozzle, source_id=source_id)
    return low_mach_handoff_from_ambient_pressure_jet(source, **kwargs)


__all__ = [
    "AmbientPressureHydrogenJet",
    "LowMachHydrogenHandoff3D",
    "ResolvedSourceConfigFields3D",
    "ambient_pressure_jet_from_energy_nozzle",
    "declared_ambient_pressure_hydrogen_jet",
    "low_mach_handoff_from_ambient_pressure_jet",
    "low_mach_handoff_from_energy_nozzle",
]
