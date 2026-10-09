"""Opt-in LH2 droplet crosswind transport and post-rainout pool coupling.

The ordinary DEGALI source path deliberately collapses the residual liquid at
the first all-vapour plane.  This module provides the complementary branch for
releases where droplets may leave the jet, evaporate in flight, or reach the
ground.  It keeps those three hydrogen streams in one explicit mass ledger.

The trajectory closure is intentionally modest: spherical droplets,
Schiller--Naumann drag, gravity/buoyancy, a declared mean wind vector, and a
declared d-squared evaporation coefficient.  It does not infer a coefficient
from concentration data and it does not claim to resolve turbulent droplet
dispersion. Deposited liquid can be sent either to the legacy fixed-footprint
heat balance or to the conservative axisymmetric dynamic-pool module.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import TYPE_CHECKING, Sequence

import numpy as np

from .pool_evaporation import (
    PoolEvaporationResult,
    SolidSubstrate,
    substrate_conduction_evaporation,
)
from .dynamic_pool import (
    DynamicPoolNumerics,
    DynamicPoolResult,
    simulate_axisymmetric_spreading_pool,
)

if TYPE_CHECKING:
    from .lh2_droplets import FlashingHydrogenDropletSource


GRAVITY_M_S2 = 9.80665


@dataclass(frozen=True)
class DropletClass:
    """One declared diameter class and its share of the liquid mass rate."""

    diameter_m: float
    mass_fraction: float


@dataclass(frozen=True)
class DropletTransportInput:
    """Boundary conditions shared by all classes in a spray population."""

    liquid_mass_flow_kg_s: float
    classes: Sequence[DropletClass]
    liquid_density_kg_m3: float
    initial_velocity_m_s: tuple[float, float, float]
    release_position_m: tuple[float, float, float]
    air_density_kg_m3: float
    air_kinematic_viscosity_m2_s: float
    evaporation_coefficient_m2_s: float
    schmidt_number: float
    wind_velocity_m_s: tuple[float, float, float]
    latent_heat_j_kg: float | None = None
    environmental_heat_input_w: float | None = None


@dataclass(frozen=True)
class DropletEvaporationSegment:
    """A mass-weighted section of one class's resolved evaporation path.

    ``vapour_mass_fraction`` is relative to that class's inlet liquid mass
    flow.  The segment is an Euler-trajectory record, not a fitted plume
    source or an assertion of turbulent droplet dispersion.
    """

    start_time_s: float
    end_time_s: float
    start_position_m: tuple[float, float, float]
    end_position_m: tuple[float, float, float]
    vapour_mass_fraction: float
    vapour_centroid_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not (
            math.isfinite(self.start_time_s)
            and math.isfinite(self.end_time_s)
            and self.end_time_s > self.start_time_s
        ):
            raise ValueError("evaporation segment time bounds must be finite and increasing")
        for name, position in {
            "start_position_m": self.start_position_m,
            "end_position_m": self.end_position_m,
            "vapour_centroid_m": self.vapour_centroid_m,
        }.items():
            if len(position) != 3 or not all(math.isfinite(float(value)) for value in position):
                raise ValueError(f"{name} must contain three finite values")
        if not math.isfinite(self.vapour_mass_fraction) or self.vapour_mass_fraction <= 0.0:
            raise ValueError("vapour_mass_fraction must be positive and finite")


@dataclass(frozen=True)
class DropletClassOutcome:
    """Endpoint and mass split for one injected droplet class."""

    diameter_m: float
    inlet_mass_flow_kg_s: float
    airborne_vapour_mass_flow_kg_s: float
    airborne_liquid_mass_flow_kg_s: float
    ground_liquid_mass_flow_kg_s: float
    terminal_position_m: tuple[float, float, float]
    terminal_velocity_m_s: tuple[float, float, float]
    terminal_time_s: float
    final_diameter_m: float
    maximum_reynolds_number: float
    status: str
    evaporation_segments: tuple[DropletEvaporationSegment, ...] = ()


@dataclass(frozen=True)
class DropletPopulationResult:
    """Population-integrated airborne/rainout source and conservation audit."""

    outcomes: tuple[DropletClassOutcome, ...]
    inlet_liquid_mass_flow_kg_s: float
    airborne_vapour_mass_flow_kg_s: float
    airborne_liquid_mass_flow_kg_s: float
    ground_liquid_mass_flow_kg_s: float
    mass_residual_kg_s: float
    evaporation_heat_requirement_w: float | None
    energy_residual_w: float | None
    impact_centroid_m: tuple[float, float, float] | None
    impact_rms_radius_m: float | None
    qualifications: tuple[str, ...]


@dataclass(frozen=True)
class RainoutPoolCouplingResult:
    """Post-release fixed-footprint pool result with an end-to-end H2 ledger."""

    droplets: DropletPopulationResult
    pool: PoolEvaporationResult | None
    direct_vapour_mass_kg: float
    airborne_droplet_vapour_mass_kg: float
    airborne_liquid_mass_kg: float
    deposited_liquid_mass_kg: float
    pool_vapour_mass_kg: float
    remaining_pool_liquid_mass_kg: float
    hydrogen_mass_residual_kg: float
    qualification: str
    escaped_pool_liquid_mass_kg: float = 0.0


def _drag_coefficient(reynolds: float) -> float:
    if reynolds <= 0.0:
        return 0.0
    if reynolds <= 1000.0:
        return 24.0 / reynolds * (1.0 + 0.15 * reynolds**0.687)
    return 0.44


def _finite_vector(values: Sequence[float], name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain three finite values")
    return vector


def _validate_input(boundary: DropletTransportInput, max_time_s: float) -> None:
    positive = {
        "liquid_mass_flow_kg_s": boundary.liquid_mass_flow_kg_s,
        "liquid_density_kg_m3": boundary.liquid_density_kg_m3,
        "air_density_kg_m3": boundary.air_density_kg_m3,
        "air_kinematic_viscosity_m2_s": boundary.air_kinematic_viscosity_m2_s,
        "evaporation_coefficient_m2_s": boundary.evaporation_coefficient_m2_s,
        "schmidt_number": boundary.schmidt_number,
        "max_time_s": max_time_s,
    }
    for name, value in positive.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    _finite_vector(boundary.initial_velocity_m_s, "initial_velocity_m_s")
    position = _finite_vector(boundary.release_position_m, "release_position_m")
    _finite_vector(boundary.wind_velocity_m_s, "wind_velocity_m_s")
    if position[2] <= 0.0:
        raise ValueError("release height must be positive")
    classes = tuple(boundary.classes)
    if not classes:
        raise ValueError("at least one droplet class is required")
    for item in classes:
        if not math.isfinite(item.diameter_m) or item.diameter_m <= 0.0:
            raise ValueError("droplet diameters must be finite and positive")
        if not math.isfinite(item.mass_fraction) or item.mass_fraction <= 0.0:
            raise ValueError("droplet mass fractions must be finite and positive")
    if not math.isclose(
        sum(item.mass_fraction for item in classes), 1.0,
        rel_tol=1.0e-10, abs_tol=1.0e-12,
    ):
        raise ValueError("droplet-class mass fractions must sum to one")
    if (boundary.latent_heat_j_kg is None) != (
        boundary.environmental_heat_input_w is None
    ):
        raise ValueError("provide both thermal-accounting inputs, or neither")
    if boundary.latent_heat_j_kg is not None and (
        not math.isfinite(boundary.latent_heat_j_kg)
        or boundary.latent_heat_j_kg <= 0.0
    ):
        raise ValueError("latent_heat_j_kg must be finite and positive")
    if boundary.environmental_heat_input_w is not None and not math.isfinite(
        boundary.environmental_heat_input_w
    ):
        raise ValueError("environmental_heat_input_w must be finite")


def _integrate_class(
    boundary: DropletTransportInput,
    item: DropletClass,
    *,
    max_time_s: float,
    trajectory_segment_duration_s: float | None,
) -> DropletClassOutcome:
    position = _finite_vector(boundary.release_position_m, "release_position_m")
    velocity = _finite_vector(boundary.initial_velocity_m_s, "initial_velocity_m_s")
    wind = _finite_vector(boundary.wind_velocity_m_s, "wind_velocity_m_s")
    diameter = float(item.diameter_m)
    initial_diameter = diameter
    time_s = 0.0
    maximum_reynolds = 0.0
    # The continuum point-drop equations become singular at exactly zero.
    # This limit removes less than 1e-18 of the original class mass.
    extinction_diameter = max(1.0e-12, initial_diameter * 1.0e-6)
    evaporation_segments: list[DropletEvaporationSegment] = []
    segment_start_time_s = time_s
    segment_start_position = position.copy()
    segment_mass_fraction = 0.0
    segment_centroid_numerator = np.zeros(3, dtype=float)

    while (
        time_s < max_time_s
        and position[2] > 0.0
        and diameter > extinction_diameter
    ):
        slip = wind - velocity
        slip_speed = float(np.linalg.norm(slip))
        reynolds = (
            slip_speed * diameter / boundary.air_kinematic_viscosity_m2_s
        )
        maximum_reynolds = max(maximum_reynolds, reynolds)
        drag = _drag_coefficient(reynolds)
        drag_rate = (
            3.0 * drag * boundary.air_density_kg_m3 * slip_speed
            / (4.0 * boundary.liquid_density_kg_m3 * diameter)
        )
        evaporation_rate = (
            boundary.evaporation_coefficient_m2_s / diameter
            * (
                1.0
                + 0.28 * math.sqrt(reynolds)
                * boundary.schmidt_number ** (1.0 / 3.0)
            )
        )
        # Bound changes in diameter, velocity relaxation and position.  The
        # small absolute ceiling also makes the Euler endpoint reproducible.
        dt = min(
            1.0e-3,
            0.002 * diameter / max(evaporation_rate, 1.0e-30),
            0.02 / max(drag_rate, 1.0e-30),
            max_time_s - time_s,
        )
        old_position = position.copy()
        old_velocity = velocity.copy()
        old_diameter = diameter
        old_time_s = time_s
        acceleration = drag_rate * slip
        acceleration[2] -= GRAVITY_M_S2 * (
            1.0 - boundary.air_density_kg_m3 / boundary.liquid_density_kg_m3
        )
        velocity = velocity + acceleration * dt
        position = position + velocity * dt
        diameter = max(0.0, diameter - evaporation_rate * dt)
        time_s += dt

        if position[2] <= 0.0:
            # Linear endpoint interpolation prevents the time step from
            # creating a fictitious below-ground impact location.
            fraction = old_position[2] / max(
                old_position[2] - position[2], 1.0e-30
            )
            fraction = min(max(float(fraction), 0.0), 1.0)
            position = old_position + fraction * (position - old_position)
            velocity = old_velocity + fraction * (velocity - old_velocity)
            diameter = old_diameter + fraction * (diameter - old_diameter)
            time_s -= (1.0 - fraction) * dt
            position[2] = 0.0
        evaporated_fraction = max(
            0.0, (old_diameter**3 - diameter**3) / initial_diameter**3
        )
        if trajectory_segment_duration_s is not None:
            segment_mass_fraction += evaporated_fraction
            segment_centroid_numerator += evaporated_fraction * 0.5 * (
                old_position + position
            )
            if time_s - segment_start_time_s >= trajectory_segment_duration_s:
                if segment_mass_fraction > 0.0:
                    evaporation_segments.append(DropletEvaporationSegment(
                        start_time_s=float(segment_start_time_s),
                        end_time_s=float(time_s),
                        start_position_m=tuple(float(value) for value in segment_start_position),
                        end_position_m=tuple(float(value) for value in position),
                        vapour_mass_fraction=float(segment_mass_fraction),
                        vapour_centroid_m=tuple(
                            float(value) for value in (
                                segment_centroid_numerator / segment_mass_fraction
                            )
                        ),
                    ))
                segment_start_time_s = time_s
                segment_start_position = position.copy()
                segment_mass_fraction = 0.0
                segment_centroid_numerator = np.zeros(3, dtype=float)
        if position[2] <= 0.0:
            break

    if diameter <= extinction_diameter:
        diameter = 0.0
        status = "complete_evaporation_before_ground"
    elif position[2] <= 0.0:
        status = "ground_impact"
    else:
        status = "trajectory_time_limit"

    inlet = boundary.liquid_mass_flow_kg_s * item.mass_fraction
    remaining_fraction = (diameter / initial_diameter) ** 3
    remaining = inlet * remaining_fraction
    ground = remaining if status == "ground_impact" else 0.0
    airborne_liquid = remaining if status == "trajectory_time_limit" else 0.0
    vapour = inlet - ground - airborne_liquid
    if trajectory_segment_duration_s is not None and segment_mass_fraction > 0.0:
        evaporation_segments.append(DropletEvaporationSegment(
            start_time_s=float(segment_start_time_s),
            end_time_s=float(time_s),
            start_position_m=tuple(float(value) for value in segment_start_position),
            end_position_m=tuple(float(value) for value in position),
            vapour_mass_fraction=float(segment_mass_fraction),
            vapour_centroid_m=tuple(
                float(value) for value in (
                    segment_centroid_numerator / segment_mass_fraction
                )
            ),
        ))
    if trajectory_segment_duration_s is not None and evaporation_segments:
        recorded_fraction = sum(item.vapour_mass_fraction for item in evaporation_segments)
        required_fraction = vapour / inlet
        correction = required_fraction - recorded_fraction
        if abs(correction) > 1.0e-15:
            final = evaporation_segments[-1]
            evaporation_segments[-1] = replace(
                final,
                vapour_mass_fraction=final.vapour_mass_fraction + correction,
            )
    return DropletClassOutcome(
        diameter_m=initial_diameter,
        inlet_mass_flow_kg_s=inlet,
        airborne_vapour_mass_flow_kg_s=vapour,
        airborne_liquid_mass_flow_kg_s=airborne_liquid,
        ground_liquid_mass_flow_kg_s=ground,
        terminal_position_m=tuple(float(value) for value in position),
        terminal_velocity_m_s=tuple(float(value) for value in velocity),
        terminal_time_s=float(time_s),
        final_diameter_m=float(diameter),
        maximum_reynolds_number=float(maximum_reynolds),
        status=status,
        evaporation_segments=tuple(evaporation_segments),
    )


def transport_droplet_population(
    boundary: DropletTransportInput,
    *,
    max_time_s: float = 120.0,
    trajectory_segment_duration_s: float | None = None,
) -> DropletPopulationResult:
    """Transport classes and optionally retain conservative evaporation paths.

    With ``trajectory_segment_duration_s=None`` (the default), the historical
    endpoint-only result is retained.  A positive interval records only
    mass-weighted vapour-producing trajectory sections for an explicit
    downstream atmospheric handoff.
    """

    _validate_input(boundary, max_time_s)
    if trajectory_segment_duration_s is not None and (
        not math.isfinite(trajectory_segment_duration_s)
        or trajectory_segment_duration_s <= 0.0
    ):
        raise ValueError("trajectory_segment_duration_s must be positive and finite or None")
    outcomes = tuple(
        _integrate_class(
            boundary, item, max_time_s=max_time_s,
            trajectory_segment_duration_s=trajectory_segment_duration_s,
        )
        for item in boundary.classes
    )
    vapour = sum(item.airborne_vapour_mass_flow_kg_s for item in outcomes)
    airborne_liquid = sum(
        item.airborne_liquid_mass_flow_kg_s for item in outcomes
    )
    ground = sum(item.ground_liquid_mass_flow_kg_s for item in outcomes)
    residual = boundary.liquid_mass_flow_kg_s - vapour - airborne_liquid - ground
    impacts = [item for item in outcomes if item.ground_liquid_mass_flow_kg_s > 0.0]
    if impacts:
        centroid_array = sum(
            item.ground_liquid_mass_flow_kg_s
            * np.asarray(item.terminal_position_m)
            for item in impacts
        ) / ground
        rms = math.sqrt(sum(
            item.ground_liquid_mass_flow_kg_s
            * float(np.sum(
                (np.asarray(item.terminal_position_m) - centroid_array)[:2] ** 2
            ))
            for item in impacts
        ) / ground)
        centroid = tuple(float(value) for value in centroid_array)
    else:
        centroid = None
        rms = None
    heat = (
        vapour * boundary.latent_heat_j_kg
        if boundary.latent_heat_j_kg is not None else None
    )
    energy_residual = (
        boundary.environmental_heat_input_w - heat
        if heat is not None else None
    )
    maximum_reynolds = max(item.maximum_reynolds_number for item in outcomes)
    drag_scope = (
        "schiller_naumann_validated_to_re800"
        if maximum_reynolds <= 800.0 else
        "schiller_naumann_extrapolated_re800_to_1000"
        if maximum_reynolds <= 1000.0 else
        "constant_cd_above_re1000"
    )
    return DropletPopulationResult(
        outcomes=outcomes,
        inlet_liquid_mass_flow_kg_s=boundary.liquid_mass_flow_kg_s,
        airborne_vapour_mass_flow_kg_s=float(vapour),
        airborne_liquid_mass_flow_kg_s=float(airborne_liquid),
        ground_liquid_mass_flow_kg_s=float(ground),
        mass_residual_kg_s=float(residual),
        evaporation_heat_requirement_w=heat,
        energy_residual_w=energy_residual,
        impact_centroid_m=centroid,
        impact_rms_radius_m=rms,
        qualifications=(
            "declared_steady_mean_wind_vector",
            "declared_d2_evaporation_coefficient",
            drag_scope,
            "no_turbulent_droplet_dispersion",
            "impact_points_do_not_define_pool_area",
        ),
    )


def droplet_transport_input_from_flash(
    source: "FlashingHydrogenDropletSource",
    *,
    release_position_m: tuple[float, float, float],
    jet_direction: tuple[float, float, float] = (1.0, 0.0, 0.0),
    wind_velocity_m_s: tuple[float, float, float],
    air_density_kg_m3: float,
    air_kinematic_viscosity_m2_s: float,
    evaporation_coefficient_m2_s: float,
    schmidt_number: float,
    classes: Sequence[DropletClass] | None = None,
    latent_heat_j_kg: float | None = None,
    environmental_heat_input_w: float | None = None,
) -> DropletTransportInput:
    """Map the existing post-flash source into the droplet transport boundary.

    The default is the source model's published monodisperse breakup diameter.
    A caller may instead provide a normalized measured/correlated population;
    no distribution is manufactured here.
    """

    if source.liquid_mass_flow <= 0.0 or source.droplet_diameter is None:
        raise ValueError("post-flash source contains no residual liquid droplets")
    direction = _finite_vector(jet_direction, "jet_direction")
    norm = float(np.linalg.norm(direction))
    if norm <= 0.0:
        raise ValueError("jet_direction must be non-zero")
    declared_classes = (
        (DropletClass(source.droplet_diameter, 1.0),)
        if classes is None else tuple(classes)
    )
    return DropletTransportInput(
        liquid_mass_flow_kg_s=source.liquid_mass_flow,
        classes=declared_classes,
        liquid_density_kg_m3=source.liquid_density,
        initial_velocity_m_s=tuple(
            float(value) for value in source.postflash_velocity * direction / norm
        ),
        release_position_m=release_position_m,
        air_density_kg_m3=air_density_kg_m3,
        air_kinematic_viscosity_m2_s=air_kinematic_viscosity_m2_s,
        evaporation_coefficient_m2_s=evaporation_coefficient_m2_s,
        schmidt_number=schmidt_number,
        wind_velocity_m_s=wind_velocity_m_s,
        latent_heat_j_kg=latent_heat_j_kg,
        environmental_heat_input_w=environmental_heat_input_w,
    )


def post_release_rainout_pool(
    droplets: DropletPopulationResult,
    *,
    direct_vapour_mass_flow_kg_s: float,
    release_duration_s: float,
    pool_area_m2: float,
    pool_duration_s: float,
    pool_time_step_s: float,
    substrate: SolidSubstrate,
    saturation_temperature_k: float = 20.27,
    latent_heat_j_kg: float = 4.46e5,
) -> RainoutPoolCouplingResult:
    """Feed deposited liquid to the fixed-area pool after release cessation.

    The pool footprint must be supplied.  A handful of deterministic impact
    centroids cannot determine a turbulent deposition area, so this function
    refuses to invent one.  Deposition is accumulated over the declared
    release and the existing conduction model is then run as a post-release
    pool; concurrent inflow/spreading requires a dynamic pool solver.
    """

    if not isinstance(substrate, SolidSubstrate):
        raise TypeError("fixed-area pool requires SolidSubstrate")

    for name, value in {
        "direct_vapour_mass_flow_kg_s": direct_vapour_mass_flow_kg_s,
        "release_duration_s": release_duration_s,
        "pool_area_m2": pool_area_m2,
        "pool_duration_s": pool_duration_s,
        "pool_time_step_s": pool_time_step_s,
    }.items():
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if release_duration_s <= 0.0 or pool_area_m2 <= 0.0:
        raise ValueError("release duration and pool area must be positive")
    if pool_duration_s <= 0.0 or pool_time_step_s <= 0.0:
        raise ValueError("pool duration and time step must be positive")

    direct = direct_vapour_mass_flow_kg_s * release_duration_s
    airborne = droplets.airborne_vapour_mass_flow_kg_s * release_duration_s
    airborne_liquid = (
        droplets.airborne_liquid_mass_flow_kg_s * release_duration_s
    )
    deposited = droplets.ground_liquid_mass_flow_kg_s * release_duration_s
    pool = None
    pool_vapour = 0.0
    remaining = 0.0
    if deposited > 0.0:
        pool = substrate_conduction_evaporation(
            substrate,
            area_m2=pool_area_m2,
            duration_s=pool_duration_s,
            time_step_s=pool_time_step_s,
            saturation_temperature_k=saturation_temperature_k,
            latent_heat_j_kg=latent_heat_j_kg,
            initial_liquid_mass_kg=deposited,
        )
        pool_vapour = pool.total_evaporated_mass_kg
        remaining_value = pool.steps[-1].remaining_liquid_kg
        remaining = 0.0 if remaining_value is None else remaining_value
    inlet = (
        direct_vapour_mass_flow_kg_s + droplets.inlet_liquid_mass_flow_kg_s
    ) * release_duration_s
    outlet = direct + airborne + airborne_liquid + pool_vapour + remaining
    return RainoutPoolCouplingResult(
        droplets=droplets,
        pool=pool,
        direct_vapour_mass_kg=float(direct),
        airborne_droplet_vapour_mass_kg=float(airborne),
        airborne_liquid_mass_kg=float(airborne_liquid),
        deposited_liquid_mass_kg=float(deposited),
        pool_vapour_mass_kg=float(pool_vapour),
        remaining_pool_liquid_mass_kg=float(remaining),
        hydrogen_mass_residual_kg=float(inlet - outlet),
        qualification=(
            "post-release fixed-footprint conduction pool; deposition mass is "
            "conserved, but concurrent rainout and pool spreading are outside "
            "this model"
        ),
    )


def concurrent_rainout_pool(
    droplets: DropletPopulationResult,
    *,
    direct_vapour_mass_flow_kg_s: float,
    release_duration_s: float,
    post_release_duration_s: float,
    pool_area_m2: float,
    pool_time_step_s: float,
    substrate: SolidSubstrate,
    saturation_temperature_k: float = 20.27,
    latent_heat_j_kg: float = 4.46e5,
) -> RainoutPoolCouplingResult:
    """Form and evaporate a fixed-footprint pool during continuous rainout.

    Ground-liquid deposition is applied at the class-integrated rate for the
    source duration.  The pool is then followed for the requested additional
    time.  This closes the concurrency gap in the post-release approximation,
    while pool spreading and a time-varying deposition footprint remain
    explicitly outside the model.
    """

    if not isinstance(substrate, SolidSubstrate):
        raise TypeError("fixed-area pool requires SolidSubstrate")

    for name, value in {
        "direct_vapour_mass_flow_kg_s": direct_vapour_mass_flow_kg_s,
        "release_duration_s": release_duration_s,
        "post_release_duration_s": post_release_duration_s,
        "pool_area_m2": pool_area_m2,
        "pool_time_step_s": pool_time_step_s,
    }.items():
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if release_duration_s <= 0.0 or pool_area_m2 <= 0.0:
        raise ValueError("release duration and pool area must be positive")
    if pool_time_step_s <= 0.0:
        raise ValueError("pool time step must be positive")

    direct = direct_vapour_mass_flow_kg_s * release_duration_s
    airborne = droplets.airborne_vapour_mass_flow_kg_s * release_duration_s
    airborne_liquid = (
        droplets.airborne_liquid_mass_flow_kg_s * release_duration_s
    )
    deposited = droplets.ground_liquid_mass_flow_kg_s * release_duration_s
    pool = None
    pool_vapour = 0.0
    remaining = 0.0
    if deposited > 0.0:
        pool = substrate_conduction_evaporation(
            substrate,
            area_m2=pool_area_m2,
            duration_s=release_duration_s + post_release_duration_s,
            time_step_s=pool_time_step_s,
            saturation_temperature_k=saturation_temperature_k,
            latent_heat_j_kg=latent_heat_j_kg,
            initial_liquid_mass_kg=0.0,
            liquid_inflow_rate_kg_s=droplets.ground_liquid_mass_flow_kg_s,
            inflow_duration_s=release_duration_s,
        )
        pool_vapour = pool.total_evaporated_mass_kg
        remaining_value = pool.steps[-1].remaining_liquid_kg
        remaining = 0.0 if remaining_value is None else remaining_value
    inlet = (
        direct_vapour_mass_flow_kg_s + droplets.inlet_liquid_mass_flow_kg_s
    ) * release_duration_s
    outlet = direct + airborne + airborne_liquid + pool_vapour + remaining
    return RainoutPoolCouplingResult(
        droplets=droplets,
        pool=pool,
        direct_vapour_mass_kg=float(direct),
        airborne_droplet_vapour_mass_kg=float(airborne),
        airborne_liquid_mass_kg=float(airborne_liquid),
        deposited_liquid_mass_kg=float(deposited),
        pool_vapour_mass_kg=float(pool_vapour),
        remaining_pool_liquid_mass_kg=float(remaining),
        hydrogen_mass_residual_kg=float(inlet - outlet),
        qualification=(
            "concurrent constant-rate rainout into a declared fixed footprint; "
            "mass and substrate-conduction evaporation are coupled, but pool "
            "spreading and a time-varying footprint are outside this model"
        ),
    )


def dynamic_rainout_pool(
    droplets: DropletPopulationResult,
    *,
    direct_vapour_mass_flow_kg_s: float,
    release_duration_s: float,
    post_release_duration_s: float,
    pool_time_step_s: float,
    substrate: SolidSubstrate,
    liquid_density_kg_m3: float,
    liquid_kinematic_viscosity_m2_s: float,
    numerics: DynamicPoolNumerics,
    saturation_temperature_k: float = 20.27,
    latent_heat_j_kg: float = 4.46e5,
    evaporation_momentum_closure: str = "zero_radial_momentum_vapor",
    solid_heat_flux_multiplier: float = 1.0,
    solid_heat_flux_cap_w_m2: float | None = None,
) -> RainoutPoolCouplingResult:
    """Form, spread and evaporate a rainout pool during a finite release.

    Ground-liquid deposition is applied over the footprint declared by
    ``numerics``.  The pool then spreads beyond that footprint under the
    axisymmetric shallow-layer equations.  Domain escape and numerical
    positivity corrections remain explicit in the returned conservation
    ledger.
    """

    for name, value in {
        "direct_vapour_mass_flow_kg_s": direct_vapour_mass_flow_kg_s,
        "release_duration_s": release_duration_s,
        "post_release_duration_s": post_release_duration_s,
        "pool_time_step_s": pool_time_step_s,
        "liquid_density_kg_m3": liquid_density_kg_m3,
        "liquid_kinematic_viscosity_m2_s": (
            liquid_kinematic_viscosity_m2_s
        ),
    }.items():
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if release_duration_s <= 0.0 or pool_time_step_s <= 0.0:
        raise ValueError("release duration and pool time step must be positive")
    if liquid_density_kg_m3 <= 0.0:
        raise ValueError("liquid density must be positive")

    direct = direct_vapour_mass_flow_kg_s * release_duration_s
    airborne = droplets.airborne_vapour_mass_flow_kg_s * release_duration_s
    airborne_liquid = (
        droplets.airborne_liquid_mass_flow_kg_s * release_duration_s
    )
    deposited = droplets.ground_liquid_mass_flow_kg_s * release_duration_s
    pool: DynamicPoolResult | None = None
    pool_vapour = 0.0
    remaining = 0.0
    escaped = 0.0
    if deposited > 0.0:
        pool = simulate_axisymmetric_spreading_pool(
            substrate,
            liquid_inflow_rate_kg_s=(
                droplets.ground_liquid_mass_flow_kg_s
            ),
            inflow_duration_s=release_duration_s,
            duration_s=release_duration_s + post_release_duration_s,
            output_time_step_s=pool_time_step_s,
            liquid_density_kg_m3=liquid_density_kg_m3,
            liquid_kinematic_viscosity_m2_s=(
                liquid_kinematic_viscosity_m2_s
            ),
            saturation_temperature_k=saturation_temperature_k,
            latent_heat_j_kg=latent_heat_j_kg,
            numerics=numerics,
            evaporation_momentum_closure=evaporation_momentum_closure,
            solid_heat_flux_multiplier=solid_heat_flux_multiplier,
            solid_heat_flux_cap_w_m2=solid_heat_flux_cap_w_m2,
        )
        pool_vapour = pool.total_evaporated_mass_kg
        remaining = pool.remaining_liquid_mass_kg
        escaped = pool.escaped_domain_liquid_mass_kg

    inlet = (
        direct_vapour_mass_flow_kg_s + droplets.inlet_liquid_mass_flow_kg_s
    ) * release_duration_s
    outlet = (
        direct + airborne + airborne_liquid + pool_vapour + remaining + escaped
    )
    return RainoutPoolCouplingResult(
        droplets=droplets,
        pool=pool,
        direct_vapour_mass_kg=float(direct),
        airborne_droplet_vapour_mass_kg=float(airborne),
        airborne_liquid_mass_kg=float(airborne_liquid),
        deposited_liquid_mass_kg=float(deposited),
        pool_vapour_mass_kg=float(pool_vapour),
        remaining_pool_liquid_mass_kg=float(remaining),
        hydrogen_mass_residual_kg=float(inlet - outlet),
        qualification=(
            "concurrent rainout into a declared footprint with conservative "
            "axisymmetric shallow-layer spreading and local-contact-age "
            f"substrate evaporation ({evaporation_momentum_closure}); slope, "
            "obstacles, drains, wind shear on "
            "the liquid and non-axisymmetric footprints are outside scope"
        ),
        escaped_pool_liquid_mass_kg=float(escaped),
    )


__all__ = [
    "DropletClass", "DropletTransportInput", "DropletEvaporationSegment",
    "DropletClassOutcome",
    "DropletPopulationResult", "RainoutPoolCouplingResult",
    "transport_droplet_population", "droplet_transport_input_from_flash",
    "post_release_rainout_pool", "concurrent_rainout_pool",
    "dynamic_rainout_pool",
]
