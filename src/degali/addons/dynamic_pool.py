"""Conservative axisymmetric spreading and evaporation of a rainout pool.

This module closes the fixed-footprint limitation of :mod:`pool_evaporation`.
It solves the depth-averaged radial shallow-layer equations with a first-order
Rusanov finite-volume flux.  Liquid is deposited over a declared disk or
annulus, spreads under gravity, loses momentum through declared bed/viscous
friction, and evaporates through independent semi-infinite substrate columns.

The formulation follows the shallow-layer route documented by Dienhart
(JUEL-3155, 1995).  It deliberately remains axisymmetric and horizontal: slope,
bund walls, drains, wind shear on the liquid, boiling-regime transitions and
water/ice dynamics require different boundary models and are not inferred.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .pool_evaporation import SolidSubstrate


GRAVITY_M_S2 = 9.80665


@dataclass(frozen=True)
class ConstantHeatFluxSurface:
    """Declared short-duration water/ice heat-flux boundary.

    This is an experiment-specific boundary condition, not a generic water
    property.  Ice growth, waves, detached floes and pool-front pulsation are
    not resolved by the axisymmetric layer model.
    """

    heat_flux_w_m2: float
    calibration_label: str

    def __post_init__(self) -> None:
        if not math.isfinite(self.heat_flux_w_m2) or self.heat_flux_w_m2 <= 0.0:
            raise ValueError("heat_flux_w_m2 must be finite and positive")
        if not self.calibration_label.strip():
            raise ValueError("calibration_label must identify the evidence source")


@dataclass(frozen=True)
class DynamicPoolNumerics:
    """Grid, source-footprint and front-observation settings in SI units."""

    domain_radius_m: float = 10.0
    radial_step_m: float = 0.02
    maximum_time_step_s: float = 0.01
    cfl: float = 0.2
    chezy_coefficient: float = 1.0e-3
    source_inner_radius_m: float = 0.0
    source_radius_m: float = 0.25
    surface_retention_depth_m: float = 0.0
    dry_depth_m: float = 1.0e-8
    reported_front_depth_m: float = 7.5e-4

    def __post_init__(self) -> None:
        for name in (
            "domain_radius_m", "radial_step_m", "maximum_time_step_s",
            "cfl", "source_radius_m", "dry_depth_m",
            "reported_front_depth_m",
        ):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        if self.domain_radius_m <= self.source_radius_m:
            raise ValueError("domain_radius_m must exceed source_radius_m")
        if not 0.0 < self.cfl <= 0.5:
            raise ValueError("cfl must be in (0, 0.5]")
        if (
            not math.isfinite(self.source_inner_radius_m)
            or self.source_inner_radius_m < 0.0
            or self.source_inner_radius_m >= self.source_radius_m
        ):
            raise ValueError(
                "source_inner_radius_m must be in [0, source_radius_m)"
            )
        for name in ("chezy_coefficient", "surface_retention_depth_m"):
            value = getattr(self, name)
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")


@dataclass(frozen=True)
class DynamicPoolStep:
    """One conservation-ledger row for the spreading pool."""

    elapsed_s: float
    liquid_inflow_rate_kg_s: float
    cumulative_liquid_inflow_kg: float
    liquid_mass_kg: float
    reported_radius_m: float
    reported_area_m2: float
    wet_area_m2: float
    wet_area_mean_depth_m: float
    evaporation_rate_kg_s: float
    cumulative_evaporated_mass_kg: float
    escaped_domain_liquid_mass_kg: float
    cumulative_numerical_mass_adjustment_kg: float
    mass_residual_kg: float

    @property
    def vapour_rate_kg_s(self) -> float:
        """Compatibility alias for quasi-steady pool-plume observation tools."""
        return self.evaporation_rate_kg_s

    @property
    def evaporated_mass_kg(self) -> float:
        return self.cumulative_evaporated_mass_kg

    @property
    def remaining_liquid_kg(self) -> float:
        return self.liquid_mass_kg


@dataclass(frozen=True)
class DynamicPoolResult:
    """Time history, final radial profile and audit of a dynamic LH2 pool."""

    steps: tuple[DynamicPoolStep, ...]
    radial_centres_m: tuple[float, ...]
    final_depth_m: tuple[float, ...]
    saturation_temperature_k: float
    liquid_density_kg_m3: float
    latent_heat_j_kg: float
    substrate: SolidSubstrate | ConstantHeatFluxSurface
    numerics: DynamicPoolNumerics
    evaporation_momentum_closure: str
    solid_heat_flux_multiplier: float
    solid_heat_flux_cap_w_m2: float | None
    solver_id: str = "axisymmetric_rusanov_shallow_layer_local_contact_v1"

    @property
    def total_evaporated_mass_kg(self) -> float:
        return self.steps[-1].cumulative_evaporated_mass_kg

    @property
    def total_liquid_inflow_kg(self) -> float:
        return self.steps[-1].cumulative_liquid_inflow_kg

    @property
    def mean_vapour_rate_kg_s(self) -> float:
        duration = self.steps[-1].elapsed_s
        return self.total_evaporated_mass_kg / duration if duration > 0.0 else 0.0

    @property
    def remaining_liquid_mass_kg(self) -> float:
        return self.steps[-1].liquid_mass_kg

    @property
    def escaped_domain_liquid_mass_kg(self) -> float:
        return self.steps[-1].escaped_domain_liquid_mass_kg

    @property
    def maximum_absolute_mass_residual_kg(self) -> float:
        return max(abs(step.mass_residual_kg) for step in self.steps)


def _radial_grid(
    numerics: DynamicPoolNumerics,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    count = int(math.ceil(numerics.domain_radius_m / numerics.radial_step_m))
    faces = np.arange(count + 1, dtype=float) * numerics.radial_step_m
    centres = (np.arange(count, dtype=float) + 0.5) * numerics.radial_step_m
    areas = math.pi * (faces[1:] ** 2 - faces[:-1] ** 2)
    return centres, faces, areas


def _source_shape(
    centres: np.ndarray,
    areas: np.ndarray,
    numerics: DynamicPoolNumerics,
) -> np.ndarray:
    if numerics.source_inner_radius_m == 0.0:
        shape = np.maximum(
            0.0, 1.0 - (centres / numerics.source_radius_m) ** 2
        )
    else:
        shape = (
            (centres >= numerics.source_inner_radius_m)
            & (centres <= numerics.source_radius_m)
        ).astype(float)
    integral = float(np.dot(shape, areas))
    if integral <= 0.0:
        raise ValueError(
            "source footprint is unresolved by the radial grid; refine the grid"
        )
    return shape / integral


def _stable_step(
    depth: np.ndarray,
    momentum: np.ndarray,
    numerics: DynamicPoolNumerics,
) -> float:
    mobile = np.maximum(depth - numerics.surface_retention_depth_m, 0.0)
    velocity = np.divide(
        momentum, mobile, out=np.zeros_like(depth),
        where=mobile > numerics.dry_depth_m,
    )
    signal = np.abs(velocity) + np.sqrt(GRAVITY_M_S2 * mobile)
    maximum = float(np.max(signal))
    if maximum <= 1.0e-12:
        return numerics.maximum_time_step_s
    return min(
        numerics.maximum_time_step_s,
        numerics.cfl * numerics.radial_step_m / maximum,
    )


def _flow_step(
    depth: np.ndarray,
    momentum: np.ndarray,
    *,
    dt_s: float,
    centres: np.ndarray,
    faces: np.ndarray,
    areas: np.ndarray,
    liquid_kinematic_viscosity_m2_s: float,
    numerics: DynamicPoolNumerics,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Advance radial shallow-layer flow and expose all mass corrections."""

    mobile = np.maximum(depth - numerics.surface_retention_depth_m, 0.0)
    velocity = np.divide(
        momentum, mobile, out=np.zeros_like(depth),
        where=mobile > numerics.dry_depth_m,
    )
    count = len(depth)
    volume_flux = np.zeros(count + 1)
    momentum_flux = np.zeros(count + 1)
    if count > 1:
        h_l, h_r = mobile[:-1], mobile[1:]
        m_l, m_r = momentum[:-1], momentum[1:]
        u_l, u_r = velocity[:-1], velocity[1:]
        signal = np.maximum(
            np.abs(u_l) + np.sqrt(GRAVITY_M_S2 * h_l),
            np.abs(u_r) + np.sqrt(GRAVITY_M_S2 * h_r),
        )
        volume_flux[1:count] = (
            0.5 * (m_l + m_r) - 0.5 * signal * (h_r - h_l)
        )
        momentum_flux[1:count] = 0.5 * (
            m_l * u_l + 0.5 * GRAVITY_M_S2 * h_l**2
            + m_r * u_r + 0.5 * GRAVITY_M_S2 * h_r**2
        ) - 0.5 * signal * (m_r - m_l)

    outer_signal = abs(velocity[-1]) + math.sqrt(GRAVITY_M_S2 * mobile[-1])
    volume_flux[-1] = (
        0.5 * momentum[-1] + 0.5 * outer_signal * mobile[-1]
    )
    momentum_flux[-1] = 0.5 * (
        momentum[-1] * velocity[-1]
        + 0.5 * GRAVITY_M_S2 * mobile[-1] ** 2
    ) + 0.5 * outer_signal * momentum[-1]

    dr = numerics.radial_step_m
    new_depth = depth - dt_s / (centres * dr) * (
        faces[1:] * volume_flux[1:] - faces[:-1] * volume_flux[:-1]
    )
    new_momentum = momentum - dt_s / (centres * dr) * (
        faces[1:] * momentum_flux[1:] - faces[:-1] * momentum_flux[:-1]
    )
    # Cylindrical geometric source balances the radial pressure-flux term.
    new_momentum += dt_s * 0.5 * GRAVITY_M_S2 * mobile**2 / centres
    escaped_volume = max(
        0.0, 2.0 * math.pi * faces[-1] * volume_flux[-1] * dt_s
    )

    negative = new_depth < 0.0
    numerical_adjustment = float(
        np.dot(np.maximum(-new_depth, 0.0), areas)
    )
    new_depth[negative] = 0.0
    new_momentum[negative] = 0.0

    new_mobile = np.maximum(
        new_depth - numerics.surface_retention_depth_m, 0.0
    )
    new_velocity = np.divide(
        new_momentum, new_mobile, out=np.zeros_like(new_depth),
        where=new_mobile > numerics.dry_depth_m,
    )
    damping = 1.0 + dt_s * (
        numerics.chezy_coefficient * np.abs(new_velocity)
        / np.maximum(new_mobile, numerics.dry_depth_m)
        + liquid_kinematic_viscosity_m2_s
        / np.maximum(new_mobile, numerics.dry_depth_m) ** 2
    )
    new_momentum = new_mobile * new_velocity / damping
    new_momentum[new_mobile <= numerics.dry_depth_m] = 0.0
    return new_depth, new_momentum, escaped_volume, numerical_adjustment


def simulate_axisymmetric_spreading_pool(
    substrate: SolidSubstrate | ConstantHeatFluxSurface,
    *,
    liquid_inflow_rate_kg_s: float,
    inflow_duration_s: float,
    duration_s: float,
    output_time_step_s: float,
    liquid_density_kg_m3: float,
    liquid_kinematic_viscosity_m2_s: float,
    saturation_temperature_k: float = 20.27,
    latent_heat_j_kg: float = 4.46e5,
    numerics: DynamicPoolNumerics = DynamicPoolNumerics(),
    evaporation_momentum_closure: str = "zero_radial_momentum_vapor",
    solid_heat_flux_multiplier: float = 1.0,
    solid_heat_flux_cap_w_m2: float | None = None,
) -> DynamicPoolResult:
    """Integrate concurrent deposition, spreading and substrate evaporation.

    The heat flux singularity at first wetting is integrated analytically for
    every radial cell.  Evaporation removes local liquid and proportionally
    removes its radial momentum.  The output ledger includes domain escape and
    positivity clipping, so neither can silently masquerade as evaporation.
    """

    positive = {
        "liquid_inflow_rate_kg_s": liquid_inflow_rate_kg_s,
        "inflow_duration_s": inflow_duration_s,
        "duration_s": duration_s,
        "output_time_step_s": output_time_step_s,
        "liquid_density_kg_m3": liquid_density_kg_m3,
        "saturation_temperature_k": saturation_temperature_k,
        "latent_heat_j_kg": latent_heat_j_kg,
    }
    for name, value in positive.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(liquid_kinematic_viscosity_m2_s) \
            or liquid_kinematic_viscosity_m2_s < 0.0:
        raise ValueError(
            "liquid_kinematic_viscosity_m2_s must be finite and non-negative"
        )
    if inflow_duration_s > duration_s:
        raise ValueError("inflow_duration_s cannot exceed duration_s")
    if not isinstance(substrate, (SolidSubstrate, ConstantHeatFluxSurface)):
        raise TypeError(
            "substrate must be SolidSubstrate or ConstantHeatFluxSurface"
        )
    if evaporation_momentum_closure not in {
        "zero_radial_momentum_vapor", "liquid_velocity_carryoff",
    }:
        raise ValueError("unsupported evaporation_momentum_closure")
    if not math.isfinite(solid_heat_flux_multiplier) \
            or solid_heat_flux_multiplier <= 0.0:
        raise ValueError("solid_heat_flux_multiplier must be finite and positive")
    if solid_heat_flux_cap_w_m2 is not None and (
        not math.isfinite(solid_heat_flux_cap_w_m2)
        or solid_heat_flux_cap_w_m2 <= 0.0
    ):
        raise ValueError("solid_heat_flux_cap_w_m2 must be positive or None")
    if isinstance(substrate, ConstantHeatFluxSurface) and (
        solid_heat_flux_multiplier != 1.0
        or solid_heat_flux_cap_w_m2 is not None
    ):
        raise ValueError(
            "solid heat-flux controls do not apply to ConstantHeatFluxSurface"
        )

    centres, faces, areas = _radial_grid(numerics)
    source_shape = _source_shape(centres, areas, numerics)
    depth = np.zeros_like(centres)
    momentum = np.zeros_like(centres)
    contact_age = np.zeros_like(centres)
    cumulative_inflow = 0.0
    cumulative_evaporation = 0.0
    escaped_mass = 0.0
    numerical_adjustment_mass = 0.0
    time_s = 0.0
    last_output_evaporation = 0.0
    last_output_time = 0.0
    steps: list[DynamicPoolStep] = []

    output_times = [0.0]
    while output_times[-1] < duration_s - 1.0e-12:
        output_times.append(min(duration_s, output_times[-1] + output_time_step_s))

    if isinstance(substrate, SolidSubstrate):
        heat_coefficient = (
            substrate.conductivity_w_m_k
            * max(0.0, substrate.initial_temperature_k - saturation_temperature_k)
            * solid_heat_flux_multiplier
            / (
                liquid_density_kg_m3 * latent_heat_j_kg
                * math.sqrt(math.pi * substrate.diffusivity_m2_s)
            )
        )
        heat_flux_cap_depth_rate = (
            solid_heat_flux_cap_w_m2
            / (liquid_density_kg_m3 * latent_heat_j_kg)
            if solid_heat_flux_cap_w_m2 is not None else None
        )
    else:
        heat_coefficient = 0.0
        heat_flux_cap_depth_rate = None

    for target_time in output_times:
        while time_s < target_time - 1.0e-12:
            dt = min(_stable_step(depth, momentum, numerics), target_time - time_s)
            if time_s < inflow_duration_s < time_s + dt:
                dt = inflow_duration_s - time_s
            inflow_rate = (
                liquid_inflow_rate_kg_s if time_s < inflow_duration_s else 0.0
            )
            depth, momentum, escaped_volume, numerical_volume = _flow_step(
                depth, momentum, dt_s=dt, centres=centres, faces=faces,
                areas=areas,
                liquid_kinematic_viscosity_m2_s=(
                    liquid_kinematic_viscosity_m2_s
                ),
                numerics=numerics,
            )
            numerical_adjustment_mass += numerical_volume * liquid_density_kg_m3
            depth += (
                inflow_rate / liquid_density_kg_m3 * source_shape * dt
            )

            wet = depth > numerics.dry_depth_m
            old_age = contact_age.copy()
            contact_age[wet] += dt
            evaporation_capacity = np.zeros_like(depth)
            if isinstance(substrate, ConstantHeatFluxSurface):
                evaporation_capacity[wet] = (
                    substrate.heat_flux_w_m2 * dt
                    / (liquid_density_kg_m3 * latent_heat_j_kg)
                )
            elif heat_flux_cap_depth_rate is None:
                evaporation_capacity[wet] = 2.0 * heat_coefficient * (
                    np.sqrt(contact_age[wet]) - np.sqrt(old_age[wet])
                )
            else:
                switch_age = (
                    heat_coefficient / heat_flux_cap_depth_rate
                ) ** 2
                capped_duration = np.maximum(
                    0.0,
                    np.minimum(contact_age, switch_age)
                    - np.minimum(old_age, switch_age),
                )
                conduction_start = np.maximum(old_age, switch_age)
                conduction = np.where(
                    contact_age > conduction_start,
                    2.0 * heat_coefficient * (
                        np.sqrt(contact_age) - np.sqrt(conduction_start)
                    ),
                    0.0,
                )
                evaporation_capacity[wet] = (
                    heat_flux_cap_depth_rate * capped_duration[wet]
                    + conduction[wet]
                )
            before_evaporation = depth.copy()
            evaporated_depth = np.minimum(evaporation_capacity, depth)
            depth -= evaporated_depth
            if evaporation_momentum_closure == "liquid_velocity_carryoff":
                remaining_fraction = np.divide(
                    depth, before_evaporation, out=np.zeros_like(depth),
                    where=before_evaporation > 0.0,
                )
                momentum *= remaining_fraction
            momentum[depth <= numerics.dry_depth_m] = 0.0

            evaporated_mass = float(
                np.dot(evaporated_depth, areas) * liquid_density_kg_m3
            )
            cumulative_inflow += inflow_rate * dt
            cumulative_evaporation += evaporated_mass
            escaped_mass += escaped_volume * liquid_density_kg_m3
            time_s += dt

        liquid_mass = float(np.dot(depth, areas) * liquid_density_kg_m3)
        wet = depth > 0.0
        wet_area = float(np.sum(areas[wet]))
        detected = np.flatnonzero(depth >= numerics.reported_front_depth_m)
        radius = float(faces[detected[-1] + 1]) if len(detected) else 0.0
        interval = target_time - last_output_time
        evaporation_rate = (
            (cumulative_evaporation - last_output_evaporation) / interval
            if interval > 0.0 else 0.0
        )
        residual = (
            cumulative_inflow + numerical_adjustment_mass
            - cumulative_evaporation - escaped_mass - liquid_mass
        )
        steps.append(DynamicPoolStep(
            elapsed_s=float(target_time),
            liquid_inflow_rate_kg_s=float(
                liquid_inflow_rate_kg_s
                if target_time < inflow_duration_s else 0.0
            ),
            cumulative_liquid_inflow_kg=float(cumulative_inflow),
            liquid_mass_kg=float(liquid_mass),
            reported_radius_m=radius,
            reported_area_m2=float(math.pi * radius**2),
            wet_area_m2=wet_area,
            wet_area_mean_depth_m=(
                liquid_mass / (liquid_density_kg_m3 * wet_area)
                if wet_area > 0.0 else 0.0
            ),
            evaporation_rate_kg_s=float(evaporation_rate),
            cumulative_evaporated_mass_kg=float(cumulative_evaporation),
            escaped_domain_liquid_mass_kg=float(escaped_mass),
            cumulative_numerical_mass_adjustment_kg=float(
                numerical_adjustment_mass
            ),
            mass_residual_kg=float(residual),
        ))
        last_output_evaporation = cumulative_evaporation
        last_output_time = target_time

    return DynamicPoolResult(
        steps=tuple(steps),
        radial_centres_m=tuple(float(value) for value in centres),
        final_depth_m=tuple(float(value) for value in depth),
        saturation_temperature_k=float(saturation_temperature_k),
        liquid_density_kg_m3=float(liquid_density_kg_m3),
        latent_heat_j_kg=float(latent_heat_j_kg),
        substrate=substrate,
        numerics=numerics,
        evaporation_momentum_closure=evaporation_momentum_closure,
        solid_heat_flux_multiplier=float(solid_heat_flux_multiplier),
        solid_heat_flux_cap_w_m2=solid_heat_flux_cap_w_m2,
    )


__all__ = [
    "ConstantHeatFluxSurface", "DynamicPoolNumerics", "DynamicPoolStep",
    "DynamicPoolResult",
    "simulate_axisymmetric_spreading_pool",
]
