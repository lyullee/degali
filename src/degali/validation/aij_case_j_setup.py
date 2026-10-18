"""Create a physical-property Case J configuration for a resolved FV run.

This is a coarse development setup, not a validated/adequately resolved
benchmark: the caller chooses grid spacing and simulated duration, and the
5-mm source is subgrid whenever the selected spacing exceeds its diameter.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.resolved_obstacle_flow_3d import ResolvedObstacleFlow3DConfig
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from .aij_case_j import AijCaseJDataset
from .aij_case_j_comparison import ETHYLENE_MOLAR_MASS_KG_MOL


R_UNIVERSAL_J_MOL_K = 8.31446261815324
AIR_MOLAR_MASS_KG_MOL = 28.96546e-3


@dataclass(frozen=True)
class AijCaseJResolvedSetup:
    config: ResolvedObstacleFlow3DConfig
    coordinate_translation_model_minus_experiment_m: tuple[float, float, float]
    grid_spacing_m: float
    source_ethylene_mass_rate_kg_s: float
    source_total_to_ethylene_mass_ratio: float
    physical_air_kinematic_viscosity_m2_s: float
    reynolds_from_properties: float
    reference_reynolds_number: float
    source_hole_resolved_cells_across: float
    scope_warnings: tuple[str, ...]


def _air_kinematic_viscosity_sutherland(temperature_k: float, density: float) -> float:
    reference_temperature_k = 273.15
    sutherland_k = 111.0
    dynamic_viscosity = (
        1.716e-5
        * (temperature_k / reference_temperature_k) ** 1.5
        * (reference_temperature_k + sutherland_k)
        / (temperature_k + sutherland_k)
    )
    return dynamic_viscosity / density


def build_aij_case_j_resolved_setup(
    dataset: AijCaseJDataset,
    *,
    grid_spacing_m: float = 0.02,
    duration_s: float = 0.60,
    time_step_s: float = 5.0e-4,
) -> AijCaseJResolvedSetup:
    """Map published Case J conditions to a Cartesian 3-D binary-mixture FV run."""
    dx = float(grid_spacing_m)
    duration = float(duration_s)
    dt = float(time_step_s)
    if not all(math.isfinite(value) and value > 0.0 for value in (dx, duration, dt)):
        raise ValueError("grid spacing, duration, and time step must be finite and positive")
    if duration < 2.0 * dt:
        raise ValueError("duration must contain at least two time steps")

    # Preserve the published experiment coordinates with fixed upstream,
    # lateral, and downstream margins. x increases with the approaching wind.
    exp_x_min, exp_x_max = -0.40, 0.48
    exp_y_min, exp_y_max = -0.20, 0.20
    exp_z_min, exp_z_max = 0.0, 0.40
    lengths = (exp_x_max - exp_x_min, exp_y_max - exp_y_min, exp_z_max)
    ratios = tuple(length / dx for length in lengths)
    cells = tuple(int(round(value)) for value in ratios)
    if any(not math.isclose(ratio, count, rel_tol=0.0, abs_tol=1.0e-9)
           for ratio, count in zip(ratios, cells)):
        raise ValueError("grid spacing must exactly divide the mapped Case J domain")

    translation = (-exp_x_min, -exp_y_min, -exp_z_min)
    block_x = (translation[0] - dataset.building_width_m, translation[0])
    block_y = (
        translation[1] - 0.5 * dataset.building_depth_m,
        translation[1] + 0.5 * dataset.building_depth_m,
    )
    obstacle = AxisAlignedCuboid(
        block_x[0], block_x[1], block_y[0], block_y[1],
        0.0, dataset.building_height_m, "AIJ Case J 1:1:2 building",
    )

    ambient_temperature = dataset.approach_temperature_k
    ambient_pressure = 101_325.0
    closure = BinaryIdealGasMixtureClosure(
        species_name="ethylene",
        species_molar_mass_kg_mol=ETHYLENE_MOLAR_MASS_KG_MOL,
        species_heat_capacity_j_kg_k=1_500.0,
        ambient_temperature_k=ambient_temperature,
        ambient_pressure_pa=ambient_pressure,
        gravity_m_s2=9.80665,
    )
    air_density = closure.ambient_density_kg_m3
    air_nu = _air_kinematic_viscosity_sutherland(ambient_temperature, air_density)
    source_temperature = dataset.source_temperature_k
    source_mole_fraction = dataset.source_tracer_mole_fraction
    source_mass_fraction = (
        source_mole_fraction * ETHYLENE_MOLAR_MASS_KG_MOL
        / (
            source_mole_fraction * ETHYLENE_MOLAR_MASS_KG_MOL
            + (1.0 - source_mole_fraction) * AIR_MOLAR_MASS_KG_MOL
        )
    )
    source_ratio = 1.0 / source_mass_fraction
    source_molar_density = ambient_pressure / (R_UNIVERSAL_J_MOL_K * source_temperature)
    ethylene_mass_rate = (
        dataset.source_mixture_volume_flow_m3_s
        * source_mole_fraction
        * source_molar_density
        * ETHYLENE_MOLAR_MASS_KG_MOL
    )

    velocity_z, velocity_u = dataset.inlet_velocity_profile
    temperature_z, temperature_k = dataset.inlet_temperature_profile
    if temperature_z[0] > 0.0:
        temperature_z = (0.0,) + temperature_z
        temperature_k = (temperature_k[0],) + temperature_k
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(0.12, 0.5 * lengths[1]),
        flat_to_flat_m=max(4.0 * dx, 0.04),
        z_min_m=0.0,
        z_max_m=min(0.08, lengths[2]),
        species_mole_fraction=0.0,
        temperature_k=ambient_temperature,
        label="zero-species binary-state initialization marker",
    )
    source_x_experiment = (
        dataset.source_downstream_offset_building_heights
        * dataset.building_height_m
    )
    source_position = (
        source_x_experiment + translation[0],
        translation[1],
        0.5 * dataset.source_hole_diameter_m,
    )
    source_sigma = max(0.5 * dx, 0.5 * dataset.source_hole_diameter_m)
    source_schedule = SourceRateSchedule(
        (0.0, duration - dt, duration),
        (ethylene_mass_rate, ethylene_mass_rate, 0.0),
        source_id="AIJ-Case-J-5mol-percent-ethylene-floor-source",
        rate_operator="piecewise_constant",
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=lengths,
        cells=cells,
        wind_speed_m_s=dataset.reference_speed_m_s,
        kinematic_viscosity_m2_s=air_nu,
        flow_steps=1,
        minimum_flow_steps=1,
        convergence_check_interval=1,
        lattice_inlet_speed=0.03,
        duration_s=duration,
        scalar_time_step_s=dt,
        scalar_diffusivity_m2_s=1.3e-5,
        thermal_diffusivity_m2_s=2.0e-5,
        source_schedule=source_schedule,
        source_total_mass_to_species_mass_ratio=source_ratio,
        source_position_m=source_position,
        source_sigma_m=source_sigma,
        obstacles=(obstacle,),
        inlet_profile_z_m=velocity_z,
        inlet_profile_u_m_s=velocity_u,
        inlet_temperature_profile_z_m=temperature_z,
        inlet_temperature_profile_k=temperature_k,
        obstacle_wall_temperature_k=dataset.mean_building_surface_temperature_k,
        ground_temperature_k=dataset.floor_temperature_k,
        smagorinsky_constant=0.0,
        wale_constant=0.5,
        subgrid_model="wale",
        ambient_temperature_k=ambient_temperature,
        ambient_pressure_pa=ambient_pressure,
        ambient_air_density_kg_m3=air_density,
        thermal_buoyancy_coupling=True,
        source_temperature_k=source_temperature,
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.9,
        scalar_turbulent_diffusivity_model="molecular_plus_subgrid_eddy_viscosity",
        thermal_turbulent_diffusivity_model="molecular_plus_subgrid_eddy_viscosity",
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud,
        binary_mixture_closure=closure,
        transient_average_start_s=0.5 * duration,
        thermodynamic_corrector_maximum_iterations=12,
        thermodynamic_corrector_minimum_iterations=2,
        require_thermodynamic_corrector_convergence=False,
        store_scalar_fields=True,
    )
    return AijCaseJResolvedSetup(
        config=config,
        coordinate_translation_model_minus_experiment_m=translation,
        grid_spacing_m=dx,
        source_ethylene_mass_rate_kg_s=ethylene_mass_rate,
        source_total_to_ethylene_mass_ratio=source_ratio,
        physical_air_kinematic_viscosity_m2_s=air_nu,
        reynolds_from_properties=(
            dataset.reference_speed_m_s * dataset.building_height_m / air_nu
        ),
        reference_reynolds_number=dataset.reference_reynolds_number,
        source_hole_resolved_cells_across=dataset.source_hole_diameter_m / dx,
        scope_warnings=(
            "coarse development mesh: building width has only {:.3g} cells".format(
                dataset.building_width_m / dx
            ),
            "5-mm source is subgrid; source spreading uses a grid-dependent Gaussian width",
            "measured building wall-temperature field is replaced by its arithmetic mean",
            "AIJ Case J reports time-averaged measurements; this finite transient is not temporal validation",
            "the species uses a calorically-perfect ethylene/air EOS with constant heat capacities",
            "the zero-species initial-cloud marker carries no gas; the physical source is the continuing schedule",
            "WALE, turbulent Schmidt 0.7, and turbulent Prandtl 0.9 are candidate closures, not Case J-calibrated values",
        ),
    )


__all__ = ["AijCaseJResolvedSetup", "build_aij_case_j_resolved_setup"]
