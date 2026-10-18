"""Run the coordinate-consistent TI021 streaming handoff preflight.

This audit exercises the actual 200 x 200 x 40 m resolved-inner domain from
the frozen full-transient plan.  It samples the physical site-y=100 m plane
(solver x=160 m), writes the five conservative flux histories directly to a
durable HDF5 file, reopens that file, and remaps it to the 640 m wide nested
outer inlet.  Concentration observations are deliberately not accessed.

The short horizon proves execution, coordinates, persistence and conservative
remapping only.  The initial puff cannot reach the handoff plane in 0.25 s, so
this audit makes no plume-arrival or dispersion-accuracy claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle
import numpy as np

from degali.addons.atmospheric_surface_layer_mixing import (
    NeutralSurfaceLayerMixing,
)
from degali.addons.binary_mixture_handoff_plane_3d import (
    BinaryMixtureHandoffPlane3DRequest,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.conservative_handoff_hdf5_3d import (
    ConservativeHandoffHdf5Writer3D,
    load_conservative_handoff_hdf5_3d,
)
from degali.addons.conservative_nested_handoff_3d import (
    remap_conservative_nested_handoff_result_3d,
)
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedInletVelocityHistory3D,
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import CircularArcWall


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-handoff-streaming-preflight.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _log_profile(height_m: float) -> float:
    if height_m <= 0.0:
        return 0.0
    return 0.205 / 0.4 * math.log((height_m + 0.005) / 0.005)


def build_configuration(
    properties: dict[str, object], *,
    cell_size_m: float = 5.0,
    time_step_s: float = 0.125,
    active_steps: int = 2,
    corrector_tolerance: float = 1.0e-6,
    circular_arc_representation: str = "finite_thickness_voxel",
    obstacle_height_m: float = 5.0,
    inner_domain_m: tuple[float, float, float] = (200.0, 200.0, 40.0),
    handoff_site_downwind_y_m: float = 100.0,
    initial_cloud_flat_to_flat_m: float = 14.0,
    initial_cloud_z_max_m: float = 13.0,
    finite_source_release_duration_s: float | None = None,
    finite_source_start_delay_s: float = 0.0,
    finite_source_mass_kg: float | None = None,
    finite_source_mass_basis: str = "species",
    finite_source_position_z_m: float = 1.0,
    finite_source_momentum_velocity_m_s: tuple[float, float, float] = (0.0, 0.0, 0.0),
    finite_source_total_mass_to_species_mass_ratio: float = 1.0,
    wind_speed_m_s: float = 3.9,
    wind_heading_relative_deg: float = 6.1,
    receptor_vertical_sampling: str = "cell_center",
    finite_source_release_shape: str = "constant",
    finite_source_sigma_xyz_m: tuple[float, float, float] | None = None,
    finite_source_kernel: str = "gaussian",
    scalar_advection_scheme: str = "muscl_van_leer",
    source_time_integration: str = "prestep",
) -> tuple[
    ResolvedObstacleFlow3DConfig,
    BinaryIdealGasMixtureClosure,
    dict[str, object],
]:
    """Build the observation-blind report-coordinate preflight."""

    dx = float(cell_size_m)
    dt = float(time_step_s)
    tolerance = float(corrector_tolerance)
    if dx <= 0.0 or dt <= 0.0 or active_steps < 1:
        raise ValueError("cell size, time step and active_steps must be positive")
    if not math.isfinite(tolerance) or not 0.0 < tolerance <= 1.0e-3:
        raise ValueError("corrector_tolerance must lie in (0, 1e-3]")
    domain = tuple(float(value) for value in inner_domain_m)
    handoff_site_y = float(handoff_site_downwind_y_m)
    handoff_x = 60.0 + handoff_site_y
    if len(domain) != 3 or any(not math.isfinite(value) or value <= 0.0 for value in domain):
        raise ValueError("inner_domain_m must contain three positive finite lengths")
    if not math.isfinite(handoff_site_y) or handoff_site_y < 0.0:
        raise ValueError("handoff site downwind coordinate must be finite and non-negative")
    counts = tuple(length / dx for length in domain)
    if any(not math.isclose(value, round(value), abs_tol=1.0e-12) for value in counts):
        raise ValueError("cell size must divide every inner-domain dimension")
    cells = tuple(int(round(value)) for value in counts)
    if handoff_x >= domain[0] or not math.isclose(handoff_x / dx, round(handoff_x / dx), abs_tol=1.0e-12):
        raise ValueError("handoff must lie on an interior x-face of the resolved domain")
    inner_crosswind_origin = 0.5 * domain[1]

    if properties.get("status") != "accepted_molecular_property_component":
        raise ValueError("R12 property audit is not accepted")
    policy = properties["solver_policy"]
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=float(policy["species_molar_mass_kg_mol"]),
        species_heat_capacity_j_kg_k=float(
            policy["species_heat_capacity_j_kg_k"]
        ),
        ambient_temperature_k=291.15,
        ambient_pressure_pa=100_500.0,
        minimum_temperature_k=150.0,
    )
    duration = active_steps * dt
    release_duration = (
        None if finite_source_release_duration_s is None
        else float(finite_source_release_duration_s)
    )
    start_delay = float(finite_source_start_delay_s)
    if not math.isfinite(start_delay) or start_delay < 0.0:
        raise ValueError("finite_source_start_delay_s must be finite and non-negative")
    if release_duration is None and start_delay > 0.0:
        raise ValueError(
            "finite_source_start_delay_s requires finite_source_release_duration_s"
        )
    if release_duration is not None and start_delay + release_duration > duration + 1.0e-12:
        raise ValueError(
            "finite source start delay plus release duration cannot exceed active duration"
        )
    if release_duration is not None and (
        not math.isfinite(release_duration)
        or release_duration <= 0.0
        or release_duration > duration + 1.0e-12
    ):
        raise ValueError(
            "finite_source_release_duration_s must lie in (0, active duration]"
        )
    if finite_source_release_shape not in {
        "constant", "front_loaded", "back_loaded", "smooth_ramp",
    }:
        raise ValueError(
            "finite_source_release_shape must be constant, front_loaded, "
            "back_loaded or smooth_ramp"
        )
    if finite_source_sigma_xyz_m is not None:
        if len(finite_source_sigma_xyz_m) != 3 or any(
            not math.isfinite(float(value)) or float(value) <= 0.0
            for value in finite_source_sigma_xyz_m
        ):
            raise ValueError("finite_source_sigma_xyz_m must contain three positive values")
    wind_speed = float(wind_speed_m_s)
    wind_heading = float(wind_heading_relative_deg)
    if not math.isfinite(wind_speed) or wind_speed <= 0.0:
        raise ValueError("wind_speed_m_s must be positive and finite")
    if not math.isfinite(wind_heading):
        raise ValueError("wind_heading_relative_deg must be finite")
    heading = math.radians(wind_heading)
    streamwise = wind_speed * math.cos(heading)
    crosswind = wind_speed * math.sin(heading)
    profile_z = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0)
    profile_u = tuple(_log_profile(value) for value in profile_z)
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(60.0, inner_crosswind_origin),
        flat_to_flat_m=float(initial_cloud_flat_to_flat_m),
        z_min_m=0.0,
        z_max_m=float(initial_cloud_z_max_m),
        species_mole_fraction=0.0 if release_duration is not None else 0.32,
        temperature_k=291.15,
        label="TI021 report-scale initial R12/air cloud",
    )
    obstacle_height = float(obstacle_height_m)
    if not math.isfinite(obstacle_height) or obstacle_height <= 0.0:
        raise ValueError("obstacle_height_m must be positive and finite")
    wall = CircularArcWall(
        center_m=(60.0, inner_crosswind_origin),
        radius_m=50.0,
        thickness_m=(
            0.0 if circular_arc_representation == "zero_thickness_surface"
            else dx
        ),
        z_min_m=0.0,
        z_max_m=obstacle_height,
        arc_start_deg=-90.0,
        arc_end_deg=90.0,
        label="TI021 report-priority 50 m radius, 5 m fence",
    )
    mixing = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=0.205,
        roughness_length_m=0.005,
        source_id="TI021 concentration-blind mixing policy v1",
    )
    release_mass_kg = float(
        closure.state_from_mole_fraction_and_temperature(
            0.32, 291.15,
        ).species_mass_concentration_kg_m3 * cloud.exact_volume_m3
    )
    source_mass_kg = (
        release_mass_kg if finite_source_mass_kg is None
        else float(finite_source_mass_kg)
    )
    if finite_source_mass_basis not in {"species", "total"}:
        raise ValueError("finite_source_mass_basis must be 'species' or 'total'")
    if finite_source_mass_basis == "total" and finite_source_mass_kg is None:
        raise ValueError("finite_source_mass_basis='total' requires finite_source_mass_kg")
    if finite_source_mass_kg is not None and release_duration is None:
        raise ValueError(
            "finite_source_mass_kg requires finite_source_release_duration_s"
        )
    source_z_m = float(finite_source_position_z_m)
    if not math.isfinite(source_z_m) or source_z_m < 0.0 or source_z_m > float(inner_domain_m[2]):
        raise ValueError("finite_source_position_z_m must lie within the domain height")
    source_momentum = tuple(float(value) for value in finite_source_momentum_velocity_m_s)
    if len(source_momentum) != 3 or any(not math.isfinite(value) for value in source_momentum):
        raise ValueError("finite_source_momentum_velocity_m_s must contain three finite values")
    if any(abs(value) > 0.0 for value in source_momentum) and release_duration is None:
        raise ValueError(
            "finite_source_momentum_velocity_m_s requires a finite source release"
        )
    source_mass_ratio = float(finite_source_total_mass_to_species_mass_ratio)
    if not math.isfinite(source_mass_ratio) or source_mass_ratio < 1.0:
        raise ValueError(
            "finite_source_total_mass_to_species_mass_ratio must be finite and at least one"
        )
    if not math.isfinite(source_mass_kg) or source_mass_kg <= 0.0:
        raise ValueError("finite_source_mass_kg must be positive and finite")
    # The resolved solver's source schedule is deliberately expressed in the
    # tracked-species mass.  Make a total-mixture CLI input explicit and
    # convert it once here; otherwise a caller can silently inject ratio times
    # too much material when using the workbook's 4,962 kg total-mixture value.
    requested_source_mass_kg = source_mass_kg
    if finite_source_mass_basis == "total":
        source_mass_kg /= source_mass_ratio
    if release_duration is None:
        source_schedule = SourceRateSchedule(
            (0.0, duration), (0.0, 0.0),
            source_id="TI021-initial-cloud-only",
        )
    else:
        release_end = start_delay + release_duration
        if finite_source_release_shape == "constant":
            local_times = (0.0, release_duration)
            local_rates = (source_mass_kg / release_duration, 0.0)
            source_id = "TI021-finite-release-diagnostic"
            rate_operator = "piecewise_constant"
        elif finite_source_release_shape == "front_loaded":
            local_times = (0.0, 0.5 * release_duration, release_duration)
            local_rates = (2.0 * source_mass_kg / release_duration, 0.0, 0.0)
            source_id = "TI021-finite-release-front-loaded-diagnostic"
            rate_operator = "piecewise_constant"
        elif finite_source_release_shape == "smooth_ramp":
            # A C0 temporal pulse with zero source rate at both endpoints.
            # The four linear segments integrate exactly to source_mass_kg:
            # [0,1,2,1,0] * M/T over quarter-duration intervals.
            local_times = (
                0.0,
                0.25 * release_duration,
                0.5 * release_duration,
                0.75 * release_duration,
                release_duration,
            )
            scale = source_mass_kg / release_duration
            local_rates = (0.0, scale, 2.0 * scale, scale, 0.0)
            source_id = "TI021-finite-release-smooth-ramp-diagnostic"
            rate_operator = "linear"
        else:
            local_times = (0.0, 0.5 * release_duration, release_duration)
            local_rates = (0.0, 4.0 * source_mass_kg / release_duration, 0.0)
            source_id = "TI021-finite-release-back-loaded-diagnostic"
            rate_operator = "linear"
        shifted_times = tuple(start_delay + value for value in local_times)
        shifted_rates = local_rates
        if start_delay > 0.0:
            schedule_times = (0.0,) + shifted_times
            schedule_rates = (0.0,) + shifted_rates
        else:
            schedule_times = shifted_times
            schedule_rates = shifted_rates
        if release_end < duration - 1.0e-12:
            schedule_times = schedule_times + (duration,)
            schedule_rates = schedule_rates + (0.0,)
        source_schedule = SourceRateSchedule(
            schedule_times, schedule_rates,
            source_id=source_id,
            rate_operator=rate_operator,
        )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=domain,
        cells=cells,
        wind_speed_m_s=wind_speed,
        kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=1,
        minimum_flow_steps=1,
        convergence_check_interval=1,
        lattice_inlet_speed=0.03,
        duration_s=duration,
        scalar_time_step_s=dt,
        scalar_diffusivity_m2_s=float(
            policy["molecular_species_diffusivity_m2_s"]
        ),
        thermal_diffusivity_m2_s=float(
            properties["pure_r12_ambient_properties"][
                "thermal_diffusivity_m2_s"
            ]
        ),
        source_schedule=source_schedule,
        source_position_m=(60.0, inner_crosswind_origin, source_z_m),
        source_momentum_velocity_m_s=source_momentum,
        source_total_mass_to_species_mass_ratio=source_mass_ratio,
        source_sigma_m=dx,
        source_sigma_xyz_m=(
            tuple(float(value) for value in finite_source_sigma_xyz_m)
            if release_duration is not None and finite_source_sigma_xyz_m is not None
            else None
        ),
        source_kernel=finite_source_kernel,
        source_temperature_k=291.15,
        receptor_vertical_sampling=receptor_vertical_sampling,
        obstacles=(wall,),
        circular_arc_representation=circular_arc_representation,
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
        active_inlet_velocity_history=ResolvedInletVelocityHistory3D(
            time_s=(0.0, duration),
            u_m_s=(streamwise, streamwise),
            v_m_s=(crosswind, crosswind),
            w_m_s=(0.0, 0.0),
            history_id=(
                f"TI021-explicit-wind-{wind_speed:g}mps-"
                f"{wind_heading:g}deg-relative"
            ),
        ),
        ambient_temperature_k=closure.ambient_temperature_k,
        ambient_pressure_pa=closure.ambient_pressure_pa,
        ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
        gravity_m_s2=closure.gravity_m_s2,
        dry_air_heat_capacity_j_kg_k=closure.carrier_heat_capacity_j_kg_k,
        hydrogen_heat_capacity_j_kg_k=closure.species_heat_capacity_j_kg_k,
        thermal_buoyancy_coupling=True,
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud,
        binary_mixture_closure=closure,
        initial_binary_cloud_momentum="quiescent",
        thermodynamic_corrector_maximum_iterations=32,
        thermodynamic_corrector_minimum_iterations=2,
        thermodynamic_corrector_density_tolerance=tolerance,
        thermodynamic_corrector_velocity_tolerance=tolerance,
        require_thermodynamic_corrector_convergence=True,
        variable_density_projection_linear_solver="amg_cg",
        subgrid_model="none",
        smagorinsky_constant=0.0,
        background_surface_layer_mixing=mixing,
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.74,
        scalar_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        thermal_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        scalar_advection_scheme=scalar_advection_scheme,
        source_time_integration=source_time_integration,
        store_scalar_fields=True,
    )
    coordinates = {
        "site_frame": "+y downwind; +x site-right",
        "solver_frame": "+x downwind; +y=-site_x; +z vertical",
        # Keep the coordinate audit aligned with the actual resolved source
        # position.  The old hard-coded z=0 entry silently disagreed with the
        # finite-source z input (default 1 m).
        "source_solver_xyz_m": [60.0, inner_crosswind_origin, source_z_m],
        "inner_crosswind_origin_m": inner_crosswind_origin,
        "handoff_site_downwind_y_m": handoff_site_y,
        "handoff_inner_solver_x_m": handoff_x,
        "source_transverse_q_equals_minus_site_x_bounds_m": [
            -inner_crosswind_origin, inner_crosswind_origin,
        ],
        "nested_outer_transverse_q_equals_minus_site_x_bounds_m": [
            -320.0, 320.0,
        ],
        "nested_outer_x_at_handoff_m": 0.0,
        "finite_source_mass_input_kg": (
            None if finite_source_mass_kg is None else requested_source_mass_kg
        ),
        "finite_source_mass_basis": finite_source_mass_basis,
        "finite_source_schedule_species_mass_kg": (
            None if finite_source_release_duration_s is None else source_mass_kg
        ),
        "finite_source_schedule_total_mass_kg": (
            None if finite_source_release_duration_s is None
            else source_mass_kg * source_mass_ratio
        ),
    }
    return config, closure, coordinates


def run_audit(
    output_directory: Path, *,
    property_audit: Path,
    replace: bool = False,
) -> tuple[dict[str, object], dict[str, Path]]:
    """Execute, reopen and conservatively remap the streaming checkpoint."""

    output = output_directory.resolve()
    property_path = property_audit.resolve()
    properties = json.loads(property_path.read_text(encoding="utf-8"))
    config, _closure, coordinates = build_configuration(properties)
    dx = config.domain_m[0] / config.cells[0]
    source_z_edges = np.linspace(0.0, config.domain_m[2], config.cells[2] + 1)
    source_q_edges = np.linspace(
        -0.5 * config.domain_m[1],
        0.5 * config.domain_m[1],
        config.cells[1] + 1,
    )
    target_z_edges = source_z_edges.copy()
    target_q_edges = np.arange(-320.0, 320.0 + dx, dx)
    output.mkdir(parents=True, exist_ok=True)
    hdf5_path = output / "ti021-inner-to-outer-handoff.h5"
    json_path = output / "ti021-handoff-streaming-preflight.json"
    plot_path = output / "ti021-handoff-streaming-preflight.png"
    if replace:
        for path in (hdf5_path, json_path, plot_path):
            if path.is_file():
                path.unlink()
    metadata = {
        "case": "Thorney Island trial 21",
        "species": "R12",
        "inner_domain_m": list(config.domain_m),
        "inner_cells": list(config.cells),
        "cell_size_m": dx,
        "handoff_inner_solver_x_m": float(coordinates["handoff_inner_solver_x_m"]),
        "transverse_coordinate": (
            f"q=-site_x=solver_y-{0.5 * config.domain_m[1]:g}m"
        ),
        "concentration_observation_fields_accessed": [],
    }
    request = BinaryMixtureHandoffPlane3DRequest(
        x_face_m=float(coordinates["handoff_inner_solver_x_m"]),
        sample_interval_s=config.scalar_time_step_s,
        include_initial=True,
    )
    started = time.perf_counter()
    with ConservativeHandoffHdf5Writer3D(
        hdf5_path,
        source_z_edges_m=source_z_edges,
        source_y_edges_m=source_q_edges,
        source_id=(
            "TI021 resolved inner streaming preflight "
            f"x={coordinates['handoff_inner_solver_x_m']:g}m"
        ),
        metadata=metadata,
    ) as writer:
        result = solve_resolved_obstacle_flow_3d(
            config,
            handoff_request=request,
            handoff_observer=writer,
        )
        writer_frame_count = writer.frame_count
    execution_runtime = time.perf_counter() - started
    loaded = load_conservative_handoff_hdf5_3d(hdf5_path)
    remapped = remap_conservative_nested_handoff_result_3d(
        loaded,
        target_z_edges_m=target_z_edges,
        target_y_edges_m=target_q_edges,
        source_id="TI021 x=160m HDF5 history remapped to 640m outer inlet",
    )
    diagnostics = result.diagnostics
    expected_times = np.arange(
        0.0,
        config.duration_s + 0.5 * config.scalar_time_step_s,
        config.scalar_time_step_s,
    )
    maximum_species_flux = float(np.max(np.abs(
        loaded.species_mass_flux_kg_m2_s
    )))
    maximum_remap_residual = float(
        remapped.diagnostics["maximum_relative_integral_residual"]
    )
    gates = {
        "actual_200x200x40m_inner_domain_executed": (
            tuple(config.domain_m) == (200.0, 200.0, 40.0)
            and tuple(config.cells) == (40, 40, 8)
        ),
        "physical_site_y100m_equals_inner_x160m": (
            math.isclose(
                coordinates["source_solver_xyz_m"][0]
                + coordinates["handoff_site_downwind_y_m"],
                coordinates["handoff_inner_solver_x_m"],
            )
        ),
        "transverse_coordinate_orientation_is_shared": (
            coordinates[
                "source_transverse_q_equals_minus_site_x_bounds_m"
            ] == [-100.0, 100.0]
            and coordinates[
                "nested_outer_transverse_q_equals_minus_site_x_bounds_m"
            ] == [-320.0, 320.0]
        ),
        "durable_hdf5_reopened_with_expected_frames": (
            writer_frame_count == expected_times.size
            and loaded.time_s.tolist() == expected_times.tolist()
            and loaded.diagnostics["frame_count"] == expected_times.size
        ),
        "source_plane_shape_matches_inner_grid": (
            loaded.species_mass_flux_kg_m2_s.shape
            == (expected_times.size, 8, 40)
        ),
        "target_plane_shape_matches_outer_grid": (
            remapped.species_mass_flux_kg_m2_s.shape
            == (expected_times.size, 8, 128)
        ),
        "all_five_fluxes_conserved_by_area_remap": (
            maximum_remap_residual <= 1.0e-12
        ),
        "solver_observer_report_matches_checkpoint": (
            diagnostics["binary_handoff_recording_enabled"]
            and diagnostics["binary_handoff_x_face_m"] == 160.0
            and diagnostics["binary_handoff_frame_count"] == expected_times.size
        ),
        "species_mass_conservative": (
            diagnostics["maximum_relative_mass_residual"] <= 1.0e-6
        ),
        "sensible_enthalpy_conservative": (
            diagnostics["maximum_relative_energy_residual"] <= 1.0e-6
        ),
        "continuity_residual_at_most_1e_7": (
            diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ] <= 1.0e-7
        ),
        "thermodynamic_corrector_converged_every_step": (
            diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
        ),
        "no_concentration_observations_accessed": (
            metadata["concentration_observation_fields_accessed"] == []
        ),
        "puff_not_falsely_claimed_at_handoff_in_short_horizon": (
            maximum_species_flux <= 1.0e-15
        ),
    }
    all_gates_pass = all(gates.values())
    record = {
        "schema": SCHEMA,
        "status": (
            "accepted_coordinate_consistent_streaming_handoff_preflight"
            if all_gates_pass
            else "withheld_streaming_handoff_preflight"
        ),
        "scope": (
            "short observation-blind execution proof for actual inner-domain "
            "coordinates, streaming persistence and conservative outer remap"
        ),
        "limitations": [
            "0.25 s is too short for the initial puff to reach the x=160 m handoff",
            "this audit does not validate plume arrival, dispersion, sensor dose or grid convergence",
            "the full 600 s no-fit inner and outer transients remain required",
        ],
        "property_audit": {
            "path": str(property_path),
            "sha256": _sha256(property_path),
        },
        "coordinates": coordinates,
        "configuration": {
            "domain_m": list(config.domain_m),
            "cells": list(config.cells),
            "cell_count": int(np.prod(config.cells)),
            "cell_size_m": dx,
            "duration_s": config.duration_s,
            "time_step_s": config.scalar_time_step_s,
            "active_steps": int(round(
                config.duration_s / config.scalar_time_step_s
            )),
            "handoff_sample_interval_s": request.sample_interval_s,
            "scalar_advection_scheme": config.scalar_advection_scheme,
            "effective_schmidt_number": config.effective_schmidt_number,
            "effective_prandtl_number": config.effective_prandtl_number,
            "observation_value_fields_accessed": [],
        },
        "checkpoint": {
            "path": str(hdf5_path),
            "sha256": _sha256(hdf5_path),
            "size_bytes": hdf5_path.stat().st_size,
            "frame_count": int(loaded.time_s.size),
            "time_s": loaded.time_s.tolist(),
            "source_plane_shape": list(
                loaded.species_mass_flux_kg_m2_s.shape
            ),
            "target_plane_shape": list(
                remapped.species_mass_flux_kg_m2_s.shape
            ),
        },
        "flux_diagnostics": {
            "maximum_absolute_species_flux_kg_m2_s": maximum_species_flux,
            "maximum_absolute_total_mass_flux_anomaly_kg_m2_s": float(
                np.max(np.abs(
                    loaded.total_mass_flux_anomaly_kg_m2_s
                ))
            ),
            "maximum_absolute_lateral_momentum_flux_anomaly_pa": float(
                np.max(np.abs(
                    loaded.lateral_momentum_flux_anomaly_pa
                ))
            ),
            "maximum_absolute_vertical_momentum_flux_anomaly_pa": float(
                np.max(np.abs(
                    loaded.vertical_momentum_flux_anomaly_pa
                ))
            ),
            "relative_integral_residuals": remapped.diagnostics[
                "relative_integral_residuals"
            ],
            "maximum_relative_integral_residual": maximum_remap_residual,
            "clipping_or_posthoc_renormalisation": False,
        },
        "solver_diagnostics": {
            "maximum_relative_species_mass_residual": diagnostics[
                "maximum_relative_mass_residual"
            ],
            "maximum_relative_sensible_enthalpy_residual": diagnostics[
                "maximum_relative_energy_residual"
            ],
            "maximum_continuity_residual_kg_m3_s": diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
            "maximum_global_mass_residual_kg_s": diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "maximum_momentum_balance_residual_kg_m_s": diagnostics[
                "maximum_momentum_balance_residual_kg_m_s"
            ],
            "corrector_nonconverged_steps": diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ],
            "corrector_iterations_per_step": diagnostics[
                "thermodynamic_corrector_iterations_per_step"
            ],
        },
        "runtime": {
            "audit_wall_clock_s": execution_runtime,
            "solver_reported_total_s": diagnostics["total_runtime_s"],
            "hdf5_size_bytes": hdf5_path.stat().st_size,
        },
        "gates": gates,
        "all_gates_pass": all_gates_pass,
        "next_required": (
            "add restartable full-state checkpoints, then run the 600 s "
            "5 m inner history and drive the 440 x 640 x 40 m outer domain"
        ),
    }
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    area = dx * dx
    integrated = {
        "mass anomaly (kg/s)": np.sum(
            loaded.total_mass_flux_anomaly_kg_m2_s, axis=(-2, -1)
        ) * area,
        "species (kg/s)": np.sum(
            loaded.species_mass_flux_kg_m2_s, axis=(-2, -1)
        ) * area,
        "lateral momentum (N)": np.sum(
            loaded.lateral_momentum_flux_anomaly_pa, axis=(-2, -1)
        ) * area,
        "vertical momentum (N)": np.sum(
            loaded.vertical_momentum_flux_anomaly_pa, axis=(-2, -1)
        ) * area,
    }
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
    axes[0].add_patch(Rectangle((0.0, -100.0), 200.0, 200.0,
                                facecolor="#f0a44b", alpha=0.18,
                                edgecolor="#b45f06", linewidth=1.5))
    axes[0].add_patch(Rectangle((160.0, -320.0), 440.0, 640.0,
                                facecolor="#5b9bd5", alpha=0.12,
                                edgecolor="#1f4e79", linewidth=1.5))
    axes[0].axvline(160.0, color="black", linestyle="--", linewidth=1.5)
    axes[0].scatter([60.0], [0.0], color="crimson", marker="*", s=90,
                    label="initial puff")
    axes[0].add_patch(Arc((60.0, 0.0), 100.0, 100.0,
                          theta1=-90.0, theta2=90.0,
                          color="black", linewidth=2.0))
    axes[0].set_xlim(0.0, 600.0)
    axes[0].set_ylim(-330.0, 330.0)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("physical downwind distance from inner origin (m)")
    axes[0].set_ylabel("q = -site x (m)")
    axes[0].set_title("coordinate-consistent nested domains")
    axes[0].legend(loc="upper right", fontsize=8)

    final_mass = loaded.total_mass_flux_anomaly_kg_m2_s[-1]
    image = axes[1].pcolormesh(
        source_q_edges, source_z_edges, final_mass,
        shading="flat", cmap="coolwarm",
    )
    axes[1].set_xlabel("q = -site x (m)")
    axes[1].set_ylabel("height z (m)")
    axes[1].set_title("final mass-flux anomaly at x=160 m")
    fig.colorbar(image, ax=axes[1], label="kg m$^{-2}$ s$^{-1}$")

    for label, values in integrated.items():
        axes[2].plot(loaded.time_s, values, marker="o", label=label)
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("integrated flux (quantity-specific SI unit)")
    axes[2].set_title("integrated streamed flux histories")
    axes[2].legend(fontsize=8)
    fig.suptitle(
        "TI021 actual-inner streaming handoff preflight — no concentration fit"
    )
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    return record, {
        "json": json_path,
        "hdf5": hdf5_path,
        "plot": plot_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--property-audit",
        type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
        ),
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-handoff-streaming-preflight-2026-10-11-v1"
        ),
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="replace only the three exact output files in the selected directory",
    )
    args = parser.parse_args()
    record, paths = run_audit(
        args.output_directory,
        property_audit=args.property_audit,
        replace=args.replace,
    )
    print(json.dumps({
        "status": record["status"],
        "all_gates_pass": record["all_gates_pass"],
        "paths": {key: str(path) for key, path in paths.items()},
        "runtime": record["runtime"],
        "maximum_relative_integral_residual": record[
            "flux_diagnostics"
        ]["maximum_relative_integral_residual"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
