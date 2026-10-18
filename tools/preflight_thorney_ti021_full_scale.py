"""Run a concentration-blind full-site TI021 coupled numerical preflight.

The horizon is intentionally one or a few small time steps.  This checks that
the report-scale initial cloud and 50 m curved wall can enter the same coupled
momentum/species/enthalpy/EOS/continuity solve before the expensive no-fit
transient and before any observation is scored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.atmospheric_surface_layer_mixing import (
    NeutralSurfaceLayerMixing,
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-size-m", type=float, default=2.0)
    parser.add_argument("--active-steps", type=int, default=1)
    parser.add_argument(
        "--time-step-s", type=float,
        help="finite-volume step; default is 0.001 times cell size",
    )
    parser.add_argument("--flow-steps", type=int, default=20)
    parser.add_argument(
        "--projection-solver",
        choices=("amg_cg", "cg", "cached_direct", "auto"),
        default="amg_cg",
    )
    parser.add_argument(
        "--surface-layer-mixing", action="store_true",
        help="apply the frozen neutral primary Sc_t=0.7, Pr_t=0.74 policy",
    )
    parser.add_argument("--effective-schmidt", type=float, default=0.7)
    parser.add_argument("--effective-prandtl", type=float, default=0.74)
    parser.add_argument(
        "--property-audit", type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
        ),
    )
    parser.add_argument("--output-directory", type=Path)
    args = parser.parse_args()
    dx = float(args.cell_size_m)
    if dx not in {1.0, 2.0}:
        raise ValueError("preflight cell size must be 1.0 or 2.0 m")
    if args.active_steps < 1 or args.flow_steps < 1:
        raise ValueError("step counts must be positive")
    output = args.output_directory or (
        ROOT / "outputs"
        / f"thorney-ti021-full-scale-preflight-dx{dx:g}m-2026-10-11-v1"
    )
    property_path = args.property_audit.resolve()
    properties = json.loads(property_path.read_text(encoding="utf-8"))
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
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(60.0, 60.0), flat_to_flat_m=14.0,
        z_min_m=0.0, z_max_m=13.0,
        species_mole_fraction=0.32, temperature_k=291.15,
        label="TI021 report-scale initial R12/air cloud",
    )
    wall = CircularArcWall(
        center_m=(60.0, 60.0), radius_m=50.0,
        thickness_m=dx, z_min_m=0.0, z_max_m=5.0,
        arc_start_deg=-90.0, arc_end_deg=90.0,
        label="TI021 report-priority 50 m radius, 5 m fence",
    )
    domain = (130.0, 140.0, 16.0)
    cells = tuple(int(round(value / dx)) for value in domain)
    if any(not math.isclose(count * dx, length) for count, length in zip(cells, domain)):
        raise ValueError("domain must contain an integer cell count")
    active_dt = (
        0.001 * dx if args.time_step_s is None else float(args.time_step_s)
    )
    if not math.isfinite(active_dt) or active_dt <= 0.0:
        raise ValueError("time step must be positive and finite")
    duration = args.active_steps * active_dt
    heading = math.radians(6.1)
    streamwise = 3.9 * math.cos(heading)
    crosswind = 3.9 * math.sin(heading)
    profile_z = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 16.0)
    profile_u = tuple(_log_profile(value) for value in profile_z)
    surface_layer_mixing = (
        NeutralSurfaceLayerMixing(
            friction_velocity_m_s=0.205,
            roughness_length_m=0.005,
            source_id="TI021 concentration-blind mixing policy v1",
        )
        if args.surface_layer_mixing else None
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=domain,
        cells=cells,
        wind_speed_m_s=3.9,
        kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=args.flow_steps,
        minimum_flow_steps=args.flow_steps,
        convergence_check_interval=max(1, min(5, args.flow_steps)),
        lattice_inlet_speed=0.03,
        duration_s=duration,
        scalar_time_step_s=active_dt,
        scalar_diffusivity_m2_s=float(
            policy["molecular_species_diffusivity_m2_s"]
        ),
        thermal_diffusivity_m2_s=float(
            properties["pure_r12_ambient_properties"][
                "thermal_diffusivity_m2_s"
            ]
        ),
        source_schedule=SourceRateSchedule(
            (0.0, duration), (0.0, 0.0),
            source_id="TI021-initial-cloud-only",
        ),
        source_position_m=(60.0, 60.0, 1.0),
        source_sigma_m=dx,
        obstacles=(wall,),
        inlet_profile_z_m=profile_z,
        inlet_profile_u_m_s=profile_u,
        active_inlet_velocity_history=ResolvedInletVelocityHistory3D(
            time_s=(0.0, duration),
            u_m_s=(streamwise, streamwise),
            v_m_s=(crosswind, crosswind),
            w_m_s=(0.0, 0.0),
            history_id="TI021-report-3p9mps-minus6p1deg",
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
        thermodynamic_corrector_maximum_iterations=32,
        thermodynamic_corrector_minimum_iterations=2,
        thermodynamic_corrector_density_tolerance=1.0e-6,
        thermodynamic_corrector_velocity_tolerance=1.0e-6,
        require_thermodynamic_corrector_convergence=True,
        variable_density_projection_linear_solver=args.projection_solver,
        subgrid_model="none",
        smagorinsky_constant=0.0,
        background_surface_layer_mixing=surface_layer_mixing,
        effective_schmidt_number=(
            args.effective_schmidt if args.surface_layer_mixing else None
        ),
        effective_prandtl_number=(
            args.effective_prandtl if args.surface_layer_mixing else None
        ),
        scalar_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
            if args.surface_layer_mixing
            else "legacy_total_kinematic_viscosity_ratio"
        ),
        thermal_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
            if args.surface_layer_mixing
            else "legacy_total_kinematic_viscosity_ratio"
        ),
        store_scalar_fields=True,
    )
    result = solve_resolved_obstacle_flow_3d(config)
    diagnostics = result.diagnostics
    voxel = diagnostics["obstacle_voxelization"][0]
    effective_height = (
        voxel["effective_voxel_bounds_m"][5]
        - voxel["effective_voxel_bounds_m"][4]
    )
    gates = {
        "report_scale_cloud_and_wall_share_domain": True,
        "physical_molecular_viscosity_retained": math.isclose(
            config.kinematic_viscosity_m2_s, 1.5e-5,
            rel_tol=0.0, abs_tol=1.0e-15,
        ),
        "no_d3q19_low_reynolds_spinup": (
            diagnostics["carrier_flow_initialization"] == "projected_inlet"
            and diagnostics["flow_steps_completed"] == 0
            and diagnostics["base_lattice_relaxation_time"] is None
        ),
        "initial_cloud_volume_exact": abs(
            diagnostics["initial_cloud_relative_geometry_volume_error"]
        ) <= 1.0e-11,
        "species_mass_conservative": (
            diagnostics["maximum_relative_mass_residual"] <= 1.0e-6
        ),
        "total_continuity_residual_at_most_1e_7": (
            diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ] <= 1.0e-7
        ),
        "corrector_converged_every_step": (
            diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
        ),
        "lattice_mach_gate_passed": (
            diagnostics["maximum_lattice_mach_number"]
            <= diagnostics["maximum_allowed_lattice_mach_number"]
        ),
        "mole_fraction_bounded": bool(
            np.nanmin(result.species_mole_fraction) >= -1.0e-12
            and np.nanmax(result.maximum_species_mole_fraction) <= 1.0 + 1.0e-12
        ),
        "initial_maximum_principle_within_corrector_tolerance": bool(
            np.nanmax(result.maximum_species_mole_fraction)
            <= (
                cloud.species_mole_fraction
                + config.thermodynamic_corrector_density_tolerance
            )
        ),
        "report_wall_height_exact_on_grid": abs(effective_height - 5.0) <= 1.0e-12,
    }
    numerical_gates = {
        key: value for key, value in gates.items()
        if key != "report_wall_height_exact_on_grid"
    }
    passed_numerics = all(numerical_gates.values())
    passed_all = all(gates.values())
    record = {
        "schema": "degali.thorney-ti021-full-scale-preflight.v1",
        "status": (
            "accepted_full_scale_single_step_preflight" if passed_all else (
                "conditional_grid_resolution_preflight"
                if passed_numerics else "withheld_full_scale_preflight"
            )
        ),
        "scope": (
            "report-scale geometry and coupled numerical execution only; this "
            "short horizon does not validate dispersion, dose or turbulence"
        ),
        "property_audit": {
            "path": str(property_path), "sha256": _sha256(property_path),
        },
        "configuration": {
            "cell_size_m": dx,
            "domain_m": list(domain),
            "cells": list(cells),
            "cell_count": int(np.prod(cells)),
            "carrier_flow_initialization": diagnostics[
                "carrier_flow_initialization"
            ],
            "d3q19_spinup_steps_requested_placeholder": args.flow_steps,
            "d3q19_spinup_steps_completed": diagnostics[
                "flow_steps_completed"
            ],
            "active_steps": args.active_steps,
            "active_time_step_s": active_dt,
            "active_duration_s": duration,
            "wind_speed_m_s": 3.9,
            "streamwise_wind_m_s": streamwise,
            "crosswind_m_s": crosswind,
            "friction_velocity_m_s": 0.205,
            "roughness_length_m": 0.005,
            "turbulence_closure": (
                "neutral_monin_obukhov_background_surface_layer"
                if args.surface_layer_mixing else "none_preflight_only"
            ),
            "effective_schmidt_number": (
                args.effective_schmidt if args.surface_layer_mixing else None
            ),
            "effective_prandtl_number": (
                args.effective_prandtl if args.surface_layer_mixing else None
            ),
            "projection_solver_requested": args.projection_solver,
        },
        "geometry": {
            "cloud_exact_volume_m3": cloud.exact_volume_m3,
            "represented_cloud_volume_m3": diagnostics[
                "initial_cloud_represented_geometry_volume_m3"
            ],
            "initial_r12_mass_kg": result.domain_mass_kg[0],
            "wall_declared_height_m": 5.0,
            "wall_effective_voxel_height_m": effective_height,
            "wall_voxelization": voxel,
        },
        "numerics": {
            "minimum_final_species_mole_fraction": float(
                np.nanmin(result.species_mole_fraction)
            ),
            "maximum_final_species_mole_fraction": float(
                np.nanmax(result.species_mole_fraction)
            ),
            "maximum_over_time_species_mole_fraction": float(
                np.nanmax(result.maximum_species_mole_fraction)
            ),
            "maximum_relative_species_mass_residual": diagnostics[
                "maximum_relative_mass_residual"
            ],
            "maximum_continuity_residual_kg_m3_s": diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
            "maximum_global_mass_residual_kg_s": diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "maximum_relative_density_departure": diagnostics[
                "maximum_relative_mixture_density_departure"
            ],
            "maximum_lattice_mach_number": diagnostics[
                "maximum_lattice_mach_number"
            ],
            "maximum_momentum_advective_courant_number": diagnostics[
                "maximum_momentum_advective_courant_number"
            ],
            "maximum_momentum_viscous_stability_number": diagnostics[
                "maximum_momentum_viscous_stability_number"
            ],
            "maximum_corrector_iterations": diagnostics[
                "thermodynamic_corrector_maximum_iterations_used"
            ],
            "corrector_iterations_per_step": diagnostics[
                "thermodynamic_corrector_iterations_per_step"
            ],
            "corrector_nonconverged_steps": diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ],
            "projection_maximum_linear_iterations_per_step": diagnostics[
                "variable_density_projection_maximum_linear_iterations_per_step"
            ],
            "projection_linear_iterations_per_step": diagnostics[
                "variable_density_projection_linear_iterations_per_step"
            ],
            "projection_solver_used": diagnostics[
                "variable_density_linear_solver"
            ],
            "projection_preconditioner": diagnostics[
                "variable_density_preconditioner"
            ],
            "projection_preconditioner_level_count": diagnostics[
                "variable_density_preconditioner_level_count"
            ],
            "projection_preconditioner_operator_complexity": diagnostics[
                "variable_density_preconditioner_operator_complexity"
            ],
            "pyamg_version": diagnostics[
                "variable_density_pyamg_version"
            ],
            "background_surface_layer_mixing": diagnostics[
                "background_surface_layer_mixing"
            ],
            "background_eddy_viscosity_range_m2_s": [
                diagnostics["background_eddy_viscosity_minimum_m2_s"],
                diagnostics["background_eddy_viscosity_maximum_m2_s"],
            ],
            "effective_species_diffusivity_range_m2_s": [
                diagnostics["applied_scalar_diffusivity_minimum_m2_s"],
                diagnostics["applied_scalar_diffusivity_maximum_m2_s"],
            ],
            "effective_thermal_diffusivity_range_m2_s": [
                diagnostics["applied_thermal_diffusivity_minimum_m2_s"],
                diagnostics["applied_thermal_diffusivity_maximum_m2_s"],
            ],
        },
        "runtime": {
            "flow_solve_s": diagnostics["flow_solve_runtime_s"],
            "transport_s": diagnostics["transport_runtime_s"],
            "total_s": diagnostics["total_runtime_s"],
        },
        "gates": gates,
        "all_numerical_gates_pass": passed_numerics,
        "all_gates_pass": passed_all,
        "next_required": (
            (
                "repeat at 1 m so the 5 m wall is exact, then freeze turbulent "
                "mixing without TI021 concentration and run grid/time/domain audits"
            ) if dx > 1.0 else (
                "freeze turbulent mixing without TI021 concentration, then run "
                "the full transient grid/time/domain audits and no-fit sensor score"
            )
        ),
    }

    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "ti021-full-scale-preflight.json"
    plot_path = output / "ti021-full-scale-preflight.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    assert result.concentration_fields_kg_m3 is not None
    initial_state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        result.concentration_fields_kg_m3[0],
        np.zeros_like(result.concentration_fields_kg_m3[0]),
    )
    k = int(np.argmin(abs(result.z_m - 0.5 * dx)))
    panels = (
        (100.0 * initial_state.species_mole_fraction[k], "initial R12 (vol%)", "viridis"),
        (result.concentration_vol_percent[k], "final R12 (vol%)", "viridis"),
        (result.velocity_w_m_s[k], "final vertical velocity (m/s)", "coolwarm"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5), constrained_layout=True)
    for axis, (field, title, cmap) in zip(axes, panels):
        image = axis.pcolormesh(
            result.x_m, result.y_m, field, shading="nearest", cmap=cmap,
        )
        axis.contour(
            result.x_m, result.y_m, result.solid_mask[k].astype(float),
            levels=[0.5], colors="black", linewidths=1.0,
        )
        axis.set_aspect("equal")
        axis.set_title(title)
        axis.set_xlabel("solver x / site downwind (m)")
        fig.colorbar(image, ax=axis, shrink=0.78)
    axes[0].set_ylabel("solver y / site crosswind (m)")
    fig.suptitle(
        f"TI021 full-site numerical preflight — dx={dx:g} m, {duration:g} s"
    )
    fig.savefig(plot_path, dpi=170)
    plt.close(fig)

    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "cell_count": record["configuration"]["cell_count"],
        "runtime": record["runtime"],
        "gates": gates,
    }, indent=2, ensure_ascii=False))
    return 0 if passed_numerics else 2


if __name__ == "__main__":
    raise SystemExit(main())
