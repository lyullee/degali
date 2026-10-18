"""Audit the full-mixture finite-volume initial state for Thorney TI021."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud,
    initialise_regular_polygon_prism_cloud,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule


ROOT = Path(__file__).resolve().parents[1]


def _centres(start: float, stop: float, spacing: float) -> np.ndarray:
    count = int(round((stop - start) / spacing))
    if not np.isclose(start + count * spacing, stop):
        raise ValueError("domain extent must be an integer number of cells")
    return start + (np.arange(count) + 0.5) * spacing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-directory", type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-initial-cloud-audit-2026-10-11-v1"
        ),
    )
    args = parser.parse_args()

    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0748985077907,
        ambient_temperature_k=291.15,
        ambient_pressure_pa=100_500.0,
        minimum_temperature_k=150.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(0.0, 0.0),
        flat_to_flat_m=14.0,
        z_min_m=0.0,
        z_max_m=13.0,
        species_mole_fraction=0.32,
        temperature_k=291.15,
        sides=12,
        label="Thorney Island Trial 021 initial R12/air cloud",
    )
    cloud_state = closure.state_from_mole_fraction_and_temperature(
        cloud.species_mole_fraction, cloud.temperature_k,
    )
    spacings = (1.0, 0.5, 0.25)
    grid_records = []
    initial_states = []
    for spacing in spacings:
        x = _centres(-10.0, 10.0, spacing)
        y = _centres(-10.0, 10.0, spacing)
        z = _centres(0.0, 14.0, spacing)
        initial = initialise_regular_polygon_prism_cloud(
            x, y, z, cloud=cloud, closure=closure,
        )
        initial_states.append((spacing, x, y, z, initial))
        grid_records.append({
            "spacing_m": spacing,
            "cells": [len(x), len(y), len(z)],
            "represented_geometry_volume_m3": (
                initial.represented_geometry_volume_m3
            ),
            "relative_geometry_volume_error": (
                initial.relative_geometry_volume_error
            ),
            "integrated_r12_mass_kg": initial.integrated_species_mass_kg,
            "relative_r12_mass_error": initial.relative_species_mass_error,
            "full_cloud_cell_count": int(np.count_nonzero(
                initial.volume_fraction >= 1.0 - 1.0e-12
            )),
            "partial_cloud_cell_count": int(np.count_nonzero(
                (initial.volume_fraction > 1.0e-12)
                & (initial.volume_fraction < 1.0 - 1.0e-12)
            )),
            "maximum_reconstructed_mole_fraction": float(np.max(
                initial.thermodynamic_state.species_mole_fraction
            )),
        })

    reported_volume_m3 = 2050.0
    reported_effective_molar_mass = 58.5e-3
    reported_relative_density = 2.02
    effective_molar_mass = float(cloud_state.mixture_molar_mass_kg_mol)
    relative_density = float(
        cloud_state.mixture_density_kg_m3 / closure.ambient_density_kg_m3
    )
    canonical_cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0), flat_to_flat_m=4.0,
        z_min_m=0.0, z_max_m=3.0,
        species_mole_fraction=0.32, temperature_k=291.15,
        label="TI021 reduced initial-cloud integration check",
    )
    canonical_result = solve_resolved_obstacle_flow_3d(
        ResolvedObstacleFlow3DConfig(
            domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
            wind_speed_m_s=0.2, kinematic_viscosity_m2_s=0.1,
            flow_steps=5, minimum_flow_steps=5,
            convergence_check_interval=1, lattice_inlet_speed=0.02,
            duration_s=0.01, scalar_time_step_s=0.01,
            scalar_diffusivity_m2_s=0.01, thermal_diffusivity_m2_s=0.01,
            source_schedule=SourceRateSchedule(
                (0.0, 0.01), (0.0, 0.0), source_id="initial-cloud-only",
            ),
            source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
            obstacles=(), ambient_temperature_k=closure.ambient_temperature_k,
            ambient_pressure_pa=closure.ambient_pressure_pa,
            ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
            gravity_m_s2=closure.gravity_m_s2,
            thermal_buoyancy_coupling=True,
            active_momentum_solver="finite_volume_projection",
            low_mach_mass_flux_correction=True,
            initial_binary_cloud=canonical_cloud,
            binary_mixture_closure=closure,
            thermodynamic_corrector_maximum_iterations=32,
            thermodynamic_corrector_minimum_iterations=2,
            thermodynamic_corrector_density_tolerance=1.0e-6,
            thermodynamic_corrector_velocity_tolerance=1.0e-6,
            require_thermodynamic_corrector_convergence=True,
            subgrid_model="none", smagorinsky_constant=0.0,
        )
    )
    canonical_diagnostics = canonical_result.diagnostics
    gates = {
        "dodecagonal_volume_matches_report_within_0p2pct": (
            abs(cloud.exact_volume_m3 / reported_volume_m3 - 1.0) <= 0.002
        ),
        "effective_molar_mass_matches_workbook_within_0p5pct": (
            abs(effective_molar_mass / reported_effective_molar_mass - 1.0)
            <= 0.005
        ),
        "relative_density_matches_report_within_0p5pct": (
            abs(relative_density / reported_relative_density - 1.0) <= 0.005
        ),
        "finite_volume_geometry_exact_on_all_grids": all(
            abs(item[4].relative_geometry_volume_error) <= 1.0e-11
            for item in initial_states
        ),
        "finite_volume_species_mass_exact_on_all_grids": all(
            abs(item[4].relative_species_mass_error) <= 1.0e-11
            for item in initial_states
        ),
        "full_binary_canonical_step_conservative": (
            canonical_diagnostics["maximum_relative_mass_residual"] <= 1.0e-10
            and canonical_diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ] <= 1.0e-8
        ),
        "full_binary_canonical_step_corrector_converged": (
            canonical_diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ] == 0
        ),
        "full_binary_canonical_step_exceeds_active_scalar_density_range": (
            canonical_diagnostics[
                "maximum_relative_mixture_density_departure"
            ] > 1.0
        ),
    }
    record = {
        "schema": "degali.thorney-ti021-initial-cloud-audit.v1",
        "status": (
            "accepted_initial_state_component"
            if all(gates.values()) else "withheld_initial_state_component"
        ),
        "scope": (
            "geometry, composition, EOS and finite-volume initialisation only; "
            "not a dispersion or obstacle-validation result"
        ),
        "source_policy": {
            "reported_cloud_volume_m3": reported_volume_m3,
            "reported_holder_description": (
                "twelve-sided, 14 m across, 14 m high, filled to 13 m"
            ),
            "across_interpretation": (
                "flat-to-flat; this resolves the ambiguous word 'across' because "
                "the resulting 2048.2 m3 independently matches the reported 2050 m3"
            ),
            "initial_r12_mole_fraction": 0.32,
            "ambient_temperature_k": closure.ambient_temperature_k,
            "ambient_pressure_pa": closure.ambient_pressure_pa,
            "r12_heat_capacity_note": (
                "606.0749 J/kg/K is the CoolProp 8.0.0 R12 value at the event "
                "state and is inactive at the ambient-temperature initial state; "
                "see the separate R12 property audit"
            ),
        },
        "reconstructed_cloud": {
            "exact_regular_dodecagonal_volume_m3": cloud.exact_volume_m3,
            "volume_difference_from_report_percent": (
                100.0 * (cloud.exact_volume_m3 / reported_volume_m3 - 1.0)
            ),
            "effective_molar_mass_kg_mol": effective_molar_mass,
            "effective_molar_mass_difference_percent": (
                100.0 * (
                    effective_molar_mass / reported_effective_molar_mass - 1.0
                )
            ),
            "relative_density": relative_density,
            "relative_density_difference_percent": (
                100.0 * (relative_density / reported_relative_density - 1.0)
            ),
            "mixture_density_kg_m3": float(cloud_state.mixture_density_kg_m3),
            "r12_mass_fraction": float(cloud_state.species_mass_fraction),
            "r12_mass_concentration_kg_m3": float(
                cloud_state.species_mass_concentration_kg_m3
            ),
            "exact_initial_r12_mass_kg": float(
                cloud_state.species_mass_concentration_kg_m3
                * cloud.exact_volume_m3
            ),
        },
        "grid_records": grid_records,
        "reduced_canonical_coupled_step": {
            "scope": (
                "one 0.01 s full-mixture finite-volume momentum/species/enthalpy/"
                "EOS/total-mass-projection integration check; not TI021 validation"
            ),
            "resolved_solver_schema": canonical_result.schema,
            "species_name": canonical_result.species_name,
            "initial_species_mass_kg": float(
                canonical_result.domain_mass_kg[0]
            ),
            "maximum_relative_species_mass_residual": (
                canonical_diagnostics["maximum_relative_mass_residual"]
            ),
            "maximum_continuity_residual_kg_m3_s": canonical_diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
            "maximum_global_mass_residual_kg_s": canonical_diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "corrector_nonconverged_steps": canonical_diagnostics[
                "thermodynamic_corrector_nonconverged_steps"
            ],
            "maximum_relative_density_departure": canonical_diagnostics[
                "maximum_relative_mixture_density_departure"
            ],
            "total_runtime_s": canonical_diagnostics["total_runtime_s"],
        },
        "gates": gates,
        "remaining_before_dispersion_score": [
            "full-scale domain and wall/inlet placement",
            "attach the accepted R12 molecular-property audit",
            "independently freeze the turbulent scalar-mixing closure",
            "grid/time/domain convergence of the coupled transient",
            "no-fit sensor dose and arc-summary scoring",
        ],
    }

    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-initial-cloud-audit.json"
    plot_path = args.output_directory / "ti021-initial-cloud-grid-comparison.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    for axis, (spacing, x, y, z, initial), grid_record in zip(
        axes, initial_states, grid_records,
    ):
        k = int(np.argmin(np.abs(z - 6.5)))
        image = axis.pcolormesh(
            x, y, initial.volume_fraction[k],
            shading="nearest", cmap="Blues", vmin=0.0, vmax=1.0,
        )
        vertices = cloud.vertices_xy_m()
        closed = np.vstack((vertices, vertices[0]))
        axis.plot(closed[:, 0], closed[:, 1], color="crimson", linewidth=1.4)
        axis.set_aspect("equal")
        axis.set_xlim(-8.5, 8.5)
        axis.set_ylim(-8.5, 8.5)
        axis.set_title(
            f"dx={spacing:g} m\n"
            f"volume error={100*grid_record['relative_geometry_volume_error']:.2e}%"
        )
        axis.set_xlabel("x relative to release (m)")
        axis.grid(alpha=0.15)
    axes[0].set_ylabel("y relative to release (m)")
    bar = fig.colorbar(image, ax=axes, shrink=0.82)
    bar.set_label("cell cloud-volume fraction")
    fig.suptitle(
        "TI021 initial 12-sided R12/air cloud — finite-volume representation"
    )
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "exact_volume_m3": cloud.exact_volume_m3,
        "reported_volume_m3": reported_volume_m3,
        "effective_molar_mass_kg_mol": effective_molar_mass,
        "relative_density": relative_density,
        "gates": gates,
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
