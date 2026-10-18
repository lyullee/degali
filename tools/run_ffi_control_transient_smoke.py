"""Run a short no-fit direct-H2 control-channel startup for the FFI facility.

This is a solver-path smoke, not an FFI sensor score: it uses the published
article A-case source rate in a zero-wind, closed-wall/open-end channel and
stops before the experimental 60--120 s averaging window.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.high_pressure_source_handoff_3d import (
    declared_ambient_pressure_hydrogen_jet,
    low_mach_handoff_from_ambient_pressure_jet,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs/ffi-control-transient-smoke-2026-10-13-v3"
H2_MOLAR_MASS_KG_MOL = 2.01588e-3
AIR_MOLAR_MASS_KG_MOL = 28.96546e-3
R_GAS = 8.314462618


def build_config(*, duration_s: float = 0.01, cell_size_m: float = 0.025):
    """Construct a deliberately short FFI A-control source/EOS smoke case."""
    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise ValueError("duration_s must be positive and finite")
    if not math.isfinite(cell_size_m) or cell_size_m <= 0.0:
        raise ValueError("cell_size_m must be positive and finite")
    domain = (5.8, 0.9, 0.8)  # FFI y, x, z mapped to solver streamwise/lateral/z
    cells = tuple(int(round(length / cell_size_m)) for length in domain)
    actual = tuple(length / count for length, count in zip(domain, cells))
    if not np.allclose(actual, actual[0], rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("cell_size_m must divide all mapped FFI channel dimensions")
    temperature_k = 293.0
    pressure_pa = 101_325.0
    source_rate_kg_s = 0.086e-3  # published LES input for paired control A
    orifice_diameter_m = 0.0046
    source_density_kg_m3 = pressure_pa * H2_MOLAR_MASS_KG_MOL / (R_GAS * temperature_k)
    source_speed_m_s = source_rate_kg_s / (
        source_density_kg_m3 * math.pi * (orifice_diameter_m / 2.0) ** 2
    )
    # An explicitly coarse regularized source width is retained as a smoke
    # setting; it is not a claim that the 4.6 mm nozzle is resolved.
    port_density = pressure_pa * H2_MOLAR_MASS_KG_MOL / (R_GAS * temperature_k)
    port_speed = source_rate_kg_s / (
        port_density * math.pi * (orifice_diameter_m / 2.0) ** 2
    )
    port = declared_ambient_pressure_hydrogen_jet(
        source_id="FFI-control-A-published-LES-source-rate-no-fit",
        source_model=(
            "ambient-pressure FFI source: published mean mass rate and port diameter; "
            "uniform exit state inferred at ambient T/P"
        ),
        mass_flow_kg_s=source_rate_kg_s,
        diameter_m=orifice_diameter_m,
        velocity_m_s=port_speed,
        density_kg_m3=port_density,
        temperature_k=temperature_k,
        ambient_pressure_pa=pressure_pa,
    )
    handoff = low_mach_handoff_from_ambient_pressure_jet(
        port,
        ambient_temperature_k=temperature_k,
        ambient_pressure_pa=pressure_pa,
        ambient_density_kg_m3=pressure_pa * AIR_MOLAR_MASS_KG_MOL / (R_GAS * temperature_k),
        ambient_axial_velocity_m_s=0.0,
        axis_unit_vector=(0.0, 0.0, -1.0),
        maximum_mach_number=0.30,
        maximum_relative_density_departure=0.10,
    )
    source = handoff.resolved_source_fields(
        nozzle_position_m=(0.5, 0.45, 0.8),
        grid_spacing_m=actual[0],
        duration_s=duration_s,
        minimum_cells_across_diameter=4.0,
        # The 4.6 mm orifice is intentionally regularised on this coarse
        # startup grid. This smoke is not nozzle-resolved validation.
        require_spatial_resolution=False,
    )
    source_sigma_m = source.source_sigma_m
    closure = BinaryIdealGasMixtureClosure(
        species_name="hydrogen",
        species_molar_mass_kg_mol=H2_MOLAR_MASS_KG_MOL,
        species_heat_capacity_j_kg_k=14_300.0,
        carrier_molar_mass_kg_mol=AIR_MOLAR_MASS_KG_MOL,
        carrier_heat_capacity_j_kg_k=1005.0,
        ambient_temperature_k=temperature_k,
        ambient_pressure_pa=pressure_pa,
    )
    initial_air = RegularPolygonPrismInitialCloud(
        center_xy_m=(2.9, 0.45), flat_to_flat_m=0.02,
        z_min_m=0.001, z_max_m=0.002,
        species_mole_fraction=0.0, temperature_k=temperature_k,
        label="ambient-air-initialization-no-hydrogen",
    )
    return ResolvedObstacleFlow3DConfig(
        domain_m=domain,
        cells=cells,
        wind_speed_m_s=0.0,
        velocity_reference_scale_m_s=2.0,
        kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=1,
        minimum_flow_steps=1,
        lattice_inlet_speed=0.03,
        duration_s=duration_s,
        scalar_time_step_s=5.0e-3,
        scalar_diffusivity_m2_s=7.0e-5,
        thermal_diffusivity_m2_s=2.2e-5,
        source_schedule=source.source_schedule,
        source_position_m=source.source_position_m,
        source_sigma_m=source_sigma_m,
        source_sigma_xyz_m=source.source_sigma_xyz_m,
        obstacles=(),
        source_momentum_velocity_m_s=source.source_momentum_velocity_m_s,
        nonstreamwise_wall_boundary="all_no_slip_tunnel",
        carrier_flow_initialization="projected_inlet",
        active_momentum_solver="finite_volume_projection",
        thermal_buoyancy_coupling=True,
        source_temperature_k=temperature_k,
        ambient_temperature_k=temperature_k,
        ambient_pressure_pa=pressure_pa,
        ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
        binary_mixture_closure=closure,
        initial_binary_cloud=initial_air,
        maximum_active_scalar_density_departure=1.0,
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.85,
        subgrid_model="smagorinsky",
        smagorinsky_constant=0.12,
        scalar_turbulent_diffusivity_model="molecular_plus_subgrid_eddy_viscosity",
        thermal_turbulent_diffusivity_model="molecular_plus_subgrid_eddy_viscosity",
        require_thermodynamic_corrector_convergence=True,
        store_scalar_fields=True,
    )


def run(output: Path, *, duration_s: float = 0.01, cell_size_m: float = 0.1) -> dict[str, object]:
    if output.exists():
        raise FileExistsError(output)
    config = build_config(duration_s=duration_s, cell_size_m=cell_size_m)
    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(config)
    runtime_s = time.perf_counter() - started
    output.mkdir(parents=True, exist_ok=False)
    field_path = output / "ffi-control-A-smoke-fields.npz"
    np.savez_compressed(
        field_path,
        x_m=result.x_m, y_m=result.y_m, z_m=result.z_m,
        solid_mask=result.solid_mask,
        hydrogen_vol_percent=result.concentration_vol_percent,
        maximum_hydrogen_vol_percent=result.maximum_concentration_vol_percent,
        u_m_s=result.velocity_u_m_s, v_m_s=result.velocity_v_m_s,
        w_m_s=result.velocity_w_m_s,
        pressure_gauge_pa=result.pressure_pa_gauge,
        temperature_k=result.temperature_k,
        density_kg_m3=result.mixture_density_kg_m3,
    )
    record = {
        "schema": "degali.ffi-control-transient-smoke.v1",
        "status": "solver_path_smoke_only_not_experimental_validation",
        "case": "unobstructed FFI control A geometry",
        "configuration": {
            "domain_solver_axes_m": list(config.domain_m),
            "ffi_axis_mapping": {"solver_x": "FFI channel y / streamwise", "solver_y": "FFI channel x / lateral", "solver_z": "FFI channel z / vertical"},
            "cells": list(config.cells),
            "cell_size_m": config.domain_m[0] / config.cells[0],
            "duration_s": duration_s,
            "scalar_time_step_s": config.scalar_time_step_s,
            "source_rate_kg_s": config.source_schedule.rate_kg_s[0],
            "source_schedule_id": config.source_schedule.source_id,
            "source_exit_speed_inferred_m_s": abs(config.source_momentum_velocity_m_s[2]),
            "low_mach_handoff_distance_m": 0.8 - config.source_position_m[2],
            "low_mach_handoff_equivalent_diameter_m": 4.0 * config.source_sigma_m,
            "low_mach_handoff_relative_density_departure_target": -0.10,
            "low_mach_source_cells_across_equivalent_diameter": (
                4.0 * config.source_sigma_m / (config.domain_m[0] / config.cells[0])
            ),
            "low_mach_source_momentum_increment_n": (
                config.source_schedule.rate_kg_s[0]
                * abs(config.source_momentum_velocity_m_s[2])
            ),
            "source_sigma_xyz_m": list(config.source_sigma_xyz_m),
            "source_nozzle_resolved": False,
            "obstacle_geometry_included": False,
            "ambient_wind_m_s": 0.0,
            "wall_boundary": config.nonstreamwise_wall_boundary,
            "active_momentum_solver": config.active_momentum_solver,
            "binary_mixture_eos": "ideal H2/dry-air with sensible enthalpy and buoyancy",
        },
        "execution": {
            "wall_runtime_s": runtime_s,
            "solver_runtime_s": result.diagnostics.get("total_runtime_s"),
            "time_steps": int(len(result.time_s) - 1),
            "final_time_s": float(result.time_s[-1]),
            "field_npz": str(field_path.resolve()),
            "max_hydrogen_vol_pct": float(np.max(result.maximum_concentration_vol_percent)),
            "max_velocity_m_s": float(np.max(np.sqrt(
                result.velocity_u_m_s**2 + result.velocity_v_m_s**2 + result.velocity_w_m_s**2
            ))),
            "min_temperature_k": float(np.min(result.temperature_k)),
            "max_temperature_k": float(np.max(result.temperature_k)),
            "max_abs_mass_residual_kg": float(np.max(np.abs(result.mass_residual_kg))),
            "max_abs_energy_residual_j": float(np.max(np.abs(result.energy_residual_j))),
            "max_continuity_residual_kg_m3_s": float(np.max(np.abs(result.continuity_residual_kg_m3_s))),
            "solver_diagnostics": result.diagnostics,
        },
        "comparison": {
            "sensor_score_performed": False,
            "matched_experimental_window": False,
            "promotion_allowed": False,
            "limitations": [
                "run stops before FFI's published 60--120 s mean window",
                "control case only; no obstacle geometry is supplied",
                "4.6 mm nozzle exit speed is inferred from reported mean rate, ideal-gas density, and circular port area",
                "near-field low-Mach handoff is an entrainment model, not a measured station",
                "handoff source is Gaussian-regularized; only four-cell equivalent-diameter resolution is required here",
                "not an outdoor or cryogenic LH2 validation",
            ],
        },
    }
    (output / "ffi-control-A-smoke.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--duration-s", type=float, default=0.01)
    parser.add_argument("--cell-size-m", type=float, default=0.025)
    args = parser.parse_args()
    result = run(args.output, duration_s=args.duration_s, cell_size_m=args.cell_size_m)
    print(json.dumps({
        "status": result["status"],
        "wall_runtime_s": result["execution"]["wall_runtime_s"],
        "steps": result["execution"]["time_steps"],
        "max_hydrogen_vol_pct": result["execution"]["max_hydrogen_vol_pct"],
        "max_abs_mass_residual_kg": result["execution"]["max_abs_mass_residual_kg"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
