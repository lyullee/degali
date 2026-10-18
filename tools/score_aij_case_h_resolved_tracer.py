"""No-fit 3-D passive-tracer component comparison for AIJ Case H.

This pilot reuses the resolved 3-D AIJ geometry/inlet and the public ethylene
source rate.  The solver stores scalar mass concentration in its hydrogen
species slot; that slot is used only as a passive conserved scalar here, and
the output is converted using ethylene's molecular weight.  It is not an H2,
TKE, or facility-qualified concentration validation.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.validation.aij_case_h import read_aij_case_h
from tools.audit_aij_resolved_3d_grid import _config


DATA = ROOT / "reference/aij_case_h"
DEFAULT_OUTPUT = ROOT / "outputs/aij-case-h-resolved-tracer-coarse-2026-10-12-v1"
X_ORIGIN_M = -0.25
Y_ORIGIN_M = -0.15
ETHYLENE_MOLAR_MASS_KG_MOL = 28.0532e-3
AIR_MOLAR_MASS_KG_MOL = 28.96546e-3
R_GAS_J_MOL_K = 8.314462618
SOURCE_ORIFICE_DIAMETER_M = 0.004


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bracket(axis: np.ndarray, value: float) -> tuple[int, int, float] | None:
    if value < axis[0] - 1.0e-12 or value > axis[-1] + 1.0e-12:
        return None
    upper = int(np.searchsorted(axis, value, side="right"))
    upper = min(max(upper, 1), len(axis) - 1)
    lower = upper - 1
    fraction = float((value - axis[lower]) / (axis[upper] - axis[lower]))
    return lower, upper, fraction


def _sample_fluid(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    field_zyx: np.ndarray,
    solid_zyx: np.ndarray,
    point_xyz: tuple[float, float, float],
) -> float | None:
    bx, by, bz = (
        _bracket(axis, coordinate) for axis, coordinate in zip((x, y, z), point_xyz)
    )
    if bx is None or by is None or bz is None:
        return None
    x0, x1, fx = bx
    y0, y1, fy = by
    z0, z1, fz = bz
    weighted = []
    for k, wz in ((z0, 1.0 - fz), (z1, fz)):
        for j, wy in ((y0, 1.0 - fy), (y1, fy)):
            for i, wx in ((x0, 1.0 - fx), (x1, fx)):
                weight = wx * wy * wz
                if weight > 1.0e-14:
                    if solid_zyx[k, j, i]:
                        return None
                    weighted.append(weight * float(field_zyx[k, j, i]))
    return float(sum(weighted))


def _ethylene_ppmv_from_mass_concentration(
    ethylene_kg_m3: float, *, temperature_k: float, pressure_pa: float,
) -> float:
    """Convert dilute ethylene mass concentration to ideal-gas ppmv."""
    if not math.isfinite(ethylene_kg_m3) or ethylene_kg_m3 < 0.0:
        raise ValueError("ethylene mass concentration must be finite and nonnegative")
    gas_moles_m3 = pressure_pa / (R_GAS_J_MOL_K * temperature_k)
    ethylene_moles_m3 = ethylene_kg_m3 / ETHYLENE_MOLAR_MASS_KG_MOL
    return 1.0e6 * ethylene_moles_m3 / (gas_moles_m3 + ethylene_moles_m3)


def _pearson(left: np.ndarray, right: np.ndarray) -> float:
    if left.size < 2 or np.std(left) == 0.0 or np.std(right) == 0.0:
        return 0.0
    return float(np.corrcoef(left, right)[0, 1])


def _source_exit_velocity_m_s(volume_flow_m3_s: float) -> float:
    """Infer a uniform mean exit speed from reported flow and 4 mm opening."""
    area = math.pi * (SOURCE_ORIFICE_DIAMETER_M / 2.0) ** 2
    return volume_flow_m3_s / area


def _source_sigma_for_run(cell_size_m: float, requested_sigma_m: float | None) -> float:
    """Resolve source regularization without silently tying sensitivity to dx."""
    cell = float(cell_size_m)
    if not math.isfinite(cell) or cell <= 0.0:
        raise ValueError("cell_size_m must be positive and finite")
    sigma = cell if requested_sigma_m is None else float(requested_sigma_m)
    if not math.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("source_sigma_m must be positive and finite")
    return sigma


def run(
    *, output_directory: Path = DEFAULT_OUTPUT, duration_s: float = 20.0,
    cell_size_m: float = 0.025, scalar_time_step_s: float = 0.002,
    source_sigma_m: float | None = None,
    average_start_s: float = 15.0, use_sgs_scalar_mixing: bool = False,
    use_binary_fv_source_momentum: bool = False,
    use_binary_fv_neutral_source: bool = False,
    lateral_domain_m: float | None = None,
    y_origin_m: float = Y_ORIGIN_M,
) -> dict[str, object]:
    dataset = read_aij_case_h(DATA)
    if not 0.0 < average_start_s < duration_s:
        raise ValueError("average_start_s must lie within the scalar run")
    resolved_source_sigma_m = _source_sigma_for_run(cell_size_m, source_sigma_m)

    # The README gives a 100%-ethylene volumetric source.  Convert at declared
    # neutral tunnel ambient and keep that conversion frozen (no concentration fit).
    source_volume_flow_m3_s = 5.83e-6
    source_density_kg_m3 = (
        101_325.0 * ETHYLENE_MOLAR_MASS_KG_MOL
        / (R_GAS_J_MOL_K * 293.0)
    )
    source_mass_rate_kg_s = source_volume_flow_m3_s * source_density_kg_m3
    source_exit_velocity_m_s = _source_exit_velocity_m_s(source_volume_flow_m3_s)
    base = _config(dataset, cell_size_m)
    if lateral_domain_m is not None:
        lateral_domain_m = float(lateral_domain_m)
        if not math.isfinite(lateral_domain_m) or lateral_domain_m <= 0.0:
            raise ValueError("lateral_domain_m must be positive and finite")
        y_origin_m = float(y_origin_m)
        if not math.isfinite(y_origin_m):
            raise ValueError("y_origin_m must be finite")
        shift = Y_ORIGIN_M - y_origin_m
        lateral_cells = int(round(lateral_domain_m / cell_size_m))
        if lateral_cells < 8 or not math.isclose(
            lateral_cells * cell_size_m, lateral_domain_m, rel_tol=0.0, abs_tol=1.0e-10
        ):
            raise ValueError("lateral_domain_m must be an integer number of cells and provide at least eight cells")
        base = replace(
            base,
            domain_m=(base.domain_m[0], lateral_domain_m, base.domain_m[2]),
            cells=(base.cells[0], lateral_cells, base.cells[2]),
            source_position_m=(
                base.source_position_m[0],
                base.source_position_m[1] + shift,
                base.source_position_m[2],
            ),
            obstacles=tuple(
                replace(
                    obstacle,
                    y_min_m=obstacle.y_min_m + shift,
                    y_max_m=obstacle.y_max_m + shift,
                )
                for obstacle in base.obstacles
            ),
        )
    flow_dt = base.lattice_inlet_speed * cell_size_m / base.wind_speed_m_s
    flow_steps = int(round(0.5572755417956656 / flow_dt))
    averaging_steps = max(2, int(round(0.18575851393188854 / flow_dt)))
    schedule = SourceRateSchedule(
        time_s=(0.0, 0.02, duration_s - 0.02, duration_s),
        rate_kg_s=(0.0, source_mass_rate_kg_s, source_mass_rate_kg_s, 0.0),
        source_id="AIJ-Case-H-published-ethylene-volume-flow-no-fit",
        rate_operator="piecewise_constant",
    )
    config_changes = {
        "duration_s": duration_s,
        "scalar_time_step_s": scalar_time_step_s,
        "source_schedule": schedule,
        "source_position_m": (0.30, 0.15, cell_size_m / 2.0),
        "source_sigma_m": resolved_source_sigma_m,
        "flow_steps": flow_steps,
        "minimum_flow_steps": flow_steps,
        "flow_average_last_steps": averaging_steps,
        "flow_statistical_convergence_tolerance": 0.01,
        "subgrid_model": "smagorinsky" if use_sgs_scalar_mixing else "none",
        "smagorinsky_constant": 0.12 if use_sgs_scalar_mixing else 0.0,
        "effective_schmidt_number": 0.7 if use_sgs_scalar_mixing else None,
        "scalar_turbulent_diffusivity_model": (
            "molecular_plus_subgrid_eddy_viscosity"
            if use_sgs_scalar_mixing else "legacy_total_kinematic_viscosity_ratio"
        ),
        "transient_average_start_s": average_start_s,
        "ambient_temperature_k": 293.0,
        "ambient_pressure_pa": 101_325.0,
        "ambient_air_density_kg_m3": 1.204,
        "thermal_buoyancy_coupling": False,
        "initial_binary_cloud": None,
        "binary_mixture_closure": None,
    }
    if use_binary_fv_source_momentum or use_binary_fv_neutral_source:
        # This is an ethylene/air ideal-gas mixture, not the solver's H2 path.
        # A zero ethylene fraction cloud provides a consistent ambient binary
        # thermodynamic state; the measured ethylene mass source then supplies
        # both species and total mixture mass to the FV mass-flux projection.
        config_changes.update({
            "source_temperature_k": 293.0,
            "ambient_air_density_kg_m3": 101_325.0 * AIR_MOLAR_MASS_KG_MOL
            / (R_GAS_J_MOL_K * 293.0),
            "thermal_buoyancy_coupling": True,
            "low_mach_mass_flux_correction": True,
            "active_momentum_solver": "finite_volume_projection",
            "source_momentum_velocity_m_s": (
                (0.0, 0.0, source_exit_velocity_m_s)
                if use_binary_fv_source_momentum else (0.0, 0.0, 0.0)
            ),
            "binary_mixture_closure": BinaryIdealGasMixtureClosure(
                species_name="ethylene",
                species_molar_mass_kg_mol=ETHYLENE_MOLAR_MASS_KG_MOL,
                species_heat_capacity_j_kg_k=1_500.0,
                ambient_temperature_k=293.0,
                ambient_pressure_pa=101_325.0,
                gravity_m_s2=9.80665,
            ),
            "initial_binary_cloud": RegularPolygonPrismInitialCloud(
                center_xy_m=(0.0, 0.0), flat_to_flat_m=0.02,
                z_min_m=0.001, z_max_m=0.002,
                species_mole_fraction=0.0, temperature_k=293.0,
                label="ambient-air-initial-state-no-ethylene",
            ),
            "subgrid_model": "smagorinsky" if use_sgs_scalar_mixing else "none",
            "smagorinsky_constant": 0.12 if use_sgs_scalar_mixing else 0.0,
            "effective_schmidt_number": 0.7 if use_sgs_scalar_mixing else None,
            "effective_prandtl_number": 0.85,
            "scalar_turbulent_diffusivity_model": (
                "molecular_plus_subgrid_eddy_viscosity"
                if use_sgs_scalar_mixing else "legacy_total_kinematic_viscosity_ratio"
            ),
            "thermal_turbulent_diffusivity_model": (
                "molecular_plus_subgrid_eddy_viscosity"
                if use_sgs_scalar_mixing else "legacy_total_kinematic_viscosity_ratio"
            ),
            "maximum_active_scalar_density_departure": 0.10,
            "require_thermodynamic_corrector_convergence": True,
        })
    config = replace(base, **config_changes)

    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(config)
    wall_clock_s = time.perf_counter() - started
    mean_scalar = result.time_averaged_concentration_kg_m3
    if mean_scalar is None:
        raise RuntimeError("solver did not return the requested scalar time average")

    rows = []
    measured = []
    predicted = []
    skipped = {"missing_measurement": 0, "outside_or_solid_stencil": 0}
    for observation in dataset.observations:
        if observation.concentration_ppm is None:
            skipped["missing_measurement"] += 1
            continue
        model_point = (
            observation.x_m - X_ORIGIN_M,
            observation.y_m - y_origin_m,
            observation.z_m,
        )
        mass_concentration = _sample_fluid(
            result.x_m, result.y_m, result.z_m, mean_scalar,
            result.solid_mask, model_point,
        )
        if mass_concentration is None:
            skipped["outside_or_solid_stencil"] += 1
            continue
        model_ppmv = _ethylene_ppmv_from_mass_concentration(
            max(0.0, mass_concentration), temperature_k=293.0, pressure_pa=101_325.0,
        )
        measured.append(float(observation.concentration_ppm))
        predicted.append(model_ppmv)
        rows.append({
            "point_id": observation.point_id,
            "x_m": observation.x_m,
            "y_m": observation.y_m,
            "z_m": observation.z_m,
            "measured_ethylene_ppmv": observation.concentration_ppm,
            "measured_c_rms_ppmv": observation.concentration_rms_ppm,
            "predicted_mean_ethylene_ppmv": model_ppmv,
            "residual_predicted_minus_measured_ppmv": model_ppmv - observation.concentration_ppm,
        })
    obs = np.asarray(measured, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    if not obs.size:
        raise RuntimeError("no AIJ concentration observations overlap the resolved model")
    residual = pred - obs
    q = source_volume_flow_m3_s
    u_h, h = dataset.reference_speed_m_s, dataset.reference_height_m
    measured_k = obs * 1.0e-6 * u_h * h * h / q
    predicted_k = pred * 1.0e-6 * u_h * h * h / q
    released_mass_kg = schedule.released_mass_kg
    record = {
        "schema": "degali.aij-case-h-resolved-tracer-component-score.v1",
        "status": "conditional_passive_tracer_component_score_not_facility_or_hydrogen_validation",
        "data": {
            "record_doi": dataset.record_doi,
            "readme_sha256": dataset.readme_sha256,
            "result_sha256": dataset.result_sha256,
            "observation_file": "reference/aij_case_h/RS_caseH.csv",
            "source_file": "reference/aij_case_h/readme_caseH.md",
            "source_volume_flow_m3_s": source_volume_flow_m3_s,
            "source_concentration": "100% ethylene per supplied Case H README",
            "sensor_mean_count": int(len(obs)),
            "sensor_rms_is_fluctuation_not_mean_measurement_uncertainty": True,
        },
        "scope": {
            "calibration_used": False,
            "comparison": "mean ethylene concentration in ppmv at overlapping AIJ Case H points",
            "model_scalar": "conservative transported mass concentration; converted to ethylene ppmv using ideal-gas molar masses",
            "not_hydrogen_validation": True,
            "not_lh2_pool_validation": True,
            "not_tke_validation": True,
            "quantitative_facility_prediction_allowed": False,
            "limitations": [
                "3-D carrier uses the existing stability-qualified effective-viscosity flow component, not matched experimental Reynolds number",
                (
                    "source is a cell-integrated Gaussian; its declared sigma is not a resolved nozzle profile"
                    if resolved_source_sigma_m < cell_size_m else
                    "coarse cells smear the 4 mm source over the declared Gaussian width"
                ),
                "steady time-mean carrier velocity does not resolve scalar/velocity covariance over the experiment's 120 s record",
                (
                    "the 4 mm nozzle velocity inferred from flow and opening diameter is not measured; its momentum is distributed over a Gaussian source kernel"
                    if config_changes.get("active_momentum_solver") == "finite_volume_projection" else
                    "the D3Q19 scalar-only carrier does not apply the inferred 4 mm nozzle-exit momentum"
                ),
                "the SGS scalar option is a Smagorinsky/Schmidt gradient-diffusion sensitivity, not a calibrated or independently validated closure",
                "the default source Gaussian width follows cell size; explicit sigma sweeps isolate source regularization but do not replace a resolved nozzle boundary",
                "molecular ethylene is treated as a passive scalar and no concentration-dependent momentum or thermal feedback is qualified",
                "the scalar mean is a 15-20 s model window rather than a 120 s experimental acquisition",
            ],
        },
        "configuration": {
            "domain_m": list(config.domain_m),
            "cells": list(config.cells),
            "cell_size_m": cell_size_m,
            "flow_steps": config.flow_steps,
            "flow_average_last_steps": config.flow_average_last_steps,
            "scalar_duration_s": duration_s,
            "scalar_time_step_s": scalar_time_step_s,
            "scalar_mean_window_s": [average_start_s, duration_s],
            "declared_scalar_diffusivity_m2_s": config.scalar_diffusivity_m2_s,
            "source_mass_rate_kg_s": source_mass_rate_kg_s,
            "integrated_released_species_mass_kg": released_mass_kg,
            "source_sigma_m": config.source_sigma_m,
            "momentum_model": config.active_momentum_solver,
            "source_momentum_velocity_m_s": list(config.source_momentum_velocity_m_s),
            "source_orifice_diameter_m": SOURCE_ORIFICE_DIAMETER_M,
            "source_exit_speed_from_reported_flow_m_s": source_exit_velocity_m_s,
            "source_exit_speed_is_inferred_not_directly_measured": True,
            "source_momentum_applied": use_binary_fv_source_momentum,
            "subgrid_model": config.subgrid_model,
            "smagorinsky_constant": config.smagorinsky_constant,
            "effective_schmidt_number": config.effective_schmidt_number,
            "scalar_turbulent_diffusivity_model": config.scalar_turbulent_diffusivity_model,
            "x_origin_in_experiment_m": X_ORIGIN_M,
            "y_origin_in_experiment_m": y_origin_m,
            "lateral_domain_override_m": lateral_domain_m,
        },
        "scores": {
            "overlap_count": int(obs.size),
            "skipped": skipped,
            "bias_ppmv": float(np.mean(residual)),
            "mae_ppmv": float(np.mean(np.abs(residual))),
            "rmse_ppmv": float(np.sqrt(np.mean(residual ** 2))),
            "measured_mean_ppmv": float(np.mean(obs)),
            "predicted_mean_ppmv": float(np.mean(pred)),
            "measured_p95_ppmv": float(np.percentile(obs, 95)),
            "predicted_p95_ppmv": float(np.percentile(pred, 95)),
            "measured_max_ppmv": float(np.max(obs)),
            "predicted_max_ppmv": float(np.max(pred)),
            "pearson_r": _pearson(obs, pred),
            "median_measured_c_rms_ppmv_descriptive_only": float(np.nanmedian(
                [float(row["measured_c_rms_ppmv"]) for row in rows
                 if row["measured_c_rms_ppmv"] is not None]
            )),
            "measured_K_mean": float(np.mean(measured_k)),
            "predicted_K_mean": float(np.mean(predicted_k)),
        },
        "numerics": {
            "source_kernel_fwhm_cell_count_xyz": result.diagnostics.get(
                "source_kernel_fwhm_cell_count_xyz"
            ),
            "source_kernel_underresolved": result.diagnostics.get(
                "source_kernel_underresolved"
            ),
            "maximum_relative_species_mass_residual": result.diagnostics.get("maximum_relative_mass_residual"),
            "maximum_relative_sensible_enthalpy_residual": result.diagnostics.get("maximum_relative_energy_residual"),
            "maximum_absolute_divergence_metric_s_1": result.diagnostics.get("maximum_absolute_divergence_s_1"),
            "maximum_absolute_raw_face_divergence_s_1": result.diagnostics.get("maximum_absolute_raw_face_divergence_s_1"),
            "maximum_variable_density_continuity_residual_kg_m3_s": result.diagnostics.get("maximum_variable_density_continuity_residual_kg_m3_s"),
            "maximum_variable_density_global_residual_kg_s": result.diagnostics.get("maximum_variable_density_global_residual_kg_s"),
            "maximum_momentum_balance_residual_kg_m_s": result.diagnostics.get("maximum_momentum_balance_residual_kg_m_s"),
            "flow_converged": result.diagnostics.get("flow_converged"),
            "flow_statistical_residual": result.diagnostics.get("flow_statistical_residual"),
            "maximum_advective_number": result.diagnostics.get("advective_number"),
            "maximum_diffusion_number": result.diagnostics.get("diffusion_number"),
            "thermodynamic_corrector_nonconverged_steps": result.diagnostics.get("thermodynamic_corrector_nonconverged_steps"),
            "applied_scalar_diffusivity_mean_m2_s": result.diagnostics.get("applied_scalar_diffusivity_mean_m2_s"),
            "applied_thermal_diffusivity_mean_m2_s": result.diagnostics.get("applied_thermal_diffusivity_mean_m2_s"),
            "thermodynamic_corrector_maximum_iterations_used": result.diagnostics.get("thermodynamic_corrector_maximum_iterations_used"),
            "active_scalar_applicability_passed": result.diagnostics.get("active_scalar_applicability_passed"),
            "active_momentum_solver": result.diagnostics.get("active_momentum_solver"),
            "variable_density_linear_solver": result.diagnostics.get("variable_density_linear_solver"),
        },
        "runtime": {
            "wall_clock_s": wall_clock_s,
            "solver_reported_total_s": result.diagnostics.get("total_runtime_s"),
            "model_time_s": duration_s,
        },
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "aij-case-h-tracer-score.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with (output_directory / "aij-case-h-tracer-point-scores.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    figure, axis = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    axis.scatter(obs, pred, c=[row["z_m"] for row in rows], cmap="viridis", s=18, alpha=0.65)
    low = min(float(np.min(obs)), float(np.min(pred)))
    high = max(float(np.max(obs)), float(np.max(pred)))
    axis.plot([low, high], [low, high], "k--", lw=1)
    axis.set(xlabel="Measured mean ethylene (ppmv)", ylabel="Predicted mean ethylene (ppmv)")
    axis.set_title(
        "AIJ Case H — 3-D resolved passive-tracer component\n"
        f"RMSE={record['scores']['rmse_ppmv']:.1f} ppmv; r={record['scores']['pearson_r']:.3f}; no fit"
    )
    axis.grid(alpha=0.2)
    figure.savefig(output_directory / "aij-case-h-tracer-parity.png", dpi=180)
    plt.close(figure)
    record["output_hashes"] = {
        name: _sha256(output_directory / name)
        for name in ("aij-case-h-tracer-point-scores.csv", "aij-case-h-tracer-parity.png")
    }
    (output_directory / "aij-case-h-tracer-score.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--duration-s", type=float, default=20.0)
    parser.add_argument("--cell-size-m", type=float, default=0.025)
    parser.add_argument("--scalar-time-step-s", type=float, default=0.002)
    parser.add_argument(
        "--source-sigma-m", type=float, default=None,
        help="Gaussian source sigma in metres (default: one cell); sensitivity only, not a resolved nozzle",
    )
    parser.add_argument("--average-start-s", type=float, default=15.0)
    parser.add_argument(
        "--lateral-domain-m", type=float, default=None,
        help="optional wider lateral domain for a boundary-sensitivity audit; no concentration fitting",
    )
    parser.add_argument(
        "--y-origin-m", type=float, default=Y_ORIGIN_M,
        help="experiment-coordinate y origin used with --lateral-domain-m",
    )
    parser.add_argument(
        "--use-sgs-scalar-mixing", action="store_true",
        help="enable uncalibrated Smagorinsky/Schmidt scalar-gradient-diffusion sensitivity",
    )
    parser.add_argument(
        "--binary-fv-source-momentum", action="store_true",
        help="use generic ethylene-air EOS, low-Mach finite-volume projection, and inferred 4 mm nozzle momentum",
    )
    parser.add_argument(
        "--binary-fv-neutral-source", action="store_true",
        help="use the same ethylene-air finite-volume path but omit nozzle momentum as a matched control",
    )
    args = parser.parse_args()
    record = run(
        output_directory=args.output_directory,
        duration_s=args.duration_s,
        cell_size_m=args.cell_size_m,
        scalar_time_step_s=args.scalar_time_step_s,
        source_sigma_m=args.source_sigma_m,
        average_start_s=args.average_start_s,
        use_sgs_scalar_mixing=args.use_sgs_scalar_mixing,
        use_binary_fv_source_momentum=args.binary_fv_source_momentum,
        use_binary_fv_neutral_source=args.binary_fv_neutral_source,
        lateral_domain_m=args.lateral_domain_m,
        y_origin_m=args.y_origin_m,
    )
    print(json.dumps({"status": record["status"], "scores": record["scores"], "runtime": record["runtime"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
