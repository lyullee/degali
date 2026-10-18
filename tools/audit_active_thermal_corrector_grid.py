"""Three-level joint grid/time audit for the active Tier-3 obstacle path.

The physical domain, obstacle, source, wind, flow spin-up and transient are
fixed.  Cell size and time step are refined together at constant lattice inlet
speed, which is the appropriate fixed-Mach refinement path for this D3Q19
prototype.  This is development evidence, not external validation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.fast_obstacle_transport import (
    HYDROGEN_MOLAR_MASS_KG_MOL,
    R_GAS_J_MOL_K,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"
DOMAIN_M = (4.8, 3.6, 2.4)
SPACINGS_M = (0.30, 0.20, 2.0 / 15.0)  # exact 1.5 refinement ratio
WIND_M_S = 0.50
LATTICE_INLET_SPEED = 0.015
FLOW_SPINUP_S = 3.6
DURATION_S = 1.8
OBSTACLE = AxisAlignedCuboid(
    1.2, 2.4, 1.2, 2.4, 0.0, 1.2, "grid-aligned finite block",
)


def _config(
    spacing: float, scalar_advection_scheme: str = "muscl_van_leer",
) -> ResolvedObstacleFlow3DConfig:
    cells = tuple(int(round(length / spacing)) for length in DOMAIN_M)
    if any(not math.isclose(count * spacing, length, abs_tol=1.0e-12)
           for count, length in zip(cells, DOMAIN_M)):
        raise RuntimeError("grid does not tile the physical domain")
    dt = LATTICE_INLET_SPEED * spacing / WIND_M_S
    flow_steps = round(FLOW_SPINUP_S / dt)
    transport_steps = round(DURATION_S / dt)
    if not math.isclose(flow_steps * dt, FLOW_SPINUP_S, abs_tol=1.0e-12):
        raise RuntimeError("flow spin-up is not aligned with the selected grid")
    if not math.isclose(transport_steps * dt, DURATION_S, abs_tol=1.0e-12):
        raise RuntimeError("transient is not aligned with the selected grid")
    return ResolvedObstacleFlow3DConfig(
        domain_m=DOMAIN_M,
        cells=cells,
        wind_speed_m_s=WIND_M_S,
        kinematic_viscosity_m2_s=0.20,
        flow_steps=flow_steps,
        minimum_flow_steps=flow_steps,
        lattice_inlet_speed=LATTICE_INLET_SPEED,
        duration_s=DURATION_S,
        scalar_time_step_s=dt,
        scalar_diffusivity_m2_s=0.030,
        thermal_diffusivity_m2_s=0.030,
        source_schedule=SourceRateSchedule(
            (0.0, 0.36, DURATION_S), (0.015, 0.0, 0.0),
        ),
        source_position_m=(0.60, 1.80, 0.30),
        source_sigma_m=0.30,
        obstacles=(OBSTACLE,),
        scalar_advection_scheme=scalar_advection_scheme,
        smagorinsky_constant=0.10,
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        maximum_active_scalar_density_departure=0.10,
        require_thermodynamic_corrector_convergence=True,
        store_scalar_fields=False,
    )


def _dose(time_s: np.ndarray, trace: np.ndarray) -> float:
    return float(np.trapezoid(trace, time_s))


def _arrival(time_s: np.ndarray, trace: np.ndarray) -> float | None:
    indices = np.flatnonzero(trace >= 4.0)
    return None if indices.size == 0 else float(time_s[indices[0]])


def _case_record(result, config: ResolvedObstacleFlow3DConfig) -> dict[str, object]:
    maximum_vol = result.maximum_concentration_vol_percent
    spacing = config.domain_m[0] / config.cells[0]
    return {
        "cell_size_m": spacing,
        "cells": list(config.cells),
        "time_step_s": config.scalar_time_step_s,
        "flow_steps": config.flow_steps,
        "transport_steps": int(round(config.duration_s / config.scalar_time_step_s)),
        "runtime_s": result.diagnostics["total_runtime_s"],
        "maximum_h2_vol_percent": float(np.nanmax(maximum_vol)),
        "four_vol_percent_volume_m3": float(
            np.count_nonzero(maximum_vol >= 4.0) * spacing**3
        ),
        "minimum_temperature_k": float(np.nanmin(result.temperature_k)),
        "maximum_lattice_mach_number": result.diagnostics[
            "maximum_lattice_mach_number"
        ],
        "maximum_relative_density_departure": result.diagnostics[
            "maximum_relative_mixture_density_departure"
        ],
        "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        "maximum_energy_residual_j": result.maximum_energy_residual_j,
        "maximum_continuity_residual_kg_m3_s": result.diagnostics[
            "maximum_variable_density_continuity_residual_kg_m3_s"
        ],
        "maximum_global_continuity_residual_kg_s": result.diagnostics[
            "maximum_variable_density_global_residual_kg_s"
        ],
        "maximum_momentum_balance_residual_kg_m_s": result.diagnostics[
            "maximum_momentum_balance_residual_kg_m_s"
        ],
        "corrector_nonconverged_steps": result.diagnostics[
            "thermodynamic_corrector_nonconverged_steps"
        ],
        "variable_density_linear_solver": result.diagnostics[
            "variable_density_linear_solver"
        ],
        "receptors": {
            label: {
                "peak_h2_vol_percent": float(np.max(trace)),
                "dose_h2_vol_percent_s": _dose(result.time_s, trace),
                "arrival_time_4vol_s": _arrival(result.time_s, trace),
            }
            for label, trace in result.receptor_traces_vol_percent.items()
        },
    }


def _overlap_average_weights(
    source_coordinates: np.ndarray, target_coordinates: np.ndarray,
) -> np.ndarray:
    """Map uniform source-cell averages to uniform target-cell averages."""
    source_spacing = float(source_coordinates[1] - source_coordinates[0])
    target_spacing = float(target_coordinates[1] - target_coordinates[0])
    source_lower = source_coordinates - 0.5 * source_spacing
    source_upper = source_coordinates + 0.5 * source_spacing
    target_lower = target_coordinates - 0.5 * target_spacing
    target_upper = target_coordinates + 0.5 * target_spacing
    overlap = np.maximum(
        0.0,
        np.minimum(target_upper[:, None], source_upper[None, :])
        - np.maximum(target_lower[:, None], source_lower[None, :]),
    )
    overlap[
        overlap
        <= 1.0e-12 * min(abs(source_spacing), abs(target_spacing))
    ] = 0.0
    weights = overlap / target_spacing
    if not np.allclose(np.sum(weights, axis=1), 1.0, atol=1.0e-12):
        raise RuntimeError("source and target finite-volume domains do not coincide")
    return weights


def _conservative_remap_cell_averages(
    field: np.ndarray,
    source_coordinates: tuple[np.ndarray, np.ndarray, np.ndarray],
    target_coordinates: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    """Conservatively remap a z-y-x finite-volume field, retaining solid NaNs."""
    source_x, source_y, source_z = source_coordinates
    target_x, target_y, target_z = target_coordinates
    wx = _overlap_average_weights(source_x, target_x)
    wy = _overlap_average_weights(source_y, target_y)
    wz = _overlap_average_weights(source_z, target_z)
    finite = np.isfinite(field)
    numerator = np.einsum(
        "ai,bj,ck,ijk->abc", wz, wy, wx,
        np.where(finite, field, 0.0), optimize=True,
    )
    denominator = np.einsum(
        "ai,bj,ck,ijk->abc", wz, wy, wx,
        finite.astype(float), optimize=True,
    )
    remapped = np.full(numerator.shape, np.nan)
    np.divide(numerator, denominator, out=remapped, where=denominator > 0.0)
    return remapped


def _sample_on_coarse_grid(results, field_getter):
    coarse = results[0]
    target = (coarse.x_m, coarse.y_m, coarse.z_m)
    samples = []
    for result in results:
        remapped = _conservative_remap_cell_averages(
            field_getter(result),
            (result.x_m, result.y_m, result.z_m), target,
        )
        samples.append(remapped.ravel())
    valid = np.logical_and.reduce([np.isfinite(sample) for sample in samples])
    return [sample[valid] for sample in samples], int(np.count_nonzero(valid))


def _field_metric(samples, refinement_ratio: float) -> dict[str, float | None]:
    coarse_medium = float(np.sqrt(np.mean((samples[0] - samples[1])**2)))
    medium_fine = float(np.sqrt(np.mean((samples[1] - samples[2])**2)))
    order = None
    if coarse_medium > 0.0 and medium_fine > 0.0:
        order = float(math.log(coarse_medium / medium_fine) / math.log(
            refinement_ratio
        ))
    return {
        "coarse_medium_rms": coarse_medium,
        "medium_fine_rms": medium_fine,
        "observed_order": order,
        "successive_difference_decreased": medium_fine < coarse_medium,
    }


def _relative_change(left: float, right: float, floor: float) -> float:
    return abs(right - left) / max(abs(left), abs(right), floor)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scheme", choices=("upwind", "muscl_van_leer", "muscl_ssprk2"),
        default="muscl_van_leer",
    )
    parser.add_argument("--effective-schmidt", type=float, default=None)
    parser.add_argument("--effective-prandtl", type=float, default=None)
    parser.add_argument(
        "--common-time-maximum", action="store_true",
        help="compare maximum fields only at the 0.036 s times shared by all grids",
    )
    parser.add_argument(
        "--source-time-integration", choices=("prestep", "strang"),
        default="prestep",
    )
    parser.add_argument(
        "--active-momentum-solver",
        choices=("d3q19_bgk", "finite_volume_projection"),
        default="d3q19_bgk",
    )
    arguments = parser.parse_args(argv)
    for name, value in (
        ("effective Schmidt", arguments.effective_schmidt),
        ("effective Prandtl", arguments.effective_prandtl),
    ):
        if value is not None and (not math.isfinite(value) or value <= 0.0):
            parser.error(f"{name} number must be finite and positive")
    receptors = (
        FastObstacleReceptor3D("windward", (1.05, 1.80, 0.60)),
        FastObstacleReceptor3D("roof", (1.80, 1.80, 1.35)),
        FastObstacleReceptor3D("side", (1.80, 1.05, 0.60)),
        FastObstacleReceptor3D("lee", (2.70, 1.80, 0.60)),
    )
    configs = []
    for spacing in SPACINGS_M:
        base = _config(spacing, arguments.scheme)
        configs.append(replace(
            base,
            effective_schmidt_number=arguments.effective_schmidt,
            effective_prandtl_number=arguments.effective_prandtl,
            source_time_integration=arguments.source_time_integration,
            active_momentum_solver=arguments.active_momentum_solver,
            store_scalar_fields=arguments.common_time_maximum,
        ))
    results = [
        solve_resolved_obstacle_flow_3d(config, receptors=receptors)
        for config in configs
    ]
    cases = [_case_record(result, config) for result, config in zip(results, configs)]
    refinement_ratio = SPACINGS_M[0] / SPACINGS_M[1]
    common_time_maximum_fields = None
    common_time_s = None
    if arguments.common_time_maximum:
        common_time_s = np.arange(0.0, DURATION_S + 1.0e-12, 0.036)
        common_time_maximum_fields = {}
        for result in results:
            if (
                result.concentration_fields_kg_m3 is None
                or result.temperature_fields_k is None
            ):
                raise RuntimeError("common-time maximum requires stored active fields")
            result_dt = float(result.time_s[1] - result.time_s[0])
            indices = np.rint(common_time_s / result_dt).astype(int)
            if np.any(indices >= result.time_s.size) or not np.allclose(
                result.time_s[indices], common_time_s,
                rtol=1.0e-10, atol=1.0e-11,
            ):
                raise RuntimeError("declared common output times do not align")
            concentration = result.concentration_fields_kg_m3[indices]
            temperature = result.temperature_fields_k[indices]
            vol_percent = np.clip(
                concentration / HYDROGEN_MOLAR_MASS_KG_MOL
                * R_GAS_J_MOL_K * temperature / result.ambient_pressure_pa
                * 100.0,
                0.0, 100.0,
            )
            common_time_maximum_fields[id(result)] = np.nanmax(
                vol_percent, axis=0,
            )
        maximum_getter = lambda result: common_time_maximum_fields[id(result)]
    else:
        maximum_getter = lambda result: result.maximum_concentration_vol_percent
    maximum_samples, point_count = _sample_on_coarse_grid(
        results, maximum_getter,
    )
    final_samples, final_point_count = _sample_on_coarse_grid(
        results, lambda result: result.concentration_vol_percent,
    )
    temperature_samples, temperature_point_count = _sample_on_coarse_grid(
        results, lambda result: result.temperature_k,
    )
    field_metrics = {
        "maximum_h2_vol_percent": _field_metric(
            maximum_samples, refinement_ratio,
        ),
        "final_h2_vol_percent": _field_metric(final_samples, refinement_ratio),
        "final_temperature_k": _field_metric(
            temperature_samples, refinement_ratio,
        ),
    }

    medium, fine = cases[-2], cases[-1]
    finest_pair = {
        "maximum_h2_relative_change": _relative_change(
            float(medium["maximum_h2_vol_percent"]),
            float(fine["maximum_h2_vol_percent"]), 0.01,
        ),
        "four_vol_percent_volume_relative_change": _relative_change(
            float(medium["four_vol_percent_volume_m3"]),
            float(fine["four_vol_percent_volume_m3"]), SPACINGS_M[-1]**3,
        ),
        "receptors": {},
    }
    for receptor in receptors:
        left = medium["receptors"][receptor.label]
        right = fine["receptors"][receptor.label]
        finest_pair["receptors"][receptor.label] = {
            "peak_relative_change": _relative_change(
                float(left["peak_h2_vol_percent"]),
                float(right["peak_h2_vol_percent"]), 0.01,
            ),
            "dose_relative_change": _relative_change(
                float(left["dose_h2_vol_percent_s"]),
                float(right["dose_h2_vol_percent_s"]), 0.01,
            ),
            "arrival_absolute_change_s": (
                None
                if left["arrival_time_4vol_s"] is None
                or right["arrival_time_4vol_s"] is None
                else abs(
                    float(left["arrival_time_4vol_s"])
                    - float(right["arrival_time_4vol_s"])
                )
            ),
            "arrival_status_changed": bool(
                (left["arrival_time_4vol_s"] is None)
                != (right["arrival_time_4vol_s"] is None)
            ),
        }

    all_correctors_converged = all(
        case["corrector_nonconverged_steps"] == 0 for case in cases
    )
    continuity_gates_passed = all(
        float(case["maximum_continuity_residual_kg_m3_s"]) < 1.0e-8
        and abs(float(case["maximum_global_continuity_residual_kg_s"])) < 1.0e-9
        for case in cases
    )
    monotonic = all(
        bool(metric["successive_difference_decreased"])
        and metric["observed_order"] is not None
        and float(metric["observed_order"]) > 0.0
        for metric in field_metrics.values()
    )
    gate_reasons = []
    if not all_correctors_converged:
        gate_reasons.append("one or more thermodynamic corrector steps did not converge")
    if not continuity_gates_passed:
        gate_reasons.append("one or more continuity residual gates failed")
    if not monotonic:
        gate_reasons.append("successive full-field differences are not monotonically decreasing")
    if finest_pair["maximum_h2_relative_change"] > 0.10:
        gate_reasons.append("finest-pair domain maximum changes by more than 10%")
    if finest_pair["four_vol_percent_volume_relative_change"] > 0.10:
        gate_reasons.append("finest-pair 4 vol% volume changes by more than 10%")
    for label, metric in finest_pair["receptors"].items():
        if metric["peak_relative_change"] > 0.10:
            gate_reasons.append(f"{label} peak changes by more than 10%")
        if metric["dose_relative_change"] > 0.10:
            gate_reasons.append(f"{label} dose changes by more than 10%")
        if metric["arrival_status_changed"]:
            gate_reasons.append(f"{label} 4 vol% arrival status changes")

    record = {
        "schema": "degali.active-thermal-corrector-grid-audit.v1",
        "status": (
            "accepted_component_grid"
            if not gate_reasons else "withheld_grid_refinement_gate_failed"
        ),
        "scope": (
            "joint grid/time refinement at fixed lattice Mach; development "
            "evidence only, not external concentration validation"
        ),
        "basis": {
            "domain_m": DOMAIN_M,
            "obstacle": asdict(OBSTACLE),
            "spacings_m": SPACINGS_M,
            "refinement_ratio": refinement_ratio,
            "obstacle_boundaries_aligned_on_all_grids": True,
            "same_physical_flow_spinup_s": FLOW_SPINUP_S,
            "same_transient_duration_s": DURATION_S,
            "same_source_history": True,
            "source_spatial_discretization": (
                "cell-integrated separable Gaussian over Cartesian cells"
            ),
            "same_lattice_inlet_speed": LATTICE_INLET_SPEED,
            "scalar_advection_scheme": arguments.scheme,
            "effective_schmidt_number": arguments.effective_schmidt,
            "effective_prandtl_number": arguments.effective_prandtl,
            "source_time_integration": arguments.source_time_integration,
            "active_momentum_solver": arguments.active_momentum_solver,
            "maximum_field_observation_operator": (
                "maximum on common 0.036 s output times"
                if arguments.common_time_maximum
                else "maximum on each grid native time steps"
            ),
            "spatial_field_comparison_operator": (
                "conservative finite-volume remap to the coarse control volumes"
            ),
            "common_maximum_time_s": (
                None if common_time_s is None else common_time_s.tolist()
            ),
            "configs": [asdict(config) for config in configs],
        },
        "common_interpolation_point_counts": {
            "maximum_h2": point_count,
            "final_h2": final_point_count,
            "temperature": temperature_point_count,
        },
        "cases_coarse_to_fine": cases,
        "field_difference_metrics": field_metrics,
        "finest_pair_changes": finest_pair,
        "numerical_gates": {
            "all_correctors_converged": all_correctors_converged,
            "continuity_residual_gates_passed": continuity_gates_passed,
            "formal_joint_grid_time_convergence_demonstrated": monotonic,
        },
        "gate_reasons": gate_reasons,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    suffix = {
        "muscl_van_leer": "_cell_integrated_source",
        "upwind": "_upwind_cell_integrated_source",
        "muscl_ssprk2": "_ssprk2_cell_integrated_source",
    }[arguments.scheme]
    if arguments.effective_schmidt is not None:
        suffix += f"_sc{arguments.effective_schmidt:g}".replace(".", "p")
    if arguments.effective_prandtl is not None:
        suffix += f"_pr{arguments.effective_prandtl:g}".replace(".", "p")
    if arguments.common_time_maximum:
        suffix += "_common_time"
    if arguments.source_time_integration != "prestep":
        suffix += "_" + arguments.source_time_integration
    if arguments.active_momentum_solver != "d3q19_bgk":
        suffix += "_" + arguments.active_momentum_solver
    json_path = OUTPUT / f"active_thermal_corrector_grid_audit{suffix}.json"
    figure_path = OUTPUT / f"active_thermal_corrector_grid_audit{suffix}.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    common_max = max(float(np.nanmax(
        result.maximum_concentration_vol_percent
    )) for result in results)
    for axis, result, case, label in zip(
        axes[0], results, cases, ("coarse", "medium", "fine"),
    ):
        center = int(np.argmin(np.abs(result.y_m - 1.80)))
        field = result.maximum_concentration_vol_percent[:, center, :]
        image = axis.pcolormesh(
            result.x_m, result.z_m, field, shading="nearest", cmap="magma",
            vmin=0.0, vmax=common_max,
        )
        if np.nanmax(field) >= 4.0:
            axis.contour(
                result.x_m, result.z_m, field, levels=[4.0], colors="white",
                linewidths=1.2,
            )
        obstacle = result.solid_mask[:, center, :]
        axis.contourf(
            result.x_m, result.z_m, obstacle.astype(float),
            levels=[0.5, 1.5], colors=["0.25"],
        )
        axis.set_title(
            f"{label}: dx={case['cell_size_m']:.3f} m, "
            f"dt={case['time_step_s']:.4f} s"
        )
        axis.set_xlabel("x (m)")
        axis.set_ylabel("z (m)")
        figure.colorbar(image, ax=axis, label="maximum H2 (vol%)")

    spacing_values = np.asarray([case["cell_size_m"] for case in cases])
    axes[1, 0].plot(
        spacing_values,
        [case["maximum_h2_vol_percent"] for case in cases], "o-",
        label="maximum H2",
    )
    axes[1, 0].set_ylabel("maximum H2 (vol%)")
    volume_axis = axes[1, 0].twinx()
    volume_axis.plot(
        spacing_values,
        [case["four_vol_percent_volume_m3"] for case in cases], "s--",
        color="#16a34a", label="4 vol% volume",
    )
    volume_axis.set_ylabel("4 vol% volume (m3)")
    axes[1, 0].set_title("Hazard metrics")

    for receptor in receptors:
        axes[1, 1].plot(
            spacing_values,
            [case["receptors"][receptor.label]["peak_h2_vol_percent"]
             for case in cases],
            "o-", label=receptor.label,
        )
    axes[1, 1].set_title("Receptor peak sensitivity")
    axes[1, 1].set_ylabel("peak H2 (vol%)")
    axes[1, 1].legend(fontsize=8)

    metric_names = list(field_metrics)
    axes[1, 2].bar(
        np.arange(3) - 0.18,
        [field_metrics[name]["coarse_medium_rms"] for name in metric_names],
        width=0.36, label="coarse-medium",
    )
    axes[1, 2].bar(
        np.arange(3) + 0.18,
        [field_metrics[name]["medium_fine_rms"] for name in metric_names],
        width=0.36, label="medium-fine",
    )
    axes[1, 2].set_xticks(np.arange(3), ("max H2", "final H2", "final T"))
    axes[1, 2].set_yscale("log")
    axes[1, 2].set_title("Common-point RMS differences")
    axes[1, 2].set_ylabel("native field units")
    axes[1, 2].legend(fontsize=8)
    for axis in axes[1, :2]:
        axis.set_xlabel("cell size (m), finer →")
        axis.invert_xaxis()
        axis.grid(alpha=0.25)
    axes[1, 2].set_xlabel("field")
    axes[1, 2].grid(axis="y", alpha=0.25)
    figure.suptitle(
        "Active Tier-3 joint grid/time audit — "
        + arguments.scheme + " — "
        + ("DEVELOPMENT" if not gate_reasons else "GATE FAILED")
    )
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "cases_coarse_to_fine": cases,
        "field_difference_metrics": field_metrics,
        "finest_pair_changes": finest_pair,
        "gate_reasons": gate_reasons,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
