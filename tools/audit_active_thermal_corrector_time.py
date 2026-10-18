"""Three-level physical-time sensitivity audit for the coupled Tier-3 path.

All cases use the same grid, physical flow spin-up horizon, release history,
and 3 s transient.  Changing the lattice inlet speed changes both the LBM and
scalar time step while the physical wind speed remains fixed.
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
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"
DX_M = 0.20
WIND_M_S = 0.50
DURATION_S = 3.0
FLOW_SPINUP_S = 7.2
LATTICE_SPEEDS = (0.030, 0.015, 0.0075)  # coarse, medium, fine


def _config(lattice_speed: float) -> ResolvedObstacleFlow3DConfig:
    dt = lattice_speed * DX_M / WIND_M_S
    transport_steps = round(DURATION_S / dt)
    flow_steps = round(FLOW_SPINUP_S / dt)
    if not math.isclose(transport_steps * dt, DURATION_S, abs_tol=1.0e-12):
        raise RuntimeError("duration is not aligned with the selected time step")
    if not math.isclose(flow_steps * dt, FLOW_SPINUP_S, abs_tol=1.0e-12):
        raise RuntimeError("flow spin-up is not aligned with the selected time step")
    return ResolvedObstacleFlow3DConfig(
        domain_m=(4.0, 2.4, 2.0),
        cells=(20, 12, 10),
        wind_speed_m_s=WIND_M_S,
        kinematic_viscosity_m2_s=0.20,
        flow_steps=flow_steps,
        minimum_flow_steps=flow_steps,
        lattice_inlet_speed=lattice_speed,
        duration_s=DURATION_S,
        scalar_time_step_s=dt,
        scalar_diffusivity_m2_s=0.030,
        thermal_diffusivity_m2_s=0.030,
        source_schedule=SourceRateSchedule(
            (0.0, 0.48, DURATION_S), (0.015, 0.0, 0.0),
        ),
        source_position_m=(1.00, 1.20, 0.30),
        source_sigma_m=0.25,
        obstacles=(
            AxisAlignedCuboid(1.40, 2.00, 0.80, 1.60, 0.0, 1.00, "block"),
        ),
        smagorinsky_constant=0.10,
        thermal_buoyancy_coupling=True,
        source_temperature_k=80.0,
        maximum_active_scalar_density_departure=0.10,
        require_thermodynamic_corrector_convergence=True,
        store_scalar_fields=False,
    )


def _rms(left: np.ndarray, right: np.ndarray, mask: np.ndarray) -> float:
    delta = np.asarray(left)[mask] - np.asarray(right)[mask]
    return float(np.sqrt(np.mean(delta**2)))


def _observed_order(coarse_medium: float, medium_fine: float) -> float | None:
    if coarse_medium <= 0.0 or medium_fine <= 0.0:
        return None
    return float(math.log(coarse_medium / medium_fine) / math.log(2.0))


def _record(result, config: ResolvedObstacleFlow3DConfig) -> dict[str, object]:
    dx = float(result.x_m[1] - result.x_m[0])
    maximum_vol = result.maximum_concentration_vol_percent
    return {
        "time_step_s": config.scalar_time_step_s,
        "flow_steps": config.flow_steps,
        "transport_steps": int(round(config.duration_s / config.scalar_time_step_s)),
        "runtime_s": result.diagnostics["total_runtime_s"],
        "maximum_h2_vol_percent": float(np.nanmax(maximum_vol)),
        "four_vol_percent_volume_m3": float(
            np.count_nonzero(maximum_vol >= 4.0) * dx**3
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
        "corrector_maximum_iterations": result.diagnostics[
            "thermodynamic_corrector_maximum_iterations_used"
        ],
        "corrector_mean_iterations": result.diagnostics[
            "thermodynamic_corrector_mean_iterations_used"
        ],
        "corrector_nonconverged_steps": result.diagnostics[
            "thermodynamic_corrector_nonconverged_steps"
        ],
        "variable_density_linear_solver": result.diagnostics[
            "variable_density_linear_solver"
        ],
        "receptor_peak_h2_vol_percent": {
            label: float(np.max(trace))
            for label, trace in result.receptor_traces_vol_percent.items()
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--effective-schmidt", type=float, default=None)
    parser.add_argument("--effective-prandtl", type=float, default=None)
    parser.add_argument(
        "--scheme",
        choices=("upwind", "muscl_van_leer", "muscl_ssprk2"),
        default="muscl_van_leer",
    )
    parser.add_argument(
        "--source-time-integration", choices=("prestep", "strang"),
        default="prestep",
    )
    parser.add_argument(
        "--gravity", type=float, default=9.80665,
        help="gravity used by the active-scalar buoyancy closure (m/s^2)",
    )
    parser.add_argument(
        "--disable-mass-flux-correction", action="store_true",
        help=(
            "diagnostic only: omit variable-density continuity projection and "
            "corrected-velocity reinjection into D3Q19"
        ),
    )
    parser.add_argument(
        "--active-provisional-face-projection",
        choices=(
            "incompressible_then_variable_density", "variable_density_only",
        ),
        default="incompressible_then_variable_density",
    )
    parser.add_argument(
        "--pressure-correction-momentum-coupling",
        choices=(
            "equilibrium_replacement", "deferred_guo_force", "transport_only",
        ),
        default="equilibrium_replacement",
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
    if not math.isfinite(arguments.gravity) or arguments.gravity < 0.0:
        parser.error("gravity must be finite and non-negative")
    receptors = (
        FastObstacleReceptor3D("windward", (1.20, 1.20, 0.50)),
        FastObstacleReceptor3D("roof", (1.70, 1.20, 1.20)),
        FastObstacleReceptor3D("side", (1.70, 0.60, 0.50)),
        FastObstacleReceptor3D("lee", (2.40, 1.20, 0.50)),
    )
    base_configs = (
        [
            replace(
                _config(0.015),
                scalar_time_step_s=speed * DX_M / WIND_M_S,
            )
            for speed in LATTICE_SPEEDS
        ]
        if arguments.active_momentum_solver == "finite_volume_projection"
        else [_config(speed) for speed in LATTICE_SPEEDS]
    )
    configs = [replace(
        config,
        effective_schmidt_number=arguments.effective_schmidt,
        effective_prandtl_number=arguments.effective_prandtl,
        scalar_advection_scheme=arguments.scheme,
        source_time_integration=arguments.source_time_integration,
        gravity_m_s2=arguments.gravity,
        low_mach_mass_flux_correction=(
            not arguments.disable_mass_flux_correction
        ),
        require_thermodynamic_corrector_convergence=(
            not arguments.disable_mass_flux_correction
        ),
        active_provisional_face_projection=(
            arguments.active_provisional_face_projection
        ),
        pressure_correction_momentum_coupling=(
            arguments.pressure_correction_momentum_coupling
        ),
        active_momentum_solver=arguments.active_momentum_solver,
    ) for config in base_configs]
    results = [
        solve_resolved_obstacle_flow_3d(config, receptors=receptors)
        for config in configs
    ]
    coarse, medium, fine = results
    fluid = ~fine.solid_mask
    metrics: dict[str, dict[str, float | None]] = {}
    for name, values in (
        ("maximum_h2_vol_percent", [
            result.maximum_concentration_vol_percent for result in results
        ]),
        ("final_h2_vol_percent", [
            result.concentration_vol_percent for result in results
        ]),
        ("final_temperature_k", [result.temperature_k for result in results]),
    ):
        coarse_medium = _rms(values[0], values[1], fluid)
        medium_fine = _rms(values[1], values[2], fluid)
        metrics[name] = {
            "coarse_medium_rms": coarse_medium,
            "medium_fine_rms": medium_fine,
            "observed_order": _observed_order(coarse_medium, medium_fine),
        }

    cases = [_record(result, config) for result, config in zip(results, configs)]
    all_correctors_converged = all(
        case["corrector_nonconverged_steps"] == 0 for case in cases
    )
    residual_gates_passed = all(
        float(case["maximum_continuity_residual_kg_m3_s"]) < 1.0e-8
        and abs(float(case["maximum_global_continuity_residual_kg_s"])) < 1.0e-9
        for case in cases
    )
    successive_field_differences_decrease = all(
        float(metric["medium_fine_rms"])
        < float(metric["coarse_medium_rms"])
        for metric in metrics.values()
    )
    record = {
        "schema": "degali.active-thermal-corrector-time-audit.v1",
        "status": (
            "withheld_time_refinement_not_monotonic"
            if not successive_field_differences_decrease
            else "accepted_component_time"
        ),
        "basis": {
            "same_grid": True,
            "same_physical_flow_spinup_s": FLOW_SPINUP_S,
            "same_transient_duration_s": DURATION_S,
            "same_source_history": True,
            "refinement_ratio": 2.0,
            "effective_schmidt_number": arguments.effective_schmidt,
            "effective_prandtl_number": arguments.effective_prandtl,
            "scalar_advection_scheme": arguments.scheme,
            "source_time_integration": arguments.source_time_integration,
            "gravity_m_s2": arguments.gravity,
            "low_mach_mass_flux_correction": (
                not arguments.disable_mass_flux_correction
            ),
            "active_provisional_face_projection": (
                arguments.active_provisional_face_projection
            ),
            "pressure_correction_momentum_coupling": (
                arguments.pressure_correction_momentum_coupling
            ),
            "active_momentum_solver": arguments.active_momentum_solver,
            "configs": [asdict(config) for config in configs],
        },
        "numerical_gates": {
            "all_correctors_converged": all_correctors_converged,
            "continuity_residual_gates_passed": residual_gates_passed,
            "successive_field_differences_decrease": (
                successive_field_differences_decrease
            ),
            "formal_time_convergence_demonstrated": bool(
                successive_field_differences_decrease
                and all(
                    metric["observed_order"] is not None
                    and float(metric["observed_order"]) > 0.0
                    for metric in metrics.values()
                )
            ),
        },
        "cases_coarse_to_fine": cases,
        "field_difference_metrics": metrics,
        "interpretation": (
            "The audit measures time-step sensitivity of the current coupled "
            "implementation. It is not external validation and an observed "
            "order alone cannot qualify the active-scalar thermodynamics. "
            "A non-decreasing successive field norm withholds a formal "
            "time-convergence claim even when the absolute differences are small."
        ),
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    suffix = ""
    if arguments.effective_schmidt is not None:
        suffix += f"_sc{arguments.effective_schmidt:g}".replace(".", "p")
    if arguments.effective_prandtl is not None:
        suffix += f"_pr{arguments.effective_prandtl:g}".replace(".", "p")
    if arguments.scheme != "muscl_van_leer":
        suffix += "_" + arguments.scheme
    if arguments.source_time_integration != "prestep":
        suffix += "_" + arguments.source_time_integration
    if not math.isclose(arguments.gravity, 9.80665, abs_tol=1.0e-12):
        suffix += f"_g{arguments.gravity:g}".replace(".", "p")
    if arguments.disable_mass_flux_correction:
        suffix += "_no_mass_flux"
    if arguments.active_provisional_face_projection == "variable_density_only":
        suffix += "_variable_density_only"
    if arguments.pressure_correction_momentum_coupling != "equilibrium_replacement":
        suffix += "_" + arguments.pressure_correction_momentum_coupling
    if arguments.active_momentum_solver != "d3q19_bgk":
        suffix += "_" + arguments.active_momentum_solver
    json_path = OUTPUT / f"active_thermal_corrector_time_audit{suffix}.json"
    figure_path = OUTPUT / f"active_thermal_corrector_time_audit{suffix}.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    center = int(np.argmin(np.abs(fine.y_m - 1.20)))
    x, z = fine.x_m, fine.z_m
    figure, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    common_max = max(float(np.nanmax(
        result.maximum_concentration_vol_percent[:, center, :]
    )) for result in results)
    for axis, result, case, label in zip(
        axes[0], results, cases, ("coarse", "medium", "fine"),
    ):
        field = result.maximum_concentration_vol_percent[:, center, :]
        image = axis.pcolormesh(
            x, z, field, shading="nearest", cmap="magma", vmin=0.0,
            vmax=common_max,
        )
        axis.contour(x, z, field, levels=[4.0], colors="white", linewidths=1.2)
        obstacle = result.solid_mask[:, center, :]
        axis.contourf(x, z, obstacle.astype(float), levels=[0.5, 1.5], colors=["0.25"])
        axis.set_title(f"{label}: dt={case['time_step_s']:.3f} s")
        axis.set_xlabel("x (m)")
        axis.set_ylabel("z (m)")
        figure.colorbar(image, ax=axis, label="maximum H2 (vol%)")

    time_steps = np.asarray([case["time_step_s"] for case in cases])
    axes[1, 0].loglog(
        time_steps,
        [case["maximum_h2_vol_percent"] for case in cases],
        "o-",
    )
    axes[1, 0].invert_xaxis()
    axes[1, 0].set_title("Domain maximum")
    axes[1, 0].set_xlabel("time step (s)")
    axes[1, 0].set_ylabel("maximum H2 (vol%)")
    for label in receptors:
        axes[1, 1].plot(
            time_steps,
            [case["receptor_peak_h2_vol_percent"][label.label] for case in cases],
            "o-", label=label.label,
        )
    axes[1, 1].invert_xaxis()
    axes[1, 1].set_title("Receptor peak sensitivity")
    axes[1, 1].set_xlabel("time step (s)")
    axes[1, 1].set_ylabel("peak H2 (vol%)")
    axes[1, 1].legend(fontsize=8)
    metric_names = list(metrics)
    axes[1, 2].bar(
        np.arange(len(metric_names)) - 0.18,
        [metrics[name]["coarse_medium_rms"] for name in metric_names],
        width=0.36, label="coarse-medium",
    )
    axes[1, 2].bar(
        np.arange(len(metric_names)) + 0.18,
        [metrics[name]["medium_fine_rms"] for name in metric_names],
        width=0.36, label="medium-fine",
    )
    axes[1, 2].set_xticks(np.arange(len(metric_names)), (
        "max H2", "final H2", "final T",
    ))
    axes[1, 2].set_yscale("log")
    axes[1, 2].set_title("Successive-field RMS differences")
    axes[1, 2].set_ylabel("native field units")
    axes[1, 2].legend(fontsize=8)
    figure.suptitle(
        "Coupled Tier-3 time-step sensitivity — " + record["status"]
    )
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "numerical_gates": record["numerical_gates"],
        "cases_coarse_to_fine": cases,
        "field_difference_metrics": metrics,
    }, indent=2))
    return 0 if all_correctors_converged and residual_gates_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
