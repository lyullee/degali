"""Compare neutral and active-thermal resolved obstacle transport.

This is a development witness, not validation evidence.  It demonstrates
that the same local source and obstacle field can be advanced with or without
conservative sensible enthalpy and density-derived buoyancy feedback.
"""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.field_observation import (
    DRY_AIR_MOLAR_MASS_KG_MOL,
    HYDROGEN_MOLAR_MASS_KG_MOL,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"


def _config() -> ResolvedObstacleFlow3DConfig:
    dx = 0.20
    wind = 0.5
    flow_dt = 0.030 * dx / wind
    duration = 250.0 * flow_dt
    return ResolvedObstacleFlow3DConfig(
        domain_m=(4.0, 2.4, 2.0),
        cells=(20, 12, 10),
        wind_speed_m_s=wind,
        kinematic_viscosity_m2_s=0.20,
        flow_steps=600,
        minimum_flow_steps=600,
        lattice_inlet_speed=0.030,
        duration_s=duration,
        scalar_time_step_s=flow_dt,
        scalar_diffusivity_m2_s=0.030,
        thermal_diffusivity_m2_s=0.030,
        source_schedule=SourceRateSchedule(
            (0.0, 40.0 * flow_dt, duration),
            (0.015, 0.0, 0.0),
        ),
        source_position_m=(1.00, 1.20, 0.30),
        source_sigma_m=0.25,
        obstacles=(
            AxisAlignedCuboid(1.40, 2.00, 0.80, 1.60, 0.0, 1.00, "block"),
        ),
        smagorinsky_constant=0.10,
        source_temperature_k=80.0,
        maximum_active_scalar_density_departure=0.10,
        store_scalar_fields=True,
    )


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _neutral_vol_percent(result, concentration: np.ndarray) -> np.ndarray:
    hydrogen_moles = concentration / HYDROGEN_MOLAR_MASS_KG_MOL
    air_moles = result.diagnostics["ambient_air_density_kg_m3"] / (
        DRY_AIR_MOLAR_MASS_KG_MOL
    )
    return np.clip(100.0 * hydrogen_moles / (hydrogen_moles + air_moles), 0.0, 100.0)


def _hazard_volume(
    result, field_vol_percent: np.ndarray, threshold_vol_percent: float = 4.0,
) -> float:
    dx = float(result.x_m[1] - result.x_m[0])
    return float(np.count_nonzero(field_vol_percent >= threshold_vol_percent) * dx**3)


def _masked_centerplane(result, field: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    index = int(np.argmin(np.abs(result.y_m - 1.20)))
    values = np.asarray(field[:, index, :], dtype=float).copy()
    values[result.solid_mask[:, index, :]] = np.nan
    return result.x_m, result.z_m, values


def _plot_field(axis, result, field, title, label, *, cmap="viridis", symmetric=False):
    x, z, values = _masked_centerplane(result, field)
    if symmetric:
        limit = max(float(np.nanmax(np.abs(values))), 1.0e-12)
        image = axis.pcolormesh(x, z, values, shading="nearest", cmap=cmap,
                                vmin=-limit, vmax=limit)
    else:
        image = axis.pcolormesh(x, z, values, shading="nearest", cmap=cmap)
    if label == "H2 (vol%)" and np.nanmax(values) >= 4.0:
        axis.contour(x, z, values, levels=[4.0], colors="white", linewidths=1.3)
    obstacle = result.solid_mask[:, int(np.argmin(np.abs(result.y_m - 1.20))), :]
    axis.contourf(x, z, obstacle.astype(float), levels=[0.5, 1.5], colors=["0.25"])
    axis.set_title(title)
    axis.set_xlabel("x (m)")
    axis.set_ylabel("z (m)")
    plt.colorbar(image, ax=axis, label=label)


def main() -> int:
    base = _config()
    receptors = (
        FastObstacleReceptor3D("windward", (1.20, 1.20, 0.50)),
        FastObstacleReceptor3D("roof", (1.70, 1.20, 1.20)),
        FastObstacleReceptor3D("side", (1.70, 0.60, 0.50)),
        FastObstacleReceptor3D("lee", (2.40, 1.20, 0.50)),
    )
    neutral = solve_resolved_obstacle_flow_3d(
        replace(base, thermal_buoyancy_coupling=False), receptors=receptors,
    )
    active = solve_resolved_obstacle_flow_3d(
        replace(base, thermal_buoyancy_coupling=True), receptors=receptors,
    )
    assert active.temperature_k is not None
    assert active.mixture_density_kg_m3 is not None
    assert active.buoyancy_acceleration_m_s2 is not None
    neutral_maximum_vol_percent = _neutral_vol_percent(
        neutral, neutral.maximum_concentration_kg_m3,
    )
    active_maximum_vol_percent = active.maximum_concentration_vol_percent

    record = {
        "schema": "degali.obstacle-thermal-buoyancy-development-comparison.v2",
        "status": "withheld",
        "purpose": "development witness; not external validation",
        "config": _json_safe(asdict(base)),
        "neutral": {
            "runtime_s": neutral.diagnostics["total_runtime_s"],
            "maximum_h2_vol_percent": float(np.nanmax(neutral_maximum_vol_percent)),
            "four_vol_percent_volume_m3": _hazard_volume(
                neutral, neutral_maximum_vol_percent,
            ),
            "maximum_mass_residual_kg": neutral.maximum_mass_residual_kg,
        },
        "active_thermal": {
            "runtime_s": active.diagnostics["total_runtime_s"],
            "maximum_h2_vol_percent": float(np.nanmax(
                active_maximum_vol_percent
            )),
            "four_vol_percent_volume_m3": _hazard_volume(
                active, active_maximum_vol_percent,
            ),
            "minimum_temperature_k": float(np.nanmin(active.temperature_k)),
            "minimum_mixture_density_kg_m3": float(np.nanmin(
                active.mixture_density_kg_m3
            )),
            "maximum_mixture_density_kg_m3": float(np.nanmax(
                active.mixture_density_kg_m3
            )),
            "maximum_mass_residual_kg": active.maximum_mass_residual_kg,
            "maximum_energy_residual_j": active.maximum_energy_residual_j,
            "diagnostics": _json_safe(active.diagnostics),
        },
        "comparison": {
            "maximum_absolute_h2_difference_vol_percent_point": float(np.nanmax(
                np.abs(active_maximum_vol_percent - neutral_maximum_vol_percent)
            )),
            "four_vol_percent_volume_difference_m3": (
                _hazard_volume(active, active_maximum_vol_percent)
                - _hazard_volume(neutral, neutral_maximum_vol_percent)
            ),
            "receptor_peak_vol_percent": {
                label: float(np.max(trace))
                for label, trace in active.receptor_traces_vol_percent.items()
            },
        },
        "limitations": [
            "variable-density mass-flux continuity is enforced, but carrier thermodynamics remains an active-scalar approximation",
            "constant heat capacities omit real-gas and phase-change thermodynamics",
            "flow Reynolds number and boundaries are not externally qualified",
        ],
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "obstacle_thermal_buoyancy_comparison.json"
    figure_path = OUTPUT / "obstacle_thermal_buoyancy_comparison.png"
    json_path.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")

    figure, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    _plot_field(
        axes[0, 0], neutral, neutral_maximum_vol_percent,
        "Neutral one-way: maximum H2", "H2 (vol%)", cmap="magma",
    )
    _plot_field(
        axes[0, 1], active, active_maximum_vol_percent,
        "Active thermal: maximum H2", "H2 (vol%)", cmap="magma",
    )
    _plot_field(
        axes[0, 2], active,
        active_maximum_vol_percent - neutral_maximum_vol_percent,
        "Active minus neutral maximum H2", "difference (vol%-point)",
        cmap="coolwarm", symmetric=True,
    )
    _plot_field(
        axes[1, 0], active, active.temperature_k,
        "Active thermal: final temperature", "temperature (K)", cmap="coolwarm",
    )
    _plot_field(
        axes[1, 1], active, active.buoyancy_acceleration_m_s2,
        "Active thermal: final vertical buoyancy", "a_z (m/s2)",
        cmap="coolwarm", symmetric=True,
    )
    for label, trace in active.receptor_traces_vol_percent.items():
        axes[1, 2].plot(active.time_s, trace, label=label)
    axes[1, 2].axhline(4.0, color="black", linestyle="--", linewidth=1.0,
                       label="4 vol%")
    axes[1, 2].set_title("Active-thermal receptor histories")
    axes[1, 2].set_xlabel("time (s)")
    axes[1, 2].set_ylabel("H2 (vol%)")
    axes[1, 2].legend(fontsize=8)
    figure.suptitle(
        "DEGALI finite-block thermal-buoyancy development witness — WITHHELD",
        fontsize=15,
    )
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "neutral_runtime_s": record["neutral"]["runtime_s"],
        "active_runtime_s": record["active_thermal"]["runtime_s"],
        "minimum_temperature_k": record["active_thermal"]["minimum_temperature_k"],
        "maximum_density_departure": active.diagnostics[
            "maximum_relative_mixture_density_departure"
        ],
        "active_scalar_applicability_passed": active.diagnostics[
            "active_scalar_applicability_passed"
        ],
        "thermodynamic_corrector_nonconverged_steps": active.diagnostics[
            "thermodynamic_corrector_nonconverged_steps"
        ],
        "maximum_variable_density_continuity_residual_kg_m3_s": (
            active.diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ]
        ),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
