"""Canonical verification of the 3-D variable-density mass-flux projection."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.variable_density_projection_3d import (
    project_variable_density_mass_flux_3d,
)


OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"
SHAPE = (8, 10, 14)
SPACING = (0.20, 0.20, 0.20)
DT = 0.04


def _zero_flux(shape=SHAPE):
    nz, ny, nx = shape
    return (
        np.zeros((nz, ny, nx + 1)),
        np.zeros((nz, ny + 1, nx)),
        np.zeros((nz + 1, ny, nx)),
    )


def _run_cases():
    density = np.full(SHAPE, 1.20)
    fluid = np.zeros(SHAPE, dtype=bool)
    cases = []

    rng = np.random.default_rng(20261010)
    fx, fy, fz = _zero_flux()
    fx[:, :, 1:-1] = rng.normal(0.0, 0.08, fx[:, :, 1:-1].shape)
    fy[:, 1:-1, :] = rng.normal(0.0, 0.08, fy[:, 1:-1, :].shape)
    fz[1:-1, :, :] = rng.normal(0.0, 0.08, fz[1:-1, :, :].shape)
    cases.append(("uniform_random_flux", fx, fy, fz, density, density, 0.0, fluid))

    fx, fy, fz = _zero_flux()
    cases.append((
        "uniform_density_decrease", fx, fy, fz,
        density, np.full(SHAPE, 1.19), 0.0, fluid,
    ))

    fx, fy, fz = _zero_flux()
    source = np.zeros(SHAPE)
    source[1:3, 2:5, 2:5] = 0.08
    obstacle = np.zeros(SHAPE, dtype=bool)
    obstacle[:6, 4:7, 6:9] = True
    cases.append((
        "localized_source_with_obstacle", fx, fy, fz,
        density, density, source, obstacle,
    ))

    fx, fy, fz = _zero_flux()
    front = density.copy()
    front[:, :, :6] -= 0.018
    front[:, :, 6:10] += 0.006
    cases.append((
        "nonuniform_density_front", fx, fy, fz,
        density, front, 0.0, fluid,
    ))

    output = []
    started = time.perf_counter()
    for name, fx, fy, fz, old, new, source, solid in cases:
        case_started = time.perf_counter()
        result = project_variable_density_mass_flux_3d(
            fx, fy, fz,
            density_old_kg_m3=old,
            density_new_kg_m3=new,
            mass_source_rate_kg_m3_s=source,
            solid_mask=solid,
            spacing_m=SPACING,
            time_step_s=DT,
        )
        output.append({
            "name": name,
            "result": result,
            "density_old": old,
            "density_new": new,
            "source": np.broadcast_to(np.asarray(source), SHAPE),
            "solid": solid,
            "runtime_s": time.perf_counter() - case_started,
        })
    return output, time.perf_counter() - started


def main() -> int:
    cases, total_runtime = _run_cases()
    local_tolerance = 1.0e-8
    global_tolerance = 1.0e-9
    gate_reasons = []
    records = []
    for case in cases:
        result = case["result"]
        if result.maximum_absolute_residual_after_kg_m3_s > local_tolerance:
            gate_reasons.append(
                f"{case['name']} local residual exceeds {local_tolerance:.1e}"
            )
        if abs(result.global_mass_balance_residual_kg_s) > global_tolerance:
            gate_reasons.append(
                f"{case['name']} global residual exceeds {global_tolerance:.1e}"
            )
        records.append({
            "name": case["name"],
            "maximum_absolute_residual_before_kg_m3_s": (
                result.maximum_absolute_residual_before_kg_m3_s
            ),
            "maximum_absolute_residual_after_kg_m3_s": (
                result.maximum_absolute_residual_after_kg_m3_s
            ),
            "domain_mass_change_rate_kg_s": result.domain_mass_change_rate_kg_s,
            "integrated_source_rate_kg_s": result.integrated_source_rate_kg_s,
            "net_boundary_outflow_rate_kg_s": result.net_boundary_outflow_rate_kg_s,
            "global_mass_balance_residual_kg_s": (
                result.global_mass_balance_residual_kg_s
            ),
            "solver_method": result.solver_method,
            "solver_iterations": result.solver_iterations,
            "runtime_s": case["runtime_s"],
        })
    record = {
        "schema": "degali.variable-density-mass-flux-canonical-verification.v1",
        "status": "accepted" if not gate_reasons else "withheld",
        "scope": "numerical continuity operator only; not physical dispersion validation",
        "grid_shape_zyx": SHAPE,
        "spacing_m_xyz": SPACING,
        "time_step_s": DT,
        "local_residual_tolerance_kg_m3_s": local_tolerance,
        "global_residual_tolerance_kg_s": global_tolerance,
        "cases": records,
        "total_runtime_s": total_runtime,
        "gate_reasons": gate_reasons,
        "next_gate": (
            "coupled predictor-corrector for EOS density, species, enthalpy and momentum"
        ),
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "variable_density_mass_flux_projection_verification.json"
    figure_path = OUTPUT / "variable_density_mass_flux_projection_verification.png"
    json_path.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")

    names = [case["name"].replace("_", "\n") for case in cases]
    before = [
        case["result"].maximum_absolute_residual_before_kg_m3_s for case in cases
    ]
    after = [
        max(case["result"].maximum_absolute_residual_after_kg_m3_s, 1.0e-18)
        for case in cases
    ]
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    positions = np.arange(len(cases))
    axes[0, 0].bar(positions - 0.18, before, 0.36, label="before")
    axes[0, 0].bar(positions + 0.18, after, 0.36, label="after")
    axes[0, 0].set_yscale("log")
    axes[0, 0].axhline(local_tolerance, color="black", linestyle="--",
                       label="acceptance limit")
    axes[0, 0].set_xticks(positions, names)
    axes[0, 0].set_ylabel("max continuity residual (kg/m3/s)")
    axes[0, 0].set_title("Cellwise continuity closure")
    axes[0, 0].legend(fontsize=8)

    width = 0.24
    axes[0, 1].bar(
        positions - width,
        [case["result"].domain_mass_change_rate_kg_s for case in cases],
        width, label="dM/dt",
    )
    axes[0, 1].bar(
        positions,
        [case["result"].net_boundary_outflow_rate_kg_s for case in cases],
        width, label="boundary outflow",
    )
    axes[0, 1].bar(
        positions + width,
        [-case["result"].integrated_source_rate_kg_s for case in cases],
        width, label="- source",
    )
    axes[0, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[0, 1].set_xticks(positions, names)
    axes[0, 1].set_ylabel("mass-rate term (kg/s)")
    axes[0, 1].set_title("Global budget terms: sum = 0")
    axes[0, 1].legend(fontsize=8)

    front_case = cases[-1]
    mid_z, mid_y = SHAPE[0] // 2, SHAPE[1] // 2
    x = (np.arange(SHAPE[2]) + 0.5) * SPACING[0]
    axes[1, 0].plot(x, front_case["density_old"][mid_z, mid_y], label="old")
    axes[1, 0].plot(x, front_case["density_new"][mid_z, mid_y], label="new")
    axes[1, 0].set_xlabel("x (m)")
    axes[1, 0].set_ylabel("density (kg/m3)")
    axes[1, 0].set_title("Manufactured nonuniform density front")
    axes[1, 0].legend()

    obstacle_case = cases[2]
    result = obstacle_case["result"]
    j = SHAPE[1] // 2
    x_cell = (np.arange(SHAPE[2]) + 0.5) * SPACING[0]
    z_cell = (np.arange(SHAPE[0]) + 0.5) * SPACING[2]
    cell_fx = 0.5 * (
        result.mass_flux_x_kg_m2_s[:, j, :-1]
        + result.mass_flux_x_kg_m2_s[:, j, 1:]
    )
    cell_fz = 0.5 * (
        result.mass_flux_z_kg_m2_s[:-1, j, :]
        + result.mass_flux_z_kg_m2_s[1:, j, :]
    )
    speed = np.hypot(cell_fx, cell_fz)
    image = axes[1, 1].pcolormesh(
        x_cell, z_cell, speed, shading="nearest", cmap="viridis",
    )
    axes[1, 1].quiver(
        x_cell, z_cell, cell_fx, cell_fz, color="white", scale=0.16,
        width=0.0035,
    )
    solid_slice = obstacle_case["solid"][:, j, :]
    axes[1, 1].contourf(
        x_cell, z_cell, solid_slice.astype(float),
        levels=[0.5, 1.5], colors=["0.25"],
    )
    axes[1, 1].set_xlabel("x (m)")
    axes[1, 1].set_ylabel("z (m)")
    axes[1, 1].set_title("Localized source mass flux around impermeable block")
    plt.colorbar(image, ax=axes[1, 1], label="mass-flux magnitude (kg/m2/s)")
    figure.suptitle(
        f"Variable-density mass-flux projection — {record['status'].upper()}",
        fontsize=15,
    )
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "maximum_local_residual_after_kg_m3_s": max(after),
        "maximum_global_residual_kg_s": max(
            abs(case["result"].global_mass_balance_residual_kg_s)
            for case in cases
        ),
        "total_runtime_s": total_runtime,
        "gate_reasons": gate_reasons,
    }, indent=2))
    return 0 if not gate_reasons else 2


if __name__ == "__main__":
    raise SystemExit(main())
