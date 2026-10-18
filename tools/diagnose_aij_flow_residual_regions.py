"""Localise D3Q19 residuals for the AIJ open-boundary candidates."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import (
    _advance_d3q19,
    _equilibrium_d3q19,
    _inlet_profile,
    _moments_d3q19,
    _solid_mask,
)
from degali.validation.aij_case_h import read_aij_case_h
from tools.audit_aij_resolved_3d_grid import DATA, OUTPUT, _config


CELL_SIZE_M = 0.020
BOUNDARY_CASES = (
    ("equilibrium_copy", "Equilibrium/copy"),
    ("nonequilibrium_extrapolation", "NEE velocity/pressure"),
)


def _run(config) -> dict[str, object]:
    x, _y, z, solid = _solid_mask(config)
    dx = config.domain_m[0] / config.cells[0]
    physical_dt = config.lattice_inlet_speed * dx / config.wind_speed_m_s
    _inlet_physical, inlet_lattice = _inlet_profile(config, z)
    lattice_viscosity = config.kinematic_viscosity_m2_s * physical_dt / dx**2
    nz, ny, nx = solid.shape
    rho = np.ones_like(solid, dtype=float)
    ux = np.broadcast_to(inlet_lattice[:, None, None], (nz, ny, nx)).copy()
    uy = np.zeros_like(ux)
    uz = np.zeros_like(ux)
    ux[solid] = 0.0
    distribution = _equilibrium_d3q19(rho, ux, uy, uz)
    previous = (ux.copy(), uy.copy(), uz.copy())
    times = []
    x_profiles = []
    peak_x = []
    peak_value = []
    started = time.perf_counter()
    for step in range(1, config.flow_steps + 1):
        distribution, _tau, _eddy = _advance_d3q19(
            distribution,
            solid=solid,
            inlet_u_lattice=inlet_lattice,
            lattice_viscosity=lattice_viscosity,
            flow_boundary_scheme=config.flow_boundary_scheme,
            outlet_lattice_density=config.outlet_lattice_density,
            subgrid_model=config.subgrid_model,
            smagorinsky_constant=config.smagorinsky_constant,
            wale_constant=config.wale_constant,
            maximum_relaxation_time=config.maximum_relaxation_time,
        )
        if step % config.convergence_check_interval != 0 and step != config.flow_steps:
            continue
        _rho, current_u, current_v, current_w = _moments_d3q19(distribution, solid)
        change = np.sqrt(
            (current_u - previous[0])**2
            + (current_v - previous[1])**2
            + (current_w - previous[2])**2
        )
        reference = max(
            float(np.max(np.sqrt(
                current_u**2 + current_v**2 + current_w**2
            ))),
            float(np.max(inlet_lattice)),
            1.0e-15,
        )
        normalized = change / reference
        normalized[solid] = np.nan
        profile = np.nanmax(normalized, axis=(0, 1))
        index = int(np.nanargmax(profile))
        times.append(step * physical_dt)
        x_profiles.append(profile)
        peak_x.append(float(x[index]))
        peak_value.append(float(profile[index]))
        previous = (current_u.copy(), current_v.copy(), current_w.copy())
    profiles = np.asarray(x_profiles)
    return {
        "scheme": config.flow_boundary_scheme,
        "x_m": x.tolist(),
        "time_s": times,
        "x_residual_profiles": profiles.tolist(),
        "peak_x_m": peak_x,
        "peak_residual": peak_value,
        "final_inlet_zone_max": float(np.nanmax(profiles[-1, :3])),
        "final_obstacle_wake_zone_max": float(np.nanmax(
            profiles[-1, (x >= 0.13) & (x <= 0.65)]
        )),
        "final_outlet_zone_max": float(np.nanmax(profiles[-1, -3:])),
        "runtime_s": time.perf_counter() - started,
    }


def main() -> int:
    dataset = read_aij_case_h(DATA)
    base = replace(
        _config(dataset, CELL_SIZE_M),
        subgrid_model="none",
        smagorinsky_constant=0.0,
    )
    cases = [
        _run(replace(base, flow_boundary_scheme=scheme))
        for scheme, _label in BOUNDARY_CASES
    ]
    record = {
        "schema": "degali.aij-case-h-flow-residual-regions.v1",
        "status": "diagnostic_only",
        "cell_size_m": CELL_SIZE_M,
        "cells": list(base.cells),
        "flow_steps": base.flow_steps,
        "convergence_check_interval": base.convergence_check_interval,
        "normalization": "maximum velocity change over each check interval / current maximum speed",
        "cases": cases,
        "quantitative_hydrogen_validation": False,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "aij_case_h_flow_residual_regions.json"
    figure_path = OUTPUT / "aij_case_h_flow_residual_regions.png"
    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), constrained_layout=True)
    for column, (case, (_scheme, label)) in enumerate(zip(cases, BOUNDARY_CASES)):
        x = np.asarray(case["x_m"])
        times = np.asarray(case["time_s"])
        profiles = np.asarray(case["x_residual_profiles"])
        mesh = axes[0, column].pcolormesh(
            x, times, np.maximum(profiles, 1.0e-8),
            shading="nearest", cmap="magma", norm="log",
        )
        axes[0, column].axvspan(0.15, 0.25, color="#38bdf8", alpha=0.25)
        axes[0, column].set(
            title=label,
            xlabel="x (m)",
            ylabel="Physical horizon (s)",
        )
        figure.colorbar(mesh, ax=axes[0, column], label="Normalised residual")
        axes[1, column].semilogy(
            times, case["peak_residual"], "o-", color="#dc2626",
            label="domain maximum",
        )
        axes[1, column].plot(
            times, np.asarray(case["peak_x_m"]) + 1.0e-8,
            "s--", color="#2563eb", label="peak x (m)",
        )
        axes[1, column].set(
            title="Peak magnitude and location",
            xlabel="Physical horizon (s)",
        )
        axes[1, column].grid(alpha=0.25)
        axes[1, column].legend(fontsize=8)
    figure.suptitle("AIJ Case H — spatial localisation of flow residual", fontsize=14)
    figure.savefig(figure_path, dpi=220)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "final_regions": {
            case["scheme"]: {
                "inlet": case["final_inlet_zone_max"],
                "obstacle_wake": case["final_obstacle_wake_zone_max"],
                "outlet": case["final_outlet_zone_max"],
                "peak_x_m": case["peak_x_m"][-1],
            }
            for case in cases
        },
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
