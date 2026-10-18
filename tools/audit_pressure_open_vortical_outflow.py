"""Audit finite-amplitude vortex passage through the D3Q19 pressure outlet.

A pressure-balanced isothermal Gaussian vortex is convected by a uniform
carrier through the short-domain x outlet.  A longer domain supplies a
no-return reference over the audit horizon.  The comparison measures velocity,
vorticity and density contamination in a buffered common interior, plus the
vorticity left behind after the vortex has exited.  This is a nonlinear
vortical boundary-component test; it is not a hydrogen-plume validation.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

from degali.addons.resolved_obstacle_flow_3d import (
    _CS2_3D,
    _advance_d3q19,
    _equilibrium_d3q19,
    _moments_d3q19,
)


SCHEMA = "degali.pressure-open-vortical-outflow-audit.v1"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-vortical-outflow-audit-2026-10-10.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-vortical-outflow-audit-2026-10-10.png"
)
SHORT_NX = 96
REFERENCE_NX = 448
NY = 64
NZ = 4
VORTEX_CENTER_X = 32.0
VORTEX_CENTER_Y = (NY - 1.0) / 2.0
VORTEX_SIGMA_CELLS = 8.0
CARRIER_SPEED = 0.08
PEAK_SWIRL_SPEED = 0.015
LATTICE_VISCOSITY = 0.01
STEPS = 1150
SAMPLE_STRIDE = 5
BUFFER_CELLS = 8
INTERACTION_SIGMA_MULTIPLE = 3.0
NORMALISED_INTERIOR_ERROR_LIMIT = 0.05
RESIDUAL_ENSTROPHY_RATIO_LIMIT = 0.05
MEAN_DENSITY_ERROR_LIMIT = 1.0e-4
CANDIDATE_SCHEME = "characteristic_nscbc_convective"
PREDECESSOR_SCHEME = "characteristic_multidimensional"


def _initial_distribution(nx: int) -> np.ndarray:
    """Return an isothermal pressure-balanced Gaussian vortex."""
    x = np.arange(nx, dtype=float)[None, None, :]
    y = np.arange(NY, dtype=float)[None, :, None]
    dx = x - VORTEX_CENTER_X
    dy = y - VORTEX_CENTER_Y
    radius_squared = dx**2 + dy**2
    envelope = np.exp(
        -0.5 * radius_squared / VORTEX_SIGMA_CELLS**2
    )
    # u_theta = amplitude * r/sigma * exp(-r^2/(2 sigma^2)); choosing
    # amplitude=peak*sqrt(e) makes PEAK_SWIRL_SPEED the exact continuum peak.
    amplitude = PEAK_SWIRL_SPEED * math.sqrt(math.e)
    swirl_factor = amplitude / VORTEX_SIGMA_CELLS * envelope
    ux_2d = CARRIER_SPEED - swirl_factor * dy
    uy_2d = swirl_factor * dx
    # Isothermal radial balance: c_s^2 d(ln rho)/dr = u_theta^2/r.
    density_2d = np.exp(
        -(amplitude**2 / (2.0 * _CS2_3D))
        * np.exp(-radius_squared / VORTEX_SIGMA_CELLS**2)
    )
    shape = (NZ, NY, nx)
    density = np.broadcast_to(density_2d, shape).copy()
    ux = np.broadcast_to(ux_2d, shape).copy()
    uy = np.broadcast_to(uy_2d, shape).copy()
    uz = np.zeros(shape, dtype=float)
    return _equilibrium_d3q19(density, ux, uy, uz)


def _state(distribution: np.ndarray) -> tuple[np.ndarray, ...]:
    solid = np.zeros(distribution.shape[1:], dtype=bool)
    return _moments_d3q19(distribution, solid)


def _vorticity_z(ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    mean_ux = np.mean(ux, axis=0)
    mean_uy = np.mean(uy, axis=0)
    return np.gradient(mean_uy, axis=1) - np.gradient(mean_ux, axis=0)


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def _normalisations() -> dict[str, float]:
    density, ux, uy, uz = _state(_initial_distribution(REFERENCE_NX))
    common = (
        slice(None), slice(BUFFER_CELLS, -BUFFER_CELLS),
        slice(BUFFER_CELLS, SHORT_NX - BUFFER_CELLS),
    )
    velocity = np.sqrt(
        (ux[common] - CARRIER_SPEED) ** 2
        + uy[common] ** 2 + uz[common] ** 2
    )
    vorticity = _vorticity_z(ux, uy)[
        BUFFER_CELLS:-BUFFER_CELLS,
        BUFFER_CELLS:SHORT_NX - BUFFER_CELLS,
    ]
    return {
        "velocity_rms": _rms(velocity),
        "density_rms": _rms(density[common] - 1.0),
        "vorticity_rms": _rms(vorticity),
        "enstrophy": float(np.sum(vorticity**2)),
    }


def classify_vortical_outflow(
    *,
    maximum_velocity_error: float,
    maximum_vorticity_error: float,
    maximum_density_error: float,
    residual_enstrophy_ratio: float,
    final_mean_density_error: float,
) -> dict[str, Any]:
    values = {
        "maximum_normalised_velocity_l2": maximum_velocity_error,
        "maximum_normalised_vorticity_l2": maximum_vorticity_error,
        "maximum_normalised_density_l2": maximum_density_error,
        "residual_enstrophy_ratio": residual_enstrophy_ratio,
        "final_mean_density_error": final_mean_density_error,
    }
    if any(not math.isfinite(float(value)) for value in values.values()):
        raise ValueError("vortical-outflow metrics must be finite")
    component_passes = {
        "velocity": maximum_velocity_error <= NORMALISED_INTERIOR_ERROR_LIMIT,
        "vorticity": maximum_vorticity_error <= NORMALISED_INTERIOR_ERROR_LIMIT,
        "density": maximum_density_error <= NORMALISED_INTERIOR_ERROR_LIMIT,
        "residual_enstrophy": (
            residual_enstrophy_ratio <= RESIDUAL_ENSTROPHY_RATIO_LIMIT
        ),
        "mean_density": final_mean_density_error <= MEAN_DENSITY_ERROR_LIMIT,
    }
    passed = all(component_passes.values())
    return {
        **values,
        "normalised_interior_error_limit": NORMALISED_INTERIOR_ERROR_LIMIT,
        "residual_enstrophy_ratio_limit": RESIDUAL_ENSTROPHY_RATIO_LIMIT,
        "mean_density_error_limit": MEAN_DENSITY_ERROR_LIMIT,
        "component_passes": component_passes,
        "passed": passed,
    }


def _advance(distribution: np.ndarray, scheme: str) -> np.ndarray:
    solid = np.zeros(distribution.shape[1:], dtype=bool)
    advanced, _, _ = _advance_d3q19(
        distribution,
        solid=solid,
        inlet_u_lattice=np.full(NZ, CARRIER_SPEED),
        lattice_viscosity=LATTICE_VISCOSITY,
        flow_boundary_scheme=scheme,
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary="all_symmetry_free_slip",
        subgrid_model="none",
        smagorinsky_constant=0.0,
        wale_constant=0.0,
        maximum_relaxation_time=1.95,
    )
    return advanced


def run_audit(steps: int = STEPS) -> dict[str, Any]:
    if not isinstance(steps, int) or steps < 1:
        raise ValueError("steps must be a positive integer")
    started = time.perf_counter()
    reference = _initial_distribution(REFERENCE_NX)
    short = {
        PREDECESSOR_SCHEME: _initial_distribution(SHORT_NX),
        CANDIDATE_SCHEME: _initial_distribution(SHORT_NX),
    }
    norms = _normalisations()
    common_short = (
        slice(None), slice(BUFFER_CELLS, -BUFFER_CELLS),
        slice(BUFFER_CELLS, SHORT_NX - BUFFER_CELLS),
    )
    common_ref = (
        slice(None), slice(BUFFER_CELLS, -BUFFER_CELLS),
        slice(BUFFER_CELLS, SHORT_NX - BUFFER_CELLS),
    )
    vorticity_slice = (
        slice(BUFFER_CELLS, -BUFFER_CELLS),
        slice(BUFFER_CELLS, SHORT_NX - BUFFER_CELLS),
    )
    interaction_start = int(math.floor(
        (SHORT_NX - 1.0 - VORTEX_CENTER_X
         - INTERACTION_SIGMA_MULTIPLE * VORTEX_SIGMA_CELLS)
        / CARRIER_SPEED
    ))
    traces: dict[str, dict[str, list[float]]] = {
        scheme: {
            "step": [],
            "normalised_velocity_l2": [],
            "normalised_vorticity_l2": [],
            "normalised_density_l2": [],
            "residual_enstrophy_ratio": [],
            "mean_density_error": [],
            "probe_transverse_velocity_rms": [],
        }
        for scheme in short
    }
    for step in range(steps + 1):
        if step % SAMPLE_STRIDE == 0 or step == steps:
            ref_rho, ref_ux, ref_uy, ref_uz = _state(reference)
            ref_vorticity = _vorticity_z(ref_ux, ref_uy)[vorticity_slice]
            for scheme, distribution in short.items():
                rho, ux, uy, uz = _state(distribution)
                velocity_difference = np.sqrt(
                    (ux[common_short] - ref_ux[common_ref]) ** 2
                    + (uy[common_short] - ref_uy[common_ref]) ** 2
                    + (uz[common_short] - ref_uz[common_ref]) ** 2
                )
                vorticity = _vorticity_z(ux, uy)[vorticity_slice]
                trace = traces[scheme]
                trace["step"].append(step)
                trace["normalised_velocity_l2"].append(
                    _rms(velocity_difference) / norms["velocity_rms"]
                )
                trace["normalised_vorticity_l2"].append(
                    _rms(vorticity - ref_vorticity)
                    / norms["vorticity_rms"]
                )
                trace["normalised_density_l2"].append(
                    _rms(rho[common_short] - ref_rho[common_ref])
                    / norms["density_rms"]
                )
                trace["residual_enstrophy_ratio"].append(
                    float(np.sum(vorticity**2)) / norms["enstrophy"]
                )
                trace["mean_density_error"].append(
                    abs(float(np.mean(rho)) - 1.0)
                )
                trace["probe_transverse_velocity_rms"].append(
                    _rms(uy[:, BUFFER_CELLS:-BUFFER_CELLS, 60])
                )
        if step == steps:
            break
        reference = _advance(reference, CANDIDATE_SCHEME)
        for scheme in short:
            short[scheme] = _advance(short[scheme], scheme)

    schemes: dict[str, Any] = {}
    for scheme, trace in traces.items():
        sample_steps = np.asarray(trace["step"])
        interaction = sample_steps >= interaction_start
        if not np.any(interaction):
            # Unit/smoke runs may intentionally stop before boundary contact.
            interaction = np.ones_like(sample_steps, dtype=bool)
        classification = classify_vortical_outflow(
            maximum_velocity_error=float(np.max(
                np.asarray(trace["normalised_velocity_l2"])[interaction]
            )),
            maximum_vorticity_error=float(np.max(
                np.asarray(trace["normalised_vorticity_l2"])[interaction]
            )),
            maximum_density_error=float(np.max(
                np.asarray(trace["normalised_density_l2"])[interaction]
            )),
            residual_enstrophy_ratio=float(
                trace["residual_enstrophy_ratio"][-1]
            ),
            final_mean_density_error=float(
                trace["mean_density_error"][-1]
            ),
        )
        schemes[scheme] = {**classification, "trace": trace}
    passed = schemes[CANDIDATE_SCHEME]["passed"]
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "nonlinear_vortical_outflow_component_gate_passed"
            if passed else "nonlinear_vortical_outflow_component_gate_failed"
        ),
        "passed": passed,
        "promotion_allowed": False,
        "configuration": {
            "short_domain_cells_xyz": [SHORT_NX, NY, NZ],
            "reference_domain_cells_xyz": [REFERENCE_NX, NY, NZ],
            "vortex_center_cells": [VORTEX_CENTER_X, VORTEX_CENTER_Y],
            "vortex_sigma_cells": VORTEX_SIGMA_CELLS,
            "carrier_lattice_speed": CARRIER_SPEED,
            "peak_swirl_lattice_speed": PEAK_SWIRL_SPEED,
            "carrier_mach_number": CARRIER_SPEED / math.sqrt(_CS2_3D),
            "lattice_viscosity": LATTICE_VISCOSITY,
            "vortex_reynolds_number": (
                PEAK_SWIRL_SPEED * 2.0 * VORTEX_SIGMA_CELLS
                / LATTICE_VISCOSITY
            ),
            "time_steps": steps,
            "sample_stride": SAMPLE_STRIDE,
            "buffer_cells": BUFFER_CELLS,
            "interaction_start_step": interaction_start,
            "candidate_boundary_scheme": CANDIDATE_SCHEME,
            "predecessor_boundary_scheme": PREDECESSOR_SCHEME,
            "transverse_boundary": "symmetry/free-slip",
        },
        "initial_normalisations": norms,
        "schemes": schemes,
        "scope": (
            "finite_amplitude_isothermal_gaussian_vortex_outflow_component_"
            "test_not_hydrogen_plume_or_external_validation"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    metrics = (
        ("normalised_velocity_l2", "velocity L2 / initial"),
        ("normalised_vorticity_l2", "vorticity L2 / initial"),
        ("normalised_density_l2", "density L2 / initial"),
        ("residual_enstrophy_ratio", "common-domain enstrophy / initial"),
    )
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), sharex=True)
    for axis, (metric, label) in zip(axes.flat, metrics):
        for scheme, result in record["schemes"].items():
            trace = result["trace"]
            axis.plot(trace["step"], trace[metric], label=scheme)
        axis.axvline(
            record["configuration"]["interaction_start_step"],
            color="0.4", linestyle="--", linewidth=1.0,
        )
        if metric != "residual_enstrophy_ratio":
            axis.axhline(
                NORMALISED_INTERIOR_ERROR_LIMIT,
                color="tab:red", linestyle=":", linewidth=1.0,
            )
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
    axes[-1, 0].set_xlabel("lattice step")
    axes[-1, 1].set_xlabel("lattice step")
    axes[0, 0].legend(fontsize=8)
    figure.suptitle(
        "Finite-amplitude vortex passage through pressure-open outlet"
    )
    figure.tight_layout()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_audit(record: dict[str, Any], output: str | Path) -> None:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    parser.add_argument("--steps", type=int, default=STEPS)
    args = parser.parse_args(argv)
    record = run_audit(args.steps)
    write_audit(record, args.output)
    plot_audit(record, args.plot)
    candidate = record["schemes"][CANDIDATE_SCHEME]
    print(
        f"wrote {args.output} status={record['status']} "
        f"velocity={candidate['maximum_normalised_velocity_l2']:.6f} "
        f"vorticity={candidate['maximum_normalised_vorticity_l2']:.6f} "
        f"density={candidate['maximum_normalised_density_l2']:.6f} "
        f"runtime_s={record['runtime_s']:.3f}"
    )
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
