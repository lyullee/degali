"""Thirty-degree oblique acoustic reflection audit for the D3Q19 side outlet.

The tangential x direction is periodic and z is symmetry isolated.  A small
Fourier-Gaussian acoustic packet travels toward y-max at 30 degrees from the
boundary normal.  A longer y-domain supplies the no-return reference.  The
complex tangential Fourier coefficient at the probe separates the incident
and reflected packet without fitting a boundary parameter.
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
)


SCHEMA = "degali.pressure-open-oblique-reflection-audit.v1"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-oblique-reflection-audit-2026-10-10.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-oblique-reflection-audit-2026-10-10.png"
)
SHORT_NY = 200
REFERENCE_NY = 400
NX = 32
NZ = 4
SOURCE_Y = 50.0
PROBE_Y = 130
PULSE_SIGMA_CELLS = 16.0
PULSE_AMPLITUDE = 1.0e-4
INCIDENCE_ANGLE_DEGREES = 30.0
TANGENTIAL_MODE = 1
LATTICE_VISCOSITY = 0.02
STEPS = 560
REFLECTION_AMPLITUDE_RATIO_LIMIT = 0.05


def _wavenumbers() -> tuple[float, float]:
    angle = math.radians(INCIDENCE_ANGLE_DEGREES)
    tangential = 2.0 * math.pi * TANGENTIAL_MODE / NX
    normal = tangential / math.tan(angle)
    return tangential, normal


def _initial_distribution(ny: int) -> np.ndarray:
    tangential_wavenumber, normal_wavenumber = _wavenumbers()
    x = np.arange(NX, dtype=float)
    y = np.arange(ny, dtype=float)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    envelope = np.exp(
        -0.5 * ((yy - SOURCE_Y) / PULSE_SIGMA_CELLS) ** 2
    )
    phase = (
        tangential_wavenumber * xx
        + normal_wavenumber * (yy - SOURCE_Y)
    )
    density_perturbation = PULSE_AMPLITUDE * envelope * np.cos(phase)
    density_plane = 1.0 + density_perturbation
    density = np.broadcast_to(
        density_plane[None, :, :], (NZ, ny, NX),
    ).copy()
    sound_speed = math.sqrt(_CS2_3D)
    angle = math.radians(INCIDENCE_ANGLE_DEGREES)
    velocity_x_plane = (
        sound_speed * math.sin(angle)
        * density_perturbation / density_plane
    )
    velocity_y_plane = (
        sound_speed * math.cos(angle)
        * density_perturbation / density_plane
    )
    velocity_x = np.broadcast_to(
        velocity_x_plane[None, :, :], density.shape,
    ).copy()
    velocity_y = np.broadcast_to(
        velocity_y_plane[None, :, :], density.shape,
    ).copy()
    return _equilibrium_d3q19(
        density, velocity_x, velocity_y, np.zeros_like(density),
    )


def _probe_mode(distribution: np.ndarray) -> complex:
    density = np.sum(distribution, axis=0)
    perturbation = np.mean(density[:, PROBE_Y, :] - 1.0, axis=0)
    tangential_wavenumber, _ = _wavenumbers()
    phase = np.exp(
        -1j * tangential_wavenumber * np.arange(NX, dtype=float)
    )
    return complex(2.0 * np.mean(perturbation * phase))


def _trace(
    ny: int,
    boundary_scheme: str = "characteristic_nonreflecting",
) -> np.ndarray:
    distribution = _initial_distribution(ny)
    solid = np.zeros((NZ, ny, NX), dtype=bool)
    trace = np.empty(STEPS + 1, dtype=complex)
    for step in range(STEPS + 1):
        trace[step] = _probe_mode(distribution)
        if step == STEPS:
            break
        distribution, _, _ = _advance_d3q19(
            distribution,
            solid=solid,
            inlet_u_lattice=np.zeros(NZ),
            lattice_viscosity=LATTICE_VISCOSITY,
            flow_boundary_scheme=boundary_scheme,
            outlet_lattice_density=1.0,
            nonstreamwise_wall_boundary="all_pressure_open",
            subgrid_model="none",
            smagorinsky_constant=0.0,
            wale_constant=0.0,
            maximum_relaxation_time=1.95,
            pressure_open_nonstream_faces=("y_min", "y_max"),
            apply_streamwise_boundaries=False,
        )
    return trace


def classify_oblique_reflection(
    short_trace: np.ndarray,
    reference_trace: np.ndarray,
) -> dict[str, Any]:
    short = np.asarray(short_trace, dtype=complex)
    reference = np.asarray(reference_trace, dtype=complex)
    if short.shape != reference.shape or short.ndim != 1:
        raise ValueError("short and reference traces must be equal-length vectors")
    if np.any(~np.isfinite(short)) or np.any(~np.isfinite(reference)):
        raise ValueError("reflection traces must be finite")
    sound_speed = math.sqrt(_CS2_3D)
    normal_speed = sound_speed * math.cos(
        math.radians(INCIDENCE_ANGLE_DEGREES)
    )
    direct_center = (PROBE_Y - SOURCE_Y) / normal_speed
    reflected_center = (
        (SHORT_NY - 1.0 - SOURCE_Y)
        + (SHORT_NY - 1.0 - PROBE_Y)
    ) / normal_speed
    half_window = 3.0 * PULSE_SIGMA_CELLS / normal_speed
    direct_slice = slice(
        max(0, int(math.floor(direct_center - half_window))),
        min(short.size, int(math.ceil(direct_center + half_window)) + 1),
    )
    reflected_slice = slice(
        max(0, int(math.floor(reflected_center - half_window))),
        min(short.size, int(math.ceil(reflected_center + half_window)) + 1),
    )
    incident_amplitude = float(np.max(np.abs(reference[direct_slice])))
    reflected_amplitude = float(np.max(np.abs(
        short[reflected_slice] - reference[reflected_slice]
    )))
    ratio = (
        math.inf if incident_amplitude <= 0.0
        else reflected_amplitude / incident_amplitude
    )
    passed = bool(
        math.isfinite(ratio)
        and ratio <= REFLECTION_AMPLITUDE_RATIO_LIMIT
    )
    return {
        "incident_window_steps": [direct_slice.start, direct_slice.stop - 1],
        "reflection_window_steps": [
            reflected_slice.start, reflected_slice.stop - 1,
        ],
        "incident_mode_amplitude": incident_amplitude,
        "reflected_mode_amplitude": reflected_amplitude,
        "reflection_amplitude_ratio": ratio,
        "reflection_amplitude_ratio_limit": (
            REFLECTION_AMPLITUDE_RATIO_LIMIT
        ),
        "passed": passed,
        "status": (
            "canonical_30deg_oblique_reflection_gate_passed"
            if passed else "canonical_30deg_oblique_reflection_gate_failed"
        ),
        "promotion_allowed": False,
    }


def run_audit() -> dict[str, Any]:
    started = time.perf_counter()
    short = _trace(SHORT_NY)
    reference = _trace(REFERENCE_NY)
    classification = classify_oblique_reflection(short, reference)
    tangential_wavenumber, normal_wavenumber = _wavenumbers()
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **classification,
        "configuration": {
            "short_domain_cells_xyz": [NX, SHORT_NY, NZ],
            "reference_domain_cells_xyz": [NX, REFERENCE_NY, NZ],
            "source_y_cell": SOURCE_Y,
            "probe_y_cell": PROBE_Y,
            "pulse_sigma_cells": PULSE_SIGMA_CELLS,
            "pulse_density_amplitude": PULSE_AMPLITUDE,
            "incidence_angle_degrees_from_normal": (
                INCIDENCE_ANGLE_DEGREES
            ),
            "tangential_mode": TANGENTIAL_MODE,
            "tangential_wavenumber_rad_per_cell": tangential_wavenumber,
            "normal_wavenumber_rad_per_cell": normal_wavenumber,
            "lattice_viscosity": LATTICE_VISCOSITY,
            "time_steps": STEPS,
            "tangential_x_boundary": "periodic",
            "vertical_z_boundary": "specular symmetry",
            "normal_y_boundary": "isothermal characteristic open",
        },
        "short_probe_mode_real": short.real.tolist(),
        "short_probe_mode_imaginary": short.imag.tolist(),
        "reference_probe_mode_real": reference.real.tolist(),
        "reference_probe_mode_imaginary": reference.imag.tolist(),
        "scope": (
            "single_angle_single_wavenumber_linear_component_test_not_a_plume_domain_test"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    short = (
        np.asarray(record["short_probe_mode_real"])
        + 1j * np.asarray(record["short_probe_mode_imaginary"])
    )
    reference = (
        np.asarray(record["reference_probe_mode_real"])
        + 1j * np.asarray(record["reference_probe_mode_imaginary"])
    )
    difference = short - reference
    steps = np.arange(short.size)
    window = record["reflection_window_steps"]
    figure, axes = plt.subplots(
        2, 1, figsize=(9.2, 6.7), sharex=True, constrained_layout=True,
    )
    axes[0].plot(
        steps, np.abs(reference), color="#2166ac",
        label="long-domain mode magnitude",
    )
    axes[0].plot(
        steps, np.abs(short), color="#1b7837",
        label="short characteristic mode magnitude",
    )
    axes[0].axvspan(window[0], window[1], alpha=0.15, color="#fdae61")
    axes[0].set_ylabel("Tangential Fourier amplitude")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(steps, difference.real, label="real", color="#762a83")
    axes[1].plot(
        steps, difference.imag, label="imaginary", color="#e08214",
    )
    axes[1].plot(
        steps, np.abs(difference), label="magnitude", color="black",
        linewidth=1.2,
    )
    axes[1].axvspan(window[0], window[1], alpha=0.15, color="#fdae61")
    axes[1].set_xlabel("Lattice time step")
    axes[1].set_ylabel("Short - reference")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    figure.suptitle(
        "30° oblique D3Q19 characteristic-boundary reflection: "
        f"{float(record['reflection_amplitude_ratio']):.3%} "
        f"({record['status']})",
        fontsize=11,
    )
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    args = parser.parse_args(argv)
    record = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot_audit(record, args.plot)
    print(json.dumps({
        key: record[key]
        for key in (
            "status", "reflection_amplitude_ratio",
            "reflection_amplitude_ratio_limit", "runtime_s",
        )
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
