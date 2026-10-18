"""Angle/wavelength sweep for the D3Q19 characteristic side outlet.

Each case reuses the same Fourier-Gaussian acoustic packet and long-domain
reference as the single 30-degree audit.  The declared quantitative scope is
limited to waves resolved by at least eight lattice cells per wavelength.
Under-resolved cases are still reported, but cannot make the gate pass or fail.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np


TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import audit_pressure_open_boundary_oblique_reflection as oblique


SCHEMA = "degali.pressure-open-boundary-sweep.v1"
DEFAULT_OUTPUT = (
    TOOLS_DIR.parents[0]
    / "outputs/pressure-open-angle-wavelength-sweep-2026-10-10.json"
)
DEFAULT_PLOT = (
    TOOLS_DIR.parents[0]
    / "outputs/pressure-open-angle-wavelength-sweep-2026-10-10.png"
)
ANGLES_DEGREES = (15.0, 30.0, 45.0, 60.0)
TANGENTIAL_MODES = (1, 2)
MINIMUM_RESOLVED_WAVELENGTH_CELLS = 8.0
REFLECTION_AMPLITUDE_RATIO_LIMIT = 0.05


def _case_definition(
    angle_degrees: float,
    tangential_mode: int,
    boundary_scheme: str = "characteristic_nonreflecting",
) -> dict[str, Any]:
    if not 0.0 < angle_degrees < 90.0:
        raise ValueError("angle_degrees must be between 0 and 90")
    if tangential_mode < 1:
        raise ValueError("tangential_mode must be positive")
    wavelength = (
        oblique.NX * math.sin(math.radians(angle_degrees))
        / tangential_mode
    )
    sound_speed = math.sqrt(oblique._CS2_3D)
    normal_speed = sound_speed * math.cos(math.radians(angle_degrees))
    reflected_center = (
        (oblique.SHORT_NY - 1.0 - oblique.SOURCE_Y)
        + (oblique.SHORT_NY - 1.0 - oblique.PROBE_Y)
    ) / normal_speed
    half_window = 3.0 * oblique.PULSE_SIGMA_CELLS / normal_speed
    steps = int(math.ceil(reflected_center + half_window)) + 4
    return {
        "angle_degrees_from_normal": float(angle_degrees),
        "tangential_mode": int(tangential_mode),
        "boundary_scheme": boundary_scheme,
        "wavelength_cells": float(wavelength),
        "resolved_for_gate": bool(
            wavelength + 1.0e-12 >= MINIMUM_RESOLVED_WAVELENGTH_CELLS
        ),
        "time_steps": steps,
    }


def build_case_matrix(
    boundary_scheme: str = "characteristic_nonreflecting",
) -> list[dict[str, Any]]:
    return [
        _case_definition(angle, mode, boundary_scheme)
        for angle in ANGLES_DEGREES
        for mode in TANGENTIAL_MODES
    ]


def _run_case(case: dict[str, Any]) -> dict[str, Any]:
    started = time.perf_counter()
    # Workers are separate processes, so changing these module-level audit
    # parameters cannot leak between cases.
    oblique.INCIDENCE_ANGLE_DEGREES = float(case["angle_degrees_from_normal"])
    oblique.TANGENTIAL_MODE = int(case["tangential_mode"])
    oblique.STEPS = int(case["time_steps"])
    boundary_scheme = str(case["boundary_scheme"])
    short = oblique._trace(oblique.SHORT_NY, boundary_scheme)
    reference = oblique._trace(oblique.REFERENCE_NY, boundary_scheme)
    classification = oblique.classify_oblique_reflection(short, reference)
    ratio = float(classification["reflection_amplitude_ratio"])
    resolved = bool(case["resolved_for_gate"])
    return {
        **case,
        "incident_mode_amplitude": float(
            classification["incident_mode_amplitude"]
        ),
        "reflected_mode_amplitude": float(
            classification["reflected_mode_amplitude"]
        ),
        "reflection_amplitude_ratio": ratio,
        "reflection_amplitude_ratio_limit": REFLECTION_AMPLITUDE_RATIO_LIMIT,
        "case_gate_passed": bool(
            resolved and math.isfinite(ratio)
            and ratio <= REFLECTION_AMPLITUDE_RATIO_LIMIT
        ),
        "gate_disposition": (
            "pass" if resolved and ratio <= REFLECTION_AMPLITUDE_RATIO_LIMIT
            else "fail" if resolved
            else "reported_but_excluded_underresolved"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def classify_sweep(cases: list[dict[str, Any]]) -> dict[str, Any]:
    in_scope = [case for case in cases if case["resolved_for_gate"]]
    excluded = [case for case in cases if not case["resolved_for_gate"]]
    passed = bool(in_scope) and all(case["case_gate_passed"] for case in in_scope)
    worst = max(
        in_scope,
        key=lambda case: float(case["reflection_amplitude_ratio"]),
        default=None,
    )
    return {
        "passed": passed,
        "status": (
            "resolved_angle_wavelength_sweep_gate_passed"
            if passed else "resolved_angle_wavelength_sweep_gate_failed"
        ),
        "promotion_allowed": False,
        "in_scope_case_count": len(in_scope),
        "excluded_underresolved_case_count": len(excluded),
        "worst_in_scope_case": worst,
    }


def run_audit(
    workers: int | None = None,
    boundary_scheme: str = "characteristic_nonreflecting",
) -> dict[str, Any]:
    started = time.perf_counter()
    definitions = build_case_matrix(boundary_scheme)
    worker_count = min(len(definitions), workers or 4)
    if worker_count == 1:
        cases = [_run_case(case) for case in definitions]
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            cases = list(executor.map(_run_case, definitions))
    classification = classify_sweep(cases)
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **classification,
        "configuration": {
            "angles_degrees_from_normal": list(ANGLES_DEGREES),
            "tangential_modes": list(TANGENTIAL_MODES),
            "minimum_resolved_wavelength_cells": (
                MINIMUM_RESOLVED_WAVELENGTH_CELLS
            ),
            "reflection_amplitude_ratio_limit": (
                REFLECTION_AMPLITUDE_RATIO_LIMIT
            ),
            "boundary_scheme": boundary_scheme,
            "short_domain_cells_xyz": [
                oblique.NX, oblique.SHORT_NY, oblique.NZ,
            ],
            "reference_domain_cells_xyz": [
                oblique.NX, oblique.REFERENCE_NY, oblique.NZ,
            ],
            "boundary_under_test": "y-min/y-max isothermal characteristic open",
            "tangential_boundary": "x periodic",
            "vertical_boundary": "z specular symmetry",
            "worker_processes": worker_count,
        },
        "cases": cases,
        "scope": (
            "linear_acoustic_component_gate_not_corner_or_plume_domain_validation"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    angles = list(record["configuration"]["angles_degrees_from_normal"])
    modes = list(record["configuration"]["tangential_modes"])
    ratios = np.full((len(modes), len(angles)), np.nan)
    wavelengths = np.full_like(ratios, np.nan)
    resolved = np.zeros_like(ratios, dtype=bool)
    for case in record["cases"]:
        row = modes.index(int(case["tangential_mode"]))
        column = angles.index(float(case["angle_degrees_from_normal"]))
        ratios[row, column] = float(case["reflection_amplitude_ratio"])
        wavelengths[row, column] = float(case["wavelength_cells"])
        resolved[row, column] = bool(case["resolved_for_gate"])

    figure, axis = plt.subplots(figsize=(9.4, 4.8), constrained_layout=True)
    image = axis.imshow(
        100.0 * ratios, cmap="viridis", vmin=0.0,
        vmax=max(5.0, float(np.nanmax(100.0 * ratios))), aspect="auto",
    )
    for row in range(len(modes)):
        for column in range(len(angles)):
            label = (
                f"{100.0 * ratios[row, column]:.2f}%\n"
                f"lambda={wavelengths[row, column]:.1f} dx"
            )
            if not resolved[row, column]:
                label += "\nexcluded"
            color = "white" if ratios[row, column] > 0.03 else "black"
            axis.text(column, row, label, ha="center", va="center", color=color)
    axis.set_xticks(range(len(angles)), [f"{angle:g}°" for angle in angles])
    axis.set_yticks(range(len(modes)), [f"mode {mode}" for mode in modes])
    axis.set_xlabel("Incidence angle from boundary normal")
    axis.set_ylabel("Tangential Fourier mode")
    axis.set_title(
        "Characteristic pressure-open reflection sweep\n"
        "5% gate; wavelength < 8 cells reported but excluded"
    )
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Reflected / incident amplitude (%)")
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--boundary-scheme",
        choices=(
            "characteristic_nonreflecting",
            "characteristic_multidimensional",
        ),
        default="characteristic_nonreflecting",
    )
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    record = run_audit(
        workers=args.workers,
        boundary_scheme=args.boundary_scheme,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot_audit(record, args.plot)
    worst = record["worst_in_scope_case"] or {}
    print(json.dumps({
        "status": record["status"],
        "in_scope_case_count": record["in_scope_case_count"],
        "excluded_underresolved_case_count": (
            record["excluded_underresolved_case_count"]
        ),
        "worst_in_scope_angle_degrees": worst.get(
            "angle_degrees_from_normal"
        ),
        "worst_in_scope_tangential_mode": worst.get("tangential_mode"),
        "worst_in_scope_reflection_amplitude_ratio": worst.get(
            "reflection_amplitude_ratio"
        ),
        "runtime_s": record["runtime_s"],
    }, indent=2, allow_nan=False))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
