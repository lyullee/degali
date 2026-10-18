"""Run and point-score a non-promotion AIJ Case J development comparison."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.validation.aij_case_j import read_aij_case_j
from degali.validation.aij_case_j_comparison import compare_aij_case_j
from degali.validation.aij_case_j_setup import build_aij_case_j_resolved_setup


DATA = ROOT / "reference" / "obstacle_wake" / "aij_case_j_zenodo_15447888"


def run(*, grid_spacing_m: float, duration_s: float, time_step_s: float) -> dict[str, object]:
    dataset = read_aij_case_j(DATA)
    setup = build_aij_case_j_resolved_setup(
        dataset,
        grid_spacing_m=grid_spacing_m,
        duration_s=duration_s,
        time_step_s=time_step_s,
    )
    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(setup.config)
    wall_runtime_s = time.perf_counter() - started
    comparison = compare_aij_case_j(
        dataset,
        result,
        coordinate_translation_model_minus_experiment_m=(
            setup.coordinate_translation_model_minus_experiment_m
        ),
    )
    return {
        "schema": "degali.aij-case-j-resolved-development-comparison.v1",
        "status": "withheld_development_comparison_not_validation",
        "benchmark": {
            "record_doi": dataset.record_doi,
            "approach_sha256": dataset.approach_sha256,
            "surface_temperature_sha256": dataset.surface_temperature_sha256,
            "observations_sha256": dataset.observations_sha256,
            "observation_count": len(dataset.observations),
            "published_observations_are_time_averaged": True,
            "model_statistic": comparison.model_statistic,
            "statistic_match": False,
            "statistic_mismatch_note": (
                "experimental 120-s statistical means are not equivalent to a short finite transient average"
            ),
            "tracer": "ethylene",
            "hydrogen_validation": False,
        },
        "setup": {
            "grid_spacing_m": setup.grid_spacing_m,
            "cells": list(setup.config.cells),
            "domain_m": list(setup.config.domain_m),
            "model_minus_experiment_translation_m": list(
                setup.coordinate_translation_model_minus_experiment_m
            ),
            "source_ethylene_mass_rate_kg_s": setup.source_ethylene_mass_rate_kg_s,
            "source_total_to_ethylene_mass_ratio": setup.source_total_to_ethylene_mass_ratio,
            "physical_air_kinematic_viscosity_m2_s": setup.physical_air_kinematic_viscosity_m2_s,
            "property_derived_reynolds_number": setup.reynolds_from_properties,
            "reference_reynolds_number": setup.reference_reynolds_number,
            "source_hole_resolved_cells_across": setup.source_hole_resolved_cells_across,
            "scope_warnings": list(setup.scope_warnings),
        },
        "comparison": {
            "status": comparison.status,
            "point_count": comparison.metrics["point_count"],
            "skipped_point_ids": list(comparison.skipped_point_ids),
            "skip_reasons": comparison.skip_reasons,
            "metrics": comparison.metrics,
            "concentration_unit": "ethylene vol% (absolute errors are vol%-points)",
            "calibration_used": False,
        },
        "runtime": {
            "wall_runtime_s": wall_runtime_s,
            "solver_reported_runtime_s": result.diagnostics["total_runtime_s"],
            "simulated_duration_s": result.time_s[-1].item(),
            "transient_average_start_s": result.time_average_start_s,
            "transient_average_sample_count": result.time_average_sample_count,
        },
        "conservation": {
            "maximum_tracked_species_mass_residual_kg": result.maximum_mass_residual_kg,
            "maximum_relative_mass_residual": result.diagnostics[
                "maximum_relative_mass_residual"
            ],
            "maximum_total_mass_projection_residual_kg_s": result.diagnostics[
                "maximum_variable_density_global_residual_kg_s"
            ],
            "maximum_energy_residual_j": result.maximum_energy_residual_j,
            "maximum_relative_energy_residual": result.diagnostics[
                "maximum_relative_energy_residual"
            ],
        },
        "solver_diagnostics": result.diagnostics,
        "pointwise": [asdict(point) for point in comparison.points],
    }


def _figure(record: dict[str, object], path: Path) -> None:
    points = record["pointwise"]
    if not points:
        return
    figure, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    comparisons = (
        ("u", "measured_u_m_s", "predicted_u_m_s", "U (m/s)"),
        ("v", "measured_v_m_s", "predicted_v_m_s", "V (m/s)"),
        ("w", "measured_w_m_s", "predicted_w_m_s", "W (m/s)"),
        (
            "ethylene", "measured_ethylene_vol_percent",
            "predicted_ethylene_vol_percent", "Ethylene (vol%)",
        ),
    )
    for axis, (_key, measured, predicted, label) in zip(axes.flat, comparisons):
        x = np.asarray([item[measured] for item in points], dtype=float)
        y = np.asarray([item[predicted] for item in points], dtype=float)
        axis.scatter(x, y, s=15, alpha=0.75)
        lower = min(float(np.min(x)), float(np.min(y)))
        upper = max(float(np.max(x)), float(np.max(y)))
        if upper == lower:
            upper = lower + 1.0
        axis.plot((lower, upper), (lower, upper), "k--", linewidth=1)
        axis.set(xlabel=f"Measured {label}", ylabel=f"Model {label}")
        axis.grid(alpha=0.25)
    figure.suptitle(
        "AIJ Case J — finite-window transient mean (not validation)",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-spacing-m", type=float, default=0.02)
    parser.add_argument("--duration-s", type=float, default=0.60)
    parser.add_argument("--time-step-s", type=float, default=5.0e-4)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "outputs" / "aij-case-j-resolved-development-2026-10-12-v1",
    )
    args = parser.parse_args()
    record = run(
        grid_spacing_m=args.grid_spacing_m,
        duration_s=args.duration_s,
        time_step_s=args.time_step_s,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "aij-case-j-resolved-development-comparison.json"
    figure_path = args.output_dir / "aij-case-j-resolved-development-scatter.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _figure(record, figure_path)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "points": record["comparison"]["point_count"],
        "wall_runtime_s": record["runtime"]["wall_runtime_s"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
