"""Ablate the sensor-blind crosswind routing of TI021 wall-impact material."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
from tools.run_ti021_measured_wind_candidate import _group_metrics
from tools.sweep_ti021_vertical_trajectory import _base_parameters, _execute


def _sensor_arrays(report: dict) -> tuple[np.ndarray, np.ndarray]:
    predicted = np.asarray(
        [row["predicted_dose_vol_percent_s"] for row in report["sensor_rows"]],
        dtype=float,
    )
    observed = np.asarray(
        [row["observed_dose_vol_percent_s"] for row in report["sensor_rows"]],
        dtype=float,
    )
    return predicted, observed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline",
        type=Path,
        default=ROOT
        / "outputs/ti021-cavity-frame-cavity-currentcode5-600s-2026-10-17.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs/ti021-crosswind-path-ablation-2026-10-17.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    common = {
        "wall_impulse_elevated_conserved_fraction": 0.20,
        "wall_impulse_elevated_conserved_center_m": 3.5,
        "wall_impulse_elevated_conserved_sigma_m": 1.5,
        "wall_impulse_elevated_conserved_vertical_diffusivity_m2_s": 0.01,
        "cavity_backflow_frame": "cavity",
    }
    variants = {
        "inherited_plus_0p3_outward_0p1": {},
        "zero_uniform_outward_0p1": {
            "wall_impulse_elevated_conserved_crosswind_velocity_m_s": 0.0,
        },
        "minus_0p15_outward_0p1": {
            "wall_impulse_elevated_conserved_crosswind_velocity_m_s": -0.15,
        },
        "zero_uniform_outward_0p3": {
            "wall_impulse_elevated_conserved_crosswind_velocity_m_s": 0.0,
            "wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s": 0.3,
        },
        "zero_uniform_no_outward": {
            "wall_impulse_elevated_conserved_crosswind_velocity_m_s": 0.0,
            "wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s": 0.0,
        },
        "impact_head_plus_0p3_outward_0p1": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
        },
        "impact_head_plus_0p3_no_outward": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s": 0.0,
        },
        "impact_head_zero_uniform_no_outward": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_crosswind_velocity_m_s": 0.0,
            "wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s": 0.0,
        },
        "impact_head_fraction_0p08": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.08,
        },
        "impact_head_fraction_0p10": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.10,
        },
        "impact_head_fraction_0p12": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.12,
        },
        "impact_head_fraction_0p15": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.15,
        },
        "impact_head_fraction_0p18": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.18,
        },
        "impact_head_fraction_0p25": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.25,
        },
        "impact_head_fraction_0p30": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.30,
        },
        "impact_head_fraction_0p35": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.35,
        },
        "impact_head_fraction_0p40": {
            "wall_impulse_elevated_conserved_center_mode": "impact_head",
            "wall_impulse_elevated_conserved_fraction": 0.40,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    jobs = []
    result_paths: dict[str, Path] = {}
    for name, changes in variants.items():
        params = _base_parameters(baseline)
        params.update(common)
        params.update(changes)
        path = args.output.with_name(f"{args.output.stem}-{name}-600s.json")
        result_paths[name] = path
        if not (args.skip_existing and path.exists()):
            jobs.append((str(path), params))

    if args.workers == 1:
        for job in jobs:
            _execute(job)
    elif jobs:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(_execute, jobs))

    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    baseline_predicted, observed = _sensor_arrays(baseline)
    baseline_groups = _group_metrics(
        baseline_predicted, observed, dataset.dose_sensors
    )
    positive = observed > 0.0
    records = []
    for name, path in result_paths.items():
        report = json.loads(path.read_text(encoding="utf-8"))
        predicted, candidate_observed = _sensor_arrays(report)
        if not np.array_equal(candidate_observed, observed):
            raise RuntimeError(f"observation mismatch in {name}")
        groups = _group_metrics(predicted, observed, dataset.dose_sensors)
        comparison = {
            key: {
                "baseline_mape_percent": baseline_groups[key]["mape_percent"],
                "candidate_mape_percent": groups[key]["mape_percent"],
                "delta_percentage_points": (
                    groups[key]["mape_percent"]
                    - baseline_groups[key]["mape_percent"]
                ),
                "nondegraded": (
                    groups[key]["mape_percent"]
                    <= baseline_groups[key]["mape_percent"] + 1.0e-12
                ),
            }
            for key in baseline_groups
        }
        diagnostics = report["diagnostics"]
        records.append(
            {
                "name": name,
                "result": str(path.resolve()),
                "dose_mape_percent": float(
                    np.mean(np.abs(predicted[positive] / observed[positive] - 1.0))
                    * 100.0
                ),
                "runtime_including_output_s": report["runtime_including_output_s"],
                "maximum_relative_species_mass_residual": diagnostics[
                    "maximum_relative_species_mass_residual"
                ],
                "mixture_mass_residual_kg": diagnostics[
                    "mixture_mass_residual_kg"
                ],
                "momentum_conservation_gate": diagnostics[
                    "momentum_conservation_gate"
                ],
                "group_comparison": comparison,
                "all_groups_nondegraded": all(
                    item["nondegraded"] for item in comparison.values()
                ),
            }
        )
    summary = {
        "schema": "degali.ti021-crosswind-path-ablation.v1",
        "basis": (
            "sensor-blind wall-impact path ablation; common 20% conserved upper "
            "split with fixed inventory and vertical kernel"
        ),
        "baseline": str(args.baseline.resolve()),
        "baseline_dose_mape_percent": baseline["dose_mape_percent"],
        "variants": records,
    }
    args.output.write_text(
        json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8"
    )
    print(
        json.dumps(
            [
                {
                    "name": row["name"],
                    "mape": row["dose_mape_percent"],
                    "runtime": row["runtime_including_output_s"],
                    "all_groups_nondegraded": row["all_groups_nondegraded"],
                }
                for row in records
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
