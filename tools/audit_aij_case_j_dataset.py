"""Build a hash-pinned structural audit for the official AIJ Case J data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degali.validation.aij_case_j import read_aij_case_j


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reference" / "obstacle_wake" / "aij_case_j_zenodo_15447888"
COMPARISON = ROOT / "outputs" / "aij-case-j-resolved-development-0p60s-dx0p02-dt0p001-2026-10-12-v2" / "aij-case-j-resolved-development-comparison.json"


def build() -> dict[str, object]:
    dataset = read_aij_case_j(DATA)
    concentrations = [item.concentration_vol_percent for item in dataset.observations]
    temperatures = [item.temperature_k for item in dataset.observations]
    comparison = json.loads(COMPARISON.read_text(encoding="utf-8")) if COMPARISON.is_file() else None
    comparison_metrics = comparison.get("comparison", {}).get("metrics", {}) if comparison else {}
    comparison_details = {
        "status": "scored_development_component_comparison_validation_withheld" if comparison else "withheld_not_run",
        "pointwise_scoring_implementation_available": True,
        "calibration_used": False,
        "hydrogen_validation": False,
    }
    if comparison:
        pointwise = comparison.get("pointwise", [])
        concentration_rows = [row for row in pointwise if "measured_ethylene_vol_percent" in row]
        concentration_by_x = []
        for x_m in sorted({row["x_experiment_m"] for row in concentration_rows}):
            station = [row for row in concentration_rows if row["x_experiment_m"] == x_m]
            measured = [row["measured_ethylene_vol_percent"] for row in station]
            predicted = [row["predicted_ethylene_vol_percent"] for row in station]
            errors = [model - observed for model, observed in zip(predicted, measured)]
            concentration_by_x.append({
                "x_experiment_m": x_m,
                "point_count": len(station),
                "measured_max_vol_percent": max(measured),
                "predicted_max_vol_percent": max(predicted),
                "bias_vol_percentage_points": sum(errors) / len(errors),
                "mae_vol_percentage_points": sum(abs(error) for error in errors) / len(errors),
                "rmse_vol_percentage_points": (sum(error * error for error in errors) / len(errors)) ** 0.5,
            })
        comparison_details.update({
            "comparison_artifact": str(COMPARISON.relative_to(ROOT)).replace("\\", "/"),
            "model_status": comparison.get("status"),
            "point_count": comparison.get("comparison", {}).get("point_count"),
            "observation_count": comparison.get("benchmark", {}).get("observation_count"),
            "skipped_point_ids": comparison.get("comparison", {}).get("skipped_point_ids", []),
            "ethylene_vol_percent_mae_vol_percentage_points": comparison_metrics.get("ethylene_vol_percent_mae"),
            "ethylene_vol_percent_rmse_vol_percentage_points": comparison_metrics.get("ethylene_vol_percent_rmse"),
            "ethylene_vol_percent_correlation": comparison_metrics.get("ethylene_vol_percent_correlation"),
            "finite_transient_average_start_s": comparison.get("runtime", {}).get("transient_average_start_s"),
            "finite_transient_average_samples": comparison.get("runtime", {}).get("transient_average_sample_count"),
            "experimental_statistic_match": comparison.get("benchmark", {}).get("statistic_match"),
            "concentration_by_downstream_station": concentration_by_x,
            "reason": "Development score only: the 0.6-s finite transient mean is not statistically equivalent to the experiment's 120-s mean; coarse grid, subgrid source, and spatially averaged wall temperature also withhold validation.",
        })
    else:
        comparison_details["reason"] = (
            "Case J is an ethylene wind-tunnel component benchmark. A matched-condition "
            "DEGALI score is not present; the scorer requires passive ethylene concentration."
        )
    return {
        "schema": "degali.aij-case-j-dataset-audit.v1",
        "status": comparison_details["status"],
        "record_doi": dataset.record_doi,
        "dataset_scope": {
            "tracer": dataset.tracer,
            "building_m": [dataset.building_width_m, dataset.building_depth_m, dataset.building_height_m],
            "unstable_stratification": True,
            "bulk_richardson_number": dataset.bulk_richardson_number,
            "published_data_are_transient": False,
        },
        "row_counts": {
            "approach_flow": len(dataset.approach_points),
            "building_surface_temperature": len(dataset.surface_points),
            "wake_observations": len(dataset.observations),
        },
        "sha256": {
            "approach_flow": dataset.approach_sha256,
            "building_surface_temperature": dataset.surface_temperature_sha256,
            "wake_observations": dataset.observations_sha256,
        },
        "wake_observation_ranges": {
            "x_m": [min(item.x_m for item in dataset.observations), max(item.x_m for item in dataset.observations)],
            "y_m": [min(item.y_m for item in dataset.observations), max(item.y_m for item in dataset.observations)],
            "z_m": [min(item.z_m for item in dataset.observations), max(item.z_m for item in dataset.observations)],
            "temperature_k": [min(temperatures), max(temperatures)],
            "ethylene_vol_percent": [min(concentrations), max(concentrations)],
        },
        "vol_percent_conversion": "C(ppm) / 10000; ideal gas mole fraction basis",
        "model_comparison": comparison_details,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.output} rows={record['row_counts']} comparison={record['model_comparison']['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
