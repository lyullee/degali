"""Assemble sealed LH2 downstream thermal, mechanical and buoyancy evidence."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"
FILES = {
    "cold": REF / "downstream_cold_envelope_2026-09-06.json",
    "mechanical": REF / "mechanical_thermal_scale_2026-09-06.json",
    "buoyancy": REF / "trial10_vertical_momentum_budget_2026-09-05.json",
    "tke": REF / "joint_tke_normal_budget_2026-09-06.json",
}
PROTOCOL = ROOT / "docs/prereg-downstream-residual-map.md"
DISTANCE_TOLERANCE_M = 0.011


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nearest(rows: list[dict], x: float, field: str, tolerance: float = DISTANCE_TOLERANCE_M):
    row = min(rows, key=lambda value: abs(float(value[field]) - x))
    difference = float(row[field]) - x
    if abs(difference) > tolerance:
        return None
    return row, difference


def thermal_classification(observations: dict, mechanical_temperature_scale: float) -> dict:
    ratios = [
        float(value["positive_enthalpy_gap_over_local_mean_ke"])
        for value in observations.values()
        if float(value["conditional_specific_enthalpy_gap_J_kg"]) > 0
    ]
    minimum_gap = float(observations["minimum"]["temperature_gap_K"])
    median_gap = float(observations["median"]["temperature_gap_K"])
    return {
        "mean_ke_primary_explanation_rejected": bool(ratios and min(ratios) > 10.0),
        "minimum_positive_enthalpy_to_mean_ke_ratio": min(ratios) if ratios else None,
        "accounted_mechanical_temperature_scale_below_1K": mechanical_temperature_scale < 1.0,
        "minimum_median_residual_sign_change": minimum_gap * median_gap < 0.0,
    }


def assemble(cold: dict, mechanical: dict, buoyancy: dict, tke: dict) -> dict:
    point_rows = [
        row for row in mechanical["point_scales"]
        if row["variant"] == "pressure_loss_density" and row["trial"] in (10, 23)
    ]
    budgets = {
        row["trial"]: row for row in mechanical["budgets"]
        if row["variant"] == "pressure_loss_density" and row["trial"] in (10, 23)
    }
    momentum = buoyancy["points"]
    rows = []
    for point in point_rows:
        trial, x = int(point["trial"]), float(point["x"])
        section_hit = nearest(budgets[trial]["sections"], x, "x_m")
        if section_hit is None:
            raise ValueError(f"no matching mechanical section for trial {trial}, x={x}")
        section, offset = section_hit
        observations = point["observations"]
        temperature_scale = float(section["handoff_plus_accounted_input_uniform_sensible_temperature_scale_K"])
        row = {
            "trial": trial,
            "x_m": x,
            "channel": point["channel"],
            "predicted_peak_temperature_K": float(point["exact_peak_temperature_K"]),
            "temperature_residual_K": {
                basis: float(value["temperature_gap_K"])
                for basis, value in observations.items()
            },
            "conditional_specific_enthalpy_gap_J_kg": {
                basis: float(value["conditional_specific_enthalpy_gap_J_kg"])
                for basis, value in observations.items()
            },
            "positive_enthalpy_gap_over_local_mean_ke": {
                basis: float(value["positive_enthalpy_gap_over_local_mean_ke"])
                for basis, value in observations.items()
            },
            "mechanical_section_offset_m": offset,
            "accounted_mechanical_input_fraction_of_initial_thermal": (
                float(section["cumulative_accounted_ambient_mechanical_input_W"])
                / abs(float(budgets[trial]["initial_thermal_flux_W"]))
            ),
            "accounted_mechanical_uniform_temperature_scale_K": temperature_scale,
            "classification": thermal_classification(observations, temperature_scale),
        }
        if trial == 10:
            force_hit = nearest(momentum, x, "x_m")
            if force_hit is not None:
                force, force_offset = force_hit
                row["trial10_buoyancy"] = {
                    "offset_m": force_offset,
                    "centre_z_m": float(force["candidate"]["z_m"]),
                    "vertical_momentum_N": float(force["candidate"]["vertical_momentum_N"]),
                    "buoyancy_N_per_m": float(force["candidate"]["buoyancy_N_per_m"]),
                    "cumulative_buoyancy_N": force["candidate"].get("cumulative_buoyancy_N"),
                }
        rows.append(row)

    tke_bounds = {}
    for trial in (10, 23):
        source = next(row for row in tke["rows"] if row["trial"] == trial)
        evaluation = next(row for row in source["evaluations"] if row["order"] == 16)
        normal = next(row for row in evaluation["normal_work_bounds"] if row["tke_fraction_cap"] == 0.12)
        kbound = next(row for row in evaluation["tke_bounds"] if row["normal_work_fraction_cap"] == 0.05)
        tke_bounds[str(trial)] = {
            "mean_kinetic_flux_W": float(evaluation["mean_kinetic_flux"]),
            "shear_bound_flux_W": float(evaluation["shear_bound_flux"]),
            "minimum_normal_work_fraction_at_tke_fraction_0.12": float(normal["normal_work_fraction_minimum"]),
            "minimum_tke_fraction_at_normal_work_fraction_0.05": float(kbound["tke_fraction_minimum"]),
            "not_an_initial_condition": True,
        }

    return {
        "rows": rows,
        "row_count": len(rows),
        "trials": [10, 23],
        "tke_normal_work_feasibility": tke_bounds,
        "trajectory_reintegrated": False,
        "coefficient_selected": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    data = {key: json.loads(path.read_text(encoding="utf-8")) for key, path in FILES.items()}
    result = assemble(**data)
    result.update({
        "completed": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": str(PROTOCOL.relative_to(ROOT)).replace("\\", "/"),
        "inputs_sha256": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
            for path in (*FILES.values(), PROTOCOL, Path(__file__))
        },
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(f"saved {args.output}: {result['row_count']} rows")


if __name__ == "__main__":
    main()

