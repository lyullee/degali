"""Compare concentration-blind CG and AMG-CG TI021 preflight records."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(path: Path, expected_solver: str) -> dict[str, object]:
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("schema") != "degali.thorney-ti021-full-scale-preflight.v1":
        raise ValueError("unsupported TI021 preflight schema")
    numerics = record.get("numerics", {})
    if numerics.get("projection_solver_used") != expected_solver:
        raise ValueError(
            f"expected {expected_solver!r} record, got "
            f"{numerics.get('projection_solver_used')!r}"
        )
    return record


def compare_records(
    cg: dict[str, object], amg: dict[str, object],
) -> dict[str, object]:
    cg_config = cg["configuration"]
    amg_config = amg["configuration"]
    comparison_keys = (
        "cell_size_m", "domain_m", "cells", "cell_count", "active_steps",
        "active_time_step_s", "active_duration_s", "wind_speed_m_s",
        "streamwise_wind_m_s", "crosswind_m_s",
        "carrier_flow_initialization", "turbulence_closure",
    )
    matched_configuration = all(
        cg_config[key] == amg_config[key] for key in comparison_keys
    )
    cg_numerics = cg["numerics"]
    amg_numerics = amg["numerics"]
    cg_runtime = float(cg["runtime"]["total_s"])
    amg_runtime = float(amg["runtime"]["total_s"])
    cg_iterations = int(
        cg_numerics["projection_maximum_linear_iterations_per_step"]
    )
    amg_iterations = int(
        amg_numerics["projection_maximum_linear_iterations_per_step"]
    )
    gates = {
        "configuration_matched": matched_configuration,
        "both_preflight_numerical_gates_pass": bool(
            cg["all_numerical_gates_pass"]
            and amg["all_numerical_gates_pass"]
        ),
        "maximum_mole_fraction_matches_1e_10": math.isclose(
            float(cg_numerics["maximum_final_species_mole_fraction"]),
            float(amg_numerics["maximum_final_species_mole_fraction"]),
            rel_tol=0.0,
            abs_tol=1.0e-10,
        ),
        "density_departure_matches_1e_10": math.isclose(
            float(cg_numerics["maximum_relative_density_departure"]),
            float(amg_numerics["maximum_relative_density_departure"]),
            rel_tol=0.0,
            abs_tol=1.0e-10,
        ),
        "amg_continuity_residual_at_most_1e_7": (
            float(amg_numerics["maximum_continuity_residual_kg_m3_s"])
            <= 1.0e-7
        ),
        "amg_mass_residual_at_most_1e_6": (
            float(amg_numerics["maximum_relative_species_mass_residual"])
            <= 1.0e-6
        ),
        "amg_iteration_reduction_at_least_10x": (
            cg_iterations / max(amg_iterations, 1) >= 10.0
        ),
        "amg_runtime_speedup_at_least_1p5x": (
            cg_runtime / amg_runtime >= 1.5
        ),
    }
    return {
        "schema": "degali.thorney-ti021-projection-acceleration-audit.v1",
        "status": (
            "accepted_projection_acceleration_component"
            if all(gates.values())
            else "withheld_projection_acceleration_component"
        ),
        "scope": (
            "linear-solver equivalence and runtime only; this does not validate "
            "turbulence, dispersion, dose or the long transient"
        ),
        "configuration": {
            key: cg_config[key] for key in comparison_keys
        },
        "comparison": {
            "cg": {
                "runtime_s": cg_runtime,
                "projection_iterations": cg_iterations,
                "maximum_continuity_residual_kg_m3_s": cg_numerics[
                    "maximum_continuity_residual_kg_m3_s"
                ],
            },
            "amg_cg": {
                "runtime_s": amg_runtime,
                "projection_iterations": amg_iterations,
                "maximum_continuity_residual_kg_m3_s": amg_numerics[
                    "maximum_continuity_residual_kg_m3_s"
                ],
                "preconditioner": amg_numerics["projection_preconditioner"],
                "preconditioner_level_count": amg_numerics[
                    "projection_preconditioner_level_count"
                ],
                "preconditioner_operator_complexity": amg_numerics[
                    "projection_preconditioner_operator_complexity"
                ],
                "pyamg_version": amg_numerics["pyamg_version"],
            },
            "runtime_speedup": cg_runtime / amg_runtime,
            "projection_iteration_reduction_factor": (
                cg_iterations / max(amg_iterations, 1)
            ),
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cg-record", type=Path, required=True)
    parser.add_argument("--amg-record", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    cg_path = args.cg_record.resolve()
    amg_path = args.amg_record.resolve()
    cg = _load(cg_path, "cg")
    amg = _load(amg_path, "amg_cg")
    record = compare_records(cg, amg)
    record["provenance"] = {
        "cg_record": {"path": str(cg_path), "sha256": _sha256(cg_path)},
        "amg_record": {"path": str(amg_path), "sha256": _sha256(amg_path)},
    }
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-projection-acceleration-audit.json"
    plot_path = args.output_directory / "ti021-projection-acceleration-audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    comparison = record["comparison"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=True)
    axes[0].bar(
        ("Jacobi-CG", "AMG-CG"),
        (comparison["cg"]["runtime_s"], comparison["amg_cg"]["runtime_s"]),
        color=("0.55", "tab:blue"),
    )
    axes[0].set_ylabel("wall time (s)")
    axes[0].set_title(f"runtime: {comparison['runtime_speedup']:.2f}x faster")
    axes[1].bar(
        ("Jacobi-CG", "AMG-CG"),
        (
            comparison["cg"]["projection_iterations"],
            comparison["amg_cg"]["projection_iterations"],
        ),
        color=("0.55", "tab:blue"),
    )
    axes[1].set_ylabel("summed CG iterations / coupled step")
    axes[1].set_title(
        f"iterations: {comparison['projection_iteration_reduction_factor']:.1f}x fewer"
    )
    fig.suptitle("TI021 dx=2 m projection acceleration (same numerical case)")
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "comparison": comparison,
        "gates": record["gates"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
