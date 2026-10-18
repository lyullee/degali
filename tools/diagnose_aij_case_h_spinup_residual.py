"""Record the effect of extending AIJ Case H resolved-flow spin-up.

This is a solver-convergence diagnostic, not an experimental validation run.
It preserves the established medium-grid AIJ configuration and lengthens only
the D3Q19 flow integration so that a residual plateau is not mistaken for a
short-horizon startup transient.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

from audit_aij_resolved_3d_grid import _config
from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.validation.aij_case_h import read_aij_case_h


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flow-steps", type=int, default=4500)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/aij-resolved-spinup-residual-2026-10-12.json"),
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.flow_steps < 1500:
        raise SystemExit("--flow-steps must be at least 1500 for comparison")

    dataset = read_aij_case_h(Path("reference/aij_case_h"))
    base = _config(dataset, 0.020)
    config = replace(
        base,
        flow_steps=args.flow_steps,
        minimum_flow_steps=args.flow_steps,
    )
    result = solve_resolved_obstacle_flow_3d(config)
    residual_history = [
        {"step": int(step), "relative_velocity_update": float(residual)}
        for step, residual in result.flow_residual_history
    ]
    residual_by_step = {item["step"]: item["relative_velocity_update"] for item in residual_history}
    late_start_step = min(2000, args.flow_steps)
    late_start = residual_by_step.get(late_start_step)
    final_residual = float(result.diagnostics["final_flow_residual"])
    late_change = (
        None if late_start is None
        else abs(final_residual - late_start) / max(abs(late_start), abs(final_residual), 1.0e-15)
    )
    record = {
        "schema": "degali.aij-case-h-resolved-spinup-diagnostic.v1",
        "status": "diagnostic_only_not_validation_promotion",
        "dataset": {
            "record_doi": dataset.record_doi,
            "approach_flow_sha256": dataset.approach_flow_sha256,
            "result_sha256": dataset.result_sha256,
        },
        "configuration": {
            "cell_size_m": 0.020,
            "cells": list(config.cells),
            "flow_steps": args.flow_steps,
            "minimum_flow_steps": args.flow_steps,
            "physical_flow_time_s": (
                float(args.flow_steps) * config.lattice_inlet_speed
                * (config.domain_m[0] / config.cells[0])
                / config.wind_speed_m_s
            ),
        },
        "residual_summary": {
            "residual_at_1500_steps": residual_by_step.get(1500),
            "residual_at_2000_steps": residual_by_step.get(2000),
            "final_residual": final_residual,
            "late_window_relative_change_from_2000_steps": late_change,
            "flow_converged": bool(result.diagnostics["flow_converged"]),
            "maximum_absolute_divergence_s_1": float(
                result.diagnostics["maximum_absolute_divergence_s_1"]
            ),
            "minimum_lee_streamwise_velocity_m_s": float(
                result.diagnostics["minimum_lee_streamwise_velocity_m_s"]
            ),
            "recirculating_lee_cell_count": int(
                result.diagnostics["recirculating_lee_cell_count"]
            ),
        },
        "runtime_s": float(result.diagnostics["total_runtime_s"]),
        "residual_history": residual_history,
        "interpretation": (
            "The instantaneous velocity-update residual is not a continuity "
            "residual. A plateau alone does not distinguish physical wake "
            "unsteadiness from numerical non-convergence; use the separately "
            "declared time-averaged/statistical convergence gate for mean-flow "
            "validation."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(args.output.resolve()),
        "flow_steps": args.flow_steps,
        **record["residual_summary"],
        "runtime_s": record["runtime_s"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
