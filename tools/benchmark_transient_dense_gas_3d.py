"""Grid/time convergence benchmark for the 3-D transient dense-gas solver."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.transient_dense_gas_3d import (
    DenseGas3DReceptor,
    TransientDenseGas3DConfig,
    WindHistory,
    solve_transient_dense_gas_3d,
)


def _case(cells: tuple[int, int, int], time_step_s: float) -> TransientDenseGas3DConfig:
    return TransientDenseGas3DConfig(
        domain_m=(20.0, 20.0, 10.0), cells=cells,
        duration_s=1.0, time_step_s=time_step_s,
        source_schedule=SourceRateSchedule(
            (0.0, 0.3, 0.6, 1.0), (1.5, 1.0, 0.6, 0.0),
            source_id="transient-3d-error-study",
        ),
        wind_history=WindHistory(
            (0.0, 0.4, 0.8, 1.0), (1.2, 2.0, 1.5, 1.0),
            (270.0, 270.0, 180.0, 180.0),
        ),
        source_position_m=(3.0, 10.0, 5.0), source_sigma_m=1.2,
        source_density_kg_m3=10.0, source_h2_mass_fraction=0.9,
        ambient_density_kg_m3=1.2, diffusivity_m2_s=0.08,
        buoyancy_length_m=1.0, wind_relaxation_time_s=0.25,
    )


def _receptors() -> tuple[DenseGas3DReceptor, ...]:
    return tuple(DenseGas3DReceptor(label, position) for label, position in {
        "near_axis": (4.0, 10.0, 5.0),
        "mid_axis": (7.0, 10.0, 4.0),
        "crosswind": (7.0, 12.0, 5.0),
        "low": (10.0, 10.0, 3.0),
    }.items())


def _trace_metrics(reference, candidate, receptor: DenseGas3DReceptor) -> dict[str, float]:
    ref = reference.receptor_traces_kg_m3[receptor.label]
    value = np.interp(reference.time_s, candidate.time_s,
                      candidate.receptor_traces_kg_m3[receptor.label])
    peak_ref = float(np.max(ref))
    peak_value = float(np.max(value))
    dose_ref = float(np.trapezoid(ref, reference.time_s))
    dose_value = float(np.trapezoid(value, reference.time_s))
    threshold = 0.1 * peak_ref
    reference_arrival = float(
        reference.time_s[np.argmax(ref >= threshold)]
        if np.any(ref >= threshold) else math.nan
    )
    candidate_arrival = float(
        reference.time_s[np.argmax(value >= threshold)]
        if np.any(value >= threshold) else math.nan
    )
    return {
        "mae_kg_m3": float(np.mean(np.abs(value - ref))),
        "rmse_kg_m3": float(np.sqrt(np.mean((value - ref) ** 2))),
        "max_absolute_error_kg_m3": float(np.max(np.abs(value - ref))),
        "reference_peak_kg_m3": peak_ref,
        "peak_relative_error_pct": (
            100.0 * (peak_value - peak_ref) / peak_ref
            if peak_ref > 1.0e-15 else math.nan
        ),
        "dose_relative_error_pct": (
            100.0 * (dose_value - dose_ref) / dose_ref
            if dose_ref > 1.0e-15 else math.nan
        ),
        "arrival_error_s": candidate_arrival - reference_arrival,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    receptors = _receptors()
    levels = (
        ("coarse", (12, 12, 8), 0.02),
        ("medium", (24, 24, 16), 0.01),
        ("fine", (36, 36, 24), 0.005),
        ("reference", (48, 48, 32), 0.0025),
    )
    results = {}
    for name, cells, step in levels:
        result = solve_transient_dense_gas_3d(
            _case(cells, step), receptors=receptors,
        )
        results[name] = result
    reference = results["reference"]
    report: dict[str, object] = {
        "schema": "degali.transient-dense-gas-3d-error-study.v1",
        "units": {
            "concentration": "kg H2/m3",
            "error": "absolute kg H2/m3; relative peak/dose in percent",
        },
        "reference": {"cells": [48, 48, 32], "time_step_s": 0.0025},
        "levels": {},
    }
    aggregate_errors = []
    for name, cells, step in levels:
        result = results[name]
        level = {
            "cells": list(cells),
            "time_step_s": step,
            "steps": len(result.time_s) - 1,
            "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        }
        if name != "reference":
            metrics = {
                receptor.label: _trace_metrics(reference, result, receptor)
                for receptor in receptors
            }
            level["receptor_metrics_vs_reference"] = metrics
            all_errors = []
            all_reference = []
            for receptor in receptors:
                candidate = np.interp(
                    reference.time_s, result.time_s,
                    result.receptor_traces_kg_m3[receptor.label],
                )
                all_errors.extend(candidate - reference.receptor_traces_kg_m3[receptor.label])
                all_reference.extend(reference.receptor_traces_kg_m3[receptor.label])
            errors = np.asarray(all_errors)
            aggregate_errors.append((name, float(np.sqrt(np.mean(errors ** 2)))))
            level["aggregate_mae_kg_m3"] = float(np.mean(np.abs(errors)))
            level["aggregate_rmse_kg_m3"] = float(np.sqrt(np.mean(errors ** 2)))
            level["aggregate_max_absolute_error_kg_m3"] = float(np.max(np.abs(errors)))
        report["levels"][name] = level
    report["apparent_order_from_aggregate_rmse"] = {
        f"{left}_to_{right}": math.log(left_error / right_error, 2)
        for (left, left_error), (right, right_error) in zip(
            aggregate_errors, aggregate_errors[1:]
        )
    }
    report["interpretation"] = (
        "Mass residual is the solver conservation error. Concentration metrics "
        "are numerical resolution errors against the declared high-resolution "
        "reference, not experimental validation errors. Relative errors at the "
        "low receptor can be large because the reference concentration is near zero."
    )
    (args.output_dir / "error_study.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# 3-D transient dense-gas numerical concentration error study",
        "",
        "Reference: 48×48×32 cells, dt=0.0025 s; concentration unit kg H2/m³.",
        "The reference is a numerical solution, not measured validation data.",
        "",
        "| level | grid | dt (s) | aggregate MAE | aggregate RMSE | max abs | max mass residual |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, cells, step in levels:
        level = report["levels"][name]
        lines.append(
            f"| {name} | {'×'.join(map(str, cells))} | {step:g} | "
            f"{level.get('aggregate_mae_kg_m3', 0.0):.6g} | "
            f"{level.get('aggregate_rmse_kg_m3', 0.0):.6g} | "
            f"{level.get('aggregate_max_absolute_error_kg_m3', 0.0):.6g} | "
            f"{level['maximum_mass_residual_kg']:.3e} |"
        )
    lines.extend([
        "",
        "Interpretation: the mass residual is near machine precision. The remaining "
        "concentration differences are discretization/smearing errors, not sensor "
        "or field-model errors.",
    ])
    (args.output_dir / "error_study.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
