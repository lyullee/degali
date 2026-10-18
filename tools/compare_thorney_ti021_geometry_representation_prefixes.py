"""Compare matched TI021 prefixes from voxel-wall and zero-thickness walls.

The output is a numerical geometry-representation sensitivity report, not an
independent validation. Inputs must have identical grid, time, source and
receptor setup and independently pass their materialized-prefix gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


SCHEMA = "degali.thorney-ti021-geometry-representation-prefix-comparison.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_prefix(root: Path) -> tuple[dict, dict, dict, dict, dict[str, str]]:
    prefix_path = root / "ti021-checkpoint-prefix.json"
    record_path = root / "ti021-inner-transient.json"
    field_path = root / "ti021-inner-final-field.npz"
    trace_path = root / "ti021-inner-receptor-traces.npz"
    for path in (prefix_path, record_path, field_path, trace_path):
        if not path.is_file():
            raise FileNotFoundError(f"required prefix artifact is missing: {path}")
    prefix = json.loads(prefix_path.read_text(encoding="utf-8"))
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if prefix.get("schema") != "degali.thorney-ti021-checkpoint-prefix.v1":
        raise ValueError(f"unsupported checkpoint-prefix schema in {root}")
    if prefix.get("all_gates_pass") is not True:
        raise ValueError(f"checkpoint prefix has failed gates: {root}")
    source_journal = Path(prefix["source_checkpoint_journal"])
    if not source_journal.is_file():
        raise FileNotFoundError(f"source checkpoint journal is missing: {source_journal}")
    journal = json.loads(source_journal.read_text(encoding="utf-8"))
    configuration = journal["configuration"]
    mode = configuration.get("circular_arc_representation", "finite_thickness_voxel")
    if mode not in {"finite_thickness_voxel", "zero_thickness_surface"}:
        raise ValueError(f"unsupported obstacle representation {mode!r}")
    with np.load(field_path, allow_pickle=False) as field_file:
        field = {key: field_file[key].copy() for key in field_file.files}
    with np.load(trace_path, allow_pickle=False) as trace_file:
        traces = {key: trace_file[key].copy() for key in trace_file.files}
    hashes = {
        "prefix_manifest_sha256": _sha256(prefix_path),
        "transient_record_sha256": _sha256(record_path),
        "field_npz_sha256": _sha256(field_path),
        "receptor_traces_sha256": _sha256(trace_path),
        "source_journal_sha256_at_prefix_extraction": str(
            prefix["source_checkpoint_journal_sha256_at_extraction"]
        ),
    }
    return prefix, record, {"mode": mode, **configuration}, {"field": field, "traces": traces}, hashes


def _first_threshold_arrivals(times: np.ndarray, values: np.ndarray, threshold: float) -> list[float | None]:
    arrivals: list[float | None] = []
    for series in values:
        indices = np.flatnonzero(series >= threshold)
        arrivals.append(float(times[indices[0]]) if indices.size else None)
    return arrivals


def _scope_warning(stop_time_s: float) -> str:
    return (
        f"This compares two wall regularisations over the same {stop_time_s:g}-s "
        "computed window. It is not an experimental validation or evidence that "
        "the obstacle representations are equivalent outside this window."
    )


def compare_prefixes(voxel_root: Path, surface_root: Path) -> dict:
    voxel_prefix, voxel_record, voxel_config, voxel_data, voxel_hashes = _load_prefix(voxel_root)
    surface_prefix, surface_record, surface_config, surface_data, surface_hashes = _load_prefix(surface_root)
    if voxel_config["mode"] != "finite_thickness_voxel":
        raise ValueError("first input must be a finite_thickness_voxel prefix")
    if surface_config["mode"] != "zero_thickness_surface":
        raise ValueError("second input must be a zero_thickness_surface prefix")
    common_config_keys = (
        "duration_s", "time_step_s", "cell_size_m", "corrector_tolerance",
        "cells", "workbook_sha256", "property_audit_sha256",
    )
    for key in common_config_keys:
        if voxel_config.get(key) != surface_config.get(key):
            raise ValueError(f"prefix configurations differ for {key}")
    voxel_setup = voxel_record["configuration"]["receptors"]
    surface_setup = surface_record["configuration"]["receptors"]
    if voxel_setup != surface_setup:
        raise ValueError("prefix receptor coordinates or ordering differ")

    vf, sf = voxel_data["field"], surface_data["field"]
    vt, st = voxel_data["traces"], surface_data["traces"]
    for axis in ("x_m", "y_m", "z_m"):
        if not np.array_equal(vf[axis], sf[axis]):
            raise ValueError(f"field coordinate axis {axis} differs")
    if not np.array_equal(vt["time_s"], st["time_s"]):
        raise ValueError("receptor trace time axes differ")
    if not np.array_equal(vt["labels"], st["labels"]):
        raise ValueError("receptor labels or ordering differ")

    voxel_solid = np.asarray(vf["solid_mask"], dtype=bool)
    surface_solid = np.asarray(sf["solid_mask"], dtype=bool)
    common_fluid = ~(voxel_solid | surface_solid)
    if not np.any(common_fluid):
        raise ValueError("the two geometries have no common fluid cells")
    voxel_final = np.asarray(vf["concentration_vol_percent"], dtype=float)
    surface_final = np.asarray(sf["concentration_vol_percent"], dtype=float)
    voxel_max = np.asarray(vf["maximum_concentration_vol_percent"], dtype=float)
    surface_max = np.asarray(sf["maximum_concentration_vol_percent"], dtype=float)
    final_delta = surface_final[common_fluid] - voxel_final[common_fluid]
    max_delta = surface_max[common_fluid] - voxel_max[common_fluid]
    voxel_trace = np.asarray(vt["concentration_vol_percent"], dtype=float)
    surface_trace = np.asarray(st["concentration_vol_percent"], dtype=float)
    trace_delta = surface_trace - voxel_trace
    cell_volume = float(voxel_config["cell_size_m"]) ** 3

    voxel_arrivals = _first_threshold_arrivals(vt["time_s"], voxel_trace, 4.0)
    surface_arrivals = _first_threshold_arrivals(st["time_s"], surface_trace, 4.0)
    comparable_arrivals = [
        abs(a - b) for a, b in zip(voxel_arrivals, surface_arrivals)
        if a is not None and b is not None
    ]
    return {
        "schema": SCHEMA,
        "status": "matched_computed_window_sensitivity_only",
        "scope_warning": _scope_warning(float(vt["time_s"][-1])),
        "shared_configuration": {key: voxel_config[key] for key in common_config_keys},
        "stop_time_s": float(vt["time_s"][-1]),
        "receptor_count": int(vt["labels"].size),
        "frame_count": int(vt["time_s"].size),
        "geometry_cell_counts": {
            "finite_thickness_voxel_solid_cells": int(voxel_solid.sum()),
            "zero_thickness_surface_solid_cells": int(surface_solid.sum()),
            "common_fluid_cells_for_field_error": int(common_fluid.sum()),
        },
        "maximum_concentration_vol_percent": {
            "finite_thickness_voxel": float(np.max(voxel_final[~voxel_solid])),
            "zero_thickness_surface": float(np.max(surface_final[~surface_solid])),
        },
        "common_fluid_field_difference_vol_percentage_points": {
            "endpoint_mae": float(np.mean(np.abs(final_delta))),
            "endpoint_rmse": float(np.sqrt(np.mean(final_delta**2))),
            "endpoint_max_abs": float(np.max(np.abs(final_delta))),
            "maximum_envelope_rmse": float(np.sqrt(np.mean(max_delta**2))),
        },
        "four_vol_percent_volume_m3": {
            "finite_thickness_voxel_endpoint": float(
                np.count_nonzero((voxel_final >= 4.0) & ~voxel_solid) * cell_volume
            ),
            "zero_thickness_surface_endpoint": float(
                np.count_nonzero((surface_final >= 4.0) & ~surface_solid) * cell_volume
            ),
            "finite_thickness_voxel_ever": float(
                np.count_nonzero((voxel_max >= 4.0) & ~voxel_solid) * cell_volume
            ),
            "zero_thickness_surface_ever": float(
                np.count_nonzero((surface_max >= 4.0) & ~surface_solid) * cell_volume
            ),
        },
        "receptor_difference_vol_percentage_points": {
            "pooled_rmse": float(np.sqrt(np.mean(trace_delta**2))),
            "median_per_sensor_rmse": float(np.median(np.sqrt(np.mean(trace_delta**2, axis=1)))),
            "maximum_per_sensor_rmse": float(np.max(np.sqrt(np.mean(trace_delta**2, axis=1)))),
            "sensors_reaching_4_vol_percent": {
                "finite_thickness_voxel": sum(item is not None for item in voxel_arrivals),
                "zero_thickness_surface": sum(item is not None for item in surface_arrivals),
            },
            "median_absolute_arrival_difference_s": (
                None if not comparable_arrivals else float(np.median(comparable_arrivals))
            ),
        },
        "input_artifact_sha256": {
            "finite_thickness_voxel": voxel_hashes,
            "zero_thickness_surface": surface_hashes,
        },
        "prefix_journal_hashes_at_extraction": {
            "finite_thickness_voxel": voxel_prefix["source_checkpoint_journal_sha256_at_extraction"],
            "zero_thickness_surface": surface_prefix["source_checkpoint_journal_sha256_at_extraction"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--voxel-prefix", type=Path, required=True)
    parser.add_argument("--surface-prefix", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path,
        default=Path("outputs/thorney-ti021-geometry-representation-prefix-comparison-2026-10-12.json"),
    )
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {args.output}")
    record = compare_prefixes(args.voxel_prefix, args.surface_prefix)
    record["runner_script_sha256"] = _sha256(Path(__file__).resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(args.output.resolve()),
        "status": record["status"],
        "stop_time_s": record["stop_time_s"],
        "common_fluid_field_rmse": record["common_fluid_field_difference_vol_percentage_points"]["endpoint_rmse"],
        "pooled_receptor_rmse": record["receptor_difference_vol_percentage_points"]["pooled_rmse"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
