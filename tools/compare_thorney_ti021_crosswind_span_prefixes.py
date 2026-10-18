"""Compare valid TI021 nested runs over a common prefix and sensor geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

try:
    from tools.score_thorney_ti021_nested_145_sensors import (
        _interpolate_rows,
        _overlap_geometry_diagnostic,
        _trace_metrics,
    )
except ModuleNotFoundError:
    from score_thorney_ti021_nested_145_sensors import (
        _interpolate_rows,
        _overlap_geometry_diagnostic,
        _trace_metrics,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_pair(inner_dir: Path, outer_dir: Path, duration_s: float) -> dict[str, Any]:
    inner_record_path = inner_dir / "ti021-inner-transient.json"
    outer_record_path = outer_dir / "ti021-outer-transient.json"
    inner_trace_path = inner_dir / "ti021-inner-receptor-traces.npz"
    outer_trace_path = outer_dir / "ti021-outer-receptor-traces.npz"
    paths = (inner_record_path, outer_record_path, inner_trace_path, outer_trace_path)
    if not all(path.is_file() for path in paths):
        raise FileNotFoundError(f"incomplete inner/outer run pair: {inner_dir}, {outer_dir}")
    inner_record = json.loads(inner_record_path.read_text(encoding="utf-8"))
    outer_record = json.loads(outer_record_path.read_text(encoding="utf-8"))
    if inner_record.get("status") not in {
        "accepted_completed_inner_transient", "accepted_inner_transient_qualification",
    }:
        raise ValueError(f"inner run not accepted: {inner_record.get('status')}")
    if outer_record.get("status") != "accepted_completed_outer_transient":
        raise ValueError(f"outer run not accepted: {outer_record.get('status')}")
    if inner_record["configuration"].get("observation_values_used_to_configure_run") != []:
        raise ValueError("inner run is not observation-blind")
    if outer_record["configuration"].get("observation_values_used_to_configure_run") != []:
        raise ValueError("outer run is not observation-blind")

    with np.load(inner_trace_path, allow_pickle=False) as archive:
        inner = {key: np.asarray(archive[key]).copy() for key in archive.files}
    with np.load(outer_trace_path, allow_pickle=False) as archive:
        outer = {key: np.asarray(archive[key]).copy() for key in archive.files}
    time_s = inner["time_s"]
    if duration_s <= time_s[0] or duration_s > time_s[-1] + 1.0e-9:
        raise ValueError(f"requested duration {duration_s} s is outside inner trace coverage")
    mask = time_s <= duration_s + 1.0e-9
    time_s = time_s[mask]
    inner_metrics = _trace_metrics(time_s, inner["concentration_vol_percent"][:, mask])
    outer_on_inner_clock = _interpolate_rows(
        outer["time_s"], outer["concentration_vol_percent"], time_s,
    )
    outer_metrics = _trace_metrics(time_s, outer_on_inner_clock)
    inner_index = {str(label): index for index, label in enumerate(inner["labels"])}
    outer_index = {str(label): index for index, label in enumerate(outer["labels"])}
    overlap = sorted(set(inner_index) & set(outer_index))
    geometry = _overlap_geometry_diagnostic(inner_record, outer_record, np.asarray(overlap))
    inner_receptors = {
        item["label"]: item for item in inner_record["configuration"]["receptors"]
    }
    outer_receptors = {
        item["label"]: item for item in outer_record["configuration"]["receptors"]
    }
    interior = set(geometry["interior_overlap_sensor_labels"])
    return {
        "inner_record": inner_record,
        "outer_record": outer_record,
        "inner_record_path": inner_record_path,
        "outer_record_path": outer_record_path,
        "inner_trace_path": inner_trace_path,
        "outer_trace_path": outer_trace_path,
        "inner_index": inner_index,
        "outer_index": outer_index,
        "inner_metrics": inner_metrics,
        "outer_metrics": outer_metrics,
        "inner_receptors": inner_receptors,
        "outer_receptors": outer_receptors,
        "interior_labels": interior,
        "geometry": geometry,
    }


def _relative_summary(
    pair: dict[str, Any], labels: list[str], *, peak_floor: float, dose_floor: float,
) -> dict[str, Any]:
    im, om = pair["inner_metrics"], pair["outer_metrics"]
    ii, oi = pair["inner_index"], pair["outer_index"]
    peak = np.asarray([
        abs(float(im["peak"][ii[label]]) - float(om["peak"][oi[label]]))
        / max(abs(float(im["peak"][ii[label]])), peak_floor)
        for label in labels
    ])
    dose = np.asarray([
        abs(float(im["dose"][ii[label]]) - float(om["dose"][oi[label]]))
        / max(abs(float(im["dose"][ii[label]])), dose_floor)
        for label in labels
    ])

    def summarize(values: np.ndarray) -> dict[str, float | None]:
        if not values.size:
            return {"maximum_percent": None, "p95_percent": None, "median_percent": None}
        return {
            "maximum_percent": 100.0 * float(np.max(values)),
            "p95_percent": 100.0 * float(np.percentile(values, 95)),
            "median_percent": 100.0 * float(np.median(values)),
        }

    return {
        "sensor_count": len(labels),
        "labels": labels,
        "peak_change": summarize(peak),
        "dose_change": summarize(dose),
    }


def compare_prefixes(
    *, baseline_inner_dir: Path, baseline_outer_dir: Path,
    candidate_inner_dir: Path, candidate_outer_dir: Path,
    duration_s: float = 120.0, peak_floor: float = 0.1, dose_floor: float = 0.1,
) -> dict[str, Any]:
    """Compare coupling metrics on exactly common interior sensor labels."""
    baseline = _load_pair(baseline_inner_dir, baseline_outer_dir, duration_s)
    candidate = _load_pair(candidate_inner_dir, candidate_outer_dir, duration_s)
    if not baseline["geometry"]["all_shared_sensor_coordinates_match"]:
        raise ValueError("baseline inner/outer sampler-coordinate gate failed")
    if not candidate["geometry"]["all_shared_sensor_coordinates_match"]:
        raise ValueError("candidate inner/outer sampler-coordinate gate failed")
    labels = sorted(baseline["interior_labels"] & candidate["interior_labels"])
    if not labels:
        raise ValueError("run pairs have no common interior overlap sensors")
    for label in labels:
        for run in (baseline, candidate):
            ix, ox = run["inner_receptors"][label], run["outer_receptors"][label]
            if not np.allclose(ix["site_local_xyz_m"], ox["site_local_xyz_m"], rtol=0.0, atol=1.0e-8):
                raise ValueError(f"inner/outer physical receptor coordinates differ for {label}")
        a = baseline["inner_receptors"][label]["site_local_xyz_m"]
        b = candidate["inner_receptors"][label]["site_local_xyz_m"]
        if not np.allclose(a, b, rtol=0.0, atol=1.0e-8):
            raise ValueError(f"baseline/candidate physical receptor coordinates differ for {label}")
    active = [
        label for label in labels
        if max(
            baseline["inner_metrics"]["peak"][baseline["inner_index"][label]],
            baseline["outer_metrics"]["peak"][baseline["outer_index"][label]],
            candidate["inner_metrics"]["peak"][candidate["inner_index"][label]],
            candidate["outer_metrics"]["peak"][candidate["outer_index"][label]],
        ) >= 0.1
    ]

    def run_record(pair: dict[str, Any]) -> dict[str, Any]:
        return {
            "inner_record": {"path": str(pair["inner_record_path"]), "sha256": _sha256(pair["inner_record_path"])},
            "outer_record": {"path": str(pair["outer_record_path"]), "sha256": _sha256(pair["outer_record_path"])},
            "inner_trace": {"path": str(pair["inner_trace_path"]), "sha256": _sha256(pair["inner_trace_path"])},
            "outer_trace": {"path": str(pair["outer_trace_path"]), "sha256": _sha256(pair["outer_trace_path"])},
            "inner_domain_m": pair["inner_record"]["configuration"]["domain_m"],
            "outer_domain_m": pair["outer_record"]["configuration"]["domain_m"],
            "coordinate_gate": pair["geometry"]["all_shared_sensor_coordinates_match"],
            "interior_sensor_count": len(pair["interior_labels"]),
            "interior_location_count": pair["geometry"]["interior_overlap_horizontal_location_count"],
        }

    return {
        "schema": "degali.thorney-ti021-crosswind-span-prefix-comparison.v1",
        "status": "diagnostic_only_no_full_event_or_experimental_validation_claim",
        "duration_s": float(duration_s),
        "concentration_unit": "vol%",
        "dose_unit": "vol% s",
        "relative_error_denominator": {
            "peak_floor_vol_percent": float(peak_floor),
            "dose_floor_vol_percent_s": float(dose_floor),
            "formula": "abs(outer-inner)/max(abs(inner), floor)",
        },
        "sensor_selection": "intersection of interior overlap labels in both runs; active if any pair peak >= 0.1 vol%",
        "common_interior_sensor_count": len(labels),
        "active_sensor_count": len(active),
        "baseline": run_record(baseline),
        "candidate": run_record(candidate),
        "baseline_all_common": _relative_summary(baseline, labels, peak_floor=peak_floor, dose_floor=dose_floor),
        "candidate_all_common": _relative_summary(candidate, labels, peak_floor=peak_floor, dose_floor=dose_floor),
        "baseline_active": _relative_summary(baseline, active, peak_floor=peak_floor, dose_floor=dose_floor),
        "candidate_active": _relative_summary(candidate, active, peak_floor=peak_floor, dose_floor=dose_floor),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-inner-directory", type=Path, required=True)
    parser.add_argument("--baseline-outer-directory", type=Path, required=True)
    parser.add_argument("--candidate-inner-directory", type=Path, required=True)
    parser.add_argument("--candidate-outer-directory", type=Path, required=True)
    parser.add_argument("--duration-s", type=float, default=120.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare_prefixes(
        baseline_inner_dir=args.baseline_inner_directory,
        baseline_outer_dir=args.baseline_outer_directory,
        candidate_inner_dir=args.candidate_inner_directory,
        candidate_outer_dir=args.candidate_outer_directory,
        duration_s=args.duration_s,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "common_interior_sensor_count": result["common_interior_sensor_count"],
        "active_sensor_count": result["active_sensor_count"],
        "baseline_active": result["baseline_active"],
        "candidate_active": result["candidate_active"],
    }, indent=2))


if __name__ == "__main__":
    main()
