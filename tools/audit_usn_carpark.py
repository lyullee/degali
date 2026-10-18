"""Audit the public USN ventilated-car-park H2 dynamic data.

This tool reproduces the publication's 4 vol.% first-arrival/duration
operator for blowdown Experiments 19, 20, 21 and 23 and evaluates a no-fit
well-mixed mass-balance baseline.  It is deliberately component-level:
neither the confined GH2 experiment nor the 0-D comparator is promoted as an
outdoor LH2 validation score for DEGALI.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from degali.validation.usn_carpark import (
    DATASET_DOI,
    compare_well_mixed_baseline,
    integrated_positive_source_mass_g,
    read_usn_carpark_mat,
    threshold_metrics,
)


SCHEMA = "degali.usn-carpark-dynamic-audit.v1"
EXPECTED_MD5 = {
    19: "372c4ddee15a6fdfeac3d29c99659060",
    20: "5cc966bb36c57a3a210c0c5d404b43e6",
    21: "dbee45f5c855e149db1e963918071375",
    23: "d66804a544ede24b013479b9843df657",
}
PUBLISHED_TABLE3 = {
    19: {"first_arrival_s": 16.0, "flammable_duration_s": 285.0},
    20: {"first_arrival_s": 18.0, "flammable_duration_s": 336.0},
    21: {"first_arrival_s": 32.0, "flammable_duration_s": 195.0},
    23: {"first_arrival_s": 30.0, "flammable_duration_s": 83.0},
}
PUBLISHED_TABLE2_CALCULATED_PEAK_PERCENT = {
    19: 12.6,
    20: 15.2,
    21: 8.8,
    23: 7.5,
}


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _supporting_file(root: Path, name: str, role: str) -> dict[str, object] | None:
    path = root / name
    if not path.is_file():
        return None
    return {
        "role": role,
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": _digest(path, "sha256"),
    }


def _finite_mean(values: list[float]) -> float:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("USN audit summary contains a non-finite value")
    return float(np.mean(array))


def audit_usn_carpark(
    root: str | Path, *, experiments: tuple[int, ...] = (19, 20, 21, 23),
) -> dict[str, Any]:
    directory = Path(root)
    if not directory.is_dir():
        return {
            "schema": SCHEMA,
            "status": "not_available",
            "dataset_doi": DATASET_DOI,
            "root": str(directory.resolve()),
            "events": [],
            "promotion_allowed": False,
            "gate_codes": ["dataset_not_available", "not_lh2_field_validation"],
        }
    if not experiments or len(set(experiments)) != len(experiments):
        raise ValueError("USN experiments must be a non-empty unique tuple")

    events: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    for number in experiments:
        path = directory / f"Exp{number:02d}.mat"
        try:
            if not path.is_file():
                raise FileNotFoundError(path)
            md5 = _digest(path, "md5")
            if number in EXPECTED_MD5 and md5 != EXPECTED_MD5[number]:
                raise ValueError(
                    f"published MD5 mismatch: expected {EXPECTED_MD5[number]}, got {md5}"
                )
            run = read_usn_carpark_mat(path)
            metric = threshold_metrics(run)
            comparison = compare_well_mixed_baseline(run)
            paper = PUBLISHED_TABLE3.get(number)
            reproduction = None
            if paper is not None:
                reproduction = {
                    "published_first_arrival_s": paper["first_arrival_s"],
                    "reproduced_first_arrival_s": metric.first_arrival_after_source_s,
                    "absolute_first_arrival_difference_s": abs(
                        float(metric.first_arrival_after_source_s)
                        - paper["first_arrival_s"]
                    ),
                    "published_flammable_duration_s": paper["flammable_duration_s"],
                    "reproduced_flammable_duration_s": metric.exceedance_duration_s,
                    "absolute_duration_difference_s": abs(
                        metric.exceedance_duration_s
                        - paper["flammable_duration_s"]
                    ),
                }
            events.append({
                "experiment": number,
                "file": path.name,
                "bytes": path.stat().st_size,
                "md5": md5,
                "published_md5_verified": number not in EXPECTED_MD5 or md5 == EXPECTED_MD5[number],
                "sha256": run.raw_sha256,
                "source_t0_s": run.source_t0_s,
                "source_sample_count": int(run.source_time_s.size),
                "sensor_count": len(run.sensors),
                "canbus_sensor_count": len(run.canbus_sensors),
                "wifi_sensor_count": len(run.wifi_sensors),
                "room_cloud_sensor_count": len(run.room_cloud_sensors),
                "sensor_without_numeric_position_count": sum(
                    sensor.position_raw_mm is None for sensor in run.sensors
                ),
                "ventilation_flow_mean_m3_h": float(np.mean(run.ventilation_flow_m3_h)),
                "positive_source_mass_0_900_s_g": integrated_positive_source_mass_g(run),
                "threshold_operator": asdict(metric),
                "published_table3_reproduction": reproduction,
                "well_mixed_baseline": {
                    **asdict(comparison),
                    "published_table2_calculated_peak_vol_percent": (
                        PUBLISHED_TABLE2_CALCULATED_PEAK_PERCENT.get(number)
                    ),
                    "model_class": "no-fit 0-D ideal well-mixed ventilation mass balance",
                    "is_degali_free_field_solver": False,
                },
            })
        except (OSError, RuntimeError, TypeError, ValueError) as error:
            errors.append({"experiment": str(number), "error": str(error)})

    arrival_differences = [
        float(event["published_table3_reproduction"]["absolute_first_arrival_difference_s"])
        for event in events if event["published_table3_reproduction"] is not None
    ]
    duration_differences = [
        float(event["published_table3_reproduction"]["absolute_duration_difference_s"])
        for event in events if event["published_table3_reproduction"] is not None
    ]
    peak_ratios = [
        float(event["well_mixed_baseline"]["peak_ratio_model_to_observed"])
        for event in events
    ]
    reproduced = (
        len(events) == len(experiments)
        and len(arrival_differences) == len(experiments)
        and max(arrival_differences, default=math.inf) <= 0.6
        and max(duration_differences, default=math.inf) <= 6.0
    )
    status = "component_evidence_ready" if reproduced and not errors else "withheld"
    support = [
        _supporting_file(directory, "File handling guide.pdf", "public_file_handling_guide"),
        _supporting_file(directory, "USNReadmetemplate.txt", "public_dataset_readme"),
        _supporting_file(directory, "Lach_Gaathaug_2021_ICHS64.pdf", "public_experiment_paper"),
    ]
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "dataset_doi": DATASET_DOI,
        "root": str(directory.resolve()),
        "license": "CC BY 4.0",
        "experiments_requested": list(experiments),
        "events": events,
        "errors": errors,
        "supporting_files": [item for item in support if item is not None],
        "reproduction": {
            "paper_table3_operator_reproduced": reproduced,
            "maximum_first_arrival_difference_s": max(arrival_differences, default=None),
            "maximum_flammable_duration_difference_s": max(duration_differences, default=None),
        },
        "comparison_summary": {
            "event_count": len(events),
            "mean_peak_ratio_well_mixed_to_observed": (
                _finite_mean(peak_ratios) if peak_ratios else None
            ),
            "interpretation": (
                "The no-fit well-mixed baseline overpredicts the local room-array peak and "
                "flammable duration; the spatial/dynamic mismatch is evidence for retaining "
                "resolved sensor operators rather than a single enclosure concentration."
            ),
        },
        "scope": "gh2_confined_dynamic_transport_component_only",
        "promotion_allowed": False,
        "quantitative_lh2_field_validation_allowed": False,
        "gate_codes": [
            "component_evidence_only",
            "not_lh2_field_validation",
            "confined_ventilation_not_supported_by_degali_free_field_core",
            "event_calibration_certificate_missing",
        ],
    }


def write_audit(record: dict[str, Any], output: str | Path) -> None:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--experiment", type=int, action="append")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    experiments = tuple(args.experiment or (19, 20, 21, 23))
    record = audit_usn_carpark(args.root, experiments=experiments)
    if args.output is None:
        print(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False))
    else:
        write_audit(record, args.output)
        print(
            f"wrote {args.output} status={record['status']} "
            f"events={len(record['events'])}"
        )
    return 0 if record["status"] in {"component_evidence_ready", "not_available"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
