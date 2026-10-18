"""Build a fail-closed FieldEvidenceManifest intake record for FFI obstacle cases.

The public USN/FFI files contain a wide source/sensor CSV and a separate
obstacle/sensor registry.  They are useful evidence, but they do not by
themselves establish a weather channel or a common source/weather/sensor
clock.  This tool fingerprints the available files and emits an intake record
without fabricating a five-channel ``FieldEvidenceManifest`` or allowing
promotion.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from typing import Any


SCHEMA = "degali.ffi-field-evidence-intake.v1"
_SENSOR_HEADER_RE = re.compile(r"^sensor \d+ h2 concentration")
_MASS_FLOW_HEADER_RE = re.compile(r"^mass flow meter \d+ \[g/s\]$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _csv_summary(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        rows = sum(1 for _ in reader)
    columns = [str(value).strip() for value in header]
    sensor_columns = [value for value in columns if _SENSOR_HEADER_RE.match(value)]
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "row_count": rows,
        "columns": columns,
        "sensor_column_count": len(sensor_columns),
        "source_time_column_available": "flow time [s]" in columns,
        "sensor_time_column_available": "h2 sensor time" in columns,
        "mass_flow_column_available": any(
            _MASS_FLOW_HEADER_RE.fullmatch(value) for value in columns
        ),
    }


def _file_summary(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def _find_csv(data_root: Path, filename: str) -> Path:
    matches = sorted(data_root.rglob(filename))
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one extracted CSV named {filename!r}; found {len(matches)}"
        )
    return matches[0]


def build_intake(
    data_root: str | Path,
    registry_path: str | Path,
    *,
    article_case: str,
    output: str | Path,
) -> dict[str, Any]:
    root = Path(data_root).resolve()
    registry_file = Path(registry_path).resolve()
    output_file = Path(output).resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    if not registry_file.is_file():
        raise FileNotFoundError(registry_file)
    registry = json.loads(registry_file.read_text(encoding="utf-8"))
    cases = [
        case for case in registry.get("case_crosswalk", [])
        if case.get("article_case") == article_case
    ]
    if len(cases) != 1:
        raise ValueError(f"article case {article_case!r} is not unique in the registry")
    case = cases[0]
    csv_paths = [
        _find_csv(root, filename) for filename in case["dataset_files"]
    ]
    csv_summaries = [_csv_summary(path) for path in csv_paths]
    obstacle = registry.get("obstacle", {})
    time_operator = registry.get("time_operator", {})
    promotion = registry.get("promotion", {})
    sensor_instrumentation = registry.get("sensor_instrumentation", {})
    sensor_count = sensor_instrumentation.get(
        "sensor_count",
        time_operator.get("sensor_count", len(registry.get("sensors", []))),
    )

    channels: dict[str, dict[str, Any]] = {
        "source_boundary": {
            "status": "conditional",
            "artifacts": csv_summaries,
            "source_boundary_id": None,
            "available_evidence": [
                "FFI mass-flow history is present in each selected raw CSV",
                "registry article LES source-rate crosswalk is present",
            ],
            "missing_requirements": [
                "event-specific source-boundary ID in the raw source record",
                "complete release boundary history with end-state reconciliation",
            ],
        },
        "weather": {
            "status": "withheld",
            "artifacts": [],
            "missing_requirements": [
                "event-linked wind speed, direction, measurement height and stability history",
            ],
        },
        "obstacle_geometry": {
            "status": "withheld",
            "artifacts": [_file_summary(registry_file)],
            "missing_requirements": [
                "machine-registered beam assembly dimensions and coordinate anchor",
            ],
            "registry_geometry_allowed": bool(obstacle.get("quantitative_geometry_allowed", False)),
        },
        "receptor_observations": {
            "status": "conditional",
            "artifacts": csv_summaries + [_file_summary(registry_file)],
            "sensor_set_id": f"ffi-usn-xen5320-{sensor_count}",
            "available_evidence": [
                "wide H2 sensor concentration histories are present",
                "registry supplies sensor coordinates and heights",
            ],
            "missing_requirements": [
                "explicit joined receptor table with sensor ID, x/y/z and timestamp",
                "event-specific calibration certificate",
            ],
        },
        "common_clock": {
            "status": "withheld",
            "artifacts": [],
            "missing_requirements": [
                "hash-pinned source/weather/sensor synchronization transform",
            ],
            "dataset_source_and_sensor_clocks_are_distinct": bool(
                time_operator.get("dataset_source_and_sensor_clocks_are_distinct", True)
            ),
        },
    }
    missing_channels = [
        name for name, channel in channels.items()
        if channel["status"] == "withheld"
    ]
    gate_codes = list(dict.fromkeys([
        "ffi_field_evidence_manifest_withheld",
        *promotion.get("gate_codes", []),
        "weather_channel_missing",
        "common_clock_unverified",
    ]))
    reasons = [
        "the intake is not a five-channel FieldEvidenceManifest and promotion_allowed is false",
        "raw FFI source and receptor histories require explicit event/geometry/clock joining",
    ]
    if missing_channels:
        reasons.append("required channels withheld: " + ", ".join(missing_channels))
    record: dict[str, Any] = {
        "schema": SCHEMA,
        "event_id": f"USN-FFI-{article_case}",
        "article_case": article_case,
        "status": "withheld",
        "promotion_allowed": False,
        "operator_id": "legacy-unspecified",
        "manifest_candidate": {
            "source_boundary_id": None,
            "weather_id": None,
            "obstacle_geometry_id": (
                f"ffi-usn-{article_case}-registry-conditional"
            ),
            "receptor_geometry_id": f"ffi-usn-{article_case}-sensor-registry",
            "common_clock_id": None,
            "sensor_set_id": channels["receptor_observations"]["sensor_set_id"],
        },
        "selected_trials": [
            {
                "trial_id": trial_id,
                "dataset_file": summary["path"],
                "dataset_sha256": summary["sha256"],
                "row_count": summary["row_count"],
            }
            for trial_id, summary in zip(case["article_trials"], csv_summaries)
        ],
        "channels": channels,
        "gate_codes": gate_codes,
        "reasons": reasons,
        "registry": _file_summary(registry_file),
    }
    if output_file.exists():
        raise FileExistsError(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("reference/hydrogen/open_channel_usn/extracted"))
    parser.add_argument("--registry", type=Path, default=Path("reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"))
    parser.add_argument("--article-case", choices=("O1", "O2"), default="O1")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = build_intake(
        args.data_root, args.registry,
        article_case=args.article_case, output=args.output,
    )
    print(json.dumps({
        "status": record["status"],
        "article_case": record["article_case"],
        "gate_codes": record["gate_codes"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
