"""Audit the public H2SAFE helium-surrogate indoor-release package.

The audit verifies the official archive/document digests, proves that the
extracted files equal their ZIP members, joins every trace sensor to the
published coordinate registry, and emits one fail-closed evidence manifest
per event.  Missing concentration units, release-clock synchronization and
event calibration certificates keep quantitative concentration validation
withheld.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
from typing import Any
from zipfile import ZipFile

import numpy as np

from degali.validation.h2safe import (
    CONCENTRATION_UNIT_STATUS,
    DATASET_DOI,
    DATASET_URL,
    LAB1_VENT_BOUNDARIES,
    LAB2_RECIRCULATION_FANS,
    LAB2_SUPPLY_BOUNDARIES,
    PUBLISHED_VERTICAL_AXIS,
    RELEASES,
    SENSOR_CALIBRATION_STATUS,
    SOURCE_CLOCK_ALIGNMENT_STATUS,
    SPECIES,
    H2SafeRun,
    read_h2safe_ascii_stl_bounds,
    read_h2safe_run,
    read_h2safe_step_product_names,
)


SCHEMA = "degali.h2safe-dataset-audit.v1"
MANIFEST_SCHEMA = "degali.h2safe-field-evidence-manifest.v1"
OFFICIAL_ARCHIVE_SHA256 = (
    "f1b8bf747b7787f00152ec416940c80a143ea044fc18ee8a69ac9da204d9158e"
)
OFFICIAL_HVAC_SHA256 = (
    "e99119f5f9aadac6e0bd8269b3fc7989a6026bc0a4c68221ea2a54b1eeee58bb"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _zip_member_sha256(archive: ZipFile, name: str) -> str:
    digest = hashlib.sha256()
    with archive.open(name) as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _geometry_files(extracted_root: Path, lab_id: str) -> tuple[Path, Path]:
    lab = extracted_root / lab_id
    step = next(iter(sorted(lab.glob("*.STEP"))), None)
    stl = next(iter(sorted(lab.glob("*.STL"))), None)
    if step is None or stl is None:
        raise FileNotFoundError(f"H2SAFE {lab_id} STEP/STL geometry is incomplete")
    return step, stl


def _bounds(run: H2SafeRun) -> dict[str, list[float]]:
    points = np.asarray([sensor.position_m for sensor in run.registry], dtype=float)
    return {
        "minimum_m": [float(value) for value in np.min(points, axis=0)],
        "maximum_m": [float(value) for value in np.max(points, axis=0)],
    }


def _hvac_summary(lab_id: str) -> dict[str, Any]:
    if lab_id == "Lab1_Small":
        boundaries = LAB1_VENT_BOUNDARIES
        fans = ()
    elif lab_id == "Lab2_Large":
        boundaries = LAB2_SUPPLY_BOUNDARIES
        fans = LAB2_RECIRCULATION_FANS
    else:  # pragma: no cover - event metadata controls this branch
        raise ValueError(f"unsupported H2SAFE HVAC lab {lab_id}")
    inlet = sum(
        item.volumetric_flow_m3_s for item in boundaries if item.role == "inlet"
    )
    outlet = sum(
        item.volumetric_flow_m3_s for item in boundaries if item.role == "outlet"
    )
    relative_imbalance = None
    if inlet > 0.0 and outlet > 0.0:
        relative_imbalance = abs(inlet - outlet) / max(inlet, outlet)
    return {
        "boundaries": [
            {**asdict(item), "volumetric_flow_m3_s": item.volumetric_flow_m3_s}
            for item in boundaries
        ],
        "reported_inlet_flow_m3_s": inlet,
        "reported_outlet_flow_m3_s": outlet,
        "reported_relative_flow_imbalance": relative_imbalance,
        "recirculation_fans": [asdict(item) for item in fans],
        "mass_balance_status": (
            "conditional_reported_inlet_outlet_imbalance"
            if relative_imbalance is not None and relative_imbalance > 0.01
            else "withheld_external_outlet_not_reported"
            if outlet == 0.0
            else "accepted_reported_flow_balance"
        ),
        "boundary_face_mapping_status": "withheld_not_derived_from_step_topology",
    }


def _event_manifest(
    run: H2SafeRun, *, hvac_path: Path, step_path: Path, stl_path: Path,
) -> dict[str, Any]:
    release = run.release
    lab_token = release.lab_id.lower().replace("_", "-")
    event_token = f"{lab_token}-{release.test_id}"
    missing_requirements = [
        "concentration_unit_declaration",
        "relative_trace_to_release_clock_synchronization",
        "exact_release_duration_bounds_for_clock_text_events",
        "event_sensor_calibration_certificate",
        "helium_to_hydrogen_species_equivalence_model",
        "vent_boundary_face_mapping_and_normals",
    ]
    gate_codes = [
        "helium_surrogate_component_only",
        "concentration_unit_missing",
        "release_clock_alignment_missing",
        "release_duration_bounds_missing",
        "sensor_calibration_certificate_missing",
        "species_equivalence_model_missing",
        "vent_boundary_faces_unresolved",
        "not_lh2_field_validation",
    ]
    if release.lab_id == "Lab2_Large":
        missing_requirements.append("external_hvac_outlet_flow_boundary")
        gate_codes.append("external_hvac_outlet_flow_missing")
    return {
        "schema": MANIFEST_SCHEMA,
        "manifest_id": f"h2safe-{event_token}-manifest-v1",
        "event_id": run.event_id,
        "dataset_id": f"H2SAFE-{DATASET_DOI}",
        "species": SPECIES,
        "source_boundary_id": f"h2safe-{event_token}-documented-release-v1",
        "weather_id": f"h2safe-{lab_token}-forced-hvac-boundary-v1",
        "obstacle_geometry_id": f"h2safe-{lab_token}-step-stl-v1",
        "receptor_geometry_id": f"h2safe-{lab_token}-sensor-registry-v1",
        "sensor_set_id": f"h2safe-{lab_token}-published-sensors-v1",
        "operator_id": "degali-h2safe-audit-v1",
        "temporal_operator_id": "native-one-second-relative-clock-no-interpolation-v1",
        "common_clock_id": f"h2safe-{event_token}-relative-seconds-unaligned-v1",
        "time_synchronization": {
            "status": "withheld",
            "trace_clock": "CSV relative seconds",
            "release_clock": release.documented_clock_window,
            "alignment_status": SOURCE_CLOCK_ALIGNMENT_STATUS,
        },
        "concentration": {
            "native_unit_status": CONCENTRATION_UNIT_STATUS,
            "vol_percent_conversion_allowed": False,
            "native_maximum": run.native_maximum,
        },
        "sensor_calibration_status": SENSOR_CALIBRATION_STATUS,
        "artifacts": {
            "source": {"path": str(hvac_path.resolve()), "sha256": _sha256(hvac_path)},
            "weather_hvac": {"path": str(hvac_path.resolve()), "sha256": _sha256(hvac_path)},
            "obstacle_step": {"path": str(step_path.resolve()), "sha256": _sha256(step_path)},
            "obstacle_stl": {"path": str(stl_path.resolve()), "sha256": _sha256(stl_path)},
            "receptor_observations": {
                "path": str(run.source_path), "sha256": run.source_sha256,
            },
            "sensor_registry": {
                "path": str(run.registry_path), "sha256": run.registry_sha256,
            },
            "common_clock": {
                "path": str(run.source_path), "sha256": run.source_sha256,
            },
        },
        "status": "conditional",
        "missing_requirements": missing_requirements,
        "qualitative_component_use_allowed": True,
        "quantitative_concentration_validation_allowed": False,
        "quantitative_hydrogen_validation_allowed": False,
        "promotion_allowed": False,
        "gate_codes": gate_codes,
    }


def _run_record(
    run: H2SafeRun, *, hvac_path: Path, step_path: Path, stl_path: Path,
) -> dict[str, Any]:
    values = run.concentration_raw
    stl_bounds = read_h2safe_ascii_stl_bounds(stl_path)
    product_names = read_h2safe_step_product_names(step_path)
    ventilation_names = tuple(
        name for name in product_names
        if any(token in name.lower() for token in (
            "inlet", "outlet", "exhaust", "supply", "duct", "vent",
            "fan", "snorkel", "recirculation", "room",
        ))
    )
    if not ventilation_names:
        raise ValueError("published STEP has no identifiable room/HVAC components")
    source_inside = stl_bounds.contains(run.release.release_position_m)
    sensors_inside = all(
        stl_bounds.contains(sensor.position_m) for sensor in run.registry
    )
    if not source_inside or not sensors_inside:
        raise ValueError(
            "published source/sensor coordinates fall outside the STL coordinate bounds"
        )
    return {
        "event_id": run.event_id,
        "lab_id": run.release.lab_id,
        "test_id": run.release.test_id,
        "release": asdict(run.release),
        "sample_count": int(run.time_s.size),
        "time_start_s": float(run.time_s[0]),
        "time_end_s": float(run.time_s[-1]),
        "sample_interval_s": run.sample_interval_s,
        "clock_strictly_increasing": bool(np.all(np.diff(run.time_s) > 0.0)),
        "registry_sensor_count": len(run.registry),
        "trace_sensor_count": len(run.sensor_ids),
        "finite_trace_sensor_count": len(run.finite_sensor_ids),
        "fully_missing_sensor_ids": list(run.fully_missing_sensor_ids),
        "missing_cell_count": int(np.sum(~np.isfinite(values))),
        "native_concentration_maximum": run.native_maximum,
        "native_concentration_unit_status": CONCENTRATION_UNIT_STATUS,
        "sensor_coordinate_bounds": _bounds(run),
        "coordinate_axis_mapping": {
            "published_columns": ["X", "Y", "Z"],
            "published_vertical_axis": PUBLISHED_VERTICAL_AXIS,
            "solver_mapping": {"x": "X", "y": "Z", "z": "Y"},
            "release_position_solver_m": list(run.release.solver_position_m),
            "status": "accepted",
        },
        "stl_coordinate_bounds": {
            "minimum_m": list(stl_bounds.minimum_m),
            "maximum_m": list(stl_bounds.maximum_m),
            "vertex_record_count": stl_bounds.vertex_record_count,
            "sha256": stl_bounds.sha256,
        },
        "coordinate_frame_check": {
            "source_inside_stl_bounds": source_inside,
            "all_registry_sensors_inside_stl_bounds": sensors_inside,
            "tolerance_m": 1.0e-5,
            "status": "accepted",
        },
        "step_semantic_inventory": {
            "product_count": len(product_names),
            "ventilation_or_room_products": list(ventilation_names),
            "boundary_face_mapping_status": "withheld_not_derived_from_step_topology",
        },
        "hvac_boundary_summary": _hvac_summary(run.release.lab_id),
        "validation_scope": run.validation_scope,
        "quantitative_hydrogen_validation_allowed": False,
        "vol_percent_conversion_allowed": False,
        "manifest": _event_manifest(
            run, hvac_path=hvac_path, step_path=step_path, stl_path=stl_path,
        ),
    }


def audit_h2safe_dataset(
    raw_root: str | Path,
    *,
    expected_archive_sha256: str | None = OFFICIAL_ARCHIVE_SHA256,
    expected_hvac_sha256: str | None = OFFICIAL_HVAC_SHA256,
) -> dict[str, Any]:
    root = Path(raw_root)
    archive_path = root / "Dataset.zip"
    hvac_path = root / "HVAC_information.docx"
    extracted_root = root / "extracted"
    if not root.is_dir():
        return {
            "schema": SCHEMA, "status": "not_available", "dataset_doi": DATASET_DOI,
            "root": str(root.resolve()), "events": [], "promotion_allowed": False,
            "gate_codes": ["dataset_not_available", "not_lh2_field_validation"],
        }

    errors: list[str] = []
    archive_sha = _sha256(archive_path) if archive_path.is_file() else None
    hvac_sha = _sha256(hvac_path) if hvac_path.is_file() else None
    if archive_sha is None:
        errors.append("Dataset.zip is missing")
    elif expected_archive_sha256 is not None and archive_sha != expected_archive_sha256:
        errors.append("Dataset.zip SHA-256 differs from the pinned public resource")
    if hvac_sha is None:
        errors.append("HVAC_information.docx is missing")
    elif expected_hvac_sha256 is not None and hvac_sha != expected_hvac_sha256:
        errors.append("HVAC_information.docx SHA-256 differs from the pinned public resource")

    archive_members: list[dict[str, Any]] = []
    if archive_path.is_file() and extracted_root.is_dir():
        with ZipFile(archive_path) as archive:
            for info in archive.infolist():
                name = PurePosixPath(info.filename)
                if info.is_dir():
                    continue
                if name.is_absolute() or ".." in name.parts:
                    errors.append(f"unsafe archive member: {info.filename}")
                    continue
                extracted = extracted_root.joinpath(*name.parts)
                member_sha = _zip_member_sha256(archive, info.filename)
                extracted_sha = _sha256(extracted) if extracted.is_file() else None
                matches = extracted_sha == member_sha
                if not matches:
                    errors.append(f"archive/extracted mismatch: {info.filename}")
                archive_members.append({
                    "path": info.filename,
                    "bytes": info.file_size,
                    "sha256": member_sha,
                    "extracted_path": str(extracted.resolve()),
                    "extracted_matches_archive": matches,
                })

    events: list[dict[str, Any]] = []
    if hvac_path.is_file() and extracted_root.is_dir():
        for lab_id, test_id in RELEASES:
            try:
                run = read_h2safe_run(extracted_root, lab_id=lab_id, test_id=test_id)
                step_path, stl_path = _geometry_files(extracted_root, lab_id)
                events.append(_run_record(
                    run, hvac_path=hvac_path, step_path=step_path, stl_path=stl_path,
                ))
            except (OSError, TypeError, ValueError) as error:
                errors.append(f"{lab_id}/{test_id}: {error}")

    structurally_complete = len(events) == len(RELEASES) and not errors
    status = "conditional_component_evidence_helium_surrogate" if structurally_complete else "withheld"
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "dataset_doi": DATASET_DOI,
        "dataset_url": DATASET_URL,
        "root": str(root.resolve()),
        "license": "public catalog; consult linked license record",
        "source_integrity": {
            "archive": {
                "path": str(archive_path.resolve()), "sha256": archive_sha,
                "expected_sha256": expected_archive_sha256,
                "verified": archive_sha is not None and (
                    expected_archive_sha256 is None or archive_sha == expected_archive_sha256
                ),
            },
            "hvac_document": {
                "path": str(hvac_path.resolve()), "sha256": hvac_sha,
                "expected_sha256": expected_hvac_sha256,
                "verified": hvac_sha is not None and (
                    expected_hvac_sha256 is None or hvac_sha == expected_hvac_sha256
                ),
            },
            "archive_members": archive_members,
        },
        "event_count": len(events),
        "events": events,
        "errors": errors,
        "evidence_boundary": {
            "species": SPECIES,
            "scope": "indoor forced-ventilation transport component evidence",
            "concentration_unit_status": CONCENTRATION_UNIT_STATUS,
            "source_clock_alignment_status": SOURCE_CLOCK_ALIGNMENT_STATUS,
            "sensor_calibration_status": SENSOR_CALIBRATION_STATUS,
            "direct_vol_percent_comparison_allowed": False,
            "quantitative_hydrogen_validation_allowed": False,
            "quantitative_lh2_field_validation_allowed": False,
        },
        "promotion_allowed": False,
        "gate_codes": [
            "public_files_verified" if structurally_complete else "public_files_withheld",
            "helium_surrogate_component_only",
            "concentration_unit_missing",
            "release_clock_alignment_missing",
            "release_duration_bounds_missing",
            "sensor_calibration_certificate_missing",
            "vent_boundary_faces_unresolved",
            "not_lh2_field_validation",
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
    parser.add_argument("raw_root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--accept-unpinned-version", action="store_true")
    args = parser.parse_args(argv)
    record = audit_h2safe_dataset(
        args.raw_root,
        expected_archive_sha256=(None if args.accept_unpinned_version else OFFICIAL_ARCHIVE_SHA256),
        expected_hvac_sha256=(None if args.accept_unpinned_version else OFFICIAL_HVAC_SHA256),
    )
    if args.output is None:
        print(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False))
    else:
        write_audit(record, args.output)
        print(f"wrote {args.output} status={record['status']} events={record['event_count']}")
    return 0 if record["status"] in {
        "conditional_component_evidence_helium_surrogate", "not_available",
    } else 1


if __name__ == "__main__":
    raise SystemExit(main())
