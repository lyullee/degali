"""Fail-closed audit for the public ELVHYS WP4.2 vent-restriction cases.

This package is a useful cryogenic-H2 enclosure/ventilation component dataset,
not a free-field obstacle-wake benchmark. The audit pins a small paired subset
and keeps source-flow, geometry, calibration and clock limitations explicit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.field_evidence_manifest import (  # noqa: E402
    FieldEvidenceManifest,
    FieldEvidenceManifestArtifact,
)


DOI = "10.18710/JXJP0H"
DATASET_URL = "https://dataverse.no/dataset.xhtml?persistentId=doi%3A10.18710%2FJXJP0H"
DATAFILE_IDS = {
    "00_ELVHYSWP4_2_README.txt": 254380,
    "ELE402HSEMETA.csv": 254386,
    "ELVHYS_WP4.2_TCS_Sensor_Details.pdf": 254561,
    "test029_CONC.csv": 254306,
    "test029_FLMT.csv": 254364,
    "test029_MISC.csv": 254240,
    "test029_PRES.csv": 254217,
    "test031_CONC.csv": 254412,
    "test031_FLMT.csv": 254255,
    "test031_MISC.csv": 254390,
    "test031_PRES.csv": 254246,
}
EXPECTED_SHA256 = {
    "00_ELVHYSWP4_2_README.txt": "c6166236c7087ba4a457990470a69c125df436ac20e67b4b35d94e3a3df3fc42",
    "ELE402HSEMETA.csv": "84d71b9ee69b7814174a15a98243c7681da8cc1e9a8af7b7ebc1307e28e8d488",
    "ELVHYS_WP4.2_TCS_Sensor_Details.pdf": "bedb3115c22e7f080ab10beef6ed0f3269a0997929e56a6455e80e3acd2911b6",
    "test029_CONC.csv": "1457c8a9d13adb698d7c14541c5dbeca2d054139cba1fffdbd262e1cd04345ec",
    "test029_FLMT.csv": "34e71f75d3c2ee295bac9fc6a3e10b1cd0530a71b5b1fe5a1c3d7e9f704f48dc",
    "test029_MISC.csv": "45bd5c4ef399938d12ca4b3302fc8775415e1cafc84cd37efcecc21e5ae0d724",
    "test029_PRES.csv": "3142b9ea12847b8befe8d9d305787f3c25832ea9f49928adbbb6227f640ff00c",
    "test031_CONC.csv": "f13510cf8ce5dc8690114af95d39b1db412d24bd2daa168da7765649780b9be4",
    "test031_FLMT.csv": "cd3e53bcbcbe08ecd52de21baa6bd1aa1d79391c914d08f900a2b4f6c537e4be",
    "test031_MISC.csv": "9c1af0fe044e439c80d12c0287c171114acf6a82c9cb56da4328acd4f8af6f37",
    "test031_PRES.csv": "ebfe71240f150c695187cfa5d902eeb386b1936a5a47ea7b329f17986af06346",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _csv_summary(path: Path, *, concentration: bool = False) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames or []
        times: list[float] = []
        time_tokens: list[str] = []
        minima: dict[str, float] = {}
        maxima: dict[str, float] = {}
        row_count = 0
        for row in reader:
            row_count += 1
            time = float(row["Time"])
            if times and time <= times[-1]:
                raise ValueError(f"{path.name}: Time must increase strictly")
            times.append(time)
            time_tokens.append(row["Time"])
            if concentration:
                for key in columns:
                    if key.startswith("H2"):
                        value = float(row[key])
                        minima[key] = min(minima.get(key, value), value)
                        maxima[key] = max(maxima.get(key, value), value)
    if row_count < 2:
        raise ValueError(f"{path.name}: expected at least two time samples")
    intervals = [b - a for a, b in zip(times, times[1:])]
    h2_channels = [key for key in columns if key.startswith("H2")]
    return {
        "file": path.name,
        "sha256": _sha256(path),
        "rows": row_count,
        "columns": columns,
        "time_s": {"start": times[0], "stop": times[-1]},
        "time_vector_sha256": hashlib.sha256(
            ("\n".join(time_tokens) + "\n").encode("ascii")
        ).hexdigest(),
        "sample_interval_s": {"minimum": min(intervals), "maximum": max(intervals)},
        "h2_sensor_count": len(h2_channels),
        "h2_sensor_range_vol_percent": (
            {key: [minima[key], maxima[key]] for key in h2_channels}
            if concentration else None
        ),
    }


def _sensor_registration(readme: str, observed_channels: list[str]) -> dict[str, Any]:
    lines = readme.splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Sensor\tX (mm)")),
        None,
    )
    if header_index is None:
        raise ValueError("README sensor coordinate table was not found")
    locations: dict[str, dict[str, Any]] = {}
    for line in lines[header_index + 1:]:
        if not line.strip():
            if locations:
                break
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        try:
            xyz = [int(fields[index]) for index in (1, 2, 3)]
        except ValueError:
            continue
        locations[fields[0].strip()] = {
            "coordinates_mm": xyz,
            "sensor_types": [part.strip() for part in fields[4].split(",") if part.strip()],
        }

    def channel_for_location(location: str) -> str | None:
        if location == "Bottom1_Front":
            return "H2Bottom1"
        if location.startswith("Bottom"):
            return "H2" + location
        if location == "Side1_Top":
            return "H2Side1"
        if location.startswith("Side"):
            return "H2" + location
        if location in {"VentBottom", "VentTop", "ChimneyBottom", "ChimneyTop"}:
            return "H2" + location
        return None

    expected_channels = {
        channel_for_location(name)
        for name, value in locations.items()
        if "C" in value["sensor_types"] and channel_for_location(name) is not None
    }
    expected_channels.discard(None)
    observed = {channel for channel in observed_channels if channel.startswith("H2")}
    mappings = []
    for channel in sorted(observed):
        location = channel.removeprefix("H2")
        if channel == "H2Bottom1":
            location = "Bottom1_Front"
        elif channel == "H2Side1":
            location = "Side1_Top"
        record = locations.get(location)
        mappings.append({
            "channel": channel,
            "location_id": location if record is not None else None,
            "coordinates_m": (
                [value / 1000.0 for value in record["coordinates_mm"]]
                if record is not None else None
            ),
            "sensor_types_declared": record["sensor_types"] if record is not None else [],
            "concentration_sensor_declared": bool(
                record is not None and "C" in record["sensor_types"]
            ),
        })
    return {
        "coordinate_table_rows": len(locations),
        "observed_h2_channels": sorted(observed),
        "expected_concentration_channels_from_registry": sorted(expected_channels),
        "observed_channel_without_concentration_registry": sorted(observed - expected_channels),
        "registry_concentration_channel_missing_from_csv": sorted(expected_channels - observed),
        "all_observed_channels_have_coordinates": all(item["coordinates_m"] is not None for item in mappings),
        "sensor_type_registry_matches_csv": observed == expected_channels,
        "channel_mappings": mappings,
    }


def audit(data_root: Path) -> dict[str, Any]:
    files: dict[str, Any] = {}
    for name, file_id in DATAFILE_IDS.items():
        path = data_root / name
        if not path.is_file():
            raise FileNotFoundError(f"missing DataverseNO file {path}")
        digest = _sha256(path)
        expected = EXPECTED_SHA256[name]
        if digest != expected:
            raise ValueError(f"SHA-256 mismatch for {name}: {digest}")
        files[name] = {"datafile_id": file_id, "bytes": path.stat().st_size, "sha256": digest}

    readme = (data_root / "00_ELVHYSWP4_2_README.txt").read_text(encoding="utf-8-sig")
    required_declarations = (
        "Hydrogen sensors are calibrated at 5 points between 0 - 100 vol%H2",
        "PT1ReleaseRig: Pressure in LH2 release rig (barg)",
        "PT2Nozzle: Pressure at LH2 nozzle (barg)",
        "FanFlowMeter: Flow rate of ventilation fan (L/min)",
        "+- 1 vol% H2",
    )
    declaration_checks = {text: text in readme for text in required_declarations}
    if not all(declaration_checks.values()):
        raise ValueError("README is missing a required sensor/source declaration")
    sensor_registration = _sensor_registration(readme, [])

    geometry_path = data_root / "ELVHYS_WP4.2_geometry_registry_v1.json"
    geometry = {
        "schema": "degali.elvhys-tcs-geometry-registry.v1",
        "dataset_doi": DOI,
        "source_artifacts": {
            name: {"datafile_id": DATAFILE_IDS[name], "sha256": files[name]["sha256"]}
            for name in (
                "00_ELVHYSWP4_2_README.txt",
                "ELVHYS_WP4.2_TCS_Sensor_Details.pdf",
            )
        },
        "facility": {
            "nominal_internal_volume_m3": 1.0,
            "wall_material": "stainless steel",
            "wall_thickness_m": 0.003,
            "rear_passive_vent_count": 4,
            "rear_passive_vent_dimensions_reported_m": [0.4, 0.2],
            "vent_area_callout_in_sensor_drawing_m2": 0.02,
            "vent_area_callout_scope": "drawing does not say whether this is per vent or aggregate",
            "blast_panel_dimensions_reported_m": [0.4, 0.8],
            "optional_top_chimney_dimensions_reported_m": [0.075, 0.8],
            "source_horizontal_coordinate_m": [0.02, 0.5, 0.25],
            "source_coordinate_frame": "bottom-right-front origin; x width, y length from front, z upward; see sensor-details PDF",
        },
        "pressure_peaking_test_vent_configurations": {
            "29": {"arrangement": "four VAV plates", "reported_total_open_area_m2": 0.0416},
            "30": {"arrangement": "two blanks and two VAV plates at top", "reported_total_open_area_m2": 0.0208},
            "31": {"arrangement": "three blanks and one top-left VAV plate", "reported_total_open_area_m2": 0.0104},
        },
        "geometry_quantitative_allowed": False,
        "unresolved_items": [
            "The README's vent-panel dimensions and drawing's 0.02 m2 passive-vent callout are not reconciled; the drawing does not state the callout scope.",
            "Test-specific VAV effective areas are reported, but the physical plate aperture profiles and machine-registered opening coordinates are not provided.",
            "This is a confined enclosure with vent restrictions, not a machine-registered open-field obstacle-wake geometry.",
        ],
    }
    geometry_path.write_text(
        json.dumps(geometry, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    geometry_sha256 = _sha256(geometry_path)

    tests: dict[str, Any] = {}
    for test_id, vent_description, expected_dt in (
        (29, "4 VAV plates", 0.05),
        (31, "3 blanks, 1 VAV plate (top left)", 0.002),
    ):
        prefix = f"test{test_id:03d}"
        paths = {
            kind: data_root / f"{prefix}_{kind}.csv"
            for kind in ("CONC", "FLMT", "MISC", "PRES")
        }
        conc = _csv_summary(paths["CONC"], concentration=True)
        sensor_registration = _sensor_registration(readme, conc["columns"])
        flow = _csv_summary(paths["FLMT"])
        weather = _csv_summary(paths["MISC"])
        pressure = _csv_summary(paths["PRES"])
        clock_hashes = {
            summary["time_vector_sha256"]
            for summary in (conc, flow, weather, pressure)
        }
        common_clock_verified = len(clock_hashes) == 1 and len({
            summary["rows"] for summary in (conc, flow, weather, pressure)
        }) == 1
        if not common_clock_verified:
            raise ValueError(f"test {test_id} channel time vectors do not match")
        if abs(conc["sample_interval_s"]["minimum"] - expected_dt) > 1e-9:
            raise ValueError(f"unexpected concentration sample rate for test {test_id}")
        if "PT1ReleaseRig" not in pressure["columns"] or "PT2Nozzle" not in pressure["columns"]:
            raise ValueError(f"test {test_id} is missing pressure-based release channels")
        file_paths = {
            "source_boundary": paths["PRES"],
            "weather": paths["MISC"],
            "obstacle_geometry": geometry_path,
            "receptor_observations": paths["CONC"],
            "common_clock": paths["MISC"],
        }
        artifacts = tuple(
            FieldEvidenceManifestArtifact(
                channel=channel,
                path=str(path.resolve()),
                sha256=_sha256(path),
            )
            for channel, path in file_paths.items()
        )
        manifest = FieldEvidenceManifest(
            manifest_id=f"elvhys-wp4-2-test{test_id:03d}-manifest-v1",
            event_id=f"ELVHYS-WP4.2-test{test_id:03d}",
            audit_root=str(data_root.resolve()),
            channel_artifacts=artifacts,
            dataset_id=f"doi:{DOI}:test{test_id:03d}",
            observed_dataset_path=str(paths["CONC"].resolve()),
            observed_dataset_sha256=_sha256(paths["CONC"]),
            observed_row_count=conc["rows"],
            source_boundary_id=(
                f"test{test_id:03d}:LH2-pressure-history-no-measured-mass-flow"
            ),
            weather_id=f"test{test_id:03d}:time-series-ambient-conditions",
            obstacle_geometry_id=(
                f"test{test_id:03d}:{vent_description}:reported-vent-restriction-"
                "not-machine-registered"
            ),
            receptor_geometry_id=(
                f"test{test_id:03d}:coordinate-label-join-with-sensor-type-registry-mismatch"
            ),
            temporal_operator_id=(
                "H2-sample-delay-corrected-in-published-data-processing"
            ),
            common_clock_id=(
                f"test{test_id:03d}:exact-shared-raw-Time-vector:"
                f"{expected_dt:g}s-sampling"
            ),
            scope="lh2_obstacle_transport",
            sensor_set_id=f"elvhys-test{test_id:03d}-xensor-h2-v1",
            operator_id="degali-local-audit-2026-10-12",
            sensor_registry_artifact=FieldEvidenceManifestArtifact(
                channel="sensor_registry",
                path=str((data_root / "00_ELVHYSWP4_2_README.txt").resolve()),
                sha256=_sha256(data_root / "00_ELVHYSWP4_2_README.txt"),
            ),
            sensor_calibration_status="specification_only",
        )
        manifest.verify_files()
        tests[str(test_id)] = {
            "declared_configuration": vent_description,
            "nominal_release": "LH2, 25.4 mm nozzle, 2 barg",
            "concentration": conc,
            "ventilation_flow": flow,
            "release_pressure": pressure,
            "weather": weather,
            "common_clock_verified": common_clock_verified,
            "geometry_quantitative_allowed": geometry["geometry_quantitative_allowed"],
            "sensor_registration": sensor_registration,
            "field_evidence_manifest": manifest.as_record(),
            "release_mass_flow_measured": any(
                "mass" in column.lower() and "flow" in column.lower()
                for column in flow["columns"]
            ),
        }

    limitations = [
        "The test cases concern a 1 m3 transfer-connection enclosure and restricted vents, not an open-field obstacle wake.",
        "Flow-meter channels are ventilation flow, not measured hydrogen release mass flow; release pressure channels are available.",
        "The README reports nominal vent configurations but the detailed registered obstacle/vent geometry and event-specific calibration certificates are not in this subset.",
        "No DEGALI concentration score is emitted: LH2 source thermodynamics and geometry must first be mapped without fitting concentration observations.",
    ]
    if not sensor_registration["sensor_type_registry_matches_csv"]:
        limitations.append(
            "Concentration CSV sensor labels do not fully agree with the README/PDF sensor-type registry; unmatched IDs are listed explicitly."
        )
    return {
        "schema": "degali.elvhys-tcs-component-evidence-audit.v1",
        "dataset": {"doi": DOI, "url": DATASET_URL, "license": "CC0-1.0"},
        "status": "conditional_cryogenic_h2_enclosure_vent_restriction_component",
        "promotion_allowed": False,
        "files": files,
        "geometry_registry": {
            "path": str(geometry_path.resolve()),
            "sha256": geometry_sha256,
            "quantitative_geometry_allowed": geometry["geometry_quantitative_allowed"],
            "unresolved_items": geometry["unresolved_items"],
        },
        "declaration_checks": declaration_checks,
        "sensor_registration": sensor_registration,
        "tests": tests,
        "limitations": limitations,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = audit(args.data_root.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "promotion_allowed": False, "tests": list(record["tests"])}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
