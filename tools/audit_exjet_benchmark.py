"""Fail-closed readiness audit for the public ExJet obstacle benchmark.

The paper is useful evidence, but paper plots are not silently treated as raw
observations.  This tool accepts machine-readable geometry and observation
files only when they are supplied explicitly and reports a withheld status
otherwise.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REQUIRED_OBSERVATION_COLUMNS = {
    "station_x_m",
    "radial_coordinate_m",
    "h2_molar_fraction",
    "quantity_kind",
    "coordinate_uncertainty_m",
    "concentration_uncertainty_mole_fraction",
}
REQUIRED_QUANTITIES = {
    "centreline_hydrogen_molar_fraction",
    "radial_hydrogen_molar_fraction",
}
EXPECTED_STATIONS_M = (1.4, 2.5, 4.4)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_observations(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or ())
        missing = sorted(REQUIRED_OBSERVATION_COLUMNS - columns)
        if missing:
            return [], missing
        return list(reader), []


def _read_json_object(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def audit(
    benchmark_path: Path,
    observations: Path | None,
    geometry: Path | None,
    common_clock: Path | None,
    source_term: Path | None,
    source_pdf: Path | None = None,
) -> dict[str, Any]:
    benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
    requirements = benchmark["degali_use"]["required_before_scoring"]
    spatial_missing: list[str] = []
    transient_missing: list[str] = []
    if observations is None:
        spatial_missing.append("digitised_or_author_supplied_observation_coordinates")
    elif not observations.is_file():
        spatial_missing.append("observation_file_not_found")
    if geometry is None:
        spatial_missing.append("machine_registered_cylinder_array_geometry")
    elif not geometry.is_file():
        spatial_missing.append("geometry_file_not_found")
    if source_term is None or not source_term.is_file():
        spatial_missing.append("source_term_provenance_and_uncertainty")
    elif not (_read_json_object(source_term) or {}).get("provenance"):
        spatial_missing.append("source_term_provenance_and_uncertainty")

    observation_count = 0
    quantities: list[str] = []
    stations: list[float] = []
    observation_sha256 = None
    observation_expected_sha256 = benchmark.get("observations", {}).get(
        "paper_plot_digitisation_sha256"
    )
    observation_integrity_verified = False
    geometry_sha256 = None
    source_pdf_sha256 = None
    source_pdf_expected_sha256 = benchmark.get("source", {}).get("local_pdf_sha256")
    if source_pdf is not None:
        if not source_pdf.is_file():
            missing.append("source_pdf_not_found")
        else:
            source_pdf_sha256 = _sha256(source_pdf)
            if (
                source_pdf_expected_sha256
                and source_pdf_sha256.lower() != str(source_pdf_expected_sha256).lower()
            ):
                missing.append("source_pdf_sha256_mismatch")
    common_clock_sha256 = _sha256(common_clock) if common_clock and common_clock.is_file() else None
    source_term_sha256 = _sha256(source_term) if source_term and source_term.is_file() else None
    observations_have_time = False
    if observations is not None and observations.is_file():
        rows, missing_columns = _read_observations(observations)
        spatial_missing.extend(f"observation_column:{column}" for column in missing_columns)
        observation_count = len(rows)
        quantities = sorted({row.get("quantity_kind", "") for row in rows})
        observations_have_time = bool(rows) and all(
            row.get("observation_time_s", "").strip() for row in rows
        )
        try:
            stations = sorted({float(row["station_x_m"]) for row in rows})
            concentrations = [float(row["h2_molar_fraction"]) for row in rows]
            coordinate_uncertainties = [
                float(row["coordinate_uncertainty_m"]) for row in rows
            ]
            concentration_uncertainties = [
                float(row["concentration_uncertainty_mole_fraction"])
                for row in rows
            ]
        except (KeyError, TypeError, ValueError):
            spatial_missing.append("observation_numeric_values")
            stations = []
            concentrations = []
            coordinate_uncertainties = []
            concentration_uncertainties = []
        if concentrations and any(value < 0.0 or value > 1.0 for value in concentrations):
            spatial_missing.append("observation_h2_molar_fraction_out_of_range")
        if any(
            not math.isfinite(value) or value <= 0.0
            for value in coordinate_uncertainties + concentration_uncertainties
        ):
            spatial_missing.append("observation_uncertainties_invalid")
        for expected in EXPECTED_STATIONS_M:
            if not any(abs(actual - expected) <= 1.0e-6 for actual in stations):
                spatial_missing.append(f"observation_station:{expected:g}m")
        if not rows:
            spatial_missing.append("observation_rows")
        observation_sha256 = _sha256(observations)
        observation_integrity_verified = bool(
            observation_expected_sha256
            and observation_sha256.lower() == str(observation_expected_sha256).lower()
        )
        if observation_expected_sha256 and not observation_integrity_verified:
            spatial_missing.append("observation_file_sha256_mismatch")
    if geometry is not None and geometry.is_file():
        geometry_record = _read_json_object(geometry)
        if (
            not geometry_record
            or geometry_record.get("shape") != "cylinder_array"
            or geometry_record.get("registration_status")
            != "machine_registered_exact"
        ):
            spatial_missing.append("machine_registered_cylinder_array_geometry")
        geometry_sha256 = _sha256(geometry)

    absent_quantities = sorted(REQUIRED_QUANTITIES - set(quantities))
    spatial_missing.extend(
        f"observation_quantity:{quantity}" for quantity in absent_quantities
    )
    clock_record = (
        _read_json_object(common_clock)
        if common_clock is not None and common_clock.is_file() else None
    )
    clock_valid = bool(
        clock_record
        and clock_record.get("release_time_s") is not None
        and isinstance(clock_record.get("time_basis_id"), str)
        and clock_record["time_basis_id"].strip()
    )
    if not observations_have_time:
        transient_missing.append("time_resolved_observation_history")
    if not clock_valid:
        transient_missing.append("common_release_time_definition")
    missing = sorted(set(spatial_missing + transient_missing))
    spatial_ready = not spatial_missing
    transient_ready = spatial_ready and not transient_missing
    if transient_ready:
        status = "accepted_for_spatial_and_transient_quantitative_scoring"
    elif spatial_ready:
        status = "spatial_profile_ready_transient_withheld"
    else:
        status = "withheld_spatial_profile_inputs"
    return {
        "schema": "degali.exjet-benchmark-audit.v1",
        "benchmark_id": benchmark["benchmark_id"],
        "status": status,
        "quantitative_promotion_allowed": transient_ready,
        "spatial_profile_comparison_allowed": spatial_ready,
        "transient_validation_allowed": transient_ready,
        "benchmark_required_fields": requirements,
        "benchmark_spatial_required_fields": benchmark["degali_use"].get(
            "required_before_spatial_profile_scoring", []
        ),
        "benchmark_transient_required_fields": benchmark["degali_use"].get(
            "required_before_transient_validation", []
        ),
        "missing_requirements": sorted(set(missing)),
        "missing_spatial_requirements": sorted(set(spatial_missing)),
        "missing_transient_requirements": sorted(set(transient_missing)),
        "observation_count": observation_count,
        "observation_quantities": quantities,
        "observation_stations_m": stations,
        "observations_have_time_resolved_history": observations_have_time,
        "observation_sha256": observation_sha256,
        "observation_expected_sha256": observation_expected_sha256,
        "observation_integrity_verified": observation_integrity_verified,
        "geometry_sha256": geometry_sha256,
        "source_pdf_sha256": source_pdf_sha256,
        "source_pdf_expected_sha256": source_pdf_expected_sha256,
        "source_pdf_integrity_verified": bool(
            source_pdf_sha256
            and source_pdf_expected_sha256
            and source_pdf_sha256.lower() == str(source_pdf_expected_sha256).lower()
        ),
        "common_clock_sha256": common_clock_sha256,
        "source_term_sha256": source_term_sha256,
        "source_term_provenance_required": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=Path("reference/obstacle_wake/exjet_obstructed_free_jet_benchmark.json"))
    parser.add_argument("--observations", type=Path)
    parser.add_argument("--geometry", type=Path)
    parser.add_argument("--common-clock", type=Path)
    parser.add_argument("--source-term", type=Path)
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = audit(
        args.benchmark, args.observations, args.geometry, args.common_clock,
        args.source_term, args.source_pdf,
    )
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
