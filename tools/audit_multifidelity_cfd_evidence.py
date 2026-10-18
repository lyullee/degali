"""Audit the public MultHyFuel D2.2 CFD/experiment comparison boundary.

The report contains a useful obstructed high-pressure GH2 benchmark and a
multi-code comparison, but not a machine-readable receptor table.  This tool
hash-pins the report, verifies the relevant textual anchors and records why
the published plots can provide contextual CFD evidence but cannot support a
new DEGALI-versus-CFD ranking without matched numerical exports.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any


SCHEMA = "degali.multifidelity-cfd-evidence-audit.v1"
DOCUMENT_URL = (
    "https://multhyfuel.eu/images/event-documents/deliverables/"
    "MHYF_WP2_D22_Assessment_of_dispersion_for_high-pressure_H2_20221111.pdf"
)
EXPECTED_SHA256 = "30c5131e30c9158daaf14559fdcb7e35b736a338e51ed8b6f4df65e7b2a7303f"
EXPECTED_PAGE_COUNT = 105


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reader(path: Path):
    try:
        from pypdf import PdfReader
    except ImportError as error:  # pragma: no cover - environment boundary
        raise RuntimeError(
            "MultHyFuel audit requires pypdf; install DEGALI research dependencies"
        ) from error
    return PdfReader(path)


def _page_text(reader, page_number: int) -> str:
    return " ".join((reader.pages[page_number - 1].extract_text() or "").split())


def audit_multifidelity_cfd_evidence(pdf: str | Path) -> dict[str, Any]:
    path = Path(pdf)
    if not path.is_file():
        return {
            "schema": SCHEMA,
            "status": "not_available",
            "document_url": DOCUMENT_URL,
            "path": str(path.resolve()),
            "ranking_allowed": False,
            "gate_codes": ["document_not_available", "matched_numeric_export_missing"],
        }

    digest = _sha256(path)
    errors: list[str] = []
    if digest != EXPECTED_SHA256:
        errors.append(
            f"SHA-256 mismatch: expected {EXPECTED_SHA256}, got {digest}"
        )
    try:
        reader = _reader(path)
        page_count = len(reader.pages)
        if page_count != EXPECTED_PAGE_COUNT:
            errors.append(
                f"page-count mismatch: expected {EXPECTED_PAGE_COUNT}, got {page_count}"
            )
        anchors = {
            "obstructed_setup": (
                20,
                ("12 mm diameter nozzle", "0.75 m above the ground", "obstacle array"),
            ),
            "receptor_coordinates": (
                21,
                ("1.4, 2.5 and 4.4 m downstream", "Sensor Location", "0.25", "0.50"),
            ),
            "selected_models": (
                23,
                ("FLACS v10.4", "CFX v19.0", "OpenFOAM v1912+", "KFX 3.6.8"),
            ),
            "obstructed_concentration_result": (
                39,
                ("FRED results give the closest agreement", "Figure 16"),
            ),
            "obstructed_radial_result": (
                40,
                ("Figure 17", "obstacle array", "radial profiles"),
            ),
            "source_term_conclusion": (
                92,
                ("specification of the source term is critical", "choice of CFD model"),
            ),
        }
        anchor_results = {}
        for name, (page, phrases) in anchors.items():
            text = _page_text(reader, page) if page <= page_count else ""
            missing = [phrase for phrase in phrases if phrase.lower() not in text.lower()]
            anchor_results[name] = {
                "pdf_page": page,
                "verified": not missing,
                "missing_phrases": missing,
            }
            if missing:
                errors.append(f"page {page} missing anchor(s): {', '.join(missing)}")
    except (OSError, RuntimeError, ValueError) as error:
        page_count = None
        anchor_results = {}
        errors.append(str(error))

    status = "context_evidence_ready" if not errors else "withheld"
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "document_url": DOCUMENT_URL,
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": digest,
        "page_count": page_count,
        "anchors": anchor_results,
        "errors": errors,
        "benchmark": {
            "name": "ExJet obstructed high-pressure GH2 free jet",
            "source": {
                "gas": "H2",
                "phase": "ambient-temperature gas",
                "nozzle_diameter_m": 0.012,
                "release_height_m": 0.75,
                "upstream_pressure_barg": 38.0,
                "upstream_pressure_uncertainty_bar": 0.1,
                "upstream_temperature_c": 8.0,
                "upstream_temperature_uncertainty_c": 1.0,
                "ambient_temperature_c": 15.0,
                "exact_measured_mass_flow_table_in_report": False,
            },
            "obstacle_array": {
                "cylinder_length_m": 1.5,
                "cylinder_diameter_m": 0.1,
                "approximate_vertical_spacing_m": 0.35,
                "array_length_m": 3.5,
                "geometry_source": "report schematic; no machine-readable CAD in D2.2",
            },
            "receptor_coordinates_relative_to_release_m": {
                "x": [1.4, 2.5, 4.4],
                "y_z_pairs": [
                    [0.0, -0.05], [0.0, 0.20], [0.0, 0.45],
                    [0.0, 0.70], [0.25, -0.05], [0.50, -0.05],
                ],
            },
            "published_models": [
                "CFX 19.0", "FLACS 10.4/10.6", "OpenFOAM 1812",
                "OpenFOAM 1912+", "FRED AEROPLUME",
            ],
            "published_obstacle_result": (
                "FRED is reported closest for centreline concentration even though it does "
                "not resolve the obstacles; OpenFOAM 1812 is reported closest for "
                "centreline velocity decay. CFD radial profiles reveal obstacle-induced "
                "asymmetry. These are report-level conclusions, not newly computed scores."
            ),
        },
        "comparison_policy": {
            "primary_reference": "experiment",
            "cfd_role": "contextual multi-fidelity comparator",
            "published_figure_digitised": False,
            "raw_numeric_receptor_table_available": False,
            "degali_prediction_generated": False,
            "ranking_allowed": False,
            "reason": (
                "No matched machine-readable experimental/CFD receptor export is supplied, "
                "and DEGALI's current free-field core does not resolve the cylinder array."
            ),
        },
        "scope": "gh2_obstructed_jet_cfd_context_only",
        "promotion_allowed": False,
        "quantitative_lh2_field_validation_allowed": False,
        "gate_codes": [
            "not_lh2_field_validation",
            "published_plot_without_raw_numeric_export",
            "exact_mass_flow_not_tabulated",
            "obstacle_geometry_not_machine_readable",
            "degali_free_field_core_does_not_resolve_obstacle_array",
            "matched_operator_model_export_missing",
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
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    record = audit_multifidelity_cfd_evidence(args.pdf)
    if args.output is None:
        print(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False))
    else:
        write_audit(record, args.output)
        print(f"wrote {args.output} status={record['status']}")
    return 0 if record["status"] in {"context_evidence_ready", "not_available"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
