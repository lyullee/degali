"""Build an author-neutral IJHE upload bundle from verified local artifacts.

The bundle intentionally refuses to imply final submission readiness while
author metadata and matched field evidence remain unresolved. It copies only
manuscript-facing source/derived files; third-party raw data are never copied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


FILES = {
    "Manuscript.docx": ROOT / "outputs" / "ijhe-manuscript-draft-2026-10-09.docx",
    "Manuscript.pdf": ROOT / "outputs" / "ijhe-manuscript-draft-2026-10-09.pdf",
    "Manuscript_source.md": ROOT / "MANUSCRIPT_LH2_VALIDATION.md",
    "Supplementary_Information.md": ROOT / "IJHE_SUPPLEMENTARY_INFORMATION.md",
    "Supplementary_Information.docx": ROOT / "outputs" / "ijhe-supplementary-information-2026-10-09.docx",
    "Supplementary_Information.pdf": ROOT / "outputs" / "ijhe-supplementary-information-2026-10-09.pdf",
    "Highlights.docx": ROOT / "outputs" / "ijhe-highlights-2026-10-09.docx",
    "Highlights.txt": ROOT / "IJHE_HIGHLIGHTS.md",
    "Cover_letter_draft.md": ROOT / "IJHE_COVER_LETTER_DRAFT.md",
    "Graphical_abstract.svg": ROOT / "IJHE_GRAPHICAL_ABSTRACT.svg",
    "Graphical_abstract.pdf": ROOT / "outputs" / "ijhe-graphical-abstract-2026-10-09.pdf",
    "Figures_combined_review.pdf": ROOT / "outputs" / "ijhe-figures-combined-2026-10-09-v2.pdf",
    "Submission_manifest.json": ROOT / "outputs" / "ijhe-submission-manifest-2026-10-09.json",
    "Submission_readiness.md": ROOT / "docs" / "ijhe-submission-readiness.md",
    "Author_intake.md": ROOT / "docs" / "ijhe-author-submission-intake.md",
    "Closeout_worksheet.md": ROOT / "docs" / "ijhe-closeout-worksheet.md",
    "Figure_table_register.md": ROOT / "docs" / "ijhe-figure-table-register.md",
    "Final_gap_audit.md": ROOT / "docs" / "ijhe-final-gap-audit-2026-10-09.md",
    "External_evidence_search.md": ROOT / "docs" / "ijhe-external-evidence-search-2026-10-09.md",
    "ELVHYS_conditional_evidence.md": ROOT / "docs" / "ijhe-elvhys-conditional-evidence-2026-10-09.md",
    "PRESLHY_Trial10_conditional_evidence.md": ROOT / "docs" / "ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md",
    "PRESLHY_Trial10_field_evidence_manifest.json": ROOT / "outputs" / "preslhy-e35-trial10-evidence-2026-10-09-v2" / "field-evidence-manifest.json",
    "PRESLHY_Trial10_field_audit_execution.json": ROOT / "outputs" / "preslhy-e35-trial10-evidence-2026-10-09-v2" / "field-audit-execution.json",
    "PRESLHY_Trial10_public_retrieval_audit.json": ROOT / "outputs" / "preslhy-kitopen-trial10-retrieval-audit-2026-10-09.json",
    "PRESLHY_public_retrieval_index.json": ROOT / "outputs" / "preslhy-public-retrieval-index-2026-10-09.json",
    "PRESLHY_conditional_manifest_index.json": ROOT / "outputs" / "preslhy-conditional-manifest-index-2026-10-09.json",
    "Center_operational_data_audit.md": ROOT / "docs" / "ijhe-center-operational-data-audit-2026-10-09.md",
    "Research_requirements.txt": ROOT / "requirements-research.txt",
}


PRESLHY_MANIFEST = (
    ROOT
    / "outputs"
    / "preslhy-e35-trial10-evidence-2026-10-09-v2"
    / "field-evidence-manifest.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_preslhy_manifest(path: Path = PRESLHY_MANIFEST) -> dict[str, object]:
    """Reject a stale or incomplete conditional evidence manifest.

    The upload bundle is allowed to carry conditional evidence, but it must not
    silently copy an older manifest that predates the explicit calibration
    status and promotion gate.  This preflight keeps the derived hand-off
    bundle aligned with the evidence gate used by the manuscript audit.
    """

    if not path.is_file():
        raise FileNotFoundError(f"PRESLHY field-evidence manifest is missing: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid PRESLHY field-evidence manifest: {error}") from error
    field = payload.get("field_evidence_manifest")
    if not isinstance(field, dict):
        raise ValueError("PRESLHY manifest is missing field_evidence_manifest")
    readiness = field.get("evidence_readiness")
    if not isinstance(readiness, dict):
        raise ValueError("PRESLHY manifest is missing evidence_readiness")
    required = {
        "sensor_calibration_status": field.get("sensor_calibration_status"),
        "evidence_status": readiness.get("status"),
        "promotion_allowed": field.get("promotion_allowed"),
        "readiness_promotion_allowed": readiness.get("promotion_allowed"),
    }
    if required["sensor_calibration_status"] is None:
        raise ValueError("PRESLHY manifest lacks explicit sensor_calibration_status")
    if required["evidence_status"] not in {"conditional", "accepted", "withheld"}:
        raise ValueError(f"unexpected PRESLHY evidence status: {required['evidence_status']!r}")
    if required["promotion_allowed"] != required["readiness_promotion_allowed"]:
        raise ValueError("PRESLHY promotion gate disagrees between manifest levels")
    if not isinstance(required["promotion_allowed"], bool):
        raise ValueError("PRESLHY promotion_allowed must be boolean")
    if required["evidence_status"] == "conditional" and required["promotion_allowed"]:
        raise ValueError("conditional PRESLHY evidence cannot be promotion-allowed")
    return {
        "status": required["evidence_status"],
        "sensor_calibration_status": required["sensor_calibration_status"],
        "promotion_allowed": required["promotion_allowed"],
        "missing_requirements": readiness.get("missing_requirements", []),
    }


def build(output: Path) -> dict[str, object]:
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"output directory not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    preslhy_gate = _validate_preslhy_manifest()

    copied: list[dict[str, object]] = []
    missing: list[str] = []
    for name, source in FILES.items():
        if not source.is_file():
            missing.append(str(source.relative_to(ROOT)))
            continue
        target = output / name
        shutil.copy2(source, target)
        copied.append({
            "name": name,
            "source": str(source.relative_to(ROOT)),
            "size": target.stat().st_size,
            "sha256": sha256(target),
        })

    for source in sorted((ROOT / "outputs" / "ijhe-figure-pdfs-2026-10-09").glob("figure-*.pdf")):
        target = output / "figures" / source.name
        target.parent.mkdir(exist_ok=True)
        shutil.copy2(source, target)
        copied.append({
            "name": f"figures/{source.name}",
            "source": str(source.relative_to(ROOT)),
            "size": target.stat().st_size,
            "sha256": sha256(target),
        })

    status = {
        "bundle_status": "draft_author_input_pending",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "author_neutral": True,
        "final_upload_allowed": False,
        "preslhy_evidence_gate": preslhy_gate,
        "remaining_gates": [
            "confirm complete author list, order, affiliations, ORCIDs and corresponding author",
            "complete author contributions, funding, competing interests, exclusive-submission and preprint/prior-publication declarations",
            "obtain custodian-confirmed sensor calibration and measured three-dimensional site geometry for the matched transient event",
            "rerun strict DOCX export, portal checks and permissions review after author edits",
        ],
        "copied_files": copied,
        "missing_files": missing,
    }
    (output / "BUNDLE_STATUS.md").write_text(
        "# IJHE upload draft bundle\n\n"
        "This is an author-neutral draft bundle for internal co-author review. "
        "It is not final-upload-ready: author/declaration placeholders and the "
        "matched field-evidence gate remain open. The conditional PRESLHY Trial 10 "
        "note, multi-event manifest index and public retrieval audits are included as derived evidence; "
        "no third-party raw workbook, report PDF or proprietary source is copied here.\n\n"
        f"The bundled PRESLHY gate is `{preslhy_gate['status']}` with "
        f"sensor calibration status `{preslhy_gate['sensor_calibration_status']}` "
        f"and promotion allowed `{str(preslhy_gate['promotion_allowed']).lower()}`.\n\n"
        "## Before upload\n\n"
        + "\n".join(f"- {item}" for item in status["remaining_gates"])
        + "\n\nThe machine-readable file hash record is `bundle-manifest.json`.\n",
        encoding="utf-8",
    )
    (output / "bundle-manifest.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return status


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.output_dir)
    print(json.dumps({
        "output": str(args.output_dir),
        "status": result["bundle_status"],
        "copied": len(result["copied_files"]),
        "missing": result["missing_files"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
