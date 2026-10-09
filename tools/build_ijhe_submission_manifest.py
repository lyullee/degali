"""Build a hash-pinned manifest for the IJHE upload-facing source package."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from audit_ijhe_submission import audit


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.ijhe-submission-manifest.v1"

SOURCE_FILES = (
    "MANUSCRIPT_LH2_VALIDATION.md",
    "IJHE_SUPPLEMENTARY_INFORMATION.md",
    "IJHE_HIGHLIGHTS.md",
    "IJHE_COVER_LETTER_DRAFT.md",
    "IJHE_GRAPHICAL_ABSTRACT.svg",
    "docs/ijhe-submission-readiness.md",
    "docs/ijhe-figure-table-register.md",
    "docs/ijhe-author-submission-intake.md",
    "docs/ijhe-author-metadata-source-2026-10-09.md",
    "docs/ijhe-closeout-worksheet.md",
    "docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md",
    "docs/ijhe-figure-derivative-map-2026-10-09.md",
    "docs/ijhe-word-artifact-boundary.md",
    "docs/ijhe-final-gap-audit-2026-10-09.md",
    "docs/ijhe-external-evidence-search-2026-10-09.md",
    "docs/ijhe-field-evidence-request-2026-10-09.md",
    "docs/ijhe-center-operational-data-audit-2026-10-09.md",
    "docs/ijhe-smedis-transfer-boundary-2026-10-09.md",
    "docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md",
    "docs/ijhe-elvhys-conditional-evidence-2026-10-09.md",
    "tools/audit_ijhe_submission.py",
    "tools/audit_ijhe_figure_artifacts.py",
    "tools/plot_ijhe_operator_comparison.py",
    "tools/plot_transient_3d_error_study.py",
    "tools/render_ijhe_graphical_abstract.ps1",
    "tools/export_ijhe_word_pdfs.ps1",
    "tools/build_ijhe_submission_manifest.py",
    "tools/build_ijhe_docx.py",
    "tools/build_ijhe_highlights_docx.py",
    "tools/build_ijhe_upload_bundle.py",
    "tools/build_preslhy_trial_evidence_bundle.py",
    "tools/build_preslhy_manifest_index.py",
    "tools/build_preslhy_public_retrieval_index.py",
    "tools/audit_ijhe_docx_draft.py",
    "src/degali/addons/field_evidence_manifest.py",
    "src/degali/addons/__init__.py",
    "src/degali/cli.py",
    "CITATION.cff",
    "pyproject.toml",
    "requirements-research.txt",
    "tmp/ijhe-graphical-abstract-print.html",
)

DERIVED_FILES = (
    "outputs/ijhe-figure-artifact-audit-2026-10-09.json",
    "outputs/smedis-transfer-field-audit-2026-10-09.json",
    "outputs/ijhe-operator-comparison-2026-10-08.svg",
    "outputs/transient-3d-error-study-2026-10-08/transient-3d-convergence.svg",
    "outputs/ijhe-manuscript-draft-2026-10-09.docx",
    "outputs/ijhe-manuscript-draft-2026-10-09.pdf",
    "outputs/ijhe-supplementary-information-2026-10-09.docx",
    "outputs/ijhe-supplementary-information-2026-10-09.pdf",
    "outputs/ijhe-highlights-2026-10-09.docx",
    "outputs/ijhe-graphical-abstract-2026-10-09.pdf",
    "outputs/ijhe-figures-combined-2026-10-09-v2.pdf",
    "outputs/ijhe-main-figure-images-2026-10-09-clean/figure-1.png",
    "outputs/ijhe-main-figure-images-2026-10-09-clean/figure-2.png",
    "outputs/ijhe-main-figure-images-2026-10-09-clean/figure-3.png",
    "outputs/ijhe-main-figure-images-2026-10-09-clean/figure-4.png",
    "outputs/ijhe-main-figure-images-2026-10-09-clean/figure-5.png",
    "outputs/ijhe-main-figure-images-2026-10-09-clean/figure-6.png",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-1.pdf",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-2.pdf",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-3.pdf",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-4.pdf",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-5.pdf",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-6.pdf",
    "outputs/ijhe-figure-pdfs-2026-10-09/figure-7.pdf",
    "outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/bundle_metadata.json",
    "outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/field-audit-execution.json",
    "outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/field-evidence-manifest.json",
    "outputs/preslhy-kitopen-trial10-retrieval-audit-2026-10-09.json",
    "outputs/preslhy-public-retrieval-index-2026-10-09.json",
    "outputs/preslhy-conditional-manifest-index-2026-10-09.json",
    "tmp/ijhe-lh2-paper-figures/preslhy-concentration-residual-map.svg",
    "tmp/ijhe-lh2-paper-figures/ffi-arc-residual-map.svg",
    "tmp/ijhe-time-aligned/e35_lfl_uncertainty_envelope.svg",
    "tmp/ijhe-decision-scenarios/decision_lfl_distances.svg",
)


def _file_record(relative: str) -> dict[str, object]:
    path = ROOT / relative
    if not path.is_file():
        return {"path": relative, "present": False}
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": relative,
        "present": True,
        "bytes": path.stat().st_size,
        "sha256": digest,
    }


def build(output: Path) -> dict[str, object]:
    readiness = audit()
    payload: dict[str, object] = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "package_status": readiness["status"],
        "audit_summary": readiness["summary"],
        "source_files": [_file_record(relative) for relative in SOURCE_FILES],
        "derived_figure_files": [_file_record(relative) for relative in DERIVED_FILES],
        "redistribution_boundary": (
            "Derived figures and hashes may be retained for submission control; "
            "third-party raw workbooks, report PDFs and original Fortran source "
            "are not redistributed by this manifest."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = build(args.output)
    missing = sum(not item["present"] for item in payload["source_files"])
    print(json.dumps({"output": str(args.output), "source_missing": missing, "status": payload["package_status"]}))
    return 0 if missing == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
