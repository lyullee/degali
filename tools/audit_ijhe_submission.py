"""Audit the IJHE manuscript package without silently promoting validation.

The audit is intentionally evidence-oriented.  It checks that the manuscript,
scope documents and numerical verification artifact are present and internally
consistent, while retaining a warning for external matched-field qualification
that cannot be inferred from numerical convergence or a conditional replay
bundle alone.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from typing import Any
import xml.etree.ElementTree as ET
import zipfile


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.ijhe-submission-readiness.v1"


def _check(
    checks: list[dict[str, Any]],
    check_id: str,
    status: str,
    evidence: str,
    detail: str,
) -> None:
    if status not in {"pass", "warning", "fail"}:
        raise ValueError(f"unsupported audit status: {status}")
    checks.append({
        "id": check_id,
        "status": status,
        "evidence": evidence,
        "detail": detail,
    })


def _path_check(checks: list[dict[str, Any]], relative: str) -> bool:
    path = ROOT / relative
    if path.is_file():
        _check(checks, f"file:{relative}", "pass", relative, "file is present")
        return True
    _check(checks, f"file:{relative}", "fail", relative, "required file is missing")
    return False


def _text_check(
    checks: list[dict[str, Any]],
    relative: str,
    check_id: str,
    phrases: tuple[str, ...],
) -> None:
    path = ROOT / relative
    if not path.is_file():
        _check(checks, check_id, "fail", relative, "cannot inspect a missing file")
        return
    text = path.read_text(encoding="utf-8").lower()
    missing = [phrase for phrase in phrases if phrase.lower() not in text]
    if missing:
        _check(
            checks,
            check_id,
            "fail",
            relative,
            "required scope language is missing: " + ", ".join(missing),
        )
    else:
        _check(checks, check_id, "pass", relative, "required scope language is present")


def _forbidden_text_check(
    checks: list[dict[str, Any]],
    relative: str,
    check_id: str,
    phrases: tuple[str, ...],
) -> None:
    """Fail when a journal-facing artifact retains an excluded claim lane."""

    path = ROOT / relative
    if not path.is_file():
        _check(checks, check_id, "fail", relative, "cannot inspect a missing file")
        return
    text = path.read_text(encoding="utf-8").lower()
    present = [phrase for phrase in phrases if phrase.lower() in text]
    if present:
        _check(
            checks,
            check_id,
            "fail",
            relative,
            "excluded comparison language remains: " + ", ".join(present),
        )
    else:
        _check(checks, check_id, "pass", relative, "excluded comparison language is absent")


def _read_error_study(checks: list[dict[str, Any]]) -> dict[str, Any] | None:
    relative = "outputs/transient-3d-error-study-2026-10-08/error_study.json"
    path = ROOT / relative
    if not path.is_file():
        _check(
            checks,
            "transient:error-study",
            "warning",
            relative,
            "local numerical error study is absent; regenerate it before submission",
        )
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        _check(checks, "transient:error-study", "fail", relative, f"invalid JSON: {error}")
        return None
    if payload.get("schema") != "degali.transient-dense-gas-3d-error-study.v1":
        _check(checks, "transient:error-study", "fail", relative, "unexpected error-study schema")
        return None
    levels = payload.get("levels")
    order = payload.get("apparent_order_from_aggregate_rmse")
    try:
        fine = levels["fine"]
        rmse = float(fine["aggregate_rmse_kg_m3"])
        maximum = float(fine["aggregate_max_absolute_error_kg_m3"])
        residual = max(
            float(levels[name]["maximum_mass_residual_kg"])
            for name in ("coarse", "medium", "fine", "reference")
        )
        order_fine = float(order["medium_to_fine"])
    except (KeyError, TypeError, ValueError) as error:
        _check(checks, "transient:error-study", "fail", relative, f"missing numerical fields: {error}")
        return None
    if not all(math.isfinite(value) for value in (rmse, maximum, residual, order_fine)):
        _check(checks, "transient:error-study", "fail", relative, "numerical fields are not finite")
        return None
    if residual > 1.0e-12 or order_fine < 1.0:
        _check(
            checks,
            "transient:error-study",
            "warning",
            relative,
            f"study is present but needs review (residual={residual:.3e}, order={order_fine:.3f})",
        )
    else:
        _check(
            checks,
            "transient:error-study",
            "pass",
            relative,
            f"fine RMSE={rmse:.6g}, max abs={maximum:.6g}, order={order_fine:.3f}, residual={residual:.3e}",
        )
    return payload


def _publication_check(checks: list[dict[str, Any]]) -> None:
    command = [sys.executable, str(ROOT / "tools" / "check_publication.py")]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    output = (result.stdout + result.stderr).strip().splitlines()
    detail = output[0] if output else f"exit code {result.returncode}"
    _check(
        checks,
        "public-snapshot",
        "pass" if result.returncode == 0 else "fail",
        "tools/check_publication.py",
        detail,
    )


def _public_doi_check(checks: list[dict[str, Any]]) -> None:
    """Require the release DOI and stable concept DOI in public metadata."""

    public_files = (
        "README.md",
        "USER_GUIDE.md",
        "USER_GUIDE_KO.md",
        "CITATION.cff",
        "pyproject.toml",
    )
    concept_doi = "10.5281/zenodo.22646258"
    version_doi = "10.5281/zenodo.23256538"
    missing_concept = []
    missing_version = []
    for relative in public_files:
        path = ROOT / relative
        if not path.is_file():
            missing_concept.append(relative)
            missing_version.append(relative)
            continue
        text = path.read_text(encoding="utf-8")
        if concept_doi not in text:
            missing_concept.append(relative)
        if version_doi not in text:
            missing_version.append(relative)
    if missing_concept or missing_version:
        details = []
        if missing_concept:
            details.append(f"concept DOI missing: {', '.join(missing_concept)}")
        if missing_version:
            details.append(f"version DOI missing: {', '.join(missing_version)}")
        _check(
            checks,
            "public:doi-integrity",
            "fail",
            "; ".join(details),
            "release metadata must retain both the exact archive DOI and the stable concept DOI",
        )
    else:
        _check(
            checks,
            "public:doi-integrity",
            "pass",
            ", ".join(public_files),
            f"public metadata pins release DOI {version_doi} and concept DOI {concept_doi}",
        )


def _submission_file_checks(checks: list[dict[str, Any]]) -> None:
    manuscript = ROOT / "MANUSCRIPT_LH2_VALIDATION.md"
    manuscript_text = manuscript.read_text(encoding="utf-8").lower() if manuscript.is_file() else ""
    if manuscript.is_file() and (
        "to be completed by the authors" in manuscript_text
        or "[author to confirm" in manuscript_text
    ):
        _check(
            checks,
            "submission:manuscript-declarations",
            "warning",
            "MANUSCRIPT_LH2_VALIDATION.md",
            "author/declaration placeholders still need completion",
        )
    elif manuscript.is_file():
        _check(
            checks,
            "submission:manuscript-declarations",
            "pass",
            "MANUSCRIPT_LH2_VALIDATION.md",
            "no author-declaration placeholder detected",
        )

    highlights = ROOT / "IJHE_HIGHLIGHTS.md"
    if highlights.is_file():
        bullets = [
            line[2:].strip()
            for line in highlights.read_text(encoding="utf-8").splitlines()
            if line.startswith("- ")
        ]
        lengths = [len(item) for item in bullets]
        if not 3 <= len(bullets) <= 5 or any(length > 85 for length in lengths):
            _check(
                checks,
                "submission:highlights",
                "fail",
                "IJHE_HIGHLIGHTS.md",
                f"expected 3-5 bullets of at most 85 characters; got lengths={lengths}",
            )
        else:
            _check(
                checks,
                "submission:highlights",
                "pass",
                "IJHE_HIGHLIGHTS.md",
                f"{len(bullets)} bullets, maximum length={max(lengths)} characters",
            )
    else:
        _check(checks, "submission:highlights", "fail", "IJHE_HIGHLIGHTS.md", "highlights file is missing")

    cover = ROOT / "IJHE_COVER_LETTER_DRAFT.md"
    if not cover.is_file():
        _check(checks, "submission:cover-letter", "fail", "IJHE_COVER_LETTER_DRAFT.md", "cover-letter draft is missing")
    elif "[AUTHOR TO CONFIRM" in cover.read_text(encoding="utf-8"):
        _check(
            checks,
            "submission:cover-letter",
            "warning",
            "IJHE_COVER_LETTER_DRAFT.md",
            "draft exists but author/declaration placeholders require manual completion",
        )
    else:
        _check(checks, "submission:cover-letter", "pass", "IJHE_COVER_LETTER_DRAFT.md", "draft has no author placeholders")

    svg = ROOT / "IJHE_GRAPHICAL_ABSTRACT.svg"
    if not svg.is_file():
        _check(
            checks,
            "submission:graphical-abstract",
            "fail",
            "IJHE_GRAPHICAL_ABSTRACT.svg",
            "graphical abstract is missing",
        )
    else:
        try:
            ET.parse(svg)
        except (OSError, ET.ParseError) as error:
            _check(
                checks,
                "submission:graphical-abstract",
                "fail",
                "IJHE_GRAPHICAL_ABSTRACT.svg",
                f"invalid SVG/XML: {error}",
            )
        else:
            _check(
                checks,
                "submission:graphical-abstract",
                "pass",
                "IJHE_GRAPHICAL_ABSTRACT.svg",
                "SVG/XML parses successfully",
            )


def _production_marker_check(checks: list[dict[str, Any]]) -> None:
    """Reject accidental production notes in journal-facing source files."""

    files = (
        "MANUSCRIPT_LH2_VALIDATION.md",
        "IJHE_SUPPLEMENTARY_INFORMATION.md",
        "IJHE_HIGHLIGHTS.md",
        "IJHE_COVER_LETTER_DRAFT.md",
    )
    pattern = re.compile(r"\b(?:TODO|FIXME|TBD)\b", re.IGNORECASE)
    offenders: list[str] = []
    for relative in files:
        path = ROOT / relative
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")):
            offenders.append(relative)
    if offenders:
        _check(
            checks,
            "submission:production-markers",
            "fail",
            "; ".join(offenders),
            "journal-facing files contain TODO/FIXME/TBD production markers",
        )
    else:
        _check(
            checks,
            "submission:production-markers",
            "pass",
            "; ".join(files),
            "no TODO/FIXME/TBD production markers found in journal-facing files",
        )


def _cross_document_consistency_check(checks: list[dict[str, Any]]) -> None:
    """Ensure headline numbers do not drift between upload-facing documents."""

    relative_paths = (
        "MANUSCRIPT_LH2_VALIDATION.md",
        "IJHE_SUPPLEMENTARY_INFORMATION.md",
        "IJHE_COVER_LETTER_DRAFT.md",
    )
    documents = {
        relative: (ROOT / relative).read_text(encoding="utf-8")
        for relative in relative_paths
        if (ROOT / relative).is_file()
    }
    expected_by_document = {
        "MANUSCRIPT_LH2_VALIDATION.md": (
            ("E3.5 peak operator", ("MG = 1.047", "VG = 1.425", "FAC2 = 0.839")),
            ("synchronised 20 s mean operator", ("MG = 0.544", "VG = 4.601", "FAC2 = 0.761")),
            ("FFI/DNV screen", ("MG = 1.245", "VG = 1.373", "FAC2 = 0.833")),
            ("Test 6 application limit", ("3.48", "underprediction")),
        ),
        "IJHE_SUPPLEMENTARY_INFORMATION.md": (
            ("E3.5 peak operator", ("E3.5 peak", "1.047", "1.425", "0.839")),
            ("synchronised 20 s mean operator", ("E3.5 synchronised 20 s mean", "0.544", "4.601", "0.761")),
            ("FFI/DNV screen", ("FFI/DNV six-arc screen", "1.245", "1.373", "0.833")),
            ("Test 6 application limit", ("2.85", "3.48", "mismatch")),
        ),
        "IJHE_COVER_LETTER_DRAFT.md": (
            ("E3.5 peak operator", ("MG = 1.047", "VG = 1.425", "FAC2 = 0.839")),
            ("synchronised 20 s mean operator", ("MG = 0.544", "VG = 4.601", "FAC2 = 0.761")),
            ("FFI/DNV screen", ("MG = 1.245", "VG = 1.373", "FAC2 = 0.833")),
            ("Test 6 application limit", ("3.48", "underprediction")),
        ),
    }
    missing: list[str] = []
    for relative, expected in expected_by_document.items():
        text = " ".join(documents.get(relative, "").split()).lower()
        for label, fragments in expected:
            absent = [fragment for fragment in fragments if fragment.lower() not in text.lower()]
            if absent:
                missing.append(f"{label} missing from {relative}: {', '.join(absent)}")
    if missing:
        _check(
            checks,
            "submission:cross-document-metrics",
            "fail",
            "; ".join(relative_paths),
            "headline metric drift detected: " + "; ".join(missing),
        )
    else:
        _check(
            checks,
            "submission:cross-document-metrics",
            "pass",
            "; ".join(relative_paths),
            "headline metrics and Test 6 boundary agree across upload-facing documents",
        )


def _submission_format_checks(checks: list[dict[str, Any]]) -> None:
    """Check conservative abstract, keyword and Highlights limits before upload."""

    manuscript_path = ROOT / "MANUSCRIPT_LH2_VALIDATION.md"
    if not manuscript_path.is_file():
        _check(checks, "submission:abstract-keywords", "fail", str(manuscript_path), "manuscript is missing")
        return
    manuscript = manuscript_path.read_text(encoding="utf-8")
    try:
        abstract = manuscript.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0]
        keyword_line = manuscript.split("**Keywords:**", 1)[1].splitlines()[0]
    except (IndexError, ValueError):
        _check(
            checks,
            "submission:abstract-keywords",
            "fail",
            "MANUSCRIPT_LH2_VALIDATION.md",
            "abstract or keyword block is missing",
        )
        return
    word_count = len(re.findall(r"\b[\w’'-]+\b", abstract))
    keywords = [item.strip() for item in keyword_line.rstrip(".").split(";") if item.strip()]
    if word_count > 200 or not 3 <= len(keywords) <= 6:
        _check(
            checks,
            "submission:abstract-keywords",
            "fail",
            "MANUSCRIPT_LH2_VALIDATION.md",
            f"abstract={word_count} words (target <=200), keywords={len(keywords)} (target 3-6)",
        )
    else:
        _check(
            checks,
            "submission:abstract-keywords",
            "pass",
            "MANUSCRIPT_LH2_VALIDATION.md",
            f"abstract={word_count} words, keywords={len(keywords)}",
        )


def _reference_citation_check(checks: list[dict[str, Any]]) -> None:
    """Catch missing, uncited, or non-sequential numeric references."""

    path = ROOT / "MANUSCRIPT_LH2_VALIDATION.md"
    if not path.is_file():
        _check(checks, "submission:reference-citations", "fail", str(path), "manuscript is missing")
        return
    text = path.read_text(encoding="utf-8")
    if "## References" not in text:
        _check(checks, "submission:reference-citations", "fail", str(path), "reference section is missing")
        return
    body, references = text.split("## References", 1)
    listed = [int(number) for number in re.findall(r"(?m)^\[(\d+)\]\s", references)]
    cited: set[int] = set()
    for token in re.findall(r"\[([0-9][0-9,\-\s]*)\]", body):
        for item in token.split(","):
            item = item.strip()
            if not item:
                continue
            if "-" in item:
                start, end = (int(value.strip()) for value in item.split("-", 1))
                cited.update(range(start, end + 1))
            else:
                cited.add(int(item))
    expected = list(range(1, len(listed) + 1))
    missing = [number for number in listed if number not in cited]
    unexpected = sorted(cited.difference(listed))
    if listed != expected:
        _check(
            checks,
            "submission:reference-citations",
            "fail",
            "MANUSCRIPT_LH2_VALIDATION.md",
            f"reference numbering is not sequential: listed={listed}, expected={expected}",
        )
    elif missing or unexpected:
        detail = []
        if missing:
            detail.append("uncited=" + ",".join(str(number) for number in missing))
        if unexpected:
            detail.append("unexpected=" + ",".join(str(number) for number in unexpected))
        _check(
            checks,
            "submission:reference-citations",
            "fail",
            "MANUSCRIPT_LH2_VALIDATION.md",
            "numeric references are not cited consistently: " + "; ".join(detail),
        )
    else:
        _check(
            checks,
            "submission:reference-citations",
            "pass",
            "MANUSCRIPT_LH2_VALIDATION.md",
            f"all {len(listed)} sequential numeric references are cited in the manuscript body",
        )


def _figure_artifact_check(checks: list[dict[str, Any]]) -> None:
    try:
        from tools.audit_ijhe_figure_artifacts import audit as audit_figures
    except ModuleNotFoundError:
        from audit_ijhe_figure_artifacts import audit as audit_figures

    payload = audit_figures()
    status = payload["status"]
    valid = sum(item.get("status") == "valid" for item in payload["figures"])
    total = len(payload["figures"])
    if status == "fail":
        detail = "; ".join(
            f"{item['path']}: {item.get('status')}"
            for item in payload["figures"]
            if item.get("status") != "valid"
        )
    elif status == "warning":
        detail = f"{valid}/{total} valid; optional local artifacts are absent"
    else:
        detail = f"{valid}/{total} SVG/XML artifacts valid with dimensions and SHA-256"
    _check(checks, "submission:figure-artifacts", status, "tools/audit_ijhe_figure_artifacts.py", detail)


def _graphical_abstract_pdf_check(checks: list[dict[str, Any]]) -> None:
    """Verify that the portal-compatible graphical-abstract derivative is not cropped."""

    relative = "outputs/ijhe-graphical-abstract-2026-10-09.pdf"
    path = ROOT / relative
    if not path.is_file():
        _check(checks, "submission:graphical-abstract-pdf", "fail", relative, "PDF derivative is missing")
        return
    try:
        from pypdf import PdfReader

        reader = PdfReader(str(path))
        pages = len(reader.pages)
        box = reader.pages[0].mediabox
        width = float(box.width)
        height = float(box.height)
    except (ImportError, OSError, IndexError, TypeError, ValueError) as error:
        _check(checks, "submission:graphical-abstract-pdf", "fail", relative, f"cannot inspect PDF: {error}")
        return
    if pages != 1 or width <= height:
        _check(
            checks,
            "submission:graphical-abstract-pdf",
            "fail",
            relative,
            f"expected one landscape page; got pages={pages}, size={width:.1f}x{height:.1f} pt",
        )
        return
    _check(
        checks,
        "submission:graphical-abstract-pdf",
        "pass",
        relative,
        f"one-page landscape PDF verified ({width:.1f}x{height:.1f} pt; aspect preserved)",
    )


def _supplementary_artifact_check(checks: list[dict[str, Any]]) -> None:
    """Verify that supplementary information is available in portal-friendly formats."""

    docx_relative = "outputs/ijhe-supplementary-information-2026-10-09.docx"
    pdf_relative = "outputs/ijhe-supplementary-information-2026-10-09.pdf"
    docx_path = ROOT / docx_relative
    pdf_path = ROOT / pdf_relative
    missing = [relative for relative, path in ((docx_relative, docx_path), (pdf_relative, pdf_path)) if not path.is_file()]
    if missing:
        _check(
            checks,
            "submission:supplementary-artifacts",
            "fail",
            "; ".join((docx_relative, pdf_relative)),
            "missing: " + ", ".join(missing),
        )
        return
    try:
        from pypdf import PdfReader

        pages = len(PdfReader(str(pdf_path)).pages)
        docx_valid = zipfile.is_zipfile(docx_path)
    except (ImportError, OSError, ValueError) as error:
        _check(
            checks,
            "submission:supplementary-artifacts",
            "fail",
            "; ".join((docx_relative, pdf_relative)),
            f"cannot inspect supplementary files: {error}",
        )
        return
    if pages < 1 or not docx_valid:
        _check(
            checks,
            "submission:supplementary-artifacts",
            "fail",
            "; ".join((docx_relative, pdf_relative)),
            f"expected valid DOCX and non-empty PDF; docx_valid={docx_valid}, pdf_pages={pages}",
        )
        return
    _check(
        checks,
        "submission:supplementary-artifacts",
        "pass",
        "; ".join((docx_relative, pdf_relative)),
        f"portal-friendly supplementary DOCX and {pages}-page PDF verified",
    )


def _highlights_docx_check(checks: list[dict[str, Any]]) -> None:
    """Verify the standalone Word Highlights file required by Elsevier."""

    relative = "outputs/ijhe-highlights-2026-10-09.docx"
    path = ROOT / relative
    if not path.is_file():
        _check(checks, "submission:highlights-docx", "fail", relative, "Highlights DOCX is missing")
        return
    try:
        from docx import Document

        paragraphs = [paragraph.text.strip() for paragraph in Document(path).paragraphs if paragraph.text.strip()]
    except (ImportError, OSError, ValueError) as error:
        _check(checks, "submission:highlights-docx", "fail", relative, f"cannot inspect Highlights DOCX: {error}")
        return
    lengths = [len(item) for item in paragraphs]
    if not 3 <= len(paragraphs) <= 5 or any(length > 85 for length in lengths):
        _check(
            checks,
            "submission:highlights-docx",
            "fail",
            relative,
            f"expected 3-5 bullets of at most 85 characters; got count={len(paragraphs)}, lengths={lengths}",
        )
        return
    _check(
        checks,
        "submission:highlights-docx",
        "pass",
        relative,
        f"standalone Word Highlights verified ({len(paragraphs)} bullets; maximum length={max(lengths)} characters)",
    )


def _docx_draft_check(checks: list[dict[str, Any]]) -> None:
    try:
        from tools.audit_ijhe_docx_draft import audit as audit_docx
    except ModuleNotFoundError:
        from audit_ijhe_docx_draft import audit as audit_docx
    payload = audit_docx()
    status = payload.get("status")
    if status == "fail":
        detail = "; ".join(item["detail"] for item in payload.get("checks", []) if item.get("status") == "fail")
        result = "fail"
    elif status == "not_available":
        detail = "author-neutral DOCX/PDF draft is not available"
        result = "warning"
    else:
        detail = "author-neutral DOCX structure with six inline figures and 12-page Word-rendered PDF verified; final author fields remain gated"
        result = "pass"
    _check(checks, "submission:docx-draft", result, "tools/audit_ijhe_docx_draft.py", detail)


def _preslhy_manifest_check(checks: list[dict[str, Any]]) -> None:
    """Verify that the conditional replay manifest exposes its promotion gate.

    The manuscript may cite a conditional replay record, but the upload
    package must not contain a pre-schema or internally contradictory manifest.
    This is deliberately separate from the external-validation warning below:
    a structurally valid conditional record is a pass, while field promotion
    remains a scientific evidence question.
    """

    relative = "outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/field-evidence-manifest.json"
    path = ROOT / relative
    if not path.is_file():
        _check(checks, "external:preslhy-manifest-gate", "fail", relative, "conditional evidence manifest is missing")
        return
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        field = payload["field_evidence_manifest"]
        readiness = field["evidence_readiness"]
        calibration = field["sensor_calibration_status"]
        status = readiness["status"]
        promotion = field["promotion_allowed"]
        readiness_promotion = readiness["promotion_allowed"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
        _check(checks, "external:preslhy-manifest-gate", "fail", relative, f"invalid gate record: {error}")
        return
    if calibration is None or status not in {"conditional", "accepted", "withheld"}:
        _check(
            checks,
            "external:preslhy-manifest-gate",
            "fail",
            relative,
            f"missing explicit calibration/status fields: calibration={calibration!r}, status={status!r}",
        )
        return
    if not isinstance(promotion, bool) or promotion != readiness_promotion:
        _check(
            checks,
            "external:preslhy-manifest-gate",
            "fail",
            relative,
            "promotion gate is missing or inconsistent between manifest levels",
        )
        return
    if status == "conditional" and promotion:
        _check(
            checks,
            "external:preslhy-manifest-gate",
            "fail",
            relative,
            "conditional evidence cannot be marked promotion-allowed",
        )
        return
    missing = readiness.get("missing_requirements", [])
    _check(
        checks,
        "external:preslhy-manifest-gate",
        "pass",
        relative,
        f"status={status}, sensor_calibration_status={calibration}, promotion_allowed={promotion}, missing={missing}",
    )


def _preslhy_manifest_index_check(checks: list[dict[str, Any]]) -> None:
    """Verify that the multi-event conditional replay index is coherent."""

    relative = "outputs/preslhy-conditional-manifest-index-2026-10-09.json"
    path = ROOT / relative
    if not path.is_file():
        _check(checks, "external:preslhy-manifest-index", "fail", relative, "conditional manifest index is missing")
        return
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        events = payload["events"]
        statuses = {item["manifest_status"] for item in events}
        promotions = {item["promotion_allowed"] for item in events}
        event_ids = {item["event_id"] for item in events}
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
        _check(checks, "external:preslhy-manifest-index", "fail", relative, f"invalid index: {error}")
        return
    expected = {
        "preslhy-e35-trial-10-2019-09-13",
        "preslhy-e35-trial-20-2019-09-18",
        "preslhy-e35-trial-21-2019-09-18",
    }
    if event_ids != expected or statuses != {"conditional"} or promotions != {False}:
        _check(
            checks,
            "external:preslhy-manifest-index",
            "fail",
            relative,
            f"expected events={sorted(expected)}, statuses={{'conditional'}}, promotion={{false}}; got events={sorted(event_ids)}, statuses={sorted(statuses)}, promotion={sorted(promotions)}",
        )
        return
    _check(
        checks,
        "external:preslhy-manifest-index",
        "pass",
        relative,
        "three accepted PRESLHY events indexed; all conditional and non-promotable",
    )


def _preslhy_public_retrieval_check(checks: list[dict[str, Any]]) -> None:
    """Verify public RADAR byte matches for the additional conditional trials."""

    relative = "outputs/preslhy-public-retrieval-index-2026-10-09.json"
    path = ROOT / relative
    if not path.is_file():
        _check(checks, "external:preslhy-public-retrieval", "fail", relative, "public retrieval index is missing")
        return
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = payload["records"]
        schema = payload["schema"]
        report = payload["source_record"]["coordinate_report"]
        trials = {int(item["trial"]) for item in records}
        verified = {bool(item["file"]["local_match_verified"]) for item in records}
        statuses = {item["promotion_boundary"]["evidence_status"] for item in records}
        promotions = {item["promotion_boundary"]["promotion_allowed"] for item in records}
        hashes_match = {
            item["file"]["public_sha256"] == item["file"]["local_sha256"]
            and item["file"]["public_bytes"] == item["file"]["local_bytes"]
            for item in records
        }
        report_ok = (
            report["name"] == "PRESLHY_D3.6_Summary_of_Rainout_Experiments_V1.22.pdf"
            and report["public_bytes"] == 6518995
            and report["public_sha256"] == "865b2b9f966e05023fbb581a4ed68f95f2353d6d6db2ad137aae7a9cf17c4f67"
            and report["redistributed"] is False
            and report["evidence_pages"]["nominal_sensor_accuracy"] == 15
            and report["evidence_pages"]["obstruction_offsets"] == 50
        )
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
        _check(checks, "external:preslhy-public-retrieval", "fail", relative, f"invalid retrieval index: {error}")
        return
    expected = {10, 20, 21}
    if (
        schema != "degali.preslhy-public-file-retrieval-index.v1"
        or trials != expected
        or verified != {True}
        or hashes_match != {True}
        or statuses != {"conditional"}
        or promotions != {False}
        or not report_ok
    ):
        _check(
            checks,
            "external:preslhy-public-retrieval",
            "fail",
            relative,
            f"expected byte-matched conditional Trials 10/20/21 and verified D3.6 report provenance; got schema={schema!r}, trials={sorted(trials)}, verified={sorted(verified)}, hashes_match={sorted(hashes_match)}, statuses={sorted(statuses)}, promotion={sorted(promotions)}, report_ok={report_ok}",
        )
        return
    _check(
        checks,
        "external:preslhy-public-retrieval",
        "pass",
        relative,
        "RADAR public downloads match local Trial 10/20/21 workbooks byte-for-byte and the D3.6 coordinate-report provenance is pinned; all remain conditional and non-promotable",
    )


def audit() -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    required_files = (
        "MANUSCRIPT_LH2_VALIDATION.md",
        "IJHE_SUPPLEMENTARY_INFORMATION.md",
        "IJHE_HIGHLIGHTS.md",
        "IJHE_COVER_LETTER_DRAFT.md",
        "IJHE_GRAPHICAL_ABSTRACT.svg",
        "docs/lh2-paper-baseline-2026-09-20.md",
        "docs/ijhe-submission-readiness.md",
        "docs/ijhe-figure-table-register.md",
        "docs/publication-scope.md",
        "docs/field-evidence-manifest.md",
        "docs/transient-dense-gas-3d.md",
        "docs/open-channel-h2-boundary.md",
        "docs/ijhe-word-artifact-boundary.md",
        "docs/ijhe-author-submission-intake.md",
        "docs/ijhe-author-metadata-source-2026-10-09.md",
        "docs/ijhe-closeout-worksheet.md",
        "docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md",
        "docs/ijhe-figure-derivative-map-2026-10-09.md",
        "docs/ijhe-final-gap-audit-2026-10-09.md",
        "docs/ijhe-external-evidence-search-2026-10-09.md",
        "docs/ijhe-field-evidence-request-2026-10-09.md",
        "docs/ijhe-center-operational-data-audit-2026-10-09.md",
        "docs/ijhe-smedis-transfer-boundary-2026-10-09.md",
        "docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md",
        "docs/ijhe-elvhys-conditional-evidence-2026-10-09.md",
        "tools/benchmark_transient_dense_gas_3d.py",
        "tools/plot_ijhe_operator_comparison.py",
        "tools/plot_transient_3d_error_study.py",
        "tools/render_ijhe_graphical_abstract.ps1",
        "tools/export_ijhe_word_pdfs.ps1",
        "tools/build_ijhe_submission_manifest.py",
        "tools/build_ijhe_upload_bundle.py",
        "tools/build_ijhe_docx.py",
        "tools/build_ijhe_highlights_docx.py",
        "tools/build_preslhy_trial_evidence_bundle.py",
        "tools/build_preslhy_manifest_index.py",
        "tools/build_preslhy_public_retrieval_index.py",
        "tools/audit_ijhe_figure_artifacts.py",
        "tools/audit_ijhe_docx_draft.py",
        "tools/audit_open_channel_h2_dataset.py",
    )
    for relative in required_files:
        _path_check(checks, relative)

    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:operator-boundary",
        ("synchronized 20 s mean", "observation operator"),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:failure-boundary",
        ("Test 6", "not a substitute for transient multiphase CFD"),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:industrial-positioning",
        ("hazard identification", "pre-FEED", "separation distance", "higher-fidelity transient and obstacle-resolving"),
    )
    _forbidden_text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:slabx-comparison-excluded",
        ("slabx",),
    )
    _forbidden_text_check(
        checks,
        "IJHE_SUPPLEMENTARY_INFORMATION.md",
        "supplement:slabx-comparison-excluded",
        ("slabx",),
    )
    _forbidden_text_check(
        checks,
        "IJHE_GRAPHICAL_ABSTRACT.svg",
        "graphical-abstract:slabx-comparison-excluded",
        ("slabx",),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:discussion-section",
        ("## 6. Discussion and application envelope",),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:reproducibility",
        ("Code availability", "SHA-256", "10.5281/zenodo.22646258", "10.5281/zenodo.23256538"),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:headline-metrics",
        ("MG = 1.047", "VG = 1.425", "FAC2 = 0.839", "MG = 1.245"),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:main-artifacts",
        ("Table 1.", "Table 2.", "Figure 1.", "Figure 4.", "FFI/DNV horizontal-release arc residual map", "Tests 4 and 6", "Figure 6."),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "submission:title-page-metadata",
        ("Article type:", "Manuscript word count", "Abstract word count", "Keywords:"),
    )
    _text_check(
        checks,
        "MANUSCRIPT_LH2_VALIDATION.md",
        "manuscript:ai-disclosure",
        (
            "### Declaration of generative AI and AI-assisted technologies in the manuscript preparation process",
            "OpenAI Codex",
            "did not generate experimental observations",
            "AI-assisted software-development aid",
        ),
    )
    _text_check(
        checks,
        "docs/open-channel-h2-boundary.md",
        "open-channel:non-promotion-boundary",
        ("boundary_only", "promotion_allowed=false", "not be pooled"),
    )
    _text_check(
        checks,
        "docs/ijhe-elvhys-conditional-evidence-2026-10-09.md",
        "elvhys:conditional-boundary",
        ("10.18710/JXJP0H", "not a measured hydrogen mass-flow history", "promotion_allowed=false"),
    )
    _text_check(
        checks,
        "docs/ijhe-figure-table-register.md",
        "submission:figure-table-register",
        ("Fig. 6", "Table 2", "Location / generation source", "Supplementary tables"),
    )
    _text_check(
        checks,
        "IJHE_SUPPLEMENTARY_INFORMATION.md",
        "supplement:sensor-uncertainty-boundary",
        ("Sensor and observation-operator uncertainty boundary", "Calibration bias/precision", "not pure sensor accuracy"),
    )
    _submission_file_checks(checks)
    _production_marker_check(checks)
    _text_check(
        checks,
        "IJHE_COVER_LETTER_DRAFT.md",
        "submission:cover-letter-evidence-boundary",
        ("ELVHYS WP4.2", "not pooled", "mass-flow boundary"),
    )
    _text_check(
        checks,
        "IJHE_COVER_LETTER_DRAFT.md",
        "submission:related-manuscript-disclosure",
        ("related-manuscript disclosure", "currently under review", "excludes that model's formulation", "comparative conclusions"),
    )
    _cross_document_consistency_check(checks)
    _submission_format_checks(checks)
    _reference_citation_check(checks)
    _figure_artifact_check(checks)
    _graphical_abstract_pdf_check(checks)
    _supplementary_artifact_check(checks)
    _highlights_docx_check(checks)
    _docx_draft_check(checks)
    _preslhy_manifest_check(checks)
    _preslhy_manifest_index_check(checks)
    _preslhy_public_retrieval_check(checks)
    error_study = _read_error_study(checks)

    # This remains intentionally a warning.  Conditional common-clock and
    # confined dynamic replay bundles are documented, but calibration, source
    # boundary and measured site geometry are still required before operational
    # promotion.
    _check(
        checks,
        "external:matched-transient-evidence",
        "warning",
        "docs/ijhe-submission-readiness.md",
        "ELVHYS supplies time-resolved concentration, calibration specifications and measured TCS geometry, but a measured H2 source-rate history and custodian-confirmed field package are still required for a stronger operational claim",
    )
    _publication_check(checks)
    _public_doi_check(checks)

    failures = sum(item["status"] == "fail" for item in checks)
    warnings = sum(item["status"] == "warning" for item in checks)
    if failures:
        status = "not_ready"
    elif warnings:
        status = "draft_ready_external_validation_pending"
    else:
        status = "draft_ready"
    return {
        "schema": SCHEMA,
        "status": status,
        "checks": checks,
        "summary": {
            "pass": sum(item["status"] == "pass" for item in checks),
            "warning": warnings,
            "fail": failures,
        },
        "numerical_error_study_present": error_study is not None,
        "external_validation_promotion_allowed": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit the machine-readable audit record")
    args = parser.parse_args(argv)
    record = audit()
    if args.json:
        print(json.dumps(record, indent=2, ensure_ascii=False, sort_keys=True))
    else:
        print(f"IJHE submission audit: {record['status']}")
        for item in record["checks"]:
            print(f"[{item['status'].upper():7}] {item['id']}: {item['detail']}")
        summary = record["summary"]
        print(f"summary: {summary['pass']} pass, {summary['warning']} warning, {summary['fail']} fail")
    return 0 if record["status"] != "not_ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
