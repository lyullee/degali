from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from tools.build_ijhe_upload_bundle import _validate_preslhy_manifest, build


@pytest.mark.skipif(
    not Path(
        "outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/field-evidence-manifest.json"
    ).is_file(),
    reason="requires the local IJHE submission and conditional-evidence artifacts",
)
def test_ijhe_upload_bundle_is_author_neutral_and_hash_pinned(tmp_path):
    result = build(tmp_path / "bundle")

    assert result["bundle_status"] == "draft_author_input_pending"
    assert result["author_neutral"] is True
    assert result["final_upload_allowed"] is False
    assert result["missing_files"] == []
    assert len(result["copied_files"]) == 39
    assert any(
        item["name"] == "External_evidence_search.md"
        for item in result["copied_files"]
    )
    assert any(
        item["name"] == "PRESLHY_Trial10_public_retrieval_audit.json"
        for item in result["copied_files"]
    )
    assert any(
        item["name"] == "PRESLHY_conditional_manifest_index.json"
        for item in result["copied_files"]
    )
    assert any(
        item["name"] == "PRESLHY_public_retrieval_index.json"
        for item in result["copied_files"]
    )
    assert result["preslhy_evidence_gate"] == {
        "status": "conditional",
        "sensor_calibration_status": "specification_only",
        "promotion_allowed": False,
        "missing_requirements": ["sensor_calibration_certificate"],
    }

    bundle = tmp_path / "bundle"
    manifest = json.loads((bundle / "bundle-manifest.json").read_text(encoding="utf-8"))
    for item in manifest["copied_files"]:
        digest = hashlib.sha256((bundle / item["name"]).read_bytes()).hexdigest()
        assert digest == item["sha256"]
    assert (bundle / "Graphical_abstract.pdf").is_file()
    assert (bundle / "Highlights.docx").is_file()
    assert (bundle / "Supplementary_Information.docx").is_file()
    assert (bundle / "Supplementary_Information.pdf").is_file()
    assert (bundle / "Center_operational_data_audit.md").is_file()
    assert (bundle / "Research_requirements.txt").is_file()
    assert (bundle / "Closeout_worksheet.md").is_file()
    assert (bundle / "PRESLHY_Trial10_conditional_evidence.md").is_file()
    assert (bundle / "PRESLHY_Trial10_field_evidence_manifest.json").is_file()
    assert (bundle / "PRESLHY_Trial10_field_audit_execution.json").is_file()
    assert (bundle / "PRESLHY_conditional_manifest_index.json").is_file()
    assert (bundle / "PRESLHY_public_retrieval_index.json").is_file()
    assert (bundle / "Multifidelity_evidence.md").is_file()
    assert (bundle / "USN_enclosure_validation_report.md").is_file()
    assert (bundle / "USN_enclosure_validation_metrics.json").is_file()
    assert (bundle / "figures" / "figure-8.pdf").is_file()
    status_text = (bundle / "BUNDLE_STATUS.md").read_text(encoding="utf-8")
    assert "third-party raw" in status_text
    assert "sensor calibration status `specification_only`" in status_text


@pytest.mark.skipif(
    not Path("outputs/ijhe-manuscript-author-approved-2026-10-10.docx").is_file(),
    reason="requires the strict author-approved IJHE artifacts",
)
def test_ijhe_upload_bundle_can_close_only_the_author_gate(tmp_path):
    result = build(tmp_path / "approved-bundle", author_approved=True)

    assert result["bundle_status"] == "author_approved_evidence_bounded"
    assert result["author_neutral"] is False
    assert result["author_approval_recorded"] is True
    assert result["final_upload_allowed"] is False
    assert result["preslhy_evidence_gate"]["status"] == "conditional"
    assert len(result["copied_files"]) == 41
    bundle = tmp_path / "approved-bundle"
    assert (bundle / "Author_approval.md").is_file()
    assert (bundle / "Cover_letter.docx").is_file()
    assert (bundle / "Manuscript.docx").read_bytes() == Path(
        "outputs/ijhe-manuscript-author-approved-2026-10-10.docx"
    ).read_bytes()
    status_text = (bundle / "BUNDLE_STATUS.md").read_text(encoding="utf-8")
    assert "author metadata and declarations in this bundle are approved" in status_text
    assert "matched field-evidence gate remains conditional" in status_text


def test_ijhe_upload_bundle_rejects_stale_manifest(tmp_path):
    stale = tmp_path / "stale-manifest.json"
    stale.write_text(
        json.dumps({
            "field_evidence_manifest": {
                "evidence_readiness": {
                    "status": "accepted",
                    "promotion_allowed": True,
                },
                "promotion_allowed": True,
            }
        }),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="sensor_calibration_status"):
        _validate_preslhy_manifest(stale)
