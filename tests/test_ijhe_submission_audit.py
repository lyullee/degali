import json

from tools.audit_ijhe_submission import audit, main


def test_ijhe_audit_keeps_external_validation_pending_without_failing_draft(tmp_path, capsys):
    record = audit()

    assert record["schema"] == "degali.ijhe-submission-readiness.v1"
    assert record["status"] in {"draft_ready_external_validation_pending", "draft_ready"}
    assert record["external_validation_promotion_allowed"] is False
    assert any(
        item["id"] == "external:matched-transient-evidence"
        and item["status"] == "warning"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:cross-document-metrics"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:abstract-keywords"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:title-page-metadata"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:reference-citations"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:figure-artifacts"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:graphical-abstract-pdf"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:supplementary-artifacts"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:highlights-docx"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "manuscript:industrial-positioning"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "manuscript:slabx-comparison-excluded"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "supplement:slabx-comparison-excluded"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "graphical-abstract:slabx-comparison-excluded"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:related-manuscript-disclosure"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "manuscript:main-artifacts"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "file:tools/build_ijhe_docx.py"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:docx-draft"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:manuscript-declarations"
        and item["status"] == "warning"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "external:preslhy-manifest-gate"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "external:preslhy-public-retrieval"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "public:doi-integrity"
        and item["status"] == "pass"
        for item in record["checks"]
    )
    assert any(
        item["id"] == "submission:production-markers"
        and item["status"] == "pass"
        for item in record["checks"]
    )

    assert main(["--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema"] == record["schema"]
