import json
from pathlib import Path

import pytest

from tools.build_ijhe_submission_manifest import build


pytestmark = pytest.mark.skipif(
    not Path("outputs/ijhe-manuscript-draft-2026-10-09.docx").is_file(),
    reason="requires the local IJHE submission artifact set",
)


def test_ijhe_submission_manifest_pins_source_files(tmp_path):
    output = tmp_path / "ijhe-manifest.json"
    payload = build(output)

    assert payload["schema"] == "degali.ijhe-submission-manifest.v1"
    assert payload["package_status"] == "draft_ready_external_validation_pending"
    assert payload["audit_summary"]["fail"] == 0
    assert all(item["present"] for item in payload["source_files"])
    assert all(len(item["sha256"]) == 64 for item in payload["source_files"])
    assert json.loads(output.read_text(encoding="utf-8"))["schema"] == payload["schema"]
