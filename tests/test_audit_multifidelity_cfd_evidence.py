from pathlib import Path
import runpy

import pytest


PDF = Path(
    "reference/datasets/multhyfuel-d22/raw/"
    "MHYF_WP2_D22_Assessment_of_dispersion_for_high-pressure_H2_20221111.pdf"
)


@pytest.mark.skipif(not PDF.is_file(), reason="requires locally controlled MultHyFuel D2.2")
def test_multifidelity_audit_keeps_published_plot_out_of_ranking():
    tool = runpy.run_path("tools/audit_multifidelity_cfd_evidence.py")
    record = tool["audit_multifidelity_cfd_evidence"](PDF)

    assert record["status"] == "context_evidence_ready"
    assert record["page_count"] == 105
    assert all(item["verified"] for item in record["anchors"].values())
    assert record["comparison_policy"]["primary_reference"] == "experiment"
    assert record["comparison_policy"]["ranking_allowed"] is False
    assert record["promotion_allowed"] is False
    assert record["quantitative_lh2_field_validation_allowed"] is False


def test_multifidelity_audit_reports_missing_file_without_claim(tmp_path):
    tool = runpy.run_path("tools/audit_multifidelity_cfd_evidence.py")
    record = tool["audit_multifidelity_cfd_evidence"](tmp_path / "missing.pdf")

    assert record["status"] == "not_available"
    assert record["ranking_allowed"] is False
