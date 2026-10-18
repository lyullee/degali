from pathlib import Path

import numpy as np

from tools.audit_ti021_time_origin import _summary, audit
from tools.score_thorney_ti021_initial_cloud_sweep import (
    _manifest_clock_alignment,
    score_run,
)


def test_summary_reports_absolute_timing_error():
    result = _summary(np.asarray([-2.0, 4.0, 6.0]))
    assert result["count"] == 3
    assert result["mean_s"] == 8.0 / 3.0
    assert result["mae_s"] == 4.0


def test_ti021_time_origin_audit_is_explicitly_diagnostic():
    root = Path(__file__).resolve().parents[1]
    result = audit(root / "ti21_input_ascii.xls")
    assert result["status"] == "diagnostic_only"
    assert result["positive_dose_arrival_row_count"] == 82
    assert result["diagnostic_only"] is True
    assert result["descriptive_linear_fit"]["effective_speed_m_s"] > 0.0
    assert result["advection_sanity_by_speed"]["3.4"][
        "residual_observed_minus_local_y_over_wind_s"
    ]["median_s"] > 0.0


def test_score_time_origin_offset_is_separate_from_raw_error():
    root = Path(__file__).resolve().parents[1]
    run = root / "outputs" / "thorney-ti021-inner-10s-dx5-workbookmass10s-sigma332p5-2026-10-15"
    result = score_run(
        root / "ti21_input_ascii.xls", run,
        diagnostic_time_origin_offset_s=17.25,
    )
    assert result["mean_absolute_error_s"] == 16.6
    assert result["time_origin_corrected_mean_absolute_error_s"] == 6.15
    assert result["time_origin_correction_status"] == "observation_operator_sensitivity_only"


def test_score_relative_peak_operator_is_explicitly_diagnostic():
    root = Path(__file__).resolve().parents[1]
    run = root / "outputs" / "thorney-ti021-inner-10s-dx5-workbookmass10s-sigma332p5-2026-10-15"
    result = score_run(
        root / "ti21_input_ascii.xls", run,
        relative_peak_fraction=0.05,
    )
    assert result["relative_peak_fraction"] == 0.05
    assert result["relative_operator_status"] == "diagnostic_only_truncated_peak_window"
    assert result["relative_model_detected_among_available_count"] >= 0


def test_manifest_clock_alignment_fails_closed_when_readiness_is_conditional(
    monkeypatch, tmp_path,
):
    class FakeManifest:
        def quantitative_readiness_record(self):
            return {
                "status": "conditional",
                "clock_alignment": {
                    "status": "accepted",
                    "offset_s": 17.25,
                    "uncertainty_s": 0.2,
                },
            }

    monkeypatch.setattr(
        "tools.score_thorney_ti021_initial_cloud_sweep.read_field_evidence_manifest_json",
        lambda path: FakeManifest(),
    )
    result = _manifest_clock_alignment(tmp_path / "manifest.json")
    assert result["status"] == "withheld"
    assert "not accepted" in result["reason"]


def test_manifest_clock_alignment_fails_closed_when_file_is_missing(tmp_path):
    result = _manifest_clock_alignment(tmp_path / "missing-manifest.json")
    assert result["status"] == "withheld"
    assert "could not be verified" in result["reason"]
