import pytest

from tools.audit_d48_source_boundary import audit, exact_screen, range_screen


def test_exact_screen_accepts_approximate_report_precision():
    assert exact_screen(284.34, 298.0)["status"] == "pass"
    assert exact_screen(250.0, 298.0)["status"] == "fail"


def test_range_screen_distinguishes_inside_near_and_failed():
    assert range_screen(95.0, 90.0, 100.0)["status"] == "pass"
    assert range_screen(105.5, 90.0, 100.0)["status"] == "near_range"
    assert range_screen(120.0, 90.0, 100.0)["status"] == "fail"


def test_audit_does_not_extrapolate_summary_to_other_trials():
    source = {"trials": [
        {"trial": 10, "pressure_loss_mass_flow_g_s": 284.34},
        {"trial": 11, "pressure_loss_mass_flow_g_s": 265.27},
        {"trial": 12, "pressure_loss_mass_flow_g_s": 105.54},
    ]}
    result = audit(source)
    assert result["decision"] == "retain_trial_specific_pressure_loss"
    assert result["stop_before_downstream"] is False
    assert result["trials_22_to_25_inferred_from_summary"] is False
    assert result["screens"]["12"]["status"] == "near_range"


def test_audit_rejects_nonphysical_flow():
    source = {"trials": [
        {"trial": 10, "pressure_loss_mass_flow_g_s": 0.0},
        {"trial": 11, "pressure_loss_mass_flow_g_s": 265.0},
        {"trial": 12, "pressure_loss_mass_flow_g_s": 95.0},
    ]}
    with pytest.raises(ValueError, match="finite and positive"):
        audit(source)
