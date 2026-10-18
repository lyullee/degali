from pathlib import Path

from tools.audit_ti021_vertical_ordering import audit


def test_vertical_ordering_audit_finds_non_monotone_ti021_group():
    root = Path(__file__).resolve().parents[1]
    result = audit(root / "ti21_input_ascii.xls")
    assert result["status"] == "diagnostic_only"
    assert result["groups_with_at_least_three_heights"] > 0
    assert result["non_monotone_group_count"] > 0
    target = next(
        group for group in result["groups"]
        if group["local_x_m"] == -14.0 and group["local_y_m"] == 14.0
    )
    assert target["nondecreasing_with_height"] is False
    assert target["nonincreasing_with_height"] is False
