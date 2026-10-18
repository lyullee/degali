from pathlib import Path

from tools.audit_dataverse_no_obstructed_release import audit


ROOT = Path(__file__).resolve().parents[1]


def test_dataverse_no_source_pressure_component_is_registered_but_not_promoted():
    result = audit(ROOT / "reference" / "obstacle_wake" / "dataverse_no_wskbij")
    assert result["status"] == "component_source_pressure_ready"
    assert result["promotion_allowed"] is False
    assert result["concentration_validation_allowed"] is False
    assert result["workbook"]["experiment_row_count"] == 51
    assert result["workbook"]["experiment_count_with_complete_pressure_rows"] == 39
    assert result["metadata_gates"]["ambient_weather_missing"] is True
