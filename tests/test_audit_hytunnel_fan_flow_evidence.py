from pathlib import Path
import runpy

import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = runpy.run_path(ROOT / "tools/audit_hytunnel_fan_flow_evidence.py")


def test_public_fan_workbook_permits_mean_profile_but_withholds_turbulence():
    record = TOOL["audit_fan_flow_evidence"]()
    plane = record["selected_first_average_plane"]

    assert record["status"] == "conditional_mean_profile_only_turbulence_withheld"
    assert record["mean_profile_candidate_allowed"] is True
    assert record["turbulence_inlet_qualification_allowed"] is False
    assert record["promotion_allowed"] is False
    assert record["dataset"]["selected_configuration_id"].endswith("sheet-01")
    assert record["dataset"]["available_spatial_statistics"] == [
        "average", "maximum", "minimum",
    ]
    assert plane["x_m"] == pytest.approx(0.25)
    assert plane["relative_y_m"] == pytest.approx(
        [-0.5, -0.25, 0.0, 0.25, 0.5]
    )
    assert plane["minimum_velocity_m_s"] == pytest.approx(2.902463414634146)
    assert plane["mean_velocity_m_s"] == pytest.approx(3.57159837398374)
    assert plane["maximum_velocity_m_s"] == pytest.approx(4.240520325203254)
    assert "velocity_time_series" in record["missing_required_turbulence_evidence"]
