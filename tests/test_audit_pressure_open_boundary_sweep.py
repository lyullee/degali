import runpy

import pytest


TOOL = runpy.run_path("tools/audit_pressure_open_boundary_sweep.py")


def test_case_matrix_spans_declared_angles_and_wavelengths():
    cases = TOOL["build_case_matrix"]()
    assert len(cases) == 8
    assert {case["angle_degrees_from_normal"] for case in cases} == {
        15.0, 30.0, 45.0, 60.0,
    }
    assert {case["tangential_mode"] for case in cases} == {1, 2}
    assert all(case["time_steps"] > 0 for case in cases)
    assert any(not case["resolved_for_gate"] for case in cases)


def test_wavelength_geometry_matches_packet_direction():
    case = TOOL["_case_definition"](30.0, 1)
    assert case["wavelength_cells"] == pytest.approx(16.0)
    assert case["resolved_for_gate"] is True
    threshold_case = TOOL["_case_definition"](30.0, 2)
    assert threshold_case["wavelength_cells"] == pytest.approx(8.0)
    assert threshold_case["resolved_for_gate"] is True


def test_sweep_gate_excludes_underresolved_case_but_does_not_hide_it():
    cases = [
        {
            "resolved_for_gate": True,
            "case_gate_passed": True,
            "reflection_amplitude_ratio": 0.04,
        },
        {
            "resolved_for_gate": False,
            "case_gate_passed": False,
            "reflection_amplitude_ratio": 0.20,
        },
    ]
    result = TOOL["classify_sweep"](cases)
    assert result["passed"] is True
    assert result["in_scope_case_count"] == 1
    assert result["excluded_underresolved_case_count"] == 1
    assert result["promotion_allowed"] is False


def test_sweep_gate_fails_when_any_resolved_case_fails():
    cases = [
        {
            "resolved_for_gate": True,
            "case_gate_passed": True,
            "reflection_amplitude_ratio": 0.04,
        },
        {
            "resolved_for_gate": True,
            "case_gate_passed": False,
            "reflection_amplitude_ratio": 0.06,
        },
    ]
    result = TOOL["classify_sweep"](cases)
    assert result["passed"] is False
    assert result["status"] == "resolved_angle_wavelength_sweep_gate_failed"
