import runpy

import numpy as np
import pytest


TOOL = runpy.run_path(
    "tools/audit_pressure_open_boundary_oblique_reflection.py"
)


def _synthetic_traces(reflected_ratio):
    size = TOOL["STEPS"] + 1
    reference = np.zeros(size, dtype=complex)
    short = np.zeros(size, dtype=complex)
    reference[140:180] = 1.0e-4 * np.exp(0.3j)
    short[:] = reference
    preliminary = TOOL["classify_oblique_reflection"](short, reference)
    start, stop = preliminary["reflection_window_steps"]
    short[(start + stop) // 2] = reflected_ratio * 1.0e-4 * np.exp(-0.4j)
    return short, reference


def test_oblique_classifier_passes_below_predeclared_limit():
    short, reference = _synthetic_traces(0.04)
    result = TOOL["classify_oblique_reflection"](short, reference)
    assert result["reflection_amplitude_ratio"] == pytest.approx(0.04)
    assert result["reflection_amplitude_ratio_limit"] == 0.05
    assert result["passed"] is True
    assert result["promotion_allowed"] is False


def test_oblique_classifier_fails_above_predeclared_limit():
    short, reference = _synthetic_traces(0.06)
    result = TOOL["classify_oblique_reflection"](short, reference)
    assert result["reflection_amplitude_ratio"] == pytest.approx(0.06)
    assert result["passed"] is False


def test_oblique_packet_wavenumbers_match_declared_angle():
    tangential, normal = TOOL["_wavenumbers"]()
    angle = np.degrees(np.arctan2(tangential, normal))
    assert angle == pytest.approx(TOOL["INCIDENCE_ANGLE_DEGREES"])
