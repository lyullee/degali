import runpy

import numpy as np
import pytest


TOOL = runpy.run_path("tools/audit_pressure_open_boundary_reflection.py")


def test_reflection_classifier_uses_predeclared_return_window_and_limit():
    size = TOOL["STEPS"] + 1
    reference = np.zeros(size)
    short = np.zeros(size)
    reference[130:150] = 1.0e-4
    short[:] = reference
    preliminary = TOOL["classify_reflection"](short, reference)
    start, stop = preliminary["reflection_window_steps"]
    short[(start + stop) // 2] = 4.0e-6

    result = TOOL["classify_reflection"](short, reference)

    assert result["reflection_amplitude_ratio"] == pytest.approx(0.04)
    assert result["reflection_amplitude_ratio_limit"] == 0.05
    assert result["passed"] is True
    assert result["promotion_allowed"] is False


def test_reflection_classifier_fails_above_limit():
    size = TOOL["STEPS"] + 1
    reference = np.zeros(size)
    short = np.zeros(size)
    reference[130:150] = 1.0e-4
    preliminary = TOOL["classify_reflection"](short, reference)
    start, stop = preliminary["reflection_window_steps"]
    short[(start + stop) // 2] = 6.0e-6

    result = TOOL["classify_reflection"](short, reference)

    assert result["reflection_amplitude_ratio"] == pytest.approx(0.06)
    assert result["passed"] is False
