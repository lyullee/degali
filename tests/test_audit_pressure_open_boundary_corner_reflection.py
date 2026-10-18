import runpy

import numpy as np
import pytest


TOOL = runpy.run_path(
    "tools/audit_pressure_open_boundary_corner_reflection.py"
)


def _synthetic_traces(reflected_ratio):
    size = TOOL["STEPS"] + 1
    reference = np.zeros(size)
    short = np.zeros(size)
    preliminary = TOOL["classify_corner_reflection"](short, reference + 1.0e-4)
    direct_start, direct_stop = preliminary["incident_window_steps"]
    reference[(direct_start + direct_stop) // 2] = 1.0e-4
    short[:] = reference
    preliminary = TOOL["classify_corner_reflection"](short, reference)
    start, stop = preliminary["reflection_window_steps"]
    short[(start + stop) // 2] += reflected_ratio * 1.0e-4
    return short, reference


def test_corner_classifier_passes_below_predeclared_limit():
    short, reference = _synthetic_traces(0.04)
    result = TOOL["classify_corner_reflection"](short, reference)
    assert result["reflection_amplitude_ratio"] == pytest.approx(0.04)
    assert result["reflection_amplitude_ratio_limit"] == 0.05
    assert result["passed"] is True


def test_corner_classifier_fails_above_predeclared_limit():
    short, reference = _synthetic_traces(0.06)
    result = TOOL["classify_corner_reflection"](short, reference)
    assert result["reflection_amplitude_ratio"] == pytest.approx(0.06)
    assert result["passed"] is False


def test_corner_packet_has_declared_diagonal_acoustic_velocity():
    distribution = TOOL["_initial_distribution"](TOOL["SHORT_SIZE"])
    density = np.sum(distribution, axis=0)
    cx = np.asarray([item[0] for item in TOOL["_equilibrium_d3q19"].__globals__["_C3"]])
    cy = np.asarray([item[1] for item in TOOL["_equilibrium_d3q19"].__globals__["_C3"]])
    ux = np.tensordot(cx, distribution, axes=(0, 0)) / density
    uy = np.tensordot(cy, distribution, axes=(0, 0)) / density
    assert ux == pytest.approx(uy, abs=1.0e-14)
