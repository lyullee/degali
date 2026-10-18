import runpy
from pathlib import Path

import numpy as np
import pytest


TOOL = runpy.run_path(
    "tools/audit_pressure_open_boundary_3d_corner_reflection.py"
)


def _synthetic_traces(reflected_ratio, active_axis_count):
    size = TOOL["STEPS"] + 1
    reference = np.zeros(size)
    short = np.zeros(size)
    preliminary = TOOL["classify_corner_reflection"](
        short, reference + 1.0e-4,
        active_axis_count=active_axis_count,
    )
    direct_start, direct_stop = preliminary["incident_window_steps"]
    reference[(direct_start + direct_stop) // 2] = 1.0e-4
    short[:] = reference
    preliminary = TOOL["classify_corner_reflection"](
        short, reference, active_axis_count=active_axis_count,
    )
    start, stop = preliminary["reflection_window_steps"]
    short[(start + stop) // 2] += reflected_ratio * 1.0e-4
    return short, reference


@pytest.mark.parametrize("active_axis_count", [2, 3])
def test_3d_corner_classifier_applies_frozen_limit(active_axis_count):
    short, reference = _synthetic_traces(0.04, active_axis_count)
    result = TOOL["classify_corner_reflection"](
        short, reference, active_axis_count=active_axis_count,
    )
    assert result["reflection_amplitude_ratio"] == pytest.approx(0.04)
    assert result["reflection_amplitude_ratio_limit"] == 0.05
    assert result["passed"] is True


@pytest.mark.parametrize(
    "orientation", [
        "x_y_scaled_control", "x_z_edge", "y_z_edge", "x_y_z_corner",
    ],
)
def test_packet_velocity_is_diagonal_only_on_declared_axes(orientation):
    active_axes = TOOL["ORIENTATIONS"][orientation]
    distribution = TOOL["_initial_distribution"](
        TOOL["SHORT_SIZE"], active_axes,
    )
    density = np.sum(distribution, axis=0)
    velocities = {}
    for component, axis in enumerate(("x", "y", "z")):
        lattice_component = np.asarray([
            direction[component] for direction in TOOL["_C3"]
        ])
        velocities[axis] = (
            np.tensordot(lattice_component, distribution, axes=(0, 0))
            / density
        )
    active_reference = velocities[active_axes[0]]
    for axis in ("x", "y", "z"):
        if axis in active_axes:
            assert velocities[axis] == pytest.approx(
                active_reference, abs=1.0e-14,
            )
        else:
            assert velocities[axis] == pytest.approx(0.0, abs=1.0e-14)


def test_orientation_shapes_keep_only_active_axes_large():
    assert TOOL["_shape"](48, ("x", "z")) == (48, 4, 48)
    assert TOOL["_shape"](48, ("y", "z")) == (48, 48, 4)
    assert TOOL["_shape"](48, ("x", "y", "z")) == (48, 48, 48)


def test_packet_wavelength_is_an_explicit_validated_input():
    default = TOOL["_initial_distribution"](48, ("x", "y"))
    alternate = TOOL["_initial_distribution"](48, ("x", "y"), 10.0)
    assert not np.allclose(default, alternate)
    with pytest.raises(ValueError, match="pulse_wavelength_cells"):
        TOOL["_initial_distribution"](48, ("x", "y"), 0.0)


def test_plot_accepts_every_declared_orientation(tmp_path):
    trace = [0.0] * (TOOL["STEPS"] + 1)
    scheme = {
        "reflection_amplitude_ratio": 0.01,
        "reflection_window_steps": [10, 20],
        "probe_trace": trace,
    }
    record = {"orientations": {}}
    for name in TOOL["ORIENTATIONS"]:
        record["orientations"][name] = {
            "reference_probe_trace": trace,
            "schemes": {
                "characteristic_nonreflecting": scheme,
                "characteristic_multidimensional": scheme,
            },
        }
    output = tmp_path / "corner.png"
    TOOL["plot_audit"](record, output)
    assert Path(output).is_file()
