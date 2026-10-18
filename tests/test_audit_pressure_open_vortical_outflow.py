from pathlib import Path
import runpy

import numpy as np
import pytest


TOOL = runpy.run_path("tools/audit_pressure_open_vortical_outflow.py")


def test_initial_vortex_has_declared_peak_and_pressure_depression():
    distribution = TOOL["_initial_distribution"](TOOL["SHORT_NX"])
    density, ux, uy, uz = TOOL["_state"](distribution)
    speed = np.sqrt(
        (ux - TOOL["CARRIER_SPEED"]) ** 2 + uy**2 + uz**2
    )
    assert float(np.max(speed)) == pytest.approx(
        TOOL["PEAK_SWIRL_SPEED"], rel=0.01,
    )
    center_y = int(np.floor(TOOL["VORTEX_CENTER_Y"]))
    center_x = int(TOOL["VORTEX_CENTER_X"])
    assert density[0, center_y, center_x] < 1.0
    assert float(np.max(np.abs(uz))) < 1.0e-14


def test_classifier_requires_every_frozen_component_gate():
    passing = TOOL["classify_vortical_outflow"](
        maximum_velocity_error=0.04,
        maximum_vorticity_error=0.04,
        maximum_density_error=0.04,
        residual_enstrophy_ratio=0.04,
        final_mean_density_error=5.0e-5,
    )
    assert passing["passed"] is True
    failing = TOOL["classify_vortical_outflow"](
        maximum_velocity_error=0.04,
        maximum_vorticity_error=0.051,
        maximum_density_error=0.04,
        residual_enstrophy_ratio=0.04,
        final_mean_density_error=5.0e-5,
    )
    assert failing["passed"] is False
    assert failing["component_passes"]["vorticity"] is False


def test_classifier_rejects_nonfinite_metrics():
    with pytest.raises(ValueError, match="finite"):
        TOOL["classify_vortical_outflow"](
            maximum_velocity_error=np.nan,
            maximum_vorticity_error=0.0,
            maximum_density_error=0.0,
            residual_enstrophy_ratio=0.0,
            final_mean_density_error=0.0,
        )


def test_short_run_and_plot_are_reproducible(tmp_path):
    record = TOOL["run_audit"](steps=2)
    assert record["configuration"]["time_steps"] == 2
    assert set(record["schemes"]) == {
        TOOL["PREDECESSOR_SCHEME"], TOOL["CANDIDATE_SCHEME"],
    }
    output = tmp_path / "vortex.png"
    TOOL["plot_audit"](record, output)
    assert Path(output).is_file()
