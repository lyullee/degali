"""Saved-state diagnostics may not silently reinterpret another experiment."""

import json
from pathlib import Path
import runpy

import pytest

pytest.importorskip("CoolProp")
ROOT = Path(__file__).resolve().parents[1]
AUDIT = runpy.run_path(str(ROOT / "tools/run_preslhy_ambient_profile_audit.py"))


@pytest.mark.parametrize("field", [
    {"downstream_thermodynamic_profile": "unknown"},
    {"mode": "flow_only"},
    {"mode": "full", "droplet_equilibrium_bound": True, "hydrogen_spin_isomer": "para"},
    {"mode": "full", "droplet_equilibrium_bound": True,
     "source_energy_ledger": "moving_pipe_phase_enthalpy_v2", "phase_ambient_closure": "unknown"},
])
def test_replay_rejects_incompatible_physics(field):
    with pytest.raises(ValueError):
        AUDIT["replay"](field, {})


def test_frozen_control_replay_reproduces_temperature_flux_and_sensor_arcs():
    field_path = ROOT / "reference/preslhy/source_phase_energy_ledger_field_complete_2026-09-05.json"
    if not field_path.exists():
        pytest.skip("frozen field artifact not distributed")
    field = json.loads(field_path.read_text(encoding="utf-8"))
    trial = next(t for t in json.loads((ROOT / "reference/preslhy/e35_reduced.json").read_text())["trials"]
                 if t["trial"] == 23)
    trajectory, check = AUDIT["replay"](field, trial)
    assert trajectory is not None
    assert check["temperature_max_abs_K"] < 1e-6
    assert check["sampled_flux_max_scaled_error"] < 1e-9
    assert check["arc_max_relative_error"] < 1e-9
