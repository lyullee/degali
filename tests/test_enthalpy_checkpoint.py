"""The enthalpy experiment cannot reuse the old energy-only convergence gate."""

import json
from pathlib import Path
import runpy
import sys

import pytest


@pytest.mark.parametrize("quadrature", [None, "exact_symmetric_nodes_v1"])
def test_enthalpy_field_rejects_energy_only_checkpoint(tmp_path, monkeypatch, capsys, quadrature):
    checkpoint = {
        "interfaces_accepted": True, "failures": {},
        "selected_trials": [10, 11, 12, 22, 23, 24, 25],
        "downstream_thermodynamic_profile": "enthalpy",
        "phase_ambient_closure": "consistent_explicit_ideal_v1",
        "mode": "full", "droplet_equilibrium_bound": True,
        "hydrogen_spin_isomer": "normal",
    }
    if quadrature is not None:
        checkpoint["energy_quadrature"] = quadrature
    path = tmp_path / "checkpoint.json"
    path.write_text(json.dumps(checkpoint), encoding="utf-8")
    script = Path(__file__).resolve().parents[1] / "tools/run_preslhy_source_ablation.py"
    monkeypatch.setattr(sys, "argv", [
        str(script), "--phase", "field", "--mode", "full", "--droplet-equilibrium-bound",
        "--consistent-phase-ambient", "--downstream-profile", "enthalpy",
        "--interface-checkpoint", str(path),
    ])
    module = runpy.run_path(str(script))
    with pytest.raises(SystemExit) as error:
        module["main"]()
    assert error.value.code == 2
    assert "checkpoint does not match" in capsys.readouterr().err
