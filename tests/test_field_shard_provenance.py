"""Reject merging field results with different source physics or isomers."""

import copy
import json
import runpy
import sys
from pathlib import Path

import pytest


def _shard(trial):
    metric = {
        "arcs": 1, "MG": 1.1, "VG": 1.2, "FAC2": 1.0,
        "geometry": {"width_ratio": 1.05, "centre_mean_error": 0.02, "centre_mae": 0.02},
        "pairs": [{"trial": trial, "observed": 1.0, "predicted": 1.1}],
        "vertical_profiles": [{"trial": trial}],
    }
    return {
        "mode": "full", "droplet_equilibrium_bound": True,
        "hydrogen_spin_isomer": "normal", "source_energy_ledger": "moving_pipe_phase_enthalpy_v2",
        "selected_trials": [trial], "failures": {}, "vertical_rows": 1,
        "interfaces": {str(trial): {"accepted": True}},
        "baseline": copy.deepcopy(metric), "candidate": copy.deepcopy(metric),
    }


def _merge(tmp_path, monkeypatch, first, second):
    paths = [tmp_path / "a.json", tmp_path / "b.json", tmp_path / "out.json"]
    for path, data in zip(paths, (first, second)):
        path.write_text(json.dumps(data), encoding="utf-8")
    script = Path(__file__).resolve().parents[1] / "tools/merge_preslhy_field_shards.py"
    monkeypatch.setattr(sys, "argv", [str(script), *map(str, paths)])
    module = runpy.run_path(str(script))
    module["main"]()
    return json.loads(paths[2].read_text(encoding="utf-8"))


@pytest.mark.parametrize("key,value", [
    ("hydrogen_spin_isomer", "para"), ("source_energy_ledger", "legacy"),
    ("phase_ambient_closure", "consistent_explicit_ideal_v1"),
    ("downstream_thermodynamic_profile", "enthalpy"),
    ("energy_quadrature", "polar_square_all_flux_v1"),
])
def test_shard_merge_rejects_changed_physics(tmp_path, monkeypatch, key, value):
    first, second = _shard(11), _shard(12)
    second[key] = value
    with pytest.raises(ValueError, match=key):
        _merge(tmp_path, monkeypatch, first, second)


def test_shard_merge_preserves_pointwise_results_and_physics(tmp_path, monkeypatch):
    result = _merge(tmp_path, monkeypatch, _shard(11), _shard(12))
    assert result["source_energy_ledger"] == "moving_pipe_phase_enthalpy_v2"
    assert result["hydrogen_spin_isomer"] == "normal"
    assert result["phase_ambient_closure"] == "legacy_air_eos"
    assert result["energy_quadrature"] == "unspecified"
    assert len(result["aggregation"]["source_sha256"]) == 2
    for kind in ("baseline", "candidate"):
        assert [row["trial"] for row in result[kind]["pairs"]] == [11, 12]
        assert len(result[kind]["vertical_profiles"]) == 2
        assert result[kind]["MG"] == pytest.approx(1.1)
