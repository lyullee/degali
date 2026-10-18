import json
from pathlib import Path

import numpy as np
import pytest

from tools.audit_thorney_ti021_stratification_diagnostic import (
    audit_stratification,
)
from tools.audit_thorney_ti021_stratification_envelope import (
    audit_stratification_envelope,
)
from tools.run_thorney_ti021_inner_transient import _sha256


def _write_case(root: Path, *, complete: bool = True):
    run = root / "run"
    segment_dir = run / "checkpoint-segments"
    segment_dir.mkdir(parents=True)
    x = np.arange(0.5, 120.0, 5.0)
    y = np.arange(2.5, 200.0, 5.0)
    z = np.arange(2.5, 20.0, 5.0)
    shape = (z.size, y.size, x.size)
    solid = np.zeros(shape, dtype=bool)
    density = np.broadcast_to(
        np.linspace(1.4, 1.1, z.size)[:, None, None], shape,
    ).copy()
    u = np.broadcast_to(np.linspace(1.0, 2.0, z.size)[:, None, None], shape).copy()
    v = np.zeros(shape)
    segment = segment_dir / "segment-0001.npz"
    np.savez_compressed(
        segment,
        x_m=x,
        y_m=y,
        z_m=z,
        solid_mask=solid,
        final_mixture_density_kg_m3=density,
        final_velocity_u_m_s=u,
        final_velocity_v_m_s=v,
    )
    duration = 1.0 if complete else 2.0
    journal = {
        "schema": "degali.thorney-ti021-inner-checkpoint.v1",
        "configuration": {
            "duration_s": duration,
            "cell_size_m": 5.0,
        },
        "segments": [{
            "segment_index": 1,
            "start_time_s": 0.0,
            "stop_time_s": 1.0,
            "segment_file": "checkpoint-segments/segment-0001.npz",
            "segment_sha256": _sha256(segment),
        }],
    }
    (run / "ti021-inner-checkpoint.json").write_text(
        json.dumps(journal), encoding="utf-8",
    )
    run_record = {
        "all_gates_pass": True,
        "configuration": {
            "observation_values_used_to_configure_run": [],
            "receptors": [{
                "label": "R1",
                "solver_xyz_m": [75.0, 100.0, 2.5],
            }],
        },
    }
    (run / "ti021-inner-transient.json").write_text(
        json.dumps(run_record), encoding="utf-8",
    )
    policy = root / "policy.json"
    policy.write_text(json.dumps({
        "schema": "degali.thorney-ti021-turbulent-mixing-policy.v1",
        "all_gates_pass": True,
        "evidence_provenance": {
            "concentration_observation_fields_accessed": [],
        },
        "closure": {
            "friction_velocity_m_s": 0.205,
            "roughness_length_m": 0.005,
            "von_karman_constant": 0.4,
        },
        "frozen_policy": {"primary_turbulent_schmidt_number": 0.7},
    }), encoding="utf-8")
    plan = root / "plan.json"
    plan.write_text(json.dumps({
        "schema": "degali.thorney-ti021-stratification-diagnostic-plan.v1",
        "all_gates_pass": True,
        "provenance": {
            "observation_value_fields_accessed": [],
            "mixing_policy_sha256": _sha256(policy),
        },
        "members": [
            {"member_id": "low", "role": "mandatory_sensitivity",
             "shear_regularization_s_inv": 0.025},
            {"member_id": "primary", "role": "primary",
             "shear_regularization_s_inv": 0.05},
            {"member_id": "high", "role": "mandatory_sensitivity",
             "shear_regularization_s_inv": 0.1},
        ],
        "execution_policy": {"selection_after_execution": "prohibited"},
    }), encoding="utf-8")
    return run, policy, plan


def test_stratification_audit_is_observation_blind_and_non_promoting(tmp_path):
    run, policy, plan = _write_case(tmp_path)
    record, paths = audit_stratification(
        run_directory=run,
        mixing_policy_path=policy,
        diagnostic_plan_path=plan,
        shear_regularization_s_inv=0.05,
        output_directory=tmp_path / "output",
    )
    assert record["all_gates_pass"]
    assert record["provenance"]["observation_value_fields_accessed"] == []
    assert record["configuration"]["damping_applied"] is False
    assert record["promotion_policy"]["closure_promotion_allowed"] is False
    assert len(record["segments"]) == 1
    assert len(record["segments"][0]["receptors"]) == 1
    assert set(record["segments"][0]["regions"]) == {
        "windward_inner_face_shell",
        "roof_shear_band",
        "side_tip_bypass",
        "exterior_lee_wake_shell",
    }
    assert paths["json"].is_file()
    assert paths["figure"].is_file()


def test_stratification_audit_rejects_incomplete_checkpoint_journal(tmp_path):
    run, policy, plan = _write_case(tmp_path, complete=False)
    with pytest.raises(ValueError, match="not complete"):
        audit_stratification(
            run_directory=run,
            mixing_policy_path=policy,
            diagnostic_plan_path=plan,
            shear_regularization_s_inv=0.05,
            output_directory=tmp_path / "output",
        )


def test_stratification_audit_rejects_segment_hash_drift(tmp_path):
    run, policy, plan = _write_case(tmp_path)
    journal_path = run / "ti021-inner-checkpoint.json"
    journal = json.loads(journal_path.read_text(encoding="utf-8"))
    journal["segments"][0]["segment_sha256"] = "0" * 64
    journal_path.write_text(json.dumps(journal), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        audit_stratification(
            run_directory=run,
            mixing_policy_path=policy,
            diagnostic_plan_path=plan,
            shear_regularization_s_inv=0.05,
            output_directory=tmp_path / "output",
        )


def test_stratification_audit_rejects_unplanned_regularization(tmp_path):
    run, policy, plan = _write_case(tmp_path)
    with pytest.raises(ValueError, match="not a frozen plan member"):
        audit_stratification(
            run_directory=run,
            mixing_policy_path=policy,
            diagnostic_plan_path=plan,
            shear_regularization_s_inv=0.075,
            output_directory=tmp_path / "output",
        )


def test_stratification_envelope_executes_every_frozen_member(tmp_path):
    run, policy, plan = _write_case(tmp_path)
    record, paths = audit_stratification_envelope(
        run_directory=run,
        mixing_policy_path=policy,
        diagnostic_plan_path=plan,
        output_directory=tmp_path / "envelope",
    )
    assert record["all_gates_pass"]
    assert len(record["member_artifacts"]) == 3
    assert record["promotion_policy"]["member_selection_allowed"] is False
    assert record["provenance"]["observation_value_fields_accessed"] == []
    assert paths["json"].is_file()
    assert paths["figure"].is_file()
