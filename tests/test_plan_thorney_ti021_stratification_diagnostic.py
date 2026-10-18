import json

import pytest

from tools.plan_thorney_ti021_stratification_diagnostic import build_plan


def test_stratification_plan_freezes_factor_two_envelope_without_observations(
    tmp_path,
):
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({
        "schema": "degali.thorney-ti021-turbulent-mixing-policy.v1",
        "all_gates_pass": True,
        "evidence_provenance": {
            "concentration_observation_fields_accessed": [],
        },
        "closure": {"friction_velocity_m_s": 0.205},
    }), encoding="utf-8")
    record = build_plan(policy)
    assert record["all_gates_pass"]
    assert record["provenance"]["observation_value_fields_accessed"] == []
    assert record["scale_basis"]["primary_shear_regularization_s_inv"] == (
        pytest.approx(0.041)
    )
    assert [item["shear_regularization_s_inv"] for item in record["members"]] == (
        pytest.approx([0.0205, 0.041, 0.082])
    )
    assert record["execution_policy"]["selection_after_execution"] == "prohibited"
    assert record["execution_policy"]["closure_promotion_allowed"] is False


def test_stratification_plan_rejects_observation_access(tmp_path):
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({
        "schema": "degali.thorney-ti021-turbulent-mixing-policy.v1",
        "all_gates_pass": True,
        "evidence_provenance": {
            "concentration_observation_fields_accessed": ["dose"],
        },
        "closure": {"friction_velocity_m_s": 0.205},
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="not observation-blind"):
        build_plan(policy)
