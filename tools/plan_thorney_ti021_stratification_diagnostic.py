"""Freeze the observation-blind TI021 stratification diagnostic envelope."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-stratification-diagnostic-plan.v1"
POLICY_SCHEMA = "degali.thorney-ti021-turbulent-mixing-policy.v1"


def build_plan(mixing_policy_path: Path) -> dict[str, object]:
    path = mixing_policy_path.resolve()
    policy = json.loads(path.read_text(encoding="utf-8"))
    if policy.get("schema") != POLICY_SCHEMA or not policy.get("all_gates_pass"):
        raise ValueError("accepted TI021 turbulent-mixing policy is required")
    accessed = policy["evidence_provenance"].get(
        "concentration_observation_fields_accessed"
    )
    if accessed != []:
        raise ValueError("TI021 turbulent-mixing policy is not observation-blind")
    friction_velocity = float(policy["closure"]["friction_velocity_m_s"])
    wall_height = 5.0
    primary = friction_velocity / wall_height
    multipliers = (0.5, 1.0, 2.0)
    members = [{
        "member_id": f"epsilon-u-star-over-H-x{multiplier:g}",
        "multiplier": multiplier,
        "shear_regularization_s_inv": primary * multiplier,
        "role": "primary" if multiplier == 1.0 else "mandatory_sensitivity",
    } for multiplier in multipliers]
    gates = {
        "mixing_policy_gates_pass": True,
        "selection_uses_no_concentration_dose_or_timing_observation": True,
        "regularization_scale_uses_only_declared_u_star_and_wall_height": True,
        "factor_two_low_primary_high_envelope_is_retained": True,
        "all_members_are_diagnostic_only": True,
        "closure_promotion_is_prohibited": True,
    }
    return {
        "schema": SCHEMA,
        "status": "accepted_observation_blind_stratification_diagnostic_plan",
        "scope": (
            "regularisation sensitivity for model-only gradient Richardson "
            "diagnosis; no flux damping and no closure promotion"
        ),
        "provenance": {
            "mixing_policy_path": str(path),
            "mixing_policy_sha256": _sha256(path),
            "observation_value_fields_accessed": [],
        },
        "scale_basis": {
            "formula": "epsilon_s = multiplier * u_star / H_wall",
            "friction_velocity_m_s": friction_velocity,
            "wall_height_m": wall_height,
            "primary_shear_regularization_s_inv": primary,
            "rationale": (
                "event momentum scale divided by report-priority obstacle "
                "height; factor-two bracket retained before diagnostic output"
            ),
        },
        "members": members,
        "execution_policy": {
            "required_members": [item["member_id"] for item in members],
            "selection_after_execution": "prohibited",
            "damping_applied": False,
            "closure_promotion_allowed": False,
            "required_regions": [
                "windward_inner_face_shell",
                "roof_shear_band",
                "side_tip_bypass",
                "exterior_lee_wake_shell",
            ],
            "required_receptors": "all coordinate-only inner receptor cells",
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mixing-policy", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    record = build_plan(args.mixing_policy)
    destination = args.output_json.resolve()
    destination.parent.mkdir(parents=True, exist_ok=False)
    destination.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
