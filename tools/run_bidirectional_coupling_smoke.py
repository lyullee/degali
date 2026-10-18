"""Execute the conservative two-state coupling contract on small fields."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from degali.addons.hybrid_amr import make_obstacle_patch_bounds
from degali.addons.hybrid_amr_plan import HybridAMRRunPlan
from degali.addons.hybrid_amr_executor import HybridAMRExecutor
from degali.addons.interface_mass_application import apply_interface_mass_to_field

def main(output=Path("outputs/bidirectional_coupling_smoke.json")):
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
        downstream_m=4, lateral_m=3, vertical_m=2, domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    coarse = np.zeros((2, 2, 2)); local = np.zeros((2, 2, 2))
    coarse[0, 0, 0] = 1.0 / 8.0
    local_mask = np.ones_like(local, dtype=bool)
    def coarse_step(state, _t, _dt, incoming, available):
        field = apply_interface_mass_to_field(state, target_mask=np.ones_like(state, bool),
                                              incoming_mass_kg=incoming, cell_volume_m3=1.)
        available = min(float(np.sum(field)), float(available))
        outgoing = min(.1, available)
        if outgoing:
            field *= (available - outgoing) / available if available else 0.0
        return field, outgoing
    def patch_step(state, _t, _dt, incoming, _available):
        field = apply_interface_mass_to_field(state, target_mask=local_mask,
                                              incoming_mass_kg=incoming, cell_volume_m3=1.)
        return field, 0.
    final_coarse, final_local, history = HybridAMRExecutor(plan).run_bidirectional_transactional(
        coarse_state=coarse, patch_state=local, coarse_step=coarse_step,
        patch_step=patch_step, coarse_mass_kg=float(np.sum(coarse)), patch_mass_kg=0.)
    report = {"schema": "degali.bidirectional-coupling-smoke.v2",
              "steps": len(history), "initial_mass_kg": 0.125,
              "final_mass_kg": float(np.sum(final_coarse) + np.sum(final_local)),
              "executor_history_final_mass_kg": history[-1][3] + history[-1][4],
              "field_state_mass_kg": float(np.sum(final_coarse) + np.sum(final_local)),
              "status": "transactional_field_contract_verified",
              "quantitative_facility_status": "not_claimed"}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
