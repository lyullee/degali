from degali.addons.hybrid_amr import make_obstacle_patch_bounds
from degali.addons.hybrid_amr_plan import HybridAMRRunPlan


def test_plan_records_local_patch_and_clock_ratio():
    patch = make_obstacle_patch_bounds(
        (10, 20, 5, 15, 0, 8), upstream_m=4, downstream_m=4,
        lateral_m=3, vertical_m=2, domain_m=(30, 25, 20),
    )
    plan = HybridAMRRunPlan(
        domain_m=(30, 25, 20), coarse_cell_size_m=5, patch_cell_size_m=1,
        patch=patch, duration_s=600, coarse_time_step_s=1, patch_time_step_s=.25,
    )
    assert plan.patch_volume_fraction < 1
    assert plan.substeps_per_coarse_step == 4
    coarse_weight, patch_weight = plan.partition_weights
    assert coarse_weight > patch_weight
    assert coarse_weight + patch_weight == 1.0
