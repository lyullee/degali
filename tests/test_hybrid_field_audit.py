import numpy as np
import pytest
from degali.addons.hybrid_field_audit import (
    aggregate_substeps_to_coarse, audit_hybrid_mass_histories, audit_solver_results,
    audit_solver_results_with_patch_substeps, audit_partitioned_solver_results,
    infer_net_internal_flux_rate, threshold_region_metrics)

def test_threshold_region_metrics_reports_volume_and_centroid():
    field = np.zeros((2, 2, 3), dtype=float)
    field[0, 0, 1] = 4.0
    field[1, 1, 2] = 8.0
    result = threshold_region_metrics(
        field, threshold_vol_percent=4.0, cell_volume_m3=0.5,
        coordinates_m=(np.array([0.0, 1.0]), np.array([0.0, 2.0]),
                       np.array([0.0, 3.0, 6.0])),
    )
    assert result["cell_count"] == 2
    assert result["volume_m3"] == 1.0
    assert result["maximum_vol_percent"] == 8.0
    assert result["centroid_m"] == [0.5, 1.0, 4.5]

def test_audit_accepts_common_clock_and_reports_drift():
    t = np.arange(4, dtype=float)
    audit = audit_hybrid_mass_histories(
        coarse_time_s=t, coarse_mass_kg=[3, 2.5, 2, 1.5],
        patch_time_s=t, patch_mass_kg=[1, 1.5, 2, 2.5])
    assert audit.quantitative_ready
    assert audit.maximum_total_mass_drift_kg < 1e-12

def test_audit_rejects_misaligned_histories():
    with pytest.raises(ValueError, match="common time axis"):
        audit_hybrid_mass_histories(
            coarse_time_s=[0, 1], coarse_mass_kg=[1, 1],
            patch_time_s=[0, 1.1], patch_mass_kg=[1, 1])

def test_audit_extracts_public_solver_histories():
    class Coarse:
        time_s = np.array([0., 1.])
        domain_mass_kg = np.array([2., 1.])
    class Patch:
        time_s = np.array([0., 1.])
        domain_mass_kg = np.array([0., 1.])
    audit = audit_solver_results(Coarse(), Patch())
    assert audit.quantitative_ready

def test_substep_aggregation_uses_exact_endpoints_without_interpolation():
    fine_t = np.arange(0., 2.01, .25)
    fine_m = fine_t * 2
    coarse_t = np.array([0., 1., 2.])
    assert np.allclose(aggregate_substeps_to_coarse(
        fine_time_s=fine_t, fine_mass_kg=fine_m, coarse_time_s=coarse_t), [0, 2, 4])

def test_substep_aggregation_rejects_missing_endpoint():
    with pytest.raises(ValueError, match="exact fine-grid endpoint"):
        aggregate_substeps_to_coarse(
            fine_time_s=[0., .3, .6], fine_mass_kg=[1., 1., 1.], coarse_time_s=[0., 1.])

def test_solver_result_audit_aggregates_patch_substeps():
    class Coarse:
        time_s = np.array([0., 1., 2.])
        domain_mass_kg = np.array([3., 2.5, 2.])
    class Patch:
        time_s = np.arange(0., 2.01, .25)
        domain_mass_kg = None
    Patch.domain_mass_kg = Patch.time_s / 2
    audit = audit_solver_results_with_patch_substeps(Coarse(), Patch())
    assert audit.quantitative_ready

def test_partitioned_audit_rejects_double_counting_weights():
    class R:
        time_s = np.array([0., 1.])
        domain_mass_kg = np.array([2., 1.])
    with pytest.raises(ValueError, match="sum to one"):
        audit_partitioned_solver_results(R(), R(), coarse_weight=.8, patch_weight=.8)

def test_partitioned_audit_preserves_weighted_mass():
    class C:
        time_s = np.array([0., 1.])
        domain_mass_kg = np.array([2., 1.])
    class P:
        time_s = np.array([0., .5, 1.])
        domain_mass_kg = np.array([2., 2.5, 3.])
    audit = audit_partitioned_solver_results(C(), P(), coarse_weight=.5, patch_weight=.5)
    assert audit.quantitative_ready

def test_infer_net_internal_flux_from_solver_balance():
    class R:
        time_s = np.array([0., 1., 2.])
        domain_mass_kg = np.array([1., 1.5, 2.])
        injected_mass_kg = np.array([0., .5, 1.])
        outflow_mass_kg = np.array([0., 0., 0.])
    assert np.allclose(infer_net_internal_flux_rate(R()), [0., 0.])
