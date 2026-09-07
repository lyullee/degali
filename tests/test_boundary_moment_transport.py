"""Algebra/field invariants, not scientific validation of row replacement."""

import numpy as np
import pytest
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.boundary_moment_transport import lift_same_fields,raw_modal_system,BoundaryMomentTransport


def test_field_lift_preserves_actual_values_and_derivatives(transport):
    p,par=transport.projection,transport.parameters
    view,new=lift_same_fields(p,par,p.basis.degree+2)
    rng=np.random.default_rng(410)
    a,b=rng.uniform(0.,1.,(2,81))
    oldprep,newprep=p.prepare(a,b),view.prepare(a,b)
    for key in ('psi','da','db'):
        for i in range(2):
            assert oldprep[key]@np.split(par,2)[i] == pytest.approx(newprep[key]@np.split(new,2)[i],abs=2e-14)
    old,newfields=p.fields(par,oldprep),view.fields(new,newprep)
    for key in ('rho','c','h','y','temperature'):
        assert newfields[key] == pytest.approx(old[key],rel=1e-12,abs=1e-10)
    assert view.basis is not p.basis


def test_raw_test_transformation_is_not_whitened_index_truncation(transport):
    p=transport.projection
    rng=np.random.default_rng(112)
    matrix=rng.normal(size=(5+p.count,7+p.count))
    right=rng.normal(size=5+p.count)
    raw,rhs,ordered=raw_modal_system(p,matrix,right)
    size=p.basis.size
    assert p.basis.whitener.T@raw[:size] == pytest.approx(matrix[5:5+size],abs=1e-13)
    assert p.basis.whitener.T@rhs[:size] == pytest.approx(right[5:5+size],abs=1e-13)
    assert [sum(p.basis.pairs[j]) for j in ordered] == sorted(sum(v) for v in p.basis.pairs)


def test_degenerate_boundary_system_is_rejected_without_pseudoinverse(transport):
    p,par=transport.projection,transport.parameters
    # Zero scalar enrichment can give dependent scalar boundary constraints.
    p,par=lift_same_fields(p,par,4)
    model=BoundaryMomentTransport(p,par,scalar_mixing=transport.mixing,
        thermal_species_ratio=1.,mechanical_work='reduced_buoyancy_work_immediate_shear_heat',
        order=8,angular_order=24,split_angles=False)
    with pytest.raises(ValueError,match='rank deficient'):
        model.assemble()
    diag=model.linear_diagnostics
    assert diag['singular_values'][-1] <= 1e-12*diag['singular_values'][0]


def test_boundary_replacement_retains_global_rows_and_reports_bad_solve(transport):
    p,par=transport.projection,transport.parameters.copy()
    par += np.random.default_rng(512).normal(0.,1e-4,len(par))
    model=BoundaryMomentTransport(p,par,scalar_mixing=transport.mixing,
        thermal_species_ratio=1.,mechanical_work='reduced_buoyancy_work_immediate_shear_heat',
        order=8,angular_order=24,split_angles=False)
    out=model.assemble()
    assert out['matrix'].shape==(5+p.count,7+p.count)
    assert len(out['omitted_raw_indices'])==8
    assert not out['all_original_weak_equations_enforced']
    assert np.array_equal(out['matrix'][:5],out['original_matrix'][:5])
    assert np.array_equal(out['right'][:5],out['original_right'][:5])
    independent_error=max(abs(out['matrix']@out['rates']-out['right'])/np.maximum(abs(out['right']),1.))
    assert out['linear_scaled_error']==pytest.approx(independent_error)
    # These nearly redundant synthetic constraints produce severe numerical
    # cancellation. A successful matrix solve must not imply acceptance.
    assert independent_error>1e-8
    assert not out['local_numerics_passed']
    assert np.all(np.isfinite(out['omitted_scaled_residual']))
    assert not out['adopted']
