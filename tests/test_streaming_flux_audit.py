import numpy as np
import pytest
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.solenoidal_transport import SolenoidalTransport
from degali.addons.streaming_flux_audit import stream_retained_operator


def test_streamed_rows_match_dense_retained_rows_without_hard_boundary_solve(transport):
    dense=SolenoidalTransport(transport.projection,transport.parameters,scalar_mixing=transport.mixing,
        thermal_species_ratio=1.,mechanical_work='reduced_buoyancy_work_immediate_shear_heat',
        flux_modes=4,order=4,angular_order=12)
    try:
        full=dense.assemble()
    except ValueError as exc:
        assert 'rank' in str(exc)
        full=dense.linear_diagnostics
    streamed=stream_retained_operator(transport,order=4,angular_order=12,split_angles=False,batch_size=3)
    n=transport.count-2
    for key in ('matrix','right'):
        assert streamed[key]==pytest.approx(full[key][:n],rel=1e-10,abs=1e-8)
    assert streamed['angles']==12


def test_fixed_witness_diagnostics_are_batch_independent(transport):
    rates=np.random.default_rng(789).normal(0.,.01,transport.count)
    witness=dict(rates=rates,mass_modes=np.array([.01,-.02,.03,.004]),momentum_modes=np.array([.02,.01,.03,.002]))
    a,b=[stream_retained_operator(transport,order=4,angular_order=12,split_angles=False,batch_size=n,witness=witness) for n in (2,6)]
    for key in ('matrix','right','advective_jacobian','retained_scaled_residual'):
        assert a[key]==pytest.approx(b[key],rel=1e-10,abs=1e-8)
    for key in ('weak_heat_scaled_error','mass_boundary_error','maximum_retained_error','minimum_shear_production'):
        assert a[key]==pytest.approx(b[key],rel=1e-10,abs=1e-8)
    assert a['peak_batch_nodes']<b['peak_batch_nodes']
