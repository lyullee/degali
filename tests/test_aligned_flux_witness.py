import numpy as np
import pytest
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.aligned_flux_witness import stream_aligned_columns,affine_vectors,modes


def test_affine_flux_increment_preserves_normalized_stress_direction(transport):
    d=transport.local(np.array([.1,.3,.7]),np.array([.2,.4,.5]))
    n=transport.count
    rng=np.random.default_rng(730)
    fm,fp=rng.normal(size=(2,3,n+1))
    rates=rng.normal(size=n); response=rng.normal(size=(n,4))
    out=affine_vectors(transport,d,fm,fp,rates,response)
    stress=out['momentum_j']-d['u'][:,None,None]*out['mass_j']
    cross=d['a'][:,None]*stress[:,1]-d['b'][:,None]*stress[:,0]
    assert cross==pytest.approx(np.zeros_like(cross),abs=1e-12)
    delta=rng.normal(size=4)
    actual=affine_vectors(transport,d,fm,fp,rates+response@delta,np.zeros((n,4)))
    v=modes(transport,d)
    for key,mode in (('mass','mass_modes'),('momentum','momentum_modes')):
        expected=actual[key]+np.einsum('nik,k->ni',v[mode],delta)
        assert out[key]+np.einsum('nik,k->ni',out[key+'_j'],delta)==pytest.approx(expected,abs=1e-11)


def test_aligned_columns_streaming_and_global_conservation(transport):
    a,b=[stream_aligned_columns(transport,order=8,angular_order=8,batch_size=batch) for batch in (2,8)]
    assert a['columns']==pytest.approx(b['columns'],rel=1e-10,abs=1e-7)
    assert np.all(a['columns'][:5]==0.)
    for key in ('global_heat_correction','mass_face_correction','momentum_face_correction'):
        assert a[key]==pytest.approx(np.zeros(4),abs=1e-7)
