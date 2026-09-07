import numpy as np
import pytest
from resume_yawed_crosswind_trial10 import join_segments


def segment(start):
    f=np.ones((3,6))*2+np.arange(start,start+3)[:,None]
    return dict(arc_length=np.arange(3,dtype=float),states=f.copy(),fluxes=f,
        sources=np.ones_like(f),cumulative_sources=f-f[0])


def test_cumulative_sources_and_shared_boundary_are_preserved():
    a,b=segment(0),segment(2)
    out=join_segments(a,b)
    assert out['arc_length']==pytest.approx(np.arange(5))
    assert out['fluxes'][:,0]==pytest.approx(np.arange(2,7))
    assert out['cumulative_sources'][:,0]==pytest.approx(np.arange(5))
    assert out['maximum_relative_balance_residual']==0.


def test_boundary_mismatch_is_rejected():
    with pytest.raises(ValueError,match='boundary mismatch'):
        join_segments(segment(0),segment(3))


def test_nonzero_segment_origin_is_rejected():
    a,b=segment(0),segment(2);b['arc_length']+=1
    with pytest.raises(ValueError,match='start at zero'):
        join_segments(a,b)
