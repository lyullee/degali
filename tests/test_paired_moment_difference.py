import numpy as np
import pytest
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.paired_moment_difference import paired_moment_difference
from degali.addons.enriched_transport import shifted_advective_moments


def test_zero_direction_is_exactly_zero_on_common_grid(transport):
    out=paired_moment_difference(transport.projection,transport.parameters,np.zeros(transport.count),1e-5,order=4,angular_order=4)
    assert np.array_equal(out['derivative'],np.zeros(5+transport.projection.count))


def test_paired_actual_difference_matches_separate_large_step_values(transport):
    rates=np.zeros(transport.count); rates[1],rates[4]=.03,-.02
    step=1e-3
    out=paired_moment_difference(transport.projection,transport.parameters,rates,step,order=4,angular_order=4)
    plus=shifted_advective_moments(transport.projection,transport.parameters,rates,step,order=4,angular_order=4)
    minus=shifted_advective_moments(transport.projection,transport.parameters,rates,-step,order=4,angular_order=4)
    expected=(plus-minus)/(2*step)
    assert max(abs(out['derivative']-expected)/np.maximum(abs(expected),1.))<1e-6
    assert out['common_union_phase_grid'] and out['fields_are_binary64']
