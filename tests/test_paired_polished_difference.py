import numpy as np
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.paired_moment_difference import paired_moment_difference
from degali.addons.paired_polished_difference import paired_polished_moment_difference


def test_polished_pair_uses_same_fields_and_explicit_root_statistics(transport):
    rates=np.zeros(transport.count); rates[1],rates[4]=.03,-.02
    args=(transport.projection,transport.parameters,rates,1e-3)
    old=paired_moment_difference(*args,order=4,angular_order=4)
    new=paired_polished_moment_difference(*args,order=4,angular_order=4)
    assert max(abs(new['derivative']-old['derivative'])/np.maximum(abs(old['derivative']),1.))<1e-6
    assert new['roots_polished_on_same_table']
    for stat in new['root_statistics']:
        assert stat['calls']>0 and stat['points']>0
        assert stat['maximum_relative_residual_after']<=stat['maximum_relative_residual_before']
