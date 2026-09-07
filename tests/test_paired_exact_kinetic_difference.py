import numpy as np
import pytest
pytest.importorskip('mpmath')
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.enriched_segments import encode_enriched
from degali.addons.exact_transverse_geometry import exact_geometry_field_view
from degali.addons.paired_exact_kinetic_difference import paired_exact_kinetic_difference
from degali.addons.enriched_transport import EnrichedModalTransport


@pytest.fixture
def exact_transport(transport,monkeypatch):
    p=transport.projection
    # Supply real atmospheric parameters to the inherited toy geometry fixture.
    for name,value in dict(ustar=.2,zr=.1,rml=0.,spread_floor=False).items():
        monkeypatch.setattr(p.section.jetplume,name,value,raising=False)
    view,parameters,_=exact_geometry_field_view(p,encode_enriched(p,transport.parameters))
    return EnrichedModalTransport(view,parameters,scalar_mixing=transport.mixing,
        thermal_species_ratio=1.,mechanical_work='reduced_buoyancy_work_immediate_shear_heat')


def test_zero_actual_kinetic_difference_is_exact(exact_transport):
    m=exact_transport
    out=paired_exact_kinetic_difference(m.projection,m.parameters,np.zeros(m.count),1e-6,order=4,angular_order=4)
    assert out['derivative']==0.


def test_actual_velocity_only_kinetic_difference_matches_derivative(exact_transport):
    m=exact_transport; rates=np.zeros(m.count); rates[4]=.03
    expected=0.
    for prep in m.quadrature.rule(m.parameters,order=8,angular_order=8):
        local=m.local(prep['a'],prep['b'])
        expected+=float(prep['weights']@local['bk']@rates/m.area)
    out=paired_exact_kinetic_difference(m.projection,m.parameters,rates,1e-6,order=8,angular_order=8)
    assert out['derivative']==pytest.approx(expected,rel=1e-9,abs=1e-6)
