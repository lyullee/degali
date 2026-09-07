import numpy as np
import pytest
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.polished_phase_inverse import PolishedPhaseMassEnthalpyInverter


def test_polishing_improves_residual_without_changing_oracle_or_branch(transport):
    p=transport.projection
    x=np.linspace(.02,1.,37)
    d=p.fields(transport.parameters,p.prepare(x,x*.7))
    original=p.section.phase_inverse
    refined=PolishedPhaseMassEnthalpyInverter(p.section.thermodynamics)
    rho,y,t=refined.state(d['c'],d['h'])
    before=abs(original.enthalpy_and_slope(d['rho'],d['c'])[0]-d['h'])
    after=abs(original.enthalpy_and_slope(rho,d['c'])[0]-d['h'])
    assert np.all(after<=before)
    assert np.any(after<before)
    assert rho==pytest.approx(d['rho'],rel=1e-9,abs=1e-12)
    assert np.all((rho>0.)&(y>0.)&(y<1.))
    assert np.array_equal(refined.h_values,original.h_values)
    assert refined.statistics['improved_points']>0


def test_polishing_does_not_accept_an_unphysical_request(transport):
    refined=PolishedPhaseMassEnthalpyInverter(transport.projection.section.thermodynamics)
    with pytest.raises(ValueError):
        refined.state(np.array([-1.]),np.array([0.]))
