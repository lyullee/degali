import numpy as np
import pytest

from audit_downstream_cold_envelope import reflected_peak, scalar_shape, envelope


@pytest.mark.parametrize('ratio',[0.,.1,.9,1.,1.+1e-10,1.01,1.2,2.,5.,100.])
def test_reflected_mode_against_dense_independent_profile(ratio):
    sigma=.23
    centre=ratio*sigma
    z,q=reflected_peak(centre,sigma)
    coords=np.linspace(0,centre+6*sigma,20001)
    profile=np.exp(-.5*((coords-centre)/sigma)**2)+np.exp(-.5*((coords+centre)/sigma)**2)
    assert q>=max(profile)-2e-14
    assert q-max(profile)<1e-5
    assert 0<=z<=centre
    if ratio<=1: assert z==0
    if ratio>1.01:
        assert z==pytest.approx(centre*np.tanh(z*centre/sigma**2),abs=1e-12)


@pytest.mark.parametrize('centre,sigma',[(-1.,1.),(1.,0.),(1.,-1.),(float('nan'),1.),(1.,float('inf'))])
def test_invalid_geometry(centre,sigma):
    with pytest.raises(ValueError): reflected_peak(centre,sigma)


def test_sensor_shape_symmetry_and_peak_bound():
    state=np.array([1.,.1,.2,0.,10.,1.,.5])
    z,q=reflected_peak(.5,.3)
    assert scalar_shape(state,.2,.3,0,z)==pytest.approx(q)
    for y in np.linspace(-1,1,7):
        for height in np.linspace(0,2,20):
            signal=scalar_shape(state,.2,.3,y,height)
            assert signal<=q+1e-14
            assert signal==pytest.approx(scalar_shape(state,.2,.3,-y,-height))


class ManufacturedThermodynamics:
    def __init__(self, temperature):
        self.temperature=temperature

    def _condensed_air_state(self,rho,y):
        # Manufactured rho=1+2q. No physical-EOS claim for these test functions.
        return self.temperature((rho-1)/2), np.zeros_like(rho)


class ManufacturedModel:
    rhoa=1.
    def __init__(self, temperature):
        self.thermodynamics=ManufacturedThermodynamics(temperature)

    def section_widths(self,state):
        return .2,.3

    def centre_temperature(self,state):
        return float(self.thermodynamics.temperature(1.))


@pytest.mark.parametrize('law,minimum,location,monotone',[
    (lambda q:300-10*q,280.,2.,True),
    (lambda q:280+10*q,280.,0.,False),
    (lambda q:280+(q-.3)**2,280.,.3,False),
    (lambda q:np.zeros_like(q)+280,280.,0.,True),
])
def test_envelope_endpoints_and_internal_minimum(law,minimum,location,monotone):
    state=np.array([3.,.1,.2,0.,10.,1.,0.])
    result=envelope(ManufacturedModel(law),state)
    assert result['minimum_temperature_K']==pytest.approx(minimum,abs=1e-10)
    assert result['minimum_temperature_shape']==pytest.approx(location,abs=1e-6)
    assert result['sample_monotone_cooling']==monotone


def test_envelope_does_not_clip_invalid_reflected_composition():
    state=np.array([.8,.9,.2,0.,10.,1.,0.])
    with pytest.raises(ValueError,match='nonphysical'):
        envelope(ManufacturedModel(lambda q:300-q),state)
