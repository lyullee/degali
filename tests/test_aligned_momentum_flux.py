import numpy as np
import pytest
from degali.addons.aligned_momentum_flux import aligned_mode_vectors,aligned_mode_vectors_analytic
from degali.addons.phase_radial_quadrature import gauss_rule


PARAMS=dict(q0=3.,ambient_parallel=2.,excess_velocity=12.,velocity_exponent=.7,flux_modes=4)


def test_aligned_momentum_extension_is_divergence_free():
    a,b=np.array([.1,.3,.6,.85]),np.array([.2,.4,.7,.25])
    h=1e-6
    ap=aligned_mode_vectors(a+h,b,**PARAMS)['momentum_modes']
    am=aligned_mode_vectors(a-h,b,**PARAMS)['momentum_modes']
    bp=aligned_mode_vectors(a,b+h,**PARAMS)['momentum_modes']
    bm=aligned_mode_vectors(a,b-h,**PARAMS)['momentum_modes']
    div=((ap-am)[:,0]+(bp-bm)[:,1])/(2*h)
    assert div==pytest.approx(np.zeros_like(div),abs=2e-8)


def test_stress_remains_parallel_to_normalized_velocity_gradient():
    a,b=np.array([0.,.2,.7,1.]),np.array([0.,.3,.4,1.])
    out=aligned_mode_vectors(a,b,**PARAMS)
    u=2.+12.*np.exp(-.7*3.*(a*a+b*b))
    stress=out['momentum_modes']-u[:,None,None]*out['mass_modes']
    cross=a[:,None]*stress[:,1]-b[:,None]*stress[:,0]
    assert cross==pytest.approx(np.zeros_like(cross),abs=1e-13)
    assert np.all(out['radial_stress_modes'][0]==0.)
    assert not out['circulation_amplitudes_determined']
    assert not out['positive_diffusivity_established']


def test_both_flux_corrections_have_zero_integrated_normal_flow():
    x,w=gauss_rule(48)
    out=aligned_mode_vectors(np.ones_like(x),x,**PARAMS)
    for key in ('mass_modes','momentum_modes'):
        assert w@out[key][:,0]==pytest.approx(np.zeros(4),abs=2e-12)


def test_gaussian_closed_integrals_match_independent_radial_quadrature():
    rng=np.random.default_rng(992)
    a,b=rng.uniform(0.,1.,(2,51))
    a=np.r_[0.,1e-8,a]; b=np.r_[0.,1e-8,b]
    quad=aligned_mode_vectors(a,b,**PARAMS,order=48)
    exact=aligned_mode_vectors_analytic(a,b,**PARAMS)
    for key in ('mass_modes','momentum_modes','radial_stress_modes'):
        assert exact[key]==pytest.approx(quad[key],rel=1e-10,abs=3e-11)
