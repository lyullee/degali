"""Conservative stream-function identities and full vector production."""

import numpy as np
import pytest
from numpy.polynomial.legendre import legvander
from test_enriched_transport import phase,candidate,section,mixing,reservoir,projection,transport
from degali.addons.phase_radial_quadrature import gauss_rule
from degali.addons.solenoidal_transport import solenoidal_basis,redistributed_fluxes,SolenoidalTransport


@pytest.mark.parametrize('count',[4,6])
def test_stream_modes_are_divergence_free_and_face_mean_zero(count):
    rng=np.random.default_rng(909)
    a,b=rng.uniform(.05,.9,(2,31))
    h=1e-6
    da=(solenoidal_basis(a+h,b,count)-solenoidal_basis(a-h,b,count))/(2*h)
    db=(solenoidal_basis(a,b+h,count)-solenoidal_basis(a,b-h,count))/(2*h)
    assert da[:,0]+db[:,1]==pytest.approx(np.zeros((len(a),count)),abs=3e-8)
    x,w=gauss_rule(16)
    face=solenoidal_basis(np.ones_like(x),x,count)[:,0]
    expected=legvander(x,2*count)[:,2*np.arange(1,count+1)]*np.sqrt(4*np.arange(1,count+1)+1)
    assert face==pytest.approx(expected,abs=1e-13)
    assert w@face==pytest.approx(np.zeros(count),abs=1e-13)
    normal_axis=solenoidal_basis(np.zeros_like(x),x,count)[:,0]
    assert normal_axis==pytest.approx(np.zeros_like(normal_axis),abs=1e-13)


def test_stream_basis_respects_square_symmetry():
    a,b=np.array([.1,.3,.8]),np.array([.2,.7,.4])
    v=solenoidal_basis(a,b,4)
    assert solenoidal_basis(b,a,4)==pytest.approx(v[:,::-1,:],abs=1e-13)
    assert solenoidal_basis(-a,b,4)[:,0]==pytest.approx(-v[:,0],abs=1e-13)
    assert solenoidal_basis(-a,b,4)[:,1]==pytest.approx(v[:,1],abs=1e-13)


def test_zero_modes_recover_radial_stress_and_production(transport):
    a,b=np.array([.1,.4,.8]),np.array([.2,.5,.9])
    d=transport.local(a,b)
    n=transport.count
    fm,fp=np.zeros((3,n+1)),np.zeros((3,n+1))
    fm[:,0],fp[:,0]=-1.,-.1
    rates=np.zeros(n)
    out=redistributed_fluxes(transport,d,fm,fp,rates,np.zeros(4),np.zeros(4))
    stress=fp[:,0]-d['u']*fm[:,0]
    assert out['production']==pytest.approx(-2*d['q']*d['uq']*stress,rel=1e-13)
    assert out['radial_momentum_diffusivity']==pytest.approx(-stress/(transport.area*d['rho']*d['uq']),rel=1e-13)
    assert max(out['nonradial_stress_fraction'])<1e-14


def test_solenoidal_kinetic_work_integral_is_boundary_flux(transport):
    # A smooth polynomial velocity permits a high-accuracy independent square
    # test of P=-div(u*delta_FP-.5*u²*delta_FM) for divergence-free corrections.
    x,w=gauss_rule(24)
    a,b=np.meshgrid(x,x,indexing='ij')
    a,b=a.ravel(),b.ravel()
    v=solenoidal_basis(a,b,4)
    ma=np.array([.1,-.2,.05,.01]); pa=np.array([.3,.1,-.1,.02])
    fm=np.einsum('nik,k->ni',v,ma); fp=np.einsum('nik,k->ni',v,pa)
    u=3.-a*a-b*b
    production=-np.sum(np.column_stack([-2*a,-2*b])*(fp-u[:,None]*fm),axis=1)
    integral=np.outer(w,w).ravel()@production
    edge=solenoidal_basis(np.ones_like(x),x,4)[:,0]
    ue=2.-x*x
    boundary=2*w@(ue*(edge@pa)-.5*ue**2*(edge@ma))
    assert integral+boundary==pytest.approx(0.,abs=1e-12)


def test_redistribution_keeps_all_modal_rows_or_reports_rank_failure(transport):
    model=SolenoidalTransport(transport.projection,transport.parameters,scalar_mixing=transport.mixing,
        thermal_species_ratio=1.,mechanical_work='reduced_buoyancy_work_immediate_shear_heat',
        flux_modes=4,order=4,angular_order=24)
    try:
        out=model.assemble()
    except ValueError as exc:
        assert 'rank' in str(exc)
        diag=model.linear_diagnostics
        assert diag['matrix'].shape==(29,31)
        assert diag['singular_values'][-1]<=1e-12*diag['singular_values'][0]
    else:
        assert out['matrix'].shape==(29,31)
        assert out['all_original_weak_equations_retained']
        assert not out['scalar_isotropic_momentum_closure']
        assert not out['reynolds_tensor_realizability_established']
        assert not out['adopted']
