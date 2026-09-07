import numpy as np
import pytest
from degali.addons.stress_realizability import shear_tke_lower_bound


def test_random_positive_covariances_obey_shear_tke_bound():
    rng=np.random.default_rng(922)
    b=rng.normal(size=(100,3,3)); covariance=b@np.swapaxes(b,1,2)
    area=.02; rho=rng.uniform(.5,2.,100); lengths=np.array([.2,.1])
    stress=area*rho[:,None]*covariance[:,0,1:]/lengths
    out=shear_tke_lower_bound(stress,area,rho,lengths)
    tke=.5*np.trace(covariance,axis1=1,axis2=2)
    assert np.all(tke>=out['minimum_tke']-1e-13)
    assert out['physical_shear_covariance']==pytest.approx(covariance[:,0,1:],abs=1e-13)


def test_sharp_bound_is_attained_by_positive_semidefinite_covariance():
    stress=np.array([[.1,.2],[-.4,.3],[0.,0.]])
    out=shear_tke_lower_bound(stress,.02,np.array([1.,.5,2.]),[.3,.1])
    matrix=out['attaining_covariance']
    assert np.min(np.linalg.eigvalsh(matrix))>=-1e-13
    assert matrix[:,0,1:]==pytest.approx(out['physical_shear_covariance'],abs=1e-13)
    assert .5*np.trace(matrix,axis1=1,axis2=2)==pytest.approx(out['minimum_tke'],abs=1e-13)
    assert not out['covariance_is_physical_closure']


def test_zero_density_is_rejected_not_clipped():
    with pytest.raises(ValueError,match='positive'):
        shear_tke_lower_bound(np.array([[1.,2.]]),.02,np.array([0.]),[.3,.1])
