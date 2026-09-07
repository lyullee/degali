import numpy as np
import pytest

from degali.addons.thermochemical_ensemble import extreme
from degali.addons.thermochemical_moment_solver_v2 import whitened_extreme


@pytest.mark.parametrize('maximize',[False,True])
def test_same_objective_and_moments_as_original(maximize):
    y=np.array([0.,1.,.3,.7,.4]);h=np.array([0.,1.,2.,-.5,.6])*1e5
    v=np.array([1.,.5,2.,.8,.9]);w=np.array([.1,.2,.3,.4,0]);c=np.array([1.,4.,9.,16.,-10.])
    kwargs=dict(mean_fraction=w@y,mean_enthalpy=w@h,mean_specific_volume=w@v,values=c,maximize=maximize)
    a=extreme(y,h,v,**kwargs);b=whitened_extreme(y,h,v,**kwargs)
    assert b.objective==pytest.approx(a.objective,abs=1e-8)
    assert b.maximum_scaled_moment_error<1e-10
    assert b.scaled_duality_gap<1e-10


def test_dependent_rows_and_incompatible_target():
    args=dict(mean_fraction=.2,mean_enthalpy=-1e5,mean_specific_volume=.7,values=[2.,-3.,8.])
    b=whitened_extreme([.2]*3,[-1e5]*3,[.7]*3,**args)
    assert b.objective==pytest.approx(-3)
    with pytest.raises(ValueError):
        whitened_extreme([.2]*3,[-1e5]*3,[.7]*3,**dict(args,mean_fraction=.3))


def test_nearly_parallel_moment_rows_retain_exact_four_weights():
    y=np.array([0.,1.,.3,.7]);v=np.array([1.,.5,2.,.8]);perturb=np.array([.5,2.,-.2,3.])
    h=1e6*(2+3*y+5*v+1e-5*perturb);w=np.array([.1,.2,.3,.4]);c=np.array([1.,4.,9.,16.])
    b=whitened_extreme(y,h,v,mean_fraction=w@y,mean_enthalpy=w@h,mean_specific_volume=w@v,values=c)
    assert b.weights==pytest.approx(w,abs=1e-8)
    assert b.objective==pytest.approx(w@c,abs=1e-7)


def test_same_reference_transformation_with_extra_moment():
    y=np.array([0.,1.,.3,.7,.4]);h=np.array([0.,1.,2.,-.5,.6])*1e5
    v=np.array([1.,.5,2.,.8,.9]);w=np.array([.1,.2,.3,.4,0]);c=np.array([1.,4.,9.,16.,-10.])
    extra=np.array([.2,.3,.8,.7,.1]);kwargs=dict(mean_fraction=w@y,mean_enthalpy=w@h,
        mean_specific_volume=w@v,values=c,additional_moments=[extra],additional_targets=[w@extra])
    a=whitened_extreme(y,h,v,**kwargs)
    b=whitened_extreme(y,h+1e6*y-5e5,v,**dict(kwargs,mean_enthalpy=w@h+1e6*(w@y)-5e5))
    assert a.objective==pytest.approx(b.objective,abs=1e-8)
    assert a.maximum_scaled_moment_error<1e-10
