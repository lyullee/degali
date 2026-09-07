"""Parcel identities and independently soluble finite-support moment problems."""
import numpy as np
import pytest

from degali.addons.thermochemical_ensemble import ensemble,extreme


def test_direct_mass_volume_and_averages():
    m=np.array([2.,1.,3.]);t=np.array([60.,150.,300.]);y=np.array([.6,.2,0.])
    h=np.array([-8e5,-1e5,0.]);v=np.array([.2,.8,1.])
    a=ensemble(m,t,y,h,v);p=m*v/sum(m*v)
    assert a['mean_density']==pytest.approx(sum(m)/sum(m*v))
    assert a['reynolds_temperature']==pytest.approx(p@t)
    assert a['favre_temperature']==pytest.approx(m@t/sum(m))
    assert p@(1/v)==pytest.approx(a['mean_density'])
    assert (p@(y/v))/a['mean_density']==pytest.approx(a['favre_mass_fraction'])
    assert (p@(h/v))/a['mean_density']==pytest.approx(a['favre_enthalpy'])
    assert a['favre_temperature']!=a['reynolds_temperature']


def test_equal_density_and_single_state_limits():
    a=ensemble([1,2],[70,150],[.4,.1],[-1e5,-3e4],[.8,.8])
    assert a['favre_temperature']==a['reynolds_temperature']
    b=ensemble([0,2],[70,150],[.4,.1],[-1e5,-3e4],[.2,.8])
    assert b['favre_temperature']==b['reynolds_temperature']==150
    assert b['favre_temperature_variance']==b['favre_fraction_variance']==0


def test_mass_scaling_and_enthalpy_reference():
    args=(np.array([80.,200.]),np.array([.4,.1]),np.array([-3e5,-1e5]),np.array([.2,.8]))
    a=ensemble([1,2],*args);b=ensemble([10,20],*args)
    for key in ('favre_enthalpy','mean_density','favre_temperature','reynolds_temperature'):
        assert a[key]==pytest.approx(b[key])
    offset,constant=1e6,-5e5
    c=ensemble([1,2],args[0],args[1],args[2]+offset*args[1]+constant,args[3])
    assert c['favre_enthalpy']==pytest.approx(a['favre_enthalpy']+offset*a['favre_mass_fraction']+constant)
    assert c['reynolds_temperature']==a['reynolds_temperature']


@pytest.mark.parametrize('m,t,y,h,v',[
    ([],[],[],[],[]),([-1,2],[60,70],[.4,.2],[0,0],[1,1]),
    ([0,0],[60,70],[.4,.2],[0,0],[1,1]),([1,1],[0,70],[.4,.2],[0,0],[1,1]),
    ([1,1],[60,70],[.4,1.2],[0,0],[1,1]),([1,1],[60,70],[.4,.2],[0,0],[1,0]),
    ([1,1],[60,70],[.4,.2],[0,float('nan')],[1,1]),
])
def test_invalid_parcels(m,t,y,h,v):
    with pytest.raises(ValueError): ensemble(m,t,y,h,v)


def test_unique_four_state_solution():
    y=np.array([0.,1.,.3,.7]);h=np.array([0.,1.,2.,-.5]);v=np.array([1.,.5,2.,.8])
    w=np.array([.1,.2,.3,.4]);c=np.array([1.,4.,9.,16.])
    for maximize in (False,True):
        r=extreme(y,h,v,mean_fraction=w@y,mean_enthalpy=w@h,mean_specific_volume=w@v,values=c,maximize=maximize)
        assert r.weights==pytest.approx(w,abs=1e-9)
        assert r.objective==pytest.approx(w@c,abs=1e-9)
        assert r.maximum_scaled_moment_error<1e-10


def test_known_nonunique_extrema():
    # Repeated thermodynamic moments leave an independently known interval.
    y=np.full(3,.2);h=np.full(3,-1e5);v=np.full(3,.7);c=np.array([2.,-3.,8.])
    args=dict(mean_fraction=.2,mean_enthalpy=-1e5,mean_specific_volume=.7,values=c)
    assert extreme(y,h,v,**args).objective==pytest.approx(-3)
    assert extreme(y,h,v,**args,maximize=True).objective==pytest.approx(8)


def test_reference_covariant_extrema_and_support_extension():
    y=np.array([0.,1.,.3,.7,.4]);h=np.array([0.,1.,2.,-.5,.6])*1e5
    v=np.array([1.,.5,2.,.8,.9]);w=np.array([.1,.2,.3,.4,0]);c=np.array([1.,4.,9.,16.,-10.])
    args=dict(mean_fraction=w@y,mean_enthalpy=w@h,mean_specific_volume=w@v,values=c)
    a=extreme(y,h,v,**args)
    offset,constant=1e6,-5e5
    b=extreme(y,h+offset*y+constant,v,**dict(args,mean_enthalpy=args['mean_enthalpy']+offset*args['mean_fraction']+constant))
    assert a.objective==pytest.approx(b.objective,abs=1e-8)
    assert a.objective<=w@c+1e-10
    assert a.scaled_duality_gap<1e-10


def test_eulerian_cdf_matches_ensemble():
    t=np.array([60.,80.,200.,300.]);y=np.array([0.,1.,.3,.7]);h=np.array([0.,1.,2.,-.5])
    v=np.array([1.,.5,2.,.8]);w=np.array([.1,.2,.3,.4]);vm=w@v
    r=extreme(y,h,v,mean_fraction=w@y,mean_enthalpy=w@h,mean_specific_volume=vm,
        values=v/vm*(t<=80),maximize=True)
    e=ensemble(r.weights,t,y,h,v)
    assert r.objective==pytest.approx(sum(e['volume_weights'][t<=80]))
    assert r.objective!=pytest.approx(sum(w[t<=80]))


def test_infeasible_and_invalid_target_rejected():
    args=dict(mean_fraction=.9,mean_enthalpy=0.,mean_specific_volume=1.,values=[1.,2.])
    with pytest.raises(ValueError): extreme([.1,.2],[0.,0.],[1.,1.],**args)
    with pytest.raises(ValueError): extreme([.1,.2],[0.,0.],[1.,1.],**dict(args,mean_specific_volume=0))


def test_additional_observable_constraint():
    y=np.full(3,.2);h=np.full(3,-1e5);v=np.full(3,.7);c=np.array([2.,-3.,8.])
    a=extreme(y,h,v,mean_fraction=.2,mean_enthalpy=-1e5,mean_specific_volume=.7,
        values=c,additional_moments=[c],additional_targets=[1.])
    assert a.objective==pytest.approx(1.)
    assert a.moment_residuals.shape==(5,)
    with pytest.raises(ValueError):
        extreme(y,h,v,mean_fraction=.2,mean_enthalpy=-1e5,mean_specific_volume=.7,
            values=c,additional_moments=[c],additional_targets=[])
