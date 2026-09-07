import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.addons.observed_wind_profile import with_observed_wind
from degali.core.atmosphere import friction_velocity,psi
from degali.core.constants import VKC
from degali.core.jetplume import JetPlume


def base():
    jp=JetPlume.__new__(JetPlume)
    jp.u0,jp.z0,jp.zr,jp.rml=2.4,1.5,.01,0.
    jp.ustar=friction_velocity(jp.u0,jp.z0,jp.zr,jp.rml)
    jp.k=SimpleNamespace(delta=2.15)
    out=IndependentEnergyCrosswind.__new__(IndependentEnergyCrosswind)
    out.jetplume=jp
    return out


def test_same_inputs_recover_original_profile_without_mutation():
    old=base();before=vars(old.jetplume).copy()
    new=with_observed_wind(old,speed=old.jetplume.u0,reference_height=old.jetplume.z0)
    assert new is not old and new.jetplume is not old.jetplume
    for z in (.1,.5,1.5,3.,6.):
        assert new.jetplume._wind_profile(z,.2,.9)==old.jetplume._wind_profile(z,.2,.9)
    assert vars(old.jetplume)==before


def test_observed_speed_is_recovered_at_reference_height():
    old=base()
    new=with_observed_wind(old,speed=1.9222640246473406,reference_height=3.)
    jp=new.jetplume
    recovered=jp.ustar/VKC*(math.log((jp.z0+jp.zr)/jp.zr)-psi(jp.z0,jp.rml))
    assert recovered==pytest.approx(jp.u0,rel=1e-14)
    assert old.jetplume.z0==1.5
    assert old.jetplume.u0==2.4


@pytest.mark.parametrize('speed,height',[(-1.,3.),(math.nan,3.),(math.inf,3.),(1.,0.),(1.,math.nan)])
def test_invalid_observations_rejected(speed,height):
    with pytest.raises(ValueError): with_observed_wind(base(),speed=speed,reference_height=height)


def test_quiescent_observation_is_supported():
    model=with_observed_wind(base(),speed=0.,reference_height=3.)
    assert model.jetplume._wind_profile(.5,.1,1.)==(0.,0.)
