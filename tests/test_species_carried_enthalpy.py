import numpy as np
import pytest

from degali.addons.species_carried_enthalpy import binary_enthalpy_flux
from degali.addons.thermal_moments import eddy_enthalpy_flux


def flux(rho=1.2,dy=.01,dt=.02,gh=4e4,gy=.3,hy=-2e5):
    return binary_enthalpy_flux(density=rho,species_diffusivity=dy,thermal_diffusivity=dt,
        enthalpy_gradient=gh,mass_fraction_gradient=gy,enthalpy_composition_derivative=hy)


@pytest.mark.parametrize('dy,dt',[(0.,0.),(0.,.1),(.1,0.),(.1,.1),(.1,.07),(.1,.2)])
def test_two_equivalent_forms_and_species_sum(dy,dt):
    f=flux(dy=dy,dt=dt)
    assert f.total_enthalpy_flux==pytest.approx(f.naive_enthalpy_flux+f.composition_correction_flux)
    assert f.total_enthalpy_flux==pytest.approx(f.heat_only_flux+f.species_enthalpy_flux)
    assert f.species_mass_flux+(-f.species_mass_flux)==0


@pytest.mark.parametrize('d',[0.,1e-5,.01,1.])
def test_equal_diffusivity_recovers_existing_operator(d):
    f=flux(dy=d,dt=d)
    assert f.composition_correction_flux==0
    assert f.total_enthalpy_flux==pytest.approx(eddy_enthalpy_flux(1.2,4e4,d))


def test_isothermal_composition_gradient_carries_enthalpy_not_heat():
    f=flux(gh=-2e5*.3)
    assert f.heat_only_flux==0
    assert f.total_enthalpy_flux==pytest.approx(-2e5*f.species_mass_flux)


def test_uniform_composition_is_pure_heat_transport():
    f=flux(gy=0.)
    assert f.species_mass_flux==0
    assert f.species_enthalpy_flux==0
    assert f.total_enthalpy_flux==pytest.approx(f.naive_enthalpy_flux)


@pytest.mark.parametrize('offset',[-1e9,-2e5,0.,1e6,1e9])
def test_reference_change_is_species_flux_covariant(offset):
    a=flux()
    b=flux(gh=4e4+offset*.3,hy=-2e5+offset)
    assert b.heat_only_flux==pytest.approx(a.heat_only_flux,abs=1e-7)
    assert b.total_enthalpy_flux==pytest.approx(a.total_enthalpy_flux+offset*a.species_mass_flux,abs=1e-7)


def test_unequal_old_enthalpy_only_operator_fails_reference_covariance():
    rho,dy,dt,gy,gh,offset=1.2,.01,.02,.3,4e4,1e6
    old=eddy_enthalpy_flux(rho,gh,dt)
    changed=eddy_enthalpy_flux(rho,gh+offset*gy,dt)
    required=old+offset*(-rho*dy*gy)
    assert changed-required==pytest.approx(-rho*(dt-dy)*offset*gy)
    assert abs(changed-required)>1


def two_cell_step(offset):
    # Closed1m cells, rho=1. No reaction, advection or boundary flux.
    rho=1.;cp=2000.;base_species_offset=2e6
    t=np.array([100.,150.]);y=np.array([.8,.2])
    hy=base_species_offset+offset
    h=cp*t+hy*y+7e4
    f=flux(rho=rho,dy=.03,dt=.02,gh=h[1]-h[0],gy=y[1]-y[0],hy=hy)
    step=.1
    yn=y+step*np.array([-f.species_mass_flux,f.species_mass_flux])/rho
    hn=h+step*np.array([-f.total_enthalpy_flux,f.total_enthalpy_flux])/rho
    tn=(hn-hy*yn-7e4)/cp
    assert sum(yn)==pytest.approx(sum(y),abs=1e-14)
    assert sum(hn)==pytest.approx(sum(h),rel=1e-15)
    return tn


@pytest.mark.parametrize('offset',[-1e9,0.,1e9])
def test_finite_volume_temperature_independent_of_species_reference(offset):
    assert two_cell_step(offset)==pytest.approx([100.1,149.9],abs=1e-9)


def test_array_broadcast():
    f=flux(rho=np.ones((3,1)),dy=np.array([.01,.02]),dt=.03,
           gh=np.array([[3.],[4.],[5.]]),gy=.3,hy=-2e5)
    assert f.total_enthalpy_flux.shape==(3,2)
    for i in range(3):
        for j in range(2):
            assert f.total_enthalpy_flux[i,j]==pytest.approx(flux(rho=1,dy=[.01,.02][j],dt=.03,gh=i+3).total_enthalpy_flux)


@pytest.mark.parametrize('kwargs',[{'rho':0.},{'rho':-1.},{'dy':-.1},{'dt':-.1},
    {'gh':float('nan')},{'gy':float('inf')},{'hy':float('nan')}])
def test_invalid_inputs(kwargs):
    with pytest.raises(ValueError): flux(**kwargs)
