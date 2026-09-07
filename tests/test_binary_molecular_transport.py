import numpy as np
import pytest
from degali.addons.binary_molecular_transport import (
    TEMPERATURE_K,DP_M2_S_ATM,equimolar_reference_diffusivity,
    binary_gas_species_flux,constant_reference_spreading)


def test_printed_table_units_and_pressure_scaling():
    d=equimolar_reference_diffusivity(TEMPERATURE_K,101325.,allow_flagged_65K=True)
    assert d==pytest.approx(DP_M2_S_ATM,rel=1e-14)
    assert d[-1]==pytest.approx(.7664e-4)
    assert equimolar_reference_diffusivity(90.2,202650.)==pytest.approx(4.5e-6)


def test_midpoint_log_interpolation_and_monotonicity():
    assert equimolar_reference_diffusivity(np.sqrt(90.2*169.3),101325.)==pytest.approx(np.sqrt(9*28.94)*1e-6)
    assert np.all(np.diff(equimolar_reference_diffusivity(np.linspace(77.35,294.8,201),101325.))>0)


@pytest.mark.parametrize('t,p',[(65.25,101325.),(300.,101325.),(90.,0),(90.,np.nan),(np.nan,101325.)])
def test_no_silent_extrapolation(t,p):
    with pytest.raises(ValueError): equimolar_reference_diffusivity(t,p)


def test_binary_mass_conservation_and_species_relabelling():
    args=dict(density=2.,diffusivity=.03,mass_fraction=.2,mass_fraction_gradient=.04,
        log_temperature_gradient=.05,thermal_diffusion_factor=-.7)
    r=binary_gas_species_flux(**args)
    assert r['species1']+r['species2']==0
    assert r['fick']==pytest.approx(-.0024)
    assert r['soret']==pytest.approx(.000336)
    args.update(mass_fraction=.8,mass_fraction_gradient=-.04,thermal_diffusion_factor=.7)
    assert binary_gas_species_flux(**args)['species1']==pytest.approx(-r['species1'])


@pytest.mark.parametrize('y',[0.,1.])
def test_pure_species_no_soret(y):
    r=binary_gas_species_flux(density=1.,diffusivity=1.,mass_fraction=y,
        mass_fraction_gradient=0.,log_temperature_gradient=2.,thermal_diffusion_factor=3.)
    assert r['species1']==0


def test_zero_D_and_reversed_gradient():
    args=dict(density=1.,diffusivity=0.,mass_fraction=.5,mass_fraction_gradient=1.,log_temperature_gradient=2.,thermal_diffusion_factor=3.)
    assert binary_gas_species_flux(**args)['species1']==0
    args['diffusivity']=1
    positive=binary_gas_species_flux(**args)['species1']
    args.update(mass_fraction_gradient=-1.,log_temperature_gradient=-2.)
    assert binary_gas_species_flux(**args)['species1']==-positive


@pytest.mark.parametrize('key,value',[('density',0),('diffusivity',-1),('mass_fraction',1.1),('thermal_diffusion_factor',np.nan)])
def test_invalid_flux(key,value):
    args=dict(density=1.,diffusivity=1.,mass_fraction=.5,mass_fraction_gradient=1.,log_temperature_gradient=2.,thermal_diffusion_factor=3.)
    args[key]=value
    with pytest.raises(ValueError): binary_gas_species_flux(**args)


def test_reference_variance_identity_and_zero_time():
    r=constant_reference_spreading([0.,2.],1.,.03)
    assert r['diffusion_length'][0]==0
    assert r['relative_width_increment'][0]==0
    assert np.isinf(r['diffusivity_for_one_percent_width'][0])
    assert r['diffusion_length'][1]**2==pytest.approx(.12)
    assert r['relative_width_increment'][1]==pytest.approx(np.sqrt(1.12)-1)
    d=r['diffusivity_for_one_percent_width'][1]
    assert constant_reference_spreading(2.,1.,d)['relative_width_increment']==pytest.approx(.01)
