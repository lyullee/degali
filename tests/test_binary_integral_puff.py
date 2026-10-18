from dataclasses import replace
import math
import numpy as np
import pytest
from degali.addons.binary_integral_puff import initial_prism_puff,binary_puff_traces,axisymmetric_puff_inventory
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_puff import GaussianPuffConfig
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud


@pytest.fixture
def puff():
    closure=BinaryIdealGasMixtureClosure('R12',.120913,606.0748985077907)
    cloud=RegularPolygonPrismInitialCloud((20.,20.),14.,0.,13.,.32,291.15)
    config=GaussianPuffConfig(duration_s=2.,time_step_s=.125,wind_velocity_m_s=(3.,.4,0.),width_cutoff_sigma=math.sqrt(3))
    return initial_prism_puff(cloud,closure,config)


def test_initial_inventory_and_resting_source_are_preserved(puff):
    s=puff.states[0]
    assert s.hydrogen_mass_kg==pytest.approx(3290.1158326487234,rel=1e-10)
    assert s.centre_position_m[2]==pytest.approx(6.5)
    assert s.bulk_velocity_m_s==(0.,0.,0.)
    assert s.bulk_density_kg_m3/puff.thermodynamics.ambient_density>2.
    assert puff.maximum_relative_hydrogen_residual==0.


def test_vectorized_generic_r12_observer_matches_native_observer(puff):
    points=np.array([[22.,22.,.4],[22.,22.,2.4],[25.,23.,8.]])
    values=binary_puff_traces(puff,points)
    for i,point in enumerate(points):
        for j in (0,5,16):
            mole,_=puff.observation_at(puff.states[j].elapsed_s,point)
            assert values[i,j]==pytest.approx(100*mole,rel=2e-10,abs=1e-10)


def test_integral_to_cell_inventory_keeps_truncation_exterior_not_renormalized(puff):
    edges=np.linspace(0.,40.,17);ze=np.linspace(0.,30.,13)
    field,exterior=axisymmetric_puff_inventory(puff,10,edges,edges,ze)
    mapped=field.sum()*2.5**3
    assert exterior>0
    assert mapped+exterior==pytest.approx(puff.states[10].hydrogen_mass_kg,abs=1e-10)
    assert mapped<puff.states[10].hydrogen_mass_kg


def test_anisotropic_rotated_inventory_is_rejected(puff):
    puff.states=(replace(puff.states[0],lateral_half_width_m=10.),)
    with pytest.raises(ValueError,match='axisymmetric'):
        axisymmetric_puff_inventory(puff,0,[0.,40.],[0.,40.],[0.,30.])


def test_bad_dense_entrainment_mode_is_rejected():
    with pytest.raises(ValueError,match='dense_entrainment_model'):
        GaussianPuffConfig(duration_s=1.,dense_entrainment_model='fit_to_sensors').validate()


def test_grounded_dense_puff_does_not_numerically_lift_when_width_changes(puff):
    assert all(s.ground_contact for s in puff.states)
    assert all(s.centre_position_m[2]==pytest.approx(s.vertical_half_width_m) for s in puff.states)
    assert all(s.bulk_velocity_m_s[2]==0. for s in puff.states)
