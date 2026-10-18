from dataclasses import replace
import numpy as np
import pytest

from degali.addons.integral_obstacle_patch import (
    GaussianIntegralInlet, downstream_trace, solve_integral_obstacle_patch,
)
from degali.addons.resolved_obstacle_flow_3d import ResolvedObstacleFlow3DConfig
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import CircularArcWall


def test_integral_face_packet_keeps_tails_and_exact_source_clock():
    source = SourceRateSchedule((0., 1., 4.), (2., 0., 0.))
    inlet = GaussianIntegralInlet(source, 4., 2., .1, 0., 1., .2)
    broad = inlet.area_probabilities(np.array([-20.,20.]), np.array([0.,20.]))
    narrow = inlet.area_probabilities(np.array([0.,.2]), np.array([0.,.2]))
    assert broad.sum() == pytest.approx(1.)
    assert 0 < narrow.sum() < .1  # No renormalization of truncated plume.
    assert inlet.mass_between(0.,2.) == 0
    assert inlet.mass_between(2.,3.) == pytest.approx(2.)
    assert inlet.mass_between(3.,4.) == 0


@pytest.fixture(scope='module')
def inlet_witness():
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(16.,8.,8.),cells=(16,8,8),wind_speed_m_s=1.,
        kinematic_viscosity_m2_s=.5,flow_steps=60,minimum_flow_steps=60,
        lattice_inlet_speed=.03,duration_s=2.,scalar_time_step_s=.1,
        scalar_diffusivity_m2_s=.03,
        source_schedule=SourceRateSchedule((0.,2.),(0.,0.)),
        source_position_m=(1.,4.,1.),source_sigma_m=1.,obstacles=(),
    )
    return solve_integral_obstacle_patch(config,
        incoming_flux=lambda t,dt,y,z: np.full((len(z)-1,len(y)-1),1e-4))


def test_face_packet_actually_enters_field_and_closes_independent_ledgers(inlet_witness):
    result = inlet_witness
    assert result.inlet_mass_kg[-1] == pytest.approx(.0128)
    assert result.domain_mass_kg[-1] > .01
    assert np.max(result.final_kg_m3) > 0
    assert np.max(np.abs(result.mass_residual_kg)) < 1e-15
    assert result.diagnostics['exterior_time_marched_grid_cells'] == 0


def test_downstream_characteristic_has_no_prearrival_signal(inlet_witness):
    trace = downstream_trace(inlet_witness,point_m=(20.,4.,1.),
                             wind_speed_m_s=1.,diffusivity_m2_s=.03)
    assert np.all(trace == 0)  # 4 s delay exceeds this 2 s run.


def test_invalid_face_flux_is_rejected():
    inlet = GaussianIntegralInlet(SourceRateSchedule((0.,2.),(1.,0.)),1.,1.,.1,0.,1.,.2)
    with pytest.raises(ValueError,match='strictly increasing'):
        inlet.area_probabilities(np.array([1.,0.]),np.array([0.,1.]))


def test_initial_inventory_and_embedded_fence_are_in_conservation_ledger():
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(32.,32.,16.), cells=(16,16,8), wind_speed_m_s=1.,
        kinematic_viscosity_m2_s=.5,flow_steps=60,minimum_flow_steps=60,
        lattice_inlet_speed=.03,duration_s=.2,scalar_time_step_s=.1,
        scalar_diffusivity_m2_s=.03,
        source_schedule=SourceRateSchedule((0.,.2),(0.,0.)),
        source_position_m=(16.,16.,1.),source_sigma_m=1.,
        obstacles=(CircularArcWall((16.,16.),8.,0.,0.,4.,-90.,90.),),
        circular_arc_representation='zero_thickness_surface',
    )
    initial=np.zeros((8,16,16))
    initial[0:2,7:9,7:9]=.02
    result=solve_integral_obstacle_patch(
        config,incoming_flux=lambda t,dt,y,z:np.zeros((len(z)-1,len(y)-1)),
        initial_species_mass_kg_m3=initial)
    assert result.domain_mass_kg[0] == pytest.approx(initial.sum()*8)
    assert result.diagnostics['embedded_surface_blocked_face_count']>0
    assert result.inlet_mass_kg[-1]==0
    assert np.max(np.abs(result.mass_residual_kg))<1e-12
    assert np.all(result.boundary_rate_kg_s['z_min']==0)


def test_open_ground_is_rejected_before_expensive_carrier_solve():
    config=ResolvedObstacleFlow3DConfig(domain_m=(8.,8.,8.),cells=(8,8,8),
        wind_speed_m_s=1.,kinematic_viscosity_m2_s=.5,flow_steps=60,minimum_flow_steps=60,
        lattice_inlet_speed=.03,duration_s=1.,scalar_time_step_s=.1,scalar_diffusivity_m2_s=.03,
        source_schedule=SourceRateSchedule((0.,1.),(0.,0.)),source_position_m=(1.,4.,1.),
        source_sigma_m=1.,obstacles=(),nonstreamwise_wall_boundary='all_pressure_open')
    with pytest.raises(ValueError,match='open ground carrier'):
        solve_integral_obstacle_patch(config,incoming_flux=lambda *a:0.)


def test_ground_free_slip_keeps_open_top_faces_but_zero_ground_faces():
    from degali.addons.resolved_obstacle_flow_3d import _raw_faces_from_cell_velocity
    field=np.ones((8,8,8));solid=np.zeros_like(field,dtype=bool)
    u,v,w,_=_raw_faces_from_cell_velocity(u=field,v=field,w=field,solid=solid,
        inlet_physical=field[:,:,0],spacing_m=1.,
        nonstreamwise_wall_boundary='ground_free_slip_lateral_top_pressure_open')
    assert np.all(w[0]==0)
    assert np.all(w[-1]==1)
    assert np.all(v[:,0,:]==1)
