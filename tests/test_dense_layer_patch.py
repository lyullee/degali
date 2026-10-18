import numpy as np
import pytest
from degali.addons.binary_integral_puff import BinaryPuffThermodynamics
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.dense_layer_patch import (
    solve_dense_layer_patch,_flux_pair,_transport_multilayer_relative,
    _repair_multilayer_species_bounds,
    _diffuse_horizontal_through_open_faces,
)


def run_layer(*,wall=False,duration=.3,wind=(0.,0.),uniform=False,**extra):
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    q=np.zeros((4,12,12));rho=thermo._density_from_temperature(291.15,.6)
    if uniform:q[0]=rho
    else:q[0,5:7,5:7]=rho
    q[1]=.6*q[0]
    wx=np.zeros((12,13));wy=np.zeros((13,12))
    if wall:wx[:,6]=100.
    return solve_dense_layer_patch(initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),
        cell_size_m=1.,duration_s=duration,wind_xy_m_s=wind,wall_x_height_m=wx,
        wall_y_height_m=wy,receptor_points_m=np.array([[5.5,5.5,.4],[5.5,5.5,2.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,**extra)


@pytest.mark.parametrize('profile,extra', [
    ('dense_power', {'dense_vertical_exponent':2.}),
    ('top_hat', {'cavity_vertical_exponent':12.,'cavity_geometry':'full_disk',
                 'cavity_center_xy_m':(1.5,1.5),'cavity_radius_m':1.}),
    ('exponential_tail', {'vertical_tail_scale_m':1.,'vertical_tail_extension_m':3.}),
])
def test_conservative_observer_reconstructs_species_inventory_without_changing_solver(profile,extra):
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('mixture',.0585,800.))
    density=thermo._density_from_temperature(thermo.ambient_temperature,.2)
    q=np.zeros((4,3,3));q[0]=2.*density;q[1]=.2*q[0]
    z=(np.arange(8192)+.5)*5./8192
    kwargs=dict(initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.001,wind_xy_m_s=(0.,0.),wall_x_height_m=np.zeros((3,4)),
        wall_y_height_m=np.zeros((4,3)),max_dt_s=.001,entrainment_coefficient=0.,
        drag_coefficient=0.,gravity_scale=0.,vertical_observation_profile=profile,
        receptor_points_m=np.column_stack((np.full_like(z,1.5),np.full_like(z,1.5),z)),**extra)
    legacy=solve_dense_layer_patch(**kwargs)
    fixed=solve_dense_layer_patch(**kwargs,vertical_observation_mode='conservative_volume')
    pure_density=thermo._density_from_temperature(thermo.ambient_temperature,1.)
    reconstructed=np.mean(fixed.species_traces_vol_percent[:,0]/100.)*5.*pure_density
    assert reconstructed==pytest.approx(q[1,1,1],rel=.003)
    np.testing.assert_array_equal(fixed.final_conserved,legacy.final_conserved)
    assert fixed.diagnostics['vertical_observation_maximum_relative_volume_residual']<1e-12


def test_conservative_observer_marker_share_is_partitioned_not_double_counted():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('mixture',.0585,800.))
    density=thermo._density_from_temperature(thermo.ambient_temperature,.1)
    q=np.zeros((4,3,3));q[0]=2.*density;q[1]=.1*q[0]
    z=(np.arange(8192)+.5)*10./8192
    fixed=solve_dense_layer_patch(initial_conserved=q,thermo=thermo,
        origin_xy_m=(0.,0.),cell_size_m=1.,duration_s=.001,wind_xy_m_s=(0.,0.),
        wall_x_height_m=np.zeros((3,4)),wall_y_height_m=np.zeros((4,3)),
        receptor_points_m=np.column_stack((np.full_like(z,1.5),np.full_like(z,1.5),z)),
        max_dt_s=.001,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        initial_elevated_species_fraction=.5,wall_impulse_elevated_center_m=2.,
        wall_impulse_elevated_sigma_m=.5,vertical_observation_mode='conservative_volume')
    reconstructed=np.mean(fixed.species_traces_vol_percent[:,0]/100.)*10.
    reconstructed*=thermo._density_from_temperature(thermo.ambient_temperature,1.)
    assert reconstructed==pytest.approx(q[1,1,1],rel=.003)


def test_constant_wind_history_is_bitwise_equivalent_to_constant_input():
    base=run_layer(duration=.3,wind=(1.,.2),uniform=False)
    history=run_layer(
        duration=.3,wind=(1.,.2),uniform=False,
        wind_history_time_s=np.array([0.,.15,.3]),
        wind_history_xy_m_s=np.array([[1.,.2],[1.,.2],[1.,.2]]),
    )
    assert np.array_equal(history.time_s,base.time_s)
    assert np.array_equal(history.final_conserved,base.final_conserved)
    assert np.array_equal(
        history.species_traces_vol_percent,base.species_traces_vol_percent)
    assert history.diagnostics['wind_history_enabled'] is True
    assert history.diagnostics['wind_history_sample_count']==3


def test_time_varying_wind_forcing_changes_drag_path_and_keeps_ledgers_closed():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    q=np.zeros((4,12,12));rho=thermo._density_from_temperature(291.15,.6)
    q[0,5:7,5:7]=rho;q[1]=.6*q[0]
    kwargs=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((12,13)),
        wall_y_height_m=np.zeros((13,12)),
        receptor_points_m=np.array([[5.5,5.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=.5,
    )
    steady=solve_dense_layer_patch(**kwargs)
    varying=solve_dense_layer_patch(
        **kwargs,
        wind_history_time_s=np.array([0.,.15,.3]),
        wind_history_xy_m_s=np.array([[1.,0.],[0.,1.],[-1.,0.]]),
    )
    assert not np.allclose(varying.final_conserved[2:4],steady.final_conserved[2:4])
    assert varying.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(varying.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert varying.diagnostics['momentum_conservation_gate'] is True


def test_gravity_spreading_has_density_feedback_and_conserves_species():
    r=run_layer()
    assert r.final_conserved[1,5,4]>0
    assert r.final_conserved[2,5,4]<0
    assert r.final_conserved[2,5,7]>0
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(r.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert r.diagnostics['density_momentum_feedback'] is True


def test_continuous_layer_observer_distinguishes_point_heights():
    r=run_layer(uniform=True,duration=.02)
    assert r.species_traces_vol_percent[0,0]>0
    assert r.species_traces_vol_percent[1,0]==0


def test_exponential_vertical_observation_tail_is_diagnostic_only():
    r=solve_dense_layer_patch(
        initial_conserved=np.pad(np.array([[[1.0]],[[.6]],[[0.]],[[0.]]]),
                                 ((0,0),(2,2),(2,2))),
        thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.)),
        origin_xy_m=(0.,0.),cell_size_m=1.,duration_s=.02,wind_xy_m_s=(0.,0.),
        wall_x_height_m=np.zeros((5,6)),wall_y_height_m=np.zeros((6,5)),
        receptor_points_m=np.array([[2.5,2.5,.4],[2.5,2.5,2.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='exponential_tail',
        vertical_tail_scale_m=1.,vertical_tail_extension_m=3.)
    assert r.species_traces_vol_percent[0,0]>0
    assert r.species_traces_vol_percent[1,0]>0
    assert r.diagnostics['vertical_observation_profile']=='exponential_tail'
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_diffusive_vertical_observation_tail_starts_delayed_and_preserves_mass():
    r=solve_dense_layer_patch(
        initial_conserved=np.pad(np.array([[[1.0]],[[.6]],[[0.]],[[0.]]]),
                                 ((0,0),(2,2),(2,2))),
        thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.)),
        origin_xy_m=(0.,0.),cell_size_m=1.,duration_s=.2,wind_xy_m_s=(0.,0.),
        wall_x_height_m=np.zeros((5,6)),wall_y_height_m=np.zeros((6,5)),
        receptor_points_m=np.array([[2.5,2.5,.4],[2.5,2.5,2.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='diffusive_tail',
        vertical_tail_extension_m=3.,vertical_tail_diffusivity_m2_s=.5)
    assert r.species_traces_vol_percent[0,0]>0.
    assert r.species_traces_vol_percent[1,0]==pytest.approx(0.)
    assert r.species_traces_vol_percent[1,-1]>0.
    assert r.diagnostics['vertical_observation_profile']=='diffusive_tail'
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_mass_conserving_diffusive_tail_normalizes_resolved_layer():
    kwargs=dict(
        initial_conserved=np.pad(np.array([[[1.0]],[[.6]],[[0.]],[[0.]]]),
                                 ((0,0),(2,2),(2,2))),
        thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.)),
        origin_xy_m=(0.,0.),cell_size_m=1.,duration_s=.2,wind_xy_m_s=(0.,0.),
        wall_x_height_m=np.zeros((5,6)),wall_y_height_m=np.zeros((6,5)),
        receptor_points_m=np.array([[2.5,2.5,.4],[2.5,2.5,2.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='diffusive_tail',
        vertical_tail_extension_m=3.,vertical_tail_diffusivity_m2_s=.5)
    raw=solve_dense_layer_patch(**kwargs)
    conservative=solve_dense_layer_patch(
        **kwargs,mass_conserving_vertical_profile=True)
    assert conservative.species_traces_vol_percent[0,-1] < raw.species_traces_vol_percent[0,-1]
    assert conservative.species_traces_vol_percent[1,-1] > 0.
    assert conservative.diagnostics['mass_conserving_vertical_profile'] is True
    assert conservative.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_multilayer_vertical_transport_conserves_species_and_mixture():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,5));q[0,2:4,2:4]=rho;q[1]=.6*q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(0.,0.),wall_x_height_m=np.zeros((5,6)),
        wall_y_height_m=np.zeros((6,5)),
        receptor_points_m=np.array([[2.5,2.5,.4],[2.5,2.5,2.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='multilayer',vertical_layers=4,
        vertical_mixing_diffusivity_m2_s=.5)
    assert r.diagnostics['vertical_layers']==4
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(r.diagnostics['mixture_mass_residual_kg'])<1e-10


def test_multilayer_observer_volume_averages_partially_occupied_fixed_layer():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,3,3));q[0,1,1]=.25*rho;q[1]=.6*q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.01,wind_xy_m_s=(0.,0.),wall_x_height_m=np.zeros((3,4)),
        wall_y_height_m=np.zeros((4,3)),
        receptor_points_m=np.array([[1.5,1.5,.25]]),max_dt_s=.01,
        entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='multilayer',vertical_layers=2,
        vertical_domain_top_m=1.0,vertical_mixing_diffusivity_m2_s=1.e-12)
    pure=100.*thermo.closure.mole_fraction_from_mass_fraction(.6)
    assert r.species_traces_vol_percent[0,0]==pytest.approx(.5*pure)
    assert r.diagnostics['vertical_observer'].endswith('EOS volume occupancy')


def test_multilayer_entrainment_enters_top_boundary_not_ground_layer():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,3,3));q[0,1,1]=rho;q[1]=.6*q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.01,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((3,4)),
        wall_y_height_m=np.zeros((4,3)),
        receptor_points_m=np.array([[1.5,1.5,.25],[1.5,1.5,.75]]),
        max_dt_s=.001,entrainment_coefficient=1.,drag_coefficient=0.,gravity_scale=0.,
        vertical_observation_profile='multilayer',vertical_layers=2,
        vertical_observation_mode='conservative_volume',
        vertical_domain_top_m=1.,vertical_mixing_diffusivity_m2_s=1.e-12)
    corrected=solve_dense_layer_patch(**common)
    historical=solve_dense_layer_patch(
        **common,multilayer_entrainment_placement='legacy_ground')
    # Both start identically.  Afterwards, top-surface entrainment dilutes the
    # upper layer; the explicitly historical ground placement does the reverse.
    assert corrected.species_traces_vol_percent[:,0]==pytest.approx(
        historical.species_traces_vol_percent[:,0])
    assert corrected.species_traces_vol_percent[0,-1] > historical.species_traces_vol_percent[0,-1]
    assert corrected.species_traces_vol_percent[1,-1] < historical.species_traces_vol_percent[1,-1]
    assert corrected.diagnostics['multilayer_entrainment_placement']=='top'
    assert corrected.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(corrected.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert corrected.diagnostics['momentum_conservation_gate'] is True


def test_declared_initial_vertical_distribution_is_not_clipped_by_column_height():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,3,3));q[0,1,1]=.25*rho;q[1]=.6*q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.01,wind_xy_m_s=(0.,0.),wall_x_height_m=np.zeros((3,4)),
        wall_y_height_m=np.zeros((4,3)),
        receptor_points_m=np.array([[1.5,1.5,2.25]]),max_dt_s=.01,
        entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='multilayer',vertical_layers=4,
        vertical_domain_top_m=4.0,vertical_mixing_diffusivity_m2_s=1.e-12)
    collapsed=solve_dense_layer_patch(**common)
    elevated=solve_dense_layer_patch(
        **common,initial_vertical_gaussian=(2.0,.5))
    assert collapsed.species_traces_vol_percent[0,0]==0.0
    assert elevated.species_traces_vol_percent[0,0]>0.0
    assert elevated.diagnostics['maximum_relative_species_mass_residual']<1.e-12


def test_multilayer_horizontal_flux_carries_donor_elevation_downstream():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,8));q[0,2,1]=.25*rho;q[1]=.6*q[0];q[2]=q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=2.,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((5,9)),
        wall_y_height_m=np.zeros((6,8)),
        receptor_points_m=np.array([[3.5,2.5,2.5],[3.5,2.5,.5]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        vertical_observation_profile='multilayer',vertical_layers=4,
        vertical_domain_top_m=4.,vertical_mixing_diffusivity_m2_s=1.e-12,
        initial_vertical_gaussian=(2.5,.3))
    peak=r.species_traces_vol_percent.max(axis=1)
    assert peak[0]>1.0
    assert peak[0]>1.e5*peak[1]
    assert r.diagnostics['multilayer_horizontal_flux_transport'] is True
    assert r.diagnostics['maximum_multilayer_horizontal_closure_residual_kg_m2']<1.e-12


def test_multilayer_relative_transport_separates_paths_and_conserves_every_component():
    state=np.zeros((4,2,3,7),dtype=float)
    state[0,:,1,3]=1.0
    state[1,:,1,3]=.4
    state[2,:,1,3]=(.25,1.25)
    state[3,:,1,3]=(-.1,.2)
    initial=state.sum(axis=(1,2,3)).copy()
    u=np.zeros((2,3,7));v=np.zeros_like(u)
    u[0]=-1.0;u[1]=1.0
    moved,mass_moved,substeps=_transport_multilayer_relative(
        state,u,v,dt_s=.5,cell_size_m=1.0,
        wall_x_height_m=np.zeros((3,8)),
        wall_y_height_m=np.zeros((4,7)),z_edges_m=np.array([0.,1.,2.]))
    assert mass_moved>0.0
    assert substeps>=2
    assert moved[0,0,1,:3].sum()>0.0
    assert moved[0,1,1,4:].sum()>0.0
    assert moved.sum(axis=(1,2,3))==pytest.approx(initial,abs=1.e-12)
    assert np.min(moved[:2])>=-1.e-12
    assert np.max(moved[1]-moved[0])<=1.e-12


def test_multilayer_species_bound_repair_moves_excess_without_changing_column():
    state=np.zeros((4,3,1,1),dtype=float)
    state[0,:,0,0]=(1.0,.2,.8)
    state[1,:,0,0]=(1.1,-.1,.5)
    before=state[1].sum(axis=0).copy()
    repaired,moved=_repair_multilayer_species_bounds(state)
    assert moved==pytest.approx(.1)
    assert repaired[1].sum(axis=0)==pytest.approx(before)
    assert np.min(repaired[1])>=0.0
    assert np.max(repaired[1]-repaired[0])<=1.e-12


def test_obstacle_lift_advects_multilayer_mass_upward_conservatively():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,6,8));q[0,2,3]=.5*rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=1.5,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((6,9)),
        wall_y_height_m=np.zeros((7,8)),
        receptor_points_m=np.array([[4.5,2.5,1.5]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        vertical_observation_profile='multilayer',vertical_layers=4,
        vertical_domain_top_m=4.,vertical_mixing_diffusivity_m2_s=1.e-12,
        cavity_center_xy_m=(2.5,2.5),cavity_radius_m=1.)
    no_lift=solve_dense_layer_patch(**common)
    lifted=solve_dense_layer_patch(
        **common,obstacle_vertical_lift_fraction=.5,
        obstacle_vertical_lift_decay_m=2.)
    assert lifted.species_traces_vol_percent.max()>1.0
    assert lifted.species_traces_vol_percent.max()>1.e6*no_lift.species_traces_vol_percent.max()
    assert lifted.diagnostics['maximum_obstacle_vertical_velocity_m_s']==pytest.approx(.5)
    assert lifted.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(lifted.diagnostics['mixture_mass_residual_kg'])<1.e-10


def test_wall_overflow_storage_delays_only_crossing_mass_and_closes_ledgers():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,10));q[0,2,2]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((5,11));wx[:,4]=.5
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=3.,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((6,10)),
        receptor_points_m=np.array([[4.5,2.5,.2]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.)
    direct=solve_dense_layer_patch(**common)
    delayed=solve_dense_layer_patch(
        **common,wall_overflow_capture_fraction=1.,
        wall_overflow_return_rate_s_inv=.5)
    direct_trace=direct.species_traces_vol_percent[0]
    delayed_trace=delayed.species_traces_vol_percent[0]
    direct_arrival=direct.time_s[np.flatnonzero(direct_trace>1.)[0]]
    delayed_arrival=delayed.time_s[np.flatnonzero(delayed_trace>1.)[0]]
    direct_centroid=np.trapezoid(direct.time_s*direct_trace,direct.time_s)/np.trapezoid(direct_trace,direct.time_s)
    delayed_centroid=np.trapezoid(delayed.time_s*delayed_trace,delayed.time_s)/np.trapezoid(delayed_trace,delayed.time_s)
    assert delayed_arrival>direct_arrival
    assert delayed_centroid>direct_centroid
    assert delayed.diagnostics['wall_overflow_captured_species_kg']>0.
    assert delayed.diagnostics['wall_overflow_returned_species_kg']>0.
    assert delayed.diagnostics['wall_overflow_stored_species_kg']>0.
    assert delayed.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(delayed.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert delayed.diagnostics['momentum_conservation_gate'] is True


def test_wall_overflow_storage_carries_elevated_subset_and_age():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,10));q[0,2,2]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((5,11));wx[:,4]=.5
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=3.,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((6,10)),
        receptor_points_m=np.array([[4.5,2.5,.2],[4.5,2.5,1.2]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        gravity_scale=0.,initial_elevated_species_fraction=.5,
        wall_impulse_elevated_center_m=1.,
        wall_impulse_elevated_sigma_m=.3,
        wall_overflow_capture_fraction=1.,
        wall_overflow_return_rate_s_inv=.2)
    assert result.diagnostics['wall_overflow_captured_species_kg']>0.
    assert result.diagnostics['wall_overflow_elevated_stored_species_kg']>0.
    assert result.diagnostics['wall_overflow_elevated_age_moment_kg_s']>0.
    assert result.diagnostics['wall_impulse_elevated_marker_active_species_kg']>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_upwind_semicircle_control_surface_captures_outward_side_bypass():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,7,12));q[0,3,3]=2.*rho;q[1]=.6*q[0];q[2]=-1.5*q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=2.,wind_xy_m_s=(-1.,0.),wall_x_height_m=np.zeros((7,13)),
        wall_y_height_m=np.zeros((8,12)),
        receptor_points_m=np.array([[2.5,3.5,.2]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        cavity_center_xy_m=(5.5,3.5),cavity_radius_m=2.,
        wall_overflow_capture_fraction=1.,wall_overflow_return_rate_s_inv=.2,
        wall_overflow_side_gate_length_m=2.)
    assert result.diagnostics['wall_overflow_side_gate_face_count']>0
    assert result.diagnostics['wall_overflow_side_captured_species_kg']>0.
    assert result.diagnostics['wall_overflow_top_captured_species_kg']==pytest.approx(0.)
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_spatially_advected_wind_history_closes_transport_ledgers():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,8));q[0,2,1]=rho;q[1]=.6*q[0];q[2]=q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=1.,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((5,9)),
        wall_y_height_m=np.zeros((6,8)),receptor_points_m=np.array([[3.5,2.5,.2]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=.2,
        gravity_scale=0.,wind_history_time_s=np.array([0.,5.5,11.]),
        wind_history_xy_m_s=np.array([[1.,0.],[.7,.7],[0.,1.]]),
        wind_history_advection_reference_x_m=10.,
        wind_history_advection_speed_m_s=1.)
    assert result.diagnostics['wind_history_spatial_advection_enabled'] is True
    assert result.diagnostics['wind_history_advection_reference_x_m']==pytest.approx(10.)
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_wall_impulse_lift_uses_only_blocked_donor_flow_and_is_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,10));q[0,2,2]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((5,11));wx[:,4]=2.
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=2.,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((6,10)),
        receptor_points_m=np.array([[4.5,2.5,1.5]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        vertical_observation_profile='multilayer',vertical_layers=4,
        vertical_domain_top_m=4.,vertical_mixing_diffusivity_m2_s=1.e-12,
        height_resolved_wall_flux=True,wall_impulse_vertical_lift_fraction=.5)
    assert result.diagnostics['maximum_wall_impulse_vertical_velocity_m_s']>0.
    assert result.diagnostics['wall_impulse_lifted_species_kg']>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_wall_impulse_elevated_branch_changes_only_impacted_vertical_projection():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,10));q[0,2,2]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((5,11));wx[:,4]=2.
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=2.,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((6,10)),
        receptor_points_m=np.array(
            [[x+.5,2.5,2.5] for x in range(10)]+[[1.5,.5,2.5]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.)
    direct=solve_dense_layer_patch(**common)
    elevated=solve_dense_layer_patch(
        **common,wall_impulse_elevated_transfer_fraction=1.,
        wall_impulse_elevated_center_m=2.5,wall_impulse_elevated_sigma_m=.5)
    assert elevated.diagnostics['wall_impulse_elevated_species_transferred_kg']>0.
    assert elevated.species_traces_vol_percent[:,:10].max()>direct.species_traces_vol_percent[:,:10].max()
    assert np.array_equal(
        elevated.species_traces_vol_percent[0],direct.species_traces_vol_percent[0])
    assert elevated.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(elevated.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert elevated.diagnostics['momentum_conservation_gate'] is True


def test_wall_crossing_flux_alone_marks_elevated_branch_and_preserves_solution():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,10));q[0,2,2]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((5,11));wx[:,4]=.5
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=2.,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((6,10)),
        receptor_points_m=np.array([[x+.5,2.5,2.5] for x in range(10)]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.)
    direct=solve_dense_layer_patch(**common)
    elevated=solve_dense_layer_patch(
        **common,wall_flux_elevated_transfer_fraction=1.,
        wall_impulse_elevated_center_m=2.5,wall_impulse_elevated_sigma_m=.5)
    assert elevated.diagnostics['wall_flux_elevated_species_marked_kg']>0.
    assert elevated.diagnostics['wall_impulse_elevated_species_transferred_kg']==0.
    assert not np.array_equal(
        elevated.species_traces_vol_percent,direct.species_traces_vol_percent)
    assert elevated.diagnostics['wall_flux_elevated_relative_transport_kg']>0.
    assert np.all(np.isfinite(elevated.final_conserved))
    assert np.min(elevated.final_conserved[:2])>=-1.e-12
    assert np.max(elevated.final_conserved[1]-elevated.final_conserved[0])<=1.e-12
    assert elevated.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(elevated.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert elevated.diagnostics['momentum_conservation_gate'] is True


def test_wall_crossing_marker_and_hybrid_impact_split_close_without_duplication():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,7,12));q[0,2:5,1:4]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((7,13));wx[:,5]=1.
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=1.2,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((8,12)),
        receptor_points_m=np.array([[x+.5,3.5,z] for x in range(12)
                                    for z in (.4,2.4)]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        gravity_scale=0.,vertical_observation_mode='conservative_volume',
        vertical_observation_component_diagnostics=True,
        wall_flux_elevated_transfer_fraction=1.,
        wall_impulse_elevated_transfer_fraction=1.,
        wall_impulse_elevated_state_mode='hybrid',
        wall_impulse_elevated_conserved_fraction=.4,
        wall_impulse_elevated_return_rate_s_inv=.1,
        wall_impulse_reattached_fraction=.35,
        wall_impulse_reattached_return_rate_s_inv=.05,
        wall_impulse_elevated_center_m=2.5,
        wall_impulse_elevated_sigma_m=.5,
        cavity_center_xy_m=(5.,3.5),cavity_radius_m=3.)
    assert result.diagnostics['wall_flux_elevated_species_marked_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_captured_total_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_marker_active_species_kg']>=0.
    components=result.observation_component_traces_vol_percent
    assert components is not None
    assert np.allclose(
        sum(components.values()),result.species_traces_vol_percent,
        rtol=1.e-12,atol=1.e-12)
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_initial_elevated_species_branch_preserves_total_species_inventory():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,4,6));q[0,1:3,1:3]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.4,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((4,7)),
        wall_y_height_m=np.zeros((5,6)),receptor_points_m=np.array([[1.5,1.5,3.]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.)
    direct=solve_dense_layer_patch(**common)
    elevated=solve_dense_layer_patch(
        **common,initial_elevated_species_fraction=.5,
        wall_impulse_elevated_center_m=3.,wall_impulse_elevated_sigma_m=.5)
    assert elevated.species_traces_vol_percent[0,0]>direct.species_traces_vol_percent[0,0]
    assert elevated.diagnostics['wall_impulse_elevated_species_active_kg']>0.
    assert elevated.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(elevated.diagnostics['mixture_mass_residual_kg'])<1.e-10


def test_elevated_branch_age_changes_vertical_projection_without_mass_change():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,4,8));q[0,1:3,1:3]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=1.,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((4,9)),
        wall_y_height_m=np.zeros((5,8)),
        receptor_points_m=np.array([[2.5,1.5,.4],[2.5,1.5,3.5]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        initial_elevated_species_fraction=.5,
        wall_impulse_elevated_center_m=3.,wall_impulse_elevated_sigma_m=.5)
    fixed=solve_dense_layer_patch(**common)
    aged=solve_dense_layer_patch(
        **common,wall_impulse_elevated_settling_velocity_m_s=.5,
        wall_impulse_elevated_vertical_diffusivity_m2_s=.2,
        wall_impulse_elevated_min_center_m=.2)
    assert aged.diagnostics['wall_impulse_elevated_age_moment_kg_s']>0.
    assert not np.array_equal(
        aged.species_traces_vol_percent,fixed.species_traces_vol_percent)
    assert np.array_equal(aged.final_conserved,fixed.final_conserved)
    assert aged.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(aged.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert aged.diagnostics['momentum_conservation_gate'] is True


def test_elevated_branch_relative_transport_is_bounded_and_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,10));q[0,2,2]=2.*rho;q[1]=.6*q[0];q[2]=1.5*q[0]
    wx=np.zeros((5,11));wx[:,4]=2.
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=2.,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((6,10)),
        receptor_points_m=np.array([[x+.5,2.5,2.5] for x in range(10)]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        wall_impulse_elevated_transfer_fraction=1.,
        wall_impulse_elevated_center_m=2.5,wall_impulse_elevated_sigma_m=.5)
    carried=solve_dense_layer_patch(**common)
    slipped=solve_dense_layer_patch(
        **common,wall_impulse_elevated_transport_fraction=1.)
    assert slipped.diagnostics['wall_impulse_elevated_relative_transport_kg']>0.
    assert not np.array_equal(slipped.final_conserved,carried.final_conserved)
    assert np.min(slipped.final_conserved[1])>=-1.e-10
    assert np.max(slipped.final_conserved[1]-slipped.final_conserved[0])<=1.e-10
    assert slipped.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(slipped.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert slipped.diagnostics['momentum_conservation_gate'] is True


def test_elevated_branch_crosswind_transport_changes_path_and_is_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,7,9));q[0,2:5,2:4]=rho;q[1,3,2:4]=.6*rho;q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.8,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((7,10)),
        wall_y_height_m=np.zeros((8,9)),receptor_points_m=np.array([[4.5,3.5,2.]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        initial_elevated_species_fraction=.5,
        wall_impulse_elevated_center_m=2.,wall_impulse_elevated_sigma_m=.5)
    direct=solve_dense_layer_patch(**common)
    def centroid_y(result):
        species=result.final_conserved[1]
        y=np.arange(species.shape[0],dtype=float)+.5
        return float((species*y[:,None]).sum()/species.sum())
    deflected=solve_dense_layer_patch(
        **common,wall_impulse_elevated_crosswind_velocity_m_s=.5)
    assert centroid_y(deflected)>centroid_y(direct)
    assert deflected.diagnostics['wall_impulse_elevated_relative_transport_kg']>0.
    assert deflected.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(deflected.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert deflected.diagnostics['momentum_conservation_gate'] is True


def test_elevated_branch_crosswind_diffusion_spreads_and_conserves_marker():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,9,9));q[0,2:7,2:5]=rho;q[1,4,2:5]=.6*rho;q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.8,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((9,10)),
        wall_y_height_m=np.zeros((10,9)),receptor_points_m=np.array([[4.5,4.5,2.]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        initial_elevated_species_fraction=1.,
        wall_impulse_elevated_center_m=2.,wall_impulse_elevated_sigma_m=.5)
    direct=solve_dense_layer_patch(**common)
    spread=solve_dense_layer_patch(
        **common,wall_impulse_elevated_crosswind_diffusivity_m2_s=.4)
    def variance_y(result):
        species=result.final_conserved[1]
        y=np.arange(species.shape[0],dtype=float)+.5
        mean=float((species*y[:,None]).sum()/species.sum())
        return float((species*(y[:,None]-mean)**2).sum()/species.sum())
    assert variance_y(spread)>variance_y(direct)
    assert spread.diagnostics['wall_impulse_elevated_crosswind_diffusion_kg']>0.
    assert spread.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(spread.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert spread.diagnostics['momentum_conservation_gate'] is True


def test_elevated_branch_outward_crosswind_splits_about_obstacle_center():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,9,10));q[0,2:7,2:5]=rho
    q[1,3,2:5]=.6*rho;q[1,5,2:5]=.6*rho;q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.8,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((9,11)),
        wall_y_height_m=np.zeros((10,10)),cavity_center_xy_m=(4.5,4.5),
        receptor_points_m=np.array([[4.5,4.5,2.]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        initial_elevated_species_fraction=1.,wall_impulse_elevated_center_m=2.,
        wall_impulse_elevated_sigma_m=.5)
    direct=solve_dense_layer_patch(**common)
    split=solve_dense_layer_patch(
        **common,wall_impulse_elevated_crosswind_outward_velocity_m_s=.5,
        wall_impulse_elevated_crosswind_outward_scale_m=1.)
    decayed=solve_dense_layer_patch(
        **common,wall_impulse_elevated_crosswind_outward_velocity_m_s=.5,
        wall_impulse_elevated_crosswind_outward_scale_m=1.,
        wall_impulse_elevated_crosswind_outward_decay_s=.1)
    def moments(result):
        species=result.final_conserved[1]
        y=np.arange(species.shape[0],dtype=float)+.5
        mean=float((species*y[:,None]).sum()/species.sum())
        variance=float((species*(y[:,None]-mean)**2).sum()/species.sum())
        return mean,variance
    direct_mean,direct_variance=moments(direct)
    split_mean,split_variance=moments(split)
    _,decayed_variance=moments(decayed)
    assert abs(split_mean-direct_mean)<1.e-12
    assert split_variance>direct_variance
    assert direct_variance<decayed_variance<split_variance
    assert split.diagnostics['wall_impulse_elevated_relative_transport_kg']>0.
    assert split.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(split.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert split.diagnostics['momentum_conservation_gate'] is True


def test_elevated_branch_finite_advection_scale_slows_path_and_is_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,12));q[0,2,2:5]=rho;q[1,2,2:5]=.6*rho;q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.8,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((5,13)),
        wall_y_height_m=np.zeros((6,12)),receptor_points_m=np.array([[5.5,2.5,2.]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        initial_elevated_species_fraction=.5,
        wall_impulse_elevated_center_m=2.,wall_impulse_elevated_sigma_m=.5)
    direct=solve_dense_layer_patch(**common)
    slowed=solve_dense_layer_patch(
        **common,wall_impulse_elevated_advection_scale=.5,
        wall_impulse_elevated_advection_duration_s=.5)
    x=np.arange(q.shape[2],dtype=float)+.5
    direct_x=float((direct.final_conserved[1]*x[None,:]).sum()
                   /direct.final_conserved[1].sum())
    slowed_x=float((slowed.final_conserved[1]*x[None,:]).sum()
                   /slowed.final_conserved[1].sum())
    assert slowed_x<direct_x
    assert slowed.diagnostics['wall_impulse_elevated_relative_transport_kg']>0.
    assert slowed.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(slowed.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert slowed.diagnostics['momentum_conservation_gate'] is True


def test_height_resolved_wall_flux_uses_vertical_mass_above_wall_and_closes_ledger():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,5,5));q[0,2:4,1:4]=rho;q[1]=.6*q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,0.),wall_x_height_m=np.full((5,6),.75),
        wall_y_height_m=np.zeros((6,5)),receptor_points_m=np.array([[2.5,2.5,.4],[2.5,2.5,1.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        vertical_observation_profile='multilayer',vertical_layers=4,
        vertical_mixing_diffusivity_m2_s=.5,height_resolved_wall_flux=True)
    assert r.diagnostics['height_resolved_wall_flux'] is True
    assert r.diagnostics['height_resolved_wall_steps']>0
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_tangential_wall_bypass_is_conservative_and_records_wall_impulse():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0]
    wx=np.zeros((8,9));wy=np.zeros((9,8));wx[:,4]=5.;wy[4,:]=5.
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,.2),wall_x_height_m=wx,wall_y_height_m=wy,
        receptor_points_m=np.array([[2.5,2.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        wall_tangential_bypass_fraction=.5,cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    assert r.diagnostics['maximum_wall_bypass_mass_kg']>0.
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_tip_wall_bypass_is_conservative_and_records_wall_impulse():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0]
    wx=np.zeros((8,9));wy=np.zeros((9,8));wx[:,4]=5.;wy[4,:]=5.
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,.2),wall_x_height_m=wx,wall_y_height_m=wy,
        receptor_points_m=np.array([[2.5,2.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        wall_tip_bypass_fraction=.5,wall_tip_bypass_power=2.,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    assert r.diagnostics['maximum_wall_bypass_mass_kg']>0.
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['maximum_momentum_budget_residual_kg_m_s']<1e-9
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_reservoir_capture_annulus_is_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0]
    wx=np.zeros((8,9));wy=np.zeros((9,8));wx[:,4]=5.;wy[4,:]=5.
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,.2),wall_x_height_m=wx,wall_y_height_m=wy,
        receptor_points_m=np.array([[2.5,2.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.01,
        cavity_reservoir_return_rate_s_inv=.01,
        cavity_reservoir_return_sigma_m=2.,cavity_reservoir_return_radius_m=1.,
        cavity_reservoir_capture_radius_m=2.,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    assert r.diagnostics['cavity_reservoir_capture_radius_m']==pytest.approx(2.)
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_local_cavity_reservoir_preserves_spatial_inventory_and_conservation():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.4,wind_xy_m_s=(1.,.2),wall_x_height_m=np.zeros((8,9)),
        wall_y_height_m=np.zeros((9,8)),receptor_points_m=np.array([[2.5,3.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.1,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    mixed=solve_dense_layer_patch(**common)
    local=solve_dense_layer_patch(
        **common,cavity_reservoir_spatial_mode='local')
    assert local.diagnostics['cavity_reservoir_spatial_mode']=='local'
    assert local.diagnostics['local_cavity_reservoir_species_mass_kg']>0.
    assert not np.array_equal(local.final_conserved,mixed.final_conserved)
    assert local.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(local.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert local.diagnostics['momentum_conservation_gate'] is True


def test_local_reservoir_can_advect_conservatively_into_wake_corridor():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,7,14));q[0,3,2:5]=rho;q[1]=.6*q[0];q[2]=q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=3.,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((7,15)),
        wall_y_height_m=np.zeros((8,14)),receptor_points_m=np.array([[8.5,3.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        cavity_reservoir_capture_rate_s_inv=2.,
        cavity_reservoir_return_rate_s_inv=.01,
        cavity_reservoir_spatial_mode='local',
        cavity_reservoir_transport_fraction=1.,
        cavity_reservoir_transport_mask='wake',
        cavity_reservoir_wake_length_m=10.,
        cavity_reservoir_wake_half_width_m=1.5,
        cavity_geometry='full_disk',cavity_center_xy_m=(3.5,3.5),
        cavity_radius_m=1.5)
    assert result.diagnostics['cavity_reservoir_transport_mask']=='wake'
    assert result.diagnostics['local_cavity_reservoir_outside_cavity_mass_kg']>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_local_reservoir_transport_does_not_require_return_to_main_layer():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,7,14));q[0,3,2:5]=rho;q[1]=.6*q[0];q[2]=q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=3.,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((7,15)),
        wall_y_height_m=np.zeros((8,14)),receptor_points_m=np.array([[8.5,3.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        friction_velocity_m_s=.1,roughness_length_m=.01,
        cavity_reservoir_capture_rate_s_inv=2.,
        cavity_reservoir_return_rate_s_inv=0.,
        cavity_reservoir_spatial_mode='local',
        cavity_reservoir_transport_fraction=1.,
        cavity_reservoir_transport_height_m=.5,
        cavity_reservoir_transport_mask='wake',
        cavity_reservoir_wake_length_m=10.,
        cavity_reservoir_wake_half_width_m=1.5,
        cavity_geometry='full_disk',cavity_center_xy_m=(3.5,3.5),
        cavity_radius_m=1.5)
    assert result.diagnostics['local_cavity_reservoir_outside_cavity_mass_kg']>0.
    assert result.diagnostics['local_cavity_reservoir_species_mass_kg']>0.
    assert result.diagnostics[
        'cavity_reservoir_transport_height_m']==pytest.approx(.5)
    assert result.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_wall_impact_reservoir_captures_only_blocked_approach_conservatively():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,10));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    wx=np.zeros((8,11));wx[:,4]=10.0
    wy=np.zeros((9,10))
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_y_height_m=wy,
        receptor_points_m=np.array([[3.5,3.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=5.,
        cavity_reservoir_return_rate_s_inv=.01,
        cavity_reservoir_capture_mode='wall_impact',
        cavity_reservoir_spatial_mode='local',
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    blocked=solve_dense_layer_patch(**common,wall_x_height_m=wx)
    open_flow=solve_dense_layer_patch(
        **common,wall_x_height_m=np.zeros_like(wx))
    assert blocked.diagnostics['cavity_reservoir_capture_mode']=='wall_impact'
    assert blocked.diagnostics['local_cavity_reservoir_species_mass_kg']>0.
    assert open_flow.diagnostics['local_cavity_reservoir_species_mass_kg']==pytest.approx(0.)
    assert blocked.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(blocked.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert blocked.diagnostics['momentum_conservation_gate'] is True


def test_below_wall_core_capture_preserves_prompt_elevated_complement():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,10));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    wx=np.zeros((8,11));wx[:,4]=10.0
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((9,10)),
        receptor_points_m=np.array([[3.5,3.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=5.,
        cavity_reservoir_return_rate_s_inv=.01,
        cavity_reservoir_capture_mode='wall_impact',
        cavity_reservoir_spatial_mode='local',
        initial_elevated_species_fraction=.25,
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    column=solve_dense_layer_patch(**common)
    split=solve_dense_layer_patch(
        **common,cavity_reservoir_vertical_capture_mode='below_wall_core')
    assert split.diagnostics[
        'cavity_reservoir_vertical_capture_mode']=='below_wall_core'
    assert 0.<split.diagnostics['local_cavity_reservoir_species_mass_kg'] \
        <column.diagnostics['local_cavity_reservoir_species_mass_kg']
    assert split.diagnostics['local_cavity_reservoir_elevated_species_kg'] \
        ==pytest.approx(0.)
    assert column.diagnostics['local_cavity_reservoir_elevated_species_kg']>0.
    assert split.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(split.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert split.diagnostics['momentum_conservation_gate'] is True


def test_separate_wall_wake_branch_coexists_with_cavity_reservoir():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,10));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    wx=np.zeros((8,11));wx[:,4]=10.0
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((9,10)),
        receptor_points_m=np.array([[3.5,3.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.1,
        cavity_reservoir_spatial_mode='local',
        wall_impulse_wake_capture_rate_s_inv=2.,
        wall_impulse_wake_return_rate_s_inv=.1,
        wall_impulse_wake_transport_fraction=.3,
        wall_impulse_wake_crosswind_outward_velocity_m_s=.4,
        wall_impulse_wake_crosswind_outward_scale_m=1.,
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    assert result.diagnostics['local_cavity_reservoir_species_mass_kg']>0.
    assert result.diagnostics['wall_impulse_wake_captured_species_kg']>0.
    assert result.diagnostics['wall_impulse_wake_transport_moved_kg']>0.
    assert result.diagnostics[
        'wall_impulse_wake_crosswind_outward_velocity_m_s']==pytest.approx(.4)
    assert result.diagnostics['wall_impulse_wake_stored_species_kg']>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert result.diagnostics['momentum_conservation_gate'] is True
    weighted=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((9,10)),
        receptor_points_m=np.array([[3.5,3.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        wall_impulse_wake_capture_rate_s_inv=2.,
        wall_impulse_wake_capture_mode='species_weighted',
        wall_impulse_wake_return_rate_s_inv=.1,
        wall_impulse_wake_transport_fraction=.3,
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    assert 0.<weighted.diagnostics['wall_impulse_wake_captured_total_kg'] \
        <result.diagnostics['wall_impulse_wake_captured_total_kg']
    assert weighted.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(weighted.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert weighted.diagnostics['momentum_conservation_gate'] is True


def test_wall_wake_age_moment_broadens_vertical_observation_conservatively():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,10));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    wx=np.zeros((8,11));wx[:,4]=10.0
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.6,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((9,10)),
        receptor_points_m=np.array([[3.5,3.5,.4],[3.5,3.5,2.5]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        wall_impulse_wake_capture_rate_s_inv=2.,
        wall_impulse_wake_capture_mode='species_weighted',
        wall_impulse_wake_return_rate_s_inv=.1,
        wall_impulse_wake_observation_center_m=.5,
        wall_impulse_wake_observation_sigma_m=.2,
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    static=solve_dense_layer_patch(**common)
    diffused=solve_dense_layer_patch(
        **common,wall_impulse_wake_vertical_diffusivity_m2_s=2.)
    assert diffused.diagnostics['wall_impulse_wake_mean_species_age_s']>0.
    assert diffused.diagnostics['wall_impulse_wake_species_age_moment_kg_s']>0.
    assert diffused.species_traces_vol_percent[1,-1] \
        > static.species_traces_vol_percent[1,-1]
    assert diffused.species_traces_vol_percent[0,-1] \
        < static.species_traces_vol_percent[0,-1]
    assert diffused.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(diffused.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert diffused.diagnostics['momentum_conservation_gate'] is True


def test_conserved_prompt_elevated_branch_carries_full_state_and_age():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,10));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    wx=np.zeros((8,11));wx[:,4]=10.0
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.6,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((9,10)),
        receptor_points_m=np.array([[3.5,3.5,.4],[3.5,3.5,2.5]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        wall_impulse_elevated_transfer_fraction=1.,
        wall_impulse_elevated_state_mode='conserved',
        wall_impulse_elevated_return_rate_s_inv=.1,
        wall_impulse_elevated_crosswind_velocity_m_s=.2,
        wall_impulse_elevated_center_m=2.,wall_impulse_elevated_sigma_m=.5,
        wall_impulse_wake_capture_rate_s_inv=.2,
        wall_impulse_wake_capture_mode='species_weighted',
        wall_impulse_wake_return_rate_s_inv=.1,
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    assert result.diagnostics['wall_impulse_elevated_state_mode']=='conserved'
    assert result.diagnostics['wall_impulse_elevated_captured_total_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_species_transferred_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_stored_total_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_age_moment_kg_s']>0.
    assert result.diagnostics['wall_impulse_elevated_reservoir_transport_moved_kg']>0.
    assert result.diagnostics['wall_impulse_wake_captured_species_kg']>0.
    assert result.species_traces_vol_percent[1,-1]>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_hybrid_prompt_split_keeps_marker_and_conserved_paths_without_duplication():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,10));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    wx=np.zeros((8,11));wx[:,4]=10.0
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.6,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((9,10)),
        receptor_points_m=np.array([[3.5,3.5,.4],[3.5,3.5,2.5]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        wall_impulse_elevated_transfer_fraction=1.,
        wall_impulse_elevated_state_mode='hybrid',
        wall_impulse_elevated_conserved_fraction=.4,
        wall_impulse_elevated_return_rate_s_inv=.1,
        wall_impulse_elevated_crosswind_velocity_m_s=.2,
        wall_impulse_elevated_center_m=2.,wall_impulse_elevated_sigma_m=.5,
        wall_impulse_elevated_rise_velocity_m_s=.2,
        wall_impulse_elevated_rise_time_s=3.,
        wall_impulse_elevated_conserved_center_m=1.25,
        wall_impulse_elevated_conserved_sigma_m=.4,
        wall_impulse_elevated_conserved_settling_velocity_m_s=.01,
        wall_impulse_elevated_conserved_vertical_diffusivity_m2_s=.02,
        wall_impulse_elevated_conserved_min_center_m=.2,
        wall_impulse_elevated_conserved_rise_velocity_m_s=.3,
        wall_impulse_elevated_conserved_rise_time_s=2.,
        wall_impulse_elevated_conserved_transport_fraction=.3,
        wall_impulse_elevated_conserved_crosswind_velocity_m_s=.15,
        wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s=.04,
        wall_impulse_elevated_conserved_crosswind_outward_scale_m=20.,
        wall_impulse_elevated_conserved_center_mode='impact_head',
        wall_impulse_reattached_fraction=.35,
        wall_impulse_reattached_transport_fraction=.2,
        wall_impulse_reattached_crosswind_velocity_m_s=-.1,
        wall_impulse_reattached_crosswind_outward_velocity_m_s=.05,
        wall_impulse_reattached_return_rate_s_inv=.05,
        cavity_center_xy_m=(3.,4.),cavity_radius_m=4.)
    assert result.diagnostics['wall_impulse_elevated_state_mode']=='hybrid'
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_fraction']==pytest.approx(.4)
    assert result.diagnostics['wall_impulse_elevated_captured_total_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_stored_total_kg']>0.
    assert result.diagnostics['wall_impulse_elevated_marker_active_species_kg']>0.
    assert result.diagnostics[
        'wall_impulse_elevated_rise_velocity_m_s']==pytest.approx(.2)
    assert result.diagnostics[
        'wall_impulse_elevated_rise_time_s']==pytest.approx(3.)
    assert result.diagnostics['wall_impulse_elevated_conserved_active_species_kg']>0.
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_center_m']==pytest.approx(1.25)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_sigma_m']==pytest.approx(.4)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_settling_velocity_m_s']==pytest.approx(.01)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_vertical_diffusivity_m2_s']==pytest.approx(.02)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_min_center_m']==pytest.approx(.2)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_rise_velocity_m_s']==pytest.approx(.3)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_rise_time_s']==pytest.approx(2.)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_transport_fraction']==pytest.approx(.3)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_crosswind_velocity_m_s']==pytest.approx(.15)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s']==pytest.approx(.04)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_crosswind_outward_scale_m']==pytest.approx(20.)
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_center_mode']=='impact_head'
    assert result.diagnostics[
        'wall_impulse_elevated_conserved_mean_center_m']>0.
    assert result.diagnostics[
        'wall_impulse_reattached_fraction']==pytest.approx(.35)
    assert result.diagnostics['wall_impulse_reattached_stored_total_kg']>0.
    assert result.diagnostics['wall_impulse_reattached_active_species_kg']>0.
    assert result.diagnostics[
        'wall_impulse_reattached_reservoir_transport_moved_kg']>0.
    assert result.species_traces_vol_percent[1,-1]>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_local_reservoir_crosswind_convergence_moves_toward_wake_axis_conservatively():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,9,10));q[0,1,2:5]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=1.,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((9,11)),
        wall_y_height_m=np.zeros((10,10)),receptor_points_m=np.array([[5.5,4.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        cavity_reservoir_capture_rate_s_inv=2.,
        cavity_reservoir_return_rate_s_inv=.01,
        cavity_reservoir_spatial_mode='local',cavity_geometry='full_disk',
        cavity_center_xy_m=(4.5,4.5),cavity_radius_m=4.5)
    direct=solve_dense_layer_patch(**common)
    focused=solve_dense_layer_patch(
        **common,cavity_reservoir_crosswind_convergence_rate_s_inv=.5)
    y=np.arange(q.shape[1],dtype=float)+.5
    direct_y=float((direct.final_conserved[1]*y[:,None]).sum()
                   /direct.final_conserved[1].sum())
    focused_y=float((focused.final_conserved[1]*y[:,None]).sum()
                    /focused.final_conserved[1].sum())
    assert focused_y>direct_y
    assert focused.diagnostics[
        'cavity_reservoir_crosswind_convergence_rate_s_inv']==pytest.approx(.5)
    assert focused.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert abs(focused.diagnostics['mixture_mass_residual_kg'])<1.e-10
    assert focused.diagnostics['momentum_conservation_gate'] is True


def test_local_cavity_reservoir_can_contribute_to_sensor_observation():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=1.,wind_xy_m_s=(0.,0.),wall_x_height_m=np.zeros((8,9)),
        wall_y_height_m=np.zeros((9,8)),receptor_points_m=np.array([[2.5,3.5,.95]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.1,
        cavity_reservoir_spatial_mode='local',
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    hidden=solve_dense_layer_patch(
        **common,local_reservoir_observation_fraction=0.)
    visible=solve_dense_layer_patch(
        **common,local_reservoir_observation_fraction=1.)
    assert visible.diagnostics['local_reservoir_observation_fraction']==pytest.approx(1.)
    assert visible.species_traces_vol_percent.max()>hidden.species_traces_vol_percent.max()
    assert np.array_equal(visible.final_conserved,hidden.final_conserved)
    assert visible.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(visible.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert visible.diagnostics['momentum_conservation_gate'] is True


def test_local_reservoir_can_return_into_height_advected_reattached_branch():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,12));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.5,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((8,13)),
        wall_y_height_m=np.zeros((9,12)),
        receptor_points_m=np.array([[5.5,3.5,.7]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        friction_velocity_m_s=.2,roughness_length_m=.01,
        cavity_reservoir_capture_rate_s_inv=1.,
        cavity_reservoir_return_rate_s_inv=1.,
        cavity_reservoir_return_target='reattached',
        cavity_reservoir_spatial_mode='local',
        cavity_reservoir_transport_fraction=.1,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.,
        wall_impulse_elevated_state_mode='conserved',
        wall_impulse_elevated_transfer_fraction=.1,
        wall_impulse_reattached_fraction=.5,
        wall_impulse_reattached_transport_fraction=1.,
        wall_impulse_reattached_transport_height_m=2.,
        wall_impulse_reattached_return_rate_s_inv=0.)
    assert result.diagnostics['cavity_reservoir_return_target']=='reattached'
    assert result.diagnostics[
        'wall_impulse_reattached_transport_height_m']==pytest.approx(2.)
    assert result.diagnostics['wall_impulse_reattached_stored_total_kg']>0.
    assert result.diagnostics[
        'wall_impulse_reattached_reservoir_transport_moved_kg']>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_local_reservoir_return_can_split_main_and_reattached_carriers():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,12));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.5,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((8,13)),
        wall_y_height_m=np.zeros((9,12)),
        receptor_points_m=np.array([[5.5,3.5,.7]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        friction_velocity_m_s=.2,roughness_length_m=.01,
        cavity_reservoir_capture_rate_s_inv=1.,
        cavity_reservoir_return_rate_s_inv=1.,
        cavity_reservoir_return_target='reattached',
        cavity_reservoir_spatial_mode='local',
        cavity_reservoir_transport_fraction=.1,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.,
        wall_impulse_elevated_state_mode='conserved',
        wall_impulse_elevated_transfer_fraction=.1,
        wall_impulse_reattached_fraction=.5,
        wall_impulse_reattached_transport_fraction=1.,
        wall_impulse_reattached_transport_height_m=2.,
        wall_impulse_reattached_return_rate_s_inv=0.)
    split=solve_dense_layer_patch(
        **common,cavity_reservoir_reattached_fraction=.25)
    full=solve_dense_layer_patch(
        **common,cavity_reservoir_reattached_fraction=1.)
    assert split.diagnostics[
        'cavity_reservoir_reattached_fraction']==pytest.approx(.25)
    assert 0.<split.diagnostics['wall_impulse_reattached_captured_total_kg']<(
        full.diagnostics['wall_impulse_reattached_captured_total_kg'])
    assert split.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(split.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert split.diagnostics['momentum_conservation_gate'] is True


def test_species_only_local_reservoir_delays_scalar_without_removing_carrier():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.4,wind_xy_m_s=(1.,.2),wall_x_height_m=np.zeros((8,9)),
        wall_y_height_m=np.zeros((9,8)),receptor_points_m=np.array([[2.5,3.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.1,
        cavity_reservoir_spatial_mode='local',
        cavity_reservoir_content_mode='species_only',
        local_reservoir_observation_fraction=1.,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    result=solve_dense_layer_patch(**common)
    assert result.diagnostics['cavity_reservoir_content_mode']=='species_only'
    assert result.diagnostics['local_cavity_reservoir_species_mass_kg']>0.
    assert result.diagnostics['local_cavity_reservoir_mass_kg']==pytest.approx(0.)
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_local_reservoir_preserves_elevated_species_marker():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0];q[2]=q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.4,wind_xy_m_s=(1.,.2),wall_x_height_m=np.zeros((8,9)),
        wall_y_height_m=np.zeros((9,8)),receptor_points_m=np.array([[2.5,3.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.1,
        cavity_reservoir_spatial_mode='local',
        cavity_reservoir_transport_fraction=.5,
        initial_elevated_species_fraction=.25,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    assert result.diagnostics['local_cavity_reservoir_elevated_species_kg']>0.
    assert result.diagnostics['cavity_reservoir_transport_fraction']==pytest.approx(.5)
    assert result.diagnostics['wall_impulse_elevated_species_active_kg']>0.
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(result.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert result.diagnostics['momentum_conservation_gate'] is True


def test_wind_aligned_cavity_return_is_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,10,10));q[0,4:6,3:7]=rho;q[1]=.6*q[0]
    wx=np.zeros((10,11));wy=np.zeros((11,10))
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,.3),wall_x_height_m=wx,wall_y_height_m=wy,
        receptor_points_m=np.array([[5.5,5.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.5,
        cavity_return_streamwise_bias=.8,
        cavity_align_with_wind=True,
        cavity_center_xy_m=(5.,5.),cavity_radius_m=4.)
    assert r.diagnostics['cavity_align_with_wind'] is True
    assert r.diagnostics['cavity_return_streamwise_bias']==pytest.approx(.8)
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_delayed_reservoir_return_keeps_pending_inventory_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,8,8));q[0,3:5,1:4]=rho;q[1]=.6*q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,.2),wall_x_height_m=np.zeros((8,9)),
        wall_y_height_m=np.zeros((9,8)),receptor_points_m=np.array([[2.5,2.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_reservoir_capture_rate_s_inv=.5,
        cavity_reservoir_return_rate_s_inv=.5,
        cavity_reservoir_return_delay_s=1.,cavity_reservoir_return_sigma_m=2.,
        cavity_center_xy_m=(4.,4.),cavity_radius_m=4.)
    assert r.diagnostics['cavity_reservoir_return_delay_s']==pytest.approx(1.)
    assert r.diagnostics['cavity_reservoir_pending_mass_kg']>0.
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_directional_wake_diagnostic_is_local_and_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,4,8));q[0,1:3,1:3]=rho;q[1]=.6*q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((4,9)),
        wall_y_height_m=np.zeros((5,8)),receptor_points_m=np.array([[6.5,1.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        wake_deficit_strength=.5,wake_length_m=4.,wake_base_width_m=1.,
        wake_spread_slope=.2,wake_center_xy_m=(1.,1.),wake_radius_m=1.)
    assert r.diagnostics['wake_deficit_strength']==pytest.approx(.5)
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_convected_wake_meander_is_opt_in_and_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,6,10));q[0,2:4,2:5]=rho;q[1]=.6*q[0];q[2]=q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.4,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((6,11)),
        wall_y_height_m=np.zeros((7,10)),receptor_points_m=np.array([[6.5,2.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=2.,
        wake_deficit_strength=.5,wake_length_m=8.,wake_base_width_m=2.,
        wake_spread_slope=.1,wake_center_xy_m=(1.,3.),wake_radius_m=1.)
    baseline=solve_dense_layer_patch(**common)
    explicit_zero=solve_dense_layer_patch(**common,wake_meander_strength=0.)
    active=solve_dense_layer_patch(
        **common,wake_meander_strength=.4,wake_shedding_strouhal=.2,
        wake_meander_phase_deg=90.)
    np.testing.assert_array_equal(explicit_zero.final_conserved,baseline.final_conserved)
    assert not np.array_equal(active.final_conserved[3],baseline.final_conserved[3])
    assert active.diagnostics['wake_meander_strength']==pytest.approx(.4)
    assert active.diagnostics['wake_shedding_strouhal']==pytest.approx(.2)
    assert active.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert abs(active.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert active.diagnostics['momentum_conservation_gate'] is True


def test_boundary_current_speed_is_exported_with_species_packets():
    r=run_layer(duration=.2,wind=(1.,0.))
    assert r.boundary_current_speed_m_s.shape==r.boundary_species_rate_kg_s.shape
    assert np.all(np.isfinite(r.boundary_current_speed_m_s))
    assert np.all(r.boundary_current_speed_m_s>=0.)


@pytest.mark.parametrize('axis',[0,1])
def test_diffusive_flux_obeys_finite_wall_height_and_porosity(axis):
    q=np.zeros((4,1,2));q[:,0,0]=[4.,2.,-3.,1.]
    height=np.array([[2.,0.]])
    if axis==1:
        q=q.transpose(0,2,1);height=height.T
    ny,nx=height.shape
    wx=np.zeros((ny,nx+1));wy=np.zeros((ny+1,nx))
    def step(wall,blockage=1.):
        if axis==0:wx[:,1]=wall
        else:wy[1,:]=wall
        return _diffuse_horizontal_through_open_faces(
            q,height,dxx=1.,dyy=1.,dt_s=.1,cell_size_m=1.,
            wall_x_height_m=wx,wall_y_height_m=wy,
            wall_blockage_fraction=blockage)[0]
    opened=step(0.)
    closed=step(3.)
    half=step(1.)
    porous=step(3.,.5)
    assert np.array_equal(closed,q)
    assert np.allclose(half-q,.5*(opened-q))
    assert np.allclose(porous,half)
    for result in (opened,closed,half,porous):
        assert np.allclose(result.sum(axis=(1,2)),q.sum(axis=(1,2)))
        assert np.all(result[:2]>=0.)
        assert np.all(result[1]<=result[0])


def test_solver_diffusion_cannot_leak_under_a_solid_wall():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    q=np.zeros((4,6,8));q[0,2:4,3]=1.;q[1]=.6*q[0]
    wx=np.zeros((6,9));wx[:,4]=100.
    kwargs=dict(initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),
        cell_size_m=1.,duration_s=.2,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=np.zeros((7,8)),receptor_points_m=np.array([[4.5,2.5,.1]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        horizontal_diffusivity_m2_s=1.)
    legacy=solve_dense_layer_patch(**kwargs)
    fixed=solve_dense_layer_patch(**kwargs,horizontal_diffusion_respect_walls=True)
    assert legacy.final_conserved[1,:,4:].sum()>0.
    assert fixed.final_conserved[1,:,4:].sum()==0.
    assert fixed.diagnostics['horizontal_diffusion_blocked_species_throughput_kg']>0.
    assert fixed.diagnostics['maximum_relative_species_mass_residual']<1.e-12
    assert fixed.diagnostics['momentum_conservation_gate'] is True


def test_horizontal_diffusion_is_conservative_and_bounded():
    r=run_layer(duration=.3,wind=(1.,0.))
    # Re-run with a stable diagnostic diffusivity on the same local patch.
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,12,12));q[0,5:7,5:7]=rho;q[1]=.6*q[0]
    wx=np.zeros((12,13));wy=np.zeros((13,12))
    d=solve_dense_layer_patch(initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),
        cell_size_m=1.,duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,
        wall_y_height_m=wy,receptor_points_m=np.array([[5.5,5.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        horizontal_diffusivity_m2_s=2.)
    assert d.diagnostics['horizontal_diffusivity_m2_s']==pytest.approx(2.)
    assert d.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert np.min(d.final_conserved[:2])>=-1e-12


def test_entrainment_relaxes_to_declared_far_field_without_species_source():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,12,12));q[0,5:7,5:7]=rho;q[1]=.6*q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((12,13)),
        wall_y_height_m=np.zeros((13,12)),
        receptor_points_m=np.array([[5.5,5.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=.3,drag_coefficient=0.)
    constant=solve_dense_layer_patch(**common)
    relaxed=solve_dense_layer_patch(
        **common,entrainment_far_coefficient=.05,
        entrainment_relaxation_time_s=.1)
    assert relaxed.mixture_mass_kg[-1] < constant.mixture_mass_kg[-1]
    assert relaxed.diagnostics['entrainment_far_coefficient']==pytest.approx(.05)
    assert relaxed.diagnostics['entrainment_relaxation_time_s']==pytest.approx(.1)
    assert relaxed.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert relaxed.diagnostics['momentum_conservation_gate'] is True


def test_entrainment_relaxation_rejects_nonpositive_time_scale():
    with pytest.raises(ValueError,match='relaxation_time'):
        run_layer(
            entrainment_far_coefficient=.05,
            entrainment_relaxation_time_s=0.0)


def test_cavity_vertical_grounding_relaxes_radially_without_changing_mass():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,11,11));q[0]=5.0*rho;q[1]=.6*q[0]
    result=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.02,wind_xy_m_s=(0.,0.),wall_x_height_m=np.zeros((11,12)),
        wall_y_height_m=np.zeros((12,11)),
        receptor_points_m=np.array([[5.5,5.5,2.4],[9.5,5.5,2.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_center_xy_m=(5.5,5.5),cavity_radius_m=5.0,
        cavity_vertical_exponent=10.0,
        cavity_vertical_radial_decay_power=1.0)
    assert (result.species_traces_vol_percent[0,0]
            < result.species_traces_vol_percent[1,0])
    assert result.diagnostics['cavity_vertical_radial_decay_power']==pytest.approx(1.)
    assert result.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_exterior_entrainment_relaxes_with_distance_and_conserves_species():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,12,12));q[0,5:7,5:7]=rho;q[1]=.6*q[0]
    common=dict(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.3,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((12,13)),
        wall_y_height_m=np.zeros((13,12)),
        receptor_points_m=np.array([[5.5,5.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=.3,drag_coefficient=0.,
        cavity_center_xy_m=(0.,0.),cavity_radius_m=1.)
    constant=solve_dense_layer_patch(**common)
    relaxed=solve_dense_layer_patch(
        **common,exterior_entrainment_far_factor=.2,
        exterior_entrainment_decay_length_m=1.)
    assert relaxed.mixture_mass_kg[-1] < constant.mixture_mass_kg[-1]
    assert relaxed.diagnostics['exterior_entrainment_far_factor']==pytest.approx(.2)
    assert relaxed.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_subfence_species_flux_is_reflected_not_deleted():
    r=run_layer(wall=True)
    left=r.final_conserved[1,:,:6].sum();right=r.final_conserved[1,:,6:].sum()
    assert left==pytest.approx(right,rel=1e-12)
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12


def test_light_species_is_not_silently_treated_as_grounded_dense_gas():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('H2',.002016,14000.))
    q=np.zeros((4,2,2));q[0]=1.;q[1]=.1
    with pytest.raises(ValueError,match='positively buoyant'):
        solve_dense_layer_patch(initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
            duration_s=.1,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((2,3)),
            wall_y_height_m=np.zeros((3,2)),receptor_points_m=np.array([[.5,.5,.4]]))


def test_face_species_cannot_cross_high_wall_but_can_overtop_low_wall():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    left=np.zeros((4,1));left[0]=2*rho;left[1]=.6*left[0];left[2]=left[0]
    right=np.zeros_like(left)
    closed_l,closed_r=_flux_pair(left,right,0,thermo,np.array([5.]))
    assert closed_l[1,0]==closed_r[1,0]==0.
    opened_l,opened_r=_flux_pair(left,right,0,thermo,np.array([1.]))
    assert opened_l[1,0]>0
    assert opened_l[1,0]==opened_r[1,0]
    ground_loaded_l,ground_loaded_r=_flux_pair(
        left,right,0,thermo,np.array([1.]),wall_overtop_exponent=2.)
    assert 0.0 < ground_loaded_l[1,0] < opened_l[1,0]
    assert ground_loaded_l[1,0]==ground_loaded_r[1,0]


@pytest.mark.parametrize('power,expected', [(0., .5), (2., .125), (3., .0625)])
def test_wall_overtopping_is_the_integral_of_the_declared_vertical_profile(power,expected):
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    left=np.zeros((4,1));left[0]=2.*rho;left[1]=.6*left[0];left[2]=left[0]
    right=np.zeros_like(left)
    open_species=_flux_pair(left,right,0,thermo,np.array([0.]))[0][1,0]
    wall_species=_flux_pair(
        left,right,0,thermo,np.array([1.]),wall_overtop_exponent=power)[0][1,0]
    assert wall_species/open_species==pytest.approx(expected,rel=1.e-13)


def test_component_diagnostics_are_returned_and_close_the_total_trace():
    result=run_layer(
        duration=.02,vertical_observation_mode='conservative_volume',
        vertical_observation_component_diagnostics=True)
    components=result.observation_component_traces_vol_percent
    assert components is not None
    assert set(components)=={'core'}
    np.testing.assert_allclose(
        sum(components.values()),result.species_traces_vol_percent,
        rtol=1.e-12,atol=1.e-12)


def test_gravity_ablation_removes_pressure_not_species_inventory():
    from degali.addons.dense_layer_patch import _primitive
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    q=np.array([[[4.]],[[2.]],[[0.]],[[0.]]])
    normal=_primitive(q,thermo)
    off=_primitive(q,thermo,0.)
    assert normal[4][0,0]>0
    assert off[4][0,0]==0
    assert off[3][0,0]==0
    assert np.array_equal(normal[0],off[0])


def test_interior_interface_packet_replaces_downstream_flux_without_double_counting():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,2,4));q[0,:,0]=rho;q[1,:,0]=.6*rho
    wx=np.zeros((2,5));wy=np.zeros((3,4))
    calls=[]
    def packet(time_s,dt_s,ny,cell_size_m):
        calls.append((time_s,dt_s))
        out=np.zeros((4,ny));out[0]=.02;out[1]=.012
        return out
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,0.),wall_x_height_m=wx,wall_y_height_m=wy,
        receptor_points_m=np.array([[3.5,.5,.4]]),max_dt_s=.02,
        entrainment_coefficient=0.,drag_coefficient=0.,interface_x_m=2.,
        interface_source_callback=packet)
    assert calls
    assert r.diagnostics['interface_injected_total_kg']==pytest.approx(.02*.2*2)
    assert r.diagnostics['interface_injected_species_kg']==pytest.approx(.012*.2*2)
    assert r.diagnostics['interface_negative_mass_kg']<=1e-12
    assert abs(r.diagnostics['mixture_mass_residual_kg'])<1e-10
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['momentum_conservation_gate'] is True


def test_cavity_geometry_and_backflow_rate_are_conservative():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    q=np.zeros((4,6,8));q[0,2:4,1:4]=rho;q[1]=.6*q[0]
    r=solve_dense_layer_patch(
        initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,
        duration_s=.2,wind_xy_m_s=(1.,0.),wall_x_height_m=np.zeros((6,9)),
        wall_y_height_m=np.zeros((7,8)),receptor_points_m=np.array([[2.5,2.5,.4]]),
        max_dt_s=.02,entrainment_coefficient=0.,drag_coefficient=0.,
        cavity_geometry='downwind_semidisk',cavity_center_xy_m=(3.,3.),
        cavity_radius_m=2.,cavity_backflow_strength=.7,
        cavity_backflow_rate_s_inv=.01)
    assert r.diagnostics['cavity_geometry']=='downwind_semidisk'
    assert r.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert r.diagnostics['maximum_momentum_budget_residual_kg_m_s']<1e-10


def test_wind_aligned_cavity_backflow_is_rotation_covariant():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    rho=thermo._density_from_temperature(291.15,.6)
    cells=10
    forward_q=np.zeros((4,cells,cells))
    forward_q[0,4:6,1:4]=rho
    forward_q[1]=.6*forward_q[0]
    forward_q[2]=forward_q[0]
    reverse_q=forward_q[:,:,::-1][:,::-1,:].copy()
    reverse_q[2:]*=-1.

    common=dict(
        thermo=thermo,origin_xy_m=(0.,0.),cell_size_m=1.,duration_s=.4,
        wall_x_height_m=np.zeros((cells,cells+1)),
        wall_y_height_m=np.zeros((cells+1,cells)),max_dt_s=.01,
        entrainment_coefficient=0.,drag_coefficient=0.,gravity_scale=0.,
        cavity_geometry='downwind_semidisk',cavity_align_with_wind=True,
        cavity_center_xy_m=(5.,5.),cavity_radius_m=4.,
        cavity_backflow_strength=.8,cavity_backflow_rate_s_inv=.5)
    forward=solve_dense_layer_patch(
        **common,initial_conserved=forward_q,wind_xy_m_s=(1.,0.),
        receptor_points_m=np.array([[3.5,5.5,.4]]))
    reverse=solve_dense_layer_patch(
        **common,initial_conserved=reverse_q,wind_xy_m_s=(-1.,0.),
        receptor_points_m=np.array([[6.5,4.5,.4]]))
    reverse_in_forward_frame=(
        reverse.final_conserved[:,:,::-1][:,::-1,:].copy())
    reverse_in_forward_frame[2:]*=-1.

    assert np.array_equal(forward.final_conserved,reverse_in_forward_frame)
    assert np.array_equal(
        forward.species_traces_vol_percent,
        reverse.species_traces_vol_percent)
    assert forward.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert reverse.diagnostics['maximum_relative_species_mass_residual']<1e-12
    assert forward.diagnostics['momentum_conservation_gate'] is True
    assert reverse.diagnostics['momentum_conservation_gate'] is True


def test_cavity_geometry_rejects_unknown_sector():
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    q=np.zeros((4,2,2));q[0]=1.;q[1]=.5
    with pytest.raises(ValueError,match='unsupported cavity_geometry'):
        solve_dense_layer_patch(initial_conserved=q,thermo=thermo,origin_xy_m=(0.,0.),
            cell_size_m=1.,duration_s=.1,wind_xy_m_s=(1.,0.),
            wall_x_height_m=np.zeros((2,3)),wall_y_height_m=np.zeros((3,2)),
            receptor_points_m=np.array([[.5,.5,.4]]),cavity_geometry='bad')
