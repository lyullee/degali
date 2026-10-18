from dataclasses import replace
import math
import numpy as np
import pytest
from degali.addons.binary_integral_puff import initial_prism_puff
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.finite_puff import GaussianPuffConfig
from degali.addons.dense_prism_integral import (
    integrate_inertial_prism,prism_from_gaussian_trajectory,prism_layer_inventory,
)


@pytest.fixture
def source():
    cloud=RegularPolygonPrismInitialCloud((60.,320.),14.,0.,13.,.32,291.15)
    closure=BinaryIdealGasMixtureClosure('R12',.120913,606.0748985077907,
        ambient_temperature_k=291.15,ambient_pressure_pa=100500.)
    return cloud,closure


def test_initial_uniform_profile_inventory_and_rest_are_one_consistent_state(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=.1)
    points=np.array([[60.,320.,6.5],[66.,320.,12.9],[68.,320.,1.],[60.,320.,13.1]])
    assert result.traces(points)[:,0]==pytest.approx([32.,32.,0.,0.])
    s=result.states[0]
    assert s.species_mass_kg==pytest.approx(3290.1158326487234,rel=1e-10)
    assert s.radial_velocity_m_s==0.
    assert s.apothem_m==7.
    assert s.height_m==pytest.approx(13.)
    assert s.bulk_velocity_xy_m_s==(0.,0.)


def test_cosine_dense_vertical_profile_is_smooth_and_mass_mean_preserving(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=.1,vertical_profile='cosine_dense',
        vertical_profile_amplitude=.4)
    points=np.array([[60.,320.,0.],[60.,320.,6.5],[60.,320.,13.]])
    values=result.traces(points)[:,0]
    assert values[0]>values[1]>values[2]
    assert values[1]==pytest.approx(32.,abs=1e-12)
    assert result.diagnostics['vertical_profile']=='cosine_dense'
    assert result.diagnostics['vertical_profile_amplitude']==.4
    _,ledger=prism_layer_inventory(result,0,np.arange(20.,102.5,2.5),np.arange(280.,362.5,2.5))
    assert ledger['vertical_profile_transferred'] is False
    with pytest.raises(ValueError,match='admissible mass fraction'):
        integrate_inertial_prism(
            cloud,closure,duration_s=.1,vertical_profile='cosine_dense',
            vertical_profile_amplitude=.6)


def test_compact_legacy_bridge_changes_geometry_not_total_or_width_trajectory(source):
    cloud,closure=source
    puff=initial_prism_puff(cloud,closure,GaussianPuffConfig(duration_s=.2,time_step_s=.01,
        width_cutoff_sigma=math.sqrt(3),wind_velocity_m_s=(3.9,0.,0.)))
    compact=prism_from_gaussian_trajectory(puff,cloud)
    assert compact.traces(np.array([[60.,320.,6.5]]))[0,0]==pytest.approx(32.)
    for a,b in zip(puff.states,compact.states):
        assert a.total_mass_kg==b.total_mass_kg
        assert a.hydrogen_mass_kg==b.species_mass_kg
        assert a.centre_position_m[:2]==b.centre_xy_m


def test_exact_polygon_inventory_and_affine_radial_momentum_survive_handoff(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=1.)
    xe=np.arange(20.,102.5,2.5);ye=np.arange(280.,362.5,2.5)
    q,ledger=prism_layer_inventory(result,-1,xe,ye)
    s=result.states[-1]
    assert q[1].sum()*2.5**2==pytest.approx(s.species_mass_kg,rel=1e-12)
    assert q[0].sum()*2.5**2==pytest.approx(s.total_mass_kg,rel=1e-12)
    assert q[2:].sum(axis=(1,2))*2.5**2==pytest.approx(np.array(s.bulk_velocity_xy_m_s)*s.total_mass_kg,rel=1e-12,abs=1e-9)
    radial_p=q[2]-s.bulk_velocity_xy_m_s[0]*q[0]
    assert radial_p.min()<0
    assert radial_p.max()>0
    assert np.max(np.abs(ledger['exterior_conserved']))<1e-8


def test_truncated_prism_mass_is_reported_not_renormalized(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=.1)
    q,ledger=prism_layer_inventory(result,0,[60.,70.],[310.,320.])
    assert 0<ledger['mapped_conserved'][1]<result.states[0].species_mass_kg
    assert ledger['exterior_conserved'][1]>0
    assert ledger['mapped_conserved'][1]+ledger['exterior_conserved'][1]==pytest.approx(result.states[0].species_mass_kg)


def test_without_entrainment_gravity_conversion_conserves_radial_energy(source):
    cloud,closure=source
    r=integrate_inertial_prism(cloud,closure,duration_s=10.,dt_s=.0125,
        wind_xy_m_s=(0.,0.),entrainment_shear=0.,entrainment_ambient=0.,drag_coefficient=0.)
    assert r.diagnostics['maximum_relative_radial_energy_ledger_residual']<1e-7
    assert r.diagnostics['maximum_relative_species_residual']==0.
    assert r.states[-1].apothem_m>r.states[0].apothem_m
    assert r.states[-1].height_m<r.states[0].height_m
    assert r.states[-1].total_mass_kg==r.states[0].total_mass_kg
    assert r.states[-1].centre_xy_m==r.states[0].centre_xy_m


def test_vertical_inertia_reduces_initial_acceleration(source):
    cloud,closure=source
    on=integrate_inertial_prism(cloud,closure,duration_s=.1)
    off=integrate_inertial_prism(cloud,closure,duration_s=.1,vertical_inertia=False)
    assert on.states[-1].radial_velocity_m_s<off.states[-1].radial_velocity_m_s
    assert on.states[-1].apothem_m<off.states[-1].apothem_m


def test_face_specific_entrainment_preserves_legacy_closure_and_audits_surfaces(source):
    cloud,closure=source
    legacy=integrate_inertial_prism(cloud,closure,duration_s=1.,dt_s=.0125)
    explicit=integrate_inertial_prism(cloud,closure,duration_s=1.,dt_s=.0125,
        entrainment_top_shear=.1,entrainment_top_ambient=.1,
        entrainment_side_shear=.1,entrainment_side_ambient=.1)
    assert explicit.states[-1].total_mass_kg==pytest.approx(legacy.states[-1].total_mass_kg,rel=1e-13)
    assert explicit.states[-1].species_mass_kg==pytest.approx(legacy.states[-1].species_mass_kg,rel=1e-13)
    assert explicit.states[-1].centre_xy_m==pytest.approx(legacy.states[-1].centre_xy_m,abs=1e-12)
    surface=explicit.diagnostics['entrained_air_by_surface_kg']
    assert surface['top']>0 and surface['side']>0
    assert surface['top']+surface['side']==pytest.approx(explicit.states[-1].entrained_air_kg,rel=1e-13)
    top_off=integrate_inertial_prism(cloud,closure,duration_s=1.,dt_s=.0125,
        entrainment_top_shear=0.,entrainment_top_ambient=0.)
    assert top_off.diagnostics['entrained_air_by_surface_kg']['top']==pytest.approx(0.)
    assert top_off.diagnostics['entrained_air_by_surface_kg']['side']>0


def test_parker_bulk_suppression_is_selectable_and_reduces_dense_entrainment(source):
    cloud,closure=source
    neutral=integrate_inertial_prism(cloud,closure,duration_s=1.,dt_s=.0125)
    suppressed=integrate_inertial_prism(
        cloud,closure,duration_s=1.,dt_s=.0125,
        dense_entrainment_model='parker_bulk_suppression')
    assert suppressed.states[-1].entrained_air_kg < neutral.states[-1].entrained_air_kg
    assert suppressed.diagnostics['dense_entrainment_model']=='parker_bulk_suppression'
    assert 0.0 < suppressed.diagnostics['last_entrainment_suppression'] <= 1.0
    with pytest.raises(ValueError,match='dense_entrainment_model'):
        integrate_inertial_prism(cloud,closure,dense_entrainment_model='unknown')


def test_parker_surface_selective_suppression_keeps_the_other_surface_neutral(source):
    cloud,closure=source
    top=integrate_inertial_prism(
        cloud,closure,duration_s=1.,dt_s=.0125,
        dense_entrainment_model='parker_top_suppression')
    side=integrate_inertial_prism(
        cloud,closure,duration_s=1.,dt_s=.0125,
        dense_entrainment_model='parker_side_suppression')
    assert top.diagnostics['last_top_entrainment_suppression'] < 1.0
    assert top.diagnostics['last_side_entrainment_suppression'] == 1.0
    assert side.diagnostics['last_top_entrainment_suppression'] == 1.0
    assert side.diagnostics['last_side_entrainment_suppression'] < 1.0


def test_wake_tail_profile_is_downwind_only_and_zero_at_tail_end(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=2.,wind_xy_m_s=(3.9,0.),
        front_profile='wake_tail',wake_extension_fraction=1.0)
    state=result.states[-1]
    points=np.array([
        [state.centre_xy_m[0],state.centre_xy_m[1],1.],
        [state.centre_xy_m[0]-state.apothem_m*2.0,state.centre_xy_m[1],1.],
        [state.centre_xy_m[0]+state.apothem_m*3.0,state.centre_xy_m[1],1.],
    ])
    values=result.traces(points)[:,-1]
    assert values[0]>0.0
    assert values[1]==pytest.approx(0.0)
    assert values[2]==pytest.approx(0.0)
    assert result.diagnostics['front_profile']=='wake_tail'


def test_vertical_tail_profile_is_above_cloud_only_and_zero_at_extension_end(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=2.,wind_xy_m_s=(3.9,0.),
        front_profile='vertical_tail',vertical_tail_scale_m=1.0,
        vertical_tail_extension_m=4.0)
    state=result.states[-1]
    points=np.array([
        [state.centre_xy_m[0],state.centre_xy_m[1],state.height_m+1.0],
        [state.centre_xy_m[0],state.centre_xy_m[1],state.height_m+4.1],
        [60.,320.,6.5],
    ])
    values=result.traces(points)
    assert values[0,-1]>0.0
    assert values[1,-1]==pytest.approx(0.0)
    assert values[2,0]==pytest.approx(32.0)
    assert result.diagnostics['front_profile']=='vertical_tail'
    assert result.diagnostics['vertical_tail_scale_m']==pytest.approx(1.0)


def test_wind_anisotropic_leading_profile_contracts_only_downwind_front(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=2.,wind_xy_m_s=(3.9,0.),
        front_profile='wind_anisotropic_leading',leading_edge_scale=.5)
    state=result.states[-1]
    points=np.array([
        [state.centre_xy_m[0]+state.apothem_m*.75,state.centre_xy_m[1],1.],
        [state.centre_xy_m[0],state.centre_xy_m[1],1.],
    ])
    values=result.traces(points)[:,-1]
    assert values[0]==pytest.approx(0.)
    assert values[1]>0.
    assert result.diagnostics['front_profile']=='wind_anisotropic_leading'
    assert result.diagnostics['leading_edge_scale']==pytest.approx(.5)
    with pytest.raises(ValueError,match='leading_edge_scale'):
        integrate_inertial_prism(cloud,closure,front_profile='wind_anisotropic_leading',
            leading_edge_scale=0.)


def test_wind_anisotropic_leading_onset_preserves_near_core(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=2.,wind_xy_m_s=(3.9,0.),
        front_profile='wind_anisotropic_leading',leading_edge_scale=.5,
        leading_edge_onset_fraction=.75)
    state=result.states[-1]
    points=np.array([
        [state.centre_xy_m[0]+state.apothem_m*.5,state.centre_xy_m[1],1.],
        [state.centre_xy_m[0]+state.apothem_m*.9,state.centre_xy_m[1],1.],
    ])
    values=result.traces(points)[:,-1]
    assert values[0]>0.
    assert values[1]==pytest.approx(0.)
    assert result.diagnostics['leading_edge_onset_fraction']==pytest.approx(.75)
    with pytest.raises(ValueError,match='leading_edge_onset_fraction'):
        integrate_inertial_prism(cloud,closure,front_profile='wind_anisotropic_leading',
            leading_edge_onset_fraction=1.1)


def test_entrained_edge_profile_starts_uniform_then_conserves_species_at_full_handoff(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,front_profile='entrained_edge')
    points=np.array([[60.,320.,6.5],[66.,320.,6.5],[74.,334.,.4]])
    assert result.traces(points)[:,0]==pytest.approx([32.,32.,0.])
    assert result.traces(np.array([[74.,334.,.4]]))[:,-1][0]<32.
    xe=np.arange(-10.,160.1,2.5);ye=np.arange(200.,440.1,2.5)
    q,ledger=prism_layer_inventory(result,-1,xe,ye)
    state=result.states[-1]
    assert q[0].sum()*2.5**2==pytest.approx(state.total_mass_kg,rel=1e-12)
    assert q[1].sum()*2.5**2==pytest.approx(state.species_mass_kg,rel=1e-12)
    assert ledger['front_profile']=='entrained_edge'
    assert ledger['front_profile_species_normalization']>1.


def test_grounded_head_diagnostic_splits_heights_but_is_not_silently_mass_uniform(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,front_profile='grounded_head',
        grounded_head_fraction=.38)
    centre_heights=np.array([[60.,320.,.4],[60.,320.,2.4]])
    assert result.traces(centre_heights)[:,0]==pytest.approx([32.,32.])
    same_xy=np.array([[74.,334.,.4],[74.,334.,2.4]])
    assert result.traces(same_xy)[:,-1][0]>0.
    assert result.traces(same_xy)[:,-1][1]==pytest.approx(0.)
    with pytest.raises(ValueError,match='grounded_head_fraction'):
        integrate_inertial_prism(cloud,closure,grounded_head_fraction=0.)


def test_diffuse_edge_profile_extends_support_without_changing_initial_state(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='diffuse_edge',edge_extension_fraction=.5)
    points=np.array([[60.,320.,1.],[80.,320.,1.]])
    assert result.traces(points)[:,0]==pytest.approx([32.,0.])
    final=result.traces(points)[:,-1]
    assert final[0]>0. and final[1]>0.
    assert result.diagnostics['front_profile']=='diffuse_edge'
    with pytest.raises(ValueError,match='edge_extension_fraction'):
        integrate_inertial_prism(cloud,closure,front_profile='diffuse_edge',
            edge_extension_fraction=-.1)


def test_smooth_edge_profile_only_softens_the_front_diagnostic(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='smooth_edge',edge_smoothing_fraction=.5)
    points=np.array([[60.,320.,6.5],[130.,320.,1.]])
    assert result.traces(points)[:,0]==pytest.approx([32.,0.])
    uniform=integrate_inertial_prism(cloud,closure,duration_s=10.)
    assert result.traces(points)[:,-1][1]<uniform.traces(points)[:,-1][1]
    with pytest.raises(ValueError,match='edge_smoothing_fraction'):
        integrate_inertial_prism(cloud,closure,front_profile='smooth_edge',
            edge_smoothing_fraction=0.)


def test_longitudinal_ramp_softens_only_the_leading_face(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='longitudinal_ramp',longitudinal_ramp_fraction=.5)
    points=np.array([[60.,320.,6.5],[120.,320.,1.]])
    assert result.traces(points)[:,0]==pytest.approx([32.,0.])
    uniform=integrate_inertial_prism(cloud,closure,duration_s=10.)
    assert result.traces(points)[:,-1][1]<=uniform.traces(points)[:,-1][1]
    with pytest.raises(ValueError,match='longitudinal_ramp_fraction'):
        integrate_inertial_prism(cloud,closure,front_profile='longitudinal_ramp',
            longitudinal_ramp_fraction=0.)


def test_longitudinal_ramp_anisotropic_combines_front_and_crosswind_support(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='longitudinal_ramp_anisotropic',
        longitudinal_ramp_fraction=.5,crosswind_scale=.5)
    points=np.array([[60.,320.,6.5],[60.,335.,1.]])
    uniform=integrate_inertial_prism(cloud,closure,duration_s=10.)
    assert result.traces(points)[:,0]==pytest.approx([32.,0.])
    assert result.traces(points)[:,-1][1]<=uniform.traces(points)[:,-1][1]
    assert result.diagnostics['front_profile']=='longitudinal_ramp_anisotropic'


def test_longitudinal_ramp_gaussian_attenuates_crosswind_without_hard_cutoff(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='longitudinal_ramp_gaussian',
        longitudinal_ramp_fraction=.5,crosswind_scale=.5,
        crosswind_drift_ratio=.2)
    points=np.array([[60.,320.,6.5],[60.,335.,1.]])
    uniform=integrate_inertial_prism(cloud,closure,duration_s=10.)
    assert result.traces(points)[:,0]==pytest.approx([32.,0.])
    assert 0.<result.traces(points)[:,-1][1]<uniform.traces(points)[:,-1][1]
    assert result.diagnostics['front_profile']=='longitudinal_ramp_gaussian'
    assert result.diagnostics['crosswind_drift_ratio']==pytest.approx(.2)
    with pytest.raises(ValueError,match='crosswind_drift_ratio'):
        integrate_inertial_prism(cloud,closure,front_profile='longitudinal_ramp_gaussian',
            crosswind_drift_ratio=-.1)


def test_gaussian_source_profile_can_expose_diagnostic_vertical_tail(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='longitudinal_ramp_gaussian',crosswind_scale=.22,
        vertical_observation_profile='exponential_tail',
        vertical_observation_tail_scale_m=8.,
        vertical_observation_tail_extension_m=16.)
    value=result.traces(np.array([[60.,320.,4.]]))[:,-1][0]
    assert value>0.
    assert result.diagnostics['vertical_observation_profile']=='exponential_tail'
    assert result.diagnostics['vertical_observation_tail_scale_m']==pytest.approx(8.)


def test_gaussian_source_diffusive_tail_is_delayed_and_observation_only(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=10.,front_profile='longitudinal_ramp_gaussian',
        crosswind_scale=.22,vertical_observation_profile='diffusive_tail',
        vertical_observation_tail_extension_m=16.,
        vertical_observation_tail_diffusivity_m2_s=1.25)
    point=np.array([[result.states[-1].centre_xy_m[0],
                     result.states[-1].centre_xy_m[1],
                     result.states[-1].height_m+2.]])
    # At release there is no instantaneous signal above the compact cloud;
    # the diagnostic tail appears only after finite elapsed time.
    assert result.traces(np.array([[60.,320.,14.]]))[0,0]==pytest.approx(0.)
    assert result.traces(point)[0,-1]>0.
    assert result.diagnostics['vertical_observation_profile']=='diffusive_tail'
    assert result.diagnostics['vertical_observation_tail_diffusivity_m2_s']==pytest.approx(1.25)
    assert result.states[-1].species_mass_kg==pytest.approx(result.states[0].species_mass_kg,
                                                              rel=1e-12)


def test_gaussian_source_diffusive_tail_can_end_at_declared_transient_time(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=10.,front_profile='longitudinal_ramp_gaussian',
        crosswind_scale=.22,vertical_observation_profile='diffusive_tail',
        vertical_observation_tail_extension_m=16.,
        vertical_observation_tail_diffusivity_m2_s=1.25,
        vertical_observation_tail_until_s=5.)
    point=np.array([[result.states[-1].centre_xy_m[0],
                     result.states[-1].centre_xy_m[1],
                     result.states[-1].height_m+2.]])
    assert result.traces(point)[0,-1]==pytest.approx(0.)
    assert result.diagnostics['vertical_observation_tail_until_s']==pytest.approx(5.)
    with pytest.raises(ValueError,match='vertical_observation_tail_until_s'):
        integrate_inertial_prism(cloud,closure,vertical_observation_tail_until_s=-1.)


def test_gaussian_source_vertical_crosswind_spread_only_broadens_upper_observation(source):
    cloud,closure=source
    kwargs=dict(duration_s=10.,front_profile='longitudinal_ramp_gaussian',
                crosswind_scale=.22,vertical_observation_profile='diffusive_tail',
                vertical_observation_tail_extension_m=8.,
                vertical_observation_tail_diffusivity_m2_s=1.25)
    narrow=integrate_inertial_prism(cloud,closure,**kwargs)
    broad=integrate_inertial_prism(cloud,closure,
        vertical_crosswind_spread_m_per_m=1.0,**kwargs)
    s=broad.states[-1]
    point=np.array([[s.centre_xy_m[0],s.centre_xy_m[1]+s.apothem_m*.9,
                     s.height_m+2.]])
    assert broad.traces(point)[0,-1]>=narrow.traces(point)[0,-1]
    assert broad.states[-1].species_mass_kg==pytest.approx(
        narrow.states[-1].species_mass_kg,rel=1e-12)
    assert broad.diagnostics['vertical_crosswind_spread_m_per_m']==pytest.approx(1.)
    with pytest.raises(ValueError,match='vertical_crosswind_spread_m_per_m'):
        integrate_inertial_prism(cloud,closure,vertical_crosswind_spread_m_per_m=-.1)


def test_staged_source_release_is_conservative_and_delays_radial_spread(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=4.,dt_s=.01,release_duration_s=2.,
        entrainment_shear=0.,entrainment_ambient=0.,drag_coefficient=0.)
    assert result.diagnostics['source_release_staged'] is True
    assert result.diagnostics['release_duration_s']==pytest.approx(2.)
    assert result.diagnostics['released_source_mass_kg']==pytest.approx(
        result.diagnostics['release_target_mass_kg'],rel=1e-10)
    assert result.states[0].total_mass_kg < result.states[-1].total_mass_kg
    # No radial gravity spread is spent while the source is still being formed.
    assert result.states[1].radial_velocity_m_s==pytest.approx(0.)
    assert result.states[-1].radial_velocity_m_s>0.
    assert result.states[-1].species_mass_kg==pytest.approx(
        result.diagnostics['release_target_species_mass_kg'],rel=1e-10)
    with pytest.raises(ValueError,match='release_duration_s'):
        integrate_inertial_prism(cloud,closure,release_duration_s=0.)
    with pytest.raises(ValueError,match='release_injection_velocity_xy_m_s'):
        integrate_inertial_prism(cloud,closure,release_duration_s=2.,
                                 release_injection_velocity_xy_m_s=(1.,))


def test_staged_release_can_separate_source_and_ambient_wind_momentum(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=1.,dt_s=.01,release_duration_s=1.,
        release_injection_velocity_xy_m_s=(0.,0.),
        entrainment_shear=0.,entrainment_ambient=0.,drag_coefficient=0.)
    assert result.diagnostics['release_injection_velocity_xy_m_s']==pytest.approx((0.,0.))
    assert result.states[-1].bulk_velocity_xy_m_s==pytest.approx((0.,0.),abs=1e-12)
    assert result.states[-1].species_mass_kg==pytest.approx(
        result.diagnostics['release_target_species_mass_kg'],rel=1e-10)


def test_source_wind_scale_is_explicit_and_conservative(source):
    cloud,closure=source
    ambient=integrate_inertial_prism(cloud,closure,duration_s=1.,wind_xy_m_s=(4.,0.))
    reduced=integrate_inertial_prism(cloud,closure,duration_s=1.,wind_xy_m_s=(4.,0.),
        source_wind_scale=.75)
    assert reduced.states[-1].centre_xy_m[0] < ambient.states[-1].centre_xy_m[0]
    assert reduced.states[-1].species_mass_kg == pytest.approx(
        ambient.states[-1].species_mass_kg,rel=1e-12)
    assert reduced.diagnostics['maximum_relative_species_residual'] == pytest.approx(0.)
    assert reduced.diagnostics['mixture_mass_residual_kg'] == pytest.approx(0.,abs=1e-8)
    assert reduced.diagnostics['source_wind_scale']==pytest.approx(.75)
    assert reduced.diagnostics['source_wind_xy_m_s']==pytest.approx((3.,0.))
    assert reduced.diagnostics['ambient_wind_xy_m_s']==pytest.approx((4.,0.))
    with pytest.raises(ValueError,match='source_wind_scale'):
        integrate_inertial_prism(cloud,closure,source_wind_scale=0.)


def test_gaussian_observation_tail_extends_crosswind_support_without_handoff_mass_change(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=2.,wind_xy_m_s=(3.9,0.),
        front_profile='longitudinal_ramp_gaussian',crosswind_scale=.5,
        crosswind_observation_tail_sigma=2.)
    state=result.states[-1]
    core=integrate_inertial_prism(
        cloud,closure,duration_s=2.,wind_xy_m_s=(3.9,0.),
        front_profile='longitudinal_ramp_gaussian',crosswind_scale=.5)
    point=np.array([[state.centre_xy_m[0],state.centre_xy_m[1]+
                     1.5*state.apothem_m,1.]])
    assert core.traces(point)[0,-1]==pytest.approx(0.)
    assert result.traces(point)[0,-1]>0.
    assert result.diagnostics['crosswind_observation_tail_sigma']==pytest.approx(2.)
    with pytest.raises(ValueError,match='crosswind_observation_tail_sigma'):
        integrate_inertial_prism(cloud,closure,crosswind_observation_tail_sigma=-.1)


def test_anisotropic_support_narrows_only_crosswind_diagnostic(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,
        front_profile='anisotropic_support',crosswind_scale=.5)
    points=np.array([[60.,320.,6.5],[60.,370.,1.]])
    uniform=integrate_inertial_prism(cloud,closure,duration_s=10.)
    assert result.traces(points)[:,0]==pytest.approx([32.,0.])
    assert result.traces(points)[:,-1][1]<=uniform.traces(points)[:,-1][1]
    with pytest.raises(ValueError,match='crosswind_scale'):
        integrate_inertial_prism(cloud,closure,front_profile='anisotropic_support',
            crosswind_scale=0.)


def test_radial_froude_limit_records_dissipation_and_bounds_front_speed(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,radial_froude_limit=.5)
    state=result.states[-1]
    reduced_g=max((state.total_mass_kg-result.thermodynamics.ambient_density*state.total_mass_kg/
                   state.density_kg_m3)/state.total_mass_kg,0.)*9.80665
    assert state.radial_velocity_m_s <= .5*math.sqrt(reduced_g*state.height_m)+1e-10
    assert result.diagnostics['radial_dissipation_j']>0.
    assert result.diagnostics['maximum_relative_radial_energy_ledger_residual']<1e-7
    with pytest.raises(ValueError,match='radial_froude_limit'):
        integrate_inertial_prism(cloud,closure,radial_froude_limit=0.)


def test_radial_damping_is_energy_audited_and_reduces_front_speed(source):
    cloud,closure=source
    undamped=integrate_inertial_prism(cloud,closure,duration_s=10.)
    damped=integrate_inertial_prism(cloud,closure,duration_s=10.,
        radial_damping_coefficient=.1)
    assert damped.states[-1].radial_velocity_m_s<undamped.states[-1].radial_velocity_m_s
    assert damped.diagnostics['radial_dissipation_j']>0.
    assert damped.diagnostics['maximum_relative_radial_energy_ledger_residual']<1e-7
    with pytest.raises(ValueError,match='radial_damping_coefficient'):
        integrate_inertial_prism(cloud,closure,radial_damping_coefficient=-.1)


def test_aspect_ratio_damping_activates_as_height_collapses_and_is_audited(source):
    cloud,closure=source
    neutral=integrate_inertial_prism(
        cloud,closure,duration_s=10.,radial_damping_coefficient=1.,
        radial_damping_aspect_ratio=0.)
    aspect=integrate_inertial_prism(
        cloud,closure,duration_s=10.,radial_damping_coefficient=1.,
        radial_damping_aspect_ratio=.5,radial_damping_aspect_threshold=0.)
    assert aspect.states[-1].radial_velocity_m_s < neutral.states[-1].radial_velocity_m_s
    assert aspect.states[-1].apothem_m < neutral.states[-1].apothem_m
    assert aspect.diagnostics['radial_damping_aspect_ratio']==pytest.approx(.5)
    assert aspect.diagnostics['radial_damping_aspect_threshold']==pytest.approx(0.)
    assert aspect.diagnostics['radial_dissipation_j'] > 0.
    with pytest.raises(ValueError,match='radial_damping_aspect_ratio'):
        integrate_inertial_prism(cloud,closure,radial_damping_aspect_ratio=-.1)
    with pytest.raises(ValueError,match='radial_damping_aspect_threshold'):
        integrate_inertial_prism(cloud,closure,radial_damping_aspect_threshold=-.1)


def test_minimum_height_floor_preserves_volume_and_records_clamp_loss(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=10.,minimum_height_m=6.)
    assert result.states[-1].height_m>=6.-1e-10
    assert result.diagnostics['minimum_height_m']==6.
    assert result.diagnostics['radial_dissipation_j']>0.
    assert result.diagnostics['maximum_relative_radial_energy_ledger_residual']<1e-7
    with pytest.raises(ValueError,match='minimum_height_m'):
        integrate_inertial_prism(cloud,closure,minimum_height_m=14.)


def test_light_or_nonisothermal_source_not_silently_accepted(source):
    cloud,closure=source
    with pytest.raises(ValueError,match='dense species'):
        integrate_inertial_prism(cloud,replace(closure,species_molar_mass_kg_mol=.002016))
    with pytest.raises(ValueError,match='ambient source temperature'):
        integrate_inertial_prism(replace(cloud,temperature_k=280.),closure)


def test_interface_reports_unresolved_energy_instead_of_claiming_full_conservation(source):
    cloud,closure=source
    result=integrate_inertial_prism(cloud,closure,duration_s=1.)
    _,coarse=prism_layer_inventory(result,-1,np.arange(20.,102.5,2.5),np.arange(280.,362.5,2.5))
    _,fine=prism_layer_inventory(result,-1,np.arange(20.,101.25,1.25),np.arange(280.,361.25,1.25))
    assert coarse['subcell_horizontal_kinetic_energy_not_transferred_j']>0
    assert fine['subcell_horizontal_kinetic_energy_not_transferred_j']<coarse['subcell_horizontal_kinetic_energy_not_transferred_j']
    assert coarse['full_source_vertical_affine_spreading_kinetic_energy_not_transferred_j']>0
    assert coarse['full_energy_conserving_interface'] is False


def test_matching_front_handoff_preserves_mixture_and_species_inventory(source):
    cloud,closure=source
    result=integrate_inertial_prism(
        cloud,closure,duration_s=2.,front_profile='longitudinal_ramp_gaussian',
        longitudinal_ramp_fraction=.8,crosswind_scale=.22,
        handoff_profile='match_front_profile')
    q,audit=prism_layer_inventory(
        result,-1,np.arange(20.,102.5,2.5),np.arange(280.,362.5,2.5))
    state=result.states[-1]
    assert audit['handoff_profile']=='match_front_profile'
    assert audit['mapped_conserved'][0]==pytest.approx(state.total_mass_kg,rel=1e-10)
    assert audit['mapped_conserved'][1]==pytest.approx(state.species_mass_kg,rel=1e-10)
    assert np.all(q[1] <= q[0] + 1e-10)
    assert audit['exterior_conserved'][0]==pytest.approx(0.,abs=1e-8)
    assert audit['exterior_conserved'][1]==pytest.approx(0.,abs=1e-8)
    assert audit['front_profile_species_residual_kg']==pytest.approx(0.,abs=1e-8)
    with pytest.raises(ValueError,match='handoff_profile'):
        integrate_inertial_prism(cloud,closure,handoff_profile='unknown')
