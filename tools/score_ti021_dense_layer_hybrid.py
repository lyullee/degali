"""Dense integral -> local hydrostatic layer -> grid-free exported packets.

Experimental isothermal TI021 hybrid.  The historical R12-component observable
is retained only for exact regression; the field report's GAS/HGAS readings are
scored as released-mixture concentration.  This does not qualify hydrogen or a
resolved 3-D wake.  No sensor-specific dose factor is applied.
"""
from pathlib import Path
import argparse
import json
import math
import time
import hashlib
from dataclasses import replace
import numpy as np
from scipy.special import ndtr
from scipy.ndimage import gaussian_filter, map_coordinates
from degali.addons.binary_integral_puff import initial_prism_puff,binary_puff_traces,BinaryPuffThermodynamics
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_puff import GaussianPuffConfig,integrate_finite_gaussian_puff
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.dense_layer_patch import solve_dense_layer_patch
from degali.addons.hdf5_dense_layer_interface import Hdf5DenseLayerPacketSource
from degali.addons.dense_prism_integral import integrate_inertial_prism, prism_layer_inventory
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

ROOT=Path(__file__).resolve().parents[1]


def _handoff_vertical_geometry(source_model, state, puff):
    """Return the vertical distribution that exists at layer handoff.

    The near-field model has already advanced the source before this point.
    Reusing the initial prism bounds would undo its computed collapse and
    double-count the source initialization, so both source paths use their
    current handoff state.
    """
    if source_model == 'gaussian':
        return (
            float(state.centre_position_m[2]),
            float(state.vertical_half_width_m / puff.config.width_cutoff_sigma),
        ), None
    height = float(state.height_m)
    return (0.5 * height, height / math.sqrt(12.0)), height


def curved_fence_faces(xe,ye,centre=(60.,320.),radius=50.,height=5.):
    wx=np.zeros((len(ye)-1,len(xe)));wy=np.zeros((len(ye),len(xe)-1))
    for i,x in enumerate(xe):
        a=x-centre[0]
        if not 0<=a<=radius:continue
        for y in (centre[1]-math.sqrt(max(radius**2-a*a,0)),centre[1]+math.sqrt(max(radius**2-a*a,0))):
            wx[(ye[:-1]<=y)&(ye[1:]>=y),i]=height
    for j,y in enumerate(ye):
        a=y-centre[1]
        if abs(a)>radius:continue
        x=centre[0]+math.sqrt(max(radius**2-a*a,0))
        wy[j,(xe[:-1]<=x)&(xe[1:]>=x)]=height
    return wx,wy


def exterior_trace(result,point,wind,h,handoff_time,transport='ambient_wind',
                   density_current_strength=1.0,
                   exterior_lateral_diffusivity_factor=1.0):
    """Continue outward face packets to an exterior receptor.

    ``ambient_wind`` is the historical one-way continuation.  The optional
    ``density_current`` mode adds the conserved boundary gravity-current
    speed to the down-range component of each packet and uses its resolved
    hydrostatic height for the associated lateral spreading.  It is still a
    reduced exterior closure, not a recoupled 3-D outer grid: no packet mass
    is multiplied or renormalized for the observation.
    """
    if transport not in {'ambient_wind','density_current'}:
        raise ValueError('unsupported exterior transport')
    if not math.isfinite(density_current_strength) or not 0.0 <= density_current_strength <= 1.0:
        raise ValueError('density_current_strength must be finite in [0,1]')
    if (not math.isfinite(exterior_lateral_diffusivity_factor)
            or exterior_lateral_diffusivity_factor <= 0.0):
        raise ValueError('exterior_lateral_diffusivity_factor must be finite and positive')
    t=result.time_s; rates=result.boundary_species_rate_kg_s
    boundary=result.boundary_positions_m
    height_all=result.boundary_height_m
    speed_all=result.boundary_current_speed_m_s
    if transport=='ambient_wind':
        # Retain the historical x-directed packet route exactly for the
        # reference comparison.
        distance=point[0]-boundary[:,0]
        valid=distance>0
        if wind[0] <= 0.0:
            raise ValueError('ambient-wind exterior continuation requires positive x velocity')
        carrier_speed=np.full(int(np.count_nonzero(valid)),float(wind[0]))
        delay=distance[valid]/carrier_speed
        dy=point[1]-boundary[valid,1]-wind[1]*delay
        d=exterior_lateral_diffusivity_factor*.4*.205*max(point[2],.4)/.7
        sigma=np.sqrt(2*d*delay)
        ky=(ndtr((h/2-dy)/sigma)-ndtr((-h/2-dy)/sigma))/h
        height=height_all[:,valid]
    else:
        # A packet may leave any face and spread as a shallow density current.
        # The resolved hydrostatic gravity-current speed is projected along
        # the receptor ray; this preserves a positive travel-time gate and
        # avoids inventing an upstream source when the projected carrier flow
        # cannot reach the receptor.
        delta=np.asarray(point[:2],float)[None,:]-boundary[:,:2]
        distance=np.linalg.norm(delta,axis=1)
        unit=np.divide(delta,distance[:,None],out=np.zeros_like(delta),where=distance[:,None]>1e-30)
        wind_projected=unit@np.asarray(wind,float)
        # Boundary current speed is recorded at every layer step.  Use the
        # time-average per face so each packet retains the same one-dimensional
        # delay representation as the historical exterior operator.
        speed_face=np.nanmean(speed_all,axis=0)
        effective=wind_projected+density_current_strength*speed_face
        valid=(distance>0.)&(effective>1e-9)
        carrier_speed=effective[valid]
        delay=np.divide(distance[valid],carrier_speed)
        normal=np.column_stack((-unit[valid,1],unit[valid,0]))
        cross=np.sum(normal*(delta[valid]-np.asarray(wind)[None,:]*delay[:,None]),axis=1)
        h_eff=np.maximum(height_all[:,valid],1e-6)
        d=(exterior_lateral_diffusivity_factor*.4*.205*max(point[2],.4)/.7
           + .5*speed_face[valid]*h_eff)
        sigma=np.sqrt(np.maximum(2*d*delay[None,:],1e-30))
        ky=(ndtr((h/2-cross[None,:])/sigma)-ndtr((-h/2-cross[None,:])/sigma))/h
        height=h_eff
    kz=np.divide(ndtr((height-point[2])/sigma)-ndtr(-point[2]/sigma)
                 +ndtr((height+point[2])/sigma)-ndtr(point[2]/sigma),height,
                 out=np.zeros_like(height),where=height>1e-20)
    q=rates[:,valid]*ky*kz/carrier_speed
    trace=np.zeros(len(t))
    for k,lag in enumerate(delay):
        idx=np.searchsorted(t,t-lag,side='right')-1
        ok=(idx>=0)&(idx<len(q))
        trace[ok]+=q[idx[ok],k]
    return trace


def run(output,*,h=5.,dt=.25,source_dt=.025,layer_wind='neutral_log',
        entrainment='constant',entrainment_coefficient=.1,
        entrainment_far_coefficient=None,entrainment_relaxation_time_s=60.0,
        drag_coefficient=.2,
        domain_mode='local',source_mass_scale=1.0,
        source_spread_scale=1.0,
        source_longitudinal_scale=1.0,
        source_crosswind_scale=1.0,
        source_initial_velocity_scale=0.0,
        source_release_duration_s=None,
        puff_shear_entrainment_coefficient=.1,
        puff_ambient_turbulence_coefficient=.1,
        puff_top_shear_entrainment_coefficient=None,
        puff_top_ambient_turbulence_coefficient=None,
        puff_side_shear_entrainment_coefficient=None,
        puff_side_ambient_turbulence_coefficient=None,
        puff_gravity_spreading_coefficient=.75,
        puff_form_drag_coefficient=.2,
        handoff_radius_m=50.0,
        vertical_observation_profile='top_hat',dense_vertical_exponent=0.0,vertical_tail_extension_m=0.0,
        vertical_grounding_exponent=0.0,
        vertical_elevated_fraction=0.0,
        wall_overtop_exponent=0.0,
        wall_blockage_fraction=1.0,
        vertical_tail_diffusivity_m2_s=1.0,
        mass_conserving_vertical_profile=False,vertical_layers=1,
        vertical_observation_mode='conservative_volume',vertical_observation_bins=64,
        vertical_observation_component_diagnostics=False,
        initial_vertical_gaussian=True,
        vertical_mixing_diffusivity_m2_s=0.0,height_resolved_wall_flux=False,
         multilayer_entrainment_placement='top',
        multilayer_horizontal_shear_fraction=0.0,
        multilayer_obstacle_crosswind_shear_rate_s_inv=0.0,
        multilayer_obstacle_crosswind_decay_length_m=70.0,
        multilayer_relative_transport_interval_s=0.0,
        obstacle_vertical_lift_fraction=0.0,
        obstacle_vertical_lift_decay_m=30.0,
        wall_impulse_vertical_lift_fraction=0.0,
        wall_impulse_elevated_transfer_fraction=0.0,
        wall_impulse_elevated_state_mode='species_marker',
        wall_impulse_elevated_conserved_fraction=1.0,
        wall_flux_elevated_transfer_fraction=0.0,
        wall_impulse_elevated_transport_fraction=0.0,
        wall_impulse_elevated_crosswind_velocity_m_s=0.0,
        wall_impulse_elevated_crosswind_diffusivity_m2_s=0.0,
        wall_impulse_elevated_crosswind_outward_velocity_m_s=0.0,
        wall_impulse_elevated_crosswind_outward_scale_m=50.0,
        wall_impulse_elevated_crosswind_outward_decay_s=0.0,
        wall_impulse_elevated_crosswind_spinup_s=0.0,
        wall_impulse_elevated_transport_interval_s=0.0,
        wall_impulse_elevated_advection_scale=1.0,
        wall_impulse_elevated_advection_duration_s=0.0,
        wall_impulse_elevated_return_rate_s_inv=0.0,
        wall_impulse_elevated_center_m=6.0,
        wall_impulse_elevated_sigma_m=2.0,
        wall_impulse_elevated_observation_ramp_m=0.0,
        wall_impulse_elevated_settling_velocity_m_s=0.0,
        wall_impulse_elevated_vertical_diffusivity_m2_s=0.0,
        wall_impulse_elevated_min_center_m=0.0,
        wall_impulse_elevated_rise_velocity_m_s=0.0,
        wall_impulse_elevated_rise_time_s=1.0,
        wall_impulse_elevated_conserved_center_m=None,
        wall_impulse_elevated_conserved_sigma_m=None,
        wall_impulse_elevated_conserved_settling_velocity_m_s=None,
        wall_impulse_elevated_conserved_vertical_diffusivity_m2_s=None,
        wall_impulse_elevated_conserved_min_center_m=None,
        wall_impulse_elevated_conserved_rise_velocity_m_s=0.0,
        wall_impulse_elevated_conserved_rise_time_s=1.0,
        wall_impulse_elevated_conserved_transport_fraction=None,
        wall_impulse_elevated_conserved_crosswind_velocity_m_s=None,
        wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s=None,
        wall_impulse_elevated_conserved_crosswind_outward_scale_m=None,
        wall_impulse_elevated_conserved_center_mode='fixed',
        wall_impulse_reattached_fraction=0.0,
        wall_impulse_reattached_transport_fraction=0.0,
        wall_impulse_reattached_transport_height_m=None,
        wall_impulse_reattached_crosswind_velocity_m_s=0.0,
        wall_impulse_reattached_crosswind_outward_velocity_m_s=0.0,
        wall_impulse_reattached_crosswind_outward_scale_m=50.0,
        wall_impulse_reattached_return_rate_s_inv=0.0,
        wall_impulse_reattached_center_m=0.7,
        wall_impulse_reattached_sigma_m=0.35,
        wall_impulse_reattached_settling_velocity_m_s=0.002,
        wall_impulse_reattached_vertical_diffusivity_m2_s=0.0,
        wall_impulse_reattached_min_center_m=0.0,
        wall_impulse_wake_capture_rate_s_inv=0.0,
        wall_impulse_wake_capture_mode='bulk',
        wall_impulse_wake_return_rate_s_inv=0.0,
        wall_impulse_wake_transport_fraction=0.0,
        wall_impulse_wake_crosswind_velocity_m_s=0.0,
        wall_impulse_wake_crosswind_outward_velocity_m_s=0.0,
        wall_impulse_wake_crosswind_outward_scale_m=50.0,
        wall_impulse_wake_observation_fraction=1.0,
        wall_impulse_wake_observation_center_m=None,
        wall_impulse_wake_observation_sigma_m=2.0,
        wall_impulse_wake_settling_velocity_m_s=0.0,
        wall_impulse_wake_vertical_diffusivity_m2_s=0.0,
        wall_impulse_wake_min_center_m=0.0,
        initial_elevated_species_fraction=0.0,
        dense_settling_velocity_m_s=0.0,
        receptor_interpolation='nearest',
        wall_tangential_bypass_fraction=0.0,
        wall_tip_bypass_fraction=0.0,wall_tip_bypass_power=2.0,
        wall_overflow_capture_fraction=0.0,
        wall_overflow_return_rate_s_inv=0.0,
        wall_overflow_side_gate_length_m=0.0,
        horizontal_diffusivity_m2_s=0.0,
        horizontal_diffusion_respect_walls=False,
        crosswind_diffusivity_factor=1.0,
        wake_deficit_strength=0.0,wake_length_m=1.0,wake_base_width_m=1.0,
        wake_spread_slope=0.0,wake_recirculation_strength=0.0,
        wake_meander_strength=0.0,wake_shedding_strouhal=0.2,
        wake_meander_phase_deg=0.0,
       core_meander_strength=0.0,core_meander_length_m=100.0,
       core_meander_width_m=50.0,
        cavity_shelter_strength=0.0,
       cavity_entrainment_factor=1.0,
       cavity_core_entrainment_reduction=0.0,
       cavity_vortex_rate_s_inv=0.0,cavity_vortex_strength=0.0,
       cavity_reservoir_capture_rate_s_inv=0.005,cavity_reservoir_return_rate_s_inv=0.005,
       cavity_reservoir_return_target='main',
       cavity_reservoir_reattached_fraction=1.0,
       cavity_reservoir_capture_mode='cavity_rate',
       cavity_reservoir_vertical_capture_mode='column',
       cavity_reservoir_spatial_mode='well_mixed',
       cavity_reservoir_content_mode='conserved',
       cavity_reservoir_transport_fraction=0.0,
       cavity_reservoir_transport_height_m=None,
       cavity_reservoir_transport_mask='cavity',
       cavity_reservoir_wake_length_m=0.0,
       cavity_reservoir_wake_half_width_m=None,
       cavity_reservoir_crosswind_convergence_rate_s_inv=0.0,
       cavity_reservoir_crosswind_target_offset_m=0.0,
       local_reservoir_observation_fraction=0.0,
       local_reservoir_observation_center_m=None,
       local_reservoir_observation_sigma_m=2.0,
       cavity_reservoir_return_sigma_m=16.7,
       cavity_reservoir_return_radius_m=10.0,
       cavity_reservoir_return_offset_m=0.0,
       cavity_reservoir_return_crosswind_scale=1.0,
       cavity_reservoir_face_return_fraction=0.0,
       cavity_reservoir_capture_radius_m=0.0,
       cavity_reservoir_start_time_s=0.0,
       cavity_reservoir_return_delay_s=0.0,
       cavity_return_streamwise_bias=0.0,
       cavity_vertical_exponent=0.0,
       cavity_vertical_radial_decay_power=0.0,
       cavity_vertical_mixing_factor=0.0,
       exterior_entrainment_far_factor=1.0,
       exterior_entrainment_decay_length_m=100.0,
       cavity_backflow_strength=0.9,
       cavity_backflow_rate_s_inv=0.01,
       cavity_backflow_frame='cavity',
       cavity_stagnation_strength=0.0,
       exterior_reentry_fraction=0.0,
       exterior_reentry_delay_s=0.0,
       exterior_reentry_sigma_m=15.0,
       cavity_geometry='downwind_semidisk',
       cavity_radius_m=50.0,
       cavity_align_with_wind=False,
        exterior_lateral_diffusivity_factor=1.0,
        exterior_transport='ambient_wind',density_current_strength=1.0,
        handoff_interface_path=None,interface_x_m=None,
        wind_history_path=None,wind_history_component_mode='measured',
        wind_history_direction_offset_deg=0.0,
        wind_history_advection_reference_x_m=None,
        wind_history_advection_speed_m_s=None,
        wind_history_scale_friction_velocity=True,
        wind_speed_m_s=3.9,layer_wind_speed_m_s=None,friction_velocity_m_s=None,
        wind_angle_deg=6.1,wall_height_m=5.0,gravity_scale=1.0,
        wall_radius_m=50.0,
        source_model='gaussian',source_front_profile='uniform',
        source_initial_radius_m=14.0,
        observation_species_mode='r12_component',
        source_radial_damping_coefficient=0.0,source_radial_damping_aspect_ratio=0.0,
        source_minimum_height_m=0.0,source_wind_scale=1.0):
    started=time.perf_counter()
    if not math.isfinite(source_mass_scale) or source_mass_scale<=0.0:
        raise ValueError('source_mass_scale must be positive and finite')
    if not isinstance(wind_history_scale_friction_velocity,bool):
        raise TypeError('wind_history_scale_friction_velocity must be boolean')
    if friction_velocity_m_s is None:
        friction_velocity_m_s=0.205
    if (not math.isfinite(float(friction_velocity_m_s))
            or float(friction_velocity_m_s)<=0.0):
        raise ValueError('friction_velocity_m_s must be positive and finite')
    if wind_history_component_mode not in {'measured','direction_only','speed_only'}:
        raise ValueError("wind_history_component_mode must be 'measured', 'direction_only', or 'speed_only'")
    if not math.isfinite(wind_history_direction_offset_deg):
        raise ValueError('wind_history_direction_offset_deg must be finite')
    wind_history_time_s=None
    wind_history_xy_m_s=None
    wind_history_sha256=None
    wind_history_vector_scale=None
    if wind_history_path is not None:
        wind_history_path=Path(wind_history_path).resolve()
        with np.load(wind_history_path) as history:
            if 'time_s' not in history or 'wind_xy_m_s' not in history:
                raise ValueError('wind history must contain time_s and wind_xy_m_s')
            wind_history_time_s=np.asarray(history['time_s'],dtype=float).copy()
            wind_history_xy_m_s=np.asarray(history['wind_xy_m_s'],dtype=float).copy()
        if (wind_history_time_s.ndim!=1 or wind_history_xy_m_s.shape!=(wind_history_time_s.size,2)
                or wind_history_time_s.size<2):
            raise ValueError('invalid wind history array shapes')
        if (not np.all(np.isfinite(wind_history_time_s))
                or not np.all(np.isfinite(wind_history_xy_m_s))
                or np.any(np.diff(wind_history_time_s)<=0.0)):
            raise ValueError('wind history must be finite and strictly increasing')
        wind_history_sha256=hashlib.sha256(wind_history_path.read_bytes()).hexdigest()
    if not math.isfinite(source_spread_scale) or source_spread_scale<=0.0:
        raise ValueError('source_spread_scale must be positive and finite')
    for value,name in ((source_longitudinal_scale,'source_longitudinal_scale'),
                       (source_crosswind_scale,'source_crosswind_scale')):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f'{name} must be positive and finite')
    if not math.isfinite(source_initial_velocity_scale) or source_initial_velocity_scale<0.0:
        raise ValueError('source_initial_velocity_scale must be finite and non-negative')
    if source_release_duration_s is not None and (
            not math.isfinite(source_release_duration_s) or source_release_duration_s <= 0.0):
        raise ValueError('source_release_duration_s must be finite and positive when provided')
    if not math.isfinite(handoff_radius_m) or handoff_radius_m <= 0.0 or handoff_radius_m > 50.0:
        raise ValueError('handoff_radius_m must be finite in (0, 50]')
    if not math.isfinite(entrainment_coefficient) or entrainment_coefficient<0.0:
        raise ValueError('entrainment_coefficient must be finite and non-negative')
    if entrainment_far_coefficient is not None and (
            not math.isfinite(entrainment_far_coefficient)
            or entrainment_far_coefficient < 0.0):
        raise ValueError(
            'entrainment_far_coefficient must be finite and non-negative')
    if (not math.isfinite(entrainment_relaxation_time_s)
            or entrainment_relaxation_time_s <= 0.0):
        raise ValueError('entrainment_relaxation_time_s must be positive and finite')
    if not math.isfinite(drag_coefficient) or drag_coefficient<0.0:
        raise ValueError('drag_coefficient must be finite and non-negative')
    for value,name in ((puff_shear_entrainment_coefficient,'puff_shear_entrainment_coefficient'),
                       (puff_ambient_turbulence_coefficient,'puff_ambient_turbulence_coefficient'),
                       (puff_gravity_spreading_coefficient,'puff_gravity_spreading_coefficient'),
                       (puff_form_drag_coefficient,'puff_form_drag_coefficient')):
        if not math.isfinite(value) or value<0.0:
            raise ValueError(f'{name} must be finite and non-negative')
    if vertical_observation_profile not in {'top_hat','dense_power','exponential_tail','diffusive_tail','handoff_gaussian','multilayer'}:
        raise ValueError('unsupported vertical_observation_profile')
    if (not math.isfinite(vertical_tail_extension_m)
            or vertical_tail_extension_m<0.0
            or not math.isfinite(vertical_tail_diffusivity_m2_s)
            or vertical_tail_diffusivity_m2_s<=0.0):
        raise ValueError('invalid vertical tail controls')
    if not isinstance(mass_conserving_vertical_profile,bool):
        raise TypeError('mass_conserving_vertical_profile must be boolean')
    if isinstance(vertical_layers,bool) or not isinstance(vertical_layers,int) or vertical_layers<1:
        raise ValueError('vertical_layers must be a positive integer')
    if not math.isfinite(vertical_mixing_diffusivity_m2_s) or vertical_mixing_diffusivity_m2_s<0:
        raise ValueError('vertical_mixing_diffusivity_m2_s must be finite and non-negative')
    if not isinstance(height_resolved_wall_flux,bool):
        raise TypeError('height_resolved_wall_flux must be boolean')
    if receptor_interpolation not in {'nearest','bilinear'}:
        raise ValueError("receptor_interpolation must be 'nearest' or 'bilinear'")
    if not math.isfinite(wall_tangential_bypass_fraction) or not 0.0 <= wall_tangential_bypass_fraction <= 1.0:
        raise ValueError('wall_tangential_bypass_fraction must lie in [0,1]')
    if not math.isfinite(horizontal_diffusivity_m2_s) or horizontal_diffusivity_m2_s<0:
        raise ValueError('horizontal_diffusivity_m2_s must be finite and non-negative')
    if not math.isfinite(crosswind_diffusivity_factor) or crosswind_diffusivity_factor <= 0.0:
        raise ValueError('crosswind_diffusivity_factor must be positive and finite')
    if (not math.isfinite(exterior_lateral_diffusivity_factor)
            or exterior_lateral_diffusivity_factor <= 0.0):
        raise ValueError('exterior_lateral_diffusivity_factor must be finite and positive')
    for value,name in ((wake_deficit_strength,'wake_deficit_strength'),
                       (wake_length_m,'wake_length_m'),
                       (wake_base_width_m,'wake_base_width_m'),
                       (wake_spread_slope,'wake_spread_slope'),
                       (wake_recirculation_strength,'wake_recirculation_strength')):
        if not math.isfinite(value) or value<0:
            raise ValueError(f'{name} must be finite and non-negative')
    if not math.isfinite(cavity_shelter_strength) or not 0.0 <= cavity_shelter_strength <= 1.0:
        raise ValueError('cavity_shelter_strength must be finite in [0,1]')
    if not math.isfinite(cavity_vertical_exponent) or cavity_vertical_exponent < 0.0:
        raise ValueError('cavity_vertical_exponent must be finite and non-negative')
    if not math.isfinite(vertical_grounding_exponent) or vertical_grounding_exponent < 0.0:
        raise ValueError('vertical_grounding_exponent must be finite and non-negative')
    if not math.isfinite(wall_overtop_exponent) or wall_overtop_exponent < 0.0:
        raise ValueError('wall_overtop_exponent must be finite and non-negative')
    if (not math.isfinite(wall_blockage_fraction)
            or not 0.0 <= wall_blockage_fraction <= 1.0):
        raise ValueError('wall_blockage_fraction must lie in [0,1]')
    if not math.isfinite(vertical_elevated_fraction) or not 0.0 <= vertical_elevated_fraction <= 1.0:
        raise ValueError('vertical_elevated_fraction must lie in [0,1]')
    if not math.isfinite(cavity_backflow_strength) or not 0.0 <= cavity_backflow_strength <= 1.0:
        raise ValueError('cavity_backflow_strength must lie in [0,1]')
    if not math.isfinite(cavity_backflow_rate_s_inv) or cavity_backflow_rate_s_inv < 0.0:
        raise ValueError('cavity_backflow_rate_s_inv must be finite and non-negative')
    if cavity_backflow_frame not in {'cavity','legacy_x'}:
        raise ValueError("cavity_backflow_frame must be 'cavity' or 'legacy_x'")
    if not math.isfinite(cavity_reservoir_return_delay_s) or cavity_reservoir_return_delay_s < 0.0:
        raise ValueError('cavity_reservoir_return_delay_s must be finite and non-negative')
    if not math.isfinite(cavity_stagnation_strength) or not 0.0 <= cavity_stagnation_strength <= 1.0:
        raise ValueError('cavity_stagnation_strength must lie in [0,1]')
    if cavity_geometry not in {'full_disk','downwind_semidisk','upwind_semidisk'}:
        raise ValueError('unsupported cavity_geometry')
    if not math.isfinite(cavity_radius_m) or cavity_radius_m <= 0.0:
        raise ValueError('cavity_radius_m must be finite and positive')
    if exterior_transport not in {'ambient_wind','density_current'}:
        raise ValueError('unsupported exterior transport')
    if not math.isfinite(density_current_strength) or not 0.0 <= density_current_strength <= 1.0:
        raise ValueError('density_current_strength must be finite in [0,1]')
    if not math.isfinite(wind_speed_m_s) or wind_speed_m_s <= 0.0:
        raise ValueError('wind_speed_m_s must be finite and positive')
    if layer_wind_speed_m_s is None:
        layer_wind_speed_m_s=float(wind_speed_m_s)
    if not math.isfinite(layer_wind_speed_m_s) or layer_wind_speed_m_s <= 0.0:
        raise ValueError('layer_wind_speed_m_s must be finite and positive')
    if not math.isfinite(wind_angle_deg):
        raise ValueError('wind_angle_deg must be finite')
    if not math.isfinite(wall_height_m) or wall_height_m <= 0.0:
        raise ValueError('wall_height_m must be finite and positive')
    if not math.isfinite(wall_radius_m) or wall_radius_m <= 0.0:
        raise ValueError('wall_radius_m must be finite and positive')
    if not math.isfinite(gravity_scale) or gravity_scale < 0.0:
        raise ValueError('gravity_scale must be finite and non-negative')
    if source_model not in {'gaussian','finite_prism'}:
        raise ValueError('unsupported source_model')
    if observation_species_mode not in {'r12_component','released_mixture'}:
        raise ValueError('unsupported observation_species_mode')
    if source_front_profile not in {'uniform','grounded_head','entrained_edge','diffuse_edge','smooth_edge','vertical_tail','wake_tail','anisotropic_support','wind_anisotropic_leading','longitudinal_ramp','longitudinal_ramp_gaussian'}:
        raise ValueError('unsupported source_front_profile')
    if not math.isfinite(source_initial_radius_m) or source_initial_radius_m <= 0.0:
        raise ValueError('source_initial_radius_m must be finite and positive')
    for value,name in ((source_radial_damping_coefficient,'source_radial_damping_coefficient'),
                       (source_radial_damping_aspect_ratio,'source_radial_damping_aspect_ratio'),
                       (source_minimum_height_m,'source_minimum_height_m'),
                       (source_wind_scale,'source_wind_scale')):
        if not math.isfinite(value) or value<0.0 or (name=='source_wind_scale' and value<=0.0):
            raise ValueError(f'{name} must be finite and valid')
    workbook=ROOT/'ti21_input_ascii.xls';dataset=read_thorney_ti021_workbook(workbook)
    if observation_species_mode=='released_mixture':
        r12_mass_fraction=(dataset.initial_r12_mole_fraction*.120913
                           /dataset.effective_molecular_weight_kg_mol)
        released_mixture_cp=(r12_mass_fraction*606.0748985077907
                             +(1.0-r12_mass_fraction)*1040.0)
        closure=BinaryIdealGasMixtureClosure(
            'R12/N2 release mixture',dataset.effective_molecular_weight_kg_mol,
            released_mixture_cp,ambient_temperature_k=291.15,
            ambient_pressure_pa=100500.)
        source_mole_fraction=1.0
        gas_label='R12/N2 release mixture'
    else:
        closure=BinaryIdealGasMixtureClosure('R12',.120913,606.0748985077907,
            ambient_temperature_k=291.15,ambient_pressure_pa=100500.)
        source_mole_fraction=dataset.initial_r12_mole_fraction
        gas_label='R12 component'
    thermo=BinaryPuffThermodynamics(closure)
    angle=np.deg2rad(wind_angle_deg)
    wind=np.array([wind_speed_m_s*np.cos(angle),wind_speed_m_s*np.sin(angle)])
    layer_wind_vec=np.array([layer_wind_speed_m_s*np.cos(angle),layer_wind_speed_m_s*np.sin(angle)])
    if wind_history_xy_m_s is not None:
        if wind_history_direction_offset_deg!=0.0:
            angle=math.radians(float(wind_history_direction_offset_deg))
            cosine=math.cos(angle);sine=math.sin(angle)
            hx=wind_history_xy_m_s[:,0].copy()
            hy=wind_history_xy_m_s[:,1].copy()
            wind_history_xy_m_s[:,0]=cosine*hx-sine*hy
            wind_history_xy_m_s[:,1]=sine*hx+cosine*hy
        history_speed=np.linalg.norm(wind_history_xy_m_s,axis=1)
        history_mean_speed=float(np.mean(history_speed))
        if not math.isfinite(history_mean_speed) or history_mean_speed<=0.0:
            raise ValueError('wind history mean speed must be positive and finite')
        if wind_history_component_mode=='direction_only':
            wind_history_xy_m_s=(wind_history_xy_m_s/history_speed[:,None])*history_mean_speed
        elif wind_history_component_mode=='speed_only':
            fixed_direction=layer_wind_vec/np.linalg.norm(layer_wind_vec)
            wind_history_xy_m_s=history_speed[:,None]*fixed_direction[None,:]
        # TI021 meteorology is the 10 m ambient record, while this reduced
        # layer model transports with a separately declared effective speed.
        # Preserve measured relative speed/direction changes without replacing
        # that independently frozen effective-speed closure.
        wind_history_vector_scale=float(layer_wind_speed_m_s/history_mean_speed)
        wind_history_xy_m_s=wind_history_xy_m_s*wind_history_vector_scale
    source_result=None
    source_cloud=None
    if source_model=='gaussian':
        puff=initial_prism_puff(RegularPolygonPrismInitialCloud((60.,320.),source_initial_radius_m,0.,13.,source_mole_fraction,291.15),
            closure,GaussianPuffConfig(duration_s=60.,time_step_s=source_dt,
                wind_velocity_m_s=(*wind,0.),width_cutoff_sigma=math.sqrt(3),
                shear_entrainment_coefficient=puff_shear_entrainment_coefficient,
                ambient_turbulence_coefficient=puff_ambient_turbulence_coefficient,
                gravity_spreading_coefficient=puff_gravity_spreading_coefficient,
                form_drag_coefficient=puff_form_drag_coefficient,
                dense_entrainment_model=entrainment),
            initial_velocity_m_s=(float(wind[0]*source_initial_velocity_scale),
                                 float(wind[1]*source_initial_velocity_scale),0.0))
    else:
        source_cloud=RegularPolygonPrismInitialCloud((60.,320.),source_initial_radius_m,0.,13.,source_mole_fraction,291.15)
        source_result=integrate_inertial_prism(
            source_cloud,
            closure, duration_s=60., dt_s=source_dt,
            wind_xy_m_s=tuple(float(v) for v in wind),
            source_wind_scale=source_wind_scale,
            entrainment_shear=puff_shear_entrainment_coefficient,
            entrainment_ambient=puff_ambient_turbulence_coefficient,
            entrainment_top_shear=puff_top_shear_entrainment_coefficient,
            entrainment_top_ambient=puff_top_ambient_turbulence_coefficient,
            entrainment_side_shear=puff_side_shear_entrainment_coefficient,
            entrainment_side_ambient=puff_side_ambient_turbulence_coefficient,
            drag_coefficient=puff_form_drag_coefficient,
            front_profile=source_front_profile,
            handoff_profile=(
                'match_front_profile'
                if source_front_profile in {
                    'longitudinal_ramp', 'longitudinal_ramp_anisotropic',
                    'longitudinal_ramp_gaussian',
                }
                else 'uniform'
            ),
            radial_damping_coefficient=source_radial_damping_coefficient,
            radial_damping_aspect_ratio=source_radial_damping_aspect_ratio,
            minimum_height_m=source_minimum_height_m,
            release_injection_velocity_xy_m_s=(
                tuple(float(v * source_initial_velocity_scale) for v in wind)
            ),
            release_duration_s=source_release_duration_s)
        puff=None
    if source_model=='gaussian' and not math.isclose(source_spread_scale,1.0):
        # Preserve the source inventory while broadening only the horizontal
        # Gaussian moment used by the short-prefix and handoff observers.  The
        # vertical source geometry remains the declared 0--13 m prism.
        widths=puff.handoff.transverse_widths_m
        handoff=replace(
            puff.handoff,
            transverse_widths_m=(widths[0]*source_spread_scale,widths[1]),
            longitudinal_half_length_m=(
                puff.handoff.longitudinal_half_length_m*source_spread_scale),
        )
        puff=integrate_finite_gaussian_puff(handoff,puff.thermodynamics,puff.config)
    if source_model=='gaussian':
        contact=next(i for i,s in enumerate(puff.states) if
            np.linalg.norm(np.array(s.centre_position_m[:2])-np.array([60.,320.]))
            +s.longitudinal_half_length_m>=handoff_radius_m)
        state=puff.states[contact];handoff_time=state.elapsed_s
    else:
        contact=next(i for i,s in enumerate(source_result.states) if
            np.linalg.norm(np.array(s.centre_xy_m)-np.array([60.,320.]))
            +s.apothem_m>=handoff_radius_m)
        state=source_result.states[contact];handoff_time=state.elapsed_s
    if domain_mode not in {'local','sensor_envelope'}:
        raise ValueError('unsupported domain_mode')
    # The local mode is the historical hybrid score.  The sensor-envelope
    # mode is a diagnostic: it removes the analytic exterior packet route so
    # that its contribution can be tested without changing source physics.
    if domain_mode=='local':
        origin=(20.,250.);extent=(110.,140.)
    else:
        positive=dataset.dose_sensors
        sensor_xy=np.array([(60.+s.local_y_m,320.-s.local_x_m)
                            for s in positive])
        # Include the source/fence neighbourhood as well as every receptor;
        # otherwise the enlarged patch would clip the initial cloud or wall.
        anchors=np.vstack((sensor_xy,np.array([
            [10.,270.],[10.,370.],[110.,270.],[110.,370.]])))
        lower=np.floor(np.min(anchors,axis=0)/h)*h
        upper=np.ceil(np.max(anchors,axis=0)/h)*h
        origin=(float(lower[0]),float(lower[1]))
        extent=(float(upper[0]-lower[0]),float(upper[1]-lower[1]))
    if any(not np.isclose(length/h,round(length/h)) for length in extent):raise ValueError('h must tile local domain')
    nx,ny=(round(length/h) for length in extent)
    xe=origin[0]+np.arange(nx+1)*h;ye=origin[1]+np.arange(ny+1)*h
    if source_model=='gaussian':
        sigma=state.longitudinal_half_length_m/puff.config.width_cutoff_sigma
        if not np.isclose(state.longitudinal_half_length_m,state.lateral_half_width_m):raise ValueError('axisymmetry required')
        prob=np.diff(ndtr((ye-state.centre_position_m[1])/sigma))[:,None]*np.diff(ndtr((xe-state.centre_position_m[0])/sigma))[None,:]
        q=np.empty((4,ny,nx));q[0]=source_mass_scale*state.total_mass_kg*prob/h**2;q[1]=source_mass_scale*state.hydrogen_mass_kg*prob/h**2
        q[2]=q[0]*state.bulk_velocity_m_s[0];q[3]=q[0]*state.bulk_velocity_m_s[1]
        unmapped_species=source_mass_scale*state.hydrogen_mass_kg-float(q[1].sum()*h*h)
    else:
        q,inventory=prism_layer_inventory(source_result,contact,xe,ye,radial_momentum=True)
        q*=float(source_mass_scale)
        unmapped_species=float(inventory['exterior_conserved'][1])*float(source_mass_scale)
        if not math.isclose(source_spread_scale,1.0):
            if source_spread_scale <= 0.0:
                raise ValueError('source_spread_scale must be positive')
            # Broaden the finite-prism footprint with a zero-mass Gaussian
            # overlap while renormalizing each conserved component.  The
            # added variance is the declared prism-radius variance change;
            # no receptor dose is used in this operation.
            sigma_m=7.0*math.sqrt(max(source_spread_scale**2-1.0,0.0))
            if sigma_m > 0.0:
                sigma_cells=sigma_m/h
                before=q.copy()
                for component in range(4):
                    filtered=gaussian_filter(q[component],sigma=sigma_cells,
                                             mode='constant',cval=0.0)
                    target=float(before[component].sum())
                    current=float(filtered.sum())
                    if current > 1.0e-30:
                        filtered*=target/current
                    q[component]=filtered
        if (not math.isclose(source_longitudinal_scale,1.0)
                or not math.isclose(source_crosswind_scale,1.0)):
            # Rescale the handoff footprint about the registered cloud centre
            # along the wind/crosswind axes.  At the TI021 six-degree angle the
            # axis rotation is negligible at 5 m, so x/y interpolation is the
            # bounded equivalent.  Every conserved component is renormalized
            # to its pre-transform inventory; no receptor values enter here.
            cx=(float(state.centre_xy_m[0])-origin[0])/h-.5
            cy=(float(state.centre_xy_m[1])-origin[1])/h-.5
            yy,xx=np.indices(q.shape[1:],dtype=float)
            sample_y=cy+(yy-cy)/float(source_crosswind_scale)
            sample_x=cx+(xx-cx)/float(source_longitudinal_scale)
            coords=np.stack((sample_y,sample_x))
            before=q.copy()
            for component in range(4):
                transformed=map_coordinates(q[component],coords,order=1,
                                             mode='constant',cval=0.0)
                target=float(before[component].sum())
                current=float(transformed.sum())
                if current>1.0e-30:
                    transformed*=target/current
                q[component]=transformed
    source_injection=None
    source_injection_metadata={'enabled':False}
    if source_model=='finite_prism' and source_release_duration_s is not None:
        # Complete the declared source history after the near-field handoff.
        # The previous implementation retained only the partial handoff state;
        # this callback transfers the remaining source inventory on the layer
        # clock without using receptor observations.
        source_y=float(closure.mass_fraction_from_mole_fraction(source_cloud.species_mole_fraction))
        source_rho=thermo._density_from_temperature(source_cloud.temperature_k,source_y)
        full_mass=float(source_cloud.exact_volume_m3*source_rho)*float(source_mass_scale)
        release_elapsed=min(max(float(handoff_time),0.0),float(source_release_duration_s))
        remaining_mass=max(full_mass*(1.0-release_elapsed/float(source_release_duration_s)),0.0)
        remaining_duration=max(float(source_release_duration_s)-release_elapsed,0.0)
        centre=np.asarray(state.centre_xy_m,float)
        sigma=max(float(state.apothem_m)/2.0,0.5)
        vel=np.asarray(state.bulk_velocity_xy_m_s,float)
        # Build the footprint on the actual layer grid after its origin is known.
        _xx,_yy=np.meshgrid(origin[0]+(np.arange(nx)+.5)*h,
                            origin[1]+(np.arange(ny)+.5)*h)
        weights=np.exp(-0.5*((_xx-centre[0])**2+(_yy-centre[1])**2)/sigma**2)
        weights/=max(float(weights.sum()*h*h),1.0e-30)
        def source_injection(t_s,dt_s):
            out=np.zeros_like(q)
            if remaining_mass<=0.0 or remaining_duration<=0.0:
                return out
            a0=min(max(float(t_s),0.0),remaining_duration)
            a1=min(max(float(t_s)+float(dt_s),0.0),remaining_duration)
            dm=remaining_mass*max(a1-a0,0.0)/remaining_duration
            if dm<=0.0:return out
            out[0]=dm*weights
            out[1]=out[0]*source_y
            out[2]=out[0]*vel[0]
            out[3]=out[0]*vel[1]
            return out
        source_injection_metadata={
            'enabled':True,'remaining_mass_kg':remaining_mass,
            'remaining_species_kg':remaining_mass*source_y,
            'remaining_duration_s':remaining_duration,
            'handoff_time_s':float(handoff_time),'footprint_sigma_m':sigma}
    wx,wy=curved_fence_faces(xe,ye,radius=wall_radius_m,height=wall_height_m)
    positions=np.array([(60.+s.local_y_m,320.-s.local_x_m,s.height_m) for s in dataset.dose_sensors])
    interface_source=None
    if handoff_interface_path is not None:
        if interface_x_m is None:
            raise ValueError('interface_x_m is required with handoff_interface_path')
        # The resolved inner plane is centred at crosswind coordinate 0
        # (metadata ``inner_crosswind_origin_m=320``), whereas the layer patch
        # uses the global site y coordinate.
        interface_y_edges=(origin[1]+np.arange(ny+1,dtype=float)*h)-320.0
        interface_source=Hdf5DenseLayerPacketSource(
            handoff_interface_path,target_y_edges_m=interface_y_edges)
    handoff_vertical_gaussian, handoff_vertical_top = _handoff_vertical_geometry(
        source_model, state, puff)
    resolved_initial_vertical_gaussian = (
        handoff_vertical_gaussian
        if initial_vertical_gaussian is True else initial_vertical_gaussian
    ) if (initial_vertical_gaussian and
          (vertical_layers > 1
           or vertical_observation_profile == 'handoff_gaussian')) else None
    result=solve_dense_layer_patch(initial_conserved=q,thermo=thermo,origin_xy_m=origin,
        cell_size_m=h,duration_s=600.-handoff_time,wind_xy_m_s=layer_wind_vec,wall_x_height_m=wx,
        wall_y_height_m=wy,receptor_points_m=positions,max_dt_s=dt,
        wind_history_time_s=wind_history_time_s,
        wind_history_xy_m_s=wind_history_xy_m_s,
        wind_history_time_offset_s=handoff_time,
        wind_history_advection_reference_x_m=wind_history_advection_reference_x_m,
        wind_history_advection_speed_m_s=wind_history_advection_speed_m_s,
        wind_history_scale_friction_velocity=wind_history_scale_friction_velocity,
        friction_velocity_m_s=(
            float(friction_velocity_m_s) if layer_wind=='neutral_log' else None),
        dense_entrainment_model=entrainment,
        entrainment_coefficient=entrainment_coefficient,
        entrainment_far_coefficient=entrainment_far_coefficient,
        entrainment_relaxation_time_s=entrainment_relaxation_time_s,
        drag_coefficient=drag_coefficient,
        vertical_observation_profile=vertical_observation_profile,
        dense_vertical_exponent=dense_vertical_exponent,
        vertical_grounding_exponent=vertical_grounding_exponent,
        vertical_elevated_fraction=vertical_elevated_fraction,
        wall_overtop_exponent=wall_overtop_exponent,
        wall_blockage_fraction=wall_blockage_fraction,
        vertical_tail_extension_m=vertical_tail_extension_m,
        vertical_tail_diffusivity_m2_s=vertical_tail_diffusivity_m2_s,
        mass_conserving_vertical_profile=mass_conserving_vertical_profile,
        vertical_observation_mode=vertical_observation_mode,
        vertical_observation_bins=vertical_observation_bins,
        vertical_observation_component_diagnostics=(
            vertical_observation_component_diagnostics),
        vertical_layers=vertical_layers,
        vertical_mixing_diffusivity_m2_s=vertical_mixing_diffusivity_m2_s,
        multilayer_entrainment_placement=multilayer_entrainment_placement,
        multilayer_horizontal_shear_fraction=(
            multilayer_horizontal_shear_fraction),
        multilayer_obstacle_crosswind_shear_rate_s_inv=(
            multilayer_obstacle_crosswind_shear_rate_s_inv),
        multilayer_obstacle_crosswind_decay_length_m=(
            multilayer_obstacle_crosswind_decay_length_m),
        multilayer_relative_transport_interval_s=(
            multilayer_relative_transport_interval_s),
        obstacle_vertical_lift_fraction=obstacle_vertical_lift_fraction,
        obstacle_vertical_lift_decay_m=obstacle_vertical_lift_decay_m,
        wall_impulse_vertical_lift_fraction=wall_impulse_vertical_lift_fraction,
        wall_impulse_elevated_transfer_fraction=wall_impulse_elevated_transfer_fraction,
        wall_impulse_elevated_state_mode=wall_impulse_elevated_state_mode,
        wall_impulse_elevated_conserved_fraction=(
            wall_impulse_elevated_conserved_fraction),
        wall_flux_elevated_transfer_fraction=wall_flux_elevated_transfer_fraction,
        wall_impulse_elevated_transport_fraction=wall_impulse_elevated_transport_fraction,
        wall_impulse_elevated_crosswind_velocity_m_s=(
            wall_impulse_elevated_crosswind_velocity_m_s),
        wall_impulse_elevated_crosswind_diffusivity_m2_s=(
            wall_impulse_elevated_crosswind_diffusivity_m2_s),
        wall_impulse_elevated_crosswind_outward_velocity_m_s=(
            wall_impulse_elevated_crosswind_outward_velocity_m_s),
        wall_impulse_elevated_crosswind_outward_scale_m=(
            wall_impulse_elevated_crosswind_outward_scale_m),
        wall_impulse_elevated_crosswind_outward_decay_s=(
            wall_impulse_elevated_crosswind_outward_decay_s),
        wall_impulse_elevated_crosswind_spinup_s=(
            wall_impulse_elevated_crosswind_spinup_s),
        wall_impulse_elevated_transport_interval_s=(
            wall_impulse_elevated_transport_interval_s),
        wall_impulse_elevated_advection_scale=(
            wall_impulse_elevated_advection_scale),
        wall_impulse_elevated_advection_duration_s=(
            wall_impulse_elevated_advection_duration_s),
        wall_impulse_elevated_return_rate_s_inv=wall_impulse_elevated_return_rate_s_inv,
        wall_impulse_elevated_center_m=wall_impulse_elevated_center_m,
        wall_impulse_elevated_sigma_m=wall_impulse_elevated_sigma_m,
        wall_impulse_elevated_observation_ramp_m=(
            wall_impulse_elevated_observation_ramp_m),
        wall_impulse_elevated_settling_velocity_m_s=(
            wall_impulse_elevated_settling_velocity_m_s),
        wall_impulse_elevated_vertical_diffusivity_m2_s=(
            wall_impulse_elevated_vertical_diffusivity_m2_s),
        wall_impulse_elevated_min_center_m=(
            wall_impulse_elevated_min_center_m),
        wall_impulse_elevated_rise_velocity_m_s=(
            wall_impulse_elevated_rise_velocity_m_s),
        wall_impulse_elevated_rise_time_s=(
            wall_impulse_elevated_rise_time_s),
        wall_impulse_elevated_conserved_center_m=(
            wall_impulse_elevated_conserved_center_m),
        wall_impulse_elevated_conserved_sigma_m=(
            wall_impulse_elevated_conserved_sigma_m),
        wall_impulse_elevated_conserved_settling_velocity_m_s=(
            wall_impulse_elevated_conserved_settling_velocity_m_s),
        wall_impulse_elevated_conserved_vertical_diffusivity_m2_s=(
            wall_impulse_elevated_conserved_vertical_diffusivity_m2_s),
        wall_impulse_elevated_conserved_min_center_m=(
            wall_impulse_elevated_conserved_min_center_m),
        wall_impulse_elevated_conserved_rise_velocity_m_s=(
            wall_impulse_elevated_conserved_rise_velocity_m_s),
        wall_impulse_elevated_conserved_rise_time_s=(
            wall_impulse_elevated_conserved_rise_time_s),
        wall_impulse_elevated_conserved_transport_fraction=(
            wall_impulse_elevated_conserved_transport_fraction),
        wall_impulse_elevated_conserved_crosswind_velocity_m_s=(
            wall_impulse_elevated_conserved_crosswind_velocity_m_s),
        wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s=(
            wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s),
        wall_impulse_elevated_conserved_crosswind_outward_scale_m=(
            wall_impulse_elevated_conserved_crosswind_outward_scale_m),
        wall_impulse_elevated_conserved_center_mode=(
            wall_impulse_elevated_conserved_center_mode),
        wall_impulse_reattached_fraction=wall_impulse_reattached_fraction,
        wall_impulse_reattached_transport_fraction=(
            wall_impulse_reattached_transport_fraction),
        wall_impulse_reattached_transport_height_m=(
            wall_impulse_reattached_transport_height_m),
        wall_impulse_reattached_crosswind_velocity_m_s=(
            wall_impulse_reattached_crosswind_velocity_m_s),
        wall_impulse_reattached_crosswind_outward_velocity_m_s=(
            wall_impulse_reattached_crosswind_outward_velocity_m_s),
        wall_impulse_reattached_crosswind_outward_scale_m=(
            wall_impulse_reattached_crosswind_outward_scale_m),
        wall_impulse_reattached_return_rate_s_inv=(
            wall_impulse_reattached_return_rate_s_inv),
        wall_impulse_reattached_center_m=wall_impulse_reattached_center_m,
        wall_impulse_reattached_sigma_m=wall_impulse_reattached_sigma_m,
        wall_impulse_reattached_settling_velocity_m_s=(
            wall_impulse_reattached_settling_velocity_m_s),
        wall_impulse_reattached_vertical_diffusivity_m2_s=(
            wall_impulse_reattached_vertical_diffusivity_m2_s),
        wall_impulse_reattached_min_center_m=(
            wall_impulse_reattached_min_center_m),
        wall_impulse_wake_capture_rate_s_inv=(
            wall_impulse_wake_capture_rate_s_inv),
        wall_impulse_wake_capture_mode=wall_impulse_wake_capture_mode,
        wall_impulse_wake_return_rate_s_inv=(
            wall_impulse_wake_return_rate_s_inv),
        wall_impulse_wake_transport_fraction=(
            wall_impulse_wake_transport_fraction),
        wall_impulse_wake_crosswind_velocity_m_s=(
            wall_impulse_wake_crosswind_velocity_m_s),
        wall_impulse_wake_crosswind_outward_velocity_m_s=(
            wall_impulse_wake_crosswind_outward_velocity_m_s),
        wall_impulse_wake_crosswind_outward_scale_m=(
            wall_impulse_wake_crosswind_outward_scale_m),
        wall_impulse_wake_observation_fraction=(
            wall_impulse_wake_observation_fraction),
        wall_impulse_wake_observation_center_m=(
            wall_impulse_wake_observation_center_m),
        wall_impulse_wake_observation_sigma_m=(
            wall_impulse_wake_observation_sigma_m),
        wall_impulse_wake_settling_velocity_m_s=(
            wall_impulse_wake_settling_velocity_m_s),
        wall_impulse_wake_vertical_diffusivity_m2_s=(
            wall_impulse_wake_vertical_diffusivity_m2_s),
        wall_impulse_wake_min_center_m=wall_impulse_wake_min_center_m,
        initial_elevated_species_fraction=initial_elevated_species_fraction,
        dense_settling_velocity_m_s=dense_settling_velocity_m_s,
        height_resolved_wall_flux=height_resolved_wall_flux,
        receptor_interpolation=receptor_interpolation,
        wall_tangential_bypass_fraction=wall_tangential_bypass_fraction,
        wall_tip_bypass_fraction=wall_tip_bypass_fraction,
        wall_tip_bypass_power=wall_tip_bypass_power,
        wall_overflow_capture_fraction=wall_overflow_capture_fraction,
        wall_overflow_return_rate_s_inv=wall_overflow_return_rate_s_inv,
        wall_overflow_side_gate_length_m=wall_overflow_side_gate_length_m,
        initial_vertical_gaussian=resolved_initial_vertical_gaussian,
        vertical_domain_top_m=(
            handoff_vertical_top if source_model != 'gaussian' else None
        ),
        horizontal_diffusivity_m2_s=horizontal_diffusivity_m2_s,
        horizontal_diffusion_respect_walls=horizontal_diffusion_respect_walls,
        crosswind_diffusivity_factor=crosswind_diffusivity_factor,
        wake_deficit_strength=wake_deficit_strength,wake_length_m=wake_length_m,
        wake_base_width_m=wake_base_width_m,wake_spread_slope=wake_spread_slope,
        wake_center_xy_m=(60.,320.),wake_radius_m=50.,
        wake_recirculation_strength=wake_recirculation_strength,
        wake_meander_strength=wake_meander_strength,
        wake_shedding_strouhal=wake_shedding_strouhal,
        wake_meander_phase_deg=wake_meander_phase_deg,
        core_meander_strength=core_meander_strength,
        core_meander_length_m=core_meander_length_m,
        core_meander_width_m=core_meander_width_m,
        cavity_shelter_strength=cavity_shelter_strength,
        cavity_entrainment_factor=cavity_entrainment_factor,
        cavity_core_entrainment_reduction=cavity_core_entrainment_reduction,
        cavity_vortex_rate_s_inv=cavity_vortex_rate_s_inv,
        cavity_vortex_strength=cavity_vortex_strength,
        cavity_reservoir_capture_rate_s_inv=cavity_reservoir_capture_rate_s_inv,
        cavity_reservoir_return_rate_s_inv=cavity_reservoir_return_rate_s_inv,
        cavity_reservoir_return_target=cavity_reservoir_return_target,
        cavity_reservoir_reattached_fraction=(
            cavity_reservoir_reattached_fraction),
        cavity_reservoir_capture_mode=cavity_reservoir_capture_mode,
        cavity_reservoir_vertical_capture_mode=(
            cavity_reservoir_vertical_capture_mode),
        cavity_reservoir_spatial_mode=cavity_reservoir_spatial_mode,
        cavity_reservoir_content_mode=cavity_reservoir_content_mode,
        cavity_reservoir_transport_fraction=(
            cavity_reservoir_transport_fraction),
        cavity_reservoir_transport_height_m=(
            cavity_reservoir_transport_height_m),
        cavity_reservoir_transport_mask=cavity_reservoir_transport_mask,
        cavity_reservoir_wake_length_m=cavity_reservoir_wake_length_m,
        cavity_reservoir_wake_half_width_m=(
            cavity_reservoir_wake_half_width_m),
        cavity_reservoir_crosswind_convergence_rate_s_inv=(
            cavity_reservoir_crosswind_convergence_rate_s_inv),
        cavity_reservoir_crosswind_target_offset_m=(
            cavity_reservoir_crosswind_target_offset_m),
        local_reservoir_observation_fraction=(
            local_reservoir_observation_fraction),
        local_reservoir_observation_center_m=(
            local_reservoir_observation_center_m),
        local_reservoir_observation_sigma_m=(
            local_reservoir_observation_sigma_m),
        cavity_reservoir_return_sigma_m=cavity_reservoir_return_sigma_m,
        cavity_reservoir_return_radius_m=cavity_reservoir_return_radius_m,
        cavity_reservoir_return_offset_m=cavity_reservoir_return_offset_m,
        cavity_reservoir_return_crosswind_scale=cavity_reservoir_return_crosswind_scale,
        cavity_reservoir_face_return_fraction=cavity_reservoir_face_return_fraction,
        cavity_reservoir_capture_radius_m=cavity_reservoir_capture_radius_m,
        cavity_reservoir_start_time_s=cavity_reservoir_start_time_s,
        cavity_reservoir_return_delay_s=cavity_reservoir_return_delay_s,
        cavity_return_streamwise_bias=cavity_return_streamwise_bias,
        cavity_vertical_exponent=cavity_vertical_exponent,
        cavity_vertical_radial_decay_power=cavity_vertical_radial_decay_power,
        cavity_vertical_mixing_factor=cavity_vertical_mixing_factor,
        exterior_entrainment_far_factor=exterior_entrainment_far_factor,
        exterior_entrainment_decay_length_m=exterior_entrainment_decay_length_m,
        cavity_backflow_strength=cavity_backflow_strength,
        cavity_backflow_rate_s_inv=cavity_backflow_rate_s_inv,
        cavity_backflow_frame=cavity_backflow_frame,
        cavity_stagnation_strength=cavity_stagnation_strength,
        exterior_reentry_fraction=exterior_reentry_fraction,
        exterior_reentry_delay_s=exterior_reentry_delay_s,
        exterior_reentry_sigma_m=exterior_reentry_sigma_m,
        interface_x_m=interface_x_m,
        interface_source_callback=interface_source,
        interface_time_offset_s=handoff_time,
        source_injection_callback=source_injection,
        cavity_geometry=cavity_geometry,
        cavity_center_xy_m=(60.,320.),cavity_radius_m=cavity_radius_m,
        cavity_align_with_wind=cavity_align_with_wind,
        gravity_scale=gravity_scale)
    interface_diagnostics=None
    if interface_source is not None:
        interface_diagnostics={
            'packet_schema':'degali.marked-layer-donor-packet.v1',
            'component_order':['mixture_mass','species_mass','momentum_x','momentum_y'],
            'time_integration':'exact integral of piecewise linear donor rates',
            'donor_history_start_s':float(interface_source.time_s[0]),
            'donor_history_end_s':float(interface_source.time_s[-1]),
            'integrated_species_mass_kg':interface_source.integrated_species_mass_kg,
            'integrated_marked_mixture_mass_kg':interface_source.integrated_total_mass_kg,
            'integrated_horizontal_momentum_kg_m_s':interface_source.integrated_horizontal_momentum_kg_m_s.tolist(),
            'discarded_vertical_impulse_kg_m_s':interface_source.discarded_vertical_impulse_kg_m_s,
            'absolute_discarded_vertical_impulse_kg_m_s':interface_source.absolute_discarded_vertical_impulse_kg_m_s,
            'omitted_lateral_species_mass_kg':interface_source.omitted_species_mass_kg,
            'vertical_momentum_retained':False,
            'sensible_enthalpy_retained':False,
            'full_3d_momentum_handoff_gate':False,
        }
        interface_source.close()
    traces=result.species_traces_vol_percent.copy()
    inside=(positions[:,:2]>=origin).all(axis=1)&(positions[:,:2]<=np.array(origin)+extent).all(axis=1)
    conversion=(100*8.314462618*closure.ambient_temperature_k
                /(closure.ambient_pressure_pa
                  *closure.species_molar_mass_kg_mol))
    for i in np.flatnonzero(~inside):
        traces[i]=exterior_trace(result,positions[i],layer_wind_vec,h,handoff_time,
                                 transport=exterior_transport,
                                 density_current_strength=density_current_strength,
                                 exterior_lateral_diffusivity_factor=exterior_lateral_diffusivity_factor)*conversion
    # Linear source scaling is explicitly a diagnostic sensitivity, not a
    # calibrated production correction. It keeps the short source prefix and
    # handed-off inventory on the same scale for a closed 600 s comparison.
    if source_model=='gaussian':
        before=source_mass_scale*binary_puff_traces(puff,positions,state_count=contact)
        before_times=np.array([s.elapsed_s for s in puff.states[:contact]])
    else:
        before=source_result.traces(positions)[:,:contact]
        before_times=np.array([s.elapsed_s for s in source_result.states[:contact]])
    traces=np.concatenate((before,traces),axis=1)
    times=np.concatenate((before_times,result.time_s+handoff_time))
    component_labels=[]
    component_concentration=None
    if result.observation_component_traces_vol_percent is not None:
        components={
            'source_before_handoff':np.concatenate((
                before,np.zeros((positions.shape[0],result.time_s.size))),axis=1)
        }
        for label,values in result.observation_component_traces_vol_percent.items():
            values=np.asarray(values,float)
            if values.shape != (positions.shape[0],result.time_s.size):
                raise RuntimeError('vertical observation component trace shape mismatch')
            components[label]=np.concatenate((
                np.zeros((positions.shape[0],before_times.size)),values),axis=1)
        exterior=np.zeros_like(traces)
        exterior[np.flatnonzero(~inside),before_times.size:]=traces[
            np.flatnonzero(~inside),before_times.size:]
        components['exterior_packet_continuation']=exterior
        component_labels=list(components)
        component_concentration=np.stack(
            [components[label] for label in component_labels],axis=0)
        if not np.allclose(component_concentration.sum(axis=0),traces,
                           rtol=1.e-12,atol=1.e-12):
            raise RuntimeError('full-event observation component attribution does not close')
    prediction_runtime=time.perf_counter()-started
    observed=np.array([s.dose_percent_s for s in dataset.dose_sensors]);positive=observed>0
    dose=np.trapezoid(traces,times,axis=1)
    report={'schema':'degali.ti021-dense-layer-hybrid.v1','gas':gas_label,'duration_s':600.,
        'status':('experimental_multilayer_hybrid_not_field_qualified'
                  if vertical_layers>1 else
                  'experimental_single_layer_hybrid_not_field_qualified'),
        'dose_mape_percent':float(np.abs(dose[positive]/observed[positive]-1).mean()*100),
        'dose_mae_vol_percent_s':float(np.abs(dose[positive]-observed[positive]).mean()),
        'normalized_mean_bias':float((dose[positive]-observed[positive]).sum()/observed[positive].sum()),
        'sensor_count':145,'positive_observed_sensor_count':int(positive.sum()),
        'prediction_runtime_s':prediction_runtime,'handoff_time_s':handoff_time,
        'unmapped_species_kg':unmapped_species,'unmapped_species_fraction':unmapped_species/(
            source_mass_scale*(state.hydrogen_mass_kg if source_model=='gaussian' else state.species_mass_kg)),
        'patch_cells':[nx,ny],'horizontal_cell_size_m':h,'max_dt_s':dt,
        'domain_mode':domain_mode,
        'vertical_observation_profile':vertical_observation_profile,
        'dense_vertical_exponent':float(dense_vertical_exponent),
        'vertical_grounding_exponent':float(vertical_grounding_exponent),
        'vertical_elevated_fraction':float(vertical_elevated_fraction),
        'wall_overtop_exponent':float(wall_overtop_exponent),
        'wall_blockage_fraction':float(wall_blockage_fraction),
        'vertical_tail_extension_m':float(vertical_tail_extension_m),
        'vertical_tail_diffusivity_m2_s':float(vertical_tail_diffusivity_m2_s),
        'mass_conserving_vertical_profile':mass_conserving_vertical_profile,
        'vertical_observation_mode':vertical_observation_mode,
        'vertical_observation_bins':vertical_observation_bins,
        'vertical_observation_component_diagnostics':bool(
            vertical_observation_component_diagnostics),
        'vertical_observation_component_labels':component_labels,
        'vertical_layers':int(vertical_layers),
        'initial_vertical_gaussian':bool(initial_vertical_gaussian and vertical_layers>1),
        'initial_vertical_gaussian_center_m':(
            None if resolved_initial_vertical_gaussian is None
            else float(resolved_initial_vertical_gaussian[0])),
        'initial_vertical_gaussian_sigma_m':(
            None if resolved_initial_vertical_gaussian is None
            else float(resolved_initial_vertical_gaussian[1])),
        'handoff_vertical_domain_top_m':(
            None if handoff_vertical_top is None else float(handoff_vertical_top)),
        'multilayer_entrainment_placement':multilayer_entrainment_placement,
        'vertical_mixing_diffusivity_m2_s':float(vertical_mixing_diffusivity_m2_s),
        'multilayer_horizontal_shear_fraction':float(
            multilayer_horizontal_shear_fraction),
        'multilayer_obstacle_crosswind_shear_rate_s_inv':float(
            multilayer_obstacle_crosswind_shear_rate_s_inv),
        'multilayer_obstacle_crosswind_decay_length_m':float(
            multilayer_obstacle_crosswind_decay_length_m),
        'multilayer_relative_transport_interval_s':float(
            multilayer_relative_transport_interval_s),
        'obstacle_vertical_lift_fraction':float(obstacle_vertical_lift_fraction),
        'obstacle_vertical_lift_decay_m':float(obstacle_vertical_lift_decay_m),
        'wall_impulse_vertical_lift_fraction':float(
            wall_impulse_vertical_lift_fraction),
        'wall_impulse_elevated_transfer_fraction':float(
            wall_impulse_elevated_transfer_fraction),
        'wall_impulse_elevated_state_mode':wall_impulse_elevated_state_mode,
        'wall_impulse_elevated_conserved_fraction':float(
            wall_impulse_elevated_conserved_fraction),
        'wall_flux_elevated_transfer_fraction':float(
            wall_flux_elevated_transfer_fraction),
        'wall_impulse_elevated_transport_fraction':float(
            wall_impulse_elevated_transport_fraction),
        'wall_impulse_elevated_crosswind_velocity_m_s':float(
            wall_impulse_elevated_crosswind_velocity_m_s),
        'wall_impulse_elevated_crosswind_diffusivity_m2_s':float(
            wall_impulse_elevated_crosswind_diffusivity_m2_s),
        'wall_impulse_elevated_crosswind_outward_velocity_m_s':float(
            wall_impulse_elevated_crosswind_outward_velocity_m_s),
        'wall_impulse_elevated_crosswind_outward_scale_m':float(
            wall_impulse_elevated_crosswind_outward_scale_m),
        'wall_impulse_elevated_crosswind_outward_decay_s':float(
            wall_impulse_elevated_crosswind_outward_decay_s),
        'wall_impulse_elevated_crosswind_spinup_s':float(
            wall_impulse_elevated_crosswind_spinup_s),
        'wall_impulse_elevated_transport_interval_s':float(
            wall_impulse_elevated_transport_interval_s),
        'wall_impulse_elevated_advection_scale':float(
            wall_impulse_elevated_advection_scale),
        'wall_impulse_elevated_advection_duration_s':float(
            wall_impulse_elevated_advection_duration_s),
        'wall_impulse_elevated_return_rate_s_inv':float(
            wall_impulse_elevated_return_rate_s_inv),
        'wall_impulse_elevated_center_m':float(wall_impulse_elevated_center_m),
        'wall_impulse_elevated_sigma_m':float(wall_impulse_elevated_sigma_m),
        'wall_impulse_elevated_observation_ramp_m':float(
            wall_impulse_elevated_observation_ramp_m),
        'wall_impulse_elevated_settling_velocity_m_s':float(
            wall_impulse_elevated_settling_velocity_m_s),
        'wall_impulse_elevated_vertical_diffusivity_m2_s':float(
            wall_impulse_elevated_vertical_diffusivity_m2_s),
        'wall_impulse_elevated_min_center_m':float(
            wall_impulse_elevated_min_center_m),
        'wall_impulse_elevated_rise_velocity_m_s':float(
            wall_impulse_elevated_rise_velocity_m_s),
        'wall_impulse_elevated_rise_time_s':float(
            wall_impulse_elevated_rise_time_s),
        'wall_impulse_elevated_conserved_center_m':(
            None if wall_impulse_elevated_conserved_center_m is None
            else float(wall_impulse_elevated_conserved_center_m)),
        'wall_impulse_elevated_conserved_sigma_m':(
            None if wall_impulse_elevated_conserved_sigma_m is None
            else float(wall_impulse_elevated_conserved_sigma_m)),
        'wall_impulse_elevated_conserved_settling_velocity_m_s':(
            None if wall_impulse_elevated_conserved_settling_velocity_m_s is None
            else float(wall_impulse_elevated_conserved_settling_velocity_m_s)),
        'wall_impulse_elevated_conserved_vertical_diffusivity_m2_s':(
            None if wall_impulse_elevated_conserved_vertical_diffusivity_m2_s is None
            else float(wall_impulse_elevated_conserved_vertical_diffusivity_m2_s)),
        'wall_impulse_elevated_conserved_min_center_m':(
            None if wall_impulse_elevated_conserved_min_center_m is None
            else float(wall_impulse_elevated_conserved_min_center_m)),
        'wall_impulse_elevated_conserved_rise_velocity_m_s':float(
            wall_impulse_elevated_conserved_rise_velocity_m_s),
        'wall_impulse_elevated_conserved_rise_time_s':float(
            wall_impulse_elevated_conserved_rise_time_s),
        'wall_impulse_elevated_conserved_transport_fraction':(
            None if wall_impulse_elevated_conserved_transport_fraction is None
            else float(wall_impulse_elevated_conserved_transport_fraction)),
        'wall_impulse_elevated_conserved_crosswind_velocity_m_s':(
            None if wall_impulse_elevated_conserved_crosswind_velocity_m_s is None
            else float(wall_impulse_elevated_conserved_crosswind_velocity_m_s)),
        'wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s':(
            None if wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s is None
            else float(
                wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s)),
        'wall_impulse_elevated_conserved_crosswind_outward_scale_m':(
            None if wall_impulse_elevated_conserved_crosswind_outward_scale_m is None
            else float(wall_impulse_elevated_conserved_crosswind_outward_scale_m)),
        'wall_impulse_elevated_conserved_center_mode':(
            wall_impulse_elevated_conserved_center_mode),
        'wall_impulse_reattached_fraction':float(
            wall_impulse_reattached_fraction),
        'wall_impulse_reattached_transport_fraction':float(
            wall_impulse_reattached_transport_fraction),
        'wall_impulse_reattached_transport_height_m':(
            None if wall_impulse_reattached_transport_height_m is None
            else float(wall_impulse_reattached_transport_height_m)),
        'wall_impulse_reattached_crosswind_velocity_m_s':float(
            wall_impulse_reattached_crosswind_velocity_m_s),
        'wall_impulse_reattached_crosswind_outward_velocity_m_s':float(
            wall_impulse_reattached_crosswind_outward_velocity_m_s),
        'wall_impulse_reattached_crosswind_outward_scale_m':float(
            wall_impulse_reattached_crosswind_outward_scale_m),
        'wall_impulse_reattached_return_rate_s_inv':float(
            wall_impulse_reattached_return_rate_s_inv),
        'wall_impulse_reattached_center_m':float(
            wall_impulse_reattached_center_m),
        'wall_impulse_reattached_sigma_m':float(
            wall_impulse_reattached_sigma_m),
        'wall_impulse_reattached_settling_velocity_m_s':float(
            wall_impulse_reattached_settling_velocity_m_s),
        'wall_impulse_reattached_vertical_diffusivity_m2_s':float(
            wall_impulse_reattached_vertical_diffusivity_m2_s),
        'wall_impulse_reattached_min_center_m':float(
            wall_impulse_reattached_min_center_m),
        'wall_impulse_wake_capture_rate_s_inv':float(
            wall_impulse_wake_capture_rate_s_inv),
        'wall_impulse_wake_capture_mode':wall_impulse_wake_capture_mode,
        'wall_impulse_wake_return_rate_s_inv':float(
            wall_impulse_wake_return_rate_s_inv),
        'wall_impulse_wake_transport_fraction':float(
            wall_impulse_wake_transport_fraction),
        'wall_impulse_wake_crosswind_velocity_m_s':float(
            wall_impulse_wake_crosswind_velocity_m_s),
        'wall_impulse_wake_crosswind_outward_velocity_m_s':float(
            wall_impulse_wake_crosswind_outward_velocity_m_s),
        'wall_impulse_wake_crosswind_outward_scale_m':float(
            wall_impulse_wake_crosswind_outward_scale_m),
        'wall_impulse_wake_observation_fraction':float(
            wall_impulse_wake_observation_fraction),
        'wall_impulse_wake_observation_center_m':(
            None if wall_impulse_wake_observation_center_m is None
            else float(wall_impulse_wake_observation_center_m)),
        'wall_impulse_wake_observation_sigma_m':float(
            wall_impulse_wake_observation_sigma_m),
        'wall_impulse_wake_settling_velocity_m_s':float(
            wall_impulse_wake_settling_velocity_m_s),
        'wall_impulse_wake_vertical_diffusivity_m2_s':float(
            wall_impulse_wake_vertical_diffusivity_m2_s),
        'wall_impulse_wake_min_center_m':float(
            wall_impulse_wake_min_center_m),
        'initial_elevated_species_fraction':float(initial_elevated_species_fraction),
        'dense_settling_velocity_m_s':float(dense_settling_velocity_m_s),
        'height_resolved_wall_flux':bool(height_resolved_wall_flux),
        'receptor_interpolation':receptor_interpolation,
        'wall_tangential_bypass_fraction':float(wall_tangential_bypass_fraction),
        'wall_tip_bypass_fraction':float(wall_tip_bypass_fraction),
        'wall_tip_bypass_power':float(wall_tip_bypass_power),
        'wall_overflow_capture_fraction':float(wall_overflow_capture_fraction),
        'wall_overflow_return_rate_s_inv':float(
            wall_overflow_return_rate_s_inv),
        'wall_overflow_side_gate_length_m':float(
            wall_overflow_side_gate_length_m),
        'horizontal_diffusivity_m2_s':float(horizontal_diffusivity_m2_s),
        'horizontal_diffusion_respect_walls':horizontal_diffusion_respect_walls,
        'crosswind_diffusivity_factor':float(crosswind_diffusivity_factor),
        'wake_deficit_strength':float(wake_deficit_strength),
        'wake_length_m':float(wake_length_m),
        'wake_base_width_m':float(wake_base_width_m),
        'wake_spread_slope':float(wake_spread_slope),
        'wake_recirculation_strength':float(wake_recirculation_strength),
        'wake_meander_strength':float(wake_meander_strength),
        'wake_shedding_strouhal':float(wake_shedding_strouhal),
        'wake_meander_phase_deg':float(wake_meander_phase_deg),
        'core_meander_strength':float(core_meander_strength),
        'core_meander_length_m':float(core_meander_length_m),
        'core_meander_width_m':float(core_meander_width_m),
        'cavity_shelter_strength':float(cavity_shelter_strength),
        'cavity_entrainment_factor':float(cavity_entrainment_factor),
        'cavity_core_entrainment_reduction':float(cavity_core_entrainment_reduction),
        'cavity_vortex_rate_s_inv':float(cavity_vortex_rate_s_inv),
        'cavity_vortex_strength':float(cavity_vortex_strength),
        'cavity_reservoir_capture_rate_s_inv':float(cavity_reservoir_capture_rate_s_inv),
        'cavity_reservoir_return_rate_s_inv':float(cavity_reservoir_return_rate_s_inv),
        'cavity_reservoir_return_target':cavity_reservoir_return_target,
        'cavity_reservoir_reattached_fraction':float(
            cavity_reservoir_reattached_fraction),
        'cavity_reservoir_capture_mode':cavity_reservoir_capture_mode,
        'cavity_reservoir_vertical_capture_mode':(
            cavity_reservoir_vertical_capture_mode),
        'cavity_reservoir_spatial_mode':cavity_reservoir_spatial_mode,
        'cavity_reservoir_content_mode':cavity_reservoir_content_mode,
        'cavity_reservoir_transport_fraction':float(
            cavity_reservoir_transport_fraction),
        'cavity_reservoir_transport_height_m':(
            None if cavity_reservoir_transport_height_m is None
            else float(cavity_reservoir_transport_height_m)),
        'cavity_reservoir_transport_mask':cavity_reservoir_transport_mask,
        'cavity_reservoir_wake_length_m':float(
            cavity_reservoir_wake_length_m),
        'cavity_reservoir_wake_half_width_m':(
            None if cavity_reservoir_wake_half_width_m is None
            else float(cavity_reservoir_wake_half_width_m)),
        'cavity_reservoir_crosswind_convergence_rate_s_inv':float(
            cavity_reservoir_crosswind_convergence_rate_s_inv),
        'cavity_reservoir_crosswind_target_offset_m':float(
            cavity_reservoir_crosswind_target_offset_m),
        'local_reservoir_observation_fraction':float(
            local_reservoir_observation_fraction),
        'local_reservoir_observation_center_m':(
            None if local_reservoir_observation_center_m is None
            else float(local_reservoir_observation_center_m)),
        'local_reservoir_observation_sigma_m':float(
            local_reservoir_observation_sigma_m),
        'cavity_reservoir_return_sigma_m':None if cavity_reservoir_return_sigma_m is None else float(cavity_reservoir_return_sigma_m),
        'cavity_reservoir_return_radius_m':float(cavity_reservoir_return_radius_m),
        'cavity_reservoir_return_offset_m':float(cavity_reservoir_return_offset_m),
        'cavity_reservoir_return_crosswind_scale':float(cavity_reservoir_return_crosswind_scale),
        'cavity_reservoir_face_return_fraction':float(cavity_reservoir_face_return_fraction),
        'cavity_reservoir_capture_radius_m':float(cavity_reservoir_capture_radius_m),
        'cavity_reservoir_start_time_s':float(cavity_reservoir_start_time_s),
        'cavity_reservoir_return_delay_s':float(cavity_reservoir_return_delay_s),
        'cavity_return_streamwise_bias':float(cavity_return_streamwise_bias),
        'cavity_vertical_exponent':float(cavity_vertical_exponent),
        'cavity_vertical_radial_decay_power':float(
            cavity_vertical_radial_decay_power),
        'exterior_entrainment_far_factor':float(
            exterior_entrainment_far_factor),
        'exterior_entrainment_decay_length_m':float(
            exterior_entrainment_decay_length_m),
        'cavity_vertical_mixing_factor':float(cavity_vertical_mixing_factor),
        'cavity_backflow_strength':float(cavity_backflow_strength),
        'cavity_backflow_rate_s_inv':float(cavity_backflow_rate_s_inv),
        'cavity_backflow_frame':cavity_backflow_frame,
        'cavity_stagnation_strength':float(cavity_stagnation_strength),
        'exterior_reentry_fraction':float(exterior_reentry_fraction),
        'exterior_reentry_delay_s':float(exterior_reentry_delay_s),
        'exterior_reentry_sigma_m':float(exterior_reentry_sigma_m),
        'cavity_geometry':cavity_geometry,
        'cavity_align_with_wind':bool(cavity_align_with_wind),
        'cavity_radius_m':float(cavity_radius_m),
        'exterior_transport':exterior_transport,
        'exterior_lateral_diffusivity_factor':float(exterior_lateral_diffusivity_factor),
        'density_current_strength':float(density_current_strength),
        'handoff_interface_path':None if handoff_interface_path is None else str(Path(handoff_interface_path).resolve()),
        'interface_x_m':None if interface_x_m is None else float(interface_x_m),
        'wind_history_path':None if wind_history_path is None else str(wind_history_path),
        'wind_history_sha256':wind_history_sha256,
        'wind_history_sample_count':0 if wind_history_time_s is None else int(wind_history_time_s.size),
        'wind_history_component_mode':wind_history_component_mode,
        'wind_history_direction_offset_deg':float(
            wind_history_direction_offset_deg),
        'wind_history_vector_scale_to_layer_speed':wind_history_vector_scale,
        'wind_history_advection_reference_x_m':(
            None if wind_history_advection_reference_x_m is None
            else float(wind_history_advection_reference_x_m)),
        'wind_history_advection_speed_m_s':(
            None if wind_history_advection_speed_m_s is None
            else float(wind_history_advection_speed_m_s)),
        'wind_history_scale_friction_velocity':bool(wind_history_scale_friction_velocity),
        'wind_history_applied_after_handoff_only':wind_history_path is not None,
        'source_time_step_s':source_dt,'layer_wind':layer_wind,
        'source_mass_scale':float(source_mass_scale),'source_model':source_model,
        'source_front_profile':source_front_profile,
        'source_initial_radius_m':float(source_initial_radius_m),
        'source_initial_diameter_m':float(source_initial_radius_m),
        'source_initial_apothem_m':float(source_initial_radius_m)/2.0,
        'source_geometry_input_semantics':'historical parameter name; value is full 12-sided-prism diameter',
        'observation_species_mode':observation_species_mode,
        'source_species_mole_fraction':float(source_mole_fraction),
        'handoff_profile':(
            'match_front_profile'
            if source_front_profile in {
                'longitudinal_ramp', 'longitudinal_ramp_anisotropic',
                'longitudinal_ramp_gaussian',
            }
            else 'uniform'
        ),
        'source_radial_damping_coefficient':float(source_radial_damping_coefficient),
        'source_radial_damping_aspect_ratio':float(source_radial_damping_aspect_ratio),
        'source_minimum_height_m':float(source_minimum_height_m),
        'source_wind_scale':float(source_wind_scale),
        'source_spread_scale':float(source_spread_scale),
        'source_longitudinal_scale':float(source_longitudinal_scale),
        'source_crosswind_scale':float(source_crosswind_scale),
        'source_initial_velocity_scale':float(source_initial_velocity_scale),
        'source_release_duration_s':None if source_release_duration_s is None else float(source_release_duration_s),
        'source_injection':source_injection_metadata,
        'puff_shear_entrainment_coefficient':float(puff_shear_entrainment_coefficient),
        'puff_ambient_turbulence_coefficient':float(puff_ambient_turbulence_coefficient),
        'puff_top_shear_entrainment_coefficient':None if puff_top_shear_entrainment_coefficient is None else float(puff_top_shear_entrainment_coefficient),
        'puff_top_ambient_turbulence_coefficient':None if puff_top_ambient_turbulence_coefficient is None else float(puff_top_ambient_turbulence_coefficient),
        'puff_side_shear_entrainment_coefficient':None if puff_side_shear_entrainment_coefficient is None else float(puff_side_shear_entrainment_coefficient),
        'puff_side_ambient_turbulence_coefficient':None if puff_side_ambient_turbulence_coefficient is None else float(puff_side_ambient_turbulence_coefficient),
        'puff_gravity_spreading_coefficient':float(puff_gravity_spreading_coefficient),
        'puff_form_drag_coefficient':float(puff_form_drag_coefficient),
        'handoff_radius_m':float(handoff_radius_m),
        'wind_speed_m_s':float(wind_speed_m_s),'layer_wind_speed_m_s':float(layer_wind_speed_m_s),
        'friction_velocity_m_s':float(friction_velocity_m_s),'wind_angle_deg':float(wind_angle_deg),
        'wall_height_m':float(wall_height_m),
        'wall_radius_m':float(wall_radius_m),
        'gravity_scale':float(gravity_scale),
        'dense_entrainment_model':entrainment,
        'entrainment_coefficient':float(entrainment_coefficient),
        'entrainment_far_coefficient':(
            None if entrainment_far_coefficient is None
            else float(entrainment_far_coefficient)),
        'entrainment_relaxation_time_s':float(entrainment_relaxation_time_s),
        'drag_coefficient':float(drag_coefficient),
        'source_momentum_and_mixture_mass_transferred':True,'energy_model':'isothermal; source enthalpy not advanced',
        'handoff_temperature_departure_k':(state.bulk_temperature_k-thermo.ambient_temperature
            if source_model=='gaussian' else None),
        'exterior_time_marched_grid_cells':0,'all_boundary_species_packets_continued':True,
        'experimental_concentration_time_series_available':False,'metric':'dose MAPE; not instantaneous H2 concentration MAPE',
        'configuration_used_observed_dose_peak_arrival_values':False,
        'workbook_sha256':hashlib.sha256(workbook.read_bytes()).hexdigest(),
        'diagnostics':result.diagnostics,
        'interface_diagnostics':interface_diagnostics,
        'remaining_physics':['3D vertical bypass/separation/wake turbulence not resolved',
            *(['legacy mass-fraction observation reshaping does not conserve reconstructed species volume']
              if vertical_layers==1 and vertical_observation_mode=='legacy_mass_fraction' else []),
            *(['legacy horizontal turbulent diffusion ignores finite-height wall faces']
              if horizontal_diffusivity_m2_s>0.0 and not horizontal_diffusion_respect_walls else []),
            *([] if vertical_layers>1 else ['uniform vertical layer approximation; no continuous stratified profile']),
            *(['vertical momentum equation is not resolved; horizontal interlayer shear uses a reduced conservative closure']
              if (vertical_layers>1 and (
                  multilayer_horizontal_shear_fraction>0.0
                  or multilayer_obstacle_crosswind_shear_rate_s_inv!=0.0))
              else (['vertical momentum equation and interlayer shear are not resolved']
                    if vertical_layers>1 else [])),
            *([] if exterior_transport=='density_current' else
              ['one-way exterior continuation; patch re-entry is not recoupled']),
            *(['reduced exterior density-current route is not recoupled to the patch']
              if exterior_transport=='density_current' else []),
            'Gaussian initial footprint contact neglects earlier tails hitting the wall',
            'unmapped exterior source tails recorded but not propagated',
            *([] if handoff_interface_path is not None else
              ['interior 3-D handoff packet interface is disabled'])],
        'sensor_rows':[{'label':f'TI021-S{i+1:03d}','observed_dose_vol_percent_s':float(observed[i]),
            'predicted_dose_vol_percent_s':float(dose[i]),'route':'layer_patch' if inside[i] else 'all_face_packets'} for i in range(145)]}
    report['implementation_sha256']={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        for path in ('src/degali/addons/dense_layer_patch.py','src/degali/addons/conservative_column_observation.py','src/degali/addons/binary_integral_puff.py',
                     'src/degali/addons/finite_puff.py','src/degali/addons/hdf5_dense_layer_interface.py',
                     'tools/score_ti021_dense_layer_hybrid.py')}
    output.parent.mkdir(parents=True,exist_ok=True)
    npz_payload=dict(time_s=times,concentration_vol_percent=traces,
        observed_dose_vol_percent_s=observed,predicted_dose_vol_percent_s=dose,
        layer_species_mass_kg=result.species_mass_kg,
        layer_mixture_mass_kg=result.mixture_mass_kg)
    if component_concentration is not None:
        npz_payload['observation_component_labels']=np.asarray(component_labels)
        npz_payload['observation_component_concentration_vol_percent']=(
            component_concentration)
    np.savez_compressed(output.with_suffix('.npz'),**npz_payload)
    report['runtime_including_output_s']=time.perf_counter()-started
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='sensor_rows'},indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cell-size-m',type=float,default=5.);p.add_argument('--time-step-s',type=float,default=.25)
    p.add_argument('--source-time-step-s',type=float,default=.025)
    p.add_argument('--layer-wind',choices=['uniform','neutral_log'],default='neutral_log')
    p.add_argument('--entrainment',choices=['constant','parker_bulk_suppression'],default='constant')
    p.add_argument('--entrainment-coefficient',type=float,default=.1)
    p.add_argument('--drag-coefficient',type=float,default=.2)
    p.add_argument('--domain-mode',choices=['local','sensor_envelope'],default='local')
    p.add_argument('--source-mass-scale',type=float,default=1.0)
    p.add_argument('--source-model',choices=['gaussian','finite_prism'],default='gaussian')
    p.add_argument('--source-initial-radius-m',type=float,default=14.0,
                   help='legacy name: full initial 12-sided-prism diameter in metres')
    p.add_argument('--observation-species-mode',
                   choices=['r12_component','released_mixture'],
                   default='r12_component')
    p.add_argument('--source-front-profile',choices=['uniform','grounded_head','entrained_edge','diffuse_edge','smooth_edge','vertical_tail','wake_tail','anisotropic_support','wind_anisotropic_leading','longitudinal_ramp','longitudinal_ramp_gaussian'],default='uniform')
    p.add_argument('--source-radial-damping-coefficient',type=float,default=0.0)
    p.add_argument('--source-radial-damping-aspect-ratio',type=float,default=0.0)
    p.add_argument('--source-minimum-height-m',type=float,default=0.0)
    p.add_argument('--source-wind-scale',type=float,default=1.0)
    p.add_argument('--source-spread-scale',type=float,default=1.0)
    p.add_argument('--source-longitudinal-scale',type=float,default=1.0)
    p.add_argument('--source-crosswind-scale',type=float,default=1.0)
    p.add_argument('--source-initial-velocity-scale',type=float,default=0.0)
    p.add_argument('--source-release-duration-s',type=float,default=None)
    p.add_argument('--puff-shear-entrainment-coefficient',type=float,default=.1)
    p.add_argument('--puff-ambient-turbulence-coefficient',type=float,default=.1)
    p.add_argument('--puff-top-shear-entrainment-coefficient',type=float,default=None)
    p.add_argument('--puff-top-ambient-turbulence-coefficient',type=float,default=None)
    p.add_argument('--puff-side-shear-entrainment-coefficient',type=float,default=None)
    p.add_argument('--puff-side-ambient-turbulence-coefficient',type=float,default=None)
    p.add_argument('--puff-gravity-spreading-coefficient',type=float,default=.75)
    p.add_argument('--puff-form-drag-coefficient',type=float,default=.2)
    p.add_argument('--handoff-radius-m',type=float,default=50.0)
    p.add_argument('--vertical-observation-profile',choices=['top_hat','dense_power','exponential_tail','diffusive_tail','handoff_gaussian','multilayer'],default='top_hat')
    p.add_argument('--dense-vertical-exponent',type=float,default=0.0)
    p.add_argument('--vertical-grounding-exponent',type=float,default=0.0)
    p.add_argument('--vertical-elevated-fraction',type=float,default=0.0)
    p.add_argument('--wall-overtop-exponent',type=float,default=0.0)
    p.add_argument('--wall-blockage-fraction',type=float,default=1.0)
    p.add_argument('--vertical-tail-extension-m',type=float,default=0.0)
    p.add_argument('--vertical-tail-diffusivity-m2-s',type=float,default=1.0)
    p.add_argument('--mass-conserving-vertical-profile',action='store_true')
    p.add_argument('--vertical-observation-mode',
                   choices=['legacy_mass_fraction','conservative_volume'],
                   default='conservative_volume')
    p.add_argument('--vertical-observation-bins',type=int,default=64)
    p.add_argument('--vertical-observation-component-diagnostics',action='store_true')
    p.add_argument('--vertical-layers',type=int,default=1)
    p.add_argument('--no-initial-vertical-gaussian',dest='initial_vertical_gaussian',action='store_false')
    p.add_argument('--vertical-mixing-diffusivity-m2-s',type=float,default=0.0)
    p.add_argument('--multilayer-entrainment-placement',
                   choices=['top','legacy_ground'],default='top')
    p.add_argument('--multilayer-horizontal-shear-fraction',type=float,default=0.0)
    p.add_argument('--multilayer-obstacle-crosswind-shear-rate-s-inv',type=float,default=0.0)
    p.add_argument('--multilayer-obstacle-crosswind-decay-length-m',type=float,default=70.0)
    p.add_argument('--multilayer-relative-transport-interval-s',type=float,default=0.0)
    p.add_argument('--height-resolved-wall-flux',action='store_true')
    p.add_argument('--wall-tangential-bypass-fraction',type=float,default=0.0)
    p.add_argument('--wall-tip-bypass-fraction',type=float,default=0.0)
    p.add_argument('--wall-tip-bypass-power',type=float,default=2.0)
    p.add_argument('--wall-overflow-capture-fraction',type=float,default=0.0)
    p.add_argument('--wall-overflow-return-rate-s-inv',type=float,default=0.0)
    p.add_argument('--wall-overflow-side-gate-length-m',type=float,default=0.0)
    p.add_argument('--horizontal-diffusivity-m2-s',type=float,default=0.0)
    p.add_argument('--horizontal-diffusion-respect-walls',action='store_true')
    p.add_argument('--crosswind-diffusivity-factor',type=float,default=1.0)
    p.add_argument('--wake-deficit-strength',type=float,default=0.0)
    p.add_argument('--wake-length-m',type=float,default=1.0)
    p.add_argument('--wake-base-width-m',type=float,default=1.0)
    p.add_argument('--wake-spread-slope',type=float,default=0.0)
    p.add_argument('--wake-recirculation-strength',type=float,default=0.0)
    p.add_argument('--wake-meander-strength',type=float,default=0.0)
    p.add_argument('--wake-shedding-strouhal',type=float,default=0.2)
    p.add_argument('--wake-meander-phase-deg',type=float,default=0.0)
    p.add_argument('--core-meander-strength',type=float,default=0.0)
    p.add_argument('--core-meander-length-m',type=float,default=100.0)
    p.add_argument('--core-meander-width-m',type=float,default=50.0)
    p.add_argument('--cavity-shelter-strength',type=float,default=0.0)
    p.add_argument('--cavity-core-entrainment-reduction',type=float,default=0.0)
    p.add_argument('--cavity-entrainment-factor',type=float,default=1.0)
    p.add_argument('--cavity-vortex-rate-s-inv',type=float,default=0.0)
    p.add_argument('--cavity-vortex-strength',type=float,default=0.0)
    p.add_argument('--cavity-reservoir-capture-rate-s-inv',type=float,default=0.005)
    p.add_argument('--cavity-reservoir-return-rate-s-inv',type=float,default=0.005)
    p.add_argument('--cavity-reservoir-capture-mode',
                   type=str,
                   choices=('cavity_rate','wall_impact'),default='cavity_rate')
    p.add_argument('--cavity-reservoir-vertical-capture-mode',
                   choices=('column','below_wall_core'),default='column')
    p.add_argument('--cavity-reservoir-return-sigma-m',type=float,default=16.7)
    p.add_argument('--cavity-reservoir-return-radius-m',type=float,default=10.0)
    p.add_argument('--cavity-reservoir-return-offset-m',type=float,default=0.0)
    p.add_argument('--cavity-reservoir-return-crosswind-scale',type=float,default=1.0)
    p.add_argument('--cavity-reservoir-face-return-fraction',type=float,default=0.0)
    p.add_argument('--cavity-reservoir-capture-radius-m',type=float,default=0.0)
    p.add_argument('--cavity-reservoir-start-time-s',type=float,default=0.0)
    p.add_argument('--cavity-return-streamwise-bias',type=float,default=0.0)
    p.add_argument('--cavity-vertical-exponent',type=float,default=0.0)
    p.add_argument('--cavity-vertical-mixing-factor',type=float,default=0.0)
    p.add_argument('--cavity-backflow-strength',type=float,default=0.9)
    p.add_argument('--cavity-backflow-rate-s-inv',type=float,default=0.01)
    p.add_argument('--cavity-backflow-frame',choices=['cavity','legacy_x'],default='cavity')
    p.add_argument('--cavity-stagnation-strength',type=float,default=0.0)
    p.add_argument('--cavity-geometry',type=str,choices=['full_disk','downwind_semidisk','upwind_semidisk'],default='downwind_semidisk')
    p.add_argument('--cavity-radius-m',type=float,default=50.0)
    p.add_argument('--cavity-align-with-wind',action='store_true')
    p.add_argument('--exterior-lateral-diffusivity-factor',type=float,default=1.0)
    p.add_argument('--exterior-transport',choices=['ambient_wind','density_current'],default='ambient_wind')
    p.add_argument('--density-current-strength',type=float,default=1.0)
    p.add_argument('--handoff-interface-path',type=Path,default=None)
    p.add_argument('--interface-x-m',type=float,default=None)
    p.add_argument('--wind-speed-m-s',type=float,default=3.9);p.add_argument('--wind-angle-deg',type=float,default=6.1)
    p.add_argument('--layer-wind-speed-m-s',type=float,default=None)
    p.add_argument('--wall-height-m',type=float,default=5.0)
    p.add_argument('--wall-radius-m',type=float,default=50.0)
    p.add_argument('--gravity-scale',type=float,default=1.0)
    a=p.parse_args();run(a.output,h=a.cell_size_m,dt=a.time_step_s,source_dt=a.source_time_step_s,layer_wind=a.layer_wind,entrainment=a.entrainment,entrainment_coefficient=a.entrainment_coefficient,drag_coefficient=a.drag_coefficient,domain_mode=a.domain_mode,source_mass_scale=a.source_mass_scale,source_spread_scale=a.source_spread_scale,source_longitudinal_scale=a.source_longitudinal_scale,source_crosswind_scale=a.source_crosswind_scale,source_initial_velocity_scale=a.source_initial_velocity_scale,source_release_duration_s=a.source_release_duration_s,puff_shear_entrainment_coefficient=a.puff_shear_entrainment_coefficient,puff_ambient_turbulence_coefficient=a.puff_ambient_turbulence_coefficient,puff_top_shear_entrainment_coefficient=a.puff_top_shear_entrainment_coefficient,puff_top_ambient_turbulence_coefficient=a.puff_top_ambient_turbulence_coefficient,puff_side_shear_entrainment_coefficient=a.puff_side_shear_entrainment_coefficient,puff_side_ambient_turbulence_coefficient=a.puff_side_ambient_turbulence_coefficient,puff_gravity_spreading_coefficient=a.puff_gravity_spreading_coefficient,puff_form_drag_coefficient=a.puff_form_drag_coefficient,handoff_radius_m=a.handoff_radius_m,vertical_observation_profile=a.vertical_observation_profile,dense_vertical_exponent=a.dense_vertical_exponent,vertical_grounding_exponent=a.vertical_grounding_exponent,vertical_elevated_fraction=a.vertical_elevated_fraction,wall_overtop_exponent=a.wall_overtop_exponent,wall_blockage_fraction=a.wall_blockage_fraction,vertical_tail_extension_m=a.vertical_tail_extension_m,vertical_tail_diffusivity_m2_s=a.vertical_tail_diffusivity_m2_s,mass_conserving_vertical_profile=a.mass_conserving_vertical_profile,vertical_observation_mode=a.vertical_observation_mode,vertical_observation_bins=a.vertical_observation_bins,vertical_observation_component_diagnostics=a.vertical_observation_component_diagnostics,vertical_layers=a.vertical_layers,initial_vertical_gaussian=a.initial_vertical_gaussian,vertical_mixing_diffusivity_m2_s=a.vertical_mixing_diffusivity_m2_s,multilayer_entrainment_placement=a.multilayer_entrainment_placement,multilayer_horizontal_shear_fraction=a.multilayer_horizontal_shear_fraction,multilayer_obstacle_crosswind_shear_rate_s_inv=a.multilayer_obstacle_crosswind_shear_rate_s_inv,multilayer_obstacle_crosswind_decay_length_m=a.multilayer_obstacle_crosswind_decay_length_m,multilayer_relative_transport_interval_s=a.multilayer_relative_transport_interval_s,height_resolved_wall_flux=a.height_resolved_wall_flux,wall_tangential_bypass_fraction=a.wall_tangential_bypass_fraction,wall_tip_bypass_fraction=a.wall_tip_bypass_fraction,wall_tip_bypass_power=a.wall_tip_bypass_power,wall_overflow_capture_fraction=a.wall_overflow_capture_fraction,wall_overflow_return_rate_s_inv=a.wall_overflow_return_rate_s_inv,wall_overflow_side_gate_length_m=a.wall_overflow_side_gate_length_m,horizontal_diffusivity_m2_s=a.horizontal_diffusivity_m2_s,horizontal_diffusion_respect_walls=a.horizontal_diffusion_respect_walls,crosswind_diffusivity_factor=a.crosswind_diffusivity_factor,wake_deficit_strength=a.wake_deficit_strength,wake_length_m=a.wake_length_m,wake_base_width_m=a.wake_base_width_m,wake_spread_slope=a.wake_spread_slope,wake_recirculation_strength=a.wake_recirculation_strength,wake_meander_strength=a.wake_meander_strength,wake_shedding_strouhal=a.wake_shedding_strouhal,wake_meander_phase_deg=a.wake_meander_phase_deg,cavity_shelter_strength=a.cavity_shelter_strength,cavity_entrainment_factor=a.cavity_entrainment_factor,cavity_core_entrainment_reduction=a.cavity_core_entrainment_reduction,cavity_vertical_exponent=a.cavity_vertical_exponent,cavity_vertical_mixing_factor=a.cavity_vertical_mixing_factor,cavity_backflow_strength=a.cavity_backflow_strength,cavity_backflow_rate_s_inv=a.cavity_backflow_rate_s_inv,cavity_backflow_frame=a.cavity_backflow_frame,cavity_stagnation_strength=a.cavity_stagnation_strength,cavity_radius_m=a.cavity_radius_m,cavity_vortex_rate_s_inv=a.cavity_vortex_rate_s_inv,cavity_vortex_strength=a.cavity_vortex_strength,cavity_reservoir_capture_rate_s_inv=a.cavity_reservoir_capture_rate_s_inv,cavity_reservoir_return_rate_s_inv=a.cavity_reservoir_return_rate_s_inv,cavity_reservoir_capture_mode=a.cavity_reservoir_capture_mode,cavity_reservoir_vertical_capture_mode=a.cavity_reservoir_vertical_capture_mode,cavity_reservoir_return_sigma_m=a.cavity_reservoir_return_sigma_m,cavity_reservoir_return_radius_m=a.cavity_reservoir_return_radius_m,cavity_reservoir_return_offset_m=a.cavity_reservoir_return_offset_m,cavity_reservoir_return_crosswind_scale=a.cavity_reservoir_return_crosswind_scale,cavity_reservoir_face_return_fraction=a.cavity_reservoir_face_return_fraction,cavity_reservoir_capture_radius_m=a.cavity_reservoir_capture_radius_m,cavity_reservoir_start_time_s=a.cavity_reservoir_start_time_s,cavity_return_streamwise_bias=a.cavity_return_streamwise_bias,cavity_align_with_wind=a.cavity_align_with_wind,exterior_lateral_diffusivity_factor=a.exterior_lateral_diffusivity_factor,exterior_transport=a.exterior_transport,density_current_strength=a.density_current_strength,handoff_interface_path=a.handoff_interface_path,interface_x_m=a.interface_x_m,wind_speed_m_s=a.wind_speed_m_s,layer_wind_speed_m_s=a.layer_wind_speed_m_s,wind_angle_deg=a.wind_angle_deg,wall_height_m=a.wall_height_m,wall_radius_m=a.wall_radius_m,gravity_scale=a.gravity_scale,source_model=a.source_model,source_front_profile=a.source_front_profile,source_initial_radius_m=a.source_initial_radius_m,observation_species_mode=a.observation_species_mode,source_radial_damping_coefficient=a.source_radial_damping_coefficient,source_radial_damping_aspect_ratio=a.source_radial_damping_aspect_ratio,source_minimum_height_m=a.source_minimum_height_m,source_wind_scale=a.source_wind_scale)
