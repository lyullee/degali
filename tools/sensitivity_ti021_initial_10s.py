"""Short-window OFAT sensitivities; never compare a 10s dose to a full dose.

Current integral/local hybrid, plus a separate compact-holder local-layer
structural control. Initial horizon is exactly ten physical seconds. Observed
data provide arrival summaries, NOT concentration traces in that interval.
The earliness score is a one-sided, censored timing consistency diagnostic.
"""
from dataclasses import replace
from pathlib import Path
import argparse
import hashlib
import json
import math
import time
import numpy as np
from scipy.special import ndtr

from degali.addons.binary_integral_puff import initial_prism_puff,binary_puff_traces,BinaryPuffThermodynamics
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_puff import GaussianPuffConfig
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud,initialise_regular_polygon_prism_cloud
from degali.addons.dense_layer_patch import solve_dense_layer_patch,_primitive
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
try:
    from tools.score_ti021_dense_layer_hybrid import curved_fence_faces
except ModuleNotFoundError:
    from score_ti021_dense_layer_hybrid import curved_fence_faces

ROOT=Path(__file__).resolve().parents[1]


def censored_timing_score(times,traces,arrivals,*,threshold=.1,endpoint=10.):
    """No detection by 10s is censored, not a correct later arrival prediction."""
    t=np.asarray(times,float);c=np.asarray(traces,float);a=np.asarray(arrivals,float)
    valid=np.isfinite(a)
    hit=c>=threshold
    first=np.argmax(hit,axis=1)
    prediction=np.where(hit.any(axis=1),t[first],np.nan)
    early=valid&np.isfinite(prediction)&(prediction<a)&(prediction<=endpoint)
    earliness=np.where(early,np.minimum(a,endpoint)-np.nan_to_num(prediction),0.)
    # Only compare within the observed pre-arrival interval; the large measured
    # arrival later than 10s must not inflate a short-window timing penalty.
    earliness=np.maximum(earliness,0.)
    dt=np.diff(t);mid=.5*(t[1:]+t[:-1])
    above=np.maximum((c[:,1:]+c[:,:-1])*.5-threshold,0.)
    mask=valid[:,None]&(mid[None,:]<np.minimum(a,endpoint)[:,None])
    excess=np.sum(above*mask*dt,axis=1)
    return {'threshold_vol_percent':threshold,'known_arrival_sensor_count':int(valid.sum()),
        'early_detected_sensor_count':int(np.sum(early&(earliness>0))),
        'sum_censored_earliness_s':float(earliness.sum()),
        'mean_censored_earliness_s_over_known_arrivals':float(earliness[valid].mean()) if valid.any() else None,
        'prearrival_threshold_excess_vol_percent_s':float(excess.sum()),
        'first_detection_time_s':[float(v) if np.isfinite(v) else None for v in prediction],
        'not_a_full_arrival_mae':True,'concentration_measurements_in_window_available':False}


def bounded_arrival_window_score(times,traces,arrivals,*,threshold=.1,endpoint=None):
    """Two-sided, horizon-bounded timing diagnostic.

    Unlike :func:`censored_timing_score`, this reports a lower-bound absolute
    timing error for sensors whose observed arrival lies inside the simulated
    horizon but whose model signal is still absent at the horizon.  It must
    not be called a full-event MAE: missing detections contribute only the
    known minimum delay ``endpoint - observed_arrival``.
    """
    t=np.asarray(times,float);c=np.asarray(traces,float);a=np.asarray(arrivals,float)
    if t.ndim!=1 or c.ndim!=2 or c.shape[1]!=len(t) or a.ndim!=1 or c.shape[0]!=len(a):
        raise ValueError('times, traces and arrivals must be aligned')
    if np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0.) or np.any(~np.isfinite(c)):
        raise ValueError('times and traces must be finite and increasing')
    horizon=float(t[-1] if endpoint is None else endpoint)
    if not math.isfinite(horizon) or horizon<=0. or horizon>float(t[-1])+1.0e-12:
        raise ValueError('endpoint must be finite, positive and within times')
    if not math.isfinite(threshold) or threshold<0.:
        raise ValueError('threshold must be finite and non-negative')
    valid=np.isfinite(a)&(a<=horizon)
    hit=(c>=threshold)&(t[None,:]<=horizon+1.0e-12)
    first=np.argmax(hit,axis=1)
    prediction=np.where(hit.any(axis=1),t[first],np.nan)
    in_window=valid&np.isfinite(prediction)&(prediction<=horizon+1.0e-12)
    missing=valid&~in_window
    abs_error=np.where(in_window,np.abs(prediction-a),
                       np.where(missing,np.maximum(horizon-a,0.),0.))
    return {
        'threshold_vol_percent':float(threshold),
        'horizon_s':horizon,
        'observed_arrivals_in_window_count':int(valid.sum()),
        'predicted_arrivals_in_window_count':int(in_window.sum()),
        'missing_by_horizon_count':int(missing.sum()),
        'absolute_timing_error_lower_bound_s':float(abs_error[valid].sum()),
        'mean_absolute_timing_error_lower_bound_s':float(abs_error[valid].mean()) if valid.any() else None,
        'first_detection_time_s':[float(v) if np.isfinite(v) else None for v in prediction],
        'not_a_full_arrival_mae':True,
        'concentration_measurements_in_window_available':False,
    }


def relative_peak_timing_score(times, traces, arrivals, *,
                               relative_threshold=.05, endpoint=10.):
    """Apply the TI021 5%-of-local-peak arrival operator to model traces.

    The workbook does not provide the measured sensor peaks, so this function
    uses each model trace's own peak only to make the *operator* compatible.
    It remains a censored timing diagnostic, not a concentration validation or
    an arrival-MAE estimate.
    """
    t=np.asarray(times,float);c=np.asarray(traces,float);a=np.asarray(arrivals,float)
    if t.ndim!=1 or c.ndim!=2 or c.shape[1]!=len(t) or np.any(~np.isfinite(t)) or np.any(~np.isfinite(c)):
        raise ValueError('times and traces must be finite and aligned')
    if not 0.0 < relative_threshold < 1.0 or endpoint <= 0.0:
        raise ValueError('relative threshold and endpoint must be positive')
    peaks=np.max(c,axis=1)
    threshold=relative_threshold*peaks
    hit=(c>=threshold[:,None])&(peaks[:,None]>0.0)
    first=np.argmax(hit,axis=1)
    prediction=np.where(hit.any(axis=1),t[first],np.nan)
    valid=np.isfinite(a)
    early=valid&np.isfinite(prediction)&(prediction<a)&(prediction<=endpoint)
    earliness=np.where(early,np.minimum(a,endpoint)-np.nan_to_num(prediction),0.)
    earliness=np.maximum(earliness,0.)
    dt=np.diff(t);mid=.5*(t[1:]+t[:-1])
    above=np.maximum((c[:,1:]+c[:,:-1])*.5-threshold[:,None],0.)
    mask=valid[:,None]&(mid[None,:]<np.minimum(a,endpoint)[:,None])
    excess=np.sum(above*mask*dt,axis=1)
    return {'relative_peak_threshold':float(relative_threshold),
        'known_arrival_sensor_count':int(valid.sum()),
        'early_detected_sensor_count':int(np.sum(early&(earliness>0))),
        'sum_censored_earliness_s':float(earliness.sum()),
        'mean_censored_earliness_s_over_known_arrivals':float(earliness[valid].mean()) if valid.any() else None,
        'prearrival_relative_excess_vol_percent_s':float(excess.sum()),
        'first_detection_time_s':[float(v) if np.isfinite(v) else None for v in prediction],
        'model_local_peak_vol_percent':[float(v) for v in peaks],
        'not_a_full_arrival_mae':True,
        'concentration_measurements_in_window_available':False,
        'measured_local_peak_available':False}


def compact_cylinder_traces(puff,points):
    """Observation-only structural control, matched inventory and 2nd moments."""
    out=np.zeros((len(points),len(puff.states)))
    for i,s in enumerate(puff.states):
        sigma=s.longitudinal_half_length_m/puff.config.width_cutoff_sigma
        radius=2*sigma;depth=2*s.vertical_half_width_m
        inside=(np.sum((points[:,:2]-np.array(s.centre_position_m[:2]))**2,axis=1)<=radius**2)&(points[:,2]<=depth)
        volume=np.pi*radius**2*depth
        density=s.hydrogen_mass_kg/volume
        mole=100*density*8.314462618*s.bulk_temperature_k/(puff.thermodynamics.closure.ambient_pressure_pa*puff.thermodynamics.fuel_molecular_weight)
        if mole>100+1e-9:raise RuntimeError('compact cylinder observation exceeds pure species')
        out[:,i]=np.where(inside,mole,0.)
    return out


def sensor_average(times,traces,window_s=.6):
    """Causal boxcar with zero pre-event gas, exact piecewise-linear integral."""
    t=np.asarray(times,float);c=np.asarray(traces,float)
    if window_s<=0:raise ValueError('positive averaging window required')
    cumulative=np.concatenate((np.zeros((len(c),1)),np.cumsum(.5*(c[:,1:]+c[:,:-1])*np.diff(t),axis=1)),axis=1)
    earlier=t-window_s
    previous=np.array([np.interp(earlier,t,row,left=0.) for row in cumulative])
    return (cumulative-previous)/window_s


def make_cases():
    cases=[('hybrid_baseline','hybrid',{},{}),
        ('gravity_off','hybrid',{'gravity_spreading_coefficient':0.},{}),
        ('gravity_half','hybrid',{'gravity_spreading_coefficient':.375},{}),
        ('gravity_double','hybrid',{'gravity_spreading_coefficient':1.5},{}),
        ('shear_entrainment_off','hybrid',{'shear_entrainment_coefficient':0.},{}),
        ('shear_entrainment_half','hybrid',{'shear_entrainment_coefficient':.05},{}),
        ('shear_entrainment_double','hybrid',{'shear_entrainment_coefficient':.2},{}),
        ('ambient_entrainment_off','hybrid',{'ambient_turbulence_coefficient':0.},{}),
        ('ambient_entrainment_half','hybrid',{'ambient_turbulence_coefficient':.05},{}),
        ('ambient_entrainment_double','hybrid',{'ambient_turbulence_coefficient':.2},{}),
        ('all_entrainment_off','hybrid',{'shear_entrainment_coefficient':0.,'ambient_turbulence_coefficient':0.},{}),
        ('form_drag_off','hybrid',{'form_drag_coefficient':0.},{}),
        ('form_drag_half','hybrid',{'form_drag_coefficient':.1},{}),
        ('form_drag_double','hybrid',{'form_drag_coefficient':.4},{}),
        ('ri_entrainment_suppression','hybrid',{'dense_entrainment_model':'parker_bulk_suppression'},{}),
        ('source_dt_half','hybrid',{'time_step_s':.0125},{}),
        ('compact_observer_only','puff_observer',{},{}),
        ('compact_observer_gravity_half','puff_observer',{'gravity_spreading_coefficient':.375},{}),
        ('compact_observer_gravity_off','puff_observer',{'gravity_spreading_coefficient':0.},{}),
        ('puff_without_patch','puff_observer_gaussian',{},{}),
        ('layer_fence_off','hybrid',{}, {'fence_off':True}),
        ('layer_entrainment_off_after_handoff','hybrid',{}, {'entrainment_coefficient':0.}),
        ('layer_drag_off_after_handoff','hybrid',{}, {'drag_coefficient':0.}),
        ('layer_gravity_off_after_handoff','hybrid',{}, {'gravity_scale':0.}),
        ('layer_uniform_wind_after_handoff','hybrid',{}, {'friction_velocity_m_s':None}),
        ('holder_layer_baseline','holder_layer',{},{}),
        ('holder_layer_gravity_off','holder_layer',{}, {'gravity_scale':0.}),
        ('holder_layer_gravity_half','holder_layer',{}, {'gravity_scale':.5}),
        ('holder_layer_gravity_double','holder_layer',{}, {'gravity_scale':2.}),
        ('holder_layer_entrainment_off','holder_layer',{}, {'entrainment_coefficient':0.}),
        ('holder_layer_entrainment_half','holder_layer',{}, {'entrainment_coefficient':.05}),
        ('holder_layer_entrainment_double','holder_layer',{}, {'entrainment_coefficient':.2}),
        ('holder_layer_drag_off','holder_layer',{}, {'drag_coefficient':0.}),
        ('holder_layer_drag_half','holder_layer',{}, {'drag_coefficient':.1}),
        ('holder_layer_drag_double','holder_layer',{}, {'drag_coefficient':.4}),
        ('holder_layer_ri_suppression','holder_layer',{}, {'dense_entrainment_model':'parker_bulk_suppression'}),
        ('holder_layer_fence_off','holder_layer',{}, {'fence_off':True}),
        ('holder_layer_uniform_wind','holder_layer',{}, {'friction_velocity_m_s':None})]
    return cases


def run(output,*,h=2.5):
    start=time.perf_counter()
    dataset=read_thorney_ti021_workbook(ROOT/'ti21_input_ascii.xls')
    closure=BinaryIdealGasMixtureClosure('R12',.120913,606.0748985077907,
        ambient_temperature_k=291.15,ambient_pressure_pa=100500.)
    thermo=BinaryPuffThermodynamics(closure)
    cloud=RegularPolygonPrismInitialCloud((60.,320.),14.,0.,13.,.32,291.15)
    angle=np.deg2rad(6.1);wind=np.array([3.9*np.cos(angle),3.9*np.sin(angle)])
    source=GaussianPuffConfig(duration_s=10.,time_step_s=.025,wind_velocity_m_s=(*wind,0.),width_cutoff_sigma=math.sqrt(3))
    points=np.array([(60.+s.local_y_m,320.-s.local_x_m,s.height_m) for s in dataset.dose_sensors])
    arrivals=np.array([s.arrival_time_s if s.arrival_time_s is not None else np.nan for s in dataset.dose_sensors])
    origin=(20.,250.);extent=(110.,140.)
    if any(not np.isclose(v/h,round(v/h)) for v in extent):raise ValueError('h must tile local extent')
    nx,ny=(round(v/h) for v in extent)
    xe=origin[0]+np.arange(nx+1)*h;ye=origin[1]+np.arange(ny+1)*h
    wx,wy=curved_fence_faces(xe,ye)
    # One 13 m vertical slice exactly holds the initial prism. The horizontal
    # geometry comes from the existing exact polygon/cell intersection map.
    init=initialise_regular_polygon_prism_cloud((xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2,
        np.array([6.5,19.5]),cloud=cloud,closure=closure)
    initial_state=closure.state_from_mole_fraction_and_temperature(.32,291.15)
    holder=np.zeros((4,ny,nx))
    holder[0]=init.volume_fraction[0]*13*float(initial_state.mixture_density_kg_m3)
    holder[1]=init.species_mass_concentration_kg_m3[0]*13
    records=[];arrays={};puff_cache={}
    for name,kind,source_changes,layer_changes in make_cases():
        case_start=time.perf_counter();layer=None;contact_time=None;outside=0.
        cfg=replace(source,**source_changes)
        key=repr(cfg)
        if key not in puff_cache:puff_cache[key]=initial_prism_puff(cloud,closure,cfg)
        puff=puff_cache[key]
        t=np.array([s.elapsed_s for s in puff.states]);trace=binary_puff_traces(puff,points)
        source_final=puff.states[-1]
        if kind=='puff_observer':trace=compact_cylinder_traces(puff,points)
        elif kind=='holder_layer':q=holder.copy();contact_time=0.
        elif kind=='hybrid':
            contact=next((i for i,s in enumerate(puff.states) if
                np.linalg.norm(np.array(s.centre_position_m[:2])-np.array([60.,320.]))+s.longitudinal_half_length_m>=50.),None)
            if contact is not None and t[contact]<10-1e-10:
                s=puff.states[contact];contact_time=t[contact]
                sigma=s.longitudinal_half_length_m/cfg.width_cutoff_sigma
                prob=np.diff(ndtr((ye-s.centre_position_m[1])/sigma))[:,None]*np.diff(ndtr((xe-s.centre_position_m[0])/sigma))[None,:]
                q=np.empty((4,ny,nx));q[0]=s.total_mass_kg*prob/h**2;q[1]=s.hydrogen_mass_kg*prob/h**2
                q[2]=q[0]*s.bulk_velocity_m_s[0];q[3]=q[0]*s.bulk_velocity_m_s[1]
                outside=s.hydrogen_mass_kg-float(q[1].sum()*h*h)
        if contact_time is not None:
            settings=dict(max_dt_s=.05,entrainment_coefficient=.1,drag_coefficient=.2,
                          friction_velocity_m_s=.205,roughness_length_m=.005)
            settings.update(layer_changes);fence_off=settings.pop('fence_off',False)
            layer=solve_dense_layer_patch(initial_conserved=q,thermo=thermo,origin_xy_m=origin,
                cell_size_m=h,duration_s=10.-contact_time,wind_xy_m_s=wind,
                wall_x_height_m=np.zeros_like(wx) if fence_off else wx,
                wall_y_height_m=np.zeros_like(wy) if fence_off else wy,
                receptor_points_m=points,**settings)
            if kind=='holder_layer':t=layer.time_s;trace=layer.species_traces_vol_percent
            else:
                trace=np.concatenate((trace[:,:contact],layer.species_traces_vol_percent),axis=1)
                t=np.concatenate((t[:contact],layer.time_s+contact_time))
        common=np.linspace(0.,10.,401)
        values=np.array([np.interp(common,t,row) for row in trace])
        scores={str(threshold):censored_timing_score(common,values,arrivals,threshold=threshold) for threshold in (.05,.1,.5)}
        averaged=sensor_average(common,values)
        averaged_score=censored_timing_score(common,averaged,arrivals)
        diagnostic=layer.diagnostics if layer is not None else {}
        record={'case':name,'protocol':kind,'source_changes':source_changes,'layer_changes':layer_changes,
            'source_coefficients':{k:getattr(cfg,k) for k in ('gravity_spreading_coefficient','shear_entrainment_coefficient',
                'ambient_turbulence_coefficient','form_drag_coefficient','dense_entrainment_model','time_step_s')},
            'runtime_s':time.perf_counter()-case_start,'physical_end_time_s':float(t[-1]),
            'integral_to_layer_time_s':float(contact_time) if contact_time is not None else None,
            'source_10s_centre_m':list(source_final.centre_position_m),
            'source_10s_half_widths_m':[source_final.longitudinal_half_length_m,source_final.lateral_half_width_m,source_final.vertical_half_width_m],
            'source_10s_density_kg_m3':source_final.bulk_density_kg_m3,
            'source_10s_entrained_air_kg':source_final.cumulative_entrained_air_kg,
            'source_initial_mass_residual':puff.maximum_relative_mass_residual,
            'source_species_mass_residual':puff.maximum_relative_hydrogen_residual,
            'scoring':scores,'causal_0p6s_averaged_scoring':averaged_score,
            'initial_detected_sensor_count_0p1':int(np.sum(values[:,0]>=.1)),
            'peak_vol_percent':float(np.max(values)),
            'sum_0to10_sensor_dose_vol_percent_s':float(np.trapezoid(values,common,axis=1).sum()),
            'layer_diagnostics':diagnostic,'handoff_unmapped_species_kg':outside,
            'prefix_is_right_censored':True,
            'concentration_error_mape':None,'full_event_dose_error_mape':None,
            'sensor_rows':[{'sensor':i+1,'height_m':float(points[i,2]),
                'observed_arrival_s':float(arrivals[i]) if np.isfinite(arrivals[i]) else None,
                'prediction_first_0p1vol_percent_s':scores['0.1']['first_detection_time_s'][i],
                'dose_0to10_vol_percent_s':float(np.trapezoid(values[i],common)),
                'concentration_at_0s_vol_percent':float(values[i,0]),
                'concentration_at_10s_vol_percent':float(values[i,-1])} for i in range(len(points))]}
        records.append(record);arrays[name]=values
    baselines={r['protocol']:r for r in records if r['case'] in ('hybrid_baseline','holder_layer_baseline')}
    for r in records:
        baseline=baselines.get(r['protocol'],baselines['hybrid'])
        r['comparison_baseline']=baseline['case']
        r['early_count_change']=r['scoring']['0.1']['early_detected_sensor_count']-baseline['scoring']['0.1']['early_detected_sensor_count']
        r['censored_earliness_mean_change_s']=r['scoring']['0.1']['mean_censored_earliness_s_over_known_arrivals']-baseline['scoring']['0.1']['mean_censored_earliness_s_over_known_arrivals']
    report={'schema':'degali.ti021-initial10s-ofat.v1','gas':'R12','physical_horizon_s':10.,
        'status':'short_window_diagnostic_not_parameter_calibration_or_validation',
        'case_count':len(records),'sensor_count':145,'known_arrival_count':int(np.isfinite(arrivals).sum()),
        'earliest_observed_arrival_s':float(np.nanmin(arrivals)),
        'observed_arrivals_before10s':int(np.sum(arrivals<10)),
        'observed_arrivals_at10s':int(np.sum(arrivals==10)),
        'instantaneous_concentration_mape_available':False,'full_event_dose_mape_computed':False,
        'threshold_protocol':'0.1 vol% nominal resolution proxy; .05 and .5 sensitivity; no proven experimental trigger equivalence',
        'censored_timing_meaning':'compare premature detection before min(observed arrival,10s); non-detection is censored, not a correct later arrival',
        'sensor_response_filter_applied':'raw primary plus causal .6s boxcar sensitivity','experimental_sensor_averaging_time_s':.6,
        'source_gaussian_vs_compact_is_structural_not_single_physics_ablation':True,
        'all_case_wall_clock_scope':'cache shared identical source integrations; no precomputed carrier; 10s max per case',
        'original_workbook_sha256':dataset.workbook_sha256,
        'cases':records}
    baseline_puff=puff_cache[repr(source)]
    initial_point=np.array([baseline_puff.states[0].centre_position_m])
    report['source_profile_audit']={
        'reported_initial_uniform_concentration_vol_percent':32.,
        'gaussian_initial_centre_concentration_vol_percent':float(binary_puff_traces(baseline_puff,initial_point,state_count=1)[0,0]),
        'compact_initial_centre_concentration_vol_percent':float(compact_cylinder_traces(baseline_puff,initial_point)[0,0]),
        'initial_full_depth_over_radius':13/cloud.apothem_m,
        'initial_bulk_density_ratio':baseline_puff.states[0].bulk_density_kg_m3/thermo.ambient_density,
        'initial_gravity_width_growth_speed_m_s':.75*math.sqrt(2*9.80665*(baseline_puff.states[0].bulk_density_kg_m3/thermo.ambient_density-1.)*baseline_puff.states[0].vertical_half_width_m),
        'source_gaussian_does_not_reproduce_uniform_holder_concentration':True,
        'initial_layer_is_not_a_slender_shallow_cloud':True}
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),time_s=np.linspace(0,10,401),**arrays)
    report['total_wall_time_including_output_s']=time.perf_counter()-start
    report['implementation_sha256']={str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
        Path('src/degali/addons/finite_puff.py'),Path('src/degali/addons/dense_layer_patch.py'),
        Path('src/degali/addons/binary_integral_puff.py'),Path('tools/sensitivity_ti021_initial_10s.py'))}
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(f"Cases={len(records)}, total wall={report['total_wall_time_including_output_s']:.3f}s")
    for r in records:
        s=r['scoring']['0.1'];print(f"{r['case']:42s} early={s['early_detected_sensor_count']:2d} mean={s['mean_censored_earliness_s_over_known_arrivals']:.3f}s S001={r['sensor_rows'][0]['prediction_first_0p1vol_percent_s']} S004={r['sensor_rows'][3]['prediction_first_0p1vol_percent_s']} air={r['source_10s_entrained_air_kg']:.1f}kg")
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cell-size-m',type=float,default=2.5)
    a=p.parse_args();run(a.output,h=a.cell_size_m)
