"""Same-data TI021 dose MAPE for the runtime-bounded integral/local path.

This is deliberately a new passive-transport diagnostic. It is not a cryogenic
H2 run or a density-buoyant R12 qualification. All source/wind/geometry/mixing
choices are set before the observed doses are used for scoring.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import math
import numpy as np
from scipy.stats import pearsonr, spearmanr

from degali.addons.integral_obstacle_patch import (
    solve_integral_obstacle_patch, downstream_trace, all_boundary_trace,
)
from degali.addons.resolved_obstacle_flow_3d import ResolvedObstacleFlow3DConfig
from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud, initialise_regular_polygon_prism_cloud,
)
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import CircularArcWall
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
from degali.addons.binary_integral_puff import initial_prism_puff, binary_puff_traces, axisymmetric_puff_inventory
from degali.addons.finite_puff import GaussianPuffConfig
from dataclasses import replace

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'outputs/thorney-ti021-wide-nested-145-sensor-score-600s-x450y640-centered-2026-10-13-v2/ti021-nested-145-sensor-score.json'


def _boundary_packet_moments(result):
    """Return conservative space/time moments of positive exported packets."""
    h=float(result.diagnostics['cell_size_m'])
    dt=np.diff(result.time_s)
    t_mid=.5*(result.time_s[:-1]+result.time_s[1:])
    x,y,z=result.x_m,result.y_m,result.z_m
    zy_z,zy_y=np.meshgrid(z,y,indexing='ij')
    zx_z,zx_x=np.meshgrid(z,x,indexing='ij')
    yx_y,yx_x=np.meshgrid(y,x,indexing='ij')
    coordinates={
        'x_min':np.column_stack((np.full(zy_y.size,x[0]-h/2),zy_y.ravel(),zy_z.ravel())),
        'x_max':np.column_stack((np.full(zy_y.size,x[-1]+h/2),zy_y.ravel(),zy_z.ravel())),
        'y_min':np.column_stack((zx_x.ravel(),np.full(zx_x.size,y[0]-h/2),zx_z.ravel())),
        'y_max':np.column_stack((zx_x.ravel(),np.full(zx_x.size,y[-1]+h/2),zx_z.ravel())),
        'z_min':np.column_stack((yx_x.ravel(),yx_y.ravel(),np.full(yx_x.size,z[0]-h/2))),
        'z_max':np.column_stack((yx_x.ravel(),yx_y.ravel(),np.full(yx_x.size,z[-1]+h/2))),
    }
    summary={}
    for face,raw in result.boundary_rate_kg_s.items():
        mass=np.maximum(np.asarray(raw,float),0.)*dt[:,None,None]
        flat=mass.reshape(len(dt),-1)
        total=float(flat.sum())
        if total<=0.:
            summary[face]={'positive_exported_mass_kg':0.}
            continue
        cell_mass=flat.sum(axis=0)
        time_mass=flat.sum(axis=1)
        mean_time=float(np.dot(time_mass,t_mid)/total)
        mean_xyz=np.sum(coordinates[face]*cell_mass[:,None],axis=0)/total
        std_xyz=np.sqrt(np.sum(
            (coordinates[face]-mean_xyz[None,:])**2*cell_mass[:,None],axis=0
        )/total)
        vertical_partitions={}
        for name,mask in (
                ('at_or_below_5m',coordinates[face][:,2]<=5.0),
                ('above_5m',coordinates[face][:,2]>5.0)):
            partition_mass=float(cell_mass[mask].sum())
            if partition_mass<=0.0:
                vertical_partitions[name]={'positive_exported_mass_kg':0.0}
                continue
            partition_xyz=np.sum(
                coordinates[face][mask]*cell_mass[mask,None],axis=0
            )/partition_mass
            partition_time_mass=flat[:,mask].sum(axis=1)
            partition_time=float(
                np.dot(partition_time_mass,t_mid)/partition_mass)
            vertical_partitions[name]={
                'positive_exported_mass_kg':partition_mass,
                'mass_fraction_of_face':partition_mass/total,
                'mean_patch_exit_time_s':partition_time,
                'mass_centroid_xyz_m':[float(v) for v in partition_xyz],
            }
        summary[face]={
            'positive_exported_mass_kg':total,
            'mean_patch_exit_time_s':mean_time,
            'exit_time_standard_deviation_s':float(np.sqrt(np.dot(time_mass,(t_mid-mean_time)**2)/total)),
            'mass_centroid_xyz_m':[float(v) for v in mean_xyz],
            'mass_standard_deviation_xyz_m':[float(v) for v in std_xyz],
            'mass_fraction_above_5m':float(
                cell_mass[coordinates[face][:,2]>5.0].sum()/total),
            'vertical_partitions':vertical_partitions,
        }
    return summary


def run(*, output:Path, cell_size=5., flow_steps=1000, source_model='initial_prism',
        include_obstacle=True):
    started=time.perf_counter()
    workbook=ROOT/'ti21_input_ascii.xls'
    dataset=read_thorney_ti021_workbook(workbook)
    # Fixed report-coordinate patch enclosing the 50 m radius fence and the
    # holder. The distant radial sensor arcs remain outside the time-marched grid.
    origin=np.array([40.,250.,0.]); extent=(90.,140.,40.)
    h=float(cell_size)
    if any(not np.isclose(l/h,round(l/h)) for l in extent):
        raise ValueError('cell size must tile the patch')
    cells=tuple(round(l/h) for l in extent)
    heading=np.deg2rad(6.1)
    speed=3.9*np.cos(heading); lateral=3.9*np.sin(heading)
    u_star=.205; kappa=.4; z0=.005; schmidt=.7
    profile_z=(0.,.5,1.,2.,5.,10.,20.,40.)
    profile_u=tuple(u_star/kappa*np.log((z+z0)/z0)*np.cos(heading) for z in profile_z)
    zero=SourceRateSchedule((0.,600.),(0.,0.))
    wall=CircularArcWall((20.,70.),50.,0.,0.,5.,-90.,90.,'TI021 50m/5m fence')
    config=ResolvedObstacleFlow3DConfig(
        domain_m=extent,cells=cells,wind_speed_m_s=float(speed),
        # Mid-height neutral eddy viscosity, declared from weather rather
        # than selected using observed concentrations. Molecular-only D3Q19
        # at this coarse high-Re fence would not represent atmospheric mixing.
        kinematic_viscosity_m2_s=kappa*u_star*(extent[2]/2),
        flow_steps=flow_steps,minimum_flow_steps=flow_steps,
        flow_average_last_steps=flow_steps//5,lattice_inlet_speed=.03,
        duration_s=600.,scalar_time_step_s=.25,scalar_diffusivity_m2_s=9.788649302046554e-6,
        source_schedule=zero,source_position_m=(20.,70.,1.),source_sigma_m=h,
        obstacles=((wall,) if include_obstacle else ()),
        circular_arc_representation='zero_thickness_surface',
        inlet_profile_z_m=profile_z,inlet_profile_u_m_s=profile_u,
        ambient_temperature_k=291.15,ambient_pressure_pa=100500.,
        smagorinsky_constant=.1,effective_schmidt_number=schmidt,
        scalar_turbulent_diffusivity_model='molecular_plus_subgrid_eddy_viscosity',
        nonstreamwise_wall_boundary='ground_free_slip_lateral_top_pressure_open',
    )
    x=(np.arange(cells[0])+.5)*h
    y=(np.arange(cells[1])+.5)*h
    z=(np.arange(cells[2])+.5)*h
    r12_mass_fraction=(dataset.initial_r12_mole_fraction*.120913
                       /dataset.effective_molecular_weight_kg_mol)
    released_mixture_cp=(r12_mass_fraction*606.0748985077907
                         +(1.0-r12_mass_fraction)*1040.0)
    closure=BinaryIdealGasMixtureClosure(
        species_name='R12/N2 release mixture',
        species_molar_mass_kg_mol=dataset.effective_molecular_weight_kg_mol,
        species_heat_capacity_j_kg_k=released_mixture_cp,
        ambient_temperature_k=291.15,ambient_pressure_pa=100500.,minimum_temperature_k=150.)
    cloud=RegularPolygonPrismInitialCloud((20.,70.),14.,0.,13.,1.0,291.15)
    initial=initialise_regular_polygon_prism_cloud(x,y,z,cloud=cloud,closure=closure)
    positions=np.array([(60.+s.local_y_m,320.-s.local_x_m,s.height_m) for s in dataset.dose_sensors])
    handoff_time=0.; exterior_handoff_mass=0.; initial_field=initial.species_mass_concentration_kg_m3
    puff=None; contact_index=None
    if source_model=='integral_until_fence':
        puff=initial_prism_puff(RegularPolygonPrismInitialCloud((60.,320.),14.,0.,13.,1.0,291.15),
            closure,GaussianPuffConfig(duration_s=600.,time_step_s=.25,
                wind_velocity_m_s=(float(speed),float(lateral),0.),width_cutoff_sigma=math.sqrt(3)))
        # Declared footprint contact, not observed arrival or dose. Gaussian
        # tails have no finite front, so this is explicitly a geometric gate.
        contact_index=next(i for i,s in enumerate(puff.states) if
            np.linalg.norm(np.array(s.centre_position_m[:2])-np.array([60.,320.]))
            +max(s.longitudinal_half_length_m,s.lateral_half_width_m)>=50.)
        handoff_time=puff.states[contact_index].elapsed_s
        initial_field,exterior_handoff_mass=axisymmetric_puff_inventory(puff,contact_index,
            origin[0]+np.arange(cells[0]+1)*h,origin[1]+np.arange(cells[1]+1)*h,
            np.arange(cells[2]+1)*h)
        config=replace(config,duration_s=600.-handoff_time,
            source_schedule=SourceRateSchedule((0.,600.-handoff_time),(0.,0.)))
    elif source_model!='initial_prism':
        raise ValueError('unsupported source model')
    local=(positions>=origin).all(axis=1)&(positions<=origin+extent).all(axis=1)
    indices=np.flatnonzero(local)
    receptors=tuple(FastObstacleReceptor3D(f'TI021-S{i+1:03d}',tuple(positions[i])) for i in indices)
    additional=np.broadcast_to(kappa*u_star*z[:,None,None]/schmidt,(cells[2],cells[1],cells[0])).copy()
    result=solve_integral_obstacle_patch(
        config,incoming_flux=lambda t,dt,ye,ze: np.zeros((len(ze)-1,len(ye)-1)),
        receptors=receptors,origin_m=tuple(origin),
        initial_species_mass_kg_m3=initial_field,
        additional_diffusivity_m2_s=additional,
        carrier_velocity_offset_m_s=(0.,float(lateral),0.),
    )
    # GAS/HGAS report the released R12/N2 mixture concentration. Treat that
    # registered mixture as one pseudo-species, matching the production layer
    # scorer instead of comparing an R12 component concentration with it.
    conversion=(100*8.314462618*291.15
                /(100500*dataset.effective_molecular_weight_kg_mol))
    traces=np.zeros((145,len(result.time_s)))
    route=[]
    for i,point in enumerate(positions):
        label=f'TI021-S{i+1:03d}'
        if local[i]:
            traces[i]=result.receptor_traces_kg_m3[label]*conversion
            route.append('local_patch')
        else:
            # Continue every positive x/y/z boundary packet.  The old route
            # used x_max only and silently set lateral/top sensors to zero.
            # This operator keeps the height-resolved packet inventory and
            # applies one declared exterior Gaussian kernel; it is still a
            # diagnostic continuation, not a recoupled outer CFD grid.
            d=kappa*u_star*max(float(point[2]),h/2)/schmidt
            traces[i]=all_boundary_trace(
                result, point_m=tuple(point),
                wind_xy_m_s=(float(speed), float(lateral)),
                horizontal_diffusivity_m2_s=d,
                vertical_diffusivity_m2_s=d,
            )*conversion
            route.append('all_face_boundary_continuation')
    observation_times=result.time_s
    if puff is not None:
        before=binary_puff_traces(puff,positions)[:,:contact_index]
        before_times=np.array([s.elapsed_s for s in puff.states[:contact_index]])
        observation_times=np.concatenate((before_times,result.time_s+handoff_time))
        traces=np.concatenate((before,traces),axis=1)
    prediction_wall=time.perf_counter()-started
    # The observed values are only consulted after the transient is complete.
    observed=np.array([s.dose_percent_s for s in dataset.dose_sensors])
    predicted=np.trapezoid(traces,observation_times,axis=1)
    positive=observed>0
    ape=np.abs(predicted[positive]/observed[positive]-1)*100
    old=json.loads(OLD.read_text(encoding='utf-8'))
    old_mape=100*old['dose_sensor_score']['mean_absolute_percentage_error']
    groups={}
    for key in set(route):
        mask=positive&np.array([v==key for v in route])
        groups[key]={'positive_observed_sensor_count':int(mask.sum()),
            'dose_mape_percent':float(np.mean(np.abs(predicted[mask]/observed[mask]-1))*100) if mask.any() else None}
    report={
        'schema':'degali.ti021-integral-patch-30s-dose-score.v1',
        'gas':'R12/N2 release mixture','duration_s':600.,
        'observation_species_mode':'released_mixture',
        'metric':'sensor dose MAPE on positive observed released-mixture doses',
        'unit':'vol% s; equivalently average vol% over a common 600 s horizon',
        'not_an_instantaneous_concentration_time_series_mape':True,
        'sensor_count':145,'positive_observed_sensor_count':int(positive.sum()),
        'old_dose_mape_percent':old_mape,'new_dose_mape_percent':float(ape.mean()),
        'mape_change_percentage_points':float(ape.mean()-old_mape),
        'dose_mae_vol_percent_s':float(np.mean(np.abs(predicted[positive]-observed[positive]))),
        'equivalent_600s_average_concentration_mae_vol_percent':float(np.mean(np.abs(predicted[positive]-observed[positive]))/600),
        'normalized_mean_bias':float(np.sum(predicted[positive]-observed[positive])/observed[positive].sum()),
        'pearson_r':float(pearsonr(observed[positive],predicted[positive]).statistic),
        'spearman_rho':float(spearmanr(observed[positive],predicted[positive]).statistic),
        'routes':groups,'cell_size_m':h,'patch_cells':list(cells),
        'prediction_runtime_s':prediction_wall,'prediction_30s_target_passed':prediction_wall<=30,
        'patch_diagnostics':result.diagnostics,
        'boundary_packet_moments':_boundary_packet_moments(result),
        'score_status':('diagnostic_passive_model_not_experimentally_qualified'
                        if result.diagnostics['carrier_flow_converged']
                        else 'diagnostic_carrier_not_converged'),
        'initial_cloud_volume_m3':initial.represented_geometry_volume_m3,
        'initial_cloud_species_mass_kg':initial.integrated_species_mass_kg,
        'configuration_used_observed_dose_peak_arrival_values':False,
        'source_model':source_model,'integral_to_patch_handoff_time_s':handoff_time,
        'include_obstacle':bool(include_obstacle),
        'unmapped_exterior_species_mass_at_handoff_kg':exterior_handoff_mass,
        'handoff_species_ledger_residual_kg':float(initial_field.sum()*h**3+exterior_handoff_mass-initial.integrated_species_mass_kg),
        'handoff_density_momentum_enthalpy_transferred':False,
        'carrier_effective_viscosity_m2_s':config.kinematic_viscosity_m2_s,
        'carrier_viscosity_reference_height_m':extent[2]/2,
        'workbook_path':str(workbook),'workbook_sha256':hashlib.sha256(workbook.read_bytes()).hexdigest(),
        'baseline_score_path':str(OLD),
        'changed_physics':[
            'passive released R12/N2 mixture; dense-gas gravity spreading and evolving density omitted',
            'resolved time-averaged local carrier; no transient momentum/pressure feedback',
            'neutral background scalar mixing and SGS retained',
            'outer continuation uses constant per-receptor diffusivity, not a density-dependent integral model',
            'impermeable free-slip ground paired with scalar ground; not a no-slip surface-layer validation',
            'all positive x/y/z boundary packets exported to a one-way exterior continuation',
            ('initial quiescent dense integral puff until footprint contact; species-only cell-integrated handoff'
             if puff is not None else 'source is ambient-carried immediately rather than initially quiescent'),
            'carrier viscosity uses a bulk 20 m-height neutral closure, not the original molecular/variable-density momentum path',
        ],
        'sensor_rows':[{'label':f'TI021-S{i+1:03d}','position_m':point.tolist(),
            'route':route[i],'observed_dose_vol_percent_s':float(observed[i]),
            'predicted_dose_vol_percent_s':float(predicted[i]),
            'absolute_percentage_error':float(abs(predicted[i]/observed[i]-1)*100) if observed[i]>0 else None}
            for i,point in enumerate(positions)],
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),time_s=observation_times,
        labels=np.array([f'TI021-S{i+1:03d}' for i in range(145)]),
        predicted_concentration_vol_percent=traces,
        observed_dose_vol_percent_s=observed,predicted_dose_vol_percent_s=predicted,
        mass_residual_kg=result.mass_residual_kg)
    report['runtime_including_scoring_and_output_s']=time.perf_counter()-started
    report['end_to_end_30s_target_passed']=report['runtime_including_scoring_and_output_s']<=30
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='sensor_rows'},indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cell-size-m',type=float,default=5.)
    p.add_argument('--flow-steps',type=int,default=1000)
    p.add_argument('--source-model',choices=['initial_prism','integral_until_fence'],default='initial_prism')
    p.add_argument('--no-obstacle',action='store_true',
                   help='remove the 50 m arc wall while keeping every other input unchanged')
    a=p.parse_args()
    run(output=a.output,cell_size=a.cell_size_m,flow_steps=a.flow_steps,
        source_model=a.source_model,include_obstacle=not a.no_obstacle)
