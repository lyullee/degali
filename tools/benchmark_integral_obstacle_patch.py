"""Cold-start 600 s integral/local-grid/integral runtime and closure audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np

from degali.addons.integral_obstacle_patch import (
    GaussianIntegralInlet, downstream_trace, solve_integral_obstacle_patch,
)
from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D
from degali.addons.resolved_obstacle_flow_3d import ResolvedObstacleFlow3DConfig
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


def run(*, duration=600., spacing=1., step=.2, output: Path, flow_steps=1000, obstacle=True):
    started = time.perf_counter()
    if not all(np.isfinite(value) and value>0 for value in (duration,spacing,step)):
        raise ValueError("duration, spacing and step must be finite and positive")
    if flow_steps < 10 or any(not np.isclose(size/spacing,round(size/spacing)) for size in (20.,12.,10.)):
        raise ValueError("grid must tile the patch and flow needs at least ten steps")
    # Global exterior x=0..100; only x=30..50,y=6..18,z=0..10 is gridded.
    origin = (30., 6., 0.)
    source = SourceRateSchedule((0., min(120., duration/2), duration), (.02, 0., 0.))
    inlet = GaussianIntegralInlet(source, 25., 2., .03, 12., 1., 1.5)
    zero = SourceRateSchedule((0., duration), (0., 0.))
    cells = tuple(int(round(size/spacing)) for size in (20., 12., 10.))
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(20.,12.,10.), cells=cells, wind_speed_m_s=2.,
        kinematic_viscosity_m2_s=.5, flow_steps=flow_steps,
        minimum_flow_steps=flow_steps, lattice_inlet_speed=.03,
        flow_average_last_steps=max(2, flow_steps//5),
        duration_s=duration, scalar_time_step_s=step, scalar_diffusivity_m2_s=.03,
        source_schedule=zero, source_position_m=(1.,6.,1.), source_sigma_m=1.5,
        obstacles=(AxisAlignedCuboid(6.,10.,4.,8.,0.,4.,'block'),) if obstacle else (),
        smagorinsky_constant=.1, effective_schmidt_number=.7,
        scalar_turbulent_diffusivity_model="molecular_plus_subgrid_eddy_viscosity",
        scalar_advection_scheme="muscl_van_leer",
        nonstreamwise_wall_boundary=(
            "ground_no_slip_lateral_top_free_slip" if obstacle
            else "all_symmetry_free_slip"
        ),
    )
    probabilities = inlet.area_probabilities(
        origin[1]+np.arange(cells[1]+1)*spacing, np.arange(cells[2]+1)*spacing,
    )
    fraction = float(probabilities.sum())
    def supply(start, dt, y_edges, z_edges):
        return inlet.mass_between(start, start+dt)*probabilities/(spacing**2*dt)
    receptors = tuple(FastObstacleReceptor3D(label, point) for label, point in (
        ('front',(35.,12.,2.)), ('roof',(38.,12.,5.)),
        ('side',(38.,9.,2.)), ('wake',(43.,12.,2.)),
    ))
    result = solve_integral_obstacle_patch(config, incoming_flux=supply,
                                         receptors=receptors, origin_m=origin)
    factor = 100*8.314462618*293.15/(101325*.00201588)
    traces = {key: value*factor for key,value in result.receptor_traces_kg_m3.items()}
    traces['downstream'] = factor*downstream_trace(
        result, point_m=(70.,12.,2.), wind_speed_m_s=2., diffusivity_m2_s=.03,
    )
    # Exterior bypass is independent, positive mass. It is not forced into the
    # patch to improve the score or make the packet ledger appear closed.
    transmitted = np.array([inlet.mass_between(0.,t) for t in result.time_s])
    released = np.array([source.mass_between(0.,t) for t in result.time_s])
    upstream_inventory = released-transmitted
    bypass = transmitted*(1-fraction)
    global_residual = (upstream_inventory+bypass+result.domain_mass_kg
                       +result.outlet_mass_kg+result.other_boundary_mass_kg-released)
    total_wall = time.perf_counter()-started
    arrival = {}
    for key, values in traces.items():
        hits = np.flatnonzero(values >= .1)
        arrival[key] = None if len(hits)==0 else float(result.time_s[hits[0]])
    report = {
        'schema':'degali.integral-local-patch-benchmark.v1',
        'duration_s':duration, 'cold_start_compute_runtime_s':total_wall,
        'runtime_target_s':30., 'runtime_target_passed':total_wall<=30.,
        'global_extent_m':[100.,24.,20.], 'patch_origin_m':list(origin),
        'patch_extent_m':list(config.domain_m), 'patch_cell_count':int(np.prod(cells)),
        'uniform_same_resolution_cell_count':int(100*24*20/spacing**3),
        'exterior_time_marched_grid_cells':0, 'source_mass_kg':source.released_mass_kg,
        'fraction_of_inlet_mass_entering_patch':fraction,
        'maximum_global_species_residual_kg':float(np.max(np.abs(global_residual))),
        'maximum_global_relative_species_residual':float(np.max(np.abs(global_residual))/max(released[-1],1e-30)),
        'patch_diagnostics':result.diagnostics,
        'sensor_peak_vol_percent':{k:float(v.max()) for k,v in traces.items()},
        'sensor_arrival_at_0p1_vol_percent_s':arrival,
        'sensor_dose_vol_percent_s':{k:float(np.trapezoid(v,result.time_s)) for k,v in traces.items()},
        'obstacle_present':obstacle,
        'four_vol_percent_ever_patch_volume_m3':float(np.count_nonzero((result.maximum_kg_m3*factor>=4)&~result.solid_mask)*spacing**3),
        'accuracy_comparison':'not yet paired with independent experiment or CFD',
        'limitations':[
            'passive isothermal H2; cryogenic buoyancy not represented',
            'fixed wind and fixed resolved carrier; no transient momentum feedback',
            'exterior uses analytic Gaussian transport, not the DEGADIS dense-gas integral closure',
            'lateral/top patch exits are inventoried but not continued to downstream receptors',
            'patch outlet face packets are continued analytically; inlet-to-outlet reverse pressure feedback omitted',
        ],
    }
    if not obstacle:
        expected = {}
        for r in receptors:
            rx,ry,rz = r.position_m
            delay = (rx-5.)/2.
            sigma = np.sqrt(1.5**2+2*.03*delay)
            ky = np.exp(-.5*((ry-12.)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
            kz = (np.exp(-.5*((rz-1.)/sigma)**2)+np.exp(-.5*((rz+1.)/sigma)**2))/(np.sqrt(2*np.pi)*sigma)
            expected[r.label] = float(.02/2.*ky*kz*factor)
        report['no_obstacle_analytic_peak_vol_percent'] = expected
        report['no_obstacle_peak_relative_error_percent'] = {
            key:float(100*abs(traces[key].max()-value)/value) for key,value in expected.items()
        }
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'), time_s=result.time_s,
        x_m=result.x_m,y_m=result.y_m,z_m=result.z_m,solid_mask=result.solid_mask,
        maximum_h2_vol_percent=result.maximum_kg_m3*factor,
        final_h2_vol_percent=result.final_kg_m3*factor,
        sensor_labels=np.asarray(list(traces)),sensor_vol_percent=np.asarray(list(traces.values())),
        species_mass_residual_kg=result.mass_residual_kg,
        outlet_rate_kg_s=result.outlet_rate_kg_s)
    report['runtime_including_output_s'] = time.perf_counter()-started
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--duration-s',type=float,default=600.)
    parser.add_argument('--spacing-m',type=float,default=1.)
    parser.add_argument('--time-step-s',type=float,default=.2)
    parser.add_argument('--flow-steps',type=int,default=1000)
    parser.add_argument('--no-obstacle',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(run(duration=args.duration_s,spacing=args.spacing_m,
        step=args.time_step_s,flow_steps=args.flow_steps,obstacle=not args.no_obstacle,
        output=args.output),indent=2))
