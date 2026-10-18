"""Read-only instrumentation of the failed fast TI021 execution.

Wrappers observe the existing solver calls without changing fields, boundary
operators, timestep, source, or scoring. Production solver files are not fixed
by this diagnostic. The report separates inventory loss from misplaced mass.
"""
from __future__ import annotations
import argparse
import contextlib
import io
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np

import degali.addons.integral_obstacle_patch as transport
import degali.addons.resolved_obstacle_flow_3d as flow
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
try:
    from tools.score_ti021_integral_patch_30s import run
except ModuleNotFoundError:
    from score_ti021_integral_patch_30s import run


def diagnose(*,steps:int,output:Path,sealed_ground=False,ground_free_slip=False):
    seed_solver=transport.solve_resolved_obstacle_flow_3d
    projector=transport._projected_faces_from_cell_velocity
    scalar_step=transport._transport_step
    outer_operator=flow._apply_outer_boundaries_d3q19
    raw_faces=flow._raw_faces_from_cell_velocity
    evidence={}
    ledger={face:0. for face in ('x_min','x_max','y_min','y_max','z_min','z_max')}
    elapsed=0.; snapshots=[]; frames={}
    selected=np.array([0.,1.,3.,5.,10.,20.,40.,80.,160.,320.,600.])
    next_snapshot=0
    def seed_observer(config,**kwargs):
        seed=seed_solver(config,**kwargs)
        evidence['configuration']={
            'nonstreamwise_wall_boundary':config.nonstreamwise_wall_boundary,
            'flow_boundary_scheme':config.flow_boundary_scheme,
            'cell_size_m':config.domain_m[0]/config.cells[0],
            'initial_momentum':'fast runner immediately uses ambient carrier',
            'scalar_open_faces':['x_min','x_max','y_min','y_max','z_max'],
        }
        evidence['carrier']={k:seed.diagnostics.get(k) for k in (
            'flow_converged','final_flow_residual','maximum_lattice_mach_number',
            'maximum_relative_lattice_density_variation',
            'minimum_lee_streamwise_velocity_m_s','recirculating_lee_cell_count')}
        evidence['z_centres_m']=seed.z_m.tolist()
        return seed
    def ground_boundary_probe(distribution,*args,**kwargs):
        if (sealed_ground or ground_free_slip) and kwargs.get('nonstreamwise_wall_boundary')=='all_pressure_open':
            kwargs['pressure_open_nonstream_faces']=('y_min','y_max','z_max')
        result=outer_operator(distribution,*args,**kwargs)
        if sealed_ground and not ground_free_slip:
            # Diagnostic intervention only: retain open top/sides and use
            # the solver's existing no-slip bounce-back rule at the ground.
            for q,(_,_,cz) in enumerate(flow._C3):
                if cz>0:
                    distribution[q,0,:,:]=distribution[flow._OPPOSITE3[q],0,:,:].copy()
        return result
    def ground_faces_probe(**kwargs):
        u,v,w,div=raw_faces(**kwargs)
        if sealed_ground or ground_free_slip:
            w[0]=0.
            h=kwargs['spacing_m']
            div=(np.diff(u,axis=2)+np.diff(v,axis=1)+np.diff(w,axis=0))/h
            div[kwargs['solid']]=0.
        return u,v,w,div
    def project_observer(**kwargs):
        projected=projector(**kwargs)
        u,v,w=projected[:3]
        def stats(a):
            return {'min':float(a.min()),'max':float(a.max()),'mean':float(a.mean())}
        evidence['actual_transport_faces']={
            'u_m_s':stats(u),'v_m_s':stats(v),'w_m_s':stats(w),
            'ground_normal_velocity_m_s':stats(w[0]),
            'roof_normal_velocity_m_s':stats(w[-1]),
            'maximum_abs_ground_velocity_m_s':float(np.max(np.abs(w[0]))),
            'maximum_abs_ground_scalar_flux_kg_m2_s':0.,
            'post_projection_divergence_s_1':projected[5],
        }
        return projected
    def scalar_observer(c,**kwargs):
        nonlocal elapsed,next_snapshot
        dt=kwargs['dt']; h=float(kwargs['spacing'][0]); vol=h**3
        z=np.array(evidence['z_centres_m'])
        def inventory(field,time):
            layers=field.sum(axis=(1,2))*vol
            total=float(layers.sum())
            low=float(layers[z<5].sum())
            return {'time_s':float(time),'domain_species_mass_kg':total,
                    'lowest_5m_species_mass_kg':low,
                    'lowest_5m_fraction':low/total if total>0 else 0.,
                    'mass_centroid_height_m':float(layers@z/total) if total>0 else None,
                    'layer_species_mass_kg':layers.tolist(),
                    'cumulative_exit_mass_kg':dict(ledger)}
        if next_snapshot==0:
            snapshots.append(inventory(c,0.));frames['t0']=c.copy();next_snapshot+=1
        callback=kwargs.get('face_flux_callback')
        def flux_observer(fx,fy,fz):
            ledger['x_min']-=float(fx[:,:,0].sum()*h*h*dt)
            ledger['x_max']+=float(fx[:,:,-1].sum()*h*h*dt)
            ledger['y_min']-=float(fy[:,0,:].sum()*h*h*dt)
            ledger['y_max']+=float(fy[:,-1,:].sum()*h*h*dt)
            ledger['z_min']-=float(fz[0].sum()*h*h*dt)
            ledger['z_max']+=float(fz[-1].sum()*h*h*dt)
            if callback is not None:
                callback(fx,fy,fz)
        kwargs['face_flux_callback']=flux_observer
        advanced,boundary_net=scalar_step(c,**kwargs)
        elapsed+=dt
        if next_snapshot<len(selected) and elapsed>=selected[next_snapshot]-1e-9:
            requested=selected[next_snapshot]
            snapshots.append(inventory(advanced,elapsed))
            if requested in (5.,20.,40.):frames[f't{int(requested)}']=advanced.copy()
            next_snapshot+=1
        return advanced,boundary_net
    run_output=output.with_name(output.stem+'-run.json')
    with patch.object(transport,'solve_resolved_obstacle_flow_3d',seed_observer), \
         patch.object(transport,'_projected_faces_from_cell_velocity',project_observer), \
         patch.object(transport,'_transport_step',scalar_observer), \
         patch.object(flow,'_apply_outer_boundaries_d3q19',ground_boundary_probe), \
         patch.object(flow,'_raw_faces_from_cell_velocity',ground_faces_probe), \
         contextlib.redirect_stdout(io.StringIO()):
        scored=run(output=run_output,flow_steps=steps)
    evidence['snapshots']=snapshots
    initial=scored['initial_cloud_species_mass_kg']
    evidence['final_exit_mass_kg']=ledger
    evidence['final_exit_fraction_of_initial']={k:v/initial for k,v in ledger.items()}
    evidence['dose_score']=scored['new_dose_mape_percent']
    evidence['initial_species_mass_kg']=initial
    evidence['workbook_sha256']=scored['workbook_sha256']
    evidence['run_json']=str(run_output.resolve())
    closure=BinaryIdealGasMixtureClosure(
        species_name='R12',species_molar_mass_kg_mol=.120913,
        species_heat_capacity_j_kg_k=606.0748985077907)
    x=.32
    mixture_mw=x*.120913+(1-x)*closure.carrier_molar_mass_kg_mol
    density=100500*mixture_mw/(8.314462618*291.15)
    ratio=density/closure.ambient_density_kg_m3
    acceleration=9.80665*(ratio-1)/ratio
    evidence['missing_density_physics_scale']={
        'initial_mole_fraction':x,'initial_mixture_density_kg_m3':density,
        'ambient_air_density_kg_m3':closure.ambient_density_kg_m3,
        'density_ratio':ratio,
        'nonboussinesq_vertical_acceleration_scale_m_s2':acceleration,
        'gravity_height_velocity_scale_m_s':float(np.sqrt(acceleration*13)),
        'initial_height_collapse_timescale_s':float(np.sqrt(2*13/acceleration)),
        'gravity_to_wind_inertia_ratio':acceleration*13/3.9**2,
        'force_applied_by_fast_model_m_s2':0.,
    }
    evidence['instrumentation_changed_solver_states_or_coefficients']=bool(sealed_ground or ground_free_slip)
    evidence['diagnostic_intervention']=(
        'impermeable_free_slip_ground_only' if ground_free_slip else
        ('sealed_no_slip_ground_only' if sealed_ground else 'none'))
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(evidence,indent=2,allow_nan=False),encoding='utf-8')
    np.savez_compressed(output.with_suffix('.npz'),**frames)
    print(json.dumps({k:v for k,v in evidence.items() if k not in ('snapshots',)},indent=2))
    return evidence


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--flow-steps',type=int,default=1000)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--sealed-ground',action='store_true')
    p.add_argument('--ground-free-slip',action='store_true');a=p.parse_args()
    diagnose(steps=a.flow_steps,output=a.output,sealed_ground=a.sealed_ground,
             ground_free_slip=a.ground_free_slip)
