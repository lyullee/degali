"""Separate actual stored-field mechanical/thermal budgets from point cold gaps."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import math
from pathlib import Path

import numpy as np

from stage2_matched_source import ROOT,REF,read,digest,write_new
from audit_stage2_flow_effect import replay_document
from audit_downstream_cold_envelope import reflected_peak
from audit_mixed_air_liquid import independent_split
from run_preslhy_ambient_profile_audit import replay
from degali.addons.axisymmetric_jet import _air_phase_property_table,R_UNIVERSAL
from stage1 import verify


def velocity_moment(model,state,power):
    """Integral rho*u**power over the ORIGINAL finite Gaussian section."""
    if power not in (1,2,3): raise ValueError('moment power must be1,2or3')
    rho,_,area,theta,u,_,_=model._physical(state)
    wind=model._wind(state)*math.cos(theta)
    lam=model.velocity_shape_exponent
    return area*sum(math.comb(power,j)*wind**(power-j)*u**j*
        (model.rhoa*model.k.profile_integral(j*lam)+(rho-model.rhoa)*model.k.profile_integral(1+j*lam))
        for j in range(power+1))


def legacy_forward(th,t,y):
    """Old independent-phase EOS at given T,Y; NOT the new mixed-liquid EOS."""
    if not math.isfinite(t) or not 14.1<=t<=300. or not math.isfinite(y) or not 0<y<=1:
        raise ValueError('old diagnostic domain14.1--300K,0<Y<=1 required')
    if th._dry_argon_mass_fraction!=0 or not th.consistent_phase_ambient or th._phase_enthalpy_tables is None:
        raise ValueError('requires the frozen explicit N2/O2/water phase caloric configuration')
    tab=_air_phase_property_table()
    names=('nitrogen','oxygen','water')
    grid=tab['temperature']
    mw=np.array([.0280134,.0319988,th.water_molecular_weight])
    dry=(1-y)/(1+th.ambient_absolute_humidity)
    mass=dry*np.array([th._dry_nitrogen_mass_fraction,th._dry_oxygen_mass_fraction,th.ambient_absolute_humidity])
    totals=mass/mw
    inert=max(y,1e-14)/th.fuel_molecular_weight
    k=np.array([np.interp(t,grid,tab[s+'_saturation']) for s in names])/th.ambient_pressure
    gas,cond=independent_split(inert,totals,k)
    cm=cond*mw
    latent=np.array([np.interp(t,grid,tab[s+'_latent']) for s in names])
    rho_c=np.array([np.interp(t,grid,tab[s+'_density']) for s in names])
    specific_volume=(inert+float(sum(gas)))*R_UNIVERSAL*t/th.ambient_pressure+float(sum(cm/rho_c))
    h=float(th._mixture_enthalpy(t,y))-float(sum(cm*latent))
    return dict(density_kg_m3=1/specific_volume,specific_enthalpy_J_kg=h,
                condensed_kg_kg=cm.tolist())


def point_scale(trajectory,old):
    model=trajectory.model
    th=model.thermodynamics
    state=trajectory.state_at(old['x'])
    _,sz=model.section_widths(state)
    peak_z,q=reflected_peak(float(state[6]),sz)
    rho=model.rhoa+(state[0]-model.rhoa)*q
    y=state[0]*state[1]*q/rho
    t,hv=th._condensed_air_state_exact(np.array([rho]),np.array([y]))
    t,h=float(t[0]),float(hv[0]/rho)
    checked=legacy_forward(th,t,y)
    density_error=abs(checked['density_kg_m3']/rho-1)
    h_error=abs(checked['specific_enthalpy_J_kg']-h)
    if density_error>1e-7 or h_error>1e-5: raise ValueError('forward exact old EOS replay failed')
    maximum_unreflected_mean_speed=model._wind(state)*math.cos(state[3])+state[4]
    maximum_unreflected_mean_ke=.5*maximum_unreflected_mean_speed**2
    gaps={}
    for basis in ('minimum','p05','median'):
        observed=old[f'observed_{basis}_K']
        cold=legacy_forward(th,observed,y)
        difference=h-cold['specific_enthalpy_J_kg']
        gaps[basis]=dict(observed_K=observed,temperature_gap_K=t-observed,
            conditional_specific_enthalpy_gap_J_kg=difference,
            positive_enthalpy_gap_over_local_mean_ke=max(difference,0.)/maximum_unreflected_mean_ke,
            conditional_observed_state_density_kg_m3=cold['density_kg_m3'],
            conditional_observed_condensed_kg_kg=cold['condensed_kg_kg'])
    return dict(channel=old['channel'],x=old['x'],exact_peak_temperature_K=t,peak_shape_z_m=peak_z,
        peak_Y_H2=y,maximum_unreflected_mean_speed_m_s=maximum_unreflected_mean_speed,
        maximum_unreflected_mean_ke_J_kg=maximum_unreflected_mean_ke,
        forward_replay_density_relative_error=density_error,forward_replay_h_error_J_kg=h_error,
        observations=gaps)


def audit():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    verify()
    paths=[REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
        REF/'stage2_coriolis_2026-09-06/field_density_10-23_step0.02.json',
        REF/'stage1_single_image_2026-09-06.json',REF/'e35_reduced.json',
        ROOT/'docs/prereg-mechanical-thermal-scale.md',Path(__file__).resolve(),
        ROOT/'tools/audit_downstream_cold_envelope.py',ROOT/'tools/audit_mixed_air_liquid.py',
        ROOT/'tools/audit_stage2_flow_effect.py',ROOT/'tools/run_preslhy_ambient_profile_audit.py']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    fields=dict(pressure_loss_density=read(paths[0]),coriolis_density=replay_document(read(paths[1])))
    frozen=read(paths[2]); trials={t['trial']:t for t in read(paths[3])['trials']}
    budgets,point_rows,checks=[],[],{}
    for variant,field in fields.items():
        for n in (10,23):
            trajectory,checks[f'{variant}_{n}']=replay(field,trials[n])
            model=trajectory.model
            data=field['interfaces'][str(n)]['downstream']
            states=np.array(data['states']); fluxes=np.array(data['fluxes'])
            cumulative=np.array(data['cumulative_sources'])
            ke=np.array([.5*velocity_moment(model,state,3) for state in states])
            thermal=fluxes[:,4]-ke
            budget_error=thermal-thermal[0]-(cumulative[:,4]-(ke-ke[0]))
            nearest={0:{'kind':'initial'},len(states)-1:{'kind':'last'}}
            local=[r for r in frozen['temperature_rows'] if r['trial']==n]
            for x in sorted({r['x'] for r in local}):
                index=int(np.argmin(abs(states[:,5]-x)))
                nearest[index]=dict(kind='nearest_saved_observation_section',requested_x_m=x)
            sections=[]
            for index,selection in sorted(nearest.items()):
                state=states[index]
                vel,rho,y,t,hv=model.profiles(state)
                _,weights=model._quadrature(model.quadrature_points)
                weights=weights*state[2]
                actual_ke=float(np.sum(.5*rho*vel**3*weights))
                actual_h=float(np.sum((hv-model.thermodynamics._ambient_enthalpy*rho)*vel*weights))
                er=max(abs(actual_ke-ke[index])/max(ke[index],1),
                       abs(actual_h-thermal[index])/max(abs(thermal[index]),1))
                mass=velocity_moment(model,state,1)
                momentum=velocity_moment(model,state,2)
                er=max(er,abs(mass/fluxes[index,0]-1),abs(momentum/np.hypot(*fluxes[index,2:4])-1))
                if er>1e-9: raise ValueError(f'actual stored energy split replay failed:{er}')
                th=model.thermodynamics
                cp=(th._mixture_enthalpy(t+.005,y)-th._mixture_enthalpy(t-.005,y))/.01
                cp2=(th._mixture_enthalpy(t+.0025,y)-th._mixture_enthalpy(t-.0025,y))/.005
                cdot=float(np.sum(rho*vel*cp*weights)); cdot2=float(np.sum(rho*vel*cp2*weights))
                if np.any(cp<=0): raise ValueError('nonpositive component sensibleCp')
                sections.append(dict(index=index,x_m=float(state[5]),**selection,
                    distance_from_requested_m=None if 'requested_x_m' not in selection else float(state[5]-selection['requested_x_m']),
                    mass_flux_kg_s=mass,hydrogen_flux_kg_s=float(fluxes[index,1]),
                    thermal_flux_W=float(thermal[index]),mean_kinetic_flux_W=float(ke[index]),
                    cumulative_accounted_ambient_mechanical_input_W=float(cumulative[index,4]),
                    thermal_change_from_handoff_W=float(thermal[index]-thermal[0]),
                    sensible_heat_capacity_flux_W_K=cdot,
                    sensible_cp_step_relative_change=abs(cdot2/cdot-1),
                    local_ke_uniform_sensible_temperature_scale_K=float(ke[index]/cdot),
                    handoff_plus_accounted_input_uniform_sensible_temperature_scale_K=float((ke[0]+cumulative[index,4])/cdot),
                    independent_integral_relative_error=er))
            budgets.append(dict(variant=variant,trial=n,stored_node_count=len(states),sections=sections,
                maximum_energy_identity_absolute_residual_W=float(max(abs(budget_error))),
                maximum_energy_identity_scaled_residual=float(max(abs(budget_error))/max(abs(thermal[0]),1)),
                initial_kinetic_flux_W=float(ke[0]),initial_thermal_flux_W=float(thermal[0]),
                largest_accounted_mechanical_budget_fraction_of_initial_thermal=float(max(ke[0]+cumulative[:,4])/abs(thermal[0])),
                no_prehandoff_tke_or_pressure_work_upper_bound=True))
            for old in local:
                if abs(old['y'])>1e-9 or abs(old['z']-trials[n]['release_height_m'])>1e-9: continue
                point_rows.append(dict(variant=variant,trial=n,**point_scale(trajectory,old)))
            print('Mechanical/thermal scale completed',variant,n,flush=True)
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('input changed')
    verify()
    result=dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        field_replay_checks=checks,budgets=budgets,point_scales=point_rows,
        new_dispersion_field=False,candidate_promoted=False,
        interpretation='Energy-scale diagnostic, not a closure or universal turbulence bound. Point gap at fixed modelY is not integrated demand and not measured composition.')
    write_new(args.output,result)
    print('Saved',args.output,flush=True)


if __name__=='__main__': audit()
