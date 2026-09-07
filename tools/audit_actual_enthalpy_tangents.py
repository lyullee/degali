"""Ratio-free thermal/species flux response on actual frozen radial profiles."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import math
from pathlib import Path

import numpy as np

from stage2_matched_source import ROOT,REF,TRIALS,read,digest,write_new
from run_preslhy_ambient_profile_audit import replay
from audit_mechanical_thermal_scale import legacy_forward
from degali.addons.species_carried_enthalpy import binary_enthalpy_flux
from stage1 import verify


def caloric_tangent(th,t,y,scale=1.):
    dt=.002*scale
    dy=min(1e-5,.005*min(y,1-y))*scale
    if dy<=0: raise ValueError('tangent requires an interior mass fraction')
    ht=(legacy_forward(th,t+dt,y)['specific_enthalpy_J_kg']-
        legacy_forward(th,t-dt,y)['specific_enthalpy_J_kg'])/(2*dt)
    hy=(legacy_forward(th,t,y+dy)['specific_enthalpy_J_kg']-
        legacy_forward(th,t,y-dy)['specific_enthalpy_J_kg'])/(2*dy)
    return ht,hy


def profile_q(model,state,q):
    q=np.atleast_1d(np.asarray(q,dtype=float))
    shape=np.exp(-q)
    rho=model.rhoa+(state[0]-model.rhoa)*shape
    c=state[0]*state[1]*shape
    y=c/rho
    t,hv=model.thermodynamics._condensed_air_state_exact(rho,y)
    return rho,c,y,t,hv/rho


def sample(model,state,q):
    th=model.thermodynamics
    rho,c,y,t,h=(float(v[0]) for v in profile_q(model,state,q))
    forward=legacy_forward(th,t,y)
    rho_error=abs(forward['density_kg_m3']/rho-1)
    h_error=abs(forward['specific_enthalpy_J_kg']-h)
    if rho_error>1e-7 or h_error>1e-5: raise ValueError('old exact forward replay failed')
    htc,hyc=caloric_tangent(th,t,y,1.)
    ht,hy=caloric_tangent(th,t,y,.5)
    derivatives=[]
    for dq in (1e-4,5e-5):
        rp,cp,yp,tp,hp=profile_q(model,state,np.array([q-dq,q+dq]))
        derivatives.append(np.array([(a[1]-a[0])/(2*dq) for a in (yp,tp,hp)]))
    coarse,fine=derivatives
    yq,tq,hq=fine
    analytic_yq=-c*model.rhoa/rho**2
    chain=ht*tq+hy*analytic_yq
    norm=max(abs(hq),abs(ht*tq),abs(hy*analytic_yq),1.)
    chain_error=abs(hq-chain)/norm
    refinement=max(abs(ht-htc)/max(abs(ht),1),abs(hy-hyc)/max(abs(hy),1),
        float(max(abs(fine-coarse)/np.maximum(abs(fine),1))))
    y_derivative_error=abs(yq-analytic_yq)/max(abs(analytic_yq),1e-8)
    _,sz=model.section_widths(state)
    qz=math.sqrt(2*q)/sz
    gh,gy=hq*qz,analytic_yq*qz
    thermal_response=-rho*(gh-hy*gy)
    naive_response=-rho*gh
    # Evaluate polynomial coefficients at0and1, not proposed physical diffusivities.
    flux0=binary_enthalpy_flux(density=rho,species_diffusivity=1.,thermal_diffusivity=0.,
        enthalpy_gradient=gh,mass_fraction_gradient=gy,enthalpy_composition_derivative=hy)
    flux1=binary_enthalpy_flux(density=rho,species_diffusivity=1.,thermal_diffusivity=1.,
        enthalpy_gradient=gh,mass_fraction_gradient=gy,enthalpy_composition_derivative=hy)
    response_scale=max(abs(thermal_response),abs(naive_response),1.)
    coefficient_error=abs(float(flux1.total_enthalpy_flux-flux0.total_enthalpy_flux)-thermal_response)/response_scale
    same_d_error=abs(float(flux1.total_enthalpy_flux)-naive_response)/response_scale
    offset=1e6
    gauge_response=-rho*((gh+offset*gy)-(hy+offset)*gy)
    gauge_error=abs(gauge_response-thermal_response)/response_scale
    if max(coefficient_error,same_d_error,gauge_error)>1e-10:
        raise ValueError('new flux coefficient or reference covariance failed')
    condensed=np.array(forward['condensed_kg_kg'])
    phase='dry_air_condensed' if sum(condensed[:2])>1e-9 else ('water_condensed' if condensed[2]>1e-9 else 'all_gas')
    gas_hy_error=None
    if phase=='all_gas':
        expected=th._mixture_enthalpy(t,1.)-th._mixture_enthalpy(t,0.)
        gas_hy_error=abs(hy-expected)/max(abs(expected),1)
        if gas_hy_error>1e-6: raise ValueError('single-gas composition tangent differs from component enthalpies')
    resolved=ht>0 and max(chain_error,refinement,y_derivative_error)<=1e-3
    return dict(q=q,T_K=t,Y_H2=y,density_kg_m3=rho,specific_h_J_kg=h,
        phase_classification=phase,condensed_kg_kg=condensed,forward_rho_relative_error=rho_error,
        forward_h_error_J_kg=h_error,h_T_J_kg_K=ht,h_Y_J_kg=hy,Y_q=analytic_yq,T_q=tq,h_q=hq,
        thermal_h_q=ht*tq,composition_h_q=hy*analytic_yq,chain_scaled_error=chain_error,
        derivative_refinement=refinement,analytic_Y_q_relative_error=y_derivative_error,
        status='resolved' if resolved else 'nondifferentiable_or_unresolved',
        thermal_gradient_fraction=None if abs(hq)<1 else ht*tq/hq,
        composition_gradient_fraction=None if abs(hq)<1 else hy*analytic_yq/hq,
        derivative_JH_over_DY_wrt_DT_over_DY_J_m4=thermal_response,
        naive_derivative_JH_over_DY_wrt_DT_over_DY_J_m4=naive_response,
        coefficient_scaled_error=coefficient_error,same_D_scaled_error=same_d_error,
        species_reference_response_scaled_error=gauge_error,all_gas_composition_tangent_relative_error=gas_hy_error)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    verify()
    paths=[REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
        REF/'stage1_single_image_2026-09-06.json',REF/'e35_reduced.json',
        ROOT/'docs/prereg-actual-enthalpy-tangents.md',Path(__file__).resolve(),
        ROOT/'src/degali/addons/species_carried_enthalpy.py',
        ROOT/'tools/audit_mechanical_thermal_scale.py',ROOT/'tools/run_preslhy_ambient_profile_audit.py']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    field,frozen=read(paths[0]),read(paths[1]); trials={t['trial']:t for t in read(paths[2])['trials']}
    rows,checks=[],{}
    for n in TRIALS:
        trajectory,checks[str(n)]=replay(field,trials[n])
        states=np.array(field['interfaces'][str(n)]['downstream']['states'])
        selections=[dict(kind='initial',index=0)]
        for row in frozen['concentration_pairs']:
            if row['trial']==n:
                selections.append(dict(kind='nearest_stored_concentration_section',requested_x_m=row['x'],
                    index=int(np.argmin(abs(states[:,5]-row['x'])))))
        for chosen in selections:
            state=states[chosen['index']]
            for q in (0.,.125,.25,.5,1.,2.,4.):
                row=dict(trial=n,x_m=float(state[5]),**chosen)
                try: row.update(sample(trajectory.model,state,q))
                except (ValueError,RuntimeError,FloatingPointError) as err:
                    row.update(q=q,status='evaluation_failed',error=str(err))
                rows.append(row)
        print('Actual enthalpy tangents trial',n,'done',flush=True)
    resolved=[r for r in rows if r['status']=='resolved']
    groups={}
    for phase in ('all_gas','water_condensed','dry_air_condensed'):
        group=[r for r in resolved if r['phase_classification']==phase]
        fractions=[r['composition_gradient_fraction'] for r in group if r['composition_gradient_fraction'] is not None]
        groups[phase]=dict(n=len(group),composition_gradient_fraction_range=None if not fractions else [min(fractions),max(fractions)],
            composition_gradient_fraction_median=None if not fractions else float(np.median(fractions)))
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('input changed')
    verify()
    result=dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,rows=rows,replay_checks=checks,
        summary=dict(total_samples=len(rows),resolved=len(resolved),
            unresolved=sum(r['status']=='nondifferentiable_or_unresolved' for r in rows),
            evaluation_failed=sum(r['status']=='evaluation_failed' for r in rows),by_phase=groups,
            maximum_resolved_chain_error=max((r['chain_scaled_error'] for r in resolved),default=None)),
        physical_diffusivity_selected=False,new_field=False,candidate_promoted=False,
        limitation='Local fixed-equilibrium single-velocity tangents and unit response, not a multiphase transport validation or dispersion improvement.')
    write_new(args.output,result)
    print(result['summary'],flush=True)
    print('Saved',args.output,flush=True)


if __name__=='__main__': main()
