"""Frozen-key local screening of common nonideal liquid G/h; no new field."""
from dataclasses import asdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from degali.addons import liquid_excess_potential as ep
from degali.addons import liquid_phase_potential as lp
from degali.addons.axisymmetric_jet import _air_phase_property_table
from degali.addons.mixed_air_liquid import MW_N,MW_O
from audit_liquid_phase_potential import local_state as ideal_local
from run_preslhy_ambient_profile_audit import thermodynamics
from stage2_matched_source import ROOT,REF,read,digest,write_new
from stage1 import verify


def binary(data):
    rows=[]
    for row in data['rows']:
        t,x=row['T_K'],row['x_N2']
        variants={}
        for method in ('ideal','quadratic','pchip'):
            gamma=(1.,1.) if method=='ideal' else ep.activity_coefficients(t,x,interpolation=method)
            def terms(p):
                return [z*g*lp.equilibrium_partial_pressure(sp,t,p)
                    for sp,z,g in zip(('Nitrogen','Oxygen'),(x,1-x),gamma)]
            p=brentq(lambda p:sum(terms(p))-p,1.,lp.MAX_P,xtol=1e-8)
            variants[method]=dict(P_atm=p/101325,y_N2=terms(p)[0]/p,
                relative_pressure_error=p/(101325*row['P_atm'])-1,
                vapor_fraction_error=terms(p)[0]/p-row['y_N2'])
        rows.append(dict(**row,variants=variants))
    summary={m:dict(mean_absolute_pressure_error=float(np.mean([abs(r['variants'][m]['relative_pressure_error']) for r in rows])),
        maximum_absolute_pressure_error=max(abs(r['variants'][m]['relative_pressure_error']) for r in rows),
        mean_absolute_vapor_fraction_error=float(np.mean([abs(r['variants'][m]['vapor_fraction_error']) for r in rows]))) for m in ('ideal','quadratic','pchip')}
    return dict(rows=rows,summary=summary,same_source_correlation_not_independent_validation=True,
        gas_model='same ideal-gas common reference in all variants; no vapor virial correction')


def local_state(thermo,t,y,method):
    """Existing gas/water thermal ledger with a coupled nonideal dry-air liquid.

    This remains a local hybrid diagnostic; not a full common H2/air/H2O EOS.
    The only added heat is liquid_moles*hE, with SAME-G activity partitioning.
    """
    p=thermo.ambient_pressure
    table=_air_phase_property_table()
    water_latent=float(np.interp(t,table['temperature'],table['water_latent']))
    water_k=float(np.interp(t,table['temperature'],table['water_saturation']))/p
    water_rho=float(np.interp(t,table['temperature'],table['water_density']))
    mw=np.array([MW_N,MW_O,thermo.water_molecular_weight])
    dry=(1-y)/(1+thermo.ambient_absolute_humidity)
    totals=dry*np.array([thermo._dry_nitrogen_mass_fraction,thermo._dry_oxygen_mass_fraction,thermo.ambient_absolute_humidity])/mw
    inert=y/thermo.fuel_molecular_weight
    flash=ep.flash(t,p,inert+totals[2],*totals[:2],interpolation=method)
    f=flash.amounts
    if totals[2]/f.gas_mol<=water_k*(1+1e-12):
        water_gas=totals[2]
    else:
        flash=ep.flash(t,p,inert,*totals[:2],interpolation=method,reservoir_vapor_fraction=water_k)
        f=flash.amounts
        water_gas=f.reservoir_vapor_mol
        if water_gas>totals[2]*(1+1e-12): raise ValueError('water inventory exceeded')
    ice_mass=(totals[2]-water_gas)*mw[2]
    liquid_moles=np.array([f.liquid_nitrogen_mol,f.liquid_oxygen_mol])
    liquid_pure=[lp.liquid(sp,t,p) for sp in ('Nitrogen','Oxygen')]
    decrement=[lp.ideal_gas_to_liquid_enthalpy(sp,t,p) for sp in ('Nitrogen','Oxygen')]
    excess_heat=0.
    liquid_x=None
    if sum(liquid_moles)>0:
        liquid_x=liquid_moles[0]/sum(liquid_moles)
        excess_heat=sum(liquid_moles)*ep.excess(t,liquid_x,interpolation=method).enthalpy_J_mol
    h=thermo._mixture_enthalpy(t,y)-np.dot(liquid_moles,decrement)-ice_mass*water_latent+excess_heat
    liquid_volume=sum(n*v.volume_m3_mol for n,v in zip(liquid_moles,liquid_pure))+ice_mass/water_rho
    gas_moles=inert+f.gas_nitrogen_mol+f.gas_oxygen_mol+water_gas
    rho=1/(liquid_volume+gas_moles*lp.R*t/p)
    return dict(T_K=t,specific_enthalpy_J_kg=float(h),density_kg_m3=float(rho),
        liquid_N2_kg_kg=float(liquid_moles[0]*mw[0]),liquid_O2_kg_kg=float(liquid_moles[1]*mw[1]),
        liquid_x_N2=None if liquid_x is None else float(liquid_x),
        ice_water_kg_kg=float(ice_mass),excess_heat_at_candidate_amounts_J_kg=float(excess_heat),
        flash=asdict(flash))


def local_screen(previous,trials):
    rows=[];thermos={}
    for raw in previous['rows']:
        t,y,n=raw['original_exact_T_K'],raw['fraction'],raw['trial']
        row={k:raw[k] for k in ('kind','trial','x_m','original_exact_T_K','fraction','channel','Gaussian_shape') if k in raw}
        if not ep.MIN_T<=t<=ep.MAX_T:
            row['status']='below_domain' if t<ep.MIN_T else 'above_domain'
            rows.append(row);continue
        if n not in thermos: thermos[n]=thermodynamics(trials[n],consistent=True)
        th=thermos[n]
        old=ideal_local(th,t,y)
        row.update(status='evaluated',common_ideal=old,variants={})
        for method in ('quadratic','pchip'):
            new=local_state(th,t,y,method)
            change=new['specific_enthalpy_J_kg']-old['specific_enthalpy_J_kg']
            def residual(temp): return local_state(th,temp,y,method)['specific_enthalpy_J_kg']-old['specific_enthalpy_J_kg']
            endpoints=[residual(z) for z in (ep.MIN_T,ep.MAX_T)]
            fixed_h=None
            if endpoints[0]*endpoints[1]<=0:
                temp=brentq(residual,ep.MIN_T,ep.MAX_T,xtol=1e-9)
                fixed_h=dict(**local_state(th,temp,y,method),temperature_shift_K=temp-t,
                             enthalpy_residual_J_kg=residual(temp))
            # Purely caloric size at unchanged ideal amounts, distinguished
            # from the coupled equilibrium response above.
            f=old['flash'];ln=f['liquid_nitrogen_mol'];lo=f['liquid_oxygen_mol']
            direct=0. if ln+lo==0 else (ln+lo)*ep.excess(t,ln/(ln+lo),interpolation=method).enthalpy_J_mol
            row['variants'][method]=dict(candidate=new,enthalpy_change_vs_common_ideal_J_kg=change,
                amount_redistribution_enthalpy_change_J_kg=change-new['excess_heat_at_candidate_amounts_J_kg'],
                excess_heat_at_original_ideal_amounts_J_kg=direct,
                density_relative_change=new['density_kg_m3']/old['density_kg_m3']-1,
                fixed_common_ideal_enthalpy=fixed_h,endpoint_enthalpy_residuals_J_kg=endpoints)
        rows.append(row)
    evaluated=[r for r in rows if r['status']=='evaluated']
    summaries={}
    for method in ('quadratic','pchip'):
        vals=[r['variants'][method] for r in evaluated]
        fixed=[v['fixed_common_ideal_enthalpy'] for v in vals if v['fixed_common_ideal_enthalpy'] is not None]
        summaries[method]=dict(
            maximum_abs_h_change_J_kg=max(abs(v['enthalpy_change_vs_common_ideal_J_kg']) for v in vals),
            maximum_abs_caloric_at_original_amounts_J_kg=max(abs(v['excess_heat_at_original_ideal_amounts_J_kg']) for v in vals),
            maximum_abs_density_relative_change=max(abs(v['density_relative_change']) for v in vals),
            fixed_h_solved=len(fixed),fixed_h_outside_domain=len(vals)-len(fixed),
            maximum_abs_fixed_h_T_shift_K=max(abs(v['temperature_shift_K']) for v in fixed))
    summaries['interpolation_difference']=dict(maximum_abs_h_J_kg=max(abs(r['variants']['quadratic']['candidate']['specific_enthalpy_J_kg']-r['variants']['pchip']['candidate']['specific_enthalpy_J_kg']) for r in evaluated))
    return dict(rows=rows,total=len(rows),evaluated=len(evaluated),
        below_domain=sum(r['status']=='below_domain' for r in rows),above_domain=sum(r['status']=='above_domain' for r in rows),
        summary=summaries,new_field=False,old_gas_water_ledger_retained=True)


def potential_screen():
    rows=[]
    for method in ('quadratic','pchip'):
        for t in np.linspace(65.,77.5,51):
            for x in np.linspace(0,1,101):
                e=ep.excess(float(t),float(x),interpolation=method)
                rows.append(dict(interpolation=method,T_K=float(t),x_N2=float(x),**asdict(e)))
    return dict(sample_count=len(rows),summary={m:{name:dict(minimum=min(r[name] for r in rows if r['interpolation']==m),maximum=max(r[name] for r in rows if r['interpolation']==m)) for name in ('gibbs_J_mol','enthalpy_J_mol','heat_capacity_J_mol_K')} for m in ('quadratic','pchip')},
        note='Sampled grid only; hE/CpE derived from interpolation, not measured caloric validation. No claimed global or sub65K bound.')


def main():
    output=REF/'liquid_excess_potential_2026-09-06.json'
    if output.exists(): raise FileExistsError(output)
    paths=[REF/'mixed_air_liquid_property_screen_2026-09-06.json',REF/'liquid_phase_potential_2026-09-06.json',
        ROOT/'reference/air_mixture/NBS3921_selected_binary_vle.json',ROOT/'reference/air_mixture/NBS_Report_3921_1955.pdf',
        REF/'e35_reduced.json',ROOT/'src/degali/addons/liquid_excess_potential.py',
        ROOT/'tests/test_liquid_excess_potential.py',ROOT/'docs/prereg-liquid-excess-potential.md',
        ROOT/'tools/audit_liquid_phase_potential.py',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    previous=read(paths[0]);common=read(paths[1]);data=read(paths[2])
    if digest(paths[3])!=data['sha256']: raise ValueError('source PDF changed')
    for doc,key in ((previous,'sha256'),(common,'input_sha256')):
        for path,value in doc[key].items():
            if digest(ROOT/path)!=value: raise ValueError('previous source changed:'+path)
    verify()
    benchmark=binary(data);print('Binary:',benchmark['summary'],flush=True)
    potentials=potential_screen();print('Excess:',potentials['summary'],flush=True)
    local=local_screen(previous,{r['trial']:r for r in read(paths[4])['trials']})
    print('Local:',local['summary'],'counts',local['total'],local['evaluated'],local['below_domain'],local['above_domain'],flush=True)
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('audit input changed')
    verify()
    write_new(output,dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        binary_benchmark=benchmark,potential_grid=potentials,local_screen=local,
        stage1_unchanged=True,new_integrated_dispersion_field=False,candidate_promoted=False))
    print('Saved',output,flush=True)


if __name__=='__main__': main()
