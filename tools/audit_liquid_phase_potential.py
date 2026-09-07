"""Common liquid G property/local audit, NOT an integrated plume prediction."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import math

import CoolProp as CP
import numpy as np
from scipy.optimize import brentq

from degali.addons import liquid_phase_potential as lp
from degali.addons.axisymmetric_jet import _air_phase_property_table
from degali.addons.mixed_air_liquid import ideal_flash, MW_N, MW_O
from stage2_matched_source import ROOT, REF, read, digest, write_new
from run_preslhy_ambient_profile_audit import thermodynamics
from stage1 import verify


def pressure_expansion_check():
    rows = []
    for species,tmin in lp.TRIPLE_T.items():
        eos = CP.AbstractState('HEOS',species)
        for t in np.linspace(tmin+1,lp.MAX_T,61):
            t = float(t)
            sat = lp.saturated_reference(species,t)
            for p in (50000.,101325.,200000.):
                if p <= sat.pressure_Pa*1.0001:
                    continue
                # Only compressed-liquid PT states; never force a vapor state
                # onto a liquid root. Stay1K above the triple vicinity.
                eos.update(CP.PT_INPUTS,p,t)
                if eos.phase() != CP.iphase_liquid:
                    raise ValueError('unexpected non-liquid comparison state')
                c = lp.liquid(species,t,p)
                rows.append(dict(species=species,T_K=t,P_Pa=p,
                    pressure_offset_Pa=p-sat.pressure_Pa,
                    enthalpy_difference_J_kg=(c.enthalpy_J_mol-eos.hmolar())/eos.molar_mass(),
                    gibbs_difference_J_mol=c.gibbs_J_mol-eos.gibbsmolar(),
                    volume_relative_difference=c.volume_m3_mol*eos.rhomolar()-1,
                    Cp_relative_difference=c.heat_capacity_J_mol_K/eos.cpmolar()-1))
    keys = ('enthalpy_difference_J_kg','gibbs_difference_J_mol','volume_relative_difference','Cp_relative_difference')
    return dict(rows=rows,maximum_absolute={key:max(abs(r[key]) for r in rows) for key in keys},
        limitation='Compressed pure liquids only; does not bound metastable-reference error or mixed-solid stability.')


def binary_benchmark(data):
    rows = []
    for row in data['rows']:
        t,x = row['T_K'],row['x_N2']
        def terms(p):
            return [z*lp.equilibrium_partial_pressure(sp,t,p)
                    for sp,z in zip(('Nitrogen','Oxygen'),(x,1-x))]
        p = brentq(lambda p:sum(terms(p))-p,1.,lp.MAX_P,xtol=1e-8)
        rows.append(dict(**row,candidate_bubble_P_atm=p/101325,
            candidate_y_N2=terms(p)[0]/p,relative_pressure_error=p/(101325*row['P_atm'])-1))
    return dict(rows=rows,n=len(rows),mean_absolute_relative_pressure_error=float(np.mean([abs(r['relative_pressure_error']) for r in rows])),
        maximum_absolute_relative_pressure_error=max(abs(r['relative_pressure_error']) for r in rows),
        limitation='Ideal gas/ideal liquid property comparison, not independent dispersion validation. No NBS activity/virial correlation added.')


def local_state(thermo,t,y):
    """Isolate new N2/O2 potential inside the existing gas/water ledger.

    Existing gas enthalpy interpolation, water phase and molecular weights
    are retained explicitly. This local comparison is not a complete common
    H2/N2/O2/H2O Gibbs closure or a new conservative field.
    """
    if not lp.TRIPLE_T['Nitrogen'] <= t <= lp.MAX_T:
        raise ValueError('local liquid audit outside common temperature domain')
    p = thermo.ambient_pressure
    table = _air_phase_property_table()
    old_latent = np.array([np.interp(t,table['temperature'],table[f'{sp}_latent'])
                           for sp in ('nitrogen','oxygen','water')])
    water_k = float(np.interp(t,table['temperature'],table['water_saturation']))/p
    water_rho = float(np.interp(t,table['temperature'],table['water_density']))
    mw = np.array([MW_N,MW_O,thermo.water_molecular_weight])
    dry = (1-y)/(1+thermo.ambient_absolute_humidity)
    mass = dry*np.array([thermo._dry_nitrogen_mass_fraction,thermo._dry_oxygen_mass_fraction,
                        thermo.ambient_absolute_humidity])
    totals = mass/mw
    inert = y/thermo.fuel_molecular_weight
    k = [lp.equilibrium_partial_pressure(sp,t,p)/p for sp in ('Nitrogen','Oxygen')]
    f = ideal_flash(inert+totals[2],*totals[:2],*k)
    if totals[2]/f.gas_mol <= water_k*(1+1e-12):
        gas = np.array([f.gas_nitrogen_mol,f.gas_oxygen_mol,totals[2]])
    else:
        f = ideal_flash(inert,*totals[:2],*k,reservoir_vapor_fraction=water_k)
        if f.reservoir_vapor_mol > totals[2]*(1+1e-12):
            raise ValueError('water reservoir inventory exceeded')
        gas = np.array([f.gas_nitrogen_mol,f.gas_oxygen_mol,f.reservoir_vapor_mol])
    condensed = totals-gas
    condensed_mass = condensed*mw
    # Consistent molar h decrements use the SAME audit mole convention as
    # the amount solver. Differences from HEOS molecular weights are recorded.
    delta_h = [lp.ideal_gas_to_liquid_enthalpy(sp,t,p) for sp in ('Nitrogen','Oxygen')]
    v_l = sum(n*lp.liquid(sp,t,p).volume_m3_mol
              for sp,n in zip(('Nitrogen','Oxygen'),condensed[:2]))+condensed_mass[2]/water_rho
    v_g = (inert+sum(gas))*lp.R*t/p
    h = thermo._mixture_enthalpy(t,y)-float(np.dot(condensed[:2],delta_h))-condensed_mass[2]*old_latent[2]
    # Exact HEOS ideal-air h increment vs frozen interpolation; separate it
    # from the new caloric correction rather than hide it in a field score.
    grid,tables = thermo._phase_enthalpy_tables
    interpolation_error = 0.
    for sp,m in zip(('Hydrogen','Nitrogen','Oxygen'),(y,mass[0],mass[1])):
        h0,_,_,mw_eos = lp.ideal_standard(sp,t)
        ha = lp.ideal_standard(sp,thermo.ambient_temperature)[0]
        old = np.interp(t,grid,tables[sp])-np.interp(thermo.ambient_temperature,grid,tables[sp])
        interpolation_error += m*((h0-ha)/mw_eos-old)
    return dict(T_K=t,specific_enthalpy_J_kg=float(h),density_kg_m3=float(1/(v_l+v_g)),
        liquid_N2_kg_kg=float(condensed_mass[0]),liquid_O2_kg_kg=float(condensed_mass[1]),
        ice_water_kg_kg=float(condensed_mass[2]),gas_H2_mole_fraction=float(inert/(inert+sum(gas))),
        caloric_change_at_candidate_amounts_J_kg=float(np.dot(condensed_mass[:2],old_latent[:2])-np.dot(condensed[:2],delta_h)),
        frozen_gas_interpolation_difference_J_kg=float(interpolation_error),flash=asdict(f))


def local_audit(previous,trials):
    thermos = {}
    rows = []
    for raw in previous['rows']:
        t = raw['original_exact_T_K']
        r = {key:raw[key] for key in ('kind','trial','x_m','original_exact_T_K','fraction')}
        for key in ('channel','Gaussian_shape'):
            if key in raw: r[key] = raw[key]
        if not lp.TRIPLE_T['Nitrogen'] <= t <= lp.MAX_T:
            r['status'] = 'below_liquid_domain' if t < lp.TRIPLE_T['Nitrogen'] else 'above_liquid_audit_domain'
            rows.append(r)
            continue
        n,y = raw['trial'],raw['fraction']
        if n not in thermos:
            thermos[n] = thermodynamics(trials[n],consistent=True)
        thermo = thermos[n]
        old,screen = raw['old_pure_phase'],raw['new_ideal_solution']
        new = local_state(thermo,t,y)
        target_h = old['specific_enthalpy_J_kg']
        def residual(temp): return local_state(thermo,temp,y)['specific_enthalpy_J_kg']-target_h
        th = None
        if residual(lp.TRIPLE_T['Nitrogen'])*residual(lp.MAX_T) <= 0:
            tn = brentq(residual,lp.TRIPLE_T['Nitrogen'],lp.MAX_T,xtol=1e-9)
            th = dict(**local_state(thermo,tn,y),temperature_shift_K=tn-t,enthalpy_residual_J_kg=residual(tn))
        r.update(status='evaluated',candidate=new,fixed_old_enthalpy=th,
            delta_h_vs_old_pure_J_kg=new['specific_enthalpy_J_kg']-target_h,
            delta_h_vs_previous_ideal_screen_J_kg=new['specific_enthalpy_J_kg']-screen['specific_enthalpy_J_kg'],
            relative_density_change_vs_old_pure=new['density_kg_m3']/old['density_kg_m3']-1)
        rows.append(r)
    evaluated = [r for r in rows if r['status']=='evaluated']
    return dict(rows=rows,summary=dict(total=len(rows),evaluated=len(evaluated),
        below_domain=sum(r['status']=='below_liquid_domain' for r in rows),
        above_domain=sum(r['status']=='above_liquid_audit_domain' for r in rows),
        maximum_absolute_h_change_vs_old_pure_J_kg=max(abs(r['delta_h_vs_old_pure_J_kg']) for r in evaluated),
        maximum_absolute_h_change_vs_previous_mixed_screen_J_kg=max(abs(r['delta_h_vs_previous_ideal_screen_J_kg']) for r in evaluated),
        maximum_absolute_caloric_only_change_J_kg=max(abs(r['candidate']['caloric_change_at_candidate_amounts_J_kg']) for r in evaluated),
        maximum_absolute_gas_interpolation_difference_J_kg=max(abs(r['candidate']['frozen_gas_interpolation_difference_J_kg']) for r in evaluated),
        maximum_absolute_density_relative_change=max(abs(r['relative_density_change_vs_old_pure']) for r in evaluated),
        maximum_fixed_h_T_shift_K=max(abs(r['fixed_old_enthalpy']['temperature_shift_K']) for r in evaluated if r['fixed_old_enthalpy'] is not None)),
        limitation='Reused local states, not reintegrated fields. Legacy gas interpolation/water ledger and molecular weights retained to isolate air-liquid potential. Subtriple/mixed-solid domain remains unsupported.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    verify()
    old_path = REF/'mixed_air_liquid_property_screen_2026-09-06.json'
    previous = read(old_path)
    for name,expected in previous['sha256'].items():
        if digest(ROOT/name) != expected: raise ValueError(f'completed screen input hash mismatch:{name}')
    paths = [old_path,ROOT/'reference/air_mixture/NBS3921_selected_binary_vle.json',
        REF/'e35_reduced.json',ROOT/'src/degali/addons/liquid_phase_potential.py',
        ROOT/'tests/test_liquid_phase_potential.py',ROOT/'docs/prereg-liquid-phase-potential.md',
        Path(__file__).resolve(),ROOT/'tools/run_preslhy_ambient_profile_audit.py']
    hashes = {str(p.relative_to(ROOT)):digest(p) for p in paths}
    pressure = pressure_expansion_check()
    binary = binary_benchmark(read(paths[1]))
    local = local_audit(previous,{t['trial']:t for t in read(paths[2])['trials']})
    metadata = {sp:dict(backend_R=lp.saturated_reference(sp,70).backend_gas_constant,
        backend_molecular_weight=lp.saturated_reference(sp,70).molecular_weight_kg_mol,
        audit_molecular_weight=mw) for sp,mw in zip(('Nitrogen','Oxygen'),(MW_N,MW_O))}
    assert hashes == {str(p.relative_to(ROOT)):digest(p) for p in paths}
    verify()
    write_new(args.output,dict(schema='stage2-common-liquid-potential-v1',completed=True,
        created_utc=datetime.now(timezone.utc).isoformat(),input_sha256=hashes,
        reused_screen_inputs_verified=len(previous['sha256']),metadata=metadata,
        pressure_expansion=pressure,binary_benchmark=binary,local_effect=local,
        stage1_unchanged=True,new_integrated_dispersion_field=False,candidate_promoted=False))
    print('Pressure approximation:',pressure['maximum_absolute'])
    print('Binary mean/max:',binary['mean_absolute_relative_pressure_error'],binary['maximum_absolute_relative_pressure_error'])
    print('Local:',local['summary'])
    print('Saved:',args.output)


if __name__ == '__main__': main()
