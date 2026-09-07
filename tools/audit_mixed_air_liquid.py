"""Bounded phase-property screen; deliberately NOT a new dispersion field."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import math

import numpy as np
from scipy.optimize import brentq

from stage2_matched_source import ROOT, REF, TRIALS, read, digest, write_new
from run_preslhy_ambient_profile_audit import replay
from degali.addons.axisymmetric_jet import _air_phase_property_table
from degali.addons.mixed_air_liquid import ideal_flash, liquid_properties, nbs_activity_coefficients, R, MW_N, MW_O, MIN_T, MAX_T
from stage1 import verify


def independent_split(inert, totals, k):
    """Independent pure-condensate active set; N2,O2,H2O on a mole basis."""
    for mask in range(8):
        active = np.array([bool(mask & (1 << i)) for i in range(3)])
        denominator = 1-float(sum(k[active]))
        if denominator <= 0:
            continue
        v = (inert+float(sum(totals[~active])))/denominator
        gas = np.where(active, k*v, totals)
        if np.any(gas > totals+1e-12*max(v,1)):
            continue
        if np.any(gas[~active]/v > k[~active]*(1+1e-12)):
            continue
        return gas, totals-gas
    raise ValueError('no pure-condensate split')


def property_state(thermo, temperature, fraction, *, mixed):
    """Forward T,Y closure using IDENTICAL component tables and ice handling."""
    if not MIN_T <= temperature <= MAX_T:
        raise ValueError('bounded liquid diagnostic domain exceeded')
    if thermo._dry_argon_mass_fraction != 0 or not thermo.consistent_phase_ambient:
        raise ValueError('screen only supports frozen explicit N2/O2/H2O atmosphere')
    table = _air_phase_property_table()
    names = ('nitrogen', 'oxygen', 'water')
    grid = table['temperature']
    k = np.array([np.interp(temperature, grid, table[f'{s}_saturation']) for s in names])/thermo.ambient_pressure
    latent = np.array([np.interp(temperature, grid, table[f'{s}_latent']) for s in names])
    liquid_rho = np.array([np.interp(temperature, grid, table[f'{s}_density']) for s in names])
    mw = np.array([MW_N, MW_O, thermo.water_molecular_weight])
    dry = (1-fraction)/(1+thermo.ambient_absolute_humidity)
    masses = dry*np.array([thermo._dry_nitrogen_mass_fraction, thermo._dry_oxygen_mass_fraction,
                          thermo.ambient_absolute_humidity])
    totals = masses/mw
    inert = fraction/thermo.fuel_molecular_weight
    if mixed:
        # Try all water gaseous, then an ice reservoir. Verify its inventory.
        f = ideal_flash(inert+totals[2], *totals[:2], *k[:2])
        if totals[2]/f.gas_mol <= k[2]*(1+1e-12):
            gas = np.array([f.gas_nitrogen_mol, f.gas_oxygen_mol, totals[2]])
        else:
            f = ideal_flash(inert, *totals[:2], *k[:2], reservoir_vapor_fraction=float(k[2]))
            if f.reservoir_vapor_mol > totals[2]*(1+1e-12):
                raise ValueError('ice reservoir does not have enough water')
            gas = np.array([f.gas_nitrogen_mol, f.gas_oxygen_mol, f.reservoir_vapor_mol])
        condensed = totals-gas
        flash_check = asdict(f)
    else:
        gas, condensed = independent_split(inert, totals, k)
        flash_check = None
    condensed_mass = condensed*mw
    v_gas = (inert+float(sum(gas)))*R*temperature/thermo.ambient_pressure
    v_condensed = float(sum(condensed_mass/liquid_rho))
    h = thermo._mixture_enthalpy(temperature, fraction)-float(sum(condensed_mass*latent))
    return dict(temperature_K=float(temperature), fraction_H2=float(fraction), density_kg_m3=1/(v_gas+v_condensed),
        specific_enthalpy_J_kg=float(h), gas_moles_per_kg=(inert+float(sum(gas))),
        condensed_N2_kg_kg=float(condensed_mass[0]), condensed_O2_kg_kg=float(condensed_mass[1]),
        condensed_water_kg_kg=float(condensed_mass[2]), gas_H2_mole_fraction=inert/(inert+float(sum(gas))),
        flash_check=flash_check)


def nbs_comparison(data):
    rows = []
    for row in data['rows']:
        t, x = row['T_K'], row['x_N2']
        prop = liquid_properties(t)
        pn, po = [prop[s]['pressure']/101325. for s in ('Nitrogen','Oxygen')]
        ideal_p = x*pn+(1-x)*po
        pure_p = pn+po
        gn, go = nbs_activity_coefficients(t, x)
        correlated_p = x*gn*pn+(1-x)*go*po
        # NBS eq printed p10-11 includes a binary-gas fugacity correction.
        # Reproduce this binary benchmark ONLY, never apply its B' to H2 gas.
        bn = float(np.interp(t, [65,70,77.5], [-.0646,-.0539,-.0385]))
        bo = float(np.interp(t, [65,70,77.5], [-.0839,-.0669,-.0491]))
        def nbs_terms(p):
            return x*gn*pn*math.exp(-bn*(p-pn)), (1-x)*go*po*math.exp(-bo*(p-po))
        nbs_p = brentq(lambda p: sum(nbs_terms(p))-p, .5*ideal_p, 2*correlated_p)
        rows.append(dict(**row, ideal_bubble_P_atm=ideal_p, ideal_y_N2=x*pn/ideal_p,
            independent_pure_condensates_P_atm=pure_p,
            published_activity_binary_fugacity_P_atm=nbs_p,
            published_activity_binary_fugacity_y_N2=nbs_terms(nbs_p)[0]/nbs_p,
            ideal_pressure_error_fraction=ideal_p/row['P_atm']-1,
            pure_pressure_error_fraction=pure_p/row['P_atm']-1,
            published_correlation_pressure_error_fraction=nbs_p/row['P_atm']-1))
    summary = {label: dict(n=len(rows), mean_absolute_relative_pressure_error=float(np.mean([
        abs(row[key]) for row in rows])), maximum_absolute_relative_pressure_error=max(abs(row[key]) for row in rows))
        for label, key in [('ideal_solution','ideal_pressure_error_fraction'),
                           ('independent_pure_condensates','pure_pressure_error_fraction'),
                           ('published_same_data_correlation','published_correlation_pressure_error_fraction')]}
    return dict(rows=rows, summary=summary, note='binary solution property test, not H2-plume validation; published correlation reuses this source data')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    verify()
    paths = [ROOT/'reference/air_mixture/NBS3921_selected_binary_vle.json', ROOT/'reference/air_mixture/NBS_Report_3921_1955.pdf',
        REF/'phase_ambient_consistency_field_complete_2026-09-05.json', REF/'stage1_single_image_2026-09-06.json',
        REF/'e35_reduced.json', ROOT/'docs/prereg-mixed-air-liquid-screen.md',
        ROOT/'src/degali/addons/mixed_air_liquid.py', Path(__file__).resolve()]
    hashes = {str(p.relative_to(ROOT)).replace('\\','/'): digest(p) for p in paths}
    nbs = read(paths[0])
    if digest(paths[1]) != nbs['sha256']:
        raise ValueError('NBS transcription source hash mismatch')
    benchmark = nbs_comparison(nbs)
    field, frozen = read(paths[2]), read(paths[3])
    trials = {t['trial']: t for t in read(paths[4])['trials']}
    rows, checks = [], {}
    for n in TRIALS:
        trajectory, checks[str(n)] = replay(field, trials[n])
        model, thermo = trajectory.model, trajectory.model.thermodynamics
        xs = sorted({float(trajectory.result.states[0,5]), *[r['x'] for r in frozen['concentration_pairs'] if r['trial'] == n]})
        samples = []
        for x in xs:
            state = trajectory.state_at(x)
            for shape in (1.,.9,.75,.5,.25,.1,.01):
                rho = model.rhoa+(state[0]-model.rhoa)*shape
                fuel = state[0]*state[1]*shape/rho
                samples.append(dict(kind='radial_section', trial=n, x_m=x, Gaussian_shape=shape, density=rho, fraction=fuel))
        for row in frozen['temperature_rows']:
            if row['trial'] != n:
                continue
            state = trajectory.state_at(row['x'])
            sy, sz = model.section_widths(state)
            shape = math.exp(-.5*(row['y']/sy)**2)*(math.exp(-.5*((row['z']-state[6])/sz)**2)
                +math.exp(-.5*((row['z']+state[6])/sz)**2))
            rho = model.rhoa+(state[0]-model.rhoa)*shape
            samples.append(dict(kind='thermocouple', trial=n, channel=row['channel'], x_m=row['x'],
                                Gaussian_shape=shape, density=rho, fraction=state[0]*state[1]*shape/rho))
        densities, fractions = [np.array([s[k] for s in samples]) for k in ('density','fraction')]
        exact_t, exact_rho_h = thermo._condensed_air_state_exact(densities, fractions)
        interpolated_t, interpolated_h = thermo._condensed_air_state(densities, fractions)
        for i, sample in enumerate(samples):
            t, y = float(exact_t[i]), sample['fraction']
            record = dict(**sample, original_exact_T_K=t, original_interpolated_T_K=float(interpolated_t[i]),
                lookup_T_error_K=float(interpolated_t[i]-t))
            if not MIN_T <= t <= MAX_T:
                record['status'] = 'below_liquid_domain' if t < MIN_T else 'above_screen_domain'
                rows.append(record)
                continue
            old = property_state(thermo,t,y,mixed=False)
            new = property_state(thermo,t,y,mixed=True)
            # This audits the independent forward reconstruction, not a new
            # tolerance for the original interpolation-based plume solve.
            rhoerr = abs(old['density_kg_m3']/sample['density']-1.)
            herr = abs(old['specific_enthalpy_J_kg']-float(exact_rho_h[i]/sample['density']))
            if rhoerr > 1e-7 or herr > 1e-5:
                raise ValueError('forward pure-phase reconstruction differs from exact old EOS')
            h0 = old['specific_enthalpy_J_kg']
            def objective(temp):
                return property_state(thermo,temp,y,mixed=True)['specific_enthalpy_J_kg']-h0
            same_h = None
            if objective(MIN_T)*objective(MAX_T) <= 0:
                nt = brentq(objective, MIN_T, MAX_T, xtol=1e-9)
                same_h = property_state(thermo,nt,y,mixed=True)
                same_h['temperature_shift_K'] = nt-t
                same_h['enthalpy_residual_J_kg'] = objective(nt)
            record.update(status='evaluated', old_pure_phase=old, new_ideal_solution=new,
                fixed_enthalpy_state=same_h, forward_replay_density_relative_error=rhoerr,
                forward_replay_enthalpy_error_J_kg=herr,
                delta_h_at_same_T_J_kg=new['specific_enthalpy_J_kg']-h0,
                relative_density_change_at_same_T=new['density_kg_m3']/old['density_kg_m3']-1.,
                delta_condensed_air_kg_kg=sum(new[k]-old[k] for k in ('condensed_N2_kg_kg','condensed_O2_kg_kg')))
            rows.append(record)
        print('Mixed liquid property screen trial',n,'done',flush=True)
    evaluated = [r for r in rows if r['status']=='evaluated']
    summary = dict(sample_count=len(rows), evaluated=len(evaluated),
        below_liquid_domain=sum(r['status']=='below_liquid_domain' for r in rows),
        above_screen_domain=sum(r['status']=='above_screen_domain' for r in rows),
        changed_condensed_air_states=sum(abs(r['delta_condensed_air_kg_kg'])>1e-10 for r in evaluated),
        maximum_added_condensed_air_kg_kg=max(r['delta_condensed_air_kg_kg'] for r in evaluated),
        maximum_absolute_enthalpy_change_J_kg=max(abs(r['delta_h_at_same_T_J_kg']) for r in evaluated),
        maximum_absolute_density_relative_change=max(abs(r['relative_density_change_at_same_T']) for r in evaluated),
        maximum_fixed_enthalpy_temperature_change_K=max(abs(r['fixed_enthalpy_state']['temperature_shift_K'])
            for r in evaluated if r['fixed_enthalpy_state'] is not None),
        maximum_lookup_temperature_error_K=max(abs(r['lookup_T_error_K']) for r in rows))
    if {name:digest(ROOT/name) for name in hashes} != hashes:
        raise ValueError('screen input changed')
    verify()
    write_new(args.output, dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(), sha256=hashes,
        benchmark=benchmark, summary=summary, rows=rows, original_replay_checks=checks, promoted=False,
        limitation='local bounded ideal-liquid phase test; no integrated field prediction, no below63.151K continuation'))
    import json
    print(json.dumps(dict(benchmark=benchmark['summary'],field_property_screen=summary),indent=2),flush=True)
    print('Saved:',args.output,flush=True)


if __name__ == '__main__':
    main()
