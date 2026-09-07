"""Reproducible frozen-phase audit and gamma-O2 property candidate screen."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import math

import numpy as np
import CoolProp as CP
from CoolProp.CoolProp import PropsSI

from stage2_matched_source import ROOT, read, digest, write_new
from stage1 import verify
from degali.addons.axisymmetric_jet import _air_phase_property_table, R_UNIVERSAL
from degali.addons.cryogenic_air import (
    _SOLID_VAPOUR_TABLES, _SOLID_SUBLIMATION_ENTHALPY, _SOLID_DENSITY,
    _TRIPLE_TEMPERATURE, _solid_vapour_pressure, _MOLECULAR_WEIGHT,
)
from degali.addons.oxygen_phase_potential import (
    MIN_T, MAX_T, MAX_P, solid, ideal_gas, pure_vapor_pressure,
    equilibrium_partial_pressure, ideal_gas_to_solid_enthalpy,
    nbs1977_vapor_pressure, triple_reference, FUSION_J_MOL,
    volume_and_derivatives,
)


def frozen_audit():
    output = {}
    grid = _air_phase_property_table()
    for species in ['Nitrogen','Oxygen']:
        mw, rho, old_l = _MOLECULAR_WEIGHT[species], _SOLID_DENSITY[species], _SOLID_SUBLIMATION_ENTHALPY[species]
        intervals = []
        table = _SOLID_VAPOUR_TABLES[species]
        for (t0,p0),(t1,p1) in zip(table[:-1],table[1:]):
            t = (t0+t1)/2
            p = _solid_vapour_pressure(species,t)
            slope = (math.log(p1)-math.log(p0))/(1/t0-1/t1)
            ideal_l = R_UNIVERSAL*slope/mw
            finite_v_l = ideal_l*(1-p*mw/(rho*R_UNIVERSAL*t))
            intervals.append(dict(lower_K=t0,upper_K=t1,old_latent_J_kg=old_l,
                slope_implied_ideal_latent_J_kg=ideal_l,
                slope_implied_with_old_solid_volume_J_kg=finite_v_l,
                relative_caloric_mismatch=old_l/finite_v_l-1))
        tt = _TRIPLE_TEMPERATURE[species]
        ps = _solid_vapour_pressure(species,tt)
        pl = PropsSI('P','T',tt,'Q',0,species)
        lv = PropsSI('H','T',tt,'Q',1,species)-PropsSI('H','T',tt,'Q',0,species)
        j = np.searchsorted(grid['temperature'],tt)
        grid_bounds = list(map(float,grid['temperature'][j-1:j+1]))
        liquid = []
        for t in sorted(set([tt,63.151,65.,70.,77.5,90.,100.,120.])):
            if t < tt:
                continue
            gas = CP.AbstractState('HEOS',species)
            gas.update(CP.DmolarT_INPUTS,1.,t)
            h0 = gas.hmass_idealgas()
            hv,hl = PropsSI('H','T',t,'Q',1,species),PropsSI('H','T',t,'Q',0,species)
            liquid.append(dict(T_K=t,old_real_vapor_latent_J_kg=hv-hl,
                ideal_reference_condensation_decrement_J_kg=h0-hl,
                old_condensed_enthalpy_bias_J_kg=h0-hv,
                pure_saturation_pressure_Pa=PropsSI('P','T',t,'Q',0,species)))
        output[species] = dict(pressure_intervals=intervals,
            maximum_absolute_interval_caloric_mismatch=max(abs(r['relative_caloric_mismatch']) for r in intervals),
            triple=dict(T_K=tt,old_solid_extrapolated_pressure_Pa=ps,liquid_pressure_Pa=pl,
                relative_pressure_jump_on_heating=pl/ps-1,
                implied_old_fusion_enthalpy_J_kg=old_l-lv,
                grid_transition_bounds_K=grid_bounds,grid_transition_width_K=grid_bounds[1]-grid_bounds[0],
                grid_has_explicit_triple_knot=bool(tt in grid['temperature'])),
            liquid_caloric_reference=liquid)
    # Independent check of the published formula, not a new fitted N2 EOS.
    output['Nitrogen']['nbs1955_analytic_solid'] = dict(A_log10_mmHg=7.65894,B_K=359.093,
        implied_ideal_latent_J_kg=R_UNIVERSAL*math.log(10)*359.093/_MOLECULAR_WEIGHT['Nitrogen'],
        omitted_62K_pressure_mmHg=73.6,
        frozen_62K_pressure_mmHg=_solid_vapour_pressure('Nitrogen',62)/133.322368,
        note='No calibrated replacement or solid caloric validity follows from differentiating a rounded pressure correlation.')
    return output


def candidate_audit(data):
    ref = triple_reference()
    rows = []
    for raw in data['table3']:
        t = raw['T_K']
        p = pure_vapor_pressure(t)
        h = ideal_gas_to_solid_enthalpy(t,p)
        old = _SOLID_SUBLIMATION_ENTHALPY['Oxygen']*_MOLECULAR_WEIGHT['Oxygen']
        published_p = nbs1977_vapor_pressure(t)
        rows.append(dict(**raw,candidate_pure_ideal_vapor_pressure_Pa=p,
            published_pressure_Pa=published_p,candidate_pressure_departure_fraction=p/published_p-1,
            candidate_ideal_sublimation_J_mol=h,old_constant_sublimation_J_mol=old,
            candidate_latent_departure_fraction=h/raw['boundary_latent_J_mol']-1,
            old_latent_departure_fraction=old/raw['boundary_latent_J_mol']-1))
    grid_rows = []
    for t in np.linspace(MIN_T,MAX_T,201):
        t = float(t)
        p = pure_vapor_pressure(t)
        s = solid(t,p)
        v,dv,_ = volume_and_derivatives(t)
        # Cp-CSAT = T*(dV/dT)_P*dP_sat/dT. Evaluate using the published
        # pressure equation on its own43.801--54.359K domain only.
        caloric_approximation = None
        if t <= 54.359:
            npres = nbs1977_vapor_pressure(t)
            dpdt = npres*(1096.562485/t**2-2.025578307/t)
            caloric_approximation = t*dv*dpdt/s.heat_capacity_J_mol_K
        dt = 1e-3
        clapeyron_error = None
        if MIN_T+dt <= t <= MAX_T-dt:
            g = ideal_gas(t,p)
            derivative = (pure_vapor_pressure(t+dt)-pure_vapor_pressure(t-dt))/(2*dt)
            implied_h = t*(g.volume_m3_mol-v)*derivative
            clapeyron_error = implied_h/(g.enthalpy_J_mol-s.enthalpy_J_mol)-1
        grid_rows.append(dict(T_K=t,pure_pressure_Pa=p,
            Cp_minus_CSAT_estimate_fraction=caloric_approximation,
            clapeyron_relative_residual=clapeyron_error,
            old_vs_new_ideal_solid_h_difference_J_kg=(ideal_gas_to_solid_enthalpy(t,101325)-_SOLID_SUBLIMATION_ENTHALPY['Oxygen']*ref['molecular_weight'])/ref['molecular_weight'],
            atmospheric_solid=asdict(solid(t,101325)),
            ideal_equilibrium_partial_pressure_at_1atm_Pa=equilibrium_partial_pressure(t,101325)))
    max_approx = max(abs(r['Cp_minus_CSAT_estimate_fraction']) for r in grid_rows if r['Cp_minus_CSAT_estimate_fraction'] is not None)
    max_cc = max(abs(r['clapeyron_relative_residual']) for r in grid_rows if r['clapeyron_relative_residual'] is not None)
    assert max_approx < .005, 'Cp approximation exceeds cited caloric error; do not use candidate'
    assert max_cc < 1e-7, 'single-potential Clapeyron regression failed'
    measured_fusion = FUSION_J_MOL/ref['molecular_weight']
    old_fusion = _SOLID_SUBLIMATION_ENTHALPY['Oxygen']-(PropsSI('H','T',MAX_T,'Q',1,'Oxygen')-PropsSI('H','T',MAX_T,'Q',0,'Oxygen'))
    summary = dict(table3_rows=len(rows),
        maximum_absolute_candidate_latent_departure=max(abs(r['candidate_latent_departure_fraction']) for r in rows),
        maximum_absolute_old_latent_departure=max(abs(r['old_latent_departure_fraction']) for r in rows),
        maximum_absolute_candidate_pressure_departure=max(abs(r['candidate_pressure_departure_fraction']) for r in rows),
        maximum_relative_Cp_CSAT_approximation_estimate=max_approx,
        maximum_clapeyron_relative_residual=max_cc,
        measured_fusion_J_kg=measured_fusion,implied_old_fusion_J_kg=old_fusion,
        implied_old_fusion_relative_departure=old_fusion/measured_fusion-1,
        maximum_absolute_old_to_new_solid_enthalpy_change_J_kg=max(abs(r['old_vs_new_ideal_solid_h_difference_J_kg']) for r in grid_rows))
    return dict(reference=ref,domain=dict(minimum_K=MIN_T,maximum_K=MAX_T,maximum_pressure_Pa=MAX_P),
        rows=rows,grid=grid_rows,summary=summary,
        interpretation='Shared-source correlation reproduction and single-potential consistency. No mixture solid equilibrium or new dispersion field score.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    verify()
    paths = [ROOT/p for p in [
        'reference/air_mixture/NBSIR77-859_oxygen_selected.json',
        'reference/air_mixture/NBSIR77-859_slush_1977.pdf',
        'reference/air_mixture/NBS_Circular_564_1955.pdf',
        'src/degali/addons/cryogenic_air.py','src/degali/addons/axisymmetric_jet.py',
        'src/degali/addons/oxygen_phase_potential.py','tests/test_oxygen_phase_potential.py',
        'tools/audit_phase_caloric_consistency.py','docs/prereg-phase-caloric-consistency.md']]
    hashes = {str(p.relative_to(ROOT)):digest(p) for p in paths}
    data = read(paths[0])
    result = dict(schema='stage2-phase-caloric-consistency-v1',created_utc=datetime.now(timezone.utc).isoformat(),
        frozen=frozen_audit(),candidate_gamma_oxygen=candidate_audit(data),input_sha256=hashes,
        candidate_promoted=False,stage1_unchanged=True,new_integrated_dispersion_field=False)
    assert hashes == {str(p.relative_to(ROOT)):digest(p) for p in paths}
    verify()
    write_new(args.output,result)
    print(result['candidate_gamma_oxygen']['summary'])


if __name__ == '__main__':
    main()
