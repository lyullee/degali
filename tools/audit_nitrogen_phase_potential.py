"""Corrected original nitrogen caloric benchmark; no dispersion score."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import math

import numpy as np
import CoolProp as CP
from scipy.optimize import brentq

from stage1 import verify
from stage2_matched_source import ROOT, read, digest, write_new
from degali.addons import nitrogen_phase_potential as n2
from degali.addons.cryogenic_air import (
    _SOLID_SUBLIMATION_ENTHALPY, _SOLID_DENSITY, _MOLECULAR_WEIGHT,
)


def real_vapor_equilibrium(t):
    """Diagnostic pure HEOS vapor on candidate solid; never a mixture EOS."""
    gas = CP.AbstractState('HEOS','Nitrogen')
    gas.specify_phase(CP.iphase_gas)
    pideal = n2.pure_vapor_pressure(t)
    def residual(logp):
        p = math.exp(logp)
        gas.update(CP.PT_INPUTS,p,t)
        return gas.gibbsmolar()-n2.solid(t,p).gibbs_J_mol
    root = brentq(residual,math.log(pideal*.5),math.log(pideal*2),xtol=2e-13)
    p = math.exp(root)
    residual(root)
    solid = n2.solid(t,p)
    return dict(pressure_Pa=p,latent_J_mol=gas.hmolar()-solid.enthalpy_J_mol,
        chemical_potential_residual_J_mol=residual(root),
        gas_Z=gas.compressibility_factor())


def audit(data):
    assert tuple(v*4.184 for v in data['beta_Cp_coefficients_cal_mol_K']) == n2.CP_COEFFICIENTS
    ref = n2.triple_reference()
    old_latent = _SOLID_SUBLIMATION_ENTHALPY['Nitrogen']*_MOLECULAR_WEIGHT['Nitrogen']
    rows = []
    for raw in data['solid_table_XVI_and_corrections']:
        t = raw['T_K']
        if not n2.MIN_T <= t <= n2.MAX_T:
            rows.append(dict(**raw,status='outside_current_reference_domain'))
            continue
        p = n2.pure_vapor_pressure(t)
        h = n2.ideal_gas_to_solid_enthalpy(t,p)
        real = real_vapor_equilibrium(t)
        target_p,target_h = raw['P_mmHg']*133.322368,raw['latent_cal_mol']*4.184
        rows.append(dict(**raw,status='evaluated',ideal_pressure_Pa=p,ideal_latent_J_mol=h,
            real_vapor_same_solid=real,old_constant_latent_J_mol=old_latent,
            old_latent_relative_departure=old_latent/target_h-1,
            ideal_latent_relative_departure=h/target_h-1,
            real_latent_relative_departure=real['latent_J_mol']/target_h-1,
            ideal_pressure_relative_departure=p/target_p-1,
            real_pressure_relative_departure=real['pressure_Pa']/target_p-1))
    local = []
    for t in [59.,60.,61.,62.,63.]:
        s,g = n2.solid(t,101325.),n2.ideal_gas(t,100.)
        local.append(dict(T_K=t,candidate_solid=asdict(s),
            ideal_reference_latent_J_kg=(g.enthalpy_J_mol-s.enthalpy_J_mol)/ref['molecular_weight'],
            candidate_minus_old_solid_h_J_kg=(s.enthalpy_J_mol-(g.enthalpy_J_mol-old_latent))/ref['molecular_weight'],
            candidate_minus_old_solid_Cp_J_kg_K=(s.heat_capacity_J_mol_K-g.heat_capacity_J_mol_K)/ref['molecular_weight'],
            candidate_density_kg_m3=ref['molecular_weight']/s.volume_m3_mol,
            old_density_kg_m3=_SOLID_DENSITY['Nitrogen']))
    density = [dict(**raw,candidate_mean_volume_density_kg_m3=ref['molecular_weight']/n2.VOLUME_M3_MOL,
        candidate_relative_departure=(ref['molecular_weight']/n2.VOLUME_M3_MOL)/raw['rho_kg_m3']-1)
        for raw in data['beta_density_measurements_as_reported']]
    cc = []
    for t in np.linspace(n2.MIN_T+.002,n2.MAX_T-.002,301):
        t = float(t)
        p = n2.pure_vapor_pressure(t)
        s,g = n2.solid(t,p),n2.ideal_gas(t,p)
        dpdt = (n2.pure_vapor_pressure(t+.001)-n2.pure_vapor_pressure(t-.001))/.002
        cc.append(t*(g.volume_m3_mol-s.volume_m3_mol)*dpdt/(g.enthalpy_J_mol-s.enthalpy_J_mol)-1)
    assert max(map(abs,cc)) < 1e-7
    complete = [r for r in rows if r['status']=='evaluated']
    summary = {key: max(abs(r[key]) for r in complete) for key in [
        'old_latent_relative_departure','ideal_latent_relative_departure',
        'real_latent_relative_departure','ideal_pressure_relative_departure','real_pressure_relative_departure']}
    summary.update(evaluated_rows=len(complete),excluded_rows=len(rows)-len(complete),
        corrected_rows=sum(r['corrected'] for r in complete),
        maximum_clapeyron_residual=max(map(abs,cc)),
        measured_fusion_J_kg=n2.FUSION_J_MOL/ref['molecular_weight'],
        maximum_59_63K_solid_h_change_J_kg=max(abs(r['candidate_minus_old_solid_h_J_kg']) for r in local),
        maximum_59_63K_solid_Cp_change_J_kg_K=max(abs(r['candidate_minus_old_solid_Cp_J_kg_K']) for r in local))
    return dict(reference=ref,rows=rows,local_pure_material_changes=local,
        density_approximation_check=density,summary=summary,
        interpretation='Shared-source caloric correlation reproduction and pure-material consistency. No new integrated dispersion predictions. No claim of accurate mixed solids or constant solid density.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    verify()
    paths = [ROOT/p for p in [
        'reference/air_mixture/GeorgiaTech_A663_nitrogen_selected.json',
        'reference/air_mixture/GeorgiaTech_A663_1963.pdf',
        'src/degali/addons/nitrogen_phase_potential.py',
        'src/degali/addons/oxygen_phase_potential.py',
        'src/degali/addons/cryogenic_air.py',
        'tests/test_nitrogen_phase_potential.py',
        'tools/audit_nitrogen_phase_potential.py',
        'docs/prereg-nitrogen-phase-potential.md']]
    hashes = {str(p.relative_to(ROOT)):digest(p) for p in paths}
    result = dict(schema='stage2-beta-nitrogen-potential-v1',created_utc=datetime.now(timezone.utc).isoformat(),
        candidate=audit(read(paths[0])),input_sha256=hashes,
        candidate_promoted=False,new_integrated_dispersion_field=False,stage1_unchanged=True)
    assert hashes == {str(p.relative_to(ROOT)):digest(p) for p in paths}
    verify()
    write_new(args.output,result)
    print(result['candidate']['summary'])


if __name__ == '__main__': main()
