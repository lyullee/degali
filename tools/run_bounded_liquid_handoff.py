"""Actual10/23 same-source boundary feasibility for bounded liquid air."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime,timezone
import math
from pathlib import Path

import numpy as np

from degali.addons.bounded_liquid_air import BoundedLiquidAir
from degali.addons.energy_crosswind import IndependentEnergyCrosswind
from degali.core.jetplume import JetIntegralFluxes
from run_preslhy_ambient_profile_audit import replay
from stage2_matched_source import ROOT,REF,read,digest,write_new
from stage1 import verify


def clone_with_liquid(original):
    if type(original) is not IndependentEnergyCrosswind:
        raise ValueError('do not silently discard another original model option')
    thermo = BoundedLiquidAir(original.thermodynamics)
    model = IndependentEnergyCrosswind(original.jetplume,thermo,
        quadrature_points=original.quadrature_points,energy_transport=original.energy_transport,
        houf_width_mapping=original.houf_width_mapping,ground_interaction=original.ground_interaction)
    model.velocity_shape_exponent = original.velocity_shape_exponent
    return model


def original_target_check(original,initial,entry,target_record):
    tv = target_record['target_fluxes']
    target = JetIntegralFluxes(total_mass=tv['total_mass'],contaminant_mass=tv['hydrogen'],
        momentum_x=tv['momentum_x'],momentum_z=tv['momentum_z'],energy=tv['energy'])
    values = original._as_array(target)
    residual = np.max(abs(original._as_array(original.integral_fluxes(initial))-values)/np.maximum(abs(values),1))
    t_error = abs(original.centre_temperature(initial)-target_record['target_temperature_K'])
    w_error = abs(math.sqrt(2*math.log(2)*initial[2])/target_record['target_halfwidth_m']-1)
    if (residual>1e-8 or abs(t_error-entry['temperature_residual_k'])>1e-7
            or abs(w_error-entry['width_residual'])>1e-10
            or target_record['station_m']!=entry['station_m']):
        raise ValueError('target record does not reproduce original same-source boundary')
    return target,dict(maximum_target_flux_replay_error=float(residual),
        original_temperature_residual_K=t_error,original_width_residual=w_error)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    verify()
    paths = [REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
        REF/'gaussian_enthalpy_profile_allflux_interface_2026-09-05.json',REF/'e35_reduced.json',
        ROOT/'src/degali/addons/bounded_liquid_air.py',ROOT/'src/degali/addons/liquid_phase_potential.py',
        ROOT/'src/degali/addons/mixed_air_liquid.py',ROOT/'src/degali/addons/oxygen_phase_potential.py',
        ROOT/'tests/test_bounded_liquid_air.py',ROOT/'docs/prereg-bounded-liquid-handoff.md',
        ROOT/'tools/run_preslhy_ambient_profile_audit.py',Path(__file__).resolve()]
    hashes = {str(p.relative_to(ROOT)):digest(p) for p in paths}
    field,targets = read(paths[0]),read(paths[1])
    trials = {r['trial']:r for r in read(paths[2])['trials']}
    rows = {}
    for n in (10,23):
        partial = args.output.with_name(args.output.stem+f'_trial{n}.json')
        if partial.exists(): raise FileExistsError(partial)
        print('Bounded liquid handoff:',n,'replaying original boundary',flush=True)
        trajectory,replay_check = replay(field,trials[n])
        original = trajectory.model
        initial = trajectory.result.states[0].copy()
        target,target_check = original_target_check(original,initial,field['interfaces'][str(n)],targets['interfaces'][str(n)])
        candidate = clone_with_liquid(original)
        row = dict(trial=n,original_replay_check=replay_check,target_check=target_check,
            original_initial_state=initial,target_fluxes=asdict(target),
            target_temperature_K=targets['interfaces'][str(n)]['target_temperature_K'],
            target_halfwidth_m=targets['interfaces'][str(n)]['target_halfwidth_m'],
            source_unchanged=True,accepted=False,failure_reasons=[])
        try:
            row['same_state_new_temperature_K']=candidate.centre_temperature(initial)
            before = candidate._as_array(candidate.integral_fluxes(initial))
            row['same_state_five_flux_changes']=before-original._as_array(target)
            print('Bounded liquid handoff:',n,'projecting conserved fluxes',flush=True)
            projection = candidate.project(target,initial)
            candidate.quadrature_points=projection.quadrature_points
            state = projection.state
            t = candidate.centre_temperature(state)
            local_t=candidate.profiles(state)[3]
            width=math.sqrt(2*math.log(2)*state[2])
            dt=abs(t-row['target_temperature_K'])
            dw=abs(width/row['target_halfwidth_m']-1)
            reasons=[]
            if not projection.success: reasons.append('five-flux projection/quadrature failed unchanged gates')
            if dt>2.: reasons.append('centre-temperature mismatch exceeds original2K')
            if dw>.05: reasons.append('hydrogen halfwidth mismatch exceeds original5%')
            row.update(projection=asdict(projection),projected_temperature_K=t,
                minimum_quadrature_temperature_K=float(np.min(local_t)),projected_halfwidth_m=width,
                temperature_residual_K=dt,halfwidth_residual=dw,accepted=not reasons,failure_reasons=reasons)
        except (ValueError,RuntimeError,FloatingPointError) as err:
            row['failure_reasons'].append(str(err))
        row.update(rejected_optimizer_or_evaluation_domain_probes=candidate.thermodynamics.domain_failures,
            last_domain_probe=candidate.thermodynamics.last_domain_failure,
            finished_utc=datetime.now(timezone.utc).isoformat())
        write_new(partial,dict(row=row,input_sha256=hashes))
        rows[str(n)]=row
        print('Bounded liquid handoff:',n,'accepted=',row['accepted'],'reasons=',row['failure_reasons'],flush=True)
    assert hashes=={str(p.relative_to(ROOT)):digest(p) for p in paths}
    verify()
    write_new(args.output,dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=hashes,selected_trials=[10,23],interfaces=rows,
        all_interfaces_accepted=all(r['accepted'] for r in rows.values()),
        candidate_promoted=False,new_integrated_dispersion_field=False,stage1_unchanged=True))
    print('Saved:',args.output,flush=True)


if __name__=='__main__': main()
