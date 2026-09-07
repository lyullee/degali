"""Reuse frozen shear integrals to bound jointly omitted normal work/TKE."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from degali.addons.normal_stress_budget import normal_work_interval, minimum_tke_for_normal_work_cap


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite joint covariance budget evidence')
    ref = ROOT/'reference/preslhy'
    origin = ref/'initial_shear_tke_bounds_2026-09-06.json'
    upstream = json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or len(upstream['rows']) != 7 or not all(r['resolved'] for r in upstream['rows']):
        raise ValueError('seven resolved fixed initial shear integrals required')
    for name, digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'frozen initial dependency changed: {name}')
    paths = [ROOT/name for name in upstream['sha256']] + [origin, Path(__file__).resolve()]
    paths += [ROOT/name for name in ('src/degali/addons/normal_stress_budget.py',
        'tests/test_normal_stress_budget.py', 'docs/prereg-joint-tke-normal-stress-budget.md')]
    def hashes():
        return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    seeds = {r['trial']: r for r in read('edge_conservative_refit_combined_2026-09-05.json')['rows']}
    measured = {r['trial']: r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    boundary = read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json')
    reduced = read('e35_reduced.json')['trials']
    data = dict(phase='joint_initial_tke_normal_stress_budgets', completed=False, rows=[], sha256=initial,
        started_utc=datetime.now(timezone.utc).isoformat(), adopted=False, field_scored=False,
        budgets_are_illustrations_not_initial_values=True, axial_derivatives_not_bounded=True,
        pressure_compensation_not_solved=True)
    def checkpoint():
        partial.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    checkpoint()
    caches = {}
    for inherited in upstream['rows']:
        trial = inherited['trial']
        row = dict(trial=trial, resolved=False)
        data['rows'].append(row)
        try:
            # Reconstruction checks positivity analytically, without repeating
            # the previously completed million-node shear integration.
            p, replay = construct_projection(seeds[trial], boundary, reduced, measured, caches)
            m = p.mixing
            min_u = m.wind*math.cos(m.state[3])+m.state[4]*math.exp(-p.section.velocity_shape_exponent*2*p.q0)
            if m.state[4] <= 0. or min_u <= 0.:
                raise ValueError('positive advective measure not established throughout the square')
            row.update(minimum_axial_velocity=min_u, source_reconstruction=replay, evaluations=[])
            for ev in inherited['evaluations']:
                b, mean = ev['minimum_tke_axial_flux'], ev['mean_axial_kinetic_flux']
                normal = []
                for cap in (.12, .20, .50):
                    interval = normal_work_interval(b, cap*mean)
                    normal.append(dict(tke_fraction_cap=cap, feasible=interval['feasible'],
                        normal_work_fraction_minimum=interval['minimum']/mean if interval['feasible'] else None))
                kinetic = [dict(normal_work_fraction_cap=cap,
                    tke_fraction_minimum=minimum_tke_for_normal_work_cap(b, cap*mean)/mean) for cap in (.01, .05, .10)]
                row['evaluations'].append(dict(order=ev['order'], mean_kinetic_flux=mean,
                    shear_bound_flux=b, normal_work_bounds=normal, tke_bounds=kinetic))
            coarse, fine = row['evaluations']
            errors = []
            for key, field in (('normal_work_bounds', 'normal_work_fraction_minimum'), ('tke_bounds', 'tke_fraction_minimum')):
                for a, b in zip(coarse[key], fine[key]):
                    if a[field] is not None and b[field] is not None:
                        errors.append(abs(a[field]-b[field])/max(abs(b[field]), 1.))
                    elif a[field] != b[field]:
                        raise ValueError('coarse/fine feasibility classification differs')
            row['maximum_bound_refinement'] = max(errors)
            row['resolved'] = bool(max(errors) <= 1e-3)
            print('joint budget', trial, row['resolved'], fine['tke_bounds'][0], flush=True)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row['failure'] = str(exc)
        checkpoint()
    if hashes() != initial:
        raise RuntimeError('joint covariance budget dependencies changed')
    data.update(completed=True, resolved_trials=[r['trial'] for r in data['rows'] if r['resolved']],
                finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('joint budgets completed', data['resolved_trials'], flush=True)


if __name__ == '__main__':
    main()
