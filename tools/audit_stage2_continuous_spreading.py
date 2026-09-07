"""Verify and aggregate the pre-registered, fixed-source spreading screen."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from stage2_matched_source import ROOT, REF, TRIALS, read, digest, write_new
from audit_stage1_existing_evidence import paired_rows, concentration_metrics, geometry_metrics
from run_preslhy_ambient_profile_audit import replay, statistics
from degali.addons.continuous_ambient_spreading import with_continuous_ambient_spreading
from degali.validation.nearfield import IndependentEnergyTrajectory
from stage1 import verify


ROW_SPECS = {
    'concentration_pairs': (('trial', 'x'), ('observed',), ('predicted',)),
    'vertical_profiles': (('trial', 'x'), ('measured_centre', 'measured_sigma_z'),
                          ('modelled_centre', 'modelled_sigma_z')),
    'temperature_rows': (('trial', 'channel', 'x', 'y', 'z'),
                         ('observed_minimum_K', 'observed_p05_K', 'observed_median_K'), ('matched_K',)),
}


def verified_groups(runs, frozen):
    selected = [n for run in runs for n in run['selected_trials']]
    if len(selected) != 7 or set(selected) != set(TRIALS):
        raise ValueError('expected seven unique trials')
    for run in runs:
        if (not run['completed'] or run['failures'] or not run['fixed_scores']['coverage_complete']
                or not run['identical_initial_states'] or run['source_choice'] != 'pressure_loss'
                or run['model'] != 'continuous_ambient_spreading_density' or run['maximum_step_m'] != .02):
            raise ValueError('incomplete or confounded candidate')
        if {name: digest(ROOT/name) for name in run['sha256']} != run['sha256']:
            raise ValueError('candidate input/code hash changed')
    groups = {key: [r for run in runs for r in run['fixed_scores'][key]] for key in ROW_SPECS}
    for key, (keys, observed, _) in ROW_SPECS.items():
        paired_rows(groups[key], frozen[key], keys, observed)
    if tuple(map(len, groups.values())) != (38, 17, 41):
        raise ValueError('wrong fixed populations')
    return groups


def resolution(coarse, fine):
    if (coarse['selected_trials'] != fine['selected_trials'] or not fine['completed']
            or fine['failures'] or not fine['fixed_scores']['coverage_complete']
            or fine['maximum_step_m'] != .01 or fine['sha256'] != coarse['sha256']):
        raise ValueError('fine-step comparison changed more than resolution or is incomplete')
    result = {}
    for group, (keys, observed, predicted) in ROW_SPECS.items():
        a, b = (run['fixed_scores'][group] for run in (coarse, fine))
        paired_rows(a, b, keys, observed)
        index = {tuple(row[k] for k in keys): row for row in b}
        result[group] = {}
        for column in predicted:
            delta = np.array([index[tuple(row[k] for k in keys)][column]-row[column] for row in a])
            reference = np.array([row[column] for row in a])
            result[group][column] = dict(n=len(delta), maximum_absolute_change=float(max(abs(delta))),
                maximum_relative_change=float(max(abs(delta)/np.maximum(abs(reference), 1e-12))))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    verify()
    paths = [REF / f'continuous_spreading_pressure_{s}_2026-09-06.json' for s in
             ('10-23_step0.02', '11-12-22-24-25_step0.02', '10-23_step0.01')]
    paths += [REF/'stage1_single_image_2026-09-06.json',
              REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
              REF/'e35_reduced.json', Path(__file__).resolve()]
    hashes = {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in paths}
    runs = [read(p) for p in paths[:3]]
    frozen, control = read(paths[3]), read(paths[4])
    groups = verified_groups(runs[:2], frozen)
    for key, control_key in [('concentration_pairs', 'pairs'), ('vertical_profiles', 'vertical_profiles')]:
        keys, observed, _ = ROW_SPECS[key]
        paired_rows(groups[key], control['candidate'][control_key], keys, observed)
    def metrics(pairs, vertical):
        return dict(**concentration_metrics(pairs), geometry=geometry_metrics(vertical))
    comparison = dict(control=metrics(control['candidate']['pairs'], control['candidate']['vertical_profiles']),
                      continuous_spreading=metrics(groups['concentration_pairs'], groups['vertical_profiles']))
    temperatures = {basis: {label: statistics(groups['temperature_rows'], column, f'observed_{basis}_K')
                           for label, column in [('control', 'CONTROL_K'), ('continuous_spreading', 'matched_K')]}
                    for basis in ('minimum', 'p05', 'median')}
    trials = {t['trial']: t for t in read(paths[5])['trials']}
    checks = {}
    for run in runs:
        for n in run['selected_trials']:
            original, original_check = replay(control, trials[n])
            model = with_continuous_ambient_spreading(original.model)
            stored = run['runs'][str(n)]['downstream']
            states = np.asarray(stored['states'])
            if not np.array_equal(states[0], original.result.states[0]):
                raise ValueError('candidate did not start at the frozen state')
            trajectory = IndependentEnergyTrajectory(SimpleNamespace(model=model), SimpleNamespace(states=states))
            ts, _ = model.thermodynamics._condensed_air_state(states[:, 0], states[:, 1])
            terr = float(max(abs(ts-np.asarray(stored['temperatures']))))
            indices = np.unique(np.linspace(0, len(states)-1, 7, dtype=int))
            ferr = max(float(max(abs(model._as_array(model.integral_fluxes(states[i]))-stored['fluxes'][i])
                                 /np.maximum(abs(np.asarray(stored['fluxes'][i])), 1.))) for i in indices)
            cerrs, sensor_errors = [], []
            for row in run['fixed_scores']['concentration_pairs']:
                if row['trial'] == n:
                    sensors = [s for s in trials[n]['sensors'] if round(s['x'], 3) == row['x']]
                    predicted = max(trajectory.concentration_at(row['x'], s['y'], s['z']) for s in sensors)
                    cerrs.append(abs(predicted/row['predicted']-1.))
            for row in run['fixed_scores']['temperature_rows']:
                if row['trial'] == n:
                    sensor_errors.append(abs(trajectory.temperature_at(row['x'], row['y'], row['z'])-row['matched_K']))
            cerr, serr = max(cerrs), max(sensor_errors, default=0.)
            if terr > 1e-6 or ferr > 1e-9 or cerr > 1e-9 or serr > 1e-6:
                raise ValueError('candidate replay mismatch')
            checks[f'{n}_step{run["maximum_step_m"]}'] = dict(temperature_max_abs_K=terr,
                sampled_flux_max_scaled_error=ferr, arc_max_relative_error=cerr,
                sensor_max_abs_K=serr, original_replay=original_check,
                balance_residual=stored['maximum_relative_balance_residual'])
            print('Replayed candidate', n, run['maximum_step_m'], flush=True)
    if {name: digest(ROOT/name) for name in hashes} != hashes:
        raise ValueError('input changed during audit')
    verify()
    result = dict(completed=True, created_utc=datetime.now(timezone.utc).isoformat(), sha256=hashes,
        groups=groups, metrics=comparison, temperature=temperatures,
        resolution=resolution(runs[0], runs[2]), replay_checks=checks, promoted=False,
        interpretation='fixed-source, coefficient-free consistency screen; sub-metre atmospheric law remains unvalidated')
    write_new(args.output, result)
    import json
    print(json.dumps({k: result[k] for k in ('metrics', 'temperature', 'resolution')}, indent=2), flush=True)
    print('Saved:', args.output, flush=True)


if __name__ == '__main__':
    main()
