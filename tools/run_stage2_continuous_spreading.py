"""Test sub-metre spreading consistency with exactly frozen source states."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from stage2_matched_source import ROOT, REF, TRIALS, read, digest, write_new, fixed_scores
from run_preslhy_ambient_profile_audit import replay
from degali.addons.continuous_ambient_spreading import with_continuous_ambient_spreading
from stage1 import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--trials', nargs='+', type=int, default=[10, 23])
    parser.add_argument('--step', type=float, choices=(.01, .02), default=.02)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if not args.trials or len(set(args.trials)) != len(args.trials) or not set(args.trials).issubset(TRIALS):
        parser.error('trials must be a unique subset of the frozen seven')
    numbers = [n for n in TRIALS if n in args.trials]
    verify()
    source_path = REF / 'phase_ambient_consistency_field_complete_2026-09-05.json'
    module = ROOT / 'src/degali/addons/continuous_ambient_spreading.py'
    protocol = ROOT / 'docs/prereg-continuous-ambient-spreading.md'
    paths = [source_path, module, protocol, Path(__file__).resolve(),
        ROOT / 'tools/stage2_matched_source.py', ROOT / 'tools/run_preslhy_ambient_profile_audit.py']
    hashes = {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in paths}
    old_field = read(source_path)
    if (old_field['selected_trials'] != list(TRIALS) or not old_field['interfaces_accepted']
            or old_field['failures']):
        raise ValueError('original source does not have seven passed interfaces')
    trials = {t['trial']: t for t in read(REF / 'e35_reduced.json')['trials']}
    data = dict(completed=False, started_utc=datetime.now(timezone.utc).isoformat(),
        selected_trials=numbers, model='continuous_ambient_spreading_density',
        source_choice='pressure_loss', maximum_step_m=args.step,
        identical_initial_states=True, promoted=False, sha256=hashes, runs={}, failures={})
    interfaces = {}
    for n in numbers:
        print('Continuous spreading trial', n, 'replay original state', flush=True)
        try:
            original, check = replay(old_field, trials[n])
            entry = old_field['interfaces'][str(n)]
            if not entry['accepted']:
                raise ValueError('original individual interface failed')
            initial = original.result.states[0].copy()
            candidate = with_continuous_ambient_spreading(original.model)
            old_flux = original.model._as_array(original.model.integral_fluxes(initial))
            new_flux = candidate._as_array(candidate.integral_fluxes(initial))
            initial_error = float(np.max(abs(old_flux-new_flux)/np.maximum(abs(old_flux), 1.)))
            if (initial_error > 1e-12
                    or candidate.centre_temperature(initial) != original.model.centre_temperature(initial)
                    or candidate.section_widths(initial) != original.model.section_widths(initial)):
                raise ValueError('spreading option changed the frozen initial source')
            arc = entry['downstream']['arc_length']
            result = candidate.solve(initial, maximum_distance=arc[-1]-arc[0],
                maximum_step=args.step, relative_tolerance=1e-5)
            interfaces[n] = SimpleNamespace(model=candidate, accepted=True, downstream_result=result)
            data['runs'][n] = dict(original_replay_check=check,
                initial_flux_error=initial_error, original_initial_state=initial,
                downstream=vars(result), velocity_width_ratio=entry['lambda'],
                energy_quadrature_points=candidate.quadrature_points)
            print('Continuous spreading trial', n, 'done; balance',
                  result.maximum_relative_balance_residual, flush=True)
        except (ValueError, RuntimeError, FloatingPointError) as error:
            data['failures'][str(n)] = str(error)
            print('Continuous spreading trial', n, 'failed:', error, flush=True)
    data['fixed_scores'] = fixed_scores(SimpleNamespace(handoffs=interfaces, failures=data['failures']), numbers)
    if {name: digest(ROOT/name) for name in hashes} != hashes:
        raise ValueError('input/code changed during candidate integration')
    verify()
    data.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat())
    write_new(args.output, data)
    import json
    print(json.dumps({k: v for k, v in data['fixed_scores'].items()
        if k not in ('concentration_pairs', 'vertical_profiles', 'temperature_rows')}, indent=2), flush=True)
    print('Saved:', args.output, flush=True)


if __name__ == '__main__':
    main()
