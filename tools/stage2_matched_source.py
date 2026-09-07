"""Reproducible matched-flow/profile experiments, without editing Stage1."""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'tools')]
REF = ROOT / 'reference/preslhy'
TRIALS = (10, 11, 12, 22, 23, 24, 25)
SOURCE = REF / 'measured_pipe_source_2026-09-05.json'
REDUCED = REF / 'e35_reduced.json'
PROTOCOL = ROOT / 'docs/prereg-stage2-matched-source.md'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite_json(value):
    if isinstance(value, np.ndarray):
        return finite_json(value.tolist())
    if isinstance(value, np.generic):
        return finite_json(value.item())
    if isinstance(value, dict):
        return {str(k): finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_new(path, data):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(finite_json(data), stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def make_source_document(original, reduced, choice):
    """Change only the flow input, retain all measured T/P and provenance."""
    if choice not in ('coriolis', 'pressure_loss'):
        raise ValueError('unknown flow choice')
    rows = {int(t['trial']): t for t in reduced['trials']}
    if len(rows) != len(reduced['trials']):
        raise ValueError('duplicate reduced trial')
    numbers = [int(r['trial']) for r in original['trials']]
    if len(set(numbers)) != len(numbers) or set(numbers) != set(TRIALS):
        raise ValueError('source must contain the seven unique frozen trials')
    out = deepcopy(original)
    out['stage2_analysis_input'] = {
        'flow_choice': choice,
        'not_a_replacement_measurement': True,
        'loader_slot': 'pressure_loss_mass_flow_g_s contains the selected flow, see flow_choice',
        'original_protocol': original.get('protocol'),
    }
    out['protocol'] = 'docs/prereg-stage2-matched-source.md'
    for row in out['trials']:
        trial = rows[int(row['trial'])]
        corio = float(trial['flow_mean_gs'])
        original_corio = float(row['coriolis_mass_flow_g_s']['mean'])
        # The sealed reduced input stores g/s at three decimal places;
        # preserve exactly that BASE flow, while checking the raw reduction.
        if (not math.isfinite(original_corio) or corio != round(original_corio, 3)
                or trial['window'] != row['window']):
            raise ValueError('Coriolis reductions use inconsistent windows or values')
        pressure = float(row['pressure_loss_mass_flow_g_s'])
        if not all(math.isfinite(v) and v > 0 for v in (corio, pressure)):
            raise ValueError('both source flow estimates must be finite and positive')
        selected = corio if choice == 'coriolis' else pressure
        row['original_pressure_loss_mass_flow_g_s'] = pressure
        row['stage2_selected_mass_flow_g_s'] = selected
        row['stage2_coriolis_rounding_difference_g_s'] = corio - original_corio
        row['stage2_flow_choice'] = choice
        row['pressure_loss_mass_flow_g_s'] = selected
        if 'hem_nozzle' in row:
            row['historical_hem_nozzle_not_used'] = row.pop('hem_nozzle')
    return out


def prepare(directory, choice):
    from stage1 import verify
    verify()
    original, reduced = read(SOURCE), read(REDUCED)
    out = make_source_document(original, reduced, choice)
    inputs = [SOURCE, REDUCED, PROTOCOL, Path(__file__),
              REF / 'stage1_delivery_manifest_2026-09-06.json']
    hashes = {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in inputs}
    directory.mkdir(parents=True, exist_ok=False)
    source_path = directory / 'selected_source.json'
    write_new(source_path, out)
    write_new(directory / 'inputs.json', dict(
        completed=True, created_utc=datetime.now(timezone.utc).isoformat(),
        flow_choice=choice, sha256=hashes, source_sha256=digest(source_path)))
    print('Prepared isolated source input:', source_path, flush=True)


def verify_inputs(directory):
    from stage1 import verify
    verify()
    config = read(directory / 'inputs.json')
    for name, expected in config['sha256'].items():
        if digest(ROOT / name) != expected:
            raise ValueError(f'experiment input changed: {name}')
    if digest(directory / 'selected_source.json') != config['source_sha256']:
        raise ValueError('selected source changed')
    return config


def settings(directory, profile):
    if profile not in ('density', 'enthalpy'):
        raise ValueError('unknown profile')
    return dict(measured_source_table=directory / 'selected_source.json',
        measured_source_mode='full', measured_lh2_equilibrium_bound=True,
        hydrogen_spin_isomer='normal', fit_velocity_spreading=False,
        consistent_phase_ambient=True, downstream_thermodynamic_profile=profile,
        progress=lambda message: print(message, flush=True))


def require_checkpoint(checkpoint, config, profile):
    if (not checkpoint.get('completed') or not checkpoint.get('interfaces_accepted')
            or checkpoint.get('failures') or checkpoint.get('selected_trials') != list(TRIALS)
            or checkpoint.get('profile') != profile
            or checkpoint.get('source_sha256') != config['source_sha256']
            or checkpoint.get('flow_choice') != config['flow_choice']
            or checkpoint.get('fit_velocity_spreading') is not False):
        raise ValueError('missing or incompatible passed seven-case interface checkpoint')


def fixed_scores(result, numbers):
    from audit_stage1_existing_evidence import concentration_metrics, geometry_metrics
    from run_preslhy_ambient_profile_audit import statistics
    from degali.validation.nearfield import IndependentEnergyTrajectory
    frozen = read(REF / 'stage1_single_image_2026-09-06.json')
    trials = {t['trial']: t for t in read(REDUCED)['trials']}
    trajectories, failures = {}, dict(result.failures)
    for n in numbers:
        interface = result.handoffs.get(n)
        downstream = getattr(interface, 'downstream_result', None)
        if interface is None or not interface.accepted or downstream is None:
            failures[str(n)] = 'no accepted complete trajectory'
            continue
        if (not np.all(np.isfinite(downstream.states))
                or downstream.maximum_relative_balance_residual > 1e-5):
            failures[str(n)] = 'nonfinite states or failed conservation gate'
            continue
        trajectories[n] = IndependentEnergyTrajectory(interface, downstream)
    pairs, vertical, temperatures, missing = [], [], [], []
    for group, output in [('concentration_pairs', pairs), ('vertical_profiles', vertical),
                          ('temperature_rows', temperatures)]:
        for old in frozen[group]:
            n = old['trial']
            if n not in numbers:
                continue
            model = trajectories.get(n)
            state = None if model is None else model.state_at(old['x'])
            if state is None:
                missing.append(dict(group=group, trial=n, x=old['x'], channel=old.get('channel')))
                continue
            row = dict(old)
            if group == 'concentration_pairs':
                sensors = [s for s in trials[n]['sensors'] if round(s['x'], 3) == old['x']]
                value = max(model.concentration_at(old['x'], s['y'], s['z']) for s in sensors)
                if not math.isfinite(value) or value <= 0:
                    raise ValueError('invalid fixed concentration prediction')
                row['predicted'] = value
            elif group == 'vertical_profiles':
                row['modelled_centre'] = float(state[6])
                row['modelled_sigma_z'] = model.model.section_widths(state)[1]
            else:
                row['matched_K'] = model.temperature_at(old['x'], old['y'], old['z'])
                if not math.isfinite(row['matched_K']):
                    raise ValueError('nonfinite sensor temperature')
            output.append(row)
    complete = not failures and not missing
    out = dict(coverage_complete=complete, failures=failures, missing=missing,
        selected_trials=list(numbers), full_seven_case_claim=list(numbers) == list(TRIALS),
        concentration_pairs=pairs, vertical_profiles=vertical, temperature_rows=temperatures)
    if complete:
        out['metrics'] = dict(**concentration_metrics(pairs), geometry=geometry_metrics(vertical))
        out['temperature_summaries'] = {
            basis: statistics(temperatures, 'matched_K', f'observed_{basis}_K')
            for basis in ('minimum', 'p05', 'median')}
    return out


def run(directory, profile, phase, numbers, step):
    from degali.validation.nearfield import (
        independent_energy_interfaces_from_reduced, independent_energy_from_reduced)
    from run_preslhy_source_ablation import _interfaces, _field
    config = verify_inputs(directory)
    checkpoint_path = directory / f'interface_{profile}.json'
    name = f'field_{profile}_' + '-'.join(map(str, numbers)) + f'_step{step:g}.json'
    output = checkpoint_path if phase == 'interface' else directory / name
    if output.exists():
        raise FileExistsError(output)
    data = dict(completed=False, started_utc=datetime.now(timezone.utc).isoformat(),
        phase=phase, profile=profile, flow_choice=config['flow_choice'],
        source_sha256=config['source_sha256'], fit_velocity_spreading=False,
        promoted=False, protocol_sha256=digest(PROTOCOL))
    if phase == 'interface':
        result = independent_energy_interfaces_from_reduced(REDUCED, **settings(directory, profile))
        data.update(selected_trials=result.selected_trials, interfaces=_interfaces(result),
                    interfaces_accepted=result.all_interfaces_accepted, failures=result.failures)
    else:
        checkpoint = read(checkpoint_path)
        require_checkpoint(checkpoint, config, profile)
        result = independent_energy_from_reduced(REDUCED, **settings(directory, profile),
            trial_filter=numbers, maximum_step=step)
        data.update(field=_field(result), fixed_scores=fixed_scores(result, numbers),
                    maximum_step_m=step, interface_checkpoint_sha256=digest(checkpoint_path))
        data['field']['promoted'] = False
    verify_inputs(directory)
    data.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat())
    write_new(output, data)
    print('Saved:', output, flush=True)
    if phase == 'field':
        print(json.dumps(finite_json({k: v for k, v in data['fixed_scores'].items()
            if k not in ('concentration_pairs', 'vertical_profiles', 'temperature_rows')}), indent=2), flush=True)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('prepare', 'interface', 'field'))
    parser.add_argument('directory', type=Path)
    parser.add_argument('--flow', choices=('coriolis', 'pressure_loss'), default='coriolis')
    parser.add_argument('--profile', choices=('density', 'enthalpy'), default='enthalpy')
    parser.add_argument('--trials', nargs='+', type=int, default=[10, 23])
    parser.add_argument('--step', type=float, choices=(.01, .02), default=.02)
    args = parser.parse_args()
    if (not args.trials or len(set(args.trials)) != len(args.trials)
            or not set(args.trials).issubset(TRIALS)):
        parser.error('trials must be a unique subset of the fixed seven')
    numbers = [n for n in TRIALS if n in args.trials]
    directory = args.directory.resolve()
    if args.phase == 'prepare':
        prepare(directory, args.flow)
    else:
        run(directory, args.profile, args.phase, numbers, args.step)


if __name__ == '__main__':
    main()
