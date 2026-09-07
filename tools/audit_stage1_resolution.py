"""Two-case output/integration-step sensitivity; never replace full-field scores."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import numpy as np
from audit_stage1_existing_evidence import ROOT, REF, FIELD, digest, read, relative, serial, concentration_metrics, geometry_metrics
from audit_stage1_temperature import baseline_arguments
from run_preslhy_ambient_profile_audit import replay, statistics
from degali.validation.nearfield import hydrogen_jet, Trajectory, REACH


def differences(values, field_names):
    output = dict(rows=len(values))
    for name in field_names:
        changes = np.array([r['fine'][name]-r['coarse'][name] for r in values])
        scale = np.array([max(abs(r['coarse'][name]), 1e-12) for r in values])
        output[name] = dict(maximum_absolute_change=float(max(abs(changes))),
                            mean_absolute_change=float(np.mean(abs(changes))),
                            maximum_relative_change=float(max(abs(changes)/scale)))
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite resolution diagnosis')
    inventory_path = REF/'stage1_existing_evidence_2026-09-06.json'
    inventory = read(inventory_path)
    temperature_path = REF/'stage1_temperature_2026-09-06.json'
    temperature = read(temperature_path)
    paths = {ROOT/name for name in temperature['sha256']} | {temperature_path, inventory_path,
              Path(__file__).resolve(), ROOT/'docs/prereg-stage1-resolution-diagnostic.md'}
    for name, expected in temperature['sha256'].items():
        if digest(ROOT/name) != expected:
            raise ValueError(f'upstream input/code changed: {name}')
    for n in (10, 23):
        paths.add(REF/f'stage1_resolution_{n}_2026-09-06.json')
    hashes = {relative(path): digest(path) for path in sorted(paths)}
    field = read(FIELD)
    trials = {r['trial']: r for r in read(REF/'e35_reduced.json')['trials']}
    output = dict(phase='stage1_two_case_step_sensitivity', completed=False, sha256=hashes,
                  started_utc=datetime.now(timezone.utc).isoformat(), cases=[])
    for number in (10, 23):
        trial = trials[number]
        fine_field = read(REF/f'stage1_resolution_{number}_2026-09-06.json')
        if fine_field['selected_trials'] != [number] or not fine_field['interfaces_accepted'] or fine_field['failures']:
            raise ValueError('fine candidate failed source/downstream integration')
        candidate, replay_checks = replay(fine_field, trial)
        jp, initial = hydrogen_jet(**baseline_arguments(trial))
        base = Trajectory(jp.th.table, jp.run(initial, distmx=.1, smax=REACH).rows)
        if not base.ok:
            raise ValueError('fine BASE run failed')
        for variant, model, old_kind in (('BASE', base, 'baseline'), ('CANDIDATE', candidate, 'candidate')):
            arc_rows, vertical_rows, temperature_rows = [], [], []
            for old in field[old_kind]['pairs']:
                if old['trial'] != number:
                    continue
                if model.at(old['x']) is None:
                    raise ValueError('fine model misses a fixed concentration arc')
                sensors = [s for s in trial['sensors'] if round(s['x'], 3) == old['x']]
                value = max(model.concentration_at(old['x'], s['y'], s['z']) for s in sensors)
                arc_rows.append(dict(coarse=old, fine=dict(old, predicted=value)))
            for old in field[old_kind]['vertical_profiles']:
                if old['trial'] != number:
                    continue
                state = model.at(old['x'])
                if state is None:
                    raise ValueError('fine model misses a fixed vertical profile')
                vertical_rows.append(dict(coarse=old, fine=dict(old, modelled_centre=state.z, modelled_sigma_z=state.sz)))
            for old in temperature['temperature_rows']:
                if old['trial'] != number:
                    continue
                if model.at(old['x']) is None:
                    raise ValueError('fine model misses a fixed thermocouple')
                value = model.temperature_at(old['x'], old['y'], old['z'])
                coarse = dict(old, predicted_K=old[f'{variant}_K'])
                temperature_rows.append(dict(coarse=coarse, fine=dict(coarse, predicted_K=value)))
            case = dict(trial=number, variant=variant,
                coarse_step_m=.2 if variant == 'BASE' else .02,
                fine_step_m=.1 if variant == 'BASE' else .01,
                concentration_changes=differences(arc_rows, ('predicted',)),
                vertical_changes=differences(vertical_rows, ('modelled_centre', 'modelled_sigma_z')),
                temperature_changes=differences(temperature_rows, ('predicted_K',)),
                concentration_rows=arc_rows, vertical_rows=vertical_rows, temperature_rows=temperature_rows,
                metrics={level: dict(concentration=concentration_metrics([r[level] for r in arc_rows]),
                    geometry=geometry_metrics([r[level] for r in vertical_rows]),
                    temperature={basis: statistics([r[level] for r in temperature_rows], 'predicted_K', f'observed_{basis}_K')
                                 for basis in ('minimum', 'p05', 'median')}) for level in ('coarse', 'fine')})
            if variant == 'CANDIDATE':
                case['fine_replay_checks'] = replay_checks
            output['cases'].append(case)
            print(number, variant, case['concentration_changes'], case['temperature_changes'], flush=True)
    if {relative(path): digest(path) for path in sorted(paths)} != hashes:
        raise ValueError('input/code changed during resolution diagnosis')
    output.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
        seven_case_convergence_proven=False, original_scores_replaced=False, new_physics_fitted=False)
    args.output.write_text(json.dumps(output, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('saved', args.output, flush=True)


if __name__ == '__main__':
    main()
