"""Bounded, coefficient-free correction of the BASE reporting adapter."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import numpy as np

from audit_stage1_existing_evidence import ROOT, REF, FIELD, read, digest, relative, serial, concentration_metrics, geometry_metrics
from audit_stage1_temperature import baseline_arguments
from audit_stage1_resolution import differences
from degali.validation.nearfield import hydrogen_jet, Trajectory, REACH
from degali.validation.reported_jet_trajectory import ReportedJetTrajectory
from run_preslhy_ambient_profile_audit import statistics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite single-image evidence')
    upstream_path = REF/'stage1_resolution_diagnostic_2026-09-06.json'
    upstream = read(upstream_path)
    paths = {ROOT/name for name in upstream['sha256']} | {upstream_path, Path(__file__).resolve(),
        ROOT/'src/degali/validation/reported_jet_trajectory.py',
        ROOT/'tests/test_reported_jet_trajectory.py', ROOT/'docs/prereg-stage1-single-ground-image.md'}
    for name, expected in upstream['sha256'].items():
        if digest(ROOT/name) != expected:
            raise ValueError(f'upstream changed: {name}')
    hashes = {relative(path): digest(path) for path in sorted(paths)}
    temperature = read(REF/'stage1_temperature_2026-09-06.json')
    field = read(FIELD)
    trials = {r['trial']: r for r in read(REF/'e35_reduced.json')['trials']}
    pairs, vertical, temperatures, checks, sensitivity = [], [], [], [], []
    for old_run in temperature['baseline_runs']:
        number = old_run['trial']
        trial = trials[number]
        jp, initial = hydrogen_jet(**baseline_arguments(trial))
        rows = np.asarray(old_run['output_rows'])
        legacy = Trajectory(jp.th.table, rows)
        fixed = ReportedJetTrajectory(jp.th.table, rows)
        fine = None
        if number in (10, 23):
            fine = ReportedJetTrajectory(jp.th.table, jp.run(initial, distmx=.1, smax=REACH).rows)
        errors, arc_changes, temp_changes = [], [], []
        for old in field['baseline']['pairs']:
            if old['trial'] != number:
                continue
            sensors = [s for s in trial['sensors'] if round(s['x'], 3) == old['x']]
            def predict(model):
                return max(model.concentration_at(old['x'], s['y'], s['z']) for s in sensors)
            errors.append(abs(predict(legacy)-old['predicted']))
            current = dict(old, predicted=predict(fixed))
            pairs.append(current)
            if fine is not None:
                arc_changes.append(dict(coarse=current, fine=dict(current, predicted=predict(fine))))
        for old in field['baseline']['vertical_profiles']:
            if old['trial'] != number:
                continue
            state = fixed.at(old['x'])
            errors.extend([abs(state.z-old['modelled_centre']), abs(state.sz-old['modelled_sigma_z'])])
            vertical.append(dict(old, modelled_centre=state.z, modelled_sigma_z=state.sz))
        for old in temperature['temperature_rows']:
            if old['trial'] != number:
                continue
            q = (old['x'], old['y'], old['z'])
            errors.append(abs(legacy.temperature_at(*q)-old['BASE_K']))
            current = dict(old, BASE_SINGLE_IMAGE_K=fixed.temperature_at(*q))
            temperatures.append(current)
            if fine is not None:
                temp_changes.append(dict(coarse=current, fine=dict(current, BASE_SINGLE_IMAGE_K=fine.temperature_at(*q))))
        if max(errors) > 1e-8:
            raise ValueError('uncorrected reconstruction/unchanged geometry differs from frozen BASE')
        checks.append(dict(trial=number, maximum_old_replay_and_geometry_error=max(errors),
                           reported_rows_reused=True, minimum_elevation_m=float(rows[:, 1].min())))
        if fine is not None:
            sensitivity.append(dict(trial=number, coarse_step_m=.2, fine_step_m=.1,
                concentration=differences(arc_changes, ('predicted',)),
                temperature=differences(temp_changes, ('BASE_SINGLE_IMAGE_K',))))
        print('single-image trial', number, 'old replay/geometry error', max(errors), flush=True)
    if (len(pairs), len(vertical), len(temperatures)) != (38, 17, 41):
        raise ValueError('fixed population changed')
    if {relative(path): digest(path) for path in sorted(paths)} != hashes:
        raise ValueError('inputs changed during audit')
    output = dict(completed=True, phase='stage1_BASE_single_ground_image', sha256=hashes,
        created_utc=datetime.now(timezone.utc).isoformat(), model='BASE_SINGLE_IMAGE',
        concentration_pairs=pairs, vertical_profiles=vertical, temperature_rows=temperatures,
        metrics=dict(**concentration_metrics(pairs), geometry=geometry_metrics(vertical)),
        temperature_summaries={basis: statistics(temperatures, 'BASE_SINGLE_IMAGE_K', f'observed_{basis}_K')
                              for basis in ('minimum', 'p05', 'median')},
        replay_checks=checks, corrected_adapter_step_sensitivity=sensitivity,
        core_changed=False, ODE_reintegrated_all_seven=False, touchdown_supported=False,
        no_new_fitted_parameters=True, historical_default_unchanged=True)
    args.output.write_text(json.dumps(output, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(metrics=output['metrics'], temperature=output['temperature_summaries']), indent=2), flush=True)


if __name__ == '__main__':
    main()
