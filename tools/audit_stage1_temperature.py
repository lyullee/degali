"""Fresh corrected-BASE runs and the same41 sensor temperature comparison."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time
import numpy as np
from audit_stage1_existing_evidence import ROOT, REF, FIELD, CONTROL, TEMPERATURE, digest, read, relative, serial
from run_preslhy_ambient_profile_audit import replay, statistics
from degali.validation.nearfield import hydrogen_jet, Trajectory, STEP, REACH, _liquid_source_temperature


def baseline_arguments(trial):
    return dict(rate=trial['flow_mean_gs']/1000., diameter=trial['orifice_mm']/1000.,
        wind=trial['wind_ms'], height=trial['release_height_m'], ambient_temperature=trial['T_C']+273.15,
        relative_humidity=trial['RH_pct'], storage_pressure_barg=trial['tanker_barg'],
        storage_temperature=_liquid_source_temperature(trial, 'tank_saturation'),
        wind_reference_height=trial['wind_ref_m'], corrections=True, ground_effect=False,
        source_table_consistency=True, source_momentum_consistency=True)


def checked_baseline(trial, field):
    settings = baseline_arguments(trial)
    jp, initial = hydrogen_jet(**settings)
    started = time.perf_counter()
    rows = jp.run(initial, distmx=STEP, smax=REACH).rows
    trajectory = Trajectory(jp.th.table, rows)
    if not trajectory.ok or not np.all(np.isfinite(rows)):
        raise ValueError('baseline run failed or returned nonfinite rows')
    number = trial['trial']
    pairs, vertical, errors = [], [], []
    for expected in field['baseline']['pairs']:
        if expected['trial'] != number:
            continue
        x = expected['x']
        if trajectory.at(x) is None:
            raise ValueError(f'BASE fails to reach concentration arc {number}/{x}')
        sensors = [s for s in trial['sensors'] if round(s['x'], 3) == x]
        value = max(trajectory.concentration_at(x, s['y'], s['z']) for s in sensors)
        errors.append(abs(value-expected['predicted'])/max(abs(expected['predicted']), 1.))
        pairs.append(dict(expected, predicted=value, plume_centre=trajectory.at(x).z))
    for expected in field['baseline']['vertical_profiles']:
        if expected['trial'] != number:
            continue
        state = trajectory.at(expected['x'])
        if state is None:
            raise ValueError('BASE fails to reach a vertical profile')
        for value, key in ((state.z, 'modelled_centre'), (state.sz, 'modelled_sigma_z')):
            errors.append(abs(value-expected[key])/max(abs(expected[key]), 1.))
        vertical.append(dict(expected, modelled_centre=state.z, modelled_sigma_z=state.sz))
    if not pairs or max(errors) > 1e-8:
        raise ValueError(f'BASE differs from frozen concentration/geometry: {max(errors)}')
    out = dict(trial=number, settings=settings, output_step_m=STEP, arclength_limit_m=REACH,
        output_rows=rows, output_columns=['x_m', 'centre_z_m', 'C_H2_kg_m3', 'sigma_y_m', 'sigma_z_m'],
        concentration_pairs=pairs, vertical_profiles=vertical,
        maximum_reproduction_error=max(errors), x_range_m=[float(min(rows[:, 0])), float(max(rows[:, 0]))],
        elapsed_seconds=time.perf_counter()-started)
    return trajectory, out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite stage1 temperature evidence')
    inventory_path = REF/'stage1_existing_evidence_2026-09-06.json'
    inventory = read(inventory_path)
    paths = {ROOT/name for name in inventory['sha256']}
    for name, expected in inventory['sha256'].items():
        if digest(ROOT/name) != expected:
            raise ValueError(f'current stage1 manifest changed: {name}')
    paths.update([inventory_path, Path(__file__).resolve(), ROOT/'tools/audit_stage1_existing_evidence.py',
                  ROOT/'tools/run_preslhy_ambient_profile_audit.py'])
    hashes = {relative(path): digest(path) for path in sorted(paths)}
    field, old_temperature = read(FIELD), read(TEMPERATURE)
    trials = {t['trial']: t for t in read(REF/'e35_reduced.json')['trials']}
    data = dict(phase='stage1_paired_temperature_and_BASE_reproduction', completed=False,
                started_utc=datetime.now(timezone.utc).isoformat(), sha256=hashes,
                baseline_runs=[], temperature_rows=[], candidate_replay_checks={},
                candidate_reintegrated_by_this_tool=False, CONTROL_reintegrated=False)
    def checkpoint():
        partial.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    checkpoint()
    for number in field['selected_trials']:
        print(f'BASE trial {number}: fresh integration', flush=True)
        base, result = checked_baseline(trials[number], field)
        data['baseline_runs'].append(result)
        print(f'BASE trial {number}: {len(result["concentration_pairs"])} arcs, '
              f'error={result["maximum_reproduction_error"]:.3e}', flush=True)
        checkpoint()
        if number not in (10, 23):
            continue
        candidate, check = replay(field, trials[number])
        if candidate is None:
            raise ValueError('required candidate temperature trajectory is missing')
        data['candidate_replay_checks'][str(number)] = check
        for old in old_temperature['temperature_rows']:
            if old['trial'] != number:
                continue
            if base.at(old['x']) is None or candidate.state_at(old['x']) is None:
                raise ValueError('fixed sensor is outside a BASE/candidate trajectory; do not use ambient fallback')
            q = (old['x'], old['y'], old['z'])
            tb, tc = base.temperature_at(*q), candidate.temperature_at(*q)
            if not np.isfinite(tb) or not np.isfinite(tc) or abs(tc-old['candidate_K']) > 1e-6:
                raise ValueError('candidate temperature does not reproduce frozen sensor prediction')
            row = {k: v for k, v in old.items() if k not in ('control_K', 'candidate_K')}
            row.update(BASE_K=tb, CONTROL_K=old['control_K'], CANDIDATE_K=tc)
            data['temperature_rows'].append(row)
        checkpoint()
    rows = data['temperature_rows']
    sensor_keys = [(r['trial'], r['channel'], r['serial'], r['x'], r['y'], r['z']) for r in rows]
    if len(rows) != 41 or len(set(sensor_keys)) != 41:
        raise ValueError('fixed41 sensor population is not complete/unique')
    data['temperature_summaries'] = {}
    for label, subset in [('paired41', rows)]+[(f'trial_{n}', [r for r in rows if r['trial'] == n]) for n in (10, 23)]:
        data['temperature_summaries'][label] = {basis: {kind: statistics(subset, f'{kind}_K', f'observed_{basis}_K')
            for kind in ('BASE', 'CONTROL', 'CANDIDATE')} for basis in ('minimum', 'p05', 'median')}
    if {relative(path): digest(path) for path in sorted(paths)} != hashes:
        raise ValueError('stage1 source/input changed during evaluation')
    data.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
        baseline_all_seven_reproduced=True, field_scoring_population_unchanged=True,
        interpretation='same-sensor diagnostic, not independent holdout or validated steady-mean sensor dynamics')
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(temperature=data['temperature_summaries'], saved=str(args.output)), indent=2), flush=True)


if __name__ == '__main__':
    main()
