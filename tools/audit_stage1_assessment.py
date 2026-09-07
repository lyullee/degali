"""Assemble Stage1 evidence; do not integrate, retune, or overwrite history."""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import numpy as np

from audit_stage1_existing_evidence import ROOT, REF, FIELD, read, digest, relative, serial, paired_rows

TRIALS = [10, 11, 12, 22, 23, 24, 25]
BASE_COLUMNS = ['x_m', 'centre_z_m', 'reported_imaged_C_H2_kg_m3',
    'sigma_y_m', 'sigma_z_m', 'theta_rad', 'excess_centre_velocity_m_s',
    'section_averaged_wind_m_s', 'reported_imaged_density_kg_m3',
    'reported_imaged_temperature_K', 'reported_imaged_H2_mole_fraction', 'arc_length_m']


def numerical_match(old, new, tolerance=1e-8):
    a, b = np.asarray(old, dtype=float), np.asarray(new, dtype=float)
    if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError('reproduction shape/nonfinite mismatch')
    error = float(np.max(abs(a-b)/np.maximum(1., abs(a))))
    if error > tolerance:
        raise ValueError(f'reproduction scaled error {error} exceeds {tolerance}')
    return dict(shape=list(a.shape), maximum_scaled_error=error,
                maximum_absolute_error=float(np.max(abs(a-b))))


def promotion_gates(base, candidate, interfaces_accepted):
    return dict(all_interfaces_accepted=bool(interfaces_accepted),
        concentration_bias_improved=abs(math.log(candidate['MG'])) < abs(math.log(base['MG'])),
        concentration_scatter_improved=candidate['VG'] < base['VG'],
        FAC2_not_worse=candidate['FAC2'] >= base['FAC2'],
        width_ratio_improved=abs(candidate['geometry']['width_ratio']-1) < abs(base['geometry']['width_ratio']-1),
        centre_MAE_not_worse=candidate['geometry']['centre_mae'] <= base['geometry']['centre_mae'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite assessment')
    evidence_path = REF/'stage1_existing_evidence_2026-09-06.json'
    temperature_path = REF/'stage1_temperature_2026-09-06.json'
    resolution_path = REF/'stage1_resolution_diagnostic_2026-09-06.json'
    correction_path = REF/'stage1_single_image_2026-09-06.json'
    evidence, temperature, resolution, correction = map(read,
        (evidence_path, temperature_path, resolution_path, correction_path))
    paths = {evidence_path, temperature_path, resolution_path, correction_path, Path(__file__).resolve()}
    for record in (evidence, temperature, resolution, correction):
        if not record['completed']:
            raise ValueError('required audit incomplete')
        for name, expected in record['sha256'].items():
            path = ROOT/name
            if digest(path) != expected:
                raise ValueError(f'audit dependency changed: {name}')
            paths.add(path)
    field = read(FIELD)
    trials = {r['trial']: r for r in read(REF/'e35_reduced.json')['trials']}
    sources = {r['trial']: r for r in read(REF/'measured_pipe_source_2026-09-05.json')['trials']}
    reproductions = []
    for number in (10, 23):
        path = REF/f'stage1_reproduction_{number}_2026-09-06.json'
        paths.add(path)
        fresh = read(path)
        if fresh['selected_trials'] != [number] or not fresh['interfaces_accepted'] or fresh['failures']:
            raise ValueError('representative reproduction failed')
        for name, expected in fresh['provenance_sha256'].items():
            if digest(ROOT/name) != expected:
                raise ValueError('fresh field provenance mismatch')
            paths.add(ROOT/name)
        checkpoint = REF/'gaussian_enthalpy_profile_allflux_interface_2026-09-05.json'
        if fresh['interface_checkpoint_sha256'] != digest(checkpoint):
            raise ValueError('fresh field checkpoint changed')
        checks = {}
        for variant in ('baseline', 'candidate'):
            for collection, observed, measured in (
                ('pairs', ('observed',), ('predicted', 'plume_centre')),
                ('vertical_profiles', ('measured_centre', 'measured_sigma_z'), ('modelled_centre', 'modelled_sigma_z'))):
                old = [r for r in field[variant][collection] if r['trial'] == number]
                new = fresh[variant][collection]
                keys = paired_rows(old, new, ('trial', 'x'), observed)
                a, b = ({(r['trial'], r['x']): r for r in group} for group in (old, new))
                checks[f'{variant}_{collection}'] = numerical_match(
                    [[a[k][v] for v in measured] for k in keys], [[b[k][v] for v in measured] for k in keys])
        old_interface, new_interface = field['interfaces'][str(number)], fresh['interfaces'][str(number)]
        for key in ('station_m', 'lambda', 'width_residual', 'temperature_residual_k', 'max_flux_residual', 'energy_quadrature_points'):
            checks[f'interface_{key}'] = numerical_match(old_interface[key], new_interface[key])
        if not new_interface['accepted'] or new_interface['failure_reasons']:
            raise ValueError('new interface failed existing acceptance conditions')
        old_state, new_state = old_interface['downstream'], new_interface['downstream']
        for key in old_state:
            checks[f'downstream_{key}'] = numerical_match(old_state[key], new_state[key])
        xs = np.asarray(new_state['states'])[:, 5]
        required_x = [r['x'] for r in field['candidate']['pairs'] if r['trial'] == number]
        required_x += [r['x'] for r in field['candidate']['vertical_profiles'] if r['trial'] == number]
        required_x += [r['x'] for r in temperature['temperature_rows'] if r['trial'] == number]
        if min(xs) > min(required_x) or max(xs) < max(required_x):
            raise ValueError('representative does not cover its fixed receptors')
        reproductions.append(dict(trial=number, checks=checks, state_count=len(xs),
            x_range_m=[float(min(xs)), float(max(xs))], required_x_range_m=[min(required_x), max(required_x)],
            maximum_relative_balance_residual=new_state['maximum_relative_balance_residual']))
    if [r['trial'] for r in temperature['baseline_runs']] != TRIALS:
        raise ValueError('fresh BASE population changed')
    if not temperature['baseline_all_seven_reproduced']:
        raise ValueError('BASE reproduction incomplete')
    if any(r['maximum_reproduction_error'] > 1e-8 for r in temperature['baseline_runs']):
        raise ValueError('BASE reproduction error')
    for run in temperature['baseline_runs']:
        if any(len(row) != len(BASE_COLUMNS) for row in run['output_rows']):
            raise ValueError('unexpected BASE row schema')
    metrics = dict(evidence['metrics'], BASE_SINGLE_IMAGE=correction['metrics'])
    gates = promotion_gates(metrics['BASE'], metrics['CANDIDATE'], field['interfaces_accepted'])
    corrected_gates = promotion_gates(metrics['BASE_SINGLE_IMAGE'], metrics['CANDIDATE'], field['interfaces_accepted'])
    if all(gates.values()) != field['promoted']:
        raise ValueError('recomputed original promotion decision differs')
    conditions = []
    for run in temperature['baseline_runs']:
        n, settings = run['trial'], run['settings']
        s, t = sources[n], trials[n]
        nozzle = s['topology'] == 'nozzle'
        conditions.append(dict(trial=n, topology=s['topology'], baseline_kwargs=settings,
            BASE_output_step_m=run['output_step_m'], BASE_arclength_limit_m=run['arclength_limit_m'],
            candidate_mass_flow_kg_s=s['pressure_loss_mass_flow_g_s']/1000.,
            candidate_to_BASE_flow_ratio=s['pressure_loss_mass_flow_g_s']/t['flow_mean_gs'],
            candidate_source_temperature_K=s['tc3_temperature_k']['median'] if nozzle else settings['storage_temperature'],
            candidate_source_upstream_pressure_barg=s['pt2_barg']['median'] if nozzle else None,
            candidate_temperature_basis='TC3_median' if nozzle else 'prior_tank_saturation_not_TC3',
            candidate_maximum_step_m=.02, candidate_relative_tolerance=1e-5))
    data = dict(completed=True, phase='stage1_final_assessment', created_utc=datetime.now(timezone.utc).isoformat(),
        version='0.1.0.dev0 with Stage1 hash snapshot, not a release certification',
        operational_path='BASE_SINGLE_IMAGE for positive-elevation non-touchdown JetPlume rows',
        historical_default_changed=False, candidate_promoted=all(gates.values()),
        original_promotion_gates=gates, candidate_vs_corrected_BASE_diagnostic_gates=corrected_gates,
        fixed_populations=dict(trials=TRIALS, concentration_arcs=38, vertical_profiles=17, temperature_receptors=41),
        metrics=metrics, temperatures=temperature['temperature_summaries'],
        corrected_BASE_temperatures=correction['temperature_summaries'],
        representative_reproductions=reproductions, baseline_all_seven_reproduced=True,
        baseline_maximum_reproduction_error=max(r['maximum_reproduction_error'] for r in temperature['baseline_runs']),
        executed_source_conditions=conditions, runtime=evidence['runtime'],
        reporting_schema_correction=dict(applies_to=relative(temperature_path),
            original_output_columns_named_only_first_five=True, complete_columns=BASE_COLUMNS,
            third_column_already_includes_ground_image=True, original_evidence_not_rewritten=True),
        historical_provenance_mismatches=evidence['code_or_input_mismatches'],
        resolution_diagnostic_file=relative(resolution_path), single_image_correction_file=relative(correction_path),
        limitations=['Same previously investigated campaign, not independent holdout.',
            'Only trials10/23 freshly reintegrated for candidate; other5 reuse verified evidence.',
            'BASE and candidate source reconstruction differ; not an isolated equation comparison.',
            'Two step levels on two cases do not establish convergence for all seven.',
            'Temperature extrema/p05/median diagnostics do not validate steady sensor dynamics.',
            'New enriched thermal/TKE/normal-stress modules are excluded from field adoption.',
            'No full Linux oracle, safety design certification, or general LH2 plant validation.'],
        sha256={relative(path): digest(path) for path in sorted(paths)})
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(completed=True, hashes=len(data['sha256']), gates=gates,
        reproduction_max_error=max(c['maximum_scaled_error'] for r in reproductions for c in r['checks'].values()),
        saved=str(args.output)), indent=2), flush=True)


if __name__ == '__main__':
    main()
