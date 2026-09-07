"""Stage1: verify retained evidence and recompute metrics, without model edits."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT/'reference/preslhy'
FIELD = REF/'gaussian_enthalpy_profile_allflux_field_complete_2026-09-05.json'
CONTROL = REF/'phase_ambient_consistency_field_complete_2026-09-05.json'
TEMPERATURE = REF/'gaussian_enthalpy_profile_temperature_audit_2026-09-05.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT)).replace('\\', '/')


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def concentration_metrics(rows):
    observed = np.array([r['observed'] for r in rows])
    predicted = np.array([r['predicted'] for r in rows])
    if not len(rows) or np.any(observed <= 0.) or np.any(predicted <= 0.):
        raise ValueError('positive paired concentration values required')
    error = np.log(observed/predicted)
    return dict(arcs=len(rows), MG=float(np.exp(np.mean(error))),
        VG=float(np.exp(np.mean(error**2))), FAC2=float(np.mean(abs(error) <= np.log(2.))))


def geometry_metrics(rows):
    center = np.array([r['modelled_centre']-r['measured_centre'] for r in rows])
    return dict(width_ratio=float(np.mean([r['modelled_sigma_z']/r['measured_sigma_z'] for r in rows])),
        centre_mean_error=float(np.mean(center)), centre_mae=float(np.mean(abs(center))))


def paired_rows(first, second, keys, observed_keys):
    a = {tuple(r[k] for k in keys): r for r in first}
    b = {tuple(r[k] for k in keys): r for r in second}
    if len(a) != len(first) or len(b) != len(second) or a.keys() != b.keys():
        raise ValueError('duplicate or unmatched observation keys')
    for key in a:
        for name in observed_keys:
            if abs(a[key][name]-b[key][name]) > 1e-12*max(abs(a[key][name]), 1.):
                raise ValueError(f'observations differ at {key}: {name}')
    return sorted(a)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite stage1 evidence')
    paths, checks, seen = set(), [], set()
    def trace(path):
        path = path.resolve()
        if path in seen:
            return
        seen.add(path); paths.add(path)
        data = read(path)
        for name, expected in data.get('aggregation', {}).get('source_sha256', {}).items():
            child = path.parent/name
            if not child.resolve().is_relative_to(ROOT):
                raise ValueError('aggregation escapes project')
            actual = digest(child)
            checks.append(dict(record=relative(path), dependency=relative(child), kind='aggregation',
                               matches=actual == expected, expected=expected, actual=actual))
            if actual != expected:
                raise ValueError('frozen field aggregation changed')
            trace(child)
        for name, expected in data.get('provenance_sha256', data.get('sha256', {})).items():
            child = Path(name)
            if not child.is_absolute():
                child = ROOT/child
            child = child.resolve()
            if not child.is_relative_to(ROOT):
                raise ValueError('recorded source dependency escapes project')
            paths.add(child)
            actual = digest(child) if child.is_file() else None
            checks.append(dict(record=relative(path), dependency=relative(child), kind='recorded_input_or_code',
                               matches=actual == expected, expected=expected, actual=actual))
    for path in (FIELD, CONTROL, TEMPERATURE,
                 REF/'gaussian_enthalpy_profile_allflux_interface_2026-09-05.json'):
        trace(path)
    field, control, temperature = read(FIELD), read(CONTROL), read(TEMPERATURE)
    keys = paired_rows(field['baseline']['pairs'], field['candidate']['pairs'],
                       ('trial', 'x'), ('observed',))
    vertical_keys = paired_rows(field['baseline']['vertical_profiles'], field['candidate']['vertical_profiles'],
                                ('trial', 'x'), ('measured_centre', 'measured_sigma_z'))
    paired_rows(field['candidate']['pairs'], control['candidate']['pairs'], ('trial', 'x'), ('observed',))
    paired_rows(field['candidate']['vertical_profiles'], control['candidate']['vertical_profiles'],
                ('trial', 'x'), ('measured_centre', 'measured_sigma_z'))
    if len(keys) != 38 or len(vertical_keys) != 17 or field['selected_trials'] != [10, 11, 12, 22, 23, 24, 25]:
        raise ValueError('fixed38/17/7 population changed')
    metrics = {}
    for label, item in (('BASE', field['baseline']), ('CONTROL', control['candidate']),
                         ('CANDIDATE', field['candidate'])):
        current = dict(**concentration_metrics(item['pairs']), geometry=geometry_metrics(item['vertical_profiles']))
        errors = [abs(current[k]-item[k])/max(1., abs(item[k])) for k in ('MG', 'VG', 'FAC2')]
        errors += [abs(current['geometry'][k]-v)/max(1., abs(v)) for k, v in item['geometry'].items()]
        if max(errors) > 1e-10:
            raise ValueError(f'{label} saved aggregate does not match actual paired rows')
        current['maximum_aggregate_recomputation_error'] = max(errors)
        current['by_trial'] = {str(t): concentration_metrics([r for r in item['pairs'] if r['trial'] == t])
                               for t in field['selected_trials']}
        metrics[label] = current
    raw_checks = []
    for number, entry in temperature['raw_workbooks'].items():
        path = REF/'raw'/entry['file']
        paths.add(path)
        actual = digest(path)
        raw_checks.append(dict(trial=int(number), file=relative(path), expected=entry['sha256'], actual=actual,
                               matches=actual == entry['sha256'], window=entry['data_window_zero_based_half_open']))
        if actual != entry['sha256']:
            raise ValueError('raw temperature workbook no longer matches fixed reduction')
    if len(temperature['temperature_rows']) != 41:
        raise ValueError('fixed41 temperature population changed')
    # Freeze a complete current source snapshot, not merely the historical partial manifest.
    paths.update((ROOT/'src').rglob('*.py'))
    paths.update([Path(__file__).resolve(), ROOT/'pyproject.toml',
                  ROOT/'docs/prereg-first-milestone-reproduction.md',
                  ROOT/'docs/first-milestone-scope-2026-09-06.md'])
    data = dict(completed=True, phase='stage1_existing_evidence',
        created_utc=datetime.now(timezone.utc).isoformat(), code_or_input_mismatches=[r for r in checks if not r['matches']],
        provenance_checks=checks, raw_workbook_checks=raw_checks, metrics=metrics,
        concentration_keys=keys, vertical_keys=vertical_keys,
        temperature_control_is_BASE=False, temperature_BASE_requires_new_evaluation=True,
        historical_reproduction_not_new_holdout=True, latest_candidate_promoted=False,
        runtime=dict(python=sys.version, platform=platform.platform(),
            packages={name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'CoolProp', 'pytest')}),
        sha256={relative(path): digest(path) for path in sorted(paths)},
        limitations=['Historical recorded manifests cover only their listed dependencies.',
                     'Historical runner/code hash differences are retained, not silently accepted as identical.',
                     'Current snapshot plus representative reintegration is not seven fresh candidate ODE runs.'])
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(metrics=metrics, mismatch_count=len(data['code_or_input_mismatches']),
                          hashes=len(data['sha256']), saved=str(args.output)), indent=2), flush=True)


if __name__ == '__main__':
    main()
