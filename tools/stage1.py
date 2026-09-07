"""Portable entry point for the frozen Stage1 assessment and corrected BASE.

No installation, network, deletion, or overwrite. Use this project's Python.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT/'reference/preslhy'
ASSESSMENT = REF/'stage1_assessment_2026-09-06.json'
MANIFEST = REF/'stage1_delivery_manifest_2026-09-06.json'
REPORT = ROOT/'docs/stage1-results-2026-09-06.md'
GUIDE = ROOT/'docs/stage1-reproduction-2026-09-06.md'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_path(root, name):
    path = (root/name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'manifest dependency escapes project: {name}')
    return path


def verify_hashes(record, root=ROOT):
    failures = []
    if not record.get('completed') or not record.get('sha256'):
        raise ValueError('incomplete/empty evidence manifest')
    for name, expected in record['sha256'].items():
        path = checked_path(root, name)
        if not path.is_file() or digest(path) != expected:
            failures.append(name)
    if failures:
        raise ValueError('missing/changed files: '+', '.join(failures))
    return len(record['sha256'])


def verify():
    manifest = read(MANIFEST)
    count = verify_hashes(manifest)
    expected = read(ASSESSMENT)['runtime']['packages']
    differences = {}
    for name, version in expected.items():
        try:
            actual = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            actual = None
        if actual != version:
            differences[name] = dict(expected=version, actual=actual)
    if sys.version.split()[0] != manifest['python_version']:
        differences['python'] = dict(expected=manifest['python_version'], actual=sys.version.split()[0])
    print(json.dumps(dict(files_verified=count, runtime_matches=not differences,
        runtime_differences=differences, candidate_promoted=False,
        operational_path='BASE_SINGLE_IMAGE, non-touchdown only',
        historical_runner_hash_mismatch_retained=True), indent=2), flush=True)
    if differences:
        raise ValueError('runtime differs from snapshot; use the documented environment before reproduction')
    return manifest


def write_new(path, data):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def seal():
    if MANIFEST.exists():
        raise ValueError('refusing to replace sealed Stage1 bundle')
    assessment = read(ASSESSMENT)
    verify_hashes(assessment)
    paths = {checked_path(ROOT, p) for p in assessment['sha256']}
    paths.update([ASSESSMENT, REPORT, GUIDE, Path(__file__).resolve(),
        ROOT/'tests/test_stage1_tools.py', REF/'stage1_delivery_tests_2026-09-06.xml',
        REF/'stage1_delivery_smoke_2026-09-06/baseline_single_image.json'])
    paths.update((ROOT/'src').rglob('*.py'))
    # Tests and a fresh run must exist and pass before calling this command.
    from xml.etree import ElementTree
    xml = ElementTree.parse(REF/'stage1_delivery_tests_2026-09-06.xml')
    suites = list(xml.getroot().iter('testsuite'))
    if not suites or any(int(s.get('failures', 0))+int(s.get('errors', 0)) for s in suites):
        raise ValueError('delivery tests failed/missing')
    smoke = read(REF/'stage1_delivery_smoke_2026-09-06/baseline_single_image.json')
    if not smoke.get('completed') or smoke['maximum_saved_correction_reproduction_error'] > 1e-8:
        raise ValueError('fresh corrected-BASE smoke reproduction failed')
    manifest = dict(completed=True, stage='first_milestone_delivery',
        created_utc=datetime.now(timezone.utc).isoformat(), python_version=sys.version.split()[0],
        assessment=str(ASSESSMENT.relative_to(ROOT)).replace('\\', '/'),
        report=str(REPORT.relative_to(ROOT)).replace('\\', '/'),
        delivery_test_count=sum(int(s.get('tests', 0)) for s in suites),
        candidate_promoted=assessment['candidate_promoted'],
        sha256={str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in sorted(paths)})
    write_new(MANIFEST, manifest)
    print('sealed', len(manifest['sha256']), 'files:', MANIFEST, flush=True)


def fresh_baseline(directory):
    import numpy as np
    from audit_stage1_assessment import BASE_COLUMNS
    from audit_stage1_existing_evidence import concentration_metrics, geometry_metrics
    from audit_stage1_temperature import baseline_arguments
    from degali.validation.nearfield import hydrogen_jet, STEP, REACH
    from degali.validation.reported_jet_trajectory import ReportedJetTrajectory
    from run_preslhy_ambient_profile_audit import statistics
    assessment = read(ASSESSMENT)
    corrected = read(REF/'stage1_single_image_2026-09-06.json')
    temperatures = read(REF/'stage1_temperature_2026-09-06.json')
    trials = {t['trial']: t for t in read(REF/'e35_reduced.json')['trials']}
    runs, pairs, vertical, sensor_rows, errors = [], [], [], [], []
    for n in assessment['fixed_populations']['trials']:
        settings = baseline_arguments(trials[n])
        jp, initial = hydrogen_jet(**settings)
        rows = jp.run(initial, distmx=STEP, smax=REACH).rows
        model = ReportedJetTrajectory(jp.th.table, rows)
        runs.append(dict(trial=n, settings=settings, output_step_m=STEP, arclength_limit_m=REACH,
                         output_columns=BASE_COLUMNS, output_rows=rows.tolist()))
        for old in corrected['concentration_pairs']:
            if old['trial'] != n:
                continue
            sensors = [s for s in trials[n]['sensors'] if round(s['x'], 3) == old['x']]
            value = max(model.concentration_at(old['x'], s['y'], s['z']) for s in sensors)
            pairs.append(dict(old, predicted=value))
            errors.append(abs(value-old['predicted'])/max(1., abs(old['predicted'])))
        for old in corrected['vertical_profiles']:
            if old['trial'] != n:
                continue
            state = model.at(old['x'])
            if state is None:
                raise ValueError('missing vertical receptor')
            vertical.append(dict(old, modelled_centre=state.z, modelled_sigma_z=state.sz))
            errors.extend([abs(state.z-old['modelled_centre']), abs(state.sz-old['modelled_sigma_z'])])
        for old in corrected['temperature_rows']:
            if old['trial'] != n:
                continue
            value = model.temperature_at(old['x'], old['y'], old['z'])
            sensor_rows.append(dict(old, BASE_SINGLE_IMAGE_K=value))
            errors.append(abs(value-old['BASE_SINGLE_IMAGE_K']))
        print('fresh BASE_SINGLE_IMAGE trial', n, 'completed', flush=True)
    if (len(pairs), len(vertical), len(sensor_rows)) != (38, 17, 41) or max(errors) > 1e-8:
        raise ValueError('fresh corrected BASE differs from its fixed reference')
    verify_hashes(assessment)
    output = dict(completed=True, model='BASE_SINGLE_IMAGE', all_seven_freshly_integrated=True,
        created_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
        assessment_sha256=digest(ASSESSMENT), runner_sha256=digest(Path(__file__)),
        baseline_runs=runs, concentration_pairs=pairs, vertical_profiles=vertical,
        temperature_rows=sensor_rows, metrics=dict(**concentration_metrics(pairs), geometry=geometry_metrics(vertical)),
        temperature_summaries={basis: statistics(sensor_rows, 'BASE_SINGLE_IMAGE_K', f'observed_{basis}_K')
                              for basis in ('minimum', 'p05', 'median')},
        maximum_saved_correction_reproduction_error=float(max(errors)))
    write_new(directory/'baseline_single_image.json', output)
    print('saved', directory/'baseline_single_image.json', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('verify', 'baseline', 'candidate', 'seal'))
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--maximum-step', type=float, choices=(.01, .02), default=.02,
                        help='candidate only; baseline keeps the frozen .2 m output interval')
    args = parser.parse_args()
    if args.action == 'seal':
        seal()
        return
    if args.action == 'verify':
        verify()
        return
    if args.output_dir is None:
        parser.error('baseline/candidate requires a NEW --output-dir')
    if MANIFEST.exists():
        verify()
    else:
        # Used once for the verified delivery smoke test before sealing.
        verify_hashes(read(ASSESSMENT))
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    sys.path[:0] = [str(ROOT/'src'), str(ROOT/'tools')]
    if args.action == 'baseline':
        fresh_baseline(directory)
        return
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT/'src'), str(ROOT/'tools')]))
    for n in (10, 23):
        command = [sys.executable, str(ROOT/'tools/run_preslhy_source_ablation.py'),
            '--phase', 'field', '--mode', 'full', '--droplet-equilibrium-bound',
            '--consistent-phase-ambient', '--downstream-profile', 'enthalpy',
            '--interface-checkpoint', str(REF/'gaussian_enthalpy_profile_allflux_interface_2026-09-05.json'),
            '--trials', str(n), '--maximum-step', str(args.maximum_step),
            '--output', str(directory/f'candidate_{n}.json')]
        subprocess.run(command, cwd=ROOT, env=env, check=True)
    verify_hashes(read(ASSESSMENT))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, FileExistsError, FileNotFoundError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
