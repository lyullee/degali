"""Frozen manufactured source-energy reallocation, no trial Q selected."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_finite_tke_manufactured import ROOT, manufactured_model, serial
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments
from degali.addons.prescribed_tke_source import retract_source_with_prescribed_tke


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite source reallocation evidence')
    origin = ROOT/'reference/preslhy/finite_tke_manufactured_2026-09-06.json'
    upstream = json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or not upstream['numerical_passed']:
        raise ValueError('completed finite-TKE numerical operator required first')
    for name, digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'upstream frozen dependency changed: {name}')
    paths = [ROOT/name for name in upstream['sha256']] + [origin, Path(__file__).resolve()]
    paths += [ROOT/name for name in ('src/degali/addons/prescribed_tke_source.py',
        'tests/test_prescribed_tke_source.py', 'docs/prereg-prescribed-tke-source-retraction.md')]
    def hashes():
        return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    data = dict(phase='manufactured_prescribed_tke_source_reallocation', completed=False,
        numerical_passed=False, physical_initialization_passed=False, field_scored=False,
        adopted=False, q_was_fitted=False, sha256=initial, started_utc=datetime.now(timezone.utc).isoformat())
    def checkpoint():
        partial.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    def progress(row):
        data['progress'] = row
        checkpoint()
        print('source reallocation', row, flush=True)
    checkpoint()
    try:
        model, geometry = manufactured_model()
        p, par = model.base.projection, model.base.parameters
        target = FaceSplitSquareMoments(p).moments(par, order=8, angular_order=8)
        original_target = p.target.copy()
        data.update(geometry=geometry, original_mean_thermal_target=target,
                    explicit_q_parameters=model.parameters.copy())
        result = retract_source_with_prescribed_tke(p, par, model.parameters,
            total_moment_target=target, callback=progress)
        data['result'] = result
        data['original_projection_unchanged'] = bool(np.array_equal(original_target, p.target))
        data['parameters_changed_maximum'] = float(max(abs(result['parameters']-par)))
        data['mean_thermal_energy_change'] = float(result['actual_mean_thermal_moments'][3]-target[3])
        data['total_source_energy_change'] = float(result['actual_total_moments'][3]-target[3])
        data['numerical_passed'] = bool(result['numerical_passed'] and data['original_projection_unchanged'])
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        data['failure'] = str(exc)
    if hashes() != initial:
        raise RuntimeError('frozen source reallocation dependencies changed')
    data.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('source reallocation completed', data['numerical_passed'], data.get('failure'), flush=True)


if __name__ == '__main__':
    main()
