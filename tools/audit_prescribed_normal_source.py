"""Exact-geometry manufactured Q/normal-stress source cases, not LH2 inputs."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_finite_tke_manufactured import ROOT, manufactured_model, serial
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments
from degali.addons.prescribed_normal_source import retract_source_with_normal_stress


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite normal-source evidence')
    origin = ROOT/'reference/preslhy/prescribed_tke_source_2026-09-06.json'
    upstream = json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed'] or not upstream['numerical_passed']:
        raise ValueError('completed Q-only source verification required')
    for name, digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'frozen source dependency changed: {name}')
    paths = [ROOT/name for name in upstream['sha256']] + [origin, Path(__file__).resolve()]
    paths += [ROOT/name for name in ('src/degali/addons/prescribed_normal_source.py',
        'tests/test_prescribed_normal_source.py', 'docs/prereg-prescribed-normal-source.md')]
    def hashes():
        return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    data = dict(phase='manufactured_prescribed_normal_source', completed=False, rows=[], sha256=initial,
        started_utc=datetime.now(timezone.utc).isoformat(), physical_initialization_passed=False,
        adopted=False, field_scored=False, variance_ratio_was_fitted=False,
        pressure_assumption='ambient_pressure_no_compensation')
    def checkpoint():
        partial.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    checkpoint()
    model, geometry = manufactured_model()
    p, par = model.base.projection, model.base.parameters
    target = FaceSplitSquareMoments(p).moments(par, order=8, angular_order=8)
    old_target = p.target.copy()
    data.update(geometry=geometry, original_total_target=target)
    for center in (.02, 2.):
        q_parameters = np.r_[math.log(center), np.zeros(model.base.size)]
        row = dict(q_center=center, q_parameters=q_parameters, axial_variance_fraction=2/3, numerical_passed=False)
        data['rows'].append(row)
        def progress(item):
            row['progress'] = item
            checkpoint()
            print('normal source', center, item, flush=True)
        try:
            result = retract_source_with_normal_stress(p, par, q_parameters, total_moment_target=target,
                axial_variance_fraction=2/3, pressure_assumption=data['pressure_assumption'], callback=progress)
            row.update(result=result, numerical_passed=result['numerical_passed'],
                       maximum_shape_change=float(max(abs(result['parameters']-par))))
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row['failure'] = str(exc)
        checkpoint()
        print('normal source completed case', center, row['numerical_passed'], row.get('failure'), flush=True)
    if hashes() != initial or not np.array_equal(old_target, p.target):
        raise RuntimeError('frozen dependency or original source target changed')
    data.update(completed=True, original_projection_unchanged=True,
                finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('normal source audit completed', [(r['q_center'], r['numerical_passed']) for r in data['rows']], flush=True)


if __name__ == '__main__':
    main()
