"""Frozen six-DOF Q/normal source audit; retain both manufactured cases."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_finite_tke_manufactured import ROOT, manufactured_model, serial
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments
from degali.addons.mobile_tke_normal_source import MobileTkeNormalSource


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite mobile-source evidence')
    origin = ROOT/'reference/preslhy/prescribed_normal_source_2026-09-06.json'
    upstream = json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed fixed-center normal-source audit required, including failures')
    for name, digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'frozen normal-source dependency changed: {name}')
    paths = [ROOT/name for name in upstream['sha256']] + [origin, Path(__file__).resolve()]
    paths += [ROOT/name for name in ('src/degali/addons/mobile_tke_normal_source.py',
        'tests/test_mobile_tke_normal_source.py', 'docs/prereg-mobile-tke-normal-source.md')]
    def hashes():
        return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    data = dict(phase='manufactured_mobile_tke_normal_source', completed=False, rows=[], sha256=initial,
        started_utc=datetime.now(timezone.utc).isoformat(), physical_initialization_passed=False,
        adopted=False, field_scored=False, q_was_fitted=False, variance_ratio_was_fitted=False,
        pressure_assumption='ambient_pressure_no_compensation')
    def checkpoint():
        partial.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    checkpoint()
    model, geometry = manufactured_model()
    p, par = model.base.projection, model.base.parameters
    target = FaceSplitSquareMoments(p).moments(par, order=8, angular_order=8)
    original_state, old_target, old_parameters = p.mixing.state.copy(), p.target.copy(), par.copy()
    data.update(original_geometry=geometry, original_total_target=target, original_state=original_state)
    for center in (.02, 2.):
        q_parameters = np.r_[math.log(center), np.zeros(model.base.size)]
        row = dict(q_center=center, q_parameters=q_parameters, axial_variance_fraction=2/3, numerical_passed=False)
        data['rows'].append(row)
        source = MobileTkeNormalSource(model.base, q_parameters, axial_variance_fraction=2/3,
                                      pressure_assumption=data['pressure_assumption'])
        def progress(item):
            row['progress'] = item
            checkpoint()
            print('mobile source', center, item, flush=True)
        try:
            result = source.solve(total_moment_target=target, callback=progress)
            row.update(result=result, numerical_passed=result['numerical_passed'],
                center_area_velocity_ratios=np.exp(result['changes'][:4]),
                maximum_shape_change=float(max(abs(result['parameters']-par))))
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row['failure'] = str(exc)
        checkpoint()
        print('mobile source completed case', center, row['numerical_passed'], row.get('failure'), flush=True)
    if (hashes() != initial or not np.array_equal(old_target, p.target)
            or not np.array_equal(original_state, p.mixing.state) or not np.array_equal(old_parameters, par)):
        raise RuntimeError('frozen dependency or original source changed')
    data.update(completed=True, original_projection_unchanged=True,
                finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('mobile source audit completed', [(r['q_center'], r['numerical_passed']) for r in data['rows']], flush=True)


if __name__ == '__main__':
    main()
