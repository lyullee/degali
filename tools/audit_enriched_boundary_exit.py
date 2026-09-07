"""Post-screen independent reconstruction of rejected boundary endpoints."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK, independent_rays
from degali.addons.enriched_transport import EnrichedModalTransport, PrescribedRadialMixing
from degali.addons.enriched_segments import EnrichedShortSegment


def witnesses(model, rates):
    """Independent radial integration; keep physical residual components."""
    at = np.r_[1., rates]
    p, m = model.projection, model.projection.mixing
    rows = []
    angles = np.linspace(0., math.pi/4, 65)
    for angle, cuts in zip(angles, model.quadrature.partitions(model.parameters, angles)):
        extra = model.mixing.knots[(model.mixing.knots > 0.) & (model.mixing.knots < cuts[-1])]
        cuts = np.sort(np.r_[cuts, extra])
        cuts = cuts[np.r_[True, np.diff(cuts) > cuts[-1]*2e-13]]
        _, _, _, edge, _ = model._ray(angle, cuts, 16)
        t = math.tan(angle)
        d = model.local(np.array([1.]), np.array([t]))
        fm, fp = edge@at
        u, rho, y, h = [d[k][0] for k in ('u', 'rho', 'y', 'specific_h')]
        q = np.array([p.q0*(1+t*t)])
        chi, oldfm = model.mixing.evaluate(q)[0]
        old = m.local(q)
        scale = np.array([max(abs(old['y'][0]*oldfm), 1e-12),
            max(abs(old['h'][0]/old['rho'][0]*oldfm), 1.), max(abs(old['u'][0]*oldfm), 1.)])
        adv = np.array([y*fm, h*fm])
        diff = -model.area*rho*chi*np.array([d['grad_y'][0, 0], d['grad_h'][0, 0]])/(2*p.q0)
        heat_external = .5*m.wind**2*fm-(u*fp-.5*u*u*fm)
        residual = np.r_[adv+diff-np.array([0., heat_external]), fp-m.wind*math.cos(m.state[3])*fm]
        rows.append(dict(angle_rad=angle, a=1., b=t, physical_y_m=m.sy,
            physical_n_m=m.sn*t, mass=fm, momentum=fp, scalar_advection=adv,
            scalar_diffusion=diff, heat_external=heat_external, residuals=residual,
            scales=scale, normalized=residual/scale,
            chi_momentum=-(fp-u*fm)/(model.area*rho*d['uq'][0])))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite boundary-exit evidence')
    ref = ROOT/'reference/preslhy'
    origin = ref/'enriched_segments_a_2026-09-06.json'
    data = json.loads(origin.read_text(encoding='utf-8'))
    if not data['completed']:
        raise ValueError('completed original screen required')
    for name, digest in data['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths = [ROOT/n for n in data['sha256']]+[origin, Path(__file__).resolve(),
        ROOT/'docs/prereg-enriched-boundary-exit-diagnosis.md']
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    frozen_hashes = hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    source = next(r for r in read('edge_conservative_refit_combined_2026-09-05.json')['rows'] if r['trial']==10)
    initial = next(r for r in read('coupled_initialization_combined_2026-09-06.json')['rows'] if r['trial']==10)
    original = next(r for r in data['rows'] if r['trial']==10)
    measured = {r['trial']:r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    p, replay = construct_projection(source, read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json'),
        read('e35_reduced.json')['trials'], measured, {})
    mixing = PrescribedRadialMixing.from_weak_baseline(p, order=32)
    kwargs = dict(scalar_mixing=mixing, thermal_species_ratio=1., mechanical_work=WORK)
    driver = EnrichedShortSegment(p, np.array(initial['parameters']), mixing_update='equilibrium_velocity_width', **kwargs)
    results = dict(phase='postscreen_trial10_boundary_exit', completed=False, trial=10, rows=[],
        sha256=frozen_hashes, source_reconstruction=replay, original_passed=False,
        field_scored=False, promoted=False, started_utc=datetime.now(timezone.utc).isoformat())
    def checkpoint():
        partial.write_text(json.dumps(results, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    checkpoint()
    for segment in original['segments']:
        if segment['divisions'] not in (4, 8):
            continue
        if segment['completed'] or segment['failure_stage'] not in ('midpoint','endpoint'):
            raise ValueError('expected a saved rejection after accepted state')
        value = np.array(segment['parameters'])
        ds = original['length_m']/segment['divisions']
        first = driver.evaluate(value)
        midpoint = value+.5*ds*first['rates']
        middle = driver.evaluate(midpoint, enforce=False)
        at_midpoint = segment['failure_stage'] == 'midpoint'
        endpoint = midpoint if at_midpoint else value+ds*middle['rates']
        row = dict(divisions=segment['divisions'], failure_stage=segment['failure_stage'],
            start_m=segment['reached_m'], endpoint_m=segment['reached_m']+ds*(.5 if at_midpoint else 1.),
            start=value, midpoint=midpoint, endpoint=endpoint, start_edges=first['edge_defects'], midpoint_edges=middle['edge_defects'])
        results['rows'].append(row)
        checkpoint()
        evaluations = [driver.evaluate(endpoint, order=o, angular_order=a, exact=exact, enforce=False)
            for o,a,exact in ((8,48,False),(4,8,True),(8,16,True))]
        row['evaluations'] = evaluations
        fine = evaluations[-1]
        row['rate_refinement'] = max(abs(fine['rates']-evaluations[-2]['rates'])/np.maximum(abs(fine['rates']),1.))
        current, par, supplied, _ = driver.context(endpoint)
        model = EnrichedModalTransport(current, par, scalar_mixing=supplied, thermal_species_ratio=1., mechanical_work=WORK)
        row['independent_rays'] = independent_rays(model, fine['rates'])
        row['witnesses'] = witnesses(model, fine['rates'])
        row['worst_heat_witness'] = max(row['witnesses'], key=lambda r:abs(r['normalized'][1]))
        print('trial10 reconstructed',segment['divisions'],'endpoint_m',row['endpoint_m'],
            'fine_edges',fine['edge_defects'],'rays',row['independent_rays']['edge_defects'],flush=True)
        checkpoint()
    # Compare at the same distance, but a rejected explicit internal stage is
    # lower order than an endpoint. This is not a standard Richardson test.
    a,b = results['rows']
    results['same_endpoint_distance'] = bool(abs(a['endpoint_m']-b['endpoint_m'])<1e-14)
    results['comparison_is_richardson_endpoint_refinement'] = False
    if results['same_endpoint_distance']:
        results['relative_endpoint_change'] = max(abs(a['endpoint']-b['endpoint'])/np.maximum(abs(b['endpoint']),1.))
    alt = EnrichedShortSegment(p, np.array(initial['parameters']), mixing_update='constant_geometric_diffusivity', **kwargs)
    results['fixed_state_alternative'] = alt.evaluate(b['endpoint'], order=8, angular_order=16, exact=True, enforce=False)
    if hashes() != frozen_hashes:
        raise RuntimeError('boundary-exit dependencies changed during the diagnostic')
    results['completed'] = True
    args.output.write_text(json.dumps(results, default=serial, indent=2, allow_nan=False)+'\n',encoding='utf-8')
    print('completed independent boundary-exit diagnosis',flush=True)


if __name__ == '__main__':
    main()
