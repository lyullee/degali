"""Fixed-section cold-envelope diagnostic; no fitted or transported new state."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from stage2_matched_source import ROOT, REF, read, digest, write_new
from audit_stage2_flow_effect import replay_document
from run_preslhy_ambient_profile_audit import replay
from stage1 import verify


def reflected_peak(centre, sigma):
    """Exact nonnegative mode and amplitude of equal reflected Gaussians."""
    if not math.isfinite(centre) or centre < 0 or not math.isfinite(sigma) or sigma <= 0:
        raise ValueError('finite nonnegative centre and positive sigma required')
    a = centre / sigma
    if a <= 1:
        z = 0.
    else:
        aa = a*a
        # Near the bifurcation the symmetric mode differs below float accuracy.
        if aa - 1 < 1e-12:
            z = 0.
        else:
            u = brentq(lambda u: math.tanh(aa*u)-u, 1e-12, 1., xtol=5e-15)
            z = centre*u
    q = math.exp(-.5*((z-centre)/sigma)**2) + math.exp(-.5*((z+centre)/sigma)**2)
    return z, q


def scalar_shape(state, sy, sz, y, z):
    return math.exp(-.5*(y/sy)**2) * (math.exp(-.5*((z-state[6])/sz)**2)
                                            + math.exp(-.5*((z+state[6])/sz)**2))


def envelope(model, state):
    sy, sz = model.section_widths(state)
    peak_z, qmax = reflected_peak(float(state[6]), sz)
    th = model.thermodynamics
    def temperatures(q):
        q = np.atleast_1d(q)
        rho = model.rhoa + (state[0]-model.rhoa)*q
        mass = state[0]*state[1]*q
        if np.any(rho <= 0) or np.any(mass < 0) or np.any(mass > rho):
            raise ValueError('fixed reflected section includes nonphysical composition')
        t, _ = th._condensed_air_state(rho, mass/rho)
        if not np.all(np.isfinite(t)):
            raise ValueError('nonfinite envelope temperature')
        return t
    coarse = temperatures(np.linspace(0, qmax, 257))
    grid = np.linspace(0, qmax, 1025)
    fine = temperatures(grid)
    candidates = [(float(fine[0]), 0.), (float(fine[-1]), qmax)]
    for i in range(1, len(grid)-1):
        if fine[i] <= fine[i-1] and fine[i] < fine[i+1]:
            r = minimize_scalar(lambda q: float(temperatures(q)[0]), bounds=(grid[i-1],grid[i+1]),
                                method='bounded', options={'xatol': 1e-12})
            if not r.success: raise ValueError('local cold-envelope minimization failed')
            candidates.append((float(r.fun), float(r.x)))
    tmin, qmin = min(candidates)
    centre_q = scalar_shape(state, sy, sz, 0., float(state[6]))
    return dict(centre_z_m=float(state[6]), sigma_y_m=sy, sigma_z_m=sz,
                peak_shape_z_m=peak_z, maximum_reflected_shape=qmax,
                minimum_temperature_K=tmin, minimum_temperature_shape=qmin,
                sample_monotone_cooling=bool(np.all(np.diff(fine) <= 1e-7)),
                grid_refinement_K=abs(float(np.min(coarse))-float(np.min(fine))),
                local_refinement_K=abs(float(np.min(fine))-tmin),
                internal_centre_temperature_K=model.centre_temperature(state),
                reflected_centre_temperature_K=float(temperatures(centre_q)[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    verify()
    paths = [REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
             REF/'stage2_coriolis_2026-09-06/field_density_10-23_step0.02.json',
             REF/'stage1_single_image_2026-09-06.json', REF/'e35_reduced.json',
             Path(__file__).resolve(), ROOT/'docs/prereg-downstream-cold-envelope.md',
             ROOT/'tools/audit_stage2_flow_effect.py', ROOT/'tools/run_preslhy_ambient_profile_audit.py']
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    fields = dict(pressure_loss_density=read(paths[0]), coriolis_density=replay_document(read(paths[1])))
    frozen, reduced = read(paths[2]), read(paths[3])
    trials = {t['trial']:t for t in reduced['trials']}
    checks, sections, rows, failures = {}, [], [], []
    for n in (10,23):
        for variant, field in fields.items():
            trajectory, checks[f'{variant}_{n}'] = replay(field,trials[n])
            local = [r for r in frozen['temperature_rows'] if r['trial']==n]
            for x in sorted({r['x'] for r in local}):
                state = trajectory.state_at(x)
                if state is None: raise ValueError('fixed sensor not covered')
                try:
                    section = envelope(trajectory.model,state)
                except ValueError as err:
                    failures.append(dict(trial=n,variant=variant,x=x,error=str(err)))
                    continue
                sections.append(dict(trial=n,variant=variant,x=x,**section))
                for old in [r for r in local if r['x']==x]:
                    prediction = trajectory.temperature_at(x,old['y'],old['z'])
                    if variant=='pressure_loss_density':
                        stored = old['CONTROL_K']
                    else:
                        stored = next(r['matched_K'] for r in read(paths[1])['fixed_scores']['temperature_rows']
                                      if r['trial']==n and r['channel']==old['channel'])
                    if abs(prediction-stored)>1e-6: raise ValueError('fixed sensor replay changed')
                    q = scalar_shape(state,section['sigma_y_m'],section['sigma_z_m'],old['y'],old['z'])
                    if q>section['maximum_reflected_shape']+1e-12 or prediction<section['minimum_temperature_K']-1e-6:
                        raise ValueError('sensor lies below computed fixed-section cold envelope')
                    row = dict(trial=n,variant=variant,channel=old['channel'],x=x,y=old['y'],z=old['z'],
                               release_height_m=trials[n]['release_height_m'],sensor_temperature_K=prediction,
                               fixed_section_minimum_K=section['minimum_temperature_K'])
                    for basis in ('minimum','p05','median'):
                        observation = old[f'observed_{basis}_K']
                        row[f'observed_{basis}_K']=observation
                        row[f'sensor_minus_{basis}_K']=prediction-observation
                        row[f'cold_envelope_gap_{basis}_K']=max(section['minimum_temperature_K']-observation,0.)
                    rows.append(row)
    summary = {}
    for variant in fields:
        summary[variant]={}
        for n in (10,23):
            for group in ('all','release_axis','off_axis'):
                selected = [r for r in rows if r['variant']==variant and r['trial']==n and
                            (group=='all' or ((abs(r['y'])<1e-9 and abs(r['z']-r['release_height_m'])<1e-9)==(group=='release_axis')))]
                summary[variant][f'{n}_{group}']={basis:dict(n=len(selected),
                    sensor_mae_K=float(np.mean([abs(r[f'sensor_minus_{basis}_K']) for r in selected])),
                    cold_envelope_gap_mean_K=float(np.mean([r[f'cold_envelope_gap_{basis}_K'] for r in selected])),
                    cold_envelope_gap_max_K=max(r[f'cold_envelope_gap_{basis}_K'] for r in selected),
                    n_gap_over_1K=sum(r[f'cold_envelope_gap_{basis}_K']>1 for r in selected))
                    for basis in ('minimum','p05','median')}
    if {str(p.relative_to(ROOT)):digest(p) for p in paths}!=hashes: raise ValueError('input changed')
    verify()
    out=dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
             replay_checks=checks,sections=sections,rows=rows,summary=summary,failures=failures,
             expected_rows=82,covered_rows=len(rows),coverage_complete=len(rows)==82 and not failures,
             reintegrated=False,new_model=False,candidate_promoted=False,
             limitation='Fixed-section spatial diagnostic, not a formal global EOS bound or causal attribution; extrema are not steady means.')
    write_new(args.output,out)
    print(summary,flush=True)
    print('Saved',args.output,'coverage',out['coverage_complete'],flush=True)


if __name__=='__main__': main()
