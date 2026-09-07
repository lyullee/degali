"""Aggregate fixed populations and localize the matched-flow effect."""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import math
from pathlib import Path

import numpy as np

from stage2_matched_source import ROOT, REF, TRIALS, read, digest, write_new, verify_inputs
from audit_stage1_existing_evidence import concentration_metrics, geometry_metrics, paired_rows
from run_preslhy_ambient_profile_audit import replay, statistics


def replay_document(run):
    if (not run['completed'] or run['profile'] != 'density'
            or run['flow_choice'] != 'coriolis' or run['fit_velocity_spreading']):
        raise ValueError('unexpected executed model configuration')
    field = deepcopy(run['field'])
    # These values come from the executed settings in stage2_matched_source,
    # not from guessing a thermodynamic convention based on stored scores.
    field.update(mode='full', droplet_equilibrium_bound=True,
        hydrogen_spin_isomer='normal', source_energy_ledger='moving_pipe_phase_enthalpy_v2',
        phase_ambient_closure='consistent_explicit_ideal_v1', downstream_thermodynamic_profile='density')
    return field


def metric_bundle(pairs, vertical):
    return dict(**concentration_metrics(pairs), geometry=geometry_metrics(vertical))


def aggregate(runs, frozen, control):
    selected = [n for r in runs for n in r['fixed_scores']['selected_trials']]
    if len(selected) != 7 or set(selected) != set(TRIALS):
        raise ValueError('field groups do not cover the seven unique trials')
    groups = {}
    for key in ('concentration_pairs', 'vertical_profiles', 'temperature_rows'):
        groups[key] = [row for r in runs for row in r['fixed_scores'][key]]
    if any(not r['completed'] or not r['fixed_scores']['coverage_complete'] for r in runs):
        raise ValueError('a field run is incomplete or has missing observations')
    pairs, vertical, temperatures = (groups[k] for k in groups)
    paired_rows(pairs, frozen['concentration_pairs'], ('trial', 'x'), ('observed',))
    paired_rows(vertical, frozen['vertical_profiles'], ('trial', 'x'),
                ('measured_centre', 'measured_sigma_z'))
    paired_rows(temperatures, frozen['temperature_rows'], ('trial', 'channel', 'x', 'y', 'z'),
                ('observed_minimum_K', 'observed_p05_K', 'observed_median_K'))
    paired_rows(pairs, control['candidate']['pairs'], ('trial', 'x'), ('observed',))
    paired_rows(vertical, control['candidate']['vertical_profiles'], ('trial', 'x'),
                ('measured_centre', 'measured_sigma_z'))
    if tuple(map(len, (pairs, vertical, temperatures))) != (38, 17, 41):
        raise ValueError('wrong fixed population cardinalities')
    metrics = {
        'pressure_loss_density': metric_bundle(control['candidate']['pairs'], control['candidate']['vertical_profiles']),
        'coriolis_density': metric_bundle(pairs, vertical),
        'BASE_SINGLE_IMAGE': metric_bundle(frozen['concentration_pairs'], frozen['vertical_profiles']),
    }
    per_trial = {}
    for n in TRIALS:
        per_trial[str(n)] = {}
        for label, ps, vs in (
            ('pressure_loss_density', control['candidate']['pairs'], control['candidate']['vertical_profiles']),
            ('coriolis_density', pairs, vertical)):
            pp, vv = [r for r in ps if r['trial'] == n], [r for r in vs if r['trial'] == n]
            per_trial[str(n)][label] = metric_bundle(pp, vv)
    temps = {}
    for label, subset in [('all41', temperatures)] + [
            (f'trial_{n}', [r for r in temperatures if r['trial'] == n]) for n in (10, 23)]:
        temps[label] = {basis: {
            kind: statistics(subset, key, f'observed_{basis}_K')
            for kind, key in [('pressure_loss_density', 'CONTROL_K'), ('coriolis_density', 'matched_K'),
                              ('BASE_SINGLE_IMAGE', 'BASE_SINGLE_IMAGE_K')]}
            for basis in ('minimum', 'p05', 'median')}
    return dict(groups=groups, metrics=metrics, per_trial=per_trial, temperature=temps)


def section_record(trajectory, x):
    model = trajectory.model
    state = trajectory.state_at(x)
    if state is None:
        raise ValueError('requested diagnostic station is not covered')
    sy, sz, ua, _, _, _, _, sya, sza, dsya, dsza = model._geometry(state)
    sources = model.source_terms(state)
    vel, rho, _, _, h = model.profiles(state)
    _, w = model._quadrature(model.quadrature_points)
    ke = float(np.sum(.5 * rho * vel**3 * w) * state[2])
    thermal = float(np.sum(h * vel * w) * state[2])
    # A diagnostic of the existing x>1m derivative switch. It changes no RHS.
    jp = model.jetplume
    ct = math.cos(state[3])
    exact_da = sya / x * jp.betay * ct
    exact_dz = sza / x * (jp.betaz + 2 * jp.gammaz * math.log(x)) * ct
    factor = model.rhoa * model.k.profile_integral(0.) * ua
    supplied_spread_mass = factor * (sza * dsya + sya * dsza)
    differentiated_spread_mass = factor * (sza * exact_da + sya * exact_dz)
    return dict(x_m=float(x), centre_z_m=float(state[6]), centre_T_K=model.centre_temperature(state),
        centre_density_kg_m3=float(state[0]), centre_H2_mass_fraction=float(state[1]),
        sigma_y_m=sy, sigma_z_m=sz, angle_rad=float(state[3]),
        centre_axial_velocity_m_s=float(state[4] + ua * ct),
        ambient_axial_velocity_m_s=ua * ct,
        fluxes=model._as_array(model.integral_fluxes(state)),
        thermal_flux_W=thermal, mean_kinetic_flux_W=ke,
        buoyancy_N_m=model.buoyancy_force(state), sources=sources,
        ambient_spread_mass_implemented_kg_s_m=supplied_spread_mass,
        ambient_spread_mass_from_width_derivative_kg_s_m=differentiated_spread_mass,
        missing_spread_over_total_mass_source=(differentiated_spread_mass-supplied_spread_mass) / sources[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    directory = args.directory.resolve()
    if args.output.exists():
        raise FileExistsError(args.output)
    verify_inputs(directory)
    paths = [directory / name for name in ('field_density_10-23_step0.02.json',
        'field_density_11-12-22-24-25_step0.02.json')]
    paths += [REF / 'stage1_single_image_2026-09-06.json',
        REF / 'phase_ambient_consistency_field_complete_2026-09-05.json',
        REF / 'e35_reduced.json', Path(__file__).resolve()]
    hashes = {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in paths}
    runs = [read(p) for p in paths[:2]]
    frozen, control = read(paths[2]), read(paths[3])
    data = aggregate(runs, frozen, control)
    trials = {t['trial']: t for t in read(paths[4])['trials']}
    replay_checks, sections, sensor_geometry = {}, [], []
    for n in TRIALS:
        run = next(r for r in runs if n in r['fixed_scores']['selected_trials'])
        new_field = replay_document(run)
        models = {}
        for label, field in [('pressure_loss_density', control), ('coriolis_density', new_field)]:
            if not math.isclose(field['interfaces'][str(n)]['lambda'], 1.16, abs_tol=1e-13):
                raise ValueError('velocity/scalar width changes would confound the flow-only comparison')
            models[label], replay_checks[f'{label}_{n}'] = replay(field, trials[n])
        start = max(float(m.result.states[0, 5]) for m in models.values())
        xs = sorted({start, *[r['x'] for r in data['groups']['concentration_pairs'] if r['trial'] == n]})
        for x in xs:
            sections.append(dict(trial=n, x_m=x, **{label: section_record(m, x) for label, m in models.items()}))
        for old in data['groups']['temperature_rows']:
            if old['trial'] != n:
                continue
            values = {}
            for label, model in models.items():
                state = model.state_at(old['x'])
                sy, sz = model.model.section_widths(state)
                predicted = model.temperature_at(old['x'], old['y'], old['z'])
                key = 'CONTROL_K' if label == 'pressure_loss_density' else 'matched_K'
                if abs(predicted-old[key]) > 1e-6:
                    raise ValueError('sensor-temperature replay changed')
                values[label] = dict(predicted_K=predicted, centre_T_K=model.model.centre_temperature(state),
                    centre_z_m=float(state[6]), sigma_y_m=sy, sigma_z_m=sz,
                    lateral_sigma=float(old['y']/sy), vertical_sigma=float((old['z']-state[6])/sz))
            sensor_geometry.append(dict(trial=n, channel=old['channel'], x_m=old['x'],
                observed_minimum_K=old['observed_minimum_K'], observed_median_K=old['observed_median_K'], **values))
        print('Verified flow comparison trial', n, flush=True)
    if {name: digest(ROOT/name) for name in hashes} != hashes:
        raise ValueError('input changed during audit')
    verify_inputs(directory)
    data.update(completed=True, created_utc=datetime.now(timezone.utc).isoformat(), sha256=hashes,
        replay_checks=replay_checks, sections=sections, sensor_geometry=sensor_geometry,
        promoted=False, interpretation='same research source method and profiles; only flow estimate changed')
    write_new(args.output, data)
    import json
    print(json.dumps(data['metrics'], indent=2), flush=True)
    print('Saved:', args.output, flush=True)


if __name__ == '__main__':
    main()
