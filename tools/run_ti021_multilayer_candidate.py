"""Run a conservative fixed-layer TI021 candidate against the frozen baseline."""
import argparse
import inspect
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
from tools.run_ti021_measured_wind_candidate import _group_metrics
from tools.score_ti021_dense_layer_hybrid import run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--layers', type=int, default=13)
    parser.add_argument('--vertical-diffusivity-m2-s', type=float, default=.28)
    parser.add_argument(
        '--multilayer-entrainment-placement',
        choices=('top', 'legacy_ground'), default='top')
    parser.add_argument('--height-resolved-wall-flux', action='store_true')
    parser.add_argument('--initial-vertical-gaussian', action='store_true')
    parser.add_argument('--obstacle-vertical-lift-fraction', type=float, default=0.0)
    parser.add_argument('--obstacle-vertical-lift-decay-m', type=float, default=30.0)
    parser.add_argument('--wall-impulse-vertical-lift-fraction', type=float, default=0.0)
    parser.add_argument('--dense-settling-velocity-m-s', type=float, default=0.0)
    parser.add_argument('--baseline', type=Path, default=(
        ROOT / 'outputs/ti021-baseline-after-align-default-2026-10-16.json'))
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text(encoding='utf-8'))
    accepted = inspect.signature(run).parameters
    params = {key: value for key, value in baseline.items() if key in accepted}
    params.update(
        h=baseline['horizontal_cell_size_m'],
        dt=baseline['max_dt_s'],
        source_dt=baseline['source_time_step_s'],
        entrainment=baseline['dense_entrainment_model'],
        vertical_observation_profile='multilayer',
        vertical_layers=args.layers,
        vertical_mixing_diffusivity_m2_s=args.vertical_diffusivity_m2_s,
        multilayer_entrainment_placement=args.multilayer_entrainment_placement,
        height_resolved_wall_flux=args.height_resolved_wall_flux,
        initial_vertical_gaussian=args.initial_vertical_gaussian,
        obstacle_vertical_lift_fraction=args.obstacle_vertical_lift_fraction,
        obstacle_vertical_lift_decay_m=args.obstacle_vertical_lift_decay_m,
        wall_impulse_vertical_lift_fraction=args.wall_impulse_vertical_lift_fraction,
        dense_settling_velocity_m_s=args.dense_settling_velocity_m_s,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    config_path = args.output.with_suffix('.configuration.json')
    config_path.write_text(json.dumps(params, indent=2), encoding='utf-8')

    started = time.perf_counter()
    result = run(args.output, **params)
    predicted = np.array([
        row['predicted_dose_vol_percent_s'] for row in result['sensor_rows']])
    old = np.array([
        row['predicted_dose_vol_percent_s'] for row in baseline['sensor_rows']])
    observed = np.array([
        row['observed_dose_vol_percent_s'] for row in baseline['sensor_rows']])
    labels_match = ([row['label'] for row in result['sensor_rows']]
                    == [row['label'] for row in baseline['sensor_rows']])
    dataset = read_thorney_ti021_workbook(ROOT / 'ti21_input_ascii.xls')
    baseline_groups = _group_metrics(old, observed, dataset.dose_sensors)
    candidate_groups = _group_metrics(predicted, observed, dataset.dose_sensors)
    group_comparison = {}
    for key in baseline_groups:
        before = baseline_groups[key]['mape_percent']
        after = candidate_groups[key]['mape_percent']
        group_comparison[key] = {
            'baseline_mape_percent': before,
            'candidate_mape_percent': after,
            'delta_percentage_points': after-before,
            'nondegraded': after <= before+1.e-12,
            'positive_count': candidate_groups[key]['positive_count'],
        }
    positive = observed > 0.0
    baseline_mape = float(np.mean(np.abs(old[positive]/observed[positive]-1.0))*100.0)
    candidate_mape = float(np.mean(
        np.abs(predicted[positive]/observed[positive]-1.0))*100.0)
    diagnostics = result['diagnostics']
    conservation_ok = bool(
        abs(diagnostics.get('maximum_relative_species_mass_residual', np.inf)) <= 1.e-8
        and abs(diagnostics.get('mixture_mass_residual_kg', np.inf)) <= 1.e-6
        and diagnostics.get('momentum_conservation_gate') is True)
    promotion_allowed = bool(
        labels_match and candidate_mape < baseline_mape
        and all(item['nondegraded'] for item in group_comparison.values())
        and conservation_ok)
    result['paired_baseline'] = {
        'source': str(args.baseline.resolve()),
        'identical_sensor_order': labels_match,
        'baseline_all_positive_mape_percent': baseline_mape,
        'candidate_all_positive_mape_percent': candidate_mape,
        'delta_percentage_points': candidate_mape-baseline_mape,
        'configuration_path': str(config_path.resolve()),
        'group_comparison': group_comparison,
        'conservation_gate_passed': conservation_ok,
        'promotion_allowed': promotion_allowed,
        'target_10_percent_reached': candidate_mape <= 10.0,
    }
    result['promotion_allowed'] = promotion_allowed
    result['candidate_basis'] = (
        'fixed vertical layers with conservative diffusion and EOS volume-occupancy observation')
    result['runtime_with_audit_s'] = time.perf_counter()-started
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result['paired_baseline'], indent=2))


if __name__ == '__main__':
    main()
