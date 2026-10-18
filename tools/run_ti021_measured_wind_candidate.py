"""Evaluate digitized TI021 meteorology against the frozen 600 s baseline.

The wind history is independent of gas receptor values.  It is applied only
after the integral-to-layer handoff, and promotion requires an overall MAPE
improvement without worsening any populated height or downwind-distance group.
"""
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
from tools.score_ti021_dense_layer_hybrid import run


def _mape(predicted, observed, mask):
    selected = mask & (observed > 0.0)
    if not np.any(selected):
        return None, 0
    value = np.mean(np.abs(predicted[selected] / observed[selected] - 1.0)) * 100.0
    return float(value), int(np.count_nonzero(selected))


def _group_metrics(predicted, observed, sensors):
    heights = np.array([sensor.height_m for sensor in sensors], dtype=float)
    downwind = np.array([sensor.local_y_m for sensor in sensors], dtype=float)
    groups = {}
    for height in sorted(set(heights)):
        value, count = _mape(predicted, observed, np.isclose(heights, height))
        if count:
            groups[f'height_{height:g}_m'] = {'mape_percent': value, 'positive_count': count}
    for name, lower, upper in (
            ('downwind_lt_150_m', -np.inf, 150.0),
            ('downwind_150_to_300_m', 150.0, 300.0),
            ('downwind_ge_300_m', 300.0, np.inf)):
        mask = (downwind >= lower) & (downwind < upper)
        value, count = _mape(predicted, observed, mask)
        if count:
            groups[name] = {'mape_percent': value, 'positive_count': count}
    return groups


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wind-history', type=Path, required=True)
    parser.add_argument('--component-mode', choices=['measured', 'direction_only', 'speed_only'],
                        default='measured')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--advection-reference-x-m', type=float, default=None)
    parser.add_argument('--advection-speed-m-s', type=float, default=None)
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
        wind_history_path=args.wind_history.resolve(),
        wind_history_component_mode=args.component_mode,
        wind_history_advection_reference_x_m=args.advection_reference_x_m,
        wind_history_advection_speed_m_s=args.advection_speed_m_s,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    config_path = args.output.with_suffix('.configuration.json')
    config_path.write_text(json.dumps({
        key: str(value) if isinstance(value, Path) else value
        for key, value in params.items()
    }, indent=2), encoding='utf-8')

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
    workbook = ROOT / 'ti21_input_ascii.xls'
    dataset = read_thorney_ti021_workbook(workbook)
    baseline_groups = _group_metrics(old, observed, dataset.dose_sensors)
    candidate_groups = _group_metrics(predicted, observed, dataset.dose_sensors)
    observation_modes_match = (
        baseline.get('observation_species_mode', 'r12_component')
        == result.get('observation_species_mode', 'r12_component')
    )
    group_comparison = {}
    if observation_modes_match:
        for key in baseline_groups:
            old_value = baseline_groups[key]['mape_percent']
            new_value = candidate_groups[key]['mape_percent']
            group_comparison[key] = {
                'baseline_mape_percent': old_value,
                'candidate_mape_percent': new_value,
                'delta_percentage_points': new_value - old_value,
                'nondegraded': new_value <= old_value + 1e-12,
                'positive_count': candidate_groups[key]['positive_count'],
            }
    positive = observed > 0.0
    baseline_mape = float(np.mean(np.abs(old[positive] / observed[positive] - 1.0)) * 100.0)
    candidate_mape = float(np.mean(
        np.abs(predicted[positive] / observed[positive] - 1.0)) * 100.0)
    conservation = result['diagnostics']
    conservation_ok = bool(
        abs(conservation.get('maximum_relative_species_mass_residual', np.inf)) <= 1e-8
        and abs(conservation.get('mixture_mass_residual_kg', np.inf)) <= 1e-6
        and conservation.get('momentum_conservation_gate') is True)
    promotion_allowed = bool(
        labels_match
        and observation_modes_match
        and candidate_mape < baseline_mape
        and all(item['nondegraded'] for item in group_comparison.values())
        and conservation_ok)
    result['paired_baseline'] = {
        'source': str(args.baseline.resolve()),
        'identical_sensor_order': labels_match,
        'observation_species_modes_match': observation_modes_match,
        'baseline_all_positive_mape_percent': (
            baseline_mape if observation_modes_match else None),
        'candidate_all_positive_mape_percent': candidate_mape,
        'delta_percentage_points': (
            candidate_mape - baseline_mape if observation_modes_match else None),
        'configuration_path': str(config_path.resolve()),
        'group_comparison': group_comparison,
        'conservation_gate_passed': conservation_ok,
        'promotion_allowed': promotion_allowed,
        'target_10_percent_reached': candidate_mape <= 10.0,
    }
    result['promotion_allowed'] = promotion_allowed
    result['candidate_basis'] = (
        'digitized independent 3-minute-running-mean meteorology; no gas receptor values used')
    result['runtime_with_audit_s'] = time.perf_counter() - started
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result['paired_baseline'], indent=2))


if __name__ == '__main__':
    main()
