from copy import deepcopy
import math

import pytest

from stage2_matched_source import TRIALS, make_source_document, require_checkpoint, finite_json


def documents():
    original = {'protocol': 'original', 'trials': [
        {'trial': n, 'pressure_loss_mass_flow_g_s': 200. + n,
         'coriolis_mass_flow_g_s': {'mean': 100. + n},
         'tc3_temperature_k': {'median': 20. + n / 100.},
         'pt2_barg': {'median': .1 * n}, 'hem_nozzle': {'old': n}, 'window': [10, 20]}
        for n in TRIALS]}
    reduced = {'trials': [{'trial': n, 'flow_mean_gs': 100. + n, 'window': [10, 20]} for n in TRIALS]}
    return original, reduced


def test_flow_ablation_retains_temperature_pressure_and_original_record():
    original, reduced = documents()
    before = deepcopy(original)
    out = make_source_document(original, reduced, 'coriolis')
    assert original == before
    for old, new in zip(original['trials'], out['trials']):
        assert new['pressure_loss_mass_flow_g_s'] == old['coriolis_mass_flow_g_s']['mean']
        assert new['original_pressure_loss_mass_flow_g_s'] == old['pressure_loss_mass_flow_g_s']
        assert new['tc3_temperature_k'] == old['tc3_temperature_k']
        assert new['pt2_barg'] == old['pt2_barg']
        assert 'hem_nozzle' not in new
        assert new['historical_hem_nozzle_not_used'] == old['hem_nozzle']


def test_pressure_control_retains_original_flow():
    original, reduced = documents()
    out = make_source_document(original, reduced, 'pressure_loss')
    assert [r['pressure_loss_mass_flow_g_s'] for r in out['trials']] == [
        r['pressure_loss_mass_flow_g_s'] for r in original['trials']]


@pytest.mark.parametrize('defect', ['duplicate', 'missing', 'changed_mean', 'nonfinite', 'window'])
def test_inconsistent_source_inputs_are_rejected(defect):
    original, reduced = documents()
    if defect == 'duplicate':
        original['trials'].append(deepcopy(original['trials'][0]))
    elif defect == 'missing':
        original['trials'].pop()
    elif defect == 'changed_mean':
        reduced['trials'][0]['flow_mean_gs'] += 1
    elif defect == 'window':
        reduced['trials'][0]['window'] = [11, 21]
    else:
        original['trials'][0]['pressure_loss_mass_flow_g_s'] = math.nan
    with pytest.raises(ValueError):
        make_source_document(original, reduced, 'coriolis')


def test_preserves_sealed_three_decimal_flow_and_reports_rounding():
    original, reduced = documents()
    original['trials'][0]['coriolis_mass_flow_g_s']['mean'] += .0004
    out = make_source_document(original, reduced, 'coriolis')
    assert out['trials'][0]['stage2_selected_mass_flow_g_s'] == 110.
    assert out['trials'][0]['stage2_coriolis_rounding_difference_g_s'] == pytest.approx(-.0004)


@pytest.mark.parametrize('change', [
    {'source_sha256': 'different'}, {'profile': 'density'},
    {'interfaces_accepted': False}, {'failures': {'10': 'bad'}},
    {'fit_velocity_spreading': True}, {'selected_trials': [10, 23]},
])
def test_field_cannot_reuse_unmatched_or_partial_interface_gate(change):
    config = {'source_sha256': 'input', 'flow_choice': 'coriolis'}
    checkpoint = dict(completed=True, interfaces_accepted=True, failures={},
        selected_trials=list(TRIALS), profile='enthalpy', fit_velocity_spreading=False, **config)
    require_checkpoint(checkpoint, config, 'enthalpy')
    checkpoint.update(change)
    with pytest.raises(ValueError):
        require_checkpoint(checkpoint, config, 'enthalpy')


def test_failed_diagnostic_nan_is_explicit_null_in_strict_json():
    assert finite_json({'r2': math.nan, 'ok': 1.}) == {'r2': None, 'ok': 1.}
