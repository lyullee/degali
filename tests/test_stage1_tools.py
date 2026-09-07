import importlib.util
from pathlib import Path
import sys

import pytest

TOOLS = Path(__file__).resolve().parents[1]/'tools'
sys.path.insert(0, str(TOOLS))
from audit_stage1_assessment import numerical_match, promotion_gates, BASE_COLUMNS
import stage1


def test_zero_reproduction_and_explicit_mismatch():
    assert numerical_match([[1, 2]], [[1, 2]])['maximum_scaled_error'] == 0
    with pytest.raises(ValueError):
        numerical_match([1], [1.01])
    with pytest.raises(ValueError):
        numerical_match([1], [[1]])
    with pytest.raises(ValueError):
        numerical_match([float('nan')], [float('nan')])


def test_joint_promotion_never_hides_adverse_center():
    base = dict(MG=1.1, VG=1.2, FAC2=.9, geometry=dict(width_ratio=1.1, centre_mae=.05))
    candidate = dict(MG=1.01, VG=1.1, FAC2=.95, geometry=dict(width_ratio=1.01, centre_mae=.06))
    gates = promotion_gates(base, candidate, True)
    assert not gates['centre_MAE_not_worse']
    assert not all(gates.values())
    assert not promotion_gates(base, candidate, False)['all_interfaces_accepted']


def test_manifest_checks_mutation_missing_and_escape(tmp_path):
    path = tmp_path/'input.txt'
    path.write_text('original', encoding='utf-8')
    record = dict(completed=True, sha256={'input.txt': stage1.digest(path)})
    assert stage1.verify_hashes(record, tmp_path) == 1
    path.write_text('changed', encoding='utf-8')
    with pytest.raises(ValueError):
        stage1.verify_hashes(record, tmp_path)
    with pytest.raises(ValueError):
        stage1.verify_hashes(dict(completed=True, sha256={'missing': 'x'}), tmp_path)
    with pytest.raises(ValueError):
        stage1.checked_path(tmp_path, '../outside.txt')
    with pytest.raises(ValueError):
        stage1.verify_hashes(dict(completed=False, sha256={}), tmp_path)


def test_write_new_never_overwrites(tmp_path):
    path = tmp_path/'result.json'
    stage1.write_new(path, {'complete': True})
    with pytest.raises(FileExistsError):
        stage1.write_new(path, {'complete': False})
    assert stage1.read(path) == {'complete': True}


def test_complete_reported_schema_is_twelve_columns():
    assert len(BASE_COLUMNS) == 12
    assert BASE_COLUMNS[2] == 'reported_imaged_C_H2_kg_m3'
    assert BASE_COLUMNS[-1] == 'arc_length_m'
