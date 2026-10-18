import numpy as np
import pytest

from tools.audit_ti021_internal_promotion import _trace_evidence, _wall_diffusion_evidence
from tools.audit_ti021_internal_promotion import _vertical_observation_evidence


def test_legacy_observer_cannot_pass_on_unchanged_solver_mass_ledger():
    record={'vertical_observation_mode':'legacy_mass_fraction','diagnostics':{
        'vertical_observation_mode':'legacy_mass_fraction',
        'vertical_observation_maximum_relative_volume_residual':0.}}
    assert not _vertical_observation_evidence(record)['passed']
    record['vertical_observation_mode']='conservative_volume'
    record['diagnostics']['vertical_observation_mode']='conservative_volume'
    assert _vertical_observation_evidence(record)['passed']
    record['diagnostics']['vertical_observation_maximum_relative_volume_residual']=.01
    assert not _vertical_observation_evidence(record)['passed']


@pytest.mark.parametrize('end,reported,expected',[(600.,600.,True),(30.,30.,False),(600.,300.,False)])
def test_full_event_gate_reintegrates_traces_and_rejects_prefix_or_wrong_dose(tmp_path,end,reported,expected):
    path=tmp_path/'candidate.json'
    record={'duration_s':600.,'sensor_rows':[{'label':'TI021-S001',
        'observed_dose_vol_percent_s':600.,'predicted_dose_vol_percent_s':reported}]}
    np.savez(path.with_suffix('.npz'),time_s=np.array([0.,end]),
        concentration_vol_percent=np.ones((1,2)),
        predicted_dose_vol_percent_s=np.array([reported]),
        observed_dose_vol_percent_s=np.array([600.]))
    assert _trace_evidence(path,record)['passed'] is expected


def test_legacy_wall_blind_diffusion_cannot_be_promoted_on_mass_balance_alone():
    record={'wall_height_m':6.,'horizontal_diffusivity_m2_s':2.5}
    assert _wall_diffusion_evidence(record)['passed'] is False
    record['horizontal_diffusion_respect_walls']=True
    assert _wall_diffusion_evidence(record)['passed'] is True
    record['horizontal_diffusion_respect_walls']=False
    record['horizontal_diffusivity_m2_s']=0.
    assert _wall_diffusion_evidence(record)['passed'] is True
