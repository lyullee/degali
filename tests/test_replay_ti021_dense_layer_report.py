import inspect

import pytest

from tools.replay_ti021_dense_layer_report import (
    ALIASES,
    _parse_overrides,
    replay_parameters,
)
from tools.score_ti021_dense_layer_hybrid import run


def test_replay_parameter_map_covers_the_complete_run_signature():
    record={}
    expected={}
    for name,parameter in inspect.signature(run).parameters.items():
        if name=='output':
            continue
        key=ALIASES.get(name,name)
        value=(parameter.default if parameter.default is not inspect.Parameter.empty
               else f'value-for-{name}')
        record[key]=value
        expected[name]=value
    recovered=replay_parameters(record,component_diagnostics=True)
    expected['vertical_observation_component_diagnostics']=True
    assert recovered==expected


def test_replay_parameter_override_is_explicit_and_signature_checked():
    recovered = replay_parameters(
        {"horizontal_cell_size_m": 5.0},
        overrides={"h": 2.5, "wall_height_m": 5.0},
    )
    assert recovered["h"] == 2.5
    assert recovered["wall_height_m"] == 5.0
    with pytest.raises(ValueError, match="not a run"):
        replay_parameters({}, overrides={"sensor_gain": 2.0})


def test_parse_overrides_uses_json_types_and_rejects_duplicates():
    assert _parse_overrides(["wall_height_m=5.0", "height_resolved_wall_flux=true"]) == {
        "wall_height_m": 5.0,
        "height_resolved_wall_flux": True,
    }
    with pytest.raises(ValueError, match="duplicate"):
        _parse_overrides(["wall_height_m=5", "wall_height_m=6"])
