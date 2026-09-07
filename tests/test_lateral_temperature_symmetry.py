import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('lateral_symmetry',Path(__file__).resolve().parents[1]/'tools/audit_lateral_temperature_symmetry.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_triangle_lower_bound_and_achieving_interval():
    p=np.array([-1,0,1,2,3,4])
    r=module.symmetric_errors(0,3,p)
    assert np.all(r['pair_absolute_error_sum']>=r['pair_absolute_error_sum_floor'])
    assert r['pair_MAE_excess'][1:5]==pytest.approx([0,0,0,0])
    assert r['pair_MAE_floor']==pytest.approx(np.full(6,1.5))
    assert r['pair_squared_error_sum_floor']==pytest.approx(np.full(6,4.5))


def test_swap_and_temperature_reference_invariant():
    a=module.symmetric_errors(60,90,80)
    b=module.symmetric_errors(90+273.15,60+273.15,80+273.15)
    for k in a: assert a[k]==pytest.approx(b[k])


def test_cardinal_wind_from_and_rotation():
    assert module.wind_components(2,255)==pytest.approx([2,0],abs=1e-14)
    assert module.wind_components(2,345)==pytest.approx([0,2],abs=1e-14)
    assert module.wind_components(2,75)==pytest.approx([-2,0],abs=1e-14)
    assert module.wind_components(2,338)==pytest.approx(module.wind_components(2,338+37,75+37))


def test_vector_speed_identity_and_wraparound():
    uv=module.wind_components([1,2,3],[359,1,721])
    assert np.sum(uv**2,axis=1)==pytest.approx([1,4,9])
    assert module.wind_components(2,1)==pytest.approx(module.wind_components(2,361))


@pytest.mark.parametrize('speed,angle',[(-1,20),(1,np.nan),(np.inf,20)])
def test_bad_wind(speed,angle):
    with pytest.raises(ValueError): module.wind_components(speed,angle)


def test_nonfinite_pair_rejected():
    with pytest.raises(ValueError): module.symmetric_errors(1,2,np.nan)
