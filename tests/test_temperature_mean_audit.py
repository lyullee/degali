import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location('temperature_mean_audit',Path(__file__).resolve().parents[1]/'tools/audit_temperature_mean.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_known_statistics_and_endpoint_weighting():
    r = module.statistics([0,0,9],[10,11,12],1)
    assert r['mean_K']==3
    assert r['median_K']==0
    assert r['trapezoid_time_mean_K']==2.25
    assert r['fraction_at_or_below_CONTROL']==pytest.approx(2/3)
    assert r['standard_deviation_sample_K']==pytest.approx(np.sqrt(27))
    assert r['moment_skewness']==pytest.approx(1/np.sqrt(2))


def test_constant_series():
    r=module.statistics([70]*4,[0,1,2,3],70)
    assert r['mean_K']==r['median_K']==r['p95_K']==70
    assert r['moment_skewness']==r['standard_deviation_sample_K']==0


@pytest.mark.parametrize('a,t', [([1,2],[0,2]),([1,np.nan],[0,1]),([1],[0]),([1,2],[0]),([1,2],[1,1])])
def test_invalid_series(a,t):
    with pytest.raises(ValueError):
        module.statistics(a,t,1)


def test_timestamp_and_invalid_type():
    assert module.seconds('13:25:28')==48328
    with pytest.raises(ValueError): module.seconds(3)
