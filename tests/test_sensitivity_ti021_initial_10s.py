import numpy as np
import pytest
from tools.sensitivity_ti021_initial_10s import (
    censored_timing_score, bounded_arrival_window_score,
    relative_peak_timing_score, sensor_average, make_cases,
)


def test_late_observation_does_not_inflate_censored_earliness():
    t=np.array([0.,5.,10.]);c=np.array([[0.,1.,1.],[0.,0.,0.]])
    score=censored_timing_score(t,c,np.array([100.,12.]))
    assert score['early_detected_sensor_count']==1
    assert score['sum_censored_earliness_s']==5.
    assert score['mean_censored_earliness_s_over_known_arrivals']==2.5
    assert score['first_detection_time_s']==[5.,None]
    assert score['not_a_full_arrival_mae'] is True


def test_no_observed_arrival_is_not_invented_as_zero_concentration():
    s=censored_timing_score([0.,10.],np.ones((1,2)),np.array([np.nan]))
    assert s['early_detected_sensor_count']==0
    assert s['known_arrival_sensor_count']==0


def test_bounded_window_score_penalizes_missing_in_window_arrivals_as_lower_bound():
    score=bounded_arrival_window_score(
        [0.,5.,10.],np.array([[0.,1.,1.],[0.,0.,0.],[0.,0.,1.]]),
        np.array([5.,8.,12.]),threshold=.1,endpoint=10.)
    assert score['observed_arrivals_in_window_count']==2
    assert score['predicted_arrivals_in_window_count']==1
    assert score['missing_by_horizon_count']==1
    assert score['absolute_timing_error_lower_bound_s']==pytest.approx(2.)
    assert score['mean_absolute_timing_error_lower_bound_s']==pytest.approx(1.)
    assert score['not_a_full_arrival_mae'] is True


def test_relative_peak_operator_uses_each_model_trace_peak():
    score=relative_peak_timing_score(
        [0., 1., 2., 3.],
        np.array([[0., .021, .4, .4], [0., .2, .1, .1]]),
        np.array([10., 10.]),
    )
    assert score['first_detection_time_s']==[1., 1.]
    assert score['model_local_peak_vol_percent']==[.4, .2]
    assert score['measured_local_peak_available'] is False


def test_sensor_average_uses_no_future_values_and_zero_preevent():
    t=np.linspace(0,1,101);c=np.ones((1,len(t)))
    v=sensor_average(t,c)
    assert v[0,0]==0
    assert v[0,30]==pytest.approx(.5)
    assert v[0,60]==pytest.approx(1.)
    assert v.max()<=1.000000001


def test_cases_have_unique_names_and_off_cases_are_explicit():
    cases=make_cases()
    assert len({c[0] for c in cases})==len(cases)
    assert any(c[0]=='gravity_off' and c[2]['gravity_spreading_coefficient']==0 for c in cases)
