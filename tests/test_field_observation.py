import numpy as np
import pytest

from degali.addons.field_contracts import BoundedValue, SensorModel
from degali.addons.field_observation import (
    apply_sensor_model,
    h2_mole_fraction_from_mass_concentration,
    trailing_time_average,
)


def test_mass_concentration_conversion_and_sensor_lag_are_explicit():
    sensor = SensorModel(
        (1.0, 0.0, 0.5),
        response_time_s=BoundedValue(1.0),
        gain=BoundedValue(1.1),
        bias_mole_fraction=BoundedValue(0.01),
        averaging_time_s=1.0,
    )
    time = np.array([0.0, 1.0, 2.0, 3.0])
    concentration = np.array([0.0, 0.01, 0.01, 0.01])
    trace = apply_sensor_model(
        time, concentration, sensor, ambient_air_density_kg_m3=1.2
    )
    assert trace.true_mole_fraction[0] == pytest.approx(0.0)
    assert trace.true_mole_fraction[-1] < 1.0
    assert trace.response_mole_fraction[-1] < trace.true_mole_fraction[-1]
    assert trace.indicated_mole_fraction[-1] > trace.response_mole_fraction[-1]


def test_trailing_average_handles_nonuniform_times_without_resampling():
    averaged = trailing_time_average(
        [0.0, 1.0, 3.0], [0.0, 2.0, 2.0], window_s=2.0
    )
    assert averaged == pytest.approx([0.0, 1.0, 2.0])
    mole_fraction = h2_mole_fraction_from_mass_concentration(
        [0.0, 0.02], ambient_air_density_kg_m3=1.2
    )
    assert mole_fraction[0] == pytest.approx(0.0)
    assert 0.0 < mole_fraction[1] < 1.0
