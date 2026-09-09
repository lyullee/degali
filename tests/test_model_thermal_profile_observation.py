import math

import numpy as np
import pytest

from tools.audit_model_thermal_profile_observation import (
    aggregate,
    distribution,
    observed_profile_summary,
)


def test_distribution_preserves_unfitted_percentiles():
    result = distribution([1.0, 2.0, 3.0, 4.0])
    assert result["median"] == pytest.approx(2.5)
    assert result["p25"] == pytest.approx(1.75)
    assert result["p75"] == pytest.approx(3.25)


def test_observed_summary_applies_only_the_frozen_centre_gate():
    z = np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    shape = np.exp(-0.5 * (z / 0.2) ** 2)
    profiles = np.array([4.9 * shape, 5.0 * shape, 10.0 * shape])
    result = observed_profile_summary(z, profiles, centre_gate=5.0)
    assert result["usable_profiles"] == 2
    assert result["centre_amplitude"]["median"] == pytest.approx(7.5)
    assert result["variance"]["minimum"] == pytest.approx(
        result["variance"]["maximum"]
    )


def test_observed_summary_rejects_missing_signal():
    z = np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    with pytest.raises(ValueError, match="no observed profile"):
        observed_profile_summary(z, np.ones((3, 5)), centre_gate=5.0)
    invalid = np.ones((3, 5))
    invalid[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite and nonnegative"):
        observed_profile_summary(z, invalid, centre_gate=0.0)


def test_aggregate_keeps_amplitude_and_variance_as_separate_gates():
    def row(thermal_centre, thermal_variance, h2_centre, h2_variance):
        def entry(error, within):
            return {
                "absolute_log_ratio": error,
                "within_observed_iqr": within,
            }

        return {
            "comparison": {
                "thermal_deficit": {
                    "centre_amplitude": entry(thermal_centre, True),
                    "variance": entry(thermal_variance, False),
                },
                "hydrogen": {
                    "centre_amplitude": entry(h2_centre, False),
                    "variance": entry(h2_variance, True),
                },
            }
        }

    result = aggregate([
        row(math.log(2.0), math.log(1.1), math.log(1.2), math.log(1.3)),
        row(math.log(1.5), math.log(1.2), math.log(1.1), math.log(1.4)),
    ])
    assert result["thermal_deficit"]["centre_amplitude"]["within_observed_iqr"] == 2
    assert result["thermal_deficit"]["variance"]["within_observed_iqr"] == 0
    assert result["thermal_deficit"]["centre_amplitude"][
        "median_absolute_log_ratio"
    ] > result["thermal_deficit"]["variance"]["median_absolute_log_ratio"]
