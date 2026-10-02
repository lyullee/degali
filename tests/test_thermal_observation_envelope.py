import pytest

from degali.validation.thermal_profile_moments import empirical_scalar_envelope


def test_empirical_scalar_envelope_has_half_tie_rank_and_sample_spread():
    result = empirical_scalar_envelope([1.0, 2.0, 2.0, 5.0], 2.0)
    assert result.count == 4
    assert result.prediction_midrank == pytest.approx(0.5)
    assert result.mean == pytest.approx(2.5)
    assert result.median == pytest.approx(2.0)
    assert result.sample_standard_deviation == pytest.approx(1.7320508075688772)
    assert result.within_p05_p95


@pytest.mark.parametrize(
    "samples,prediction",
    [([], 1.0), ([1.0], 1.0), ([1.0, float("nan")], 1.0), ([1.0, 2.0], float("nan"))],
)
def test_empirical_scalar_envelope_rejects_invalid_inputs(samples, prediction):
    with pytest.raises(ValueError, match="finite samples"):
        empirical_scalar_envelope(samples, prediction)
