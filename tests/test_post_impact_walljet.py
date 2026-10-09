import numpy as np
import pytest

from degali.addons.post_impact_walljet import (
    biased_radial_walljet_field,
    bifurcated_walljet_field,
    radial_walljet_field,
    radial_walljet_mixing,
    redistribute_vertical_profile,
    source_state_profile_weights,
)


def _field():
    x = np.linspace(-150.0, 150.0, 301)
    y = np.linspace(-150.0, 150.0, 301)
    xx, yy = np.meshgrid(x, y)
    field = 0.25 * np.exp(-((xx - 30.0) / 18.0) ** 2 - (yy / 8.0) ** 2)
    return field, x, y


def _integral(field, x, y):
    return float(np.trapezoid(np.trapezoid(field, x, axis=1), y, axis=0))


@pytest.mark.parametrize("operator", [radial_walljet_field, biased_radial_walljet_field, bifurcated_walljet_field])
def test_walljet_angular_operators_preserve_scalar_inventory(operator):
    field, x, y = _field()
    kwargs = {"x": x, "y": y, "half_angle_deg": 55.0}
    if operator is biased_radial_walljet_field:
        kwargs.update(bias_kappa=1.5, bias_center_deg=20.0)
    output = operator(field, **kwargs)

    assert np.all(np.isfinite(output))
    assert np.all((output >= 0.0) & (output <= 1.0))
    assert _integral(output, x, y) == pytest.approx(_integral(field, x, y), rel=1e-10)


def test_bifurcated_walljet_is_symmetric_without_vector_bias():
    field, x, y = _field()
    output = bifurcated_walljet_field(field, x=x, y=y, half_angle_deg=50.0)
    assert np.allclose(output, output[::-1, :], atol=2.0e-3)


def test_radial_walljet_mixing_preserves_each_downwind_column():
    field, x, y = _field()
    output = radial_walljet_mixing(
        field, x=x, y=y, height_m=1.0, wind_speed_m_s=5.0,
        jet_speed_m_s=20.0, start_m=5.0,
    )
    for index, station in enumerate(x):
        if station <= 5.0:
            continue
        assert np.trapezoid(output[:, index], y) == pytest.approx(
            np.trapezoid(field[:, index], y), rel=1e-10, abs=1e-12,
        )


def test_vertical_profile_weights_and_redistribution_are_conservative():
    heights = (0.1, 1.0, 2.0, 4.0)
    x = np.array([10.0, 20.0, 30.0])
    weights = source_state_profile_weights(
        heights, x, [0.0, 10.0, 30.0], [2.0, 3.0, 4.0], [0.2, 0.5, 1.0],
    )
    thickness = np.array([0.9, 0.95, 1.5, 2.0])
    assert np.allclose(np.sum(thickness[:, None] * weights, axis=0), 1.0)

    fields = {height: np.full((2, 3), index + 1.0) for index, height in enumerate(heights)}
    redistributed = redistribute_vertical_profile(
        fields, heights, x_m=x,
        source_x_m=[0.0, 10.0, 30.0],
        source_depth_m=[2.0, 3.0, 4.0],
        source_centre_m=[0.2, 0.5, 1.0],
    )
    before = sum(thickness[index] * fields[height] for index, height in enumerate(heights))
    after = sum(thickness[index] * redistributed[height] for index, height in enumerate(heights))
    assert np.allclose(before, after)
