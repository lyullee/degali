import numpy as np
import pytest

from degali.addons.fast_obstacle_transport import (
    FastObstacleReceptor,
    FastObstacleTransportConfig,
    solve_fast_obstacle_transport,
)
from degali.addons.semi_fv_obstacle import (
    RectangularObstacle2D,
    SourceRateSchedule,
)


def _config(*, obstacle=True, store_fields=False, time_step_s=0.005):
    return FastObstacleTransportConfig(
        length_m=12.0,
        height_m=6.0,
        nx=80,
        nz=40,
        duration_s=5.0,
        time_step_s=time_step_s,
        wind_speed_m_s=1.2,
        diffusivity_m2_s=0.05,
        source_schedule=SourceRateSchedule(
            (0.0, 3.0, 5.0),
            (0.002, 0.0, 0.0),
            source_id="finite-atmospheric-h2-per-metre-span",
        ),
        source_x_m=1.0,
        source_height_m=0.5,
        source_sigma_m=0.3,
        obstacle=(
            RectangularObstacle2D(5.0, 2.0, 2.0, label="block")
            if obstacle else None
        ),
        store_fields=store_fields,
    )


def test_open_field_streamfunction_recovers_uniform_divergence_free_wind():
    result = solve_fast_obstacle_transport(_config(obstacle=False))

    assert result.velocity_u_m_s == pytest.approx(1.2, abs=1.0e-11)
    assert result.velocity_w_m_s == pytest.approx(0.0, abs=1.0e-11)
    assert result.diagnostics["maximum_absolute_divergence_s_1"] < 1.0e-11
    assert not result.diagnostics["nonlocal_scalar_routing"]


def test_obstacle_transport_is_local_impermeable_and_mass_conservative():
    result = solve_fast_obstacle_transport(_config(store_fields=True))

    assert result.maximum_mass_residual_per_span_kg_m < 1.0e-12
    assert result.injected_mass_per_span_kg_m[-1] == pytest.approx(0.006)
    assert (
        result.domain_mass_per_span_kg_m[-1]
        + result.outflow_mass_per_span_kg_m[-1]
    ) == pytest.approx(0.006)
    assert result.diagnostics["maximum_absolute_divergence_s_1"] < 1.0e-11
    assert result.concentration_fields_kg_m3 is not None
    assert np.all(result.concentration_fields_kg_m3[:, result.solid_mask] == 0.0)
    assert np.all(np.isnan(result.concentration_kg_m3[result.solid_mask]))

    # Every face between solid and fluid is explicitly impermeable.
    solid = result.solid_mask
    assert np.all(result.u_face_m_s[:, 1:-1][solid[:, :-1] | solid[:, 1:]] == 0.0)
    assert np.all(result.w_face_m_s[1:-1, :][solid[:-1, :] | solid[1:, :]] == 0.0)


def test_windward_impact_is_not_teleported_to_the_first_roof_cell():
    receptors = (
        FastObstacleReceptor("windward", 3.8, 0.5),
        FastObstacleReceptor("roof", 5.0, 2.2),
        FastObstacleReceptor("lee", 7.0, 0.5),
    )
    result = solve_fast_obstacle_transport(_config(), receptors=receptors)
    peaks = {
        label: float(np.max(trace))
        for label, trace in result.receptor_traces_kg_m3.items()
    }

    assert peaks["windward"] > peaks["roof"] > peaks["lee"]
    maximum_j, maximum_i = np.unravel_index(
        np.nanargmax(result.concentration_kg_m3),
        result.concentration_kg_m3.shape,
    )
    # The local maximum remains on the windward side, rather than appearing
    # at the first cell over the obstacle roof as in legacy non-local routing.
    assert result.x_m[maximum_i] < 4.0
    assert result.z_m[maximum_j] < 1.0


def test_fast_obstacle_result_exposes_vol_percent_runtime_and_scope():
    result = solve_fast_obstacle_transport(_config())

    assert np.nanmax(result.concentration_vol_percent) > 0.0
    assert result.diagnostics["total_runtime_s"] > 0.0
    assert result.diagnostics["applicability"] == "conditional"
    assert result.diagnostics["quantitative_prediction_allowed"] is False
    assert "no no-slip separation" in result.diagnostics["warnings"][0]


def test_receptor_inside_solid_and_elevated_obstacle_are_rejected():
    with pytest.raises(ValueError, match="inside the obstacle"):
        solve_fast_obstacle_transport(
            _config(), receptors=(FastObstacleReceptor("solid", 5.0, 1.0),)
        )

    with pytest.raises(ValueError, match="ground-mounted"):
        FastObstacleTransportConfig(
            **{
                **_config().__dict__,
                "obstacle": RectangularObstacle2D(
                    5.0, 2.0, 1.0, base_height_m=1.0, label="elevated",
                ),
            }
        )


def test_explicit_stability_limit_is_enforced():
    with pytest.raises(ValueError, match="stability limit"):
        solve_fast_obstacle_transport(_config(time_step_s=0.2))
