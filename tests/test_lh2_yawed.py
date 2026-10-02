import math

import pytest

from degali.addons.axisymmetric_jet import AxisymmetricJetSource
from degali.addons.transient_receptor import WindHistory
from degali.lh2 import (
    ApplicabilityError,
    run_lh2_finite_release_research,
    run_lh2_yawed_crosswind_research,
)


def _source():
    return AxisymmetricJetSource(
        diameter=0.001, velocity=500.0, density=0.5,
        temperature=45.0, theta=0.0, y=0.5,
    )


def test_yawed_path_rejects_reverse_axial_branch():
    with pytest.raises(ValueError, match="reverse axial yaw"):
        run_lh2_yawed_crosswind_research(
            _source(), wind=2.5, wind_angle=0.0, release_angle=math.pi,
        )


@pytest.mark.slow
def test_yawed_path_is_research_only_and_returns_global_trajectory():
    result = run_lh2_yawed_crosswind_research(
        _source(), wind=2.5, wind_angle=0.0, release_angle=0.5,
        ambient_temperature=295.0, relative_humidity=0.0,
        maximum_nearfield_distance=0.02, radial_points=41,
        nearfield_maximum_step=0.001,
        nearfield_relative_tolerance=2.0e-6,
        maximum_distance=0.1, maximum_step=0.02,
    )
    assert result.trajectory.states.shape[1] == 10
    assert result.validated is False
    assert result.trajectory.states[-1, 7] > result.trajectory.states[0, 7]


@pytest.mark.slow
def test_finite_release_runs_jet_plume_transition_puff_and_receptor_trace():
    wind_history = WindHistory(
        time_s=[0.0, 0.002, 0.102], speed_m_s=[2.5, 2.5, 2.5],
        direction_from_deg=[270.0, 180.0, 180.0],
    )
    result = run_lh2_finite_release_research(
        _source(), source_duration_s=0.002, puff_duration_s=0.1,
        wind=2.5, ambient_temperature=295.0, relative_humidity=0.0,
        maximum_nearfield_distance=0.02, radial_points=41,
        nearfield_maximum_step=0.001,
        nearfield_relative_tolerance=2.0e-6,
        crosswind_maximum_distance=2.0,
        crosswind_maximum_step=0.02, puff_time_step=0.01,
        puff_wind_history=wind_history,
    )
    assert result.handoff.status == "transition_ready"
    assert result.handoff.transition_arc_length_m > 0.0
    assert result.puff.maximum_relative_mass_residual < 1.0e-10
    assert result.puff.maximum_relative_hydrogen_residual < 1.0e-10
    assert result.puff.states[-1].centre_position_m[0] \
        > result.transition_position_m[0]
    assert result.puff.states[-1].centre_position_m[1] > 0.0
    trace = result.receptor_trace(result.transition_position_m)
    assert trace.time_s[0] == pytest.approx(0.002)
    assert 0.0 < trace.peak_mole_fraction <= 1.0
    assert (
        "pretransition_wind_direction_span_exceeds_steady_limit"
        in result.steady.applicability_failures
    )


def test_finite_release_strict_scope_rejects_variable_pretransition_wind_early():
    wind_history = WindHistory(
        time_s=[0.0, 0.001, 0.2],
        speed_m_s=[2.0, 4.0, 2.5],
        direction_from_deg=[270.0, 180.0, 180.0],
    )
    with pytest.raises(ApplicabilityError, match="pretransition_wind"):
        run_lh2_finite_release_research(
            _source(), source_duration_s=0.002, puff_duration_s=0.1,
            wind=2.5, puff_wind_history=wind_history, strict_scope=True,
        )


@pytest.mark.slow
def test_five_second_field_scale_release_reaches_conservative_puff():
    source = AxisymmetricJetSource(
        diameter=0.0254, velocity=1000.0, density=0.5,
        temperature=45.0, theta=0.0, y=1.5,
    )
    result = run_lh2_finite_release_research(
        source, source_duration_s=5.0, puff_duration_s=0.2, wind=2.5,
        maximum_nearfield_distance=3.0, minimum_mass_fraction=1.0e-5,
        radial_points=41, nearfield_maximum_step=0.01,
        nearfield_relative_tolerance=2.0e-6,
        crosswind_maximum_distance=200.0,
        crosswind_maximum_step=0.1, puff_time_step=0.02,
    )
    assert 10.0 < result.handoff.transition_arc_length_m < 40.0
    assert result.handoff.hydrogen_mass_kg == pytest.approx(
        source.fuel_mass_flow * 5.0
    )
    assert result.puff.maximum_relative_mass_residual < 1.0e-10
    assert result.puff.maximum_relative_hydrogen_residual < 1.0e-10
    assert result.accepted
